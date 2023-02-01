# @created: 25 Nov 21
# @modified: 14 May 22
# @authors: Yoann Pradat
#
# Run a specified survival model on a specified completed table of a specified set of features.
# The fitting of the model hyperparameters is performed once a first cross-validation and then the best hyperparameters
# are used during repeated cross validation. The "survival model" may include a first step of feature selection in which
# case the robustness of this selection will be assessed through the different cross-validation repeats.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(yaml))

source("workflow/functions/metrics.R")
source("workflow/functions/feature_selection.R")
source("workflow/functions/survival_models.R")

# functions ============================================================================================================

load_indices <- function(filepath, repeat_mode){
  df_indices <- load_table(filepath)
  full_foldid <- as.numeric(df_indices %>% filter(Run=="Full") %>% arrange(Index) %>% pull(var="Split"))
  runs_not_full <- setdiff(unique(df_indices$Run), "Full")

  train_indices_all <- list()
  test_indices_all <- list()
  train_foldid_all <- list()

  for (run in runs_not_full){
    df_indices_run <- df_indices %>% filter(Run==run)

    if (repeat_mode=="boot"){
      train_indices <- df_indices_run %>% filter(Split=="Train") %>% pull(var="Index")
      test_indices <- df_indices_run %>% filter(Split=="Test") %>% pull(var="Index")

      train_indices_all[[paste0(run)]] <- train_indices
      test_indices_all[[paste0(run)]] <- test_indices
      train_foldid_all[[paste0(run)]] <- NULL
    } else if (repeat_mode=="xval"){
      for (split in unique(df_indices_run$Split)){
        train_indices <- df_indices_run %>% filter(Split!=split) %>% pull(var="Index")
        test_indices <- df_indices_run %>% filter(Split==split) %>% pull(var="Index")

        train_indices_all[[paste0(run,"/", split)]] <- train_indices
        test_indices_all[[paste0(run,"/", split)]] <- test_indices

        col_split_foldid <- paste0("Split_Nested_", split)
        train_foldid <- as.numeric(df_indices_run %>% filter(Split!=split) %>% pull(var=col_split_foldid))
        train_foldid_all[[paste0(run,"/", split)]] <- train_foldid 
      }
    }
  }

  list(train=train_indices_all, test=test_indices_all, train_foldid=train_foldid_all, full_foldid=full_foldid)
}


run_selection_and_model <- function(surv, df, name_selection, name_model, train_indices=NULL,
                                    train_foldid=NULL){

  warn_selection <- NULL
  df_fs <- withCallingHandlers(
    feature_selection(surv=surv, df=df, name_selection=name_selection, name_model=name_model,
                       train_indices=train_indices, train_foldid=train_foldid)
    , warning=function(w) {
      warn_selection <<- append(warn_selection, conditionMessage(w))
      invokeRestart("muffleWarning")
  })

  covs_all <- colnames(df)
  covs_in <- colnames(df_fs)
  covs_out <- setdiff(covs_all, covs_in) 

  out <- train_survival(surv=surv, df=df_fs, name_model=name_model, train_indices=train_indices,
                        train_foldid=train_foldid)

  # add info about covs in and out
  df_cov <- out$df_cov
  df_cov <- df_cov %>% mutate(Selected=1)

  if (length(covs_out)>0){
    df_cov_out <- tibble(Covariate=covs_out, Coefficient=0, Selected=0)
    df_cov <- bind_rows(df_cov, df_cov_out)
  }

  # add info about warning feature selection
  out$warn_selection <- warn_selection

  out$df_cov <- df_cov

  out
}


bind_results <- function(x,y,...){
  df_met <- bind_rows(x$df_met, y$df_met)
  df_cov <- bind_rows(x$df_cov, y$df_cov)
  list(df_met=df_met, df_cov=df_cov)
}


setup_parallel_computations <- function(n_cores, require_quality_metrics=F){
  cl <- makePSOCKcluster(n_cores)
  registerDoParallel(cl)
  tmp <- clusterEvalQ(cl, library(doParallel))
  tmp <- clusterEvalQ(cl, library(survival))
  tmp <- clusterEvalQ(cl, library(dplyr))
  tmp <- clusterEvalQ(cl, library(stringr))
  tmp <- clusterEvalQ(cl, library(tibble))
  clusterExport(cl, c("run_selection_and_model"))
  clusterExport(cl, c("feature_selection", "feature_selection_fdr", "feature_selection_lasso"))
  clusterExport(cl, c("train_survival", "train_coxph_formula", "train_coxph_ridge", "train_coxph_glmnet"))
  clusterExport(cl, c("cv_coxph_ridge", "get_formula_coxph_ridge"))

  if (require_quality_metrics){
    clusterExport(cl, c("cph_survest", "base_survest", "c_index_ipcw", "brier_score_ipcw"))
    clusterExport(cl, c("assess_cox_model_quality", "assess_cox_model_quality_train_test"))
        
    # Giving access to C++-coded function to all nodes is a bit special
    tmp <- clusterEvalQ(cl, library(Rcpp))
    clusterExport(cl, c("code_compute_concordance", "code_compute_brier_score"))
    tmp <- clusterEvalQ(cl, cppFunction(code_compute_concordance))
    tmp <- clusterEvalQ(cl, cppFunction(code_compute_brier_score))
  }

  cl
}


train_survival_jackknife <- function(df, surv, name_selection, name_model, n_cores, train_foldid=NULL, verbose=T){
  # generate jackknife train indices
  n_jackknife <- nrow(df)
  train_indices_all <- lapply(1:n_jackknife, function(j) setdiff(1:nrow(df), j))

  # prepare cluster
  if (verbose) cat(paste("-exporting functions to cluster..."))
  cl <- setup_parallel_computations(n_cores)
  if (verbose) cat("done!\n")

  if (verbose) cat(paste("-running jackknife iterations on cluster using", n_cores, "cores ..."))

  results <- foreach(j=1:n_jackknife, 
                     .combine = bind_rows, 
                     .inorder = F,
                     .multicombine = F, 
                     .noexport = c())  %dopar% {

    if (!is.null(train_foldid)){
      train_foldid_j <- train_foldid[-j]
    } else {
      train_foldid_j <- NULL
    }

    out <- run_selection_and_model(surv=surv,
                                   df=df,
                                   name_selection=name_selection,
                                   name_model=name_model,
                                   train_indices=train_indices_all[[j]],
                                   train_foldid=train_foldid_j)

    warn <- out$warn
    warn_selection <- out$warn_selection
    err <- out$err
    if (!is.null(warn)) warn <- paste(warn, collapse=";")
    if (is.null(warn)) warn <- NA
    if (!is.null(warn_selection)) warn_selection <- paste(warn_selection, collapse=";")
    if (is.null(warn_selection)) warn_selection <- NA
    if (is.null(err)) err <- NA
    df_cov <- out$df_cov %>% mutate(Run=paste0("Jack", j), Selection=name_selection, Model=name_model,
                                    Warning=warn, Error=err, Warning_Selection=warn_selection)

    df_cov
  }

  stopCluster(cl)
  if (verbose) cat("done!\n")

  results
}


train_test_survival_repeat <- function(surv, df, name_selection, name_model, taus, ipcw, indices, n_cores, verbose=T){
  train_indices_all <- indices$train
  test_indices_all <- indices$test
  train_foldid_all <- indices$train_foldid
  run_repeats <- names(train_indices_all)
  n_repeats <- length(train_indices_all)

  # prepare cluster
  if (verbose) cat(paste("-exporting functions to cluster..."))
  cl <- setup_parallel_computations(n_cores, require_quality_metrics=T)
  if (verbose) cat("done!\n")

  # parallel for loop
  if (verbose) cat(paste("-running repeated iterations on cluster using", n_cores, "cores ..."))

  results <- foreach(b = 1:n_repeats, 
                     .combine = 'bind_results', 
                     .inorder = F,
                     .multicombine = F, 
                     .noexport = c("compute_concordance", "compute_brier_score"))  %dopar% {
  
    run <- run_repeats[b]
    out <- run_selection_and_model(surv=surv, 
                                   df=df, 
                                   name_selection=name_selection, 
                                   name_model=name_model,
                                   train_indices=train_indices_all[[run]],
                                   train_foldid=train_foldid_all[[run]])


    out$df_met <- assess_cox_model_quality_train_test(surv=surv, df=df, linear_predictors=out$linear_predictors,
                                                      taus=taus, ipcw=ipcw, train_indices=train_indices_all[[run]],
                                                      test_indices=test_indices_all[[run]])
    warn <- out$warn
    warn_selection <- out$warn_selection
    err <- out$err
    if (!is.null(warn)) warn <- paste(warn, collapse=";")
    if (is.null(warn)) warn <- NA
    if (!is.null(warn_selection)) warn_selection <- paste(warn_selection, collapse=";")
    if (is.null(warn_selection)) warn_selection <- NA
    if (is.null(err)) err <- NA
    df_met <- out$df_met %>% mutate(Run=run, Selection=name_selection, Model=name_model, Warning=warn, Error=err,
                                    Warning_Selection=warn_selection)
    df_cov <- out$df_cov %>% mutate(Run=run, Selection=name_selection, Model=name_model, Warning=warn, Error=err,
                                    Warning_Selection=warn_selection)

    list(df_met=df_met, df_cov=df_cov)
  }

  stopCluster(cl)
  if (verbose) cat("done!\n")

  results
}


main <- function(args){
  # load data and covariates tables
  df_dat <- load_table(args$input_dat)
  df_cov <- load_table(args$input_cov)
  cov_codes <- df_cov$Code

  # load indices
  indices <- load_indices(args$indices, args$repeat_mode)
  
  # checks
  stopifnot(sum(is.na(df_dat))==0)

  # split target and covariates
  surv <- Surv(df_dat$Survival_Time, df_dat$Survival_Status)
  df <- df_dat[,cov_codes]

  # time horizons for restricted c-index computations
  tau_min <- 30
  tau_max <- ceil(quantile(df_dat$Survival_Time, 0.95)/30)*30
  taus <- seq(from=tau_min, to=tau_max, by=30)

  # compute ipcw weights
  # is it correct to use the full cohort?
  if (args$verbose) cat(paste("-computing ipcw weights using all data ...\n"))

  # fix the folds for cv.glmnet
  ipcw <- ipcw_cox(surv, df, taus, foldid=indices$full_foldid)

  # train/test on bootstrapped samples of the data
  results_repeat <- train_test_survival_repeat(surv=surv, df=df, name_selection=args$name_selection,
                                               name_model=args$name_model, taus=taus, ipcw=ipcw, 
                                               indices=indices, n_cores=args$n_cores, verbose=args$verbose)

  # train on the full data to compare standard errors estimated by bootstrap to standard errors estimated from
  # standard statistical theory
  results_full <- run_selection_and_model(surv=surv, df=df, name_selection=args$name_selection,
                                          name_model=args$name_model, train_indices=NULL,
                                          train_foldid=indices$full_foldid)

  results_full$df_met <- assess_cox_model_quality(surv=surv, df=df, linear_predictors=results_full$linear_predictors,
                                                  taus=taus, ipcw=ipcw)

  warn <- results_full$warn
  warn_selection <- results_full$warn_selection
  err <- results_full$err
  if (!is.null(warn)) warn <- paste(warn, collapse=";")
  if (is.null(warn)) warn <- NA
  if (!is.null(warn_selection)) warn_selection <- paste(warn_selection, collapse=";")
  if (is.null(warn_selection)) warn_selection <- NA
  if (is.null(err)) err <- NA

  results_full$df_met <- results_full$df_met %>%
    mutate(Split="Train", Run="Full", Selection=args$name_selection, Model=args$name_model, Warning=warn, Error=err,
           Warning_Selection=warn_selection)

  results_full$df_cov <- results_full$df_cov %>% 
    mutate(Run="Full", Selection=args$name_selection, Model=args$name_model, Warning=warn, Error=err,
           Warning_Selection=warn_selection)

  # train on jackknife samples of the data
  if (args$jackknife){
    results_jack <- train_survival_jackknife(surv=surv, df=df, name_selection=args$name_selection,
                                             name_model=args$name_model, n_cores=args$n_cores,
                                             verbose=args$verbose, train_foldid=indices$full_foldid)

    # append full data cov metrics to repeat data cov metrics
    results <- list(df_met=bind_rows(results_repeat$df_met, results_full$df_met),
                    df_cov=bind_rows(results_repeat$df_cov, results_full$df_cov, results_jack))
  } else {
    # append full data cov metrics to repeat data cov metrics
    results <- list(df_met=bind_rows(results_repeat$df_met, results_full$df_met),
                    df_cov=bind_rows(results_repeat$df_cov, results_full$df_cov))
  }


  # add info about the subset of features used and the imputation version
  results$df_met <- results$df_met %>% mutate(Features=args$name_features, Table=args$name_table)
  results$df_cov <- results$df_cov %>% mutate(Features=args$name_features, Table=args$name_table)

  # save tables of metrics and covariate coefficients
  dir.create(dirname(args$output_met), showWarnings=F, recursive=T)
  write.table(results$df_met, gsub(".gz", "", args$output_met), sep="\t", row.names=F, quote=F)
  system(paste("gzip",gsub(".gz", "", args$output_met)))
  if (args$verbose) cat(paste("-metrics saved at", args$output_met, "\n"))

  dir.create(dirname(args$output_cov), showWarnings=F, recursive=T)
  write.table(results$df_cov, gsub(".gz", "", args$output_cov), sep="\t", row.names=F, quote=F)
  system(paste("gzip",gsub(".gz", "", args$output_cov)))
  if (args$verbose) cat(paste("-covariates coefficients saved at", args$output_cov, "\n"))
}


# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  default_folder_data <- "../../../results/survival_analysis/data/sub_features/BRCA_dna_and_rna/cln_biol_1_dna_3"
  default_folder_models <- "../../../results/survival_analysis/models/BRCA_dna_and_rna/cln_biol_1_dna_3"

  parser <- ArgumentParser(description='Prepare tables with sub features selected.')
  parser$add_argument("--name_features", type="character", default="cln_biol_1_dna_3",
                      help="Name of the features selection.")
  parser$add_argument("--name_table", type="character", default="imputed_5", help="Name of the completed table.")
  parser$add_argument("--name_selection", type="character", default="10pct_fdr",
                      help="Name of the feature selection method.")
  parser$add_argument("--name_model", type="character", default="coxph_standard", help="Name of the model.")
  parser$add_argument("--input_dat", type="character", help="Path to input data table.",
                    default=file.path(default_folder_data, "processed", "data.imputed_5.tsv.gz"))
  parser$add_argument("--input_cov", type="character", help="Path to input covariates table.",
                    default=file.path(default_folder_data, "processed", "covs.final.tsv.gz"))
  parser$add_argument("--indices", type="character", help="Path to table of cross-validation indices.",
                    default=file.path(default_folder_models, "xval_indices.tsv.gz"))
  parser$add_argument("--repeat_mode", type="character", help="Choose between 'xval' and 'boot'",
                    default='xval')
  parser$add_argument("--jackknife", action="store_true", default=FALSE, 
                      help="If specified, jackknife estimates are computed.")
  parser$add_argument("--verbose", action="store_true", default=F, help="Use to print intermediate messages.")
  parser$add_argument("--n_cores", type="integer", default=4, help="Number of cores to be used.")
  parser$add_argument("--output_met", type="character", help="Path to output table of model quality metrics.",
                    default=file.path(default_folder_models, "10pct_fdr_coxph_standard", "mets.imputed_5.tsv.gz"))
  parser$add_argument("--output_cov", type="character",
                    help="Path to output table of fitted coefficients and/or related quantities.",
                    default=file.path(default_folder_models, "10pct_fdr_coxph_standard", "covs.imputed_5.tsv.gz"))
  parser$add_argument("--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
