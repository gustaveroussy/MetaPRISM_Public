# @created: 25 May 22
# @modified: 27 May 22
# @authors: Yoann Pradat
#
# Run a logistic regression model to assess the associated between platinum drug received and detection of SBS31 and
# SBS35.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(forester))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))

# functions ============================================================================================================

train_logistic <- function(df, col_y, cols_x, ci_conf=0.95, ...){
  # regular logistic regression with penalization
  # it translates in a generalized linear model with Bernoulli distribution and logit link

  formula <- paste(col_y, " ~ ", paste(cols_x, collapse="+"), ...)

  warn <- err <- NULL
  fit_lr <- withCallingHandlers(
    tryCatch(fit_lr <- glm(eval(parse(text=formula)), data=df)
    , error=function(e) {
      err <<- conditionMessage(e)
      NULL
    }), warning=function(w) {
      warn <<- append(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

  if (!is.null(err)){
    df_cov <- tibble(Covariate=cols_x, Estimate=NA)
  } else {
    # retrieve coeffient estimated and stds
    df_cov <- summary(fit_lr)$coefficients %>% as.data.frame() %>% rownames_to_column(var="Covariate")

    # compute adjusted odds ratios estimates and confidence intervals
    z_norm <- qnorm(max(1-ci_conf/2, ci_conf + (1-ci_conf)/2))
    df_cov["Estimate Low"] <- df_cov["Estimate"] - z_norm * df_cov["Std. Error"]
    df_cov["Estimate High"] <- df_cov["Estimate"] + z_norm * df_cov["Std. Error"]
  }

  list(df_cov=df_cov, fit_lr=fit_lr, warn=warn, err=err)
}


format_table <- function(dfs_cov, tab_size=2){
  tab_1 <- paste0(rep(" ", tab_size), collapse="")
  tab_2 <- paste0(rep(" ", tab_size*2), collapse="")

  # aggregate dfs_cov
  df_cov <- tibble()

  for (name in names(dfs_cov)){
    df_cov <- bind_rows(df_cov, tibble(Covariate=name))
    df_cov_name <- dfs_cov[[name]]
    df_cov_name <- df_cov_name %>% filter(Covariate!="(Intercept)")
    df_cov_name[["Covariate"]] <- paste0(tab_1, df_cov_name[["Covariate"]])
    df_cov <- bind_rows(df_cov, df_cov_name)
  }

  df_cov
}


draw_plot <- function(df_plot, output, ci_conf=0.95, ...){
  # draw forest plot
  # xmin <- max(-20, floor(1.1*min(df_plot$`Estimate Low`, na.rm=T)))
  # xmax <- min(20, ceiling(1.1*max(df_plot$`Estimate High`, na.rm=T)))
  xmin <- max(-20)
  xmax <- min(30)

  arrow_labels <- c("Negative association", "Positive association")
  coeff_name <- "Log odds ratio"

  forester(left_side_data=df_plot[,1,drop=F],
           estimate=df_plot$Estimate,
           ci_low=df_plot$`Estimate Low`,
           ci_high=df_plot$`Estimate High`,
           display=FALSE,
           null_line_at=0,
           xlim=c(xmin, xmax),
           x_scale_linear=T,
           estimate_precision=2,
           arrows=T,
           arrow_labels=arrow_labels,
           estimate_col_name=paste0(coeff_name, " (", round(ci_conf*100, 1), "% CI)"),
           render_as="pdf",
           lower_header_row=T,
           file_path=output, ...)
}


main <- function(args){
  # load data
  df_cln <- load_table(args$cln)
  df_sig <- load_table(args$sig)
  df_sig <- df_sig %>% column_to_rownames(var="Signature") %>% t() %>% 
    as.data.frame() %>% rownames_to_column(var="Sample_Id_DNA_T")

  # select patients for which mutational signature analysis was performed
  df_cln <- left_join(df_cln, df_sig, by="Sample_Id_DNA_T")
  df_cln <- df_cln %>% filter(!is.na(SBS1), !is.na(Drugs_Before_Biopsy))

  # platinum drugs and signatures
  drugs_platinum <- c("CARBOPLATINE", "CISPLATINE", "OXALIPLATINE")
  sigs_platinum <- c("SBS31", "SBS35")

  # make a column for each platinum drugs
  col_drug <- "Drugs_Before_Biopsy"
  for (drug in drugs_platinum){
    df_cln[[drug]] <- as.numeric(grepl(drug, df_cln[[col_drug]]))
  }


  # select predictors and outcome
  dfs_cov <- list()
  ci_conf <- 0.95

  for (sig in sigs_platinum){
    out <- train_logistic(df=df_cln, col_y=sig, cols_x=drugs_platinum, ci_conf=ci_conf)
    dfs_cov[[sig]] <- out$df_cov
  }

  # draw plot and save
  df_plot <- format_table(dfs_cov)
  dir.create(file.path(dirname(args$output)), showWarnings=F, recursive=T)
  draw_plot(df_plot, output=args$output, ci_conf=ci_conf)
  draw_plot(df_plot, output=args$output_paper, ci_conf=ci_conf)
  cat(paste("-plot saved at", args$output, "\n"))
  cat(paste("-plot saved at", args$output_paper, "\n"))
}


# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  default_folder <- "../../../results/mutational_signatures/projection_known_signatures/MutationalPatterns"
  basename_sig <- "counts_signatures_cosmic_sbs_96_v3.2_sbs_96_min_mut_sparse_sigprofilerjulia"

  parser <- ArgumentParser(description='Prepare tables with sub features selected.')
  parser$add_argument("--cln", type="character", help="Name of the features selection.",
                      default="../../../data/prism/clinical/curated/cln_prism_in_design_curated.tsv")
  parser$add_argument("--sig", type="character", help="Name of the features selection.",
                      default=file.path(default_folder, paste(basename_sig, "prism.tsv", sep="_")))
  parser$add_argument("--output", type="character", help="Path to output plot.",
                      default="../../../results/mutational_signatures/associations/logistic_regression_platinum.pdf")
  parser$add_argument("--output_paper", type="character", help="Path to output plot.",
                      default="../../../results/figures_paper/F2b_bot.pdf")
  parser$add_argument("--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
