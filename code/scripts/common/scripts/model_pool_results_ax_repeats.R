# created: Jan 04 2022
# modified: May 14 2022
# author: Yoann Pradat
# 
#     CentraleSupelec
#     MICS laboratory
#     9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France
# 
#     Institut Gustave Roussy
#     Prism Center
#     114 rue Edouard Vaillant, Villejuif, 94800 France
# 
# Pool coefficient/standard error and quality metrics estimates across bootstraps.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))

# ======================================================================================================================

#' Quantile computation with interpolation on the normal quantile scale.   
#'
#' For a non-integer order statistic this function interpolates between the surrounding
#' order statistics using the normal quantile scale.  See equation
#' 5.8 of Davison and Hinkley (1997)
#'
#' @param t vector of values from which quantiles are computed
#' @param alpha vector of quantiles
#' @return quantile estimate
#'
#' @references Davison, A., & Hinkley, D. (1997). Bootstrap Methods and their Application. Chapter 5.
quantile_with_norm_inter <- function(t,alpha)
{
  t <- t[is.finite(t)]
  if (length(t)==0){
    return(rep(NA, length.out=length(alpha)))
  } else {
    R <- length(t)
    rk <- (R+1)*alpha
    if (!all(rk>1 & rk<R))
      warning("extreme order statistics used as endpoints")
    k <- trunc(rk)
    inds <- seq_along(k)
    out <- inds
    kvs <- k[k>0 & k<R]
    tstar <- sort(t, partial = sort(union(c(1, R), c(kvs, kvs+1))))
    ints <- (k == rk)
    if (any(ints)) out[inds[ints]] <- tstar[k[inds[ints]]]
    out[k == 0] <- tstar[1L]
    out[k == R] <- tstar[R]
    not <- function(v) xor(rep(TRUE,length(v)),v)
    temp <- inds[not(ints) & k != 0 & k != R]
    temp1 <- qnorm(alpha[temp])
    temp2 <- qnorm(k[temp]/(R+1))
    temp3 <- qnorm((k[temp]+1)/(R+1))
    tk <- tstar[k[temp]]
    tk1 <- tstar[k[temp]+1L]
    out[temp] <- tk + (temp1-temp2)/(temp3-temp2)*(tk1 - tk)
    return(out)
  }
}

#'  Basic bootstrap confidence interval
#'
#' @param t0 parameter estimate on the full data.
#' @param tb vector of parameter estimates across bootstraps.
#' @param conf confidence level of the interval
ci_basic <- function(t0, tb, conf){
  alpha <- (1+c(conf,-conf))/2
  qq <- quantile_with_norm_inter(tb,alpha)
  ci_1 <- 2*t0 - qq[1]
  ci_2 <- 2*t0 - qq[2]
  list(ci_1=ci_1, ci_2=ci_2)
}


#' Compute percentile interval
#'
#' @param tb vector of estimates of the parameter of interest across bootstrap iterations.
#' @param conf confidence level of confidence interval
#'
#' @references Efron, B. & Tibshirani, R. J. (1994). An introduction to the bootstrap. Chapman & Hall/CRC. Chapter 13.
ci_percentile <- function(tb, conf){
  alpha <- (1+c(conf,-conf))/2
  qq <- quantile_with_norm_inter(tb,alpha)
  list(ci_1=qq[2], ci_2=qq[1])
}


#' Compute bootstrap-t interval
#'
#' Studentized version of the basic bootstrap confidence interval
#'
#' @param t0 parameter estimate on the full data.
#' @param tb vector of parameter estimates across bootstraps.
#' @param var_t0 parameter variance estimate on the full data.
#' @param var_tb vector of parameter variance estimates across bootstraps.
#' @param conf confidence level of the interval
#'
#' @references Efron, B.& Tibshirani, R. J. (1994). An introduction to the bootstrap. Chapman & Hall/CRC. Chapter 12.
ci_student <- function(t0, tb, var_t0, var_tb, conf){
  z <- (tb-t0)/sqrt(var_tb)
  alpha <- (1+c(conf,-conf))/2
  qq <- quantile_with_norm_inter(z,alpha)
  ci_1 <- t0-sqrt(var_t0)*qq[1]
  ci_2 <- t0-sqrt(var_t0)*qq[2]
  list(ci_1=ci_1, ci_2=ci_2)
}


#' Compute BCa interval
#'
#' @param t0 Estimate of the parameter of interest on the full dataset.
#' @param tb vector of estimates of the parameter of interest across bootstrap iterations.
#' @param tj Vector of estimates of the parameter of interest across jackknife iterations.
#' @param conf confidence level of the interval
#' @param use_t0_over_tj_mean If TRUE, t0 is used in place of tj_mean in the estimation of the accelation parameter.
#'   In the `boot` R package, computing BCa interval using usual jackknife for estimating the acceleration parameter
#'   leads to the use of tobs=t0 in the formula. See the internal function `usual.jack`.
#'   In the equation (14.15) of the book by Efron and Tibshirani, the mean of jackknife estimates, tj_mean, is used.
#'
#' @references Efron, B. & Tibshirani, R. J. (1994). An introduction to the bootstrap. Chapman &  Hall/CRC. Chapter 14.
#' Davison, A., & Hinkley, D. (1997). Bootstrap Methods and their Application. Chapter 5.
ci_bca <- function(t0, tb, tj, conf, use_t0_over_tj_mean=FALSE){
  # bias-correction parameter
  z_0 <- qnorm(sum(tb < t0)/(length(tb)+1))
  if (!is.finite(z_0)){
    z_0 <- 0.5
  }

  # acceleration parameter
  if (use_t0_over_tj_mean){
    tuse <- mean(tj)
  } else {
    tuse <- t0
  }
  L <- (tuse-tj)
  if (all(L==0)){
    a <- 0
  } else {
    a <- sum(L^3)/(6*(sum(L^2)^(1.5)))
  }

  # ajdusted quantiles
  z_alpha <- qnorm((1+c(-conf,conf))/2)
  adj_alpha <-  pnorm(z_0 + (z_0+z_alpha)/(1-a*(z_0+z_alpha)))
  qq <- quantile_with_norm_inter(tb, adj_alpha)

  list(ci_1=qq[1], ci_2=qq[2])
}


#' Bootstrap confidence interval
#'
#' @param t0 parameter estimate on the full data.
#' @param tb vector of parameter estimates across bootstraps.
#' @param conf confidence level of the interval
#' @param type type of confidence interval. Possible values are 'basic', 'percentile', 'student' and 'bca'.
#' @param tj (optional) vector of parameter estimates across jackknifes. Used only for type='bca'.
#' @param var_t0 (optional) parameter variance estimate on the full data.
#' @param var_tb (optional) vector of parameter variance estimates across bootstraps.
#'
#' @references Efron, B., Tibshirani, R., & Tibshirani, R. J. (1994). An introduction to the bootstrap. Chapman &
#' Hall/CRC. Chapters 12,13,14.
#' Davison, A., & Hinkley, D. (1997). Bootstrap Methods and their Application. Chapter 5.
ci_bootstrap <- function(t0, tb, conf, type=c("basic", "percentile", "student", "bca"), tj=NULL, var_t0=NULL,
                         var_tj=NULL){
  type <- match.arg(type)
  switch(type,
         basic = ci_basic(t0, tb, conf),
         percentile = ci_percentile(tb, conf),
         student = ci_student(t0, tb, var_t0, var_tb, conf),
         bca = ci_bca(t0, tb, tj, conf))
}


poolstr <- function(x){
  x_l <- x[!is.na(x)]
  x_uni_all <- sort(unique((x_l)))
  if (length(x_uni_all)==0){
    return(NA)
  } else {
    x_cnt_all <- c()
    for (x_uni in x_uni_all){
      x_cnt_all <- c(x_cnt_all, sum(x_l==x_uni))
    }
    x_uni_cnt <- paste(x_uni_all, paste0("[x", x_cnt_all, "]"))
    return(paste(x_uni_cnt, collapse=" | "))
  }
}



get_point_estimates <- function(df, cols_gby, cols_est, cols_str=NULL, drop_warn=F, drop_error=T){
  dfs_mean <- list()
  dfs_info <- list()
  
  if ("Grid_Run" %in% colnames(df)){
    mask <- (df$Run != "Full") & (df$Grid_Run != "Train")
  } else {
    mask <- df$Run != "Full"
  }

  mask_e <- is.na(df$Error)
  mask_w <- is.na(df$Warning)
  if (drop_error=="true" & drop_warn=="true"){
    mask_we <- mask_e & mask_w
    if (mean(mask_we < 0.1)){
      mask_we <- mask_e
      cat(paste("-WARNING: removing models with warnings or errors would remove over 90% of models. Remove only models",
                "with errors.\n"))
    }
  } else if (drop_error=="true"){
    mask_we <- mask_e
  } else if (drop_warn=="true"){
    mask_we <- mask_w
    if (mean(mask_we < 0.1)){
      mask_we <- rep(T, times=nrow(df)) 
      cat("-WARNING: removing models with warnings would remove over 90% of models. Keep all models.\n")
    }
  } else {
    mask_we <- rep(T, times=nrow(df))
  }

  for (col_est in cols_est){
    dfs_mean[[col_est]] <- df %>% filter(mask & mask_we) %>%
      group_by_at(vars(all_of(cols_gby))) %>%
      group_modify(~ data.frame(Mean=mean(.x[[col_est]], na.rm=T))) %>%
      rename(!!col_est:=Mean)
  }

  for (col_str in cols_str){
    dfs_info[[col_str]] <- df %>% filter(mask) %>%
      group_by_at(vars(all_of(cols_gby))) %>%
      group_modify(~ data.frame(Info=poolstr(.x[[col_str]]), N_Info=sum(!is.na(.x[[col_str]])))) %>%
      rename(!!col_str:=Info, !!paste0("N_", col_str):=N_Info)
  }

  if ("Selected" %in% colnames(df)){
    # add number of selections 
    dfs_info[["Selected"]] <- df %>% filter(mask & mask_we) %>% group_by_at(vars(all_of(cols_gby))) %>% 
        group_modify(~ data.frame(N_Selected=sum(.x[["Selected"]])))
  }

  # add total number of repeats
  dfs_info[["Total"]] <- df %>% filter(mask) %>% group_by_at(vars(all_of(cols_gby))) %>% 
    summarize(N_Repeats=n(), .groups="drop")
  
  Reduce(function(df1, df2) left_join(df1, df2, by=cols_gby), c(dfs_mean, dfs_info))
}


get_confidence_intervals <- function(df, cols_gby, cols_ci, types, repeat_mode, conf=0.95){
  dfs_ci <- list()
  if (repeat_mode=="xval"){
    prefix <- "CV"
  } else if (repeat_mode=="boot") {
    prefix <- "Boot"
  }

  for (col_ci in cols_ci){
    for (type in types){
      if (type=="bca" & sum(grepl("Jack", df$Run))==0){
        warning(paste("-cannot compute BCa confidence interval for", col_ci, "as no jackknife estimates available"))
      } else {
        name <- paste0(col_ci, "_", type)
        dfs_ci[[name]] <- df %>%
          group_by_at(vars(all_of(cols_gby))) %>%
          group_modify(~ data.frame(ci_bootstrap(t0=.x[grepl("Full", .x$Run),col_ci,drop=T],
                                                 tb=.x[grepl(prefix, .x$Run),col_ci,drop=T],
                                                 tj=.x[grepl("Jack", .x$Run),col_ci,drop=T],
                                                 conf=conf,
                                                 type=type))) %>%
          rename(!!paste0(col_ci, "_", "CI_Low_", str_to_title(type)):=ci_1,
                 !!paste0(col_ci, "_", "CI_High_", str_to_title(type)):=ci_2)
      }
    }
  }

  Reduce(function(df1, df2) left_join(df1, df2, by=cols_gby), dfs_ci)
}


main <- function(args){
  df_met <- read_tsv(args$input_met, progress=F, show_col_types=F)
  df_cov <- read_tsv(args$input_cov, progress=F, show_col_types=F)

  cols_gby_met <- c("Model", "Features", "Selection", "Time", "Split")
  cols_gby_met <- intersect(colnames(df_met), cols_gby_met)
  cols_str_met <- c("Warning", "Error", "Warning_Selection")
  cols_str_met <- intersect(colnames(df_met), cols_str_met)
  cols_est_met <- setdiff(colnames(df_met), Reduce(union, list(cols_gby_met, c("Run", "Grid_Run"), cols_str_met)))
  cols_ci_met <- cols_est_met

  cols_gby_cov <- c("Model", "Features", "Selection", "Covariate")
  cols_gby_cov <- intersect(colnames(df_cov), cols_gby_cov)
  cols_str_cov <- c("Warning", "Error", "Warning_Selection")
  cols_str_cov <- intersect(colnames(df_cov), cols_str_cov)
  cols_est_cov <- intersect(colnames(df_cov), c("Coefficient", "Standard_Error", "Selected"))
  cols_ci_cov <- c("Coefficient")

  dfs <- list(df_met, df_cov)
  cols_gbys <- list(cols_gby_met, cols_gby_cov)
  cols_ests <- list(cols_est_met, cols_est_cov)
  cols_cis <- list(cols_ci_met, cols_ci_cov)
  cols_strs <- list(cols_str_met, cols_str_cov)

  # compute point estimates by averaging estimates across bootstraps 
  dfs_point_est <- list()
  for (i in 1:length(dfs)){
    dfs_point_est[[i]] <- get_point_estimates(df=dfs[[i]], cols_gby=cols_gbys[[i]], cols_est=cols_ests[[i]], 
                                              cols_str=cols_strs[[i]], drop_warn=args$drop_warn, 
                                              drop_error=args$drop_error)
  }

  # compute confidence intervals using bootstraps iterations
  # several types of confidence intervals may be computed
  # the following are computed where possible
  #  - basic
  #  - percentile
  #  - bca (only possible if jackknife estimates were computed)

  # NB: for quality metrics, estimates on the full data set are only available for the train data
  # Applying bootstrap statistical theory on test quality metrics may in fact be incorrect given that these metrics are
  # not proper statistics in the sense that they are not functions of the observations on which the model was trained
  # but of observations outside the model training sample. For compatibility with the function for estimating confidence
  # intervals, an artificial "full data set" value will be computed for test quality metrics by averaging over bootstrap
  # estimates.
  dfs[[1]] <- bind_rows(dfs[[1]], dfs_point_est[[1]] %>% filter(Split=="test") %>% mutate(Run="Full"))

  dfs_ci_est <- list()
  for (i in 1:length(dfs)){
    dfs_ci_est[[i]] <- get_confidence_intervals(df=dfs[[i]], cols_gby=cols_gbys[[i]], cols_ci=cols_cis[[i]],
                                                types=args$ci_types, repeat_mode=args$repeat_mode, conf=args$ci_conf)
  }

  # merge point estimates and confidence interval estimates and save
  df_met_pool <- left_join(dfs_point_est[[1]], dfs_ci_est[[1]], by=cols_gbys[[1]])
  df_cov_pool <- left_join(dfs_point_est[[2]], dfs_ci_est[[2]], by=cols_gbys[[2]])

  write.table(df_met_pool, gsub(".gz", "", args$output_met), sep="\t", row.names=F, quote=F)
  system(paste("gzip",gsub(".gz", "", args$output_met)))
  cat(paste("-metrics saved at", args$output_met, "\n"))

  write.table(df_cov_pool, gsub(".gz", "", args$output_cov), sep="\t", row.names=F, quote=F)
  system(paste("gzip",gsub(".gz", "", args$output_cov)))
  cat(paste("-covariates coefficients saved at", args$output_cov, "\n"))
}


# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  results_folder <- "../../../results/survival_analysis"
  default_folder <- paste0(results_folder, "/models/BRCA_dna_and_rna/cln_biol_2/lasso_coxph_standard")

  parser <- ArgumentParser(description='Pool model quality metrics and coefficients estimates across bootstraps.')
  parser$add_argument("--repeat_mode", type="character", help="Choose between xval and boot.", default="xval")
  parser$add_argument("--input_met", type="character", help="Path to input table of pooled model quality metrics.",
                      default=file.path(default_folder, "mets.pooled_ax_imp.tsv.gz"))
  parser$add_argument("--input_cov", type="character", help="Path to input table of pooled model coefficients estimates.",
                      default=file.path(default_folder, "covs.pooled_ax_imp.tsv.gz"))
  parser$add_argument("--drop_warn", type="character", default="true",
                      help="If true, coefficients estimates from model(s) (if imputations) with warning are dropped.")
  parser$add_argument("--drop_error", type="character", default="true",
                      help="If true, coefficients estimates from model(s) (if imputations) with error are dropped.")
  parser$add_argument("--ci_types", nargs="+", type="character", default=c("basic", "percentile"),
                      help="Choose one of 'percent', 'student' or 'bca'. For 'bca' costly computations are required.")
  parser$add_argument("--ci_conf", type="double", default=0.95,
                      help="Confidence level of intervals.")
  parser$add_argument('--output_met', type="character", help='Path to output table of pooled model quality metrics.',
                      default=file.path(default_folder, "mets.pooled_ax_rep.tsv.gz"))
  parser$add_argument('--output_cov', type="character", help='Path to output table of pooled model coefficients estimates.',
                      default=file.path(default_folder, "covs.pooled_ax_rep.tsv.gz"))
  parser$add_argument("--log", type="character", help="Path where log file will be saved.")
  args = parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
