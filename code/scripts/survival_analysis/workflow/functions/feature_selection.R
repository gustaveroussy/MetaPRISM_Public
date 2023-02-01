# @created Apr 13 2022
# @modified May 12 2022
# @author: Yoann Pradat
# 
#     CentraleSupelec
#     MICS laboratory
#     9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France
#
# Functions to perform feature selection.


# Libraries
suppressMessages(library(dplyr))
suppressMessages(library(survival))
suppressMessages(library(glmnet))
suppressMessages(library(rms))
suppressMessages(library(parallel))
suppressMessages(library(doParallel))
suppressMessages(library(foreach))
suppressMessages(library(stringr))

source("workflow/functions/survival_models.R")

#### # MAIN FUNCTIONS ==================================================================================================

feature_selection <- function(surv, df, name_selection, name_model=NULL, train_indices=NULL, train_foldid=NULL){
  if (is.null(train_indices)) train_indices <- 1:nrow(df)

  # it may happen that in a cross-validation split, one column is uniquely-valued in train-indices
  # if that happens, that variable is removed
  covs_all <- colnames(df)
  covs_bad <- c()
  for (cov in covs_all){
    if (length(unique(df[train_indices,cov,drop=T]))==1){
      covs_bad <- c(covs_bad, cov)
    }
  }

  if (length(covs_bad) > 0){
    cat(paste("-WARNING: the following covariates are uniquely-valued in the train_indices provided\n\t"))
    cat(paste(covs_bad, collapse="\n\t"),"\n")
    df_fs <- df[,setdiff(covs_all, covs_bad)]
  } else {
    df_fs <- df
  }

  if (tolower(name_selection)=="none"){
    df_fs <- df_fs
  } else if (grepl("pct_fdr$", name_selection)){
    fdr_level <- as.numeric(unlist(str_split(name_selection, "pct_fdr"))[1])/100
    df_fs <- feature_selection_fdr(surv, df_fs, train_indices, fdr_level)
  } else if (tolower(name_selection)=="lasso"){
    df_fs <- feature_selection_lasso(surv, df_fs, train_indices, train_foldid)
  } else {
    stop(paste("-feature selection", name_selection, "not implemented"))
  }

  df_fs
}

#### # 1. FDR FEATURE SELECTION ========================================================================================

feature_selection_fdr <- function(surv, df, train_indices, fdr_level){
  covs_all <- colnames(df)

  # get p-value of association of each covariate with survival
  covs_pvals <- c()
  for (cov in covs_all){
    # it may happen that in a cross-validation split, one column is uniquely-valued in train-indices
    # if that happens, that variable is removed
    if (length(unique(df[train_indices,cov]==0))==1){
      pval <- 1
    } else {
      formula <- paste("surv[train_indices] ~ ", paste(cov, collapse="+"))
      fit_cph <- survival::coxph(formula = eval(parse(text=formula)),
                                 data    = df[train_indices,])
      fit_coefs <- summary(fit_cph)$coefficients
      pval <- fit_coefs[,"Pr(>|z|)"]
      if (is.na(pval)){
        pval <- 1
      }
    }
    covs_pvals <- c(covs_pvals, pval)
  }

  # correct for multiple testing
  covs_qvals <- p.adjust(covs_pvals, method="fdr")

  # retain only covariates that pass the fdr level
  covs_in <- covs_all[covs_qvals <= fdr_level]

  if (length(covs_in)==0){
    warning(paste("The FDR feature selection method at level", fdr_level, "has selected 0 variable!",
                  "As a consequence, no feature selection is applied."))
    covs_in <- covs_all
  } else {
    cat(paste("-INFO:", paste0(length(covs_in), "/", length(covs_all)), "features were selected using FDR method at",
              "level", fdr_level, "\n"))
  }

  df[,covs_in]
}


#### # 2. LASSO FEATURE SELECTION ======================================================================================

feature_selection_lasso <- function(surv, df, train_indices, train_foldid){
  covs_all <- colnames(df)

  # run lasso
  out <- train_coxph_glmnet(surv, df, train_indices=train_indices, train_foldid=train_foldid, alpha=1)

  df_cov <- out$df_cov
  covs_in <- df_cov %>% filter(Coefficient!=0) %>% pull(var="Covariate")

  if (length(covs_in)==0){
    warning(paste("The lasso feature selection method has selected 0 variable!",
                  "As a consequence, no feature selection is applied."))
    covs_in <- covs_all
  } else {
    cat(paste("-INFO:", paste0(length(covs_in), "/", length(covs_all)), "features were selected using lasso method\n"))
  }

  df[,covs_in]
}
