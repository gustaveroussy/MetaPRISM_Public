# @created Apr 13 2022
# @modified May 06 2022
# @author: Yoann Pradat
# 
#     CentraleSupelec
#     MICS laboratory
#     9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France
#
# Linear cox models as implemented in the R libraries survival, glmnet.


# Libraries
suppressMessages(library(dplyr))
suppressMessages(library(survival))
suppressMessages(library(glmnet))
suppressMessages(library(rms))
suppressMessages(library(parallel))
suppressMessages(library(doParallel))
suppressMessages(library(foreach))
suppressMessages(library(stringr))

source("workflow/functions/metrics.R")

#### # MAIN FUNCTIONS ==================================================================================================

train_survival <- function(surv, df, name_model, train_indices=NULL, train_foldid=NULL){
  if (is.null(train_indices)) train_indices <- 1:nrow(df)

  # glmnet cannot be used on a matrix with only 1 column. In this case, disregard penalization, and
  # use standard survival::coxph function
  if (name_model=="coxph_standard" | (name_model %in% c("coxph_lasso", "coxph_ridge") & ncol(df)==1)){
    formula <- paste("surv[train_indices] ~ ", paste(colnames(df), collapse="+"))
    return(train_coxph_formula(surv, df, train_indices, formula))
  } else if (name_model=="coxph_lasso"){
    return(train_coxph_glmnet(surv, df, train_indices, train_foldid, alpha=1))
  } else if (name_model=="coxph_ridge"){
    #return(train_coxph_ridge(surv, df, train_indices))
    return(train_coxph_glmnet(surv, df, train_indices, train_foldid, alpha=0))
  } else {
    stop(paste("-model", name_model, "not implemented"))
  }
}

#### # UTILS ===========================================================================================================

get_formula_coxph_ridge <- function(covs, theta, lim=30){
  # There is width.cutoff of 500 characters somewhere in the code of coxph function that makes the code break if we input
  # one long ridge. A work-around for now is to break it into several
  # If the code still raises the error "missing value where TRUE/FALSE is needed", try changing 40 to lim or 20 or try
  # shortening your variables names
  ncovs <- length(covs)
  ridge_formula <- paste0("ridge(",paste(covs[1:min(lim, length(covs))],collapse=","), ",theta=", theta, ")")
  q <- as.integer((ncovs-1)/lim)
  if (q>=1){
    for (k in 1:q){
      ridge_formula <- paste0(ridge_formula, "+", 
                paste0("ridge(",paste(covs[(lim*k+1):min(lim*(k+1), ncovs)], collapse=","), ",theta=", theta, ")")) 
    }
  }

  paste("surv[train_indices] ~ ", ridge_formula)
}


cv_coxph_ridge <- function(surv, df, tau, ipcw, n_folds, thetas, seed=NULL){
  if (!is.null(seed)) set.seed(seed)
  cv_folds <- sample(rep(seq(n_folds),length=nrow(df)))

  df_met <- foreach (i=1:length(thetas), .combine='rbind') %do% {
    theta <- thetas[[i]]
    
    df_met <- foreach(k=1:n_folds,
                      .combine = 'rbind', 
                      .multicombine = F)  %do% {
      train_indices <- which(!(cv_folds == k))
      test_indices <- which(cv_folds == k)
      formula <- get_formula_coxph_ridge(covs=colnames(df), theta=theta)

      results <- train_test_coxph_formula(surv, df, tau, ipcw, train_indices, test_indices, formula=formula)
      results$df_met %>% mutate(Fold=k)
    }

    df_met %>% mutate(Theta=theta)
  }

  # compute theta_min
  theta_min <- df_met %>% filter(Split=="Test") %>% group_by(Theta) %>% summarise(C_ipcw_avg=mean(C_ipcw)) %>%
    filter(C_ipcw_avg==max(C_ipcw_avg)) %>% pull(Theta)

  list(df_met=df_met, theta_min=theta_min)
}

#### # 1. RMS::CPH, COXPH STANDARD/RIDGE ===============================================================================

# survival package by Terry Therneau implements Cox Proportional Hazards model 
# rms pacakage by Frank Harrell provides additional functions.

train_coxph_formula <- function(surv, df, train_indices, formula, method="efron"){
  # regular coxph without penalty or with l2 penalty (ridge)
  # the log hazard is linear in the covariates
  #
  #    log(\lambda(t)) = log(\lambda_0(t)) + \beta'X
  #
  # the model is fit by maximization of partial likelihood 

  # one can alternatively use rms::cph but it fails if the formula contains ridge penalization
  warn <- err <- NULL
  fit_cph <- withCallingHandlers(
    tryCatch(fit_cph <- survival::coxph(formula = eval(parse(text=formula)),
                                        data    = df[train_indices,],
                                        ties    = method)
    ,error=function(e) {
      err <<- conditionMessage(e)
      NULL
    }), warning=function(w) {
      warn <<- append(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

  if (!is.null(err)){
    beta_hat <- rep(NA, ncol(df))
    beta_names <- colnames(df)
    linear_predictors <- rep(NA, nrow(df))

    # table of covariates estimates
    df_cov <- tibble(Covariate=beta_names,
                     Coefficient=beta_hat,
                     Standard_Error=NA)
  } else {
    # In case the input matrix is singular, some coefficients may be NA. We arbitrarily set them to 0
    # as is done in survival::coxph .....
    beta_hat <- fit_cph$coefficients
    beta_names <- names(beta_hat)
    names(beta_hat) <- NULL

    if (all(grepl("^ridge", beta_names))) beta_names <- str_extract(beta_names, "(?<=ridge\\()\\w+(?=\\))")

    if (sum(is.na(beta_hat)) > 0){
        for (i in 1:sum(is.na(beta_hat))){
            cat(paste("-warning!", names(beta_hat[is.na(beta_hat)])[i], "has coefficient NA. We set it to 0.\n"))
        }
        beta_hat[is.na(beta_hat)] <- 0
    }

    # simple and ipcw estimates of c-index and Brier score at different time horizons
    linear_predictors <- as.numeric(predict(fit_cph, newdata=as.data.frame(df), type="lp"))

    # table of covariates estimates
    df_cov <- tibble(Covariate=beta_names,
                     Coefficient=beta_hat,
                     Standard_Error=sqrt(diag(fit_cph$var)))
  }


  list(df_cov=df_cov, linear_predictors=linear_predictors, warn=warn, err=err)
}


train_coxph_ridge <- function(surv, df, train_indices, theta=NULL, method="efron"){
  if (is.null(theta)){
    thetas <- c(seq(2,10,2) %o% 10^(-1:3))
    cv_results <- cv_coxph_ridge(surv, df, quantile(surv[,1], 0.95), ipcw, n_folds=5, thetas=thetas, seed=123)
    theta <- cv_results$theta_min
  }
  formula <- get_formula_coxph_ridge(covs=colnames(df), theta=theta)
  results <- train_test_coxph_formula(surv, df, train_indices, formula, method=method)
  results$df_cov <- results$df_cov %>% mutate(Theta=theta)

  results
}

#### # 2. GLMNET::CV.GLMNET, COXPH REGULARIZED L1/L2 ===================================================================

train_coxph_glmnet <- function(surv, df, train_indices, train_foldid, alpha=1, ...){
  # regular coxph with l1/l2 penalty (ridge)
  # the log hazard is linear in the covariates
  #
  #    log(\lambda(t)) = log(\lambda_0(t)) + \beta'X
  #
  # the model is fit by maximization of penalized partial likelihood 

  if (!is.null(train_foldid)){
    nfolds <- NULL
  } else {
    nfolds <- 5
  }

  warn <- err <- NULL
  fit_cph <- withCallingHandlers(
    tryCatch(fit_cph <- glmnet::cv.glmnet(x                = as.matrix(df[train_indices,]),
                                          y                = surv[train_indices,],
                                          foldid           = train_foldid,
                                          family           = "cox",
                                          type.measure     = "C",
                                          standardize      = F,
                                          alpha            = alpha,
                                          ...)
    , error=function(e) {
      err <<- conditionMessage(e)
      NULL
    }), warning=function(w) {
      warn <<- append(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

  if (!is.null(err)){
    fit_cph <- withCallingHandlers(
      tryCatch(fit_cph <- glmnet::cv.glmnet(x                = as.matrix(df[train_indices,]),
                                            y                = surv[train_indices,],
                                            foldid           = train_foldid,
                                            family           = "cox",
                                            type.measure     = "deviance",
                                            standardize      = F,
                                            alpha            = alpha,
                                            ...)
      , error=function(e) {
        err <<- paste(err, conditionMessage(e), sep=" & ")
        NULL
      }), warning=function(w) {
        err <<- paste(err, "", sep=" & ")
        warn <<- append(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
  }


  if (!is.null(err)){
    beta_hat <- rep(NA, ncol(df))
    beta_names <- colnames(df)
    linear_predictors <- rep(NA, nrow(df))
  } else {
    # In case the input matrix is singular, some coefficients may be NA. We arbitrarily set them to 0
    # as is done in survival::coxph .....

    beta_hat <- as.vector(coef(fit_cph, s=fit_cph$lambda.min))
    beta_names <- rownames(coef(fit_cph, s=fit_cph$lambda.min))
    colnames(beta_hat) <- NULL

    if (sum(is.na(beta_hat)) > 0){
        for (i in 1:sum(is.na(beta_hat))){
            cat(paste("-warning!", names(beta_hat[is.na(beta_hat)])[i], "has coefficient NA. We set it to 0.\n"))
        }
        beta_hat[is.na(beta_hat)] <- 0
    }

    # simple and ipcw estimates of c-index and Brier score at different time horizons
    linear_predictors <- as.numeric(predict(fit_cph, newx=as.matrix(df), s=fit_cph$lambda.min, type="link"))
  }

  # table of covariates estimates
  df_cov <- tibble(Covariate=beta_names,
                   Coefficient=beta_hat)
  df_cov <- df_cov %>% mutate(Alpha=alpha, Lambda=fit_cph$lambda.min)

  list(df_cov=df_cov, linear_predictors=linear_predictors, warn=warn, err=err)
}
