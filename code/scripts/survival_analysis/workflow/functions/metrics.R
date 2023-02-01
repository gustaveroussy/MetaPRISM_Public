# @created Jun 22 2020
# @modified Jun 22 2020
# @author: Yoann Pradat
# 
#     CentraleSupelec
#     MICS laboratory
#     9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France
#
# Metrics for assessing the quality of survival models.

# Libraries
suppressMessages(library(Rcpp))
suppressMessages(library(glmnet))
suppressMessages(library(rms))

#### MAIN FUNCTIONS ====================================================================================================

assess_cox_model_quality <- function(surv, df, linear_predictors, taus, ipcw, indices=NULL){
  if (is.null(indices)) indices <- seq(1, nrow(df))

  # table of model quality metrics
  df_met <- tibble(Time=numeric(), C_simple=numeric(), C_ipcw=numeric(), F_comp=numeric())

  # retrieve ipcw weights for indices
  ipcw_ind <- list(weight.i  = ipcw[["weight.i"]][indices],
                   weight.ij = ipcw[["weight.ij"]][indices, indices],
                   weight.it = ipcw[["weight.it"]][indices,],
                   fit.cph   = ipcw[["fit.cph"]])


  # Simple and IPCW estimates of c-index and Brier score
  df_met <- tibble()

  for (i_tau in 1:length(taus)){
    tau <- taus[i_tau]

    # simple and ipcw c-index
    if (!all(is.na(linear_predictors))){
      c_stats_ind <- c_index_ipcw(surv = surv[indices,],
                                  risk = linear_predictors[indices],
                                  ipcw = ipcw_ind,
                                  tau  = tau)

      # compute estimates of survival functions at tau for train
      prob_ind <- cph_survest(lp = linear_predictors[indices],
                              surv = surv[indices,],
                              times = tau)

      # ipcw brier
      brier_score_ind <- brier_score_ipcw(surv  = surv[indices,],
                                          prob  = prob_ind,
                                          ipcw  = ipcw_ind,
                                          i.tau = i_tau,
                                          tau   = tau)

    } else {
      c_stats_ind <- list("s.index"=NA, "c.index"=NA)
      brier_score_ind <- NA
    }
    
    # write all metrics in a table
    df_met_ind <- tibble(Time=tau, F_comp=c_stats_ind[["f.comp"]],
                         C_simple=c_stats_ind[["s.index"]], C_ipcw=c_stats_ind[["c.index"]], 
                         Brier=brier_score_ind)

    df_met <- bind_rows(df_met, df_met_ind)
  }

  df_met
}


assess_cox_model_quality_train_test <- function(surv, df, linear_predictors, taus, ipcw, train_indices, test_indices){
  df_met_train <- assess_cox_model_quality(surv, df, linear_predictors, taus, ipcw, train_indices)
  df_met_test <- assess_cox_model_quality(surv, df, linear_predictors, taus, ipcw, test_indices)

  df_met_train <- df_met_train %>% mutate(Split="Train")
  df_met_test <- df_met_test %>% mutate(Split="Test")

  bind_rows(df_met_train, df_met_test)
}


#### UTIL FUNCTIONS ====================================================================================================

base_survest <- function (response, lp, times.eval = NULL, centered = FALSE){
  # copied from https://github.com/fbertran/c060/blob/main/R/peperr_glmnet.R
  # baseline survival/ hazard Breslow estimator
  # function essentially based on gbm::basehaz.gbm

  if (is.null(times.eval)) times.eval <- sort(unique(response[,1]))
  
  t.unique <- sort(unique(response[,1][response[,2] == 1]))
  alpha    <- length(t.unique)

  for (i in 1:length(t.unique)) {
    alpha[i] <- sum(response[,1][response[,2] == 1] == t.unique[i])/sum(exp(lp[response[,1] >=  t.unique[i]]))
  }

  obj <- approx(t.unique, cumsum(alpha), yleft=0, xout = times.eval, rule=2)

  if (centered) obj$y <- obj$y * exp(mean(lp))
  obj$z <- exp(-obj$y)

  names(obj) <- c("times","cum_base_haz","base_surv")
  obj
}


cph_survest <- function(lp, surv, times){
  bsurv <- base_survest(surv, lp, sort(unique(times)))

  # reorder to match order of times specifies
  orders <- match(times, bsurv$times)
  cum_base_haz <- bsurv$cum_base_haz[orders]

  exp(exp(lp) %*% -t(cum_base_haz))
}


#### MODELS FOR COMPUTING IPCW =========================================================================================

# We prefer pec::cindex over survival::survConcordance as the latter depends on the censoring distribution (bias)
# Read "Estimating a time-dependent concordance index for survival prediction models with covariate dependent 
# censoring", T.A Gerds, M.W Katter et al. 2012

# If censoring is not independent of predictors, Kaplan-Meier estimator for estimating censoring distribution
# leads to a biased estimate of the c-index.
# We choose to estimate it with a Cox regression model on censoring times

ipcw_cox <- function(surv, df, taus, foldid=NULL, regularized=T, verbose=T){
  # reverse status so as to build a survival model for censoring.
  surv[, "status"] <- 1 - surv[,"status"]

  # if no or too little censoring events, it is not possible to estimate the censoring distribution
  n_censoring <- sum(surv[,"status"])
  not_enough_censoring <- n_censoring==0

  # for some unclear reasons, it is required to explicit the formula
  formula <- eval(parse(text=paste("surv ~ ", paste(colnames(df), collapse="+"))))
    warn <- err <- NULL

  # fit cox proportional hazards model to the censoring times
  if (regularized & !not_enough_censoring){
    if (verbose)
      cat("-fitting regularized cox model (cv.glmnet) for ipcw weights estimation ... ")

      fit.cph <- withCallingHandlers(
        tryCatch(fit.cph <- glmnet::cv.glmnet(x = as.matrix(df),
                                              y = surv,
                                              type.measure = "C",
                                              foldid = foldid,
                                              family = "cox")
        ,error=function(e) {
          err <<- conditionMessage(e)
          NULL
        }), warning=function(w) {
          warn <<- append(warn, conditionMessage(w))
          invokeRestart("muffleWarning")
        })

    if (verbose)
      cat("done!\n")

    if(!is.null(err)){
      not_enough_censoring <- grepl("probably too many censored observations", err)
    } else {
      linear.predictors <- as.numeric(predict(fit.cph, newx=data.matrix(df), s=fit.cph$lambda.min, type="link"))
    }

  } else if (!not_enough_censoring) {
    if (verbose)
      cat("-fitting unregularized cox model (rms::cph) for ipcw weights estimation ... ")

    fit.cph <- withCallingHandlers(
      tryCatch(fit.cph <- rms::cph(formula = formula,
                                   data = df,
                                   method = "efron",   #### Efron approximation method for handling ties
                                   surv = T,           #### Set to TRUE for "speeding up computation of surv probas"
                                   x = T,              #### add x in the return list
                                   y = T)              #### add y in the return listx)
      ,error=function(e) {
        err <<- conditionMessage(e)
        NULL
      }), warning=function(w) {
        warn <<- append(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      })

    if (verbose)
      cat("done!\n")

    linear.predictors <- as.numeric(predict(fit.cph, newx=data.matrix(df)))
  } else {
    fit.cph <- NULL
  }

  # print error/warning
  if (!is.null(warn) & verbose)
    cat(paste("-model warning:", warn, "\n"))

  if (!is.null(err) & verbose)
    cat(paste("-model error:", err, "\n"))

  if (not_enough_censoring | is.null(fit.cph)){
    cat(paste("-all weights are set to 1 as the model is null and/or there is not enough censoring events\n"))
    weight.ij <- matrix(1, nrow(df), nrow(df))
    weight.it <- matrix(1, nrow(df), length(taus))
  } else {
    # quantify the quality of the model using simple c-index (all weights equal to 1)
    ipcw <- list(weight.ij=matrix(1,nrow(df),nrow(df)), weight.it=matrix(1,nrow(df),length(taus)))
    c.stats.train <- c_index_ipcw(surv = surv,
                                  risk = linear.predictors,
                                  ipcw = ipcw,
                                  tau  = taus[length(taus)])

    if (verbose)
      cat(paste("-ipcw cox model (uncorrected) c-index:", c.stats.train$c.index, "\n"))

    # in pec R package, the weights  w_i = P(C > t_i | X = x_i) for i = 1, ..., n
    # are estimated using
    #
    #    weight.i <- rms::survest(
    #        fit   = fit.cph,
    #        times = surv[,"time"],    #### ipcw in pec package uses surv[,"time"] - min(diff(c(0,surv["time"])))/2
    #                                  #### the right term is 0 whenever min(surv[,"time"]) is 0 or if there are ties
    #                                  #### note: in pec, surv["time"] and df.cox are sorted (!)
    #        what  = "parallel"
    #    )
    #
    # while the weights w_ij = P(C > t_j | X = x_i) for i = 1, ..., n, j = 1, ..., n are
    # estimated as done below. Theoretically, the diagonal of weights.ij should be rigourously equal
    # to weight.i but that is not the case at tied times (!).
    #
    # Moreover, it is not specified in the docs how Cox baseline's hazard is estimated (Breslow estimator?).

    # compute w_ij = P(C > t_j | X = x_i) for i = 1, ..., n, j = 1, ..., n
    # weight.ij <- rms::survest(
    #     fit    = fit.cph,
    #     times  = surv[,"time"],
    #     se.fit = F
    # )$surv

    # compute w_it = P(C > t | X = x_i) for i = 1, ..., n, t in taus
    # weight.it <- rms::survest(
    #     fit    = fit.cph,
    #     times  = taus,
    #     se.fit = F
    # )$surv

    # compute w_ij = P(C > t_j | X = x_i) for i = 1, ..., n, j = 1, ..., n
    weight.ij <- cph_survest(lp=linear.predictors,
                             surv=surv,
                             times=surv[,"time"])

    # compute w_it = P(C > t | X = x_i) for i = 1, ..., n, t in taus
    weight.it <- cph_survest(lp=linear.predictors,
                             surv=surv,
                             times=taus)
  }


  return(list(weight.ij = weight.ij, weight.it = weight.it, fit.cph = fit.cph, warn=warn, err=err))
}


#### # UNBIASED (through IPWC) C-INDEX =================================================================================

# Use a custom C-implementation of the IPCW estimate of the c-index.

code_compute_concordance <- ("
    NumericVector compute_concordance(NumericVector time, IntegerVector status, NumericVector risk, 
                                      NumericMatrix weight_ij, double t){
        int i, j;
        int n = time.size(); 
        
        double N_con=0;
        double N_cmp=0;
        double N_tot=0;
        double W_con=0;
        double W_cmp=0;
        double w_ij=0;
        NumericVector result(5);

        for(i=0; i < n; i++){
            for(j=0; j < n; j++){
                // weight_ij(i,j) is P(C > t_j | X = x_i) for i = 1, ..., n, j = 1, ..., n
                // however we want P(C > t_i | X = x_j) in the formula of the ipcw c-index
                w_ij = weight_ij(j,i);

                if(time(i) < time(j) & time(i) <= t){
                    N_tot = N_tot + 1;

                    if (status(i) == 1){
                        N_cmp = N_cmp + 1;
                        W_cmp = W_cmp + 1/(weight_ij(i,i) * w_ij);

                        if(risk(i) > risk(j)){
                            N_con = N_con + 1;
                            W_con = W_con + 1/(weight_ij(i,i) * w_ij);
                        }
                        if(risk(i) == risk(j)){
                            N_con = N_con + 0.5;
                            W_con = W_con + 0.5/(weight_ij(i,i) * w_ij);
                        }
                    }
                }
                else{
                    if(time(i) == time(j) & time(i) <= t){
                        N_tot = N_tot + 1;

                        if (status(i) == 1 & status(j) == 0){
                            N_cmp = N_cmp + 1;
                            W_cmp = W_cmp + 1/(weight_ij(i,i) * w_ij);

                            if(risk(i) > risk(j)){
                                N_con = N_con + 1;
                                W_con = W_con + 1/(weight_ij(i,i) * w_ij);
                            }
                            if(risk(i) == risk(j)){
                                N_con = N_con + 0.5;
                                W_con = W_con + 0.5/(weight_ij(i,i) * w_ij) ;
                            }
                        }

                        /* for now, the case below is considered as an incomparable pair
                          if (status(i)==1 & status(j)==1){
                            // the pair (i,j) will be seen in both directions. add 0.5 so that in total it will count
                            // for 1.
                            N_cmp = N_cmp + 0.5;
                            W_cmp = W_cmp + 0.5/(weight_ij(i,i) * w_ij);

                            if(risk(i) == risk(j)){
                                N_con = N_con + 0.5;
                                W_con = W_con + 0.5/(weight_ij(i,i) * w_ij) ;
                            } else {
                                N_con = N_con + 0.25;
                                W_con = W_con + 0.25/(weight_ij(i,i) * w_ij) ;
                            }
                        } */
                    }
                }
            }
        }

        result(0) = W_con;
        result(1) = W_cmp;
        result(2) = N_con;
        result(3) = N_cmp;
        result(4) = N_tot;
        return(result);
    }
")

cppFunction(code_compute_concordance, showOutput=F, verbose=F)

c_index_ipcw <- function(surv, risk, ipcw, tau){
    c.conc <- compute_concordance(
        surv[,"time"], 
        surv[,"status"], 
        risk, 
        ipcw[["weight.ij"]], 
        tau
    )
    c.index <- c.conc[1]/c.conc[2]  #### IPCW estimate
    s.index <- c.conc[3]/c.conc[4]  #### simple estimate
    f.comp  <- c.conc[4]/c.conc[5]  #### ratio comparable pairs at t / possible pairs at t

    return(list(c.index = c.index, s.index = s.index, f.comp = f.comp))
}

#### # UNBIASED (through IPWC) BRIER SCORE =============================================================================

code_compute_brier_score <- ("
    double compute_brier_score(NumericVector time, IntegerVector status, NumericVector prob_it, NumericVector weight_it, 
                               NumericMatrix weight_ij, double t){
        int i, j;
        int n = time.size(); 
        
        double B=0;

        for(i=0; i < n; i++){
            if(time(i) <= t & status(i) == 1){
                B = B + prob_it(i) * prob_it(i) / (n * weight_ij(i,i));
            }
        
            if(time(i) > t){
                B = B + (1-prob_it(i)) * (1-prob_it(i)) / (n * weight_it(i));
            }
        }

        return(B);
    }
")

cppFunction(code_compute_brier_score)


brier_score_ipcw <- function(surv, prob, ipcw, i.tau, tau){
    brier <- compute_brier_score(
        surv[,"time"], 
        surv[,"status"], 
        prob, 
        ipcw[["weight.it"]][,i.tau], 
        ipcw[["weight.ij"]], 
        tau
    )

    return(brier)
}
