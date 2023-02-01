
#' Compute pvalues of multiple Fisher tests
#'
#' From 2 matrices of event counts in variables and fixed samples, compute a Fisher test for each 2x2 count table
#' induced by the each pair of cell in table A and table B.
#'
#' @param df_count_a a \code{data.frame}.
#' @param df_count_b a \code{data.frame}.
#' @param margins_a a numeric vector defining the sample size of each column of \code{df_count_a}.
#' @param margins_b a numeric vector defining the sample size of each column of \code{df_count_b}.
#' @param conf_int if set to TRUE, the returned object is a list with names 'pvals', 'conf_int_l', 'conf_int_h'.
#' @param conf_lvl Use only if conf_int set to TRUE
#' @param progress set to TRUE to see progress bar.
#' @param n_cores number of cores to be used for parallel computations.
#' @return a \code{data.frame} with the same row names and column names as \code{df_count_a} and \code{df_count_b}.
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom Exact exact.test
#' @importFrom dplyr bind_rows
#'
#' @author Yoann Pradat
#' @keywords internal
get_pvals_fisher <- function(df_count_a, df_count_b, margins_a, margins_b, conf_int=F, conf_lvl=0.95, progress=F,
                             n_cores=1){
  stopifnot(identical(rownames(df_count_a), rownames(df_count_b)))
  stopifnot(identical(colnames(df_count_a), colnames(df_count_b)))
  stopifnot(identical(colnames(df_count_a), rownames(margins_a)))
  stopifnot(identical(colnames(df_count_b), rownames(margins_b)))

  n_row <- nrow(df_count_a)
  n_col <- ncol(df_count_a)
  df_pval <- data.frame(row.names=rownames(df_count_a))
  df_est <- data.frame(row.names=rownames(df_count_a))

  if (conf_int){
    df_cil <- data.frame(row.names=rownames(df_count_a))
    df_cih <- data.frame(row.names=rownames(df_count_a))
  }
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  for (j in 1:n_col){
    margin_a <- margins_a[j,]
    margin_b <- margins_b[j,]

    outs <- foreach(i = 1:n_row, .packages=c("Exact"), .combine=bind_rows) %dopar% {
      c_a <- df_count_a[i,j]
      c_b <- df_count_b[i,j]
      data <- matrix(c(c_a, margin_a - c_a, c_b, margin_b - c_b), byrow=T, nrow=2)

      if (margin_a > 0 & margin_b > 0){
        out <- exact.test(data=data, alternative="two.sided",
                          np.interval=FALSE,
                          npNumbers=100,
                          method="Boschloo", 
                          model="Binomial",
                          tsmethod="central",
                          conf.int=conf_int, conf.level=conf_lvl,
                          cond.row=TRUE, to.plot=FALSE)
        est <- out$estimate
        if (out$estimate >= 0){
          pval <- out$p.value
        } else {
          pval <- -out$p.value
        }
        if (conf_int){
          cil <- out$conf.int[[1]]
          cih <- out$conf.int[[2]]
        } else {
          cil <- NULL
          cih <- NULL
        }
      } else {
        est <- NA
        pval <- -999
        cil <- NULL
        cih <- NULL
      }

      list(pval=pval, est=est, cil=cil, cih=cih)
    }
    df_pval <- cbind(df_pval, data.frame(outs$pval))
    df_est <- cbind(df_est, data.frame(outs$est))
    if (conf_int){
      df_cil <- cbind(df_cil, data.frame(outs$cil))
      df_cih <- cbind(df_cih, data.frame(outs$cih))
    }

    if (progress){
      progress(j, n_col)
    }
  }
  stopCluster(cl)

  colnames(df_pval) <- colnames(df_count_a)
  if (conf_int){
    colnames(df_cil) <- colnames(df_count_a)
    colnames(df_cih) <- colnames(df_count_a)
  }

  if (conf_int){
    return(list(pval=df_pval, est=df_est, cil=df_cil, cih=df_cih))
  } else {
    return(list(pval=df_pval, est=df_est))
  }
}


#' Compute pvalues of multiple stratified Fisher tests.
#'
#' From 2 matrices of event counts in variables and fixed samples, compute a stratified Fisher test
#' for each 2x2 count table induced by the each pair of rows in table A and table B.
#'
#' @param df_count_a a \code{data.frame}
#' @param df_count_b a \code{data.frame}
#' @param margins_a a numeric vector defining the sample size of each column of \code{df_count_a}
#' @param margins_b a numeric vector defining the sample size of each column of \code{df_count_b}
#' @param n_cores number of cores to be used for parallel computations.
#' @return a \code{data.frame} with the same row names and column names as \code{df_count_a} and \code{df_count_b}.
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom stats mantelhaen.test
#'
#' @author Yoann Pradat
#' @keywords internal
get_pvals_cmh <- function(df_count_a, df_count_b, margins_a, margins_b, n_cores=1){
  stopifnot(identical(rownames(df_count_a), rownames(df_count_b)))
  stopifnot(identical(colnames(df_count_a), colnames(df_count_b)))
  stopifnot(identical(colnames(df_count_a), rownames(margins_a)))
  stopifnot(identical(colnames(df_count_b), rownames(margins_b)))

  # can only compare across common strats
  common_strats <- (margins_a > 0) & (margins_b > 0)
  dfc_count_a <- df_count_a[,common_strats,drop=F]
  dfc_count_b <- df_count_b[,common_strats,drop=F]
  marginsc_a <- margins_a[common_strats,,drop=F]
  marginsc_b <- margins_b[common_strats,,drop=F]

  n_row <- nrow(dfc_count_a)
  n_strat <- ncol(dfc_count_a)
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  pvals <- foreach(i = 1:n_row, .packages=c("Exact"), .combine=c) %dopar% {
    c_a <- dfc_count_a[i,]
    c_b <- dfc_count_b[i,]
    
    if (sum(c_a)==0 & sum(c_b)==0){
      pval <- -999
    } else {
      x  <- sapply(1:n_strat, function(j) 
                   {matrix(c(c_a[j], marginsc_a[j,] - c_a[j], c_b[j], marginsc_b[j,] - c_b[j]), byrow=T, nrow=2)},
                   simplify="array")

      table2x2 <- function(j){
        c(c_a[j], c_b[j], marginsc_a[j,] - c_a[j],  marginsc_b[j,] - c_b[j])
      }

      x <- array(Reduce(c, sapply(1:n_strat, table2x2)), dim=c(2,2,n_strat))

      # NOTE: more work is needed here.
      # The exact CMH test is a conditional test that assumes fixed row and column margins.
      # Is the asymptotic CMH test (chi-squared distributed statistic) a conditional test?
      # Let us choose exact=F for now.
      #
      # correct=T applies Yates' correction.
      #
      # References https://doi.org/10.3390/stats1010008 and https://onlinelibrary.wiley.com/doi/abs/10.1111/anzs.12215
      # may help in choosing a test that is conditional only on fixed row margins.
      out <- mantelhaen.test(x=x, alternative="two.sided", correct=T, exact=F)

      pval <- out$p.value
    }

    pval
  }

  stopCluster(cl)

  df_pvals <- data.frame("p.value"=pvals, row.names=rownames(df_count_a))
  df_pvals
}


#' Add missing rows and columns.
#'
#' From a table and a list of row names and column names, add rows or columns with all 0s if they are not already in the
#' table. After missing rows or columns were added, both rows and columns are sorted.
#'
#' @param df a \code{data.frame}
#' @param row_names a character vector
#' @param col_names a character vector
#' @param row_col column to be used as rownames for \code{df}.
#' @return a \code{data.frame} with missing rows and columns added.
#'
#' @importFrom tibble column_to_rownames
#'
#' @author Yoann Pradat
#' @keywords internal
prepare_count_table <- function(df, row_names, col_names, row_col){
  df <- df %>% column_to_rownames(var=row_col)
  row_names_missing <- setdiff(row_names, rownames(df))
  for (row_name_missing in row_names_missing){
    df_row <- data.frame(lapply(1:ncol(df), function(x){0}))
    rownames(df_row) <- row_name_missing 
    colnames(df_row) <- colnames(df)
    df <- rbind(df, df_row)
  }
  
  col_names_missing <- setdiff(col_names, colnames(df))
  for (col_name_missing in col_names_missing){
    df[,col_name_missing] <- 0
  }

  df[order(rownames(df)),order(colnames(df)),drop=F]
}


#' Add missing row names in margins.
#'
#' After missing rows were added, rows are sorted.
#'
#' @param df a \code{data.frame}
#' @param row_names
#'
#' @author Yoann Pradat
#' @keywords internal
prepare_margins <- function(df, row_names){
  row_names_missing <- setdiff(row_names, rownames(df))
  for (row_name_missing in row_names_missing){
    df[row_name_missing,] <- 0
  }

  df[order(rownames(df)),,drop=F]
}


#' Compute pvalues of multiple 2x2 tables.
#'
#' From a list of 2 data.frames and their margins, build 2 tables with identical rows and columns and run Fisher tests
#' on each of the 2x2 contigency table built from each pair of cells of both tables.
#'
#' @param dfs a 2-level list of data.frames. \code{cohort_a} and \code{cohort_b} must be 1st level names while 2nd-level
#'   names must contain 'count' (the count matrix) and 'count_col' (the vector of columns margins).
#' @param cohort_a Name of the cohort A.
#' @param cohort_b Name of the cohort B.
#' @param row_col Name of the column to be used as rownames.
#' @param row_names Values of \code{row_col} to run tests on.
#' @param col_names Subcohort names of each of cohort A and B to run tests on.
#' @param test Name of the statistical test employed. Use "fisher" for Fisher-Booschloo, and "cmh" for CMH stratified
#'   fisher test.
#' @param n_cores Number of cores to be used for parallel computations.
#'
#' @author Yoann Pradat
#' @export
compute_pvals_table <- function(dfs, cohort_a, cohort_b, row_col, row_names, col_names, test="fisher", n_cores=1){
  df_count_a <- prepare_count_table(df=dfs[[cohort_a]]$count, row_col=row_col, row_names=row_names, col_names=col_names)
  df_count_b <- prepare_count_table(df=dfs[[cohort_b]]$count, row_col=row_col, row_names=row_names, col_names=col_names)
  margins_a <- prepare_margins(df=dfs[[cohort_a]]$count_col, row_names=col_names)
  margins_b <- prepare_margins(df=dfs[[cohort_b]]$count_col, row_names=col_names)

  if (test=="fisher"){
    get_pvals_fisher(df_count_a, df_count_b, margins_a, margins_b, n_cores=n_cores)$pval
  } else if (test=="cmh"){
    get_pvals_cmh(df_count_a, df_count_b, margins_a, margins_b, n_cores=n_cores)
  } else {
    stop(paste0("Unsupported value '", test, "' for 'test' argument."))
  }
}


#' Perform multiple testing correction from a set of pvalues
#'
#' From a list of multiple tables of pvalues, perform one multiple testing comparison correction considering all tests of 
#' all tables at the same time.  Corrected pvalues (qvalues) are returned in the same format as the input list of
#' pvalue tables.
#'
#' @param dfs_pvals a list of dataframes. Names must contain the values in \code{cohorts}.
#' @param cohorts A character vector.
#' @param method The multiple testing correction method.
#' @return a list of dataframes. Names and format are identical to that of \code{dfs_pvals}.
#'
#' @importFrom stats p.adjust
#'
#' @author Yoann Pradat
#' @export
correct_multiple_testing <- function(dfs_pvals, cohorts, method="fdr"){
  # aggregate test p-values to be corrected together
  pvals <- c()
  for (cohort in cohorts){
    pvals <- c(pvals, as.vector(as.matrix(dfs_pvals[[cohort]])))
  }

  # do not correct where no test was performed
  mask_notest <- pvals==-999
  qvals <- abs(pvals)

  # apply correction
  qvals[!mask_notest] <- p.adjust(qvals[!mask_notest], method=method)
  qvals[mask_notest] <- 1

  # restore sign
  qvals[pvals < 0] <- -qvals[pvals < 0]

  # reshape
  dfs_qvals <- list()
  shift <- 1
  for (cohort in cohorts){
    n_row <- nrow(dfs_pvals[[cohort]])
    n_col <- ncol(dfs_pvals[[cohort]])
    size <- n_row*n_col
    df_qvals <- as.data.frame(matrix(qvals[shift:(shift+size-1)], nrow=n_row, byrow=F))
    rownames(df_qvals) <- rownames(dfs_pvals[[cohort]])
    colnames(df_qvals) <- colnames(dfs_pvals[[cohort]])

    shift <- shift + size
    dfs_qvals[[cohort]] <- df_qvals
  }

  dfs_qvals
}


#' Perform multiple testing correction from a set of pvalues
#'
#' From a list of multiple tables of pvalues, perform one multiple testing comparison correction considering all tests of 
#' all tables at the same time.  Corrected pvalues (qvalues) are returned in the same format as the input list of
#' pvalue tables.
#'
#' @param dfs_pvals a list of dataframes. Names must contain the values in \code{cohorts}.
#' @param cohorts A character vector.
#' @param method The multiple testing correction method.
#' @return a list of dataframes. Names and format are identical to that of \code{dfs_pvals}.
#'
#' @importFrom stats p.adjust
#'
#' @author Yoann Pradat
#' @export
add_pvals_tables <- function(dfs_plot, cohorts_a, cohorts_b, test, col_evt, in_plot_evt, tt_keep, n_cores, suffix=""){
  stopifnot(length(cohorts_a)==length(cohorts_b))
  names_comp <- paste(cohorts_a, cohorts_b, sep="_vs_")
  n_tests <- length(cohorts_a)
  dfs_pvals_evt <- list()

  for (i in 1:n_tests){
    cohort_a <- cohorts_a[i]
    cohort_b <- cohorts_b[i]
    name_comp <- names_comp[i]
    dfs_pvals_evt[[name_comp]] <- compute_pvals_table(dfs_plot, cohort_a, cohort_b, row_col=col_evt, test=test,
                                                      row_names=in_plot_evt, col_names=tt_keep, n_cores=n_cores)
  }

  # correct for multiple testing
  dfs_qvals_evt <- correct_multiple_testing(dfs_pvals=dfs_pvals_evt, cohorts=names_comp)

  pvals_name <- paste0("pvals", suffix)
  qvals_name <- paste0("qvals", suffix)

  for (name_comp in names_comp){
    # transform -999 (no test) to 1 (not significant)
    dfs_pvals_evt[[name_comp]][dfs_pvals_evt[[name_comp]] == -999] <- 1
    dfs_qvals_evt[[name_comp]][dfs_qvals_evt[[name_comp]] == -999] <- 1

    dfs_plot[[name_comp]][[pvals_name]] <- as_tibble(dfs_pvals_evt[[name_comp]] %>% rownames_to_column(var=col_evt))
    dfs_plot[[name_comp]][[qvals_name]] <- as_tibble(dfs_qvals_evt[[name_comp]] %>% rownames_to_column(var=col_evt))
  }

  dfs_plot
}
