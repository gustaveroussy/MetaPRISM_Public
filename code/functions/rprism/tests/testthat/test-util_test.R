test_that("get_pvals_fisher works", {
  skip_on_cran()

  col_names <- c("ACC","BLCA", "BRCA")
  row_names <- c("AKT1", "ALK", "AR")
  df_count_a <- data.frame(matrix(c(0,0,2,0,0,0,0,0,0), nrow=3, byrow=T), row.names=row_names)
  df_count_b <- data.frame(matrix(c(0,1,7,0,0,0,0,0,0), nrow=3, byrow=T), row.names=row_names)
  colnames(df_count_a) <- col_names
  colnames(df_count_b) <- col_names
  margins_a <- data.frame(Count=c(12, 45, 61), row.names=col_names)
  margins_b <- data.frame(Count=c(91, 411, 1064), row.names=col_names)

  out <- get_pvals_fisher(df_count_a, df_count_b, margins_a, margins_b, conf_int=F, conf_lvl=0.95, progress=F, n_cores=1)
  expect_true(is(out$pval, "data.frame"))
  expect_true(is(out$est, "data.frame"))

  out <- get_pvals_fisher(df_count_a, df_count_b, margins_a, margins_b, conf_int=F, conf_lvl=0.95, progress=F, n_cores=3)
  expect_true(is(out$pval, "data.frame"))
  expect_true(is(out$est, "data.frame"))
})
