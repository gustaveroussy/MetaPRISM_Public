test_that("load_dsg() prism works", {
  skip_on_cran()

  cat("\n")
  df_design <- load_dsg(study="prism")
  expect_true(is(df_design, "data.frame"))
})

