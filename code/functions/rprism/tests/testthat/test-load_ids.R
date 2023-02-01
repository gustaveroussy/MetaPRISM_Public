# MET500

test_that("load_ids() met500 all works", {
  skip_on_cran()

  cat("\n")
  df_ids <- load_ids(study="met500")
  expect_true(is(df_ids, "data.frame"))

})


# PRISM

test_that("load_ids() prism all works", {
  skip_on_cran()

  cat("\n")
  df_ids <- load_ids(study="prism")
  expect_true(is(df_ids, "data.frame"))

})


# TCGA

test_that("load_ids() tcga all works", {
  skip_on_cran()

  cat("\n")
  df_ids <- load_ids(study="tcga")
  expect_true(is(df_ids, "data.frame"))

})
