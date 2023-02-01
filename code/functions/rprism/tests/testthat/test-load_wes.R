test_that("load_wes_mut() met500 works", {
  skip_on_cran()
  cat("\n")
  df_maf <- load_wes_mut(study="met500", mode="somatic_maf")
  expect_true(is(df_maf, "data.frame"))

  cat("\n")
  df_maf <- load_wes_mut(study="met500", mode="somatic_filters")
  expect_true(is(df_maf, "data.frame"))
})

test_that("load_wes_mut() prism works", {
  skip_on_cran()
  cat("\n")
  df_maf <- load_wes_mut(study="prism", mode="somatic_maf")
  expect_true(is(df_maf, "data.frame"))

  cat("\n")
  df_maf <- load_wes_mut(study="prism", mode="somatic_filters")
  expect_true(is(df_maf, "data.frame"))
})

test_that("load_wes_mut() prism works", {
  skip_on_cran()

  cat("\n")
  df_maf <- load_wes_mut(study="tcga", mode="somatic_civic")
  expect_true(is(df_maf, "data.frame"))

  cat("\n")
  df_maf <- load_wes_mut(study="tcga", mode="somatic_oncokb")
  expect_true(is(df_maf, "data.frame"))
})
