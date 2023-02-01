# MET500

test_that("load_bio() met500 all works", {
  skip_on_cran()

  cat("\n")
  df_bio <- load_bio(study="met500")
  expect_true(is(df_bio, "data.frame"))

  cat("\n")
  df_bio <- load_bio(study="met500", mode="all")
  expect_true(is(df_bio, "data.frame"))
})

test_that("load_bio() met500 Subject_Id works", {
  skip_on_cran()

  cat("\n")
  df_bio <- load_bio(study="met500", identifiers=c("MO_1537", "MO_1196"), identifiers_name="Subject_Id")
  expect_true(is(df_bio, "data.frame"))
})



# PRISM

test_that("load_bio() prism all works", {
  skip_on_cran()

  cat("\n")
  df_bio <- load_bio(study="prism")
  expect_true(is(df_bio, "data.frame"))
})

test_that("load_bio() prism Subject_Id works", {
  skip_on_cran()

  cat("\n")
  df_bio <- load_bio(study="prism", identifiers=c("201410666GH", "200905404EN"), identifiers_name="Subject_Id")
  expect_true(is(df_bio, "data.frame"))
})


# TCGA

test_that("load_bio() tcga all works", {
  skip_on_cran()

  cat("\n")
  df_bio <- load_bio(study="tcga")
  expect_true(is(df_bio, "data.frame"))
})


test_that("load_bio() tcga Subject_Id works", {
  skip_on_cran()

  cat("\n")
  df_bio <- load_bio(study="tcga", identifiers=c("TCGA-02-0047", "TCCA-02-0055"), identifiers_name="Subject_Id")
  expect_true(is(df_bio, "data.frame"))
})
