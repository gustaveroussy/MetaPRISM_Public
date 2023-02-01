# MET500

test_that("load_cln() met500 all works", {
  skip_on_cran()

  cat("\n")
  df_cln <- load_cln(study="met500")
  expect_true(is(df_cln, "data.frame"))

  cat("\n")
  df_cln <- load_cln(study="met500", mode="all")
  expect_true(is(df_cln, "data.frame"))
})

test_that("load_cln() met500 Subject_Id works", {
  skip_on_cran()

  cat("\n")
  df_cln <- load_cln(study="met500", identifiers=c("MO_1537", "MO_1196"), identifiers_name="Subject_Id")
  expect_true(is(df_cln, "data.frame"))
})



# PRISM

test_that("load_cln() prism all works", {
  skip_on_cran()

  cat("\n")
  df_cln <- load_cln(study="prism", mode="in_design")
  expect_true(is(df_cln, "data.frame"))
})

test_that("load_cln() prism Subject_Id works", {
  skip_on_cran()

  cat("\n")
  df_cln <- load_cln(study="prism", identifiers=c("201410666GH", "200905404EN"), identifiers_name="Subject_Id")
  expect_true(is(df_cln, "data.frame"))
})


# TCGA

test_that("load_cln() tcga all works", {
  skip_on_cran()

  cat("\n")
  df_cln <- load_cln(study="tcga")
  expect_true(is(df_cln, "data.frame"))
})


test_that("load_cln() tcga Subject_Id works", {
  skip_on_cran()

  cat("\n")
  df_cln <- load_cln(study="tcga", identifiers=c("TCGA-02-0047", "TCGA-CJ-4876"), identifiers_name="Subject_Id")
  expect_true(is(df_cln, "data.frame"))
})
