test_that("load_rna_gex() met500 works", {
  skip_on_cran()

  cat("\n")
  df_rna <- load_rna_gex(study="met500", level="genes", metric="TPM", test_mode=TRUE)
  expect_true(is(df_rna, "data.frame"))

  cat("\n")
  df_rna <- load_rna_gex(study="met500",  level="genes", metric="counts", test_mode=TRUE,
                         identifiers=c("MO_1031", "MO_1107"), identifiers_name="Subject_Id")
  expect_true(is(df_rna, "data.frame"))

})


test_that("load_rna_fus() met500 works", {
  skip_on_cran()

  cat("\n")
  df_rna <- load_rna_fus(study="met500", mode="arriba")
  expect_true(is(df_rna, "data.frame"))

  cat("\n")
  df_rna <- load_rna_fus(study="met500", mode="pizzly")
  expect_true(is(df_rna, "data.frame"))
})


test_that("load_rna_gex() prism works", {
  skip_on_cran()

  cat("\n")
  df_rna <- load_rna_gex(study="prism",  level="genes", metric="TPM", test_mode=TRUE)
  expect_true(is(df_rna, "data.frame"))

  cat("\n")
  df_rna <- load_rna_gex(study="prism",  level="genes", metric="counts", test_mode=TRUE,
                         identifiers=c("201410666GH", "200905404EN"), identifiers_name="Subject_Id")
  expect_true(is(df_rna, "data.frame"))
})


test_that("load_rna_fus() prism works", {
  skip_on_cran()

  cat("\n")
  df_rna <- load_rna_fus(study="prism", mode="ericscript")
  expect_true(is(df_rna, "data.frame"))

  cat("\n")
  df_rna <- load_rna_fus(study="prism", mode="starfusion")
  expect_true(is(df_rna, "data.frame"))
})


test_that("load_rna_gex() tcga works", {
  skip_on_cran()

  cat("\n")
  df_rna <- load_rna_gex(study="tcga",  metric="TPM", test_mode=TRUE)
  expect_true(is(df_rna, "data.frame"))

  cat("\n")
  df_rna <- load_rna_gex(study="tcga",  metric="counts", test_mode=TRUE)
  expect_true(is(df_rna, "data.frame"))

  cat("\n")
  df_rna <- load_rna_gex(study="tcga",  metric="TPM", test_mode=TRUE,
                     identifiers=c("TCGA-02-0047-01A-01R-1849-01", "TCGA-3X-AAV9-01A-72R-A41I-07"),
                     identifiers_name="Sample_Id")
  expect_true(is(df_rna, "data.frame"))
  expect_equal(ncol(df_rna), 3)
})


test_that("load_rna_fus() tcga works", {
  skip_on_cran()

  cat("\n")
  df_rna <- load_rna_fus(study="tcga", mode="aggregated")
  expect_true(is(df_rna, "data.frame"))
})


test_that("load_rna_fus() tcga_6_samples works", {
  skip_on_cran()

  cat("\n")
  df_rna <- load_rna_fus(study="tcga_6_samples", mode="fusioncatcher")
  expect_true(is(df_rna, "data.frame"))
})
