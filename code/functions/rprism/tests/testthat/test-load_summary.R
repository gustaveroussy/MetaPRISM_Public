
# RNA

test_that("load_summary_rna_fus() met500 all works", {
  skip_on_cran()

  df_summary <- load_summary_rna_fus(study="met500", mode="arriba")
  expect_true(is(df_summary, "data.frame"))

  df_summary <- load_summary_rna_fus(study="met500", mode="pizzly")
  expect_true(is(df_summary, "data.frame"))
})


test_that("load_summary_rna_gex() met500 all works", {
  skip_on_cran()

  df_summary <- load_summary_rna_gex(study="met500", level="genes", metric="counts")
  expect_true(is(df_summary, "data.frame"))

  df_summary <- load_summary_rna_gex(study="met500", level="genes", metric="TPM")
  expect_true(is(df_summary, "data.frame"))
})


test_that("load_summary_rna_gex() prism all works", {
  skip_on_cran()

  df_summary <- load_summary_rna_gex(study="prism", level="genes", metric="counts")
  expect_true(is(df_summary, "data.frame"))

  df_summary <- load_summary_rna_gex(study="prism", level="genes", metric="TPM")
  expect_true(is(df_summary, "data.frame"))

  df_summary <- load_summary_rna_gex(study="prism", level="transcripts", metric="counts")
  expect_true(is(df_summary, "data.frame"))

  df_summary <- load_summary_rna_gex(study="prism", level="transcripts", metric="TPM")
  expect_true(is(df_summary, "data.frame"))
})


test_that("load_summary_rna_fus() prism all works", {
  skip_on_cran()

  df_summary <- load_summary_rna_fus(study="prism", mode="ericscript")
  expect_true(is(df_summary, "data.frame"))

  df_summary <- load_summary_rna_fus(study="prism", mode="starfusion")
  expect_true(is(df_summary, "data.frame"))
})


test_that("load_summary_rna_gex() tcga all works", {
  skip_on_cran()

  df_summary <- load_summary_rna_gex(study="tcga", level="genes", metric="counts")
  expect_true(is(df_summary, "data.frame"))

  df_summary <- load_summary_rna_gex(study="tcga", level="genes", metric="TPM")
  expect_true(is(df_summary, "data.frame"))
})

# WES

test_that("load_summary_wes_mut() met500 all works", {
  skip_on_cran()

  df_summary <- load_summary_wes_mut(study="met500", mode="oncotator_unfiltered")
  expect_true(is(df_summary, "data.frame"))

  df_summary <- load_summary_wes_mut(study="met500", mode="oncotator_filtered")
  expect_true(is(df_summary, "data.frame"))

  df_summary <- load_summary_wes_mut(study="met500", mode="annovar_unfiltered")
  expect_true(is(df_summary, "data.frame"))

  df_summary <- load_summary_wes_mut(study="met500", mode="annovar_pathogenic")
  expect_true(is(df_summary, "data.frame"))
})


test_that("load_summary_wes_mut() prism all works", {
  skip_on_cran()

  df_summary <- load_summary_wes_mut(study="prism", mode="oncotator_unfiltered")
  expect_true(is(df_summary, "data.frame"))

  df_summary <- load_summary_wes_mut(study="prism", mode="oncotator_filtered")
  expect_true(is(df_summary, "data.frame"))

  df_summary <- load_summary_wes_mut(study="prism", mode="annovar_unfiltered")
  expect_true(is(df_summary, "data.frame"))

  df_summary <- load_summary_wes_mut(study="prism", mode="annovar_pathogenic")
  expect_true(is(df_summary, "data.frame"))
})


test_that("load_summary_wes_mut() tcga all works", {
  skip_on_cran()

  df_summary <- load_summary_wes_mut(study="tcga", mode="mc3_filtered")
  expect_true(is(df_summary, "data.frame"))
})
