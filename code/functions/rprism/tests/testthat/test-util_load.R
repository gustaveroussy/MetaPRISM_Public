
test_that("load_from_data works", {
  skip_on_cran()

  df_bio <- load_from_data("./tcga/clinical/raw_bio_files/ICD-O-3_morphology_table.xlsx")
  expect_true(is(df_bio, "data.frame"))

  df_bio <- load_from_data("./prism/clinical/raw_cln_files/general/Anapath_python_trim.txt", delim="\t")
  expect_true(is(df_bio, "data.frame"))

  df_bio <- load_from_data("./tcga/clinical/raw_cln_files/tcga/TCGA-gdc_cases_days.csv", delim=";")
  expect_true(is(df_bio, "data.frame"))
})
