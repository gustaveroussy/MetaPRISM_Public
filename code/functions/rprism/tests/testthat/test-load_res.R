# GENE LIST

test_that("load_resource() works", {
  skip_on_cran()

  # OncoKB
  cat("\n")
  df_gene <- load_resource(database="oncokb", name="gene")
  expect_true(is(df_gene, "data.frame"))


  # CIViC
  cat("\n")
  df_gene <- load_resource(database="civic", name="gene")
  expect_true(is(df_gene, "data.frame"))

  # curated
  cat("\n")
  df_gene <- load_resource(database="curated", name="gene")
  expect_true(is(df_gene, "data.frame"))

  # CGC
  cat("\n")
  df_gene <- load_resource(database="cosmic", name="gene")
  expect_true(is(df_gene, "data.frame"))

  # Gencode
  df_gene <- load_resource(database="gencode", name="gencode_v27")
  expect_true(is(df_gene, "data.frame"))

})
