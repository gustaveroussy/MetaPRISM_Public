test_that("load_colors()", {
  skip_on_cran()

  global2colors <- load_colors(sheet="Global")
  expect_true(is.list(global2colors))

  BS2colors <- load_colors(sheet="Biopsy_Site")
  expect_true(is.list(BS2colors))

  df_global_colors <- load_colors(sheet="Global", as_tibble=T)
  expect_true(is(df_global_colors, "data.frame"))
})
