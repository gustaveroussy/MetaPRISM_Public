test_that("merge with rownames", {
  # no divergent info
  df_x <- data.frame("A"=c(0,1), "B"=c("J","K"), "C"=c(-1, -4), row.names=c("A", "B"))
  df_y <- data.frame("A"=c(0,3), "C"=c(-1, 4), "D"=c("hello", "you"), row.names=c("A", "C"))

  cols_x <- colnames(df_x)
  cols_y <- colnames(df_y)
  cols_i <- intersect(cols_x, cols_y)

  # inner col
  df_m <- merge_on_rows(df_x, df_y, how_cols="inner")
  expect_true(setequal(colnames(df_m), cols_i))

  # x col
  df_m <- merge_on_rows(df_x, df_y, how_cols="x")
  expect_true(setequal(colnames(df_m), cols_x))

  # y col
  df_m <- merge_on_rows(df_x, df_y, how_cols="y")
  expect_true(setequal(colnames(df_m), cols_y))

  # divergent info (first row of each df)
  df_x <- data.frame("A"=c(0,1), "B"=c("J","K"), "C"=c(-1, -4), row.names=c("A", "B"))
  df_y <- data.frame("A"=c(2,3), "C"=c(-1, 4), "D"=c("hello", "you"), row.names=c("A", "C"))

  tryCatch(
    df_m <- merge_on_rows(df_x, df_y, how_cols="inner"),
    error = function(e){
      expect_true(is(e, "simpleError"))
    }
  )
})
