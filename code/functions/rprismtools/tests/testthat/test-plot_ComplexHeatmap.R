test_that("draw_upset_plot works", {
  df <- data.frame(Group=sample(c("A","B","C","D"), size=100, replace=T),
                   Index=sample(seq(1,20), size=100, replace=T))

  m <- make_upset_m(df, field_set="Group", field_identifier="Index")
  filepath <- "./out_tests/test_draw_upset_plot.pdf"

  pdf(file=filepath, width=6, height=3)
  draw_upset_plot(m, width_set_size=2, height_top_annot=2)
  dev.off()

  expect_true(file.exists(filepath))
})

