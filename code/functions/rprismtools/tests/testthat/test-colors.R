test_that("get_color_labels works", {
  labels <- c("A", "A", "B", "C", "D", "B")
  colors <- get_label_colors(labels)

  expect_identical(colors[1], colors[2])
  expect_identical(colors[3], colors[6])
})

test_that("rect_plot_colors works", {
  colors <- list(A="blue", B="green", C="red")

  filepath <- "./out_tests/test_rect_plot_colors_1.pdf"

  pdf(file=filepath, width=6, height=3)
  rect_plot_colors(col=1, colors=colors)
  dev.off()

  expect_true(file.exists(filepath))

  filepath <- "./out_tests/test_rect_plot_colors_2.pdf"

  pdf(file=filepath, width=6, height=3)
  rect_plot_colors(col=2, colors=colors)
  dev.off()

  expect_true(file.exists(filepath))
})
