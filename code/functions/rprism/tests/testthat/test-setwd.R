test_that("working directory change works", {
  skip_on_cran()
  current_wd <- setwd_to_data()
  expect_true(grepl("data$", getwd()))
  setwd(current_wd)

  current_wd <- setwd_to_results()
  expect_true(grepl("results$", getwd()))
  setwd(current_wd)

  current_wd <- setwd_to_logs()
  expect_true(grepl("logs$", getwd()))
  setwd(current_wd)
})
