test_that("rare_residuals accepts univariate ts and derives time", {
  skip_if_not_installed("dbscan")
  x <- AirPassengers
  res <- rare_residuals(x, method = "dbscan", dbscan_args = list(minPts = 5))
  expect_type(res, "list")
  expect_true(all(c("data", "seasonality_shift") %in% names(res)))
  d <- res$data
  expect_s3_class(d, "data.frame")
  expect_true(all(c("time", "value", "residual") %in% names(d)))
  expect_equal(nrow(d), length(x))
  expect_true(is.numeric(d$time))
  # The derived time should start and end at the ts time range
  expect_equal(d$time[1], as.numeric(stats::time(x)[1]))
  expect_equal(d$time[nrow(d)], as.numeric(stats::time(x)[length(x)]))
  # Check that seasonality shift list exists (seasonal inferred)
  expect_false(is.null(res$seasonality_shift))
})
