test_that("mcc computes expected confusion matrix and bootstrap summary", {
  skip_if_not_installed("boot")

  x <- c(1, 1, 1, 0, 0, 0, 1, 0)
  y <- c(1, 1, 0, 1, 0, 0, 1, 0)

  set.seed(123)
  res <- mcc(x, y, positive_class = 1, bootstrap_reps = 200, confidence = 0.90)

  expect_equal(unname(res$confusion_matrix), matrix(c(3, 1, 1, 3), nrow = 2, byrow = TRUE))
  expect_equal(res$mcc, (3 * 3 - 1 * 1) / sqrt(4 * 4 * 4 * 4))
  expect_true(is.numeric(res$mcc_bootstrap_pv))
  expect_true(is.finite(res$mcc_bootstrap_pv))
  expect_length(res$mcc_bootstrap_ci, 2)
  expect_true(all(is.finite(res$mcc_bootstrap_ci)))
})

test_that("mcc supports time-series bootstrap with compact encoding", {
  skip_if_not_installed("boot")

  x <- rep(c(1, 1, 0, 0), 10)
  y <- c(x[-1], x[1])

  set.seed(321)
  res <- mcc(x, y, positive_class = 1, ts = TRUE, l = 4, bootstrap_reps = 50)

  expect_type(res, "list")
  expect_equal(dim(res$confusion_matrix), c(2, 2))
  expect_true(is.finite(res$mcc))
  expect_true(is.numeric(res$mcc_bootstrap_pv) || is.na(res$mcc_bootstrap_pv))
  expect_length(res$mcc_bootstrap_ci, 2)
})
