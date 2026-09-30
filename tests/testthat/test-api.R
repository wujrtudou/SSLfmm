test_that("classification metrics are coherent", {
  m <- classification_performance(c(1, 1, 2, 2), c(1, 2, 2, 2))
  expect_equal(unname(m$metrics["accuracy"]), .75)
  expect_equal(sum(m$confusion_matrix), 4)
})


test_that("macro metrics retain a true class that is never predicted", {
  truth <- c(rep(1, 53), rep(2, 22))
  predicted <- rep(1, length(truth))

  m <- classification_performance(truth, predicted)

  expect_equal(unname(m$metrics["accuracy"]), 53 / 75)
  expect_equal(unname(m$metrics["balanced_accuracy"]), 0.5)
  expect_equal(unname(m$metrics["macro_recall"]), 0.5)
  expect_equal(unname(m$metrics["macro_precision"]), (53 / 75) / 2)
  expect_equal(unname(m$metrics["macro_f1"]), 0.4140625)

  class2 <- m$per_class[m$per_class$class == "2", , drop = FALSE]
  expect_equal(class2$precision, 0)
  expect_equal(class2$recall, 0)
  expect_equal(class2$f1, 0)
})


test_that("public API is intentionally small", {
  expected <- c("classification_performance", "fit_sslfmm", "initialize_sslfmm",
                "rmix", "simulate_mixed_missingness",
                "simulate_sslfmm")
  expect_setequal(getNamespaceExports("SSLfmm"), expected)
})
