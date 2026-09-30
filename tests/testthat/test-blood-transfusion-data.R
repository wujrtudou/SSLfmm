test_that("blood_transfusion package data have the expected structure", {
  e <- new.env(parent = emptyenv())

  utils::data(
    "blood_transfusion",
    package = "SSLfmm",
    envir = e
  )

  expect_true(
    exists("blood_transfusion", envir = e, inherits = FALSE)
  )

  d <- get(
    "blood_transfusion",
    envir = e,
    inherits = FALSE
  )

  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), 748L)

  expect_identical(
    names(d),
    c(
      "id",
      "truth",
      "observed",
      "missing_indicator",
      "Recency",
      "Frequency",
      "Time"
    )
  )

  expect_equal(
    sum(is.na(d$observed)),
    166L
  )

  src <- as.character(d$missing_indicator)

  expect_equal(
    sum(src == "observed"),
    582L
  )

  expect_equal(
    sum(src == "mar"),
    91L
  )

  expect_equal(
    sum(src == "mcar"),
    75L
  )

  expect_true(
    all(
      is.na(d$observed) ==
        (src != "observed")
    )
  )

  expect_false(
    anyNA(
      d[c("Recency", "Frequency", "Time")]
    )
  )
})
