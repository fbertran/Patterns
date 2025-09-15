test_that("all datasets load", {
  ds <- data(package = "Patterns")$results[, "Item"]
  expect_true(length(ds) > 0)
  for (nm in ds) {
    expect_silent(data(list = nm, package = "Patterns"))
    obj <- get(nm, envir = .GlobalEnv, inherits = TRUE)
    expect_false(is.null(obj))
    # clean up
    rm(list = nm, envir = .GlobalEnv)
  }
})
