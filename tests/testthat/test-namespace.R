test_that("exported objects exist", {
  exports <- getNamespaceExports("Patterns")
  ns <- asNamespace("Patterns")
  expect_true(length(exports) > 0)
  for (x in setdiff(exports,"show")) {
    expect_true(exists(x, envir = ns, inherits = FALSE), info = x)
  }
})

test_that("exported S4 classes exist", {
  for (cl in c("omics_array","omics_network","omics_predict")) {
    expect_true(isClass(cl), info = cl)
  }
})


        