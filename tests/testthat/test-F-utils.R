test_that("replaceBand / replaceUp / replaceDown behave as expected", {
  a <- matrix(1:9, 3, 3)
  b <- matrix(0, 3, 3)

  expect_equal(replaceBand(a, b, 0),
               {res <- a; diag(res) <- 0; res})

  expect_equal(replaceBand(a, b, 1),
               {res <- a; res[abs(row(a)-col(a))<=1] <- 0; res})

  expect_equal(replaceDown(a, matrix(1,3,3), 0),
               {res <- a; res[lower.tri(res, diag = TRUE)] <- 1; res})

  expect_equal(replaceUp(a, matrix(1,3,3), 0),
               {res <- a; res[upper.tri(res, diag = TRUE)] <- 1; res})

  expect_equal(replaceUp(a, b, 1),
               {res <- a; res[(row(a)-col(a))<=1] <- 0; res})
})
