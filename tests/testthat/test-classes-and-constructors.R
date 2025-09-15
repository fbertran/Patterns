test_that("as.omics_array builds a valid object", {
  set.seed(1)
  time <- 1:4
  subject <- 6
  G <- 10
  M <- matrix(rnorm(G * length(time) * subject), nrow = G)
  OA <- as.omics_array(M, time = time, subject = subject)
  expect_s4_class(OA, "omics_array")
  expect_true(nrow(OA@omicsarray) == G)
  expect_true(length(OA@time) == length(time))
  expect_true(OA@subject == subject)
})

test_that("network_random returns an omics_network of expected size", {
  nb <- 6L
  time_label <- c(1,1,2,2,3,3)
  exp <- 1
  init <- 1
  regul <- rep(1, nb)
  min_expr <- 0.1
  max_expr <- 0.2
  casc.level <- 0.3
  Net <- network_random(nb, time_label, exp, init, regul, min_expr, max_expr, casc.level)
  expect_s4_class(Net, "omics_network")
  expect_true(is.matrix(Net@omics_network))
  expect_identical(dim(Net@omics_network), c(nb, nb))
})
