make_net <- function(M) {
  new("omics_network",
      omics_network = M,
      name = paste0("g", seq_len(nrow(M))),
      F = array(0, c(1,1,1)),
      convF = matrix(0,1,1),
      convO = 0,
      time_pt = 1)
}

make_net2 <- function(M) {
  new("omics_network",
      omics_network = M,
      name = paste0("g", seq_len(nrow(M))),
      F = array(0, c(1,1,1)),
      convF = matrix(0,1,1),
      convO = 0,
      time_pt = c(1,2))
}

test_that("compare() reports perfect metrics on identical networks", {
  M <- matrix(0, 4, 4)
  M[1,2] <- 0.9; M[3,4] <- 0.7
  N1 <- make_net(M)
  N2 <- make_net(M)
  crits <- Patterns::compare(N1, N2, nv = 0.5)
  expect_type(crits, "double")
  expect_equal(names(crits),
               c("Sensitivity","PPV","F-score","F-score 1/2","F-score 2"))
  expect_equal(unname(crits), rep(1, 5), tolerance = 1e-12)
})

test_that("cutoff() returns named components on tiny network", {
  M <- matrix(0, 5, 5)
  set.seed(4)
  M[upper.tri(M)] <- runif(sum(upper.tri(M)))
  N <- make_net(M)
  co <- suppressWarnings(cutoff(N, sequence = seq(0.05, 0.15, by = 0.05)))
  expect_true(is.list(co))
  expect_true(all(c("p.value", "p.value.inter", "sequence") %in% names(co)))
  expect_length(co$p.value, length(co$sequence))
})

