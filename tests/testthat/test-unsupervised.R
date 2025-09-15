test_that("unsupervised_clustering_* basic flow on tiny data", {
  set.seed(2)
  time <- 1:3
  subject <- 6
  G <- 100
  M <- matrix(exp(rnorm(G * length(time) * subject)), nrow = G)
  OA <- as.omics_array(M, time = time, subject = subject)
  mc <- unsupervised_clustering_auto_m_c(OA)
  expect_true(is.list(mc))
  expect_true("m" %in% names(mc))
  OA2 <- unsupervised_clustering(OA, 3, mc$m, screen = NULL, heatmap = FALSE, new.window = FALSE)
  expect_s4_class(OA2, "omics_array")
})
