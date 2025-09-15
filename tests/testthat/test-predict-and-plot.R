test_that("predict() on a tiny omics_array returns omics_predict and plots", {
  set.seed(3)
  # Build a 6-gene, 3-time, 2-subject omics_array with 3 groups
  time <- 1:3
  subject <- 2
  G <- 6
  M <- matrix(rnorm(G * length(time) * subject), nrow = G)
  groups <- rep(1:3, each = 2)
  OA <- as.omics_array(M, time = time, subject = subject, group = groups)

  # Tiny network with 6 nodes and three time labels
  nb <- 6
  time_label <- c(1,1,2,2,3,3)
  Net <- network_random(
    nb = nb,
    time_label = time_label,
    exp = 1, init = 1,
    regul = rep(1, nb),
    min_expr = 0.1, max_expr = 0.2,
    casc.level = 0.3
  )

  # Activation times per group (3 groups)
  act_time_group <- c(1,2,3)

  OP <- predict(OA, Net, act_time_group = act_time_group, nv = 0.5)
  expect_s4_class(OP, "omics_predict")

  # Check inner omics arrays exist and have shape
  expect_s4_class(OP@omicsarray_predict, "omics_array")
  expect_identical(dim(OP@omicsarray_predict@omicsarray), dim(OA@omicsarray))
  expect_s4_class(OP@omicsarray_changed, "omics_array")
  expect_s4_class(OP@omicsarray_unchanged, "omics_array")
})
