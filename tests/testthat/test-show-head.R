test_that("show() and head() work on omics_array and omics_network", {
  time <- 1:2; subject <- 2; G <- 4
  M <- matrix(seq_len(G * length(time) * subject), nrow = G)
  OA <- as.omics_array(M, time = time, subject = subject)

  # show/print do not error
  expect_silent(capture.output(show(OA)))
  h <- head(OA)
  expect_true(is.list(h))
  expect_true(all(c("omicsarray","name","gene_ID","time","subject","group","start_time") %in% names(h)))

  # omics_network show
  Mnet <- matrix(0, G, G); Mnet[1,2] <- 0.5
  N <- new("omics_network",
           omics_network = Mnet,
           name = paste0("g", 1:G),
           F = array(0, c(1,1,1)),
           convF = matrix(0,1,1),
           convO = 0,
           time_pt = 1)
  expect_silent(capture.output(show(N)))
})
