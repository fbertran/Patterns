test_that("plotF runs silently for Fshape and F", {
  Fs <- CascadeFshape(3, 2)
  Fi <- CascadeFinit(3, 2)
  expect_silent(plotF(Fs, choice = "Fshape"))
  expect_silent(plotF(Fi, choice = "F", nround = 0))
})
