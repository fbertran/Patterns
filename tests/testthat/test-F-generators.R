blk <- function(arr, i, j, ngrp) {
  sqF <- dim(arr)[1]
  arr[1:sqF, 1:sqF, (i-1)*ngrp + j, drop = FALSE]
}

test_that("CascadeFshape and CascadeFinit return correct dimensions", {
  sqF <- 3L; ngrp <- 2L
  Fs <- CascadeFshape(sqF, ngrp)
  Fi <- CascadeFinit(sqF, ngrp)
  expect_identical(dim(Fs), c(sqF, sqF, ngrp*ngrp))
  expect_identical(dim(Fi), c(sqF, sqF, ngrp*ngrp))
  # diagonal blocks are "0" for shape and 0 for init
  expect_true(all(blk(Fs,1,1,ngrp) == "0"))
  expect_true(all(blk(Fi,1,1,ngrp) == 0))
  expect_true(all(blk(Fs,2,2,ngrp) == "0"))
  expect_true(all(blk(Fi,2,2,ngrp) == 0))
})

test_that("IndicFshape / IndicFinit respect Indic mask", {
  sqF <- 2L; ngrp <- 2L
  Indic <- matrix(1, ngrp, ngrp); Indic[1,1] <- 0
  Fs <- IndicFshape(sqF, ngrp, Indic)
  Fi <- IndicFinit(sqF, ngrp, Indic)
  # (1,1) block must be zeroed
  expect_true(all(blk(Fs,1,1,ngrp) == "0"))
  expect_true(all(blk(Fi,1,1,ngrp) == 0))
})
