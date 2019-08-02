context("Test MSE_functions")

DLMextra(TRUE)
library(DLMextra)
library(testthat)

Obj <- updateMSE(DLMextra::testMSE)

testthat::test_that("Converge", {
  testthat::expect_error(Converge(Obj), NA)
})


testthat::test_that("checkMSE", {
  testthat::expect_error(checkMSE(Obj), NA)
})


testthat::test_that("Sub by MP", {
  MPs1 <- Obj@MPs[1:3]
  MPs2 <- Obj@MPs[4:Obj@nMPs]
  testthat::expect_error(t1 <- Sub(Obj, MPs=MPs1), NA)
  testthat::expect_error(t2 <- Sub(Obj, MPs=MPs2), NA)
})

testthat::test_that("Sub by sim", {
  nsim <- Obj@nsim
  sims1 <- 1:ceiling((nsim/2))
  sims2 <- (max(sims1)+1):nsim
  testthat::expect_error(t1 <<- Sub(Obj, sim=sims1), NA)
  testthat::expect_error(t2 <<- Sub(Obj, sim=sims2), NA)
})

testthat::test_that("joinMSE", {
  testthat::expect_error(newMSE <<- joinMSE(list(t1, t2)), NA)
})

testthat::test_that("joinMSE returns same object", {
  testthat::expect_true(all(summary(newMSE, silent=TRUE) == summary(Obj, silent=TRUE)))
})


# testthat::test_that("DOM", {
#   testthat::expect_error(DOM(Obj), NA)
# })

