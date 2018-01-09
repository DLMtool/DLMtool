
# Test runMSE with all OMs (in both DLMtool and DLMdata) and all MPs 


context("runMSE with all available OMs and all MPs")

DLMextra(TRUE)
library(DLMtool)
library(DLMextra)
library(testthat)
setup()

nsim <- 9

MPs <- c(avail("Input"), avail("Output"))
OMs <- avail("OM")

test_that("runMSE works with all OMs for all MP classes with interval = 3", {
  for (om in seq_along(OMs)) {
    OM <- get(OMs[om])
    print(OMs[om])
    OM@seed <- ceiling(runif(1, 1, 1000))
    OM@nsim <- nsim
    expect_error(runMSEnomsg(OM, MPs=MPs, interval=3), NA, info=paste(OMs[om], "seed = ", OM@seed))
  }
})


test_that("runMSE works with all OMs for all MP classes with interval = 1 ", {
  for (om in seq_along(OMs)) {
    OM <- get(OMs[om])
    print(OMs[om])
    OM@seed <- ceiling(runif(1, 1, 1000))
    OM@nsim <- nsim
    expect_error(runMSEnomsg(OM, MPs=MPs, interval=1), NA, info=paste(OMs[om], "seed = ", OM@seed))
  }
})


