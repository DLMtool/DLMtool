
context("runMSE with all available OMs and all MPs")

GetMoreData(TRUE)
library(DLMdata)
setup()

MPs <- c(avail("Input"), avail("Output"))
OMs <- avail("OM")

test_that("runMSE works with all OMs for all MP classes with interval = 3", {
  for (om in seq_along(OMs)) {
    OM <- get(OMs[om])
    OM@seed <- ceiling(runif(1, 1, 1000))
    expect_error(runMSE(OM, MPs=MPs, interval=3), NA, info=paste(om, "seed = ", OM@seed))
  }
})


test_that("runMSE works with all OMs for all MP classes with interval = 1 ", {
  for (om in seq_along(OMs)) {
    OM <- get(OMs[om])
    OM@seed <- ceiling(runif(1, 1, 1000))
    expect_error(runMSE(OM, MPs=MPs, interval=1), NA, info=paste(om, "seed = ", OM@seed))
  }
})