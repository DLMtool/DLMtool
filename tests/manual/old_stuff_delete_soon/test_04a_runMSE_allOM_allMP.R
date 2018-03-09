
# Test runMSE with all OMs (in both DLMtool and DLMdata) and all MPs 


context("runMSE with all available OMs and all MPs")

DLMextra(TRUE)
library(DLMtool)
library(DLMextra)
library(testthat)


nsim <- 9

MPs <- c('MP')
OMs <- avail("OM")

for (om in seq_along(OMs)) {
  OM <- get(OMs[om])
  OM@seed <- ceiling(runif(1, 1, 1000))
  OM@nsim <- nsim
  OM@interval <- 4
  info <-paste(OM@Name, "seed = ", OM@seed)
  test_that(paste0("runMSE interval = 4:", OM@Name), {
    testthat::expect_is(runMSEnomsg(OM, MPs=MPs, silent=TRUE), 'MSE', info=info)
  })
}

for (om in seq_along(OMs)) {
  OM <- get(OMs[om])
  OM@seed <- ceiling(runif(1, 1, 1000))
  OM@nsim <- nsim
  OM@interval <- 1
  info <-paste(OM@Name, "seed = ", OM@seed)
  test_that(paste0("runMSE interval = 1:", OM@Name), {
    testthat::expect_is(runMSEnomsg(OM, MPs=MPs, silent=TRUE), 'MSE', info=info)
  })
}

