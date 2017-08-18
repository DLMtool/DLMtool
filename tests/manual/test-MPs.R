
setup()


context("Test all MPs with testOM")


nsim <- 12

MPs <- c(avail("Input"), avail("Output"))
OM <- testOM

test_that("runMSE works with all MPs with interval = 3", {
  for (x in seq_along(MPs)) {
    OM@seed <- ceiling(runif(1, 1, 1000))
    OM@nsim <- nsim
    MP <- MPs[x]
    expect_error(runMSEnomsg(OM, MPs=MP, interval=3), NA, info=paste(MP, "seed = ", OM@seed))
    # expect_warning(runMSEnomsg(OM, MPs=MP, interval=3), NA, info=paste(MP, "seed = ", OM@seed))
  }
})


test_that("runMSE works with all MPs with interval = 3", {
  for (x in seq_along(MPs)) {
    OM@seed <- ceiling(runif(1, 1, 1000))
    OM@nsim <- nsim
    MP <- MPs[x]
    expect_error(runMSEnomsg(OM, MPs=MP, interval=1), NA, info=paste(MP, "seed = ", OM@seed))
    # expect_warning(runMSEnomsg(OM, MPs=MP, interval=1), NA, info=paste(MP, "seed = ", OM@seed))
  }
})

