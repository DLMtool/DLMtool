
DLMtool::setup()

Ntest <- NA # set to NA to run all tests 
nsim <- 3

context("Test runMSE works with built-in objects with MP of each type")

MPs <- c("AvC", "matlenlim", "slotlim", "curE", "curE75", "MRnoreal", "MRreal")

stocks <- avail("Stock")
fleets <- avail("Fleet")
obs <- avail("Obs")
imps <- avail("Imp")

all <- as.matrix(expand.grid(stocks, fleets, obs, imps))

if (is.na(Ntest)) Ntest <- nrow(all)

all <- all[sample(1:nrow(all), size=Ntest),]

for (x in 1:Ntest) {
  OM <- new("OM", get(all[x,1]), get(all[x,2]), get(all[x,3]), get(all[x,4]))
  OM@seed <- ceiling(runif(1, 1, 1000))
  OM@nsim <- nsim
  OM@interval <- ceiling(runif(1, 1, 5))
  info <- paste(OM@Name, "seed =", OM@seed, "interval =", OM@interval)
  testthat::test_that(paste0("runMSE: ",info), {
    testthat::expect_is(runMSE(OM, MPs=MPs, parallel=FALSE, silent=TRUE), 'MSE', info=info)
  })
}


