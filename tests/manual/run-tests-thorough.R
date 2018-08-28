library(DLMtool)

DLMtool::setup()

Ntest <- NA # set to NA to run all tests 
nsim <- 3

testthat::context("Test runMSE with all MPs and all built-in objects")

MPs <- avail("MP") # NA # c("AvC", "matlenlim", "slotlim", "curE", "curE75", "MRnoreal", "MRreal")

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
  for (MP in MPs) {
    info <- paste(OM@Name, "seed =", OM@seed, "interval =", OM@interval, "MP = ", MP)
    testthat::test_that(paste0("runMSE: ",info), {
      testthat::expect_is(runMSE(OM, MPs=MP, parallel=FALSE, silent=TRUE), 'MSE', info=info)
    })
  }

}
