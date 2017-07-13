
setup()

Ntest <- NA # set to NA to run all tests 
nsim <- 9

context("Test runMSE works with built-in objects with MP of each type")

MPs <- c("AvC", "matlenlim", "slotlim", "curE", "curE75", "MRnoreal", "MRreal")

stocks <- avail("Stock")
fleets <- avail("Fleet")
obs <- avail("Obs")
imps <- avail("Imp")

all <- as.matrix(expand.grid(stocks, fleets, obs, imps))

if (is.na(Ntest)) Ntest <- nrow(all)

all <- all[sample(1:nrow(all), size=Ntest),]


test_that("runMSE works with all built-in Stock, Fleet, Obs, & Imp objects", {
  for (x in 1:Ntest) {
    OM <- new("OM", get(all[x,1]), get(all[x,2]), get(all[x,3]), get(all[x,4]))
    OM@seed <- ceiling(runif(1, 1, 1000))
    OM@nsim <- nsim
    interval <- ceiling(runif(1, 1, 5))
    info <- paste(OM@Name, "seed =", OM@seed, "interval =", interval)
    expect_error(runMSEnomsg(OM, MPs=MPs, interval=interval), NA, info=info)
  }
})

  
  
 