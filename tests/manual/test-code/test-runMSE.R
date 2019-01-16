
library(DLMtool)
DLMtool::setup()

Ntest <- 12 # set to NA to run all tests 
nsim <- 6

testthat::context("Test runMSE works with built-in objects with MP of each type")

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




# SRrel = 2 
for (x in 1:Ntest) {
  OM <- new("OM", get(all[x,1]), get(all[x,2]), get(all[x,3]), get(all[x,4]))
  OM@seed <- ceiling(runif(1, 1, 1000))
  OM@nsim <- nsim
  OM@interval <- ceiling(runif(1, 1, 5))
  OM@SRrel <- 2
  info <- paste(OM@Name, "seed =", OM@seed, "interval =", OM@interval)
  testthat::test_that(paste0("runMSE with Ricker SRR: ",info), {
    testthat::expect_is(runMSE(OM, MPs=MPs, parallel=FALSE, silent=TRUE), 'MSE', info=info)
  })
}



# Parallel 
# SRrel = 2 
for (x in 1:Ntest) {
  OM <- new("OM", get(all[x,1]), get(all[x,2]), get(all[x,3]), get(all[x,4]))
  OM@seed <- ceiling(runif(1, 1, 1000))
  OM@nsim <- 48
  OM@interval <- ceiling(runif(1, 1, 5))
  OM@SRrel <- 2
  info <- paste(OM@Name, "seed =", OM@seed, "interval =", OM@interval)
  testthat::test_that(paste0("runMSE with parallel: ",info), {
    testthat::expect_is(runMSE(OM, MPs=MPs, parallel=TRUE, silent=TRUE), 'MSE', info=info)
  })
}

# CheckMPs works and run all MPs - BH SRR
OM <- new("OM", get(all[1,1]), get(all[1,2]), get(all[1,3]), get(all[1,4]))
OM@seed <- ceiling(runif(1, 1, 1000))
OM@nsim <- 6
OM@interval <- ceiling(runif(1, 1, 5))
OM@SRrel <- 1
info <- paste(OM@Name, "seed =", OM@seed, "interval =", OM@interval)
testthat::test_that(paste0("runMSE with all MPs and BH SRR: ",info), {
  testthat::expect_is(runMSE(OM, MPs=NA, parallel=FALSE, silent=TRUE), 'MSE', info=info)
})


# CheckMPs works and run all MPs - Ricker SRR
OM <- new("OM", get(all[1,1]), get(all[1,2]), get(all[1,3]), get(all[1,4]))
OM@seed <- ceiling(runif(1, 1, 1000))
OM@nsim <- 6
OM@interval <- ceiling(runif(1, 1, 5))
OM@SRrel <- 2
info <- paste(OM@Name, "seed =", OM@seed, "interval =", OM@interval)
testthat::test_that(paste0("runMSE with all MPs and Ricker SRR: ",info), {
  testthat::expect_is(runMSE(OM, MPs=NA, parallel=FALSE, silent=TRUE), 'MSE', info=info)
})


# Check runMSE works in parallel with historical simulations 
OM@nsim <- 288
info <- paste(OM@Name, "seed =", OM@seed)
testthat::test_that(paste0("runMSE with Hist=TRUE and parallel=TRUE: ",info), {
  testthat::expect_is(runMSE(OM, Hist=TRUE, parallel = TRUE), 'Hist', info=info)
})



# # Historical MPA 
# for (x in 1:Ntest) {
#   OM <- new("OM", get(all[x,1]), get(all[x,2]), get(all[x,3]), get(all[x,4]))
#   OM@seed <- ceiling(runif(1, 1, 1000))
#   OM@nsim <- 6
#   OM@interval <- ceiling(runif(1, 1, 5))
#   OM2 <- OM 
#   yrs <- 2:OM2@nyears
#   OM2@MPA <- matrix(c(yrs, rep(0,length(yrs)), rep(1, length(yrs))), ncol=3)
#   mse1 <- runMSE(OM, MPs="MRreal", silent=TRUE)
#   mse2 <- runMSE(OM2, MPs="MRreal", silent=TRUE)
#   info <- paste(OM@Name, "seed =", OM@seed, "interval =", OM@interval)
#   
#   pof1 <- round(summary(mse1, silent=TRUE)$POF,1)
#   pof2 <- round(summary(mse2, silent=TRUE)$POF,1)
#   
#   testthat::test_that(paste0("Historical MPA decrease POF: ",info), {
#     testthat::expect_true(pof1 >= pof2)
#   })
# }



