## TOP -----
library(DLMtool); library(dplyr)


OM <- new("OM", Stock=Albacore, 
          Fleet=Generic_FlatE, 
          Obs=Perfect_Info, 
          Imp=Perfect_Imp)
OM@qinc <- c(0,0)
OM@nsim <- 3

# Test without bio-economic model - no BioEco parameters
# - effort control - curE - OK 
# - TAC - AvC - OK
# - spatial or size limit - matlenlim - OK
# - mixed - to test still

Plot <- function(MSE, mm=1, sim=1) {
  
  par(mfrow=c(3,2))
  
  plot(MSE@SSB[sim, mm, ]/MSE@OM$SSB0[sim], xlab="Year", ylab="Depletion", type="l",
       bty="l", lwd=2, ylim=c(0, 1))
  abline(h=MSE@OM$SSBMSY_SSB0[sim], lty=2)
  
  profit <- MSE@Misc$Revenue[sim, mm,] - MSE@Misc$Cost[sim, mm,]
  maxPM <- max(abs(profit))
  ylim <- c(-maxPM, maxPM)
  plot(profit, xlab="Year", ylab="Profit", type="l",
       bty="l", lwd=2, ylim=ylim)
  abline(h=0, lty=2)
  
  
  plot(MSE@Misc$Cost[sim, mm,], xlab="Year", ylab="Cost (solid) Revenue (dashed)", type="l",
       bty="l", lwd=2, ylim=c(0, max(c(MSE@Misc$Cost[sim, mm,],MSE@Misc$Revenue[sim, mm,]))))
  # abline(h=1)
  
  lines(MSE@Misc$Revenue[sim, mm,], lwd=2, lty=2)
  
  # 
  # plot(MSE@Effort[sim, mm,], xlab="Year", ylab="Effort", type="l",
  #      bty="l", lwd=2, ylim=c(0, max(MSE@Effort[sim, mm,])))
  # abline(h=)
  # 
  # plot(MSE@C[sim, mm,], xlab="Year", ylab="Catch", type="l",
  #      bty="l", lwd=2, ylim=c(0, max(MSE@C[sim, mm,])))
  
  
  plot(profit, MSE@Effort[sim,mm,],
       type="l", xlab='Profit', ylab="Effort",
       bty="l")
  
  # plot(profit, MSE@C[sim,mm,],
  #      type="l", xlab='Profit', ylab="Catch",
  #      bty="l")
  
  plot(MSE@C[sim,mm,],MSE@Effort[sim,mm,],
       type="l", xlab='Catch', ylab="Effort",
       bty="l")
  
  # y1 <- 1:(OM@proyears-1)
  # y2 <- y1+1
  # ChangeEffort <- (MSE@Effort[sim,mm,y2]/ MSE@Effort[sim,mm,y1] - 1)
  # 
  # plot(ChangeEffort, MSE@Misc$PMargin[sim, mm,y1], type="l", 
  #      bty="n", xlab="Change in Effort",
  #      ylab="Profit Margin")
  
  
}


# with a bio-economic model
# OM <- tinyErr(OM)
# OM@Perr <- c(0.1,0.1)
OM@proyears <- 100
OM@h <- c(0.7, 0.7)
OM@D <- c(0.7, 0.7)
OM@RevCurr <- c(1,1)
OM@CostCurr <- c(1,1)
OM@Response <- c(0.1, 0.1) # add checks
OM@RevInc <- OM@CostInc <- c(0, 0)# add checks if first two populated 

MSE <- runMSE(OM, MPs=c('Itarget1', 'NMref', 'AvC'))

sim <- 1
Plot(MSE, 1, sim=sim)
Plot(MSE, 2, sim=sim)
Plot(MSE, 3, sim=sim)

par(mfrow=c(1,2))
matplot(t(MSE@B_BMSY[sim,,]), type="l") # >% round(2) %>% t()
matplot(t(MSE@C[sim,,]), type="l")


MSE@C[sim,,] %>% round(2) %>% t()


sim <- 1; mm <-2 
MSE@SSB[sim, mm, ]/MSE@OM$SSB0[sim]



data.frame(Profit=MSE@Misc$Revenue[sim,1,] - MSE@Misc$Cost[sim,1,],
           Effort=MSE@Effort[sim,1,]) %>% round(2)


DF <- data.frame(Cost =MSE@Misc$Cost[sim,1,],
           Revenue = MSE@Misc$Revenue[sim,1,],
           Profit = MSE@Misc$Revenue[sim,1,] - MSE@Misc$Cost[sim,1,],
           PM = 1 - MSE@Misc$Cost[sim,1,]/ MSE@Misc$Revenue[sim,1,])
DF %>% round(2)

par(mfrow=c(1,2))
plot(DF$Profit, DF$Cost, type="b")

plot(DF$PM, DF$Cost, type="b")


# Check Response calculations and explanation 
# write up as profit and profit margin and make equivalent

# units of today's effort 
0.5 * -.6

MSE@Effort[1,1, ] %>% round(2)

# effort by % change of each year

# profit v profit margin


y1 <- 1:(OM@proyears-1)
y2 <- y1+1
ChangeEffort <- (MSE@Effort[sim,mm,y2]/ MSE@Effort[sim,mm,y1] - 1)


# bio-economic model with existing TAE


# with just an existing TAE - latent effort
OM@LatentEff <- c(0.05,0.05)
MSE2 <- runMSE(OM, MPs=MPs)

matplot(t(MSE2@C[1,,]),type="l")
matplot(t(MSE2@Effort[1,,]),type="l")
# fix - curE should be actual not latent
# - NMref doesn't work


# Bio-economic model


OM@CostCurr <- c(1,1)
OM@RevCurr <- c(0.95,1.05)

OM@CostInc <- c(0,0)
OM@RevInc <- c(0,0)

OM@Response <- c(0.05,0.05)
OM@LatentEff<- c(0.3, 0.4)

OM@nsim <- 5
OM@interval <- 1

# No TAC or TAE - cost = revenue - depleted

nothing <- function(x, Data, ...) {
  new("Rec")
}
class(nothing) <- "MP"

MPs <- c('nothing', 'matlenlim', "ITM")

OM@LatentEff <- numeric(0)
OM@CostCurr <- c(1,1)
OM@RevCurr <- c(0.7,0.75)

OM@Response <- c(0.05, 0.05)

OM@h <- c(0.7, 0.7)
OM@D <- c(0.6, 0.6)
# OM <- tinyErr(OM)
OM@proyears <- 200
OM@nsim <- 2
OM@Vmaxlen <- c(1,1)
# OM@qcv <- c(0,0)
OM@qinc <- c(0,0)
MSE <- runMSE(OM, MPs=MPs)



sim <- 2; mm <- 3





plot(MSE@B_BMSY[sim, mm, ], xlab="Year", ylab="B/BMSY", type="l",
     bty="l", lwd=2, ylim=c(0, max(max(MSE@B_BMSY[sim, mm, ]), 1.5)))
abline(h=1, lty=2)

# plot(MSE@C[sim, mm, ]/MSE@OM$RefY[sim], xlab="Year", ylab="Catch", type="l",
#      bty="l", lwd=2, ylim=c(0, 1.1))

# plot(MSE@Effort[sim, mm,], xlab="Year", ylab="Effort", type="l",
     # bty="l", lwd=2, ylim=c(0, max(max(MSE@Effort[sim, mm,]),1) ))






# 
# plot(MSE@Misc$Revenue[sim, mm,]-MSE@Misc$Cost[sim, mm,], xlab="Year", ylab="Profit", type="l",
#      bty="l", lwd=2)




DF <- data.frame(TAC=MSE@TAC[sim,mm,], Catch=MSE@C[sim,mm,])

DF %>% round(2) %>%
  head(20)










cbind((OM@Response[1] * MSE@Misc$PMargin[sim, mm,y1])/MSE@Effort[sim,mm,y1],
ChangeEffort)


plot(MSE@SSB[sim, mm, ]/MSE@OM$SSB0[sim], MSE@C[sim, mm, ]/MSE@OM$RefY[sim],
     xlab="SSB/SSB0", ylab="Catch", type="l", bty="l", lwd=2,
     xlim=c(0, 1.1), ylim=c(0,1))


plot(MSE@C[sim, mm, ]/MSE@OM$RefY[sim], xlab="Year", ylab="Catch", type="l",
     bty="l", lwd=2, ylim=c(0, 1.1))




plot(MSE@Misc$PMargin[1,mm,], MSE@Effort[1,mm,], type="b")
plot(MSE@C[1,mm,], MSE@Effort[1,mm,], type="b")
plot(MSE@Misc$PMargin[1,mm,], MSE@C[1,mm,], type="b")
plot(MSE@TAC[1,mm,], MSE@C[1,mm,], type="b")





# TAC with no existing TAE 

OM@D <- c(0.1, 0.1)
OM <- tinyErr(OM)
MSE <- runMSE(OM, MPs=c("AvC", "Itarget1", 'matlenlim'))

mm <- 3
sim <- 1 
plot(MSE@Misc$PMargin[sim,mm,], MSE@Effort[sim,mm,], type="b")
plot(MSE@TAC[sim,mm,], MSE@C[sim,mm,], type="b")


# TAC with an existing TAE 
MSE@Misc$LatEffort[sim,mm,]








def.args <- DLMtool:::dev.mode(); for (nm in names(def.args)) assign(nm, def.args[[nm]])





