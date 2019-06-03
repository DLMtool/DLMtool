## TOP -----
library(DLMtool); library(dplyr)
library(ggplot2); library(tidyr)
library(cowplot)

OM <- new("OM", Stock=Albacore, 
          Fleet=Generic_FlatE, 
          Obs=Perfect_Info, 
          Imp=Perfect_Imp)
OM@qinc <- c(0,0)
OM@nsim <- 3
OM@proyears <- 100
OM@h <- c(0.7, 0.7)
OM@D <- c(0.4, 0.4)


OM2 <- OM 
OM2@cpars$RevCurr <- runif(OM@nsim, 1.1, 1.3)
OM2@cpars$CostCurr <- runif(OM@nsim,0.95, 1.05)
OM2@cpars$Response <- runif(OM@nsim, 0.1, 0.1)
OM2@cpars$CostInc <- runif(OM@nsim, 2,3)  
OM2@cpars$RevInc <- runif(OM@nsim, 2,3)

# OM2@cpars$LatentEff <- runif(OM@nsim, 0.2, 0.2) # get rid of

Eff1 <- function(x, Data, ...) {
  # Rec <- Itarget1(x, Data, ...)
  Rec@Effort <- 1.5
  Rec
}
class(Eff1) <- "MP"

MPs <- c('NMref', 'Itarget1', "ItargetE1", 'curE', "Eff1")


MSE1 <- runMSE(OM, MPs=MPs)
MSE2 <- runMSE(OM2, MPs=MPs)

sim <- 2

MSEobj <- MSE2

Depletion <- as.vector(MSEobj@SSB[sim, , ])/MSEobj@OM$SSB0[sim]
BMSY_B0 <- MSEobj@OM$SSBMSY_SSB0[sim]
MPs <- rep(MSEobj@MPs, MSEobj@proyears)
Years <- rep(1:MSEobj@proyears, each=MSEobj@nMPs)

Catch <- as.vector(MSEobj@C[sim, , ])/MSEobj@OM$RefY[sim]
TAC <- as.vector(MSEobj@TAC[sim,,])/MSEobj@OM$RefY[sim]

Effort <- as.vector(MSEobj@Effort[sim, , ])
TAE <- as.vector(MSEobj@Misc$TAE[sim,,])

Cost <- as.vector(MSEobj@Misc$Cost[sim, , ])
Revenue <- as.vector(MSEobj@Misc$Revenue[sim, , ])
Profit <- Revenue - Cost

DF <- data.frame(MP=MPs, Year=Years, Depletion=Depletion, BMSY_B0=BMSY_B0,
                 Catch=Catch, TAC=TAC,
                 Effort=Effort, TAE=TAE,
                 Cost=Cost, Revenue=Revenue, Profit=Profit)

MP1 <- MSEobj@MPs[1]
MP1 <- MSEobj@MPs[2]
MP1 <- MSEobj@MPs[3]
MP1 <- MSEobj@MPs[4]
MP1 <- MSEobj@MPs[5]

pDF <- DF %>% filter(MP==MP1)
LineSize <- 1 
p1 <- ggplot(pDF, aes(x=Year, y=Depletion)) + geom_line(size=LineSize) +
  theme_classic() + expand_limits(y=c(0,1)) +
  geom_abline(slope=0, intercept = BMSY_B0, lty=2)

p2 <- ggplot(pDF, aes(x=Year, y=Catch)) + geom_line(size=LineSize) +
  theme_classic() + expand_limits(y=0) +
  geom_line(aes(x=Year, y=TAC))

p3 <- ggplot(pDF, aes(x=Year, y=Effort)) + geom_line(size=LineSize) +
  theme_classic() + expand_limits(y=0) +
  geom_line(aes(x=Year, y=TAE))

pDF2 <- pDF %>% gather("key", "value", 9:10)

p4 <- ggplot(pDF2, aes(x=Year, y=value, color=key)) + geom_line(size=LineSize) +
  geom_line(aes(x=Year, y=Profit, color="Profit"), size=LineSize) + 
  theme_classic() + expand_limits(y=0) +
  geom_abline(slope=0, intercept = 0, lty=2)

p <- plot_grid( p1, p2, p3, p4, ncol=2, labels="auto")

title <- ggdraw() + draw_label(MP1, fontface='bold')
plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1))

# plot Catch & Revenue together
# plot Effort & Cost together

# Summarize performance in economic terms 
# Net present value 
# discounting ...



DF2 <- DF %>% tidyr::gather('key', 'value', c(3,5,7,9,10,11))

head(DF2)
DF3 <- DF2 %>% filter(MP=="curE")


ggplot(DF3, aes(x=Year, y=value)) + facet_wrap(~key, scales="free", ncol=6) +
  geom_line()



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





