## TOP -----
library(DLMtool); library(dplyr)


OM <- new("OM", Stock=Albacore, 
          Fleet=Generic_Fleet, 
          Obs=Perfect_Info, 
          Imp=Perfect_Imp)

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

par(mfrow=c(3,2))

plot(MSE@SSB[sim, mm, ]/MSE@OM$SSB0[sim], xlab="Year", ylab="Depletion", type="l",
     bty="l", lwd=2, ylim=c(0, 1))
abline(h=MSE@OM$SSBMSY_SSB0[sim], lty=2)

plot(MSE@B_BMSY[sim, mm, ], xlab="Year", ylab="B/BMSY", type="l",
     bty="l", lwd=2, ylim=c(0, max(max(MSE@B_BMSY[sim, mm, ]), 1.5)))
abline(h=1, lty=2)

# plot(MSE@C[sim, mm, ]/MSE@OM$RefY[sim], xlab="Year", ylab="Catch", type="l",
#      bty="l", lwd=2, ylim=c(0, 1.1))

# plot(MSE@Effort[sim, mm,], xlab="Year", ylab="Effort", type="l",
     # bty="l", lwd=2, ylim=c(0, max(max(MSE@Effort[sim, mm,]),1) ))


plot(MSE@Misc$Cost[sim, mm,], xlab="Year", ylab="Cost", type="l",
     bty="l", lwd=2, ylim=c(0, max(MSE@Misc$Cost[sim, mm,])))

plot(MSE@Misc$Revenue[sim, mm,], xlab="Year", ylab="Revenue", type="l",
     bty="l", lwd=2, ylim=c(0, max(MSE@Misc$Revenue[sim, mm,])))


plot(MSE@Misc$Revenue[sim, mm,]-MSE@Misc$Cost[sim, mm,], xlab="Year", ylab="Profit", type="l",
     bty="l", lwd=2)


maxPM <- min(1.5,max(abs(MSE@Misc$PMargin[sim, mm,])))
ylim <- c(-maxPM, maxPM)
plot(MSE@Misc$PMargin[sim, mm,], xlab="Year", ylab="Profit Margin", type="l",
     bty="l", lwd=2, ylim=ylim)
abline(h=0, lty=2)

DF <- data.frame(TAC=MSE@TAC[sim,mm,], Catch=MSE@C[sim,mm,])

DF %>% round(2) %>%
  head(20)








# Check Response calculations and explanation 
y1 <- 1:(OM@proyears-1)
y2 <- y1+1
ChangeEffort <- (MSE@Effort[sim,mm,y2]/ MSE@Effort[sim,mm,y1] - 1)

par(mfrow=c(1,1))
plot(ChangeEffort, MSE@Misc$PMargin[sim, mm,y1], type="l", 
     ylim=c(-1, 1),
     xlim=c(0,OM@Response[1]*2))

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





