
# Scripts to check for inconsistencies in the internal population dynamics 

library(DLMtool)

# Bag limit dev 

OM <- tinyErr(testOM)
OM@interval <- 1
MSE <- runMSE(OM, MP="FMSYref")
Pplot(MSE)


MPs='LstepCC1'
MSE <- runMSE(MPs=MPs)

DDe

LstepCC1


# Test Effort calculation with fixed TAC 

# set equilbrium conditions 
OM <- testOM 
OM <- tinyErr(OM)
OM@D <- c(0.4, 0.5)
OM@h <- c(0.6, 0.6)
OM@Vmaxlen <- c(1,1)
OM@EffLower <- c(0,1,1,1)
OM@EffUpper <- OM@EffLower
OM@nyears <- 100
OM@nsim <- 5
OM@proyears <- 30
OM@interval <- 10

curCatch <- function(x, Data, ...) {
  Rec <- new("Rec")
  Rec@TAC <- Data@Cat[x, Data@Year == Data@LHYear]
  Rec
}
class(curCatch) <- 'MP'

curCatch_MP <- function(x, Data, ...) {
  Rec <- new("Rec")
  Rec@Spatial <- c(0,1)
  Rec@Allocate <- 1
  Rec@TAC <- Data@Cat[x, Data@Year == Data@LHYear]
  Rec
}
class(curCatch_MP) <- 'MP'

# check effort if Area 1 is closed - should go up with current catch 

OM@Cbiascv <- 0
OM@Cobs <- c(0,0)
OM@qinc <- c(0,0)
OM@qcv <- c(0,0)
OM@Prob_staying <- c(0.5,0.5)
OM@Size_area_1 <- OM@Frac_area_1 <- c(0.6,0.6)
MSE <- runMSE(OM, MPs=c('curCatch', 'curCatch_MP'))

MSE@Effort[1,,]
MSE@C[1,,1:10]


plot(MSE@Effort[1,1,], type="l")

plot(MSE@B_BMSY[1,1,], type="l")
plot(MSE@F_FMSY[1,1,], type="l")


HistCatches <- apply(MSE@CB_hist, c(1,3), sum)
FutureCatch <- MSE@C
HistCatches[,OM@nyears]
FutureCatch[,1,1:5]

y <- 1:10
qinc <- 1
(1/((1 + qinc/100)^y))





AllCatch <- cbind(HistCatches, FutureCatch[,2,])
AllCatch[,OM@nyears:(OM@nyears+1)]

matplot(t(AllCatch), type="l")

Cobs[,OM@nyears]

# Check refY is calculated correctly 

OM <- testOM

# OM <- tinyErr(OM)
OM@D <- c(0.4, 0.5)
OM@h <- c(0.6, 0.6)
OM@Vmaxlen <- c(1,1)
OM@nsim <- 5
OM@proyears <- 100

RefYMP <- function(x, Data, ...) {
  Rec <- new("Rec")
  Rec@TAC <- Data@OM$RefY[x]
  Rec
}
class(RefYMP) <- 'MP'

MSYMP <- function(x, Data, ...) {
  Rec <- new("Rec")
  Rec@TAC <- Data@OM$MSY[x]
  Rec
}
class(MSYMP) <- 'MP'

MSE <- runMSE(OM, MP=c('RefYMP', 'MSYMP', 'FMSYref'))

par(mfrow=c(2,3))
RefC <- MSE@C/MSE@OM$RefY
matplot(t(RefC[,1,]), type="l", ylim=c(0,1))
matplot(t(RefC[,2,]), type="l", ylim=c(0,1))
matplot(t(RefC[,3,]), type="l", ylim=c(0,1))

matplot(t(MSE@B_BMSY[,1,]), type="l", ylim=c(0, 2))
matplot(t(MSE@B_BMSY[,2,]), type="l", ylim=c(0, 2))
matplot(t(MSE@B_BMSY[,3,]), type="l", ylim=c(0, 2))


MSE@OM$RefY
MSE@OM$MSY






MSE <- runMSE(OM, MP='MSYMP')
par(mfrow=c(2,1))
RefC <- MSE@C/MSE@OM$RefY
matplot(t(RefC[,1,]), type="l", ylim=c(0,1))
matplot(t(MSE@B_BMSY[,1,]), type="l", ylim=c(0, 2))