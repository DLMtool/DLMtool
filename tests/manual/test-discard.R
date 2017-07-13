setup()

nsim <- 9

context("Test discard mortality")


nsim <- 5 
OM <- updateMSE(testOM)

OM@nsim <- nsim
OM@isRel <- 'FALSE'
OM@L50 
OM@Linf 

# Base case 
OM@L5 <- c(30, 45)
OM@LFS <- c(80, 85)
OM@Vmaxlen <- c(1, 1)
plotSelect(OM,1)

# retention curve 
OM2 <- OM 
OM2@LR5 <- c(90,95)
OM2@LFR <- c(100, 100)
plotSelect(OM2, 1)

# discard mortality 
OM3 <- OM2
OM3@Fdisc <- c(0.2, 0.2)
plotSelect(OM3,1)

# flat discarding rate across all sizes 
OM4 <- OM3
OM4@DR <- c(0.2, 0.2)
plotSelect(OM4,1)


# assume all returned fish die 
OM5 <- OM4
OM5@Fdisc <- c(1,1)
plotSelect(OM5,1)


OM6 <- OM3 
OM6@Fdisc <- c(1,1)


# Test MSE 
OM1 <- testOM
# OM1 <- makePerf(testOM)

OM1@nsim  <- 5
# OM1@D <- c(0.2, 0.2)

OM1@isRel <- 'FALSE'
OM1@L50 <- c(80, 80)
OM1@Linf <- c(130, 130)
OM1@t0 <- c(0,0)
# vulnerability 
OM1@L5 <- c(10, 10)
OM1@LFS <- c(20, 20)
OM1@Vmaxlen <- c(1,1)

# retention 
OM1@LR5 <- c(90, 90)
OM1@LFR <- c(91, 91)

m1 <- runMSE(OM1, interval=1)

OM2 <- OM1 
OM2@Fdisc <- c(1,1)
m2 <- runMSE(OM2, interval=1)

summary(m1)
summary(m2)

Tplot(m1)
Tplot(m2)



