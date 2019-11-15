library(DLMtool)

OM <- testOM 
OM@proyears <- 200
OM@nsim <- 20
OM@qinc <- c(0,0)
OM@qcv <- c(0,0)
Hist <- runMSE(OM, Hist=TRUE)

Catch <- Hist@TSdata$Catch
Hist_Effort <- Hist@TSdata$Find # not exactly right, but approx
Hist_Effort <- Hist_Effort/Hist_Effort[,OM@nyears]

rho <- PPC <- CPE <- NA
lpars <- log(c(0.01, 10, 10))

for (sim in 1:OM@nsim) {
  optfun <- function(lpars, Hist_Effort, Catch, sim) {
    rho <- exp(lpars[1])
    PPC <- exp(lpars[2])
    CPE <- exp(lpars[3])
    ssq <- NULL
    predE <<- rep(Hist_Effort[sim,1], length(Hist_Effort[sim,]))
    for (x in 1:(ncol(Hist_Effort)-1)) {
      # predE[x+1] <- Hist_Effort[sim, x] + rho*(Catch[sim,x]*PPC - Hist_Effort[sim, x]*CPE)
      predE[x+1] <- predE[x] + rho*(Catch[sim,x]*PPC - Hist_Effort[sim, x]*CPE)
      ssq[x] <- (predE[x] - Hist_Effort[sim,x+1])^2
    }
    # predE <<- predE/predE[length(predE)]
    # ssq <- (predE - Hist_Effort[sim,])^2
    sum(ssq)
  }
  
  doOpt <- optim(lpars, optfun, Hist_Effort=Hist_Effort, Catch=Catch, sim=sim)
  
  pars<-exp(doOpt$par)
  rho[sim] <- pars[1]
  PPC[sim] <- pars[2]
  CPE[sim] <- pars[3]
  
}


OM@cpars$RevCurr <- Catch[,50]*PPC
OM@cpars$CostCurr <- Hist_Effort[,50]*CPE
OM@cpars$Response <- rho

NMref <- function(x, Data, ...) {
  new("Rec")
}
class(NMref) <- "MP"


MSE <- runMSE(OM, MPs=c('curE', 'NMref'))

cur <- cbind(Hist_Effort, MSE@Effort[,1,])
sq <- cbind(Hist_Effort, MSE@Effort[,2,])

sim <- 1


# Plot Catch, Effort, and Profit
mm <- 2
totC <- c(Catch[sim,], MSE@C[sim, mm, ])
eff <- cbind(Hist_Effort, MSE@Effort[,mm,])
Cost <- eff[sim,] * CPE[sim]
Rev <- totC * PPC[sim]
Profit <- Rev - Cost 

par(mfrow=c(1,2))
plot(totC, type="l")

plot(eff[sim,], type="l", ylim=c(0, 2))
par(new=TRUE)
plot(Profit, type="l", ylim=c(-200, 200), col="blue")
abline(h=0, lty=2)
axis(side=4)


t1 <- eff[sim,51:(OM@nyears+OM@proyears)]
deltaE <- t1[2:length(t1)]-t1[1:(length(t1)-1)]
plot(Profit[51:(length(Profit)-1)], deltaE)
abline(h=0, lty=2)
abline(v=0, lty=2)

plot(Profit, type="l")
abline(h=0, lty=2)

par(mfrow=c(1,2))
plot(Profit, type="l", col="blue")
abline(h=0, lty=2)
plot(eff[sim,]/max(eff[sim,]), col="blue", type="l")


0.98 + rho[sim]*(Profit[51])



predE2 <- rep(Hist_Effort[sim,1], length(Hist_Effort[sim,]))
for (i in 2:length(Hist_Effort[sim,])) {
  predE2[i] <- predE2[i-1] + rho2 * (Catch[sim, i-1]*PPC2 - Hist_Effort[sim,i-1]*CPE2)
  # predE2[i] <- Hist_Effort[sim,i-1] + rho2 * (Catch[sim, i-1]*PPC2 - Hist_Effort[sim,i-1]*CPE2)
}

# plot Effort and predicted Effort 
plot(Hist_Effort[sim,])
lines(predE)
lines(predE2, col="blue")

  
  
  
0.7 + rho*(0)









matplot(t(Hist_Effort), type="l")
matplot(t(Catch), type="l")

Data <- SimulatedData

OM@cpars$Data <- Data
OM@nsim <- 4

MSE <- runMSE(OM)




simcatch <- apply(CBret, c(1,3), sum)
sim <- 3

plot(simcatch[sim,], col="blue")
lines(Data@Cat[1,], col="red")
