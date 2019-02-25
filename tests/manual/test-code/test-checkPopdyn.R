
testthat::context("Quick test of pop dyn")

library(DLMtool)

OM <- tinyErr(testOM, silent = TRUE)
OM@Vmaxlen <- c(1,1)
OM@SRrel <- 1

testpopdyn <- function(OM) {
  Hist <- runMSE(OM, Hist=TRUE, silent=TRUE)
  Data <- Hist@Data
  set.seed(OM@seed)
  # Stock Parameters & assign to function environment
  StockPars <- SampleStockPars(OM, OM@nsim, OM@nyears, OM@proyears, list(), msg=FALSE)
  
  # Fleet Parameters & assign to function environment
  FleetPars <- SampleFleetPars(SubOM(OM, "Fleet"), Stock=StockPars, OM@nsim, 
                               OM@nyears, OM@proyears, cpars=list())
  ObsPars <- SampleObsPars(OM, nsim, cpars=list())
  for (X in 1:length(ObsPars)) assign(names(ObsPars)[X], ObsPars[[X]])
  
  dnormal<-function(lens,lfs,sl,sr){
    cond<-lens<=lfs
    sel<-rep(NA,length(lens))
    sel[cond]<-2.0^-((lens[cond]-lfs)/sl*(lens[cond]-lfs)/sl)
    sel[!cond]<-2.0^-((lens[!cond]-lfs)/sr*(lens[!cond]-lfs)/sr)
    sel
  }
  BHSRR <- function(SBcurr, SB0, R0, steepness) {
    (4 * R0 * steepness * SBcurr)/(SB0/R0 * R0 * (1-steepness) + (5*steepness-1)*SBcurr)
  }
  
  simpop <- function(logapicF, OM, FleetPars, StockPars, sim=1, opt=1) {
    for (X in 1:length(StockPars)) assign(names(StockPars)[X], StockPars[[X]])
    for (X in 1:length(FleetPars)) assign(names(FleetPars)[X], FleetPars[[X]])
    
    Len <- Len_age[sim,,1]
    Wght <- a * Len^b
    MAA <- 1/(1 + exp(-log(19) * ((Len - L50[sim])/(L95[sim]-L50[sim])))) 
    
    # selectivity-at-length - fishery
    sl <- (LFS[1,sim] - L5[1,sim]) /((-log(0.05,2))^0.5)
    sr <- (Linf[sim] - LFS[1,sim]) / ((-log(Vmaxlen[1,sim],2))^0.5) # selectivity parameters are constant for all years
    SAA <- dnormal(Len, LFS[1,sim], sl, sr)
    
    M_array <- rep(M[sim],maxage)
    FAA <- exp(logapicF) * SAA
    ZAA <- (FAA * Find[sim,1]) + M_array 
    ages <- 1:maxage
    N <- VB <- matrix(NA, OM@nyears, maxage)
    Rec <- SB <- rep(NA, OM@nyears)
    N[1,1] <- R0[sim] 
    N[1,2:maxage] <- R0[sim] * exp(-cumsum(M_array[2:maxage]))
    
    SB[1] <- sum(N[1,] * MAA * Wght)
    SB0 <- sum(SB[1])
    Rec[1] <- R0[sim]
    for (yr in 2:OM@nyears) {
      Rec[yr] <- BHSRR(SB[yr-1], SB0, R0[sim], steepness=hs[sim]) 
      N[yr,1] <- Rec[yr] 
      ZAA <- (FAA * Find[sim,yr-1]) + M_array 
      N[yr,2:maxage] <- N[yr-1, 1:(maxage-1)] * exp(-ZAA[1:(maxage-1)])
      SB[yr] <- sum(N[yr,]* MAA * Wght)
    }
    dep <- SB[yr]/SB0
    
    if (opt==1) return((dep-StockPars$D[sim])^2)
    return(list(N=N, SAA=SAA, LenCV=LenCV, Len=Len))
  }
  chk <- rep(NA, OM@nsim)
  for (sim in 1:OM@nsim) {
    opt <- optimize(simpop, interval=log(c(0.01, 0.9)), OM=OM, FleetPars=FleetPars, 
                    StockPars=StockPars, sim=sim, opt=1)
    simple <- simpop(opt$minimum, OM=OM, FleetPars=FleetPars, 
                     StockPars=StockPars, sim=sim, opt=2)
    
    chk[sim] <- prod(round(Hist@AtAge$Nage[sim,,]/t(simple$N) ,0))
  }
  chk
}

testthat::expect_equal(testpopdyn(OM), rep(1, OM@nsim))

