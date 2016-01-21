###############################################
# Original LBSPR with regeneration assumption #
###############################################
LBSPR <- function(x, DLM_data, yrsmth=5, perc=pstar,reps=reps) {
  dependencies="DLM_data@CAL, DLM_data@CAL_bins, DLM_data@vbLinf, 
	DLM_data@vbK, DLM_data@Mort, DLM_data@L50, DLM_data@L95, DLM_data@wlb, DLM_data@wla,
	DLM_data@Abun, DLM_data@CV_Abun"

  # Save other stuff for smoothing estimates
  TotYears <- nrow(DLM_data@CAL[1,,]) # How many years of length data exist
  if (length(DLM_data@Misc[[x]]) == 0) { # Misc List is empty
    # Create Empty List Object
	MiscList <- rep(list(0), 5) # Create empty list
	MiscList[[1]] <- rep(NA, TotYears) # SPR ests
	MiscList[[2]] <- rep(NA, TotYears) # Smoothed SPR ests
	MiscList[[3]] <- rep(NA, TotYears) # FM ests
	MiscList[[4]] <- rep(NA, TotYears) # Smoothed FM ests
	MiscList[[5]] <- list()
  }
  if (length(DLM_data@Misc[[x]]) != 0) MiscList <- DLM_data@Misc[[x]]
  
  # Add Extra Row when needed 
  if (length(MiscList[[1]]) < TotYears) {
    Diff <- TotYears - length(DLM_data@Misc[[x]][[1]])
    MiscList[[1]] <- append(MiscList[[1]], rep(NA,Diff)) 
	MiscList[[2]] <- append(MiscList[[2]], rep(NA,Diff)) 
	MiscList[[3]] <- append(MiscList[[3]], rep(NA,Diff)) 
	MiscList[[4]] <- append(MiscList[[4]], rep(NA,Diff)) 
  }

  NEmpty <- sum(is.na(MiscList[[1]])) # Number of empty spots
  IsEmpty <- which(is.na(MiscList[[1]]))

 StockPars <- NULL
 StockPars$NGTG <- 41
 StockPars$GTGLinfBy <- NA
 StockPars$Linf <- DLM_data@vbLinf[x]
 StockPars$CVLinf <- 0.1 # NEED TO ADD THIS TO INPUT VARIABLES
 StockPars$MaxSD <- 2 
 StockPars$MK  <- DLM_data@Mort[x] / DLM_data@vbK[x]
 StockPars$L50 <- DLM_data@L50[x] 
 StockPars$L95 <- DLM_data@L95[x] 
 StockPars$Walpha <- DLM_data@wla[x]
 StockPars$Wbeta <- DLM_data@wlb[x]
 StockPars$FecB  <- DLM_data@wlb[x]
 StockPars$Steepness <- 0.9 # DLM_data@steep[x] not required
 StockPars$Mpow <- 0
 StockPars$R0  <- 1000
 StockPars$CentLinf <- NULL
 StockPars$CentMpar <- NULL
 StockPars$CentKpar <- NULL
 StockPars$Mslope <- 0 
 
 # Run Assessment for every year where data exists 
 if (NEmpty !=0) {
   for (xYr in IsEmpty[1]:max(IsEmpty)) {
    
     # Average length data over years if more than *yrsmth*
     if (xYr < yrsmth) ind <- xYr
     if (xYr >= yrsmth) { 
       # ind <- (length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  	 ind <- (xYr-(yrsmth-1)):xYr
     }	 
     if (length(ind) > 1) LenDat <- apply(DLM_data@CAL[x, ind,], 2, mean, na.rm=TRUE)
     if (length(ind) == 1) LenDat <- DLM_data@CAL[x, ind,]
     
     binWidth <- DLM_data@CAL_bins[2] - DLM_data@CAL_bins[1]
     CAL_binsmid <- seq(from=0.5*binWidth, by=binWidth, length=length(DLM_data@CAL_bins)-1)
     
     SizeBins <- NULL
     SizeBins$Linc <- binWidth
     SizeBins$ToSize <- max(CAL_binsmid)
     
     SL50Start <- CAL_binsmid[which.max(LenDat)]
     DeltaStart <- 0.1 * SL50Start
     FMStart <- 1 
     Starts <- log(c(SL50Start/StockPars$Linf, DeltaStart/StockPars$Linf, FMStart))
     Lower <- log(c(0.1, 0.1, 0.001))
     Upper <- log(c(0.9, 0.9, 20))
     # run optimization
     Opt <- nlminb(Starts, OptRoutine, LenDat=LenDat, Stock=StockPars, Mids=CAL_binsmid)
     # need to add penalty for selectivity so it stays 'reasonable'
     
     # N <- sum(LenDat)
     # FleetPars <- NULL
     # FleetPars$SL50 <- exp(Opt$par)[1] * StockPars$Linf
     # FleetPars$SL95 <- FleetPars$SL50  + exp(Opt$par)[2] * StockPars$Linf
     # FleetPars$MLLKnife <- NA
     # FleetPars$FM <- exp(Opt$par)[3]
	 
	 # Estimated parameters
     estFM <- exp(Opt$par[1])	
     estSL50 <- exp(Opt$par[2]) * StockPars$Linf
     estSL95 <- estSL50 + (exp(Opt$par[3]) * estSL50	)
	 
     runMod <- LBSPRFunc(MK=StockPars$MK, Linf=StockPars$Linf, CVLinf=StockPars$CVLinf, 
	 L50=StockPars$L50, L95=StockPars$L95, Beta=StockPars$FecB, FM=estFM, 
	 SL50=estSL50, SL95=estSL95, Mids=CAL_binsmid, P=0.005, Nage=100)
     
     # tt <- barplot(LenDat, names.arg=CAL_binsmid) 
     # lines(tt, runMod$ExpLenCatchFished * N, lwd=3)
	 # title(xYr)
      
     EstFM <-  min(5, estFM)
     EstSPR <- runMod$SPR	
     MiscList[[1]][xYr] <- EstSPR 
     MiscList[[3]][xYr] <- EstFM
     
     # if (xYr == length(IsEmpty)) {
      # # SPR v FM 
      # FMVec <- seq(from=0, to=5, length.out=100)
      # SaveSPR <- sapply(1:length(FMVec), function (xx) {
      # FleetPars$FM <- FMVec[xx]
      # EqSimMod_LB(StockPars, FleetPars, SizeBins, FitControl=NULL)$SPR
      # }	)
  	  # MiscList[[5]] <- cbind(FMVec, SaveSPR)
    
     # }
	 # if (xYr > 1) 
	   # points(MiscList[[1]][1:xYr])
	   # lines(Kalman(RawEsts=MiscList[[1]][1:xYr]), col=xYr)
    }
 }
  # Smoothed estimates - SPR
  MiscList[[2]] <- Kalman(RawEsts=MiscList[[1]]) 
  MiscList[[2]][MiscList[[2]] <0] <- 0.05
  MiscList[[2]][MiscList[[2]] > 1] <- 0.99
  
  
  # Smoothed estimates - FM
  MiscList[[4]] <- Kalman(RawEsts=MiscList[[3]])
  
  
  # MiscList[[6]] <- DLM_data # Store current data for checking
  return(MiscList)
}

OptRoutine <- function(Pars, LenDat, StockPars, Mids) {
   MK <- StockPars$MK 
   Linf <- StockPars$Linf
   CVLinf <- StockPars$CVLinf
   L50 <- StockPars$L50
   L95 <- StockPars$L95
   Beta <- StockPars$FecB
   
   P <- 0.01
   Nage <- 100
   
   FM <- exp(Pars[1])
   SL50 <- exp(Pars[2]) * Linf
   SL95 <- SL50 + (exp(Pars[3]) * SL50)
   runMod <- LBSPRFunc(MK, Linf, CVLinf, L50, L95, Beta, FM, SL50, SL95, Mids, P, Nage)
   LenProb <- LenDat/sum(LenDat)
   ind <- (LenProb > 0)
   return(-sum(LenDat[ind] * log(runMod$PropLen[ind]/LenProb[ind])))
} 

LBSPRFunc <- function(MK=1.5, Linf=100, CVLinf=0.1, L50=55, L95=60, Beta=3, FM=1, SL50=30, SL95=35, 
	Mids, P=0.01, Nage=100) {
  LenMids <- Mids
  By <- Mids[2] - Mids[1] 
  LenBins <- seq(from=0, by=By, length.out=length(LenMids)+1)	
  x <- seq(from=0, to=1, length.out=Nage) # relative age vector
  EL <- (1-P^(x/MK)) * Linf # length at relative age 
  rLens <- EL/Linf # relative length 
  SDL <- EL * CVLinf # standard deviation of length-at-age
  
  Nlen <- length(LenMids) 
  Prob <- matrix(NA, nrow=Nage, ncol=Nlen)
  Prob[,1] <- pnorm((LenBins[2] - EL)/SDL, 0, 1) # probablility of length-at-age
  for (i in 2:(Nlen-1)) {
    Prob[,i] <- pnorm((LenBins[i+1] - EL)/SDL, 0, 1) - pnorm((LenBins[i] - EL)/SDL, 0, 1)
  }
  Prob[,Nlen] <- 1 - pnorm((LenBins[Nlen] - EL)/SDL, 0, 1)

  SL <- 1/(1+exp(-log(19)*(LenMids-SL50)/(SL95-SL50))) # Selectivity at length
  Sx <- apply(t(Prob) * SL, 2, sum) # Selectivity at relative age 
  MSX <- cumsum(Sx) / seq_along(Sx) # Mean cumulative selectivity for each age 
  Ns <- (1-rLens)^(MK+(MK*FM)*MSX) # number at relative age in population
  
  Cx <- t(t(Prob) * SL) # Conditional catch length-at-age probablilities
  # for (X in seq_along(x)) {
    # if (sum(Cx[X,]) > 0 )Cx[X,] <- Cx[X,]/sum(Cx[X,])
    # if (sum(Cx[X,]) == 0 )Cx[X,] <- 0  
  # }
  
  Nc <- apply(Ns * Cx, 2, sum)
  
  Ml <- 1/(1+exp(-log(19)*(LenMids-L50)/(L95-L50))) # Maturity at length
  Ma <-  apply(t(Prob) * Ml, 2, sum) # Maturity at relative age 
  
  N0 <- (1-rLens)^MK # Unfished numbers-at-age 
  SPR <- sum(Ma * Ns * rLens^Beta)/sum(Ma * N0 * rLens^Beta)
  
  Output <- NULL 
  Output$SPR <- SPR 
  Output$LenMids <- LenMids
  Output$PropLen <- Nc/sum(Nc)
  return(Output)
}  



#########################
# LB-SPR Model with GTG #
######################### 
EqSimMod_LB <- function(Stock, Fleet, SizeBins, FitControl=NULL)  {
  if (length(FitControl) > 0) { # Adjust life-history parameters for fitness
    AdjustM <- FitControl$AdjustM
	if (AdjustM) {
	  OutPath <- FitControl$OutPath
	  PredictPars <- FitControl$PredictPars
      Stock <- AdjustedParsFun(Stock, Fleet, SizeBins, AdjustM=AdjustM, OutPath=OutPath, PredictPars=PredictPars)
    }	  
  }

  # Stock$House-keeping stuff 
  NGTG <- Stock$NGTG 
  GTGLinfBy <- Stock$GTGLinfBy 
  if (!exists("GTGLinfBy")) GTGLinfBy <- NA
  if (is.null(GTGLinfBy)) GTGLinfBy <- NA
  Linf <- Stock$Linf
  CVLinf <- Stock$CVLinf 
  MaxSD <- Stock$MaxSD 
  MK <- Stock$MK 
  L50 <- Stock$L50 
  L95 <- Stock$L95 
  Walpha <- Stock$Walpha 
  Wbeta <- Stock$Wbeta 
  FecB <- Stock$FecB 
  Steepness <- Stock$Steepness 
  Mpow <- Stock$Mpow
  R0 <- Stock$R0 
  CentLinf <- Stock$CentLinf 
  CentMpar <- Stock$CentMpar 
  CentKpar <- Stock$CentKpar
  Mslope <- Stock$Mslope
  
  # These parameters for accounting for differential in fitness across groups 
  if (!exists("CentLinf")) CentLinf <- Linf 
  if (is.null(CentLinf)) CentLinf <- Stock$Linf 
  if (!exists("Mpar")) Mpar <- 0.2 # Only used to calculate M/K ratio
  if (is.null(Mpar)) Mpar <- 0.2
  if (!exists("CentMpar")) CentMpar <- Mpar 
  if (is.null(CentMpar)) CentMpar <- Mpar 
  if (!exists("CentKpar")) CentKpar <- Kpar 
  if (is.null(CentKpar)) CentKpar <- CentMpar / MK
  if (!exists("Mslope")) Mslope <- 0 
  if (is.null(Mslope)) Mslope <- 0 
  	
  SL50 <- Fleet$SL50
  SL95 <- Fleet$SL95 
  MLLKnife <- Fleet$MLLKnife
  FM <- Fleet$FM 
  
  Linc <- SizeBins$Linc 
  ToSize <- SizeBins$ToSize

  # Error Catches #
  if (!(exists("NGTG") | exists("GTGLinfBy"))) stop("NGTG or GTGLinfBy must be specified")
  if (!exists("R0")) R0 <- 1E6
  
  SDLinf <- CVLinf * Linf # Standard Deviation of Length-at-Age # Assumed constant CV here
  
  # Set up Linfs for the different GTGs
  if (exists("NGTG") & !exists("GTGLinfBy")) {
    DiffLinfs <- seq(from=CentLinf-MaxSD*SDLinf, to=CentLinf+MaxSD*SDLinf, length=NGTG)
	GTGLinfBy <- DiffLinfs[2]-DiffLinfs[1]
  } else  if (!exists("NGTG") & exists("GTGLinfBy")) {
    DiffLinfs <- seq(from=CentLinf-MaxSD*SDLinf, to=CentLinf+MaxSD*SDLinf, by=GTGLinfBy)
	NGTG <- length(DiffLinfs)
  } else if (exists("NGTG") & exists("GTGLinfBy")) {
    if (!is.na(GTGLinfBy)) {
	  DiffLinfs <- seq(from=CentLinf-MaxSD*SDLinf, to=CentLinf+MaxSD*SDLinf, by=GTGLinfBy)
	  NGTG <- length(DiffLinfs)
	} 
	if (is.na(GTGLinfBy)) {
	  DiffLinfs <- seq(from=CentLinf-MaxSD*SDLinf, to=CentLinf+MaxSD*SDLinf, length=NGTG)
	  GTGLinfBy <- DiffLinfs[2]-DiffLinfs[1]
	}  
  } 
  # Distribute Recruits across GTGS 
  RecProbs <- dnorm(DiffLinfs, CentLinf, sd=SDLinf)/sum(dnorm(DiffLinfs, CentLinf, sd=SDLinf)) 
  
  # Length Bins 
  if (is.null(ToSize)) ToSize <- max(DiffLinfs, Linf + 3 * SDLinf)
  LenMids <- seq(from=0.5*Linc, by=Linc, to=ToSize)
  LenBins <- seq(from=0, by=Linc, length.out=length(LenMids)+1)
  
  Weight <- Walpha * LenMids^Wbeta
  
  # Maturity and Fecundity for each GTG 
  L50GTG <- L50/Linf * DiffLinfs # Maturity at same relative size
  L95GTG <- L95/Linf * DiffLinfs # Assumes maturity age-dependant 
  DeltaGTG <- L95GTG - L50GTG
  MatLenGTG <- sapply(seq_along(DiffLinfs), function (X) 1.0/(1+exp(-log(19)*(LenMids-L50GTG[X])/DeltaGTG[X])))
  FecLenGTG <- MatLenGTG * LenMids^FecB # Fecundity across GTGs - no scaling parameter atm

  VulLen <- 1.0/(1+exp(-log(19)*(LenBins-(SL50+0.5*Linc))/((SL95+0.5*Linc)-(SL50+0.5*Linc)))) # Selectivity-at-Length
  if (!is.na(MLLKnife)) { # Knife-edge selectivity
    VulLen[LenBins <= MLLKnife] <- 0
	VulLen[LenBins > MLLKnife] <- 1
	SL95 <- SL50 <- NA 
  }

  # Add dome-shaped selectivity curve 
  # Add F-mortality below MLL
  SelLen <- VulLen # Selectivity is equal to vulnerability currently
  
  # Life-History Ratios 
  ModMK <- CentMpar/CentKpar
  MKL <- ModMK * (CentLinf/(LenBins+0.5*Linc))^Mpow # M/K ratio for each length class
  # Matrix of MK for each GTG
  MKMat <- sapply(seq_along(DiffLinfs), function(X) MKL + Mslope*(DiffLinfs[X] - CentLinf))

  # ModMK + Mslope*(DiffLinfs-CentLinf)  #
  # # Debugging
   # matplot(MKMat)

  FK <- FM * ModMK # F/K ratio 
  FKL <- FK * SelLen # F/K ratio for each length class   
  # FkL[Legal == 0] <- FkL[Legal == 0] * DiscardMortFrac 
  ZKLMat <- MKMat + FKL # Z/K ratio (total mortality) for each GTG
    
  # Set Up Empty Matrices 
  NPRFished <- NPRUnfished <- matrix(0, nrow=length(LenBins), ncol=NGTG) # number-per-recruit at length
  NatLUnFishedPop <- NatLFishedPop <- NatLUnFishedCatch <- NatLFishedCatch <- FecGTGUnfished <- matrix(0, nrow=length(LenMids), ncol=NGTG) # number per GTG in each length class 
  NPRFished[1, ] <- NPRUnfished[1, ] <- RecProbs * R0# Distribute Recruits into first length class
  for (L in 2:length(LenBins)) { # Calc number at each size class
    NPRUnfished[L, ] <- NPRUnfished[L-1, ] * ((DiffLinfs-LenBins[L])/(DiffLinfs-LenBins[L-1]))^MKMat[L-1, ]
    NPRFished[L, ] <- NPRFished[L-1, ] * ((DiffLinfs-LenBins[L])/(DiffLinfs-LenBins[L-1]))^ZKLMat[L-1, ]
	ind <- DiffLinfs  < LenBins[L]
	NPRFished[L, ind] <- 0
	NPRUnfished[L, ind] <- 0
  } 
  NPRUnfished[is.nan(NPRUnfished)] <- 0
  NPRFished[is.nan(NPRFished)] <- 0
  NPRUnfished[NPRUnfished < 0] <- 0
  NPRFished[NPRFished < 0] <- 0
  
  for (L in 1:length(LenMids)) { # integrate over time in each size class
    NatLUnFishedPop[L, ] <- (NPRUnfished[L,] - NPRUnfished[L+1,])/MKMat[L, ]
    NatLFishedPop[L, ] <- (NPRFished[L,] - NPRFished[L+1,])/ZKLMat[L, ]  
	FecGTGUnfished[L, ] <- NatLUnFishedPop[L, ] * FecLenGTG[L, ]
  }
 
  VulLen2 <- 1.0/(1+exp(-log(19)*(LenMids-SL50)/(SL95-SL50))) # Selectivity-at-Length
  # print(LenMids)
  # print(c(SL50, SL95))
  # plot(LenMids, VulLen2)
  
  if (!is.na(MLLKnife))  { # Knife-edge selectivity
    VulLen2[LenMids <= MLLKnife] <- 0
	VulLen2[LenMids > MLLKnife] <- 1
	SL95 <- SL50 <- NA 
  }
  
  # points(LenMids, VulLen2, col="red")
  
  # print(cbind(LenMids, VulLen2))
  NatLUnFishedCatch <- NatLUnFishedPop * VulLen2 # Unfished Vul Pop
  NatLFishedCatch <- NatLFishedPop * VulLen2 # Catch Vul Pop
  
  # plot(LenMids, apply(NatLFishedCatch, 1, sum), type="p")
  # matplot(LenMids, (NatLFishedCatch), type="l")
  
  # Expected Length Structure - standardised 
  ExpectedLenCatchFished <- apply(NatLFishedCatch, 1, sum)/sum(apply(NatLFishedCatch, 1, sum))
  ExpectedLenPopFished <- apply(NatLFishedPop, 1, sum)/sum(apply(NatLFishedPop, 1, sum))
  ExpectedLenCatchUnfished <- apply(NatLUnFishedCatch, 1, sum)/sum(apply(NatLUnFishedCatch, 1, sum))
  ExpectedLenPopUnfished <- apply(NatLUnFishedPop, 1, sum)/sum(apply(NatLUnFishedPop, 1, sum))
  
  # Calc SPR
  EPR0 <- sum(NatLUnFishedPop * FecLenGTG) # Eggs-per-recruit Unfished
  EPRf <- sum(NatLFishedPop * FecLenGTG) # Eggs-per-recruit Fished
  SPR <- EPRf/EPR0 
  
  # Equilibrium Relative Recruitment
  recK <- (4*Steepness)/(1-Steepness) # Goodyear compensation ratio 
  reca <- recK/EPR0
  recb <- (reca * EPR0 - 1)/(R0*EPR0)
  RelRec <- max(0, (reca * EPRf-1)/(recb*EPRf))
  # RelRec/R0 - relative recruitment 
  YPR <- sum(NatLFishedPop  * Weight * VulLen2) * FM 
  Yield <- YPR * RelRec
    
  # Calc Unfished Fitness 
  Fit <- apply(FecGTGUnfished, 2, sum, na.rm=TRUE) # Total Fecundity per Group
  FitPR <- Fit/RecProbs # Fitness per-recruit
  FitPR <- FitPR/median(FitPR)
  ## Debugging
  # plot(FitPR, ylim=c(0,2)) # Should be relatively flat for equal fitness across GTG
    
  ObjFun <- sum((FitPR - median(FitPR, na.rm=TRUE))^2, na.rm=TRUE) # This needs to be minimised to make fitness approximately equal across GTG - by adjusting Mslope 
  Pen <- 0; if (min(MKMat) <= 0 ) Pen <- (1/abs(min(MKMat)))^2 * 1E12 # Penalty for optimising Mslope   
  ObjFun <- ObjFun + Pen
  # print(cbind(Mslope, ObjFun, Pen))

  # Calculate spawning-per-recruit at each size class
  SPRatsize <- cumsum(rowSums(NatLUnFishedPop * FecLenGTG))
  SPRatsize <- SPRatsize/max(SPRatsize)

  Output <- NULL 
  Output$SPR <- SPR
  Output$Yield <- Yield 
  Output$YPR <- YPR
  Output$ExpLenCatchFished <- ExpectedLenCatchFished
  Output$ExpLenPopFished <- ExpectedLenPopFished
  Output$ExpLenCatchUnfished <- ExpectedLenCatchUnfished
  Output$ExpLenPopUnfished <- ExpectedLenPopUnfished
  Output$NatLFishedPop <- NatLFishedPop
  Output$NatLUnFishedCatch <- NatLUnFishedCatch
  Output$NatLUnFishedPop <- NatLUnFishedPop
  Output$NatLFishedCatch <- NatLFishedCatch
  Output$LenBins <- LenBins
  Output$LenMids <- LenMids
  Output$NGTG <- NGTG
  Output$GTGdL <- DiffLinfs[2] - DiffLinfs[1]
  Output$DiffLinfs <- DiffLinfs
  Output$RecProbs <- RecProbs
  Output$Weight <- Weight
  Output$Winf <- Walpha * Linf^Wbeta
  Output$FecLen <- FecLenGTG 
  Output$MatLen <- MatLenGTG 
  Output$SelLen <- SelLen
  Output$MKL <- MKL
  Output$MKMat <- MKMat 
  Output$FKL <- FKL 
  Output$ZKLMat <- ZKLMat 
  Output$ObjFun <- ObjFun 
  Output$Pen <- Pen
  Output$FitPR <- FitPR
  Output$Diff <- range(FitPR)[2] - range(FitPR)[1]
  Output$L50GTG <- L50GTG 
  Output$L95GTG <- L95GTG
  Output$SPRatsize <- SPRatsize
  Output$RelRec <- RelRec
  return(Output)
}


LBSPR_GTG <- function(x, DLM_data, yrsmth=5, perc=pstar,reps=reps) {
  dependencies="DLM_data@CAL, DLM_data@CAL_bins, DLM_data@vbLinf, 
	DLM_data@vbK, DLM_data@Mort, DLM_data@L50, DLM_data@L95, DLM_data@wlb, DLM_data@wla,
	DLM_data@Abun, DLM_data@CV_Abun"

  # Save other stuff for smoothing estimates
  TotYears <- nrow(DLM_data@CAL[1,,]) # How many years of length data exist
  if (length(DLM_data@Misc[[x]]) == 0) { # Misc List is empty
    # Create Empty List Object
	MiscList <- rep(list(0), 5) # Create empty list
	MiscList[[1]] <- rep(NA, TotYears) # SPR ests
	MiscList[[2]] <- rep(NA, TotYears) # Smoothed SPR ests
	MiscList[[3]] <- rep(NA, TotYears) # FM ests
	MiscList[[4]] <- rep(NA, TotYears) # Smoothed FM ests
	MiscList[[5]] <- list()
  }
  if (length(DLM_data@Misc[[x]]) != 0) MiscList <- DLM_data@Misc[[x]]
  
  # Add Extra Row when needed 
  if (length(MiscList[[1]]) < TotYears) {
    Diff <- TotYears - length(DLM_data@Misc[[x]][[1]])
    MiscList[[1]] <- append(MiscList[[1]], rep(NA,Diff)) 
	MiscList[[2]] <- append(MiscList[[2]], rep(NA,Diff)) 
	MiscList[[3]] <- append(MiscList[[3]], rep(NA,Diff)) 
	MiscList[[4]] <- append(MiscList[[4]], rep(NA,Diff)) 
  }

  NEmpty <- sum(is.na(MiscList[[1]])) # Number of empty spots
  IsEmpty <- which(is.na(MiscList[[1]]))

 StockPars <- NULL
 StockPars$NGTG <- 41
 StockPars$GTGLinfBy <- NA
 StockPars$Linf <- DLM_data@vbLinf[x]
 StockPars$CVLinf <- 0.1 # NEED TO ADD THIS TO INPUT VARIABLES
 StockPars$MaxSD <- 2 
 StockPars$MK  <- DLM_data@Mort[x] / DLM_data@vbK[x]
 StockPars$L50 <- DLM_data@L50[x] 
 StockPars$L95 <- DLM_data@L95[x] 
 StockPars$Walpha <- DLM_data@wla[x]
 StockPars$Wbeta <- DLM_data@wlb[x]
 StockPars$FecB  <- DLM_data@wlb[x]
 StockPars$Steepness <- 0.9 # DLM_data@steep[x] not required
 StockPars$Mpow <- 0
 StockPars$R0  <- 1000
 StockPars$CentLinf <- NULL
 StockPars$CentMpar <- NULL
 StockPars$CentKpar <- NULL
 StockPars$Mslope <- 0 
 
 # Run Assessment for every year where data exists 
 if (NEmpty !=0) {
   for (xYr in IsEmpty[1]:max(IsEmpty)) {
    
     # Average length data over years if more than *yrsmth*
     if (xYr < yrsmth) ind <- xYr
     if (xYr >= yrsmth) { 
       # ind <- (length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  	 ind <- (xYr-(yrsmth-1)):xYr
     }	 
     if (length(ind) > 1) LenDat <- apply(DLM_data@CAL[x, ind,], 2, mean, na.rm=TRUE)
     if (length(ind) == 1) LenDat <- DLM_data@CAL[x, ind,]
     
     binWidth <- DLM_data@CAL_bins[2] - DLM_data@CAL_bins[1]
     CAL_binsmid <- seq(from=0.5*binWidth, by=binWidth, length=length(DLM_data@CAL_bins)-1)
     
     SizeBins <- NULL
     SizeBins$Linc <- binWidth
     SizeBins$ToSize <- max(CAL_binsmid)
     
     SL50Start <- CAL_binsmid[which.max(LenDat)]
     DeltaStart <- 0.1 * SL50Start
     FMStart <- 1 
     Starts <- log(c(SL50Start/StockPars$Linf, DeltaStart/StockPars$Linf, FMStart))
     Lower <- log(c(0.1, 0.1, 0.001))
     Upper <- log(c(0.9, 0.9, 20))
     # run optimization
     Opt <- nlminb(Starts, OptRoutine, LenDat=LenDat, Stock=StockPars, SizeBins=SizeBins, lower=Lower, upper=Upper)
     # need to add penalty for selectivity so it stays 'reasonable'
     
     N <- sum(LenDat)
     FleetPars <- NULL
     FleetPars$SL50 <- exp(Opt$par)[1] * StockPars$Linf
     FleetPars$SL95 <- FleetPars$SL50  + exp(Opt$par)[2] * StockPars$Linf
     FleetPars$MLLKnife <- NA
     FleetPars$FM <- exp(Opt$par)[3]
     runMod <- EqSimMod_LB(StockPars, FleetPars, SizeBins, FitControl=NULL)
     
     # tt <- barplot(LenDat, names.arg=CAL_binsmid) 
     # lines(tt, runMod$ExpLenCatchFished * N, lwd=3)
	 # title(xYr)
      
     EstFM <-  min(5, FleetPars$FM)
     EstSPR <- runMod$SPR	
     MiscList[[1]][xYr] <- EstSPR 
     MiscList[[3]][xYr] <- EstFM
     
     # if (xYr == length(IsEmpty)) {
      # # SPR v FM 
      # FMVec <- seq(from=0, to=5, length.out=100)
      # SaveSPR <- sapply(1:length(FMVec), function (xx) {
      # FleetPars$FM <- FMVec[xx]
      # EqSimMod_LB(StockPars, FleetPars, SizeBins, FitControl=NULL)$SPR
      # }	)
  	  # MiscList[[5]] <- cbind(FMVec, SaveSPR)
    
     # }
    }
 }
  # Smoothed estimates - SPR
  MiscList[[2]] <- Kalman(RawEsts=MiscList[[1]]) 
  MiscList[[2]][MiscList[[2]] <0] <- 0.05
  MiscList[[2]][MiscList[[2]] > 1] <- 0.99
  # Smoothed estimates - FM
  MiscList[[4]] <- Kalman(RawEsts=MiscList[[3]])
  
  
  # MiscList[[6]] <- DLM_data # Store current data for checking
  return(MiscList)
}

OptRoutine_GTG <- function(FleetPars, LenDat, Stock, SizeBins, FitControl=NULL) {
  Fleet <- NULL
  Fleet$SL50 <- exp(FleetPars)[1] * Stock$Linf
  Fleet$SL95 <- Fleet$SL50  + (exp(FleetPars)[2] * Stock$Linf)
  Fleet$MLLKnife <- NA
  Fleet$FM <- exp(FleetPars)[3]
  runMod <- EqSimMod_LB(Stock, Fleet, SizeBins, FitControl=NULL)
  LenProb <- LenDat/sum(LenDat)
  ind <- (LenProb > 0)
  runMod$ExpLenCatchFished[runMod$ExpLenCatchFished == 0] <- tiny
  NLL <- -sum(LenDat[ind] * log(runMod$ExpLenCatchFished[ind]/LenProb[ind]))
  return(NLL)
}


LBSPR_TAC <- function(x, DLM_data, yrsmth=5, perc=pstar,reps=reps) {

  MiscList <- LBSPR(x, DLM_data, yrsmth=5, perc=pstar,reps=reps)
  
  XX <- 1:4 
  YY <- MiscList[[2]][(length(MiscList[[2]]) - (max(XX)-1)):length(MiscList[[2]])]
  
  EstSPR <- YY[length(YY)]
  
  TgSPR <- 0.4
  h <- DLM_data@steep[x]
  SPRLim <- -(2*(h-1))/(3*h+1) # SPR that results in 0.5 R0
  
  phi1 <- 6
  phi2 <- 1
  
  MaxDw <- -0.3
  MaxUp <- 0.3
  
  minSlope <- 0.01
  
  Slope <- coef(lm(YY~XX))[2]  
  # if (abs(Slope) < minSlope) Slope <- 0 
  Dist <- EstSPR - TgSPR 
  
  # Control Rule #
  Mod <- 0 
  Buff <- 0.1
  Buffer <- c(TgSPR - Buff,  TgSPR + Buff)
  inBuff <- FALSE
  belowTG <- FALSE 
  aboveTG <- FALSE
  slopeUp <- FALSE
  slopeDw <- FALSE 
  belowLim <- FALSE
  if (Dist < 0) belowTG <- TRUE 
  if (Dist > 0) aboveTG <- TRUE 
  if (EstSPR > min(Buffer) & EstSPR < max(Buffer)) inBuff <- TRUE
  if (Slope <= 0) slopeDw <- TRUE
  if (Slope > 0) slopeUp <- TRUE
  if (EstSPR < SPRLim) belowLim <- TRUE
   
  # If within buffer zone - only slope
  if (inBuff) Mod <- phi1 * Slope
  if (slopeUp & aboveTG) Mod <- phi1 * Slope +  phi2 * Dist
  if (slopeUp & belowTG) Mod <- phi1 * Slope 
  
  if (slopeDw & aboveTG) Mod <- phi1 * Slope 
  if (slopeDw & belowTG) Mod <- phi1 * Slope +  phi2 * Dist
  
  if (belowLim) Mod <- MaxDw
  
  Mod[Mod > MaxUp] <- MaxUp
  Mod[Mod < MaxDw] <- MaxDw
  Mod <- Mod + 1 
  

  # EstSPR <- seq(from=0, to=1, by=0.01) 

  # Mod <- phi1 * (EstSPR/TgSPR - 1)^3  + phi2 *(EstSPR/TgSPR - 1)
   
  # Mod[Mod > 1.2] <- 1.2
  # if (Slope <= 0 & Mod >= 1) Mod <- Mod * 0.8
  # if (Slope <= -0.05 & Mod >= 1) Mod <- 0.95
  # if (Slope < 0 & Mod < 1) Mod <- Mod * 0.8
  # if (Slope >= 0.05 & Mod >= 1) Mod <- Mod * 1.2 
  
  # Make new slope to target rule 
  
  # cbind(EstSPR, Mod)
  
  # plot(EstSPR, Mod, ylim=c(0,max(Mod)))
  TAC <- DLM_data@MPrec[x] * Mod
  TAC <- TACfilter(TAC)
  
  # SaveSPR <- MiscList[[5]][,2]
  # FMVec <-  MiscList[[5]][,1]
  # TgSPR <- 0.4 # Target SPR
  # FTarg <- FMVec[max(which(SaveSPR >= TgSPR))] * DLM_data@Mort[x]

  # DLM_data@MPrec
  # Abun <- max(DLM_data@Abun[x], tiny)
  # Ac <- trlnorm(reps,Abun,DLM_data@CV_Abun[x])
 
  # EstSPR <- MiscList[[2]][length(MiscList[[2]])]
  # TAC <- TACfilter(Ac*FTarg)
  
  # print(c(EstSPR, FTarg, Abun, TAC))
  
  Out <- list()
  Out[[1]] <- TAC 
  Out[[2]] <- MiscList
  
  # if (EstSPR <=0.6 & EstSPR >= 0.2) Out[[1]] <- 0.5 * TAC
  # if (EstSPR < 0.2) Out[[1]] <- TAC * 0.1
  return(Out) 
}
class(LBSPR_TAC)<-"DLM_output"


Kalman <- function(RawEsts, R=1, Q=0.1, Int=100) {
  # Kalman smoother and Rauch-Tung-Striebel smoother #http://read.pudn.com/downloads88/ebook/336360/Kalman%20Filtering%20Theory%20and%20Practice,%20Using%20MATLAB/CHAPTER4/RTSvsKF.m__.htm
  # R  # Variance of sampling noise
  # Q  # Variance of random walk increments
  # Int # Covariance of initial uncertainty
  Ppred <-  rep(Int, length(RawEsts))
  Pcorr <- xcorr <- xpred <- rep(0, length(RawEsts))
  # Kalman Filter
  for (X in 1:length(Ppred)) {
    if (X !=1) {
	  Ppred[X] <- Pcorr[X-1] + Q
	  xpred[X] <- xcorr[X-1]
	}
	W <- Ppred[X]/(Ppred[X] + R)
	xcorr[X] <- xpred[X] + W * (RawEsts[X] - xpred[X]) # Kalman filter estimate
	Pcorr[X] <- Ppred[X] - W * Ppred[X]
  }
  # Smoother 
  xsmooth <- xcorr
  for (X in (length(Pcorr)-1):1) {
    A <- Pcorr[X]/Ppred[X+1]
	xsmooth[X] <- xsmooth[X] + A*(xsmooth[X+1] - xpred[X+1]) 
  }
  return(xsmooth)

}

LBSPR_Eff <-function(x,DLM_data){
  MiscList <- LBSPR(x, DLM_data)
  EstSPR <- MiscList[[2]][length(MiscList[[2]])]
  TgSPR <- 0.4 
  phi1 <- 0.2
  phi2 <- 0.05
  Mod <- phi1 * (EstSPR/TgSPR - 1)^3  + phi2 *(EstSPR/TgSPR - 1)
  Allocate<-1
  Effort<-(1+Mod)
  Spatial<-c(0,1)
  Vuln<-rep(NA,DLM_data@MaxAge)
  c(Allocate, Effort, Spatial, Vuln)
}
class(LBSPR_Eff)<-"DLM_input"

# LBSPR_TAC2 <- function(x, DLM_data, yrsmth=5, perc=pstar,reps=reps) {

  # MiscList <- LBSPR(x, DLM_data, yrsmth=5, perc=pstar,reps=reps)
  
  # XX <- 1:4 
  # YY <- (MiscList[[2]][(length(MiscList[[2]]) - 3):length(MiscList[[2]])])
  
  # EstSPR <- MiscList[[2]][length(MiscList[[2]])]
  
  # TgSPR <- 0.5
  
  # # Traj 
  # B <- YY[length(YY)] - YY[1]
  # Slope <- atan(B)
  
  # # Tg 
  # EstSPR <- YY[length(YY)]
  # Bp <- EstSPR- TgSPR
  # Dist <- atan(Bp)
  
  # phi1 <- 0.1
  # phi2 <- 0.5
  # Mod <- phi1*Dist + phi2* Slope + 1 
  
  # Mod[Mod > 1.3] <- 1.3 
  # Mod[Mod < 0.7] <- 0.7 
  
  # # plot(EstSPR, Mod, ylim=c(0,max(Mod)))
  # TAC <- DLM_data@MPrec[x] * Mod
  # TAC <- TACfilter(TAC)
  
  # # SaveSPR <- MiscList[[5]][,2]
  # # FMVec <-  MiscList[[5]][,1]
  # # TgSPR <- 0.4 # Target SPR
  # # FTarg <- FMVec[max(which(SaveSPR >= TgSPR))] * DLM_data@Mort[x]

  # # DLM_data@MPrec
  # # Abun <- max(DLM_data@Abun[x], tiny)
  # # Ac <- trlnorm(reps,Abun,DLM_data@CV_Abun[x])
 
  # # EstSPR <- MiscList[[2]][length(MiscList[[2]])]
  # # TAC <- TACfilter(Ac*FTarg)
  
  # # print(c(EstSPR, FTarg, Abun, TAC))
  
  # Out <- list()
  # Out[[1]] <- TAC 
  # Out[[2]] <- MiscList
  
  # # if (EstSPR <=0.6 & EstSPR >= 0.2) Out[[1]] <- 0.5 * TAC
  # # if (EstSPR < 0.2) Out[[1]] <- TAC * 0.1
  # return(Out) 
# }
# class(LBSPR_TAC2)<-"DLM_output"

# LBSPR_TAC3 <- function(x, DLM_data, yrsmth=5, perc=pstar,reps=reps) {

  # MiscList <- LBSPR(x, DLM_data, yrsmth=5, perc=pstar,reps=reps)
  
  # XX <- 1:4 
  # YY <- (MiscList[[2]][(length(MiscList[[2]]) - 3):length(MiscList[[2]])])
  
  # EstSPR <- MiscList[[2]][length(MiscList[[2]])]
  
  # TgSPR <- 0.5
  
  # # Traj 
  # B <- YY[length(YY)] - YY[1]
  # Slope <- atan(B)
  
  # # Tg 
  # EstSPR <- YY[length(YY)]
  # Bp <- EstSPR- TgSPR
  # Dist <- atan(Bp)
  
  # phi1 <- 0.001
  # phi2 <- 0.5
  # Mod <- phi1*Dist + phi2* Slope + 1 
  
  # Mod[Mod > 1.3] <- 1.3 
  # Mod[Mod < 0.7] <- 0.7 
  
  # # plot(EstSPR, Mod, ylim=c(0,max(Mod)))
  # TAC <- DLM_data@MPrec[x] * Mod
  # TAC <- TACfilter(TAC)
  
  # # SaveSPR <- MiscList[[5]][,2]
  # # FMVec <-  MiscList[[5]][,1]
  # # TgSPR <- 0.4 # Target SPR
  # # FTarg <- FMVec[max(which(SaveSPR >= TgSPR))] * DLM_data@Mort[x]

  # # DLM_data@MPrec
  # # Abun <- max(DLM_data@Abun[x], tiny)
  # # Ac <- trlnorm(reps,Abun,DLM_data@CV_Abun[x])
 
  # # EstSPR <- MiscList[[2]][length(MiscList[[2]])]
  # # TAC <- TACfilter(Ac*FTarg)
  
  # # print(c(EstSPR, FTarg, Abun, TAC))
  
  # Out <- list()
  # Out[[1]] <- TAC 
  # Out[[2]] <- MiscList
  
  # # if (EstSPR <=0.6 & EstSPR >= 0.2) Out[[1]] <- 0.5 * TAC
  # # if (EstSPR < 0.2) Out[[1]] <- TAC * 0.1
  # return(Out) 
# }
# class(LBSPR_TAC3)<-"DLM_output"

# LBSPR_TAC4 <- function(x, DLM_data, yrsmth=5, perc=pstar,reps=reps) {

  # MiscList <- LBSPR(x, DLM_data, yrsmth=5, perc=pstar,reps=reps)
  
  # XX <- 1:4 
  # YY <- (MiscList[[2]][(length(MiscList[[2]]) - 3):length(MiscList[[2]])])
  
  # EstSPR <- MiscList[[2]][length(MiscList[[2]])]
  
  # TgSPR <- 0.5
  
  # # Traj 
  # B <- YY[length(YY)] - YY[1]
  # Slope <- atan(B)
  
  # # Tg 
  # EstSPR <- YY[length(YY)]
  # Bp <- EstSPR- TgSPR
  # Dist <- atan(Bp)
  
  # phi1 <- 0.3
  # phi2 <- 0.5
  # Mod <- phi1*Dist + phi2* Slope + 1 
  
  # Mod[Mod > 1.3] <- 1.3 
  # Mod[Mod < 0.7] <- 0.7 
  
  # # plot(EstSPR, Mod, ylim=c(0,max(Mod)))
  # TAC <- DLM_data@MPrec[x] * Mod
  # TAC <- TACfilter(TAC)
  
  # # SaveSPR <- MiscList[[5]][,2]
  # # FMVec <-  MiscList[[5]][,1]
  # # TgSPR <- 0.4 # Target SPR
  # # FTarg <- FMVec[max(which(SaveSPR >= TgSPR))] * DLM_data@Mort[x]

  # # DLM_data@MPrec
  # # Abun <- max(DLM_data@Abun[x], tiny)
  # # Ac <- trlnorm(reps,Abun,DLM_data@CV_Abun[x])
 
  # # EstSPR <- MiscList[[2]][length(MiscList[[2]])]
  # # TAC <- TACfilter(Ac*FTarg)
  
  # # print(c(EstSPR, FTarg, Abun, TAC))
  
  # Out <- list()
  # Out[[1]] <- TAC 
  # Out[[2]] <- MiscList
  
  # # if (EstSPR <=0.6 & EstSPR >= 0.2) Out[[1]] <- 0.5 * TAC
  # # if (EstSPR < 0.2) Out[[1]] <- TAC * 0.1
  # return(Out) 
# }
# class(LBSPR_TAC4)<-"DLM_output"

# # Old one 
# LBSPR1 <- function(x, DLM_data, yrsmth=5, perc=pstar,reps=reps) {
  # dependencies="DLM_data@CAL, DLM_data@CAL_bins, DLM_data@vbLinf, 
	# DLM_data@vbK, DLM_data@Mort, DLM_data@L50, DLM_data@L95, DLM_data@wlb, DLM_data@wla,
	# DLM_data@Abun, DLM_data@CV_Abun"
 
 # Stock <- NULL
 # Stock$NGTG <- 41
 # Stock$GTGLinfBy <- NA
 # Stock$Linf <- DLM_data@vbLinf[x]
 # Stock$CVLinf <- 0.1 # NEED TO ADD THIS TO INPUT VARIABLES
 # Stock$MaxSD <- 2 
 # Stock$MK  <- DLM_data@Mort[x] / DLM_data@vbK[x]
 # Stock$L50 <- DLM_data@L50[x] 
 # Stock$L95 <- DLM_data@L95[x] 
 # Stock$Walpha <- DLM_data@wla[x]
 # Stock$Wbeta <- DLM_data@wlb[x]
 # Stock$FecB  <- DLM_data@wlb[x]
 # Stock$Steepness <- 0.9 # DLM_data@steep[x] not required
 # Stock$Mpow <- 0
 # Stock$R0  <- 1000
 # Stock$CentLinf <- NULL
 # Stock$CentMpar <- NULL
 # Stock$CentKpar <- NULL
 # Stock$Mslope <- 0 
 
 # ind <- (length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
 # if (length(ind) > 1) LenDat <- apply(DLM_data@CAL[x, ind,], 2, mean, na.rm=TRUE)
 # if (length(ind) == 1) LenDat <- DLM_data@CAL[x, ind,]

 # binWidth <- DLM_data@CAL_bins[2] - DLM_data@CAL_bins[1]
 # CAL_binsmid <- seq(from=0.5*binWidth, by=binWidth, length=length(DLM_data@CAL_bins)-1)
 
 # SizeBins <- NULL
 # SizeBins$Linc <- binWidth
 # SizeBins$ToSize <- max(CAL_binsmid)
 
 # SL50Start <- CAL_binsmid[which.max(LenDat)]
 # DeltaStart <- 0.1 * SL50Start
 # FMStart <- 1 
 # Starts <- log(c(SL50Start/Stock$Linf, DeltaStart/Stock$Linf, FMStart))
 # Lower <- log(c(0.1, 0.1, 0.001))
 # Upper <- log(c(0.9, 0.9, 20))
 # # run optimization
 # Opt <- nlminb(Starts, OptRoutine, LenDat=LenDat, Stock=Stock, SizeBins=SizeBins, lower=Lower, upper=Upper)
 # # need to add penalty for selectivity so it stays 'reasonable'
 
 # # tt <- barplot(LenDat, names.arg=CAL_binsmid) 
 # N <- sum(LenDat)
 # Fleet <- NULL
 # Fleet$SL50 <- exp(Opt$par)[1] * Stock$Linf
 # Fleet$SL95 <- Fleet$SL50  + exp(Opt$par)[2] * Stock$Linf
 # Fleet$MLLKnife <- NA
 # Fleet$FM <- exp(Opt$par)[3]
 # runMod <- EqSimMod_LB(Stock, Fleet, SizeBins, FitControl=NULL)
  # # lines(tt, runMod$ExpLenCatchFished * N, lwd=3)
  
  # EstFM <-  min(5, Fleet$FM)
  # EstSPR <- runMod$SPR
 
  # # SPR v FM 
  # FMVec <- seq(from=0, to=5, length.out=100)
  # SaveSPR <- sapply(1:length(FMVec), function (xx) {
    # Fleet$FM <- FMVec[xx]
    # EqSimMod_LB(Stock, Fleet, SizeBins, FitControl=NULL)$SPR
  # }	)
 
  # TgSPR <- 0.5 # Target SPR
  # FTarg <- FMVec[max(which(SaveSPR >= TgSPR))] * DLM_data@Mort[x] 
  # Abun <- max(DLM_data@Abun[x], tiny)
  # Ac <- trlnorm(reps,Abun,DLM_data@CV_Abun[x])
  
  # TAC <- TACfilter(Ac*FTarg)
  # if (EstSPR > 0.6) return(TAC)
  # if (EstSPR <=0.6 & EstSPR >= 0.2) return(DLM_data@MPrec[x])
  # if (EstSPR < 0.2) return(TAC * 0.5)
 

 
# }
# class(LBSPR1)<-"DLM_output"
