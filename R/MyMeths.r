##################################
# LBSPR - Hordyk et al ICES 2015 #
##################################

LBSPRSim <- function(StockPars, FleetPars, SizeBins=NULL, P=0.001, Nage=201) {

  MK <- StockPars$MK 
  Linf <- StockPars$Linf
  CVLinf <- StockPars$CVLinf 
  L50 <- StockPars$L50 
  L95 <- StockPars$L95 
  Beta <- StockPars$FecB 
  MaxSD <- StockPars$MaxSD
  
  # Assumed constant CV here
  SDLinf <- CVLinf * Linf # Standard Deviation of Length-at-Age 
  if (is.null(SizeBins)) {
    SizeBins$Linc <- 5
	SizeBins$ToSize <- Linf + MaxSD * SDLinf
  }
  if (is.null(SizeBins$ToSize)) SizeBins$ToSize <- Linf + MaxSD * SDLinf
  Linc <- SizeBins$Linc 
  ToSize <- SizeBins$ToSize
  
  FM <- FleetPars$FM 
  SL50 <- FleetPars$SL50 
  SL95 <- FleetPars$SL95 
  
  LenBins <- seq(from=0, by=Linc, to=ToSize)	
  LenMids <- seq(from=0.5*Linc, by=Linc, length.out=length(LenBins)-1)
  x <- seq(from=0, to=1, length.out=Nage) # relative age vector
  EL <- (1-P^(x/MK)) * Linf # length at relative age 
  rLens <- EL/Linf # relative length 
  SDL <- EL * CVLinf # standard deviation of length-at-age
  
  Nlen <- length(LenMids) 
  Prob <- matrix(NA, nrow=Nage, ncol=Nlen)
  Prob[,1] <- pnorm((LenBins[2] - EL)/SDL, 0, 1) # probablility of length-at-age
  for (i in 2:(Nlen-1)) {
    Prob[,i] <- pnorm((LenBins[i+1] - EL)/SDL, 0, 1) - 
		pnorm((LenBins[i] - EL)/SDL, 0, 1)
  }
  Prob[,Nlen] <- 1 - pnorm((LenBins[Nlen] - EL)/SDL, 0, 1)
  
  # Truncate normal dist at MaxSD 
  mat <- array(1, dim=dim(Prob))
  for (X in 1:Nage) {
    ind <- which(abs((LenMids - EL[X]) /SDL[X]) >= MaxSD)
    mat[X,ind] <- 0
  }
  
  Prob <- Prob * mat

  SL <- 1/(1+exp(-log(19)*(LenMids-SL50)/(SL95-SL50))) # Selectivity at length
  Sx <- apply(t(Prob) * SL, 2, sum) # Selectivity at relative age 
  MSX <- cumsum(Sx) / seq_along(Sx) # Mean cumulative selectivity for each age 
  Ns <- (1-rLens)^(MK+(MK*FM)*MSX) # number at relative age in population
  
  Cx <- t(t(Prob) * SL) # Conditional catch length-at-age probablilities  
  Nc <- apply(Ns * Cx, 2, sum) # 
  Pop <- apply(Ns * Prob, 2, sum)
  
  Ml <- 1/(1+exp(-log(19)*(LenMids-L50)/(L95-L50))) # Maturity at length
  Ma <-  apply(t(Prob) * Ml, 2, sum) # Maturity at relative age 
  
  N0 <- (1-rLens)^MK # Unfished numbers-at-age 
  SPR <- sum(Ma * Ns * rLens^Beta)/sum(Ma * N0 * rLens^Beta)
  
  Output <- NULL 
  Output$SPR <- SPR 
  Output$LenMids <- LenMids
  Output$PropLen <- Nc/sum(Nc)
  Output$Pop <- Pop
  
  Output$LCatchFished <- Nc/sum(Nc)
  Output$LPopFished <- Pop
  Output$LCatchUnfished <- apply(N0 * Cx, 2, sum)
  return(Output)
}  

##########################
# Optimisation Functions #
##########################
OptFun <- function(tryFleetPars, LenDat, StockPars, SizeBins=NULL, 
	mod=c("GTG", "LBSPR")) {
  Fleet <- NULL
  Fleet$SL50 <- exp(tryFleetPars[1]) * StockPars$Linf
  Fleet$SL95 <- Fleet$SL50  + (exp(tryFleetPars[2]) * StockPars$Linf)
  Fleet$MLLKnife <- NA
  Fleet$FM <- exp(tryFleetPars[3])
  
  if (mod == "GTG") runMod <-  GTGLBSPRSim(StockPars, Fleet, SizeBins)
  if (mod == "LBSPR") runMod <- LBSPRSim(StockPars, Fleet, SizeBins)
  
  LenDat <- LenDat + 1E-15 # add tiny constant for zero catches
  LenProb <- LenDat/sum(LenDat)
  predProb <- runMod$LCatchFished 
  predProb <- predProb + 1E-15 # add tiny constant for zero catches
  NLL <- -sum(LenDat * log(predProb/LenProb))
  
  if(!is.finite(NLL)) return(1E9)
  
  # add penalty for SL50 
  trySL50 <- exp(tryFleetPars[1])
  PenVal <- NLL
  Pen <- dbeta(trySL50, shape1=5, shape2=0.01) * PenVal
  if (Pen == 0) Pen <- PenVal * trySL50
  
  # plot(xx, dbeta(xx, shape1=5, shape2=0.01) )
  
  NLL <- NLL+Pen 

  return(NLL)
}

DoOpt <- function(StockPars, LenDat, SizeBins=NULL, mod=c("GTG", "LBSPR")) {
  
  SDLinf <- StockPars$CVLinf * StockPars$Linf
  if (is.null(SizeBins)) {
    SizeBins$Linc <- 5
	SizeBins$ToSize <- StockPars$Linf + StockPars$MaxSD * SDLinf
  }
  if (is.null(SizeBins$ToSize)) 
	SizeBins$ToSize <- StockPars$Linf + StockPars$MaxSD * SDLinf
  
  Linc <- SizeBins$Linc 
  ToSize <- SizeBins$ToSize
 
  LenBins <- seq(from=0, by=Linc, to=ToSize)	
  LenMids <- seq(from=0.5*Linc, by=Linc, length.out=length(LenBins)-1)
  
  sSL50 <- LenMids[which.max(LenDat)]/StockPars$Linf # Starting guesses
  sDel <- 0.2 * LenMids[which.max(LenDat)]/StockPars$Linf
  sFM <- 0.5 
  Start <- log(c(sSL50, sDel, sFM))
  
  # opt <- nlminb(Start, OptFun, LenDat=LenDat, StockPars=StockPars, 
	# SizeBins=SizeBins, mod=mod) 
	
  opt2 <- nlm(OptFun, Start, steptol=1e-4,gradtol=1e-4, LenDat=LenDat, 
	StockPars=StockPars, SizeBins=SizeBins, mod=mod)
  opt <- NULL
  opt$objective <- opt2$minimum
  opt$par <- opt2$estimate
	
  ModFailed <- FALSE 	
  if (opt$objective	== 1E9) ModFailed <- TRUE 
	# ,control= list(iter.max=300, eval.max=400, abs.tol=1E-20))
  # barplot(LenDat, names.arg=LenMids) 
  
  newFleet <- NULL 
  newFleet$FM <- exp(opt$par[3])
  newFleet$SL50 <- exp(opt$par[1]) * StockPars$Linf 
  newFleet$SL95 <- newFleet$SL50 + exp(opt$par[2]) * StockPars$Linf

  if (mod == "GTG") runMod <-  GTGLBSPRSim(StockPars, newFleet, SizeBins)
  if (mod == "LBSPR") runMod <- LBSPRSim(StockPars, newFleet, SizeBins)
  
  Out <- NULL 
  Out$Ests <- c(FM=newFleet$FM, SL50=newFleet$SL50, SL95=newFleet$SL95, 
	SPR=runMod$SPR)
  Out$PredLen <- runMod$LCatchFished * sum(LenDat)
  Out$ModFailed <- ModFailed
  return(Out)
}


# Run LBSPR Model for time-series of catch length composition data 
LBSPR <- function(x, DLM_data, yrsmth=1, perc=pstar,reps=reps) {
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
 StockPars$MK <- DLM_data@Mort[x] / DLM_data@vbK[x]
 StockPars$Linf <- DLM_data@vbLinf[x]
 StockPars$CVLinf <- 0.1 # NEED TO ADD THIS TO INPUT VARIABLES
 StockPars$L50 <- DLM_data@L50[x] 
 StockPars$L95 <- DLM_data@L95[x]
 StockPars$FecB <- DLM_data@wlb[x]
 StockPars$MaxSD <- 2
 
 # yrsmth not implemented here
 LenMatrix <- DLM_data@CAL[x, IsEmpty,]
 
 binWidth <- DLM_data@CAL_bins[2] - DLM_data@CAL_bins[1]
 CAL_binsmid <- seq(from=0.5*binWidth, by=binWidth, length=length(DLM_data@CAL_bins)-1)
     
 SizeBins <- NULL
 SizeBins$Linc <- binWidth
 SizeBins$ToSize <- max(CAL_binsmid) + 0.5*binWidth

 st <- Sys.time()
 if(sfIsRunning()){
    AllOpt <- sfSapply(1:length(IsEmpty), function (X) 
		DoOpt(StockPars, LenDat=LenMatrix[X,], SizeBins=SizeBins, mod="LBSPR"))
  } else {
    AllOpt <- sapply(1:length(IsEmpty), function (X)
		DoOpt(StockPars, LenDat=LenMatrix[X,], SizeBins=SizeBins, mod="LBSPR"))
  } 
  
 EstFM <- sapply(AllOpt[1,], "[[", 1)
 estSL50 <- sapply(AllOpt[1,], "[[", 2)
 estSL95 <- sapply(AllOpt[1,], "[[", 3)
 EstSPR <- sapply(AllOpt[1,], "[[", 4)
 EstFM[EstFM > 5] <- 5 
 
 Fails <- which(sapply(AllOpt[3,], "[[", 1))
 if (length(Fails) > 0) {
   EstFM[Fails] <- NA 
   EstSPR[Fails] <- NA
 }
 while(sum(is.na(EstFM)) > 0) { # if model failed, make same as last time 
   EstFM[is.na(EstFM)] <- EstFM[which(is.na(EstFM))-1]
   EstSPR[is.na(EstSPR)] <- EstSPR[which(is.na(EstSPR))-1]
 }
 
 MiscList[[1]][IsEmpty] <- EstSPR # Save estimate of SPR for smoothing
 MiscList[[3]][IsEmpty] <- EstFM # Save estimate of F/M for smoothing 
 
  # st <- Sys.time()
 # # Run Assessment for every year where data exists 
 # if (NEmpty !=0) {
   # for (xYr in IsEmpty[1]:max(IsEmpty)) {
     # # Average length data over years if more than *yrsmth*
     # if (xYr < yrsmth) ind <- xYr
     # if (xYr >= yrsmth) ind <- (xYr-(yrsmth-1)):xYr
     # if (length(ind) > 1) LenDat <- apply(DLM_data@CAL[x, ind,], 2, mean, na.rm=TRUE)
     # if (length(ind) == 1) LenDat <- DLM_data@CAL[x, ind,]
     
     # binWidth <- DLM_data@CAL_bins[2] - DLM_data@CAL_bins[1]
     # CAL_binsmid <- seq(from=0.5*binWidth, by=binWidth, length=length(DLM_data@CAL_bins)-1)
     
     # SizeBins <- NULL
     # SizeBins$Linc <- binWidth
     # SizeBins$ToSize <- max(CAL_binsmid) + 0.5*binWidth
     
	 # Opt <- DoOpt(StockPars, LenDat=LenDat, SizeBins=SizeBins, mod="LBSPR")
	 
	 # if (Opt$ModFailed) {
       # estFM <- NA
       # estSL50 <- NA
       # estSL95 <- NA
	   # EstFM <- MiscList[[1]][xYr-1] # If it fails, assume same as last year
	   # EstSPR <- MiscList[[3]][xYr-1]
	 # } else {
	   # # Estimated parameters
       # estFM <- Opt$Ests[1]
       # estSL50 <- Opt$Ests[2]
       # estSL95 <- Opt$Ests[3]
	   # EstFM <-  min(5, estFM)
       # EstSPR <- Opt$Ests[4]
	 # }
     # MiscList[[1]][xYr] <- EstSPR # Save estimate of SPR for smoothing
     # MiscList[[3]][xYr] <- EstFM # Save estimate of F/M for smoothing 
    # }
 # }
   # Sys.time() - st 
   
  
  # Smoothed estimates - SPR
  MiscList[[2]] <- Kalman(RawEsts=MiscList[[1]]) 
  MiscList[[2]][MiscList[[2]] <0] <- 0.05
  MiscList[[2]][MiscList[[2]] > 1] <- 0.99
 
  # Smoothed estimates - FM
  MiscList[[4]] <- Kalman(RawEsts=MiscList[[3]])
  
  return(MiscList)
}

# Kalman filter and Rauch-Tung-Striebel smoother
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


###############
# DLMtool MPs #
###############
# Note - these can be a bit slow - especially in the first year that the 
# methods are run - as they apply the assessment to every historical year.

LBSPR_ItTAC <- function(x, DLM_data, yrsmth=1, perc=pstar,reps=reps) {
 dependencies="DLM_data@CAL, DLM_data@CAL_bins, DLM_data@vbLinf, 
	DLM_data@vbK, DLM_data@Mort, LM_data@vbK, DLM_data@L50, DLM_data@L95, 
	DLM_data@wlb" 
  
  st <- Sys.time()	
  MiscList <- LBSPR(x, DLM_data, yrsmth=yrsmth, perc=pstar,reps=reps)
  nd <- Sys.time()
  nd - st   
  
  
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
 
  TAC <- DLM_data@MPrec[x] * Mod
  TAC <- TACfilter(TAC)
 
  Out <- list()
  Out[[1]] <- TAC 
  Out[[2]] <- MiscList
 
  return(Out) 
}
class(LBSPR_ItTAC)<-"DLM_output"

LBSPR_ItEff <- function(x, DLM_data, yrsmth=1, perc=pstar,reps=reps) {
 dependencies="DLM_data@CAL, DLM_data@CAL_bins, DLM_data@vbLinf, 
	DLM_data@vbK, DLM_data@Mort, LM_data@vbK, DLM_data@L50, DLM_data@L95, 
	DLM_data@wlb"
  MiscList <- LBSPR(x, DLM_data, yrsmth=yrsmth, perc=pstar,reps=reps)
  
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
  
  Allocate <- 1
  Effort <- DLM_data@MPrec[x] * Mod
  MiscList[[5]] <- append(MiscList[[5]], Effort)
  Spatial <- c(1,1)
  Vuln <- rep(NA,2)
  out <- c(Allocate, Effort, Spatial, Vuln)
   
  Out <- list()
  Out[[1]] <- out 
  Out[[2]] <- MiscList
 
  return(Out) 
}
class(LBSPR_ItEff)<-"DLM_input"

  

