
#' Length-based SPR model with HCR that iteratively adjusts TAC
#' 
#' Iteratively adjusts TAC based on distance between estimated and target SPR
#' (40\%), and slope of recent SPR estimates.
#' 
#' @usage LBSPR_ItTAC(x, Data, yrsmth=1,reps=5, ...)
#' @param x Simulation number
#' @param Data Data object
#' @param yrsmth Number of years to smooth length data - not currently used
#' @param reps Number of repetitions
#' @param ... ignored
#' @author A. Hordyk
#' @export LBSPR_ItTAC
LBSPR_ItTAC <- function(x, Data, yrsmth=1,reps=5, ...) {
 
 dependencies="Data@CAL, Data@CAL_bins, Data@vbLinf, Data@vbK, Data@Mort, LM_data@vbK, 
   Data@L50, Data@L95, Data@wlb, Data@MPrec" 
  
  if (is.na(Data@MPrec[x])) {
    # message("No previous TAC recommendation")
	return(NA) # 
  }
  
  # Run the LBSPR model 
  MiscList <- LBSPR(x, Data, yrsmth=yrsmth,reps=reps, ...)
  
  if(all(is.na(MiscList[[1]]))) return(NA)
  if(all(is.na(MiscList[[1]][,2,]))) return(NA)
  XX <- 1:4 
  
  YY <- MiscList[[1]][(length(MiscList[[1]][,2,1]) - (max(XX)-1)):length(MiscList[[1]][,2,1]), 2,]
  
  if (reps ==1) EstSPR <- YY[length(YY)]
  if (reps > 1) EstSPR <- YY[nrow(YY),]
  
  TgSPR <- 0.4
  Steep <- Data@steep[x]
  if (is.na(Steep)) Steep <- 0.6
  h <- trlnorm(reps, Steep, Data@CV_steep[x])
  SPRLim <- -(2*(h-1))/(3*h+1) # SPR that results in 0.5 R0
  
  phi1 <- 6
  phi2 <- 1
  
  MaxDw <- -0.3
  MaxUp <- 0.3
  
  minSlope <- 0.01
  if (reps > 1) colnames(YY) <- 1:reps
  if (reps > 1) {
    lms <- sapply(1:reps, function(tt) {
      y <- YY[,tt]
      lm(y~x, data=data.frame(x=XX, y=y))$coefficients
    })
    Slope <- lms[2,]
  } else {
    y <- YY
    lms <- lm(y~x, data=data.frame(x=XX, y=y))$coefficients
	Slope <- lms[2]
  }
  
  # if (abs(Slope) < minSlope) Slope <- 0 
  Dist <- EstSPR - TgSPR 
  
  # Control Rule #
  Mod <- rep(0, reps)
  Buff <- 0.1
  Buffer <- c(TgSPR - Buff,  TgSPR + Buff)
  inBuff <- rep(FALSE, reps)
  belowTG <- rep(FALSE, reps)
  aboveTG <- rep(FALSE, reps)
  slopeUp <- rep(FALSE, reps)
  slopeDw <- rep(FALSE, reps)
  belowLim <- rep(FALSE, reps)
  
  belowTG[Dist < 0] <- TRUE
  aboveTG[Dist > 0] <- TRUE
  inBuff[(EstSPR > min(Buffer) & EstSPR < max(Buffer))] <- TRUE
  slopeDw[Slope <=0] <- TRUE
  slopeUp[Slope > 0] <- TRUE
  belowLim[EstSPR < SPRLim] <- TRUE
  
  # If within buffer zone - only slope
  Mod[inBuff] <- phi1 * Slope[inBuff]
  Mod[slopeUp & aboveTG] <- phi1 * Slope[slopeUp & aboveTG] +  
    phi2 * Dist[slopeUp & aboveTG]
  Mod[slopeUp & aboveTG] <- phi1 * Slope[slopeUp & aboveTG] +  
    phi2 * Dist[slopeUp & aboveTG]
 
  Mod[slopeDw & aboveTG] <- phi1 * Slope[slopeDw & aboveTG]
  Mod[slopeDw & belowTG] <- phi1 * Slope[slopeDw & belowTG] +  
    phi2 * Dist[slopeDw & belowTG]
    
  Mod[belowLim]  <- MaxDw 
  Mod[Mod > MaxUp] <- MaxUp
  Mod[Mod < MaxDw] <- MaxDw
  Mod <- Mod + 1 
 
  TAC <- Data@MPrec[x] * Mod
  TAC <- TACfilter(TAC)
 
  Out <- list()
  Out[[1]] <- TAC 
  Out[[2]] <- MiscList
 
  return(Out) 
}
class(LBSPR_ItTAC)<-"Output"

#' Length-based SPR model with HCR that iteratively adjusts Effort
#' 
#' Iteratively adjusts Effort based on distance between estimated and target
#' SPR (40\%), and slope of recent SPR estimates.
#' 
#' 
#' @usage LBSPR_ItEff(x, Data, yrsmth=1,reps=5, ...)
#' @param x Simulation number
#' @param Data Data object
#' @param yrsmth Number of years to smooth length data - not currently used
#' @param reps Number of repetitions. Not currently used 
#' @param ... ignored
#' @author A. Hordyk
#' @importClassesFrom LBSPR LB_pars LB_lengths
#' @export LBSPR_ItEff
LBSPR_ItEff <- function(x, Data, yrsmth=1, reps=5, ...) {
 dependencies="Data@CAL, Data@CAL_bins, Data@vbLinf, 
	Data@vbK, Data@Mort, Data@vbK, Data@L50, Data@L95, 
	Data@wlb"
  MiscList <- LBSPR(x, Data, yrsmth=yrsmth,reps=reps, ...)
  if(all(is.na(MiscList[[1]]))) return(rep(NA, 6))
  if(all(is.na(MiscList[[1]][,2,]))) return(rep(NA, 6))
  
  XX <- 1:4 
  YY <- MiscList[[1]][(length(MiscList[[1]][,2,1]) - (max(XX)-1)):length(MiscList[[1]][,2,1]), 2,]
  
  reps <- 1 # only allow one rep 
  if (reps ==1) EstSPR <- YY[length(YY)]
  if (reps > 1) EstSPR <- YY[nrow(YY),]
  
  TgSPR <- 0.4
  Steep <- Data@steep[x]
  if (is.na(Steep)) Steep <- 0.6
  h <- Steep
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
  
  rec <- new("InputRec")
  rec@Effort <- Data@MPeff[x] * Mod
  rec@Misc <- MiscList
  rec
}
class(LBSPR_ItEff)<-"Input"

#' Length-based SPR model with HCR that iteratively adjusts Selectivity
#' 
#' Management Procedure which adjusts size-at-selection based on estimated SPR.
#' Entirely untested, and included at to demonstrate MPs of this type.
#' 
#' 
#' @usage LBSPR_ItSel(x, Data, yrsmth=1,reps=5, ...)
#' @param x Simulation number
#' @param Data Data object
#' @param yrsmth Number of years to smooth length data - not currently used
#' @param reps Number of repetitions. Not currently used 
#' @param ... ignored
#' @author A. Hordyk
#' @export LBSPR_ItSel
LBSPR_ItSel <- function(x, Data, yrsmth=1, reps=5, ...) {

 dependencies="Data@CAL, Data@CAL_bins, Data@vbLinf, 
	Data@vbK, Data@Mort, Data@vbK, Data@L50, Data@L95, 
	Data@wlb, Data@LenCV"
  MiscList <- LBSPR(x, Data, yrsmth=yrsmth,reps=reps)
  if(all(is.na(MiscList[[1]]))) return(rep(NA, 6))
  if(all(is.na(MiscList[[1]][,2,]))) return(rep(NA, 6))
  XX <- 1:4 
  YY <- MiscList[[1]][(length(MiscList[[1]][,2,1]) - (max(XX)-1)):length(MiscList[[1]][,2,1]), 2,]
  
  reps <- 1 # force reps to be one 
  if (reps ==1) EstSPR <- YY[length(YY)]
  if (reps > 1) EstSPR <- YY[nrow(YY),]
  
  TgSPR <- 0.4
  Steep <- Data@steep[x]
  if (is.na(Steep)) Steep <- 0.6
  h <- Steep
  SPRLim <- -(2*(h-1))/(3*h+1) # SPR that results in 0.5 R0
 
  
  rec <- new("InputRec") # create recommendation object

  if (EstSPR < TgSPR) {
    rec@LR5 <- Data@L50[x] * 1.05
    rec@LFR <- Data@L50[x] * 1.1
  }
  if (EstSPR < SPRLim) {
    rec@LR5 <- Data@L50[x] * 1.2
    rec@LFR <- Data@L50[x] * 1.25
 
  }
  if (EstSPR >= TgSPR) {
    rec@LR5 <- Data@L50[x] * 0.85
    rec@LFR <- Data@L50[x] * 0.9

  }
   

  rec@Misc <- MiscList
  rec
}
class(LBSPR_ItSel)<-"Input"

#' Apply the Length-based SPR model to DLMtool Data Object
#' 
#' 
#' @param x Simulation number
#' @param Data Data object
#' @param yrsmth Number of years to smooth length data - not currently used
#' @param reps Number of repetitions
#' @param lstyrs Last number of years to run model
#' @author A. Hordyk
#' @export LBSPR
LBSPR <- function(x, Data, yrsmth=1, reps=1, lstyrs=10) {
  if (length(Data@LHYear)<1) stop("LHYear must be set to last year of data", call.=FALSE)
  if (Data@LHYear <1) stop("LHYear must be set to last year of data", call.=FALSE)
  if (all(is.na(Data@CAL))) stop("No length data", call.=FALSE)
  if (length(Data@Misc) == 0) Data@Misc <- vector("list", 1)
  
  TotYears <- nrow(Data@CAL[1,,]) # How many years of length data exist
 
  if (is.null(TotYears)) TotYears <- length(Data@CAL[1,,])
  if(is.null(dim(Data@CAL[1,,]))) TotYears <- 1

  # if (lstyrs > TotYears) lstyrs <- TotYears
  # # Only apply model to lstyrs of data - only applies to the first time the model is run in MSE 
  years <- seq_along(Data@Year)
  LHYear <- which(Data@Year == Data@LHYear)
  CurrYear <- years[length(Data@Year)]
  
  #index <- NULL
  # if (!is.null(lstyrs)) {
    # DD <- dim(Data@CAL)
 	# TotYears <- nrow(Data@CAL[1,(DD[2]-lstyrs+1):DD[2],]) # How many years of length data exist
    # if (is.null(TotYears)) TotYears <- length(Data@CAL[1,(DD[2]-lstyrs+1):DD[2],]) 
	# index <- (DD[2]-lstyrs+1):DD[2]
  # }
  
  if (length(Data@Misc[[x]]) == 0) { # Misc List is empty
	# Create Empty Object
    MiscList <- rep(list(0), 2) # Create empty list
	MiscList[[1]] <- array(NA, dim=c(TotYears, 4, reps))
	colnames(MiscList[[1]]) <- c("Raw Est SPR", "Smooth Est SPR", "Raw Est F/M", "Smooth Est F/M")
	MiscList[[2]] <- list()
  }
  if (length(Data@Misc[[x]]) != 0) MiscList <- Data@Misc[[x]]
  
  # Add Extra Row when needed 
  if (nrow(MiscList[[1]][,1,1, drop=FALSE]) < TotYears) {
    Diff <- TotYears - nrow(Data@Misc[[x]][[1]])	
	newmat <- array(NA, dim=c(Diff, 4, reps))
	colnames(newmat) <- colnames(MiscList[[1]]) 
	MiscList[[1]] <- abind::abind(MiscList[[1]], newmat, along=1) 
  }

  NEmpty <- sum(is.na(MiscList[[1]][,1,1])) # Number of empty spots
  IsEmpty <- which(is.na(MiscList[[1]][,1,1]))
   
  # yrsmth not implemented here - TO BE ADDED 
  if (length(IsEmpty) > 1) {
    LenMatrix <- t(Data@CAL[x, IsEmpty, ])
  } else { 
    LenMatrix <- (Data@CAL[x, IsEmpty, ])
  }
  
  if (TotYears == 1) LenMatrix <- (as.matrix(LenMatrix))
  binWidth <- Data@CAL_bins[2] - Data@CAL_bins[1]
  CAL_binsmid <- seq(from=Data@CAL_bins[1]+0.5*binWidth, by=binWidth, length=length(Data@CAL_bins)-1) 
  LenMatrix <- cbind(CAL_binsmid, LenMatrix)
  
  if (CurrYear == LHYear) {
    ll <- length(IsEmpty)
    LenMatrix <- LenMatrix[,c(1, (ll-lstyrs+2):(ll+1))] 
    IsEmpty <- (ll-lstyrs+1):(ll)
    MiscList[[1]][1:(min(IsEmpty)-1),1,] <- 0 
    MiscList[[1]][1:(min(IsEmpty)-1),3,] <- 0 
  }
 
  Wb <- Data@wlb[x]
  if (is.na(Wb)) Wb <- 3
  WbCV <- Data@CV_wlb[x]
  if (is.na(WbCV)) WbCV <- 0.1 
  options(warn=-1)
  out <- sapply(1:reps, function(X) {
    LBpars <- new("LB_pars", verbose=FALSE)
    LBpars@MK <- (trlnorm(1,Data@Mort[x],Data@CV_Mort[x])/trlnorm(1,Data@vbK[x],Data@CV_vbK[x]))
    LBpars@Linf <- trlnorm(1, Data@vbLinf[x], Data@CV_vbLinf[x])
    LBpars@CVLinf <- rep(Data@LenCV[x], 1)
    LBpars@L50 <- trlnorm(1, Data@L50[x],  Data@CV_L50[x])
    LBpars@L95 <- trlnorm(1, Data@L95[x],  Data@CV_L50[x])
    LBpars@FecB <- trlnorm(1, Wb, WbCV)
    
    if (LBpars@L50 >= LBpars@Linf) LBpars@L50 <- 0.9 * LBpars@Linf
    if (LBpars@L95 >= LBpars@Linf) LBpars@L95 <- 0.95 * LBpars@Linf
    if (LBpars@L50 >= LBpars@L95) LBpars@L95 <- 1.05 * LBpars@L50
    
    # matplot(LenDat@LMids, LenDat@LData, type="l")
    # matplot(LenDat@LMids, LBfit@pLCatch * apply(LenDat@LData, 2, sum), type="l", add=TRUE, lwd=2)
    
    LenDat <- new("LB_lengths", LenMatrix, LBpars, dataType="freq")
    LenDat@Years <- 1:length(LenDat@Years)
  
    LBfit <- LBSPR::LBSPRfit(LBpars, LenDat, verbose=FALSE)
    # LBfit <- LBSPR::LBSPRfit(LBpars, LenDat, verbose=FALSE, Control=list(modtype="absel"))
    data.frame(rawSPR=LBfit@SPR, rawFM=LBfit@FM)
  })
  options(warn=0)
	
  # Raw estimates 
  MiscList[[1]][IsEmpty,1,] <- matrix(unlist(out[1,]), nrow=length(IsEmpty), ncol=reps)
  MiscList[[1]][IsEmpty,3,] <- matrix(unlist(out[2,]), nrow=length(IsEmpty), ncol=reps) 
    
  # Smoothed estimates - SPR 
  if (reps > 1) {
    MiscList[[1]][,2,] <- apply(MiscList[[1]][,1,], 2, LBSPR::FilterSmooth) 
    # Smoothed estimates - FM
    MiscList[[1]][,4,] <- apply(MiscList[[1]][,3,], 2, LBSPR::FilterSmooth)
  } else {
    MiscList[[1]][,2,] <- LBSPR::FilterSmooth(MiscList[[1]][,1,])
	MiscList[[1]][,4,] <- LBSPR::FilterSmooth(MiscList[[1]][,3,])
  }
    
  return(MiscList)
}

