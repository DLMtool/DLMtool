

makeData <- function(Biomass, CBret, Cret, N, SSB, VBiomass, StockPars, 
                     FleetPars, ObsPars, ImpPars, RefPoints,
                     ErrList, OM, SampCpars, initD, silent=FALSE) {
  
  if(!silent) message("Simulating observed data")
  
  Name <- OM@Name
  nyears <- OM@nyears
  proyears <- OM@proyears
  nsim <- OM@nsim 
  nareas <- StockPars$nareas
  reps <- OM@reps
  
  Data <- new("Data", stock = "MSE")  # create a blank DLM data object
  if (reps == 1) Data <- OneRep(Data)  # make stochastic variables certain for only one rep
  Data <- replic8(Data, nsim)  # make nsim sized slots in the DLM data object
  
  Data@Name <- Name
  Data@Year <- 1:nyears
  
  # --- Observed catch ----
  # Simulated observed retained catch (biomass)
  Cobs <- ErrList$Cbiasa[, 1:nyears] * ErrList$Cerr[, 1:nyears] * apply(CBret, c(1, 3), sum)  
  Data@Cat <- Cobs 
  
  # --- Index of total abundance ----
  # Index of abundance from total biomass - beginning of year before fishing
  # apply hyperstability / hyperdepletion
  II <- (apply(Biomass, c(1, 3), sum)^ObsPars$betas) * ErrList$Ierr[, 1:nyears]  
  II <- II/apply(II, 1, mean)  # normalize
  Data@Ind <- II # index of total abundance
  
  # --- Index of recruitment ----
  Data@Rec <- apply(N[, 1, , ], c(1, 2), sum) * ErrList$Recerr[, 1:nyears] 
  Data@t <- rep(nyears, nsim) # number of years of data
  
  # --- Average catch ----
  Data@AvC <- apply(Cobs, 1, mean) # average catch over all years
  
  # --- Depletion ----
  # observed depletion
  Depletion <- apply(SSB[,,nyears,],1,sum)/RefPoints$SSB0 # current depletion
  Data@Dt <- ObsPars$Dbias * Depletion * 
    rlnorm(nsim, mconv(1, ObsPars$Derr), sdconv(1, ObsPars$Derr))
  
  Data@Dep <- ObsPars$Dbias * Depletion * 
    rlnorm(nsim, mconv(1, ObsPars$Derr), sdconv(1, ObsPars$Derr))  
  
  # --- Life-history parameters ----
  Data@vbLinf <- StockPars$Linfarray[,nyears] * ObsPars$Linfbias # observed vB Linf
  Data@vbK <- StockPars$Karray[,nyears] * ObsPars$Kbias # observed vB K
  Data@vbt0 <- StockPars$t0array[,nyears] * ObsPars$t0bias # observed vB t0
  Data@Mort <- StockPars$Marray[,nyears] * ObsPars$Mbias # natural mortality
  Data@L50 <- StockPars$L50array[,nyears] * ObsPars$lenMbias # observed length at 50% maturity
  Data@L95 <- StockPars$L95array[,nyears] * ObsPars$lenMbias # observed length at 95% maturity
  Data@L95[Data@L95 > 0.9 * Data@vbLinf] <- 0.9 * Data@vbLinf[Data@L95 > 0.9 * Data@vbLinf]  # Set a hard limit on ratio of L95 to Linf
  Data@L50[Data@L50 > 0.9 * Data@L95] <- 0.9 * Data@L95[Data@L50 > 0.9 * Data@L95]  # Set a hard limit on ratio of L95 to Linf
  Data@LenCV <- StockPars$LenCV # variablity in length-at-age - no error at this time
  Data@sigmaR <-  StockPars$procsd # observed sigmaR - assumed no obs error
  Data@MaxAge <- StockPars$maxage # maximum age - no error - used for setting up matrices only
  
  # Observed steepness values 
  hs <- StockPars$hs
  if (!is.null(OM@cpars[['hsim']])) {
    hsim <- SampCpars$hsim
    hbias <- hsim/hs  # back calculate the simulated bias
    if (OM@hbiascv == 0) hbias <- rep(1, nsim) 
    ObsPars$hbias <- hbias 
  } else {
    hsim <- rep(NA, nsim)  
    cond <- hs > 0.6
    hsim[cond] <- 0.2 + rbeta(sum(hs > 0.6), alphaconv((hs[cond] - 0.2)/0.8, (1 - (hs[cond] - 0.2)/0.8) * OM@hbiascv), 
                              betaconv((hs[cond] - 0.2)/0.8,  (1 - (hs[cond] - 0.2)/0.8) * OM@hbiascv)) * 0.8
    hsim[!cond] <- 0.2 + rbeta(sum(hs <= 0.6), alphaconv((hs[!cond] - 0.2)/0.8,  (hs[!cond] - 0.2)/0.8 * OM@hbiascv), 
                               betaconv((hs[!cond] - 0.2)/0.8, (hs[!cond] - 0.2)/0.8 * OM@hbiascv)) * 0.8
    hbias <- hsim/hs  # back calculate the simulated bias
    if (OM@hbiascv == 0) hbias <- rep(1, nsim) 
    ObsPars$hbias <- hbias
  }
  Data@steep <- hs * ObsPars$hbias # observed steepness
  
  # --- Reference points ----
  # Simulate observation error in BMSY/B0 
  ntest <- 20  # number of trials  
  BMSY_B0bias <- array(rlnorm(nsim * ntest, 
                              mconv(1, OM@BMSY_B0biascv), sdconv(1, OM@BMSY_B0biascv)), 
                       dim = c(nsim, ntest))  # trial samples of BMSY relative to unfished  
  test <- array(RefPoints$SSBMSY_SSB0 * BMSY_B0bias, dim = c(nsim, ntest))  # the simulated observed BMSY_B0 
  indy <- array(rep(1:ntest, each = nsim), c(nsim, ntest))  # index
  indy[test > max(0.9, max(RefPoints$SSBMSY_SSB0))] <- NA  # interval censor
  BMSY_B0bias <- BMSY_B0bias[cbind(1:nsim, apply(indy, 1, min, na.rm = T))]  # sample such that BMSY_B0<90%
  ObsPars$BMSY_B0bias <- BMSY_B0bias
  
  Data@FMSY_M <- RefPoints$FMSY_M * ObsPars$FMSY_Mbias # observed FMSY/M
  Data@BMSY_B0 <- RefPoints$SSBMSY_SSB0 * ObsPars$BMSY_B0bias # observed BMSY/B0
  Data@Cref <- RefPoints$MSY * ObsPars$Crefbias # Catch reference - MSY with error
  Data@Bref <- RefPoints$VBMSY * ObsPars$Brefbias # Vuln biomass ref - VBMSY with error
  
  # Generate values for reference SBMSY/SB0
  # should be calculated from unfished - won't be correct if initD is set
  I3 <- apply(Biomass, c(1, 3), sum)^ObsPars$betas  # apply hyperstability / hyperdepletion
  I3 <- I3/apply(I3, 1, mean)  # normalize index to mean 1
  if (!is.null(initD)) {
    b1 <- apply(Biomass, c(1, 3), sum)
    b2 <- matrix(RefPoints$BMSY, nrow=nsim, ncol=nyears)
    ind <- apply(abs(b1/ b2 - 1), 1, which.min) # find years closest to BMSY
    Iref <- diag(I3[1:nsim,ind])  # return the real target abundance index closest to BMSY
  } else {
    Iref <- apply(I3[, 1:5], 1, mean) * RefPoints$SSBMSY_SSB0  # return the real target abundance index corresponding to BMSY
  }
  Data@Iref <- Iref * ObsPars$Irefbias # index reference with error
  
  # --- Abundance ----
  # Calculate vulnerable and spawning biomass abundance --
  M_array <- array(0.5*StockPars$M_ageArray[,,nyears], dim=c(nsim, StockPars$maxage, nareas))
  A <- apply(VBiomass[, , nyears, ] * exp(-M_array), 1, sum) # Abundance (mid-year before fishing)
  Asp <- apply(SSB[, , nyears, ] * exp(-M_array), 1, sum)  # Spawning abundance (mid-year before fishing)
  OFLreal <- A * (1-exp(-RefPoints$FMSY))  # the true simulated Over Fishing Limit
  
  Data@Abun <- A * ObsPars$Abias * 
    rlnorm(nsim, mconv(1, ObsPars$Aerr), sdconv(1, ObsPars$Aerr)) # observed vulnerable abundance
  Data@SpAbun <- Asp * ObsPars$Abias * 
    rlnorm(nsim, mconv(1, ObsPars$Aerr), sdconv(1, ObsPars$Aerr)) # spawing abundance
  
  # --- Catch-at-age ----
  Data@CAA <- simCAA(nsim, nyears, StockPars$maxage, Cret, ObsPars$CAA_ESS, ObsPars$CAA_nsamp) 
  
  # --- Catch-at-length ----
  vn <- apply(N, c(1,2,3), sum) * FleetPars$retA[,,1:nyears] # numbers at age in population that would be retained
  vn <- aperm(vn, c(1,3, 2))
  
  CALdat <- simCAL(nsim, nyears, StockPars$maxage, ObsPars$CAL_ESS, 
                   ObsPars$CAL_nsamp, StockPars$nCALbins, StockPars$CAL_binsmid, 
                   vn, FleetPars$retL, StockPars$Linfarray, 
                   StockPars$Karray, StockPars$t0array, StockPars$LenCV)
  
  Data@CAL_bins <- StockPars$CAL_bins
  Data@CAL <- CALdat$CAL # observed catch-at-length
  Data@ML <- CALdat$ML # mean length
  Data@Lc <- CALdat$Lc # modal length 
  Data@Lbar <- CALdat$Lbar # mean length above Lc 
  
  Data@LFC <- CALdat$LFC * ObsPars$LFCbias # length at first capture
  Data@LFS <- FleetPars$LFS[nyears,] * ObsPars$LFSbias # length at full selection
  
  # --- Previous Management Recommendations ----
  Data@MPrec <- apply(CBret, c(1, 3), sum) # catch in last year
  Data@MPeff <- rep(1, nsim) # effort in last year = 1 
  
  # --- Store OM Parameters ----
  # put all the operating model parameters in one table
  ind <- which(lapply(StockPars, length) == nsim)
  stock <- as.data.frame(StockPars[ind])
  stock$Fdisc <- NULL
  ind <- which(lapply(FleetPars, length) == nsim)
  fleet <- as.data.frame(FleetPars[ind])
  
  ind <- which(lapply(ImpPars, length) == nsim)
  imp <- as.data.frame(ImpPars[ind])
  refs <- RefPoints %>% select('MSY', 'FMSY', 'SSBMSY_SSB0', 'BMSY_B0',
                               'UMSY', 'FMSY_M', 'RefY', 'Blow', 'MGT', 'SSB0')
  
  OMtable <- data.frame(stock, fleet, imp, refs, ageM=StockPars$ageM[,nyears], 
                     L5=FleetPars$L5[nyears, ], LFS=FleetPars$LFS[nyears, ], 
                     Vmaxlen=FleetPars$Vmaxlen[nyears, ],
                     LR5=FleetPars$LR5[nyears,], LFR=FleetPars$LFR[nyears,], 
                     Rmaxlen=FleetPars$Rmaxlen[nyears,], 
                     DR=FleetPars$DR[nyears,], OFLreal, maxF=OM@maxF)
                 
  OMtable <- OMtable[,order(names(OMtable))]
  Data@OM <- OMtable
  
  # --- Store Obs Parameters ----
  ObsTable <- as.data.frame(ObsPars)
  ObsTable <- ObsTable[,order(names(ObsTable))]
  Data@Obs <- ObsTable # put all the observation error model parameters in one table
  
  # --- Misc ----
  Data@Units <- "unitless"
  Data@Ref_type <- "Simulated OFL"
  Data@wla <- rep(StockPars$a, nsim)
  Data@wlb <- rep(StockPars$b, nsim)
  Data@nareas <- nareas
  Data@Ref <- OFLreal 
  Data@LHYear <- nyears  # Last historical year is nyears (for fixed MPs)
  Data@Misc <- vector("list", nsim)
  
  Data
}


updateData <- function(Data, OM, MPCalcs, Effort, Biomass, Biomass_P, CB_Pret, 
                       N_P, SSB, SSB_P, VBiomass, VBiomass_P, RefPoints, ErrList, 
                       FMSY_P, retA_P, 
                       retL_P, StockPars, FleetPars, ObsPars, 
                       upyrs, interval, y=2, 
                       mm=1, Misc, SampCpars) {
  
  yind <- upyrs[match(y, upyrs) - 1]:(upyrs[match(y, upyrs)] - 1) # index
  
  nyears <- OM@nyears
  proyears <- OM@proyears
  nsim <- OM@nsim 
  nareas <- StockPars$nareas
  reps <- OM@reps
  
  Data@Year <- 1:(nyears + y - 1)
  Data@t <- rep(nyears + y, nsim)
  
  # --- Simulate catches ---- 
  CBtemp <- CB_Pret[, , yind, , drop=FALSE] # retained catch-at-age
  CNtemp <- retA_P[,,yind+nyears, drop=FALSE] * 
    apply(N_P[,,yind,, drop=FALSE], c(1,2,3), sum) # retained age structure
  CBtemp[is.na(CBtemp)] <- tiny
  CBtemp[!is.finite(CBtemp)] <- tiny
  CNtemp[is.na(CNtemp)] <- tiny
  CNtemp[!is.finite(CNtemp)] <- tiny
  CNtemp <- aperm(CNtemp, c(1,3,2))
  
  # --- Observed catch ----
  # Simulated observed retained catch (biomass)
  Cobs <- ErrList$Cbiasa[, nyears + yind] * ErrList$Cerr[, nyears + yind] * 
    apply(CBtemp, c(1, 3), sum, na.rm = TRUE)
  Data@Cat <- cbind(Data@Cat, Cobs) 
  
  # --- Index of total abundance ----
  I2 <- (cbind(apply(Biomass, c(1, 3), sum), 
               apply(Biomass_P, c(1, 3), sum)[, 1:(y - 1)])^ObsPars$betas) * 
    ErrList$Ierr[, 1:(nyears + (y - 1))]
  
  I2[is.na(I2)] <- tiny
  I2 <- I2/apply(I2, 1, mean)
  Data@Ind <- I2
  
  # --- Update real indices (if they exist) ----
  if (!all(is.na(Data@RInd))) { # Real indices exist

    b1 <- cbind(apply(Biomass, c(1, 3), sum), apply(Biomass_P, c(1, 3), sum)[, 1:(y - 1)])
    b1 <- (exp(lcs(b1))^filter(ErrList$stats.df, Index=="Biomass")$beta) * ErrList$RIerr$BIerr[, 1:(nyears+y-1)]
    
    b2 <- cbind(apply(VBiomass, c(1, 3), sum), apply(VBiomass_P, c(1, 3), sum)[, 1:(y - 1)])
    b2 <- (exp(lcs(b2))^filter(ErrList$stats.df, Index=="VBiomass")$beta) * ErrList$RIerr$VBIerr[, 1:(nyears+y-1)]
    
    b3 <- cbind(apply(SSB, c(1, 3), sum), apply(SSB_P, c(1, 3), sum)[, 1:(y - 1)])
    b3 <- (exp(lcs(b3))^filter(ErrList$stats.df, Index=="SpBiomass")$beta) * ErrList$RIerr$SBIerr[, 1:(nyears+y-1)]

    sim.indices <- list(Biomass=b1, VBiomass=b2, SpBiomass=b3)
    old <- Data@RInd
    Data@RInd <- array(NA, dim=c(nsim, length(Data@Type), nyears+y-1))
    for (type in Data@Type) {
      yrs <- min(nyears, length(old[1, match(type,Data@Type),]))
      bio <- (sim.indices[[type]]/apply(sim.indices[[type]][,1:yrs], 1, mean)) *
        mean(old[1, match(type,Data@Type),1:yrs])
      Data@RInd[,match(type,Data@Type), ] <- array(bio, dim=c(nsim, 1, ncol(bio)))
      
      # Replace generated data with real data if more than nyears are provided
      Data.temp <- SampCpars$Data@RInd[1,match(type,Data@Type),] 
      Data.temp <- Data.temp[!is.na(Data.temp)]
      Data.yr <- length(Data.temp)
      Data.yr <- min(Data.yr, nyears-1+y)
      temp <- array(SampCpars$Data@RInd[1,match(type,Data@Type),1:Data.yr], dim=c(1, Data.yr, nsim))
      temp <- aperm(temp, c(3,1,2))
      Data@RInd[,match(type,Data@Type), 1:Data.yr] <- temp
    }
  }
  
  # --- Index of recruitment ----
  Recobs <- ErrList$Recerr[, nyears + yind] * apply(array(N_P[, 1, yind, ], 
                                                          c(nsim, interval[mm], nareas)),
                                                    c(1, 2), sum)
  Data@Rec <- cbind(Data@Rec, Recobs)
  
  # --- Average catch ----
  Data@AvC <- apply(Data@Cat, 1, mean)
  
  # --- Depletion ----
  Depletion <- apply(SSB_P[, , y, ], 1, sum)/RefPoints$SSB0 
  Depletion[Depletion < tiny] <- tiny
  Data@Dt <- ObsPars$Dbias * Depletion * rlnorm(nsim, mconv(1, ObsPars$Derr), sdconv(1, ObsPars$Derr))
  Data@Dep <- ObsPars$Dbias * Depletion * rlnorm(nsim, mconv(1, ObsPars$Derr), sdconv(1, ObsPars$Derr))
  
  # --- Update life-history parameter estimates for current year ----
  Data@vbLinf <- StockPars$Linfarray[,nyears+y] * ObsPars$Linfbias # observed vB Linf
  Data@vbK <- StockPars$Karray[,nyears+y] * ObsPars$Kbias # observed vB K
  Data@vbt0 <- StockPars$t0array[,nyears+y] * ObsPars$t0bias # observed vB t0
  Data@Mort <- StockPars$Marray[,nyears+y] * ObsPars$Mbias # natural mortality
  Data@L50 <- StockPars$L50array[,nyears+y] * ObsPars$lenMbias # observed length at 50% maturity
  Data@L95 <- StockPars$L95array[,nyears+y] * ObsPars$lenMbias # observed length at 95% maturity
  Data@L95[Data@L95 > 0.9 * Data@vbLinf] <- 0.9 * Data@vbLinf[Data@L95 > 0.9 * Data@vbLinf]  # Set a hard limit on ratio of L95 to Linf
  Data@L50[Data@L50 > 0.9 * Data@L95] <- 0.9 * Data@L95[Data@L50 > 0.9 * Data@L95]  # Set a hard limit on ratio of L95 to Linf
  
  
  # --- Abundance ----
  # Calculate vulnerable and spawning biomass abundance --
  M_array <- array(0.5*StockPars$M_ageArray[,,nyears+y], dim=c(nsim, StockPars$maxage, nareas))
  A <- apply(VBiomass_P[, , y, ] * exp(-M_array), 1, sum) # Abundance (mid-year before fishing)
  Asp <- apply(SSB_P[, , y, ] * exp(-M_array), 1, sum)  # Spawning abundance (mid-year before fishing)
  Data@Abun <- A * ObsPars$Abias * rlnorm(nsim, mconv(1, ObsPars$Aerr), sdconv(1, ObsPars$Aerr))
  Data@SpAbun <- Asp * ObsPars$Abias * rlnorm(nsim, mconv(1, ObsPars$Aerr), sdconv(1, ObsPars$Aerr))
  Data@Ref <- A * (1 - exp(-FMSY_P[,mm,y])) 
  
  # --- Catch-at-age ----
  # previous CAA
  oldCAA <- Data@CAA
  Data@CAA <- array(0, dim = c(nsim, nyears + y - 1, StockPars$maxage))
  Data@CAA[, 1:(nyears + y - interval[mm] - 1), ] <- oldCAA[, 1:(nyears + y - interval[mm] - 1), ] 
  # update CAA
  CAA <- simCAA(nsim, yrs=length(yind), StockPars$maxage, Cret=CNtemp, ObsPars$CAA_ESS, ObsPars$CAA_nsamp)
  Data@CAA[, nyears + yind, ] <- CAA
  
  

  # --- Catch-at-length ----
  oldCAL <- Data@CAL
  Data@CAL <- array(0, dim = c(nsim, nyears + y - 1, StockPars$nCALbins))
  Data@CAL[, 1:(nyears + y - interval[mm] - 1), ] <- oldCAL[, 1:(nyears + y - interval[mm] - 1), ]
  
  CAL <- array(NA, dim = c(nsim, interval[mm], StockPars$nCALbins))  
  vn <- (apply(N_P[,,,], c(1,2,3), sum) * retA_P[,,(nyears+1):(nyears+proyears)]) # numbers at age that would be retained
  vn <- aperm(vn, c(1,3,2))
  
  CALdat <- simCAL(nsim, nyears=length(yind), StockPars$maxage, ObsPars$CAL_ESS, 
                   ObsPars$CAL_nsamp, StockPars$nCALbins, StockPars$CAL_binsmid, 
                   vn=vn[,yind,, drop=FALSE], retL=retL_P[,,nyears+yind, drop=FALSE],
                   Linfarray=StockPars$Linfarray[,nyears + yind, drop=FALSE],  
                   Karray=StockPars$Karray[,nyears + yind, drop=FALSE], 
                   t0array=StockPars$t0array[,nyears + yind,drop=FALSE],
                   LenCV=StockPars$LenCV)

  Data@CAL[, nyears + yind, ] <- CALdat$CAL # observed catch-at-length
  Data@ML <- cbind(Data@ML, CALdat$ML) # mean length
  Data@Lc <- cbind(Data@Lc, CALdat$Lc) # modal length 
  Data@Lbar <- cbind(Data@Lbar, CALdat$Lbar) # mean length above Lc 
  
  Data@LFC <- CALdat$LFC * ObsPars$LFCbias # length at first capture
  Data@LFS <- FleetPars$LFS[nyears+y,] * ObsPars$LFSbias # length at full selection

  # --- Previous Management Recommendations ----
  Data@MPrec <- MPCalcs$TACrec # last MP  TAC recommendation
  Data@MPeff <- Effort[, mm, y-1] # last recommended effort
  
  Data@Misc <- Misc
  
  Data
}


