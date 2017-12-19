


applyMP <- function(Data, MPs = NA, reps = 100) {
  if (class(Data) != "Data") stop("First argument must be object of class 'Data'", call.=FALSE)
  Data <- updateMSE(Data)
  nsims <- length(Data@Mort)
  nMPs <- length(MPs)

  
  if (.hasSlot(Data, "nareas")) {
    nareas <- Data@nareas   
  } else {
    nareas <- 2 
  }
  returnList <- list() # a list nMPs long containing MPs recommendations
  recList <- list() # a list containing nsim recommendations from a single MP 
  
  if (!sfIsRunning() | (nMPs < 8 & nsims < 8)) {
    for (mp in 1:nMPs) {
      temp <- sapply(1:nsims, MPs[mp], Data = Data, reps = reps)  
      slots <- slotNames(temp[[1]])
      for (X in slots) { # sequence along recommendation slots 
        if (X == "Misc") { # convert to a list nsim by nareas
          rec <- lapply(temp, slot, name=X)
        } else {
          rec <- do.call("cbind", lapply(temp, slot, name=X)) # unlist(lapply(temp, slot, name=X))
        }
        if (X == "Spatial") { # convert to a matrix nsim by nareas
          rec <- matrix(rec, nareas, nsims, byrow=FALSE)   
        }
        recList[[X]] <- rec
        for (x in 1:nsims) Data@Misc[[x]] <- recList$Misc[[x]]
        recList$Misc <- NULL
      }
      returnList[[mp]] <- recList
 
      
    }
  } else {
    for (mp in 1:nMPs) {
      temp <- sfSapply(1:nsims, MPs[mp], Data = Data, reps = reps)  
      slots <- slotNames(temp[[1]])
      for (X in slots) { # sequence along recommendation slots 
        if (X == "Misc") { # convert to a list nsim by nareas
          rec <- lapply(temp, slot, name=X)
        } else {
          rec <- do.call("cbind", lapply(temp, slot, name=X)) # unlist(lapply(temp, slot, name=X))
        }
        if (X == "Spatial") { # convert to a matrix nsim by nareas
          rec <- matrix(rec, nareas, nsims, byrow=FALSE)  
        }
        
        recList[[X]] <- rec
        for (x in 1:nsims) Data@Misc[[x]] <- recList$Misc[[x]]
        recList$Misc <- NULL
      }
      returnList[[mp]] <- recList

    }
  }
  
  Data@MPs <- MPs
  
  list(returnList, Data)
}


CalcMPDynamics <- function(MPRecs, y, nyears, proyears, nsim,
                           LastEffort, LastSpatial, LastAllocat, LastCatch,
                           TACused, maxF,
                           LR5_P, LFR_P, Rmaxlen_P, retL_P, retA_P,
                           L5_P, LFS_P, Vmaxlen_P, SLarray_P, V_P,
                           Fdisc_P, DR_P,
                           M_ageArray, FM_P, FM_Pret, Z_P, CB_P, CB_Pret,
                           TAC_f, E_f, SizeLim_f,
                           VBiomass_P, Biomass_P, FinF, Spat_targ,
                           CAL_binsmid, Linf, Len_age, maxage, nareas, Asize, nCALbins,
                           qs, qvar, qinc) {
  # Change in Effort 
  if (length(MPRecs$Effort) == 0) { # no effort recommendation
    if (y==1) Ei <- LastEffort  * E_f[,y] # effort is unchanged but has implementation error
    if (y>1) Ei <- LastEffort / E_f[,y-1]  * E_f[,y] # effort is unchanged but has implementation error
  } else if (length(MPRecs$Effort) != nsim) {
    stop("Effort recommmendation is not 'nsim' long.\n Does MP return Effort recommendation under all conditions?")
  } else {
    Ei <- MPRecs$Effort * E_f[,y] # effort adjustment with implementation error
  }
  
  # Spatial 
  if (all(is.na(MPRecs$Spatial))) { # no spatial recommendation 
    Si <- LastSpatial # spatial is unchanged 
  } else if (any(is.na(MPRecs$Spatial))) {
    stop("Spatial recommmendation has some NAs.\n Does MP return Spatial recommendation under all conditions?")
  } else {
    Si <- MPRecs$Spatial # change spatial fishing
  }
  
  # Allocation 
  if (length(MPRecs$Allocate) == 0) { # no allocation recommendation
    Ai <- LastAllocat # allocation is unchanged 
  } else if (length(MPRecs$Allocate) != nsim) {
    stop("Allocate recommmendation is not 'nsim' long.\n Does MP return Allocate recommendation under all conditions?")
  } else {
    Ai <- MPRecs$Allocate # change in spatial allocation
  }
  Ai <- as.numeric(Ai)
  
  # Retention Curve
  RetentFlag <- FALSE # should retention curve be updated for future years?
  # LR5 
  if (length(MPRecs$LR5) == 0) { # no  recommendation
    LR5_P[(y + nyears):(nyears+proyears),] <- matrix(LR5_P[y + nyears-1,], 
                                                     nrow=(length((y + nyears):(nyears+proyears))),
                                                     ncol=nsim, byrow=TRUE) # unchanged 
    
  } else if (length(MPRecs$LR5) != nsim) {
    stop("LR5 recommmendation is not 'nsim' long.\n Does MP return LR5 recommendation under all conditions?")
  } else {
    LR5_P[(y + nyears):(nyears+proyears),] <- matrix(MPRecs$LR5 * SizeLim_f[,y], 
                                                     nrow=(length((y + nyears):(nyears+proyears))),
                                                     ncol=nsim, byrow=TRUE) # recommendation with implementation error
    RetentFlag <- TRUE
  }
  # LFR 
  if (length(MPRecs$LFR) == 0) { # no  recommendation
    LFR_P[(y + nyears):(nyears+proyears),] <- matrix(LFR_P[y + nyears-1,], 
                                                     nrow=(length((y + nyears):(nyears+proyears))),
                                                     ncol=nsim, byrow=TRUE) # unchanged 
  } else if (length(MPRecs$LFR) != nsim) {
    stop("LFR recommmendation is not 'nsim' long.\n Does MP return LFR recommendation under all conditions?")
  } else {
    LFR_P[(y + nyears):(nyears+proyears),] <- matrix(MPRecs$LFR * SizeLim_f[,y], 
                                                     nrow=(length((y + nyears):(nyears+proyears))),
                                                     ncol=nsim, byrow=TRUE) # recommendation with implementation error
    RetentFlag <- TRUE
  }
  # Rmaxlen 
  if (length(MPRecs$Rmaxlen) == 0) { # no  recommendation
    Rmaxlen_P[(y + nyears):(nyears+proyears),] <- matrix(Rmaxlen_P[y + nyears-1,], 
                                                         nrow=(length((y + nyears):(nyears+proyears))),
                                                         ncol=nsim, byrow=TRUE)   # unchanged 
    
  } else if (length(Rmaxlen) != nsim) {
    stop("Rmaxlen recommmendation is not 'nsim' long.\n Does MP return Rmaxlen recommendation under all conditions?")
  } else {
    Rmaxlen_P[(y + nyears):(nyears+proyears),] <- matrix(MPRecs$Rmaxlen, 
                                                         nrow=(length((y + nyears):(nyears+proyears))),
                                                         ncol=nsim, byrow=TRUE) # recommendation
    RetentFlag <- TRUE
  }
  
  
  # HS - harvest slot 
  if (length(MPRecs$HS) == 0) { # no  recommendation
    HS <- rep(1E5, nsim) # no harvest slot 
  } else if (length(MPRecs$HS) != nsim) {
    stop("HS recommmendation is not 'nsim' long.\n Does MP return HS recommendation under all conditions?")
  } else {
    HS <- MPRecs$HS  * SizeLim_f[,y] # recommendation
    RetentFlag <- TRUE
  }
  
  # Selectivity Curve
  SelectFlag <- FALSE # has selectivity been updated?
  # L5 
  if (length(MPRecs$L5) == 0) { # no  recommendation
    L5_P[(y + nyears):(nyears+proyears),] <- matrix(L5_P[y + nyears-1,], 
                                                    nrow=(length((y + nyears):(nyears+proyears))),
                                                    ncol=nsim, byrow=TRUE) # unchanged 
    
  } else if (length(MPRecs$L5) != nsim) {
    stop("L5 recommmendation is not 'nsim' long.\n Does MP return L5 recommendation under all conditions?")
  } else {
    L5_P[(y + nyears):(nyears+proyears),] <- matrix(MPRecs$L5 * SizeLim_f[,y], 
                                                    nrow=(length((y + nyears):(nyears+proyears))),
                                                    ncol=nsim, byrow=TRUE) # recommendation with implementation error
    SelectFlag <- TRUE
  }
  # LFS
  if (length(MPRecs$LFS) == 0) { # no  recommendation
    LFS_P[(y + nyears):(nyears+proyears),] <- matrix(LFS_P[y + nyears-1,], 
                                                     nrow=(length((y + nyears):(nyears+proyears))),
                                                     ncol=nsim, byrow=TRUE) # unchanged 
  } else if (length(MPRecs$LFS) != nsim) {
    stop("LFS recommmendation is not 'nsim' long.\n Does MP return LFS recommendation under all conditions?")
  } else {
    LFS_P[(y + nyears):(nyears+proyears),] <- matrix(MPRecs$LFS * SizeLim_f[,y], 
                                                     nrow=(length((y + nyears):(nyears+proyears))),
                                                     ncol=nsim, byrow=TRUE) # recommendation with implementation error
    SelectFlag <- TRUE
  }
  # Vmaxlen 
  if (length(MPRecs$Rmaxlen) == 0) { # no  recommendation
    Vmaxlen_P[(y + nyears):(nyears+proyears),] <- matrix(Vmaxlen_P[y + nyears-1,], 
                                                         nrow=(length((y + nyears):(nyears+proyears))),
                                                         ncol=nsim, byrow=TRUE)   # unchanged 
    
  } else if (length(MPRecs$Rmaxlen) != nsim) {
    stop("Rmaxlen recommmendation is not 'nsim' long.\n Does MP return Rmaxlen recommendation under all conditions?")
  } else {
    Vmaxlen_P[(y + nyears):(nyears+proyears),] <- matrix(MPRecs$Vmaxlen, 
                                                         nrow=(length((y + nyears):(nyears+proyears))),
                                                         ncol=nsim, byrow=TRUE) # recommendation
    SelectFlag <- TRUE
  }
  
  # Discard Mortality 
  if (length(MPRecs$Fdisc) >0) { # Fdisc has changed
    if (length(MPRecs$Fdisc) != nsim) stop("Fdisc recommmendation is not 'nsim' long.\n Does MP return Fdisc recommendation under all conditions?")
    Fdisc_P <- MPRecs$Fdisc
  }
  
  # Discard Ratio 
  if (length(MPRecs$DR)>0) { # DR has changed
    if (length(MPRecs$DR) != nsim) stop("DR recommmendation is not 'nsim' long.\n Does MP return DR recommendation under all conditions?")
    DR_P[(y+nyears):(nyears+proyears),] <- matrix(MPRecs$DR, nrow=length((y+nyears):(nyears+proyears)), ncol=nsim, byrow=TRUE) 
  }
  
  # Update Selectivity and Retention Curve 
  if (SelectFlag | RetentFlag) {
    yr <- y+nyears 
    allyrs <- (y+nyears):(nyears+proyears)  # update vulnerabilty for all future years
    
    srs <- (Linf - LFS_P[yr,]) / ((-log(Vmaxlen_P[yr,],2))^0.5) # descending limb
    sls <- (LFS_P[yr,] - L5_P[yr,]) / ((-log(0.05,2))^0.5) # ascending limb
    
    CAL_binsmidMat <- matrix(CAL_binsmid, nrow=nsim, ncol=length(CAL_binsmid), byrow=TRUE)
    selLen <- t(sapply(1:nsim, getsel, lens=CAL_binsmidMat, lfs=LFS_P[yr,], sls=sls, srs=srs))
    
    for (yy in allyrs) {
      # calculate new selectivity at age curve 
      V_P[ , , yy] <- t(sapply(1:nsim, getsel, lens=Len_age[,,yy], lfs=LFS_P[yy,], sls=sls, srs=srs))
      
      # calculate new selectivity at length curve 
      SLarray_P[,, yy] <- selLen  
    }
    
    # sim <- 2
    # plot(CAL_binsmid, selLen[sim,], type="b")
    # lines(c(L5_P[yr,sim], L5_P[yr,sim]), c(0, 0.05), lty=2)
    # lines(c(LFS_P[yr,sim], LFS_P[yr,sim]), c(0, 1), lty=2)
    # lines(c(Linf[sim], Linf[sim]), c(0, Vmaxlen_P[yr,sim]), lty=2)
    
    # calculate new retention curve
    yr <- y+nyears 
    allyrs <- (y+nyears):(nyears+proyears)  # update vulnerabilty for all future years
    
    srs <- (Linf - LFR_P[yr,]) / ((-log(Rmaxlen_P[yr,],2))^0.5) # selectivity parameters are constant for all years
    sls <- (LFR_P[yr,] - LR5_P[yr,]) / ((-log(0.05,2))^0.5)
    
    CAL_binsmidMat <- matrix(CAL_binsmid, nrow=nsim, ncol=length(CAL_binsmid), byrow=TRUE)
    relLen <- t(sapply(1:nsim, getsel, lens=CAL_binsmidMat, lfs=LFR_P[yr,], sls=sls, srs=srs))
    
    for (yy in allyrs) {
      # calculate new retention at age curve 
      retA_P[ , , yy] <- t(sapply(1:nsim, getsel, lens=Len_age[,,yy], lfs=LFR_P[yy,], sls=sls, srs=srs))
      
      # calculate new retention at length curve 
      retL_P[,, yy] <- relLen  
    }
    
    # upper harvest slot 
    aboveHS <- Len_age[,,allyrs]>HS
    tretA_P <- retA_P[,,allyrs]
    tretA_P[aboveHS] <- 0
    retA_P[,,allyrs] <- tretA_P
    for (ss in 1:nsim) {
      index <- which(CAL_binsmid >= HS[ss])
      retL_P[ss, index, allyrs] <- 0
    }	
    
    dr <- aperm(abind::abind(rep(list(DR_P), maxage), along=3), c(2,3,1))
    retA_P[,,allyrs] <- (1-dr[,,yr]) * retA_P[,,yr]
    dr <- aperm(abind::abind(rep(list(DR_P), nCALbins), along=3), c(2,3,1))
    retL_P[,,allyrs] <- (1-dr[,,yr]) * retL_P[,,yr]
    
    # update realized vulnerablity curve with retention and dead discarded fish 
    Fdisc_array1 <- array(Fdisc_P, dim=c(nsim, maxage, length(allyrs)))
    
    V_P[,,allyrs] <- V_P[,,allyrs] * (retA_P[,,allyrs] + (1-retA_P[,,allyrs])*Fdisc_array1)
    
    Fdisc_array2 <- array(Fdisc_P, dim=c(nsim, nCALbins, length(allyrs)))
    SLarray_P[,,allyrs]  <- SLarray_P[,,allyrs] * (retL_P[,,allyrs]+ (1-retL_P[,,allyrs])*Fdisc_array2)
    
    # Realised Retention curves
    retA_P[,,allyrs] <- retA_P[,,allyrs] * V_P[,,allyrs]
    retL_P[,,allyrs] <- retL_P[,,allyrs] * SLarray_P[,,allyrs] 
  }
  

  # indices 
  SAYRL <- as.matrix(expand.grid(1:nsim, 1:maxage, nyears, 1:nareas))  # Final historical year
  SAYRt <- as.matrix(expand.grid(1:nsim, 1:maxage, y + nyears, 1:nareas))  # Trajectory year
  SAYR <- as.matrix(expand.grid(1:nsim, 1:maxage, y, 1:nareas))
  SYt <- SAYRt[, c(1, 3)]
  SAYt <- SAYRt[, 1:3]
  SR <- SAYR[, c(1, 4)]
  SA1 <- SAYR[, 1:2]
  S1 <- SAYR[, 1]
  SY1 <- SAYR[, c(1, 3)]
  SAY1 <- SAYR[, 1:3]
  SYA <- as.matrix(expand.grid(1:nsim, 1, 1:maxage))  # Projection year
  SY <- SYA[, 1:2]
  SA <- SYA[, c(1, 3)]
  SAY <- SYA[, c(1, 3, 2)]
  S <- SYA[, 1]
  
  # update vulnerable biomass for selectivitity curve 
  VBiomass_P[,,y,] <- Biomass_P[, , y, ] * V_P[SAYt] # update vulnerable biomass
  newVB <- apply(Biomass_P[, , y, ] * V_P[SAYt], c(1, 3), sum)  # calculate total vuln biomass by area 
  fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ, 1, sum)  # spatial preference according to spatial vulnerable biomass
  
  d1 <- t(Si) * fishdist  # distribution of fishing effort
  fracE <- apply(d1, 1, sum) # fraction of current effort in open areas
  fracE2 <- d1 * (fracE + (1-fracE) * Ai)/fracE # re-distribution of fishing effort 
  
  fishdist <- fracE2 # fishing effort by area
  
  # Apply TAC recommendation
  if (all(is.na(TACused))) { # no TAC has been set
    
    # fishing mortality with effort control recommendation 
    
    FM_P[SAYR] <- (FinF[S1] * Ei[S1] * V_P[SAYt] * t(Si)[SR] * fishdist[SR] * 
                     qvar[SY1] * (qs[S1]*(1 + qinc[S1]/100)^y))/Asize[SR]
    
    # retained fishing mortality with effort control recommendation
    FM_Pret[SAYR] <- (FinF[S1] * Ei[S1] * retA_P[SAYt] * t(Si)[SR] * fishdist[SR] *
                        qvar[SY1] * qs[S1]*(1 + qinc[S1]/100)^y)/Asize[SR]
    
    Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt] # calculate total mortality 
    
    CB_P[SAYR] <- FM_P[SAYR]/Z_P[SAYR] * Biomass_P[SAYR] * (1 - exp(-Z_P[SAYR]))
    CB_Pret[SAYR] <- FM_Pret[SAYR]/Z_P[SAYR] * Biomass_P[SAYR] * (1 - exp(-Z_P[SAYR]))
    
    Effort <- Ei 
    
  } else { # A TAC has been set
    TACused[is.na(TACused)] <- LastCatch[is.na(TACused)] # if MP returns NA - TAC is set to catch from last year
    TACrec <- TACused             # TAC recommendation
    TACusedE<- TAC_f[,y]*TACused   # TAC taken after implementation error
    
    availB <- apply(newVB * t(Si), 1, sum)
    
    maxC <- (1 - exp(-maxF)) * availB # maximum catch given maxF
    TACusedE[TACusedE > maxC] <- maxC[TACusedE > maxC] # apply maxF limit - catch can't be higher than maxF * vulnerable biomass
    
    CB_P[SAYR] <- (Biomass_P[SAYR] * V_P[SAYt] * fishdist[SR])/Asize[SR] # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space
    # calculate distribution of retained effort 
    CB_Pret[SAYR] <- (Biomass_P[SAYR] * retA_P[SAYt] * fishdist[SR])/Asize[SR]  # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space
    
    retained <- apply(CB_Pret[,,y,], 1, sum)
    actualremovals <- apply(CB_P[,,y,], 1, sum)
    ratio <- actualremovals/retained # ratio of actual removals to retained catch 
    
    temp <- CB_Pret[, , y, ]/apply(CB_Pret[, , y, ], 1, sum) # distribution of retained fish
    CB_Pret[, , y, ] <- TACusedE * temp  # retained catch 
    
    temp <- CB_P[, , y, ]/apply(CB_P[, , y, ], 1, sum) # distribution of removals
    CB_P[,,y,] <- TACusedE *  ratio * temp # scale up total removals 
    
    temp <- CB_P[SAYR]/(Biomass_P[SAYR] * exp(-M_ageArray[SAYt]/2))  # Pope's approximation
    temp[temp > (1 - exp(-maxF))] <- 1 - exp(-maxF)
    
    FM_P[SAYR] <- -log(1 - temp)
    
    temp <- CB_Pret[SAYR]/(Biomass_P[SAYR] * exp(-M_ageArray[SAYt]/2))  # Pope's approximation
    temp[temp > (1 - exp(-maxF))] <- 1 - exp(-maxF)
    
    FM_Pret[SAYR] <- -log(1 - temp)
    
    Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt] # calculate total mortality 
    
    
    # Fishing mortality
    # Fs <- (-log(1 - apply(CB_P[, , y, ], 1, sum)/(apply(CB_P[, , y, ], 1, sum) + 
    #                                                 apply(VBiomass_P[, , y, ], 1, sum))))
    Fs <- (-log(1 - apply(CB_P[, , y, ], 1, sum)/apply(VBiomass_P[, , y, ], 1, sum)))
    
    Effort <-  Fs/(qs * qvar[,y] * ((1 + qinc/100)^y))
    
    # Make sure Effort doesn't exceed regulated effort 
    if (length(MPRecs$Effort) >0 ) { # an effort regulation also exists
      aboveE <- which(Effort > Ei)
      if (length(aboveE)>0) {
        Effort[aboveE] <- Ei[aboveE]
        SAYR <- as.matrix(expand.grid(aboveE, 1:maxage, y, 1:nareas))
        SAYRt <- as.matrix(expand.grid(aboveE, 1:maxage, y + nyears, 1:nareas))  # Trajectory year
        SYt <- SAYRt[, c(1, 3)]
        SAYt <- SAYRt[, 1:3]
        SR <- SAYR[, c(1, 4)]
        S1 <- SAYR[, 1]
        SY1 <- SAYR[, c(1, 3)]
        FM_P[SAYR] <- (FinF[S1] * Ei[S1] * V_P[SAYt] * Si[SR] * fishdist[SR] * qvar[SY1] * 
                         (qs[S1]*(1 + qinc[S1]/100)^y))/Asize[SR]
        
        # retained fishing mortality with input control recommendation
        FM_Pret[SAYR] <- (FinF[S1] * Ei[S1] * retA_P[SAYt] * Si[SR] * fishdist[SR] * 
                            qvar[SY1] * qs[S1]*(1 + qinc[S1]/100)^y)/Asize[SR]
        
        Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt] # calculate total mortality 
        
        CB_P[SAYR] <- FM_P[SAYR]/Z_P[SAYR] * Biomass_P[SAYR] * (1 - exp(-Z_P[SAYR]))
        CB_Pret[SAYR] <- FM_retain[SAYR]/Z_P[SAYR] * Biomass_P[SAYR] * (1 - exp(-Z_P[SAYR]))
      }
      
    }
    
  }
  
  # Returns
  out <- list()
  out$TACrec <- TACused
  out$V_P <- V_P
  out$SLarray_P <- SLarray_P
  out$retA_P <- retA_P
  out$retL_P <- retL_P
  out$Fdisc_P <- Fdisc_P
  out$VBiomass_ <- VBiomass_P
  out$Z_P <- Z_P
  out$FM_P <- FM_P
  out$FM_Pret <- FM_Pret
  out$CB_P <- CB_P
  out$CB_Pret <- CB_Pret
  out$Si <- Si
  out$Ai <- Ai
  out$Ei <- Ei
  out$Effort <- Effort
  out
}
