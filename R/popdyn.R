
CalcMPDynamics <- function(MPRecs, y, nyears, proyears, nsim,
                           LastEi, LastSpatial, LastAllocat, LastCatch,
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
    if (y==1) Ei <- LastEi * E_f[,y] # effort is unchanged but has implementation error
    if (y>1) Ei <- LastEi / E_f[,y-1]  * E_f[,y] # effort is unchanged but has implementation error
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
  if (all(dim(Si) != c(nareas, nsim))) stop("Spatial recommmendation not nareas long")
  
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
    
  } else if (length(MPRecs$Rmaxlen) != nsim) {
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
    srs[!is.finite(srs)] <- Inf
    sls <- (LFS_P[yr,] - L5_P[yr,]) / ((-log(0.05,2))^0.5) # ascending limb
    
    CAL_binsmidMat <- matrix(CAL_binsmid, nrow=nsim, ncol=length(CAL_binsmid), byrow=TRUE)
    selLen <- t(sapply(1:nsim, getsel, lens=CAL_binsmidMat, lfs=LFS_P[yr,], sls=sls, srs=srs))
    
    for (yy in allyrs) {
      # calculate new selectivity at age curve 
      V_P[ , , yy] <- t(sapply(1:nsim, getsel, lens=Len_age[,,yy], lfs=LFS_P[yy,], sls=sls, srs=srs))
      
      # calculate new selectivity at length curve 
      SLarray_P[,, yy] <- selLen  
    }
    
    # sim <- 158
    # plot(CAL_binsmid, selLen[sim,], type="b")
    # lines(c(L5_P[yr,sim], L5_P[yr,sim]), c(0, 0.05), lty=2)
    # lines(c(LFS_P[yr,sim], LFS_P[yr,sim]), c(0, 1), lty=2)
    # lines(c(Linf[sim], Linf[sim]), c(0, Vmaxlen_P[yr,sim]), lty=2)
    
    # calculate new retention curve
    yr <- y+nyears 
    allyrs <- (y+nyears):(nyears+proyears)  # update vulnerabilty for all future years
    
    srs <- (Linf - LFR_P[yr,]) / ((-log(Rmaxlen_P[yr,],2))^0.5) # selectivity parameters are constant for all years
    srs[!is.finite(srs)] <- Inf
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
    aboveHS <- Len_age[,,allyrs, drop=FALSE]>array(HS, dim=c(nsim, maxage, length(allyrs)))
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
    
    V_P[,,allyrs] <- V_P[,,allyrs, drop=FALSE] * (retA_P[,,allyrs, drop=FALSE] + (1-retA_P[,,allyrs, drop=FALSE])*Fdisc_array1)
    
    Fdisc_array2 <- array(Fdisc_P, dim=c(nsim, nCALbins, length(allyrs)))
    SLarray_P[,,allyrs]  <- SLarray_P[,,allyrs, drop=FALSE] * (retL_P[,,allyrs, drop=FALSE]+ (1-retL_P[,,allyrs, drop=FALSE])*Fdisc_array2)
    
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
  VBiomass_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt] # update vulnerable biomass
  
  # Calculate fishing distribution if all areas were open 
  newVB <- apply(VBiomass_P[,,y,], c(1,3), sum) # calculate total vuln biomass by area 
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
    
    # M_age_area <- array(M_ageArray[,,y], dim=c(nsim, maxage, nareas))
    # Fs <- suppressWarnings(-log(1 - apply(CB_P[, , y, ], 1, sum)/apply(VBiomass_P[, , y, ]*exp(-(0.5*M_age_area)), 1, sum))) # Pope's approx
    # Fs[!is.finite(Fs)] <- 2  # NaN for very high Fs
    
    Effort <- FinF *Ei * apply(fracE2, 1, sum) # (Fs/qs)/ FinF # Ei  # (Fs/qs)/ FinF change in catchability not included in effort calc: * qvar[,y] * ((1 + qinc/100)^y)) 
    
  } else { # A TAC has been set
    TACused[is.na(TACused)] <- LastCatch[is.na(TACused)] # if MP returns NA - TAC is set to catch from last year
    TACrec <- TACused             # TAC recommendation
    TACusedE<- TAC_f[,y]*TACused   # TAC taken after implementation error
    
    availB <- apply(newVB * t(Si), 1, sum)
    
    # maxC <- (1 - exp(-maxF)) * availB # maximum catch given maxF
    # TACusedE[TACusedE > maxC] <- maxC[TACusedE > maxC] # apply maxF limit - catch can't be higher than maxF * vulnerable biomass
    
    CB_P[SAYR] <- (Biomass_P[SAYR] * V_P[SAYt] * fishdist[SR])/Asize[SR] # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space
    # calculate distribution of retained effort 
    CB_Pret[SAYR] <- (Biomass_P[SAYR] * retA_P[SAYt] * fishdist[SR])/Asize[SR]  # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space
    
    retained <- apply(CB_Pret[,,y,], 1, sum)
    actualremovals <- apply(CB_P[,,y,], 1, sum)
    ratio <- actualremovals/retained # ratio of actual removals to retained catch 
    ratio[!is.finite(ratio)] <- 0 
    ratio[ratio>1E5] <- 1E5
    temp <- CB_Pret[, , y, ]/apply(CB_Pret[, , y, ], 1, sum) # distribution of retained fish
    CB_Pret[, , y, ] <- TACusedE * temp  # retained catch 
    
    temp <- CB_P[, , y, ]/apply(CB_P[, , y, ], 1, sum) # distribution of removals
    
    CB_P[,,y,] <- TACusedE *  ratio * temp # scale up total removals 
    
    chk <- apply(CB_P[,,y,], 1, sum) > availB # total removals can't be more than available biomass
    if (sum(chk)>0) {
      c_temp <- apply(CB_P[chk,,y,, drop=FALSE], 1, sum)
      ratio_temp <- (availB[chk]/c_temp) * 0.99
      if (sum(chk)>1) CB_P[chk,,y, ] <- CB_P[chk,,y,] * array(ratio_temp, dim=c(sum(chk), maxage, nareas))
      if (sum(chk)==1) CB_P[chk,,y, ] <- CB_P[chk,,y,] * array(ratio_temp, dim=c(maxage, nareas))
    }
  
    # total removals
    # t1 <- apply(CB_P[,,y,],1, sum)
  
    temp <- CB_P[SAYR]/(Biomass_P[SAYR] * exp(-M_ageArray[SAYt]/2))  # Pope's approximation
    temp[temp > (1 - exp(-maxF))] <- 1 - exp(-maxF) # apply maxF constraint
    FM_P[SAYR] <- -log(1 - temp)
    Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt] # calculate total mortality 
    # update removals with maxF constraint
    CB_P[SAYR] <- FM_P[SAYR]/Z_P[SAYR] * Biomass_P[SAYR] * (1 - exp(-Z_P[SAYR])) 

    # t2 <- apply(CB_P[,,y,],1, sum)
   
   
    # repeated because of approximation error in Pope's approximation - an issue if CB_P ~ AvailB
    chk <- apply(CB_P[,,y,], 1, sum) > availB # total removals can't be more than available biomass
    
    if (sum(chk)>0) {
      c_temp <- apply(CB_P[chk,,y,, drop=FALSE], 1, sum)
      ratio_temp <- (availB[chk]/c_temp) * 0.99
      if (sum(chk)>1) CB_P[chk,,y, ] <- CB_P[chk,,y,] * array(ratio_temp, dim=c(sum(chk), maxage, nareas))
      if (sum(chk)==1) CB_P[chk,,y, ] <- CB_P[chk,,y,] * array(ratio_temp, dim=c(maxage, nareas))
    }
  
    # retained catch
    temp <- CB_Pret[SAYR]/(Biomass_P[SAYR] * exp(-M_ageArray[SAYt]/2))  # Pope's approximation
    temp[temp > (1 - exp(-maxF))] <- 1 - exp(-maxF) # apply maxF constraint
    FM_Pret[SAYR] <- -log(1 - temp)
    # update catch with maxF constraint
    CB_Pret[SAYR] <- FM_Pret[SAYR]/Z_P[SAYR] * Biomass_P[SAYR] * (1 - exp(-Z_P[SAYR])) 
    
    M_age_area <- array(M_ageArray[,,y], dim=c(nsim, maxage, nareas))
  
    Fs <- suppressWarnings(-log(1 - apply(CB_P[, , y, ], 1, sum)/apply(VBiomass_P[, , y, ]*exp(-(0.5*M_age_area)), 1, sum))) # Pope's approx
    Fs[!is.finite(Fs)] <- 2  # NaN for very high Fs
   
    Effort <- Fs/(FinF * qs*qvar[,y]* (1 + qinc/100)^y) * apply(fracE2, 1, sum)  

    # Make sure Effort doesn't exceed regulated effort
    if (length(MPRecs$Effort) >0 | all(LastEi != 1)) { # an effort regulation also exists
      aboveE <- which(Effort > Ei)
      if (length(aboveE)>0) {
        Effort[aboveE] <- Ei[aboveE] * FinF[aboveE] * apply(fracE2, 1, sum)[aboveE]
        SAYR <- as.matrix(expand.grid(aboveE, 1:maxage, y, 1:nareas))
        SAYRt <- as.matrix(expand.grid(aboveE, 1:maxage, y + nyears, 1:nareas))  # Trajectory year
        SYt <- SAYRt[, c(1, 3)]
        SAYt <- SAYRt[, 1:3]
        SR <- SAYR[, c(1, 4)]
        S1 <- SAYR[, 1]
        SY1 <- SAYR[, c(1, 3)]
        FM_P[SAYR] <- (FinF[S1] * Ei[S1] * V_P[SAYt] * t(Si)[SR] * fishdist[SR] * qvar[SY1] * 
                         (qs[S1]*(1 + qinc[S1]/100)^y))/Asize[SR]
        
        # retained fishing mortality with input control recommendation
        FM_Pret[SAYR] <- (FinF[S1] * Ei[S1] * retA_P[SAYt] * t(Si)[SR] * fishdist[SR] * 
                            qvar[SY1] * qs[S1]*(1 + qinc[S1]/100)^y)/Asize[SR]
        
        Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt] # calculate total mortality 
        
        CB_P[SAYR] <- FM_P[SAYR]/Z_P[SAYR] * Biomass_P[SAYR] * (1 - exp(-Z_P[SAYR]))
        CB_Pret[SAYR] <- FM_Pret[SAYR]/Z_P[SAYR] * Biomass_P[SAYR] * (1 - exp(-Z_P[SAYR]))
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


# getMSY <- function(x, MatAge, LenAge, WtAge, MatureAge,  VAge, maxage, R0, SRrel, hs) {
#   
#   opt <- optimize(MSYCalcs, log(c(0.001, 10)), MatAge=MatAge[x,], LenAge=LenAge[x,], 
#                   WtAge=WtAge[x,], MatureAge=MatureAge[x,], 
#                   VAge=VAge[x,], maxage=maxage, R0=R0[x], SRrel=SRrel[x], hs=hs[x], opt=1) 
#   
#   
#   runMod <- MSYCalcs(logapicF=opt$minimum, MatAge=MatAge[x,], LenAge=LenAge[x,], 
#                      WtAge=WtAge[x,], MatureAge=MatureAge[x,], 
#                      VAge=VAge[x,], maxage=maxage, R0=R0[x], SRrel=SRrel[x], hs=hs[x], opt=2) 
#   
#   runMod
# }
# 

MSYCalcs <- function(logapicF, MatAge, WtAge, MatureAge, VAge, maxage, R0, SRrel, hs, opt=1) {
  # Box 3.1 Walters & Martell 2004
  apicF <- exp(logapicF)
  lx <- l0 <- rep(1, maxage)
  for (a in 2:maxage) {
    l0[a] <- l0[a-1] * exp(-MatAge[a-1])
    lx[a] <- lx[a-1] * exp(-(MatAge[a-1] + apicF*VAge[a-1]))
  }
  Egg0 <- sum(l0 * WtAge * MatureAge) # unfished egg production (assuming fecundity proportional to weight)
  EggF <- sum(lx * WtAge * MatureAge) # fished egg production (assuming fecundity proportional to weight)

  vB0 <- sum(l0 * WtAge * VAge)
  vBF <- sum(lx * WtAge * VAge)

  SB0 <- sum(l0 * WtAge * MatureAge) # same as eggs atm
  SBF <- sum(lx * WtAge * MatureAge)

  B0 <- sum(l0 * WtAge) 
  BF <- sum(lx * WtAge)

  hs[hs>0.999] <- 0.999
  recK <- (4*hs)/(1-hs) # Goodyear compensation ratio
  reca <- recK/Egg0
  if (SRrel ==1) {
    recb <- (reca * Egg0 - 1)/(R0*Egg0) # BH SRR
    RelRec <- (reca * EggF-1)/(recb*EggF)
  }
  if (SRrel ==2) {
    bR <- (log(5*hs)/(0.8*SB0))
    aR <- exp(bR*SB0)/(SB0/R0)
    RelRec <- (log(aR*EggF/R0))/(bR*EggF/R0)
  }

  RelRec[RelRec<0] <- 0
  
  Fa <- apicF*VAge
  Za <- Fa + MatAge
  relyield <- Fa/Za * lx * (1-exp(-Za)) * WtAge
  YPR <- sum(relyield)
  Yield <- YPR * RelRec
  
  if (opt == 1)  return(-Yield)
  if (opt == 2) {
    out <- c(Yield=Yield,
             F=-log(1 - (Yield/(vBF*RelRec+Yield))),
             SB = SBF * RelRec,
             SB_SB0 = (SBF * RelRec)/(SB0 * R0),
             B_B0 = (BF * RelRec + Yield)/(B0 * R0),
             B = BF * RelRec + Yield,
             VB = vBF * RelRec + Yield,
             VB_VB0 = (vBF * RelRec + Yield)/(vB0 * R0),
             RelRec=RelRec,
             SB0 = SB0 * R0,
             B0=B0 * R0,
             apicF=apicF)

    return(out)
  }


}

optMSY_eq <- function(x, M_ageArray, Wt_age, Mat_age, V, maxage, R0, SRrel, hs, yr=1) {
  bounds <- c(0.0000001, 5)
  doopt <- optimise(MSYCalcs, log(bounds), MatAge=M_ageArray[x,,yr], WtAge=Wt_age[x,,yr], 
                    MatureAge=Mat_age[x,,yr], VAge=V[x,,yr], maxage, R0=R0[x], SRrel=SRrel[x], hs=hs[x], opt=1)
  
  apicFMSY <- exp(doopt$minimum)
  apicFMSY2 <- apicFMSY
  
  MSYs <- MSYCalcs(log(apicFMSY), MatAge=M_ageArray[x,,yr], WtAge=Wt_age[x,,yr], 
           MatureAge=Mat_age[x,,yr], VAge=V[x,,yr], maxage, R0=R0[x], SRrel=SRrel[x], hs=hs[x], opt=2)
  if (MSYs[1] < 1) {
    count <- 0; stop <- FALSE
    while (apicFMSY > 0.95 * max(bounds) & count < 50 & !stop) {
      count <- count + 1
      bounds <- c(0.0000001, max(bounds)-0.1)
      if (bounds[1] < bounds[2]) {
        doopt <- optimise(MSYCalcs, log(bounds), MatAge=M_ageArray[x,,yr], WtAge=Wt_age[x,,yr], 
                          MatureAge=Mat_age[x,,yr], VAge=V[x,,yr], maxage, R0=R0[x], SRrel=SRrel[x], hs=hs[x], opt=1)
        apicFMSY <- exp(doopt$minimum)
      } else {
        stop <- TRUE
      }
    }
    if (count >=50 | stop) apicFMSY <- apicFMSY2
    MSYs <- MSYCalcs(log(apicFMSY), MatAge=M_ageArray[x,,yr], WtAge=Wt_age[x,,yr], 
                     MatureAge=Mat_age[x,,yr], VAge=V[x,,yr], maxage, R0=R0[x], SRrel=SRrel[x], hs=hs[x], opt=2)
  }
  return(MSYs)
  
}


#' optimize for catchability (q)
#' 
#' Function optimizes catchability (q, where F=qE) required to get to user-specified stock
#' depletion
#'
#' @param x Integer, the simulation number
#' @param D A numeric vector nsim long of sampled depletion
#' @param SSB0 A numeric vector nsim long of total unfished spawning biomass
#' @param nareas The number of spatial areas
#' @param maxage The maximum age
#' @param N Array of the numbers-at-age in population. Dimensions are nsim, maxage, nyears, nareas. 
#' Only values from the first year (i.e `N[,,1,]`) are used, which is the current N-at-age.
#' @param pyears The number of years to project forward. Equal to 'nyears' for optimizing for q.
#' @param M_ageArray An array (dimensions nsim, maxage, nyears+proyears) with the natural mortality-at-age and year 
#' @param Mat_age An array (dimensions nsim, maxage, proyears+nyears) with the proportion mature for each age-class
#' @param Asize A matrix (dimensions nsim, nareas) with size of each area
#' @param Wt_age An array (dimensions nsim, maxage, nyears+proyears) with the weight-at-age and year 
#' @param V An array (dimensions nsim, maxage, nyears+proyears) with the vulnerability-at-age and year
#' @param retA An array (dimensions nsim, maxage, nyears+proyears) with the probability retained-at-age and year
#' @param Perr A matrix (dimensions nsim, nyears+proyears) with the recruitment deviations
#' @param mov An array (dimensions nsim, nareas, nareas) with the movement matrix
#' @param SRrel A numeric vector nsim long specifying the recruitment curve to use
#' @param Find A matrix (dimensions nsim, nyears) with the historical fishing effort 
#' @param Spat_targ A numeric vector nsim long with the spatial targeting
#' @param hs A numeric vector nsim long with the steepness values for each simulation
#' @param R0a A matrix (dimensions nsim, nareas) with the unfished recruitment by area
#' @param SSBpR A matrix (dimensions nsim, nareas) with the unfished spawning-per-recruit by area
#' @param aR A numeric vector nareas long with the Ricker SRR a values
#' @param bR A numeric vector nareas long with the Ricker SRR b values
#' @param bounds A numeric vector of length 2 with bounds for the optimizer
#' @param maxF A numeric value specifying the maximum fishing mortality for any single age class
#' @param MPA A matrix of spatial closures by year
#' @param useCPP logical - use the CPP code? For testing purposes only
#' @author A. Hordyk
#' @keywords internal
getq3 <- function(x, D, SSB0, nareas, maxage, N, pyears, M_ageArray, Mat_age, Asize, Wt_age,
                  V, retA, Perr, mov, SRrel, Find, Spat_targ, hs, R0a, SSBpR, aR, bR, 
                  bounds = c(1e-05, 15), maxF, MPA, useCPP=TRUE) {
  
  opt <- optimize(optQ, log(bounds), depc=D[x], SSB0c=SSB0[x], nareas, maxage, Ncurr=N[x,,1,], 
                  pyears, M_age=M_ageArray[x,,], MatAge=Mat_age[x,,], Asize_c=Asize[x,], WtAge=Wt_age[x,,],
                  Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,], movc=mov[x,,,], SRrelc=SRrel[x], 
                  Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
                  SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], maxF=maxF, MPA=MPA, useCPP=useCPP)
  return(exp(opt$minimum))
}


#' Optimize q for a single simulation 
#'
#' @param logQ log q
#' @param depc Depletion value
#' @param SSB0c Unfished spawning biomass
#' @param nareas Number of areas
#' @param maxage Maximum age
#' @param Ncurr Current N-at-age
#' @param pyears Number of years to project population dynamics
#' @param M_age M-at-age
#' @param Asize_c Numeric vector (length nareas) with size of each area
#' @param MatAge Maturity-at-age
#' @param WtAge Weight-at-age
#' @param Vuln Vulnerability-at-age
#' @param Retc Retention-at-age
#' @param Prec Recruitment error by year
#' @param movc movement matrix
#' @param SRrelc SR parameter
#' @param Effind Historical fishing effort
#' @param Spat_targc Spatial targetting
#' @param hc Steepness
#' @param R0c Unfished recruitment by area
#' @param SSBpRc Unfished spawning biomass per recruit by area
#' @param aRc Ricker aR
#' @param bRc Ricker bR
#' @param maxF maximum F
#' @param MPA A matrix of spatial closures by year
#' @param useCPP Logical. Use the CPP code?
#' @author A. Hordyk
#' @keywords internal

optQ <- function(logQ, depc, SSB0c, nareas, maxage, Ncurr, pyears, M_age, Asize_c,
                 MatAge, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc, 
                 R0c, SSBpRc, aRc, bRc, maxF, MPA, useCPP) {
  if (!useCPP) {
    # simpop <- popdyn(nareas, maxage, Ncurr, pyears, M_age, Asize_c,
    #                  MatAge, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc, 
    #                  R0c=R0c, SSBpRc=SSBpRc, aRc=aRc, bRc=bRc, Qc=exp(logQ), maxF=maxF, MPA=MPA, control=1) 
    # ssb <- sum(simpop$SBarray[,pyears,]) # doesn't currently work with age-based movement
    
  } else {
    simpop <- popdynCPP(nareas, maxage, Ncurr, pyears, M_age, Asize_c,
                        MatAge, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc, 
                        R0c=R0c, SSBpRc=SSBpRc, aRc=aRc, bRc=bRc, Qc=exp(logQ), Fapic=0, 
                        maxF=maxF, MPA=MPA, control=1,  SSB0c=SSB0c) 
  
    ssb <- sum(simpop[[4]][,pyears,])
  }
  
  (log(depc) - log(ssb/SSB0c))^2
  
}


# #' Population dynamics model
# #'
# #' @param nareas Integer. The number of spatial areas
# #' @param maxage Integer. The maximum age
# #' @param Ncurr Numeric matrix (dimensions maxage, nareas) with the current N-at-age
# #' @param pyears Integer. Number of years to project the model forward
# #' @param M_age Numeric matrix (dimensions maxage, pyears) with natural mortality at age
# #' @param Asize_c Numeric vector (length nareas) with size of each area
# #' @param MatAge Numeric matrix (dimensions maxage, nyears+proyears) with proportion mature for each age-class
# #' @param WtAge Numeric matrix (dimensions maxage, pyears) with weight-at-age 
# #' @param Vuln Numeric matrix (dimensions maxage, pyears) with proportion vulnerable-at-age
# #' @param Retc Numeric matrix (dimensions maxage, pyears) with proportion retained-at-age
# #' @param Prec Numeric vector (length pyears) with recruitment error
# #' @param movc Numeric matrix (dimensions nareas, nareas) with movement matrix
# #' @param SRrelc Integer. Stock-recruitment curve
# #' @param Effind Numeric vector (length pyears) with the fishing effort by year
# #' @param Spat_targc Integer. Value of spatial targetting
# #' @param hc Numeric. Steepness of stock-recruit relationship
# #' @param R0c Numeric vector of length nareas with unfished recruitment by area
# #' @param SSBpRc Numeric vector of length nareas with unfished spawning per recruit by area
# #' @param aRc Numeric. Ricker SRR a value
# #' @param bRc Numeric. Ricker SRR b value
# #' @param Qc Numeric. Catchability coefficient
# #' @param Fapic Numeric. Apical F value
# #' @param maxF A numeric value specifying the maximum fishing mortality for any single age class
# #' @param MPA A matrix of spatial closures by year
# #' @param control Integer. 1 to use q and effort to calculate F, 2 to use Fapic (apical F) and 
# #' vulnerablity to calculate F.
# #' 
# #' @author A. Hordyk
# #'
# #' @return A named list of length 8 containing with arrays (dimensions: maxage, pyears, nareas)
# #' containing numbers-at-age, biomass-at-age, spawning stock numbers, spawning biomass, 
# #' vulnerable biomass, fishing mortality, retained fishing mortality, and total mortality
# # #' @export
# #'
# popdyn <- function(nareas, maxage, Ncurr, pyears, M_age, Asize_c,
#                    MatAge, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc, 
#                    R0c, SSBpRc, aRc, bRc, Qc, Fapic=NULL, maxF, MPA, control=1) {
#   Narray <- array(NA, dim=c(maxage, pyears, nareas))
#   Barray <- array(NA, dim=c(maxage, pyears, nareas))
#   SSNarray <- array(NA, dim=c(maxage, pyears, nareas))
#   SBarray <- array(NA, dim=c(maxage, pyears, nareas))
#   VBarray <- array(NA, dim=c(maxage, pyears, nareas))
#   Marray <- array(NA, dim=c(maxage, pyears, nareas))
#   FMarray <- array(NA, dim=c(maxage, pyears, nareas))
#   FMretarray <- array(NA, dim=c(maxage, pyears, nareas))
#   Zarray <- array(NA, dim=c(maxage, pyears, nareas))
#   
#   Narray[,1,] <- Ncurr
#   Barray[,1,] <- Narray[,1,] * WtAge[,1]
#   SSNarray[,1,] <- Ncurr * MatAge[,1] # spawning stock numbers
#   SBarray[,1,] <- Narray[,1,] * WtAge[,1] * MatAge[,1] # spawning biomass
#   VBarray[,1,] <- Narray[,1,] * WtAge[,1] * Vuln[,1] # vulnerable biomass
#   Marray[,1,] <- M_age[,1] # M-at-age
#   
#   SAYR <- as.matrix(expand.grid(1:maxage, 1, 1:nareas))  # Set up some array indexes age (A) year (Y) region/area (R)
#   
#   # Distribution of fishing effort 
#   VBa <- colSums(VBarray[,1,]) # total vuln biomass in each area 
#   
#   # fishdist <- VBa^Spat_targc/mean(VBa^Spat_targc)
#   fishdist <- VBa^Spat_targc/sum(VBa^Spat_targc)
#   
#   Asize_mat <- matrix(Asize_c, nrow=maxage, ncol=nareas, byrow=TRUE)
#                     
#   if (control == 1) {
#     FMarray[SAYR] <- (Effind[SAYR[,2]] * Qc * Vuln[SAYR[,1:2]] * fishdist[SAYR[,3]])/Asize_mat
#     FMretarray[SAYR] <- (Effind[SAYR[,2]] * Qc * Retc[SAYR[,1:2]] * fishdist[SAYR[,3]])/Asize_mat
#   }
#   if (control == 2) {
#     FMarray[SAYR] <- (Fapic * Vuln[SAYR[,1:2]] * fishdist[SAYR[,3]])/Asize_mat
#     FMretarray[SAYR] <- (Fapic * Retc[SAYR[,1:2]] * fishdist[SAYR[,3]])/Asize_mat
#   }
#   
#   FMarray[,1,][FMarray[,1,] > (1 - exp(-maxF))] <- 1 - exp(-maxF)
#   FMretarray[,1,][FMretarray[,1,] > (1 - exp(-maxF))] <- 1 - exp(-maxF)
#  
#   Zarray[,1,] <- Marray[,1,] + FMarray[,1,]
# 
#   for (y in 1:(pyears-1)) {
#     
#     NextYrN <- popdynOneTS(nareas, maxage, SSBcurr=colSums(SBarray[,y,]), Ncurr=Narray[,y,], 
#                            Zcurr=Zarray[,y,], PerrYr=Prec[y+maxage+1], hc, R0c, SSBpRc, aRc, bRc, 
#                            movc, SRrelc)
#     
#     Narray[,y+1,] <- NextYrN
#     Barray[,y+1,] <- Narray[,y+1,] * WtAge[,y+1]
#     SSNarray[,y+1,] <- Narray[,y+1,] * MatAge[,y+1] # spawning stock numbers
#     SBarray[,y+1,] <- Narray[,y+1,] * WtAge[,y+1] * MatAge[,y+1] # spawning biomass
#     VBarray[,y+1,] <- Narray[,y+1,] * WtAge[,y+1] * Vuln[,y+1] # vulnerable biomass
#     Marray[, y+1, ] <- M_age[,y+1]
#     
#     # Distribution of fishing effort 
#     VBa <- colSums(VBarray[,y+1,]) # total vuln biomass in each area 
#     # fishdist <- VBa^Spat_targc/mean(VBa^Spat_targc)
#     fishdist <- VBa^Spat_targc/sum(VBa^Spat_targc)
#     
#     d1 <- t(matrix(MPA[y,])) * fishdist  # distribution of fishing effort
#     fracE <- apply(d1, 1, sum) # fraction of current effort in open areas
#     fracE2 <- d1 * (fracE + (1-fracE))/fracE # re-distribution of fishing effort 
#     fishdist <- fracE2 # fishing effort by area
#     
#     
#     SAYR <- as.matrix(expand.grid(1:maxage, y+1, 1:nareas))  # Set up some array indexes age (A) year (Y) region/area (R)
#     if (control ==1) {
#       FMarray[SAYR] <- (Effind[SAYR[,2]] * Qc * Vuln[SAYR[,1:2]] * fishdist[SAYR[,3]])/Asize_mat
#       FMretarray[SAYR] <- (Effind[SAYR[,2]] * Qc * Retc[SAYR[,1:2]] * fishdist[SAYR[,3]])/Asize_mat
#     }
#     if (control ==2) {
#       FMarray[SAYR] <- (Fapic * Vuln[SAYR[,1:2]] * fishdist[SAYR[,3]])/Asize_mat
#       FMretarray[SAYR] <- (Fapic * Retc[SAYR[,1:2]] * fishdist[SAYR[,3]])/Asize_mat
#     }   
#     FMarray[SAYR][FMarray[SAYR] > (1 - exp(-maxF))] <- 1 - exp(-maxF)
#     FMretarray[SAYR][FMretarray[SAYR] > (1 - exp(-maxF))] <- 1 - exp(-maxF)
#     Zarray[,y+1,] <- Marray[,y+1,] + FMarray[,y+1,]
#     
#   }
#   
#   out <- list()
#   out$Narray <- Narray
#   out$Barray <- Barray
#   out$SSNarray <- SSNarray
#   out$SBarray <- SBarray
#   out$VBarray <- VBarray
#   out$FMarray <- FMarray
#   out$FMretarray <- FMretarray
#   out$Zarray <- Zarray
# 
#   out
# }
# 

# #' Population dynamics model for one annual time-step
# #'
# #' Project population forward one time-step given current numbers-at-age and total mortality
# #'
# #' @param nareas The number of spatial areas
# #' @param maxage The maximum age
# #' @param SSBcurr A numeric vector of length nareas with the current spawning biomass in each area
# #' @param Ncurr A numeric matrix (maxage, nareas) with current numbers-at-age in each area
# #' @param Zcurr A numeric matrix (maxage, nareas) with total mortality-at-age in each area
# #' @param PerrYr A numeric value with recruitment deviation for current year
# #' @param hs Steepness of SRR
# #' @param R0c Numeric vector with unfished recruitment by area
# #' @param SSBpRc Numeric vector with unfished spawning stock per recruit by area
# #' @param aRc Numeric vector with Ricker SRR a parameter by area
# #' @param bRc Numeric vector with Ricker SRR b parameter by area
# #' @param movc Numeric matrix (nareas by nareas) with the movement matrix
# #' @param SRrelc Integer indicating the stock-recruitment relationship to use (1 for Beverton-Holt, 2 for Ricker)
# #' @author A. Hordyk
# #'
# # #' @export
# #' @keywords internal
# popdynOneTS <- function(nareas, maxage, SSBcurr, Ncurr, Zcurr,
#                    PerrYr, hc, R0c, SSBpRc, aRc, bRc, movc, SRrelc)  {
# 
#   # set up some indices for indexed calculation
# 
#   indMov <- as.matrix(expand.grid(1:maxage,1:nareas, 1:nareas))  # Movement master index
#   indMov2 <- indMov[, c(1, 2)]  # Movement from index
#   indMov3 <- indMov[, c(2, 3)]  # Movement to index
# 
#   Nnext <- array(NA, dim=c(maxage, nareas))
# 
#   # Recruitment assuming regional R0 and stock wide steepness
#   if (SRrelc[1] == 1) {
#     Nnext[1,  ] <- PerrYr *  (4 * R0c * hc * SSBcurr)/(SSBpRc * R0c * (1-hc) + (5*hc-1)*SSBcurr)
#   } else {
#     # most transparent form of the Ricker uses alpha and beta params
#     Nnext[1,  ] <- PerrYr * aRc * SSBcurr * exp(-bRc * SSBcurr)
#   }
# 
#   # Mortality
#   Nnext[2:maxage, ] <- Ncurr[1:(maxage - 1),  ] * exp(-Zcurr[1:(maxage - 1), ])  # Total mortality
# 
#   # Movement of stock
#   temp <- array(Nnext[indMov2] * movc[indMov3], dim = c(maxage,nareas, nareas))  # Move individuals
#   Nnext <- apply(temp, c(1, 3), sum)
# 
#   # Numbers-at-age at beginning of next year
#   return(Nnext)
# 
# }
# 
# 
# #' Simulate population dynamics for historical years
# #'
# #' @param x Integer, the simulation number 
# #' @param nareas The number of spatial areas
# #' @param maxage The maximum age
# #' @param N Array of the numbers-at-age in population. Dimensions are nsim, maxage, nyears, nareas. 
# #' Only values from the first year (i.e `N[,,1,]`) are used, which is the current N-at-age.
# #' @param pyears The number of years to project forward. Equal to 'nyears' for optimizing for q.
# #' @param M_ageArray An array (dimensions nsim, maxage, nyears+proyears) with the natural mortality-at-age and year
# #' @param Asize A matrix (dimensions nsim, nareas) of size of areas 
# #' @param Mat_age A matrix (dimensions nsim, maxage) with the proportion mature for each age-class
# #' @param Wt_age An array (dimensions nsim, maxage, nyears+proyears) with the weight-at-age and year 
# #' @param V An array (dimensions nsim, maxage, nyears+proyears) with the vulnerability-at-age and year
# #' @param retA An array (dimensions nsim, maxage, nyears+proyears) with the probability retained-at-age and year
# #' @param Perr A matrix (dimensions nsim, nyears+proyears) with the recruitment deviations
# #' @param mov An array (dimensions nsim, nareas, nareas) with the movement matrix
# #' @param SRrel A numeric vector nsim long specifying the recruitment curve to use
# #' @param Find A matrix (dimensions nsim, nyears) with the historical fishing effort 
# #' @param Spat_targ A numeric vector nsim long with the spatial targeting
# #' @param hs A numeric vector nsim long with the steepness values for each simulation
# #' @param R0a A matrix (dimensions nsim, nareas) with the unfished recruitment by area
# #' @param SSBpR A matrix (dimensions nsim, nareas) with the unfished spawning-per-recruit by area
# #' @param aR A numeric vector nsim long with the Ricker SRR a values
# #' @param bR A numeric vector nsim long with the Ricker SRR b values
# #' @param qs A numeric vector nsim long with catchability coefficients
# #' @param MPA A matrix of spatial closures by year
# #' @param maxF A numeric value specifying the maximum fishing mortality for any single age class
# #' @param useCPP logical - use the CPP code? For testing purposes only 
# #' @param SSB0 SSB0
# #' @author A. Hordyk
# #' @keywords internal
# #' @export
# simYears <- function(x, nareas, maxage, N, pyears, M_ageArray, Asize, Mat_age, Wt_age,
#                      V, retA, Perr, mov, SRrel, Find, Spat_targ, hs, R0a, SSBpR, aR, bR, qs, 
#                      MPA, maxF, useCPP=TRUE, SSB0) {
#   if(!useCPP) {
#     # popdyn(nareas, maxage, Ncurr=N[x,,1,], pyears,  
#     #        M_age=M_ageArray[x,,], Asize_c=Asize[x,], MatAge=Mat_age[x,,], WtAge=Wt_age[x,,],
#     #        Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,], movc=mov[x,,,], SRrelc=SRrel[x], 
#     #        Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
#     #        SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], Qc=qs[x], MPA=MPA, maxF=maxF, control=1)
#     # doesn't currently work with age-based movement
#   } else {
#     popdynCPP(nareas, maxage, Ncurr=N[x,,1,], pyears,  
#            M_age=M_ageArray[x,,], Asize_c=Asize[x,], MatAge=Mat_age[x,,], WtAge=Wt_age[x,,],
#            Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,], movc=mov[x,,,], SRrelc=SRrel[x], 
#            Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
#            SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], Qc=qs[x], Fapic=0, MPA=MPA, maxF=maxF, 
#            control=1, SSB0c=SSB0[x])
#   }
#   
# }



# #' Calculate FMSY and related metrics using Rcpp code
# #'
# #' @param x Integer, the simulation number
# #' @param Asize A matrix (nsim by nareas) with size of areas
# #' @param nareas The number of spatial areas
# #' @param maxage The maximum age
# #' @param N Array of the numbers-at-age in population. Dimensions are nsim, maxage, nyears, nareas.
# #' Only values from the first year (i.e `N[,,1,]`) are used, which is the current N-at-age.
# #' @param pyears The number of years to project forward. Equal to 'nyears' for optimizing for q.
# #' @param M_ageArray An array (dimensions nsim, maxage, nyears+proyears) with the natural mortality-at-age and year
# #' @param Mat_age A matrix (dimensions nsim, maxage) with the proportion mature for each age-class
# #' @param Wt_age An array (dimensions nsim, maxage, nyears+proyears) with the weight-at-age and year
# #' @param V An array (dimensions nsim, maxage, nyears+proyears) with the vulnerability-at-age and year
# #' @param retA An array (dimensions nsim, maxage, nyears+proyears) with the probability retained-at-age and year
# #' @param Perr A matrix (dimensions nsim, nyears+proyears) with the recruitment deviations
# #' @param mov An array (dimensions nsim, nareas, nareas) with the movement matrix
# #' @param SRrel A numeric vector nsim long specifying the recruitment curve to use
# #' @param Find A matrix (dimensions nsim, nyears) with the historical fishing effort
# #' @param Spat_targ A numeric vector nsim long with the spatial targeting
# #' @param hs A numeric vector nsim long with the steepness values for each simulation
# #' @param R0a A matrix (dimensions nsim, nareas) with the unfished recruitment by area
# #' @param SSBpR A matrix (dimensions nsim, nareas) with the unfished spawning-per-recruit by area
# #' @param aR A numeric vector nsim long with the Ricker SRR a values
# #' @param bR A numeric vector nsim long with the Ricker SRR b values
# #' @param SSB0 Unfished spawning biomass
# #' @param B0 Unfished total biomass
# #' @param MPA A matrix of spatial closures by year
# #' @param maxF A numeric value specifying the maximum fishing mortality for any single age class
# #' @param useCPP logical - use the CPP code? For testing purposes only
# #'
# #' @author A. Hordyk
# #'
# getFMSY3 <- function(x, Asize, nareas, maxage, N, pyears, M_ageArray, Mat_age, Wt_age,
#                     V, retA, Perr, mov, SRrel, Find, Spat_targ, hs, R0a, SSBpR, aR, bR,
#                     SSB0, B0, MPA, maxF, useCPP=TRUE) {
# 
#  opt <- optimize(optMSY, log(c(0.001, 10)), Asize_c=Asize[x,], nareas, maxage, Ncurr=N[x,,1,],
#                  pyears, M_age=M_ageArray[x,,], MatAge=Mat_age[x,,],
#                  WtAge=Wt_age[x,,], Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,],
#                  movc=mov[x,,,], SRrelc=SRrel[x],
#                  Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,],
#                  SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], MPA=MPA, maxF=maxF, useCPP=useCPP,
#                  SSB0c=SSB0[x])
# 
#  MSY <- -opt$objective
# 
#  if (!useCPP) {
#    # simpop <- popdyn(nareas, maxage, Ncurr=N[x,,1,],
#    #                  pyears, M_age=M_ageArray[x,,], Asize_c=Asize[x,],
#    #                  MatAge=Mat_age[x,,],
#    #                  WtAge=Wt_age[x,,], Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,],
#    #                  movc=mov[x,,,], SRrelc=SRrel[x],
#    #                  Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,],
#    #                  SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], Fapic=exp(opt$minimum), MPA=MPA, maxF=maxF, control=2)
#    # 
#    # # calculate B0 and SSB0 with current conditions
#    # simpopF0 <- popdyn(nareas, maxage, Ncurr=N[x,,1,],
#    #                    pyears, M_age=M_ageArray[x,,], Asize_c=Asize[x,],
#    #                    MatAge=Mat_age[x,,],
#    #                    WtAge=Wt_age[x,,], Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,],
#    #                    movc=mov[x,,,], SRrelc=SRrel[x],
#    #                    Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,],
#    #                    SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], Fapic=0, MPA=MPA, maxF=maxF, control=2)
# 
#  } else {
#    simpop <- popdynCPP(nareas, maxage, Ncurr=N[x,,1,],
#                        pyears, M_age=M_ageArray[x,,], Asize_c=Asize[x,],
#                        MatAge=Mat_age[x,,],
#                        WtAge=Wt_age[x,,], Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,],
#                        movc=mov[x,,,], SRrelc=SRrel[x],
#                        Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,],
#                        SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], Qc=0, Fapic=exp(opt$minimum), 
#                        MPA=MPA, maxF=maxF, control=2, SSB0c = SSB0[x])
#    # calculate B0 and SSB0 with current conditions
#    simpopF0 <- popdynCPP(nareas, maxage, Ncurr=N[x,,1,],
#                          pyears, M_age=M_ageArray[x,,], Asize_c=Asize[x,],
#                          MatAge=Mat_age[x,,],
#                          WtAge=Wt_age[x,,], Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,],
#                          movc=mov[x,,,], SRrelc=SRrel[x],
#                          Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,],
#                          SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], Qc=0, Fapic=0, MPA=MPA, maxF=maxF, 
#                          control=2, SSB0c = SSB0[x])
#  }
# 
# 
#  ## Cn <- simpop[[7]]/simpop[[8]] * simpop[[1]] * (1-exp(-simpop[[8]])) # retained catch
#  Cn <- simpop[[6]]/simpop[[8]] * simpop[[1]] * (1-exp(-simpop[[8]])) # removals
#  Cb <- Cn[,pyears,] * Wt_age[x,,pyears]
# 
#  B <- sum(simpop[[2]][,pyears,] + Cb)
# 
#  SSB_MSY <- sum(simpop[[4]][,pyears,])
# 
#  V_BMSY <- sum(simpop[[5]][,pyears,])
#  F_MSYv <- -log(1 - (MSY/(V_BMSY+MSY)))
# 
# 
#  SSB0_curr <- sum(simpopF0[[4]][,pyears,])
#  B0_curr <- sum(simpopF0[[2]][,pyears,])
#  SSBMSY_SSB0 <- sum(simpop[[4]][,pyears,])/SSB0_curr
#  BMSY_B0 <- sum(simpop[[2]][,pyears,])/B0_curr
#  # SSBMSY_SSB0 <- sum(simpop[[4]][,pyears,])/SSB0[x]
#  # BMSY_B0 <- sum(simpop[[2]][,pyears,])/B0[x]
# 
# 
#  return(c(MSY = MSY, FMSY = F_MSYv, SSB = SSB_MSY, SSBMSY_SSB0=SSBMSY_SSB0,
#           BMSY_B0=BMSY_B0, B = B, VB=V_BMSY+MSY))
# 
# }
# 
# 
# 
# 
#' Optimize yield for a single simulation
#' 
#' @param logFa log apical fishing mortality
#' @param Asize_c A vector of length areas with relative size of areas
#' @param nareas Number of area
#' @param maxage Maximum age
#' @param Ncurr Current N-at-age
#' @param pyears Number of projection years
#' @param M_age M-at-age
#' @param MatAge Maturity-at-age
#' @param WtAge Weight-at-age
#' @param Vuln Vulnerablity-at-age
#' @param Retc Retention-at-age
#' @param Prec Recruitment error
#' @param movc Movement matrix
#' @param SRrelc SR Relationship
#' @param Effind Historical effort
#' @param Spat_targc Spatial targeting
#' @param hc Steepness
#' @param R0c Unfished recruitment by area
#' @param SSBpRc Unfished spawning stock per recruit by area
#' @param aRc Ricker aR
#' @param bRc Ricker bR
#' @param Qc Catchability 
#' @param MPA A matrix of spatial closures by year
#' @param maxF A numeric value specifying the maximum fishing mortality for any single age class
#' @param useCPP logical - use the CPP code? For testing purposes only
#' @param SSB0c SSB0
#' @keywords internal
#'
#' @author A. Hordyk
#' 
optMSY <- function(logFa, Asize_c, nareas, maxage, Ncurr, pyears, M_age,
                 MatAge, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc,
                 R0c, SSBpRc, aRc, bRc, Qc, MPA, maxF, useCPP=TRUE, SSB0c) {

  FMSYc <- exp(logFa)
  if(!useCPP) {
    # simpop <- popdyn(nareas, maxage, Ncurr, pyears, M_age, Asize_c,
    #                  MatAge, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc,
    #                  R0c, SSBpRc, aRc, bRc, Qc, Fapic=FMSYc, MPA=MPA, maxF=maxF, control=2)
    # doesn't work with age-based movement
    
  } else {
    simpop <- popdynCPP(nareas, maxage, Ncurr, pyears, M_age, Asize_c,
                     MatAge, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc,
                     R0c, SSBpRc, aRc, bRc, Qc=0, Fapic=FMSYc, MPA=MPA, maxF=maxF, control=2,
                     SSB0c=SSB0c)
  }

  # Yield
  # Cn <- simpop[[7]]/simpop[[8]] * simpop[[1]] * (1-exp(-simpop[[8]])) # retained catch
  Cn <- simpop[[6]]/simpop[[8]] * simpop[[1]] * (1-exp(-simpop[[8]])) # removals
  # Cb <- Cn[,pyears,] * WtAge[,pyears]
  # -sum(Cb)
  Cb <- Cn[,(pyears-4):pyears,] * array(WtAge[,(pyears-4):pyears], dim=dim(Cn[,(pyears-4):pyears,]))
 
  -mean(apply(Cb,2,sum))
  

}



#' Calculate Reference Yield 
#'
#' @param x Integer, the simulation number
#' @param Asize A matrix (dimensions nsim by nareas) with relative size of areas
#' @param nareas The number of spatial areas
#' @param maxage The maximum age
#' @param N Array of the numbers-at-age in population. Dimensions are nsim, maxage, nyears, nareas. 
#' Only values from the first year are used, which is the current N-at-age.
#' @param pyears The number of years to project forward. Equal to 'nyears' for optimizing for q.
#' @param M_ageArray An array (dimensions nsim, maxage, nyears+proyears) with the natural mortality-at-age and year 
#' @param Mat_age An array (dimensions nsim, maxage, nyears+proyears) with the proportion mature for each age-class
#' @param Wt_age An array (dimensions nsim, maxage, nyears+proyears) with the weight-at-age and year 
#' @param V An array (dimensions nsim, maxage, nyears+proyears) with the vulnerability-at-age and year
#' @param retA An array (dimensions nsim, maxage, nyears+proyears) with the probability retained-at-age and year
#' @param Perr A matrix (dimensions nsim, nyears+proyears) with the recruitment deviations
#' @param mov An array (dimensions nsim, nareas, nareas) with the movement matrix
#' @param SRrel A numeric vector nsim long specifying the recruitment curve to use
#' @param Find A matrix (dimensions nsim, nyears) with the historical fishing effort 
#' @param Spat_targ A numeric vector nsim long with the spatial targeting
#' @param hs A numeric vector nsim long with the steepness values for each simulation
#' @param R0a A matrix (dimensions nsim, nareas) with the unfished recruitment by area
#' @param SSBpR A matrix (dimensions nsim, nareas) with the unfished spawning-per-recruit by area
#' @param aR A numeric vector nareas long with the Ricker SRR a values
#' @param bR A numeric vector nareas long with the Ricker SRR b values
#' @param MPA A matrix of spatial closures by year
#' @param maxF A numeric value specifying the maximum fishing mortality for any single age class
#' @param useCPP logical - use the CPP code? For testing purposes only
#' @param SSB0 SSB0
#' @author A. Hordyk
#' @export
#' @keywords internal
getFref3 <- function(x, Asize, nareas, maxage, N, pyears, M_ageArray, Mat_age, Wt_age,
                     V, retA, Perr, mov, SRrel, Find, Spat_targ, hs, R0a, SSBpR, aR, bR, 
                     MPA, maxF, useCPP=TRUE, SSB0) {
  
  opt <- optimize(optMSY, log(c(0.001, 10)), Asize_c=Asize[x,], nareas, maxage, Ncurr=N[x,,1,], 
                  pyears, M_age=M_ageArray[x,,], MatAge=Mat_age[x,,], 
                  WtAge=Wt_age[x,,], Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,], 
                  movc=mov[x,,,], SRrelc=SRrel[x], 
                  Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
                  SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], MPA=MPA, maxF=maxF, useCPP=useCPP,
                  SSB0c=SSB0[x])
  
  -opt$objective
  
}
 

# Input Control Functions Wrapper function for input control methods


#' Runs input control MPs on a Data object.
#' 
#' Function runs a MP (or MPs) of class 'Input' and returns a list: input
#' control recommendation(s) in element 1 and Data object in element 2.
#' 
#' 
#' @usage runInMP(Data, MPs = NA, reps = 100)
#' @param Data A object of class Data
#' @param MPs A vector of MPs of class 'Input'
#' @param reps Number of stochastic repititions - often not used in input
#' control MPs.
#' @author A. Hordyk
#' @export 
runInMP <- function(Data, MPs = NA, reps = 100) {
  
  nsims <- length(Data@Mort)
  if (.hasSlot(Data, "nareas")) {
    nareas <- Data@nareas   
  } else {
    nareas <- 2 
  }
  
  nMPs <- length(MPs)
  
  returnList <- list() # a list nMPs long containing MPs recommendations
  recList <- list() # a list containing nsim recommendations from a single MP 
  
  if (!sfIsRunning() | (nMPs < 8 & nsims < 8)) {
    for (ff in 1:nMPs) {
      temp <- sapply(1:nsims, MPs[ff], Data = Data, reps = reps)
      slots <- slotNames(temp[[1]])
      for (X in slots) { # sequence along recommendation slots 
        if (X == "Misc") { # convert to a list nsim by nareas
          rec <- lapply(temp, slot, name=X)
        } else {
          rec <- unlist(lapply(temp, slot, name=X))
        }
        if (X == "Spatial") { # convert to a matrix nsim by nareas
          rec <- matrix(rec, nsims, nareas, byrow=TRUE)  
        }
        
        recList[[X]] <- rec
        for (x in 1:nsims) Data@Misc[[x]] <- recList$Misc[[x]]
        recList$Misc <- NULL
      }
      returnList[[ff]] <- recList
    }
  } else {
    sfExport(list = c("Data"))
    for (ff in 1:nMPs) {
      temp <- sfSapply(1:nsims, MPs[ff], Data = Data, reps = reps)
      slots <- slotNames(temp[[1]])
      for (X in slots) { # sequence along recommendation slots 
        if (X == "Misc") { # convert to a list nsim by nareas
          rec <- lapply(temp, slot, name=X)
        } else {
          rec <- unlist(lapply(temp, slot, name=X))
        }
        if (X == "Spatial") { # convert to a matrix nsim by nareas
          rec <- matrix(rec, nsims, nareas, byrow=TRUE)  
        }
        
        recList[[X]] <- rec
        for (x in 1:nsims) Data@Misc[[x]] <- recList$Misc[[x]]
        recList$Misc <- NULL
      }
      returnList[[ff]] <- recList
    }   
  }
  
  return(list(returnList, Data))
}
  

projectEq <- function(x, Asize, nareas, maxage, N, pyears, M_ageArray, Mat_age, Wt_age,
                      V, retA, Perr, mov, SRrel, Find, Spat_targ, hs, R0a, SSBpR, aR, bR,
                      SSB0, B0, MPA, maxF, Nyrs, R0) {
  
  simpop <- popdynCPP(nareas, maxage, Ncurr=N[x,,1,],
                      pyears, M_age=M_ageArray[x,,], Asize_c=Asize[x,],
                      MatAge=Mat_age[x,,],
                      WtAge=Wt_age[x,,], Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,],
                      movc=mov[x,,,], SRrelc=SRrel[x],
                      Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,],
                      SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], Qc=0, Fapic=0, MPA=MPA,
                      maxF=maxF, control=3, SSB0c=SSB0[x])
  
  simpop[[1]][,Nyrs,]
  
}


# calcMSYRicker <- function(MSYyr, M_ageArray, Wt_age, retA, V, Perr_y, maxage,
#                           nareas, Mat_age, nsim, Asize, N, Spat_targ, hs,
#                           SRrel, mov, Find, R0a, SSBpR, aR, bR, SSB0, 
#                           B0, maxF=maxF, cur.yr) {
#   # Note: MSY and refY are calculated from total removals not total catch (different when Fdisc>0 and there is discarding)
#   # Make arrays for future conditions assuming current conditions
#   M_ageArrayp <- array(M_ageArray[,,cur.yr], dim=c(dim(M_ageArray)[1:2], MSYyr))
#   Wt_agep <- array(Wt_age[,,cur.yr], dim=c(dim(Wt_age)[1:2], MSYyr))
#   retAp <- array(retA[,,cur.yr], dim=c(dim(retA)[1:2], MSYyr))
#   Vp <- array(V[,,cur.yr], dim=c(dim(V)[1:2], MSYyr))
#   Perrp <- array(1, dim=c(dim(Perr_y)[1], MSYyr+maxage))
#   noMPA <- matrix(1, nrow=MSYyr, ncol=nareas)
#   Mat_agep <-abind::abind(rep(list(Mat_age[,,cur.yr]), MSYyr), along=3)
#   # optimize for MSY reference points
#   if (snowfall::sfIsRunning()) {
#     MSYrefs <- snowfall::sfSapply(1:nsim, getFMSY3, Asize, nareas=nareas, 
#                                   maxage=maxage, N=N, pyears=MSYyr,
#                                   M_ageArray=M_ageArrayp, Mat_age=Mat_agep, 
#                                   Wt_age=Wt_agep, V=Vp, retA=retAp,
#                                   Perr=Perrp, mov=mov, SRrel=SRrel, 
#                                   Find=Find, Spat_targ=Spat_targ, hs=hs,
#                                   R0a=R0a, SSBpR=SSBpR, aR=aR, bR=bR, SSB0=SSB0, 
#                                   B0=B0, MPA=noMPA, maxF=maxF)  
#   } else {
#     MSYrefs <- sapply(1:nsim, getFMSY3, Asize, nareas=nareas, maxage=maxage,
#                       N=N, pyears=MSYyr, M_ageArray=M_ageArrayp, Mat_age=Mat_agep, 
#                       Wt_age=Wt_agep, V=Vp, retA=retAp,Perr=Perrp, mov=mov, 
#                       SRrel=SRrel, Find=Find, Spat_targ=Spat_targ, hs=hs,
#                       R0a=R0a, SSBpR=SSBpR, aR=aR, bR=bR, SSB0=SSB0, B0=B0, 
#                       MPA=noMPA, maxF=maxF) 
#   }
#   MSYrefs
# }
  
# #' Apply output control recommendations and calculate population dynamics  
# #'
# #' @param y Projection year
# #' @param Asize relative size of areas (matrix nsim by nareas)
# #' @param TACused TAC recommendation
# #' @param TAC_f Implementation error on TAC
# #' @param lastCatch Catch from last year
# #' @param availB Total available biomass
# #' @param maxF Maximum fishing mortality
# #' @param Biomass_P Numeric array (nsim, maxage, proyears, nareas) with Biomass at age
# #' @param VBiomass_P Numeric array (nsim, maxage, proyears, nareas) with Vulnerable Biomass at age
# #' @param CB_P Numeric array (nsim, maxage, proyears, nareas) with Catch Biomass at age
# #' @param CB_Pret Numeric array (nsim, maxage, proyears, nareas) with Retained catch biomass at age
# #' @param FM_P Numeric array (nsim, maxage, proyears, nareas) with fishing mortality at age
# #' @param Z_P Numeric array (nsim, maxage, proyears, nareas) with total mortality at age
# #' @param Spat_targ Spatial targetting
# #' @param V_P Numeric array(nsim, maxage, nyears+proyears) with vulnerability at age
# #' @param retA_P Numeric array(nsim, maxage, nyears+proyears) with retention at age
# #' @param M_ageArray Numeric array (nsim, maxage, nyears+proyears) Natural mortality at age
# #' @param qs Catchability coefficient
# #' @param nyears Number of historical years
# #' @param nsim Number of simulations 
# #' @param maxage Maximum age
# #' @param nareas Number of areas
# #'
# # #' @export
# #'
# #' @author A. Hordyk
# #' 
# CalcOutput <- function(y, Asize, TACused, TAC_f, lastCatch, availB, maxF, Biomass_P, VBiomass_P, CB_P, CB_Pret,
#                        FM_P, Z_P, Spat_targ, V_P, retA_P, M_ageArray, qs, nyears, nsim, maxage, nareas) {
#   SAYRL <- as.matrix(expand.grid(1:nsim, 1:maxage, nyears, 1:nareas))  # Final historical year
#   SAYRt <- as.matrix(expand.grid(1:nsim, 1:maxage, y + nyears, 1:nareas))  # Trajectory year
#   SAYR <- as.matrix(expand.grid(1:nsim, 1:maxage, y, 1:nareas))
#   SYt <- SAYRt[, c(1, 3)]
#   SAYt <- SAYRt[, 1:3]
#   SR <- SAYR[, c(1, 4)]
#   SA1 <- SAYR[, 1:2]
#   S1 <- SAYR[, 1]
#   SY1 <- SAYR[, c(1, 3)]
#   SAY1 <- SAYR[, 1:3]
#   SYA <- as.matrix(expand.grid(1:nsim, 1, 1:maxage))  # Projection year
#   SY <- SYA[, 1:2]
#   SA <- SYA[, c(1, 3)]
#   SAY <- SYA[, c(1, 3, 2)]
#   S <- SYA[, 1]
#   
#   TACused[is.na(TACused)] <- lastCatch[is.na(TACused)] # if MP returns NA - TAC is set to catch from last year
#   
#   TACrec <- TACused             # TAC recommendation
#   TACusedE<- TAC_f[,y]*TACused   # TAC taken after implementation error
#   
#   maxC <- (1 - exp(-maxF)) * availB # maximum catch given maxF
#   TACusedE[TACusedE > maxC] <- maxC[TACusedE > maxC] # apply maxF limit - catch can't be higher than maxF * vulnerable biomass
#   
#   # fishdist <- (apply(VBiomass_P[, , y, ], c(1, 3), sum)^Spat_targ)/
#     # apply(apply(VBiomass_P[, , y, ], c(1, 3), sum)^Spat_targ, 1, mean)  # spatial preference according to spatial biomass
# 
#   fishdist <- (apply(VBiomass_P[, , y, ], c(1, 3), sum)^Spat_targ)/
#     apply(apply(VBiomass_P[, , y, ], c(1, 3), sum)^Spat_targ, 1, sum)  # spatial preference according to spatial biomass
#   
#   
#  
#   # If there is discard mortality, actual removals are higher than TACused
#   # calculate distribution of all effort
#   CB_P[SAYR] <- (Biomass_P[SAYR] * V_P[SAYt] * fishdist[SR])/Asize[SR] # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space
#   # calculate distribution of retained effort 
#   CB_Pret[SAYR] <- (Biomass_P[SAYR] * retA_P[SAYt] * fishdist[SR])/Asize[SR]  # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space
#   
#   retained <- apply(CB_Pret[,,y,], 1, sum)
#   actualremovals <- apply(CB_P[,,y,], 1, sum)
#   
#   ratio <- actualremovals/retained # ratio of actual removals to retained catch 
# 
#   temp <- CB_Pret[, , y, ]/apply(CB_Pret[, , y, ], 1, sum) # distribution of retained fish
#   CB_Pret[, , y, ] <- TACusedE * temp  # retained catch 
#   
#   temp <- CB_P[, , y, ]/apply(CB_P[, , y, ], 1, sum) # distribution of removals
#   CB_P[,,y,] <- TACusedE *  ratio * temp # scale up total removals 
# 
#   temp <- CB_P[SAYR]/(Biomass_P[SAYR] * exp(-M_ageArray[SAYt]/2))  # Pope's approximation
#   temp[temp > (1 - exp(-maxF))] <- 1 - exp(-maxF)
# 
#   FM_P[SAYR] <- -log(1 - temp)
#  
#   # calcFs <- lapply(1:nsim, getFs, y=y, Vuln=V_P, CB=CB_P, Bio=Biomass_P, Mage=M_ageArray, Fdist=fishdist,
#   #        maxage=maxage, nareas=nareas, nyears=nyears) # numerically calculate Fs
#   # 
#   # 
#   # FM_P[,,y,] <- aperm(array(unlist(calcFs, use.names=FALSE), dim=c(maxage, nareas, nsim)), c(3, 1, 2))
#   # FM_P[,,y,][FM_P[,,y,] > (1-exp(-maxF))]  <- 1 - exp(-maxF)
#   
#   Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt]
#   
#   Effort <- (-log(1 - apply(CB_P[, , y, ], 1, sum)/(apply(CB_P[, , y, ], 1, sum) + 
#                                                       apply(VBiomass_P[, , y, ], 1, sum))))/qs	 
#   out <- list()
#   out$Z_P <- Z_P 
#   out$FM_P <- FM_P 
#   out$CB_P <- CB_P
#   out$CB_Pret <- CB_Pret
#   out$TACused <- TACused 
#   out$TACrec <- TACrec
#   out$Effort <- Effort
#   out 
# }


# #' Internal function to calculate F-at-age given catch and biomass
# #'
# #' @param x Simulation
# #' @param y year
# #' @param Vuln Vulnerabilty
# #' @param CB Catch biomass
# #' @param Bio Biomass
# #' @param Mage M-at-age
# #' @param Fdist Fishing distribution
# #' @param maxage Maximum age
# #' @param nareas Number of areas
# #' @param nyears Number of historical years
# #' @keywords internal 
# #'
# #' @export
# #'
# #' @author A. Hordyk
# getFs <- function(x, y, Vuln, CB, Bio, Mage, Fdist, maxage, nareas, nyears) {
#   
#   doopt <- optimize(optF, interval=log(c(0.01, 10)), Vuln[x,,nyears+y], CB[x,,y,],
#                     Bio[x,,y,], Mage[x,,y+nyears], Fdist[x,], maxage,nareas)
#     
#   ind <- as.matrix(expand.grid(x, 1:maxage, 1:nareas))                
#   ind2 <- as.matrix(expand.grid(1, 1:maxage, 1:nareas))  
#   FM <- array(NA, dim=c(1, maxage, nareas))
#   FM[ind2] <- exp(doopt$minimum) * Vuln[ind] * Fdist[ind[,c(1,3)]]
#   FM
# }

# 
# #' Internal function to optimize for F
# #'
# #' @param fapic Apical fishing mortality
# #' @param vuln Vulnerability
# #' @param catch Catch
# #' @param bio Biomass
# #' @param mort Natural mortality
# #' @param fdist Fishing distribution
# #' @param maxage Maximum age
# #' @param nareas Number of areas
# #'
# #' @export
# #'
# #' @author A. Hordyk
# optF <- function(fapic, vuln, catch, bio, mort, fdist, maxage, nareas) {
#   FM <- array(NA, dim=c(maxage, nareas))
#   ind <- as.matrix(expand.grid(1:maxage, 1:nareas))
#   FM[ind] <- exp(fapic) * vuln[ind[,1]] * fdist[ind[,2]]
#   
#   # FM[ind] <- (exp(fapic) * vuln[ind[,1]] * fdist[ind[,2]]) / area_size[ind[,2]]
#   
#   Z <- FM + mort 
# 
#   pCatch <- FM/Z * bio* (1-exp(-Z))
#   (log(sum(pCatch)) - log(sum(catch)))^2
# 
# }



# #' Apply input control recommendations and calculate population dynamics  
# #'
# #' Internal function
# #' 
# #' @param y Simulation year
# #' @param Asize Matrix (nsim by nareas) with relative size of areas
# #' @param nyears Number of historical 
# #' @param proyears Number of projection years
# #' @param InputRecs Input control recommendations
# #' @param nsim Number of simulations
# #' @param nareas Number of areas
# #' @param LR5_P Length at 5 percent retention
# #' @param LFR_P Length at full retention
# #' @param Rmaxlen_P Retention of maximum length
# #' @param maxage Maximum age
# #' @param retA_P Retention at age
# #' @param retL_P Retention at length
# #' @param V_P Realized vulnerability at age
# #' @param V2 Gear vulnerability at age
# #' @param pSLarray Realized vulnerability at length
# #' @param SLarray2 Gear vulnerability at length
# #' @param DR Discard ratio
# #' @param maxlen maximum length
# #' @param Len_age Length-at-age
# #' @param CAL_binsmid Length-bin mid-points
# #' @param Fdisc Fraction of discarded fish that die
# #' @param nCALbins Number of length bins
# #' @param E_f Implementation error on effort recommendation
# #' @param SizeLim_f Implementation error on size limit
# #' @param VBiomass_P Vulnerable biomass-at-age
# #' @param Biomass_P Biomass-at-age
# #' @param Spat_targ Spatial targetting
# #' @param FinF Final fishing effort
# #' @param qvar Annual ariability in catchability
# #' @param qs Catchability
# #' @param qinc Numeric vector (nsim) increased
# #' @param CB_P Numeric array (nsim, maxage, proyears, nareas) Catch biomass at age
# #' @param CB_Pret Numeric array (nsim, maxage, proyears, nareas) Retained catch biomass at age
# #' @param FM_P Numeric array (nsim, maxage, proyears, nareas) Fishing mortality at age
# #' @param FM_retain Numeric array (nsim, maxage, proyears, nareas) Retained fishing mortality at age
# #' @param Z_P Numeric array (nsim, maxage, proyears, nareas) Total mortality at age
# #' @param M_ageArray Numeric array (nsim, maxage, nyears+proyears) Natural mortality at age
# #' @param LastEffort Numeric vector (nsim) with fishing effort from last year
# #' @param LastSpatial Numeric matrix (nsim, nareas) with spatial closures from last year
# #' @param LastAllocat Numeric vector (nsim) with allocation from last year
# #'
# #' @keywords internal
# #' @export
# #'
# #' @author A. Hordyk
# #' 
# CalcInput <- function(y, Linf, Asize, nyears, proyears, InputRecs, nsim, nareas, LR5_P, LFR_P,
#                       Rmaxlen_P, maxage, retA_P, retL_P, V_P, V2, pSLarray,
#                       SLarray2, DR, maxlen, Len_age, CAL_binsmid, Fdisc, 
#                       nCALbins, E_f, SizeLim_f, VBiomass_P, Biomass_P, Spat_targ,
#                       FinF, qvar, qs, qinc, CB_P, CB_Pret, FM_P, FM_retain, Z_P,
#                       M_ageArray, LastEffort, LastSpatial, LastAllocat) {
#   
#   SAYRL <- as.matrix(expand.grid(1:nsim, 1:maxage, nyears, 1:nareas))  # Final historical year
#   SAYRt <- as.matrix(expand.grid(1:nsim, 1:maxage, y + nyears, 1:nareas))  # Trajectory year
#   SAYR <- as.matrix(expand.grid(1:nsim, 1:maxage, y, 1:nareas))
#   SYt <- SAYRt[, c(1, 3)]
#   SAYt <- SAYRt[, 1:3]
#   SR <- SAYR[, c(1, 4)]
#   SA1 <- SAYR[, 1:2]
#   S1 <- SAYR[, 1]
#   SY1 <- SAYR[, c(1, 3)]
#   SAY1 <- SAYR[, 1:3]
#   SYA <- as.matrix(expand.grid(1:nsim, 1, 1:maxage))  # Projection year
#   SY <- SYA[, 1:2]
#   SA <- SYA[, c(1, 3)]
#   SAY <- SYA[, c(1, 3, 2)]
#   S <- SYA[, 1]
#   
#   # Change in Effort 
#   if (length(InputRecs$Effort) == 0) { # no effort recommendation
#     if (y==1) Ei <- LastEffort  * E_f[,y] # effort is unchanged but has implementation error
#     if (y>1) Ei <- LastEffort / E_f[,y-1]  * E_f[,y] # effort is unchanged but has implementation error
#   } else if (length(InputRecs$Effort) != nsim) {
#     stop("Effort recommmendation is not 'nsim' long.\n Does MP return Effort recommendation under all conditions?")
#   } else {
#     Ei <- InputRecs$Effort * E_f[,y] # effort adjustment with implementation error
#   }
#   
#   # Spatial 
#   if (all(is.na(InputRecs$Spatial))) { # no spatial recommendation 
#     Si <- LastSpatial # matrix(1, nsim, nareas) # spatial is unchanged - modify this if spatial closure in historical years  
#   } else if (any(is.na(InputRecs$Spatial))) {
#     stop("Spatial recommmendation has some NAs.\n Does MP return Spatial recommendation under all conditions?")
#   } else {
#     Si <-InputRecs$Spatial # change spatial fishing
#   }
#   
#   # Allocation 
#   if (length(InputRecs$Allocate) == 0) { # no allocation recommendation
#     Ai <- LastAllocat # rep(0, nsim) # allocation is unchanged 
#   } else if (length(InputRecs$Allocate) != nsim) {
#     stop("Allocate recommmendation is not 'nsim' long.\n Does MP return Allocate recommendation under all conditions?")
#   } else {
#     Ai <- InputRecs$Allocate # change in spatial allocation
#   }
#   # Retention Curve 
#   RetentFlag <- FALSE
#   # LR5 
#   if (length(InputRecs$LR5) == 0) { # no  recommendation
#     LR5_P[(y + nyears):(nyears+proyears),] <- matrix(LR5_P[y + nyears-1,], 
#                                                      nrow=(length((y + nyears):(nyears+proyears))),
#                                                      ncol=nsim, byrow=TRUE) # unchanged 
# 
#   } else if (length(InputRecs$LR5) != nsim) {
#     stop("LR5 recommmendation is not 'nsim' long.\n Does MP return LR5 recommendation under all conditions?")
#   } else {
#     LR5_P[(y + nyears):(nyears+proyears),] <- matrix(InputRecs$LR5 * SizeLim_f[,y], 
#                                                      nrow=(length((y + nyears):(nyears+proyears))),
#                                                      ncol=nsim, byrow=TRUE) # recommendation with implementation error
#     RetentFlag <- TRUE
#   }
#   # LFR 
#   if (length(InputRecs$LFR) == 0) { # no  recommendation
#     LFR_P[(y + nyears):(nyears+proyears),] <- matrix(LFR_P[y + nyears-1,], 
#                                                      nrow=(length((y + nyears):(nyears+proyears))),
#                                                      ncol=nsim, byrow=TRUE) # unchanged 
#   } else if (length(InputRecs$LFR) != nsim) {
#     stop("LFR recommmendation is not 'nsim' long.\n Does MP return LFR recommendation under all conditions?")
#   } else {
#     LFR_P[(y + nyears):(nyears+proyears),] <- matrix(InputRecs$LFR * SizeLim_f[,y], 
#                                                      nrow=(length((y + nyears):(nyears+proyears))),
#                                                      ncol=nsim, byrow=TRUE) # recommendation with implementation error
#     RetentFlag <- TRUE
#   }
#   # Rmaxlen 
#   if (length(InputRecs$Rmaxlen) == 0) { # no  recommendation
#     Rmaxlen_P[(y + nyears):(nyears+proyears),] <- matrix(Rmaxlen_P[y + nyears-1,], 
#                                                          nrow=(length((y + nyears):(nyears+proyears))),
#                                                          ncol=nsim, byrow=TRUE)   # unchanged 
#   
#   } else if (length(Rmaxlen) != nsim) {
#     stop("Rmaxlen recommmendation is not 'nsim' long.\n Does MP return Rmaxlen recommendation under all conditions?")
#   } else {
#     Rmaxlen_P[(y + nyears):(nyears+proyears),] <- matrix(InputRecs$Rmaxlen, 
#                                                          nrow=(length((y + nyears):(nyears+proyears))),
#                                                          ncol=nsim, byrow=TRUE) # recommendation
#     RetentFlag <- TRUE
#   }
#   # HS - harvest slot 
#  
#   if (length(InputRecs$HS) == 0) { # no  recommendation
#     HS <- rep(1E5, nsim) # no harvest slot 
#   } else if (length(InputRecs$HS) != nsim) {
#     stop("HS recommmendation is not 'nsim' long.\n Does MP return HS recommendation under all conditions?")
#   } else {
#     HS <- InputRecs$HS  * SizeLim_f[,y] # recommendation
#     RetentFlag <- TRUE
#   }
#   # Change in retention - update vulnerability and retention curves 
#   if (RetentFlag) {
#     yr <- y+nyears 
#     allyrs <- (y+nyears):(nyears+proyears)  # update vulnerabilty for all future years
#   
#     srs <- (Linf - LFR_P[yr,]) / ((-log(Rmaxlen_P[yr,],2))^0.5) # selectivity parameters are constant for all years
#     sls <- (LFR_P[yr,] - LR5_P[yr,]) / ((-log(0.05,2))^0.5)
#     
#     CAL_binsmidMat <- matrix(CAL_binsmid, nrow=nsim, ncol=length(CAL_binsmid), byrow=TRUE)
#     relLen <- t(sapply(1:nsim, getsel, lens=CAL_binsmidMat, lfs=LFR_P[yr,], sls=sls, srs=srs))
# 
#     for (yy in allyrs) {
#       # calculate new retention at age curve 
#       retA_P[ , , yy] <- t(sapply(1:nsim, getsel, lens=Len_age[,,yy], lfs=LFR_P[yy,], sls=sls, srs=srs))
#       
#       # calculate new retention at length curve 
#       retL_P[,, yy] <- relLen  
#     }
#    
#     # upper harvest slot 
#     aboveHS <- Len_age[,,allyrs]>HS
#     tretA_P <- retA_P[,,allyrs]
#     tretA_P[aboveHS] <- 0
#     retA_P[,,allyrs] <- tretA_P
#     for (ss in 1:nsim) {
#       index <- which(CAL_binsmid >= HS[ss])
#       retL_P[ss, index, allyrs] <- 0
#     }	
#     
#     dr <- aperm(abind::abind(rep(list(DR), maxage), along=3), c(2,3,1))
#     retA_P[,,allyrs] <- (1-dr[,,yr]) * retA_P[,,yr]
#     dr <- aperm(abind::abind(rep(list(DR), nCALbins), along=3), c(2,3,1))
#     retL_P[,,allyrs] <- (1-dr[,,yr]) * retL_P[,,yr]
#     
#     # update realized vulnerablity curve with retention and dead discarded fish 
#     Fdisc_array1 <- array(Fdisc, dim=c(nsim, maxage, length(allyrs)))
#     
#     V_P[,,allyrs] <- V2[,,allyrs] * (retA_P[,,allyrs] + (1-retA_P[,,allyrs])*Fdisc_array1)
#     
#     Fdisc_array2 <- array(Fdisc, dim=c(nsim, nCALbins, length(allyrs)))
#     pSLarray[,,allyrs]  <- SLarray2[,,allyrs] * (retL_P[,,allyrs]+ (1-retL_P[,,allyrs])*Fdisc_array2)
#     
#     # Realised Retention curves
#     retA_P[,,allyrs] <- retA_P[,,allyrs] * V_P[,,allyrs]
#     retL_P[,,allyrs] <- retL_P[,,allyrs] * pSLarray[,,allyrs] 
#      
#   }
#   
#   
#   newVB <- apply(Biomass_P[, , y, ] * V_P[SAYt], c(1, 3), sum)  # calculate total vuln biomass by area 
#   # fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ, 1, mean)  # spatial preference according to spatial vulnerable biomass
#   fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ, 1, sum)  # spatial preference according to spatial vulnerable biomass
#   Emult <- 1 + ((2/apply(fishdist * Si, 1, sum)) - 1) *   Ai  # allocate effort to new area according to fraction allocation Ai
#  
#   # fishing mortality with input control recommendation 
#   FM_P[SAYR] <- (FinF[S1] * Ei[S1] * V_P[SAYt] * Si[SR] * fishdist[SR] * Emult[S1] * qvar[SY1] * (qs[S1]*(1 + qinc[S1]/100)^y))/Asize[SR]
#   
#   # retained fishing mortality with input control recommendation
#   FM_retain[SAYR] <- (FinF[S1] * Ei[S1] * retA_P[SAYt] * Si[SR] * fishdist[SR] * Emult[S1] * qvar[SY1] * qs[S1]*(1 + qinc[S1]/100)^y)/Asize[SR]
#   
#   VBiomass_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt]  # update vulnerable biomass 
#   Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt] # calculate total mortality 
#   
#   CB_P[SAYR] <- FM_P[SAYR]/Z_P[SAYR] * Biomass_P[SAYR] * (1 - exp(-Z_P[SAYR]))
#   CB_Pret[SAYR] <- FM_retain[SAYR]/Z_P[SAYR] * Biomass_P[SAYR] * (1 - exp(-Z_P[SAYR]))
#   
#   out <- list() 
#   out$Z_P <- Z_P 
#   out$FM_P <- FM_P 
#   out$FM_retain <- FM_retain
#   out$CB_P <- CB_P
#   out$CB_Pret <- CB_Pret
#   out$Effort <- Ei 
#   out$retA_P <- retA_P
#   out$retL_P <- retL_P
#   out$V_P <- V_P 
#   out$pSLarray <- pSLarray
#   out$Si <- Si
#   out$Ai <- Ai
#   out
#   
# }
