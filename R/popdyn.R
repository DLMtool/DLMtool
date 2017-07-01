

popdyn <- function(nareas, maxage, Ncurr, pyears, FMc, M_age, 
                   MatAge, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc, 
                   R0c, SSBpRc, aRc, bRc, Qc, Fapic=NULL, control=0) {
  Narray <- array(NA, dim=c(maxage, pyears, nareas))
  Barray <- array(NA, dim=c(maxage, pyears, nareas))
  SSNarray <- array(NA, dim=c(maxage, pyears, nareas))
  SBarray <- array(NA, dim=c(maxage, pyears, nareas))
  VBarray <- array(NA, dim=c(maxage, pyears, nareas))
  Marray <- array(NA, dim=c(maxage, pyears, nareas))
  FMarray <- array(NA, dim=c(maxage, pyears, nareas))
  FMretarray <- array(NA, dim=c(maxage, pyears, nareas))
  Zarray <- array(NA, dim=c(maxage, pyears, nareas))
  
  Narray[,1,] <- Ncurr
  Barray[,1,] <- Narray[,1,] * WtAge[,1]
  SSNarray[,1,] <- Ncurr * MatAge # spawning stock numbers
  SBarray[,1,] <- Narray[,1,] * WtAge[,1] * MatAge # spawning biomass
  VBarray[,1,] <- Narray[,1,] * WtAge[,1] * Vuln[,1] # vulnerable biomass
  Marray[,1,] <- M_age[,1] # M-at-age
  
  SAYR <- as.matrix(expand.grid(1:maxage, 1, 1:nareas))  # Set up some array indexes age (A) year (Y) region/area (R)
  
  # Distribution of fishing effort 
  VBa <- colSums(VBarray[,1,]) # total vuln biomass in each area 
  fishdist <- VBa^Spat_targc/mean(VBa^Spat_targc)
  if (control == 0) FMarray <- FMc # Fishing mortality at age and area
  if (control == 1) {
    FMarray[SAYR] <- Effind[SAYR[,2]] * Qc * Vuln[SAYR[,1:2]] * fishdist[SAYR[,3]]
    FMretarray[SAYR] <- Effind[SAYR[,2]] * Qc * Retc[SAYR[,1:2]] * fishdist[SAYR[,3]]
  }
  if (control == 2) {
    FMarray[SAYR] <- Fapic * Vuln[SAYR[,1:2]] * fishdist[SAYR[,3]]
    FMretarray[SAYR] <- Fapic * Retc[SAYR[,1:2]] * fishdist[SAYR[,3]]
  }
  Zarray[,1,] <- Marray[,1,] + FMarray[,1,]

  for (y in 1:(pyears-1)) {
    
    NextYrN <- popdynOneTS(nareas, maxage, SSBcurr=colSums(SBarray[,y,]), Ncurr=Narray[,y,], 
                           Zcurr=Zarray[,y,], PerrYr=Prec[y+maxage-1], hc, R0c, SSBpRc, aRc, bRc, 
                           movc, SRrelc)
    
    Narray[,y+1,] <- NextYrN
    Barray[,y+1,] <- Narray[,y+1,] * WtAge[,y+1]
    SSNarray[,y+1,] <- Narray[,y+1,] * MatAge # spawning stock numbers
    SBarray[,y+1,] <- Narray[,y+1,] * WtAge[,y+1] * MatAge # spawning biomass
    VBarray[,y+1,] <- Narray[,y+1,] * WtAge[,y+1] * Vuln[,y+1] # vulnerable biomass
    Marray[, y+1, ] <- M_age[,y]
    
    # Distribution of fishing effort 
    VBa <- colSums(VBarray[,y+1,]) # total vuln biomass in each area 
    fishdist <- VBa^Spat_targc/mean(VBa^Spat_targc)
    
    SAYR <- as.matrix(expand.grid(1:maxage, y+1, 1:nareas))  # Set up some array indexes age (A) year (Y) region/area (R)
    if (control ==1) {
      FMarray[SAYR] <- Effind[SAYR[,2]] * Qc * Vuln[SAYR[,1:2]] * fishdist[SAYR[,3]]
      FMretarray[SAYR] <- Effind[SAYR[,2]] * Qc * Retc[SAYR[,1:2]] * fishdist[SAYR[,3]]
    }
    if (control ==2) {
      FMarray[SAYR] <- Fapic * Vuln[SAYR[,1:2]] * fishdist[SAYR[,3]]
      FMretarray[SAYR] <- Fapic * Retc[SAYR[,1:2]] * fishdist[SAYR[,3]]
    }   

    Zarray[,y+1,] <- Marray[,y+1,] + FMarray[,y+1,]
    
  }
  
  out <- list()
  out$Narray <- Narray
  out$Barray <- Barray
  out$SSNarray <- SSNarray
  out$SBarray <- SBarray
  out$VBarray <- VBarray
  out$FMarray <- FMarray
  out$FMretarray <- FMretarray
  out$Zarray <- Zarray
  out
}


popdynOneTS <- function(nareas, maxage, SSBcurr, Ncurr, Zcurr, 
                   PerrYr, hc, R0c, SSBpRc, aRc, bRc, movc, SRrelc)  {
  
  # set up some indices for indexed calculation

  indMov <- as.matrix(expand.grid(1:nareas, 1:nareas, 1:maxage)[3:1])  # Movement master index
  indMov2 <- indMov[, c(1, 2)]  # Movement from index
  indMov3 <- indMov[, c(2, 3)]  # Movement to index
  
  Nnext <- array(NA, dim=c(maxage, nareas))

  # Recruitment assuming regional R0 and stock wide steepness
  if (SRrelc[1] == 1) {
    Nnext[1,  ] <- PerrYr *  (4 * R0c * hc * SSBcurr)/(SSBpRc * R0c * (1-hc) + (5*hc-1)*SSBcurr)                                                                                          
  } else {
    # most transpaRcent form of the Ricker uses alpha and beta paRcams
    Nnext[1,  ] <- PerrYr * aRc * SSBcurr * exp(-bRc * SSBcurr)
  } 
  
  # Mortality 
  Nnext[2:maxage, ] <- Ncurr[1:(maxage - 1),  ] * exp(-Zcurr[1:(maxage - 1), ])  # Total mortality
  
  # Movement of stock 
  temp <- array(Nnext[indMov2] * movc[indMov3], dim = c(nareas, nareas, maxage))  # Move individuals
  Nnext <- apply(temp, c(3, 1), sum)
  
  # Numbers-at-age at beginning of next yeaRc
  return(Nnext)

}

  

getq3 <- function(x, dep, SSB0, nareas, maxage, N, pyears, FM, M_ageArray, Mat_age, Wt_age,
                  V, retA, Perr, mov, SRrel, Find, Spat_targ, hs, R0a, SSBpR, aR, bR, qs, 
                  bounds = c(1e-05, 15)) {
  
  
  opt <- optimize(optQ, log(bounds), depc=dep[x], SSB0c=SSB0[x], nareas, maxage, Ncurr=N[x,,1,], 
                  pyears, FMc=FM[x,,,], 
                  M_age=M_ageArray[x,,], MatAge=Mat_age[x,], WtAge=Wt_age[x,,],
                  Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,], movc=mov[x,,], SRrelc=SRrel[x], 
                  Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
                  SSBpRc=SSBpR[x,], aRc=aR[x], bRc=bR[x])
  return(exp(opt$minimum))
}


optQ <- function(logQ, depc, SSB0c, nareas, maxage, Ncurr, pyears, FMc, M_age, 
                 MatAge, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc, 
                 R0c, SSBpRc, aRc, bRc) {
  simpop <- popdyn(nareas, maxage, Ncurr, pyears, FMc, M_age, 
                   MatAge, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc, 
                   R0c, SSBpRc, aRc, bRc, Qc=exp(logQ), control=1)
  
  
  ssb <- sum(simpop$SBarray[,pyears,])
  
  (log(depc) - log(ssb/SSB0c))^2
  
}


simYears <- function(x, nareas, maxage, N, pyears, FM, M_ageArray, Mat_age, Wt_age,
                     V, retA, Perr, mov, SRrel, Find, Spat_targ, hs, R0a, SSBpR, aR, bR, qs) {
  popdyn(nareas, maxage, Ncurr=N[x,,1,], pyears, FMc=FM[x,,,], 
         M_age=M_ageArray[x,,], MatAge=Mat_age[x,], WtAge=Wt_age[x,,],
         Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,], movc=mov[x,,], SRrelc=SRrel[x], 
         Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
         SSBpRc=SSBpR[x,], aRc=aR[x], bRc=bR[x], Qc=qs[x], control=1)
  
}

getFMSY3 <- function(x, nareas, maxage, N, pyears, FM, M_ageArray, Mat_age, Wt_age,
                     V, retA, Perr, mov, SRrel, Find, Spat_targ, hs, R0a, SSBpR, aR, bR, SSB0, B0) {
  
  opt <- optimize(optMSY, log(c(0.001, 5)), nareas, maxage, Ncurr=N[x,,1,], 
                  pyears, FMc=FM[x,,,], M_age=M_ageArray[x,,], MatAge=Mat_age[x,], 
                  WtAge=Wt_age[x,,], Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,], 
                  movc=mov[x,,], SRrelc=SRrel[x], 
                  Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
                  SSBpRc=SSBpR[x,], aRc=aR[x], bRc=bR[x])
                  
  MSY <- -opt$objective 
  simpop <- popdyn(nareas, maxage, Ncurr=N[x,,1,], 
                   pyears, FMc=FM[x,,,], M_age=M_ageArray[x,,], MatAge=Mat_age[x,], 
                   WtAge=Wt_age[x,,], Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,], 
                   movc=mov[x,,], SRrelc=SRrel[x], 
                   Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
                   SSBpRc=SSBpR[x,], aRc=aR[x], bRc=bR[x], Fapic=exp(opt$minimum), control=2)
  
  SSB_MSY <- sum(simpop$SBarray[,pyears,])
  V_BMSY <- sum(simpop$VBarray[,pyears,])
  F_MSYv <- -log(1 - (MSY/(V_BMSY+MSY)))
  
  SSBMSY_SSB0 <- sum(simpop$SBarray[,pyears,])/SSB0[x]
  B <- sum(simpop$Barray[,pyears,])
  BMSY_B0 <- sum(simpop$Barray[,pyears,])/B0[x]
  
  return(c(MSY = MSY, FMSY = F_MSYv, SSB = SSB_MSY,
           SSBMSY_SSB0=SSBMSY_SSB0, BMSY_B0=BMSY_B0, 
           B = B, VB=V_BMSY))
  
}

optMSY <- function(logFa, nareas, maxage, Ncurr, pyears, FMc, M_age, 
                 MatAge, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc, 
                 R0c, SSBpRc, aRc, bRc, Qc) {
  
  FMSYc <- exp(logFa)

  simpop <- popdyn(nareas, maxage, Ncurr, pyears, FMc, M_age, 
                   MatAge, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc, 
                   R0c, SSBpRc, aRc, bRc, Qc, Fapic=FMSYc, control=2)
  
  
  # Yield 
  Cn <- simpop$FMretarray/simpop$Zarray * simpop$Narray * (1-exp(-simpop$Zarray))
  Cb <- Cn[,pyears,] * WtAge[,pyears]
  -sum(Cb)
}


getFref3 <- function(x, nareas, maxage, N, pyears, FM, M_ageArray, Mat_age, Wt_age,
                     V, retA, Perr, mov, SRrel, Find, Spat_targ, hs, R0a, SSBpR, aR, bR) {
  
  opt <- optimize(optMSY, log(c(0.001, 5)), nareas, maxage, Ncurr=N[x,,1,], 
                  pyears, FMc=FM[x,,,], M_age=M_ageArray[x,,], MatAge=Mat_age[x,], 
                  WtAge=Wt_age[x,,], Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,], 
                  movc=mov[x,,], SRrelc=SRrel[x], 
                  Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
                  SSBpRc=SSBpR[x,], aRc=aR[x], bRc=bR[x])
  
  -opt$objective
  
}
 
  
  
CalcOutput <- function(y, TAC, TAC_f, lastCatch, availB, maxF, Biomass_P, VBiomass_P, CB_P, CB_Pret,
                       FM_P, Z_P, Spat_targ, V_P, retA_P, M_ageArray, qs, nyears, nsim, maxage, nareas) {
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
  
  TAC[is.na(TAC)] <- lastCatch[is.na(TAC)] # if MP returns NA - TAC is set to catch from last year
  
  TACrec <- TAC             # TAC recommendation
  TACused<- TAC_f[,y]*TAC   # TAC taken after implementation error
  
  maxC <- (1 - exp(-maxF)) * availB # maximum catch given maxF
  TACused[TACused > maxC] <- maxC[TACused > maxC] # apply maxF limit - catch can't be higher than maxF * vulnerable biomass 
  
  fishdist <- (apply(VBiomass_P[, , y, ], c(1, 3), sum)^Spat_targ)/
    apply(apply(VBiomass_P[, , y, ], c(1, 3), sum)^Spat_targ, 1, mean)  # spatial preference according to spatial biomass
  
  # If there is discard mortality, actual removals are higher than TACused
  # calculate distribution of all effort
  CB_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt] * fishdist[SR] # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space
  # calculate distribution of retained effort 
  CB_Pret[SAYR] <- Biomass_P[SAYR] * retA_P[SAYt] * fishdist[SR]  # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space
  
  retained <- apply(CB_Pret[,,y,], 1, sum)
  actualremovals <- apply(CB_P[,,y,], 1, sum)
  

  ratio <- actualremovals/retained # ratio of actual removals to retained catch 
  
  temp <- CB_Pret[, , y, ]/apply(CB_Pret[, , y, ], 1, sum) # distribution of retained fish
  CB_Pret[, , y, ] <- TACused * temp  # retained catch 
  
  temp <- CB_P[, , y, ]/apply(CB_P[, , y, ], 1, sum) # distribution of removals
  CB_P[,,y,] <- TACused *  ratio * temp # scale up total removals 

  temp <- CB_P[SAYR]/(Biomass_P[SAYR] * exp(-M_ageArray[SAYt]/2))  # Pope's approximation
  temp[temp > (1 - exp(-maxF))] <- 1 - exp(-maxF)
  FM_P[SAYR] <- -log(1 - temp)
  Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt]
  
  Effort <- (-log(1 - apply(CB_P[, , y, ], 1, sum)/(apply(CB_P[, , y, ], 1, sum) + 
                                                      apply(VBiomass_P[, , y, ], 1, sum))))/qs	 
  out <- list()
  out$Z_P <- Z_P 
  out$FM_P <- FM_P 
  out$CB_P <- CB_P
  out$CB_Pret <- CB_Pret
  out$TACused <- TACused 
  out$TACrec <- TACrec
  out$Effort <- Effort
  out 
}


CalcInput <- function(y, nyears, proyears, InputRecs, nsim, nareas, LR5_P, LFR_P, Rmaxlen_P, maxage,
                      retA_P, retL_P, V_P, V2, pSLarray, SLarray2, DR, maxlen, Len_age, CAL_binsmid, Fdisc, nCALbins, E_f, SizeLim_f,
                      VBiomass_P, Biomass_P, Spat_targ, FinF, qvar, qs, qinc, CB_P, CB_Pret, FM_P, FM_retain, Z_P, M_ageArray) {
  
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
  
  # Change in Effort 
  if (length(InputRecs$Effort) == 0) { # no effort recommendation
    Ei <- rep(1, nsim) # effort is unchanged 
  } else if (length(InputRecs$Effort) != nsim) {
    stop("Effort recommmendation is not 'nsim' long.\n Does MP return Effort recommendation under all conditions?")
  } else {
    Ei <- InputRecs$Effort * E_f[,y] # effort adjustment with implementation error
  }
  
  # Spatial 
  if(all(is.na(InputRecs$Spatial))) { # no spatial recommendation 
    Si <- matrix(1, nsim, nareas) # spatial is unchanged - modify this if spatial closure in historical years  
  } else if (any(is.na(InputRecs$Spatial))) {
    stop("Spatial recommmendation has some NAs.\n Does MP return Spatial recommendation under all conditions?")
  } else {
    Si <-InputRecs$Spatial # change spatial fishing
  }
  
  # Allocation 
  if (length(InputRecs$Allocation) == 0) { # no allocation recommendation
    Ai <- rep(0, nsim) # allocation is unchanged 
  } else if (length(InputRecs$Allocation) != nsim) {
    stop("Allocation recommmendation is not 'nsim' long.\n Does MP return Allocation recommendation under all conditions?")
  } else {
    Ai <- InputRecs$Allocation # change in spatial allocation
  }
  # Retention Curve 
  RetentFlag <- FALSE
  # LR5 
  if (length(InputRecs$LR5) == 0) { # no  recommendation
    LR5_P[(y + nyears):(nyears+proyears),] <- LR5_P[y + nyears-1,] # unchanged 
  } else if (length(InputRecs$LR5) != nsim) {
    stop("LR5 recommmendation is not 'nsim' long.\n Does MP return LR5 recommendation under all conditions?")
  } else {
    LR5_P[(y + nyears):(nyears+proyears),] <- InputRecs$LR5 * SizeLim_f[,y] # recommendation with implementation error
    RetentFlag <- TRUE
  }
  # LFR 
  if (length(InputRecs$LFR) == 0) { # no  recommendation
    LFR_P[(y + nyears):(nyears+proyears),] <- LFR_P[y + nyears-1,] # unchanged 
  } else if (length(InputRecs$LFR) != nsim) {
    stop("LFR recommmendation is not 'nsim' long.\n Does MP return LFR recommendation under all conditions?")
  } else {
    LFR_P[(y + nyears):(nyears+proyears),] <- InputRecs$LFR * SizeLim_f[,y] # recommendation with implementation error
    RetentFlag <- TRUE
  }
  # Rmaxlen 
  if (length(InputRecs$Rmaxlen) == 0) { # no  recommendation
    Rmaxlen_P[(y + nyears):(nyears+proyears),] <- Rmaxlen_P[y + nyears-1,] # unchanged 
  } else if (length(Rmaxlen) != nsim) {
    stop("Rmaxlen recommmendation is not 'nsim' long.\n Does MP return Rmaxlen recommendation under all conditions?")
  } else {
    Rmaxlen_P[(y + nyears):(nyears+proyears),] <- InputRecs$Rmaxlen # recommendation
    RetentFlag <- TRUE
  }
  # HS - harvest slot 
  if (length(InputRecs$HS) == 0) { # no  recommendation
    HS <- rep(1E5, nsim) # no harvest slot 
  } else if (length(HS) != nsim) {
    stop("HS recommmendation is not 'nsim' long.\n Does MP return HS recommendation under all conditions?")
  } else {
    HS <- InputRecs$HS  * SizeLim_f[,y] # recommendation
    RetentFlag <- TRUE
  }
  # Change in retention - update vulnerability and retention curves 
  if (RetentFlag) {
    s1 <- sapply(1:nsim, function(i) optimize(getSlope1, interval = c(0, 1e+05), 
                                              LFS = LFR_P[y+nyears,i], 
                                              L0.05 = LR5_P[y + nyears,i])$minimum)
    s2 <- sapply(1:nsim, function(i) optimize(getSlope2, interval = c(0, 1e+05), 
                                              LFS = LFR_P[y+nyears,i], s1=s1[i], 
                                              maxlen=maxlen[i], 
                                              MaxSel= Rmaxlen_P[y + nyears,i])$minimum)
    # calculate new retention at age curve 
    retA_P[ , , y+nyears] <- t(sapply(1:nsim, function(i) 
      TwoSidedFun(LFR_P[y+nyears,i], s1[i], s2[i], lens=Len_age[i,,y+nyears])))
    # calculate new retention at length curve 
    retL_P[,, y+nyears] <- t(sapply(1:nsim, function(i) 
      TwoSidedFun(LFR_P[y+nyears,i], s1[i], s2[i], lens=CAL_binsmid))) 
    
    # upper harvest slot 
    retA_P[Len_age[, , (y + nyears)] >= HS] <- 0
    for (ss in 1:nsim) {
      index <- which(CAL_binsmid >= HS[ss])
      retL_P[ss, index, (y+nyears):(nyears+proyears)] <- 0 
    }	
    
    # correct retention curve - retention at age/length must <= selectivity (you can't retain fish you don't catch!)
    retA_P[,,y+nyears] <- (1-DR[y+nyears, ]) * matrix(mapply(pmin, retA_P[,,y+nyears], V_P[,,y+nyears]), nsim, maxage)
    retL_P[,,y+nyears] <- (1-DR[y+nyears, ]) * matrix(mapply(pmin, retL_P[,,y+nyears], pSLarray[,,y+nyears]), nsim, nCALbins)
    
    # update realized vulnerablity curve with retention and dead discarded fish 
    V_P[,,y+nyears] <- matrix(mapply(pmax, retA_P[,,y+nyears] + 
                                       (abs(retA_P[,,y+nyears]-V2[,,y+nyears])*Fdisc), retA_P[,,y+nyears]), nsim, maxage)
    pSLarray[,,y+nyears] <- matrix(mapply(pmax, retL_P[,,y+nyears] + 
                                            (abs(retL_P[,,y+nyears] - SLarray2[,,y+nyears])*Fdisc), retL_P[,,y+nyears]), nsim, nCALbins)
    
  }
  
  newVB <- apply(Biomass_P[, , y, ] * V_P[SAYt], c(1, 3), sum)  # calculate total vuln biomass by area 
  fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ, 1, mean)  # spatial preference according to spatial vulnerable biomass
  Emult <- 1 + ((2/apply(fishdist * Si, 1, sum)) - 1) *   Ai  # allocate effort to new area according to fraction allocation Ai
  # fishing mortality with input control recommendation 
  FM_P[SAYR] <- FinF[S1] * Ei[S1] * V_P[SAYt] * Si[SR] * fishdist[SR] * Emult[S1] * qvar[SY1] * qs[S1]^(1 + qinc[S1]/100)^y
  # retained fishing mortality with input control recommendation
  FM_retain[SAYR] <- FinF[S1] * Ei[S1] * retA_P[SAYt] * Si[SR] * fishdist[SR] * Emult[S1] * qvar[SY1] * qs[S1]^(1 + qinc[S1]/100)^y
  
  VBiomass_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt]  # update vulnerable biomass 
  Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt] # calculate total mortality 
  
  CB_P[SAYR] <- FM_P[SAYR]/Z_P[SAYR] * VBiomass_P[SAYR] * (1 - exp(-Z_P[SAYR]))
  CB_Pret[SAYR] <- FM_retain[SAYR]/Z_P[SAYR] * VBiomass_P[SAYR] * (1 - exp(-Z_P[SAYR]))
  
  out <- list() 
  out$Z_P <- Z_P 
  out$FM_P <- FM_P 
  out$CB_P <- CB_P
  out$CB_Pret <- CB_Pret
  out$Effort <- Ei 
  out$retA_P <- retA_P
  out$retL_P <- retL_P
  out$V_P <- V_P 
  out$pSLarray <- pSLarray
  out
  
}
