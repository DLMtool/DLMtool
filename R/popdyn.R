



popdyn <- function(nareas, maxage, Ncurr, pyears, FMc, M_age, 
                   MatAge, WtAge, Vuln, Ret, Prec, movc, SRrelc, Effind, Spat_targc, hc, 
                   R0c, SSBpRc, aRc, bRc, Qc) {
  Narray <- array(NA, dim=c(maxage, pyears, nareas))
  SSNarray <- array(NA, dim=c(maxage, pyears, nareas))
  SBarray <- array(NA, dim=c(maxage, pyears, nareas))
  VBarray <- array(NA, dim=c(maxage, pyears, nareas))
  Marray <- array(NA, dim=c(maxage, pyears, nareas))
  FMarray <- array(NA, dim=c(maxage, pyears, nareas))
  Zarray <- array(NA, dim=c(maxage, pyears, nareas))
  
  Narray[,1,] <- Ncurr
  SSNarray[,1,] <- Ncurr * MatAge # spawning stock numbers
  SBarray[,1,] <- Narray[,1,] * WtAge[,1] * MatAge # spawning biomass
  VBarray[,1,] <- Narray[,1,] * WtAge[,1] * Vuln[,1] # vulnerable biomass
  Marray[,1,] <- M_age[,1] # M-at-age
  FMarray[,1,] <- FMc[,1,] # Fishing mortality at age and area
  SAYR <- as.matrix(expand.grid(1:maxage, 1, 1:nareas))  # Set up some array indexes age (A) year (Y) region/area (R)
  if (all(is.na(FMarray[,1,]))) {
    FMarray[SAYR] <- Effind[SAYR[,2]] * Qc * Vuln[SAYR[,1:2]] * fishdist[SAYR[,3]]
  }
  Zarray[,1,] <- Marray[,1,] + FMarray[,1,]

  for (y in 1:(pyears-1)) {
    
    NextYrN <- popdynOneTS(nareas, maxage, SSBcurr=colSums(SBarray[,y,]), Ncurr= Narray[,y,], 
                           Zcurr=Zarray[,y,], PerrYr=Prec[y], hc, R0c, SSBpRc, aRc, bRc, 
                           movc, SRrelc)
    
    Narray[,y+1,] <- NextYrN
    SSNarray[,y+1,] <- Narray[,y+1,] * MatAge # spawning stock numbers
    SBarray[,y+1,] <- Narray[,y+1,] * WtAge[,y+1] * MatAge # spawning biomass
    VBarray[,y+1,] <- Narray[,y+1,] * WtAge[,y+1] * Vuln[,y+1] # vulnerable biomass
    Marray[, y+1, ] <- M_age[,y]
    
    # Distribution of fishing effort 
    VBa <- colSums(VBarray[,y+1,]) # total vuln biomass in each area 
    fishdist <- VBa^Spat_targc/mean(VBa^Spat_targc)
    
    FMarray[,y+1,] <- FMc[,y+1,] # Fishing mortality at age and area
    
    SAYR <- as.matrix(expand.grid(1:maxage, y+1, 1:nareas))  # Set up some array indexes age (A) year (Y) region/area (R)
    if (all(is.na(FMarray[,y+1,]))) {
      FMarray[SAYR] <- Effind[SAYR[,2]] * Qc * Vuln[SAYR[,1:2]] * fishdist[SAYR[,3]]
    }
    Zarray[,y+1,] <- Marray[,y+1,] + FMarray[,y+1,]
    
  }
  
  out <- list()
  out$Narray <- Narray
  out$SSNarray <- SSNarray
  out$SBarray <- SBarray
  out$VBarray <- VBarray
  out$FMarray <- FMarray
  out$Zarray <- Zarray
  out
}


popdynOneTS <- function(nareas, maxage, SSBcurr, Ncurr, Zcurr, 
                   PerrYr, hc, R0c, SSBpRc, aRc, bRc, movc, SRrel)  {
  
  # set up some indices for indexed calculation

  indMov <- as.matrix(expand.grid(1:nareas, 1:nareas, 1:maxage)[3:1])  # Movement master index
  indMov2 <- indMov[, c(1, 2)]  # Movement from index
  indMov3 <- indMov[, c(2, 3)]  # Movement to index
  
  Nnext <- array(NA, dim=c(maxage, nareas))

  # Recruitment assuming regional R0 and stock wide steepness
  if (SRrel[1] == 1) {
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
                  V, Ret, Perr, mov, SRrel, Find, Spat_targ, hs, R0a, SSBpR, aR, bR, Qs, 
                  bounds = c(1e-05, 15)) {
  
  
  opt <- optimize(optQ, log(bounds), depc=dep[x], SSB0c=SSB0[x], nareas, maxage, Ncurr=N[x,,1,], pyears, FMc, 
                  M_age=M_ageArray[x,,], MatAge=Mat_age[x,], WtAge=Wt_age[x,,],
                  Vuln=V[x,,], Ret, Prec=Perr[x,], movc=mov[x,,], SRrelc=SRrel[x], 
                  Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
                  SSBpRc=SSBpR[x,], aRc=aR[x], bRc=bR[x])
  return(exp(opt$minimum))
}


optQ <- function(logQ, depc, SSB0c, nareas, maxage, Ncurr, pyears, FMc, M_age, 
                 MatAge, WtAge, Vuln, Ret, Prec, movc, SRrelc, Effind, Spat_targc, hc, 
                 R0c, SSBpRc, aRc, bRc) {
  simpop <- popdyn(nareas, maxage, Ncurr, pyears, FMc, M_age, 
                   MatAge, WtAge, Vuln, Ret, Prec, movc, SRrelc, Effind, Spat_targc, hc, 
                   R0c, SSBpRc, aRc, bRc, Qc=exp(logQ))
  
  
  ssb <- sum(simpop$SBarray[,pyears,])
  (log(depc) - log(ssb/SSB0c))^2
  
}


simYears <- function(x, nareas, maxage, N, pyears, FM, M_ageArray, Mat_age, Wt_age,
                     V, Ret, Perr, mov, SRrel, Find, Spat_targ, hs, R0a, SSBpR, aR, bR, Qs) {
  popdyn(nareas, maxage, Ncurr=N[x,,1,], pyears, FMc, 
         M_age=M_ageArray[x,,], MatAge=Mat_age[x,], WtAge=Wt_age[x,,],
         Vuln=V[x,,], Ret, Prec=Perr[x,], movc=mov[x,,], SRrelc=SRrel[x], 
         Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
         SSBpRc=SSBpR[x,], aRc=aR[x], bRc=bR[x], Qc=Qs[x])
  
}

 
  
  

