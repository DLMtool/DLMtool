#' optimize for catchability (q)
#' 
#' Function optimizes catchability (q, where F=qE) required to get to user-specified stock
#' depletion
#'
#' @param x Integer, the simulation number
#' @param dep A numeric vector nsim long of sampled depletion
#' @param SSB0 A numeric vector nsim long of total unfished spawning biomass
#' @param nareas The number of spatial areas
#' @param maxage The maximum age
#' @param N Array of the numbers-at-age in population. Dimensions are nsim, maxage, nyears, nareas. 
#' Only values from the first year (i.e N[,,1,]) are used, which is the current N-at-age.
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
#' @param useCPP logical - use the CPP code? For testing purposes only
#'
#' @author A. Hordyk
#' @export
#'
getq3 <- function(x, dep, SSB0, nareas, maxage, N, pyears, M_ageArray, Mat_age, Asize, Wt_age,
                  V, retA, Perr, mov, SRrel, Find, Spat_targ, hs, R0a, SSBpR, aR, bR, 
                  bounds = c(1e-05, 15), maxF, useCPP=TRUE) {
  
  opt <- optimize(optQ, log(bounds), depc=dep[x], SSB0c=SSB0[x], nareas, maxage, Ncurr=N[x,,1,], 
                  pyears, M_age=M_ageArray[x,,], MatAge=Mat_age[x,,], Asize_c=Asize[x,], WtAge=Wt_age[x,,],
                  Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,], movc=mov[x,,], SRrelc=SRrel[x], 
                  Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
                  SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], maxF, useCPP=useCPP)
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
#' @param useCPP Logical. Use the CPP code?
#'
#' @author A. Hordyk
#' @export

optQ <- function(logQ, depc, SSB0c, nareas, maxage, Ncurr, pyears, M_age, Asize_c,
                 MatAge, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc, 
                 R0c, SSBpRc, aRc, bRc, maxF, useCPP) {
  if (!useCPP) {
    simpop <- popdyn(nareas, maxage, Ncurr, pyears, M_age, Asize_c,
                     MatAge, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc, 
                     R0c=R0c, SSBpRc=SSBpRc, aRc=aRc, bRc=bRc, Qc=exp(logQ), maxF=maxF, control=1) 
    ssb <- sum(simpop$SBarray[,pyears,])
    
  } else {
    simpop <- popdynCPP(nareas, maxage, Ncurr, pyears, M_age, Asize_c,
                        MatAge, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc, 
                        R0c=R0c, SSBpRc=SSBpRc, aRc=aRc, bRc=bRc, Qc=exp(logQ), Fapic=0, 
                        maxF=maxF, control=1) 
    
    ssb <- sum(simpop[[4]][,pyears,])
  }
  
  (log(depc) - log(ssb/SSB0c))^2
  
}


#' Population dynamics model
#'
#' @param nareas Integer. The number of spatial areas
#' @param maxage Integer. The maximum age
#' @param Ncurr Numeric matrix (dimensions maxage, nareas) with the current N-at-age
#' @param pyears Integer. Number of years to project the model forward
#' @param M_age Numeric matrix (dimensions maxage, pyears) with natural mortality at age
#' @param Asize_c Numeric vector (length nareas) with size of each area
#' @param MatAge Numeric matrix (dimensions maxage, nyears+proyears) with proportion mature for each age-class
#' @param WtAge Numeric matrix (dimensions maxage, pyears) with weight-at-age 
#' @param Vuln Numeric matrix (dimensions maxage, pyears) with proportion vulnerable-at-age
#' @param Retc Numeric matrix (dimensions maxage, pyears) with proportion retained-at-age
#' @param Prec Numeric vector (length pyears) with recruitment error
#' @param movc Numeric matrix (dimensions nareas, nareas) with movement matrix
#' @param SRrelc Integer. Stock-recruitment curve
#' @param Effind Numeric vector (length pyears) with the fishing effort by year
#' @param Spat_targc Integer. Value of spatial targetting
#' @param hc Numeric. Steepness of stock-recruit relationship
#' @param R0c Numeric vector of length nareas with unfished recruitment by area
#' @param SSBpRc Numeric vector of length nareas with unfished spawning per recruit by area
#' @param aRc Numeric. Ricker SRR a value
#' @param bRc Numeric. Ricker SRR b value
#' @param Qc Numeric. Catchability coefficient
#' @param Fapic Numeric. Apical F value
#' @param maxF A numeric value specifying the maximum fishing mortality for any single age class
#' @param control Integer. 1 to use q and effort to calculate F, 2 to use Fapic (apical F) and 
#' vulnerablity to calculate F.
#' 
#' @author A. Hordyk
#'
#' @return A named list of length 8 containing with arrays (dimensions: maxage, pyears, nareas)
#' containing numbers-at-age, biomass-at-age, spawning stock numbers, spawning biomass, 
#' vulnerable biomass, fishing mortality, retained fishing mortality, and total mortality
#' @export
#'
popdyn <- function(nareas, maxage, Ncurr, pyears, M_age, Asize_c,
                   MatAge, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc, 
                   R0c, SSBpRc, aRc, bRc, Qc, Fapic=NULL, maxF, control=1) {
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
  SSNarray[,1,] <- Ncurr * MatAge[,1] # spawning stock numbers
  SBarray[,1,] <- Narray[,1,] * WtAge[,1] * MatAge[,1] # spawning biomass
  VBarray[,1,] <- Narray[,1,] * WtAge[,1] * Vuln[,1] # vulnerable biomass
  Marray[,1,] <- M_age[,1] # M-at-age
  
  SAYR <- as.matrix(expand.grid(1:maxage, 1, 1:nareas))  # Set up some array indexes age (A) year (Y) region/area (R)
  
  # Distribution of fishing effort 
  VBa <- colSums(VBarray[,1,]) # total vuln biomass in each area 
  
  # fishdist <- VBa^Spat_targc/mean(VBa^Spat_targc)
  fishdist <- VBa^Spat_targc/sum(VBa^Spat_targc)
  
  Asize_mat <- matrix(Asize_c, nrow=maxage, ncol=nareas, byrow=TRUE)
                    
  if (control == 1) {
    FMarray[SAYR] <- (Effind[SAYR[,2]] * Qc * Vuln[SAYR[,1:2]] * fishdist[SAYR[,3]])/Asize_mat
    FMretarray[SAYR] <- (Effind[SAYR[,2]] * Qc * Retc[SAYR[,1:2]] * fishdist[SAYR[,3]])/Asize_mat
  }
  if (control == 2) {
    FMarray[SAYR] <- (Fapic * Vuln[SAYR[,1:2]] * fishdist[SAYR[,3]])/Asize_mat
    FMretarray[SAYR] <- (Fapic * Retc[SAYR[,1:2]] * fishdist[SAYR[,3]])/Asize_mat
  }
  
  FMarray[,1,][FMarray[,1,] > (1 - exp(-maxF))] <- 1 - exp(-maxF)
  FMretarray[,1,][FMretarray[,1,] > (1 - exp(-maxF))] <- 1 - exp(-maxF)
 
  Zarray[,1,] <- Marray[,1,] + FMarray[,1,]

  for (y in 1:(pyears-1)) {
    
    NextYrN <- popdynOneTS(nareas, maxage, SSBcurr=colSums(SBarray[,y,]), Ncurr=Narray[,y,], 
                           Zcurr=Zarray[,y,], PerrYr=Prec[y+maxage+1], hc, R0c, SSBpRc, aRc, bRc, 
                           movc, SRrelc)
    
    Narray[,y+1,] <- NextYrN
    Barray[,y+1,] <- Narray[,y+1,] * WtAge[,y+1]
    SSNarray[,y+1,] <- Narray[,y+1,] * MatAge[,y+1] # spawning stock numbers
    SBarray[,y+1,] <- Narray[,y+1,] * WtAge[,y+1] * MatAge[,y+1] # spawning biomass
    VBarray[,y+1,] <- Narray[,y+1,] * WtAge[,y+1] * Vuln[,y+1] # vulnerable biomass
    Marray[, y+1, ] <- M_age[,y+1]
    
    # Distribution of fishing effort 
    VBa <- colSums(VBarray[,y+1,]) # total vuln biomass in each area 
    # fishdist <- VBa^Spat_targc/mean(VBa^Spat_targc)
    fishdist <- VBa^Spat_targc/sum(VBa^Spat_targc)
    
    SAYR <- as.matrix(expand.grid(1:maxage, y+1, 1:nareas))  # Set up some array indexes age (A) year (Y) region/area (R)
    if (control ==1) {
      FMarray[SAYR] <- (Effind[SAYR[,2]] * Qc * Vuln[SAYR[,1:2]] * fishdist[SAYR[,3]])/Asize_mat
      FMretarray[SAYR] <- (Effind[SAYR[,2]] * Qc * Retc[SAYR[,1:2]] * fishdist[SAYR[,3]])/Asize_mat
    }
    if (control ==2) {
      FMarray[SAYR] <- (Fapic * Vuln[SAYR[,1:2]] * fishdist[SAYR[,3]])/Asize_mat
      FMretarray[SAYR] <- (Fapic * Retc[SAYR[,1:2]] * fishdist[SAYR[,3]])/Asize_mat
    }   
    FMarray[SAYR][FMarray[SAYR] > (1 - exp(-maxF))] <- 1 - exp(-maxF)
    FMretarray[SAYR][FMretarray[SAYR] > (1 - exp(-maxF))] <- 1 - exp(-maxF)
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


#' Population dynamics model for one annual time-step
#'
#' Project population forward one time-step given current numbers-at-age and total mortality
#'
#' @param nareas The number of spatial areas
#' @param maxage The maximum age 
#' @param SSBcurr A numeric vector of length nareas with the current spawning biomass in each area
#' @param Ncurr A numeric matrix (maxage, nareas) with current numbers-at-age in each area
#' @param Zcurr A numeric matrix (maxage, nareas) with total mortality-at-age in each area
#' @param PerrYr A numeric value with recruitment deviation for current year 
#' @param hs Steepness of SRR
#' @param R0c Numeric vector with unfished recruitment by area
#' @param SSBpRc Numeric vector with unfished spawning stock per recruit by area 
#' @param aRc Numeric vector with Ricker SRR a parameter by area
#' @param bRc Numeric vector with Ricker SRR b parameter by area
#' @param movc Numeric matrix (nareas by nareas) with the movement matrix
#' @param SRrelc Integer indicating the stock-recruitment relationship to use (1 for Beverton-Holt, 2 for Ricker)
#' @author A. Hordyk
#' 
#' @export
#' @keywords internal
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
    # most transparent form of the Ricker uses alpha and beta params
    Nnext[1,  ] <- PerrYr * aRc * SSBcurr * exp(-bRc * SSBcurr)
  } 
  
  # Mortality 
  Nnext[2:maxage, ] <- Ncurr[1:(maxage - 1),  ] * exp(-Zcurr[1:(maxage - 1), ])  # Total mortality
  
  # Movement of stock 
  temp <- array(Nnext[indMov2] * movc[indMov3], dim = c(nareas, nareas, maxage))  # Move individuals
  Nnext <- apply(temp, c(3, 1), sum)
  
  # Numbers-at-age at beginning of next year
  return(Nnext)

}


#' Simulate population dynamics for historical years
#'
#' @param x Integer, the simulation number 
#' @param nareas The number of spatial areas
#' @param maxage The maximum age
#' @param N Array of the numbers-at-age in population. Dimensions are nsim, maxage, nyears, nareas. 
#' Only values from the first year (i.e N[,,1,]) are used, which is the current N-at-age.
#' @param pyears The number of years to project forward. Equal to 'nyears' for optimizing for q.
#' @param M_ageArray An array (dimensions nsim, maxage, nyears+proyears) with the natural mortality-at-age and year
#' @param Asize A matrix (dimensions nsim, nareas) of size of areas 
#' @param Mat_age A matrix (dimensions nsim, maxage) with the proportion mature for each age-class
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
#' @param aR A numeric vector nsim long with the Ricker SRR a values
#' @param bR A numeric vector nsim long with the Ricker SRR b values
#' @param qs A numeric vector nsim long with catchability coefficients
#' @param maxF A numeric value specifying the maximum fishing mortality for any single age class
#' @param useCPP logical - use the CPP code? For testing purposes only 
#' 
#' @author A. Hordyk
#' 
#' 
#' @export
simYears <- function(x, nareas, maxage, N, pyears, M_ageArray, Asize, Mat_age, Wt_age,
                     V, retA, Perr, mov, SRrel, Find, Spat_targ, hs, R0a, SSBpR, aR, bR, qs, 
                     maxF, useCPP=TRUE) {
  if(!useCPP) {
    popdyn(nareas, maxage, Ncurr=N[x,,1,], pyears,  
           M_age=M_ageArray[x,,], Asize_c=Asize[x,], MatAge=Mat_age[x,,], WtAge=Wt_age[x,,],
           Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,], movc=mov[x,,], SRrelc=SRrel[x], 
           Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
           SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], Qc=qs[x], maxF=maxF, control=1)
  } else {
    popdynCPP(nareas, maxage, Ncurr=N[x,,1,], pyears,  
           M_age=M_ageArray[x,,], Asize_c=Asize[x,], MatAge=Mat_age[x,,], WtAge=Wt_age[x,,],
           Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,], movc=mov[x,,], SRrelc=SRrel[x], 
           Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
           SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], Qc=qs[x], Fapic=0, maxF=maxF, control=1)
  }
  
}



#' Calculate FMSY and related metrics using Rcpp code
#'
#' @param x Integer, the simulation number
#' @param Asize A matrix (nsim by nareas) with size of areas
#' @param nareas The number of spatial areas
#' @param maxage The maximum age
#' @param N Array of the numbers-at-age in population. Dimensions are nsim, maxage, nyears, nareas. 
#' Only values from the first year (i.e N[,,1,]) are used, which is the current N-at-age.
#' @param pyears The number of years to project forward. Equal to 'nyears' for optimizing for q.
#' @param M_ageArray An array (dimensions nsim, maxage, nyears+proyears) with the natural mortality-at-age and year 
#' @param Mat_age A matrix (dimensions nsim, maxage) with the proportion mature for each age-class
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
#' @param aR A numeric vector nsim long with the Ricker SRR a values
#' @param bR A numeric vector nsim long with the Ricker SRR b values
#' @param SSB0 Unfished spawning biomass
#' @param B0 Unfished total biomass
#' @param maxF A numeric value specifying the maximum fishing mortality for any single age class
#' @param useCPP logical - use the CPP code? For testing purposes only
#'
#' @author A. Hordyk
#' @export
#'
getFMSY3 <- function(x, Asize, nareas, maxage, N, pyears, M_ageArray, Mat_age, Wt_age,
                     V, retA, Perr, mov, SRrel, Find, Spat_targ, hs, R0a, SSBpR, aR, bR, 
                     SSB0, B0, maxF, useCPP=TRUE) {
  
  opt <- optimize(optMSY, log(c(0.001, 10)), Asize_c=Asize[x,], nareas, maxage, Ncurr=N[x,,1,], 
                  pyears, M_age=M_ageArray[x,,], MatAge=Mat_age[x,,], 
                  WtAge=Wt_age[x,,], Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,], 
                  movc=mov[x,,], SRrelc=SRrel[x], 
                  Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
                  SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], maxF=maxF, useCPP=useCPP)
  
  MSY <- -opt$objective 
  
  if (!useCPP) {
    simpop <- popdyn(nareas, maxage, Ncurr=N[x,,1,], 
                     pyears, M_age=M_ageArray[x,,], Asize_c=Asize[x,],
                     MatAge=Mat_age[x,,], 
                     WtAge=Wt_age[x,,], Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,], 
                     movc=mov[x,,], SRrelc=SRrel[x], 
                     Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
                     SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], Fapic=exp(opt$minimum), maxF=maxF, control=2) 
    
    # calculate B0 and SSB0 with current conditions
    simpopF0 <- popdyn(nareas, maxage, Ncurr=N[x,,1,], 
                       pyears, M_age=M_ageArray[x,,], Asize_c=Asize[x,],
                       MatAge=Mat_age[x,,], 
                       WtAge=Wt_age[x,,], Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,], 
                       movc=mov[x,,], SRrelc=SRrel[x], 
                       Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
                       SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], Fapic=0, maxF=maxF, control=2) 
    
  } else {
    simpop <- popdynCPP(nareas, maxage, Ncurr=N[x,,1,], 
                        pyears, M_age=M_ageArray[x,,], Asize_c=Asize[x,],
                        MatAge=Mat_age[x,,], 
                        WtAge=Wt_age[x,,], Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,], 
                        movc=mov[x,,], SRrelc=SRrel[x], 
                        Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
                        SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], Qc=0, Fapic=exp(opt$minimum), maxF=maxF, control=2)
    # calculate B0 and SSB0 with current conditions
    simpopF0 <- popdynCPP(nareas, maxage, Ncurr=N[x,,1,], 
                          pyears, M_age=M_ageArray[x,,], Asize_c=Asize[x,],
                          MatAge=Mat_age[x,,], 
                          WtAge=Wt_age[x,,], Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,], 
                          movc=mov[x,,], SRrelc=SRrel[x], 
                          Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
                          SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], Qc=0, Fapic=0, maxF=maxF, control=2)
  }
  
  ## Cn <- simpop[[7]]/simpop[[8]] * simpop[[1]] * (1-exp(-simpop[[8]])) # retained catch 
  Cn <- simpop[[6]]/simpop[[8]] * simpop[[1]] * (1-exp(-simpop[[8]])) # removals
  Cb <- Cn[,pyears,] * Wt_age[x,,pyears]
  B <- sum(simpop[[2]][,pyears,] + Cb)
  
  SSB_MSY <- sum(simpop[[4]][,pyears,])
  
  V_BMSY <- sum(simpop[[5]][,pyears,])
  F_MSYv <- -log(1 - (MSY/(V_BMSY+MSY)))
  
  
  SSB0_curr <- sum(simpopF0[[4]][,pyears,])
  B0_curr <- sum(simpopF0[[2]][,pyears,])
  SSBMSY_SSB0 <- sum(simpop[[4]][,pyears,])/SSB0_curr
  BMSY_B0 <- sum(simpop[[2]][,pyears,])/B0_curr
  # SSBMSY_SSB0 <- sum(simpop[[4]][,pyears,])/SSB0[x]
  # BMSY_B0 <- sum(simpop[[2]][,pyears,])/B0[x]
  
  
  return(c(MSY = MSY, FMSY = F_MSYv, SSB = SSB_MSY, SSBMSY_SSB0=SSBMSY_SSB0, 
           BMSY_B0=BMSY_B0, B = B, VB=V_BMSY+MSY))
  
}




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
#' @param maxF A numeric value specifying the maximum fishing mortality for any single age class
#' @param useCPP logical - use the CPP code? For testing purposes only
#'
#' @export
#'
#' @author A. Hordyk
#' 
optMSY <- function(logFa, Asize_c, nareas, maxage, Ncurr, pyears, M_age, 
                 MatAge, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc, 
                 R0c, SSBpRc, aRc, bRc, Qc, maxF, useCPP=TRUE) {
  
  FMSYc <- exp(logFa)
  if(!useCPP) {
    simpop <- popdyn(nareas, maxage, Ncurr, pyears, M_age, Asize_c, 
                     MatAge, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc, 
                     R0c, SSBpRc, aRc, bRc, Qc, Fapic=FMSYc, maxF, control=2)
  } else {
    simpop <- popdynCPP(nareas, maxage, Ncurr, pyears, M_age, Asize_c,
                     MatAge, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc, 
                     R0c, SSBpRc, aRc, bRc, Qc=0, Fapic=FMSYc, maxF, control=2)
  }

  # Yield 
  # Cn <- simpop[[7]]/simpop[[8]] * simpop[[1]] * (1-exp(-simpop[[8]])) # retained catch
  Cn <- simpop[[6]]/simpop[[8]] * simpop[[1]] * (1-exp(-simpop[[8]])) # removals
  Cb <- Cn[,pyears,] * WtAge[,pyears]
  -sum(Cb)
}



#' Calculate Reference Yield 
#'
#' @param x Integer, the simulation number
#' @param Asize A matrix (dimensions nsim by nareas) with relative size of areas
#' @param nareas The number of spatial areas
#' @param maxage The maximum age
#' @param N Array of the numbers-at-age in population. Dimensions are nsim, maxage, nyears, nareas. 
#' Only values from the first year (i.e N[,,1,]) are used, which is the current N-at-age.
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
#' @param maxF A numeric value specifying the maximum fishing mortality for any single age class
#' @param useCPP logical - use the CPP code? For testing purposes only
#'
#' @author A. Hordyk
#' @export
#'
#' @author A. Hordyk
#' 
getFref3 <- function(x, Asize, nareas, maxage, N, pyears, M_ageArray, Mat_age, Wt_age,
                     V, retA, Perr, mov, SRrel, Find, Spat_targ, hs, R0a, SSBpR, aR, bR, 
                     maxF, useCPP=TRUE) {
  
  opt <- optimize(optMSY, log(c(0.001, 10)), Asize_c=Asize[x,], nareas, maxage, Ncurr=N[x,,1,], 
                  pyears, M_age=M_ageArray[x,,], MatAge=Mat_age[x,,], 
                  WtAge=Wt_age[x,,], Vuln=V[x,,], Retc=retA[x,,], Prec=Perr[x,], 
                  movc=mov[x,,], SRrelc=SRrel[x], 
                  Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
                  SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], maxF=maxF, useCPP=useCPP)
  
  -opt$objective
  
}
 
  
  
#' Apply output control recommendations and calculate population dynamics  
#'
#' @param y Projection year
#' @param Asize relative size of areas (matrix nsim by nareas)
#' @param TACused TAC recommendation
#' @param TAC_f Implementation error on TAC
#' @param lastCatch Catch from last year
#' @param availB Total available biomass
#' @param maxF Maximum fishing mortality
#' @param Biomass_P Numeric array (nsim, maxage, proyears, nareas) with Biomass at age
#' @param VBiomass_P Numeric array (nsim, maxage, proyears, nareas) with Vulnerable Biomass at age
#' @param CB_P Numeric array (nsim, maxage, proyears, nareas) with Catch Biomass at age
#' @param CB_Pret Numeric array (nsim, maxage, proyears, nareas) with Retained catch biomass at age
#' @param FM_P Numeric array (nsim, maxage, proyears, nareas) with fishing mortality at age
#' @param Z_P Numeric array (nsim, maxage, proyears, nareas) with total mortality at age
#' @param Spat_targ Spatial targetting
#' @param V_P Numeric array(nsim, maxage, nyears+proyears) with vulnerability at age
#' @param retA_P Numeric array(nsim, maxage, nyears+proyears) with retention at age
#' @param M_ageArray Numeric array (nsim, maxage, nyears+proyears) Natural mortality at age
#' @param qs Catchability coefficient
#' @param nyears Number of historical years
#' @param nsim Number of simulations 
#' @param maxage Maximum age
#' @param nareas Number of areas
#'
#' @export
#'
#' @author A. Hordyk
#' 
CalcOutput <- function(y, Asize, TACused, TAC_f, lastCatch, availB, maxF, Biomass_P, VBiomass_P, CB_P, CB_Pret,
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
  
  TACused[is.na(TACused)] <- lastCatch[is.na(TACused)] # if MP returns NA - TAC is set to catch from last year
  
  TACrec <- TACused             # TAC recommendation
  TACusedE<- TAC_f[,y]*TACused   # TAC taken after implementation error
  
  maxC <- (1 - exp(-maxF)) * availB # maximum catch given maxF
  TACusedE[TACusedE > maxC] <- maxC[TACusedE > maxC] # apply maxF limit - catch can't be higher than maxF * vulnerable biomass
  
  # fishdist <- (apply(VBiomass_P[, , y, ], c(1, 3), sum)^Spat_targ)/
    # apply(apply(VBiomass_P[, , y, ], c(1, 3), sum)^Spat_targ, 1, mean)  # spatial preference according to spatial biomass

  fishdist <- (apply(VBiomass_P[, , y, ], c(1, 3), sum)^Spat_targ)/
    apply(apply(VBiomass_P[, , y, ], c(1, 3), sum)^Spat_targ, 1, sum)  # spatial preference according to spatial biomass
  
  
 
  # If there is discard mortality, actual removals are higher than TACused
  # calculate distribution of all effort
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
 
  # calcFs <- lapply(1:nsim, getFs, y=y, Vuln=V_P, CB=CB_P, Bio=Biomass_P, Mage=M_ageArray, Fdist=fishdist,
  #        maxage=maxage, nareas=nareas, nyears=nyears) # numerically calculate Fs
  # 
  # 
  # FM_P[,,y,] <- aperm(array(unlist(calcFs, use.names=FALSE), dim=c(maxage, nareas, nsim)), c(3, 1, 2))
  # FM_P[,,y,][FM_P[,,y,] > (1-exp(-maxF))]  <- 1 - exp(-maxF)
  
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



#' Apply input control recommendations and calculate population dynamics  
#'
#' Internal function
#' 
#' @param y Simulation year
#' @param Asize Matrix (nsim by nareas) with relative size of areas
#' @param nyears Number of historical 
#' @param proyears Number of projection years
#' @param InputRecs Input control recommendations
#' @param nsim Number of simulations
#' @param nareas Number of areas
#' @param LR5_P Length at 5 percent retention
#' @param LFR_P Length at full retention
#' @param Rmaxlen_P Retention of maximum length
#' @param maxage Maximum age
#' @param retA_P Retention at age
#' @param retL_P Retention at length
#' @param V_P Realized vulnerability at age
#' @param V2 Gear vulnerability at age
#' @param pSLarray Realized vulnerability at length
#' @param SLarray2 Gear vulnerability at length
#' @param DR Discard ratio
#' @param maxlen maximum length
#' @param Len_age Length-at-age
#' @param CAL_binsmid Length-bin mid-points
#' @param Fdisc Fraction of discarded fish that die
#' @param nCALbins Number of length bins
#' @param E_f Implementation error on effort recommendation
#' @param SizeLim_f Implementation error on size limit
#' @param VBiomass_P Vulnerable biomass-at-age
#' @param Biomass_P Biomass-at-age
#' @param Spat_targ Spatial targetting
#' @param FinF Final fishing effort
#' @param qvar Annual ariability in catchability
#' @param qs Catchability
#' @param qinc Numeric vector (nsim) increased
#' @param CB_P Numeric array (nsim, maxage, proyears, nareas) Catch biomass at age
#' @param CB_Pret Numeric array (nsim, maxage, proyears, nareas) Retained catch biomass at age
#' @param FM_P Numeric array (nsim, maxage, proyears, nareas) Fishing mortality at age
#' @param FM_retain Numeric array (nsim, maxage, proyears, nareas) Retained fishing mortality at age
#' @param Z_P Numeric array (nsim, maxage, proyears, nareas) Total mortality at age
#' @param M_ageArray Numeric array (nsim, maxage, nyears+proyears) Natural mortality at age
#' @param LastEffort Numeric vector (nsim) with fishing effort from last year
#' @param LastSpatial Numeric matrix (nsim, nareas) with spatial closures from last year
#' @param LastAllocat Numeric vector (nsim) with allocation from last year
#'
#' @keywords internal
#' @export
#'
#' @author A. Hordyk
#' 
CalcInput <- function(y, Linf, Asize, nyears, proyears, InputRecs, nsim, nareas, LR5_P, LFR_P,
                      Rmaxlen_P, maxage, retA_P, retL_P, V_P, V2, pSLarray,
                      SLarray2, DR, maxlen, Len_age, CAL_binsmid, Fdisc, 
                      nCALbins, E_f, SizeLim_f, VBiomass_P, Biomass_P, Spat_targ,
                      FinF, qvar, qs, qinc, CB_P, CB_Pret, FM_P, FM_retain, Z_P,
                      M_ageArray, LastEffort, LastSpatial, LastAllocat) {
  
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
    if (y==1) Ei <- LastEffort  * E_f[,y] # effort is unchanged but has implementation error
    if (y>1) Ei <- LastEffort / E_f[,y-1]  * E_f[,y] # effort is unchanged but has implementation error
  } else if (length(InputRecs$Effort) != nsim) {
    stop("Effort recommmendation is not 'nsim' long.\n Does MP return Effort recommendation under all conditions?")
  } else {
    Ei <- InputRecs$Effort * E_f[,y] # effort adjustment with implementation error
  }
  
  # Spatial 
  if (all(is.na(InputRecs$Spatial))) { # no spatial recommendation 
    Si <- LastSpatial # matrix(1, nsim, nareas) # spatial is unchanged - modify this if spatial closure in historical years  
  } else if (any(is.na(InputRecs$Spatial))) {
    stop("Spatial recommmendation has some NAs.\n Does MP return Spatial recommendation under all conditions?")
  } else {
    Si <-InputRecs$Spatial # change spatial fishing
  }
  
  # Allocation 
  if (length(InputRecs$Allocate) == 0) { # no allocation recommendation
    Ai <- LastAllocat # rep(0, nsim) # allocation is unchanged 
  } else if (length(InputRecs$Allocate) != nsim) {
    stop("Allocate recommmendation is not 'nsim' long.\n Does MP return Allocate recommendation under all conditions?")
  } else {
    Ai <- InputRecs$Allocate # change in spatial allocation
  }
  # Retention Curve 
  RetentFlag <- FALSE
  # LR5 
  if (length(InputRecs$LR5) == 0) { # no  recommendation
    LR5_P[(y + nyears):(nyears+proyears),] <- matrix(LR5_P[y + nyears-1,], 
                                                     nrow=(length((y + nyears):(nyears+proyears))),
                                                     ncol=nsim, byrow=TRUE) # unchanged 

  } else if (length(InputRecs$LR5) != nsim) {
    stop("LR5 recommmendation is not 'nsim' long.\n Does MP return LR5 recommendation under all conditions?")
  } else {
    LR5_P[(y + nyears):(nyears+proyears),] <- matrix(InputRecs$LR5 * SizeLim_f[,y], 
                                                     nrow=(length((y + nyears):(nyears+proyears))),
                                                     ncol=nsim, byrow=TRUE) # recommendation with implementation error
    RetentFlag <- TRUE
  }
  # LFR 
  if (length(InputRecs$LFR) == 0) { # no  recommendation
    LFR_P[(y + nyears):(nyears+proyears),] <- matrix(LFR_P[y + nyears-1,], 
                                                     nrow=(length((y + nyears):(nyears+proyears))),
                                                     ncol=nsim, byrow=TRUE) # unchanged 
  } else if (length(InputRecs$LFR) != nsim) {
    stop("LFR recommmendation is not 'nsim' long.\n Does MP return LFR recommendation under all conditions?")
  } else {
    LFR_P[(y + nyears):(nyears+proyears),] <- matrix(InputRecs$LFR * SizeLim_f[,y], 
                                                     nrow=(length((y + nyears):(nyears+proyears))),
                                                     ncol=nsim, byrow=TRUE) # recommendation with implementation error
    RetentFlag <- TRUE
  }
  # Rmaxlen 
  if (length(InputRecs$Rmaxlen) == 0) { # no  recommendation
    Rmaxlen_P[(y + nyears):(nyears+proyears),] <- matrix(Rmaxlen_P[y + nyears-1,], 
                                                         nrow=(length((y + nyears):(nyears+proyears))),
                                                         ncol=nsim, byrow=TRUE)   # unchanged 
  
  } else if (length(Rmaxlen) != nsim) {
    stop("Rmaxlen recommmendation is not 'nsim' long.\n Does MP return Rmaxlen recommendation under all conditions?")
  } else {
    Rmaxlen_P[(y + nyears):(nyears+proyears),] <- matrix(InputRecs$Rmaxlen, 
                                                         nrow=(length((y + nyears):(nyears+proyears))),
                                                         ncol=nsim, byrow=TRUE) # recommendation
    RetentFlag <- TRUE
  }
  # HS - harvest slot 
 
  if (length(InputRecs$HS) == 0) { # no  recommendation
    HS <- rep(1E5, nsim) # no harvest slot 
  } else if (length(InputRecs$HS) != nsim) {
    stop("HS recommmendation is not 'nsim' long.\n Does MP return HS recommendation under all conditions?")
  } else {
    HS <- InputRecs$HS  * SizeLim_f[,y] # recommendation
    RetentFlag <- TRUE
  }
  # Change in retention - update vulnerability and retention curves 
  if (RetentFlag) {
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
    
    dr <- aperm(abind::abind(rep(list(DR), maxage), along=3), c(2,3,1))
    retA_P[,,allyrs] <- (1-dr[,,yr]) * retA_P[,,yr]
    dr <- aperm(abind::abind(rep(list(DR), nCALbins), along=3), c(2,3,1))
    retL_P[,,allyrs] <- (1-dr[,,yr]) * retL_P[,,yr]
    
    # update realized vulnerablity curve with retention and dead discarded fish 
    Fdisc_array1 <- array(Fdisc, dim=c(nsim, maxage, length(allyrs)))
    
    V_P[,,allyrs] <- V2[,,allyrs] * (retA_P[,,allyrs] + (1-retA_P[,,allyrs])*Fdisc_array1)
    
    Fdisc_array2 <- array(Fdisc, dim=c(nsim, nCALbins, length(allyrs)))
    pSLarray[,,allyrs]  <- SLarray2[,,allyrs] * (retL_P[,,allyrs]+ (1-retL_P[,,allyrs])*Fdisc_array2)
    
    # Realised Retention curves
    retA_P[,,allyrs] <- retA_P[,,allyrs] * V_P[,,allyrs]
    retL_P[,,allyrs] <- retL_P[,,allyrs] * pSLarray[,,allyrs] 
     
  }
  
  
  newVB <- apply(Biomass_P[, , y, ] * V_P[SAYt], c(1, 3), sum)  # calculate total vuln biomass by area 
  # fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ, 1, mean)  # spatial preference according to spatial vulnerable biomass
  fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ, 1, sum)  # spatial preference according to spatial vulnerable biomass
  Emult <- 1 + ((2/apply(fishdist * Si, 1, sum)) - 1) *   Ai  # allocate effort to new area according to fraction allocation Ai
 
  # fishing mortality with input control recommendation 
  FM_P[SAYR] <- (FinF[S1] * Ei[S1] * V_P[SAYt] * Si[SR] * fishdist[SR] * Emult[S1] * qvar[SY1] * (qs[S1]*(1 + qinc[S1]/100)^y))/Asize[SR]
  
  # retained fishing mortality with input control recommendation
  FM_retain[SAYR] <- (FinF[S1] * Ei[S1] * retA_P[SAYt] * Si[SR] * fishdist[SR] * Emult[S1] * qvar[SY1] * qs[S1]*(1 + qinc[S1]/100)^y)/Asize[SR]
  
  VBiomass_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt]  # update vulnerable biomass 
  Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt] # calculate total mortality 
  
  CB_P[SAYR] <- FM_P[SAYR]/Z_P[SAYR] * Biomass_P[SAYR] * (1 - exp(-Z_P[SAYR]))
  CB_Pret[SAYR] <- FM_retain[SAYR]/Z_P[SAYR] * Biomass_P[SAYR] * (1 - exp(-Z_P[SAYR]))
  
  out <- list() 
  out$Z_P <- Z_P 
  out$FM_P <- FM_P 
  out$FM_retain <- FM_retain
  out$CB_P <- CB_P
  out$CB_Pret <- CB_Pret
  out$Effort <- Ei 
  out$retA_P <- retA_P
  out$retL_P <- retL_P
  out$V_P <- V_P 
  out$pSLarray <- pSLarray
  out$Si <- Si
  out$Ai <- Ai
  out
  
}
