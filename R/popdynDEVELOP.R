 #' Calculate population dynamics from MP recommendation
# #' 
# #' An internal function to calculate the population dynamics for the next time
# #' step based on the recent MP recommendation
# #'
# #' @param MPRecs A named list of MP recommendations. The names are the same as `slotNames('Rec')`, except 
# #' for `Misc`. Each element in the list is a matrix. With the expection of `Spatial`, all elements in list 
# #' have `nrow=1` and `ncol=nsim`. `Spatial` has `nrow=nareas`. Matrices can be empty matrix, populated with all NAs 
# #' (both mean no change in management with respect to this element (e.g. `Effort`)), or populated with a recommendation.
# #' MPs must either return a recommendation or no recommendation for every simulation for a particular slot (i.e. cannot have some NA and some values). 
# #' @param y The projection year
# #' @param nyears The number of historical years
# #' @param proyears The number of projection years
# #' @param nsim The number of simulations
# #' @param Biomass_P An array with dimensions `nsim`, `maxage`, `proyears`, and `nareas` with total biomass in the projection years
# #' @param VBiomass_P An array with dimensions `nsim`, `maxage`, `proyears`, and `nareas` with vulnerable biomass in the projection years
# #' @param LastTAE A vector of length `nsim` with the most recent TAE
# #' @param LastSpatial A matrix of `nrow=nareas` and `ncol=nsim` with the most recent spatial management arrangements
# #' @param LastAllocat A vector of length `nsim` with the most recent allocation
# #' @param LastTAC A vector of length `nsim` with the most recent TAC
# #' @param TACused A vector of length `nsim` with the most recent TAC
# #' @param maxF A numeric value with maximum allowed F. From `OM@maxF`
# #' @param LR5_P A matrix with `nyears+proyears` rows and `nsim` columns with the first length at 5 percent retention.
# #' @param LFR_P A matrix with `nyears+proyears` rows and `nsim` columns with the first length at full retention.
# #' @param Rmaxlen_P A matrix with `nyears+proyears` rows and `nsim` columns with the retention at maximum length.
# #' @param retL_P An array with dimensions `nsim`, `nCALbins` and `nyears+proyears` with retention at length
# #' @param retA_P An array with dimensions `nsim`, `maxage` and `nyears+proyears` with retention at age
# #' @param L5_P A matrix with `nyears+proyears` rows and `nsim` columns with the first length at 5 percent selectivity
# #' @param LFS_P A matrix with `nyears+proyears` rows and `nsim` columns with the first length at full selectivity
# #' @param Vmaxlen_P A matrix with `nyears+proyears` rows and `nsim` columns with the selectivity at maximum length.
# #' @param SLarray_P An array with dimensions `nsim`, `nCALbins` and `nyears+proyears` with selectivity at length
# #' @param V_P An array with dimensions `nsim`, `maxage` and `nyears+proyears` with selectivity at age
# #' @param Fdisc_P  vector of length `nsim` with discard mortality. From `OM@Fdisc` but can be updated by MP (`Rec@Fdisc`)
# #' @param DR_P A matrix with `nyears+proyears` rows and `nsim` columns with the fraction discarded.
# #' @param M_ageArray An array with dimensions `nsim`, `maxage` and `nyears+proyears` with natural mortality at age
# #' @param FM_P An array with dimensions `nsim`, `maxage`, `proyears`, and `nareas` with total fishing mortality
# #' @param FM_Pret An array with dimensions `nsim`, `maxage`, `proyears`, and `nareas` with fishing mortality of the retained fish
# #' @param Z_P An array with dimensions `nsim`, `maxage`, `proyears`, and `nareas` with total mortality 
# #' @param CB_P An array with dimensions `nsim`, `maxage`, `proyears`, and `nareas` with total catch
# #' @param CB_Pret An array with dimensions `nsim`, `maxage`, `proyears`, and `nareas` with retained catch
# #' @param TAC_f A matrix with `nsim` rows and `proyears` columns with the TAC implementation error
# #' @param E_f A matrix with `nsim` rows and `proyears` columns with the effort implementation error
# #' @param SizeLim_f A matrix with `nsim` rows and `proyears` columns with the size limit implementation error
# #' @param FinF A numeric vector of length `nsim` with fishing mortality in the last historical year
# #' @param Spat_targ A numeric vector of length `nsim` with spatial targeting
# #' @param CAL_binsmid A numeric vector of length `nCALbins` with mid-points of the CAL bins
# #' @param Linf A numeric vector of length `nsim` with Linf (from `Stock@Linf`)
# #' @param Len_age An array with dimensions `nsim`, `maxage`, and `nyears+proyears` with length-at-age
# #' @param maxage A numeric value with maximum age from `Stock@maxage`
# #' @param nareas A numeric value with number of areas
# #' @param Asize A matrix with `nsim` rows and `nareas` columns with the relative size of each area
# #' @param nCALbins The number of CAL bins. Should be the same as `length(CAL_binsmid)`
# #' @param qs A numeric vector of length `nsim` with catchability coefficient
# #' @param qvar A matrix with `nsim` rows and `proyears` columns with catchability variability 
# #' @param qinc A numeric vector of length `nsim` with average annual change in catchability
# #' @param checks Logical. Run internal checks? Currently not used. 
# #'
# #' @return A named list with updated population dynamics
# #' @author A. Hordyk
# #' @export
# #'
# #' @keywords internal
# CalcMPDynamics <- function(MPRecs, y, nyears, proyears, nsim, Biomass_P,
#                            VBiomass_P,
#                            LastTAE, histTAE, LastSpatial, LastAllocat, LastTAC,
#                            TACused, maxF,
#                            LR5_P, LFR_P, Rmaxlen_P, retL_P, retA_P,
#                            L5_P, LFS_P, Vmaxlen_P, SLarray_P, V_P,
#                            Fdisc_P, DR_P,
#                            M_ageArray, FM_P, FM_Pret, Z_P, CB_P, CB_Pret,
#                            TAC_f, E_f, SizeLim_f,
#                            FinF, Spat_targ,
#                            CAL_binsmid, Linf, Len_age, maxage, nareas, Asize, nCALbins,
#                            qs, qvar, qinc, 
#                            Effort_pot,
#                            checks=FALSE) {
#   # Effort 
#   if (length(MPRecs$Effort) == 0) { # no max effort recommendation
#     if (y==1) TAE <- LastTAE * E_f[,y] # max effort is unchanged but has implementation error
#     if (y>1) TAE <- LastTAE / E_f[,y-1]  * E_f[,y] # max effort is unchanged but has implementation error
#   } else if (length(MPRecs$Effort) != nsim) {
#     stop("Effort recommmendation is not 'nsim' long.\n Does MP return Effort recommendation under all conditions?")
#   } else {
#     # a maximum effort recommendation
#     if (!all(is.na(histTAE))) {
#       TAE <- histTAE * MPRecs$Effort * E_f[,y] # adjust existing TAE adjustment with implementation error  
#     } else {
#       TAE <- MPRecs$Effort * E_f[,y] # adjust existing TAE adjustment with implementation error
#     }
#   }
#   
#   # Spatial 
#   if (all(is.na(MPRecs$Spatial))) { # no spatial recommendation 
#     Si <- LastSpatial # spatial is unchanged 
#   } else if (any(is.na(MPRecs$Spatial))) {
#     stop("Spatial recommmendation has some NAs.\n Does MP return Spatial recommendation under all conditions?")
#   } else {
#     Si <- MPRecs$Spatial # change spatial fishing
#   }
#   if (all(dim(Si) != c(nareas, nsim))) stop("Spatial recommmendation not nareas long")
#   
#   # Allocation 
#   if (length(MPRecs$Allocate) == 0) { # no allocation recommendation
#     Ai <- LastAllocat # allocation is unchanged 
#   } else if (length(MPRecs$Allocate) != nsim) {
#     stop("Allocate recommmendation is not 'nsim' long.\n Does MP return Allocate recommendation under all conditions?")
#   } else {
#     Ai <- MPRecs$Allocate # change in spatial allocation
#   }
#   Ai <- as.numeric(Ai)
#   
#   # Retention Curve
#   RetentFlag <- FALSE # should retention curve be updated for future years?
#   # LR5 
#   if (length(MPRecs$LR5) == 0) { # no  recommendation
#     LR5_P[(y + nyears):(nyears+proyears),] <- matrix(LR5_P[y + nyears-1,], 
#                                                      nrow=(length((y + nyears):(nyears+proyears))),
#                                                      ncol=nsim, byrow=TRUE) # unchanged 
#     
#   } else if (length(MPRecs$LR5) != nsim) {
#     stop("LR5 recommmendation is not 'nsim' long.\n Does MP return LR5 recommendation under all conditions?")
#   } else {
#     LR5_P[(y + nyears):(nyears+proyears),] <- matrix(MPRecs$LR5 * SizeLim_f[,y], 
#                                                      nrow=(length((y + nyears):(nyears+proyears))),
#                                                      ncol=nsim, byrow=TRUE) # recommendation with implementation error
#     RetentFlag <- TRUE
#   }
#   # LFR 
#   if (length(MPRecs$LFR) == 0) { # no  recommendation
#     LFR_P[(y + nyears):(nyears+proyears),] <- matrix(LFR_P[y + nyears-1,], 
#                                                      nrow=(length((y + nyears):(nyears+proyears))),
#                                                      ncol=nsim, byrow=TRUE) # unchanged 
#   } else if (length(MPRecs$LFR) != nsim) {
#     stop("LFR recommmendation is not 'nsim' long.\n Does MP return LFR recommendation under all conditions?")
#   } else {
#     LFR_P[(y + nyears):(nyears+proyears),] <- matrix(MPRecs$LFR * SizeLim_f[,y], 
#                                                      nrow=(length((y + nyears):(nyears+proyears))),
#                                                      ncol=nsim, byrow=TRUE) # recommendation with implementation error
#     RetentFlag <- TRUE
#   }
#   # Rmaxlen 
#   if (length(MPRecs$Rmaxlen) == 0) { # no  recommendation
#     Rmaxlen_P[(y + nyears):(nyears+proyears),] <- matrix(Rmaxlen_P[y + nyears-1,], 
#                                                          nrow=(length((y + nyears):(nyears+proyears))),
#                                                          ncol=nsim, byrow=TRUE)   # unchanged 
#     
#   } else if (length(MPRecs$Rmaxlen) != nsim) {
#     stop("Rmaxlen recommmendation is not 'nsim' long.\n Does MP return Rmaxlen recommendation under all conditions?")
#   } else {
#     Rmaxlen_P[(y + nyears):(nyears+proyears),] <- matrix(MPRecs$Rmaxlen, 
#                                                          nrow=(length((y + nyears):(nyears+proyears))),
#                                                          ncol=nsim, byrow=TRUE) # recommendation
#     RetentFlag <- TRUE
#   }
#   
#   # HS - harvest slot 
#   if (length(MPRecs$HS) == 0) { # no  recommendation
#     HS <- rep(1E5, nsim) # no harvest slot 
#   } else if (length(MPRecs$HS) != nsim) {
#     stop("HS recommmendation is not 'nsim' long.\n Does MP return HS recommendation under all conditions?")
#   } else {
#     HS <- MPRecs$HS  * SizeLim_f[,y] # recommendation
#     RetentFlag <- TRUE
#   }
#   
#   # Selectivity Curve
#   SelectFlag <- FALSE # has selectivity been updated?
#   # L5 
#   if (length(MPRecs$L5) == 0) { # no  recommendation
#     L5_P[(y + nyears):(nyears+proyears),] <- matrix(L5_P[y + nyears-1,], 
#                                                     nrow=(length((y + nyears):(nyears+proyears))),
#                                                     ncol=nsim, byrow=TRUE) # unchanged 
#     
#   } else if (length(MPRecs$L5) != nsim) {
#     stop("L5 recommmendation is not 'nsim' long.\n Does MP return L5 recommendation under all conditions?")
#   } else {
#     L5_P[(y + nyears):(nyears+proyears),] <- matrix(MPRecs$L5 * SizeLim_f[,y], 
#                                                     nrow=(length((y + nyears):(nyears+proyears))),
#                                                     ncol=nsim, byrow=TRUE) # recommendation with implementation error
#     SelectFlag <- TRUE
#   }
#   # LFS
#   if (length(MPRecs$LFS) == 0) { # no  recommendation
#     LFS_P[(y + nyears):(nyears+proyears),] <- matrix(LFS_P[y + nyears-1,], 
#                                                      nrow=(length((y + nyears):(nyears+proyears))),
#                                                      ncol=nsim, byrow=TRUE) # unchanged 
#   } else if (length(MPRecs$LFS) != nsim) {
#     stop("LFS recommmendation is not 'nsim' long.\n Does MP return LFS recommendation under all conditions?")
#   } else {
#     LFS_P[(y + nyears):(nyears+proyears),] <- matrix(MPRecs$LFS * SizeLim_f[,y], 
#                                                      nrow=(length((y + nyears):(nyears+proyears))),
#                                                      ncol=nsim, byrow=TRUE) # recommendation with implementation error
#     SelectFlag <- TRUE
#   }
#   # Vmaxlen 
#   if (length(MPRecs$Rmaxlen) == 0) { # no  recommendation
#     Vmaxlen_P[(y + nyears):(nyears+proyears),] <- matrix(Vmaxlen_P[y + nyears-1,], 
#                                                          nrow=(length((y + nyears):(nyears+proyears))),
#                                                          ncol=nsim, byrow=TRUE)   # unchanged 
#     
#   } else if (length(MPRecs$Rmaxlen) != nsim) {
#     stop("Rmaxlen recommmendation is not 'nsim' long.\n Does MP return Rmaxlen recommendation under all conditions?")
#   } else {
#     Vmaxlen_P[(y + nyears):(nyears+proyears),] <- matrix(MPRecs$Vmaxlen, 
#                                                          nrow=(length((y + nyears):(nyears+proyears))),
#                                                          ncol=nsim, byrow=TRUE) # recommendation
#     SelectFlag <- TRUE
#   }
#   
#   # Discard Mortality 
#   if (length(MPRecs$Fdisc) >0) { # Fdisc has changed
#     if (length(MPRecs$Fdisc) != nsim) stop("Fdisc recommmendation is not 'nsim' long.\n Does MP return Fdisc recommendation under all conditions?")
#     Fdisc_P <- MPRecs$Fdisc
#   }
#   
#   # Discard Ratio 
#   if (length(MPRecs$DR)>0) { # DR has changed
#     if (length(MPRecs$DR) != nsim) stop("DR recommmendation is not 'nsim' long.\n Does MP return DR recommendation under all conditions?")
#     DR_P[(y+nyears):(nyears+proyears),] <- matrix(MPRecs$DR, nrow=length((y+nyears):(nyears+proyears)), ncol=nsim, byrow=TRUE) 
#   }
#   
#   # Update Selectivity and Retention Curve 
#   if (SelectFlag | RetentFlag) {
#     yr <- y+nyears 
#     allyrs <- (y+nyears):(nyears+proyears)  # update vulnerabilty for all future years
#     
#     srs <- (Linf - LFS_P[yr,]) / ((-log(Vmaxlen_P[yr,],2))^0.5) # descending limb
#     srs[!is.finite(srs)] <- Inf
#     sls <- (LFS_P[yr,] - L5_P[yr,]) / ((-log(0.05,2))^0.5) # ascending limb
#     
#     CAL_binsmidMat <- matrix(CAL_binsmid, nrow=nsim, ncol=length(CAL_binsmid), byrow=TRUE)
#     selLen <- t(sapply(1:nsim, getsel, lens=CAL_binsmidMat, lfs=LFS_P[yr,], sls=sls, srs=srs))
#     
#     for (yy in allyrs) {
#       # calculate new selectivity at age curve 
#       V_P[ , , yy] <- t(sapply(1:nsim, getsel, lens=Len_age[,,yy], lfs=LFS_P[yy,], sls=sls, srs=srs))
#       SLarray_P[,, yy] <- selLen # calculate new selectivity at length curve   
#     }
#     
#     # sim <- 158
#     # plot(CAL_binsmid, selLen[sim,], type="b")
#     # lines(c(L5_P[yr,sim], L5_P[yr,sim]), c(0, 0.05), lty=2)
#     # lines(c(LFS_P[yr,sim], LFS_P[yr,sim]), c(0, 1), lty=2)
#     # lines(c(Linf[sim], Linf[sim]), c(0, Vmaxlen_P[yr,sim]), lty=2)
#     
#     # calculate new retention curve
#     yr <- y+nyears 
#     allyrs <- (y+nyears):(nyears+proyears)  # update vulnerabilty for all future years
#     
#     srs <- (Linf - LFR_P[yr,]) / ((-log(Rmaxlen_P[yr,],2))^0.5) # selectivity parameters are constant for all years
#     srs[!is.finite(srs)] <- Inf
#     sls <- (LFR_P[yr,] - LR5_P[yr,]) / ((-log(0.05,2))^0.5)
#     
#     CAL_binsmidMat <- matrix(CAL_binsmid, nrow=nsim, ncol=length(CAL_binsmid), byrow=TRUE)
#     relLen <- t(sapply(1:nsim, getsel, lens=CAL_binsmidMat, lfs=LFR_P[yr,], sls=sls, srs=srs))
#     
#     for (yy in allyrs) {
#       # calculate new retention at age curve 
#       retA_P[ , , yy] <- t(sapply(1:nsim, getsel, lens=Len_age[,,yy], lfs=LFR_P[yy,], sls=sls, srs=srs))
#       retL_P[,, yy] <- relLen  # calculate new retention at length curve 
#     }
#     
#     # upper harvest slot 
#     aboveHS <- Len_age[,,allyrs, drop=FALSE]>array(HS, dim=c(nsim, maxage, length(allyrs)))
#     tretA_P <- retA_P[,,allyrs]
#     tretA_P[aboveHS] <- 0
#     retA_P[,,allyrs] <- tretA_P
#     for (ss in 1:nsim) {
#       index <- which(CAL_binsmid >= HS[ss])
#       retL_P[ss, index, allyrs] <- 0
#     }	
#     
#     dr <- aperm(abind::abind(rep(list(DR_P), maxage), along=3), c(2,3,1))
#     retA_P[,,allyrs] <- (1-dr[,,yr]) * retA_P[,,yr]
#     dr <- aperm(abind::abind(rep(list(DR_P), nCALbins), along=3), c(2,3,1))
#     retL_P[,,allyrs] <- (1-dr[,,yr]) * retL_P[,,yr]
#     
#     # update realized vulnerablity curve with retention and dead discarded fish 
#     Fdisc_array1 <- array(Fdisc_P, dim=c(nsim, maxage, length(allyrs)))
#     
#     V_P[,,allyrs] <- V_P[,,allyrs, drop=FALSE] * (retA_P[,,allyrs, drop=FALSE] + (1-retA_P[,,allyrs, drop=FALSE])*Fdisc_array1)
#     
#     Fdisc_array2 <- array(Fdisc_P, dim=c(nsim, nCALbins, length(allyrs)))
#     SLarray_P[,,allyrs]  <- SLarray_P[,,allyrs, drop=FALSE] * (retL_P[,,allyrs, drop=FALSE]+ (1-retL_P[,,allyrs, drop=FALSE])*Fdisc_array2)
#     
#     # Realised Retention curves
#     retA_P[,,allyrs] <- retA_P[,,allyrs] * V_P[,,allyrs]
#     retL_P[,,allyrs] <- retL_P[,,allyrs] * SLarray_P[,,allyrs] 
#   }
#   
#   CurrentB <- Biomass_P[,,y,] # biomass at the beginning of year 
#   CurrentVB <- array(NA, dim=dim(CurrentB))
#   Catch_tot <- Catch_retain <- array(NA, dim=dim(CurrentB)) # catch this year arrays
#   FMc <- Zc <- array(NA, dim=dim(CurrentB)) # fishing and total mortality this year
#   
#   # indices 
#   SAYRL <- as.matrix(expand.grid(1:nsim, 1:maxage, nyears, 1:nareas))  # Final historical year
#   SAYRt <- as.matrix(expand.grid(1:nsim, 1:maxage, y + nyears, 1:nareas))  # Trajectory year
#   SAYR <- as.matrix(expand.grid(1:nsim, 1:maxage, y, 1:nareas))
#   SAR <- SAYR[, c(1,2,4)]
#   SAY <- SAYR[,c(1:3)]
#   
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
#   S <- SYA[, 1]
#   
#   CurrentVB[SAR] <- CurrentB[SAR] * V_P[SAYt] # update available biomass if selectivity has changed
#   
#   # Calculate fishing distribution if all areas were open 
#   newVB <- apply(CurrentVB, c(1,3), sum) # calculate total vuln biomass by area 
#   fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ, 1, sum)  # spatial preference according to spatial vulnerable biomass
#   
#   d1 <- t(Si) * fishdist  # distribution of fishing effort
#   fracE <- apply(d1, 1, sum) # fraction of current effort in open areas
#   fracE2 <- d1 * (fracE + (1-fracE) * Ai)/fracE # re-distribution of fishing effort accounting for re-allocation of effort
#   fishdist <- fracE2 # fishing effort by area
#   
#   # ---- no TAC - calculate F with bio-economic effort ----
#   if (all(is.na(TACused))) {
#     if (all(is.na(Effort_pot)) & all(is.na(TAE))) Effort_pot <- rep(1, nsim) # historical effort
#     if (all(is.na(Effort_pot))) Effort_pot <- TAE[1,]
#     # fishing mortality with bio-economic effort
#     FM_P[SAYR] <- (FinF[S1] * Effort_pot[S1] * V_P[SAYt] * t(Si)[SR] * fishdist[SR] *
#                      qvar[SY1] * (qs[S1]*(1 + qinc[S1]/100)^y))/Asize[SR]
#     
#     # retained fishing mortality with bio-economic effort
#     FM_Pret[SAYR] <- (FinF[S1] * Effort_pot[S1] * retA_P[SAYt] * t(Si)[SR] * fishdist[SR] *
#                         qvar[SY1] * qs[S1]*(1 + qinc[S1]/100)^y)/Asize[SR]
#   }
#   
#   # ---- calculate required F and effort for TAC recommendation ----
#   if (!all(is.na(TACused))) { # a TAC has been set
#     # if MP returns NA - TAC is set to TAC from last year
#     TACused[is.na(TACused)] <- LastTAC[is.na(TACused)] 
#     TACusedE <- TAC_f[,y]*TACused   # TAC taken after implementation error
#     
#     # Calculate total vulnerable biomass available mid-year accounting for any changes in selectivity &/or spatial closures
#     M_array <- array(0.5*M_ageArray[,,nyears+y], dim=c(nsim, maxage, nareas))
#     Atemp <- apply(CurrentVB * exp(-M_array), c(1,3), sum) # mid-year before fishing
#     availB <- apply(Atemp * t(Si), 1, sum) # adjust for spatial closures
#     
#     
#     optFun <- function(Fapic, Va, fishdist, Ma, Ba, TACusedE, Asize, opt=1) {
#       Fa <- exp(Fapic) * Va
#       nareas <- length(fishdist)
#       Far <- (matrix(Fa, ncol=nareas, nrow=maxage) * fishdist)/Asize
#       Zar <- Far + matrix(Ma, ncol=nareas, nrow=maxage)
#       
#       predC <- Far/Zar * (1-exp(-Zar)) * Ba
#       nll <- (sum(predC) - TACusedE)^2
#       if (opt==1) return(nll)
#       if (opt==2) return(predC)
#       if (opt==3) return(Far=Far)
#       
#     }
#     
#     # Calculate apical F when retained catch == TAC
#     opt <- sapply(1:nsim, function(x)
#       optimize(optFun, interval=log(c(0.01, 2)), 
#                Va=retA_P[x,,nyears+y], fishdist=fishdist[x,],
#                Ma=M_ageArray[x,,nyears+y], Ba=Biomass_P[x,,y,], 
#                TACusedE=TACusedE[x], Asize=Asize[x,], tol=1E-5))
#     
#     # Retained catch-at-age & area 
#     ret_catches <- lapply(1:nsim, function(x)
#       optFun(Fapic=opt[,x]$minimum, 
#              Va=retA_P[x,,nyears+y], fishdist=fishdist[x,],
#              Ma=M_ageArray[x,,nyears+y], Ba=Biomass_P[x,,y,], 
#              TACusedE=TACusedE[x], Asize=Asize[x,], opt=2))
#     
#     # Total catch-at-age & area 
#     vuln_catches <- lapply(1:nsim, function(x)
#       optFun(Fapic=opt[,x]$minimum, 
#              Va=V_P[x,,nyears+y], fishdist=fishdist[x,],
#              Ma=M_ageArray[x,,nyears+y], Ba=Biomass_P[x,,y,], 
#              TACusedE=TACusedE[x], Asize=Asize[x,], opt=2))
#     
#     Catch_retain <- array(unlist(ret_catches), dim=c(maxage, nareas, nsim)) %>% aperm(., c(3,1,2))
#     Catch_tot <- array(unlist(vuln_catches), dim=c(maxage, nareas, nsim)) %>% aperm(., c(3,1,2))
#     
#     # Calculate total removals when Catch_retain == TAC - total removal > retained when discarding
#     retained <- apply(Catch_retain, 1, sum) # retained - available biomass
#     actualremovals <- apply(Catch_tot, 1, sum) # removals - available biomass
#     ratio <- actualremovals/retained # ratio of actual removals to retained catch
#     ratio[!is.finite(ratio)] <- 0 
#     ratio[ratio>1E5] <- 1E5
#     
#     
#     # total removals can't be more than available biomass
#     chk <- apply(Catch_tot, 1, sum) > availB 
#     if (sum(chk)>0) {
#       c_temp <- apply(Catch_tot[chk,,, drop=FALSE], 1, sum)
#       ratio_temp <- (availB[chk]/c_temp) * 0.99
#       # scale total catches to 0.99 available biomass
#       if (sum(chk)>1) Catch_tot[chk,, ] <- Catch_tot[chk,,] * array(ratio_temp, dim=c(sum(chk), maxage, nareas))
#       if (sum(chk)==1) Catch_tot[chk,, ] <- Catch_tot[chk,,] * array(ratio_temp, dim=c(maxage, nareas))
#     }
#     
#     # Populate catch arrays
#     CB_P[SAYR] <- Catch_tot[SAR] 
#     CB_Pret[SAYR] <- Catch_retain[SAR]
#     
#     
#     # Retained F-at-age & area 
#     ret_F <- lapply(1:nsim, function(x)
#       optFun(Fapic=opt[,x]$minimum, 
#              Va=retA_P[x,,nyears+y], fishdist=fishdist[x,],
#              Ma=M_ageArray[x,,nyears+y], Ba=Biomass_P[x,,y,], 
#              TACusedE=TACusedE[x], Asize=Asize[x,], opt=3))
#     
#     # Total F-at-age & area 
#     vuln_F <- lapply(1:nsim, function(x)
#       optFun(Fapic=opt[,x]$minimum, 
#              Va=V_P[x,,nyears+y], fishdist=fishdist[x,],
#              Ma=M_ageArray[x,,nyears+y], Ba=Biomass_P[x,,y,], 
#              TACusedE=TACusedE[x], Asize=Asize[x,], opt=3))
#     
#     FM_Pret[,,y,] <- array(unlist(ret_F), dim=c(maxage, nareas, nsim)) %>% aperm(., c(3,1,2))
#     FM_P[,,y,] <- array(unlist(vuln_F), dim=c(maxage, nareas, nsim)) %>% aperm(., c(3,1,2))
#   }
#   
#   # Apply maxF constraint 
#   FM_P[SAYR][FM_P[SAYR] > maxF] <- maxF 
#   FM_Pret[SAYR][FM_Pret[SAYR] > maxF] <- maxF
#   Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt] # calculate total mortality
#   
#   # Update catches after maxF constraint
#   CB_P[SAYR] <- FM_P[SAYR]/Z_P[SAYR] * (1-exp(-Z_P[SAYR])) * Biomass_P[SAYR]
#   CB_Pret[SAYR] <- FM_Pret[SAYR]/Z_P[SAYR] * (1-exp(-Z_P[SAYR])) * Biomass_P[SAYR]
#   
#   # Calculate total fishing mortality & effort
#   # M_array <- array(0.5*M_ageArray[,,nyears+y], dim=c(nsim, maxage, nareas))
#   # Ftot <- suppressWarnings(-log(1-apply(CB_P[,,y,], 1, sum)/apply(CurrentVB * exp(-M_array), 1, sum)))
#   # Ftot[!is.finite(Ftot)] <- maxF
#   Ftot <- apply(FM_P[,,y,], 1, max)
#   
#   # Effort_req - effort required to catch TAC
#   # Effort_pot - potential effort this year (active fishers) from bio-economic model
#   # Effort_act - actual effort this year
#   # TAE - maximum actual effort limit
#   # Effort_act < Effort_pot if Effort_req < Effort_pot
#   
#   # Effort relative to last historical with this potential catch
#   Effort_req <- Ftot/(FinF * qs*qvar[,y]* (1 + qinc/100)^y) * apply(fracE2, 1, sum) # effort required for this catch
#   
#   # Limit effort to potential effort from bio-economic model
#   Effort_act <- Effort_req
#   if (!all(is.na(Effort_pot))) {
#     excessEff <- Effort_req>Effort_pot # simulations where required effort > potential effort
#     Effort_act[excessEff] <- Effort_pot[excessEff] # actual effort can't be more than bio-economic effort
#   }
#   
#   # Limit actual effort <= TAE 
#   if (!all(is.na(TAE))) { # a TAE exists
#     Effort_act[Effort_act>TAE] <- TAE[Effort_act>TAE]
#   }
#   Effort_act[Effort_act<=0] <- tiny
#   
#   # --- Re-calculate catch given actual effort ----
#   # fishing mortality with actual effort 
#   FM_P[SAYR] <- (FinF[S1] * Effort_act[S1] * V_P[SAYt] * t(Si)[SR] * fishdist[SR] *
#                    qvar[SY1] * (qs[S1]*(1 + qinc[S1]/100)^y))/Asize[SR]
#   
#   # retained fishing mortality with actual effort 
#   FM_Pret[SAYR] <- (FinF[S1] * Effort_act[S1] * retA_P[SAYt] * t(Si)[SR] * fishdist[SR] *
#                       qvar[SY1] * qs[S1]*(1 + qinc[S1]/100)^y)/Asize[SR]
#   
#   # Apply maxF constraint 
#   FM_P[SAYR][FM_P[SAYR] > maxF] <- maxF 
#   FM_Pret[SAYR][FM_Pret[SAYR] > maxF] <- maxF
#   Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt] # calculate total mortality
#   
#   # Update catches after maxF constraint
#   CB_P[SAYR] <- FM_P[SAYR]/Z_P[SAYR] * (1-exp(-Z_P[SAYR])) * Biomass_P[SAYR]
#   CB_Pret[SAYR] <- FM_Pret[SAYR]/Z_P[SAYR] * (1-exp(-Z_P[SAYR])) * Biomass_P[SAYR]
#   
#   # # Apply maxF constraint 
#   # FM_P[SAYR][FM_P[SAYR] > maxF] <- maxF 
#   # FM_Pret[SAYR][FM_Pret[SAYR] > maxF] <- maxF
#   # Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt] # calculate total mortality
#   # 
#   # # Update catches after maxF constraint
#   # CB_P[SAYR] <- (1-exp(-FM_P[SAYR])) * (Biomass_P[SAYR] * exp(-0.5*M_ageArray[SAYt]))
#   # CB_Pret[SAYR] <- (1-exp(-FM_Pret[SAYR])) * (Biomass_P[SAYR] * exp(-0.5*M_ageArray[SAYt]))
#   # 
#   # # Calculate total fishing mortality & effort
#   # M_array <- array(0.5*M_ageArray[,,nyears+y], dim=c(nsim, maxage, nareas))
#   # Ftot <- suppressWarnings(-log(1-apply(CB_P[,,y,], 1, sum)/apply(VBiomass_P[,,y,] * exp(-M_array), 1, sum)))
#   # Ftot[!is.finite(Ftot)] <- maxF
#   
#   # Returns
#   out <- list()
#   out$TACrec <- TACused
#   out$V_P <- V_P
#   out$SLarray_P <- SLarray_P
#   out$retA_P <- retA_P
#   out$retL_P <- retL_P
#   out$Fdisc_P <- Fdisc_P
#   out$VBiomass_ <- VBiomass_P
#   out$Z_P <- Z_P
#   out$FM_P <- FM_P
#   out$FM_Pret <- FM_Pret
#   out$CB_P <- CB_P
#   out$CB_Pret <- CB_Pret
#   out$Si <- Si
#   out$Ai <- Ai
#   out$TAE <- TAE
#   out$Effort <- Effort_act # actual effort this year
#   out$Ftot <- Ftot
#   out
# }
