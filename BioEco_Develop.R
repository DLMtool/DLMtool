## TOP -----
BioEco <- new("BioEco")

OM <- new("OM", Stock=Albacore, 
          Fleet=Generic_Fleet, 
          Obs=Perfect_Info, 
          Imp=Perfect_Imp, 
          BioEco=BioEco)


OM@CostCurr <- c(1,1)
OM@RevCurr <- c(0.95,1.05)

OM@CostInc <- c(0,0)
OM@RevInc <- c(0,0)

OM@Response <- c(0.05,0.05)
OM@LatentEff<- c(0.4, 0.4)

OM@nsim <- 5
OM@interval <- 1
def.args <- DLMtool:::dev.mode(); for (nm in names(def.args)) assign(nm, def.args[[nm]])

## ---- Setup -----


if (class(OM) != "OM") stop("You must specify an operating model")
Misc<-new('list') #Blank miscellaneous slot created
if("seed"%in%slotNames(OM)) set.seed(OM@seed) # set seed for reproducibility 

OM <- updateMSE(OM)
tiny <- 1e-15  # define tiny variable

# Backwards compatible with DLMtool v < 4
if("nsim"%in%slotNames(OM))nsim<-OM@nsim
if("proyears"%in%slotNames(OM))proyears<-OM@proyears

# Backwards compatible with DLMtool v < 4.4.2
if(length(OM@interval)>0) interval <- OM@interval
if(length(OM@pstar)>0) pstar <- OM@pstar
if(length(OM@maxF)>0) maxF <- OM@maxF
if(length(OM@reps)>0) reps <- OM@reps

OM@interval <- interval 
OM@pstar <- pstar 
OM@maxF <- maxF 
OM@reps <- reps 

OM@nsim<-nsim # number of simulations
OM@proyears<-proyears # number of projection years
nyears <- OM@nyears  # number of historical years

OM <- ChkObj(OM) # Check that all required slots in OM object contain values 

if (proyears < 2) stop('OM@proyears must be > 1', call.=FALSE)
if(!silent) message("Loading operating model")

# --- Sample OM parameters ----
# Custom Parameters
# custom parameters exist - sample and write to list
SampCpars <- list()
SampCpars <- if(length(OM@cpars)>0) SampleCpars(OM@cpars, nsim, msg=!silent)

# Stock Parameters & assign to function environment
StockPars <- SampleStockPars(OM, nsim, nyears, proyears, SampCpars, msg=!silent)
for (X in 1:length(StockPars)) assign(names(StockPars)[X], StockPars[[X]])

# Fleet Parameters & assign to function environment
FleetPars <- SampleFleetPars(SubOM(OM, "Fleet"), Stock=StockPars, nsim, 
                             nyears, proyears, cpars=SampCpars)
for (X in 1:length(FleetPars)) assign(names(FleetPars)[X], FleetPars[[X]])

# Obs Parameters & assign to function environment
ObsPars <- SampleObsPars(OM, nsim, cpars=SampCpars)
for (X in 1:length(ObsPars)) assign(names(ObsPars)[X], ObsPars[[X]])

# Imp Parameters & assign to function environment
ImpPars <- SampleImpPars(OM, nsim, cpars=SampCpars)
for (X in 1:length(ImpPars)) assign(names(ImpPars)[X], ImpPars[[X]])

# BioEco Parameters & assign to function environment
BioEcoPars <- SampleBioEcoPars(OM, nsim, cpars=SampCpars)
for (X in 1:length(BioEcoPars)) assign(names(BioEcoPars)[X], BioEcoPars[[X]])

# --- Initialize Arrays ----
N <- array(NA, dim = c(nsim, maxage, nyears, nareas))  # stock numbers array
Biomass <- array(NA, dim = c(nsim, maxage, nyears, nareas))  # stock biomass array
VBiomass <- array(NA, dim = c(nsim, maxage, nyears, nareas))  # vulnerable biomass array
SSN <- array(NA, dim = c(nsim, maxage, nyears, nareas))  # spawning stock numbers array
SSB <- array(NA, dim = c(nsim, maxage, nyears, nareas))  # spawning stock biomass array
FM <- array(NA, dim = c(nsim, maxage, nyears, nareas))  # fishing mortality rate array
FMret <- array(NA, dim = c(nsim, maxage, nyears, nareas))  # fishing mortality rate array for retained fish 
Z <- array(NA, dim = c(nsim, maxage, nyears, nareas))  # total mortality rate array
SPR <- array(NA, dim = c(nsim, maxage, nyears)) # store the Spawning Potential Ratio
Agearray <- array(rep(1:maxage, each = nsim), dim = c(nsim, maxage))  # Age array

#  --- Pre Equilibrium calcs ----
# Survival array with M-at-age
surv <- matrix(1, nsim, maxage)
surv[, 2:maxage] <- t(exp(-apply(M_ageArray[,,1], 1, cumsum)))[, 1:(maxage-1)]  # Survival array
Nfrac <- surv * Mat_age[,,1]  # predicted Numbers of mature ages in first year

# Set up array indexes sim (S) age (A) year (Y) region/area (R)
SAYR <- as.matrix(expand.grid(1:nareas, 1, 1:maxage, 1:nsim)[4:1])  
SAY <- SAYR[, 1:3]
SAR <- SAYR[, c(1,2,4)]
SA <- Sa <- SAYR[, 1:2]
SR <- SAYR[, c(1, 4)]
S <- SAYR[, 1]
SY <- SAYR[, c(1, 3)]
Sa[,2]<- maxage-Sa[,2] + 1 # This is the process error index for initial year

# Calculate initial distribution if mov provided in cpars
if(!exists('initdist', inherits = FALSE)) { # movement matrix has been provided in cpars
  # Pinitdist is created in SampleStockPars instead of initdist if 
  # movement matrix is provided in cpars - OM@cpars$mov
  if (!exists('Asize', inherits = FALSE)) {
    message('Asize not set in cpars. Assuming all areas equal size')
    Asize <- matrix(1/nareas, nrow=nsim, ncol=nareas)
  }
  
  SSN[SAYR] <- Nfrac[SA] * R0[S] * Pinitdist[SR]  # Calculate initial spawning stock numbers
  N[SAYR] <- R0[S] * surv[SA] * Pinitdist[SR]  # Calculate initial stock numbers
  Neq <- N
  SSB[SAYR] <- SSN[SAYR] * Wt_age[SAY]    # Calculate spawning stock biomass
  SSB0 <- apply(SSB[, , 1, ], 1, sum)  # Calculate unfished spawning stock biomass
  SSBpR <- matrix(SSB0/R0, nrow=nsim, ncol=nareas)  # Spawning stock biomass per recruit
  SSB0a <- apply(SSB[, , 1, ], c(1, 3), sum)  # Calculate unfished spawning stock numbers
  
  bR <- matrix(log(5 * hs)/(0.8 * SSB0a), nrow=nsim)  # Ricker SR params
  aR <- matrix(exp(bR * SSB0a)/SSBpR, nrow=nsim)  # Ricker SR params
  R0a <- matrix(R0, nrow=nsim, ncol=nareas, byrow=FALSE) * 1/nareas # initial distribution of recruits
  
  # Project unfished for Nyrs to calculate equilibrium spatial distribution
  Nyrs <- ceiling(3 * maxage) # Project unfished for 3 x maxage
  # Set up projection arrays 
  M_ageArrayp <- array(M_ageArray[,,1], dim=c(dim(M_ageArray)[1:2], Nyrs))
  Wt_agep <- array(Wt_age[,,1], dim=c(dim(Wt_age)[1:2], Nyrs))
  Mat_agep <- array(Mat_age[,,1], dim=c(dim(Mat_age)[1:2], Nyrs))
  Perr_yp <- array(1, dim=c(dim(Perr_y)[1], Nyrs+maxage)) # no process error 
  
  # update mov if needed
  dimMov <- dim(mov)
  movp <- mov
  if (dimMov[length(dimMov)] < Nyrs) {
    movp <- array(movp, dim=c(dimMov[1:(length(dimMov)-1)], Nyrs))
  }
  
  # Not used but make the arrays anyway
  retAp <- array(retA[,,1], dim=c(dim(retA)[1:2], Nyrs))
  Vp <- array(V[,,1], dim=c(dim(V)[1:2], Nyrs))
  noMPA <- matrix(1, nrow=Nyrs, ncol=nareas)
  
  runProj <- lapply(1:nsim, projectEq, Asize, nareas=nareas, maxage=maxage, N=N, pyears=Nyrs,
                    M_ageArray=M_ageArrayp, Mat_age=Mat_agep, Wt_age=Wt_agep, V=Vp, retA=retAp,
                    Perr=Perr_yp, mov=movp, SRrel=SRrel, Find=Find, Spat_targ=Spat_targ, hs=hs,
                    R0a=R0a, SSBpR=SSBpR, aR=aR, bR=bR, SSB0=SSB0, B0=B0, MPA=noMPA, maxF=maxF,
                    Nyrs)
  Neq1 <- aperm(array(as.numeric(unlist(runProj)), dim=c(maxage, nareas, nsim)), c(3,1,2))  # unpack the list 
  
  # --- Equilibrium spatial / age structure (initdist by SAR)
  initdist <- Neq1/array(apply(Neq1, c(1,2), sum), dim=c(nsim, maxage, nareas))
  
  # check arrays and calculations
  if (checks) {
    if(!all(round(apply(initdist, c(1,2), sum),1)==1)) warning('initdist does not sum to one')
    if(!(all(round(apply(Neq[,,1,], 1, sum) /  apply(Neq1, 1, sum),1) ==1))) warning('eq age structure')
    sim <- sample(1:nsim,1)
    yrval <- sample(1:Nyrs,1)
    if (!all(M_ageArrayp[sim,,yrval] == M_ageArray[sim,,1] )) warning('problem with M_ageArrayp')
    if(!all(Wt_agep[sim,,yrval] == Wt_age[sim,,1]))  warning('problem with Wt_agep')
    if(!all(Mat_agep[sim,,yrval] == Mat_age[sim,,1])) warning('problem with Mat_agep')
    
  } 
}

# Unfished recruitment by area - INITDIST OF AGE 1.
R0a <- matrix(R0, nrow=nsim, ncol=nareas, byrow=FALSE) * initdist[,1,] # 

# ---- Unfished Equilibrium calcs ----
surv <- array(1, dim=c(nsim, maxage, nyears+proyears)) # unfished survival for every year
surv[, 2:maxage, ] <- aperm(exp(-apply(M_ageArray, c(1,3), cumsum))[1:(maxage-1), ,], c(2,1,3)) # Survival array
Nfrac <- surv * Mat_age  # predicted numbers of mature ages in all years

# indices for all years
SAYR_a <- as.matrix(expand.grid(1:nareas, 1:(nyears+proyears), 1:maxage, 1:nsim)[4:1])  
SAY_a <- SAYR_a[, 1:3]
SAR_a <- SAYR_a[, c(1,2,4)]
SA_a <- SAYR_a[, 1:2]
SR_a <- SAYR_a[, c(1, 4)]
S_a <- SAYR_a[, 1]
SY_a <- SAYR_a[, c(1, 3)]

# arrays for unfished biomass for all years 
SSN_a <- array(NA, dim = c(nsim, maxage, nyears+proyears, nareas))  
N_a <- array(NA, dim = c(nsim, maxage, nyears+proyears, nareas))
Biomass_a <- array(NA, dim = c(nsim, maxage, nyears+proyears, nareas))
SSB_a <- array(NA, dim = c(nsim, maxage, nyears+proyears, nareas))

SSN_a[SAYR_a] <- Nfrac[SAY_a] * R0[S_a] * initdist[SAR_a]  # Calculate initial spawning stock numbers for all years
N_a[SAYR_a] <- R0[S_a] * surv[SAY_a] * initdist[SAR_a] # Calculate initial stock numbers for all years
Biomass_a[SAYR_a] <- N_a[SAYR_a] * Wt_age[SAY_a]  # Calculate initial stock biomass
SSB_a[SAYR_a] <- SSN_a[SAYR_a] * Wt_age[SAY_a]    # Calculate spawning stock biomass

SSN0_a <- apply(SSN_a, c(1,3), sum) # unfished spawning numbers for each year
N0_a <- apply(N_a, c(1,3), sum) # unfished numbers for each year)
SSB0_a <- apply(SSB_a, c(1,3), sum) # unfished spawning biomass for each year
SSB0a_a <- apply(SSB_a, c(1, 3,4), sum)  # Calculate unfished spawning stock biomass by area for each year
B0_a <- apply(Biomass_a, c(1,3), sum) # unfished biomass for each year
VB0_a <- apply(apply(Biomass_a, c(1,2,3), sum) * V, c(1,3), sum) # unfished vulnerable biomass for each year

UnfishedByYear <- list(SSN0=SSN0_a, N0=N0_a, SSB0=SSB0_a, B0=B0_a, VB0=VB0_a)

if (quantile(ageM[,1],0.95) > nyears + proyears) {
  if(!silent) message('Note: number of historical year `nyears` + `proyears` is less than the highest age of maturity')
}

# ---- Unfished Reference Points ----
SSBpRa <- array(SSB0_a/matrix(R0, nrow=nsim, ncol=nyears+proyears), dim = c(nsim, nyears+proyears))

UnfishedRefs <- sapply(1:nsim, CalcUnfishedRefs, ageM=ageM, N0_a=N0_a, SSN0_a=SSN0_a,
                       SSB0_a=SSB0_a, B0_a=B0_a, VB0_a=VB0_a, SSBpRa=SSBpRa, SSB0a_a=SSB0a_a) 

N0 <- UnfishedRefs[1,] %>% unlist() # average unfished numbers
SSN0 <- UnfishedRefs[2,] %>% unlist() # average spawning unfished numbers
SSB0 <- UnfishedRefs[3,] %>% unlist() # average unfished spawning biomass
B0 <- UnfishedRefs[4,] %>% unlist() # average unfished biomass
VB0 <- UnfishedRefs[5,] %>% unlist() # average unfished biomass

# average spawning stock biomass per recruit 
SSBpR <- matrix(UnfishedRefs[6,] %>% unlist(), nrow=nsim, ncol=nareas) 
SSB0a <- UnfishedRefs[7,] %>% unlist() %>% matrix(nrow=nsim, ncol=nareas, byrow = TRUE)# average unfished biomass
bR <- matrix(log(5 * hs)/(0.8 * SSB0a), nrow=nsim)  # Ricker SR params
aR <- matrix(exp(bR * SSB0a)/SSBpR, nrow=nsim)  # Ricker SR params

Misc$Unfished <- list(Refs=UnfishedRefs, ByYear=UnfishedByYear)

# --- Optimize for Initial Depletion ----
# Depletion in year 1 
initD <- SampCpars$initD # 
if (!is.null(initD)) { # initial depletion is not unfished
  if (!silent) message("Optimizing for user-specified depletion in first historical year")
  Perrmulti <- sapply(1:nsim, optDfunwrap, initD=initD, Nfrac=Nfrac[,,1], R0=R0,
                      Perr_y=Perr_y, surv=surv[,,1], Wt_age=Wt_age, SSB0=SSB0,
                      maxage=maxage)
  Perr_y[,1:maxage] <- Perr_y[, 1:maxage] * Perrmulti
}

# --- Non-equilibrium calcs ----
SSN[SAYR] <- Nfrac[SAY] * R0[S] * initdist[SAR]*Perr_y[Sa]  # Calculate initial spawning stock numbers
N[SAYR] <- R0[S] * surv[SAY] * initdist[SAR]*Perr_y[Sa]  # Calculate initial stock numbers
Biomass[SAYR] <- N[SAYR] * Wt_age[SAY]  # Calculate initial stock biomass
SSB[SAYR] <- SSN[SAYR] * Wt_age[SAY]    # Calculate spawning stock biomass
VBiomass[SAYR] <- Biomass[SAYR] * V[SAY]  # Calculate vunerable biomass

if (checks && !is.null(initD)) { # check initial depletion 
  plot(apply(SSB[,,1,], 1, sum)/SSB0, initD)
  if (!any(round(apply(SSB[,,1,], 1, sum)/SSB0, 2) == round(initD,2))) warning('problem with initial depletion')
}

# --- Historical Spatial closures ----
MPA <- matrix(1, nyears+proyears, ncol=nareas) # fraction open to fishing 
if (all(!is.na(OM@MPA)) && sum(OM@MPA) != 0) { # historical spatial closures have been specified
  yrindex <- OM@MPA[,1]
  if (max(yrindex)>nyears) stop("Invalid year index for spatial closures: must be <= nyears")
  if (min(yrindex)<1) stop("Invalid year index for spatial closures: must be > 1")
  if (ncol(OM@MPA)-1 != nareas) stop("OM@MPA must be nareas + 1")
  for (xx in seq_along(yrindex)) {
    MPA[yrindex[xx]:nrow(MPA),] <- matrix(OM@MPA[xx, 2:ncol(OM@MPA)], nrow=length(yrindex[xx]:nrow(MPA)),ncol=nareas, byrow = TRUE)
  }
}

# --- Optimize catchability (q) to fit depletion ---- 
if(!silent) message("Optimizing for user-specified depletion in last historical year")
bounds <- c(0.0001, 15) # q bounds for optimizer
qs <- sapply(1:nsim, getq3, D, SSB0, nareas, maxage, N, pyears=nyears, 
             M_ageArray, Mat_age, Asize, Wt_age, V, retA, Perr_y, mov, SRrel, Find, 
             Spat_targ, hs, R0a, SSBpR, aR, bR, bounds=bounds, MPA=MPA, maxF=maxF) # find the q that gives current stock depletion

# --- Check that q optimizer has converged ---- 
LimBound <- c(1.1, 0.9)*range(bounds)  # bounds for q (catchability). Flag if bounded optimizer hits the bounds 
probQ <- which(qs > max(LimBound) | qs < min(LimBound))
Nprob <- length(probQ)

# If q has hit bound, re-sample depletion and try again. Tries 'ntrials' times and then alerts user
if (length(probQ) > 0) {
  Err <- TRUE
  if(!silent) message(Nprob,' simulations have final biomass that is not close to sampled depletion') 
  if(!silent) message('Re-sampling depletion, recruitment error, and fishing effort')
  
  count <- 0
  OM2 <- OM 
  while (Err & count < ntrials) {
    # Re-sample Stock Parameters 
    Nprob <- length(probQ)
    OM2@nsim <- Nprob
    SampCpars2 <- list()
    
    if (length(OM2@cpars)>0) SampCpars2 <- SampleCpars(OM2@cpars, OM2@nsim, msg=FALSE) 
    
    ResampStockPars <- SampleStockPars(OM2, cpars=SampCpars2, msg=FALSE)  
    ResampStockPars$CAL_bins <- StockPars$CAL_bins
    ResampStockPars$CAL_binsmid <- StockPars$CAL_binsmid 
    
    D[probQ] <- ResampStockPars$D  # Re-sampled depletion  
    procsd[probQ] <- ResampStockPars$procsd # Re-sampled recruitment deviations
    AC[probQ] <- ResampStockPars$AC
    Perr_y[probQ,] <- ResampStockPars$Perr_y
    hs[probQ] <- ResampStockPars$hs # Re-sampled steepness
    
    # Re-sample historical fishing effort
    ResampFleetPars <- SampleFleetPars(SubOM(OM2, "Fleet"), Stock=ResampStockPars, 
                                       OM2@nsim, nyears, proyears, cpars=SampCpars2)
    Esd[probQ] <- ResampFleetPars$Esd 
    Find[probQ, ] <- ResampFleetPars$Find
    dFfinal[probQ] <- ResampFleetPars$dFfinal
    
    # Optimize for q 
    qs[probQ] <- sapply(probQ, getq3, D, SSB0, nareas, maxage, N, pyears=nyears, 
                        M_ageArray, Mat_age, Asize, Wt_age, V, retA, Perr_y, mov, SRrel, Find, 
                        Spat_targ, hs, R0a, SSBpR, aR, bR, bounds=bounds, MPA=MPA, maxF=maxF)
    
    probQ <- which(qs > max(LimBound) | qs < min(LimBound))
    count <- count + 1 
    if (length(probQ) == 0) Err <- FALSE
  }
  if (Err) { # still a problem
    tooLow <- length(which(qs > max(LimBound)))
    tooHigh <- length(which(qs < min(LimBound)))
    prErr <- length(probQ)/nsim
    if (prErr > fracD & length(probQ) >= 1) {
      if (length(tooLow) > 0) message(tooLow, " sims can't get down to the lower bound on depletion")
      if (length(tooHigh) > 0) message(tooHigh, " sims can't get to the upper bound on depletion")
      if(!silent) message("More than ", fracD*100, "% of simulations can't get to the specified level of depletion with these Operating Model parameters")
      stop("Change OM@seed and try again for a complete new sample, modify the input parameters, or increase ntrials")
    } else {
      if (length(tooLow) > 0) message(tooLow, " sims can't get down to the lower bound on depletion")
      if (length(tooHigh) > 0) message(tooHigh, " sims can't get to the upper bound on depletion")
      if(!silent) message("More than ", 100-fracD*100, "% simulations can get to the sampled depletion.\nContinuing")
    }
  }
}

# --- Simulate historical years ----
if(!silent) message("Calculating historical stock and fishing dynamics")  # Print a progress update

if(!is.null(control$unfished)) { # generate unfished historical simulations
  if(!silent) message("Simulating unfished historical period")
  Hist <- TRUE
  CalcBlow <- FALSE
  qs <- rep(0, nsim) # no fishing
}

histYrs <- sapply(1:nsim, function(x) 
  popdynCPP(nareas, maxage, Ncurr=N[x,,1,], nyears,  
            M_age=M_ageArray[x,,], Asize_c=Asize[x,], MatAge=Mat_age[x,,], WtAge=Wt_age[x,,],
            Vuln=V[x,,], Retc=retA[x,,], Prec=Perr_y[x,], movc=split.along.dim(mov[x,,,,],4), 
            SRrelc=SRrel[x], 
            Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
            SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], Qc=qs[x], Fapic=0, MPA=MPA, maxF=maxF, 
            control=1, SSB0c=SSB0[x]))

N <- aperm(array(as.numeric(unlist(histYrs[1,], use.names=FALSE)), dim=c(maxage, nyears, nareas, nsim)), c(4,1,2,3))
Biomass <- aperm(array(as.numeric(unlist(histYrs[2,], use.names=FALSE)), dim=c(maxage, nyears, nareas, nsim)), c(4,1,2,3))
SSN <- aperm(array(as.numeric(unlist(histYrs[3,], use.names=FALSE)), dim=c(maxage, nyears, nareas, nsim)), c(4,1,2,3))
SSB <- aperm(array(as.numeric(unlist(histYrs[4,], use.names=FALSE)), dim=c(maxage, nyears, nareas, nsim)), c(4,1,2,3))
VBiomass <- aperm(array(as.numeric(unlist(histYrs[5,], use.names=FALSE)), dim=c(maxage, nyears, nareas, nsim)), c(4,1,2,3))
FM <- aperm(array(as.numeric(unlist(histYrs[6,], use.names=FALSE)), dim=c(maxage, nyears, nareas, nsim)), c(4,1,2,3))
FMret <- aperm(array(as.numeric(unlist(histYrs[7,], use.names=FALSE)), dim=c(maxage, nyears, nareas, nsim)), c(4,1,2,3))
Z <-aperm(array(as.numeric(unlist(histYrs[8,], use.names=FALSE)), dim=c(maxage, nyears, nareas, nsim)), c(4,1,2,3))

Depletion <- apply(SSB[,,nyears,],1,sum)/SSB0

# Check that depletion is correct
if (checks) {
  if (prod(round(D, 2)/ round(Depletion,2)) != 1) {
    print(cbind(round(D,2), round(Depletion,2)))
    warning("Possible problem in depletion calculations")
  } 
} 

# --- Calculate MSY statistics for each year ----
MSY_y <- array(0, dim=c(nsim, nyears+proyears)) # store MSY for each sim and year
FMSY_y <- MSY_y # store FMSY for each sim, and year
SSBMSY_y <- MSY_y # store SSBMSY for each sim, and year 
BMSY_y <- MSY_y # store BMSY for each sim, and year
VBMSY_y <- MSY_y # store VBMSY for each sim, and year 

if(!silent) message("Calculating MSY reference points for each year")
# average life-history parameters over 10 years
for (y in 1:(nyears+proyears)) {
  MSYrefsYr <- sapply(1:nsim, optMSY_eq, M_ageArray, Wt_age, Mat_age, V, maxage, R0, SRrel, hs, yr.ind=y)
  MSY_y[,y] <- MSYrefsYr[1, ]
  FMSY_y[,y] <- MSYrefsYr[2,]
  SSBMSY_y[,y] <- MSYrefsYr[3,]
  BMSY_y[,y] <- MSYrefsYr[6,]
  VBMSY_y[,y] <- MSYrefsYr[7,] 
}

# --- MSY reference points ----
MSYRefPoints <- sapply(1:nsim, CalcMSYRefs, MSY_y=MSY_y, FMSY_y=FMSY_y, 
                       SSBMSY_y=SSBMSY_y, BMSY_y=BMSY_y, VBMSY_y=VBMSY_y, 
                       ageM=ageM, OM=OM)

MSY <- MSYRefPoints[1,] %>% unlist() # record the MSY results (Vulnerable)
FMSY <- MSYRefPoints[2,] %>% unlist()  # instantaneous FMSY (Vulnerable)
SSBMSY <- MSYRefPoints[3,] %>% unlist()  # Spawning Stock Biomass at MSY
BMSY <- MSYRefPoints[4,] %>% unlist() # total biomass at MSY
VBMSY <- MSYRefPoints[5,] %>% unlist() # Biomass at MSY (Vulnerable)
UMSY <- MSY/VBMSY  # exploitation rate 
FMSY_M <- FMSY/M  # ratio of true FMSY to natural mortality rate M
SSBMSY_SSB0 <- SSBMSY/SSB0 # SSBMSY relative to unfished (SSB)
BMSY_B0 <- BMSY/B0 # Biomass relative to unfished (B0)
VBMSY_VB0 <- VBMSY/VB0 # VBiomass relative to unfished (VB0)

if (!AnnualMSY) {
  warning('AnnualMSY argument is deprecated. MSY metrics are always calculated by year.\n Use `MSE@SSB` or `MSE@B` and `MSE@Misc$MSYRefs$ByYear` for alternative methods to calculate B/BMSY')
}

if (checks) {
  Btemp <- apply(SSB, c(1,3), sum)
  x <- Btemp[,nyears]/SSBMSY
  y <-D/SSBMSY_SSB0
  plot(x,y, xlim=c(0,max(x)), ylim=c(0,max(y)), xlab="SSB/SSBMSY", ylab="D/SSBMSY_SSB0")
  lines(c(-10,10),c(-10,10))
}

# --- Calculate B-low ---- 
# (SSB where it takes MGThorizon x MGT to reach Bfrac of BMSY)
# Znow<-apply(Z[,,nyears,]*N[,,nyears,],1:2,sum)/apply(N[,,nyears,],1:2,sum)
# MGTsurv<-t(exp(-apply(Znow,1,cumsum)))
# MGT<-apply(Agearray*(Mat_age[,,nyears]*MGTsurv),1,sum)/apply(Mat_age[,,nyears]*MGTsurv,1,sum)

MarrayArea <- replicate(nareas, M_ageArray[,,1:nyears])
Mnow<-apply(MarrayArea[,,nyears,]*N[,,nyears,],1:2,sum)/apply(N[,,nyears,],1:2,sum)
MGTsurv<-t(exp(-apply(Mnow,1,cumsum)))
MGT<-apply(Agearray*(Mat_age[,,nyears]*MGTsurv),1,sum)/apply(Mat_age[,,nyears]*MGTsurv,1,sum)

Blow <- rep(NA,nsim)
if(CalcBlow){
  if(!silent) message("Calculating B-low reference points")            
  MGThorizon<-floor(HZN*MGT)
  Blow <- sapply(1:nsim,getBlow, N, Asize, SSBMSY,SSBpR, MPA, SSB0, nareas, retA, MGThorizon,
                 Find,Perr_y,M_ageArray,hs,Mat_age, Wt_age,R0a,V,nyears,maxage,mov,
                 Spat_targ,SRrel,aR,bR,Bfrac, maxF) 
}

# --- Calculate Reference Yield ----
if(!silent) message("Calculating reference yield - best fixed F strategy")  
RefY <- sapply(1:nsim, getFref3, Asize, nareas, maxage, N=N[,,nyears,, drop=FALSE], pyears=proyears, 
               M_ageArray=M_ageArray[,,(nyears):(nyears+proyears)], Mat_age[,,(nyears):(nyears+proyears)], 
               Wt_age=Wt_age[,,nyears:(nyears+proyears)], 
               V=V[, , (nyears + 1):(nyears + proyears), drop=FALSE], 
               retA=retA[, , (nyears + 1):(nyears + proyears), drop=FALSE],  
               Perr=Perr_y[,(nyears):(nyears+maxage+proyears-1)], mov, SRrel, Find, 
               Spat_targ, hs, R0a, SSBpR, aR, bR, MPA=MPA, maxF=maxF, SSB0=SSB0)

RefPoints <- data.frame(MSY=MSY, FMSY=FMSY, SSBMSY=SSBMSY, SSBMSY_SSB0=SSBMSY_SSB0,
                        BMSY_B0=BMSY_B0, BMSY=BMSY, VBMSY=VBMSY, UMSY=UMSY, VBMSY_VB0=VBMSY_VB0,
                        FMSY_M=FMSY_M, N0=N0, SSB0=SSB0, B0=B0, VB0=VB0, RefY=RefY, 
                        Blow=Blow, MGT=MGT, R0=R0)

Misc$MSYRefs <- list(Refs=RefPoints, ByYear=list(MSY=MSY_y, FMSY=FMSY_y,
                                                 SSBMSY=SSBMSY_y,
                                                 BMSY=BMSY_y,
                                                 VBMSY=VBMSY_y))

# --- Calculate Historical Catch ----
# Calculate catch-at-age 
CB <- Biomass * (1 - exp(-Z)) * (FM/Z)  # Catch in biomass (removed from population)
CB[is.na(CB)] <- 0

# Calculate retained-at-age
Cret <- apply(N * (1 - exp(-Z)) * (FMret/Z), c(1, 3, 2), sum)  # Retained catch in numbers
Cret[is.na(Cret)] <- 0
CBret <- Biomass * (1 - exp(-Z)) * (FMret/Z)  # Retained catch in biomass 
CBret[is.na(CBret)] <- 0

# --- Observation errors for all years ----
ErrList <- list()
ErrList$Cbiasa <- array(ObsPars$Cbias, c(nsim, nyears + proyears))  # Catch bias array
if (!is.null(control$Cbias_yr)) { # catch bias specified with control argument 
  Cbiasa <- matrix(1, nsim, nyears+proyears)
  Cbiasa[,control$yrs] <- control$Cbias_yr
  ErrList$Cbiasa <- Cbiasa
} 
# composite of bias and observation error
ErrList$Cerr <- array(rlnorm((nyears + proyears) * nsim, 
                             mconv(1, rep(ObsPars$Csd, (nyears + proyears))), 
                             sdconv(1, rep(ObsPars$Csd, nyears + proyears))), 
                      c(nsim, nyears + proyears))  
# Index error
ErrList$Ierr <- array(rlnorm((nyears + proyears) * nsim, 
                             mconv(1, rep(Isd, nyears + proyears)), 
                             sdconv(1, rep(Isd, nyears + proyears))), 
                      c(nsim, nyears + proyears))

# Simulate error in observed recruitment index 
ErrList$Recerr <- array(rlnorm((nyears + proyears) * nsim, mconv(1, rep(Recsd, (nyears + proyears))), 
                               sdconv(1, rep(Recsd, nyears + proyears))), c(nsim, nyears + proyears))

# --- Implementation error time series ----
TAC_f <- array(rlnorm(proyears * nsim, mconv(TACFrac, TACSD),
                      sdconv(TACFrac, TACSD)), c(nsim, proyears))  # composite of TAC fraction and error
E_f <- array(rlnorm(proyears * nsim, mconv(TAEFrac, TAESD),
                    sdconv(TAEFrac, TAESD)), c(nsim, proyears))  # composite of TAE fraction and error
SizeLim_f<-array(rlnorm(proyears * nsim, mconv(SizeLimFrac, SizeLimSD),
                        sdconv(SizeLimFrac, SizeLimSD)), c(nsim, proyears))  # composite of size limit fraction and error

# --- Populate Data object with Historical Data ---- 
Data <- makeData(Biomass, CBret, Cret, N, SSB, VBiomass, StockPars, 
                 FleetPars, ObsPars, ImpPars, RefPoints,
                 ErrList, OM, SampCpars, initD, control=control,
                 silent=silent)

# --- Add Real Indices to Data object (if they exist) & calculate error stats ----
templist <- addRealInd(Data, SampCpars, ErrList, Biomass, VBiomass, SSB, nsim,
                       nyears, proyears, silent=silent)
Data <- templist$Data # update 
ErrList <- templist$ErrList # update
Misc$RInd.stats <- ErrList$stats.df # return stats

ObsPars <- Data@Obs # Obs pars updated in makeData 
OMPars <- Data@OM
OMPars$qs <- qs

# --- Return Historical Simulations and Data from last historical year ----
if (Hist) { # Stop the model after historical simulations are complete
  if(!silent) message("Returning historical simulations")
  HistObj <- new("Hist")
  Data@Misc <- list()
  HistObj@Data <- Data 
  HistObj@Obs <- ObsPars
  om <- OMPars[,order(colnames(OMPars))]
  ind <- which(!colnames(om) %in% colnames(RefPoints))
  HistObj@OM <- om[,ind]
  HistObj@AtAge <- list(Length=Len_age, Weight=Wt_age, Select=V,
                        Retention=retA,
                        Maturity=Mat_age, N.Mortality=M_ageArray,
                        Nage=apply(N, 1:3, sum),
                        SSBage=apply(SSB, 1:3, sum),
                        FM=FM
  )
  nout <- t(apply(N, c(1, 3), sum)) 
  vb <- t(apply(VBiomass, c(1, 3), sum))
  b <- t(apply(Biomass, c(1, 3), sum))
  ssb <- t(apply(SSB, c(1, 3), sum))
  Cc <- t(apply(CB, c(1,3), sum))
  Ccret <- t(apply(CBret, c(1,3), sum))
  rec <- t(apply((N)[, 1, , ], c(1,2), sum))
  TSdata <- list(VB=t(vb), SSB=t(ssb), B=t(b), Removals=t(Cc), Catch=t(Ccret),
                 Rec=t(rec), N=t(nout),
                 Find=Find, Marray=Marray, RecDev=Perr_y)
  HistObj@TSdata <- TSdata
  HistObj@Ref <- RefPoints[,order(colnames(RefPoints))]
  HistObj@SampPars <- c(StockPars, FleetPars, ObsPars, ImpPars)
  HistObj@Misc <- Misc
  return(HistObj)	
}

# --- Check MPs ---- 
if (is.na(MPs[1])) CheckMPs <- TRUE
if (CheckMPs) MPs <- MPCheck(MPs, Data, timelimit, silent)
nMP <- length(MPs)  # the total number of methods used
if (nMP < 1) stop("No valid MPs found", call.=FALSE)

# --- Add nMP dimension to MSY stats  ----
MSY_y <- array(MSY_y, dim=c(nsim, nyears+proyears, nMP)) %>% aperm(c(1,3,2)) # store MSY for each sim, MP and year
FMSY_y <- array(FMSY_y, dim=c(nsim, nyears+proyears, nMP)) %>% aperm(c(1,3,2)) # store FMSY for each sim, MP and year
SSBMSY_y <- array(SSBMSY_y, dim=c(nsim, nyears+proyears, nMP)) %>% aperm(c(1,3,2)) # store SSBMSY for each sim, MP and year 
BMSY_y <- array(BMSY_y, dim=c(nsim, nyears+proyears, nMP)) %>% aperm(c(1,3,2)) # store BMSY for each sim, MP and year
VBMSY_y <- array(VBMSY_y, dim=c(nsim, nyears+proyears, nMP)) %>% aperm(c(1,3,2)) # store VBMSY for each sim, MP and year 

# Calculate management interval for each MP
if (length(interval) != nMP) interval <- rep(interval, nMP)[1:nMP]
if (!all(interval == interval[1])) {
  message("Variable management intervals:")
  df <- data.frame(MP=MPs,interval=interval)
  message(paste(capture.output(print(df)), collapse = "\n"))
}

# ---- Set-up arrays and objects for projections ----
MSElist <- list(Data)[rep(1, nMP)]  # create a data object for each method (they have identical historical data and branch in projected years)
B_BMSYa <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected B_BMSY
F_FMSYa <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected F_FMSY
Ba <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected Biomass
SSBa <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected SSB
VBa <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected vulnerable biomass
FMa <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected fishing mortality rate
Ca <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected removed catch
CaRet <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected retained catch
TACa <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected TAC recommendation
Effort <- array(NA, dim = c(nsim, nMP, proyears))  # store the Effort
PAAout <- array(NA, dim = c(nsim, nMP, maxage))  # store the population-at-age in last projection year
CAAout <- array(NA, dim = c(nsim, nMP, maxage))  # store the catch-at-age in last projection year
CALout <- array(NA, dim = c(nsim, nMP, nCALbins))  # store the population-at-length in last projection year
# SPRa <- array(NA,dim=c(nsim,nMP,proyears)) # store the Spawning Potential Ratio

# ---- Bio-Economics ----


# --- Begin loop over MPs ----
mm <- 1 # for debugging
    if(!silent) message(mm, "/", nMP, " Running MSE for ", MPs[mm]) 
    checkNA <- NA # save number of NAs
    
    # years management is updated
    upyrs <- seq(from=1, to=proyears, by=interval[mm]) # 1 + (0:(floor(proyears/interval[mm]) - 1)) * interval[mm] 
    
    # reset selectivity & retention parameters for projections
    L5_P <- L5  
    LFS_P <- LFS
    Vmaxlen_P <- Vmaxlen
    SLarray_P <- SLarray # selectivity at length array - projections
    V_P <- V  #  selectivity at age array - projections
    LR5_P <- LR5
    LFR_P <- LFR
    Rmaxlen_P <- Rmaxlen
    retA_P <- retA # retention at age array - projections
    retL_P <- retL # retention at length array - projections
    Fdisc_P <- Fdisc # Discard mortality for projectons 
    DR_P <- DR # Discard ratio for projections
    
    # projection arrays
    N_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    Biomass_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    VBiomass_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    SSN_P <-array(NA, dim = c(nsim, maxage, proyears, nareas))
    SSB_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    FM_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    FM_Pret <- array(NA, dim = c(nsim, maxage, proyears, nareas)) # retained F 
    Z_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    CB_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    CB_Pret <- array(NA, dim = c(nsim, maxage, proyears, nareas)) # retained catch 
    
    # indexes
    SAYRL <- as.matrix(expand.grid(1:nsim, 1:maxage, nyears, 1:nareas))  # Final historical year
    SAYRt <- as.matrix(expand.grid(1:nsim, 1:maxage, 1 + nyears, 1:nareas))  # Trajectory year
    SAYR <- as.matrix(expand.grid(1:nsim, 1:maxage, 1, 1:nareas))
    SYt <- SAYRt[, c(1, 3)]
    SAYt <- SAYRt[, 1:3]
    SR <- SAYR[, c(1, 4)]
    SA1 <- SAYR[, 1:2]
    S1 <- SAYR[, 1]
    SY1 <- SAYR[, c(1, 3)]
    SAY1 <- SAYRt[, 1:3]
    SYA <- as.matrix(expand.grid(1:nsim, 1, 1:maxage))  # Projection year
    SY <- SYA[, 1:2]
    SA <- SYA[, c(1, 3)]
    SAY <- SYA[, c(1, 3, 2)]
    S <- SYA[, 1]
    
    # -- First projection year ----
    y <- 1
    if(!silent) {
      cat("."); flush.console()
    }
    # Recruitment and movement in first year 
    NextYrN <- lapply(1:nsim, function(x)
      popdynOneTScpp(nareas, maxage, SSBcurr=colSums(SSB[x,,nyears, ]), Ncurr=N[x,,nyears,],
                     Zcurr=Z[x,,nyears,], PerrYr=Perr_y[x, nyears+maxage-1], hs=hs[x],
                     R0a=R0a[x,], SSBpR=SSBpR[x,], aR=aR[x,], bR=bR[x,],
                     mov=mov[x,,,,nyears+1], SRrel=SRrel[x]))
    
    # The stock at the beginning of projection period
    N_P[,,1,] <- aperm(array(unlist(NextYrN), dim=c(maxage, nareas, nsim, 1)), c(3,1,4,2))
    Biomass_P[SAYR] <- N_P[SAYR] * Wt_age[SAY1]  # Calculate biomass
    VBiomass_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt]  # Calculate vulnerable biomass
    SSN_P[SAYR] <- N_P[SAYR] * Mat_age[SAY1]  # Calculate spawning stock numbers
    SSB_P[SAYR] <- SSN_P[SAYR] * Wt_age[SAY1]
    
    # Update abundance estimates - used for FMSY ref methods so that FMSY is applied to current abundance
    M_array <- array(0.5*M_ageArray[,,nyears+y], dim=c(nsim, maxage, nareas))
    Atemp <- apply(VBiomass_P[, , y, ] * exp(-M_array), 1, sum) # Abundance (mid-year before fishing)
    MSElist[[mm]]@OM$A <- Atemp 
    
    # -- Apply MP in initial projection year ----
    runMP <- applyMP(Data=MSElist[[mm]], MPs = MPs[mm], reps = reps, silent=TRUE)  # Apply MP
    
    MPRecs <- runMP[[1]][[1]] # MP recommendations
    Data_p <- runMP[[2]] # Data object object with saved info from MP 
    Data_p@TAC <- MPRecs$TAC
    
    # calculate pstar quantile of TAC recommendation dist 
    TACused <- apply(Data_p@TAC, 2, quantile, p = pstar, na.rm = T) 
    checkNA[y] <- sum(is.na(TACused))
    LastEi <- rep(1,nsim) # no effort adjustment
    LastSpatial <- array(MPA[nyears,], dim=c(nareas, nsim)) # 
    LastAllocat <- rep(1, nsim) # default assumption of reallocation of effort to open areas
    LastCatch <- apply(CB[,,nyears,], 1, sum)
    
    
    
    # ---- CalcMPDynamics Function ----
    # MPCalcs <- CalcMPDynamics(MPRecs, y, nyears, proyears, nsim, Biomass_P, VBiomass_P,
    #                           LastEi, LastSpatial, LastAllocat, LastCatch,
    #                           TACused, maxF,
    #                           LR5_P, LFR_P, Rmaxlen_P, retL_P, retA_P,
    #                           L5_P, LFS_P, Vmaxlen_P, SLarray_P, V_P,
    #                           Fdisc_P, DR_P,
    #                           M_ageArray, FM_P, FM_Pret, Z_P, CB_P, CB_Pret,
    #                           TAC_f, E_f, SizeLim_f,
    #                           FinF, Spat_targ,
    #                           CAL_binsmid, Linf, Len_age, maxage, nareas, Asize, nCALbins,
    #                           qs, qvar, qinc,
    #                           RevPC, CostCurr, RevInc, CostInc)
  
      # Effort 
      if (length(MPRecs$Effort) == 0) { # no effort recommendation
        if (y==1) Ei <- LastEi * E_f[,y] # effort is unchanged but has implementation error
        if (y>1) Ei <- LastEi / E_f[,y-1]  * E_f[,y] # effort is unchanged but has implementation error
      } else if (length(MPRecs$Effort) != nsim) {
        stop("Effort recommmendation is not 'nsim' long.\n Does MP return Effort recommendation under all conditions?")
      } else {
        # an effort recommendation 
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
          SLarray_P[,, yy] <- selLen # calculate new selectivity at length curve   
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
          retL_P[,, yy] <- relLen  # calculate new retention at length curve 
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
      
      CurrentB <- Biomass_P[,,y,] # biomass at the beginning of year 
      CurrentVB <- array(NA, dim=dim(CurrentB))
      Catch_tot <- Catch_retain <- array(NA, dim=dim(CurrentB)) # catch this year arrays
      FMc <- Zc <- array(NA, dim=dim(CurrentB)) # fishing and total mortality this year
      
      # indices 
      SAYRL <- as.matrix(expand.grid(1:nsim, 1:maxage, nyears, 1:nareas))  # Final historical year
      SAYRt <- as.matrix(expand.grid(1:nsim, 1:maxage, y + nyears, 1:nareas))  # Trajectory year
      SAYR <- as.matrix(expand.grid(1:nsim, 1:maxage, y, 1:nareas))
      SAR <- SAYR[, c(1,2,4)]
      SAY <- SAYR[,c(1:3)]
      
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
      S <- SYA[, 1]
      
      CurrentVB[SAR] <- CurrentB[SAR] * V_P[SAYt] # update available biomass if selectivity has changed
      
      # Calculate fishing distribution if all areas were open 
      newVB <- apply(CurrentVB, c(1,3), sum) # calculate total vuln biomass by area 
      fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ, 1, sum)  # spatial preference according to spatial vulnerable biomass
      
      d1 <- t(Si) * fishdist  # distribution of fishing effort
      fracE <- apply(d1, 1, sum) # fraction of current effort in open areas
      fracE2 <- d1 * (fracE + (1-fracE) * Ai)/fracE # re-distribution of fishing effort accounting for re-allocation of effort
      fishdist <- fracE2 # fishing effort by area
      
      
      # ---- No TAC recommendation ----
      if (all(is.na(TACused))) { # no TAC has been set
        # fishing mortality with effort control recommendation
        FM_P[SAYR] <- (FinF[S1] * Ei[S1] * V_P[SAYt] * t(Si)[SR] * fishdist[SR] *
                         qvar[SY1] * (qs[S1]*(1 + qinc[S1]/100)^y))/Asize[SR]
        
        # retained fishing mortality with effort control recommendation
        FM_Pret[SAYR] <- (FinF[S1] * Ei[S1] * retA_P[SAYt] * t(Si)[SR] * fishdist[SR] *
                            qvar[SY1] * qs[S1]*(1 + qinc[S1]/100)^y)/Asize[SR]
        
        # Effort <- Ei * apply(fracE2, 1, sum) # change in catchability not included in effort calc: * qvar[,y] * ((1 + qinc/100)^y))
      }
      
      # ---- Apply TAC recommendation ----
      if (!all(is.na(TACused))) { # a TAC has been set
        # if MP returns NA - TAC is set to catch from last year
        TACused[is.na(TACused)] <- LastCatch[is.na(TACused)] 
        TACusedE <- TAC_f[,y]*TACused   # TAC taken after implementation error
        
        # Calculate total vulnerable biomass available mid-year accounting for any changes in selectivity &/or spatial closures
        M_array <- array(0.5*M_ageArray[,,nyears+y], dim=c(nsim, maxage, nareas))
        Atemp <- apply(CurrentVB * exp(-M_array), c(1,3), sum) # mid-year before fishing
        availB <- apply(Atemp * t(Si), 1, sum) # adjust for spatial closures
        
        # Calculate spatial distribution of catch - temporary values
        # get distribution across age and fishdist across space - ignore magnitude of effort or q increase 
        Catch_tot[SAR] <- (CurrentB[SAR]* V_P[SAYt] * fishdist[SR])/Asize[SR]
        Catch_retain[SAR] <- (CurrentB[SAR] * retA_P[SAYt] * fishdist[SR])/Asize[SR] 
        
        # Calculate total removals when Catch_retain == TAC - total removal > retained when discarding
        retained <- apply(Catch_retain, 1, sum) # retained - available biomass
        actualremovals <- apply(Catch_tot, 1, sum) # removals - available biomass
        ratio <- actualremovals/retained # ratio of actual removals to retained catch
        ratio[!is.finite(ratio)] <- 0 
        ratio[ratio>1E5] <- 1E5
        
        temp <- Catch_retain/apply(Catch_retain, 1, sum) # distribution by age & area of retained fish
        Catch_retain <- TACusedE * temp  # retained catch 
        temp <- Catch_tot/apply(Catch_tot, 1, sum) # distribution of removals
        Catch_tot <- TACusedE * ratio * temp # scale up total removals
        
        # total removals can't be more than available biomass
        chk <- apply(Catch_tot, 1, sum) > availB 
        if (sum(chk)>0) {
          c_temp <- apply(Catch_tot[chk,,, drop=FALSE], 1, sum)
          ratio_temp <- (availB[chk]/c_temp) * 0.99
          # scale total catches to 0.99 available biomass
          if (sum(chk)>1) Catch_tot[chk,, ] <- Catch_tot[chk,,] * array(ratio_temp, dim=c(sum(chk), maxage, nareas))
          if (sum(chk)==1) Catch_tot[chk,, ] <- Catch_tot[chk,,] * array(ratio_temp, dim=c(maxage, nareas))
        }
        
        # Populate catch arrays
        CB_P[SAYR] <- Catch_tot[SAR] 
        CB_Pret[SAYR] <- Catch_retain[SAR]
        
        # Calculate F by age class
        FM_P[SAYR] <- CB_P[SAYR]/(Biomass_P[SAYR] * exp(-M_ageArray[SAYt]/2)) # Pope's approximation)
        FM_Pret[SAYR] <- CB_Pret[SAYR]/(Biomass_P[SAYR] * exp(-M_ageArray[SAYt]/2))  # Pope's approximation
        # check where C > VB (for high M species this can happen)
        FM_P[SAYR][FM_P[SAYR] >= 1] <- 0.99
        FM_Pret[SAYR][FM_Pret[SAYR] >= 1] <- 0.99
        
        FM_P[SAYR] <- -log(1-FM_P[SAYR]) # convert to instantanous
        FM_Pret[SAYR] <- -log(1-FM_Pret[SAYR]) # convert to instantanous
      }
      
      # Apply maxF constraint 
      FM_P[SAYR][FM_P[SAYR] > maxF] <- maxF 
      FM_Pret[SAYR][FM_Pret[SAYR] > maxF] <- maxF
      Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt] # calculate total mortality
      
      # Update catches after maxF constraint
      CB_P[SAYR] <- (1-exp(-FM_P[SAYR])) * (Biomass_P[SAYR] * exp(-0.5*M_ageArray[SAYt]))
      CB_Pret[SAYR] <- (1-exp(-FM_Pret[SAYR])) * (Biomass_P[SAYR] * exp(-0.5*M_ageArray[SAYt]))
      
      # Calculate total fishing mortality & effort
      M_array <- array(0.5*M_ageArray[,,nyears+y], dim=c(nsim, maxage, nareas))
      Ftot <- suppressWarnings(-log(1-apply(CB_P[,,y,], 1, sum)/apply(VBiomass_P[,,y,] * exp(-M_array), 1, sum)))
      Ftot[!is.finite(Ftot)] <- maxF
      
      # effort relative to last historical if entire TAC was caught
      Effort_req <- Ftot/(FinF * qs*qvar[,y]* (1 + qinc/100)^y) * apply(fracE2, 1, sum) # effort required to get this TAC
      
      # Bio-economics 
      # Calculate Profit Margin from previous year
      if (y == 1) {
        RetainCatch <- apply(CBret[,,nyears,], 1, sum) # retained catch last historical year
        RevPC <- RevCurr/RetainCatch # Revenue per unit catch
        PMargin <- (RevPC * RetainCatch) - (CostCurr) # profit margin in last historical year
        LastEffort <- rep(1, nsim)
        CurrEffort <- LastEffort + Response*PMargin # Effort in first projection year
        CurrEffort[CurrEffort<0] <- 0
        ### TO DO - ADD MAX EFFORT #####
        
      } else {
        RetainCatch <- apply(CB_Pret[,,y-1,], 1, sum) # retained catch last year
        LastEffort <- CurrEffort
        CostLast <- LastEffort * CostCurr # cost of effort last year
        PMargin <- (RevPC * RetainCatch) - (CostLast) # profit margin last year
        CurrEffort <- LastEffort + Response*PMargin # Effort in first projection year
        CurrEffort[CurrEffort<0] <- 0
        ### TO DO - ADD MAX EFFORT #####
      }
      
      Effort_act <- Effort_req # actual effort is reduced if TAC is met
      excessEff <- Effort_req>CurrEffort # simulations where required effort > actual effort
      Effort_act[excessEff] <- CurrEffort[excessEff]
        
      # --- Re-calculate catch given actual effort ----
      # fishing mortality with actual effort 
      FM_P[SAYR] <- (FinF[S1] * Effort_act[S1] * V_P[SAYt] * t(Si)[SR] * fishdist[SR] *
                       qvar[SY1] * (qs[S1]*(1 + qinc[S1]/100)^y))/Asize[SR]
      
      # retained fishing mortality with actual effort 
      FM_Pret[SAYR] <- (FinF[S1] * Effort_act[S1] * retA_P[SAYt] * t(Si)[SR] * fishdist[SR] *
                          qvar[SY1] * qs[S1]*(1 + qinc[S1]/100)^y)/Asize[SR]
      
      # Apply maxF constraint 
      FM_P[SAYR][FM_P[SAYR] > maxF] <- maxF 
      FM_Pret[SAYR][FM_Pret[SAYR] > maxF] <- maxF
      Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt] # calculate total mortality
      
      # Update catches after maxF constraint
      CB_P[SAYR] <- (1-exp(-FM_P[SAYR])) * (Biomass_P[SAYR] * exp(-0.5*M_ageArray[SAYt]))
      CB_Pret[SAYR] <- (1-exp(-FM_Pret[SAYR])) * (Biomass_P[SAYR] * exp(-0.5*M_ageArray[SAYt]))
      
      # Calculate total fishing mortality & effort
      M_array <- array(0.5*M_ageArray[,,nyears+y], dim=c(nsim, maxage, nareas))
      Ftot <- suppressWarnings(-log(1-apply(CB_P[,,y,], 1, sum)/apply(VBiomass_P[,,y,] * exp(-M_array), 1, sum)))
      Ftot[!is.finite(Ftot)] <- maxF

      # --- an effort regulation also exists ----
      if (length(MPRecs$Effort) >0 | all(Ei != 1)) { 
        #Make sure Effort doesn't exceed regulated effort
        aboveE <- which(Effort > Ei)
        if (length(aboveE)>0) {
          Effort[aboveE] <- Ei[aboveE] * apply(fracE2, 1, sum)[aboveE]
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
          CB_P[SAYR] <- (1-exp(-FM_P[SAYR])) * (Biomass_P[SAYR] * exp(-0.5*M_ageArray[SAYt]))
          CB_Pret[SAYR] <- (1-exp(-FM_Pret[SAYR])) * (Biomass_P[SAYR] * exp(-0.5*M_ageArray[SAYt]))
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
      out$Ftot <- Ftot
      out
    
    
    