
## - TOP ------
library(DLMtool); library(dplyr)

# def.args <- DLMtool:::dev.mode(); for (nm in names(def.args)) assign(nm, def.args[[nm]])


OM_annual <- readRDS("saveOM.rdata")
OM_annual@interval <- 1
OM_annual@K <- range(OM_annual@cpars$K)
OM_annual@Linf <- range(OM_annual@cpars$Linf)
OM_annual@cpars <- list()


OM_annual@K <- c(0.10, 0.10)
OM_annual@M <- c(0.2, 0.2)
OM_annual@Linf  <- c(100, 100)

OM_annual@L50 <- c(60, 60)
OM_annual@L50_95 <- c(5,5)

OM_annual@L5 <- c(18.127, 18.127)
OM_annual@LFS <- c(46, 46)
OM_annual@Vmaxlen <- c(1,1)
OM_annual@LR5 <- c(0,0)
OM_annual@LFR <- c(0.01,0.01)
OM_annual@Rmaxlen <- c(1,1)

OM_annual@D <- c(0.2, 0.5)
OM_annual@maxage <- 14
OM_annual@h <- c(0.9, 0.9)

OM_annual@nsim <- 10


OM_annual <- tinyErr(OM_annual)

OM <- ChangeTS(OM_annual)
OM <- ChangeTS(OM)

recVec <- c(0.8, 0.05, 0.025, 0.125)

# OM <- OM_annual


OM@Perr <- c(0, 0)
OM@AC <- c(0,0)
OM@Linfsd <- c(0,0)
OM@Ksd <- c(0,0)



# OM@Msd <- c(0.1, 0.1)

checks = FALSE

# recVec <- c(0.1, 0.2, 0.3, 0.4)

# recVec <- rep(0.25, 4)

# sum(recVec)

# recVec <- rep(1, OM@cpars$nts) # c(0.1, 0.2, 0.3, 0.4)

# OM <- OM_annual
# recVec <- 1

# recVec should be first time-step
# equilibrium by year - loop over time-steps


# Dev Setup ####
# development mode - assign default argument values to current workspace if they don't exist
# def.args <- DLMtool:::dev.mode(); for (nm in names(def.args)) assign(nm, def.args[[nm]])

if (class(OM) != "OM") stop("You must specify an operating model", call.=FALSE)
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
if(!silent) message(crayon::black("Loading operating model"))

# --- Sample OM parameters ----
# Check for time-step parameter
if (!is.null(OM@cpars$nts)) {
  nts <- OM@cpars$nts
  message(nts) # TO DO #######
} else {
  OM@cpars$nts <- 1
  nts <- 1
  recVec <- 1
}

recVec <- recVec/sum(recVec)
recTS <- rep(recVec, nyears+proyears) # recruitment per time-step

# ---- Time-step Parameters ----
histnTS <- nyears * nts # number of time-steps in historical period
projnTS <- proyears * nts # number of time-steps in projection period

# Custom Parameters
# custom parameters exist - sample and write to list
SampCpars <- list()
SampCpars <- if(length(OM@cpars)>0) SampleCpars(OM@cpars, nsim, msg=!silent)

# Stock Parameters & assign to function environment
StockPars <- SampleStockPars(SubOM(OM, "Stock"), nsim, histnTS, projnTS, SampCpars, 
                             msg=!silent)
for (X in 1:length(StockPars)) assign(names(StockPars)[X], StockPars[[X]])

# Fleet Parameters & assign to function environment
FleetPars <- SampleFleetPars(SubOM(OM, "Fleet"), Stock=StockPars, nsim, 
                             histnTS, projnTS, cpars=SampCpars)
for (X in 1:length(FleetPars)) assign(names(FleetPars)[X], FleetPars[[X]])

# Obs Parameters & assign to function environment
ObsPars <- SampleObsPars(OM, nsim, cpars=SampCpars)
for (X in 1:length(ObsPars)) assign(names(ObsPars)[X], ObsPars[[X]])

# Imp Paramerers & assign to function environment
ImpPars <- SampleImpPars(OM, nsim, cpars=SampCpars)
for (X in 1:length(ImpPars)) assign(names(ImpPars)[X], ImpPars[[X]])

# Bio-Economic Parameters
BioEcoPars <- c("RevCurr", "CostCurr", "Response", "CostInc", "RevInc", "LatentEff")
if (all(lapply(SampCpars[BioEcoPars], length) == 0)) {
  # no bio-economic model
  # if (!silent) message("No bio-economic model parameters found. \nTAC and TAE assumed to be caught in full")
  RevCurr <- CostCurr <- Response <- CostInc <- RevInc <- LatentEff <- rep(NA, nsim)
} else {
  if (!silent) message("Bio-economic model parameters found.")
  # Checks
  if (length(SampCpars$CostCurr) != nsim) stop("OM@cpars$CostCurr is not length OM@nsim", call.=FALSE)
  if (length(SampCpars$RevCurr) != nsim) stop("OM@cpars$RevCurr is not length OM@nsim", call.=FALSE)
  if (length(SampCpars$Response) != nsim) stop("OM@cpars$Response is not length OM@nsim", call.=FALSE)
  if (length(SampCpars$RevInc) != nsim) SampCpars$RevInc <- rep(0, nsim)
  if (length(SampCpars$CostInc) != nsim) SampCpars$CostInc <- rep(0, nsim)
  if (length(SampCpars$LatentEff) != nsim) SampCpars$LatentEff <- rep(NA, nsim)
  RevCurr <- SampCpars$RevCurr
  CostCurr <- SampCpars$CostCurr
  Response <- SampCpars$Response
  CostInc <- SampCpars$CostInc
  RevInc <- SampCpars$RevInc
  LatentEff <- SampCpars$LatentEff
}

# --- Initialize Arrays ----
N <- array(NA, dim = c(nsim, maxage, histnTS, nareas))  # stock numbers array
Biomass <- array(NA, dim = c(nsim, maxage, histnTS, nareas))  # stock biomass array
VBiomass <- array(NA, dim = c(nsim, maxage, histnTS, nareas))  # vulnerable biomass array
SSN <- array(NA, dim = c(nsim, maxage, histnTS, nareas))  # spawning stock numbers array
SSB <- array(NA, dim = c(nsim, maxage, histnTS, nareas))  # spawning stock biomass array
FM <- array(NA, dim = c(nsim, maxage, histnTS, nareas))  # fishing mortality rate array
FMret <- array(NA, dim = c(nsim, maxage, histnTS, nareas))  # fishing mortality rate array for retained fish 
Z <- array(NA, dim = c(nsim, maxage, histnTS, nareas))  # total mortality rate array
SPR <- array(NA, dim = c(nsim, maxage, histnTS)) # store the Spawning Potential Ratio
Agearray <- array(rep(1:maxage, each = nsim), dim = c(nsim, maxage))  # Age array

#  --- Pre Equilibrium calcs ----

# FIX  get spatial distribution  in popdynCPP.cpp
##### TO DO ####

# # Survival array with M-at-age
# surv <- matrix(1, nsim, maxage)
# surv[, 2:maxage] <- t(exp(-apply(M_ageArray[,,1], 1, cumsum)))[, 1:(maxage-1)]  # Survival array
# Nfrac <- surv * Mat_age[,,1]  # predicted Numbers of mature ages in first year
# 
# # Set up array indexes sim (S) age (A) year (Y) region/area (R)
# SAYR <- as.matrix(expand.grid(1:nareas, 1, 1:maxage, 1:nsim)[4:1])  
# SAY <- SAYR[, 1:3]
# SAR <- SAYR[, c(1,2,4)]
# SA <- Sa <- SAYR[, 1:2]
# SR <- SAYR[, c(1, 4)]
# S <- SAYR[, 1]
# SY <- SAYR[, c(1, 3)]
# Sa[,2]<- maxage-Sa[,2] + 1 # This is the process error index for initial year
# 
# 
# # Calculate initial distribution if mov provided in cpars
# 
# if(!exists('initdist', inherits = FALSE)) { # movement matrix has been provided in cpars
#   if (nts !=1) stop("Not developed yet for models with sub-years")
#   # Pinitdist is created in SampleStockPars instead of initdist if 
#   # movement matrix is provided in cpars - OM@cpars$mov
#   if (!exists('Asize', inherits = FALSE)) {
#     message('Asize not set in cpars. Assuming all areas equal size')
#     Asize <- matrix(1/nareas, nrow=nsim, ncol=nareas)
#   }
#   
#   SSN[SAYR] <- Nfrac[SA] * R0[S] * Pinitdist[SR]  # Calculate initial spawning stock numbers
#   N[SAYR] <- R0[S] * surv[SA] * Pinitdist[SR]  # Calculate initial stock numbers
#   Neq <- N
#   SSB[SAYR] <- SSN[SAYR] * Wt_age[SAY]    # Calculate spawning stock biomass
#   SSB0 <- apply(SSB[, , 1, ], 1, sum)  # Calculate unfished spawning stock biomass
#   SSBpR <- matrix(SSB0/R0, nrow=nsim, ncol=nareas)  # Spawning stock biomass per recruit
#   SSB0a <- apply(SSB[, , 1, ], c(1, 3), sum)  # Calculate unfished spawning stock numbers
#   
#   bR <- matrix(log(5 * hs)/(0.8 * SSB0a), nrow=nsim)  # Ricker SR params
#   aR <- matrix(exp(bR * SSB0a)/SSBpR, nrow=nsim)  # Ricker SR params
#   R0a <- matrix(R0, nrow=nsim, ncol=nareas, byrow=FALSE) * 1/nareas # initial distribution of recruits
#   
#   # Project unfished for Nyrs to calculate equilibrium spatial distribution
#   Nyrs <- ceiling(3 * maxage) # Project unfished for 3 x maxage
#   # Set up projection arrays 
#   M_ageArrayp <- array(M_ageArray[,,1], dim=c(dim(M_ageArray)[1:2], Nyrs))
#   Wt_agep <- array(Wt_age[,,1], dim=c(dim(Wt_age)[1:2], Nyrs))
#   Mat_agep <- array(Mat_age[,,1], dim=c(dim(Mat_age)[1:2], Nyrs))
#   Perr_yp <- array(1, dim=c(dim(Perr_y)[1], Nyrs+maxage)) # no process error 
#   
#   # update mov if needed
#   dimMov <- dim(mov)
#   movp <- mov
#   if (dimMov[length(dimMov)] < Nyrs) {
#     movp <- array(movp, dim=c(dimMov[1:(length(dimMov)-1)], Nyrs))
#   }
#   
#   # Not used but make the arrays anyway
#   retAp <- array(retA[,,1], dim=c(dim(retA)[1:2], Nyrs))
#   Vp <- array(V[,,1], dim=c(dim(V)[1:2], Nyrs))
#   noMPA <- matrix(1, nrow=Nyrs, ncol=nareas)
#   
#   runProj <- lapply(1:nsim, projectEq, Asize, nareas=nareas, maxage=maxage, N=N, pyears=Nyrs,
#                     M_ageArray=M_ageArrayp, Mat_age=Mat_agep, Wt_age=Wt_agep, V=Vp, retA=retAp,
#                     Perr=Perr_yp, mov=movp, SRrel=SRrel, Find=Find, Spat_targ=Spat_targ, hs=hs,
#                     R0a=R0a, SSBpR=SSBpR, aR=aR, bR=bR, SSB0=SSB0, B0=B0, MPA=noMPA, maxF=maxF,
#                     Nyrs)
#   Neq1 <- aperm(array(as.numeric(unlist(runProj)), dim=c(maxage, nareas, nsim)), c(3,1,2))  # unpack the list 
#   
#   # --- Equilibrium spatial / age structure (initdist by SAR)
#   initdist <- Neq1/array(apply(Neq1, c(1,2), sum), dim=c(nsim, maxage, nareas))
#   
#   # check arrays and calculations
#   if (checks) {
#     if(!all(round(apply(initdist, c(1,2), sum),1)==1)) warning('initdist does not sum to one')
#     if(!(all(round(apply(Neq[,,1,], 1, sum) /  apply(Neq1, 1, sum),1) ==1))) warning('eq age structure')
#     sim <- sample(1:nsim,1)
#     yrval <- sample(1:Nyrs,1)
#     if (!all(M_ageArrayp[sim,,yrval] == M_ageArray[sim,,1] )) warning('problem with M_ageArrayp')
#     if(!all(Wt_agep[sim,,yrval] == Wt_age[sim,,1]))  warning('problem with Wt_agep')
#     if(!all(Mat_agep[sim,,yrval] == Mat_age[sim,,1])) warning('problem with Mat_agep')
#     
#   } 
# }


# ---- Unfished Equilibrium calcs ----
# arrays for unfished biomass for all years 
N_a <- array(NA, dim = c(nsim, maxage, histnTS+projnTS, nareas))
SSN_a <- array(NA, dim = c(nsim, maxage, histnTS+projnTS, nareas))  
Biomass_a <- array(NA, dim = c(nsim, maxage, histnTS+projnTS, nareas))
SSB_a <- array(NA, dim = c(nsim, maxage, histnTS+projnTS, nareas))

calcInitRec <- function(x, R0, recVec, M_ageArray) {
  maxage <- dim(M_ageArray)[2]
  nts <- length(recVec)
  recTS <- rep(recVec, maxage)
  M_yr <- M_ageArray[x, 1:(nts-1), 1:nts]
  M_yr <- do.call("cbind", rep(list(M_yr), maxage))
  
  R0init <- sapply(1:nts, function(ts)
    R0[x] * recTS[nts+ts] * exp(sum(diag(M_yr[, ts:(ts+nts-1)])))
    )
  R0init
}

R0init <- sapply(1:nsim, calcInitRec, StockPars$R0, recVec, M_ageArray)
R0 <- do.call("rbind", replicate(nyears+proyears, R0init, simplify = FALSE))

calcEqYr <- function(x, totyears=nyears+proyears, R0init, M_ageArray) {
  maxage <- dim(M_ageArray)[2]
  nts <- length(recVec)
  recTS <- rep(recVec, maxage)
  
  N <- matrix(NA, maxage, totyears*nts)
  
  RecMat <- matrix(NA, nrow=maxage, ncol=maxage)
  if (nts>1) {
  
    RecMat[,1] <- rep(c(R0init[1,x], R0init[nts:2,x]), maxage/nts)
  } else {
    RecMat <- matrix(R0init[x], nrow=maxage, ncol=maxage)
  }
  
  for (ts in 2:maxage) {
    RecMat[,ts] <- c(RecMat[maxage,ts-1], RecMat[1:(maxage-1),ts-1])
  }
  
  for (Yr in 1:totyears) {
    tsind <- seq(from=Yr*nts-nts+1, to=(Yr*nts-nts)+nts, by=1)
    if (nts > 1) {
      N[1,tsind] <-  R0init[,x]
      Ms <- diag(M_ageArray[x, 1:(maxage-1), rev(tsind)])  
    } else{
      N[1,tsind] <-  R0init[x]
      Ms <- M_ageArray[x, 1:(maxage-1), rev(tsind)]
    }
    Ms <- rep(Ms, maxage/nts)
    N[2:maxage,tsind[1]] <- RecMat[2:maxage,1] *  exp(-cumsum(Ms[1:(maxage-1)]))
    
    if (nts>1) {
      for (ts in tsind[2:length(tsind)]) {
        N[2:maxage,ts] <- N[1:(maxage-1),ts-1] * exp(-M_ageArray[x, 1:(maxage-1), ts-1])
      }
    }
  }
  N
}


N0list <- lapply(1:nsim, calcEqYr, nyears+proyears, R0init, M_ageArray)
N02 <- aperm(array(unlist(N0list), dim=c(maxage , histnTS+projnTS, nsim)), c(3,1,2))

SAYR_a <- as.matrix(expand.grid(1:nareas, 1:(histnTS+projnTS), 1:maxage, 1:nsim)[4:1]) 
SY_a <- SAYR_a[,c(1,3)]
SAR_a <- SAYR_a[, c(1,2,4)]
SA_a <- SAYR_a[,c(1,2)]
SAY_a <- SAYR_a[,c(1,2,3)]
# Unfished Equilibrium
N_a[SAYR_a] <- N02[SAY_a] * initdist[SAR_a] # Calculate initial numbers for all years
SSN_a[SAYR_a] <- N_a[SAYR_a] * Mat_age[SAY_a]  # Calculate initial spawning stock numbers for all years
Biomass_a[SAYR_a] <- N_a[SAYR_a] * Wt_age[SAY_a]  # Calculate initial stock biomass
SSB_a[SAYR_a] <- SSN_a[SAYR_a] * Wt_age[SAY_a]    # Calculate spawning stock biomass


# Convert to annual reference points - calculated at the end of each year
ind <- seq(from=nts, by=nts, to=histnTS+projnTS)
N0_a <- apply(N_a[,,ind,], c(1,3), sum)  # unfished numbers at the end of each year
B0_a <- apply(Biomass_a[,,ind,], c(1,3), sum) # unfished biomass for each year   
SSN0_a <- apply(SSN_a[,,ind,], c(1,3), sum) # unfished spawning numbers for each year
SSB0_a <- apply(SSB_a[,,ind,], c(1,3), sum) # unfished spawning biomass for each year 
vb <- apply(Biomass_a, 1:3, sum) * V
VB0_a <- apply(vb[,,ind], c(1,3), sum) # unfished vulnerable biomass for each year 
SSB0a_a <- apply(SSB_a[,,ind,], c(1,3,4), sum) # Calculate unfished spawning stock biomass by area for each year

UnfishedByYear <- list(SSN0=SSN0_a, N0=N0_a, SSB0=SSB0_a, B0=B0_a, VB0=VB0_a)

if (quantile(ageM[,1],0.95) > nyears + proyears) {
  if(!silent) message('Note: number of historical year `nyears` + `proyears` is less than the highest age of maturity')
}

# ---- Unfished Reference Points ----
SSBpRa <- array(SSB0_a/matrix(apply(N_a[,1,1:nts,], 1, sum), 
                              nrow=nsim, ncol=nyears+proyears), 
                dim = c(nsim, nyears+proyears))

ageM2 <- ageM / nts
UnfishedRefs <- sapply(1:nsim, CalcUnfishedRefs, ageM=ageM2, N0_a=N0_a, SSN0_a=SSN0_a,
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
## DEVELOP ##############
#### FIX initD and above sections for nts > 1 ###########
Perr_y <- StockPars$Perr_y

#  Depletion in year 1 
initD <- SampCpars$initD # 
if (!is.null(initD)) { # initial depletion is not unfished
  if (!silent) message("Optimizing for user-specified depletion in first historical year")
  Perrmulti <- sapply(1:nsim, optDfunwrap, initD=initD, Nfrac=Nfrac[,,1], R0=R0,
                      Perr_y=Perr_y, surv=surv[,,1], Wt_age=Wt_age, SSB0=SSB0,
                      maxage=maxage)
  Perr_y[,1:maxage] <- Perr_y[, 1:maxage] * Perrmulti
}

# --- Non-equilibrium calcs ----
# Initial Population for initial time-step

SAYR <- as.matrix(expand.grid(1:nareas, 1, 1:maxage, 1:nsim)[4:1])  
SAY <- SAYR[,1:3]
Sa <- maxage-SAYR[,2] + 1 # This is the process error index for initial year

N[SAYR] <- N_a[SAYR] * Perr_y[Sa] # Calculate initial stock numbers
SSN[SAYR] <- N[SAYR] * Mat_age[SAY]  # Calculate initial spawning stock numbers
Biomass[SAYR] <- N[SAYR] * Wt_age[SAY]  # Calculate initial stock biomass
SSB[SAYR] <- SSN[SAYR] * Wt_age[SAY]    # Calculate spawning stock biomass
VBiomass[SAYR] <- Biomass[SAYR] * V[SAY]  # Calculate vunerable biomass


######### TO DO - CHECK initD ##########
if (checks && !is.null(initD)) { # check initial depletion 
  if (nts >1) {
    SBinit <- apply(SSB[,,1:nts,], c(1,3), sum)
    SBinit <- apply(SBinit, 1, mean) # average spawning biomass in first year
  } else {
    SBinit <- apply(SSB[,,1,], 1, sum) # spawning biomass in first year
  }
  plot(SBinit/SSB0, initD)
  if (!any(round(SBinit/SSB0, 2) == round(initD,2))) warning('problem with initial depletion')
}


# --- Historical Spatial closures ----
MPA <- matrix(1, histnTS+projnTS, ncol=nareas) # fraction open to fishing 
if (all(!is.na(OM@MPA)) && sum(OM@MPA) != 0) { # historical spatial closures have been specified
  yrindex <- OM@MPA[,1]
  if (max(yrindex)>nyears) stop("Invalid year index for spatial closures: must be <= nyears")
  if (min(yrindex)<1) stop("Invalid year index for spatial closures: must be > 1")
  if (ncol(OM@MPA)-1 != nareas) stop("OM@MPA must be nareas + 1")
  TSdf <- data.frame(TS=1:(histnTS+projnTS), Year=rep(1:(nyears + proyears), each=nts))
  for (xx in seq_along(yrindex)) {
    tsindex <- which(TSdf[,2]==yrindex[xx]) %>% min()
    MPA[tsindex:nrow(MPA),] <- matrix(OM@MPA[xx, 2:ncol(OM@MPA)], nrow=length(tsindex:nrow(MPA)),ncol=nareas, byrow = TRUE)
  }
}

# --- Optimize catchability (q) to fit depletion ---- 

# Unfished recruitment by area - INITDIST OF AGE 1.
SAYR <- as.matrix(expand.grid(1:nareas, 1:(histnTS+projnTS), 1, 1:nsim)[4:1])  
YS <- SAYR[,c(3,1)]
SYR <- SAYR[,c(1,3,4)]
SAR <- SAYR[,c(1,2,4)]
R0a <- array(NA, dim=c(nsim, histnTS+projnTS, nareas))
R0a[SYR] <- R0[YS] * initdist[SAR]

if(!silent) message("Optimizing for user-specified depletion in last historical year")
bounds <- c(0.0001, 15) # q bounds for optimizer
qs <- sapply(1:nsim, getq3, D, SSB0, nareas, maxage, N, pyears=histnTS, 
             M_ageArray, Mat_age, Asize, Wt_age, V, retA, Perr_y, mov, SRrel, Find, 
             Spat_targ, hs, R0a, SSBpR, aR, bR, bounds=bounds, maxF=maxF,
             MPA=MPA, nts=nts) # find the q that gives current stock depletion

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
    
    ResampStockPars <- SampleStockPars(SubOM(OM2, "Stock"), OM2@nsim, histnTS, projnTS, cpars=SampCpars2, msg=FALSE)  
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
    qs[probQ] <- sapply(probQ, getq3, D, SSB0, nareas, maxage, N, pyears=histnTS, 
                        M_ageArray, Mat_age, Asize, Wt_age, V, retA, Perr_y, mov, SRrel, Find, 
                        Spat_targ, hs, R0a, SSBpR, aR, bR, bounds=bounds, MPA=MPA, maxF=maxF,
                        nts=nts, recTS=recTS)
    
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
  popdynCPP(nareas, maxage, Ncurr=N[x,,1,], histnTS,  
            M_age=M_ageArray[x,,], Asize_c=Asize[x,], MatAge=Mat_age[x,,], WtAge=Wt_age[x,,],
            Vuln=V[x,,], Retc=retA[x,,], Prec=Perr_y[x,], movc=split.along.dim(mov[x,,,,],4), 
            SRrelc=SRrel[x], 
            Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,,], 
            SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], Qc=qs[x], Fapic=0, MPA=MPA, maxF=maxF, 
            control=1, SSB0c=SSB0[x]))

N <- aperm(array(as.numeric(unlist(histYrs[1,], use.names=FALSE)), dim=c(maxage, histnTS, nareas, nsim)), c(4,1,2,3))
Biomass <- aperm(array(as.numeric(unlist(histYrs[2,], use.names=FALSE)), dim=c(maxage, histnTS, nareas, nsim)), c(4,1,2,3))
SSN <- aperm(array(as.numeric(unlist(histYrs[3,], use.names=FALSE)), dim=c(maxage, histnTS, nareas, nsim)), c(4,1,2,3))
SSB <- aperm(array(as.numeric(unlist(histYrs[4,], use.names=FALSE)), dim=c(maxage, histnTS, nareas, nsim)), c(4,1,2,3))
VBiomass <- aperm(array(as.numeric(unlist(histYrs[5,], use.names=FALSE)), dim=c(maxage, histnTS, nareas, nsim)), c(4,1,2,3))
FM <- aperm(array(as.numeric(unlist(histYrs[6,], use.names=FALSE)), dim=c(maxage, histnTS, nareas, nsim)), c(4,1,2,3))
FMret <- aperm(array(as.numeric(unlist(histYrs[7,], use.names=FALSE)), dim=c(maxage, histnTS, nareas, nsim)), c(4,1,2,3))
Z <-aperm(array(as.numeric(unlist(histYrs[8,], use.names=FALSE)), dim=c(maxage, histnTS, nareas, nsim)), c(4,1,2,3))


Depletion <- apply(SSB[,,histnTS,],1,sum)/SSB0

# Check that depletion is correct
if (checks) {
  if (prod(round(D, 2)/ round(Depletion,2)) != 1) {
    print(cbind(round(D,2), round(Depletion,2)))
    warning("Possible problem in depletion calculations")
  } 
} 


######## TODO Check MSY calcs once forward projections and yield calcs are sorted #####
# --- Calculate MSY statistics for each year ----
MSY_y <- array(0, dim=c(nsim, nyears+proyears)) # store MSY for each sim and year
FMSY_y <- MSY_y # store FMSY for each sim, and year
SSBMSY_y <- MSY_y # store SSBMSY for each sim, and year
BMSY_y <- MSY_y # store BMSY for each sim, and year
VBMSY_y <- MSY_y # store VBMSY for each sim, and year 

if(!silent) message("Calculating MSY reference points for each year")
# Calculate MSY ref points for each year 

for (y in 1:(nyears+proyears)) {
  R0_annual <- apply(N_a[,nts,1:nts,], 1, sum)

  MSYrefsYr <- sapply(1:nsim, optMSY_eq, M_ageArray, Wt_age, Mat_age, V, 
                      StockPars$maxage, R0_annual, SRrel, hs, yr.ind=y, nts, tyears=nyears+proyears)
  MSY_y[,y] <- MSYrefsYr[1, ]
  FMSY_y[,y] <- MSYrefsYr[2,]
  SSBMSY_y[,y] <- MSYrefsYr[3,]
  BMSY_y[,y] <- MSYrefsYr[6,]
  VBMSY_y[,y] <- MSYrefsYr[7,] 
}


# --- MSY reference points ----
# MSY reference points are calculated by averaging the annual MSY refs (calculated
# above) over (nyears - average age-of-maturity):nyears - where nyears is the
# last historical year
MSYRefPoints <- sapply(1:nsim, CalcMSYRefs, MSY_y=MSY_y, FMSY_y=FMSY_y, 
                       SSBMSY_y=SSBMSY_y, BMSY_y=BMSY_y, VBMSY_y=VBMSY_y, 
                       ageM=ageM, nyears, nts)

MSY <- MSYRefPoints[1,] %>% unlist() # record the MSY results (Vulnerable)
FMSY <- MSYRefPoints[2,] %>% unlist()  # instantaneous FMSY (Vulnerable)
SSBMSY <- MSYRefPoints[3,] %>% unlist()  # Spawning Stock Biomass at MSY
BMSY <- MSYRefPoints[4,] %>% unlist() # total biomass at MSY
VBMSY <- MSYRefPoints[5,] %>% unlist() # Biomass at MSY (Vulnerable)
UMSY <- MSY/VBMSY  # exploitation rate 
FMSY_M <- FMSY/(M*nts)  # ratio of true FMSY to natural mortality rate M

SSBMSY_SSB0 <- SSBMSY/SSB0 # SSBMSY relative to unfished (SSB)
BMSY_B0 <- BMSY/B0 # Biomass relative to unfished (B0)
VBMSY_VB0 <- VBMSY/VB0 # VBiomass relative to unfished (VB0)

if (!AnnualMSY) {
  warning('AnnualMSY argument is deprecated. MSY metrics are always calculated by year.\n Use `MSE@SSB` or `MSE@B` and `MSE@Misc$MSYRefs$ByYear` for alternative methods to calculate B/BMSY')
}


if (checks) {
  Btemp <- apply(SSB, c(1,3), sum)
  x <- Btemp[,histnTS]/SSBMSY
  y <-D/SSBMSY_SSB0
  plot(x,y, xlim=c(0,max(x)), ylim=c(0,max(y)), xlab="SSB/SSBMSY", ylab="D/SSBMSY_SSB0")
  lines(c(-10,10),c(-10,10))
}


# --- Calculate B-low ---- 

MarrayArea <- replicate(nareas, M_ageArray[,,1:histnTS])
Mnow<-apply(MarrayArea[,,histnTS,]*N[,,histnTS,],1:2,sum)/apply(N[,,histnTS,],1:2,sum)
MGTsurv<-t(exp(-apply(Mnow,1,cumsum)))
MGT<-apply(Agearray*(Mat_age[,,histnTS]*MGTsurv),1,sum)/apply(Mat_age[,,histnTS]*MGTsurv,1,sum)

Blow <- rep(NA,nsim)
if(CalcBlow){
  if(!silent) message("Calculating B-low reference points")            
  MGThorizon<-floor(HZN*MGT)
  Blow <- sapply(1:nsim,getBlow, N, Asize, SSBMSY,SSBpR, MPA, SSB0, nareas, retA, MGThorizon,
                 Find,Perr_y,M_ageArray,hs,Mat_age, Wt_age,R0a,V,histnTS,maxage,mov,
                 Spat_targ,SRrel,aR,bR,Bfrac, maxF) 
}

# --- Calculate Reference Yield ----
if(!silent) message("Calculating reference yield - best fixed F strategy") 

RefY <- sapply(1:nsim, getFref3, Asize, nareas, maxage, N=N[,,histnTS-3,, drop=FALSE], pyears=projnTS, 
               M_ageArray=M_ageArray[,,(histnTS):(histnTS+projnTS)], Mat_age[,,(histnTS):(histnTS+projnTS)], 
               Wt_age=Wt_age[,,histnTS:(histnTS+projnTS)], 
               V=V[, , (histnTS + 1):(histnTS + projnTS), drop=FALSE], 
               retA=retA[, , (histnTS + 1):(histnTS + projnTS), drop=FALSE],  
               Perr=Perr_y[,(histnTS):(histnTS+maxage+projnTS-1)], mov, SRrel, Find, 
               Spat_targ, hs, R0a[,(histnTS+1):(histnTS+projnTS),], SSBpR, aR, bR, MPA=MPA, 
               maxF=maxF, SSB0=SSB0, nts=nts)

RefPoints <- data.frame(MSY=MSY, FMSY=FMSY, SSBMSY=SSBMSY, SSBMSY_SSB0=SSBMSY_SSB0,
                        BMSY_B0=BMSY_B0, BMSY=BMSY, VBMSY=VBMSY, UMSY=UMSY, VBMSY_VB0=VBMSY_VB0,
                        FMSY_M=FMSY_M, N0=N0, SSB0=SSB0, B0=B0, VB0=VB0, RefY=RefY, 
                        Blow=Blow, MGT=MGT, R0=R0[1,])

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
ErrList$Cbiasa <- array(ObsPars$Cbias, c(nsim, histnTS + projnTS))  # Catch bias array
if (!is.null(control$Cbias_yr)) { # catch bias specified with control argument 
  Cbiasa <- matrix(1, nsim, histnTS+proyears)
  Cbiasa[,control$yrs] <- control$Cbias_yr
  ErrList$Cbiasa <- Cbiasa
} 
# composite of bias and observation error
ErrList$Cerr <- array(rlnorm((histnTS + projnTS) * nsim, 
                             mconv(1, rep(ObsPars$Csd, (histnTS + projnTS))), 
                             sdconv(1, rep(ObsPars$Csd, histnTS + projnTS))), 
                      c(nsim, histnTS + projnTS))  
# Index error
ErrList$Ierr <- array(rlnorm((histnTS + projnTS) * nsim, 
                             mconv(1, rep(Isd, histnTS + projnTS)), 
                             sdconv(1, rep(Isd, histnTS + projnTS))), 
                      c(nsim, histnTS + projnTS))

# Simulate error in observed recruitment index 
ErrList$Recerr <- array(rlnorm((histnTS + projnTS) * nsim, mconv(1, rep(Recsd, (histnTS + projnTS))), 
                               sdconv(1, rep(Recsd, histnTS + projnTS))), c(nsim, histnTS + projnTS))

# --- Implementation error time series ----
TAC_f <- array(rlnorm(projnTS * nsim, mconv(TACFrac, TACSD),
                      sdconv(TACFrac, TACSD)), c(nsim, projnTS))  # composite of TAC fraction and error
E_f <- array(rlnorm(projnTS * nsim, mconv(TAEFrac, TAESD),
                    sdconv(TAEFrac, TAESD)), c(nsim, projnTS))  # composite of TAE fraction and error
SizeLim_f<-array(rlnorm(projnTS * nsim, mconv(SizeLimFrac, SizeLimSD),
                        sdconv(SizeLimFrac, SizeLimSD)), c(nsim, projnTS))  # composite of size limit fraction and error


# --- Populate Data object with Historical Data ---- 
Data <- makeData(Biomass, CBret, Cret, N, SSB, VBiomass, StockPars, 
                 FleetPars, ObsPars, ImpPars, RefPoints,
                 ErrList, OM, SampCpars, initD, control=control,
                 silent=silent)

# --- Add Real Indices to Data object (if they exist) & calculate error stats ----
templist <- addRealInd(Data, nts, SampCpars, ErrList, Biomass, VBiomass, SSB, nsim,
                       nyears, proyears, silent=silent)
Data <- templist$Data # update 
ErrList <- templist$ErrList # update
Misc$RInd.stats <- ErrList$stats.df # return stats
Misc$nts <- nts
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
FMats <- array(NA, dim = c(nsim, nMP, projnTS))  # store the projected fishing mortality rate in each time-step
Ca <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected removed catch
CaRet <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected retained catch
TACa <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected TAC recommendation
Effort <- array(NA, dim = c(nsim, nMP, proyears))  # store the Effort
PAAout <- array(NA, dim = c(nsim, nMP, maxage))  # store the population-at-age in last projection year
CAAout <- array(NA, dim = c(nsim, nMP, maxage))  # store the catch-at-age in last projection year
CALout <- array(NA, dim = c(nsim, nMP, nCALbins))  # store the population-at-length in last projection year
# SPRa <- array(NA,dim=c(nsim,nMP,proyears)) # store the Spawning Potential Ratio

Cost_out <- array(NA, dim = c(nsim, nMP, proyears))  # store Total Cost
Rev_out <- array(NA, dim = c(nsim, nMP, proyears))  # store Total Revenue
LatEffort_out<- array(NA, dim = c(nsim, nMP, proyears))  # store the Latent Effort
TAE_out <- array(NA, dim = c(nsim, nMP, proyears)) # store the TAE

# --- Begin loop over MPs ----
mm <- 1 # for debugging

for (mm in 1:nMP) {  # MSE Loop over methods
  tryMP <- tryCatch({
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
    LatentEff_MP <- LatentEff # Historical latent effort
    
    # projection arrays for each sub-year time-step
    N_P <- array(NA, dim = c(nsim, maxage, projnTS, nareas))
    Biomass_P <- array(NA, dim = c(nsim, maxage, projnTS, nareas))
    VBiomass_P <- array(NA, dim = c(nsim, maxage, projnTS, nareas))
    SSN_P <-array(NA, dim = c(nsim, maxage, projnTS, nareas))
    SSB_P <- array(NA, dim = c(nsim, maxage, projnTS, nareas))
    FM_P <- array(NA, dim = c(nsim, maxage, projnTS, nareas))
    FM_Pret <- array(NA, dim = c(nsim, maxage, projnTS, nareas)) # retained F 
    Z_P <- array(NA, dim = c(nsim, maxage, projnTS, nareas))
    CB_P <- array(NA, dim = c(nsim, maxage, projnTS, nareas))
    CB_Pret <- array(NA, dim = c(nsim, maxage, projnTS, nareas)) # retained catch 
    
    # indexes
    SAYRL <- as.matrix(expand.grid(1:nsim, 1:maxage, histnTS, 1:nareas))  # Final historical year
    SAYRt <- as.matrix(expand.grid(1:nsim, 1:maxage, 1 + histnTS, 1:nareas))  # First projection year
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
    
    y <- 1; ts <- 1
    if(!silent) {
      cat("."); flush.console()
    }
    
    # Recruitment and movement in first year 
    NextYrN <- lapply(1:nsim, function(x)
      popdynOneTScpp(nareas, maxage, SSBcurr=colSums(SSB[x,,histnTS, ]), Ncurr=N[x,,histnTS,],
                     Zcurr=Z[x,,histnTS,], PerrYr=Perr_y[x, histnTS+maxage-1], hs=hs[x],
                     R0a=R0a[x,(histnTS+1),], SSBpR=SSBpR[x,], aR=aR[x,], bR=bR[x,],
                     mov=mov[x,,,,histnTS+1], SRrel=SRrel[x]))
    
    # The stock at the beginning of projection period
    N_P[,,1,] <- aperm(array(unlist(NextYrN), dim=c(maxage, nareas, nsim, 1)), c(3,1,4,2))
    Biomass_P[SAYR] <- N_P[SAYR] * Wt_age[SAY1]  # Calculate biomass
    VBiomass_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt]  # Calculate vulnerable biomass
    SSN_P[SAYR] <- N_P[SAYR] * Mat_age[SAY1]  # Calculate spawning stock numbers
    SSB_P[SAYR] <- SSN_P[SAYR] * Wt_age[SAY1]
    
    # Update abundance estimates - used for FMSY ref methods so that FMSY is applied to current abundance
    
    ########## TO DO - abundance should be mid-year abundance - not after first sub-year time-step ##########
    
    M_array <- array(0.5*M_ageArray[,,histnTS+y], dim=c(nsim, maxage, nareas))
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
    LastSpatial <- array(MPA[histnTS,], dim=c(nareas, nsim)) # 
    LastAllocat <- rep(1, nsim) # default assumption of reallocation of effort to open areas
    LastTAC <- LastCatch <- apply(CB[,,(histnTS-nts+1):histnTS,], 1, sum) # last years catch
    
    # -- Bio-Economics ----
    # Calculate Profit from last historical year
    RevPC <- RevCurr/LastCatch # cost-per unit catch in last historical year
    PMargin <- 1 - CostCurr/(RevPC * LastCatch) # profit margin in last historical year
    Profit <- (RevPC * LastCatch) - CostCurr # profit in last historical year
    HistEffort <- rep(1, nsim) # future effort is relative to today's effort
    Effort_pot <- HistEffort + Response*Profit # potential effort in first projection year
    Effort_pot[Effort_pot<0] <- tiny # 
    
    # Latent Effort - Maximum Effort Limit
    if (!all(is.na(LatentEff_MP))) {
      LastTAE <- histTAE <- HistEffort / (1 - LatentEff_MP) # current TAE limit exists    
    } else {
      LastTAE <- histTAE <- rep(NA, nsim) # no current TAE exists  
    }
    
    # TO DO - add SeasonalEff to Rec object #####
    MPRecs$SeasonalEff <- matrix(NA, nrow=nts, ncol=nsim)
    SeasonalEff <- c(1,1, 1, 1)
    SeasonalEff <- SeasonalEff/sum(SeasonalEff)
    MPRecs$SeasonalEff <- matrix(SeasonalEff, nrow=nts, ncol=nsim)
    ##############################################################
    
    for (ts in 1:nts) {
      # loop over sub-year time-steps in first year
      MPCalcs <- CalcMPDynamics(MPRecs, y, ts, histnTS, projnTS, nsim, Biomass_P, VBiomass_P,
                                LastTAE, histTAE, LastSpatial, LastAllocat, LastTAC,
                                TACused, maxF,
                                LR5_P, LFR_P, Rmaxlen_P, retL_P, retA_P,
                                L5_P, LFS_P, Vmaxlen_P, SLarray_P, V_P,
                                Fdisc_P, DR_P,
                                M_ageArray, FM_P, FM_Pret, Z_P, CB_P, CB_Pret,
                                TAC_f, E_f, SizeLim_f,
                                FinF, Spat_targ,
                                CAL_binsmid, Linf, Len_age, maxage, nareas, Asize, nCALbins,
                                qs, qvar, qinc,
                                Effort_pot)
      
      CB_P <- MPCalcs$CB_P # removals
      CB_Pret <- MPCalcs$CB_Pret # retained catch 
      # apply(CB_Pret[,,1,], 1, sum)
      FM_P <- MPCalcs$FM_P # fishing mortality
      FM_Pret <- MPCalcs$FM_Pret # retained fishing mortality 
      Z_P <- MPCalcs$Z_P # total mortality
      retA_P <- MPCalcs$retA_P # retained-at-age
      retL_P <- MPCalcs$retL_P # retained-at-length
      V_P <- MPCalcs$V_P  # vulnerable-at-age
      SLarray_P <- MPCalcs$SLarray_P # vulnerable-at-length
      FMats[,mm,ts] <- MPCalcs$Ftot 
      
      # update population dynamics for next sub-year time-step (or next year if nts==1)
      SAYRt <- as.matrix(expand.grid(1:nsim, 1:maxage, ts + histnTS+1, 1:nareas))  # Trajectory year
      SAYt <- SAYRt[, 1:3]
      SAYtMP <- cbind(SAYt, mm)
      SYt <- SAYRt[, c(1, 3)]
      SAY1R <- as.matrix(expand.grid(1:nsim, 1:maxage, ts, 1:nareas))
      SAYR <- as.matrix(expand.grid(1:nsim, 1:maxage, ts+1, 1:nareas))
      SY <- SAYR[, c(1, 3)]
      SA <- SAYR[, 1:2]
      S1 <- SAYR[, 1]
      SAY <- SAYR[, 1:3]
      S <- SAYR[, 1]
      SR <- SAYR[, c(1, 4)]
      SA2YR <- as.matrix(expand.grid(1:nsim, 2:maxage, ts+1, 1:nareas))
      SA1YR <- as.matrix(expand.grid(1:nsim, 1:(maxage - 1), ts, 1:nareas))
      
      # --- Age & Growth ----
      NextYrN <- lapply(1:nsim, function(x)
        popdynOneTScpp(nareas, maxage, SSBcurr=colSums(SSB_P[x,,ts, ]), Ncurr=N_P[x,,ts,],
                       Zcurr=Z_P[x,,ts,], PerrYr=Perr_y[x, ts+histnTS+maxage], hs=hs[x],
                       R0a=R0a[x,histnTS+1+ts,], SSBpR=SSBpR[x,], aR=aR[x,], bR=bR[x,],
                       mov=mov[x,,,, histnTS+ts+1], SRrel=SRrel[x]))
      N_P[,,ts+1,] <- aperm(array(unlist(NextYrN), dim=c(maxage, nareas, nsim, 1)), c(3,1,4,2)) 
      Biomass_P[SAYR] <- N_P[SAYR] * Wt_age[SAYt]  # Calculate biomass
      VBiomass_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt]  # Calculate vulnerable biomass
      SSN_P[SAYR] <- N_P[SAYR] * Mat_age[SAYt]  # Calculate spawning stock numbers
      SSB_P[SAYR] <- SSN_P[SAYR] * Wt_age[SAYt]  # Calculate spawning stock biomass
    
    }
    
    TACa[, mm, y] <- MPCalcs$TACrec # recommended TAC 
    LastSpatial <- MPCalcs$Si
    LastAllocat <- MPCalcs$Ai
    LastTAE <- MPCalcs$TAE # TAE set by MP 
    LastTAC <- MPCalcs$TACrec # TAC set by MP
    
    Effort[, mm, y] <- MPCalcs$Effort  ## TO DO ####### Average effort?
    
    # ---- Bio-economics ----
    RetainCatch <- apply(CB_Pret[,,1:nts,], 1, sum) # retained catch this year
    RetainCatch[RetainCatch<=0] <- tiny
    Cost_out[,mm,y] <-  Effort[, mm, y] * CostCurr*(1+CostInc/100)^y # cost of effort this year
    Rev_out[,mm,y] <- (RevPC*(1+RevInc/100)^y * RetainCatch)
    PMargin <- 1 - Cost_out[,mm,y]/Rev_out[,mm,y] # profit margin this year
    Profit <- Rev_out[,mm,y] - Cost_out[,mm,y] # profit this year
    Effort_pot <- Effort_pot + Response*Profit # bio-economic effort next year
    Effort_pot[Effort_pot<0] <- tiny # 
    # LatEffort_out[,mm,y] <- LastTAE - Effort[, mm, y]  # store the Latent Effort
    TAE_out[,mm,y] <- LastTAE # store the TAE
    
    # --- Begin projections ----
    for (ts in (nts+1):projnTS) {
      
      if (ts %% nts == 1){
        y <- y +1
        if(!silent) {
          cat("."); flush.console()  # update message every year
        }  
      }
      
      SelectChanged <- FALSE
      if (AnnualMSY) {
        if (any(range(retA_P[,,histnTS+ts] - retA[,,histnTS+ts]) !=0)) SelectChanged <- TRUE
        if (any(range(V_P[,,histnTS+ts] - V[,,histnTS+ts]) !=0))  SelectChanged <- TRUE
      }
      
      # -- Calculate MSY stats for this year ----
      # if (AnnualMSY & SelectChanged) { #
        # y1 <- nyears + y
        # MSYrefsYr <- sapply(1:nsim, optMSY_eq, M_ageArray, Wt_age, Mat_age, 
        #                     V_P, maxage, R0, SRrel, hs, yr.ind=y1)
        # MSY_y[,mm,y] <- MSYrefsYr[1, ]
        # FMSY_y[,mm,y] <- MSYrefsYr[2,]
        # SSBMSY_y[,mm,y] <- MSYrefsYr[3,]
      
        ## TO DO ##########
      # }    
      
      # --- Age & Growth ----
      SAYRt <- as.matrix(expand.grid(1:nsim, 1:maxage, ts + histnTS, 1:nareas))  # Trajectory year
      SAYt <- SAYRt[, 1:3]
      SAYtMP <- cbind(SAYt, mm)
      SYt <- SAYRt[, c(1, 3)]
      SAY1R <- as.matrix(expand.grid(1:nsim, 1:maxage, ts-1, 1:nareas))
      SAYR <- as.matrix(expand.grid(1:nsim, 1:maxage, ts, 1:nareas))
      SY <- SAYR[, c(1, 3)]
      SA <- SAYR[, 1:2]
      S1 <- SAYR[, 1]
      SAY <- SAYR[, 1:3]
      S <- SAYR[, 1]
      SR <- SAYR[, c(1, 4)]

      NextYrN <- lapply(1:nsim, function(x)
        popdynOneTScpp(nareas, maxage, SSBcurr=colSums(SSB_P[x,,ts-1, ]), Ncurr=N_P[x,,ts-1,],
                       Zcurr=Z_P[x,,ts-1,], PerrYr=Perr_y[x, ts+histnTS+maxage-1], hs=hs[x],
                       R0a=R0a[x,histnTS+ts,], SSBpR=SSBpR[x,], aR=aR[x,], bR=bR[x,],
                       mov=mov[x,,,, histnTS+ts], SRrel=SRrel[x]))
      N_P[,,ts,] <- aperm(array(unlist(NextYrN), dim=c(maxage, nareas, nsim, 1)), c(3,1,4,2)) 
      Biomass_P[SAYR] <- N_P[SAYR] * Wt_age[SAYt]  # Calculate biomass
      VBiomass_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt]  # Calculate vulnerable biomass
      SSN_P[SAYR] <- N_P[SAYR] * Mat_age[SAYt]  # Calculate spawning stock numbers
      SSB_P[SAYR] <- SSN_P[SAYR] * Wt_age[SAYt]  # Calculate spawning stock biomass
      
      
      TACa[, mm, y] <- TACa[, mm, y-1] # TAC same as last year unless changed 
      
      upts <- seq(1, by=nts, length.out=proyears) # update at the beginning of year
      # --- An update year ----
      if (y %in% upyrs & ts %in% upts) {
        # --- Update Data object ---- 
        MSElist[[mm]] <- updateData(Data=MSElist[[mm]], OM, MPCalcs, Effort, Biomass, 
                                    Biomass_P, CB_Pret, N_P, SSB, SSB_P, VBiomass, VBiomass_P, 
                                    RefPoints, ErrList, FMSY_y, retA_P, retL_P, StockPars, 
                                    FleetPars, ObsPars, upyrs, upts, interval, y, ts,
                                    mm, Misc=Data_p@Misc, SampCpars)
        
        # Update Abundance and FMSY for FMSYref MPs
        # TO DO #####
        # M_array <- array(0.5*M_ageArray[,,nyears+y], dim=c(nsim, maxage, nareas))
        # Atemp <- apply(VBiomass_P[, , y, ] * exp(-M_array), 1, sum) # Abundance (mid-year before fishing)
        # MSElist[[mm]]@OM$A <- Atemp
        # MSElist[[mm]]@OM$FMSY <- FMSY_y[,mm,y+OM@nyears]
        
        # --- apply MP ----
        runMP <- applyMP(Data=MSElist[[mm]], MPs = MPs[mm], reps = reps, silent=TRUE)  # Apply MP
        MPRecs <- runMP[[1]][[1]] # MP recommendations
        Data_p <- runMP[[2]] # Data object object with saved info from MP 
        Data_p@TAC <- MPRecs$TAC
        # calculate pstar quantile of TAC recommendation dist 
        TACused <- apply(Data_p@TAC, 2, quantile, p = pstar, na.rm = T)
        
        # TO DO - add SeasonalEff to Rec object #####
        MPRecs$SeasonalEff <- matrix(NA, nrow=nts, ncol=nsim)
        SeasonalEff <- c(1,1, 1, 1)
        SeasonalEff <- SeasonalEff/sum(SeasonalEff)
        MPRecs$SeasonalEff <- matrix(SeasonalEff, nrow=nts, ncol=nsim)
        
        # ---- Bio-economics ----
        # TO DO #####
        
        # -- Calc stock dynamics ----
        MPCalcs <- CalcMPDynamics(MPRecs, y, ts, histnTS, projnTS, nsim, Biomass_P, VBiomass_P,
                                  LastTAE, histTAE, LastSpatial, LastAllocat, LastTAC,
                                  TACused, maxF,
                                  LR5_P, LFR_P, Rmaxlen_P, retL_P, retA_P,
                                  L5_P, LFS_P, Vmaxlen_P, SLarray_P, V_P,
                                  Fdisc_P, DR_P,
                                  M_ageArray, FM_P, FM_Pret, Z_P, CB_P, CB_Pret,
                                  TAC_f, E_f, SizeLim_f,
                                  FinF, Spat_targ,
                                  CAL_binsmid, Linf, Len_age, maxage, nareas, Asize, nCALbins,
                                  qs, qvar, qinc,
                                  Effort_pot)
      
        CB_P <- MPCalcs$CB_P # removals
        CB_Pret <- MPCalcs$CB_Pret # retained catch 
        # apply(CB_Pret[,,1,], 1, sum)
        FM_P <- MPCalcs$FM_P # fishing mortality
        FM_Pret <- MPCalcs$FM_Pret # retained fishing mortality 
        Z_P <- MPCalcs$Z_P # total mortality
        retA_P <- MPCalcs$retA_P # retained-at-age
        retL_P <- MPCalcs$retL_P # retained-at-length
        V_P <- MPCalcs$V_P  # vulnerable-at-age
        SLarray_P <- MPCalcs$SLarray_P # vulnerable-at-length
        FMats[,mm,ts] <- MPCalcs$Ftot 
      } else {
        # --- Not an update yr ----
        NoMPRecs <- MPRecs # TAC & TAE stay the same
        NoMPRecs[lapply(NoMPRecs, length) > 0 ] <- NULL
        NoMPRecs$Spatial <- NA
        # TO DO ####### SEASONAL EFF
        NoMPRecs$SeasonalEff <- matrix(SeasonalEff, nrow=nts, ncol=nsim)
        
        MPCalcs <- CalcMPDynamics(NoMPRecs, y, ts, histnTS, projnTS, nsim, Biomass_P, VBiomass_P,
                                  LastTAE, histTAE, LastSpatial, LastAllocat, LastTAC,
                                  TACused, maxF,
                                  LR5_P, LFR_P, Rmaxlen_P, retL_P, retA_P,
                                  L5_P, LFS_P, Vmaxlen_P, SLarray_P, V_P,
                                  Fdisc_P, DR_P,
                                  M_ageArray, FM_P, FM_Pret, Z_P, CB_P, CB_Pret,
                                  TAC_f, E_f, SizeLim_f,
                                  FinF, Spat_targ,
                                  CAL_binsmid, Linf, Len_age, maxage, nareas, Asize, nCALbins,
                                  qs, qvar, qinc,
                                  Effort_pot)
        
        CB_P <- MPCalcs$CB_P # removals
        CB_Pret <- MPCalcs$CB_Pret # retained catch 
        # apply(CB_Pret[,,1,], 1, sum)
        FM_P <- MPCalcs$FM_P # fishing mortality
        FM_Pret <- MPCalcs$FM_Pret # retained fishing mortality 
        Z_P <- MPCalcs$Z_P # total mortality
        retA_P <- MPCalcs$retA_P # retained-at-age
        retL_P <- MPCalcs$retL_P # retained-at-length
        V_P <- MPCalcs$V_P  # vulnerable-at-age
        SLarray_P <- MPCalcs$SLarray_P # vulnerable-at-length
        FMats[,mm,ts] <- MPCalcs$Ftot 
        
      } # end of update loop 
      checkNA[y] <- sum(is.na(TACused))
      
      TACa[, mm, y] <- MPCalcs$TACrec # recommended TAC 
      LastSpatial <- MPCalcs$Si
      LastAllocat <- MPCalcs$Ai
      LastTAE <- MPCalcs$TAE # TAE set by MP 
      LastTAC <- MPCalcs$TACrec # TAC set by MP
      
    } # end of time-step projection
    
    tsind <- seq(nts, by=nts, to=projnTS) # end of year index
    B_BMSYa[, mm, ] <- apply(SSB_P[,,tsind,], c(1, 3), sum, na.rm=TRUE)/SSBMSY_y[,mm,(OM@nyears+1):(OM@nyears+OM@proyears)]  # SSB relative to SSBMSY
    F_FMSYa[, mm, ] <- FMa[, mm, ]/FMSY_y[,mm,(OM@nyears+1):(OM@nyears+OM@proyears)]
    
    Ba[, mm, ] <- apply(Biomass_P[,,tsind,], c(1, 3), sum, na.rm=TRUE) # biomass 
    SSBa[, mm, ] <- apply(SSB_P[,,tsind,], c(1, 3), sum, na.rm=TRUE) # spawning stock biomass
    VBa[, mm, ] <- apply(VBiomass_P[,,tsind,], c(1, 3), sum, na.rm=TRUE) # vulnerable biomass
    
    # annual catch 
    if (nts > 1) {
      Cobs <- apply(CB_P, c(1,3), sum)
      Cobs <- apply(Cobs, 1, function(x) 
        tapply(x, ceiling(seq_along(x)/nts), sum)) # sum up sub-year catches
      Cobs <- t(Cobs)
      Ca[, mm, ] <- Cobs # removed
      
      Cobs <- apply(CB_Pret, c(1,3), sum)
      Cobs <- apply(Cobs, 1, function(x) 
        tapply(x, ceiling(seq_along(x)/nts), sum)) # sum up sub-year catches
      Cobs <- t(Cobs)
      CaRet[, mm, ] <- Cobs # retained catch 
    } else {
      Ca[, mm, ] <- apply(CB_P, c(1, 3), sum, na.rm=TRUE) # removed
      CaRet[, mm, ] <- apply(CB_Pret, c(1, 3), sum, na.rm=TRUE) # retained catch 
    }

    # Store Pop and Catch-at-age and at-length for last projection year 
    # use Data instead
    PAAout[ , mm, ] <- NA #  apply(N_P[ , , projnTS, ], c(1,2), sum) # population-at-age
    # CNtemp <- apply(CB_Pret, c(1,2,3), sum)/Wt_age[,,(histnTS+1):(histnTS+projnTS)]
    CAAout[ , mm, ] <- NA # CNtemp[,,proyears] # nsim, maxage # catch-at-age
    # CALdat <- MSElist[[mm]]@CAL
    CALout[ , mm, ] <- NA # CALdat[,dim(CALdat)[2],] # catch-at-length in last year
    
    if (!silent) {
      cat("\n")
      if (all(checkNA != nsim) & !all(checkNA == 0)) {
        ntot <- sum(checkNA[upyrs])
        totyrs <- sum(checkNA[upyrs] >0)
        nfrac <- round(ntot/(length(upyrs)*nsim),2)*100
        message(totyrs, ' years had TAC = NA for some simulations (', nfrac, "% of total simulations)")
        message('Used TAC_y = TAC_y-1')  
      }
    }
    
    if (!parallel) 
      if("progress"%in%names(control))
        if(control$progress) 
          shiny::incProgress(1/nMP, detail = round(mm*100/nMP))
    
    
  }) # end tryCatch  # TO DO - fix tryCatch with new version ######
} # end loop over MP

# Miscellaneous reporting
if(PPD) Misc$Data <- MSElist

## Create MSE Object #### 
MSEout <- new("MSE", Name = OM@Name, nyears, proyears, nMPs=nMP, MPs, nsim, 
              Data@OM, Obs=Data@Obs, B_BMSY=B_BMSYa, F_FMSY=F_FMSYa, B=Ba, 
              SSB=SSBa, VB=VBa, FM=FMa, CaRet, TAC=TACa, SSB_hist = SSB, CB_hist = CB, 
              FM_hist = FM, Effort = Effort, PAA=PAAout, CAA=CAAout, CAL=CALout, CALbins=CAL_binsmid,
              Misc = Misc)
# Store MSE info
attr(MSEout, "version") <- packageVersion("DLMtool")
attr(MSEout, "date") <- date()
attr(MSEout, "R.version") <- R.version

 


