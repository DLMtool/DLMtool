
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
OM_annual@LFS <- c(39.347, 39.347)
OM_annual@Vmaxlen <- c(1,1)

OM_annual@D <- c(0.2, 0.5)
OM_annual@maxage <- 14
OM_annual@h <- c(0.9999, 0.999)

OM_annual@nsim <- 10


OM_annual <- tinyErr(OM_annual)

OM <- OMts(OM_annual)

OM@Perr <- c(0.0, 0.0)
OM@AC <- c(0,0)
OM@Linfsd <- c(0,0)
OM@Ksd <- c(0,0)


OM@Msd <- c(0.1, 0.1)

checks = TRUE

recVec <- c(0.1, 0.2, 0.3, 0.4)
recVec <- c(0.8, 0.05, 0.025, 0.125)

sum(recVec)

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
} else {
  OM@cpars$nts <- 1
  nts <- 1
}

# ---- Time-step Parameters ----
histnTS <- nyears * nts # number of time-steps in historical period
projnTS <- proyears * nts # number of time-steps in projection period

# --------------------------------------- #
 # add check to sum to one etc
recVec <- recVec/sum(recVec)
recTS <- rep(recVec, nyears+proyears) # recruitment per time-step



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

R0init <- sapply(1:nsim, calcInitRec, R0, recVec, M_ageArray)




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


N[,1:nts] %>% head()

N0list <- lapply(1:nsim, calcEqYr, nyears+proyears, R0init, M_ageArray)
N02 <- aperm(array(unlist(N0list), dim=c(maxage , histnTS+projnTS, nsim)), c(3,1,2))

N02[1,,1:nts] %>% head(20)


N02[1,,5:8] %>% head(20)


stop()




# First Time-Step
# Recruitment should follow recVec




RecMat[4,1] * exp(sum(-diag(M_ageArray[x,1:3,2:4])))




N_a[x, ,1:nts,1]


N_a[x,1:nts,1,1] <- R0init[x,1:nts]



RecMat[,1:nts]

surv <- array(1, dim=c(nsim, maxage, histnTS+projnTS)) # unfished survival for every year
surv[, 2:maxage, ] <- aperm(exp(-apply(M_ageArray, c(1,3), cumsum))[1:(maxage-1), ,], c(2,1,3)) # Survival array

RecMat[,1] * surv[x,,1]




# 
# R0Mat <- matrix(R0_ts, nrow=nsim, ncol=nts, byrow=TRUE) %>% apply(1, rev) %>% (t)
# # R0Mat <- replicate(maxage, R0Mat) 
# R0Mat <- do.call(cbind, replicate(maxage, R0Mat, simplify=FALSE))


# recVec <- c(1/3, 1/3, 1/3); nts = 3
# recVec <- c(0.5, 0.5); nts = 2

recVec <- c(0.1, 0.2, 0.3, 0.4); nts = 4


tt = calcInitRec(x=1, R0=R0, recVec=recVec, M_ageArray=M_ageArray)
tt[1,]
tt[2,]
            
# recVec <- 1 ; nts =1

R0_ts <- vapply(1:nsim, calcInitRec, R0=R0, recVec=recVec, M_ageArray=M_ageArray, 
                FUN.VALUE=matrix(0, nrow=histnTS+projnTS, ncol=maxage, byrow=TRUE))
R0_ts <- aperm(R0_ts, c(3,2,1))

surv <- array(1, dim=c(nsim, maxage, histnTS+projnTS)) # unfished survival for every year
surv[, 2:maxage, ] <- aperm(exp(-apply(M_ageArray, c(1,3), cumsum))[1:(maxage-1), ,], c(2,1,3)) # Survival array

SAYR_a <- as.matrix(expand.grid(1:nareas, 1:(histnTS+projnTS), 1:maxage, 1:nsim)[4:1]) 
SY_a <- SAYR_a[,c(1,3)]
SAR_a <- SAYR_a[, c(1,2,4)]
SA_a <- SAYR_a[,c(1,2)]
SAY_a <- SAYR_a[,c(1,2,3)]
N_a[SAYR_a] <- surv[SAY_a] * R0_ts[SAY_a] * initdist[SAR_a]



apply(N_a[1,1:5,1:10,], 1:2, sum)

surv[x,1:56,1] * R0_ts[x,1,]
surv[x,1:56,2] * R0_ts[x,2,]







R0Mat[1, ] * surv[1,,1]

SAYR_a <- as.matrix(expand.grid(1:nareas, 1, 1:maxage, 1:nsim)[4:1]) 
SY_a <- SAYR_a[,c(1,3)]
SAR_a <- SAYR_a[, c(1,2,4)]

N_a[SAYR_a] <- R0_ts[SY_a] * initdist[SAR_a]













SAYR_a <- as.matrix(expand.grid(1:nareas, 1, 1:(maxage-1), 1:nsim)[4:1]) 

SAYR3_a <- as.matrix(expand.grid(1:nareas, nts, 2:maxage, 1:nsim)[4:1])

#                                    Region    Rec Ages    Fst Lst Sim
SAY1Y2A1A2R <- as.matrix(expand.grid(1:nareas, 1, 2:maxage, 1, nts, 1:nsim)[6:1]) 

N_a[SAY1Y2A1A2R[,c(1,4,3,6)]] <- N_a[SAY1Y2A1A2R[, c(1, 5, 2, 6)]] * surv[SAY1Y2A1A2R[,c(1,4,3)]]



SArYR_a <- as.matrix(expand.grid(1:nareas, 1, nts, 2:maxage, 1:nsim)[5:1]) 
SArR_a <- SArYR_a[,c(1,2,3,5)]
SAYR_a <- SArYR_a[,c(1,2,4,5)]
SAY <- SArYR_a[,c(1,2,4)]

N_a[SAYR_a] <- N_a[SArR_a] * surv[SAY]




x <- 1
matrix(R0_ts[,x], nrow=nts, ncol=nts, byrow = TRUE) * matrix(surv[x,1:nts,1], nrow=nts, ncol=nts)

# First Year 
SAYR_a <- as.matrix(expand.grid(1:nareas, 1, 1:maxage, 1:nsim)[4:1])  
SAY_a <- SAYR_a[, 1:3]
SAR_a <- SAYR_a[, c(1,2,4)]
S_a <- SAYR_a[, 1]
Y_a <- SAYR_a[,3]
SA_a <- Sa_a <- SAYR_a[, 1:2]
SR_a <- SAYR_a[, c(1, 4)]
SY_a <- SAYR_a[, c(1, 3)]
N_a[SAYR_a]  <- R0_ts[S_a] * histRec[SAY_a] * initdist[SAR_a] * surv[SAY_a]


N_a[SAYR_a] <- R0[S_a] * histRec[SAY_a] * initdist[SAR_a] * surv[SAY_a]


N_a[,1,,]

x <- 1
recVec[nts]

R0[x] * recVec[nts] * exp(sum(diag(M_ageArray[x,,1:(nts-1)])))
R0[x] * recVec[1] * exp(sum(diag(M_ageArray[x,,2:(nts)])))
R0[x] * recVec[2] * exp(sum(diag(M_ageArray[x,,3:(nts)])))


tot_ts <- dim(M_ageArray)[3]
tempMat <- matrix(NA, ncol = length(recTS), nrow =maxage)
trecVec <- rev(recTS[1:maxage])
for (ts in 1:tot_ts) {
  tempMat[,ts] <- trecVec
  trecVec <- c(trecVec[nts], (trecVec[1:(nts-1)]))
}

x <-1 
for (ts in 1:(nts-1)) {
  tempM <- M_ageArray[x,ts:(nts-1),] 
  if (!is.null(nrow(tempM))) {
    tempM <- apply(tempM, 2, sum)
  }
  tempMat[ts,] <- tempMat[ts,] * exp(tempM)
}
tempMat[1:(nts-1),]

tempMat[,1:4]


recVec2 <- myRep(recVec, maxage) 



surv <- array(1, dim=c(nsim, maxage, histnTS+projnTS)) # unfished survival for every year
surv[, 2:maxage, ] <- aperm(exp(-apply(M_ageArray, c(1,3), cumsum))[1:(maxage-1), ,], c(2,1,3)) # Survival array

surv[1,,1:4]

# Recruitment to initialize population
# Calculate initial stock numbers for all years
histRec <- calcHistRec(recTS, maxage, nts, nsim, M_ageArray) # recruitment survival to each age
SAYR_a <- as.matrix(expand.grid(1:nareas, 1:(histnTS+projnTS), 1:maxage, 1:nsim)[4:1])  
SAY_a <- SAYR_a[, 1:3]
SAR_a <- SAYR_a[, c(1,2,4)]
S_a <- SAYR_a[, 1]
Y_a <- SAYR_a[,3]
SA_a <- Sa_a <- SAYR_a[, 1:2]
SR_a <- SAYR_a[, c(1, 4)]
SY_a <- SAYR_a[, c(1, 3)]
N_a[SAYR_a] <- R0[S_a] * histRec[SAY_a] * initdist[SAR_a] * surv[SAY_a]

Marray[1,1:4]
apply(N_a[1,,1:4,], 1:2, sum) %>% round(2)

apply(N_a[1,,5:8,], 1:2, sum) %>% round(2)

apply(N_a[1,,(histnTS+projnTS-nts+1):(histnTS+projnTS),], 1:2, sum) %>% round(2)

surv[1,,1:4]

#### STOP ####
stop()



# indices for all years
SAYR_a <- as.matrix(expand.grid(1:nareas, 1:(histnTS+projnTS), 1:maxage, 1:nsim)[4:1])  
SAY_a <- SAYR_a[, 1:3]
SAR_a <- SAYR_a[, c(1,2,4)]
SA_a <- SAYR_a[, 1:2]
SR_a <- SAYR_a[, c(1, 4)]
S_a <- SAYR_a[, 1]
SY_a <- SAYR_a[, c(1, 3)]

SSN_a[SAYR_a] <- N_a[SAYR_a] * Mat_age[SAY_a]  # Calculate initial spawning stock numbers for all years
Biomass_a[SAYR_a] <- N_a[SAYR_a] * Wt_age[SAY_a]  # Calculate initial stock biomass
SSB_a[SAYR_a] <- SSN_a[SAYR_a] * Wt_age[SAY_a]    # Calculate spawning stock biomass


# Convert to annual reference points - calculated at the end of each year
ind <- seq(from=nts, by=nts, to=maxage)
YrMat <- matrix(1:(histnTS+projnTS), ncol=nts,byrow=TRUE)

N0_a = sapply(1:(nyears+proyears), function(x)
  apply(N_a[,ind, YrMat[x,],], 1, sum))  # unfished numbers for each year
B0_a <- sapply(1:(nyears+proyears), function(x) 
  apply(Biomass_a[,ind, YrMat[x,],], 1, sum))  # unfished biomass for each year and sub-year age
SSN0_a <- sapply(1:(nyears+proyears), function(x) 
  apply(SSN_a[,ind, YrMat[x,],], 1, sum)) # unfished spawning numbers for each year
SSB0_a <- sapply(1:(nyears+proyears), function(x) 
  apply(SSB_a[,ind, YrMat[x,],], 1, sum)) # unfished spawning biomass for each year 
vb <- apply(Biomass_a, 1:3, sum) * V
VB0_a <- sapply(1:(nyears+proyears), function(x) 
  apply(vb[,ind, YrMat[x,]], 1, sum)) # unfished vulnerable biomass for each year 

SSB0a_a <- sapply(1:(nyears+proyears), function(x) 
  apply(SSB_a[,ind, YrMat[x,],], c(1,4), sum)) # Calculate unfished spawning stock biomass by area for each year

SSB0a_a <- array(SSB0a_a, dim=c(nsim, nyears+proyears, nareas)) 


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
# Initial Population for first year (nts = number of time-steps in the year)
# R0 is recruitment alive at the end of initial year - ie ts = nts

for (ts in 1:nts) {
  SAYR_a <- as.matrix(expand.grid(1:nareas, ts, 1:maxage, 1:nsim)[4:1])  
  SAY_a <- SAYR_a[, 1:3]
  SAR_a <- SAYR_a[, c(1,2,4)]
  S_a <- SAYR_a[, 1]
  Y_a <- SAYR_a[,3]
  SA_a <- Sa_a <- SAYR_a[, 1:2]
  SR_a <- SAYR_a[, c(1, 4)]
  SY_a <- SAYR_a[, c(1, 3)]
  Sa_a[,2] <- maxage-SA_a[,2] + 1 # This is the process error index for initial year
  N[SAYR_a] <- R0[S_a] * histRec[SAY_a] * initdist[SAR_a] * surv[SAY_a] * Perr_y[Sa_a]
}

SAYR <- as.matrix(expand.grid(1:nareas, 1:nts, 1:maxage, 1:nsim)[4:1])  
SAY <- SAYR[, 1:3]
SAR <- SAYR[, c(1,2,4)]
SA <- Sa <- SAYR[, 1:2]
SR <- SAYR[, c(1, 4)]
S <- SAYR[, 1]
SY <- SAYR[, c(1, 3)]

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
R0a <- matrix(apply(N[,1,1:nts,], 1, sum), nrow=nsim, ncol=nareas, byrow=FALSE) * initdist[,1,] # Unfished annual recruitment

if(!silent) message("Optimizing for user-specified depletion in last historical year")
bounds <- c(0.0001, 15) # q bounds for optimizer
qs <- sapply(1:nsim, getq3, D, SSB0, nareas, maxage, N, pyears=histnTS, 
             M_ageArray, Mat_age, Asize, Wt_age, V, retA, Perr_y, mov, SRrel, Find, 
             Spat_targ, hs, R0a, SSBpR, aR, bR, bounds=bounds, maxF=maxF,
             MPA=MPA, nts=nts, recTS=recTS) # find the q that gives current stock depletion

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
            Effind=Find[x,],  Spat_targc=Spat_targ[x], hc=hs[x], R0c=R0a[x,], 
            SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], Qc=qs[x], Fapic=0, MPA=MPA, maxF=maxF, 
            control=1, SSB0c=SSB0[x], recTS=recTS))

N <- aperm(array(as.numeric(unlist(histYrs[1,], use.names=FALSE)), dim=c(maxage, histnTS, nareas, nsim)), c(4,1,2,3))
Biomass <- aperm(array(as.numeric(unlist(histYrs[2,], use.names=FALSE)), dim=c(maxage, histnTS, nareas, nsim)), c(4,1,2,3))
SSN <- aperm(array(as.numeric(unlist(histYrs[3,], use.names=FALSE)), dim=c(maxage, histnTS, nareas, nsim)), c(4,1,2,3))
SSB <- aperm(array(as.numeric(unlist(histYrs[4,], use.names=FALSE)), dim=c(maxage, histnTS, nareas, nsim)), c(4,1,2,3))
VBiomass <- aperm(array(as.numeric(unlist(histYrs[5,], use.names=FALSE)), dim=c(maxage, histnTS, nareas, nsim)), c(4,1,2,3))
FM <- aperm(array(as.numeric(unlist(histYrs[6,], use.names=FALSE)), dim=c(maxage, histnTS, nareas, nsim)), c(4,1,2,3))
FMret <- aperm(array(as.numeric(unlist(histYrs[7,], use.names=FALSE)), dim=c(maxage, histnTS, nareas, nsim)), c(4,1,2,3))
Z <-aperm(array(as.numeric(unlist(histYrs[8,], use.names=FALSE)), dim=c(maxage, histnTS, nareas, nsim)), c(4,1,2,3))


if (nts >1) {
  ind <- seq(nts, by=nts, to = maxage)
  ssb_temp <- apply(SSB[,ind,(histnTS-nts+1):histnTS,],1,sum) # SSB in final year
  Depletion <- ssb_temp/SSB0
} else {
  Depletion <- apply(SSB[,,nyears,],1,sum)/SSB0
}
# Check that depletion is correct
if (checks) {
  if (prod(round(D, 2)/ round(Depletion,2)) != 1) {
    print(cbind(round(D,2), round(Depletion,2)))
    warning("Possible problem in depletion calculations")
  } 
} 


######## Check MSY calcs once forward projections and yield calcs are sorted #####
# --- Calculate MSY statistics for each year ----
MSY_y <- array(0, dim=c(nsim, nyears+proyears)) # store MSY for each sim and year
FMSY_y <- MSY_y # store FMSY for each sim, and year
SSBMSY_y <- MSY_y # store SSBMSY for each sim, and year
BMSY_y <- MSY_y # store BMSY for each sim, and year
VBMSY_y <- MSY_y # store VBMSY for each sim, and year 

if(!silent) message("Calculating MSY reference points for each year")
# Calculate MSY ref points for each year 
for (y in 1:(nyears+proyears)) {
  MSYrefsYr <- sapply(1:nsim, optMSY_eq, M_ageArray, Wt_age, Mat_age, V, 
                      maxage, R0, SRrel, hs, yr.ind=y, nts, tyears=nyears+proyears)
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




# --- Calculate MSY statistics for each time-step ----
MSY_ts <- array(0, dim=c(nsim, histnTS+projnTS)) # store MSY for each sim and time-step
FMSY_ts <- MSY_ts # store FMSY for each sim, and time-step
SSBMSY_ts <- MSY_ts # store SSBMSY for each sim, and time-step
BMSY_ts <- MSY_ts # store BMSY for each sim, and time-step
VBMSY_ts <- MSY_ts # store VBMSY for each sim, and time-step 

if(!silent) message("Calculating MSY reference points for each time-step")
 
for (ts in 1:(histnTS+projnTS)) {
  MSYrefsTS <- sapply(1:nsim, optMSY_eq_ts, M_ageArray, Wt_age, Mat_age, V, 
                      maxage, R0, SRrel, hs, ts.ind=ts, recTS)
  MSY_ts[,ts] <- MSYrefsTS[1, ]
  FMSY_ts[,ts] <- MSYrefsTS[2,]
  SSBMSY_ts[,ts] <- MSYrefsTS[3,]
  BMSY_ts[,ts] <- MSYrefsTS[6,]
  VBMSY_ts[,ts] <- MSYrefsTS[7,] 
}

# --- MSY reference points ----
# MSY reference points are calculated by averaging the annual MSY refs (calculated
# above) over (nyears - average age-of-maturity):nyears - where nyears is the
# last historical year

MSYRefPoints_TS <- sapply(1:nsim, CalcMSYRefs_TS, MSY_y=MSY_ts, FMSY_y=FMSY_ts, 
                       SSBMSY_y=SSBMSY_ts, BMSY_y=BMSY_ts, VBMSY_y=VBMSY_ts, 
                       ageM=ageM, histnTS, nts)

MSY_2 <- MSYRefPoints_TS[1,] %>% unlist() # record the MSY results (Vulnerable)
FMSY_2 <- MSYRefPoints_TS[2,] %>% unlist()  # instantaneous FMSY (Vulnerable)
SSBMSY_2 <- MSYRefPoints_TS[3,] %>% unlist()  # Spawning Stock Biomass at MSY
BMSY_2 <- MSYRefPoints_TS[4,] %>% unlist() # total biomass at MSY
VBMSY_2 <- MSYRefPoints_TS[5,] %>% unlist() # Biomass at MSY (Vulnerable)
UMSY_2 <- MSY_2/VBMSY_2  # exploitation rate 
FMSY_M_2 <- FMSY_2/(M)  # ratio of true FMSY to natural mortality rate M

SSBMSY_SSB0_2 <- SSBMSY_2/SSB0 # SSBMSY relative to unfished (SSB)
BMSY_B0_2 <- BMSY_2/B0 # Biomass relative to unfished (B0)
VBMSY_VB0_2 <- VBMSY_2/VB0 # VBiomass relative to unfished (VB0)

if (!AnnualMSY) {
  warning('AnnualMSY argument is deprecated. MSY metrics are always calculated by year.\n Use `MSE@SSB` or `MSE@B` and `MSE@Misc$MSYRefs$ByYear` for alternative methods to calculate B/BMSY')
}

print(MSY_2/MSY)
print(BMSY_2/BMSY)


if (checks) {
  Btemp <- apply(SSB, c(1,3), sum)
  x <- Btemp[,histnTS]/SSBMSY
  y <-D/SSBMSY_SSB0
  plot(x,y, xlim=c(0,max(x)), ylim=c(0,max(y)), xlab="SSB/SSBMSY", ylab="D/SSBMSY_SSB0")
  lines(c(-10,10),c(-10,10))
}

# do MSY calcs need to be calculated by time-step?

# why aren't SSB/SSBMSY and D/SSBMSY_SSBO the same with process error?
# initialized pop?


# --- Calculate B-low ---- 

## TO DO ####



# --- Calculate Reference Yield ----
if(!silent) message("Calculating reference yield - best fixed F strategy") 

RefY <- sapply(1:nsim, getFref3, Asize, nareas, maxage, N=N[,,histnTS,, drop=FALSE], pyears=projnTS, 
               M_ageArray=M_ageArray[,,(histnTS):(histnTS+projnTS)], Mat_age[,,(histnTS):(histnTS+projnTS)], 
               Wt_age=Wt_age[,,histnTS:(histnTS+projnTS)], 
               V=V[, , (histnTS + 1):(histnTS + projnTS), drop=FALSE], 
               retA=retA[, , (histnTS + 1):(histnTS + projnTS), drop=FALSE],  
               Perr=Perr_y[,(histnTS):(histnTS+maxage+projnTS-1)], mov, SRrel, Find, 
               Spat_targ, hs, R0a, SSBpR, aR, bR, MPA=MPA, maxF=maxF, SSB0=SSB0, recTS=recTS)
 ####### UP TO HERE @########

stop()


RefPoints <- data.frame(MSY=MSY, FMSY=FMSY, SSBMSY=SSBMSY, SSBMSY_SSB0=SSBMSY_SSB0,
                        BMSY_B0=BMSY_B0, BMSY=BMSY, VBMSY=VBMSY, UMSY=UMSY, VBMSY_VB0=VBMSY_VB0,
                        FMSY_M=FMSY_M, N0=N0, SSB0=SSB0, B0=B0, VB0=VB0, RefY=RefY, 
                        Blow=Blow, MGT=MGT, R0=R0)

Misc$MSYRefs <- list(Refs=RefPoints, ByYear=list(MSY=MSY_y, FMSY=FMSY_y,
                                                 SSBMSY=SSBMSY_y,
                                                 BMSY=BMSY_y,
                                                 VBMSY=VBMSY_y))





