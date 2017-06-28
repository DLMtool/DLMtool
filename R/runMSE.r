
Names <- c("maxage", "R0", "Mexp", "Msd", "dep", "Mgrad", "SRrel", "hs", "procsd",
           "L50", "L95", "L50_95", "CAL_binsmid", "Len_age", "maxlen", "Linf", 
           "M_at_Length", "Frac_area_1", "Prob_staying", "M_ageArray", "Mat_age",
           "Wt_age", "V", "Spat_targ", "procmu", "recMulti", "Linfrand", "Krand",
           "Abias Aerr", "Brefbias", "CAA_ESS", "CAA_nsamp", "CAL_ESS", "CAL_bins", "CAL_nsamp",
           "CALcv", "Cbias", "Crefbias", "Csd", "Dbias", "Derr", "DiscMort", "EFrac", "ESD", "EffLower",
           "EffUpper", "EffYears", "FMSY_Mbias", "Frac_area_1", "Irefbias", "Isd", "K", "Kbias", "Kgrad",
           "Krand", "Ksd", "L5", "L5s", "LFCbias", "LFS", "LFSbias", "LFSs", "LatASD", "Linfbias", "Linfgrad",
           "Linfrand", "Linfsd", "M", "M_ageArray", "Mat_age", "Mbias", "Mrand", "Prob_staying", "Recsd",
           "SLarray", "SizeLimFrac", "SizeLimSD", "Size_area_1", "Spat_targ", "TACFrac", "TACSD", 
           "Vmaxlen", "Vmaxlens", "Wt_age", "ageM", "betas", "lenMbias", "nCALbins", "procmu", "qcv", "qinc",
           "recMulti", "recgrad", "t0", "t0bias", "Abias", "Aerr", "Perr", "Esd", "qvar", "Marray",
           "Linfarray", "Karray", "AC", "LenCV", "LenCVbias", "a", "b", "FinF", "FecB")


if(getRversion() >= "2.15.1") utils::globalVariables(Names)

#' Run a Management Strategy Evaluation
#' 
#' A function that runs a Management Strategy Evaluation (closed-loop
#' simulation) for a specified operating model
#' 
#' 
#' @param OM An operating model object (class 'OM')
#' @param MPs A vector of methods (character string) of class Output or
#' Input.
#' @param nsim Number of simulations. Note that in DLMtool V4.1+ 'nsim is ignored 
#' if OM object contains the slot 'nsim'. 
#' @param proyears Number of projected years. Note that in DLMtool V4.1+ 'proyears is ignored 
#' if OM object contains the slot 'proyears'. 
#' @param interval The assessment interval - how often would you like to update
#' the management system?
#' @param pstar The percentile of the sample of the management recommendation
#' for each method
#' @param maxF Maximum instantaneous fishing mortality rate that may be
#' simulated for any given age class
#' @param timelimit Maximum time taken for a method to carry out 10 reps
#' (methods are ignored that take longer)
#' @param reps Number of samples of the management recommendation for each
#' method. Note that when this is set to 1, the mean value of the data inputs
#' is used.
#' @param CheckMPs Logical to indicate if Can function should be used to check
#' if MPs can be run.
#' @param Hist Should model stop after historical simulations? Returns a list 
#' containing all historical data
#' @param ntrials Maximum of times depletion and recruitment deviations are 
#' resampled to optimize for depletion. After this the model stops if more than 
#' percent of simulations are not close to the required depletion
#' @param fracD maximum allowed proportion of simulations where depletion is not 
#' close to sampled depletion from OM before model stops with error
#' @param CalcBlow Should low biomass be calculated where this is the spawning
#' biomass at which it takes HZN mean generation times of zero fishing to reach 
#' Bfrac fraction of SSBMSY
#' @param HZN The number of mean generation times required to reach Bfrac SSBMSY
#' in the Blow calculation
#' @param Bfrac The target fraction of SSBMSY for calculating Blow
#' @return An object of class MSE
#' @author T. Carruthers and A. Hordyk
#' @export 
runMSE <- function(OM = DLMtool::testOM, MPs = c("AvC","DCAC","FMSYref","curE","matlenlim"),nsim=48,
                   proyears=50,interval=4,pstar = 0.5, maxF = 0.8, timelimit = 1, reps = 1, 
                   CheckMPs = FALSE, Hist=FALSE, ntrials=50, fracD=0.05, CalcBlow=FALSE, 
                   HZN=2, Bfrac=0.5) {
  
  
   # For debugging - assign default argument values to to current workspace if they don't exist
  if (interactive()) { 
    DFargs <- formals(runMSE)
    argNames <- names(DFargs)
    for (X in seq_along(argNames)) {
      if (!exists(argNames[X])) {
        
        tt <- try(as.numeric(DFargs[X]), silent=TRUE)
        if (class(tt) != "try-error") {
          assign(argNames[X], tt)
        } else {
          if (argNames[X] == "OM") OM <- DLMtool::testOM
          if (argNames[X] == "MPs") MPs <- c("AvC","DCAC","FMSYref","curE","matlenlim")
        }
      }
    }
  }
  
  if (class(OM) != "OM") stop("You must specify an operating model")
  
  if("seed"%in%slotNames(OM)) set.seed(OM@seed) # set seed for reproducibility 
  
  OM <- ChkObj(OM) # Check that all required slots in OM object contain values 
  tiny <- 1e-15  # define tiny variable
  
  # Backwards compatible with DLMtool v < 4
  if("nsim"%in%slotNames(OM))nsim<-OM@nsim
  if("proyears"%in%slotNames(OM))proyears<-OM@proyears
  
  OM@nsim<-nsim # number of simulations
  OM@proyears<-proyears # number of projection years
  nyears <- OM@nyears  # number of historical years
  
  ### Sampling OM parameters ###
  message("Loading operating model")
  
  # --- Sample custom parameters ----
  SampCpars <- list() # empty list 
  # custom parameters exist - sample and write to list
  if(length(OM@cpars)>0){
    ncparsim<-cparscheck(OM@cpars)   # check each list object has the same length and if not stop and error report
    SampCpars <- SampleCpars(OM@cpars, nsim) 
  }
  
  # --- Sample Stock Parameters ----
  StockPars <- SampleStockPars(OM, nsim, nyears, proyears, SampCpars)
  # Assign Stock pars to function environment
  for (X in 1:length(StockPars)) assign(names(StockPars)[X], StockPars[[X]])

  # --- Sample Fleet Parameters ----
  FleetPars <- SampleFleetPars(SubOM(OM, "Fleet"), Stock=StockPars, nsim, nyears, proyears, SampCpars)
  # Assign Fleet pars to function environment
  for (X in 1:length(FleetPars)) assign(names(FleetPars)[X], FleetPars[[X]])
  
  # --- Sample Obs Parameters ----
  ObsPars <- SampleObsPars(OM, nsim)
  # Assign Obs pars to function environment
  for (X in 1:length(ObsPars)) assign(names(ObsPars)[X], ObsPars[[X]])

  # --- Sample Imp Paramerers ----
  ImpPars <- SampleImpPars(OM, nsim)
  # Assign Imp pars to function environment
  for (X in 1:length(ImpPars)) assign(names(ImpPars)[X], ImpPars[[X]])
  
  ### End of sampling OM parameters ###
  
  # --- Calculate movement ----
  message("Optimizing for user-specified movement")  # Print a progress update
  
  if (snowfall::sfIsRunning()) {
    # if the cluster is initiated
    snowfall::sfExport(list = c("Frac_area_1", "Prob_staying"))  # export some of the new arrays and ...
    mov <- array(t(snowfall::sfSapply(1:nsim, getmov2, Frac_area_1 = Frac_area_1, 
                                      Prob_staying = Prob_staying)), dim = c(nsim, 2, 2))  # numerically determine movement probability parameters to match Prob_staying and Frac_area_1
  } else {
    # no cluster initiated
    mov <- array(t(sapply(1:nsim, getmov2, Frac_area_1 = Frac_area_1, 
                          Prob_staying = Prob_staying)), dim = c(nsim, 2, 2))  # numerically determine movement probability parameters to match Prob_staying and Frac_area_1
  }
  
  nareas <- 2  # default is a two area model
  N <- array(NA, dim = c(nsim, maxage, nyears, nareas))  # stock numbers array
  Biomass <- array(NA, dim = c(nsim, maxage, nyears, nareas))  # stock biomass array
  VBiomass <- array(NA, dim = c(nsim, maxage, nyears, nareas))  # vulnerable biomass array
  
  SSN <- array(NA, dim = c(nsim, maxage, nyears, nareas))  # spawning stock numbers array
  
  SSB <- array(NA, dim = c(nsim, maxage, nyears, nareas))  # spawning stock biomass array
  FM <- array(NA, dim = c(nsim, maxage, nyears, nareas))  # fishing mortality rate array
  Z <- array(NA, dim = c(nsim, maxage, nyears, nareas))  # total mortality rate array
  SPR <- array(NA, dim = c(nsim, maxage, nyears)) # store the Spawning Potential Ratio
  
  Agearray <- array(rep(1:maxage, each = nsim), dim = c(nsim, maxage))  # Age array
  # surv <- exp(-Marray[, 1])^(Agearray - 1)  # Survival array
  
  # Survival array with M-at-age
  surv <- matrix(1, nsim, maxage)
  surv[, 2:maxage] <- t(exp(-apply(M_ageArray[,,1], 1, cumsum)))[, 1:(maxage-1)]
  
  Nfrac <- surv * Mat_age  # predicted Numbers of mature ages
  initdist <- as.matrix(cbind(Frac_area_1, 1 - Frac_area_1))  # Get the initial spatial distribution of each simulated population
  
  R0a <- matrix(R0, nrow=nsim, ncol=nareas, byrow=FALSE) * initdist  # Unfished recruitment by area
  
  SAYR <- as.matrix(expand.grid(1:nareas, 1, 1:maxage, 1:nsim)[4:1])  # Set up some array indexes sim (S) age (A) year (Y) region/area (R)
  SAY <- SAYR[, 1:3]
  SA <- Sa<-SAYR[, 1:2]
  SR <- SAYR[, c(1, 4)]
  S <- SAYR[, 1]
  SY <- SAYR[, c(1, 3)]
  Sa[,2]<-maxage-Sa[,2]+1 # This is the process error index for initial year
  

  #  --- Equilibrium calcs ----
  SSN[SAYR] <- Nfrac[SA] * R0[S] * initdist[SR]  # Calculate initial spawning stock numbers
  N[SAYR] <- R0[S] * surv[SA] * initdist[SR]  # Calculate initial stock numbers
  
  Biomass[SAYR] <- N[SAYR] * Wt_age[SAY]  # Calculate initial stock biomass
  SSB[SAYR] <- SSN[SAYR] * Wt_age[SAY]    # Calculate spawning stock biomass
  VBiomass[SAYR] <- Biomass[SAYR] * V[SAY]  # Calculate vunerable biomass
  
  if (nsim > 1) {
    SSN0 <- apply(SSN[, , 1, ], c(1, 3), sum)  # Calculate unfished spawning stock numbers  
    SSB0 <- apply(SSB[, , 1, ], 1, sum)  # Calculate unfished spawning stock biomass
    SSBpR <- SSB0/R0  # Spawning stock biomass per recruit
    SSBpR <- matrix(SSB0/R0, nrow=nsim, ncol=nareas)  # Spawning stock biomass per recruit
    SSB0a <- apply(SSB[, , 1, ], c(1, 3), sum)  # Calculate unfished spawning stock numbers
    B0 <- apply(Biomass[, , 1, ], 1, sum)
    N0 <- apply(N[, , 1, ], 1, sum)
  } else {
    SSN0 <- apply(SSN[, , 1, ], 2, sum)  # Calculate unfished spawning stock numbers  
    SSB0 <-  sum(SSB[, , 1, ])  # Calculate unfished spawning stock biomass
    SSBpR <- SSB0/R0  # Spawning stock biomass per recruit
    SSB0a <- apply(SSB[, , 1, ], 2, sum)  # Calculate unfished spawning stock numbers
    B0 <- apply(Biomass[, , 1, ], 2, sum)
    N0 <- apply(N[, , 1, ], 2, sum)
  }
  
  bR <- matrix(log(5 * hs)/(0.8 * SSB0a), nrow=nsim)  # Ricker SR params
  aR <- matrix(exp(bR * SSB0a)/SSBpR, nrow=nsim)  # Ricker SR params
  
  message("Optimizing for user-specified depletion")  # Print a progress update

  # --- Optimize catchability (q) to fit depletion ---- 
  # if (snowfall::sfIsRunning()) {
  #   snowfall::sfExport(list = c("dep", "Find", "Perr", "Marray", "hs", "Mat_age", 
  #     "Wt_age", "R0", "V", "nyears", "maxage", "SRrel", "aR", "bR"))
  #   qs <- snowfall::sfSapply(1:nsim, getq2, dep, Find, Perr, Marray, hs, Mat_age, 
  #     Wt_age, R0, V, nyears, maxage, mov, Spat_targ, SRrel, aR, bR)  # find the q that gives current stock depletion
  #   # qs <- snowfall::sfSapply(1:nsim, getq, dep, Find, Perr, Marray, hs, Mat_age, 
  #     # Wt_age, R0, V, nyears, maxage, mov, Spat_targ, SRrel, aR, bR)  # find the q that gives current stock depletion	  
  # } else {
  #   qs <- sapply(1:nsim, getq2, dep, Find, Perr, Marray, hs, Mat_age, 
  #     Wt_age, R0, V, nyears, maxage, mov, Spat_targ, SRrel, aR, bR)  # find the q that gives current stock depletion
  #   # qs <- sapply(1:nsim, getq, dep, Find, Perr, Marray, hs, Mat_age, 
  #     # Wt_age, R0, V, nyears, maxage, mov, Spat_targ, SRrel, aR, bR)  # find the q that gives current stock depletion	  
  # }
  bounds <- c(0.0001, 15) # q bounds for optimizer
  if (snowfall::sfIsRunning()) {
    snowfall::sfExport(list = c("dep", "Find", "Perr", "M_ageArray", "hs", "Mat_age", 
                                "Wt_age", "R0", "V", "nyears", "maxage", "SRrel", "aR", "bR"))
    qs <- snowfall::sfSapply(1:nsim, getq2, dep, Find, Perr, M_ageArray, hs, Mat_age, 
                             Wt_age, R0, V, nyears, maxage, mov, Spat_targ, SRrel, aR, bR, bounds)  # find the q that gives current stock depletion
  } else {
    qs <- sapply(1:nsim, getq2, dep, Find, Perr, M_ageArray, hs, Mat_age, 
                 Wt_age, R0, V, nyears, maxage, mov, Spat_targ, SRrel, aR, bR, bounds)  # find the q that gives current stock depletion
  }
  
  # --- Check that q optimizer has converged ---- 
  LimBound <- c(1.1, 0.9)*range(bounds)  # bounds for q (catchability). Flag if bounded optimizer hits the bounds 
  probQ <- which(qs > max(LimBound) | qs < min(LimBound))
  Nprob <- length(probQ)

  # If q has hit bound, re-sample depletion and try again. Tries 'ntrials' times
  # and then alerts user
  if (length(probQ) > 0) {
    Err <- TRUE
    message(Nprob,' simulations have final biomass that is not close to sampled depletion') 
    message('Re-sampling depletion, recruitment error, and fishing effort')
    
    count <- 0
    OM2 <- OM 
    while (Err & count < ntrials) {
      # Re-sample Stock Parameters 
      Nprob <- length(probQ)
      OM2@nsim <- Nprob
      SampCpars2 <- list()
      if (length(OM2@cpars)>0) SampCpars2 <- SampleCpars(OM2@cpars, OM2@nsim, msg=FALSE) 
      
      ResampStockPars <- SampleStockPars(OM2, cpars=SampCpars2)
      
      # Re-sample depletion 
      dep[probQ] <- ResampStockPars$dep 
  
      # Re-sample recruitment deviations
      procsd[probQ] <- ResampStockPars$procsd 
      AC[probQ] <- ResampStockPars$AC
      Perr[probQ,] <- ResampStockPars$Perr
      hs[probQ] <- ResampStockPars$hs

      # Re-sample historical fishing effort
      ResampFleetPars <- SampleFleetPars(SubOM(OM2, "Fleet"), Stock=ResampStockPars, 
                                         OM2@nsim, nyears, proyears, cpars=SampCpars2)
      Esd[probQ] <- ResampFleetPars$Esd
      Find[probQ, ] <- ResampFleetPars$Find
      dFfinal[probQ] <- ResampFleetPars$dFfinal
      
      # Optimize for q 
      if (snowfall::sfIsRunning()) {
        snowfall::sfExport(list = c("dep", "Find", "Perr", "M_ageArray", "hs", 
                                    "Mat_age", "Wt_age", "R0", "V", "nyears", "maxage", "SRrel", 
                                    "aR", "bR"))
        qs[probQ] <- snowfall::sfSapply(probQ, getq2, dep, Find, Perr, M_ageArray, 
                                        hs, Mat_age, Wt_age, R0, V, nyears, maxage, mov, Spat_targ, 
                                        SRrel, aR, bR)  # find the q that gives current stock depletion
      } else {
        qs[probQ] <- sapply(probQ, getq2, dep, Find, Perr, M_ageArray, 
                            hs, Mat_age, Wt_age, R0, V, nyears, maxage, mov, Spat_targ, 
                            SRrel, aR, bR)  # find the q that gives current stock depletion
      }
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
        message("More than ", fracD*100, "% of simulations can't get to the specified level of depletion with these Operating Model parameters")
        stop("Try again for a complete new sample, modify the input parameters, or increase ")
      } else {
        if (length(tooLow) > 0) message(tooLow, " sims can't get down to the lower bound on depletion")
        if (length(tooHigh) > 0) message(tooHigh, " sims can't get to the upper bound on depletion")
        message("Less than ", fracD*100, "% simulations can't get to the sampled depletion.\nContinuing")
      }
    }
  }
  
  #  --- Non-equilibrium calcs ----
  
  SSN[SAYR] <- Nfrac[SA] * R0[S] * initdist[SR]*Perr[Sa]  # Calculate initial spawning stock numbers
  N[SAYR] <- R0[S] * surv[SA] * initdist[SR]*Perr[Sa]  # Calculate initial stock numbers
  
  Biomass[SAYR] <- N[SAYR] * Wt_age[SAY]  # Calculate initial stock biomass
  SSB[SAYR] <- SSN[SAYR] * Wt_age[SAY]    # Calculate spawning stock biomass
  VBiomass[SAYR] <- Biomass[SAYR] * V[SAY]  # Calculate vunerable biomass
  
  message("Calculating historical stock and fishing dynamics")  # Print a progress update

  # Distribute fishing effort
  if (nsim > 1) fishdist <- (apply(VBiomass[, , 1, ], c(1, 3), sum)^Spat_targ)/
    apply(apply(VBiomass[, , 1, ], c(1, 3), sum)^Spat_targ, 1, mean)  # spatial preference according to spatial biomass
  if (nsim == 1)  fishdist <- (matrix(apply(VBiomass[,,1,], 2, sum), nrow=nsim)^Spat_targ)/
    mean((matrix(apply(VBiomass[,,1,], 2, sum), nrow=nsim)^Spat_targ))
  
  FM[SAYR] <- qs[S] * Find[SY] * V[SAY] * fishdist[SR]  # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
  # Z[SAYR] <- FM[SAYR] + Marray[SY]  # Total mortality rate                 
  Z[SAYR] <- FM[SAYR] + M_ageArray[SAY]  # Total mortality rate   
  
  # --- Simulate historical years ----
  for (y in 1:(nyears - 1)) {
    # set up some indices for indexed calculation
    SAYR <- as.matrix(expand.grid(1:nareas, y, 1:maxage, 1:nsim)[4:1])  # Set up some array indexes sim (S) age (A) year (Y) region/area (R)
    SAY1R <- as.matrix(expand.grid(1:nareas, y + 1, 1:maxage, 1:nsim)[4:1])
    SAY <- SAYR[, 1:3]
    SA <- SAYR[, 1:2]
    SR <- SAYR[, c(1, 4)]
    S <- SAYR[, 1]
    SY <- SAYR[, c(1, 3)]
    SY1 <- SAY1R[, c(1, 3)]
    indMov <- as.matrix(expand.grid(1:nareas, 1:nareas, y + 1, 1:maxage, 1:nsim)[5:1])  # Movement master index
    indMov2 <- indMov[, c(1, 2, 3, 4)]  # Movement from index
    indMov3 <- indMov[, c(1, 4, 5)]  # Movement to index
    
    if (nsim == 1) MAR <- 2 
    if (nsim >  1) MAR <- c(1, 3)
    if (SRrel[1] == 1) {
      N[, 1, y + 1, ] <- Perr[, y+maxage-1] * (0.8 * R0a * hs * 
                                                 apply(SSB[, , y, ], MAR, sum))/(0.2 * SSBpR * R0a * (1 - hs) + 
                                                                                   (hs - 0.2) * apply(SSB[, , y, ], MAR, sum))  # Recruitment assuming regional R0 and stock wide steepness
    } else {
      # most transparent form of the Ricker uses alpha and beta params
      N[, 1, y + 1, ] <- Perr[, y+maxage-1] * aR * apply(SSB[, , y, ], MAR, sum) *
        exp(-bR * apply(SSB[, , y, ], MAR, sum))
    }
    
    if (nsim > 1) fishdist <- (apply(VBiomass[, , y, ], c(1, 3), sum)^Spat_targ)/apply(apply(VBiomass[, , y, ], c(1, 3), sum)^Spat_targ, 1, mean)  # spatial preference according to spatial biomass
    if (nsim == 1)  fishdist <- (matrix(apply(VBiomass[,, y,], 2, sum), nrow=nsim)^Spat_targ)/mean((matrix(apply(VBiomass[,,y,], 2, sum), nrow=nsim)^Spat_targ))							   
    FM[SAY1R] <- qs[S] * Find[SY1] * V[SAY] * fishdist[SR]  # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
    # Z[SAY1R] <- FM[SAY1R] + Marray[SY]  # Total mortality rate
    Z[SAY1R] <- FM[SAY1R] + M_ageArray[SAY]  # Total mortality rate
    N[, 2:maxage, y + 1, ] <- N[, 1:(maxage - 1), y, ] * exp(-Z[, 1:(maxage - 1), y, ])  # Total mortality
    
    temp <- array(N[indMov2] * mov[indMov3], dim = c(nareas, nareas, maxage, nsim))  # Move individuals
    N[, , y + 1, ] <- apply(temp, c(4, 3, 1), sum)
    Biomass[SAY1R] <- N[SAY1R] * Wt_age[SAY]  # Calculate biomass
    VBiomass[SAY1R] <- Biomass[SAY1R] * V[SAY]  # Calculate vulnerable biomass
    SSN[SAY1R] <- N[SAY1R] * Mat_age[SA]  # Calculate spawning stock numbers
    SSB[SAY1R] <- SSN[SAY1R] * Wt_age[SAY]  # Calculate spawning stock biomass
    
  }  # end of year
  
  # Depletion <- apply(Biomass[, , nyears, ], 1, sum)/apply(Biomass[, , 1, ], 1, sum)  #^betas   # apply hyperstability / hyperdepletion
  if (nsim > 1) Depletion <- apply(SSB[,,nyears,],1,sum)/SSB0#^betas
  if (nsim == 1) Depletion <- sum(SSB[,,nyears,])/SSB0 #^betas
  # # apply hyperstability / hyperdepletion
  
  # Check that depletion is correct
  # print(cbind(round(dep,2), round(Depletion,2)))
  # if (prod(round(dep, 2)/ round(Depletion,2)) != 1) warning("Possible problem in depletion calculations")
  
  # --- Calculate MSY references ----  
  message("Calculating MSY reference points")  # Print a progress update
  
  # if (snowfall::sfIsRunning()) {
  #   snowfall::sfExport(list = c("Marray", "hs", "Mat_age", "Wt_age", "R0", "V", "nyears", "maxage"))  # export some newly made arrays to the cluster
  #   # MSYrefs <- snowfall::sfSapply(1:nsim, getFMSY, Marray, hs, Mat_age, Wt_age, 
  #     # R0, V = V[, , nyears], maxage, nyears, proyears = 200, Spat_targ, 
  #     # mov, SRrel, aR, bR)  # optimize for MSY reference points\t
  #   # Using Rcpp code 	
  #   MSYrefs <- snowfall::sfSapply(1:nsim, getFMSY2, Marray, hs, Mat_age, Wt_age, 
  #     R0, V = V, maxage, nyears, proyears = 200, Spat_targ, 
  #     mov, SRrel, aR, bR)  # optimize for MSY reference points\t	  
  # } else {
      # MSYrefs_R <- sapply(1:nsim, getFMSY, Marray, hs, Mat_age, Wt_age,
      #   R0, V = V[, , nyears], maxage, nyears, proyears = 200, Spat_targ,
      #   mov, SRrel, aR, bR)  # optimize for MSY reference points
  #   # Using Rcpp code 	  
    # MSYrefs <- sapply(1:nsim, getFMSY2, Marray, hs, Mat_age, Wt_age,
    #   R0, V = V, maxage, nyears, proyears = 200, Spat_targ,
    #   mov, SRrel, aR, bR)  # optimize for MSY reference points
  # }
  if (snowfall::sfIsRunning()) {
    snowfall::sfExport(list = c("M_ageArray", "hs", "Mat_age", "Wt_age", "R0", "V", "nyears", "maxage"))  # export some newly made arrays to the cluster
    MSYrefs <- snowfall::sfSapply(1:nsim, getFMSY2, M_ageArray, hs, Mat_age, Wt_age, 
                                  R0, V = V, maxage, nyears, proyears = 200, Spat_targ, 
                                  mov, SRrel, aR, bR)  # optimize for MSY reference points\t	  
  } else {
    MSYrefs <- sapply(1:nsim, getFMSY2, M_ageArray, hs, Mat_age, Wt_age, 
                      R0, V = V, maxage, nyears, proyears = 200, Spat_targ, 
                      mov, SRrel, aR, bR)  # optimize for MSY reference points
  }
  
  
  ## Commented out MSYrefs calculations  
  # MSY <- MSYrefs[1, ]  # record the MSY results (Vulnerable)
  # FMSY <- MSYrefs[2, ]  # instantaneous apical FMSY  (Vulnerable)
  # VBMSY <- (MSY/(1 - exp(-FMSY)))  # Biomass at MSY (Vulnerable)
  # UMSY <- MSY/VBMSY  # exploitation rate [equivalent to 1-exp(-FMSY)]
  # SSBMSY <- MSYrefs[3, ]  # Spawing Stock Biomass at MSY
  # BMSY_B0 <- SSBMSY_SSB0 <- MSYrefs[4, ]  # SSBMSY relative to unfished (SSB)
  # FMSYb <- -log(1-(MSY/(SSBMSY+MSY))) # instantaneous FMSY (Spawning Biomass)
  
  MSY <- MSYrefs[1, ]  # record the MSY results (Vulnerable)
  FMSY <- MSYrefs[2, ]  # instantaneous FMSY (Vulnerable)
  SSBMSY <- MSYrefs[3, ]  # Spawning Stock Biomass at MSY  
  SSBMSY_SSB0 <- SSBMSY/SSB0 # MSYrefs[4, ] # SSBMSY relative to unfished (SSB) 
  BMSY_B0 <- MSYrefs[6, ]/B0 # Biomass relative to unfished (B0)
  BMSY <- MSYrefs[6,] # total biomass at MSY
  
  VBMSY <- (MSY/(1 - exp(-FMSY)))  # Biomass at MSY (Vulnerable)
  # FMSYb <- MSYrefs[8,]  # instantaneous FMSY (Spawning Biomass)
  UMSY <- MSY/VBMSY  # exploitation rate [equivalent to 1-exp(-FMSY)]
  FMSY_M <- FMSY/M  # ratio of true FMSY to natural mortality rate M
  
  
  # --- Code for deriving low biomass ---- 
  # (SSB where it takes MGThorizon x MGT to reach Bfrac of BMSY)
  
  if(CalcBlow){
    message("Calculating Blow reference points")              # Print a progress update  
    
    Znow<-apply(Z[,,nyears,]*N[,,nyears,],1:2,sum)/apply(N[,,nyears,],1:2,sum)
    MGTsurv<-t(exp(-apply(Znow,1,cumsum)))
    MGT<-apply(Agearray*(Mat_age*MGTsurv),1,sum)/apply(Mat_age*MGTsurv,1,sum)
    MGThorizon<-floor(HZN*MGT)
    SSBMSY<-MSYrefs[3,]
    
    # if(snowfall::sfIsRunning()){
    #   snowfall::sfExport(list=c("SSBMSY","MGT","Find","Perr","Marray","hs","Mat_age","Wt_age","R0","V","nyears","maxage","SRrel","aR","bR"))
    #   Blow<-sfSapply(1:nsim,getBlow,SSBMSY,MGThorizon,Find,Perr,Marray,hs,Mat_age,Wt_age,R0,V,nyears,maxage,mov,Spat_targ,SRrel,aR,bR,Bfrac) # find the q that gives current stock depletion
    # }else{
    #   Blow <- sapply(1:nsim,getBlow,SSBMSY,MGThorizon,Find,Perr,Marray,hs,Mat_age,Wt_age,R0,V,nyears,maxage,mov,Spat_targ,SRrel,aR,bR,Bfrac) # find the q that gives current stock depletion
    # }
    if(snowfall::sfIsRunning()){
      snowfall::sfExport(list=c("SSBMSY","MGT","Find","Perr","M_ageArray","hs","Mat_age","Wt_age","R0","V","nyears","maxage","SRrel","aR","bR"))
      Blow<-sfSapply(1:nsim,getBlow,SSBMSY,MGThorizon,Find,Perr,M_ageArray,hs,Mat_age,Wt_age,R0,V,nyears,maxage,mov,Spat_targ,SRrel,aR,bR,Bfrac) # find the q that gives current stock depletion
    }else{
      Blow <- sapply(1:nsim,getBlow,SSBMSY,MGThorizon,Find,Perr,M_ageArray,hs,Mat_age,Wt_age,R0,V,nyears,maxage,mov,Spat_targ,SRrel,aR,bR,Bfrac) # find the q that gives current stock depletion
    }
  }else{
    Blow<-rep(NA,nsim)
  }
  
  # --- Calculate Reference Yield ----
  message("Calculating reference yield - best fixed F strategy")  # Print a progress update
  if (snowfall::sfIsRunning()) {
    RefY <- snowfall::sfSapply(1:nsim, getFref2, M_ageArray = M_ageArray, Wt_age = Wt_age, 
                               Mat_age = Mat_age, Perr = Perr, N_s = N[, , nyears, , drop=FALSE], SSN_s = SSN[, , nyears, , drop=FALSE], 
                               Biomass_s = Biomass[, , nyears, , drop=FALSE], VBiomass_s = VBiomass[, , nyears, , drop=FALSE], 
                               SSB_s = SSB[, , nyears, , drop=FALSE], Vn = V[, , (nyears + 1):(nyears + proyears), drop=FALSE], 
                               hs = hs, R0a = R0a, nyears = nyears, proyears = proyears, nareas = nareas,
                               maxage = maxage, mov = mov, SSBpR = SSBpR, aR = aR, bR = bR, SRrel = SRrel, Spat_targ = Spat_targ)
    
  } else {
    RefY <- sapply(1:nsim, getFref2, M_ageArray = M_ageArray, Wt_age = Wt_age, 
                   Mat_age = Mat_age, Perr = Perr, N_s = N[, , nyears, , drop=FALSE], SSN_s = SSN[, , nyears, , drop=FALSE], 
                   Biomass_s = Biomass[, , nyears, , drop=FALSE], VBiomass_s = VBiomass[, , nyears, , drop=FALSE], 
                   SSB_s = SSB[, , nyears, , drop=FALSE], Vn = V[, , (nyears + 1):(nyears + proyears), drop=FALSE], 
                   hs = hs, R0a = R0a, nyears = nyears, proyears = proyears, nareas = nareas,
                   maxage = maxage, mov = mov, SSBpR = SSBpR, aR = aR, bR = bR, SRrel = SRrel, Spat_targ = Spat_targ)
  }
  
  # --- Calculate catch-at-age ----
  CN <- apply(N * (1 - exp(-Z)) * (FM/Z), c(1, 3, 2), sum)  # Catch in numbers
  CN[is.na(CN)] <- 0
  CB <- Biomass * (1 - exp(-Z)) * (FM/Z)  # Catch in biomass
  
  # --- Simulate observed catch ---- 
  Cbiasa <- array(Cbias, c(nsim, nyears + proyears))  # Bias array
  Cerr <- array(rlnorm((nyears + proyears) * nsim, mconv(1, rep(Csd, (nyears + proyears))), 
                       sdconv(1, rep(Csd, nyears + proyears))), c(nsim, nyears + proyears))  # composite of bias and observation error
  Cobs <- Cbiasa[, 1:nyears] * Cerr[, 1:nyears] * apply(CB, c(1, 3), sum)  # Simulated observed catch (biomass)
  
  # --- Simulate observed catch-at-age ----
  CAA <- array(NA, dim = c(nsim, nyears, maxage))  # Catch  at age array
  cond <- apply(CN, 1:2, sum, na.rm = T) < 1  # this is a fix for low sample sizes. If CN is zero across the board a single fish is caught in age class of model selectivity (dumb I know)
  fixind <- as.matrix(cbind(expand.grid(1:nsim, 1:nyears), rep(floor(maxage/3), nyears)))  # more fix
  CN[fixind[cond, ]] <- 1  # puts a catch in the most vulnerable age class
  
  # a multinomial observation model for catch-at-age data
  for (i in 1:nsim) 
    for (j in 1:nyears) 
      CAA[i, j, ] <- ceiling(-0.5 + rmultinom(1, CAA_ESS[i], CN[i, j, ]) * CAA_nsamp[i]/CAA_ESS[i])  # a multinomial observation model for catch-at-age data
  
  
  # --- Simulate observed catch-at-length ----
  # a multinomial observation model for catch-at-length data
  # assumed normally-distributed length-at-age truncated at 2 standard deviations from the mean
  CAL <- array(NA, dim=c(nsim,  nyears, nCALbins))
  LFC <- rep(NA, nsim)
  vn <- (apply(N[,,,], c(1,2,3), sum) * V[,,1:nyears]) # vulnerable numbers at age
  vn <- aperm(vn, c(1,3, 2))
  
  for (i in 1:nsim) { # Rcpp code 
    # CAL[i, , ] <-  genLenComp(CAL_bins, CAL_binsmid, SLarray[i,,], CAL_ESS[i], CAL_nsamp[i], 
                              # CN[i,,], Len_age[i,,], LatASD[i,,], truncSD=2)
 
    CAL[i, , ] <-  genLenComp(CAL_bins, CAL_binsmid, SLarray[i,,], CAL_ESS[i], CAL_nsamp[i], 
                              vn[i,,], Len_age[i,,], LatASD[i,,], truncSD=2) 
    LFC[i] <- CAL_binsmid[min(which(round(CAL[i,nyears, ],0) >= 1))] # get the smallest CAL observation	  
  }
  
  # --- Simulate index of abundance from total biomass ----
  Ierr <- array(rlnorm((nyears + proyears) * nsim, mconv(1, rep(Isd, nyears + proyears)), 
                       sdconv(1, rep(Isd, nyears + proyears))), c(nsim, nyears + proyears))
  II <- (apply(Biomass, c(1, 3), sum) * Ierr[, 1:nyears])^betas  # apply hyperstability / hyperdepletion
  II <- II/apply(II, 1, mean)  # normalize
  
  # --- Calculate vulnerable and spawning biomass abundance ----
  if (nsim > 1) A <- apply(VBiomass[, , nyears, ], 1, sum)  # Abundance
  if (nsim == 1) A <- sum(VBiomass[, , nyears, ])  # Abundance
  if (nsim > 1) Asp <- apply(SSB[, , nyears, ], 1, sum)  # SSB Abundance
  if (nsim == 1) Asp <- sum(SSB[, , nyears, ])  # SSB Abundance  
  
  OFLreal <- A * FMSY  # the true simulated Over Fishing Limit
  

  # --- Simulate observed values in reference SBMSY/SB0 ----
  I3 <- apply(Biomass, c(1, 3), sum)^betas  # apply hyperstability / hyperdepletion
  I3 <- I3/apply(I3, 1, mean)  # normalize index to mean 1
  # Iref <- apply(I3[, 1:5], 1, mean) * BMSY_B0  # return the real target abundance index corresponding to BMSY
  if (nsim > 1) Iref <- apply(I3[, 1:5], 1, mean) * SSBMSY_SSB0  # return the real target abundance index corresponding to BMSY
  if (nsim == 1) Iref <- mean(I3[1:5]) * SSBMSY_SSB0

  # --- Simulate observed values in steepness ----
  hsim <- rep(NA, nsim)  
  cond <- hs > 0.6
  hsim[cond] <- 0.2 + rbeta(sum(hs > 0.6), alphaconv((hs[cond] - 0.2)/0.8, (1 - (hs[cond] - 0.2)/0.8) * OM@hcv), 
                            betaconv((hs[cond] - 0.2)/0.8,  (1 - (hs[cond] - 0.2)/0.8) * OM@hcv)) * 0.8
  hsim[!cond] <- 0.2 + rbeta(sum(hs < 0.6), alphaconv((hs[!cond] - 0.2)/0.8,  (hs[!cond] - 0.2)/0.8 * OM@hcv), 
                             betaconv((hs[!cond] - 0.2)/0.8, (hs[!cond] - 0.2)/0.8 * OM@hcv)) * 0.8
  hbias <- hsim/hs  # back calculate the simulated bias
  if (OM@hcv == 0) hbias <- rep(1, nsim) 
  ObsPars$hbias <- hbias
  
  # Simulate error in observed recruitment index 
  Recerr <- array(rlnorm((nyears + proyears) * nsim, mconv(1, rep(Recsd, (nyears + proyears))), 
                                sdconv(1, rep(Recsd, nyears + proyears))), c(nsim, nyears + proyears))
  
  
  # --- Simulate observation error in BMSY/B0 ---- 
  ntest <- 20  # number of trials  
  BMSY_B0bias <- array(rlnorm(nsim * ntest, mconv(1, OM@BMSY_B0cv), sdconv(1, OM@BMSY_B0cv)), dim = c(nsim, ntest))  # trial samples of BMSY relative to unfished  
  # test <- array(BMSY_B0 * BMSY_B0bias, dim = c(nsim, ntest))  # the simulated observed BMSY_B0 
  test <- array(SSBMSY_SSB0 * BMSY_B0bias, dim = c(nsim, ntest))  # the simulated observed BMSY_B0 
  indy <- array(rep(1:ntest, each = nsim), c(nsim, ntest))  # index
  
  # indy[test > 0.9] <- NA  # interval censor
  indy[test > max(0.9, max(SSBMSY_SSB0))] <- NA  # interval censor
  
  BMSY_B0bias <- BMSY_B0bias[cbind(1:nsim, apply(indy, 1, min, na.rm = T))]  # sample such that BMSY_B0<90%
  ObsPars$BMSY_B0bias <- BMSY_B0bias
  
  # --- Implementation error time series ----
  
  TAC_f <- array(rlnorm(proyears * nsim, mconv(TACFrac, TACSD),
                        sdconv(TACFrac, TACSD)), c(nsim, proyears))  # composite of TAC fraction and error
  
  E_f <- array(rlnorm(proyears * nsim, mconv(EFrac, ESD),
                      sdconv(EFrac, ESD)), c(nsim, proyears))  # composite of TAC fraction and error
  
  SizeLim_f<-array(rlnorm(proyears * nsim, mconv(SizeLimFrac, SizeLimSD),
                          sdconv(SizeLimFrac, SizeLimSD)), c(nsim, proyears))  # composite of TAC fraction and error
  
  
  # --- Populate Data object with Historical Data ---- 
  Data <- new("Data", stock = "MSE")  # create a blank DLM data object
  if (reps == 1) Data <- OneRep(Data)  # make stochastic variables certain for only one rep
  Data <- replic8(Data, nsim)  # make nsim sized slots in the DLM data object
  Data@Name <- OM@Name
  Data@Year <- 1:nyears
  Data@Cat <- Cobs
  Data@Ind <- II
  Data@Rec <- apply(N[, 1, , ], c(1, 2), sum) * Recerr[, 1:nyears]
  Data@t <- rep(nyears, nsim)
  Data@AvC <- apply(Cobs, 1, mean)
  Data@Dt <- Dbias * Depletion * rlnorm(nsim, mconv(1, Derr), sdconv(1, Derr))
  Data@Mort <- M * Mbias
  Data@FMSY_M <- FMSY_M * FMSY_Mbias
  # Data@BMSY_B0 <- BMSY_B0 * BMSY_B0bias
  Data@BMSY_B0 <- SSBMSY_SSB0 * BMSY_B0bias
  Data@Cref <- MSY * Crefbias
  Data@Bref <- VBMSY * Brefbias
  Data@Iref <- Iref * Irefbias
  Data@LFC <- LFC * LFCbias
  Data@LFS <- LFS[nyears,] * LFSbias
  Data@CAA <- CAA
  Data@Dep <- Dbias * Depletion * rlnorm(nsim, mconv(1, Derr), sdconv(1, Derr))
  Data@Abun <- A * Abias * rlnorm(nsim, mconv(1, Aerr), sdconv(1, Aerr))
  Data@SpAbun <- Asp * Abias * rlnorm(nsim, mconv(1, Aerr), sdconv(1, Aerr))
  Data@vbK <- K * Kbias
  Data@vbt0 <- t0 * t0bias
  Data@LenCV <- LenCV * LenCVbias
  Data@vbLinf <- Linf * Linfbias
  Data@L50 <- L50 * lenMbias
  Data@L95 <- L95 * lenMbias
  Data@L95[Data@L95 > 0.9 * Data@vbLinf] <- 0.9 * Data@vbLinf[Data@L95 > 0.9 * Data@vbLinf]  # Set a hard limit on ratio of L95 to Linf
  Data@L50[Data@L50 > 0.9 * Data@L95] <- 0.9 * Data@L95[Data@L50 > 0.9 * Data@L95]  # Set a hard limit on ratio of L95 to Linf
  Data@steep <- hs * hbias
  Data@CAL_bins <- CAL_bins
  Data@CAL <- CAL
  MLbin <- (CAL_bins[1:(length(CAL_bins) - 1)] + CAL_bins[2:length(CAL_bins)])/2
  temp <- CAL * rep(MLbin, each = nsim * nyears)
  Data@ML <- apply(temp, 1:2, sum)/apply(CAL, 1:2, sum)
  Data@Lc <- array(MLbin[apply(CAL, 1:2, which.max)], dim = c(nsim, nyears))
  nuCAL <- CAL
  for (i in 1:nsim) for (j in 1:nyears) nuCAL[i, j, 1:match(max(1, Data@Lc[i, j]), MLbin)] <- NA
  temp <- nuCAL * rep(MLbin, each = nsim * nyears)
  Data@Lbar <- apply(temp, 1:2, sum, na.rm=TRUE)/apply(nuCAL, 1:2, sum, na.rm=TRUE)
  Data@MaxAge <- maxage
  Data@Units <- "unitless"
  Data@Ref <- OFLreal
  Data@Ref_type <- "Simulated OFL"
  Data@wla <- rep(a, nsim)
  Data@wlb <- rep(b, nsim)

  Data@OM <- data.frame(RefY, M, Depletion, A, SSBMSY_SSB0, FMSY_M, Mgrad, Msd, procsd, Esd, dFfinal, 
             MSY, qinc, qcv, FMSY, Linf, K, t0, hs, Linfgrad, Kgrad, Linfsd, recgrad, Ksd, 
             ageM, L5=L5[nyears, ], LFS=LFS[nyears, ], Vmaxlen=Vmaxlen[nyears, ], LFC, OFLreal, 
             Spat_targ, Frac_area_1, Prob_staying, AC, L50, L95, B0, N0, SSB0, BMSY_B0,
             TACSD,TACFrac,ESD,EFrac,SizeLimSD,SizeLimFrac,DiscMort,Blow,
             BMSY, SSBMSY, Mexp) # put all the operating model parameters in one table
  
  Data@Obs <- as.data.frame(ObsPars) # put all the observation error model parameters in one table
  
  Data@LHYear <- OM@nyears  # Last historical year is nyears (for fixed MPs)
  Data@MPrec <- Cobs[, nyears]
  Data@MPeff <- rep(1, nsim)
  Data@Misc <- vector("list", nsim)
  
  # --- Return Historical Simulations and Data from last historical year ----
  if (Hist) { # Stop the model after historical simulations are complete
    message("Returning historical simulations")
    nout <- t(apply(N, c(1, 3), sum))
    vb <- t(apply(VBiomass, c(1, 3), sum))
    b <- t(apply(Biomass, c(1, 3), sum))
    ssb <- t(apply(SSB, c(1, 3), sum))
    Cc <- t(apply(CB, c(1,3), sum))
    rec <- t(apply(N[, 1, , ], c(1,2), sum))
    
    TSdata <- list(VB=vb, SSB=ssb, Bio=b, Catch=Cc, Rec=rec, N=nout, E_f=E_f,TAC_f=TAC_f,SizeLim_f=SizeLim_f)
    AtAge <- list(Len_age=Len_age, Wt_age=Wt_age, Sl_age=V, Mat_age=Mat_age, 
                  Nage=apply(N, c(1:3), sum), SSBage=apply(SSB, c(1:3), sum), M_ageArray=M_ageArray)
    MSYs <- list(MSY=MSY, FMSY=FMSY, VBMSY=VBMSY, UMSY=UMSY, 
                 SSBMSY=SSBMSY, BMSY_B0=BMSY_B0, SSBMSY_SSB0=SSBMSY_SSB0, SSB0=SSB0, B0=B0)
    
    # updated sampled pars
    SampPars <- list(dep=dep, Esd=Esd, Find=Find, procsd=procsd, AC=AC, M=M, Msd=Msd, 
                     Mgrad=Mgrad, hs=hs, Linf=Linf, Linfsd=Linfsd, Linfgrad=Linfgrad, recgrad=recgrad,
                     K=K, Ksd=Ksd, Kgrad=Kgrad, t0=t0, L50=L50, L50_95=L50_95, Spat_targ=Spat_targ,
                     Frac_area_1=Frac_area_1, Prob_staying=Prob_staying, Size_area_1=Size_area_1, 
                     Csd=Csd, Cbias=Cbias, CAA_nsamp=CAA_nsamp, CAA_ESS=CAA_ESS, CAL_nsamp=CAL_nsamp,
                     CAL_ESS=CAL_ESS, CALcv=CALcv, betas=betas, Isd=Isd, Derr=Derr, Dbias=Dbias, 
                     Mbias=Mbias, FMSY_Mbias=FMSY_Mbias, lenMbias=lenMbias, LFCbias=LFCbias,
                     LFSbias=LFSbias, Aerr=Aerr, Abias=Abias, Kbias=Kbias, t0bias=t0bias, 
                     Linfbias=Linfbias, Irefbias=Irefbias, Crefbias=Crefbias, Brefbias=Brefbias,
                     Recsd=Recsd, qinc=qinc, qcv=qcv, L5=L5, LFS=LFS, Vmaxlen=Vmaxlen, L5s=L5s, 
                     LFSs=LFSs, Vmaxlens=Vmaxlens, Perr=Perr, R0=R0, Mat_age=Mat_age, 
                     Mrand=Mrand, Linfrand=Linfrand, Krand=Krand, maxage=maxage, V=V, 
                     Depletion=Depletion,qs=qs, TACFrac=TACFrac,TACSD=TACSD,EFrac=EFrac,
                     ESD=ESD,SizeLimFrac=SizeLimFrac,SizeLimSD=SizeLimSD,DiscMort=DiscMort) 
    
    HistData <- list(SampPars=SampPars, TSdata=TSdata, AtAge=AtAge, MSYs=MSYs, Data=Data)
    return(HistData)	
  }
  
  # assign('Data',Data,envir=.GlobalEnv) # for debugging fun
  
  # --- Run projections ---- 
  if (is.na(MPs[1])) CheckMPs <- TRUE
  if (CheckMPs) {
    message("Determining available methods")  # print an progress report
    PosMPs <- Can(Data, timelimit = timelimit)  # list all the methods that could be applied to a DLM data object 
    if (is.na(MPs[1])) {
      MPs <- PosMPs  # if the user does not supply an argument MPs run the MSE for all available methods
      message("No MPs specified: running all available")
    }
    if (!is.na(MPs[1])) {
      cant <- MPs[!MPs %in% PosMPs]
      if (length(cant) > 0) {
        message("Cannot run some MPs: ")
        print(DLMdiag(Data, "not available", funcs1=cant, timelimit = timelimit))
      }
      MPs <- MPs[MPs %in% PosMPs]  # otherwise run the MSE for all methods that are deemed possible
    }
    if (length(MPs) == 0) {
      message(Cant(Data, timelimit = timelimit))
      stop("MSE stopped: no viable methods \n\n")  # if none of the user specied methods are possible stop the run
    }
  }
  
  nMP <- length(MPs)  # the total number of methods used
  
  MSElist <- list(Data)[rep(1, nMP)]  # create a data object for each method (they have identical historical data and branch in projected years)
  
  B_BMSYa <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected B_BMSY
  F_FMSYa <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected F_FMSY
  Ba <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected Biomass
  SSBa <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected SSB
  VBa <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected vulnerable biomass
  FMa <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected fishing mortality rate
  Ca <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected catch
  TACa <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected TAC recommendation
  Effort <- array(NA, dim = c(nsim, nMP, proyears))  # store the Effort
  PAAout <- array(NA, dim = c(nsim, nMP, maxage))  # store the population-at-age in last projection year
  CAAout <- array(NA, dim = c(nsim, nMP, maxage))  # store the catch-at-age in last projection year
  CALout <- array(NA, dim = c(nsim, nMP, nCALbins))  # store the population-at-length in last projection year
  
  # SPRa <- array(NA,dim=c(nsim,nMP,proyears)) # store the Spawning Potential Ratio

  
  # ---Begin loop over MPs ----
  mm <- 2 # for debugging
  for (mm in 1:nMP) {
    # MSE Loop over methods
    pL5 <- L5  # reset selectivity parameters for projections
    pLFS <- LFS
    pVmaxlen <- Vmaxlen
    pSLarray <- SLarray # selectivity at length array
    
    message(paste(mm, "/", nMP, " Running MSE for ", MPs[mm], sep = ""))  # print a progress report
    
    # projection arrays
    N_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    Biomass_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    VBiomass_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    SSN_P <-array(NA, dim = c(nsim, maxage, proyears, nareas))
    SSB_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    FM_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    FM_nospace <- array(NA, dim = c(nsim, maxage, proyears, nareas))  # stores prospective F before reallocation to new areas
    FML <- array(NA, dim = c(nsim, nareas))  # last apical F
    Z_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    CB_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    
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
    SAY1 <- SAYR[, 1:3]
    SYA <- as.matrix(expand.grid(1:nsim, 1, 1:maxage))  # Projection year
    SY <- SYA[, 1:2]
    SA <- SYA[, c(1, 3)]
    SAY <- SYA[, c(1, 3, 2)]
    S <- SYA[, 1]
    
    V_P <- V  # Reset vulnerability array for MP 
    
    if (SRrel[1] == 1) {
      N_P[, 1, 1, ] <- Perr[, nyears+maxage-1] * (0.8 * R0a * hs * apply(SSB[, 
                                                                             , nyears, ], c(1, 3), sum))/(0.2 * SSBpR * R0a * (1 - hs) + 
                                                                                                            (hs - 0.2) * apply(SSB[, , nyears, ], c(1, 3), sum))  # Recruitment assuming regional R0 and stock wide steepness
    } else {
      # most transparent form of the Ricker uses alpha and beta params
      N_P[, 1, 1, ] <- Perr[, nyears+maxage-1] * aR * apply(SSB[, , nyears, ], c(1, 3), sum) * exp(-bR * apply(SSB[, , nyears, ], c(1, 3), sum))
    }
    indMov <- as.matrix(expand.grid(1:nareas, 1:nareas, 1, 1:maxage, 1:nsim)[5:1])
    indMov2 <- indMov[, c(1, 2, 3, 4)]
    indMov3 <- indMov[, c(1, 4, 5)]
    
    N_P[, 2:maxage, 1, ] <- N[, 1:(maxage - 1), nyears, ] * exp(-Z[, 1:(maxage - 1), nyears, ])  # Total mortality
    temp <- array(N_P[indMov2] * mov[indMov3], dim = c(nareas, nareas, maxage, nsim))  # Move individuals
    N_P[, , 1, ] <- apply(temp, c(4, 3, 1), sum)
    Biomass_P[SAYR] <- N_P[SAYR] * Wt_age[SAY1]  # Calculate biomass
    VBiomass_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt]  # Calculate vulnerable biomass
    SSN_P[SAYR] <- N_P[SAYR] * Mat_age[SA1]  # Calculate spawning stock numbers
    SSB_P[SAYR] <- SSN_P[SAYR] * Wt_age[SAY1]
    FML <- apply(FM[, , nyears, ], c(1, 3), max)
    
    y <- 1 
    if (class(match.fun(MPs[mm])) == "Output") {

      Data <- Sam(MSElist[[mm]], MPs = MPs[mm], perc = pstar, reps = reps)

      TACused <- apply(Data@TAC, 3, quantile, p = pstar, na.rm = T)
      # if MP returns NA - TAC is set to catch from last year
      TACused[is.na(TACused)] <- apply(CB, c(1,3), sum)[is.na(TACused), nyears]
      
      TACa[, mm, 1] <- TACused                                               # TAC recommendation
      TACused<- TAC_f[,1]*TACused                                            # TAC taken after implementation error
      availB <- MSElist[[mm]]@OM$A # apply(VBiomass_P[,,1,], 1, sum) # total available biomass
      
      maxC <- (1 - exp(-maxF)) * availB                                      # max catch given maxF
      # if the TAC is higher than maxC than catch is equal to maxC
      notNA <- which(!is.na(TACused) & !is.na(availB)) # robustify for MPs that return NA 
      TACused[notNA][TACused[notNA] > maxC[notNA]] <- maxC[notNA][TACused[notNA] > maxC[notNA]]
      
      fishdist <- (apply(VBiomass_P[, , 1, ], c(1, 3), sum)^Spat_targ)/
        apply(apply(VBiomass_P[, , 1, ], c(1, 3), sum)^Spat_targ, 1, mean)  # spatial preference according to spatial biomass
      
      # CB_P[SAYR] <- Biomass_P[SAYR] * (1 - exp(-V_P[SAYt] * fishdist[SR]))  # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space
      CB_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt] * fishdist[SR] # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space
      
      temp <- CB_P[, , 1, ]/apply(CB_P[, , 1, ], 1, sum)  # how catches are going to be distributed
      CB_P[, , 1, ] <- TACused * temp  # debug - to test distribution code make TAC = TAC2, should be identical
      
      # temp <- CB_P[SAYR]/(Biomass_P[SAYR] * exp(-Marray[SYt]/2))  # Pope's approximation	  
      temp <- CB_P[SAYR]/(Biomass_P[SAYR] * exp(-M_ageArray[SAYt]/2))  # Pope's approximation
 
      temp[temp > (1 - exp(-maxF))] <- 1 - exp(-maxF)
      FM_P[SAYR] <- -log(1 - temp)
      
      # Z_P[SAYR] <- FM_P[SAYR] + Marray[SYt]
      Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt]
      
      Effort[, mm, y] <- (-log(1 - apply(CB_P[, , y, ], 1, sum)/(apply(CB_P[, , y, ], 1, sum) + apply(VBiomass_P[, , y, ], 1, sum))))/qs	  
      
    } else {
      # input control
  
      runIn <- runInMP(MSElist[[mm]], MPs = MPs[mm], reps = reps)  # Apply input control MP
   
      inc <- runIn[[1]] # input control recommendations 
      Data <- runIn[[2]] # Data object object with saved info from MP 
      
      Ai <- inc[1, , 1] 
      Ei <- inc[2, , 1] # effort
      Effort[, mm, y] <- Ei *E_f[,y] # Change in Effort
      Si <- t(inc[3:4, , 1]) 
      newSel <- inc[5:6, , 1] # new selectivity parameters (L5 and LFS)
      newUppLim <- inc[7, , 1] # new upper size limit (ignored if NA) 
      newVmax <- inc[8, , 1] # new vulnerability at max length
      
      chngSel <- which(colSums(apply(newSel, 2, is.na)) == 0)  # selectivity pattern changed in which sims?
      ind <- as.matrix(expand.grid((y+nyears):(nyears+proyears), chngSel))
      if (length(chngSel) > 0) {
        pL5[ind] <- newSel[1, ind[,2]]	# update size of first capture for future years 
        pLFS[ind] <- newSel[2, ind[,2]] # update size of first full selection for future years 
        if (any(!is.na(inc[8, , 1]))) {
          ind <- which(!is.na(inc[8, , 1])) # update Vmaxlen for future years where applicable
          ind2 <- as.matrix(expand.grid((y+nyears):(nyears+proyears), ind))
          pVmaxlen[ind2] <- inc[8, ind2[,2], 1]
        }
      }
      
      Vi <- t(sapply(1:nsim, SelectFun, pL5[y + nyears, ]*SizeLim_f[,y], pLFS[y + nyears, ]*SizeLim_f[,y], 
                     pVmaxlen[y + nyears, ], Len_age[, maxage, nyears], Len_age[, , y + nyears])) # update vulnerability-at-age schedule with implementation error on L5 and LFS
      
      ind <- as.matrix(expand.grid(1:nsim, 1:length(CAL_binsmid), (y+nyears):(nyears+proyears)))
      pSLarray[ind] <- t(sapply(1:nsim, SelectFun, SL0.05=pL5[y+nyears, ]*SizeLim_f[,y], SL1=pLFS[y+nyears, ]*SizeLim_f[,y], 
                                MaxSel=pVmaxlen[y+nyears, ], maxlens=maxlen, Lens=CAL_binsmid)) # update vulnerability-at-length schedule with implementation error on L5 and LFS
      
      # Maximum Size Limit - upper size limit has been set
      if (!all(is.na(newUppLim))) {
        Vi[Len_age[, , (y + nyears)] >= newUppLim] <- 0
        for (ss in 1:nsim) {
          index <- which(CAL_binsmid >= newUppLim[ss])
          pSLarray[ss, index, (y+nyears):(nyears+proyears)] <- 0 
        }	
      }
      # Vuln flag
      Vchange <- any(!is.na(inc[5:8]))
      
      if (sum(Si != 1) == 0) {
        # if there is no spatial closure if no vulnerability schedule is
        # specified
        if (!Vchange) {
          newVB <- apply(VBiomass_P[, , y, ], c(1, 3), sum)  # vulnerability isn't changed
          fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ, 1, mean)  # spatial preference according to spatial biomass
          FM_P[SAYR] <- FinF[S1] * Ei[S1] * V_P[SAYt] * fishdist[SR] * 
            qvar[SY1] * qs[S1] * (1 + qinc[S1]/100)^y  # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
        } else {
          if (y < proyears) 
            V_P[, , (nyears + 1):(proyears + nyears)] <- Vi  # Update vulnerability schedule for all future years  
          newVB <- apply(VBiomass_P[, , y, ] * Vi[SA1], c(1, 3), 
                         sum)  # vulnerability modified
          fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ, 
                                              1, mean)  # spatial preference according to spatial biomass
          FM_P[SAYR] <- FinF[S1] * Ei[S1] * Vi[SA1] * fishdist[SR] * 
            qvar[SY1] * qs[S1] * (1 + qinc[S1]/100)^y  # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
        }
      } else {
        # A spatial closure if no vulnerability schedule is specified
        if (!Vchange) {
          newVB <- apply(VBiomass_P[, , y, ], c(1, 3), sum)  # vulnerability isn't changed
          fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ, 1, mean)  # spatial preference according to spatial biomass
          Emult <- 1 + ((2/apply(fishdist * Si, 1, sum)) - 1) * Ai  # allocate effort to new area according to fraction allocation Ai
          FM_P[SAYR] <- FinF[S1] * Ei[S1] * V_P[SAYt] * Si[SR] * fishdist[SR] * Emult[S1] * qvar[SY1] * qs[S1]^(1 +  qinc[S1]/100)^y
        } else {
          if (y < proyears) 
            V_P[, , (nyears + 1):(proyears + nyears)] <- Vi  # Update vulnerability schedule for all future years
          newVB <- apply(VBiomass_P[, , y, ] * Vi[SA1], c(1, 3), sum)  # vulnerability modified
          fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ, 1, mean)  # spatial preference according to spatial biomass
          Emult <- 1 + ((2/apply(fishdist * Si, 1, sum)) - 1) *    Ai  # allocate effort to new area according to fraction allocation Ai
          FM_P[SAYR] <- FinF[S1] * Ei[S1] * Vi[SA1] * Si[SR] * fishdist[SR] * 
            Emult[S1] * qvar[SY1] * qs[S1]^(1 + qinc[S1]/100)^y
        }  # vulnerability specified
      }  # spatial closure specified  
      
      VBiomass_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt]  # update vulnerable biomass 
      # Z_P[SAYR] <- FM_P[SAYR] + Marray[SYt] # calculate total mortality 
      Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt] # calculate total mortality 
      CB_P[SAYR] <- FM_P[SAYR]/Z_P[SAYR] * Biomass_P[SAYR] * (1 - exp(-Z_P[SAYR]))  	   
      # CB_P[SAYR] <- FM_P[SAYR]/Z_P[SAYR] * VBiomass_P[SAYR] * (1 - exp(-Z_P[SAYR]))
    }  # input control  
    
    
    # CB_P[SAYR] <- Biomass_P[SAYR]*(1-exp(-FM_P[SAYR]))
    
    
    # TACa[, mm, 1] <- apply(CB_P[, , 1, ], 1, sum)  # Adjust TAC to actual catch in the year 
    # To account for years where TAC is higher than catch
    
    upyrs <- 1 + (0:(floor(proyears/interval) - 1)) * interval  # the years in which there are updates (every three years)
    cat(".")
    flush.console()
    
    # --- Begin projection years ----
    for (y in 2:proyears) {
      cat(".")
      flush.console()
      if (class(match.fun(MPs[mm])) == "Output")  TACa[, mm, y] <- TACa[, mm, y-1] # TAC same as last year unless changed 
      SAYRt <- as.matrix(expand.grid(1:nsim, 1:maxage, y + nyears, 1:nareas))  # Trajectory year
      SAYt <- SAYRt[, 1:3]
      SAYtMP <- cbind(SAYt, mm)
      SYt <- SAYRt[, c(1, 3)]
      SAY1R <- as.matrix(expand.grid(1:nsim, 1:maxage, y - 1, 1:nareas))
      SAYR <- as.matrix(expand.grid(1:nsim, 1:maxage, y, 1:nareas))
      SY <- SAYR[, c(1, 3)]
      SA <- SAYR[, 1:2]
      S1 <- SAYR[, 1]
      
      SAY <- SAYR[, 1:3]
      S <- SAYR[, 1]
      SR <- SAYR[, c(1, 4)]
      SA2YR <- as.matrix(expand.grid(1:nsim, 2:maxage, y, 1:nareas))
      SA1YR <- as.matrix(expand.grid(1:nsim, 1:(maxage - 1), y -1, 1:nareas))
      indMov <- as.matrix(expand.grid(1:nareas, 1:nareas, y, 1:maxage, 1:nsim)[5:1])
      indMov2 <- indMov[, c(1, 2, 3, 4)]
      indMov3 <- indMov[, c(1, 4, 5)]
      
      N_P[SA2YR] <- N_P[SA1YR] * exp(-Z_P[SA1YR])  # Total mortality
      if (SRrel[1] == 1) {
        N_P[, 1, y, ] <- Perr[, y + nyears+maxage-1] * (0.8 * R0a * hs * 
                                                          apply(SSB_P[, , y - 1, ], c(1, 3), sum))/(0.2 * SSBpR * 
                                                                                                      R0a * (1 - hs) + (hs - 0.2) * apply(SSB_P[, , y - 1, ], c(1, 3), sum))  # Recruitment assuming regional R0 and stock wide steepness
      } else {
        # most transparent form of the Ricker uses alpha and beta params
        N_P[, 1, y, ] <- Perr[, y + nyears+maxage-1] * aR *
          apply(SSB_P[, , y - 1, ], c(1, 3), sum) * exp(-bR * apply(SSB_P[, , y - 1, ], c(1, 3), sum))
      }
      
      temp <- array(N_P[indMov2] * mov[indMov3], 
                    dim = c(nareas, nareas, maxage, nsim))  # Move individuals
      N_P[, , y, ] <- apply(temp, c(4, 3, 1), sum)
      
      Biomass_P[SAYR] <- N_P[SAYR] * Wt_age[SAYt]  # Calculate biomass
      VBiomass_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt]  # Calculate vulnerable biomass
      SSN_P[SAYR] <- N_P[SAYR] * Mat_age[SA]  # Calculate spawning stock numbers
      SSB_P[SAYR] <- SSN_P[SAYR] * Wt_age[SAYt]  # Calculate spawning stock biomass
      
      # --- An update year ----
      if (y %in% upyrs) {
        # rewrite the DLM object and run the TAC function
        yind <- upyrs[match(y, upyrs) - 1]:(upyrs[match(y, upyrs)] - 1)
        CNtemp <- array(N_P[, , yind, ] * exp(Z_P[, , yind, ]) * 
                          (1 - exp(-Z_P[, , yind, ])) * (FM_P[, , yind, ]/Z_P[, , yind, ]), c(nsim, maxage, interval, nareas))
        CBtemp <- array(Biomass_P[, , yind, ] * exp(Z_P[, , yind, ]) * 
                          (1 - exp(-Z_P[, , yind, ])) * (FM_P[, , yind, ]/Z_P[, , yind, ]), c(nsim, maxage, interval, nareas))
        CNtemp[is.na(CNtemp)] <- tiny
        CBtemp[is.na(CBtemp)] <- tiny
        CNtemp[!is.finite(CNtemp)] <- tiny
        CBtemp[!is.finite(CBtemp)] <- tiny
        CNtemp <- apply(CNtemp, c(1, 3, 2), sum, na.rm = T)
        
        Cobs <- Cbiasa[, nyears + yind] * Cerr[, nyears + yind] * 
          apply(CBtemp, c(1, 3), sum, na.rm = T)
        Cobs[is.na(Cobs)] <- tiny
        Recobs <- Recerr[, nyears + yind] * apply(array(N_P[, 1, yind, ], c(nsim, interval, nareas)), c(1, 2), sum)
        
        cond <- apply(CNtemp, 1:2, sum, na.rm = T) < 1  # this is a fix for low sample sizes. If CN is zero across the board a single fish is caught in age class of model selectivity (dumb I know)
        fixind <- as.matrix(cbind(expand.grid(1:nsim, 1:interval), 
                                  rep(floor(maxage/3), interval)))  # more fix
        
        # assign('fixind',fixind,envir=.GlobalEnv) # for debugging fun
        # assign('CNtemp',CNtemp,envir=.GlobalEnv) # for debugging fun
        
        CNtemp[fixind[cond, ]] <- 1  # puts a catch in the most vulnerable age class
        CNtemp[is.na(CNtemp)] <- tiny 
        
        CAA <- array(NA, dim = c(nsim, interval, maxage))  # Catch  at age array
        # for(i in 1:nsim)for(j in
        # 1:interval)CAA[i,j,]<-ceiling(-0.5+rmultinom(1,CAA_nsamp[i],CNtemp[i,j,])*CAA_nsamp[i]/CAA_ESS[i])
        # # a multinomial observation model for catch-at-age data
        for (i in 1:nsim) {
          for (j in 1:interval) {
            CAA[i, j, ] <- ceiling(-0.5 + 
                                     rmultinom(1, CAA_ESS[i], CNtemp[i, j, ]) * CAA_nsamp[i]/CAA_ESS[i])   # a multinomial observation model for catch-at-age data
            # rmultinom(1, CAA_ESS[i], CN[i, j, ]) * CAA_nsamp[i]/CAA_ESS[i])   # a multinomial observation model for catch-at-age data
          }
        }	  
        
        CAL <- array(NA, dim = c(nsim, interval, nCALbins))  # the catch at length array
        # # a multinomial observation model for catch-at-length data
        # cn <- as.matrix(CNtemp[i,,])
        cn <- as.matrix(V_P[i,,yind])
        if (interval > 1) cn <- aperm(cn, c(2,1))
        if (interval == 1) cn <- t(cn) # dodgy hack to ensure matrix is correct 
        for (i in 1:nsim) { # Rcpp code 
          CAL[i, 1:interval, ] <- genLenComp(CAL_bins, CAL_binsmid, 
                                             as.matrix(pSLarray[i,, nyears + yind]), 
                                             CAL_ESS[i], CAL_nsamp[i], 
                                             cn, as.matrix(Len_age[i,,nyears + yind]), 
                                             as.matrix(LatASD[i,, nyears + yind]), truncSD=0) 
          LFC[i] <- CAL_binsmid[min(which(round(CAL[i, interval, ],0) >= 1))] # get the smallest CAL observation	
        }	
        # for (i in 1:nsim) {
        # for (j in 1:interval) {
        # yy <- yind[j]
        # # tempCN <- ceiling(-0.5 + rmultinom(1, size = CAL_ESS[i], prob = CN[i, j, ]) * CAL_nsamp[i]/CAL_ESS[i])
        # tempCN <- ceiling(-0.5 + rmultinom(1, size = CAL_ESS[i], prob = CNtemp[i, j , ]) * CAL_nsamp[i]/CAL_ESS[i])
        # # ages <- rep(1:maxage,tempCN)+runif(sum(tempCN),-0.5,0.5) # sample
        # # expected age
        # lens <- unlist(sapply(1:maxage, function(X) 
        # rnorm(tempCN[X], Len_age[i, X, yy + nyears], LatASD[i, X, yy + nyears])))
        # lens[lens > (max(Linfarray) + 2 * max(LatASD)) | lens > 
        # max(CAL_bins)] <- max(Linfarray) + 2 * max(LatASD)  # truncate at 2 sd 
        # CAL[i, j, ] <- hist(lens, CAL_bins, plot = F)$counts  # assign to bins
        # LFC[i] <- min(c(lens, LFC[i]), na.rm = T)  # get the smallest CAL observation
        
        # }
        # }
        
        
        I2 <- cbind(apply(Biomass, c(1, 3), sum), apply(Biomass_P, c(1, 3), sum)[, 1:(y - 1)]) * 
          Ierr[, 1:(nyears + (y - 1))]^betas
        I2[is.na(I2)] <- tiny
        I2 <- I2/apply(I2, 1, mean)
        
        # Depletion <- apply(Biomass_P[, , y, ], 1, sum)/apply(Biomass[, , 1, ], 1, sum)
        Depletion <- apply(SSB_P[, , y, ], 1, sum)/SSB0 # apply(SSB[, , 1, ], 1, sum)
        Depletion[Depletion < tiny] <- tiny
        A <- apply(VBiomass_P[, , y, ], 1, sum)
        A[is.na(A)] <- tiny
        Asp <- apply(SSB_P[, , y, ], 1, sum)  # SSB Abundance
        Asp[is.na(Asp)] <- tiny
        OFLreal <- A * FMSY
        
        # assign all the new data
        MSElist[[mm]]@OM$A <- A
        MSElist[[mm]]@Year <- 1:(nyears + y - 1)
        MSElist[[mm]]@Cat <- cbind(MSElist[[mm]]@Cat, Cobs)
        MSElist[[mm]]@Ind <- I2
        MSElist[[mm]]@Rec <- cbind(MSElist[[mm]]@Rec, Recobs)
        MSElist[[mm]]@t <- rep(nyears + y, nsim)
        MSElist[[mm]]@AvC <- apply(MSElist[[mm]]@Cat, 1, mean)
        MSElist[[mm]]@Dt <- Dbias * Depletion * rlnorm(nsim, mconv(1, Derr), sdconv(1, Derr))
        oldCAA <- MSElist[[mm]]@CAA
        MSElist[[mm]]@CAA <- array(0, dim = c(nsim, nyears + y - 1, maxage))
        MSElist[[mm]]@CAA[, 1:(nyears + y - interval - 1), ] <- oldCAA
        MSElist[[mm]]@CAA[, nyears + yind, ] <- CAA
        MSElist[[mm]]@Dep <- Dbias * Depletion * rlnorm(nsim, mconv(1, Derr), sdconv(1, Derr))
        MSElist[[mm]]@Abun <- A * Abias * rlnorm(nsim, mconv(1, Aerr), sdconv(1, Aerr))
        MSElist[[mm]]@SpAbun <- Asp * Abias * rlnorm(nsim, mconv(1, Aerr), sdconv(1, Aerr))
        MSElist[[mm]]@CAL_bins <- CAL_bins
        oldCAL <- MSElist[[mm]]@CAL
        MSElist[[mm]]@CAL <- array(0, dim = c(nsim, nyears + y - 1, nCALbins))
        MSElist[[mm]]@CAL[, 1:(nyears + y - interval - 1), ] <- oldCAL
        MSElist[[mm]]@CAL[, nyears + yind, ] <- CAL[, 1:interval, ]
        
        temp <- CAL * rep(MLbin, each = nsim * interval)
        MSElist[[mm]]@ML <- cbind(MSElist[[mm]]@ML, apply(temp, 1:2, sum)/apply(CAL, 1:2, sum))
        MSElist[[mm]]@Lc <- cbind(MSElist[[mm]]@Lc, array(MLbin[apply(CAL, 1:2, which.max)], dim = c(nsim, interval)))
        nuCAL <- CAL
        for (i in 1:nsim) for (j in 1:interval) nuCAL[i, j, 1:match(max(1, MSElist[[mm]]@Lc[i, j]), MLbin)] <- NA 
        temp <- nuCAL * rep(MLbin, each = nsim * interval)
        MSElist[[mm]]@Lbar <- cbind(MSElist[[mm]]@Lbar, apply(temp,1:2, sum, na.rm=TRUE)/apply(nuCAL, 1:2, sum, na.rm=TRUE))
        
        MSElist[[mm]]@LFC <- LFC * LFCbias
        MSElist[[mm]]@LFS <- pLFS[nyears + y,] * LFSbias 
        
        MSElist[[mm]]@Ref <- OFLreal
        MSElist[[mm]]@Ref_type <- "Simulated OFL"
        MSElist[[mm]]@Misc <- Data@Misc
        
        # assign('Data',MSElist[[mm]],envir=.GlobalEnv) # for debugging fun
        
        if (class(match.fun(MPs[mm])) == "Output") {
          Data <- Sam(MSElist[[mm]], MPs = MPs[mm], perc = pstar, reps = reps)
          TACused <- apply(Data@TAC, 3, quantile, p = pstar, na.rm = TRUE)  #
          NAs <- which(is.na(TACused))
          # If MP returns NA - the TAC from last year is used. 
          TACused[is.na(TACused)] <- TACa[is.na(TACused), mm, y-1]
          if (length(NAs) > 0) {
            # robustifying TAC setting!
            TACused[NAs] <- TACa[NAs, mm, y - 1]  #
            if (!exists("store")) 
              store <- list()
            store <- append(store, c(MPs[mm], NAs))
          }
          TACa[, mm, y] <- TACused
          MSElist[[mm]]@MPrec <- TACused
          TACused<- TAC_f[,y]*TACused    # after implementation error   
          
          availB <- availB <- MSElist[[mm]]@OM$A # apply(VBiomass_P[,,y,], 1, sum) # total available biomass
          maxC <- (1 - exp(-maxF)) * availB
          # if the TAC is higher than maxC than catch is equal to maxC
          notNA <- which(!is.na(TACused) & !is.na(availB))
          TACused[TACused[notNA] > maxC[notNA]] <- maxC[TACused[notNA] > maxC[notNA]]
          # TACused[TACused > maxC] <- maxC[TACused > maxC] 	
          
          fishdist <- (apply(VBiomass_P[, , y, ], c(1, 3), sum)^Spat_targ)/apply(apply(VBiomass_P[, , y, ], c(1, 3), sum)^Spat_targ, 1, mean)  # spatial preference according to spatial biomass     
          # CB_P[SAYR] <- Biomass_P[SAYR] * (1 - exp(-V_P[SAYt] *  fishdist[SR]))  # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space          
          CB_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt] * fishdist[SR]  # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space          
          
          temp <- CB_P[, , y, ]/apply(CB_P[, , y, ], 1, sum)  # how catches are going to be distributed
          CB_P[, , y, ] <- TACused * temp  # debug - to test distribution code make TAC = TAC2, should be identical          
          # temp <- CB_P[SAYR]/(Biomass_P[SAYR] * exp(-Marray[SYt]/2))  # Pope's approximation
          temp <- CB_P[SAYR]/(Biomass_P[SAYR] * exp(-M_ageArray[SAYt]/2))  # Pope's approximation
          temp[temp > (1 - exp(-maxF))] <- 1 - exp(-maxF)
          FM_P[SAYR] <- -log(1 - temp)
          # Z_P[SAYR] <- FM_P[SAYR] + Marray[SYt]
          Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt]
          Effort[, mm, y] <- (-log(1 - apply(CB_P[, , y, ], 1, sum)/(apply(CB_P[, , y, ], 1, sum) + apply(VBiomass_P[, , y, ], 1, sum))))/qs
          
      
        } else {
          MSElist[[mm]]@MPeff <- Ei
          runIn <- runInMP(MSElist[[mm]], MPs = MPs[mm], reps = reps)  # Apply input control MP
          inc <- runIn[[1]]
          Data <- runIn[[2]]
          Ai <- inc[1, , 1]
          Ei <- inc[2, , 1]
          Effort[, mm, y] <- Ei * E_f[,y] # Change in Effort
          Si <- t(inc[3:4, , 1])
          newSel <- (inc[5:6, , 1])
          
          newUppLim <- inc[7, , 1]
          newVmax <- inc[8, , 1]
          
          chngSel <- which(colSums(apply(newSel, 2, is.na)) == 0)  # selectivity pattern changed
          ind <- as.matrix(expand.grid((y+nyears):(nyears+proyears), chngSel))		  
          if (length(chngSel) > 0) {
            pL5[ind] <- newSel[1, ind[,2]]	# update size of first capture for future years 
            pLFS[ind] <- newSel[2, ind[,2]] # update size of first full selection for future years 
            if (any(!is.na(inc[8, , 1]))) {
              ind <- which(!is.na(inc[8, , 1])) # update Vmaxlen for future years where applicable
              ind2 <- as.matrix(expand.grid((y+nyears):(nyears+proyears), ind))
              pVmaxlen[ind2] <- inc[8, ind2[,2], 1]
            }
          }
          
          Vi <- t(sapply(1:nsim, SelectFun, pL5[y + nyears, ]*SizeLim_f[,y], pLFS[y + nyears, ]*SizeLim_f[,y], 
                         pVmaxlen[y + nyears, ], Len_age[, maxage, nyears], Len_age[, , y + nyears])) # update vulnerability-at-age schedule with implementation error on L5 and LFS
          
          ind <- as.matrix(expand.grid(1:nsim, 1:length(CAL_binsmid), (y+nyears):(nyears+proyears)))
          pSLarray[ind] <- t(sapply(1:nsim, SelectFun, SL0.05=pL5[y+nyears, ]*SizeLim_f[,y], SL1=pLFS[y+nyears, ]*SizeLim_f[,y], 
                                    MaxSel=pVmaxlen[y+nyears, ], maxlens=maxlen, Lens=CAL_binsmid)) # update vulnerability-at-length schedule with implementation error on L5 and LFS
          
          # Maximum Size Limit - upper size limit has been set
          if (!all(is.na(newUppLim))) {
            Vi[Len_age[, , (y + nyears)] >= newUppLim] <- 0
            for (ss in 1:nsim) {
              index <- which(CAL_binsmid >= newUppLim[ss])
              pSLarray[ss, index, (y+nyears):(nyears+proyears)] <- 0 
            }	
          }
          
          # Vuln flag
          Vchange <- any(!is.na(inc[5:8]))
          
          if (sum(Si != 1) == 0) {
            # if there is no spatial closure if no vulnerability schedule is
            # specified
            if (!Vchange) {
              newVB <- apply(VBiomass_P[, , y, ], c(1, 3), sum)  # vulnerability isn't changed
              fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ,  1, mean)  # spatial preference according to spatial biomass
              FM_P[SAYR] <- FinF[S1] * Ei[S1] * V_P[SAYt] * fishdist[SR] * 
                qvar[SY] * qs[S1] * (1 + qinc[S1]/100)^y  # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
            } else {
              if (y < proyears) V_P[, , (y + nyears + 1):(proyears + nyears)] <- Vi  # Update vulnerability schedule for all future years
              newVB <- apply(Biomass_P[, , y, ] * Vi[SA], c(1, 3), sum)  # vulnerability modified
              fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ, 1, mean)  # spatial preference according to spatial biomass
              FM_P[SAYR] <- FinF[S1] * Ei[S1] * Vi[SA] * fishdist[SR] *  qvar[SY] * 
                qs[S1] * (1 + qinc[S1]/100)^y  # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass         
            }
          } else {
            # A spatial closure if no vulnerability schedule is specified
            if (!Vchange) {
              newVB <- apply(Biomass_P[, , y, ], c(1, 3), sum)  # vulnerability isn't changed
              fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ,1, mean)  # spatial preference according to spatial biomass
              Emult <- 1 + ((2/apply(fishdist * Si, 1, sum)) - 1) * Ai  # allocate effort to new area according to fraction allocation Ai
              FM_P[SAYR] <- FinF[S1] * Ei[S1] * V_P[SAYt] * Si[SR] * fishdist[SR] * 
                Emult[S1] * qvar[SY] * qs[S1] * (1 + qinc[S1]/100)^y        
            } else {
              if (y < proyears) V_P[, , (y + nyears + 1):(proyears + nyears)] <- Vi  # Update vulnerability schedule for all future years
              newVB <- apply(Biomass_P[, , y, ] * Vi[SA], c(1, 3), sum)  # vulnerability modified
              fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ,       1, mean)  # spatial preference according to spatial biomass
              Emult <- 1 + ((2/apply(fishdist * Si, 1, sum)) -   1) * Ai  # allocate effort to new area according to fraction allocation Ai
              FM_P[SAYR] <- FinF[S1] * Ei[S1] * Vi[SA] * Si[SR] * fishdist[SR] * Emult[S1] * 
                qvar[SY] * qs[S1] * (1 + qinc[S1]/100)^y
              
            }  #vuln not changed
          }  # spatial closure
          VBiomass_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt]  # update vulnerable biomass 
          # Z_P[SAYR] <- FM_P[SAYR] + Marray[SYt]
          Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt]
          # CB_P[SAYR]<-Biomass_P[SAYR]*(1-exp(-FM_P[SAYR]))
          CB_P[SAYR] <- FM_P[SAYR]/Z_P[SAYR] * Biomass_P[SAYR] *   (1 - exp(-Z_P[SAYR]))
        }  # input or output control 
        
        # TACused <- apply(CB_P[, , y, ], 1, sum)  # Set last years TAC to actual catch from last year
        # TACa[, mm, y] <- TACused
        tempcatch <- apply(CB_P[, , y, ], 1, sum) 
        
        MSElist[[mm]]@MPrec <- tempcatch
      } else {
        # --- Not an update yr ----
        vbio <- apply(VBiomass_P[, , y, ], c(1, 3), sum)
        fishdist <- (vbio^Spat_targ)/apply(vbio^Spat_targ, 1, mean)  # calculate distribution of effort \t  
        if (class(match.fun(MPs[mm])) == "Output") {
          # CB_P[SAYR] <- Biomass_P[SAYR] * (1 - exp(-fishdist[SR] *  V_P[SAYt]))  # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space
          CB_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt] * fishdist[SR]  # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space                 
          temp <- CB_P[, , y, ]/apply(CB_P[, , y, ], 1, sum)  # how catches are going to be distributed
          tempcatch <- TACa[, mm, y]
          tempcatch<- TAC_f[,y]*tempcatch # after implementation error       
          
          availB <- apply(VBiomass_P[,,y,], 1, sum) # total available biomass
          maxC <- (1 - exp(-maxF)) * availB
          # if the TAC is higher than maxC than catch is equal to maxC
          notNA <- which(!is.na(tempcatch) & !is.na(availB))		  
          tempcatch[notNA][tempcatch[notNA] > maxC[notNA]] <- maxC[notNA][tempcatch[notNA] > maxC[notNA]]		  
          # tempcatch[tempcatch > maxC] <- maxC[tempcatch > maxC] 		 
          
          CB_P[, , y, ] <- tempcatch * temp  # debug - to test distribution code make TAC = TAC2, should be identical
          #temp <- CB_P[SAYR]/(Biomass_P[SAYR] * exp(-Marray[SYt]/2))  # Pope's approximation
          temp <- CB_P[SAYR]/(Biomass_P[SAYR] * exp(-M_ageArray[SAYt]/2))  # Pope's approximation

          temp[temp > (1 - exp(-maxF))] <- 1 - exp(-maxF)
          FM_P[SAYR] <- -log(1 - temp)
          Effort[, mm, y] <- (-log(1 - apply(CB_P[, , y, ], 1, sum)/ (apply(CB_P[, , y, ], 1, sum) + apply(VBiomass_P[, , y, ], 1, sum))))/qs
          # Z_P[SAYR] <- FM_P[SAYR] + Marray[SYt]							 
          Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt]
        } else {
          # input control FM_P[SAYR] <- FM_P[SAY1R]*qvar[SY] *(1+qinc[S1]/100)^y
          # # add fishing efficiency changes and variability
          FM_P[SAYR] <- FM_P[SAY1R] * qvar[SY] * (1 + qinc[S1]/100)  # add fishing efficiency changes and variability
          Effort[, mm, y] <-  Ei * E_f[,y]   # Effort doesn't change in non-update year
          # Z_P[SAYR] <- FM_P[SAYR] + Marray[SYt]
          Z_P[SAYR] <- FM_P[SAYR] + M_ageArray[SAYt]
          # CB_P[SAYR]<-Biomass_P[SAYR]*(1-exp(-FM_P[SAYR]))
          CB_P[SAYR] <- FM_P[SAYR]/Z_P[SAYR] * Biomass_P[SAYR] * (1 - exp(-Z_P[SAYR]))
        }
        
      }  # not an update year
      
    }  # end of year
 
    
    B_BMSYa[, mm, ] <- apply(SSB_P, c(1, 3), sum, na.rm=TRUE)/SSBMSY  # SSB relative to SSBMSY
    # F_FMSYa[, mm, ] <- (-log(1 - apply(CB_P, c(1, 3), sum)/(apply(CB_P, c(1, 3), sum) + 
    # apply(VBiomass_P, c(1, 3), sum))))/FMSY 
    # VBiomass is calculated before catches are taken
    # gives an error message if CB_P or VBiomass_P is NA 
    FMa[, mm, ] <- -log(1 - apply(CB_P, c(1, 3), sum, na.rm=TRUE)/apply(VBiomass_P, c(1, 3), sum, na.rm=TRUE))		
    
    F_FMSYa[, mm, ] <- FMa[, mm, ]/FMSY
    
    Ba[, mm, ] <- apply(Biomass_P, c(1, 3), sum, na.rm=TRUE) # biomass 
    SSBa[, mm, ] <- apply(SSB_P, c(1, 3), sum, na.rm=TRUE) # spawning stock biomass
    VBa[, mm, ] <- apply(VBiomass_P, c(1, 3), sum, na.rm=TRUE) # vulnerable biomass
    # FMa[, mm, ] <- -log(1 - apply(CB_P, c(1, 3), sum)/(apply(CB_P, c(1, 3), sum) + 
    # apply(VBiomass_P, c(1, 3), sum)))
    # VBiomass is calculated before catches are taken 				   
    
    Ca[, mm, ] <- apply(CB_P, c(1, 3), sum, na.rm=TRUE)
    
    # Store Pop and Catch-at-age and at-length for last projection year 
    PAAout[ , mm, ] <- apply(N_P[ , , proyears, ], c(1,2), sum) # population-at-age
    CNtemp <- array(N_P * exp(Z_P) * (1 - exp(-Z_P)) * (FM_P/Z_P), c(nsim, maxage, proyears, nareas))
    CAAout[ , mm, ] <- apply(CNtemp[,,proyears,], c(1, 2), sum) # nsim, maxage # catch-at-age
    CALout[ , mm, ] <- CAL[,max(dim(CAL)[2]),] # catch-at-length in last year
    
    cat("\n")
  }  # end of mm methods 
  
 
  MSEout <- new("MSE", Name = OM@Name, nyears, proyears, nMPs=nMP, MPs, nsim, 
                Data@OM, Obs=Data@Obs, B_BMSY=B_BMSYa, F_FMSY=F_FMSYa, B=Ba, 
                SSB=SSBa, VB=VBa, FM=FMa, Ca, TAC=TACa, SSB_hist = SSB, CB_hist = CB, 
                FM_hist = FM, Effort = Effort, PAA=PAAout, CAA=CAAout, CAL=CALout, CALbins=CAL_binsmid)
  # Store MSE info
  attr(MSEout, "version") <- packageVersion("DLMtool")
  attr(MSEout, "interval") <- interval
  attr(MSEout, "maxF") <- maxF
  attr(MSEout, "timelimit") <- timelimit
  attr(MSEout, "pstar") <- pstar
  attr(MSEout, "reps") <- reps
  attr(MSEout, "date") <- date()
  attr(MSEout, "R.version") <- R.version	
  MSEout 
}

#' Internal function of runMSE for checking that the OM slot cpars slot is formatted correctly
#'
#' @param cpars a list of model parameters to be sampled (single parameters are a vector nsim long, time series are matrices nsim x nyears)
#' @return either an error and the length of the first dimension of the various cpars list items or passes and returns the number of simulations
#' @export cparscheck
#' @author T. Carruthers
cparscheck<-function(cpars){

  dim1check<-function(x){
    if(class(x)=="numeric")length(x)
    else dim(x)[1]
  }

  dims<-sapply(cpars,dim1check)
  if(length(unique(dims))!=1){
    print(dims)
    stop("The custom parameters in your operating model @cpars have varying number of simulations. For each simulation each parameter / variable should correspond with one another")
  }else{
    as.integer(dims[1])
  }

}


cparnamecheck<-function(cpars){

  Sampnames <- c("dep","Esd","Find","procsd","AC","M","Msd",
                 "Mgrad","hs","Linf","Linfsd","Linfgrad","recgrad",
                 "K","Ksd","Kgrad","t0","L50","L50_95","Spat_targ",
                 "Frac_area_1","Prob_staying","Size_area_1",
                 "Csd","Cbias","CAA_nsamp","CAA_ESS","CAL_nsamp",
                 "CAL_ESS","CALcv","betas","Isd","Derr","Dbias",
                 "Mbias","FMSY_Mbias","lenMbias","LFCbias",
                 "LFSbias","Aerr","Abias","Kbias","t0bias",
                 "Linfbias","Irefbias","Crefbias","Brefbias",
                 "Recsd","qinc","qcv","L5","LFS","Vmaxlen","L5s",
                 "LFSs","Vmaxlens","Perr","R0","Mat_age",
                 "Mrand","Linfrand","Krand","maxage","V","Depletion", # end of OM variables
                 "ageM", "age95", "V", "EffYears", "EffLower", "EffUpper","Mat_age", # start of runMSE derived variables
                 "Wt_age")

}
