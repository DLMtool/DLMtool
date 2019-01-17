
Names <- c("maxage", "R0", "Mexp", "Msd", "dep", "D", "Mgrad", "SRrel", "hs", "procsd",
           "L50", "L95", "L50_95", "CAL_binsmid", "Len_age", "maxlen", "Linf", 
           "M_at_Length", "Frac_area_1", "Prob_staying", "M_ageArray", "Mat_age",
           "Wt_age", "V", "Spat_targ", "procmu", "recMulti", "Linfrand", "Krand",
           "Abias Aerr", "Brefbias", "CAA_ESS", "CAA_nsamp", "CAL_ESS", "CAL_bins", "CAL_nsamp",
           "Cbias", "Crefbias", "Csd", "Dbias", "Derr", "TAEFrac", "TAESD", "EffLower",
           "EffUpper", "EffYears", "FMSY_Mbias", "Frac_area_1", "Irefbias", "Isd", "K", "Kbias", "Kgrad",
           "Krand", "Ksd", "L5", "L5s", "LFCbias", "LFS", "LFSbias", "LFSs", "LatASD", "Linfbias", "Linfgrad",
           "Linfrand", "Linfsd", "M", "M_ageArray", "Mat_age", "Mbias", "Mrand", "Prob_staying", "Recsd",
           "SLarray", "SizeLimFrac", "SizeLimSD", "Spat_targ", "TACFrac", "TACSD", 
           "Vmaxlen", "Vmaxlens", "Wt_age", "ageM", "betas", "lenMbias", "nCALbins", "procmu", "qcv", "qinc",
           "recMulti",  "t0", "t0bias", "Abias", "Aerr", "Perr", "Esd", "qvar", "Marray",
           "Linfarray", "Karray", "t0array", "mov",  "nareas", "AC", "LenCV", "a", "b", "FinF", 
           "Fdisc", "R50", "Rslope", "retA", "retL", "LR5", "LFR", "Rmaxlen",
           "V2", "SLarray2", "DR", "Asize", "Size_area_1", "L50array", "L95array",
           "Fdisc_array", "Fdisc_array2", "Pinitdist", "DataOut",
           'Perr_y', "Cobs", "Iobs", "Dobs", "Btbiascv", 'Btobs', "h")


if(getRversion() >= "2.15.1") utils::globalVariables(Names)




#' Run a Management Strategy Evaluation
#' 
#' A function that runs a Management Strategy Evaluation (closed-loop
#' simulation) for a specified operating model
#' 
#' 
#' @param OM An operating model object (class 'OM')
#' @param MPs A vector of methods (character string) of class MP
#' @param CheckMPs Logical to indicate if \link{Can} function should be used to check
#' if MPs can be run.
#' @param timelimit Maximum time taken for a method to carry out 10 reps
#' (methods are ignored that take longer)
#' @param Hist Should model stop after historical simulations? Returns an object of 
#' class 'Hist' containing all historical data
#' @param ntrials Maximum of times depletion and recruitment deviations are 
#' resampled to optimize for depletion. After this the model stops if more than 
#' percent of simulations are not close to the required depletion
#' @param fracD Maximum allowed proportion of simulations where depletion is not 
#' close to sampled depletion from OM before model stops with error
#' @param CalcBlow Should low biomass be calculated where this is the spawning
#' biomass at which it takes HZN mean generation times of zero fishing to reach 
#' Bfrac fraction of SSBMSY
#' @param HZN The number of mean generation times required to reach Bfrac SSBMSY
#' in the Blow calculation
#' @param Bfrac The target fraction of SSBMSY for calculating Blow
#' @param AnnualMSY Logical. Should MSY statistics be calculated for each projection year? 
#' May differ from MSY statistics from last historical year if there are changes in productivity
#' @param silent Should messages be printed out to the console?
#' @param PPD Logical. Should posterior predicted data be included in the MSE object Misc slot?
#' @param parallel Logical. Should the MSE be run using parallel processing?
#' @param save_name Character. Optional name to save parallel MSE list
#' @param checks Logical. Run tests?
#' @param control control options for testing and debugging
#' @param maxsims Maximum number of simulations per packet
#' @param name Character string for name of saved MSE packets (if \code{savePack=TRUE}) 
#' and final MSE object. If none provided, it uses the first five letters from the \code{OM} name
#' @param unique Logical. Should the name be unique? Current date and time appended to name. 
#' @param maxCrash Maximum number of consecutive crashes before the MSE stops
#' @param saveMSE Logical to indicate if final MSE object should be saved to current 
#' working directory (this is probably a good idea)
#' @param savePack Logical to indicate if packets should be save to current working directory
#' @param ... Arguments to runMSE function 
#' 
#' @templateVar url running-the-mse
#' @templateVar ref NULL 
#' @template userguide_link
#' 
#' @return An object of class \linkS4class{MSE}
#' @author T. Carruthers and A. Hordyk
#' @describeIn runMSE Default function to use.
#' @seealso \link{joinMSE} \link{checkMSE} \link{updateMSE}
#' @export
runMSE <- function(OM = DLMtool::testOM, MPs = c("AvC","DCAC","FMSYref","curE","matlenlim", "MRreal"), 
                   CheckMPs = FALSE, timelimit = 1, Hist=FALSE, ntrials=50, fracD=0.05, CalcBlow=TRUE, 
                   HZN=2, Bfrac=0.5, AnnualMSY=TRUE, silent=FALSE, PPD=FALSE, parallel=FALSE, 
                   save_name=NULL, checks=FALSE, control=NULL) {
  
  if (class(OM)!='OM') stop("OM is not class 'OM'", call. = FALSE)
  
  # Set DLMenv to be empty. Currently updated by Assess models in MSEtool
  rm(list = ls(DLMenv), envir = DLMenv)
  
  # check if custom MP names already exist in DLMtool
  tt <- suppressWarnings(try(lsf.str(envir=globalenv()), silent=TRUE))
  if (class(tt)!="try-error") {
    gl.funs <- as.vector(tt)
    pkg.funs <- as.vector(ls.str('package:DLMtool'))
    if ('package:MSEtool' %in% search()) pkg.funs <- c(pkg.funs, as.vector(ls.str('package:MSEtool')))
    if (length(gl.funs)>0) {
      gl.clss <- unlist(lapply(lapply(gl.funs, get), class))
      gl.MP <- gl.funs[gl.clss %in% 'MP']
      if (length(gl.MP)>0) {
        inc.gl <- gl.MP[gl.MP %in% MPs]
        if (length(inc.gl)>0) {
          dup.MPs <- inc.gl[inc.gl %in% pkg.funs]
          if (length(dup.MPs)>0) {
            stop("Custom MP names already in DLMtool: ", paste0(dup.MPs, " "), "\nRename Custom MPs")
          }
        }
      } 
    }
  }
  
  # Check MPs 
  if (!all(is.na(MPs))) {
    for (mm in MPs) {
      chkMP <- try(get(mm), silent=TRUE)
      if (class(chkMP) != 'MP') stop(mm, " is not a valid MP", call.=FALSE) 
    }
  }

  if (parallel) {
    if (OM@nsim<48) stop("nsim must be >=48 for parallel processing", call.=FALSE)
    if(!snowfall::sfIsRunning()) {
      # stop("Parallel processing hasn't been initialized. Use 'setup'", call. = FALSE)
      message("Parallel processing hasn't been initialized. Calling 'setup()' now")
      setup()
    }
    
    if (all(is.na(MPs))) MPs <- avail("MP")
    # Export Custom MPs # 
    cMPs <- MPs[!MPs %in% pkg.funs]
    if (length(cMPs)>0) snowfall::sfExport(list=cMPs)
    
    ncpu <- snowfall::sfCpus()
    nits <- ceiling(OM@nsim/48)
    
    itsim <- rep(48,nits)
    
    if (nits < ncpu) {
      nits <- ncpu
      itsim <- rep(ceiling(OM@nsim/ncpu), ncpu)
    }
    cnt <- 1
    while(sum(itsim) != OM@nsim | any(itsim<2)) {
      diff <-  OM@nsim - sum(itsim)
      if (diff >0) {
        itsim[cnt] <- itsim[cnt]+1
      } 
      if(diff < 0) {
        itsim[cnt] <- itsim[cnt]-1
      }
      cnt <- cnt+1 
      if (cnt > length(itsim)) cnt <- 1 
    }
    
   
    if (!silent & !Hist) message("Running MSE in parallel on ", ncpu, ' processors')
    if (!silent & Hist) message("Running historical simulations in parallel on ", ncpu, ' processors')
    
    temp <- snowfall::sfClusterApplyLB(1:nits, run_parallel, itsim=itsim, OM=OM, MPs=MPs,  
                             CheckMPs=CheckMPs, timelimit=timelimit, Hist=Hist, ntrials=ntrials, 
                             fracD=fracD, CalcBlow=CalcBlow, 
                             HZN=HZN, Bfrac=Bfrac, AnnualMSY=AnnualMSY, silent=TRUE, PPD=PPD,
                             control=control, parallel=parallel)
    #assign_DLMenv() # grabs objects from DLMenv in cores, then merges and assigns to 'home' environment
  
    if (!is.null(save_name) && is.character(save_name)) saveRDS(temp, paste0(save_name, '.rdata'))
    
    MSE1 <- joinMSE(temp) 
    if (class(MSE1) == "MSE") {
      if (!silent) message("MSE completed")
    } else if (class(MSE1) == "Hist"){
      if (!silent) message("Historical simulations completed")
    } else {
      warning("MSE completed but could not join MSE objects. Re-run with `save_name ='MyName'` to debug")
    }
  }

 
  if (!parallel) {
    if (OM@nsim > 48 & !silent & !Hist) message("Suggest using 'parallel = TRUE' for large number of simulations")
    MSE1 <- runMSE_int(OM, MPs, CheckMPs, timelimit, Hist, ntrials, fracD, CalcBlow, 
                       HZN, Bfrac, AnnualMSY, silent, PPD, checks=checks, control=control)
    
  }
  
  return(MSE1)
  
}

runMSE_int <- function(OM = DLMtool::testOM, MPs = c("AvC","DCAC","FMSYref","curE","matlenlim", "MRreal"), 
                      CheckMPs = FALSE, timelimit = 1, Hist=FALSE, ntrials=50, fracD=0.05, CalcBlow=TRUE, 
                      HZN=2, Bfrac=0.5, AnnualMSY=TRUE, silent=FALSE, PPD=FALSE, checks=FALSE,
                      control=NULL, parallel=FALSE) {
  
  # Dev Setup ####
  # development mode - assign default argument values to current workspace if they don't exist
  # def.args <- DLMtool:::dev.mode() 
  # for (nm in names(def.args)) assign(nm, def.args[[nm]])
  
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
  
  # Imp Paramerers & assign to function environment
  ImpPars <- SampleImpPars(OM, nsim, cpars=SampCpars)
  for (X in 1:length(ImpPars)) assign(names(ImpPars)[X], ImpPars[[X]])
  
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
    # Not used but make the arrays anyway
    retAp <- array(retA[,,1], dim=c(dim(retA)[1:2], Nyrs))
    Vp <- array(V[,,1], dim=c(dim(V)[1:2], Nyrs))
    noMPA <- matrix(1, nrow=Nyrs, ncol=nareas)
    
    runProj <- lapply(1:nsim, projectEq, Asize, nareas=nareas, maxage=maxage, N=N, pyears=Nyrs,
           M_ageArray=M_ageArrayp, Mat_age=Mat_agep, Wt_age=Wt_agep, V=Vp, retA=retAp,
           Perr=Perr_yp, mov=mov, SRrel=SRrel, Find=Find, Spat_targ=Spat_targ, hs=hs,
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
  N_a <-  array(NA, dim = c(nsim, maxage, nyears+proyears, nareas))
  Biomass_a <-  array(NA, dim = c(nsim, maxage, nyears+proyears, nareas))
  SSB_a <-  array(NA, dim = c(nsim, maxage, nyears+proyears, nareas))
  
  SSN_a[SAYR_a] <- Nfrac[SAY_a] * R0[S_a] * initdist[SAR_a]  # Calculate initial spawning stock numbers for all years
  N_a[SAYR_a] <- R0[S_a] * surv[SAY_a] * initdist[SAR_a] # Calculate initial stock numbers for all years
  Biomass_a[SAYR_a] <- N_a[SAYR_a] * Wt_age[SAY_a]  # Calculate initial stock biomass
  SSB_a[SAYR_a] <- SSN_a[SAYR_a] * Wt_age[SAY_a]    # Calculate spawning stock biomass

  SSN0_a <- apply(SSN_a, c(1,3), sum) # unfished spawning numbers for each year
  N0_a <- apply(N_a, c(1,3), sum) # unfished numbers for each year)
  SSB0_a <- apply(SSB_a, c(1,3), sum) # unfished spawning biomass for each year
  SSB0a <- apply(SSB_a, c(1, 3,4), sum)  # Calculate unfished spawning stock biomass by area for each year
  B0_a <- apply(Biomass_a, c(1,3), sum) # unfished biomass for each year
  
  # ---- Unfished Reference Points ----
  # Calculated average 10 years around current year
  yr.ind <- (nyears-5):(nyears+4)
  if (max(yr.ind) > nyears+proyears) {
    t.yrs <- 9-proyears
    yr.ind <- (nyears-8):(nyears+proyears)
  }
  
  N0 <- apply(N0_a[,yr.ind], 1, mean) # average unfished numbers
  SSB0 <- apply(SSB0_a[,yr.ind], 1, mean) # average unfished spawning biomass
  B0 <- apply(B0_a[,yr.ind], 1, mean) # average unfished biomass

  SSBpRa <- array(SSB0_a/matrix(R0, nrow=nsim, ncol=nyears+proyears), dim = c(nsim, nyears+proyears))
  SSBpR <- matrix(apply(SSBpRa[, yr.ind], 1, mean), nrow=nsim, ncol=nareas) # average spawning stock biomass per recruit 
  
  bR <- matrix(log(5 * hs)/(0.8 * apply(SSB0a[,1:10,], c(1,3), mean)), nrow=nsim)  # Ricker SR params
  aR <- matrix(exp(bR * apply(SSB0a[,1:10,], c(1,3), mean))/SSBpR, nrow=nsim)  # Ricker SR params
  
  # --- Optimize for Initial Depletion ----
  # Depletion in year 1 
  initD <- control$initD # not implemented yet - Iref calculation won't work if 
  # stock isn't unfished at beginning of time-series
  # initD <- runif(nsim, 0.01, 1.5) # initial depletion SB/SB0 in year 1 
  # Perr_y <- StockPars$Perr_y
  
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
              Vuln=V[x,,], Retc=retA[x,,], Prec=Perr_y[x,], movcy=split.along.dim(mov[x,,,,],4), 
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
  
  if (nsim > 1) Depletion <- apply(SSB[,,nyears,],1,sum)/SSB0
  if (nsim == 1) Depletion <- sum(SSB[,,nyears,])/SSB0 

  # Check that depletion is correct
  if (checks) {
    if (prod(round(D, 2)/ round(Depletion,2)) != 1) {
      print(cbind(round(D,2), round(Depletion,2)))
      warning("Possible problem in depletion calculations")
    } 
  } 
  
  # --- Calculate MSY references ----  
  if(!silent) message("Calculating MSY reference points") 
  
  # MSY reference points are calculated from life-history parameters averaged
  # over 10 years (around current year - i.e last historical)
  
  MSYrefs <- sapply(1:nsim, optMSY_eq, M_ageArray, Wt_age, Mat_age, V, maxage, 
                    R0, SRrel, hs, yr.ind=yr.ind) 
  
  MSY <- MSYrefs[1, ]  # record the MSY results (Vulnerable)
  FMSY <- MSYrefs[2, ]  # instantaneous FMSY (Vulnerable)
  SSBMSY <- MSYrefs[3, ]  # Spawning Stock Biomass at MSY
  SSBMSY_SSB0 <- MSYrefs[4, ] # SSBMSY relative to unfished (SSB)
  BMSY_B0 <- MSYrefs[5, ] # Biomass relative to unfished (B0)
  BMSY <- MSYrefs[6,] # total biomass at MSY
  VBMSY <- MSYrefs[7,] # Biomass at MSY (Vulnerable)
  UMSY <- MSY/VBMSY  # exploitation rate [equivalent to 1-exp(-FMSY)]
  FMSY_M <- FMSY/M  # ratio of true FMSY to natural mortality rate M

  if (checks) {
    Btemp <- apply(SSB, c(1,3), sum)
    x <- Btemp[,nyears]/SSBMSY
    y <-D/SSBMSY_SSB0
    plot(x,y, xlim=c(0,max(x)), ylim=c(0,max(y)), xlab="SSB/SSBMSY", ylab="D/SSBMSY_SSB0")
    lines(c(-10,10),c(-10,10))
  }
  
  # --- Calculate B-low ---- 
  # (SSB where it takes MGThorizon x MGT to reach Bfrac of BMSY)
  Znow<-apply(Z[,,nyears,]*N[,,nyears,],1:2,sum)/apply(N[,,nyears,],1:2,sum)
  MGTsurv<-t(exp(-apply(Znow,1,cumsum)))
  MGT<-apply(Agearray*(Mat_age[,,nyears]*MGTsurv),1,sum)/apply(Mat_age[,,nyears]*MGTsurv,1,sum)
  
  Blow <- rep(NA,nsim)
  if(CalcBlow){
    if(!silent) message("Calculating B-low reference points")            
    MGThorizon<-floor(HZN*MGT)
    Blow <- sapply(1:nsim,getBlow, N, Asize, MSYrefs[3,],SSBpR, MPA, SSB0, nareas, retA, MGThorizon,
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
                         BMSY_B0=BMSY_B0, BMSY=BMSY, VBMSY=VBMSY, UMSY=UMSY,
                         FMSY_M=FMSY_M, N0=N0, SSB0=SSB0, B0=B0, RefY=RefY, 
                         Blow=Blow, MGT=MGT, R0=R0)
  
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
                               ErrList, OM, control=NULL, silent=FALSE)
  
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
                           Maturity=Mat_age, N.Mortality=M_ageArray,
                           N=apply(N, 1:3, sum),
                           SSB=apply(SSB, 1:3, sum)
                           )
    nout <- t(apply(N, c(1, 3), sum)) 
    vb <- t(apply(VBiomass, c(1, 3), sum))
    b <- t(apply(Biomass, c(1, 3), sum))
    ssb <- t(apply(SSB, c(1, 3), sum))
    Cc <- t(apply(CB, c(1,3), sum))
    Ccret <- t(apply(CBret, c(1,3), sum))
    rec <- t(apply((N)[, 1, , ], c(1,2), sum))
    TSdata <- list(VB=vb, SSB=ssb, B=b, Removals=Cc, Catch=Ccret, Rec=rec, N=nout,
                   Find=t(Find))
    HistObj@TSdata <- TSdata
    HistObj@Ref <- RefPoints[,order(colnames(RefPoints))]
    return(HistObj)	
  }
  
  # --- Check MPs ---- 
  ## TURN THIS INTO A FUNCTION #####
  if (is.na(MPs[1])) CheckMPs <- TRUE
  if (CheckMPs) {
    if(!silent) message("Determining available methods") 
    PosMPs <- Can(Data, timelimit = timelimit)  # list all the methods that could be applied to a DLM data object 
    if (is.na(MPs[1])) {
      MPs <- PosMPs  # if the user does not supply an argument MPs run the MSE for all available methods
      if(!silent) message("No MPs specified: running all available")
    }
    if (!is.na(MPs[1])) {
      cant <- MPs[!MPs %in% PosMPs]
      if (length(cant) > 0) {
        if(!silent) message("Cannot run some MPs: ")
        if(!silent) print(DLMdiag(Data, "not available", funcs1=cant, timelimit = timelimit))
      }
      MPs <- MPs[MPs %in% PosMPs]  # otherwise run the MSE for all methods that are deemed possible
    }
    if (length(MPs) == 0) {
      if(!silent) message(Cant(Data, timelimit = timelimit))
      stop("MSE stopped: no viable methods \n\n")  # if none of the user specified methods are possible stop the run
    }
  }
  
  ok <- rep(TRUE, length(MPs))
  for (mm in seq_along(MPs)) {
    test <- try(get(MPs[mm]), silent=TRUE)
    if (!class(test) == 'MP') {
      ok[mm] <- FALSE
      if (class(test) == 'try-error') {
        message('Object ', paste(MPs[mm], ""), " does not exist - Ignoring")
      } else message('Dropping MP: ', paste(MPs[mm], ""), " - Not class 'MP'")
    }
  }

  MPs <- MPs[ok]
  nMP <- length(MPs)  # the total number of methods used
  
  if (nMP < 1) stop("No valid MPs found")
  
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
  
  MSY_P <- array(MSY, dim=c(nsim, nMP, proyears)) # store MSY for each sim, MP and year
  FMSY_P <- array(FMSY, dim=c(nsim, nMP, proyears)) # store FMSY for each sim, MP and year
  SSBMSY_P <- array(SSBMSY, dim=c(nsim, nMP, proyears))# store SSBMSY for each sim, MP and year 
  
  # --- Calculate MSY statistics for each projection year ----
  if (AnnualMSY) {
    if(!silent) message("Calculating MSY reference points for each projection year")
    # average life-history parameters over 10 years
    for (y in (nyears+1):(nyears+proyears)) {
      if(!silent) cat('.')
      y2 <- y - nyears 
      yr.ind <- (y-5):(y+4)
      if (max(yr.ind) > nyears+proyears) {
        yr.ind <- nyears + (proyears-9):(proyears)
      }
      MSYrefsYr <- sapply(1:nsim, optMSY_eq, M_ageArray, Wt_age, Mat_age, V, maxage, R0, SRrel, hs, yr.ind=yr.ind)
      MSY_P[,,y2] <- MSYrefsYr[1, ]
      FMSY_P[,,y2] <- MSYrefsYr[2,]
      SSBMSY_P[,,y2] <- MSYrefsYr[3,]
    }
    if(!silent) cat("\n")
  }

  # --- Begin loop over MPs ----
  mm <- 1 # for debugging

  for (mm in 1:nMP) {  # MSE Loop over methods
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
    # FM_nospace <- array(NA, dim = c(nsim, maxage, proyears, nareas))  # stores prospective F before reallocation to new areas
    # FML <- array(NA, dim = c(nsim, nareas))  # last apical F
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
    
    # Update abundance estimates - used for FMSY ref methods 
    M_array <- array(0.5*M_ageArray[,,nyears+y], dim=c(nsim, maxage, nareas))
    Atemp <- apply(VBiomass_P[, , y, ] * exp(-M_array), 1, sum) # Abundance (mid-year before fishing)
    MSElist[[mm]]@OM$A <- Atemp 

    # -- Apply MP in initial projection year ----
    runMP <- applyMP(MSElist[[mm]], MPs = MPs[mm], reps = reps, silent=TRUE)  # Apply MP
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
    
    # -- Calc stock dynamics ----
    MPCalcs <- CalcMPDynamics(MPRecs, y, nyears, proyears, nsim, Biomass_P, VBiomass_P,
                              LastEi, LastSpatial, LastAllocat, LastCatch,
                              TACused, maxF,
                              LR5_P, LFR_P, Rmaxlen_P, retL_P, retA_P,
                              L5_P, LFS_P, Vmaxlen_P, SLarray_P, V_P,
                              Fdisc_P, DR_P,
                              M_ageArray, FM_P, FM_Pret, Z_P, CB_P, CB_Pret,
                              TAC_f, E_f, SizeLim_f,
                              FinF, Spat_targ,
                              CAL_binsmid, Linf, Len_age, maxage, nareas, Asize, nCALbins,
                              qs, qvar, qinc)
  
    TACa[, mm, y] <- MPCalcs$TACrec # recommended TAC 
    LastSpatial <- MPCalcs$Si
    LastAllocat <- MPCalcs$Ai
    LastEi <- MPCalcs$Ei # adjustment to effort
    LastCatch <- MPCalcs$TACrec
    
    Effort[, mm, y] <- MPCalcs$Effort  
    CB_P <- MPCalcs$CB_P # removals
    CB_Pret <- MPCalcs$CB_Pret # retained catch 
    FM_P <- MPCalcs$FM_P # fishing mortality
    FM_Pret <- MPCalcs$FM_Pret # retained fishing mortality 
    Z_P <- MPCalcs$Z_P # total mortality
    
    retA_P <- MPCalcs$retA_P # retained-at-age
    
    retL_P <- MPCalcs$retL_P # retained-at-length
    V_P <- MPCalcs$V_P  # vulnerable-at-age
    SLarray_P <- MPCalcs$SLarray_P # vulnerable-at-length
    FMa[,mm,y] <- MPCalcs$Ftot 
    
    # --- Begin projection years ----
    for (y in 2:proyears) {
      if(!silent) {
        cat("."); flush.console()
      }
      
      SelectChanged <- FALSE
      if (AnnualMSY) {
        if (any(range(retA_P[,,nyears+y] - retA[,,nyears+y]) !=0)) SelectChanged <- TRUE
        if (any(range(V_P[,,nyears+y] - V[,,nyears+y]) !=0))  SelectChanged <- TRUE
      }
      
      # -- Calculate MSY stats for this year ----
      if (AnnualMSY & SelectChanged) { #
        # average 10 year life-history & selecivity used for MSY calcs
        y1 <- nyears + y
        y2 <- y  
        yr.ind <- (y1-5):(y1+4)
        if (max(yr.ind) > nyears+proyears) {
          yr.ind <- nyears + (proyears-9):(proyears)
        }
        
        MSYrefsYr <- sapply(1:nsim, optMSY_eq, M_ageArray, Wt_age, Mat_age, 
                            V_P, maxage, R0, SRrel, hs, yr.ind=yr.ind)
        MSY_P[,mm,y] <- MSYrefsYr[1, ]
        FMSY_P[,mm,y] <- MSYrefsYr[2,]
        SSBMSY_P[,mm,y] <- MSYrefsYr[3,]
        
      }
      
      TACa[, mm, y] <- TACa[, mm, y-1] # TAC same as last year unless changed 
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
      
      # --- Age & Growth ----
      NextYrN <- lapply(1:nsim, function(x)
        popdynOneTScpp(nareas, maxage, SSBcurr=colSums(SSB_P[x,,y-1, ]), Ncurr=N_P[x,,y-1,],
                       Zcurr=Z_P[x,,y-1,], PerrYr=Perr_y[x, y+nyears+maxage-1], hs=hs[x],
                       R0a=R0a[x,], SSBpR=SSBpR[x,], aR=aR[x,], bR=bR[x,],
                       mov=mov[x,,,, nyears+y], SRrel=SRrel[x]))
  
      N_P[,,y,] <- aperm(array(unlist(NextYrN), dim=c(maxage, nareas, nsim, 1)), c(3,1,4,2)) 
      Biomass_P[SAYR] <- N_P[SAYR] * Wt_age[SAYt]  # Calculate biomass
      VBiomass_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt]  # Calculate vulnerable biomass
      SSN_P[SAYR] <- N_P[SAYR] * Mat_age[SAYt]  # Calculate spawning stock numbers
      SSB_P[SAYR] <- SSN_P[SAYR] * Wt_age[SAYt]  # Calculate spawning stock biomass
    
      # --- An update year ----
      if (y %in% upyrs) {
        # --- Update Data object ---- 
        MSElist[[mm]] <- updateData(Data=MSElist[[mm]], OM, MPCalcs, Effort, Biomass, 
                                    Biomass_P, CB_Pret, N_P, SSB_P, VBiomass_P, 
                                    RefPoints, ErrList, FMSY_P, retA_P, retL_P, StockPars, 
                                    FleetPars, ObsPars, upyrs, interval, y, mm, 
                                    Misc=Data_p@Misc)

        # Update Abundance and FMSY for FMSYref MPs
        M_array <- array(0.5*M_ageArray[,,nyears+y], dim=c(nsim, maxage, nareas))
        Atemp <- apply(VBiomass_P[, , y, ] * exp(-M_array), 1, sum) # Abundance (mid-year before fishing)
        MSElist[[mm]]@OM$A <- Atemp
        MSElist[[mm]]@OM$FMSY <- FMSY_P[,mm,y]

        # --- apply MP ----
        runMP <- applyMP(MSElist[[mm]], MPs = MPs[mm], reps = reps, silent=TRUE)  # Apply MP
        MPRecs <- runMP[[1]][[1]] # MP recommendations
        Data_p <- runMP[[2]] # Data object object with saved info from MP 
        Data_p@TAC <- MPRecs$TAC
        # calculate pstar quantile of TAC recommendation dist 
        TACused <- apply(Data_p@TAC, 2, quantile, p = pstar, na.rm = T) 
        
        # -- Calc stock dynamics ----
        MPCalcs <- CalcMPDynamics(MPRecs, y, nyears, proyears, nsim, Biomass_P, VBiomass_P,
                                  LastEi, LastSpatial, LastAllocat, LastCatch,
                                  TACused, maxF,
                                  LR5_P, LFR_P, Rmaxlen_P, retL_P, retA_P,
                                  L5_P, LFS_P, Vmaxlen_P, SLarray_P, V_P,
                                  Fdisc_P, DR_P,
                                  M_ageArray, FM_P, FM_Pret, Z_P, CB_P, CB_Pret,
                                  TAC_f, E_f, SizeLim_f,
                                  FinF, Spat_targ,
                                  CAL_binsmid, Linf, Len_age, maxage, nareas, Asize, nCALbins,
                                  qs, qvar, qinc)
        
        TACa[, mm, y] <- TACused # recommended TAC 
        LastSpatial <- MPCalcs$Si
        LastAllocat <- MPCalcs$Ai
        LastEi <- MPCalcs$Ei # adjustment to effort
        Effort[, mm, y] <- MPCalcs$Effort 
        FMa[,mm,y] <- MPCalcs$Ftot 
      
        CB_P <- MPCalcs$CB_P # removals
        CB_Pret <- MPCalcs$CB_Pret # retained catch 
        LastCatch <- apply(CB_Pret[,,y,], 1, sum, na.rm=TRUE) 
        FM_P <- MPCalcs$FM_P # fishing mortality
        FM_Pret <- MPCalcs$FM_Pret # retained fishing mortality 
        Z_P <- MPCalcs$Z_P # total mortality
        retA_P <- MPCalcs$retA_P # retained-at-age
        retL_P <- MPCalcs$retL_P # retained-at-length
        V_P <- MPCalcs$V_P  # vulnerable-at-age
        SLarray_P <- MPCalcs$SLarray_P # vulnerable-at-length
        
      } else {
        # --- Not an update yr ----
        NoMPRecs <- MPRecs # TAC & TAE stay the same
        NoMPRecs[lapply(NoMPRecs, length) > 0 ] <- NULL
        NoMPRecs$Spatial <- NA
        MPCalcs <- CalcMPDynamics(NoMPRecs, y, nyears, proyears, nsim, Biomass_P, VBiomass_P,
                                  LastEi, LastSpatial, LastAllocat, LastCatch,
                                  TACused, maxF,
                                  LR5_P, LFR_P, Rmaxlen_P, retL_P, retA_P,
                                  L5_P, LFS_P, Vmaxlen_P, SLarray_P, V_P,
                                  Fdisc_P, DR_P,
                                  M_ageArray, FM_P, FM_Pret, Z_P, CB_P, CB_Pret,
                                  TAC_f, E_f, SizeLim_f,
                                  FinF, Spat_targ,
                                  CAL_binsmid, Linf, Len_age, maxage, nareas, Asize, nCALbins,
                                  qs, qvar, qinc)
      
        TACa[, mm, y] <- TACused # 
        LastSpatial <- MPCalcs$Si
        LastAllocat <- MPCalcs$Ai
        LastEi <- MPCalcs$Ei
        Effort[, mm, y] <- MPCalcs$Effort  
        CB_P <- MPCalcs$CB_P # removals
        CB_Pret <- MPCalcs$CB_Pret # retained catch 
        FMa[,mm,y] <- MPCalcs$Ftot 
        LastCatch <- apply(CB_Pret[,,y,], 1, sum, na.rm=TRUE) 
        FM_P <- MPCalcs$FM_P # fishing mortality
        FM_Pret <- MPCalcs$FM_Pret # retained fishing mortality 
        Z_P <- MPCalcs$Z_P # total mortality
        retA_P <- MPCalcs$retA_P # retained-at-age
        retL_P <- MPCalcs$retL_P # retained-at-length
        V_P <- MPCalcs$V_P  # vulnerable-at-age
        SLarray_P <- MPCalcs$SLarray_P # vulnerable-at-length
      } # end of update loop 
      checkNA[y] <- sum(is.na(TACused))
    }  # end of year loop
    
    B_BMSYa[, mm, ] <- apply(SSB_P, c(1, 3), sum, na.rm=TRUE)/SSBMSY_P[,mm,]  # SSB relative to SSBMSY
    F_FMSYa[, mm, ] <- FMa[, mm, ]/FMSY_P[,mm,]
    
    Ba[, mm, ] <- apply(Biomass_P, c(1, 3), sum, na.rm=TRUE) # biomass 
    SSBa[, mm, ] <- apply(SSB_P, c(1, 3), sum, na.rm=TRUE) # spawning stock biomass
    VBa[, mm, ] <- apply(VBiomass_P, c(1, 3), sum, na.rm=TRUE) # vulnerable biomass
    
    Ca[, mm, ] <- apply(CB_P, c(1, 3), sum, na.rm=TRUE) # removed
    CaRet[, mm, ] <- apply(CB_Pret, c(1, 3), sum, na.rm=TRUE) # retained catch 
    
    # Store Pop and Catch-at-age and at-length for last projection year 
    PAAout[ , mm, ] <- apply(N_P[ , , proyears, ], c(1,2), sum) # population-at-age
    
    CNtemp <- apply(CB_Pret, c(1,2,3), sum)/Wt_age[(nyears+1):nyears+proyears]
    CAAout[ , mm, ] <- CNtemp[,,proyears] # nsim, maxage # catch-at-age
    CALdat <- MSElist[[mm]]@CAL
    CALout[ , mm, ] <- CALdat[,dim(CALdat)[2],] # catch-at-length in last year
    
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
    
  }  # end of mm methods 
  
  # Miscellaneous reporting
  if(PPD)Misc<-MSElist

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
    if(class(x)=="numeric" | class(x)=="integer")length(x)
    else dim(x)[1]
  }
  
  dims <- sapply(cpars,dim1check)
  
  # drop invalid names 
  dims <- dims[!unlist(lapply(dims, is.null))]
  
  # Check Stock names 
  slts <- slotNames("Stock")
  ind <- which(names(cpars) %in% slts)
  if (length(ind)>0) {
    chkdim <- !unlist(lapply(lapply(cpars[ind], dim), is.null))
    if (any(chkdim)) {
      err <- names(chkdim[chkdim])
      stop("Incorrect dimensions in Stock parameter slots in OM@cpars:", paste(err, collapse=", "))
    }
  }
  

  # check if EffYears etc are in cpars
  effNames <- c("EffYears", "EffLower", "EffUpper")
  temp <- effNames %in%  names(dims)
  # all 3?
  if (any(temp) & !all(temp)) stop(paste(effNames[!temp], collapse=", "), " missing")
  # can't have Find as well
  if (any(temp & 'Find' %in% names(dims))) stop("Can't provide both Find and EffYears etc in cpars")
  # same length 
  if (all(temp)) {
    if (!all(dims[effNames]==dims[effNames][1])) stop(paste(effNames, collapse=", "), " are not equal length")
  }

  # ignore 
  if (any(effNames %in% names(dims))) dims <- dims[-match(effNames,names(dims))]  # ignore effNames
  dims <- dims[!grepl("CAL_bins", names(dims))]  # ignore CAL_bins
  dims <- dims[!grepl("maxage", names(dims))]  # ignore CAL_bins
  # dims <- dims[!grepl("Asize", names(dims))]  # ignore CAL_bins
  if (length(dims) > 0) {
    if(length(unique(dims))!=1){
      print(dims)
      stop("The custom parameters in your operating model @cpars have varying number of simulations. For each simulation each parameter / variable should correspond with one another")
    }else{
      as.integer(dims[1])
    }
  }
  
}


# cparnamecheck<-function(cpars){
# 
#   Sampnames <- c("D","Esd","Find","procsd","AC","M","Msd",
#                  "Mgrad","hs","Linf","Linfsd","Linfgrad",
#                  "K","Ksd","Kgrad","t0","L50","L50_95","Spat_targ",
#                  "Frac_area_1","Prob_staying",
#                  "Csd","Cbias","CAA_nsamp","CAA_ESS","CAL_nsamp",
#                  "CAL_ESS","betas","Isd","Derr","Dbias",
#                  "Mbias","FMSY_Mbias","lenMbias","LFCbias",
#                  "LFSbias","Aerr","Abias","Kbias","t0bias",
#                  "Linfbias","Irefbias","Crefbias","Brefbias",
#                  "Recsd","qinc","qcv","L5","LFS","Vmaxlen","L5s",
#                  "LFSs","Vmaxlens","Perr","R0","Mat_age",
#                  "Mrand","Linfrand","Krand","maxage","V","Depletion", # end of OM variables
#                  "ageM", "age95", "V", "EffYears", "EffLower", "EffUpper","Mat_age", # start of runMSE derived variables
#                  "Wt_age")
# 
# }


#' @describeIn runMSE Save out the results to a Rdata
#' file.  To increase speed and efficiency, particulary for runs with a large
#' number simulations (\code{nsim}), the simulations are split into a number of
#' packets.  The functions loops over the packets and combines the output into
#' a single MSE object. If the MSE model crashes during a run, the MSE is run
#' again until it is successfully completed. The MSE is stopped if the number
#' of consecutive crashes exceeds \code{maxCrash}.  There is an option to save
#' the packets as Rdata files to the current working directory (default is
#' FALSE). By default, the functions saves the completed MSE object as a Rdata
#' file (to the current working directory).
#' @export
runMSErobust <- function(OM = DLMtool::testOM, MPs = c("AvC", "DCAC", "FMSYref", "curE", "matlenlim", "MRreal"), 
                         timelimit = 1, CheckMPs = FALSE, Hist = FALSE, 
                         ntrials = 50, fracD = 0.05, CalcBlow = FALSE, HZN = 2, Bfrac = 0.5, AnnualMSY=TRUE,
                         maxsims = 64, name = NULL, unique=FALSE, maxCrash = 10, saveMSE = TRUE, 
                         savePack = FALSE) {
  
  if (!snowfall::sfIsRunning()) {
    message("Setting up parallel processing")
    setup()
  }
  
  if (class(OM) != "OM") stop("You must specify an operating model")
  cpars <- NULL
  ncparsim <- nsim
  if(length(OM@cpars)>0){
    # ncparsim<-cparscheck(OM@cpars)   # check each list object has the same length and if not stop and error report
    cpars <- OM@cpars
  }
  
  # Backwards compatible with DLMtool v < 4
  if("nsim"%in%slotNames(OM))nsim<-OM@nsim
  if("proyears"%in%slotNames(OM))proyears<-OM@proyears
  
  OM@proyears<-proyears
  
  packets <- new("list")  # a list of completed MSE objects
  simsplit <- split(1:nsim, ceiling(seq_along(1:nsim)/maxsims))  # split the runs
  message(nsim, " simulations \n")
  message("Splitting into ", length(simsplit), " packets")
  flush.console()  
  
  # set name for saved MSE object 
  if (is.null(name)) {
    st <- as.numeric(regexpr(":", OM@Name)) + 1
    nd <- st + 4  # as.numeric(regexpr(' ', OM@Name))-1
    name <- substr(OM@Name, st, nd)
    name <- gsub("\\s+", "", name)
  }
  if (nchar(name) < 1)  name <- "MSE"
  if (unique) name <- paste0(name, "_", format(Sys.time(), "%H%M_%m%d%y"))
  
  stElap <- rep(NA, length(simsplit)) # store elapsed time 
  
  index <- 1:nsim
  
  
  for (i in 1:length(simsplit)) {
    message("Packet ", i, " of ", length(simsplit), " started")
    error <- 1
    crash <- 0
    st <- Sys.time()
    while (error == 1 & crash <= maxCrash) {
      assign("last.warning", NULL, envir = baseenv())
      tryOM <- OM
      tryOM@seed <- OM@seed + i + crash  # change seed 
      tryOM@nsim <- length(simsplit[[i]]) # sub number of sims 
      if (length(cpars) > 0) tryOM@cpars  <- SampleCpars(cpars, nsim=tryOM@nsim, msg=FALSE)
      
      trialMSE <- try(runMSE(OM = tryOM, MPs = MPs, timelimit = timelimit, CheckMPs = CheckMPs, 
                             Hist=Hist, ntrials=ntrials, fracD=fracD, CalcBlow = CalcBlow, HZN = HZN, 
                             Bfrac = Bfrac, AnnualMSY=AnnualMSY, parallel=TRUE))	
      
      if (!Hist & class(trialMSE) != "MSE") {
        crash <- crash + 1
        print(warnings())
        message("Packet ", i, " crashed. Trying again\n")
      }
      if (Hist & class(trialMSE) != "list") {
        crash <- crash + 1
        print(warnings())
        message("Packet ", i, " crashed. Trying again\n")
      }     
      if (crash >= maxCrash) stop("\nNumber of crashes exceeded 'maxCrash'\n", call.=FALSE)
      
      if (class(trialMSE) == "MSE" | class(trialMSE) == "list") {
        packets[[i]] <- trialMSE
        fname <- paste0(name, "_P", i, ".rdata")
        if (savePack) {
          saveRDS(trialMSE, file = fname)
          message("Saving" , fname, " to ", getwd())
          flush.console()
        }	
        error <- 0
        crash <- 0
      }
    }
    elapse <- Sys.time() - st 
    stElap[i] <- elapse
    message("Packet ", i, " of ", length(simsplit), " complete\n")
    flush.console()
    eta <- round(mean(stElap, na.rm=TRUE) * length(simsplit) - sum(stElap, na.rm=TRUE), 2)
    units <- attributes(elapse)$units
    if (eta > 120 && units == "secs") {
      eta <- round(eta/60,2)
      units <- "mins"
    }
    if (i != length(simsplit)) message("\nEstimated time to completion is: ", eta, " ", units)
    flush.console()
    
  }
  if (i == 1) MSEobj <- packets[[1]]
  
  if (i > 1 & ! Hist) MSEobj <- joinMSE(MSEobjs = packets)
  if (i > 1 & Hist) MSEobj <- unlist(packets,F)
  if (saveMSE) {
    fname <- paste0(name, ".rdata")
    saveRDS(MSEobj, file = fname)
    message("Saving ", fname, " to ", getwd())
    flush.console()
  }
  MSEobj
}

