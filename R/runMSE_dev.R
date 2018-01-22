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
#' the management system? NOTE: since DLMtool V4.5 this slot is included in the 
#' OM object which will override the value used here. This slot to be deprecated in the future.
#' @param pstar The percentile of the sample of the management recommendation
#' for each method. NOTE: since DLMtool V4.5 this slot is included in the 
#' OM object which will override the value used here. This slot to be deprecated in the future.
#' @param maxF Maximum instantaneous fishing mortality rate that may be
#' simulated for any given age class. NOTE: since DLMtool V4.5 this slot is included in the 
#' OM object which will override the value used here. This slot to be deprecated in the future.
#' @param reps Number of samples of the management recommendation for each
#' method. Note that when this is set to 1, the mean value of the data inputs
#' is used. NOTE: since DLMtool V4.5 this slot is included in the 
#' OM object which will override the value used here. This slot to be deprecated in the future.
#' @param CheckMPs Logical to indicate if Can function should be used to check
#' if MPs can be run.
#' @param timelimit Maximum time taken for a method to carry out 10 reps
#' (methods are ignored that take longer)
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
#' @param annualMSY Logical. Should MSY statistics be calculated for each projection year? 
#' May differ from MSY statistics from last historical year if there are changes in productivity
#' @param silent Should messages be printed out to the console?
#' @param PPD Logical. Should posterior predicted data be included in the MSE object Misc slot?
#' @param save Should the MSE packets be saved?
#' @return An object of class MSE
#' @author T. Carruthers and A. Hordyk
#' @export 
runMSE_fast <- function(OM = DLMtool::testOM, MPs = c("AvC","DCAC","FMSYref","curE","matlenlim", "MRreal"),nsim=48,
                   proyears=50,interval=4,pstar = 0.5, maxF = 0.8,  reps = 1, 
                   CheckMPs = FALSE, timelimit = 1, Hist=FALSE, ntrials=50, fracD=0.05, CalcBlow=FALSE, 
                   HZN=2, Bfrac=0.5, annualMSY=TRUE, silent=TRUE, PPD=FALSE, save=FALSE) {
  
  if(!snowfall::sfIsRunning()) stop("Requires parallel. Use 'setup'", call. = FALSE)

  ncpu <- sfCpus()

  nsims <- ceiling(OM@nsim / min(OM@nsim, 48))
  nits <- nsims/ncpu
  if(nits<1) nits <- ncpu/nsims
  itsim <- rep(ceiling(OM@nsim/nits), nits)
  
  if(sum(itsim) != OM@nsim) {
    itsim[length(itsim)] <- OM@nsim - sum(itsim[1:(length(itsim)-1)] )
  }
  
  temp <- sfClusterApplyLB(1:nits, run_parallel, itsim, OM, MPs, nsim,
                           proyears,interval,pstar, maxF,  reps, 
                           CheckMPs, timelimit, Hist, ntrials, fracD, CalcBlow, 
                           HZN, Bfrac, annualMSY, silent, PPD)

  if (save) saveRDS(temp, paste0('MSEList.rdata'))
  MSE1 <- joinMSE(temp) 
  sfStop()
  closeAllConnections()
  if (class(MSE1) == "MSE") {
    message("Completed")
  } else {
    message("There's a problem!")
  }
  
  return(MSE1)
  
}


run_parallel <- function(i, itsim, OM, MPs, nsim,
                         proyears,interval,pstar, maxF,  reps, 
                         CheckMPs, timelimit, Hist, ntrials, fracD, CalcBlow, 
                         HZN, Bfrac, annualMSY, silent, PPD) {
  OM@nsim <- itsim[i]
  
  OM@seed <- OM@seed + i 
  runMSE(OM, MPs,nsim,
         proyears,interval,pstar, maxF,  reps, 
         CheckMPs, timelimit, Hist, ntrials, fracD, CalcBlow, 
         HZN, Bfrac, annualMSY, silent, PPD)
  
}


#' Determine optimal number of cpus
#'
#' @param thresh Recommended n cpus is what percent of the fastest time?
#' @param plot Logical. Show the plot?
#'
#' @export
#'
#' @author A. Hordyk
profile <- function(thresh=5, plot=TRUE) {
  cpus=1:parallel::detectCores()
  time <- NA
  for (n in cpus) {
    message(n, ' of ', max(cpus))
    if (n == 1) {
      sfStop()
      st <- Sys.time()
      tt <- runMSE()
      time[n] <- difftime(Sys.time(), st, units='secs')
    } else{
      
      sfInit(parallel=TRUE, cpus=n, slaveOutfile="test.txt")
      st <- Sys.time()
      tt <- runMSE_fast(save=TRUE)
      time[n] <- difftime(Sys.time(), st, units='secs')
    }
  } 
  df <- data.frame(ncpu=cpus, time=time)
  rec <- min(which(time < min(time) * (1 + thresh/100)))
  if (plot) {
    plot(df, type='b', ylab="time (seconds)", xlab= "# cpus", bty="l", lwd=2)
    points(rec, df[rec,2], cex=2, pch=16, col="blue")
  }
  return(df)
}





