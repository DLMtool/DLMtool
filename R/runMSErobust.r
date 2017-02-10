#' Run a Management Strategy Evaluation
#' 
#' Run a Management Strategy Evaluation and save out the results to a Rdata
#' file.  To increase speed and efficiency, particulary for runs with a large
#' number simulations (\code{nsim}), the simulations are split into a number of
#' packets.  The functions loops over the packets and combines the output into
#' a single MSE object. If the MSE model crashes during a run, the MSE is run
#' again until it is successfully completed. The MSE is stopped if the number
#' of consecutive crashes exceeds \code{maxCrash}.  There is an ption to save
#' the packets as Rdata files to the current working directory (default is
#' FALSE). By default, the functions saves the completed MSE object as a Rdata
#' file (to the current working directory).
#' 
#' 
#' @usage runMSErobust(OM = '1', MPs = NA, nsim = 200, proyears = 28, interval
#' = 4, pstar = 0.5, maxF = 0.8, timelimit = 1, reps = 1, custompars = 0,
#' CheckMPs = TRUE, maxsims = 64, name = NULL, unique=TRUE, maxCrash = 10, saveMSE = TRUE,
#' savePack = FALSE)
#' @param OM An operating model object (class OM)
#' @param MPs A vector of methods (character string) of class DLM_output or
#' DLM_input. If NA all available MPs are run.
#' @param nsim Number of simulations
#' @param proyears Number of projected years
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
#' @param custompars A data.table with nsim rows and nparameter columns. The
#' column names must respond to variables of the operating model or observation
#' model see the OM and Obs slots of the MSE class for correct names and
#' interpretation. This allows users to prescribe correlated parameters or
#' estimates from stock assessments.
#' @param CheckMPs Logical to indicate if Can function should be used to check
#' if MPs can be run
#' @param Hist Should model stop after historical simulations? Returns a list 
#' containing all historical data
#' @param ntrials Maximum of times depletion and recruitment deviations are 
#' resampled to optimize for depletion. After this the model stops if more than 
#' percent of simulations are not close to the required depletion
#' @param fracD maximum allowed proportion of simulations where depletion is not 
#' close to sampled depletion from OM before model stops with error
#' @param maxsims Maximum number of simulations per packet
#' @param name Character string for name of saved MSE packets (if
#' \code{savePack=TRUE}) and final MSE object. If none provided, it uses the
#' first four letters from the \code{OM} name
#' @param unique Logical. Should the name be unique? Current date and time appended to name. 
#' @param maxCrash Maximum number of consecutive crashes before the MSE stops
#' @param saveMSE Logical to indiciate if final MSEobject should be saved to
#' current working directory (this is probably a good idea)
#' @param savePack Logical to indicate if packets should be save to current
#' working directory
#' @return An object of class MSE
#' @author A. Hordyk and T. Carruthers
#' @export runMSErobust
runMSErobust <- function(OM = "1", MPs = NA, nsim = 200, proyears = 28, 
  interval = 4, pstar = 0.5, maxF = 0.8, timelimit = 1, reps = 1, custompars = NULL, 
  CheckMPs = TRUE, Hist=FALSE, ntrials=50, fracD=0.05,
  maxsims = 64, name = NULL, unique=TRUE, maxCrash = 10, saveMSE = TRUE, 
  savePack = FALSE) {
  
  packets <- new("list")  # a list of completed MSE objects
  simsplit <- split(1:nsim, ceiling(seq_along(1:nsim)/maxsims))  # split the runs
  
  if (!snowfall::sfIsRunning()) stop("You must set up parallel processing to use runMSErobust. Use: setup()", call.=FALSE)
  
  if (is.null(name)) {
    st <- as.numeric(regexpr(":", OM@Name)) + 1
    nd <- st + 3  # as.numeric(regexpr(' ', OM@Name))-1
    name <- substr(OM@Name, st, nd)
    name <- gsub("\\s+", "", name)
  }
  if (nchar(name) < 1)  name <- paste0("MSE_")
  
  if (unique) name <- paste0(name, "_", format(Sys.time(), "%H%M_%m%d%y"))
  
  for (i in 1:length(simsplit)) {
    error <- 1
    crash <- 0
    while (error == 1 & crash <= maxCrash) {
      trialMSE <- try(runMSE(OM = OM, MPs = MPs, nsim = length(simsplit[[i]]), 
        proyears = proyears, interval = interval, pstar = pstar, 
        maxF = maxF, timelimit = timelimit, reps = reps, custompars = custompars, 
        CheckMPs = CheckMPs, Hist, ntrial, fracD))
      crash <- crash + 1
      if (crash >= maxCrash) 
        stop("Crashed too many times!")
      if (class(trialMSE) == "MSE") {
        packets[[i]] <- trialMSE
        fname <- paste0(name, "pack", i, ".rdata")
        if (savePack) {
          saveRDS(trialMSE, file = fname)
          message(paste("Saving", fname, "to", getwd()))
          flush.console()
        }
        error <- 0
        crash <- 0
      }
    }
    message(paste("Packet", i, "of", length(simsplit), "complete"))
    flush.console()
  }
  if (i == 1) MSEobj <- packets[[1]]
  if (i > 1) MSEobj <- joinMSE(MSEobjs = packets)
  if (saveMSE) {
    fname <- paste0(name, ".rdata")
    saveRDS(MSEobj, file = fname)
    message(paste("Saving", fname, "to", getwd()))
    flush.console()
  }
  MSEobj
}
