#' Sample Observation Parameters
#'
#' @param Obs An object of class 'Obs' or class 'OM'
#' @param nsim Number of simulations. Ignored if 'Obs' is class 'OM'
#' @return A named list of sampled Observation parameters
#' @export
#'
SampleObsPars <- function(Obs, nsim=NULL){
  if (class(Obs) != "Obs" & class(Obs) != "OM") 
    stop("First argument must be class 'Obs' or 'OM'")
  if (class(Obs) == "OM") nsim <- Obs@nsim

  Obs <- updateMSE(Obs) # update to add missing slots with default values

  ObsOut <- list() 
  
  # === Sample observation error model parameters ====
  
  # fix some naming issues?
  ObsOut$Csd <- runif(nsim, Obs@Cobs[1], Obs@Cobs[2])  # Sampled catch observation error (lognormal sd)
  ObsOut$Cbias <- rlnorm(nsim, mconv(1, Obs@Cbiascv), sdconv(1, Obs@Cbiascv))  # Sampled catch bias (log normal sd)
  ObsOut$CAA_nsamp <- ceiling(runif(nsim, Obs@CAA_nsamp[1], Obs@CAA_nsamp[2]))  # Number of catch-at-age observations
  ObsOut$CAA_ESS <- ceiling(runif(nsim, Obs@CAA_ESS[1], Obs@CAA_ESS[2]))  # Effective sample size
  ObsOut$CAL_nsamp <- ceiling(runif(nsim, Obs@CAL_nsamp[1], Obs@CAL_nsamp[2]))  # Observation error standard deviation for single catch at age by area
  ObsOut$CAL_ESS <- ceiling(runif(nsim, Obs@CAL_ESS[1], Obs@CAL_ESS[2]))  # Effective sample size
  # ObsOut$CALcv <- runif(nsim, Obs@CALcv[1], Obs@CALcv[2])  # Observation error standard deviation for single catch at age by area
  ObsOut$betas <- exp(runif(nsim, log(Obs@beta[1]), log(Obs@beta[2])))  # the sampled hyperstability / hyperdepletion parameter beta>1 (hyperdepletion) beta<1 (hyperstability)
  ObsOut$Isd <- runif(nsim, Obs@Iobs[1], Obs@Iobs[2])  # Abundance index observation error (log normal sd)
  ObsOut$Derr <- runif(nsim, Obs@Dobs[1], Obs@Dobs[2])
  ObsOut$Dbias <- rlnorm(nsim, mconv(1, Obs@Dbiascv), sdconv(1, Obs@Dbiascv))  # sample of depletion bias
  ObsOut$Mbias <- rlnorm(nsim, mconv(1, Obs@Mbiascv), sdconv(1, Obs@Mbiascv))  # sample of M bias
  ObsOut$FMSY_Mbias <- rlnorm(nsim, mconv(1, Obs@FMSY_Mbiascv), sdconv(1, Obs@FMSY_Mbiascv))  # sample of FMSY/M bias
 
  ObsOut$lenMbias <- rlnorm(nsim, mconv(1, Obs@LenMbiascv), sdconv(1, Obs@LenMbiascv))  # sample of length at maturity bias - assume same error as age based maturity
  ObsOut$LFCbias <- rlnorm(nsim, mconv(1, Obs@LFCbiascv), sdconv(1, Obs@LFCbiascv))  # sample of length at first capture bias
  ObsOut$LFSbias <- rlnorm(nsim, mconv(1, Obs@LFSbiascv), sdconv(1, Obs@LFSbiascv))  # sample of length at full selection bias
  ObsOut$Aerr <- runif(nsim, Obs@Btobs[1], Obs@Btobs[2])
  ObsOut$Abias <- exp(runif(nsim, log(Obs@Btbiascv[1]), log(Obs@Btbiascv[2])))  #rlnorm(nsim,mconv(1,Obs@Btbiascv),sdconv(1,Obs@Btbiascv))    # sample of current abundance bias
  ObsOut$Kbias <- rlnorm(nsim, mconv(1, Obs@Kbiascv), sdconv(1, Obs@Kbiascv))  # sample of von B. K parameter bias
  ObsOut$t0bias <- rlnorm(nsim, mconv(1, Obs@t0biascv), sdconv(1, Obs@t0biascv))  # sample of von B. t0 parameter bias
  ObsOut$Linfbias <- rlnorm(nsim, mconv(1, Obs@Linfbiascv), sdconv(1, Obs@Linfbiascv))  # sample of von B. maximum length bias
  ObsOut$Irefbias <- rlnorm(nsim, mconv(1, Obs@Irefbiascv), sdconv(1, Obs@Irefbiascv))  # sample of bias in reference (target) abundance index
  ObsOut$Crefbias <- rlnorm(nsim, mconv(1, Obs@Crefbiascv), sdconv(1, Obs@Crefbiascv))  # sample of bias in reference (target) catch index
  ObsOut$Brefbias <- rlnorm(nsim, mconv(1, Obs@Brefbiascv), sdconv(1, Obs@Brefbiascv))  # sample of bias in reference (target) biomass index
  ObsOut$Recsd <- runif(nsim, Obs@Recbiascv[1], Obs@Recbiascv[2])  # Recruitment deviation  
  # ObsOut$LenCVbias <- rlnorm(nsim, mconv(1, Obs@CALcv), sdconv(1, Obs@CALcv)) # sample of bias in assumed CV of catch-at-length
  
  
  ObsOut
}