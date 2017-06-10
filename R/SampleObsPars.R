#' Sample Observation Parameters
#'
#' @param Obs An object of class 'Obs' or class 'OM'
#' @param nsim Number of simulations. Ignored if 'Stock' is class 'OM'
#'
#' @return A named list of sampled Observation parameters
#' @export
#'
SampleObsPars <- function(Obs, nsim=NULL) {
  if (class(Obs) != "Obs" & class(Obs) != "OM") 
    stop("First argument must be class 'Obs' or 'OM'")
  if (class(Obs) == "Obs") nsim <- Obs@nsim
  
  ObsOut <- list() 
  
  # === Sample observation error model parameters ====
  
  # fix some naming issues?
  ObsOut$Csd <- runif(nsim, Obs@Cobs[1], Obs@Cobs[2])  # Sampled catch observation error (lognormal sd)
  ObsOut$Cbias <- rlnorm(nsim, mconv(1, Obs@Cbiascv), sdconv(1, Obs@Cbiascv))  # Sampled catch bias (log normal sd)
  ObsOut$CAA_nsamp <- ceiling(runif(nsim, Obs@CAA_nsamp[1], Obs@CAA_nsamp[2]))  # Number of catch-at-age observations
  ObsOut$CAA_ESS <- ceiling(runif(nsim, Obs@CAA_ESS[1], Obs@CAA_ESS[2]))  # Effective sample size
  ObsOut$CAL_nsamp <- runif(nsim, Obs@CAL_nsamp[1], Obs@CAL_nsamp[2])  # Observation error standard deviation for single catch at age by area
  ObsOut$CAL_ESS <- ceiling(runif(nsim, Obs@CAL_ESS[1], Obs@CAL_ESS[2]))  # Effective sample size
  ObsOut$CALcv <- runif(nsim, Obs@CALcv[1], Obs@CALcv[2])  # Observation error standard deviation for single catch at age by area
  ObsOut$betas <- exp(runif(nsim, log(Obs@beta[1]), log(Obs@beta[2])))  # the sampled hyperstability / hyperdepletion parameter beta>1 (hyperdepletion) beta<1 (hyperstability)
  ObsOut$Isd <- runif(nsim, Obs@Iobs[1], Obs@Iobs[2])  # Abundance index observation error (log normal sd)
  ObsOut$Derr <- runif(nsim, Obs@Dcv[1], Obs@Dcv[2])
  ObsOut$Dbias <- rlnorm(nsim, mconv(1, Obs@Dbiascv), sdconv(1, Obs@Dbiascv))  # sample of depletion bias
  ObsOut$Mbias <- rlnorm(nsim, mconv(1, Obs@Mcv), sdconv(1, Obs@Mcv))  # sample of M bias
  ObsOut$FMSY_Mbias <- rlnorm(nsim, mconv(1, Obs@FMSY_Mcv), sdconv(1, Obs@FMSY_Mcv))  # sample of FMSY/M bias
 
  ObsOut$lenMbias <- rlnorm(nsim, mconv(1, Obs@LenMcv), sdconv(1, Obs@LenMcv))  # sample of length at maturity bias - assume same error as age based maturity
  ObsOut$LFCbias <- rlnorm(nsim, mconv(1, Obs@LFCcv), sdconv(1, Obs@LFCcv))  # sample of length at first capture bias
  ObsOut$LFSbias <- rlnorm(nsim, mconv(1, Obs@LFScv), sdconv(1, Obs@LFScv))  # sample of length at full selection bias
  ObsOut$Aerr <- runif(nsim, Obs@Btcv[1], Obs@Btcv[2])
  ObsOut$Abias <- exp(runif(nsim, log(Obs@Btbias[1]), log(Obs@Btbias[2])))  #rlnorm(nsim,mconv(1,Obs@Btbiascv),sdconv(1,Obs@Btbiascv))    # sample of current abundance bias
  ObsOut$Kbias <- rlnorm(nsim, mconv(1, Obs@Kcv), sdconv(1, Obs@Kcv))  # sample of von B. K parameter bias
  ObsOut$t0bias <- rlnorm(nsim, mconv(1, Obs@t0cv), sdconv(1, Obs@t0cv))  # sample of von B. t0 parameter bias
  ObsOut$Linfbias <- rlnorm(nsim, mconv(1, Obs@Linfcv), sdconv(1, Obs@Linfcv))  # sample of von B. maximum length bias
  ObsOut$Irefbias <- rlnorm(nsim, mconv(1, Obs@Irefcv), sdconv(1, Obs@Irefcv))  # sample of bias in reference (target) abundance index
  ObsOut$Crefbias <- rlnorm(nsim, mconv(1, Obs@Crefcv), sdconv(1, Obs@Crefcv))  # sample of bias in reference (target) catch index
  ObsOut$Brefbias <- rlnorm(nsim, mconv(1, Obs@Brefcv), sdconv(1, Obs@Brefcv))  # sample of bias in reference (target) biomass index
  ObsOut$Recsd <- runif(nsim, Obs@Reccv[1], Obs@Reccv[2])  # Recruitment deviation  
  
  ObsOut
}