#' Run a Management Strategy Evaluation
#' 
#' A function that runs a Management Strategy Evaluation (closed-loop
#' simulation) for a specified operating model
#' 
#' 
#' @param OM An operating model object (class OM)
#' @param MPs A vector of methods (character string) of class Output or
#' Input.
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
#' @param CheckMPs Logical to indicate if Can function should be used to check
#' if MPs can be run.
#' @param Hist Should model stop after historical simulations? Returns a list 
#' containing all historical data
#' @param ntrials Maximum of times depletion and recruitment deviations are 
#' resampled to optimize for depletion. After this the model stops if more than 
#' percent of simulations are not close to the required depletion
#' @param fracD maximum allowed proportion of simulations where depletion is not 
#' close to sampled depletion from OM before model stops with error
#' @return An object of class MSE
#' @author T. Carruthers and A. Hordyk
#' @export runMSEdev
runMSEdev <- function(OM = "1", MPs = c("AvC","DCAC","DD","FMSYref","curE"), interval=4,
  pstar = 0.5, maxF = 0.8, timelimit = 1, reps = 1, CheckMPs = FALSE,
  Hist=FALSE, ntrials=50, fracD=0.05) {

  if("seed"%in%slotNames(OM))set.seed(OM@seed)
  tiny <- 1e-15  # define tiny variable
  message("Loading operating model")
  flush.console()
  if (class(OM) != "OM") stop("You must specify an operating model")
  if(length(OM@cpars)>0){
    ncparsim<-cparscheck(OM@cpars)   # check each list object has the same length and if not stop and error report
    cpars<-OM@cpars
  }
  nsim<-OM@nsim
  proyears<-OM@proyears
  nyears <- OM@nyears  # number of  historical years
  maxage <- OM@maxage  # maximum age (no plus group)
  
  calcMax <- -log(0.01)/(min(OM@M))        # Age at which 1% of cohort survives
  maxage <- round(max(maxage, calcMax),0)  # If maximum age is lower, increase it to calcMax
    
  dep <- runif(nsim, OM@D[1], OM@D[2])      # sample from the range of user-specified depletion (Bcurrent/B0)
  Esd <- runif(nsim, OM@Esd[1], OM@Esd[2])  # interannual variability in fishing effort (log normal sd)

  EffLower <- OM@EffLower
  EffUpper <- OM@EffUpper 
  EffYears <- OM@EffYears
  Deriv <- getEffhist(Esd, nyears, EffYears = OM@EffYears, EffLower = OM@EffLower, EffUpper = OM@EffUpper)  # Historical fishing effort
  Find <- Deriv[[1]]  # Calculate fishing effort rate
  dFfinal <- Deriv[[2]]  # Final gradient in fishing effort yr-1 
  # dep[order(dFfinal)]<-dep[order(dep,decreasing=T)] # robustifies
  # matplot(t(Find), type='l') # plot(dep, dFfinal)

  # Sample operating model parameters
  # ===========================================================
  procsd <- runif(nsim, OM@Perr[1], OM@Perr[2])  # Process error standard deviation
  AC <- runif(nsim, OM@AC[1], OM@AC[2])  # auto correlation parameter for recruitment deviations recdev(t)<-AC*recdev(t-1)+(1-AC)*recdev_proposed(t)
  
  M <- runif(nsim, OM@M[1], OM@M[2])  # natural mortality rate \t
  Msd <- runif(nsim, OM@Msd[1], OM@Msd[2])  # sample inter annual variability in M from specified range
  Mgrad <- runif(nsim, OM@Mgrad[1], OM@Mgrad[2])  # sample gradient in M (M y-1)
  hs <- runif(nsim, OM@h[1], OM@h[2])  # sample of recruitment compensation (steepness - fraction of unfished recruitment at 20% of unfished biomass)
  Linf <- runif(nsim, OM@Linf[1], OM@Linf[2])  # sample of asymptotic length
  Linfsd <- runif(nsim, OM@Linfsd[1], OM@Linfsd[2])  # sample of interannual variability in Linf
  Linfgrad <- runif(nsim, OM@Linfgrad[1], OM@Linfgrad[2])  # sample of gradient in Linf (Linf y-1)
  recgrad <- runif(nsim, OM@recgrad[1], OM@recgrad[2])  # gradient in recent recruitment
  K <- runif(nsim, OM@K[1], OM@K[2])  # now predicted by a log-linear model
  Ksd <- runif(nsim, OM@Ksd[1], OM@Ksd[2])  #runif(nsim,OM@Ksd[1],OM@Ksd[2])# sd is already added in the linear model prediction
  Kgrad <- runif(nsim, OM@Kgrad[1], OM@Kgrad[2])  # gradient in Von-B K parameter (K y-1)
  t0 <- runif(nsim, OM@t0[1], OM@t0[2])  # a sample of theoretical age at length zero
  L50 <- array(runif(nsim * 50, OM@L50[1], OM@L50[2]), c(nsim, 50))  # length at 50% maturity
  L50_95 <- array(runif(nsim * 50, OM@L50_95[1], OM@L50_95[2]), c(nsim, 50))  # length at 95% maturity
  
  # checks for unrealistically high length at maturity 
  L50[L50/Linf > 0.95] <- NA
  L50 <- apply(L50, 1, function(x) x[!is.na(x)][1])
  L50_95[(L50+L50_95)/Linf > 0.99] <- NA
  L50_95 <- apply(L50_95, 1, function(x) x[!is.na(x)][1]) 
  L95 <- array(L50 + L50_95)

  Spat_targ <- runif(nsim, OM@Spat_targ[1], OM@Spat_targ[2])  # spatial targetting Ba^targetting param 
  Frac_area_1 <- runif(nsim, OM@Frac_area_1[1], OM@Frac_area_1[2])  # sampled fraction of unfished biomass in area 1 (its a two area model by default)
  Prob_staying <- runif(nsim, OM@Prob_staying[1], OM@Prob_staying[2])  # sampled probability of individuals staying in area 1 among years
  Size_area_1 <- runif(nsim, OM@Size_area_1[1], OM@Size_area_1[2])  # currently redundant parameter for the habitat area size of area 1
  
  # Sample observation error model parameters
  # ===============================================================
  Csd <- runif(nsim, OM@Cobs[1], OM@Cobs[2])  # Sampled catch observation error (lognormal sd)
  Cbias <- rlnorm(nsim, mconv(1, OM@Cbiascv), sdconv(1, OM@Cbiascv))  # Sampled catch bias (log normal sd)
  CAA_nsamp <- ceiling(runif(nsim, OM@CAA_nsamp[1], OM@CAA_nsamp[2]))  # Number of catch-at-age observations
  CAA_ESS <- ceiling(runif(nsim, OM@CAA_ESS[1], OM@CAA_ESS[2]))  # Effective sample size
  CAL_nsamp <- runif(nsim, OM@CAL_nsamp[1], OM@CAL_nsamp[2])  # Observation error standard deviation for single catch at age by area
  CAL_ESS <- ceiling(runif(nsim, OM@CAL_ESS[1], OM@CAL_ESS[2]))  # Effective sample size
  CALcv <- runif(nsim, OM@CALcv[1], OM@CALcv[2])  # Observation error standard deviation for single catch at age by area
  betas <- exp(runif(nsim, log(OM@beta[1]), log(OM@beta[2])))  # the sampled hyperstability / hyperdepletion parameter beta>1 (hyperdepletion) beta<1 (hyperstability)
  Isd <- runif(nsim, OM@Iobs[1], OM@Iobs[2])  # Abundance index observation error (log normal sd)
  Derr <- runif(nsim, OM@Dcv[1], OM@Dcv[2])
  Dbias <- rlnorm(nsim, mconv(1, OM@Dbiascv), sdconv(1, OM@Dbiascv))  # sample of depletion bias
  Mbias <- rlnorm(nsim, mconv(1, OM@Mcv), sdconv(1, OM@Mcv))  # sample of M bias
  FMSY_Mbias <- rlnorm(nsim, mconv(1, OM@FMSY_Mcv), sdconv(1, OM@FMSY_Mcv))  # sample of FMSY/M bias

  lenMbias <- rlnorm(nsim, mconv(1, OM@LenMcv), sdconv(1, OM@LenMcv))  # sample of length at maturity bias - assume same error as age based maturity
  LFCbias <- rlnorm(nsim, mconv(1, OM@LFCcv), sdconv(1, OM@LFCcv))  # sample of length at first capture bias
  LFSbias <- rlnorm(nsim, mconv(1, OM@LFScv), sdconv(1, OM@LFScv))  # sample of length at full selection bias
  Aerr <- runif(nsim, OM@Btcv[1], OM@Btcv[2])
  Abias <- exp(runif(nsim, log(OM@Btbias[1]), log(OM@Btbias[2])))  #rlnorm(nsim,mconv(1,OM@Btbiascv),sdconv(1,OM@Btbiascv))    # smaple of current abundance bias
  Kbias <- rlnorm(nsim, mconv(1, OM@Kcv), sdconv(1, OM@Kcv))  # sample of von B. K parameter bias
  t0bias <- rlnorm(nsim, mconv(1, OM@t0cv), sdconv(1, OM@t0cv))  # sample of von B. t0 parameter bias
  Linfbias <- rlnorm(nsim, mconv(1, OM@Linfcv), sdconv(1, OM@Linfcv))  # sample of von B. maximum length bias
  Irefbias <- rlnorm(nsim, mconv(1, OM@Irefcv), sdconv(1, OM@Irefcv))  # sample of bias in reference (target) abundance index
  Crefbias <- rlnorm(nsim, mconv(1, OM@Crefcv), sdconv(1, OM@Crefcv))  # sample of bias in reference (target) catch index
  Brefbias <- rlnorm(nsim, mconv(1, OM@Brefcv), sdconv(1, OM@Brefcv))  # sample of bias in reference (target) biomass index
  Recsd <- runif(nsim, OM@Reccv[1], OM@Reccv[2])  # Recruitment deviation 
  
  # Sample fishing efficiency parameters
  # =======================================================
  qinc <- runif(nsim, OM@qinc[1], OM@qinc[2])
  qcv <- runif(nsim, OM@qcv[1], OM@qcv[2])  # interannual variability in catchability
  
  # Sample selectivity parameters 
  # =======================================================  
  
  Selnyears <- length(OM@SelYears)
  # are selectivity parameters relative to size at maturity?
  chk <- class(OM@isRel)
  if (length(OM@isRel) < 1) 
    OM@isRel <- "true"
  if (chk == "character") {
    chkRel <- tolower(OM@isRel)
    if (chkRel == "true" | OM@isRel == "1") 
      multi <- L50
    if (chkRel == "false" | OM@isRel == "0") 
      multi <- 1
  }
  if (chk == "numeric") {
    if (OM@isRel == 1) 
      multi <- L50
    if (OM@isRel == 0) 
      multi <- 1
  }
  L5 <- runif(nsim, OM@L5[1], OM@L5[2]) * multi  # length at 0.05% selectivity ascending
  LFS <- runif(nsim, OM@LFS[1], OM@LFS[2]) * multi  # first length at 100% selection
  Vmaxlen <- runif(nsim, OM@Vmaxlen[1], OM@Vmaxlen[2])  # selectivity at maximum length
  L5s <- LFSs <- Vmaxlens <- NULL  # initialize 
  
  if (Selnyears > 1) {   # change of selectivity in historical years 
    # length at 0.05% selectivity ascending
	L5s <- mapply(runif, n = nsim, min = OM@L5Lower, max = OM@L5Upper) * multi
    # first length at 100% selection
    LFSs <- mapply(runif, n = nsim, min = OM@LFSLower, max = OM@LFSUpper) *  multi
	# selectivity at maximum length
	Vmaxlens <- mapply(runif, n = nsim, min = OM@VmaxLower, max = OM@VmaxUpper)
  }
  
  ### End of sampling OM parameters ###
   
  # dat<-as.data.frame(cbind(procsd,AC,M,Msd,Mgrad,hs,Linf,Linfsd,Linfgrad,recgrad,K,Ksd,Kgrad,t0,L50,L95,L5,LFS,
  # Vmaxlen,Spat_targ,Frac_area_1,Prob_staying,Size_area_1,Csd,Cbias,CAA_nsamp,CAA_ESS,CALcv,betas,
  # Isd,Derr,Dbias,Mbias,FMSY_Mbias,lenMbias,LFCbias,LFSbias,Aerr,Abias,Kbias,
  # t0bias,Linfbias,Irefbias,Crefbias,Brefbias,Recsd,qinc,qcv))
  
  # save(dat,file='F:/DLM/Operating models/Other/custompars')
   
  # Recruitment Deviations
  procmu <- -0.5 * (procsd)^2  # adjusted log normal mean
  Perr <- array(rnorm((nyears + proyears) * nsim, rep(procmu, nyears + 
    proyears), rep(procsd, nyears + proyears)), c(nsim, nyears + proyears))
  for (y in 2:(nyears + proyears)) Perr[, y] <- AC * Perr[, y - 1] + 
    Perr[, y] * (1 - AC * AC)^0.5  #2#AC*Perr[,y-1]+(1-AC)*Perr[,y] # apply a pseudo AR1 autocorrelation to rec devs (log space)
  Perr <- exp(Perr)  # normal space (mean 1 on average)
  
  # Add cycle (phase shift) to recruitment deviations - if specified
  if (is.finite(OM@Period[1]) & is.finite(OM@Amplitude[1])) {
    Shape <- "sin"  # default sine wave - alternative - 'shift' for step changes
    recMulti <- t(sapply(1:nsim, SetRecruitCycle, Period = OM@Period, 
      Amplitude = OM@Amplitude, TotYears = nyears + proyears, Shape = Shape))
    Perr <- Perr * recMulti  # Add cyclic pattern to recruitment
    message("Adding cyclic recruitment pattern")
    flush.console()
  }
  
  if (nsim > 1) {
    cumlRecDev <- apply(Perr[, 1:nyears], 1, prod)
    dep[order(cumlRecDev)] <- dep[order(dep, decreasing = F)]  # robustifies 
  }
  
  R0 <- OM@R0  # Initial recruitment
  if (length(R0) != nsim) R0 <- rep(R0, nsim*50)[1:nsim] # modified to allow for different R0 per sim 
   
  # Generate random numbers for random walk - done here so that they can be written
  # out to SampPars 
  Mrand <- matrix(exp(rnorm(nsim*(proyears+nyears), -0.5 * Msd^2, Msd)), nrow=nsim, ncol=proyears+nyears)
  Linfrand <- matrix(exp(rnorm(nsim*(proyears+nyears), -0.5 * Linfsd^2, Linfsd)), nrow=nsim, ncol=proyears+nyears)
  Krand <- matrix(exp(rnorm(nsim*(proyears+nyears), -0.5 * Ksd^2, Ksd)), nrow=nsim, ncol=proyears+nyears)
  
  # Vector of valid names for custompars list or data.frame. Names not in this list will be printed out in warning and ignored #	
  ParsNames <- c("dep","Esd","Find","procsd","AC","M","Msd", 
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
  
  
  # Sample custom parameters
  # ===================================================================
  if (length(OM@cpars) > 0) { # custom parameters exist
     
	  Names <- names(cpars)
	  # report not valid names 
	  invalid <- which(!Names %in% ParsNames)
	  if (length(invalid) > 0) {
	    outNames <- paste(Names[invalid], "")
	    for (i in seq(5, by=5, length.out=floor(length(outNames)/5))) outNames <- gsub(outNames[i], paste0(outNames[i], "\n"), outNames)
	    warning("ignoring invalid names found in custom parameters (OM@cpars) \n", outNames)	
	  }
	  # report found names
	  valid <- which(Names %in% ParsNames)
	  cpars <- cpars[valid]
	  if (length(OM@cpars) == 0) stop("No valid names found in custompars")
	  Names <- names(cpars)
	  outNames <- paste(Names, "")
	  for (i in seq(5, by=5, length.out=floor(length(outNames)/5)))
  	  outNames <- gsub(outNames[i], paste0(outNames[i], "\n"), outNames)
	    message("valid custom parameters (OM@cpars) found: \n", outNames)
      flush.console()
	  if (ncparsim < nsim) ind <- sample(1:ncparsim, nsim, replace=TRUE)
	  if (!ncparsim < nsim) ind <- sample(1:ncparsim, nsim, replace=FALSE)
	
	  usedName <- 0 	
    for (i in 1:length(cpars)) {
	    
      samps <- cpars[[i]]
	    name <- names(cpars)[i]
	    if (any(c("EffUpper", "EffLower", "EffYears", "maxage") %in% name)) {
	      assign(name, samps)
		    usedName <- usedName + 1
	    } else {
	      if (class(samps) == "numeric" | class(samps) == "integer") {
 		      assign(name, samps[ind])
		      usedName <- usedName + 1
		    }
	      if (class(samps) == "matrix") {
		      assign(name, samps[ind,, drop=FALSE])
		      usedName <- usedName + 1
		    }
		    if (class(samps) == "array") {
		      if (length(dim(samps)) == 3) {
		        assign(name, samps[ind, , ,drop=FALSE])
			      usedName <- usedName + 1 
          }
		    }
	    }	
    }
	  
	  if ("EffUpper" %in% Names & !"Find" %in% Names) {
      Deriv <- getEffhist(Esd, nyears, EffYears = EffYears, EffUpper = EffUpper, EffLower = EffLower)  # Historical fishing effort
      Find <- Deriv[[1]]  # Calculate fishing effort rate
      dFfinal <- Deriv[[2]]  # Final gradient in fishing effort yr-1 
    }
   
  }
  
  # Checks that parameters are correct dimensions - could be messed up with custompars   
  if (length(maxage) > 1) maxage <- maxage[1] # check if maxage has been passed in custompars
  OM@maxage <- maxage # update OM object with maxage that is used 
  if (dim(Perr)[2] != proyears + nyears){
     
      procmu <- -0.5 * (procsd)^2  # adjusted log normal mean
      Perr_nu<-cbind(log(Perr),array(rnorm(proyears * nsim, rep(procmu, proyears), rep(procsd, proyears)), c(nsim, proyears)))
      
      for (y in (nyears+1):(nyears+proyears)) Perr_nu[, y] <- AC * Perr_nu[, y - 1]  +  Perr_nu[, y] * (1 - AC * AC)^0.5   # apply a pseudo AR1 autocorrelation to rec devs (log space)
      Perr <- exp(Perr_nu)  # normal space (mean 1 on average)
     
  }
  
  if(exists("V",inherits=FALSE))if(dim(V)[3] != proyears + nyears)   V<-abind(V,array(V[,,nyears],c(nsim,maxage,proyears)),along=3) # extend future Vulnerabiliy according to final historical vulnerability
   
  if (any(dim(Find) != c(nsim, nyears))) stop("Find must be matrix with dimensions: nsim, nyears")
      
  SRrel <- rep(OM@SRrel, nsim)  # type of Stock-recruit relationship. 1=Beverton Holt, 2=Ricker
  
  Marray <- gettempvar(M, Msd, Mgrad, nyears + proyears, nsim, Mrand)  # M by sim and year according to gradient and inter annual variability
  Linfarray <- gettempvar(Linf, Linfsd, Linfgrad, nyears + proyears, nsim, Linfrand)  # Linf array
  Karray <- gettempvar(K, Ksd, Kgrad, nyears + proyears, nsim, Krand)  # the K array
  
  Agearray <- array(rep(1:maxage, each = nsim), dim = c(nsim, maxage))  # Age array
  Len_age <- array(NA, dim = c(nsim, maxage, nyears + proyears))  # Length at age array
  ind <- as.matrix(expand.grid(1:nsim, 1:maxage, 1:(nyears + proyears)))  # an index for calculating Length at age
  Len_age[ind] <- Linfarray[ind[, c(1, 3)]] * (1 - exp(-Karray[ind[, c(1, 3)]] * 
    (Agearray[ind[, 1:2]] - t0[ind[, 1]])))
  if (!exists("Wt_age", inherits=FALSE)|length(OM@cpars)==0) {
    Wt_age <- array(NA, dim = c(nsim, maxage, nyears + proyears))  # Weight at age array
    Wt_age[ind] <- OM@a * Len_age[ind]^OM@b  # Calculation of weight array
  }	
  
  # Calcaluate age at maturity 
  if (!exists("ageM", inherits=FALSE) |length(OM@cpars)==0) ageM <- -((log(1 - L50/Linf))/K) + t0  # calculate ageM from L50 and growth parameters (non-time-varying)
  ageM[ageM < 1] <- 1  # age at maturity must be at least 1
  L95 <- L50 + L50_95 # reassign if updated in custompars 
  if (!exists("age95", inherits=FALSE) |length(OM@cpars)==0) age95 <- -((log(1 - L95/Linf))/K) + t0
  age95[age95 < 1] <- 1.5  # must be greater than 0 and ageM
  
  ageMsd <- sapply(1:nsim, getroot, ageM, age95)
  ageMarray <- array(ageM, dim = c(nsim, maxage))  # Age at maturity array
  if (!exists("Mat_age", inherits=FALSE)|length(OM@cpars)==0) Mat_age <- 1/(1 + exp((ageMarray - (Agearray))/(ageMarray * ageMsd)))  # Maturity at age array

  
  # Catch at Length Classes
  LatASD <- Len_age * 0.1  # SD of length-at-age - this is currently fixed to cv of 10%
 
  MaxBin <- ceiling(max(Linfarray) + 3 * max(LatASD))
  binWidth <- ceiling(0.03 * MaxBin)
  CAL_bins <- seq(from = 0, to = MaxBin + binWidth, by = binWidth)
  CAL_binsmid <- seq(from = 0.5 * binWidth, by = binWidth, length = length(CAL_bins) - 1)
  nCALbins <- length(CAL_binsmid)
  
  
  # Selectivity at Length
  # ------------------------------------------------------ if (max(OM@L5)
  # > 1.5) { message('L5 set too high (maximum value of 1.5).
  # \nDefaulting to L5 = 1.5') OM@L5[OM@L5 > 1.5] <- 1.5 }
  SLarray <- array(NA, dim=c(nsim, nCALbins, nyears+proyears)) # Selectivity-at-length 
  
  if (exists("V", inherits=FALSE) & length(OM@cpars)>0) { # V has been passed in with custompars 
    # assign L5, LFS and Vmaxlen - dodgy loop 
	# should calculate length at 5% selectivity from vB 
	L5 <- matrix(NA, nrow = nyears + proyears, ncol = nsim)
    LFS <- matrix(NA, nrow = nyears + proyears, ncol = nsim)
    Vmaxlen <- matrix(NA, nrow = nyears + proyears, ncol = nsim)

	for (yr in 1:(nyears+proyears)) {
	  for (s in 1:nsim) {
	    ind <- min(which(V[s,,yr] >=0.05))
	    L5[yr, s] <- Len_age[s, ind, yr]
	    ind2 <- min(which(V[s,,yr] >=0.50))
		if (ind2 == ind) ind2 <- ind + 1
	    LFS[yr, s] <- Len_age[s, ind2, yr]
		Vmaxlen[yr, s] <- V[s, maxage, yr]
		SLarray[s,, yr] <- SelectFun(s, SL0.05=L5[yr, ], SL1=LFS[yr, ], MaxSel=Vmaxlen[yr, ], 
 	                            maxlens=Len_age[, maxage, nyears], Lens=CAL_binsmid)
	  }
	}
  }
  
  maxlen <- Len_age[, maxage, nyears] # reference length for Vmaxlen 
  # assume it is expected length at maximum age for current (nyears) year 
  
  if (!exists("V", inherits=FALSE) |length(OM@cpars)==0) { # don't run if V has been passed in with custompars 
    if (Selnyears <= 1) {    
      L5 <- matrix(L5, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
      LFS <- matrix(LFS, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
      Vmaxlen <- matrix(Vmaxlen, nrow = nyears + proyears, ncol = nsim, byrow = TRUE) 
    
      ind <- which(LFS/matrix(Linf, nrow = proyears + nyears, ncol = nsim, byrow = TRUE) > 1, arr.ind = T)
      if (length(ind) > 0) {
        message("LFS too high (LFS > Linf) in some cases. \nDefaulting to LFS = 0.9 Linf for the affected simulations")
        LFS[ind] <- Linf[ind[, 2]] * 0.9
      } 
	  
	  # Calculate selectivity-at-age  curve 
	  V <- array(NA, dim = c(nsim, maxage, nyears + proyears)) 
      s1 <- sapply(1:nsim, function(i) optimize(getSlope1, interval = c(0, 1e+05), 
        LFS = LFS[1, i], L0.05 = L5[1,i])$minimum)	
	  if (all(Vmaxlen >= 0.99)) s2 <- rep(1E5, nsim)
      if (!all(Vmaxlen >= 0.99)) 
	    s2 <- sapply(1:nsim, function(i) optimize(getSlope2, interval = c(0, 1e+05), 
	  	             LFS = LFS[1,i], s1=s1[i], maxlen=maxlen[i], 
	  				 MaxSel=Vmaxlen[1, i])$minimum)
      for (yr in 1:(nyears+proyears)) {
  	   # Calculate selectivity at age class 
	   V[ , , yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[1,i], s1[i], s2[i], lens=Len_age[i,,yr])))
	   # Calculate selectivity at length class 
	   SLarray[,, yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[1,i], s1[i], s2[i], lens=CAL_binsmid)))   
	  }	 
    }
	     
    if (Selnyears > 1) {
      # More than one break point in historical selection pattern
      L5 <- matrix(0, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
      LFS <- matrix(0, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
      Vmaxlen <- matrix(0, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
      SelYears <- OM@SelYears
	  
	  ind <- which(LFSs/ matrix(Linf, nrow=nsim, ncol=Selnyears) > 1, arr.ind = T)
      if (length(ind) > 0) {
        message("LFS too high (LFS > Linf) in some cases. \nDefaulting to LFS = 0.9 Linf for the affected simulations")
        LFSs[ind] <- Linf[ind[, 1]] * 0.9
      }     
	  
	  # Calculate selectivity-at-age  curve 
	  V <- array(NA, dim = c(nsim, maxage, nyears + proyears))     
      
	  for (X in 1:(Selnyears - 1)) {	
        bkyears <- SelYears[X]:SelYears[X + 1]
        L5[bkyears, ] <- matrix(rep((L5s[, X]), length(bkyears)), ncol = nsim, byrow = TRUE)
        LFS[bkyears, ] <- matrix(rep((LFSs[, X]), length(bkyears)), ncol = nsim, byrow = TRUE)
        Vmaxlen[bkyears, ] <- matrix(rep((Vmaxlens[, X]), length(bkyears)), ncol = nsim, byrow = TRUE)
        
	      s1 <- sapply(1:nsim, function(i) optimize(getSlope1, interval = c(0, 1e+05), 
          LFS = LFSs[i, X], L0.05 = L5s[i, X])$minimum)
	      s2 <- sapply(1:nsim, function(i) optimize(getSlope2, interval = c(0, 1e+05), 
	  	             LFS = LFSs[i, X], s1=s1[i], maxlen=maxlen[i], 
	  				 MaxSel=Vmaxlens[i, X])$minimum)	
	      for (yr in bkyears) {
  	      V[ , , yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[yr, i], s1[i], s2[i], lens=Len_age[i,,yr])))
          SLarray[,, yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[1,i], s1[i], s2[i], lens=CAL_binsmid)))   		 
		    }
      }
	  
      restYears <- max(SelYears):(nyears + proyears)
      L5[restYears, ] <- matrix(rep((L5s[, Selnyears]), length(restYears)), ncol = nsim, byrow = TRUE)
      LFS[restYears, ] <- matrix(rep((LFSs[, Selnyears]), length(restYears)), ncol = nsim, byrow = TRUE)
      Vmaxlen[restYears, ] <- matrix(rep((Vmaxlens[, Selnyears]), length(restYears)), ncol = nsim, byrow = TRUE)
    
      s1 <- sapply(1:nsim, function(i) optimize(getSlope1, interval = c(0, 1e+05), 
          LFS = LFSs[i, Selnyears], L0.05 = L5s[i, Selnyears])$minimum)
	    s2 <- sapply(1:nsim, function(i) optimize(getSlope2, interval = c(0, 1e+05), 
	  	             LFS = LFSs[i, Selnyears], s1=s1[i], maxlen=maxlen[i], 
	  				 MaxSel=Vmaxlens[i, Selnyears])$minimum)	
	    for (yr in restYears) { 
  	     V[ , , restYears] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[yr, i], s1[i], s2[i], lens=Len_age[i,,yr])))		
		    SLarray[,, yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[1,i], s1[i], s2[i], lens=CAL_binsmid))) 
	    }	 
    }
  } # end of 'if V exists'
   
  if (any((dim(V) != c(nsim, maxage, proyears+nyears)))) 
    stop("V must have dimensions: nsim (", nsim,") maxage (", maxage, 
	      ") proyears+nyears (", proyears+nyears, ") \nbut has ", 
	      dim(V)[1], " ", dim(V)[2], " ", dim(V)[3], call.=FALSE)
  
 
  Asize <- cbind(Size_area_1, 1 - Size_area_1)
  
  message("Optimizing for user-specified movement")  # Print a progress update
  flush.console()  # refresh the console
  
  if (snowfall::sfIsRunning()) {
    # if the cluster is initiated
    snowfall::sfExport(list = c("Frac_area_1", "Prob_staying"))  # export some of the new arrays and ...
    mov <- array(t(snowfall::sfSapply(1:nsim, getmov2, Frac_area_1 = Frac_area_1, 
      Prob_staying = Prob_staying)), dim = c(nsim, 2, 2))  # numerically determine movement probability parameters to match Prob_staying and Frac_area_1
    # mov <- array(t(snowfall::sfSapply(1:nsim, getmov, Frac_area_1 = Frac_area_1, 
      # Prob_staying = Prob_staying)), dim = c(nsim, 2, 2))  # numerically determine movement probability parameters to match Prob_staying and Frac_area_1
  } else {
    # no cluster initiated
	mov <- array(t(sapply(1:nsim, getmov2, Frac_area_1 = Frac_area_1, 
      Prob_staying = Prob_staying)), dim = c(nsim, 2, 2))  # numerically determine movement probability parameters to match Prob_staying and Frac_area_1
    # mov <- array(t(sapply(1:nsim, getmov, Frac_area_1 = Frac_area_1, 
      # Prob_staying = Prob_staying)), dim = c(nsim, 2, 2))  # numerically determine movement probability parameters to match Prob_staying and Frac_area_1	  
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
  surv <- exp(-Marray[, 1])^(Agearray - 1)  # Survival array
  Nfrac <- surv * Mat_age  # predicted Numbers of mature ages
  initdist <- as.matrix(cbind(Frac_area_1, 1 - Frac_area_1))  # Get the initial spatial distribution of each simulated population
  
  R0a <- matrix(R0, nrow=nsim, ncol=nareas, byrow=FALSE) * initdist  # Unfished recruitment by area
  
  SAYR <- as.matrix(expand.grid(1:nareas, 1, 1:maxage, 1:nsim)[4:1])  # Set up some array indexes sim (S) age (A) year (Y) region/area (R)
  SAY <- SAYR[, 1:3]
  SA <- SAYR[, 1:2]
  SR <- SAYR[, c(1, 4)]
  S <- SAYR[, 1]
  SY <- SAYR[, c(1, 3)]
  
  SSN[SAYR] <- Nfrac[SA] * R0[S] * initdist[SR]  # Calculate initial spawning stock numbers
  N[SAYR] <- R0[S] * surv[SA] * initdist[SR]  # Calculate initial stock numbers
  
  Biomass[SAYR] <- N[SAYR] * Wt_age[SAY]  # Calculate initial stock biomass
  SSB[SAYR] <- SSN[SAYR] * Wt_age[SAY]    # Calculate spawning stock biomass
  VBiomass[SAYR] <- Biomass[SAYR] * V[SAY]  # Calculate vunerable biomass
  
  if (nsim > 1) {
    SSN0 <- apply(SSN[, , 1, ], c(1, 3), sum)  # Calculate unfished spawning stock numbers  
    SSB0 <- apply(SSB[, , 1, ], 1, sum)  # Calculate unfished spawning stock biomass
    SSBpR <- SSB0/R0  # Spawning stock biomass per recruit
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
  flush.console()  # update console
  
  if (snowfall::sfIsRunning()) {
    snowfall::sfExport(list = c("dep", "Find", "Perr", "Marray", "hs", "Mat_age", 
      "Wt_age", "R0", "V", "nyears", "maxage", "SRrel", "aR", "bR"))
    qs <- snowfall::sfSapply(1:nsim, getq2, dep, Find, Perr, Marray, hs, Mat_age, 
      Wt_age, R0, V, nyears, maxage, mov, Spat_targ, SRrel, aR, bR)  # find the q that gives current stock depletion
    # qs <- snowfall::sfSapply(1:nsim, getq, dep, Find, Perr, Marray, hs, Mat_age, 
      # Wt_age, R0, V, nyears, maxage, mov, Spat_targ, SRrel, aR, bR)  # find the q that gives current stock depletion	  
  } else {
    qs <- sapply(1:nsim, getq2, dep, Find, Perr, Marray, hs, Mat_age, 
      Wt_age, R0, V, nyears, maxage, mov, Spat_targ, SRrel, aR, bR)  # find the q that gives current stock depletion
    # qs <- sapply(1:nsim, getq, dep, Find, Perr, Marray, hs, Mat_age, 
      # Wt_age, R0, V, nyears, maxage, mov, Spat_targ, SRrel, aR, bR)  # find the q that gives current stock depletion	  
  }

  # Check that depletion target is reached
  UpperBound <- 13  # bounds for q (catchability). Flag if bounded optimizer hits the bounds 
  LowerBound <- 0.008
  HighQ <- which(qs > UpperBound | qs < LowerBound)
  if (length(HighQ) > 0) {
    # If q has hit bound, re-sample depletion and try again. Tries 30 times
    # and then alerts user
    Err <- TRUE
    # Nsec <- 10
    Nprob <- length(HighQ)
	
    message(Nprob,' simulations have final biomass that is not close to sampled depletion') 
	message('Re-sampling depletion, recruitment error, and fishing effort')
    flush.console()
    count <- 0
    while (Err & count < ntrials) {
      count <- count + 1
      Nprob <- length(HighQ)
	  
	  # Re-sample depletion 
      dep[HighQ] <- runif(Nprob, OM@D[1], OM@D[2])
      
	  # Re-sample recruitment deviations
	  procsd[HighQ] <- runif(Nprob, OM@Perr[1], OM@Perr[2])  # Re-sample process error standard deviation 
      AC[HighQ] <- runif(Nprob, OM@AC[1], OM@AC[2])  # Re-sample auto correlation parameter for recruitment deviations recdev(t)<-AC*recdev(t-1)+(1-AC)*recdev_proposed(t)  
      procmu2 <- -0.5 * (procsd[HighQ])^2  # adjusted log normal mean
      Perr2 <- array(rnorm((nyears + proyears) * length(HighQ), rep(procmu[HighQ], 
        nyears + proyears), rep(procsd[HighQ], nyears + proyears)), 
        c(length(HighQ), nyears + proyears))	
      for (y in 2:(nyears + proyears)) Perr2[, y] <- AC[HighQ] * 
        Perr2[, y - 1] + Perr2[, y] * (1 - AC[HighQ] * AC[HighQ])^0.5
      Perr[HighQ, ] <- exp(Perr2)  # normal space (mean 1 on average)
      
	  if (exists("recMulti", inherits=FALSE))  Perr[HighQ,] <- Perr[HighQ,] * recMulti[HighQ,]
	  
	  # Re-sample historical fishing effort 
	  Esd2 <- runif(Nprob, OM@Esd[1], OM@Esd[2])
	  Esd[HighQ] <- Esd2
	  Deriv2 <- getEffhist(Esd2, nyears, EffYears = EffYears, EffLower = EffLower, EffUpper = EffUpper)  # Historical fishing effort
      Find[HighQ, ] <- Deriv2[[1]]  # Calculate fishing effort rate
      dFfinal[HighQ] <- Deriv2[[2]]  # Final gradient in fishing effort yr-1 
	  
      if (snowfall::sfIsRunning()) {
        snowfall::sfExport(list = c("dep", "Find", "Perr", "Marray", "hs", 
          "Mat_age", "Wt_age", "R0", "V", "nyears", "maxage", "SRrel", 
          "aR", "bR"))
        qs[HighQ] <- snowfall::sfSapply(HighQ, getq2, dep, Find, Perr, Marray, 
          hs, Mat_age, Wt_age, R0, V, nyears, maxage, mov, Spat_targ, 
          SRrel, aR, bR)  # find the q that gives current stock depletion
        # qs[HighQ] <- snowfall::sfSapply(HighQ, getq, dep, Find, Perr, Marray, 
          # hs, Mat_age, Wt_age, R0, V, nyears, maxage, mov, Spat_targ, 
          # SRrel, aR, bR)  # find the q that gives current stock depletion		  
      } else {
        qs[HighQ] <- sapply(HighQ, getq2, dep, Find, Perr, Marray, 
          hs, Mat_age, Wt_age, R0, V, nyears, maxage, mov, Spat_targ, 
          SRrel, aR, bR)  # find the q that gives current stock depletion
        #qs[HighQ] <- sapply(HighQ, getq, dep, Find, Perr, Marray, 
          # hs, Mat_age, Wt_age, R0, V, nyears, maxage, mov, Spat_targ, 
          # SRrel, aR, bR)  # find the q that gives current stock depletion		  
      }
	  
      HighQ <- which(qs > UpperBound | qs < LowerBound)
      if (length(HighQ) == 0) Err <- FALSE
    }
    if (Err) {
      # still a problem
      tooLow <- length(which(qs > UpperBound))
      tooHigh <- length(which(qs < LowerBound))
      prErr <- length(HighQ)/nsim
      if (prErr > fracD & length(HighQ) >= 1) {
        if (length(tooLow) > 0) 
          message(tooLow, " sims can't get down to the lower bound on depletion")
        if (length(tooHigh) > 0) 
          message(tooHigh, " sims can't get to the upper bound on depletion")
        message("More than ", fracD*100, "% of simulations can't get to the specified level of depletion with these Operating Model parameters")
        stop("Try again for a complete new sample, modify the input parameters, or increase ")
      } else {
        if (length(tooLow) > 0) 
          message(tooLow, " sims can't get down to the lower bound on depletion")
        if (length(tooHigh) > 0) 
          message(tooHigh, " sims can't get to the upper bound on depletion")
        message("Less than ", fracD*100, "% simulations can't get to the sampled depletion.\nContinuing")
      }
      
    }
  }
  
  message("Calculating historical stock and fishing dynamics")  # Print a progress update
  flush.console()  # update console
  
  if (nsim > 1) fishdist <- (apply(VBiomass[, , 1, ], c(1, 3), sum)^Spat_targ)/
                             apply(apply(VBiomass[, , 1, ], c(1, 3), sum)^Spat_targ, 1, mean)  # spatial preference according to spatial biomass
  if (nsim == 1)  fishdist <- (matrix(apply(VBiomass[,,1,], 2, sum), nrow=nsim)^Spat_targ)/
                               mean((matrix(apply(VBiomass[,,1,], 2, sum), nrow=nsim)^Spat_targ))
  	
  FM[SAYR] <- qs[S] * Find[SY] * V[SAY] * fishdist[SR]  # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
  Z[SAYR] <- FM[SAYR] + Marray[SY]  # Total mortality rate                 
  
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
      N[, 1, y + 1, ] <- Perr[, y] * (0.8 * R0a * hs * 
	               apply(SSB[, , y, ], MAR, sum))/(0.2 * SSBpR * R0a * (1 - hs) + 
                   (hs - 0.2) * apply(SSB[, , y, ], MAR, sum))  # Recruitment assuming regional R0 and stock wide steepness
    } else {
      # most transparent form of the Ricker uses alpha and beta params
      N[, 1, y + 1, ] <- Perr[, y] * aR * apply(SSB[, , y, ], MAR, sum) * 
	  exp(-bR * apply(SSB[, , y, ], MAR, sum))
    }
    
    if (nsim > 1) fishdist <- (apply(VBiomass[, , y, ], c(1, 3), sum)^Spat_targ)/
	                           apply(apply(VBiomass[, , y, ], c(1, 3), sum)^Spat_targ, 1, mean)  # spatial preference according to spatial biomass
    if (nsim == 1)  fishdist <- (matrix(apply(VBiomass[,, y,], 2, sum), nrow=nsim)^Spat_targ)/
                               mean((matrix(apply(VBiomass[,,y,], 2, sum), nrow=nsim)^Spat_targ))							   
    FM[SAY1R] <- qs[S] * Find[SY1] * V[SAY] * fishdist[SR]  # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
    Z[SAY1R] <- FM[SAY1R] + Marray[SY]  # Total mortality rate
    N[, 2:maxage, y + 1, ] <- N[, 1:(maxage - 1), y, ] * exp(-Z[, 1:(maxage - 1), y, ])  # Total mortality
    temp <- array(N[indMov2] * mov[indMov3], dim = c(nareas, nareas, maxage, nsim))  # Move individuals
    N[, , y + 1, ] <- apply(temp, c(4, 3, 1), sum)
    Biomass[SAY1R] <- N[SAY1R] * Wt_age[SAY]  # Calculate biomass
    VBiomass[SAY1R] <- Biomass[SAY1R] * V[SAY]  # Calculate vulnerable biomass
    SSN[SAY1R] <- N[SAY1R] * Mat_age[SA]  # Calculate spawning stock numbers
    SSB[SAY1R] <- SSN[SAY1R] * Wt_age[SAY]  # Calculate spawning stock biomass
	  
  }  # end of year
  
  # Depletion <- apply(Biomass[, , nyears, ], 1, sum)/apply(Biomass[, , 1, ], 1, sum)  #^betas   # apply hyperstability / hyperdepletion
  if (nsim > 1) Depletion <- (apply(SSB[,,nyears,],1,sum)/apply(SSB[,,1,],1,sum))#^betas
  if (nsim == 1) Depletion <- sum(SSB[,,nyears,])/sum(SSB[,,1,])#^betas
  # # apply hyperstability / hyperdepletion
  
  # print(paste("Depletion: ", round(cbind(dep,Depletion),2)))
  
  CN <- apply(N * (1 - exp(-Z)) * (FM/Z), c(1, 3, 2), sum)  # Catch in numbers
  CN[is.na(CN)] <- 0
  CB <- Biomass * (1 - exp(-Z)) * (FM/Z)  # Catch in biomass
  
  Cbiasa <- array(Cbias, c(nsim, nyears + proyears))  # Bias array
  Cerr <- array(rlnorm((nyears + proyears) * nsim, mconv(1, rep(Csd, 
    (nyears + proyears))), sdconv(1, rep(Csd, nyears + proyears))), 
    c(nsim, nyears + proyears))  # composite of bias and observation error
  Cobs <- Cbiasa[, 1:nyears] * Cerr[, 1:nyears] * apply(CB, c(1, 3), sum)  # Simulated observed catch (biomass)
  
  CAA <- array(NA, dim = c(nsim, nyears, maxage))  # Catch  at age array
  cond <- apply(CN, 1:2, sum, na.rm = T) < 1  # this is a fix for low sample sizes. If CN is zero across the board a single fish is caught in age class of model selectivity (dumb I know)
  fixind <- as.matrix(cbind(expand.grid(1:nsim, 1:nyears), rep(floor(maxage/3), nyears)))  # more fix
  CN[fixind[cond, ]] <- 1  # puts a catch in the most vulnerable age class
  # for(i in 1:nsim)for(j in
  # 1:nyears)CAA[i,j,]<-ceiling(-0.5+rmultinom(1,CAA_nsamp[i],CN[i,j,])*CAA_nsamp[i]/CAA_ESS[i])
  # # a multinomial observation model for catch-at-age data
  for (i in 1:nsim) for (j in 1:nyears) CAA[i, j, ] <- ceiling(-0.5 + 
    rmultinom(1, CAA_ESS[i], CN[i, j, ]) * CAA_nsamp[i]/CAA_ESS[i])  # a multinomial observation model for catch-at-age data
  
  # for (i in 1:nsim) {
    # for (j in 1:nyears) {
      # tempCN <- ceiling(-0.5 + rmultinom(1, size = CAL_ESS[i], prob = CN[i, j, ]) * CAL_nsamp[i]/CAL_ESS[i])
      # # ages <- rep(1:maxage,tempCN)+runif(sum(tempCN),-0.5,0.5) # sample
      # # expected age	  
      # lens <- unlist(sapply(1:maxage, function(X) 
	    # rnorm(tempCN[X], Len_age[i, X, j], LatASD[i, X, j])))
      # lens[lens > (max(Linfarray) + 2 * max(LatASD)) | lens > max(CAL_bins)] <- max(Linfarray) + 
        # 2 * max(LatASD)  # truncate at 2 sd 
      # CAL[i, j, ] <- hist(lens, CAL_bins, plot = F)$counts  # assign to bins
      # LFC[i] <- min(c(lens, LFC[i]), na.rm = T)  # get the smallest CAL observation
      # # CAL[i,j,] <- ceiling(rmultinom(1, size=ESS[i],
      # # prob=tempCAL)*CAL_nsamp[i]*CAL_ESS[i]-0.5) # could replace with
      # # Dirichlet distribution
    # }
  # }
  
  # # a multinomial observation model for catch-at-length data
  CAL <- array(NA, dim=c(nsim,  nyears, nCALbins))
  LFC <- rep(NA, nsim)
  for (i in 1:nsim) { # Rcpp code 
    CAL[i, , ] <-  genLenComp(CAL_bins, CAL_binsmid, SLarray[i,,], CAL_ESS[i], CAL_nsamp[i], 
      CN[i,,], Len_age[i,,], LatASD[i,,], truncSD=0) 
    LFC[i] <- CAL_binsmid[min(which(round(CAL[i,nyears, ],0) > 1))] # get the smallest CAL observation	  
  }
  
  Ierr <- array(rlnorm((nyears + proyears) * nsim, mconv(1, rep(Isd, 
    nyears + proyears)), sdconv(1, rep(Isd, nyears + proyears))), c(nsim, 
    nyears + proyears))
  II <- (apply(Biomass, c(1, 3), sum) * Ierr[, 1:nyears])^betas  # apply hyperstability / hyperdepletion
  II <- II/apply(II, 1, mean)  # normalize
  
  message("Calculating MSY reference points")  # Print a progress update
  flush.console()  # update the console

  if (snowfall::sfIsRunning()) {
    snowfall::sfExport(list = c("Marray", "hs", "Mat_age", "Wt_age", "R0", "V", "nyears", "maxage"))  # export some newly made arrays to the cluster
    # MSYrefs <- snowfall::sfSapply(1:nsim, getFMSY, Marray, hs, Mat_age, Wt_age, 
      # R0, V = V[, , nyears], maxage, nyears, proyears = 200, Spat_targ, 
      # mov, SRrel, aR, bR)  # optimize for MSY reference points\t
	## Above version didn't include proyears in projections

    # Using Rcpp code 	
    MSYrefs <- snowfall::sfSapply(1:nsim, getFMSY2, Marray, hs, Mat_age, Wt_age, 
      R0, V = V, maxage, nyears, proyears = 200, Spat_targ, 
      mov, SRrel, aR, bR)  # optimize for MSY reference points\t	  
  } else {
    # MSYrefs_R <- sapply(1:nsim, getFMSY, Marray, hs, Mat_age, Wt_age, 
      # R0, V = V[, , nyears], maxage, nyears, proyears = 200, Spat_targ, 
      # mov, SRrel, aR, bR)  # optimize for MSY reference points
    # Using Rcpp code 	  
    MSYrefs <- sapply(1:nsim, getFMSY2, Marray, hs, Mat_age, Wt_age, 
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
  SSBMSY_SSB0 <- MSYrefs[4, ] # SSBMSY relative to unfished (SSB) 
  BMSY_B0 <- MSYrefs[5, ] # Biomass relative to unfished (B0)
  
  VBMSY <- (MSY/(1 - exp(-FMSY)))  # Biomass at MSY (Vulnerable)
  FMSYb <- MSYrefs[8,]  # instantaneous FMSY (Spawning Biomass)
  UMSY <- MSY/VBMSY  # exploitation rate [equivalent to 1-exp(-FMSY)]
  FMSY_M <- FMSY/M  # ratio of true FMSY to natural mortality rate M
  
  message("Calculating reference yield - best fixed F strategy")  # Print a progress update
  flush.console()  # update the console
  if (snowfall::sfIsRunning()) {
    # Numerically optimize for F that provides highest long term yield
    # RefY_R <- snowfall::sfSapply(1:nsim, getFref, Marray = Marray, Wt_age = Wt_age, 
        # Mat_age = Mat_age, Perr = Perr, N_s = N[, , nyears, ], SSN_s = SSN[, 
        # , nyears, ], Biomass_s = Biomass[, , nyears, ], VBiomass_s = VBiomass[, , nyears, ], 
		# SSB_s = SSB[, , nyears, ], Vn = V[, , (nyears + 1):(nyears + proyears)], 
		# hs = hs, R0a = R0a, nyears = nyears, proyears = proyears, nareas = nareas, 
		# maxage = maxage, mov = mov, SSBpR = SSBpR, aR = aR, bR = bR, SRrel = SRrel, Spat_targ = Spat_targ)
    RefY <- snowfall::sfSapply(1:nsim, getFref2, Marray = Marray, Wt_age = Wt_age, 
      Mat_age = Mat_age, Perr = Perr, N_s = N[, , nyears, , drop=FALSE], SSN_s = SSN[, , nyears, , drop=FALSE], 
	  Biomass_s = Biomass[, , nyears, , drop=FALSE], VBiomass_s = VBiomass[, , nyears, , drop=FALSE], 
	  SSB_s = SSB[, , nyears, , drop=FALSE], Vn = V[, , (nyears + 1):(nyears + proyears), drop=FALSE], 
	  hs = hs, R0a = R0a, nyears = nyears, proyears = proyears, nareas = nareas,
	  maxage = maxage, mov = mov, SSBpR = SSBpR, aR = aR, bR = bR, SRrel = SRrel, Spat_targ = Spat_targ)
	
  } else {
    # RefY <- sapply(1:nsim, getFref, Marray = Marray, Wt_age = Wt_age, 
      # Mat_age = Mat_age, Perr = Perr, N_s = N[, , nyears, ], SSN_s = SSN[, , nyears, ], 
	  # Biomass_s = Biomass[, , nyears, ], VBiomass_s = VBiomass[, , nyears, ], 
	  # SSB_s = SSB[, , nyears, ], Vn = V[, , (nyears + 1):(nyears + proyears)], 
	  # hs = hs, R0a = R0a, nyears = nyears, proyears = proyears, nareas = nareas,
	  # maxage = maxage, mov = mov, SSBpR = SSBpR, aR = aR, bR = bR, SRrel = SRrel, Spat_targ = Spat_targ)
	  
    RefY <- sapply(1:nsim, getFref2, Marray = Marray, Wt_age = Wt_age, 
      Mat_age = Mat_age, Perr = Perr, N_s = N[, , nyears, , drop=FALSE], SSN_s = SSN[, , nyears, , drop=FALSE], 
	  Biomass_s = Biomass[, , nyears, , drop=FALSE], VBiomass_s = VBiomass[, , nyears, , drop=FALSE], 
	  SSB_s = SSB[, , nyears, , drop=FALSE], Vn = V[, , (nyears + 1):(nyears + proyears), drop=FALSE], 
	  hs = hs, R0a = R0a, nyears = nyears, proyears = proyears, nareas = nareas,
	  maxage = maxage, mov = mov, SSBpR = SSBpR, aR = aR, bR = bR, SRrel = SRrel, Spat_targ = Spat_targ)

  }
  
  # LFS<-Linf*(1-exp(-K*(mod-t0))) # Length at full selection
  if (nsim > 1) A <- apply(VBiomass[, , nyears, ], 1, sum)  # Abundance
  if (nsim == 1) A <- sum(VBiomass[, , nyears, ])  # Abundance
  if (nsim > 1) Asp <- apply(SSB[, , nyears, ], 1, sum)  # SSB Abundance
  if (nsim == 1) Asp <- sum(SSB[, , nyears, ])  # SSB Abundance  
  
  OFLreal <- A * FMSY  # the true simulated Over Fishing Limit
  Recerr <- array(rlnorm((nyears + proyears) * nsim, mconv(1, rep(Recsd, (nyears + proyears))), 
    sdconv(1, rep(Recsd, nyears + proyears))), c(nsim, nyears + proyears))

  ntest <- 20  # number of trials  
  BMSY_B0bias <- array(rlnorm(nsim * ntest, mconv(1, OM@BMSY_B0cv), sdconv(1, OM@BMSY_B0cv)), dim = c(nsim, ntest))  # trial samples of BMSY relative to unfished  
  # test <- array(BMSY_B0 * BMSY_B0bias, dim = c(nsim, ntest))  # the simulated observed BMSY_B0 
  test <- array(SSBMSY_SSB0 * BMSY_B0bias, dim = c(nsim, ntest))  # the simulated observed BMSY_B0 
  indy <- array(rep(1:ntest, each = nsim), c(nsim, ntest))  # index
  indy[test > 0.9] <- NA  # interval censor
  BMSY_B0bias <- BMSY_B0bias[cbind(1:nsim, apply(indy, 1, min, na.rm = T))]  # sample such that BMSY_B0<90%
  
  I3 <- apply(Biomass, c(1, 3), sum)^betas  # apply hyperstability / hyperdepletion
  I3 <- I3/apply(I3, 1, mean)  # normalize index to mean 1
  # Iref <- apply(I3[, 1:5], 1, mean) * BMSY_B0  # return the real target abundance index corresponding to BMSY
  if (nsim > 1) Iref <- apply(I3[, 1:5], 1, mean) * SSBMSY_SSB0  # return the real target abundance index corresponding to BMSY
  if (nsim == 1) Iref <- mean(I3[1:5]) * SSBMSY_SSB0

  hsim <- rep(NA, nsim)  # simulate values in steepness 
  cond <- hs > 0.6
  hsim[cond] <- 0.2 + rbeta(sum(hs > 0.6), alphaconv((hs[cond] - 0.2)/0.8, 
    (1 - (hs[cond] - 0.2)/0.8) * OM@hcv), betaconv((hs[cond] - 0.2)/0.8, 
    (1 - (hs[cond] - 0.2)/0.8) * OM@hcv)) * 0.8
  hsim[!cond] <- 0.2 + rbeta(sum(hs < 0.6), alphaconv((hs[!cond] - 0.2)/0.8, 
    (hs[!cond] - 0.2)/0.8 * OM@hcv), betaconv((hs[!cond] - 0.2)/0.8, 
    (hs[!cond] - 0.2)/0.8 * OM@hcv)) * 0.8
  hbias <- hsim/hs  # back calculate the simulated bias
  if (OM@hcv == 0) hbias <- rep(1, nsim) 
  qmu <- -0.5 * qcv^2  # Mean
  qvar <- array(exp(rnorm(proyears * nsim, rep(qmu, proyears), rep(qcv, proyears))), c(nsim, proyears))  # Variations in interannual variation
  # colnames(qvar) <- paste0('qvar', 1:proyears)
  FinF <- Find[, nyears]  # Effort in final historical year
  
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
  Data@vbLinf <- Linf * Linfbias
  Data@L50 <- L50 * lenMbias
  Data@L95 <- L95 * lenMbias
  Data@L95[Data@L95 > 0.9 * Data@vbLinf] <- 0.9 * Data@vbLinf[Data@L95 > 
    0.9 * Data@vbLinf]  # Set a hard limit on ratio of L95 to Linf
  Data@L50[Data@L50 > 0.9 * Data@L95] <- 0.9 * Data@L95[Data@L50 > 
    0.9 * Data@L95]  # Set a hard limit on ratio of L95 to Linf
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
  Data@wla <- rep(OM@a, nsim)
  Data@wlb <- rep(OM@b, nsim)
  # Data@OM <- as.data.frame(cbind(RefY, M, Depletion, A, BMSY_B0, 
  Data@OM <- as.data.frame(cbind(RefY, M, Depletion, A, SSBMSY_SSB0, 
    FMSY_M, Mgrad, Msd, procsd, Esd, dFfinal, MSY, qinc, qcv, FMSY, 
    Linf, K, t0, hs, Linfgrad, Kgrad, Linfsd, recgrad, Ksd, ageM, L5[nyears, ], 
	LFS[nyears, ], Vmaxlen[nyears, ], LFC, OFLreal, Spat_targ, 
    Frac_area_1, Prob_staying, AC, L50, L95, B0, N0, SSB0, BMSY_B0))  # put all the operating model parameters in one table
  
  names(Data@OM)[26:28] <- c("L5", "LFS", "Vmaxlen")  # These are missing labels in the line above
  
  Data@Obs <- as.data.frame(cbind(Cbias, Csd, CAA_nsamp, CAA_ESS, 
    CAL_nsamp, CAL_ESS, Isd, Dbias, Derr, Mbias, FMSY_Mbias, BMSY_B0bias, 
    lenMbias, LFCbias, LFSbias, Abias, Aerr, Kbias, t0bias, Linfbias, 
    hbias, Irefbias, Crefbias, Brefbias, betas))  # put all the observation error model parameters in one table
  
  Data@LHYear <- OM@nyears  # Last historical year is nyears (for fixed MPs)
  Data@MPrec <- Cobs[, nyears]
  Data@MPeff <- rep(1, nsim)
  Data@Misc <- vector("list", nsim)

  ## Write custompars ##
  # datout <- as.data.frame(cbind(procsd,AC,M,Msd,Mgrad,hs,Linf,Linfsd,Linfgrad,recgrad,K,Ksd,Kgrad,t0,L50,L95,L5,LFS,
  # Vmaxlen,Spat_targ,Frac_area_1,Prob_staying,Size_area_1,Csd,Cbias,CAA_nsamp,CAA_ESS,CALcv,betas,
  # Isd,Derr,Dbias,Mbias,FMSY_Mbias,lenMbias,LFCbias,LFSbias,Aerr,Abias,Kbias,
  # t0bias,Linfbias,Irefbias,Crefbias,Brefbias,Recsd,qinc,qcv))
  
  
  ## Return Historical Simulations and Data from last historical year ##
  if (Hist) { # Stop the model after historical simulations are complete
  	message("Returning historical simulations")
	nout <- t(apply(N, c(1, 3), sum))
	vb <- t(apply(VBiomass, c(1, 3), sum))
	b <- t(apply(Biomass, c(1, 3), sum))
	ssb <- t(apply(SSB, c(1, 3), sum))
    Cc <- t(apply(CB, c(1,3), sum))
	rec <- t(apply(N[, 1, , ], c(1,2), sum))
	
    TSdata <- list(VB=vb, SSB=ssb, Bio=b, Catch=Cc, Rec=rec, N=nout)
    AtAge <- list(Len_age=Len_age, Wt_age=Wt_age, 
	  Sl_age=V, Mat_age=Mat_age, Nage=apply(N, c(1:3), sum), SSBage=apply(SSB, c(1:3), sum))
    MSYs <- list(MSY=MSY, FMSY=FMSY, FMSYb=FMSYb, VBMSY=VBMSY, UMSY=UMSY, 
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
	  Mrand=Mrand, Linfrand=Linfrand, Krand=Krand, maxage=maxage, V=V, Depletion=Depletion,qs=qs) 

	HistData <- list(SampPars=SampPars, TSdata=TSdata, AtAge=AtAge, MSYs=MSYs, Data=Data)
	class(HistData) <- c("list", "hist")
	return(HistData)	
  }

  # assign('Data',Data,envir=.GlobalEnv) # for debugging fun
   
  # Run projections
  # ===========================================================================
  
  if (is.na(MPs[1])) CheckMPs <- TRUE
  if (CheckMPs) {
    message("Determining available methods")  # print an progress report
    flush.console()  # update the console
    PosMPs <- Can(Data, timelimit = timelimit)  # list all the methods that could be applied to a DLM data object 
    if (is.na(MPs[1])) {
      MPs <- PosMPs  # if the user does not supply an argument MPs run the MSE or all available methods
      message("No MPs specified: running all available")
    }
    if (!is.na(MPs[1])) MPs <- MPs[MPs %in% PosMPs]  # otherwise run the MSE for all methods that are deemed possible
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
  # SPRa <- array(NA,dim=c(nsim,nMP,proyears)) # store the Spawning Potential Ratio
  
  MPdur <- rep(NA, nMP)
  mm <- 1 # for debugging
  for (mm in 1:nMP) {
    # MSE Loop over methods
    pL5 <- L5  # reset selectivity parameters for projections
    pLFS <- LFS
    pVmaxlen <- Vmaxlen
	pSLarray <- SLarray # selectivity at length array
    
    message(paste(mm, "/", nMP, " Running MSE for ", MPs[mm], sep = ""))  # print a progress report
    flush.console()  # update the console
    
    # projection arrays
    N_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    Biomass_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    VBiomass_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    SSN_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
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
      N_P[, 1, 1, ] <- Perr[, nyears] * (0.8 * R0a * hs * apply(SSB[, 
        , nyears, ], c(1, 3), sum))/(0.2 * SSBpR * R0a * (1 - hs) + 
        (hs - 0.2) * apply(SSB[, , nyears, ], c(1, 3), sum))  # Recruitment assuming regional R0 and stock wide steepness
    } else {
      # most transparent form of the Ricker uses alpha and beta params
      N_P[, 1, 1, ] <- Perr[, nyears] * aR * apply(SSB[, , nyears, ], c(1, 3), sum) * exp(-bR * apply(SSB[, , nyears, ], c(1, 3), sum))
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
      st <- Sys.time()
      Data <- Sam(MSElist[[mm]], MPs = MPs[mm], perc = pstar, reps = reps)
      nd <- Sys.time()
      MPdur[mm] <- nd - st
      TACused <- apply(Data@TAC, 3, quantile, p = pstar, na.rm = T)
      TACa[, mm, 1] <- TACused
	  availB <- apply(VBiomass_P[,,1,], 1, sum) # total available biomass
	  maxC <- (1 - exp(-maxF)) * availB
      # if the TAC is higher than maxC than catch is equal to maxC
	  notNA <- which(!is.na(TACused) & !is.na(availB)) # robustify for MPs that return NA 
	  TACused[TACused[notNA] > maxC[notNA]] <- maxC[TACused[notNA] > maxC[notNA]]
	   
      fishdist <- (apply(VBiomass_P[, , 1, ], c(1, 3), sum)^Spat_targ)/
	    apply(apply(VBiomass_P[, , 1, ], c(1, 3), sum)^Spat_targ, 1, mean)  # spatial preference according to spatial biomass
      
      CB_P[SAYR] <- Biomass_P[SAYR] * (1 - exp(-V_P[SAYt] * fishdist[SR]))  # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space
      
      temp <- CB_P[, , 1, ]/apply(CB_P[, , 1, ], 1, sum)  # how catches are going to be distributed
      CB_P[, , 1, ] <- TACused * temp  # debug - to test distribution code make TAC = TAC2, should be identical
       
      temp <- CB_P[SAYR]/(Biomass_P[SAYR] * exp(-Marray[SYt]/2))  # Pope's approximation	  
      temp[temp > (1 - exp(-maxF))] <- 1 - exp(-maxF)
      FM_P[SAYR] <- -log(1 - temp)
	  Z_P[SAYR] <- FM_P[SAYR] + Marray[SYt]
		  
      Effort[, mm, y] <- (-log(1 - apply(CB_P[, , y, ], 1, sum)/(apply(CB_P[, , y, ], 1, sum) + 
	    apply(VBiomass_P[, , y, ], 1, sum))))/qs	  
    } else {
      # input control
      st <- Sys.time()
      runIn <- runInMP(MSElist[[mm]], MPs = MPs[mm], reps = reps)  # Apply input control MP
      nd <- Sys.time()
      MPdur[mm] <- nd - st
      
      inc <- runIn[[1]]
      Data <- runIn[[2]]
      
      Ai <- inc[1, , 1]
      Ei <- inc[2, , 1]
      Effort[, mm, y] <- Ei  # Change in Effort
      Si <- t(inc[3:4, , 1])
      newSel <- inc[5:6, , 1]
      
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
	  
      Vi <- t(sapply(1:nsim, SelectFun, pL5[y + nyears, ], pLFS[y + nyears, ], 
	    pVmaxlen[y + nyears, ], Len_age[, maxage, nyears], Len_age[, , y + nyears])) # update vulnerability-at-age schedule 
      
	  ind <- as.matrix(expand.grid(1:nsim, 1:length(CAL_binsmid), (y+nyears):(nyears+proyears)))
      pSLarray[ind] <- t(sapply(1:nsim, SelectFun, SL0.05=pL5[y+nyears, ], SL1=pLFS[y+nyears, ], 
	                             MaxSel=pVmaxlen[y+nyears, ], maxlens=maxlen, Lens=CAL_binsmid)) # update vulnerability-at-length schedule 
								 
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
          fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ, 
          1, mean)  # spatial preference according to spatial biomass
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
	  Z_P[SAYR] <- FM_P[SAYR] + Marray[SYt] # calculate total mortality 
      CB_P[SAYR] <- FM_P[SAYR]/Z_P[SAYR] * Biomass_P[SAYR] * (1 - exp(-Z_P[SAYR]))  	   
    }  # input control  
    
  
    # CB_P[SAYR] <- Biomass_P[SAYR]*(1-exp(-FM_P[SAYR]))
    
    
    # TACa[, mm, 1] <- apply(CB_P[, , 1, ], 1, sum)  # Adjust TAC to actual catch in the year 
    # To account for years where TAC is higher than catch
    
    upyrs <- 1 + (0:(floor(proyears/interval) - 1)) * interval  # the years in which there are updates (every three years)
    cat(".")
    flush.console()
    
    for (y in 2:proyears) {
      cat(".")
      flush.console()
      if (class(match.fun(MPs[mm])) == "Output")  TACa[, mm, y] <- TACa[, mm, y-1] # TAC same as last year unless changed 
      SAYRt <- as.matrix(expand.grid(1:nsim, 1:maxage, y + nyears, 
        1:nareas))  # Trajectory year
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
      SA1YR <- as.matrix(expand.grid(1:nsim, 1:(maxage - 1), y - 
        1, 1:nareas))
      indMov <- as.matrix(expand.grid(1:nareas, 1:nareas, y, 1:maxage, 
        1:nsim)[5:1])
      indMov2 <- indMov[, c(1, 2, 3, 4)]
      indMov3 <- indMov[, c(1, 4, 5)]
      
      N_P[SA2YR] <- N_P[SA1YR] * exp(-Z_P[SA1YR])  # Total mortality
      if (SRrel[1] == 1) {
        N_P[, 1, y, ] <- Perr[, y + nyears] * (0.8 * R0a * hs * 
          apply(SSB_P[, , y - 1, ], c(1, 3), sum))/(0.2 * SSBpR * 
          R0a * (1 - hs) + (hs - 0.2) * apply(SSB_P[, , y - 1, ], c(1, 3), sum))  # Recruitment assuming regional R0 and stock wide steepness
      } else {
        # most transparent form of the Ricker uses alpha and beta params
        N_P[, 1, y, ] <- Perr[, y + nyears] * aR *
  		  apply(SSB_P[, , y - 1, ], c(1, 3), sum) * exp(-bR * apply(SSB_P[, , y - 1, ], c(1, 3), sum))
      }
      
      temp <- array(N_P[indMov2] * mov[indMov3], 
	    dim = c(nareas, nareas, maxage, nsim))  # Move individuals
      N_P[, , y, ] <- apply(temp, c(4, 3, 1), sum)
      
      Biomass_P[SAYR] <- N_P[SAYR] * Wt_age[SAYt]  # Calculate biomass
      VBiomass_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt]  # Calculate vulnerable biomass
      SSN_P[SAYR] <- N_P[SAYR] * Mat_age[SA]  # Calculate spawning stock numbers
      SSB_P[SAYR] <- SSN_P[SAYR] * Wt_age[SAYt]  # Calculate spawning stock biomass
      
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
		cn <- as.matrix(CNtemp[i,,])
		if (interval == 1) cn <- t(cn) # dodgy hack to ensure matrix is correct 
        for (i in 1:nsim) { # Rcpp code 
          CAL[i, 1:interval, ] <- genLenComp(CAL_bins, CAL_binsmid, as.matrix(pSLarray[i,, nyears + yind]), CAL_ESS[i], CAL_nsamp[i], 
            cn, as.matrix(Len_age[i,,nyears + yind]), as.matrix(LatASD[i,, nyears + yind]), truncSD=0) 
          LFC[i] <- CAL_binsmid[min(which(round(CAL[i, interval, ],0) > 1))] # get the smallest CAL observation	
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

        
        I2 <- cbind(apply(Biomass, c(1, 3), sum), apply(Biomass_P, 
          c(1, 3), sum)[, 1:(y - 1)]) * Ierr[, 1:(nyears + (y - 
          1))]^betas
        I2[is.na(I2)] <- tiny
        I2 <- I2/apply(I2, 1, mean)
        
        # Depletion <- apply(Biomass_P[, , y, ], 1, sum)/apply(Biomass[, , 1, ], 1, sum)
		Depletion <- apply(SSB_P[, , y, ], 1, sum)/apply(SSB[, , 1, ], 1, sum)
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
          TACused <- apply(Data@TAC, 3, quantile, p = pstar, 
          na.rm = TRUE)  #
          NAs <- which(is.na(TACused))
          if (length(NAs) > 0) {
          # robustifying TAC setting!
          TACused[NAs] <- TACa[NAs, mm, y - 1]  #
          if (!exists("store")) 
            store <- list()
            store <- append(store, c(MPs[mm], NAs))
          }
          TACa[, mm, y] <- TACused
          MSElist[[mm]]@MPrec <- TACused
		  
		  availB <- apply(VBiomass_P[,,y,], 1, sum) # total available biomass
	      maxC <- (1 - exp(-maxF)) * availB
          # if the TAC is higher than maxC than catch is equal to maxC
		  notNA <- which(!is.na(TACused) & !is.na(availB))
	      TACused[TACused[notNA] > maxC[notNA]] <- maxC[TACused[notNA] > maxC[notNA]]
		  # TACused[TACused > maxC] <- maxC[TACused > maxC] 	
		  
		  fishdist <- (apply(VBiomass_P[, , y, ], c(1, 3), sum)^Spat_targ)/apply(apply(VBiomass_P[, , y, ], c(1, 3), sum)^Spat_targ, 1, mean)  # spatial preference according to spatial biomass     
          CB_P[SAYR] <- Biomass_P[SAYR] * (1 - exp(-V_P[SAYt] *  fishdist[SR]))  # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space          
          temp <- CB_P[, , y, ]/apply(CB_P[, , y, ], 1, sum)  # how catches are going to be distributed
          CB_P[, , y, ] <- TACused * temp  # debug - to test distribution code make TAC = TAC2, should be identical          
          temp <- CB_P[SAYR]/(Biomass_P[SAYR] * exp(-Marray[SYt]/2))  # Pope's approximation
          temp[temp > (1 - exp(-maxF))] <- 1 - exp(-maxF)
          FM_P[SAYR] <- -log(1 - temp)
          Z_P[SAYR] <- FM_P[SAYR] + Marray[SYt]
          Effort[, mm, y] <- (-log(1 - apply(CB_P[, , y, ], 1, sum)/(apply(CB_P[, , y, ], 1, sum) + apply(VBiomass_P[, , y, ], 1, sum))))/qs
           
        } else {
          MSElist[[mm]]@MPeff <- Ei
          runIn <- runInMP(MSElist[[mm]], MPs = MPs[mm], reps = reps)  # Apply input control MP
          inc <- runIn[[1]]
          Data <- runIn[[2]]
          Ai <- inc[1, , 1]
          Ei <- inc[2, , 1]
          Effort[, mm, y] <- Ei  # Change in Effort
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

        Vi <- t(sapply(1:nsim, SelectFun, pL5[y + nyears, ], pLFS[y + nyears, ], 
           pVmaxlen[y + nyears, ], Len_age[, maxage, nyears], Len_age[, , y + nyears])) # update vulnerability-at-age schedule 
        
        ind <- as.matrix(expand.grid(1:nsim, 1:length(CAL_binsmid), (y+nyears):(nyears+proyears)))
        pSLarray[ind] <- t(sapply(1:nsim, SelectFun, SL0.05=pL5[y+nyears, ], SL1=pLFS[y+nyears, ], 
                                MaxSel=pVmaxlen[y+nyears, ], maxlens=maxlen, Lens=CAL_binsmid)) # update vulnerability-at-length schedule 
   							 
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
              fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ, 
              1, mean)  # spatial preference according to spatial biomass
              FM_P[SAYR] <- FinF[S1] * Ei[S1] * V_P[SAYt] * fishdist[SR] * 
                            qvar[SY] * qs[S1] * (1 + qinc[S1]/100)^y  # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
            } else {
              if (y < proyears) 
              V_P[, , (y + nyears + 1):(proyears + nyears)] <- Vi  # Update vulnerability schedule for all future years
              newVB <- apply(VBiomass_P[, , y, ] * Vi[SA], c(1, 
              3), sum)  # vulnerability modified
              fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ, 1, mean)  # spatial preference according to spatial biomass
              FM_P[SAYR] <- FinF[S1] * Ei[S1] * Vi[SA] * fishdist[SR] *  qvar[SY] * 
			                qs[S1] * (1 + qinc[S1]/100)^y  # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass         
            }
          } else {
          # A spatial closure if no vulnerability schedule is specified
            if (!Vchange) {
              newVB <- apply(VBiomass_P[, , y, ], c(1, 3), sum)  # vulnerability isn't changed
              fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ,1, mean)  # spatial preference according to spatial biomass
              Emult <- 1 + ((2/apply(fishdist * Si, 1, sum)) - 1) * Ai  # allocate effort to new area according to fraction allocation Ai
              FM_P[SAYR] <- FinF[S1] * Ei[S1] * V_P[SAYt] * Si[SR] * fishdist[SR] * 
			                Emult[S1] * qvar[SY] * qs[S1] * (1 + qinc[S1]/100)^y        
            } else {
              if (y < proyears) 
              V_P[, , (y + nyears + 1):(proyears + nyears)] <- Vi  # Update vulnerability schedule for all future years
              newVB <- apply(VBiomass_P[, , y, ] * Vi[SA], c(1, 
              3), sum)  # vulnerability modified
              fishdist <- (newVB^Spat_targ)/apply(newVB^Spat_targ, 
              1, mean)  # spatial preference according to spatial biomass
              Emult <- 1 + ((2/apply(fishdist * Si, 1, sum)) - 
              1) * Ai  # allocate effort to new area according to fraction allocation Ai
              FM_P[SAYR] <- FinF[S1] * Ei[S1] * Vi[SA] * Si[SR] * 
              fishdist[SR] * Emult[S1] * qvar[SY] * qs[S1] * 
              (1 + qinc[S1]/100)^y
              
            }  #vuln not changed
          }  # spatial closure
          VBiomass_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt]  # update vulnerable biomass 
          Z_P[SAYR] <- FM_P[SAYR] + Marray[SYt]
          # CB_P[SAYR]<-Biomass_P[SAYR]*(1-exp(-FM_P[SAYR]))
          CB_P[SAYR] <- FM_P[SAYR]/Z_P[SAYR] * Biomass_P[SAYR] *   (1 - exp(-Z_P[SAYR]))		  
        }  # input or output control 

        # TACused <- apply(CB_P[, , y, ], 1, sum)  # Set last years TAC to actual catch from last year
        # TACa[, mm, y] <- TACused
		tempcatch <- apply(CB_P[, , y, ], 1, sum) 
		
        MSElist[[mm]]@MPrec <- tempcatch
      } else {
        # not an update yr
        vbio <- apply(VBiomass_P[, , y, ], c(1, 3), sum)
        fishdist <- (vbio^Spat_targ)/apply(vbio^Spat_targ, 1, mean)  # calculate distribution of effort \t  
        if (class(match.fun(MPs[mm])) == "Output") {
          CB_P[SAYR] <- Biomass_P[SAYR] * (1 - exp(-fishdist[SR] *  V_P[SAYt]))  # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space
          temp <- CB_P[, , y, ]/apply(CB_P[, , y, ], 1, sum)  # how catches are going to be distributed
          tempcatch <- TACa[, mm, y-1]
		 
	      availB <- apply(VBiomass_P[,,y,], 1, sum) # total available biomass
	      maxC <- (1 - exp(-maxF)) * availB
          # if the TAC is higher than maxC than catch is equal to maxC
	      notNA <- which(!is.na(tempcatch) & !is.na(availB))		  
	      tempcatch[tempcatch[notNA] > maxC[notNA]] <- maxC[tempcatch[notNA] > maxC[notNA]]		  
	      # tempcatch[tempcatch > maxC] <- maxC[tempcatch > maxC] 		 
  		  
		  CB_P[, , y, ] <- tempcatch * temp  # debug - to test distribution code make TAC = TAC2, should be identical
          temp <- CB_P[SAYR]/(Biomass_P[SAYR] * exp(-Marray[SYt]/2))  # Pope's approximation
          temp[temp > (1 - exp(-maxF))] <- 1 - exp(-maxF)
          FM_P[SAYR] <- -log(1 - temp)
          Effort[, mm, y] <- (-log(1 - apply(CB_P[, , y, ], 1, sum)/
		                     (apply(CB_P[, , y, ], 1, sum) + apply(VBiomass_P[, , y, ], 1, sum))))/qs
          Z_P[SAYR] <- FM_P[SAYR] + Marray[SYt]							 
        } else {
          # input control FM_P[SAYR] <- FM_P[SAY1R]*qvar[SY] *(1+qinc[S1]/100)^y
          # # add fishing efficiency changes and variability
          FM_P[SAYR] <- FM_P[SAY1R] * qvar[SY] * (1 + qinc[S1]/100)  # add fishing efficiency changes and variability
          Effort[, mm, y] <- Effort[, mm, y - 1]  # Effort doesn't change in non-update year
		  Z_P[SAYR] <- FM_P[SAYR] + Marray[SYt]
          # CB_P[SAYR]<-Biomass_P[SAYR]*(1-exp(-FM_P[SAYR]))
          CB_P[SAYR] <- FM_P[SAYR]/Z_P[SAYR] * Biomass_P[SAYR] * (1 - exp(-Z_P[SAYR]))
        }
      
      }  # not an update year
      
    }  # end of year
       
    B_BMSYa[, mm, ] <- apply(SSB_P, c(1, 3), sum, na.rm=TRUE)/SSBMSY  # SSB relative to SSBMSY  
    # F_FMSYa[, mm, ] <- (-log(1 - apply(CB_P, c(1, 3), sum)/(apply(CB_P, c(1, 3), sum) + 
	                    # apply(VBiomass_P, c(1, 3), sum))))/FMSY 
    # VBiomass is calculated before catches are taken
    suppressWarnings(	# gives an error message if CB_P or VBiomass_P is NA 
	FMa[, mm, ] <- -log(1 - apply(CB_P, c(1, 3), sum, na.rm=TRUE)/apply(VBiomass_P, c(1, 3), sum, na.rm=TRUE))		
	)
	F_FMSYa[, mm, ] <- FMa[, mm, ]/FMSY
	                    	
    Ba[, mm, ] <- apply(Biomass_P, c(1, 3), sum, na.rm=TRUE) # biomass 
	SSBa[, mm, ] <- apply(SSB_P, c(1, 3), sum, na.rm=TRUE) # spawning stock biomass
	VBa[, mm, ] <- apply(VBiomass_P, c(1, 3), sum, na.rm=TRUE) # vulnerable biomass
    # FMa[, mm, ] <- -log(1 - apply(CB_P, c(1, 3), sum)/(apply(CB_P, c(1, 3), sum) + 
	               # apply(VBiomass_P, c(1, 3), sum)))
    # VBiomass is calculated before catches are taken 				   
	
    Ca[, mm, ] <- apply(CB_P, c(1, 3), sum, na.rm=TRUE)
    cat("\n")
  }  # end of mm methods

  # Store MP duration
  attr(MPs, "duration") <- MPdur

  MSEout <- new("MSE", Name = OM@Name, nyears, proyears, nMPs=nMP, MPs, nsim, 
    Data@OM, Obs=Data@Obs, B_BMSY=B_BMSYa, F_FMSY=F_FMSYa, B=Ba, 
	SSB=SSBa, VB=VBa, FM=FMa, Ca, TAC=TACa, SSB_hist = SSB, CB_hist = CB, 
	FM_hist = FM, Effort = Effort)
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