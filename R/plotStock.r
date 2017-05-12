#' Wrapper for histogram function to remove annoying 
hist2 <- function(x, col, ...) {
  if (mean(x) == x[1]) {
    
    hist(x, col="white", border=FALSE, xlim=range(x), ...)
    abline(v=x)
    
  } else {
    hist(x, col=col, ...)
  }
}


#' Plot the Stock object parameters 
#' 
#' A function that plots histograms of samples from the Stock object parameters,
#' and time-series plots of `nsamp` samples of time-varying parameters. Used to 
#' visually examine the parameter values and ranges entered into the Stock object.
#' 
#' @param x An object of class Stock (or of class OM) 
#' @param nsamp Number of random samples for time-series plots
#' @param nsim Number of iterations for histograms
#' @param nyears Number of historical years
#' @param proyears Number of projection years 
#' @param col Color of histograms 
#' @param breaks Number of breaks for histograms 
#' @param lwd line width 
#' @param ask Ask before displaying next page?
#' @param incVB Show the sampled von Bertalanffy growth curves on second page?
#' @param ...  Optional additional arguments passed to \code{plot}
#' @rdname plot-Stock
#' @method plot Stock
#' @author A. Hordyk
#' @export 
plot.Stock <- function(x, nsamp=3, nsim=500, nyears=50, proyears=28, 
  col="darkgray", breaks=10, lwd=2, ask=FALSE, incVB=TRUE, ...) {
  Stock <- x 
  cpars <- NULL
  if (class(Stock) == "OM") {
    if (is.finite(Stock@nyears)) nyears <- Stock@nyears
	if (is.finite(Stock@proyears)) proyears <- Stock@proyears
	if (is.finite(Stock@nsim)) nsim <- Stock@nsim	
	if (length(Stock@cpars) > 0) {
	  cpars <- Stock@cpars
	  ncparsim <- cparscheck(cpars)
	}
	Stock <- SubOM(Stock)
  }
  its <- sample(1:nsim, nsamp)
	
  maxage <- Stock@maxage  # maximum age (no plus group)
  calcMax <- -log(0.01)/(min(Stock@M))        # Age at which 1% of cohort survives
  maxage <- round(max(maxage, calcMax),0)  # If maximum age is lower, increase it to calcMax
 
  ## Life History Parameters ##
  procsd <- runif(nsim, Stock@Perr[1], Stock@Perr[2])  # Process error standard deviation
  AC <- runif(nsim, Stock@AC[1], Stock@AC[2])  # auto correlation parameter for recruitment deviations recdev(t)<-AC*recdev(t-1)+(1-AC)*recdev_proposed(t)
  dep <- runif(nsim, Stock@D[1], Stock@D[2])
  M <- runif(nsim, Stock@M[1], Stock@M[2])  # natural mortality rate \t
  Msd <- runif(nsim, Stock@Msd[1], Stock@Msd[2])  # sample inter annual variability in M frStock specified range
  Mgrad <- runif(nsim, Stock@Mgrad[1], Stock@Mgrad[2])  # sample gradient in M (M y-1)
  hs <- runif(nsim, Stock@h[1], Stock@h[2])  # sample of recruitment cStockpensation (steepness - fraction of unfished recruitment at 20% of unfished biStockass)
  Linf <- runif(nsim, Stock@Linf[1], Stock@Linf[2])  # sample of asymptotic length
  Linfsd <- runif(nsim, Stock@Linfsd[1], Stock@Linfsd[2])  # sample of interannual variability in Linf
  Linfgrad <- runif(nsim, Stock@Linfgrad[1], Stock@Linfgrad[2])  # sample of gradient in Linf (Linf y-1)
  recgrad <- runif(nsim, Stock@recgrad[1], Stock@recgrad[2])  # gradient in recent recruitment
  K <- runif(nsim, Stock@K[1], Stock@K[2])  # now predicted by a log-linear model
  Ksd <- runif(nsim, Stock@Ksd[1], Stock@Ksd[2])  #runif(nsim,Stock@Ksd[1],Stock@Ksd[2])# sd is already added in the linear model prediction
  Kgrad <- runif(nsim, Stock@Kgrad[1], Stock@Kgrad[2])  # gradient in Von-B K parameter (K y-1)
  t0 <- runif(nsim, Stock@t0[1], Stock@t0[2])  # a sample of theoretical age at length zero
  L50 <- array(runif(nsim * 50, Stock@L50[1], Stock@L50[2]), c(nsim, 50))  # length at 50% maturity
  L50_95 <- array(runif(nsim * 50, Stock@L50_95[1], Stock@L50_95[2]), c(nsim, 50))  # length at 95% maturity
  
  # checks for unrealistically high length at maturity 
  L50[L50/Linf > 0.95] <- NA
  L50 <- apply(L50, 1, function(x) x[!is.na(x)][1])
  L50_95[(L50+L50_95)/Linf > 0.99] <- NA
  L50_95 <- apply(L50_95, 1, function(x) x[!is.na(x)][1]) 
  L95 <- array(L50 + L50_95) 
  
    # Generate randStock numbers for randStock walk 
  Mrand <- matrix(exp(rnorm(nsim*(proyears+nyears), -0.5 * Msd^2, Msd)), nrow=nsim, ncol=proyears+nyears)
  Linfrand <- matrix(exp(rnorm(nsim*(proyears+nyears), -0.5 * Linfsd^2, Linfsd)), nrow=nsim, ncol=proyears+nyears)
  Krand <- matrix(exp(rnorm(nsim*(proyears+nyears), -0.5 * Ksd^2, Ksd)), nrow=nsim, ncol=proyears+nyears)
 
  Marray <- gettempvar(M, Msd, Mgrad, nyears + proyears, nsim, Mrand)  # M by sim and year according to gradient and inter annual variability
  Linfarray <- gettempvar(Linf, Linfsd, Linfgrad, nyears + proyears, nsim, Linfrand)  # Linf array
  Karray <- gettempvar(K, Ksd, Kgrad, nyears + proyears, nsim, Krand)  # the K array
  
  Agearray <- array(rep(1:maxage, each = nsim), dim = c(nsim, maxage))  # Age array
  Len_age <- array(NA, dim = c(nsim, maxage, nyears + proyears))  # Length at age array
  ind <- as.matrix(expand.grid(1:nsim, 1:maxage, 1:(nyears + proyears)))  # an index for calculating Length at age
  Len_age[ind] <- Linfarray[ind[, c(1, 3)]] * (1 - exp(-Karray[ind[, c(1, 3)]] * 
    (Agearray[ind[, 1:2]] - t0[ind[, 1]])))

  Wt_age <- array(NA, dim = c(nsim, maxage, nyears + proyears))  # Weight at age array
  Wt_age[ind] <- Stock@a * Len_age[ind]^Stock@b  # Calculation of weight array
  
  # edit for cpars slot
  EffYears <- EffUpper <- EffUpper <- EffLower <- NULL # hack for CRAN checks
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
  
  if (length(cpars) > 0) { # custom parameters exist     
	  Names <- names(cpars)
	  # report not valid names 
	  invalid <- which(!Names %in% ParsNames)
	  if (length(invalid) > 0) {
	    outNames <- paste(Names[invalid], "")
	    for (i in seq(5, by=5, length.out=floor(length(outNames)/5))) outNames <- gsub(outNames[i], paste0(outNames[i], "\n"), outNames)
	    warning("ignoring invalid names found in custom parameters (cpars) \n", outNames)	
	  }
	  # report found names
	  valid <- which(Names %in% ParsNames)
	  cpars <- cpars[valid]
	  if (length(cpars) == 0) stop("No valid names found in custompars")
	  Names <- names(cpars)
	  outNames <- paste(Names, "")
	  for (i in seq(5, by=5, length.out=floor(length(outNames)/5)))
  	  outNames <- gsub(outNames[i], paste0(outNames[i], "\n"), outNames)
	    message("valid custom parameters (cpars) found: \n", outNames)
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

  }
   
  ncol <- 8
 
  m <- layout(matrix(c(c(1, 2, 3, 0, 4, 4, 4, 5),
					   c(6, 7, 8, 0, 9, 9, 9, 14),			   
					   c(10, 11, 12, 0, 13, 13, 13, 0),
					   c(20, 21, 22 ,0, 23, 23, 23, 0),
					   c(15, 16, 17, 0, 18, 18, 19, 19)				   
					   ), ncol=ncol, byrow = TRUE),
					   widths=c(1, 1, 1, 0.5, 1, 1, 1, 1))
									   
  # layout.show(m)
  # stop()
  op <- par(mar = c(2, 1, 3, 1), oma=c(1,2,4,1), ask=FALSE, las=1)
  
  histwrap(M, col=col, axes=FALSE, main="M", breaks=breaks)
  axis(side=1)  
  hist2(Msd, col=col, axes=FALSE, main="Msd", breaks=breaks)
  axis(side=1) 
  hist2(Mgrad, col=col, axes=FALSE, main="Mgrad", breaks=breaks)
  axis(side=1)
  
  # M traj
  matplot(t(Marray[its,]), type="l", bty="l", main="M by Year", lwd=lwd)
  
  # M/K
  hist2(M/K, col=col, axes=FALSE, main="M/K", breaks=breaks)
  axis(side=1)
  
  hist2(Linf, col=col, axes=FALSE, main="Linf", breaks=breaks)
  axis(side=1) 
  hist2(Linfsd, col=col, axes=FALSE, main="Linfsd", breaks=breaks)
  axis(side=1)  
  hist2(Linfgrad, col=col, axes=FALSE, main="Linfgrad", breaks=breaks)
  axis(side=1)
  
  # Linf traj 
  matplot(t(Linfarray[its,]), type="l", bty="l", main="Linf by Year", lwd=lwd)
  
  hist2(K, col=col, axes=FALSE, main="K", breaks=breaks)
  axis(side=1)
  hist2(Ksd, col=col, axes=FALSE, main="Ksd", breaks=breaks)
  axis(side=1) 
  hist2(Kgrad, col=col, axes=FALSE, main="Kgrad", breaks=breaks)
  axis(side=1)  
  
  # K traj 
  matplot(t(Karray[its,]), type="l", bty="l", main="K by Year", lwd=lwd)
    
  hist2(t0, col=col, axes=FALSE, main="t0", breaks=breaks)
  axis(side=1)
  
  # Recruitment Deviations
  procmu <- -0.5 * (procsd)^2  # adjusted log normal mean
  Perr <- array(rnorm((nyears + proyears) * nsim, rep(procmu, nyears + 
    proyears), rep(procsd, nyears + proyears)), c(nsim, nyears + proyears))
  for (y in 2:(nyears + proyears)) Perr[, y] <- AC * Perr[, y - 1] + 
    Perr[, y] * (1 - AC * AC)^0.5  #2#AC*Perr[,y-1]+(1-AC)*Perr[,y] # apply a pseudo AR1 autocorrelation to rec devs (log space)
  Perr <- exp(Perr)  # normal space (mean 1 on average)
  
  # Add cycle (phase shift) to recruitment deviations - if specified
  if (is.finite(Stock@Period[1]) & is.finite(Stock@Amplitude[1])) {
    Shape <- "sin"  # default sine wave - alternative - 'shift' for step changes
    recMulti <- t(sapply(1:nsim, SetRecruitCycle, Period = Stock@Period, 
      Amplitude = Stock@Amplitude, TotYears = nyears + proyears, Shape = Shape))
    Perr <- Perr * recMulti  # Add cyclic pattern to recruitment
  }
  
  # Recruitment 
  hist2(hs, col=col, axes=FALSE, main="Steepness (h)", breaks=breaks)
  axis(side=1)  
  hist2(procsd, col=col, axes=FALSE, main="procsd", breaks=breaks)
  axis(side=1) 
  hist2(AC, col=col, axes=FALSE, main="AC", breaks=breaks)
  axis(side=1)
  
  matplot(t(Perr[its,]), type="l", bty="l", main="Rec Devs by Year", lwd=lwd)
	
  SRrel <- rep(Stock@SRrel, nsim)	
  biomass <- seq(0, 1, length.out=20)

  bR <- matrix(log(5 * hs)/(0.8 * 2), nrow=nsim)  # Ricker SR params
  aR <- matrix(exp(bR * 2), nrow=nsim)  # Ricker SR params

  if (SRrel[1] == 1) {
    recs <- sapply(its, function(X) (0.8  * hs[X] * biomass)/(0.2 * (1 - hs[X]) + 
               (hs[X] - 0.2) * biomass)) # BH SRR            
  } else {
    # most transparent form of the Ricker uses alpha and beta params
    recs <- sapply(its, function(X) aR[X] * biomass * 
	  exp(-bR[X] * biomass))
  }
  matplot(biomass, recs, type="l", bty="l", main="Stock-Recruit", 
    ylim=c(0,1), xlim=c(0,1), axes=FALSE, lwd=lwd)
  axis(side=1)
  axis(side=2, labels=FALSE)
  
  
  # Depletion
  hist2(dep,  col=col, axes=FALSE, main="Depletion", breaks=breaks)
  axis(side=1)
  
  # Maturity 
  
  hist2(L50, col=col, axes=FALSE, main="L50", breaks=breaks)
  axis(side=1)
  hist2(L95, col=col, axes=FALSE, main="L95", breaks=breaks)
  axis(side=1)
  
  slope <- log(19)/(L50_95)
  Ls <- seq(0, to=max(Linf), length.out=200)
 
  Mat_len <- sapply(its, function(X) plogis(Ls, L50[X], 1/slope[X]))
  matplot(Ls, Mat_len, type="l", bty="l", main="Maturity-at-length", lwd=lwd)
  


  title(Stock@Name, outer=TRUE)
  title(paste("nyears =", nyears, "  proyears =", proyears, "  ", nsamp, "sampled iterations"), outer=TRUE, line=0)
  
  if (incVB) {
    par(mfrow=c(1,3), ask=ask, mar=c(5,4,1,1))
    # vB   
	cex.lab <- 2 
    fstYr <- Len_age[its,,1]
    curYr <- Len_age[its,,nyears]
    lstYr <- Len_age[its,,proyears+nyears]
    MaxL <- max(Len_age)
    matplot(t(fstYr), type="l", bty="l", main="First historical year", ylim=c(0, MaxL), xlab="Age", ylab="Length", cex.lab=cex.lab, lwd=lwd)
    matplot(t(curYr), type="l", bty="l", main="Last historical year", ylim=c(0, MaxL),  axes=FALSE, xlab="Age", ylab="", cex.lab=cex.lab, lwd=lwd)
    axis(side=1)
    axis(side=2, labels=FALSE)  
    matplot(t(lstYr), type="l", bty="l", main="Last projected year", ylim=c(0, MaxL), axes=FALSE, xlab="Age", ylab="", cex.lab=cex.lab, lwd=lwd)	
    axis(side=1)
    axis(side=2, labels=FALSE)  
	title("Sampled length-at-age curves", outer=TRUE, cex.main=2)
  }

  on.exit(par(op))
  
}

