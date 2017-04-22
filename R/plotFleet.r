#' Plot the Fleet object parameters 
#' 
#' @param Fleet An object of class Fleet (or of class Fleet) 
#' @param nsamp Number of random samples for time-series plots
#' @param nsim Number of iterations for histograms
#' @param nyears Number of historical years
#' @param proyears Number of projection years 
#' @param col Color of histograms 
#' @param breaks Number of breaks for histograms 
#' @param lwd line width 
#' @author A. Hordyk
#' @export 
plot.Fleet <- function(Fleet, Stock=Albacore, nsamp=3, nsim=500, nyears=50, proyears=28, 
  col="darkgray", breaks=10, lwd=2) { 
  cpars <- NULL
  if (class(Fleet) == "OM") {
    if (is.finite(Fleet@nyears)) nyears <- Fleet@nyears
	if (is.finite(Fleet@proyears)) proyears <- Fleet@proyears
	if (is.finite(Fleet@nsim)) nsim <- Fleet@nsim	
	if (length(Fleet@cpars) > 0) {
	  cpars <- Fleet@cpars
	  ncparsim <- cparscheck(cpars)
	}
	Stock <- SubOM(Fleet, "Stock")
	Fleet <- SubOM(Fleet, "Fleet")
  }
  if (class(Stock) != "Stock") stop("Must include a Stock object", call.=FALSE)
  
  its <- sample(1:nsim, nsamp)  
  Spat_targ <- runif(nsim, Fleet@Spat_targ[1], Fleet@Spat_targ[2]) 
  Esd <- runif(nsim, Fleet@Esd[1], Fleet@Esd[2])  # interannual variability in fishing effort (log normal sd)
  EffLower <- Fleet@EffLower
  EffUpper <- Fleet@EffUpper 
  EffYears <- Fleet@EffYears
  Deriv <- getEffhist(Esd, nyears, EffYears = Fleet@EffYears, EffLower = Fleet@EffLower, EffUpper = Fleet@EffUpper)  # Historical fishing effort
  Find <- Deriv[[1]]  # Calculate fishing effort rate 
  
  # Sample fishing efficiency parameters
  # =======================================================
  qinc <- runif(nsim, Fleet@qinc[1], Fleet@qinc[2])
  qcv <- runif(nsim, Fleet@qcv[1], Fleet@qcv[2])  # interannual variability in catchability  
  qmu <- -0.5 * qcv^2  # Mean
  qvar <- array(exp(rnorm(proyears * nsim, rep(qmu, proyears), rep(qcv, proyears))), c(nsim, proyears))  # Variations in interannual variation

  # Sample Biological Parameters for Size Comps
  maxage <- Stock@maxage  # maximum age (no plus group)
  calcMax <- -log(0.01)/(min(Stock@M))        # Age at which 1% of cohort survives
  maxage <- round(max(maxage, calcMax),0)  # If maximum age is lower, increase it to calcMax
     
  Linf <- runif(nsim, Stock@Linf[1], Stock@Linf[2])  # sample of asymptotic length
  Linfsd <- runif(nsim, Stock@Linfsd[1], Stock@Linfsd[2])  # sample of interannual variability in Linf

  K <- runif(nsim, Stock@K[1], Stock@K[2])  # now predicted by a log-linear model
  Ksd <- runif(nsim, Stock@Ksd[1], Stock@Ksd[2])  #runif(nsim,Stock@Ksd[1],Stock@Ksd[2])# sd is already added in the linear model prediction
  Kgrad <- runif(nsim, Stock@Kgrad[1], Stock@Kgrad[2])  # gradient in Von-B K parameter (K y-1)
  t0 <- runif(nsim, Stock@t0[1], Stock@t0[2])  # a sample of theoretical age at length zero  
  L50 <- array(runif(nsim * 50, Stock@L50[1], Stock@L50[2]), c(nsim, 50))  # length at 50% maturity
  L50[L50/Linf > 0.95] <- NA
  L50 <- apply(L50, 1, function(x) x[!is.na(x)][1])  
  

  Linfrand <- matrix(exp(rnorm(nsim*(proyears+nyears), -0.5 * Linfsd^2, Linfsd)), nrow=nsim, ncol=proyears+nyears)
  Krand <- matrix(exp(rnorm(nsim*(proyears+nyears), -0.5 * Ksd^2, Ksd)), nrow=nsim, ncol=proyears+nyears)  
  Linfgrad <- runif(nsim, Stock@Linfgrad[1], Stock@Linfgrad[2])  # sample of gradient in Linf (Linf y-1)

  Linfarray <- gettempvar(Linf, Linfsd, Linfgrad, nyears + proyears, nsim, Linfrand)  # Linf array
  Karray <- gettempvar(K, Ksd, Kgrad, nyears + proyears, nsim, Krand)  # the K array
  
  Agearray <- array(rep(1:maxage, each = nsim), dim = c(nsim, maxage))  # Age array  
  Len_age <- array(NA, dim = c(nsim, maxage, nyears + proyears))  # Length at age array
  ind <- as.matrix(expand.grid(1:nsim, 1:maxage, 1:(nyears + proyears)))  # an index for calculating Length at age
  Len_age[ind] <- Linfarray[ind[, c(1, 3)]] * (1 - exp(-Karray[ind[, c(1, 3)]] * 
    (Agearray[ind[, 1:2]] - t0[ind[, 1]])))

  LatASD <- Len_age * 0.1  # SD of length-at-age - this is currently fixed to cv of 10%
  MaxBin <- ceiling(max(Linfarray) + 3 * max(LatASD))
  binWidth <- ceiling(0.03 * MaxBin)
  CAL_bins <- seq(from = 0, to = MaxBin + binWidth, by = binWidth)
  CAL_binsmid <- seq(from = 0.5 * binWidth, by = binWidth, length = length(CAL_bins) - 1)
  nCALbins <- length(CAL_binsmid)
  
  # Sample selectivity parameters 
  # =======================================================  
  Selnyears <- length(Fleet@SelYears)
  # are selectivity parameters relative to size at maturity?
  chk <- class(Fleet@isRel)
  if (length(Fleet@isRel) < 1) 
    Fleet@isRel <- "true"
  if (chk == "character") {
    chkRel <- tolower(Fleet@isRel)
    if (chkRel == "true" | Fleet@isRel == "1") 
      multi <- L50
    if (chkRel == "false" | Fleet@isRel == "0") 
      multi <- 1
  }
  if (chk == "numeric") {
    if (Fleet@isRel == 1) 
      multi <- L50
    if (Fleet@isRel == 0) 
      multi <- 1
  }
  L5 <- runif(nsim, Fleet@L5[1], Fleet@L5[2]) * multi  # length at 0.05% selectivity ascending
  LFS <- runif(nsim, Fleet@LFS[1], Fleet@LFS[2]) * multi  # first length at 100% selection
  Vmaxlen <- runif(nsim, Fleet@Vmaxlen[1], Fleet@Vmaxlen[2])  # selectivity at maximum length
  L5s <- LFSs <- Vmaxlens <- NULL  # initialize 
  
  if (Selnyears > 1) {   # change of selectivity in historical years 
    # length at 0.05% selectivity ascending
	L5s <- mapply(runif, n = nsim, min = Fleet@L5Lower, max = Fleet@L5Upper) * multi
    # first length at 100% selection
    LFSs <- mapply(runif, n = nsim, min = Fleet@LFSLower, max = Fleet@LFSUpper) *  multi
	# selectivity at maximum length
	Vmaxlens <- mapply(runif, n = nsim, min = Fleet@VmaxLower, max = Fleet@VmaxUpper)
  }  
  	
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
                 "Mrand","Linfrand","Krand","maxage","V","Depletion", # end of Fleet variables
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
	  
	  if ("EffUpper" %in% Names & !"Find" %in% Names) {
      Deriv <- getEffhist(Esd, nyears, EffYears = EffYears, EffUpper = EffUpper, EffLower = EffLower)  # Historical fishing effort
      Find <- Deriv[[1]]  # Calculate fishing effort rate
      dFfinal <- Deriv[[2]]  # Final gradient in fishing effort yr-1 
    }
  }
  
  if(exists("V",inherits=FALSE))if(dim(V)[3] != proyears + nyears)   V<-abind(V,array(V[,,nyears],c(nsim,maxage,proyears)),along=3) # extend future Vulnerabiliy according to final historical vulnerability  
  if (any(dim(Find) != c(nsim, nyears))) stop("Find must be matrix with dimensions: nsim, nyears")  
  
  # Selectivity at Length
  # ------------------------------------------------------ if (max(Fleet@L5)
  # > 1.5) { message('L5 set too high (maximum value of 1.5).
  # \nDefaulting to L5 = 1.5') Fleet@L5[Fleet@L5 > 1.5] <- 1.5 }
  SLarray <- array(NA, dim=c(nsim, nCALbins, nyears+proyears)) # Selectivity-at-length 
  
  if (exists("V", inherits=FALSE) & length(cpars)>0) { # V has been passed in with custompars 
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
  
  if (!exists("V", inherits=FALSE) | length(cpars)==0) { # don't run if V has been passed in with custompars 
    if (Selnyears <= 1) {    
      if (class(L5) != "matrix") L5 <-  matrix(L5, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
      if (class(LFS) != "matrix") LFS <- matrix(LFS, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
      if (class(Vmaxlen) != "matrix") Vmaxlen <- matrix(Vmaxlen, nrow = nyears + proyears, ncol = nsim, byrow = TRUE) 
    
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
      SelYears <- Fleet@SelYears
	  
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


		  
  # Start plotting shenanigans 
  ncol <- 12
 
  m <- layout(matrix(c(c(rep(1,3), rep(2,3), rep(3,3), rep(4,3)),
					   c(rep(0, ncol)),
					   c(rep(5,6), rep(6,6)),
					   c(rep(0, ncol)),
					   c(rep(7,4), rep(8, 4), rep(9,4)),
					   c(rep(0, ncol)),
					   c(rep(10,4), rep(11,4), rep(12,4))
					   ), ncol=ncol, byrow = TRUE),
					   heights=c(1,0.2, 1,0.3, 1, 0.3, 1))
 # layout.show(m)				   
									   
  op <- par(mar = c(4, 3, 3, 1), oma=c(1,2,4,1), las=1) 
  hist(Esd, col=col, axes=FALSE, main="Esd", breaks=breaks)
  axis(side=1) 
  hist(qinc, col=col, axes=FALSE, main="qinc", breaks=breaks)
  axis(side=1)  
  hist(qcv, col=col, axes=FALSE, main="qcv", breaks=breaks)
  axis(side=1) 
  hist(Spat_targ, col=col, axes=FALSE, main="Spat_targ", breaks=breaks)
  axis(side=1) 
  title("not currently used", line=0)  

  
  # Effort 
  matplot(t(Find[its,]), type="l", lwd=lwd, bty="l", main="Fishing effort\ (historical)")
  
  # Future catchability
  ind <- as.matrix(expand.grid(its, 1:proyears, 1:nsamp))
  Qfuture <- matrix(NA, nrow=proyears, ncol=nsamp)
  X <- 0 
  for (sim in its) {
    X <- X + 1 
    Qfuture[,X] <- qvar[sim,1] * (1 + qinc[sim]/100)^(1:proyears)
	Qfuture[,X] <- Qfuture[,X]/Qfuture[1,X]
  }
  matplot(Qfuture, type="l", lwd=lwd, bty="l", main="Change in future\n catchability")  
  
  # Selectivity at length
  sampV <- V[its,,]
  matplot(t(Len_age[its,,1]), t(sampV[,,1]), type="l", lwd=lwd, bty="l", main="First historical\n year", xlab="Length")  
  matplot(t(Len_age[its,,nyears]), t(sampV[,,nyears]), type="l", lwd=lwd, bty="l", main="Last historical\n year", xlab="Length")  
  title(line=3, cex.main=1.5, "Selectivity-at-length", xpd=NA)
  matplot(t(Len_age[its,,nyears+proyears]), t(sampV[,,nyears+proyears]), type="l", lwd=lwd, bty="l", 
    main="Last projected\n year", xlab="Length")  
  
  
  matplot(t(sampV[,,1]), type="l", lwd=lwd, bty="l", main="First historical\n year", xlab="Age")  
  matplot(t(sampV[,,nyears]), type="l", lwd=lwd, bty="l", main="Last historical\n year", xlab="Age")  
  title(line=3, cex.main=1.5, "Selectivity-at-age", xpd=NA)
  matplot(t(sampV[,,nyears+proyears]), type="l", lwd=lwd, bty="l", main="Last projected\n year", xlab="Age")  
  
  title(Fleet@Name, outer=TRUE)
  title(paste("nyears =", nyears, "  proyears =", proyears, "  ", nsamp, "sampled iterations"), outer=TRUE, line=0)
    
  on.exit(par(op))	  
}
