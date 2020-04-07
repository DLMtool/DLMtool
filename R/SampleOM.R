
myrunif <- function(n, val1, val2) {
  min <- min(c(val1, val2))
  max <- max(c(val1, val2))
  
  if (is.na(n)) stop("First argument is NA")
  if (is.na(val1)) stop('Second argument is NA')
  if (is.na(val2)) stop('Third argument is NA')
  
  if (all(is.na(c(min, max)))) return(rep(NA,n))
  if (all(min == max)) {
    tt <- runif(n)
    return(rep(min, n))
  } else {
    return(runif(n, min, max))
  }
}


#' Sample Stock parameters
#'
#' @param Stock An object of class 'Stock' or class 'OM'
#' @param nsim Number of simulations. Ignored if 'Stock' is class 'OM'
#' @param nyears Number of historical years. Ignored if 'Stock' is class 'OM'
#' @param proyears Number of projection years. Ignored if 'Stock' is class 'OM'
#' @param cpars Optional named list of custom parameters. Ignored if 'Stock' is class 'OM'
#' @param msg logical. Warning message for M values?
#'
#' @return A named list of sampled Stock parameters
#' @keywords internal
#' @export
#'   
SampleStockPars <- function(Stock, nsim=48, nyears=80, proyears=50, cpars=NULL, msg=TRUE) {
  if (class(Stock) != "Stock" & class(Stock) != "OM") 
    stop("First argument must be class 'Stock' or 'OM'")
  Stock <- updateMSE(Stock) # update to add missing slots with default values
  if (all(is.na(Stock@LenCV))) Stock@LenCV <- c(0.1, 0.1)
  if (all(is.na(Stock@Mexp))) Stock@Mexp <- c(0, 0)
  
  # Warning alerts for deprecated slots 
  if (msg) {
    slots <- c("Linfgrad", "Kgrad", 'Mgrad')
    for (sl in slots) {
      val <- slot(Stock, sl)
      if (!all(is.na(val)) & !all(val ==0))
        warning(sl, " is no longer used and values are being ignored. Use 'cpars' to specify time-varying changes to ", sl, call.=FALSE)
    }
  }
 
  if (class(Stock) == "OM") {
    nsim <- Stock@nsim
    nyears <- Stock@nyears 
    proyears <- Stock@proyears
  }
  
  # Get custom pars if they exist
  if (class(Stock) == "OM" && length(Stock@cpars) > 0 && is.null(cpars)) cpars <- SampleCpars(Stock@cpars, nsim)  # custom parameters exist in Stock/OM object
  if (length(cpars) > 0) { # custom pars exist - assign to function environment 
    for (X in 1:length(cpars)) assign(names(cpars)[X], cpars[[X]])
  }
  
  StockOut <- list() 
  
  # == Maximum age ====
  if (!exists("maxage", inherits=FALSE)) {
    StockOut$maxage <- maxage <- Stock@maxage # maximum age (no plus group)
  } else StockOut$maxage <- maxage
  
  
  # == Virgin Recruitment ====
  if (!exists("R0", inherits=FALSE)) R0 <- Stock@R0  # Initial recruitment
  if (length(R0) != nsim) R0 <- rep(R0, nsim)[1:nsim] # modified to allow for different R0 per sim 
  StockOut$R0 <- R0
  
  # == Natural Mortality ====
  n_age <- maxage + 1 # number of age classes (including age-0)
  # natural mortality rate
  if (length(Stock@M) == 2 & !exists("M", inherits=FALSE)) M <- myrunif(nsim, Stock@M[1], Stock@M[2])  
  
  if (length(Stock@M) == n_age) { # Stock@M is vector of M-at-age 
    if (length(Stock@M2) == n_age && !exists("Mage", inherits=FALSE)) {
      mmat <- rbind(Stock@M, Stock@M2)
      if (all(mmat[1,] == mmat[2,])) {
        Mage <- matrix(mmat[1,], nsim, n_age, byrow=TRUE)
      } else {
        if (all(mmat[1,] < mmat[2,]) | all(mmat[1,] > mmat[2,])) {
          Mage <- matrix(NA, nsim, n_age)
          Mage[,1] <- myrunif(nsim, min(mmat[,1]), max(mmat[,1]))
          val <- (Mage[,1] - min(mmat[,1]))/ diff(mmat[,1])
          for (X in 2:n_age) Mage[,X] <- min(mmat[,X]) + diff(mmat[,X])*val  
        } else stop("All values in slot 'M' must be greater or less than corresponding values in slot 'M2'", call.=FALSE)
      }
      
    } else stop("slot 'M2' must be length 'maxage+1'", call.=FALSE)
  } 
  if (length(Stock@M) != n_age & length(Stock@M) != 2) stop("slot 'M' must be either length 2 or length maxage+1", call.=FALSE)
  
  if (length(Stock@M2) == n_age & !length(Stock@M) == n_age) {
    stop("Slot M2 is used (upper bound on M-at-age) and is length 'maxage+1' but Slot M (lower bound on M-at-age) is not length 'maxage+1'.")
  }
  if (!exists("Msd", inherits=FALSE)) Msd <- myrunif(nsim, Stock@Msd[1], Stock@Msd[2])  # sample inter annual variability in M frStock specified range
  # if (!exists("Mgrad", inherits=FALSE)) Mgrad <- myrunif(nsim, Stock@Mgrad[1], Stock@Mgrad[2])  # sample gradient in M (M y-1)
  if (.hasSlot(Stock, "Mexp") & !exists("Mexp", inherits=FALSE)) {
    if (all(is.numeric(Stock@Mexp) & is.finite(Stock@Mexp))) {
      Mexp <- myrunif(nsim, min(Stock@Mexp), max(Stock@Mexp)) # sample Lorenzen M-at-weight exponent     
    } else {
      Mexp <- rep(0, nsim) # assume constant M-at-age/size
    }
  } 

  if (!exists("M", inherits=FALSE)) M <- Mage[,n_age]
  if (!exists("Mexp", inherits=FALSE)) Mexp <- rep(0, nsim) # assume constant M-at-age/size if it is not specified 
  if (!all(Mexp == 0) & length(Stock@M2) == n_age) {
    stop("Values in both M2 and Mexp slots. Only one can be used")
  }

  # == Depletion ====
  if (!exists("D", inherits=FALSE)) {
    StockOut$D <- D <- myrunif(nsim, Stock@D[1], Stock@D[2])  # sample from the range of user-specified depletion (Bcurrent/B0)  
  } else {
    StockOut$D <- D 
  }
  
  # == Stock-Recruitment Relationship ====
  if (!exists("SRrel", inherits=FALSE)) {
    StockOut$SRrel <- rep(Stock@SRrel, nsim)  # type of Stock-recruit relationship. 1=Beverton Holt, 2=Ricker
  } else {
    StockOut$SRrel <- SRrel 
  }
  
  if (exists("h", inherits = FALSE)) hs <- h
  if (!exists("hs", inherits=FALSE)) {
    StockOut$hs <- hs <- myrunif(nsim, Stock@h[1], Stock@h[2])  # sample of recruitment compensation (steepness - fraction of unfished recruitment at 20% of unfished biStockass)
  } else {
    StockOut$hs <- hs
  }
  if (any(StockOut$hs > 1 | StockOut$hs < 0.2)) stop("Steepness (OM@h) must be between 0.2 and 1", call.=FALSE)
 
  # == Recruitment Deviations ====
  if (exists("Perr", inherits = FALSE)) {
    procsd <- Perr
  }
  
  
  if (!exists("Perr_y", inherits=FALSE)) {
    if (!exists("procsd", inherits=FALSE)) {
      StockOut$procsd <- procsd <- myrunif(nsim, Stock@Perr[1], Stock@Perr[2])  # Process error standard deviation
    } else {
      StockOut$procsd <- procsd
    }
    
    if (!exists("AC", inherits=FALSE)) {
      StockOut$AC <- AC <- myrunif(nsim, Stock@AC[1], Stock@AC[2]) 
      # auto correlation parameter for recruitment deviations recdev(t)<-AC*recdev(t-1)+(1-AC)*recdev_proposed(t)  
    } else {
      StockOut$AC <- AC 
      # auto correlation parameter for recruitment deviations recdev(t)<-AC*recdev(t-1)+(1-AC)*recdev_proposed(t)
    }
    
    # All recruitment Deviations
    # Add cycle (phase shift) to recruitment deviations - if specified
    if (is.finite(Stock@Period[1]) & is.finite(Stock@Amplitude[1])) {
      # Shape <- "sin"  # default sine wave - alternative - 'shift' for step changes
      Period <- myrunif(nsim, min(Stock@Period), max(Stock@Period))
      if (max(Stock@Amplitude)>1) {
        if (msg) message("Stock@Amplitude > 1. Defaulting to 1")
        Stock@Amplitude[Stock@Amplitude>1] <- 1
      }
      Amplitude <- myrunif(nsim, min(Stock@Amplitude), max(Stock@Amplitude))
      
      yrs <- 1:(nyears + proyears+n_age-1)
      recMulti <- t(sapply(1:nsim, function(x) 1+sin((runif(1, 0, 1)*max(yrs) + 2*yrs*pi)/Period[x])*Amplitude[x]))
      if (msg) message("Adding cyclic recruitment pattern")
    } else {
      recMulti <- 1 
    }
    StockOut$procmu <- procmu <- -0.5 * procsd^2  * (1 - AC)/sqrt(1 - AC^2) #  # adjusted log normal mean http://dx.doi.org/10.1139/cjfas-2016-0167
    
    Perr_y <- array(rnorm((nyears + proyears+n_age-1) * nsim, rep(procmu, nyears + proyears+n_age-1), 
                          rep(procsd, nyears + proyears+n_age-1)), c(nsim, nyears + proyears+n_age-1))
    for (y in 2:(nyears + proyears+n_age-1)) Perr_y[, y] <- AC * Perr_y[, y - 1] + Perr_y[, y] * (1 - AC * AC)^0.5  
    #2#AC*Perr[,y-1]+(1-AC)*Perr[,y] # apply a pseudo AR1 autocorrelation to rec devs (log space)
    
    StockOut$Perr_y <- Perr_y <- exp(Perr_y) * recMulti # normal space (mean 1 on average) 
    
  } else {
    StockOut$Perr_y <- Perr_y
    StockOut$procsd <- apply(Perr_y, 1, sd)
  }
 

  # if (nsim > 1) {
  #   cumlRecDev <- apply(Perr[, 1:(nyears+maxage-1)], 1, prod)
  #   dep[order(cumlRecDev)] <- dep[order(dep, decreasing = F)]  # robustifies 
  # }
  
  # == Growth parameters ====
  vars <- c("Linf", "Linfsd", "K", "Ksd", "t0")
  for (var in vars) {
    if (!exists(var, inherits=FALSE)) {
      if (all(is.na(slot(Stock, var)))) {
        val <- rep(0, nsim)
      } else {
        val <- myrunif(nsim, slot(Stock, var)[1],slot(Stock, var)[2])  
      }
      assign(var, val)
    } 

  }
    
  # == Sample Fecundity-Length Exponent ===
  # if (!exists("FecB", inherits=FALSE))   FecB <- runif(nsim, min(Stock@FecB), max(Stock@FecB))
  
  # == Sample Spatial Parameters ====
  if (!exists("Frac_area_1", inherits=FALSE)) Frac_area_1 <- myrunif(nsim, Stock@Frac_area_1[1], Stock@Frac_area_1[2])  # sampled fraction of unfished biStockass in area 1 (its a two area model by default)
  if (!exists("Prob_staying", inherits=FALSE)) Prob_staying <- myrunif(nsim, Stock@Prob_staying[1], Stock@Prob_staying[2])  # sampled probability of individuals staying in area 1 among years
  if (!exists("Size_area_1", inherits=FALSE)) Size_area_1 <- myrunif(nsim, Stock@Size_area_1[1], Stock@Size_area_1[2])  # currently redundant parameter for the habitat area size of area 1
  
  if (max(Size_area_1) == 0) stop("Size_area_1 must be > 0", call. = FALSE)
  if (max(Frac_area_1) == 0) stop("Frac_area_1 must be > 0", call. = FALSE)
  if (max(Prob_staying) == 0) stop("Prob_staying must be > 0", call. = FALSE)
  
  if (max(Size_area_1) >= 1) stop("Size_area_1 must be < 1", call. = FALSE)
  if (max(Frac_area_1) >= 1) stop("Frac_area_1 must be < 1", call. = FALSE)
  if (max(Prob_staying) >= 1) stop("Prob_staying must be < 1", call. = FALSE)
  
  StockOut$Frac_area_1 <- Frac_area_1
  StockOut$Prob_staying <- Prob_staying
  StockOut$Size_area_1 <- Size_area_1
  
  if (!exists('Asize', inherits=FALSE)) Asize <- cbind(StockOut$Size_area_1, 1 - StockOut$Size_area_1)
  
  # === Generate random numbers for random walk ====
  if (!exists("Mrand", inherits=FALSE)) Mrand <- matrix(exp(rnorm(nsim*(proyears+nyears), -0.5 * Msd^2, Msd)), nrow=nsim, ncol=proyears+nyears)
  if (!exists("Linfrand", inherits=FALSE)) Linfrand <- matrix(exp(rnorm(nsim*(proyears+nyears), -0.5 * Linfsd^2, Linfsd)), nrow=nsim, ncol=proyears+nyears)
  if (!exists("Krand", inherits=FALSE)) Krand <- matrix(exp(rnorm(nsim*(proyears+nyears), -0.5 * Ksd^2, Ksd)), nrow=nsim, ncol=proyears+nyears)
  
  StockOut$Mrand <- Mrand
  StockOut$Linfrand <- Linfrand
  StockOut$Krand <- Krand
  
  # === Generate time-varying Linf, K and t0 arrays ====
  # if (!exists("Linfarray", inherits=FALSE)) Linfarray <- gettempvar(Linf, Linfsd, Linfgrad, nyears + proyears, nsim, Linfrand)  # Linf array  
  # if (!exists("Karray", inherits=FALSE)) Karray <- gettempvar(K, Ksd, Kgrad, nyears + proyears, nsim, Krand)  # the K array
  
  if (!exists("Linfarray", inherits=FALSE)) Linfarray <- gettempvar(Linf, Linfsd, targgrad=0, nyears + proyears, nsim, Linfrand)  # Linf array  
  if (!exists("Karray", inherits=FALSE)) Karray <- gettempvar(K, Ksd, targgrad=0, nyears + proyears, nsim, Krand)  # the K array
  if (!exists("Agearray", inherits=FALSE))  Agearray <- array(rep(0:maxage, each = nsim), dim = c(nsim, n_age))  # Age array
  
  if (all(dim(Linfarray) != c(nsim, nyears+proyears))) stop("Linfarray must be dimensions: nsim, proyears+nyears (", nsim, ", ", proyears+nyears, ")")
  if (all(dim(Karray) != c(nsim, nyears+proyears))) stop("Karray must be dimensions: nsim, proyears+nyears (", nsim, ", ", proyears+nyears, ")")
  
  if (length(StockOut$maxage) > 1) StockOut$maxage <- StockOut$maxage[1] # check if maxage has been passed in custompars
  
  t0array <- matrix(t0, nrow=nsim, ncol=proyears+nyears)
  
  # == Sample CV Length-at-age ====
  if (!exists("LenCV", inherits=FALSE)) LenCV <- myrunif(nsim, min(Stock@LenCV), max(Stock@LenCV))
  
  if (msg && any(LenCV < 0.05)) 
    warning('Stock@LenCV is very low for at least some simulations (<0.05).\nLength composition data may not be generated successfully and MPs using length data may crash or be unreliable. \nLenCV is the variation in length-at-age. Very low values implies all individuals exactly follow the average growth curve')
  
  # === Create Mean Length-at-Age array ====
  if (!exists("Len_age", inherits=FALSE)) {
    Len_age <- array(NA, dim = c(nsim, n_age, nyears + proyears))  # Length at age array
    ind <- as.matrix(expand.grid(1:nsim, 1:n_age, 1:(nyears + proyears)))  # an index for calculating Length at age
    Len_age[ind] <- Linfarray[ind[, c(1, 3)]] * (1 - exp(-Karray[ind[, c(1, 3)]] * 
                                                           (Agearray[ind[, 1:2]] - t0[ind[, 1]])))
    
    if (class(Stock)=="OM" && length(Stock@cpars[['Linf']]) >0) {
      maxLinf <- max(Stock@cpars$Linf)
    } else {
      maxLinf <- max(Stock@Linf)
    }
    # linfs <- gettempvar(maxLinf, 0, max(Stock@Linfgrad), nyears + proyears, 
    #            1, matrix(1, nrow=1, ncol=proyears+nyears))
    # MaxBin <- ceiling(max(linfs) + 3 * max(linfs) * max(Stock@LenCV)) 
    
    MaxBin <- ceiling(max(Linfarray) + 2 * max(Linfarray) * max(LenCV))

  } else { # Len_age has been passed in with cpars
    if (any(dim(Len_age) != c(nsim, n_age, nyears + proyears))) 
      stop("'Len_age' must be array with dimensions: nsim, maxage+1, nyears + proyears") 
    # Estimate vB parameters for each year and each sim 
    if (!all(c("Linf", "K", "t0") %in% names(cpars))) { # don't calculate if Linf, K and t0 have also been passed in with cpars
      vB <- function(pars, ages) pars[1] * (1-exp(-pars[2]*(ages-pars[3])))
      fitVB <- function(pars, LatAge, ages) sum((vB(pars, ages) - LatAge)^2)
      starts <- c(max(Len_age), 0.2, 0)
      if(msg) message("Estimating growth parameters from length-at-age array in cpars")
      for (ss in 1:nsim) {
        if(msg) {
          cat(".")
          flush.console()
        }
        pars <- sapply(1:(nyears + proyears), function(X) optim(starts, fitVB, LatAge=Len_age[ss,,X], ages=0:maxage)$par)
        Linfarray[ss,] <- round(pars[1,],2)
        Karray[ss,] <- round(pars[2,],2)
        t0[ss]<- mean(pars[3,])
      }
      Linf <- Linfarray[, nyears]
      K <- Karray[, nyears]
      t0array <- matrix(t0, nrow=nsim, ncol=proyears+nyears)
      if (msg) cat("\n")
    }
    # MaxBin <- ceiling(max(Len_age) + 3 * max(Len_age) * max(Stock@LenCV)) 
    MaxBin <- ceiling(max(Linfarray) + 2 * max(Linfarray) * max(LenCV))
  }
  Len_age[Len_age<0] <- 0.001
  StockOut$maxlen <- maxlen <- Len_age[, n_age, nyears] # reference length for Vmaxlen 
  
  # == Generate Catch at Length Classes ====
  if (!exists("LatASD", inherits=FALSE)) LatASD <- Len_age * array(LenCV, dim=dim(Len_age)) # SD of length-at-age 
  if (any(dim(LatASD) != dim(Len_age))) stop("Dimensions of 'LatASD' must match dimensions of 'Len_age'", .call=FALSE)
  
  if (exists("CAL_bins", inherits=FALSE)) binWidth <- CAL_bins[2] - CAL_bins[1]
  if (exists("CAL_binsmid", inherits=FALSE)) binWidth <- CAL_binsmid[2] - CAL_binsmid[1]
    
  if (!exists("binWidth", inherits=FALSE)) binWidth <- ceiling(0.03 * MaxBin)
  
  if (!exists("CAL_bins", inherits=FALSE)) CAL_bins <- seq(from = 0, to = MaxBin + binWidth, by = binWidth)
  if (!exists("CAL_binsmid", inherits=FALSE)) CAL_binsmid <- seq(from = 0.5 * binWidth, by = binWidth, length = length(CAL_bins) - 1)
  if (length(CAL_bins) != length(CAL_binsmid)+1) stop("Length of 'CAL_bins' must be length(CAL_binsmid)+1", .call=FALSE)
  
  if (is.null(cpars$binWidth))
    binWidth <- CAL_binsmid[2] - CAL_binsmid[1]
  
  # Check bin width - in case both CAL_bins or CAL_binsmid AND binWidth have been passed in with cpars
  if (!all(diff(CAL_bins) == binWidth)) stop("width of CAL_bins != binWidth", call.=FALSE)
  if (!all(diff(CAL_binsmid) == binWidth)) stop("width of CAL_binsmid != binWidth", call.=FALSE)
  nCALbins <- length(CAL_binsmid)
  
  if (max(Linfarray) > max(CAL_bins)) stop("`max(CAL_bins)` must be larger than `max(Linfarray)`")
 
  # === Create Weight-at-Age array ====
  if (!exists("Wt_age", inherits=FALSE)) {
    Wt_age <- array(NA, dim = c(nsim, n_age, nyears + proyears))  # Weight at age array
    ind <- as.matrix(expand.grid(1:nsim, 1:n_age, 1:(nyears + proyears)))  # an index for calculating Weight at age 
    Wt_age[ind] <- Stock@a * Len_age[ind]^Stock@b  # Calculation of weight array
    Wa <- Stock@a
    Wb <- Stock@b 
  }	else {
    if (any(dim(Wt_age) != c(nsim, n_age, nyears + proyears))) 
      stop("'Wt_age' must be array with dimensions: nsim, maxage+1, nyears + proyears (", paste(c(nsim, maxage+1, nyears + proyears), ""), ") but has ", paste(dim(Wt_age), "")) 
    # Estimate length-weight parameters from the Wt_age data
    logL <- log(as.numeric(Len_age)+tiny)
    logW <- log(as.numeric(Wt_age)+tiny)
    mod  <- lm(logW ~ logL)
    EstVar <- summary(mod)$sigma^2
    Wa <- as.numeric(exp(coef(mod)[1]) * exp((EstVar)/2))
    Wb <- as.numeric(coef(mod)[2])
  }
  
  # == Sample Maturity Parameters ====
  if (exists("Mat_age", inherits=FALSE)){
    if (any(dim(Mat_age) != c(nsim, n_age, nyears+proyears))) stop("'Mat_age' must be array with dimensions: nsim, maxage+1, nyears+proyears") 
    
    # Calculate L50, L95, ageM and age95 
    ageM <- age95 <- L50array <- L95array <- matrix(NA, nsim, nyears+proyears)
    for (XX in 1:(nyears+proyears)) {
      # check that Mat_age < 0.5 values exist
     if (nsim == 1) {
       oksims <- which(min(Mat_age[1,,XX]) < 0.5)
     } else {
       oksims <- which(apply(Mat_age[,,XX], 1, min) < 0.5) 
     }
      if (length(oksims)<1) {
        ageM[,XX] <- 1 # set to 1 if < 1
        L50array[,XX] <- 1 # set to 1 if < 1
      } else {
        noksims <- (1:nsim)[-oksims]
        ageM[oksims,XX] <- unlist(sapply(oksims, function(x) LinInterp(Mat_age[x,, XX], y=1:n_age, 0.5)))
        ageM[noksims,XX] <- 1 # set to 1 
        L50array[oksims,XX] <- unlist(sapply(oksims, function(x) LinInterp(Mat_age[x,,XX], y=Len_age[x, , nyears], 0.5)))
        L50array[noksims,XX] <- 1 # set to 1 
      }
    
      age95[,XX] <- unlist(sapply(1:nsim, function(x) LinInterp(Mat_age[x,, XX], y=1:n_age, 0.95)))
      L95array[,XX]<- unlist(sapply(1:nsim, function(x) LinInterp(Mat_age[x,,XX], y=Len_age[x, , nyears], 0.95)))
    }
    
    L50array[!is.finite(L50array)] <- 0.8*Linfarray[!is.finite(L50array)]
    L95array[!is.finite(L95array)] <- 0.99*Linfarray[!is.finite(L95array)]
    
    L95array[L50array >= L95array] <- L50array[L50array >= L95array] * 1.01
    L50 <- L50array[,nyears]
    L95 <- L95array[,nyears]
    L50[!is.finite(L50)] <- 0.8*Linf[!is.finite(L50)]
    L95[!is.finite(L95)] <- 0.99*Linf[!is.finite(L95)]
    if (any(L50>= Linf)) {
      if (msg) message("Note: Some samples of L50 are above Linf. Defaulting to 0.8*Linf")
      L50[L50>=Linf] <- 0.8* Linf[L50>=Linf]
    }
    if (any(L95> Linf)) {
      if (msg)  message("Note: Some samples of L95 are above Linf. Defaulting to 0.99*Linf")
      L95[L95> Linf] <- 0.99* Linf[L95> Linf]
    }
    
    L50_95 <- L95 - L50
  } else {
    if (!exists("L50", inherits=FALSE)) {
      sL50 <- array(myrunif(nsim * 50, Stock@L50[1], Stock@L50[2]), c(nsim, 50))  # length at 50% maturity  
      # checks for unrealistically high length at maturity
      sL50[sL50/Linf > 0.95] <- NA
      L50 <- apply(sL50, 1, function(x) x[!is.na(x)][1])
      L50[is.na(L50)] <- 0.95 * Linf[is.na(L50)]
    }
    if (!exists("L50_95", inherits=FALSE)) {
      L50_95 <- array(myrunif(nsim * 50, Stock@L50_95[1], Stock@L50_95[2]), c(nsim, 50))  # length at 95% maturity
      if (!exists("sL50", inherits=FALSE)) sL50 <- matrix(L50, nsim, 50)
      L50_95[((sL50+L50_95)/matrix(Linf, nsim, 50)) > 0.99] <- NA
      L50_95 <- apply(L50_95, 1, function(x) x[!is.na(x)][1]) 
      L50_95[is.na(L50_95)] <- 2
    }
    
    if (!exists("L95", inherits=FALSE))   L95 <- L50 + L50_95
    
    if (any(L95> Linf)) {
      message("Note: Some samples of L95 are above Linf. Defaulting to 0.99*Linf")
      L95[L95> Linf] <- 0.99* Linf[L95> Linf]
      L50_95 <- L95 - L50 
    }
    
    # === Generate L50 by year ====
    relL50 <- matrix(L50/Linf, nrow=nsim, ncol=nyears + proyears, byrow=FALSE) # assume L50/Linf stays constant 
    L50array <- relL50 * Linfarray
    delLm <- L95 - L50 
    L95array <- L50array + matrix(delLm, nrow=nsim, ncol=nyears + proyears, byrow=FALSE)
    L95array[L95array>Linfarray] <- 0.99 *  Linfarray[L95array>Linfarray]
  }
  

  # == Calculate age at maturity ==== 
  if (exists('ageM', inherits=FALSE)) { # check dimensions 
    if (!all(dim(ageM) == c(nsim, proyears+nyears))) stop('"ageM" must be dimensions: nsim, nyears+proyers')
  }
  if (exists('age95', inherits=FALSE)) { # check dimensions 
    if (!all(dim(age95) == c(nsim, proyears+nyears))) stop('"age95" must be dimensions: nsim, nyears+proyers')
  }
  
  if (!exists("ageM", inherits=FALSE)) ageM <- -((log(1 - L50array/Linfarray))/Karray) + t0array # calculate ageM from L50 and growth parameters (time-varying)
  ageM[ageM < 1] <- 1  # age at maturity must be at least 1
  if (!exists("age95", inherits=FALSE)) age95 <- -((log(1 - L95array/Linfarray))/Karray) + t0array
  age95[age95 < 1] <- 1.5  # must be greater than 0 and ageM
  
  if (any(ageM >= maxage-1)) {
    if (msg) message("Note: Some samples of age of maturity are above 'maxage'-1. Defaulting to maxage-1")
    ageM[ageM >= (maxage-1)] <- maxage - 1 
  }
  if (any(age95 >= maxage)) {
    if (msg) message("Note: Some samples of age of 95 per cent maturity are above 'maxage'. Defaulting to maxage")
    age95[age95 >= maxage] <- maxage  
  }
  
  # == Generate Maturity-at-Age array ====
  if (!exists("Mat_age", inherits=FALSE)) {
    Mat_age <- array(NA, dim=c(nsim, n_age, nyears+proyears))
    for (XX in 1:(nyears+proyears)) {
      Mat_age[,,XX] <- 1/(1 + exp(-log(19) * ((Agearray - ageM[,XX])/(age95[,XX] - ageM[,XX])))) # Maturity at age array by year
    }
  } 
 
  # == Calculate M-at-Age from M-at-Length if provided ====
  if (exists("M_at_Length", inherits=FALSE)) {  # M-at-length data.frame has been provided in cpars
    
    MatLen <- matrix(NA, nsim, nrow(M_at_Length))
    MatLen[,1] <- runif(nsim, min(M_at_Length[1,2:3]), max(M_at_Length[1,2:3]))
    
    for (k in 1:nsim) {
      for (X in 2:nrow(M_at_Length)) {
        val <- (MatLen[k,1] - min(M_at_Length[1,2:3]))/ diff(t(M_at_Length[1,2:3]))
        MatLen[k,X] <- min(M_at_Length[X,2:3]) + diff(t(M_at_Length[X,2:3]))*val 
      }
    }
    
    # Calculate M at age
    Mage <- matrix(NA, nsim, n_age)
    for (sim in 1:nsim) {
      ind <- findInterval(Len_age[sim,,nyears], M_at_Length[,1])  
      Mage[sim, ] <- MatLen[sim, ind]  
    }
  }
  
  
  # == M-at-age has been provided in OM ====
  if (exists("Mage", inherits=FALSE)) {
    if (exists("M", inherits=FALSE) & length(cpars[["M"]])>0) 
      if (msg) message("M-at-age has been provided in OM. Overiding M from OM@cpars")
    
    temp <- gettempvar(1, Msd, targgrad=0, nyears + proyears, nsim, Mrand) # add Msd
    temp2 <- replicate(maxage, temp)
    temp2 <- aperm(temp2, c(1,3,2))
    M_ageArray <-  array(Mage, dim=c(nsim, maxage, proyears+nyears))
    M_ageArray <- temp2 * M_ageArray
    # M is calculated as mean M of mature ages
    M <- rep(NA, nsim)
    for (sim in 1:nsim) M[sim] <- mean(Mage[sim,(round(ageM[sim],0)+1):n_age])
  }
  
  # == Mean Natural mortality by simulation and year ====
  if (exists("M_ageArray", inherits=FALSE)) {
    if (!all(dim(M_ageArray) == c(nsim, n_age, proyears+nyears))) stop("'M_ageArray' must be array with dimensions: nsim, maxage+1, nyears + proyears") 
    if(msg) message("M_ageArray has been provided in OM@cpars. Ignoring OM@Mexp, OM@Msd, and OM@Mgrad")
    Mexp <- Msd <- Mgrad <- rep(0, nsim)
  }
   
  
  if (!exists("Marray", inherits=FALSE) & exists("M_ageArray", inherits=FALSE)) {
    Marray <- matrix(NA, nsim, nyears+proyears)
    for (yr in 1:(nyears+proyears)) {
      for (sim in 1:nsim) {
        Marray[sim, yr] <- mean(M_ageArray[sim, (ageM[sim,yr]+1):n_age,yr])
      }
    }
  }
  
  if (!exists("Marray", inherits=FALSE)) {
    Marray <- gettempvar(M, Msd, targgrad=0, nyears + proyears, nsim, Mrand)  # M by sim and year according to gradient and inter annual variability
  } else {
    if (any(dim(Marray) != c(nsim, nyears + proyears))) stop("'Marray' must be array with dimensions: nsim, nyears + proyears") 
  }
  
  # == Natural mortality by simulation, age and year ====
  if (!exists("M_ageArray", inherits=FALSE)) { # only calculate M_ageArray if it hasn't been specified in cpars
    M_ageArray <- array(NA, dim=c(nsim, n_age, nyears + proyears))
    if (exists("Mage", inherits=FALSE)) { # M-at-age has been provided
      temp1 <- Mage/ matrix(apply(Mage, 1, mean), nsim, n_age, byrow=FALSE)
      ind <- as.matrix(expand.grid(1:nsim, 1:n_age, 1:(nyears+proyears)))
      M_ageArray[ind] <- temp1[ind[,1:2]] * Marray[ind[,c(1,3)]]
    } else { # M-at-age calculated from Lorenzen curve 
      Winf <- Stock@a * Linf^Stock@b
      ind <- as.matrix(expand.grid(1:nsim, 1:n_age, 1:(nyears+proyears)))
      M_ageArray[ind] <- Marray[ind[,c(1,3)]] * (Wt_age[ind]/Winf[ind[,1]]) ^ Mexp[ind[,1]]  
    } 
    
    
    # == Scale M at age so that mean M of mature ages is equal to sampled M ====
    tempM_ageArray <- M_ageArray
    for (sim in 1:nsim) {
      matyrs <- (ageM[sim, nyears]+1):n_age
      if (length(matyrs) >1) {
        # scale <- Marray[sim,]/ apply(tempM_ageArray[sim,ageM[sim]:maxage,], 2, mean) 
        scale <- Marray[sim,]/ (apply(tempM_ageArray[sim,matyrs,], 2, sum)/length(matyrs)) # this is about 4 times faster
      } else if (length(matyrs)==1){
        scale <- Marray[sim,]/ tempM_ageArray[sim,(ageM[sim]+1):n_age,]  
      } 
      
      M_ageArray[sim,,] <- M_ageArray[sim,,] * matrix(scale, n_age, nyears+proyears, byrow=TRUE)
    }
    
  }
  
  # == Sample Discard Mortality ====
  if(!exists("Fdisc", inherits = FALSE)) Fdisc <- myrunif(nsim, min(Stock@Fdisc), max(Stock@Fdisc))
  StockOut$Fdisc <- Fdisc 
  
  # == 
   
  

  # Check if M-at-age is constant that Maxage makes sense
  if (all(M_ageArray[1,,1] == mean(M_ageArray[1,,1])) & all(M !=0)) { # constant M at age
    calcMax <- ceiling(-log(0.01)/(min(M)))        # Age at which 1% of cohort survives
    if (maxage < 0.95*calcMax && msg) {
      message("Note: Maximum age (", maxage, ") is lower than assuming 1% of cohort survives to maximum age (", calcMax, ")")
    }  
  }
  
  # --- Calculate movement ----
  initdist <- Pinitdist <- NULL
  if (!exists("mov", inherits=FALSE)) {
    if(msg) message("Optimizing for user-specified movement")  # Print a progress update
    nareas<-2 # default is a 2 area model
    mov1 <- array(t(sapply(1:nsim, getmov2, Frac_area_1 = Frac_area_1, 
                           Prob_staying = Prob_staying)), dim = c(nsim, nareas, nareas))
    mov<-array(NA,c(nsim,n_age,nareas,nareas))
    mind<-as.matrix(expand.grid(1:nsim,1:n_age,1:nareas,1:nareas))
    mov[mind]<-mov1[mind[,c(1,3,4)]]
    
    initdist <- array(0,c(nsim,n_age,nareas))
    initdist[,,1]<-Frac_area_1
    initdist[,,2]<- 1- Frac_area_1  
    
  }else{ # if mov is specified need to calculate age-based spatial distribution (Pinitdist to initdist)
    nareas<-dim(mov)[3]
    if(msg) message("Custom movement matrix detected: simulating movement among ",nareas," areas")
    if(is.na(dim(mov)[5])) {
      mind<-as.matrix(expand.grid(1:nsim,n_age,1:nareas,1:nareas))
    } else {
      mind<-as.matrix(expand.grid(1:nsim,n_age,1:nareas,1:nareas, 1)) # movement for 1st year
    }
    movedarray<-array(0,c(nsim,nareas,nareas))
    Pinitdist<-array(1/nareas,c(nsim,nareas))
    for(i in 1:20){ # convergence in initial distribution is assumed to occur in 20 iterations (generally overkill)
      movedarray[mind[,c(1,3,4)]]<-Pinitdist[mind[,c(1,3)]]*mov[mind] # distribution in from areas mulitplied by movement array
      Pinitdist<-apply(movedarray,c(1,3),sum) # add over to areas
    }
  }
  if(is.na(dim(mov)[5])) { # movement matrix only specified for one year
    mov <- array(mov, dim=c(dim(mov), nyears+proyears))
  }
  # check dimensions 
  if (any(dim(mov) != c(nsim,n_age,nareas,nareas, nyears+proyears)))
      stop('cpars$mov must be array with dimensions: \nc(nsim, maxage+1, nareas, nareas) \nOR \nc(nsim, maxage+1, nareas, nareas, nyears+proyears)', call.=FALSE)

  if (dim(Asize)[2]!=nareas) {
    if(msg) message('Asize is not length "nareas", assuming all areas equal size')
    Asize <- matrix(1/nareas, nrow=nsim, ncol=nareas)
  }
  
  StockOut$Mexp <- Mexp 
  StockOut$Msd <- Msd 
  # StockOut$Mgrad <- Mgrad
  
  StockOut$ageM <- ageM
  StockOut$age95 <- age95
  StockOut$Linfarray <- Linfarray
  StockOut$Karray <- Karray
  StockOut$Agearray <- Agearray
  StockOut$Marray <- Marray
  StockOut$M_ageArray <- M_ageArray
  StockOut$t0array <- t0array
  StockOut$Len_age <- Len_age
  StockOut$Linf <- Linf 
  StockOut$Linfsd <- Linfsd
  # StockOut$Linfgrad <- Linfgrad
  # StockOut$recgrad <- recgrad
  StockOut$K <- K
  StockOut$Ksd <- Ksd
  # StockOut$Kgrad <- Kgrad
  StockOut$t0 <- t0 
  StockOut$a <- Wa 
  StockOut$b <- Wb 
  StockOut$Wt_age <- Wt_age
  StockOut$L50 <- L50
  StockOut$L50array <- L50array
  StockOut$L95 <- L95
  StockOut$L50_95 <- L50_95
  StockOut$L95array <- L95array
  # StockOut$FecB <- FecB
  StockOut$Mat_age <- Mat_age
  
  StockOut$M <- M
  StockOut$LenCV <- LenCV
  StockOut$LatASD <- LatASD
  StockOut$CAL_binsmid <- CAL_binsmid
  StockOut$CAL_bins <- CAL_bins
  StockOut$nCALbins <- nCALbins
  
  StockOut$initdist <- initdist 
  StockOut$mov <- mov
  StockOut$Pinitdist <- Pinitdist
  StockOut$Asize <- Asize 
  StockOut$nareas <- nareas 
  
  return(StockOut)
}


#' Sample Fleet Parameters
#'
#' @param Fleet An object of class 'Fleet' or class 'OM' 
#' @param Stock An object of class 'Stock' or a list of sampled Stock parameters. 
#' Ignored if 'Fleet' is class 'OM'
#' @param nsim Number of simulations. Ignored if 'Fleet' is class 'OM'
#' @param nyears Number of historical years. Ignored if 'Fleet' is class 'OM'
#' @param proyears Number of projection years. Ignored if 'Fleet' is class 'OM'
#' @param cpars Optional named list of custom parameters. Ignored if 'Fleet' is class 'OM'
#' @param msg Logical. Print messages?
#' 
#' @keywords internal
#' 
#' @return A named list of sampled Fleet parameters
#' @export
#'
SampleFleetPars <- function(Fleet, Stock=NULL, nsim=NULL, nyears=NULL, 
                            proyears=NULL, cpars=NULL, msg=TRUE) {
  if (class(Fleet) != "Fleet" & class(Fleet) != "OM") 
    stop("First argument must be class 'Fleet' or 'OM'")
  
  if (class(Fleet) != "OM" & class(Stock) != "Stock" & class(Stock) != "list") 
    stop("Must provide 'Stock' object", call.=FALSE)
  
  # Get custom pars if they exist
  if (class(Fleet) == "OM" && length(Fleet@cpars) > 0 && is.null(cpars)) 
    cpars <- SampleCpars(Fleet@cpars, Fleet@nsim, msg=msg)  # custom parameters exist in Stock/OM object
  if (length(cpars) > 0) { # custom pars exist - assign to function environment 
    Names <- names(cpars)
    for (X in 1:length(Names)) assign(names(cpars)[X], cpars[[X]])
  }
  
  if (class(Fleet) == "OM") {
    nsim <- Fleet@nsim
    nyears <- Fleet@nyears 
    proyears <- Fleet@proyears
    StockPars <- SampleStockPars(Fleet, nsim, nyears, proyears, cpars, msg=msg)
    for (X in 1:length(StockPars)) assign(names(StockPars)[X], StockPars[[X]])
  }
  
  Fleet <- updateMSE(Fleet) # update to add missing slots with default values
  
  if (class(Stock) == "Stock") {
    Stock <- updateMSE(Stock) # update to add missing slots with default values
    # Sample Stock Pars - need some to calculate selectivity at age and length  
    StockPars <- SampleStockPars(Stock, nsim, nyears, proyears, cpars, msg=msg)
    for (X in 1:length(StockPars)) assign(names(StockPars)[X], StockPars[[X]])
  } 
  if (class(Stock) == "list") for (X in 1:length(Stock)) 
    assign(names(Stock)[X], Stock[[X]])
  
  Fleetout <- list()
  n_age <- maxage + 1 
  
  # --- Sample Historical Fishing Effort ----
  if (!exists("Esd", inherits = FALSE)) 
    Esd <- myrunif(nsim, Fleet@Esd[1], Fleet@Esd[2])
  if (!exists("EffLower", inherits = FALSE)) EffLower <- Fleet@EffLower
  if (!exists("EffUpper", inherits = FALSE)) EffUpper <- Fleet@EffUpper 
  if (!exists("EffYears", inherits = FALSE)) EffYears <- Fleet@EffYears
  
  if (!exists("Find", inherits = FALSE)) {
    if (any(is.na(EffLower)) || any(is.na(EffUpper)) || any(is.na(EffYears))) {
      message("NAs in EffLower, EffUpper, or EffYears")
      Find <- matrix(NA, nsim, nyears)
      Deriv <- list(Find, rep(NA, nsim))
    } else {
      Deriv <- getEffhist(Esd, nyears, EffYears = EffYears, EffLower = EffLower, 
                          EffUpper = EffUpper)  # Historical fishing effort  
      Find <- Deriv[[1]]  # Calculate fishing effort rate    
    }
  }
  
  if (!exists("dFfinal", inherits = FALSE)) {
    if (exists("Deriv", inherits = FALSE)) {
      dFfinal <- Deriv[[2]]  # Final gradient in fishing effort yr-1 
    } else {
      dFfinal <- rep(NA, nsim)
    }
  }
  if (any(dim(Find) != c(nsim, nyears))) 
    stop("Find must be matrix with dimensions: nsim (", nsim, "), nyears (", nyears, ") but is: ", paste(dim(Find), ""))
  
  Fleetout$Esd <- Esd
  Fleetout$Find <- Find
  Fleetout$dFfinal <- dFfinal
  
  # === Spatial Targetting ====
  # spatial targetting Ba^targetting param 
  if (!exists("Spat_targ", inherits = FALSE)) {
    Fleetout$Spat_targ <- Spat_targ <- myrunif(nsim, Fleet@Spat_targ[1], Fleet@Spat_targ[2])  
  } else {
    Fleetout$Spat_targ <- Spat_targ 
  }
  
  # === Sample fishing efficiency parameters ====
  # interannual variability in catchability
  if (!exists("qinc", inherits = FALSE)) 
    qinc <- myrunif(nsim, Fleet@qinc[1], Fleet@qinc[2])
  if (!exists("qcv", inherits = FALSE)) 
    qcv <- myrunif(nsim, Fleet@qcv[1], Fleet@qcv[2])  
  
  # === Simulate future variability in fishing efficiency ====
  qmu <- -0.5 * qcv^2  # Mean
  if (!exists("qvar", inherits = FALSE)) qvar <- array(exp(rnorm(proyears * nsim, rep(qmu, proyears), rep(qcv, proyears))), c(nsim, proyears))  # Variations in interannual variation
  FinF <- Find[, nyears]  # Effort in final historical year
  
  Fleetout$qinc <- qinc
  Fleetout$qcv <- qcv
  Fleetout$qvar <- qvar
  Fleetout$FinF <- FinF
  
  # ---- Selectivity Curve ----
  if (exists("V", inherits=FALSE) | 
      exists("SLarray", inherits=FALSE) |
      exists("retA", inherits=FALSE) | 
      exists("retL", inherits=FALSE)
      ) {
    Fleet@isRel <- 'FALSE'
  }
  
  # are selectivity parameters relative to size at maturity?
  chk <- class(Fleet@isRel)
  if (length(Fleet@isRel) < 1) Fleet@isRel <- "true"
  if (chk == "character") {
    chkRel <- substr(tolower(Fleet@isRel), 1,1)
    if (chkRel == "t" | Fleet@isRel == "1") multi <- L50
    if (chkRel == "f" | Fleet@isRel == "0") multi <- 1
  }
  if (chk == "numeric") {
    if (Fleet@isRel == 1) multi <- L50
    if (Fleet@isRel == 0) multi <- 1
  }
  if (chk == "logical") {
    if (Fleet@isRel) multi <- L50
    if (!Fleet@isRel) multi <- 1
  }
  
  if (exists("L5", inherits = FALSE) | exists("LFS", inherits = FALSE) | 
      exists("Vmaxlen", inherits = FALSE) | exists("V", inherits=FALSE)) {
    if (all(multi != 1))
      stop("Selectivity parameters provided in cpars must be absolute values. Is Fleet@isRel == 'FALSE'?")
  }
  
  if (!exists("L5", inherits = FALSE)) L5 <- myrunif(nsim, Fleet@L5[1], Fleet@L5[2]) * multi  # length at 5% selectivity ascending
  if (!exists("LFS", inherits = FALSE)) LFS <- myrunif(nsim, Fleet@LFS[1], Fleet@LFS[2]) * multi  # first length at 100% selection
  if (!exists("Vmaxlen", inherits = FALSE)) Vmaxlen <- myrunif(nsim, Fleet@Vmaxlen[1], Fleet@Vmaxlen[2])  # selectivity at maximum length
  
  L5_y <- matrix(L5, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
  LFS_y <- matrix(LFS, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
  Vmaxlen_y <- matrix(Vmaxlen, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
  
  # time-varying selectivity in OM
  SelYears <- Fleet@SelYears
  Selnyears <- length(Fleet@SelYears)
  if (Selnyears > 0) {
    L5s <- mapply(runif, n = nsim, min = Fleet@L5Lower, max = Fleet@L5Upper) * multi
    # first length at 100% selection
    LFSs <- mapply(runif, n = nsim, min = Fleet@LFSLower, max = Fleet@LFSUpper) *  multi
    ind <- which(LFSs/ matrix(Linf, nrow=nsim, ncol=Selnyears) > 1, arr.ind = T)
    if (length(ind) > 0) {
      message("LFS too high (LFS > Linf) in some cases. \nDefaulting to LFS = 0.9 Linf for the affected simulations")
      LFSs[ind] <- Linf[ind[, 1]] * 0.9
    } 
    
    # selectivity at maximum length
    Vmaxlens <- mapply(runif, n = nsim, min = Fleet@VmaxLower, max = Fleet@VmaxUpper)
    
    # update selectivity parameters 
    for (yy in seq_along(SelYears)) {
      if (yy < length(SelYears)) {
        yrind <- SelYears[yy]:(SelYears[yy+1]-1)  
      } else {
        yrind <- SelYears[yy]:(nyears+proyears)
      }
      L5_y[yrind,] <- matrix(L5s[,yy], nrow=length(yrind), ncol=nsim, byrow = TRUE)
      LFS_y[yrind,] <- matrix(LFSs[,yy], nrow=length(yrind), ncol=nsim, byrow = TRUE) 
      Vmaxlen_y[yrind,] <- matrix(Vmaxlens[,yy], nrow=length(yrind), ncol=nsim, byrow = TRUE)
    }
  }
  
  if (exists("SLarray", inherits = FALSE)) {
    # update selectivity parameters 
    if (exists("V", inherits = FALSE))  stop("Cannot pass both SLarray and V in cpars")
    nbins <- length(CAL_binsmid)
    if (any(dim(SLarray) != c(nsim, nbins, nyears+proyears)))
      stop("SLarray must be dimensions c(nsim, length(CAL_binsmid), nyears+proyears)")
    for (yr in 1:(nyears+proyears)) {
      b_ind <- apply(apply(SLarray[,,yr]>=0.05, 1, which), 2, min)
      L5_y[yr,] <- CAL_binsmid[b_ind]
      b_ind <- apply(SLarray[,,yr], 1, which.max)
      LFS_y[yr,] <- CAL_binsmid[b_ind]
      
      temp <- abs(replicate(nsim, CAL_binsmid) - Linf)
      b_ind <- apply(temp, 2, which.min)
      Vmaxlen_y[yr,] <- SLarray[,b_ind,yr][1,]
    }
  } else {
    if (exists("V", inherits = FALSE)) {
      # update selectivity parameters
      if (dim(V)[3] == nyears) {
        Dims <- dim(V)
        v2 <- array(V[,,nyears], dim=c(Dims[1], Dims[2], proyears))
        V <- abind::abind(V, v2, along=3)
      }
      if(any(dim(V)!= c(nsim, n_age, nyears + proyears)))
        stop('V must be dimensions: nsim, n_age, nyears + proyears')

      VB <- function(Linf, K, t0, age) Linf * (1-exp(-K*(age-t0)))
      for (yr in 1:(nyears+proyears)) {
        for (s in 1:nsim) {
          xout <- seq(1, n_age, by=0.1)
          tt <- approx(V[s,,yr], xout=xout)
          age5 <- tt$x[min(which(tt$y >=0.05))]-1
          L5_y[yr,s] <- VB(Linfarray[s,yr], Karray[s,yr], t0array[s,yr], age5)
          ageFS <- tt$x[which.max(tt$y)]-1
          if (ageFS == age5) ageFS <- age5 + 1
          LFS_y[yr, s] <- VB(Linfarray[s,yr], Karray[s,yr], t0array[s,yr], ageFS)
          Vmaxlen_y[yr, s] <- V[s, n_age, yr]
        }
      }
    }
    
    # calculate SLarray 
    nCALbins <- length(CAL_binsmid)
    CAL_binsmidMat <- matrix(CAL_binsmid, nrow=nsim, ncol=length(CAL_binsmid), byrow=TRUE)
    SLarray <- array(NA, dim=c(nsim, nCALbins, nyears+proyears)) # Selectivity-at-length 
    Vmaxlen_y[Vmaxlen_y<=0] <- tiny 
    for (yr in 1:(nyears+proyears)) {
      srs <- (Linf - LFS_y[yr,]) / ((-log(Vmaxlen_y[yr,],2))^0.5) 
      srs[!is.finite(srs)] <- Inf
      sls <- (LFS_y[yr,] - L5_y[yr, ]) /((-log(0.05,2))^0.5)
      SLarray[,, yr] <- t(sapply(1:nsim, getsel, lens=CAL_binsmidMat, lfs=LFS_y[yr, ], sls=sls, srs=srs))
    }
  }
  # Check LFS is greater than L5 
  chk <- sum(apply(L5_y > LFS_y, 2, prod) != 0)
  if (chk > 0) stop("L5 is greater than LFS in ", chk, ' simulations')
  
  if (!exists("V", inherits = FALSE)) {
    # calculate selectivity-at-age from selectivity-at-length
    VList <- lapply(1:nsim, calcV, Len_age=Len_age, LenCV=LenCV, SLarray=SLarray, 
                    n_age=n_age, nyears=nyears, proyears=proyears, CAL_binsmid=CAL_binsmid)
    V <- aperm(array(as.numeric(unlist(VList, use.names=FALSE)), dim=c(n_age, nyears+proyears, nsim)), c(3,1,2))
  }
  
  
  # ---- Retention Curve ---- 
  if(!exists("LR5", inherits = FALSE)) LR5 <- runif(nsim, min(Fleet@LR5), max(Fleet@LR5)) * multi
  if(!exists("LFR", inherits = FALSE)) LFR <- runif(nsim, min(Fleet@LFR), max(Fleet@LFR)) * multi
  if(!exists("Rmaxlen", inherits = FALSE)) Rmaxlen <- runif(nsim, min(Fleet@Rmaxlen), max(Fleet@Rmaxlen))
  if(!exists("DR", inherits = FALSE)) DR <- runif(nsim, min(Fleet@DR), max(Fleet@DR))
  
  LR5_y <- matrix(LR5, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
  LFR_y <- matrix(LFR, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
  Rmaxlen_y <- matrix(Rmaxlen, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
  DR_y <- matrix(DR, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)

  if (exists("retL", inherits=FALSE)) {
    # update retention parameters
    if (exists("retA", inherits = FALSE))  stop("Cannot pass both retL and retA in cpars")
    
    nbins <- length(CAL_binsmid)
    if (any(dim(retL) != c(nsim, nbins, nyears+proyears)))
      stop("retL must be dimensions c(nsim, length(CAL_binsmid), nyears+proyears)")
    for (yr in 1:(nyears+proyears)) {
      b_ind <- apply(apply(retL[,,yr]>=0.05, 1, which), 2, min)
      LR5_y[yr,] <- CAL_binsmid[b_ind]
      b_ind <- apply(retL[,,yr], 1, which.max)
      LFR_y[yr,] <- CAL_binsmid[b_ind]
      
      temp <- abs(replicate(nsim, CAL_binsmid) - Linf)
      b_ind <- apply(temp, 2, which.min)
      Rmaxlen_y[yr,] <- retL[,b_ind,yr][1,]
    }
  } else {
    if (exists("retA", inherits = FALSE)) {
      # update retention parameters
      if(any(dim(retA)!= c(nsim, n_age, nyears + proyears)))
        stop('retA must be dimensions: nsim, n_age, nyears + proyears')
      VB <- function(Linf, K, t0, age) Linf * (1-exp(-K*(age-t0)))
      for (yr in 1:(nyears+proyears)) {
        for (s in 1:nsim) {
          xout <- seq(1, n_age, by=0.1)
          tt <- approx(retA[s,,yr], xout=xout)
          ageR5 <- tt$x[min(which(tt$y >=0.05))]-1
          LR5_y[yr,s] <- VB(Linfarray[s,yr], Karray[s,yr], t0array[s,yr], ageR5)
          ageFR <- tt$x[which.max(tt$y)]-1
          if (ageFR == ageR5) ageFR <- ageR5 + 1
          LFR_y[yr, s] <- VB(Linfarray[s,yr], Karray[s,yr], t0array[s,yr], ageFR)
          Rmaxlen_y[yr, s] <- retA[s, n_age, yr]
        }
      }
    }
   # calculate retL
    nCALbins <- length(CAL_binsmid)
    CAL_binsmidMat <- matrix(CAL_binsmid, nrow=nsim, ncol=length(CAL_binsmid), byrow=TRUE)
    retL <- array(NA, dim=c(nsim, nCALbins, nyears+proyears)) # Retention-at-length 
    Rmaxlen_y[Rmaxlen_y<=0] <- tiny 
    if (any(LR5_y > LFR_y)) {
      if (all(LFR_y<0.001)) {
        LFR_y <- LR5_y + LFR_y
      } else {
        stop('LR5 is greater than LFR', call.=FALSE)   
      }
    }
    for (yr in 1:(nyears+proyears)) {
      srs <- (Linf - LFR_y[yr,]) / ((-log(Rmaxlen_y[yr,],2))^0.5) 
      srs[!is.finite(srs)] <- Inf
      sls <- (LFR_y[yr,] - LR5_y[yr, ]) /((-log(0.05,2))^0.5)
      retL[,, yr] <- t(sapply(1:nsim, getsel, lens=CAL_binsmidMat, lfs=LFR_y[yr, ], sls=sls, srs=srs))
    }
  }

  if (!exists("retA", inherits = FALSE)) {
    # calculate selectivity-at-age from selectivity-at-length
    retAList <- lapply(1:nsim, calcV, Len_age=Len_age, LenCV=LenCV, SLarray=retL, 
                    n_age=n_age, nyears=nyears, proyears=proyears, CAL_binsmid=CAL_binsmid)
    retA <- aperm(array(as.numeric(unlist(retAList, use.names=FALSE)), dim=c(n_age, nyears+proyears, nsim)), c(3,1,2))
  }
    
  V2 <- V
  SLarray2 <- SLarray
  
 
  # Apply general discard rate 
  dr <- aperm(abind::abind(rep(list(DR_y), n_age), along=3), c(2,3,1))
  retA <- (1-dr) * retA
  
  dr <- aperm(abind::abind(rep(list(DR_y), nCALbins), along=3), c(2,3,1))
  retL <- (1-dr) * retL
  
  # update realized vulnerability curve with retention and dead discarded fish 
  Fdisc_array1 <- array(Fdisc, dim=c(nsim, n_age, nyears+proyears))
  V <- V * (retA + (1-retA)*Fdisc_array1) # Realised selection at age
  
  Fdisc_array2 <- array(Fdisc, dim=c(nsim, nCALbins, nyears+proyears))
  SLarray <- SLarray2 * (retL + (1-retL)*Fdisc_array2) # Realised selection at length
  
  # Realised Retention curves
  retA <- retA * V2
  retL <- retL * SLarray2
  
  Fleetout$Fdisc <- Fdisc
  Fleetout$Fdisc_array1 <- Fdisc_array1
  Fleetout$Fdisc_array2 <- Fdisc_array2
  Fleetout$LR5 <- LR5_y  
  Fleetout$LFR <- LFR_y
  Fleetout$Rmaxlen <- Rmaxlen_y
  Fleetout$DR <- DR_y
  
  Fleetout$retA <- retA  # retention-at-age array - nsim, maxage, nyears+proyears
  Fleetout$retL <- retL  # retention-at-length array - nsim, nCALbins, nyears+proyears
  
  Fleetout$L5 <- L5_y  
  Fleetout$LFS <- LFS_y 
  Fleetout$Vmaxlen <- Vmaxlen_y
  Fleetout$V <- V  # realized vulnerability-at-age
  Fleetout$SLarray <- SLarray # realized vulnerability-at-length
  Fleetout$V2 <- V2 # original vulnerablity-at-age curve 
  Fleetout$SLarray2 <- SLarray2 # original vulnerablity-at-length curve 
  

  
  # check V 
  if (sum(apply(V, c(1,3), max) <0.01)) {
    maxV <- apply(V, c(1,3), max)
    fails <- which(maxV < 0.01, arr.ind = TRUE)
    sims <- unique(fails[,1])
    yrs <- unique(fails[,2])
    warning("Vulnerability (V) is <0.01 for all ages in:\nsims:", sims, "\nyears:", yrs, "\n", call. = FALSE)
    # warning('Check selectivity parameters. Is Fleet@isRel set correctly?', call.=FALSE)
  }
  
  Fleetout 
}

#' Sample Observation Parameters
#'
#' @param Obs An object of class 'Obs' or class 'OM'
#' @param nsim Number of simulations. Ignored if 'Obs' is class 'OM'
#' @param cpars Optional named list of custom parameters. Ignored if 'OM' is class 'OM'
#' 
#' @keywords internal
#' @return A named list of sampled Observation parameters
#' @export
#'
SampleObsPars <- function(Obs, nsim=NULL, cpars=NULL){
  if (class(Obs) != "Obs" & class(Obs) != "OM") 
    stop("First argument must be class 'Obs' or 'OM'")
  if (class(Obs) == "OM") nsim <- Obs@nsim
  
  # Get custom pars if they exist
  if (class(Obs) == "OM" && length(Obs@cpars) > 0 && is.null(cpars)) 
    cpars <- SampleCpars(Obs@cpars, Obs@nsim)  # custom parameters exist in OM object
  if (length(cpars) > 0) { # custom pars exist - assign to function environment 
    Names <- names(cpars)
    for (X in 1:length(Names)) assign(names(cpars)[X], cpars[[X]])
  }
  
  Obs <- updateMSE(Obs) # update to add missing slots with default values
  
  ObsOut <- list() 
  
  # === Sample observation error model parameters ====
  
  if (exists("Cobs", inherits = FALSE)) Csd <- Cobs
  if (!exists("Csd", inherits=FALSE)) {
    ObsOut$Csd <- myrunif(nsim, Obs@Cobs[1], Obs@Cobs[2])  # Sampled catch observation error (lognormal sd)
  } else {
    ObsOut$Csd <- Csd
  }

  if (!exists("Cbias", inherits=FALSE)) {
    ObsOut$Cbias <- rlnorm(nsim, mconv(1, Obs@Cbiascv), sdconv(1, Obs@Cbiascv))  # Sampled catch bias (log normal sd)
  } else {
    ObsOut$Cbias <- Cbias
  }
  if (!exists("CAA_nsamp", inherits=FALSE)) {
    ObsOut$CAA_nsamp <- ceiling(myrunif(nsim, Obs@CAA_nsamp[1], Obs@CAA_nsamp[2]))  # Number of catch-at-age observations
  } else {
    ObsOut$CAA_nsamp <- CAA_nsamp
  } 
  if (!exists("CAA_ESS", inherits=FALSE)) {
    ObsOut$CAA_ESS <- ceiling(myrunif(nsim, Obs@CAA_ESS[1], Obs@CAA_ESS[2]))  # Effective sample size
  } else {
    ObsOut$CAA_ESS <- CAA_ESS
  } 
  if (!exists("CAL_nsamp", inherits=FALSE)) {
    ObsOut$CAL_nsamp <- ceiling(myrunif(nsim, Obs@CAL_nsamp[1], Obs@CAL_nsamp[2]))  # Observation error standard deviation for single catch at age by area
  } else {
    ObsOut$CAL_nsamp <- CAL_nsamp
  }  
  if (!exists("CAL_ESS", inherits=FALSE)) {
    ObsOut$CAL_ESS <- ceiling(myrunif(nsim, Obs@CAL_ESS[1], Obs@CAL_ESS[2]))  # Effective sample size
  } else {
    ObsOut$CAL_ESS <- CAL_ESS
  }  
  
  if (exists("beta", inherits = FALSE)) betas <- beta
  if (!exists("betas", inherits=FALSE)) {
    ObsOut$betas <- exp(myrunif(nsim, log(Obs@beta[1]), log(Obs@beta[2])))  # the sampled hyperstability / hyperdepletion parameter beta>1 (hyperdepletion) beta<1 (hyperstability)
  } else {
    ObsOut$betas <- betas
  }  
  
  if (exists("Iobs", inherits = FALSE)) Isd <- Iobs
  if (!exists("Isd", inherits=FALSE)) {
    ObsOut$Isd <- myrunif(nsim, Obs@Iobs[1], Obs@Iobs[2])  # Abundance index observation error (log normal sd)
  } else {
    ObsOut$Isd <- Isd
  } 
  if (exists("Dobs", inherits = FALSE)) Derr <- Dobs
  if (!exists("Derr", inherits=FALSE)) {
    ObsOut$Derr <- myrunif(nsim, Obs@Dobs[1], Obs@Dobs[2])
  } else {
    ObsOut$Derr <- Derr
  }  
  if (!exists("Dbias", inherits=FALSE)) {
    ObsOut$Dbias <- rlnorm(nsim, mconv(1, Obs@Dbiascv), sdconv(1, Obs@Dbiascv))  # sample of depletion bias
  } else {
    ObsOut$Dbias <- Dbias
  }  
  if (!exists("Mbias", inherits=FALSE)) {
    ObsOut$Mbias <- rlnorm(nsim, mconv(1, Obs@Mbiascv), sdconv(1, Obs@Mbiascv))  # sample of M bias
  } else {
    ObsOut$Mbias <- Mbias
  }  
  if (!exists("FMSY_Mbias", inherits=FALSE)) {
    ObsOut$FMSY_Mbias <- rlnorm(nsim, mconv(1, Obs@FMSY_Mbiascv), sdconv(1, Obs@FMSY_Mbiascv))  # sample of FMSY/M bias
  } else {
    ObsOut$FMSY_Mbias <- FMSY_Mbias
  }  
  if (!exists("lenMbias", inherits=FALSE)) {
    ObsOut$lenMbias <- rlnorm(nsim, mconv(1, Obs@LenMbiascv), sdconv(1, Obs@LenMbiascv))  # sample of length at maturity bias - assume same error as age based maturity
  } else {
    ObsOut$lenMbias <- lenMbias
  } 
  if (!exists("LFCbias", inherits=FALSE)) {
    ObsOut$LFCbias <- rlnorm(nsim, mconv(1, Obs@LFCbiascv), sdconv(1, Obs@LFCbiascv))  # sample of length at first capture bias
  } else {
    ObsOut$LFCbias <- LFCbias
  } 
  if (!exists("LFSbias", inherits=FALSE)) {
    ObsOut$LFSbias <- rlnorm(nsim, mconv(1, Obs@LFSbiascv), sdconv(1, Obs@LFSbiascv))  # sample of length at full selection bias
  } else {
    ObsOut$LFSbias <- LFSbias
  }
  
  if (exists("Btobs", inherits = FALSE))  Aerr <- Btobs
  if (!exists("Aerr", inherits=FALSE)) {
    ObsOut$Aerr <- myrunif(nsim, Obs@Btobs[1], Obs@Btobs[2])
  } else {
    ObsOut$Aerr <- Aerr
  }
  if (exists("Btbiascv", inherits = FALSE)) {
    Abias <- Btbiascv
  }
  if (!exists("Abias", inherits=FALSE)) {
    ObsOut$Abias <- exp(myrunif(nsim, log(Obs@Btbiascv[1]), log(Obs@Btbiascv[2])))  #rlnorm(nsim,mconv(1,Obs@Btbiascv),sdconv(1,Obs@Btbiascv))    # sample of current abundance bias
  } else {
    ObsOut$Abias <- Abias
  }
  if (!exists("Kbias", inherits=FALSE)) {
    ObsOut$Kbias <- rlnorm(nsim, mconv(1, Obs@Kbiascv), sdconv(1, Obs@Kbiascv))  # sample of von B. K parameter bias
  } else {
    ObsOut$Kbias <- Kbias
  }
  if (!exists("t0bias", inherits=FALSE)) {
    ObsOut$t0bias <- rlnorm(nsim, mconv(1, Obs@t0biascv), sdconv(1, Obs@t0biascv))  # sample of von B. t0 parameter bias
  } else {
    ObsOut$t0bias <- t0bias
  } 
  if (!exists("Linfbias", inherits=FALSE)) {
    ObsOut$Linfbias <- rlnorm(nsim, mconv(1, Obs@Linfbiascv), sdconv(1, Obs@Linfbiascv))  # sample of von B. maximum length bias
  } else {
    ObsOut$Linfbias <- Linfbias
  } 
  if (!exists("Irefbias", inherits=FALSE)) {
    ObsOut$Irefbias <- rlnorm(nsim, mconv(1, Obs@Irefbiascv), sdconv(1, Obs@Irefbiascv))  # sample of bias in reference (target) abundance index
  } else {
    ObsOut$Irefbias <- Irefbias
  } 
  if (!exists("Crefbias", inherits=FALSE)) {
    ObsOut$Crefbias <- rlnorm(nsim, mconv(1, Obs@Crefbiascv), sdconv(1, Obs@Crefbiascv))  # sample of bias in reference (target) catch index
  } else {
    ObsOut$Crefbias <- Crefbias
  }
  if (!exists("Brefbias", inherits=FALSE)) {
    ObsOut$Brefbias <- rlnorm(nsim, mconv(1, Obs@Brefbiascv), sdconv(1, Obs@Brefbiascv))  # sample of bias in reference (target) biomass index
  } else {
    ObsOut$Brefbias <- Brefbias
  }
  if (!exists("Recsd", inherits=FALSE)) {
    ObsOut$Recsd <- myrunif(nsim, Obs@Recbiascv[1], Obs@Recbiascv[2])  # Recruitment deviation  
  } else {
    ObsOut$Recsd <- Recsd
  }  
  
  # ObsOut$CALcv <- runif(nsim, Obs@CALcv[1], Obs@CALcv[2])  # Observation error standard deviation for single catch at age by area
  # ObsOut$LenCVbias <- rlnorm(nsim, mconv(1, Obs@CALcv), sdconv(1, Obs@CALcv)) # sample of bias in assumed CV of catch-at-length
  
  ObsOut
}

#' Sample Implementation Error Parameters
#'
#' @param Imp An object of class 'Imp' or class 'OM'
#' @param nsim Number of simulations. Ignored if 'Imp' is class 'OM'
#' @param cpars Optional named list of custom parameters. Ignored if 'OM' is class 'OM'
#' @return A named list of sampled Implementation Error parameters
#' @keywords internal
#' @export
#'
SampleImpPars <- function(Imp, nsim=NULL, cpars=NULL) {
  if (class(Imp) != "Imp" & class(Imp) != "OM") 
    stop("First argument must be class 'Imp' or 'OM'")
  if (class(Imp) == "OM") nsim <- Imp@nsim
  
  # Get custom pars if they exist
  if (class(Imp) == "OM" && length(Imp@cpars) > 0 && is.null(cpars)) 
    cpars <- SampleCpars(Imp@cpars, Imp@nsim)  # custom parameters exist in OM object
  if (length(cpars) > 0) { # custom pars exist - assign to function environment 
    Names <- names(cpars)
    for (X in 1:length(Names)) assign(names(cpars)[X], cpars[[X]])
  }
  
  ImpOut <- list() 
  # === Sample implementation error parameters ====
  if (!exists("TACSD", inherits = FALSE)) {
    ImpOut$TACSD <- myrunif(nsim, Imp@TACSD[1], Imp@TACSD[2])  # Sampled TAC error (lognormal sd)
  } else {
    ImpOut$TACSD <- TACSD
  }
  if (!exists("TACFrac", inherits = FALSE)) {
    ImpOut$TACFrac <- myrunif(nsim, Imp@TACFrac[1], Imp@TACFrac[2])  # Sampled TAC fraction (log normal sd)
  } else {
    ImpOut$TACFrac <- TACFrac
  }
  if (!exists("TAESD", inherits = FALSE)) {
    ImpOut$TAESD <- myrunif(nsim, Imp@TAESD[1], Imp@TAESD[2])  # Sampled Effort error (lognormal sd)
  } else {
    ImpOut$TAESD <- TAESD
  }
  if (!exists("TAEFrac", inherits = FALSE)) {
    ImpOut$TAEFrac <- myrunif(nsim, Imp@TAEFrac[1], Imp@TAEFrac[2])  # Sampled Effort fraction (log normal sd)
  } else {
    ImpOut$TAEFrac <- TAEFrac
  }
  if (!exists("SizeLimSD", inherits = FALSE)) {
    ImpOut$SizeLimSD<-myrunif(nsim,Imp@SizeLimSD[1],Imp@SizeLimSD[2])
  } else {
    ImpOut$SizeLimSD <- SizeLimSD
  }
  if (!exists("SizeLimFrac", inherits = FALSE)) {
    ImpOut$SizeLimFrac<-myrunif(nsim,Imp@SizeLimFrac[1],Imp@SizeLimFrac[2])
  } else {
    ImpOut$SizeLimFrac <- SizeLimFrac
  }

  ImpOut
}

# #' Sample Bio-Economic Parameters
# #'
# #' @param BioEco An object of class 'BioEco' or class 'OM'
# #' @param nsim Number of simulations. Ignored if 'BioEco' is class 'OM'
# #' @param cpars Optional named list of custom parameters. Ignored if 'OM' is class 'OM'
# #' @return A named list of sampled Bio-Economic parameters
# #' @keywords internal
# #' @export
# #'
# SampleBioEcoPars <- function(BioEco, nsim=NULL, cpars=NULL) {
#   if (class(BioEco) != "BioEco" & class(BioEco) != "OM") 
#     stop("First argument must be class 'BioEco' or 'OM'")
#   if (class(BioEco) == "OM") nsim <- BioEco@nsim
#   
#   # Get custom pars if they exist
#   if (class(BioEco) == "OM" && length(BioEco@cpars) > 0 && is.null(cpars)) 
#     cpars <- SampleCpars(BioEco@cpars, BioEco@nsim)  # custom parameters exist in OM object
#   if (length(cpars) > 0) { # custom pars exist - assign to function environment 
#     Names <- names(cpars)
#     for (X in 1:length(Names)) assign(names(cpars)[X], cpars[[X]])
#   }
#   
#   BioEcoOut <- list() 
#   
#   if (!exists("CostCurr", inherits = FALSE)) {
#     BioEcoOut$CostCurr <- myrunif(nsim, BioEco@CostCurr[1], BioEco@CostCurr[2]) 
#   } else {
#     BioEcoOut$CostCurr <- CostCurr
#   }
#   if (!exists("RevCurr", inherits = FALSE)) {
#     BioEcoOut$RevCurr <- myrunif(nsim, BioEco@RevCurr[1], BioEco@RevCurr[2]) 
#   } else {
#     BioEcoOut$RevCurr <- RevCurr
#   }
#   if (!exists("CostInc", inherits = FALSE)) {
#     BioEcoOut$CostInc <- myrunif(nsim, BioEco@CostInc[1], BioEco@CostInc[2]) 
#   } else {
#     BioEcoOut$CostInc <- CostInc
#   }
#   if (!exists("RevInc", inherits = FALSE)) {
#     BioEcoOut$RevInc <- myrunif(nsim, BioEco@RevInc[1], BioEco@RevInc[2]) 
#   } else {
#     BioEcoOut$RevInc <- RevInc
#   }
#   if (!exists("Response", inherits = FALSE)) {
#     BioEcoOut$Response <- myrunif(nsim, BioEco@Response[1], BioEco@Response[2]) 
#   } else {
#     BioEcoOut$Response <- Response
#   }
#   if (!exists("LatentEff", inherits = FALSE)) {
#     if (length(BioEco@LatentEff) ==  0) {
#       BioEcoOut$LatentEff <- rep(NA, nsim)
#     } else {
#       if (any(BioEco@LatentEff<=0)) stop("LatentEff must be fraction > 0 and <= 1")
#       if (any(BioEco@LatentEff>1)) stop("LatentEff must be fraction > 0 and <= 1")
#       BioEcoOut$LatentEff <- myrunif(nsim, BioEco@LatentEff[1], BioEco@LatentEff[2])   
#     }
#   } else {
#     BioEcoOut$LatentEff <- LatentEff
#   }
#   BioEcoOut
# }


#' Valid custom parameters (cpars)
#'
#' @param type What cpars to show? 'all', 'Stock', 'Fleet', 'Obs', 'Imp', or 'internal'
#' @param valid Logical. Show valid cpars?
#'
#' @return a HTML datatable with variable name, description and type of valid cpars
#' @export
#' 
#' @examples
#' \dontrun{
#' validcpars() # all valid cpars
#' 
#' validcpars("Obs", FALSE) # invalid Obs cpars
#' }
#'
validcpars <- function(type=c("all", "Stock", "Fleet", "Obs", "Imp", "internal"),
                       valid=TRUE) {
  
  type <- match.arg(type, choices=c("all", "Stock", "Fleet", "Obs", "Imp", "internal"),
                    several.ok = TRUE )
  if ('all' %in% type) type <- c("Stock", "Fleet", "Obs", "Imp", "internal")
  
  Valid <- Slot <- Dim <- Description <- NULL
  
  # cpars_info <- DLMtool:::cpars_info
  cpars_info <- cpars_info[!duplicated(cpars_info$Slot),] # remove duplicated 'Name'
  
  cpars_info$type <- NA
  stock_ind <- match(slotNames("Stock"), cpars_info$Slot)
  fleet_ind <- match(slotNames("Fleet"), cpars_info$Slot)
  obs_ind <- match(slotNames("Obs"), cpars_info$Slot)
  imp_ind <- match(slotNames("Imp"), cpars_info$Slot)
  int_ind <- (1:nrow(cpars_info))[!1:nrow(cpars_info) %in% 
                       c(stock_ind, fleet_ind, obs_ind, imp_ind)]
  
  cpars_info$type[stock_ind] <- "Stock"
  cpars_info$type[fleet_ind] <- "Fleet"
  cpars_info$type[obs_ind] <- "Obs"
  cpars_info$type[imp_ind] <- "Imp"
  cpars_info$type[int_ind] <- "internal"
  
  dflist <- list(); count <- 0 
  for (ss in type) {
    count <- count + 1
    df <- cpars_info %>% dplyr::filter(Valid==valid, cpars_info$type %in% ss) %>%
      dplyr::select(Slot, Dim, Description, type)
    names(df) <- c("Var.", "Dim.", "Desc.", "Type")
    if (nrow(df)> 0) {
      if (nrow(df)>1) {
        dflist[[count]] <- df[,as.logical(apply(!apply(df, 2, is.na), 2, prod))]
      } else {
        dflist[[count]] <- df[,!is.na(df)]
      }
    } 
  }
  
  dfout <- do.call("rbind", dflist)
  if (is.null(dfout)) {
    if (valid) message("No valid  parameters")
    if (!valid) message("No invalid parameters")
  } 
  
  dfout$Type <- as.factor(dfout$Type)
  dfout$Var. <- as.factor(dfout$Var.)
  if (requireNamespace("DT", quietly = TRUE)) {
    return(DT::datatable(dfout, filter = 'top', options = list(
      columnDefs = list(list(searchable = FALSE, targets = c(2,3))),
      pageLength = 25, autoWidth = TRUE)))
  } else {
    message("Install package `DT` to display dataframe as HTML table")
    return(dfout)
  }
}





#' Sample custom pars
#'
#' @param cpars A named list containing custom parameters for the OM
#' @param nsim number of simulations
#' @param msg logical - print the names of the cpars? Turn off when using the function in a loop
#' @return A named list of sampled custom parameters
#' @keywords internal
#' @export
#'
SampleCpars <- function(cpars, nsim=48, msg=TRUE) {
  
  # Vector of valid names for custompars list or data.frame. Names not in this list will be printed out in warning and ignored #	
  # ParsNames <- validcpars(FALSE)
  
  # check Perr 
  #internal process error by simulation and year is now Perr_y instead of Perr 
  if ("Perr" %in% names(cpars)) {
    if (!is.null(dim(cpars[['Perr']]))) {
      cpars[['Perr_y']] <- cpars[['Perr']]
      cpars[['Perr']] <- NULL
    }
  }
  # cpars_info <- DLMtool:::cpars_info # get internal data from sysdata
  CparsInfo <- cpars_info # get internal data from sysdata
  
  
  sampCpars <- list()
  ncparsim<-cparscheck(cpars)
  if ('CAL_bins' %in% names(cpars)) {
    sampCpars$CAL_bins <- cpars$CAL_bins
  }
  if ('maxage' %in% names(cpars)) {
    sampCpars$maxage <- cpars$maxage
  }
  if ('binWidth' %in% names(cpars)) {
    sampCpars$binWidth <- cpars$binWidth
  }
  # if (is.null(ncparsim)) return(sampCpars)
  
  Names <- names(cpars)
  
  ValNames <- c(CparsInfo$Slot[which(CparsInfo$Valid>0)], CparsInfo$Legacy[which(CparsInfo$Valid>0)])
  ValNames <- ValNames[!is.na(ValNames)]
  InvalNames <- c(CparsInfo$Slot[!which(CparsInfo$Valid>0)], CparsInfo$Legacy[!which(CparsInfo$Valid>0)])
  InvalNames <- unique(InvalNames[!is.na(InvalNames)])
  
  # report invalid names 
  invalid <- Names[!Names %in% ValNames]
  if (length(invalid)>0) {
    invdf <- data.frame(name=invalid, action='ignoring', alt="", stringsAsFactors = FALSE)
    alt_inval <- invalid[invalid %in% InvalNames]
    if (length(alt_inval)>0) {
      alt <- CparsInfo$Description[match(alt_inval, CparsInfo$Slot)]
      invdf$alt[match(alt_inval, invdf$name)] <-alt
    }
    if(msg) {
      message("invalid names found in custom parameters (OM@cpars)")	
      message(paste0(capture.output(invdf), collapse = "\n"))
    }
  }
  # report found names
  valid <- which(Names %in% ValNames)
  cpars <- cpars[valid]
  if (length(valid) == 0) {
    message("No valid names found in custompars (OM@cpars). Ignoring `OM@cpars`")
    return(list())
  }
  
  Names <- names(cpars)
  outNames <- paste(Names, "")
  for (i in seq(5, by=5, length.out=floor(length(outNames)/5)))
    outNames <- gsub(outNames[i], paste0(outNames[i], "\n"), outNames)
  if(msg) message("valid custom parameters (OM@cpars) found: \n", paste0(outNames, collapse="\n"))
  
  # # report invalid names 
  # invalid <- which(!Names %in% ParsNames)
  # if (length(invalid) > 0) {
  #   outNames <- paste(Names[invalid], "")
  #   for (i in seq(5, by=5, length.out=floor(length(outNames)/5))) outNames <- gsub(outNames[i], paste0(outNames[i], "\n"), outNames)
  #   if(msg) message("ignoring invalid names found in custom parameters (OM@cpars) \n", outNames)	
  # }

  # Sample custom pars 
  if (!is.null(ncparsim)) {
    if (ncparsim < nsim) {
      ind <- sample(1:ncparsim, nsim, replace=TRUE)
    } else {
      if (ncparsim == nsim) {
        ind <- 1:nsim
      } else {
        ind <- sample(1:ncparsim, nsim, replace=FALSE)
      }
    }
   
  }
  
  # if (!ncparsim < nsim) ind <- sample(1:ncparsim, nsim, replace=FALSE)
  if ('Data' %in% names(cpars)) {
    sampCpars$Data <- cpars$Data
    cpars$Data <- NULL
    
  }
  if (length(cpars)>0) {
    for (i in 1:length(cpars)) {
      samps <- cpars[[i]]
      name <- names(cpars)[i]
      if (any(c("maxage", "M_at_Length", "CAL_binsmid", "CAL_bins", "binWidth", "AddIunits") %in% name)) {
        sampCpars[[name]] <- samps
      } else {
        if ("numeric" %in% class(samps) | "integer" %in% class(samps)) sampCpars[[name]] <- samps[ind]
        
        if ('matrix' %in% class(samps)| 'array' %in% class(samps)) {
          if (length(dim(samps)) == 2) {
            sampCpars[[name]] <- samps[ind,, drop=FALSE]   
          }  else {
            dims <- dim(samps)
            tout <- array(NA, dim=c(length(ind), dims[2:length(dims)]))
            tlist <- c(list(ind), lapply(dims[2:length(dims)], seq))
            tlist2 <- c(list(1:nsim), lapply(dims[2:length(dims)], seq))
            varind <- expand.grid(tlist) %>% as.matrix()
            varind2 <- expand.grid(tlist2) %>% as.matrix()
            tout[varind2] <- samps[varind]
            sampCpars[[name]] <- tout
          }
        }
        
        if ("data.frame" %in% class(samps))   sampCpars[[name]] <- samps 
      }
    }
  }

  sampCpars
}





