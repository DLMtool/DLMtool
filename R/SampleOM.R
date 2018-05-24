
#' Sample Stock parameters
#'
#' @param Stock An object of class 'Stock' or class 'OM'
#' @param nsim Number of simulations. Ignored if 'Stock' is class 'OM'
#' @param nyears Number of historical years. Ignored if 'Stock' is class 'OM'
#' @param proyears Number of projection years. Ignored if 'Stock' is class 'OM'
#' @param cpars Optional named list of custom parameters. Ignored if 'Stock' is class 'OM'
#' @param Msg logical. Warning message for M values?
#'
#' @return A named list of sampled Stock parameters
#' @export
#'   
SampleStockPars <- function(Stock, nsim=48, nyears=80, proyears=50, cpars=NULL, Msg=TRUE) {
  if (class(Stock) != "Stock" & class(Stock) != "OM") 
    stop("First argument must be class 'Stock' or 'OM'")
  Stock <- updateMSE(Stock) # update to add missing slots with default values
  if (all(is.na(Stock@LenCV))) Stock@LenCV <- c(0.1, 0.1)
  if (all(is.na(Stock@Mexp))) Stock@Mexp <- c(0, 0)
  
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
  if (length(Stock@M) == 2 & !exists("M", inherits=FALSE)) M <- runif(nsim, Stock@M[1], Stock@M[2])  # natural mortality rate
  
  if (length(Stock@M) == maxage) { # Stock@M is vector of M-at-age 
    if (length(Stock@M2) == maxage && !exists("Mage", inherits=FALSE)) {
      mmat <- rbind(Stock@M, Stock@M2)
      if (all(mmat[1,] < mmat[2,]) | all(mmat[1,] > mmat[2,])) {
        Mage <- matrix(NA, nsim, maxage)
        Mage[,1] <- runif(nsim, min(mmat[,1]), max(mmat[,1]))
        val <- (Mage[,1] - min(mmat[,1]))/ diff(mmat[,1])
        for (X in 2:maxage) Mage[,X] <- min(mmat[,X]) + diff(mmat[,X])*val  
      } else stop("All values in slot 'M' must be greater or less than corresponding values in slot 'M2'", call.=FALSE)
    } else stop("slot 'M2' must be length 'maxage'", call.=FALSE)
  } 
  if (length(Stock@M) != maxage & length(Stock@M) != 2) stop("slot 'M' must be either length 2 or length maxage", call.=FALSE)
  
  if (length(Stock@M2) == maxage & !length(Stock@M) == maxage) {
    stop("Slot M2 is used (upper bound on M-at-age) and is length 'maxage' but Slot M (lower bound on M-at-age) is not length 'maxage'.")
  }
  if (!exists("Msd", inherits=FALSE)) Msd <- runif(nsim, Stock@Msd[1], Stock@Msd[2])  # sample inter annual variability in M frStock specified range
  if (!exists("Mgrad", inherits=FALSE)) Mgrad <- runif(nsim, Stock@Mgrad[1], Stock@Mgrad[2])  # sample gradient in M (M y-1)
  if (.hasSlot(Stock, "Mexp") & !exists("Mexp", inherits=FALSE)) {
    if (all(is.numeric(Stock@Mexp) & is.finite(Stock@Mexp))) {
      Mexp <- runif(nsim, min(Stock@Mexp), max(Stock@Mexp)) # sample Lorenzen M-at-weight exponent     
    } else {
      Mexp <- rep(0, nsim) # assume constant M-at-age/size
    }
  } 
  
  if (!exists("M", inherits=FALSE)) M <- Mage[,maxage]
  
  if (!exists("Mexp", inherits=FALSE)) Mexp <- rep(0, nsim) # assume constant M-at-age/size if it is not specified 
  
  if (!all(Mexp == 0) & length(Stock@M2) == maxage) {
    stop("Values in both M2 and Mexp slots. Only one can be used")
  }

  # == Depletion ====
  if (!exists("D", inherits=FALSE)) {
    StockOut$D <- D <- runif(nsim, Stock@D[1], Stock@D[2])  # sample from the range of user-specified depletion (Bcurrent/B0)  
  } else {
    StockOut$D <- D 
  }
  
  # == Stock-Recruitment Relationship ====
  if (!exists("SRrel", inherits=FALSE)) {
    StockOut$SRrel <- rep(Stock@SRrel, nsim)  # type of Stock-recruit relationship. 1=Beverton Holt, 2=Ricker
  } else {
    StockOut$SRrel <- SRrel 
  }
  
  if (!exists("hs", inherits=FALSE)) {
    StockOut$hs <- hs <- runif(nsim, Stock@h[1], Stock@h[2])  # sample of recruitment compensation (steepness - fraction of unfished recruitment at 20% of unfished biStockass)
  } else {
    StockOut$hs <- hs
  }
  if (any(StockOut$hs > 1 | StockOut$hs < 0.2)) stop("Steepness (OM@h) must be between 0.2 and 1", call.=FALSE)
  
  # == Recruitment Deviations ====
  if (!exists("procsd", inherits=FALSE)) {
    StockOut$procsd <- procsd <- runif(nsim, Stock@Perr[1], Stock@Perr[2])  # Process error standard deviation
  } else {
    StockOut$procsd <- procsd
  }
  
  if (!exists("AC", inherits=FALSE)) {
    StockOut$AC <- AC <- runif(nsim, Stock@AC[1], Stock@AC[2])  # auto correlation parameter for recruitment deviations recdev(t)<-AC*recdev(t-1)+(1-AC)*recdev_proposed(t)  
  } else {
    StockOut$AC <- AC  # auto correlation parameter for recruitment deviations recdev(t)<-AC*recdev(t-1)+(1-AC)*recdev_proposed(t)
  }
  
  # All recruitment Deviations
  # Add cycle (phase shift) to recruitment deviations - if specified
  if (is.finite(Stock@Period[1]) & is.finite(Stock@Amplitude[1])) {
    # Shape <- "sin"  # default sine wave - alternative - 'shift' for step changes
    Period <- runif(nsim, min(Stock@Period), max(Stock@Period))
    if (max(Stock@Amplitude)>1) {
      message("Stock@Amplitude > 1. Defaulting to 1")
      Stock@Amplitude[Stock@Amplitude>1] <- 1
    }
    Amplitude <- runif(nsim, min(Stock@Amplitude), max(Stock@Amplitude))
    
    yrs <- 1:(nyears + proyears+maxage-1)
    recMulti <- t(sapply(1:nsim, function(x) 1+sin((runif(1, 0, 1)*max(yrs) + 2*yrs*pi)/Period[x])*Amplitude[x]))
    message("Adding cyclic recruitment pattern")
    
    # recMulti <-  t(sapply(1:nsim, SetRecruitCycle, Period, Amplitude, TotYears=length(yrs), Shape = "sin"))
    
  } else {
    recMulti <- 1 
  }
  
  StockOut$procmu <- procmu <- -0.5 * (procsd)^2  # adjusted log normal mean
  if (!exists("Perr", inherits=FALSE)) {
    Perr <- array(rnorm((nyears + proyears+maxage-1) * nsim, rep(procmu, nyears + proyears+maxage-1), 
                        rep(procsd, nyears + proyears+maxage-1)), c(nsim, nyears + proyears+maxage-1))
    for (y in 2:(nyears + proyears+maxage-1)) Perr[, y] <- AC * Perr[, y - 1] + Perr[, y] * (1 - AC * AC)^0.5  #2#AC*Perr[,y-1]+(1-AC)*Perr[,y] # apply a pseudo AR1 autocorrelation to rec devs (log space)
    StockOut$Perr <- Perr <- exp(Perr) * recMulti # normal space (mean 1 on average) 
    
    
  } else {
    StockOut$Perr <- Perr
  }
  
  # if (nsim > 1) {
  #   cumlRecDev <- apply(Perr[, 1:(nyears+maxage-1)], 1, prod)
  #   dep[order(cumlRecDev)] <- dep[order(dep, decreasing = F)]  # robustifies 
  # }
  
  # == Growth parameters ====
  vars <- c("Linf", "Linfsd", "Linfgrad", "K", "Ksd", "Kgrad", "t0")
  for (var in vars) {
    if (!exists(var, inherits=FALSE)) {
      if (all(is.na(slot(Stock, var)))) {
        val <- rep(0, nsim)
      } else {
        val <- runif(nsim, slot(Stock, var)[1],slot(Stock, var)[2])  # sample of asymptotic length
      }
      assign(var, val)
    } 

  }
    
  # == Sample Fecundity-Length Exponent ===
  # if (!exists("FecB", inherits=FALSE))   FecB <- runif(nsim, min(Stock@FecB), max(Stock@FecB))
  
  
  # == Sample Spatial Parameters ====
  if (!exists("Frac_area_1", inherits=FALSE)) Frac_area_1 <- runif(nsim, Stock@Frac_area_1[1], Stock@Frac_area_1[2])  # sampled fraction of unfished biStockass in area 1 (its a two area model by default)
  if (!exists("Prob_staying", inherits=FALSE)) Prob_staying <- runif(nsim, Stock@Prob_staying[1], Stock@Prob_staying[2])  # sampled probability of individuals staying in area 1 among years
  if (!exists("Size_area_1", inherits=FALSE)) Size_area_1 <- runif(nsim, Stock@Size_area_1[1], Stock@Size_area_1[2])  # currently redundant parameter for the habitat area size of area 1
  
  if (max(Size_area_1) == 0) stop("Size_area_1 must be > 0")
  if (max(Frac_area_1) == 0) stop("Frac_area_1 must be > 0")
  if (max(Prob_staying) == 0) stop("Prob_staying must be > 0")
  
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
  if (!exists("Linfarray", inherits=FALSE)) Linfarray <- gettempvar(Linf, Linfsd, Linfgrad, nyears + proyears, nsim, Linfrand)  # Linf array  
  if (!exists("Karray", inherits=FALSE)) Karray <- gettempvar(K, Ksd, Kgrad, nyears + proyears, nsim, Krand)  # the K array
  if (!exists("Agearray", inherits=FALSE))  Agearray <- array(rep(1:maxage, each = nsim), dim = c(nsim, maxage))  # Age array
  
  if (length(StockOut$maxage) > 1) StockOut$maxage <- StockOut$maxage[1] # check if maxage has been passed in custompars
  
  t0array <- matrix(t0, nrow=nsim, ncol=proyears+nyears)
  
  
  
  # === Create Mean Length-at-Age array ====
  if (!exists("Len_age", inherits=FALSE)) {
    Len_age <- array(NA, dim = c(nsim, maxage, nyears + proyears))  # Length at age array
    ind <- as.matrix(expand.grid(1:nsim, 1:maxage, 1:(nyears + proyears)))  # an index for calculating Length at age
    Len_age[ind] <- Linfarray[ind[, c(1, 3)]] * (1 - exp(-Karray[ind[, c(1, 3)]] * 
                                                           (Agearray[ind[, 1:2]] - t0[ind[, 1]])))
    
    if (length(Stock@cpars$Linf) >0) {
      maxLinf <- max(Stock@cpars$Linf)
    } else {
      maxLinf <- max(Stock@Linf)
    }
    linfs <- gettempvar(maxLinf, 0, max(Stock@Linfgrad), nyears + proyears, 
               1, matrix(1, nrow=1, ncol=proyears+nyears))
    MaxBin <- ceiling(max(linfs) + 3 * max(linfs) * max(Stock@LenCV)) 
    
  } else { # Len_age has been passed in with cpars
    if (any(dim(Len_age) != c(nsim, maxage, nyears + proyears))) 
      stop("'Len_age' must be array with dimensions: nsim, maxage, nyears + proyears") 
    # Estimate vB parameters for each year and each sim 
    if (!all(c("Linf", "K", "t0") %in% names(cpars))) { # don't calculate if Linf, K and t0 have also been passed in with cpars
      vB <- function(pars, ages) pars[1] * (1-exp(-pars[2]*(ages-pars[3])))
      fitVB <- function(pars, LatAge, ages) sum((vB(pars, ages) - LatAge)^2)
      starts <- c(max(Len_age), 0.2, 0)
      if(Msg) message("Estimating growth parameters from length-at-age array in cpars")
      for (ss in 1:nsim) {
        if(Msg) {
          cat(".")
          flush.console()
        }
        pars <- sapply(1:(nyears + proyears), function(X) optim(starts, fitVB, LatAge=Len_age[ss,,X], ages=1:maxage)$par)
        Linfarray[ss,] <- round(pars[1,],2)
        Karray[ss,] <- round(pars[2,],2)
        t0[ss]<- mean(pars[3,])
      }
      Linf <- Linfarray[, nyears]
      K <- Karray[, nyears]
      t0array <- matrix(t0, nrow=nsim, ncol=proyears+nyears)
      if (Msg) cat("\n")
    }
    MaxBin <- ceiling(max(Stock@cpars$Len_age) + 3 * max(Stock@cpars$Len_age) * max(Stock@LenCV)) 
  }
  
  StockOut$maxlen <- maxlen <- Len_age[, maxage, nyears] # reference length for Vmaxlen 
  
  # == Sample CV Length-at-age ====
  if (!exists("LenCV", inherits=FALSE)) LenCV <- runif(nsim, min(Stock@LenCV), max(Stock@LenCV))
  
  # == Generate Catch at Length Classes ====
  if (!exists("LatASD", inherits=FALSE)) LatASD <- Len_age * array(LenCV, dim=dim(Len_age)) # SD of length-at-age 
  if (any(dim(LatASD) != dim(Len_age))) stop("Dimensions of 'LatASD' must match dimensions of 'Len_age'", .call=FALSE)
  
  binWidth <- ceiling(0.03 * MaxBin)
  if (!exists("CAL_bins", inherits=FALSE)) CAL_bins <- seq(from = 0, to = MaxBin + binWidth, by = binWidth)
  if (!exists("CAL_binsmid", inherits=FALSE)) CAL_binsmid <- seq(from = 0.5 * binWidth, by = binWidth, length = length(CAL_bins) - 1)
  if (length(CAL_bins) != length(CAL_binsmid)+1) stop("Length of 'CAL_bins' must be length(CAL_binsmid)+1", .call=FALSE)
  
  nCALbins <- length(CAL_binsmid)
  
  # === Create Weight-at-Age array ====
  if (!exists("Wt_age", inherits=FALSE)) {
    Wt_age <- array(NA, dim = c(nsim, maxage, nyears + proyears))  # Weight at age array
    ind <- as.matrix(expand.grid(1:nsim, 1:maxage, 1:(nyears + proyears)))  # an index for calculating Weight at age 
    Wt_age[ind] <- Stock@a * Len_age[ind]^Stock@b  # Calculation of weight array
    Wa <- Stock@a
    Wb <- Stock@b 
  }	else {
    if (any(dim(Wt_age) != c(nsim, maxage, nyears + proyears))) 
      stop("'Wt_age' must be array with dimensions: nsim, maxage, nyears + proyears (", paste(c(nsim, maxage, nyears + proyears), ""), ") but has ", paste(dim(Wt_age), "")) 
    # Estimate length-weight parameters from the Wt_age data
    logL <- log(as.numeric(Len_age))
    logW <- log(as.numeric(Wt_age))
    mod  <- lm(logW ~ logL)
    EstVar <- summary(mod)$sigma^2
    Wa <- exp(coef(mod)[1]) * exp((EstVar)/2)
    Wb <- coef(mod)[2]
  }
  
  # == Sample Maturity Parameters ====
  if (exists("Mat_age", inherits=FALSE)){
    if (any(dim(Mat_age) != c(nsim, maxage, nyears+proyears))) stop("'Mat_age' must be array with dimensions: nsim, maxage, nyears+proyears") 
    
    # Calculate L50, L95, ageM and age95 
    ageM <- age95 <- L50array <- L95array <- matrix(NA, nsim, nyears+proyears)
    for (XX in 1:(nyears+proyears)) {
      ageM[,XX] <- unlist(sapply(1:nsim, function(x) LinInterp(Mat_age[x,, XX], y=1:maxage, 0.5)))
      age95[,XX] <- unlist(sapply(1:nsim, function(x) LinInterp(Mat_age[x,, XX], y=1:maxage, 0.95)))
      L50array[,XX] <- unlist(sapply(1:nsim, function(x) LinInterp(Mat_age[x,,XX], y=Len_age[x, , nyears], 0.5)))
      L95array[,XX]<- unlist(sapply(1:nsim, function(x) LinInterp(Mat_age[x,,XX], y=Len_age[x, , nyears], 0.95)))
    }
    L50 <- L50array[,nyears]
    L95 <- L95array[,nyears]
    L50[!is.finite(L50)] <- 0.8*Linf[!is.finite(L50)]
    L95[!is.finite(L95)] <- 0.99*Linf[!is.finite(L95)]
    if (any(L50>= Linf)) {
      if (Msg) message("Note: Some samples of L50 are above Linf. Defaulting to 0.8*Linf")
      L50[L50>=Linf] <- 0.8* Linf[L50>=Linf]
    }
    if (any(L95> Linf)) {
      if (Msg)  message("Note: Some samples of L95 are above Linf. Defaulting to 0.99*Linf")
      L95[L95> Linf] <- 0.99* Linf[L95> Linf]
    }
    
    L50_95 <- L95 - L50
  } else {
    if (!exists("L50", inherits=FALSE)) {
      sL50 <- array(runif(nsim * 50, Stock@L50[1], Stock@L50[2]), c(nsim, 50))  # length at 50% maturity  
      # checks for unrealistically high length at maturity
      sL50[sL50/Linf > 0.95] <- NA
      L50 <- apply(sL50, 1, function(x) x[!is.na(x)][1])
      L50[is.na(L50)] <- 0.95 * Linf[is.na(L50)]
    }
    if (!exists("L50_95", inherits=FALSE)) {
      L50_95 <- array(runif(nsim * 50, Stock@L50_95[1], Stock@L50_95[2]), c(nsim, 50))  # length at 95% maturity
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
    if (Msg) message("Note: Some samples of age of maturity are above 'maxage'-1. Defaulting to maxage-1")
    ageM[ageM >= (maxage-1)] <- maxage - 1 
  }
  if (any(age95 >= maxage)) {
    if (Msg) message("Note: Some samples of age of 95 per cent maturity are above 'maxage'. Defaulting to maxage")
    age95[age95 >= maxage] <- maxage  
  }
  
  # == Generate Maturity-at-Age array ====
  if (!exists("Mat_age", inherits=FALSE)) {
    Mat_age <- array(NA, dim=c(nsim, maxage, nyears+proyears))
    
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
    Mage <- matrix(NA, nsim, maxage)
    for (sim in 1:nsim) {
      ind <- findInterval(Len_age[sim,,nyears], M_at_Length[,1])  
      Mage[sim, ] <- MatLen[sim, ind]  
    }
  }
  
  
  # == M-at-age has been provided in OM ====
  if (exists("Mage", inherits=FALSE)) {
    if (exists("M", inherits=FALSE) & length(cpars[["M"]])>0) 
      if (Msg) message("M-at-age has been provided in OM. Overiding M from OM@cpars")
    # M is calculated as mean M of mature ages
    M <- rep(NA, nsim)
    for (sim in 1:nsim) M[sim] <- mean(Mage[sim,round(ageM[sim],0):maxage])
  }
  
  # == Mean Natural mortality by simulation and year ====
  if (exists("M_ageArray", inherits=FALSE)) {
    if (!all(dim(M_ageArray) == c(nsim, maxage, proyears+nyears))) stop("'M_ageArray' must be array with dimensions: nsim, maxage, nyears + proyears") 
    if(Msg) message("M_ageArray has been provided in OM@cpars. Ignoring OM@Mexp, OM@Msd, and OM@Mgrad")
    Mexp <- Msd <- Mgrad <- rep(0, nsim)
  }
   
  
  if (!exists("Marray", inherits=FALSE) & exists("M_ageArray", inherits=FALSE)) {
    Marray <- matrix(NA, nsim, nyears+proyears)
    for (yr in 1:(nyears+proyears)) {
      for (sim in 1:nsim) {
        Marray[sim, yr] <- mean(M_ageArray[sim, ageM[sim,yr]:maxage,yr])
      }
    }
  }
  
  if (!exists("Marray", inherits=FALSE)) {
    Marray <- gettempvar(M, Msd, Mgrad, nyears + proyears, nsim, Mrand)  # M by sim and year according to gradient and inter annual variability
  } else {
    if (any(dim(Marray) != c(nsim, nyears + proyears))) stop("'Marray' must be array with dimensions: nsim, nyears + proyears") 
  }
  
  # == Natural mortality by simulation, age and year ====
  if (!exists("M_ageArray", inherits=FALSE)) { # only calculate M_ageArray if it hasn't been specified in cpars
    
    M_ageArray <- array(NA, dim=c(nsim, maxage, nyears + proyears))
    if (exists("Mage", inherits=FALSE)) { # M-at-age has been provided
      temp1 <- Mage/ matrix(apply(Mage, 1, mean), nsim, maxage, byrow=FALSE)
      ind <- as.matrix(expand.grid(1:nsim, 1:maxage, 1:(nyears+proyears)))
      M_ageArray[ind] <- temp1[ind[,1:2]] * Marray[ind[,c(1,3)]]
    } else { # M-at-age calculated from Lorenzen curve 
      Winf <- Stock@a * Linf^Stock@b
      ind <- as.matrix(expand.grid(1:nsim, 1:maxage, 1:(nyears+proyears)))
      M_ageArray[ind] <- Marray[ind[,c(1,3)]] * (Wt_age[ind]/Winf[ind[,1]]) ^ Mexp[ind[,1]]  
    }  
    
    
    # == Scale M at age so that mean M of mature ages is equal to sampled M ====
    tempM_ageArray <- M_ageArray
    for (sim in 1:nsim) {
      matyrs <- ageM[sim, nyears]:maxage
      if (length(matyrs) >1) {
        # scale <- Marray[sim,]/ apply(tempM_ageArray[sim,ageM[sim]:maxage,], 2, mean) 
        scale <- Marray[sim,]/ (apply(tempM_ageArray[sim,matyrs,], 2, sum)/length(matyrs)) # this is about 4 times faster
      } else if (length(matyrs)==1){
        scale <- Marray[sim,]/ tempM_ageArray[sim,ageM[sim]:maxage,]  
      } 
      
      M_ageArray[sim,,] <- M_ageArray[sim,,] * matrix(scale, maxage, nyears+proyears, byrow=TRUE)
    }
    
  }
  
  # == Sample Discard Mortality ====
  if(!exists("Fdisc", inherits = FALSE)) Fdisc <- runif(nsim, min(Stock@Fdisc), max(Stock@Fdisc))
  StockOut$Fdisc <- Fdisc 
  
  # == 
   
  

  # Check if M-at-age is constant that Maxage makes sense
  if (all(M_ageArray[1,,1] == mean(M_ageArray[1,,1]))) { # constant M at age
    calcMax <- ceiling(-log(0.01)/(min(M)))        # Age at which 1% of cohort survives
    if (maxage < 0.8*calcMax && Msg) {
      message("Note: Maximum age (", maxage, ") is lower than assuming 1% of cohort survives to maximum age (", calcMax, ")")
    }  
  }
  
  # --- Calculate movement ----
  initdist <- Pinitdist <- NULL
  if (!exists("mov", inherits=FALSE)) {
    if(Msg) message("Optimizing for user-specified movement")  # Print a progress update
    nareas<-2 # default is a 2 area model
    mov1 <- array(t(sapply(1:nsim, getmov2, Frac_area_1 = Frac_area_1, 
                           Prob_staying = Prob_staying)), dim = c(nsim, nareas, nareas))
    mov<-array(NA,c(nsim,maxage,nareas,nareas))
    mind<-as.matrix(expand.grid(1:nsim,1:maxage,1:nareas,1:nareas))
    mov[mind]<-mov1[mind[,c(1,3,4)]]
    
    initdist <- array(0,c(nsim,maxage,nareas))
    initdist[,,1]<-Frac_area_1
    initdist[,,2]<- 1- Frac_area_1  
    
  }else{ # if mov is specified need to calculate age-based spatial distribution (Pinitdist to initdist)
    nareas<-dim(mov)[3]
    if(Msg) message(paste("Custom movement matrix detected, simulating movement among",nareas,"areas"))
    
    mind<-as.matrix(expand.grid(1:nsim,maxage,1:nareas,1:nareas))
    movedarray<-array(0,c(nsim,nareas,nareas))
    Pinitdist<-array(1/nareas,c(nsim,nareas))
    for(i in 1:20){ # convergence in initial distribution is assumed to occur in 20 iterations (generally overkill)
      movedarray[mind[,c(1,3,4)]]<-Pinitdist[mind[,c(1,3)]]*mov[mind] # distribution in from areas mulitplied by movement array
      Pinitdist<-apply(movedarray,c(1,3),sum) # add over to areas
      #print(initdist[1:2,]) # debugging to check convergence
    }
    
  }
  
  if (dim(Asize)[2]!=nareas) {
    if(Msg) message('Asize is not length "nareas", assuming all areas equal size')
    Asize <- matrix(1/nareas, nrow=nsim, ncol=nareas)
  }
  
  
  
  StockOut$Mexp <- Mexp 
  StockOut$Msd <- Msd 
  StockOut$Mgrad <- Mgrad
  
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
  StockOut$Linfgrad <- Linfgrad
  # StockOut$recgrad <- recgrad
  StockOut$K <- K
  StockOut$Ksd <- Ksd
  StockOut$Kgrad <- Kgrad
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
#'
#' @return A named list of sampled Fleet parameters
#' @export
#'
SampleFleetPars <- function(Fleet, Stock=NULL, nsim=NULL, nyears=NULL, proyears=NULL, cpars=NULL) {
  if (class(Fleet) != "Fleet" & class(Fleet) != "OM") 
    stop("First argument must be class 'Fleet' or 'OM'")
  
  if (class(Fleet) != "OM" & class(Stock) != "Stock" & class(Stock) != "list") 
    stop("Must provide 'Stock' object", call.=FALSE)
  
  # Get custom pars if they exist
  if (class(Fleet) == "OM" && length(Fleet@cpars) > 0 && is.null(cpars)) 
    cpars <- SampleCpars(Fleet@cpars, Fleet@nsim)  # custom parameters exist in Stock/OM object
  if (length(cpars) > 0) { # custom pars exist - assign to function environment 
    Names <- names(cpars)
    for (X in 1:length(Names)) assign(names(cpars)[X], cpars[[X]])
  }
  
  if (class(Fleet) == "OM") {
    nsim <- Fleet@nsim
    nyears <- Fleet@nyears 
    proyears <- Fleet@proyears
    StockPars <- SampleStockPars(Fleet, nsim, nyears, proyears, cpars)
    for (X in 1:length(StockPars)) assign(names(StockPars)[X], StockPars[[X]])
  }
  
  Fleet <- updateMSE(Fleet) # update to add missing slots with default values
  if (class(Stock) == "Stock") {
    Stock <- updateMSE(Stock) # update to add missing slots with default values
    # Sample Stock Pars - need some to calculate selectivity at age and length  
    StockPars <- SampleStockPars(Stock, nsim, nyears, proyears, cpars)
    for (X in 1:length(StockPars)) assign(names(StockPars)[X], StockPars[[X]])
  } 
  if (class(Stock) == "list") for (X in 1:length(Stock)) assign(names(Stock)[X], Stock[[X]])
  
  Fleetout <- list()
  
  # == Sample Historical Fishing Effort =====
  if (!exists("Esd", inherits = FALSE)) Esd <- runif(nsim, Fleet@Esd[1], Fleet@Esd[2])  # interannual variability in fishing effort (log normal sd)
  if (!exists("EffLower", inherits = FALSE)) EffLower <- Fleet@EffLower
  if (!exists("EffUpper", inherits = FALSE)) EffUpper <- Fleet@EffUpper 
  if (!exists("EffYears", inherits = FALSE)) EffYears <- Fleet@EffYears
  
  # if (max(EffYears) != nyears && max(EffYears) != 1) stop("Maximum EffYears (", max(EffYears), ") not equal to nyears (", nyears, ")")
  
  if (!exists("Find", inherits = FALSE)) {
    if (any(is.na(EffLower)) || any(is.na(EffUpper)) || any(is.na(EffYears))) {
      message("NAs in EffLower, EffUpper, or EffYears")
      Find <- matrix(NA, nsim, nyears)
      Deriv <- list(Find, rep(NA, nsim))
    } else {
      Deriv <- getEffhist(Esd, nyears, EffYears = EffYears, EffLower = EffLower, EffUpper = EffUpper)  # Historical fishing effort  
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
  if (any(dim(Find) != c(nsim, nyears))) stop("Find must be matrix with dimensions: nsim (", nsim, "), nyears (", nyears, ") but is: ", paste(dim(Find), ""))
  
  Fleetout$Esd <- Esd
  Fleetout$Find <- Find
  Fleetout$dFfinal <- dFfinal
  
  # === Spatial Targetting ====
  if (!exists("Spat_targ", inherits = FALSE))  Spat_targ <- runif(nsim, Fleet@Spat_targ[1], Fleet@Spat_targ[2])  # spatial targetting Ba^targetting param 
  
  Fleetout$Spat_targ <- Spat_targ
  
  # === Sample fishing efficiency parameters ====
  if (!exists("qinc", inherits = FALSE)) qinc <- runif(nsim, Fleet@qinc[1], Fleet@qinc[2])
  if (!exists("qcv", inherits = FALSE)) qcv <- runif(nsim, Fleet@qcv[1], Fleet@qcv[2])  # interannual variability in catchability
  
  # === Simulate future variability in fishing efficiency ====
  qmu <- -0.5 * qcv^2  # Mean
  if (!exists("qvar", inherits = FALSE)) qvar <- array(exp(rnorm(proyears * nsim, rep(qmu, proyears), rep(qcv, proyears))), c(nsim, proyears))  # Variations in interannual variation
  FinF <- Find[, nyears]  # Effort in final historical year
  
  Fleetout$qinc <- qinc
  Fleetout$qcv <- qcv
  Fleetout$qvar <- qvar
  Fleetout$FinF <- FinF
  
  # ==== Sample selectivity parameters ====
  if (exists("V", inherits=FALSE) | exists("retA", inherits=FALSE)) {
    Fleet@isRel <- 'FALSE'
  }
  
  
  Selnyears <- length(Fleet@SelYears)
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
 
  if (exists("L5", inherits = FALSE) | exists("LFS", inherits = FALSE) | 
      exists("Vmaxlen", inherits = FALSE) | exists("V", inherits=FALSE)) {
    if (all(multi != 1)) stop("Selectivity parameters provided in cpars must be absolute values. Is Fleet@isRel == 'FALSE'?")
  }
  
  if (!exists("L5", inherits = FALSE)) L5 <- runif(nsim, Fleet@L5[1], Fleet@L5[2]) * multi  # length at 0.05% selectivity ascending
  if (!exists("LFS", inherits = FALSE)) LFS <- runif(nsim, Fleet@LFS[1], Fleet@LFS[2]) * multi  # first length at 100% selection
  if (!exists("Vmaxlen", inherits = FALSE)) Vmaxlen <- runif(nsim, Fleet@Vmaxlen[1], Fleet@Vmaxlen[2])  # selectivity at maximum length
  
  Vmaxlen[Vmaxlen<=0] <- tiny
  L5s <- LFSs <- Vmaxlens <- NULL  # initialize 
  
  if (Selnyears > 1) {   # change of selectivity in historical years 
    # length at 0.05% selectivity ascending
    L5s <- mapply(runif, n = nsim, min = Fleet@L5Lower, max = Fleet@L5Upper) * multi
    # first length at 100% selection
    LFSs <- mapply(runif, n = nsim, min = Fleet@LFSLower, max = Fleet@LFSUpper) *  multi
    # selectivity at maximum length
    Vmaxlens <- mapply(runif, n = nsim, min = Fleet@VmaxLower, max = Fleet@VmaxUpper)
  } else {
    L5s <- LFSs <- Vmaxlens <- NA
  }
  Fleetout$L5s <- L5s
  Fleetout$LFSs <- LFSs
  Fleetout$Vmaxlens <- Vmaxlens
  
  # == Calculate Selectivity at Length ====
  nCALbins <- length(CAL_binsmid)
  SLarray <- array(NA, dim=c(nsim, nCALbins, nyears+proyears)) # Selectivity-at-length 
  CAL_binsmidMat <- matrix(CAL_binsmid, nrow=nsim, ncol=length(CAL_binsmid), byrow=TRUE)
  if (exists("V", inherits=FALSE)) { # V has been passed in with custompars
    if(dim(V)[3] != proyears + nyears) V<-abind::abind(V,array(V[,,nyears],c(nsim,maxage,proyears)),along=3) # extend future Vulnerabiliy according to final historical vulnerability
    # assign L5, LFS and Vmaxlen - dodgy loop 
    # could calculate length at 5% selectivity from vB
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
        # SLarray[s,, yr] <- SelectFun(s, SL0.05=L5[yr, ], SL1=LFS[yr, ], MaxSel=Vmaxlen[yr, ], 
        # maxlens=Len_age[, maxage, nyears], Lens=CAL_binsm
        
      }
      srs <- (Linf - LFS[yr,]) / ((-log(Vmaxlen[yr,drop=FALSE],2))^0.5) # selectivity parameters are constant for all years
      sls <- (LFS[yr,] - L5[yr, ]) /((-log(0.05,2))^0.5)
      SLarray[,, yr] <- t(sapply(1:nsim, getsel, lens=CAL_binsmidMat, lfs=LFS[yr, ], sls=sls, srs=srs))
      
    }
    
  }
  
  # == Calculate Selectivity at Age and Length ====
  CAL_binsmidMat <- matrix(CAL_binsmid, nrow=nsim, ncol=length(CAL_binsmid), byrow=TRUE)
  if (!exists("V", inherits=FALSE)) { # don't run if V has been passed in with custompars 
    if (Selnyears <= 1) {    
      L5 <- matrix(L5, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
      LFS <- matrix(LFS, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
      Vmaxlen <- matrix(Vmaxlen, nrow = nyears + proyears, ncol = nsim, byrow = TRUE) 
      
      # ind <- which(LFS/matrix(Linf, nrow = proyears + nyears, ncol = nsim, byrow = TRUE) > 1, arr.ind = T)
      # if (length(ind) > 0) {
      # message("LFS too high (LFS > Linf) in some cases. \nDefaulting to LFS = 0.9 Linf for the affected simulations")
      # LFS[ind] <- Linf[ind[, 2]] * 0.9
      # } 
      
      
      # Calculate selectivity-at-age  curve 
      V <- array(NA, dim = c(nsim, maxage, nyears + proyears)) 
      # s1 <- sapply(1:nsim, function(i) optimize(getSlope1, interval = c(0, 1e+05), 
      #                                           LFS = LFS[1, i], L0.05 = L5[1,i])$minimum)	
      # if (all(Vmaxlen >= 0.99)) s2 <- rep(1E5, nsim)
      # Vmaxlen[Vmaxlen ==0] <- 0.001 # fix for when Vmaxlen == 0 
      # if (!all(Vmaxlen >= 0.99)) 
      #   s2 <- sapply(1:nsim, function(i) optimize(getSlope2, interval = c(0, 1e+05), 
      #                                             LFS = LFS[1,i], s1=s1[i], maxlen=maxlen[i], 
      #                                             MaxSel=Vmaxlen[1, i])$minimum)
      
      srs <- (Linf - LFS[1,]) / ((-log(Vmaxlen[1,],2))^0.5) # selectivity parameters are constant for all years
      sls <- (LFS[1,] - L5[1, ]) /((-log(0.05,2))^0.5)
      
      
      # Calculate selectivity at length class 
      
      if (nsim>1) SelLength <- t(sapply(1:nsim, getsel, lens=CAL_binsmidMat, lfs=LFS[1, ], sls=sls, srs=srs))
      if (nsim == 1) SelLength <- getsel(1, lens=CAL_binsmidMat, lfs=LFS[1, ], sls=sls, srs=srs)
      
      for (yr in 1:(nyears+proyears)) {
        
        # Calculate selectivity at age class 
        # V[ , , yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[1,i], s1[i], s2[i], lens=Len_age[i,,yr])))
        if(nsim>1) V[ , , yr] <- t(sapply(1:nsim, getsel, lens=Len_age[,,yr], lfs=LFS[1,], sls=sls, srs=srs))
        
        if(nsim == 1) V[ , , yr] <- getsel(x=1, lens=t(matrix(Len_age[,,yr])), lfs=LFS[1,], sls=sls, srs=srs)
        # Calculate selectivity at length class 
        # SLarray[,, yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[1,i], s1[i], s2[i], lens=CAL_binsmid)))
        SLarray[,, yr] <- SelLength
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
        if (nsim>1) {
          LFS[bkyears, ] <- matrix(rep((LFSs[, X]), length(bkyears)), ncol = nsim, byrow = TRUE)
          Vmaxlen[bkyears, ] <- matrix(rep((Vmaxlens[, X]), length(bkyears)), ncol = nsim, byrow = TRUE)
          L5[bkyears, ] <- matrix(rep((L5s[, X]), length(bkyears)), ncol = nsim, byrow = TRUE)
        } else {
          LFS[bkyears, ] <- matrix(rep((LFSs[X]), length(bkyears)), ncol = nsim, byrow = TRUE)
          Vmaxlen[bkyears, ] <- matrix(rep((Vmaxlens[X]), length(bkyears)), ncol = nsim, byrow = TRUE)
          L5[bkyears, ] <- matrix(rep((L5s[X]), length(bkyears)), ncol = nsim, byrow = TRUE)
          
        }
        
        srs <- (Linf - LFS[bkyears[1],]) / ((-log(Vmaxlen[bkyears[1],],2))^0.5) #
        sls <- (LFS[bkyears[1],] - L5[bkyears[1], ]) /((-log(0.05,2))^0.5)
        
        # Calculate selectivity at length class 
        SelLength <- t(sapply(1:nsim, getsel, lens=CAL_binsmidMat, lfs=LFS[bkyears[1],], sls=sls, srs=srs))
        
        
        # s1 <- sapply(1:nsim, function(i) optimize(getSlope1, interval = c(0, 1e+05), 
        #                                           LFS = LFSs[i, X], L0.05 = L5s[i, X])$minimum)
        # s2 <- sapply(1:nsim, function(i) optimize(getSlope2, interval = c(0, 1e+05), 
        #                                           LFS = LFSs[i, X], s1=s1[i], maxlen=maxlen[i], 
        #                                           MaxSel=Vmaxlens[i, X])$minimum)	
        for (yr in bkyears) {
          # V[ , , yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[yr, i], s1[i], s2[i], lens=Len_age[i,,yr])))
          # SLarray[,, yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[1,i], s1[i], s2[i], lens=CAL_binsmid)))   
          V[ , , yr] <-  t(sapply(1:nsim, getsel, lens=Len_age[,,yr], lfs=LFS[yr,], sls=sls, srs=srs))
          SLarray[,, yr] <- SelLength 
        }
      }
      
      restYears <- max(SelYears):(nyears + proyears)
      if (nsim>1) {
        L5[restYears, ] <- matrix(rep((L5s[, Selnyears]), length(restYears)), ncol = nsim, byrow = TRUE)
        LFS[restYears, ] <- matrix(rep((LFSs[, Selnyears]), length(restYears)), ncol = nsim, byrow = TRUE)
        Vmaxlen[restYears, ] <- matrix(rep((Vmaxlens[, Selnyears]), length(restYears)), ncol = nsim, byrow = TRUE)
      } else {
        L5[restYears, ] <- matrix(rep((L5s[Selnyears]), length(restYears)), ncol = nsim, byrow = TRUE)
        LFS[restYears, ] <- matrix(rep((LFSs[Selnyears]), length(restYears)), ncol = nsim, byrow = TRUE)
        Vmaxlen[restYears, ] <- matrix(rep((Vmaxlens[Selnyears]), length(restYears)), ncol = nsim, byrow = TRUE)
      }
      
      
      # s1 <- sapply(1:nsim, function(i) optimize(getSlope1, interval = c(0, 1e+05), 
      #                                           LFS = LFSs[i, Selnyears], L0.05 = L5s[i, Selnyears])$minimum)
      # s2 <- sapply(1:nsim, function(i) optimize(getSlope2, interval = c(0, 1e+05), 
      #                                           LFS = LFSs[i, Selnyears], s1=s1[i], maxlen=maxlen[i], 
      #                                           MaxSel=Vmaxlens[i, Selnyears])$minimum)	
      
      srs <- (Linf - LFS[restYears[1],]) / ((-log(Vmaxlen[restYears[1],],2))^0.5) #
      
      sls <- (LFS[restYears[1],] - L5[restYears[1], ]) /((-log(0.05,2))^0.5)
      
      # Calculate selectivity at length class 
      SelLength <- t(sapply(1:nsim, getsel, lens=CAL_binsmidMat, lfs=LFS[restYears[1],], sls=sls, srs=srs))
      
      for (yr in restYears) { 
        # V[ , , restYears] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[yr, i], s1[i], s2[i], lens=Len_age[i,,yr])))		
        # SLarray[,, yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[1,i], s1[i], s2[i], lens=CAL_binsmid))) 
        V[ , , yr] <- t(sapply(1:nsim, getsel, lens=Len_age[,,yr], lfs=LFS[yr,], sls=sls, srs=srs))
        SLarray[,, yr] <- SelLength
        
      }	 
    }
  } # end of 'if V exists'
  
  # Check LFS is greater than L5 
  chk <- sum(apply(L5 > LFS, 2, prod) != 0)
  if (chk > 0) stop("L5 is greater than LFS in ", chk, ' simulations')
  
  
  if (any((dim(V) != c(nsim, maxage, proyears+nyears)))) 
    stop("V must have dimensions: nsim (", nsim,") maxage (", maxage, 
         ") proyears+nyears (", proyears+nyears, ") \nbut has ", 
         dim(V)[1], " ", dim(V)[2], " ", dim(V)[3], call.=FALSE)
  
  
  # == Sample Retention Parameters ====
  if(!exists("LR5", inherits = FALSE)) LR5 <- runif(nsim, min(Fleet@LR5), max(Fleet@LR5)) * multi
  if(!exists("LFR", inherits = FALSE)) LFR <- runif(nsim, min(Fleet@LFR), max(Fleet@LFR)) * multi
  if(!exists("Rmaxlen", inherits = FALSE)) Rmaxlen <- runif(nsim, min(Fleet@Rmaxlen), max(Fleet@Rmaxlen))
  if(!exists("DR", inherits = FALSE)) DR <- runif(nsim, min(Fleet@DR), max(Fleet@DR))
  
  if (any(LR5 > LFR)) stop('LR5 is greater than LFR', call.=FALSE)
  # == Calculate Retention Curve ====
  Rmaxlen[Rmaxlen<=0] <- tiny 
  LR5 <- matrix(LR5, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
  LFR <- matrix(LFR, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
  Rmaxlen <- matrix(Rmaxlen, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
  DR <- matrix(DR, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
  
  
  # s1 <- sapply(1:nsim, function(i) optimize(getSlope1, interval = c(0, 1e+05), 
  #                                           LFS = LFR[1, i], L0.05 = LR5[1,i])$minimum)
  # if (all(Rmaxlen >= 0.99)) s2 <- rep(1E5, nsim)
  # if (!all(Rmaxlen >= 0.99)) 
  #   s2 <- sapply(1:nsim, function(i) optimize(getSlope2, interval = c(0, 1e+05), 
  #                                             LFS = LFR[1,i], s1=s1[i], maxlen=maxlen[i], 
  #                                             MaxSel=Rmaxlen[1, i])$minimum)
  
  srs <- (Linf - LFR[1,]) / ((-log(Rmaxlen[1,],2))^0.5) # selectivity parameters are constant for all years
  sls <- (LFR[1,] - LR5[1,]) /((-log(0.05,2))^0.5)
  
  RetLength <- t(sapply(1:nsim, getsel, lens=CAL_binsmidMat, lfs=LFR[1,], sls=sls, srs=srs))
  
  if (!exists("retA", inherits=FALSE)) {
    retA <- array(NA, dim = c(nsim, maxage, nyears + proyears)) # retention at age
    for (yr in 1:(nyears+proyears)) {
      # Calculate retention at age class 
      # retA[ , , yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFR[yr,i], s1[i], s2[i], lens=Len_age[i,,yr])))
      if (nsim>1) retA[ , , yr] <- t(sapply(1:nsim, getsel, lens=Len_age[,,yr], lfs=LFR[1,], sls=sls, srs=srs))
      if (nsim == 1) retA[ , , yr] <- getsel(1, lens=t(matrix(Len_age[,,yr])), lfs=LFR[1,], sls=sls, srs=srs)
      
    } 
  } else {
    # check dimensions 
    if (any((dim(retA) != c(nsim, maxage, proyears+nyears)))) 
      stop("retA must have dimensions: nsim (", nsim,") maxage (", maxage, 
           ") proyears+nyears (", proyears+nyears, ") \nbut has ", 
           dim(retA)[1], " ", dim(retA)[2], " ", dim(retA)[3], call.=FALSE) 
  }
  
  # 
  if (!exists("retL", inherits=FALSE)) {
    retL <- array(NA, dim = c(nsim, nCALbins, nyears + proyears)) # retention at length
    for (yr in 1:(nyears+proyears)) {
      # Calculate retention at length class 
      # retL[,, yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFR[yr,i], s1[i], s2[i], lens=CAL_binsmid)))   
      retL[,, yr] <- RetLength
    }
  } else {
    # check dimensions 
    if (any((dim(retL) != c(nsim, nCALbins, proyears+nyears)))) 
      stop("retL must have dimensions: nsim (", nsim,") nCALbins (", nCALbins, 
           ") proyears+nyears (", proyears+nyears, ") \nbut has ", 
           dim(retL)[1], " ", dim(retL)[2], " ", dim(retL)[3], call.=FALSE) 
  }
  
  
  V2 <- V
  SLarray2 <- SLarray
  
  # Apply general discard rate 
  dr <- aperm(abind::abind(rep(list(DR), maxage), along=3), c(2,3,1))
  retA <- (1-dr) * retA
  
  dr <- aperm(abind::abind(rep(list(DR), nCALbins), along=3), c(2,3,1))
  retL <- (1-dr) * retL
  
  # update realized vulnerablity curve with retention and dead discarded fish 
  Fdisc_array1 <- array(Fdisc, dim=c(nsim, maxage, nyears+proyears))
  V <- V * (retA + (1-retA)*Fdisc_array1) # Realised selection at age
  
  Fdisc_array2 <- array(Fdisc, dim=c(nsim, nCALbins, nyears+proyears))
  SLarray <- SLarray2 * (retL + (1-retL)*Fdisc_array2) # Realised selection at length
  
  # Realised Retention curves
  retA <- retA * V2
  retL <- retL * SLarray2
  
  Fleetout$Fdisc <- Fdisc
  Fleetout$Fdisc_array1 <- Fdisc_array1
  Fleetout$Fdisc_array2 <- Fdisc_array2
  Fleetout$LR5 <- LR5  
  Fleetout$LFR <- LFR 
  Fleetout$Rmaxlen <- Rmaxlen
  Fleetout$DR <- DR
  
  Fleetout$retA <- retA  # retention-at-age array - nsim, maxage, nyears+proyears
  Fleetout$retL <- retL  # retention-at-length array - nsim, nCALbins, nyears+proyears
  
  Fleetout$L5 <- L5  
  Fleetout$LFS <- LFS 
  Fleetout$Vmaxlen <- Vmaxlen 
  Fleetout$V <- V  # realized vulnerability-at-age
  Fleetout$SLarray <- SLarray # realized vulnerability-at-length
  Fleetout$V2 <- V2 # original vulnerablity-at-age curve 
  Fleetout$SLarray2 <- SLarray2 # original vulnerablity-at-length curve 
  
  Fleetout 
}

#' Sample Observation Parameters
#'
#' @param Obs An object of class 'Obs' or class 'OM'
#' @param nsim Number of simulations. Ignored if 'Obs' is class 'OM'
#' @param cpars Optional named list of custom parameters. Ignored if 'OM' is class 'OM'
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
  
  # fix some naming issues?
  
  if (!exists("Csd", inherits=FALSE)) {
    ObsOut$Csd <- runif(nsim, Obs@Cobs[1], Obs@Cobs[2])  # Sampled catch observation error (lognormal sd)
  } else {
    ObsOut$Csd <- Csd
  }
  if (!exists("Cbias", inherits=FALSE)) {
    ObsOut$Cbias <- rlnorm(nsim, mconv(1, Obs@Cbiascv), sdconv(1, Obs@Cbiascv))  # Sampled catch bias (log normal sd)
  } else {
    ObsOut$Cbias <- Cbias
  }
  if (!exists("CAA_nsamp", inherits=FALSE)) {
    ObsOut$CAA_nsamp <- ceiling(runif(nsim, Obs@CAA_nsamp[1], Obs@CAA_nsamp[2]))  # Number of catch-at-age observations
  } else {
    ObsOut$CAA_nsamp <- CAA_nsamp
  } 
  if (!exists("CAA_ESS", inherits=FALSE)) {
    ObsOut$CAA_ESS <- ceiling(runif(nsim, Obs@CAA_ESS[1], Obs@CAA_ESS[2]))  # Effective sample size
  } else {
    ObsOut$CAA_ESS <- CAA_ESS
  } 
  if (!exists("CAL_nsamp", inherits=FALSE)) {
    ObsOut$CAL_nsamp <- ceiling(runif(nsim, Obs@CAL_nsamp[1], Obs@CAL_nsamp[2]))  # Observation error standard deviation for single catch at age by area
  } else {
    ObsOut$CAL_nsamp <- CAL_nsamp
  }  
  if (!exists("CAL_ESS", inherits=FALSE)) {
    ObsOut$CAL_ESS <- ceiling(runif(nsim, Obs@CAL_ESS[1], Obs@CAL_ESS[2]))  # Effective sample size
  } else {
    ObsOut$CAL_ESS <- CAL_ESS
  }  
  if (!exists("betas", inherits=FALSE)) {
    ObsOut$betas <- exp(runif(nsim, log(Obs@beta[1]), log(Obs@beta[2])))  # the sampled hyperstability / hyperdepletion parameter beta>1 (hyperdepletion) beta<1 (hyperstability)
  } else {
    ObsOut$betas <- betas
  }  
  if (!exists("Isd", inherits=FALSE)) {
    ObsOut$Isd <- runif(nsim, Obs@Iobs[1], Obs@Iobs[2])  # Abundance index observation error (log normal sd)
  } else {
    ObsOut$Isd <- Isd
  } 
  if (!exists("Derr", inherits=FALSE)) {
    ObsOut$Derr <- runif(nsim, Obs@Dobs[1], Obs@Dobs[2])
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
  if (!exists("Aerr", inherits=FALSE)) {
    ObsOut$Aerr <- runif(nsim, Obs@Btobs[1], Obs@Btobs[2])
  } else {
    ObsOut$Aerr <- Aerr
  }
  if (!exists("Abias", inherits=FALSE)) {
    ObsOut$Abias <- exp(runif(nsim, log(Obs@Btbiascv[1]), log(Obs@Btbiascv[2])))  #rlnorm(nsim,mconv(1,Obs@Btbiascv),sdconv(1,Obs@Btbiascv))    # sample of current abundance bias
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
    ObsOut$Recsd <- runif(nsim, Obs@Recbiascv[1], Obs@Recbiascv[2])  # Recruitment deviation  
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
#' @param nsim Number of simulations. Ignored if 'Stock' is class 'OM'
#' @param cpars Optional named list of custom parameters. Ignored if 'OM' is class 'OM'
#' @return A named list of sampled Implementation Error parameters
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
    ImpOut$TACSD <- runif(nsim, Imp@TACSD[1], Imp@TACSD[2])  # Sampled TAC error (lognormal sd)
  } else {
    ImpOut$TACSD <- TACSD
  }
  if (!exists("TACFrac", inherits = FALSE)) {
    ImpOut$TACFrac <- runif(nsim, Imp@TACFrac[1], Imp@TACFrac[2])  # Sampled TAC fraction (log normal sd)
  } else {
    ImpOut$TACFrac <- TACFrac
  }
  if (!exists("TAESD", inherits = FALSE)) {
    ImpOut$TAESD <- runif(nsim, Imp@TAESD[1], Imp@TAESD[2])  # Sampled Effort error (lognormal sd)
  } else {
    ImpOut$TAESD <- TAESD
  }
  if (!exists("TAEFrac", inherits = FALSE)) {
    ImpOut$TAEFrac <- runif(nsim, Imp@TAEFrac[1], Imp@TAEFrac[2])  # Sampled Effort fraction (log normal sd)
  } else {
    ImpOut$TAEFrac <- TAEFrac
  }
  if (!exists("SizeLimSD", inherits = FALSE)) {
    ImpOut$SizeLimSD<-runif(nsim,Imp@SizeLimSD[1],Imp@SizeLimSD[2])
  } else {
    ImpOut$SizeLimSD <- SizeLimSD
  }
  if (!exists("SizeLimFrac", inherits = FALSE)) {
    ImpOut$SizeLimFrac<-runif(nsim,Imp@SizeLimFrac[1],Imp@SizeLimFrac[2])
  } else {
    ImpOut$SizeLimFrac <- SizeLimFrac
  }
  
  
  ImpOut
}


#' Valid custom parameters (cpars)
#'
#' @param print Print the valid names for cpars?
#'
#' @return invisibly returns vector of valid cpars names
#' @export
#'
validcpars <- function(print=TRUE) {
  vnames <- sort(c("D","Esd","Find","procsd","AC","M","Msd", 
                   "Mgrad","hs","Linf","Linfsd","Linfgrad",
                   "K","Ksd","Kgrad","t0","L50", "L95", "L50_95","Spat_targ",
                   "Frac_area_1","Prob_staying","Size_area_1","mov","initdist", "Asize",
                   "Csd","Cbias","CAA_nsamp","CAA_ESS","CAL_nsamp",
                   "CAL_ESS","betas","Isd","Derr","Dbias", 
                   "Mbias","FMSY_Mbias","lenMbias","LFCbias",
                   "LFSbias","Aerr","Abias","Kbias","t0bias", 
                   "Linfbias","Irefbias","Crefbias","Brefbias",
                   "Recsd","qinc","qcv","L5","LFS","Vmaxlen","Perr","R0","Mat_age", 
                   "Mrand","Linfrand","Krand","maxage","V",  
                   "ageM", "age95", "EffYears", "EffLower", "EffUpper",
                   "Wt_age", "Len_age", "Marray", "M_at_Length", "LenCV", 
                   "CAL_binsmid", "CAL_bins", "LatASD", "dFfinal",
                   "LR5", "LFR", "Rmaxlen", "DR", "Fdisc","M_ageArray",
                   "Linfarray", "Karray")) 
  
  if (print) {
    n <- length(vnames)
    vec <- 3:7
    nc <- vec[which.min(n  %% vec)]
    
    options(warn=-1)
    temp <- matrix(vnames, ncol=nc, byrow=TRUE)
    options(warn=1)
    temp[duplicated(as.vector(temp))] <- ""
    print(temp)
  }
  invisible(vnames)
}





#' Sample custom pars
#'
#' @param cpars A named list containing custom parameters for the OM
#' @param nsim number of simulations
#' @param msg logical - print the names of the cpars? Turn off when using the function in a loop
#' @return A named list of sampled custom parameters
#' @export
#'
SampleCpars <- function(cpars, nsim=48, msg=TRUE) {
  
  # Vector of valid names for custompars list or data.frame. Names not in this list will be printed out in warning and ignored #	
  ParsNames <- validcpars(FALSE)
  
  sampCpars <- list()
  ncparsim<-cparscheck(cpars)
  Names <- names(cpars)
  # report invalid names 
  invalid <- which(!Names %in% ParsNames)
  if (length(invalid) > 0) {
    outNames <- paste(Names[invalid], "")
    for (i in seq(5, by=5, length.out=floor(length(outNames)/5))) outNames <- gsub(outNames[i], paste0(outNames[i], "\n"), outNames)
    if(msg) message("ignoring invalid names found in custom parameters (OM@cpars) \n", outNames)	
  }
  # report found names
  valid <- which(Names %in% ParsNames)
  cpars <- cpars[valid]
  if (length(valid) == 0) stop("No valid names found in custompars (OM@cpars)", call.=FALSE)
  Names <- names(cpars)
  outNames <- paste(Names, "")
  for (i in seq(5, by=5, length.out=floor(length(outNames)/5)))
    outNames <- gsub(outNames[i], paste0(outNames[i], "\n"), outNames)
  if(msg) message("valid custom parameters (OM@cpars) found: \n", outNames)
  
  # Sample custom pars 
  if (ncparsim < nsim) ind <- sample(1:ncparsim, nsim, replace=TRUE)
  if (!ncparsim < nsim) ind <- sample(1:ncparsim, nsim, replace=FALSE)
  
  for (i in 1:length(cpars)) {
    samps <- cpars[[i]]
    name <- names(cpars)[i]
    if (any(c("EffUpper", "EffLower", "EffYears", "maxage", "M_at_Length", "CAL_binsmid", "CAL_bins") %in% name)) {
      sampCpars[[name]] <- samps
    } else {
      if (class(samps) == "numeric" | class(samps) == "integer") sampCpars[[name]] <- samps[ind]
      
      if (class(samps) == "matrix") sampCpars[[name]] <- samps[ind,, drop=FALSE] 
      
      if (class(samps) == "array") {
        if (length(dim(samps)) == 3)  sampCpars[[name]] <- samps[ind, , ,drop=FALSE]
        if (length(dim(samps)) == 4)  sampCpars[[name]] <- samps[ind, , , ,drop=FALSE]
      }
      if (class(samps) == "data.frame")   sampCpars[[name]] <- samps 
    }
  }
  
  
  
  sampCpars
}
