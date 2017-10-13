
#' Sample Stock parameters
#'
#' @param Stock An object of class 'Stock' or class 'OM'
#' @param nsim Number of simulations. Ignored if 'Stock' is class 'OM'
#' @param nyears Number of historical years. Ignored if 'Stock' is class 'OM'
#' @param proyears Number of projection years. Ignored if 'Stock' is class 'OM'
#' @param cpars Optional named list of custom parameters. Ignored if 'Stock' is class 'OM'
#'
#' @return A named list of sampled Stock parameters
#' @export
#'   
SampleStockPars <- function(Stock, nsim=48, nyears=80, proyears=50, cpars=NULL) {
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
  if (length(R0) != nsim) R0 <- rep(R0, nsim*50)[1:nsim] # modified to allow for different R0 per sim 
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
  StockOut$Mexp <- Mexp 
  StockOut$Msd <- Msd 
  StockOut$Mgrad <- Mgrad
  
  # == Depletion ====
  if (!exists("dep", inherits=FALSE)) {
    StockOut$dep <- dep <- runif(nsim, Stock@D[1], Stock@D[2])  # sample from the range of user-specified depletion (Bcurrent/B0)  
  } else {
    StockOut$dep <- dep 
  }
  
 
  # == Stock-Recruitment Relationship ====
  if (!exists("SRrel", inherits=FALSE)) {
    StockOut$SRrel <- rep(Stock@SRrel, nsim)  # type of Stock-recruit relationship. 1=Beverton Holt, 2=Ricker
  } else {
    StockOut$SRrel <- SRrel 
  }
  
  if (!exists("hs", inherits=FALSE)) {
    StockOut$hs <- hs <- runif(nsim, Stock@h[1], Stock@h[2])  # sample of recruitment cStockpensation (steepness - fraction of unfished recruitment at 20% of unfished biStockass)
  } else {
    StockOut$hs <- hs
  }
  
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
  if (!exists("Linf", inherits=FALSE)) Linf <- runif(nsim, Stock@Linf[1], Stock@Linf[2])  # sample of asymptotic length
  if (!exists("Linfsd", inherits=FALSE)) Linfsd <- runif(nsim, Stock@Linfsd[1], Stock@Linfsd[2])  # sample of interannual variability in Linf
  if (!exists("Linfgrad", inherits=FALSE)) Linfgrad <- runif(nsim, Stock@Linfgrad[1], Stock@Linfgrad[2])  # sample of gradient in Linf (Linf y-1)
  # if (!exists("recgrad", inherits=FALSE)) recgrad <- runif(nsim, Stock@recgrad[1], Stock@recgrad[2])  # gradient in recent recruitment
  if (!exists("K", inherits=FALSE)) K <- runif(nsim, Stock@K[1], Stock@K[2])  # now predicted by a log-linear model
  if (!exists("Ksd", inherits=FALSE)) Ksd <- runif(nsim, Stock@Ksd[1], Stock@Ksd[2])  #runif(nsim,Stock@Ksd[1],Stock@Ksd[2])# sd is already added in the linear model prediction
  if (!exists("Kgrad", inherits=FALSE)) Kgrad <- Kgrad <- runif(nsim, Stock@Kgrad[1], Stock@Kgrad[2])  # gradient in Von-B K parameter (K y-1)
  if (!exists("t0", inherits=FALSE)) t0 <- runif(nsim, Stock@t0[1], Stock@t0[2])  # a sample of theoretical age at length zero
  
  # == Sample Maturity Parameters ====
  if (!exists("L50", inherits=FALSE)) {
    sL50 <- array(runif(nsim * 50, Stock@L50[1], Stock@L50[2]), c(nsim, 50))  # length at 50% maturity  
    # checks for unrealistically high length at maturity
    sL50[sL50/Linf > 0.95] <- NA
    L50 <- apply(sL50, 1, function(x) x[!is.na(x)][1])
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
  
  StockOut$Asize <- cbind(StockOut$Size_area_1, 1 - StockOut$Size_area_1)
  
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
  
  # === Generate L50 by year ====
  relL50 <- matrix(L50/Linf, nrow=nsim, ncol=nyears + proyears, byrow=FALSE) # assume L50/Linf stays constant 
  L50array <- relL50 * Linfarray
  delLm <- L95 - L50 
  L95array <- L50array + matrix(delLm, nrow=nsim, ncol=nyears + proyears, byrow=FALSE)
  
  
  # === Create Mean Length-at-Age array ====
  if (!exists("Len_age", inherits=FALSE)) {
    Len_age <- array(NA, dim = c(nsim, maxage, nyears + proyears))  # Length at age array
    ind <- as.matrix(expand.grid(1:nsim, 1:maxage, 1:(nyears + proyears)))  # an index for calculating Length at age
    Len_age[ind] <- Linfarray[ind[, c(1, 3)]] * (1 - exp(-Karray[ind[, c(1, 3)]] * 
                                                           (Agearray[ind[, 1:2]] - t0[ind[, 1]])))
  } else { # Len_age has been passed in with cpars
    if (any(dim(Len_age) != c(nsim, maxage, nyears + proyears))) 
      stop("'Len_age' must be array with dimensions: nsim, maxage, nyears + proyears") 
    # Estimate vB parameters for each year and each sim 
    if (!all(c("Linf", "K", "t0") %in% names(cpars))) { # don't calculate if Linf, K and t0 have also been passed in with cpars
      vB <- function(pars, ages) pars[1] * (1-exp(-pars[2]*(ages-pars[3])))
      fitVB <- function(pars, LatAge, ages) sum((vB(pars, ages) - LatAge)^2)
      starts <- c(max(Len_age), 0.2, 0)
      message("Estimating growth parameters from length-at-age array in cpars")
      for (ss in 1:nsim) {
        pars <- sapply(1:(nyears + proyears), function(X) optim(starts, fitVB, LatAge=Len_age[ss,,X], ages=1:maxage)$par)
        Linfarray[ss,] <- pars[1,]
        Karray[ss,] <- pars[2,]
        t0[ss]<- mean(pars[3,])
      }
      Linf <- Linfarray[, nyears]
      K <- Karray[, nyears]
    }
  }
  StockOut$maxlen <- maxlen <- Len_age[, maxage, nyears] # reference length for Vmaxlen 
  
  # == Sample CV Length-at-age ====
  if (!exists("LenCV", inherits=FALSE)) LenCV <- runif(nsim, min(Stock@LenCV), max(Stock@LenCV))
  
  # == Generate Catch at Length Classes ====
  if (!exists("LatASD", inherits=FALSE)) LatASD <- Len_age * array(LenCV, dim=dim(Len_age)) # SD of length-at-age 
  if (any(dim(LatASD) != dim(Len_age))) stop("Dimensions of 'LatASD' must match dimensions of 'Len_age'", .call=FALSE)
  
  MaxBin <- ceiling(max(Linfarray) + 3 * max(LatASD))
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
    message("Note: Some samples of age of maturity are above 'maxage'-1. Defaulting to maxage-1")
    ageM[ageM >= (maxage-1)] <- maxage - 1 
  }
  if (any(ageM >= maxage)) {
    message("Note: Some samples of age of maturity are above 'maxage'. Defaulting to maxage")
    age95[age95 >= maxage] <- maxage  
  }
  
  # == Generate Maturity-at-Age array ====
  if (!exists("Mat_age", inherits=FALSE)) {
    Mat_age <- array(NA, dim=c(nsim, maxage, nyears+proyears))
    
    for (XX in 1:(nyears+proyears)) {
     Mat_age[,,XX] <- 1/(1 + exp(-log(19) * ((Agearray - ageM[,XX])/(age95[,XX] - ageM[,XX])))) # Maturity at age array by year
    }
    
  } else {
    if (any(dim(Mat_age) != c(nsim, maxage, nyears+proyears))) stop("'Mat_age' must be array with dimensions: nsim, maxage, nyears+proyears") 
    
    # Calculate L50, L95, ageM and age95 
    ageM <- age95 <- L50array <- L95array <- matrix(NA, nsim, nyears+proyears)
    for (XX in 1:(nyears+proyears)) {
      ageM[,XX] <- unlist(sapply(1:nsim, function(x) LinInterp(Mat_age[x,, XX], y=1:maxage, 0.5)))
      age95[,XX] <- unlist(sapply(1:nsim, function(x) LinInterp(Mat_age[x,, XX], y=1:maxage, 0.95)))
      L50array[,XX] <- unlist(sapply(1:nsim, function(x) LinInterp(Mat_age[x,], y=Len_age[x, , nyears], 0.5)))
      L95array[,XX]<- unlist(sapply(1:nsim, function(x) LinInterp(Mat_age[x,], y=Len_age[x, , nyears], 0.95)))
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
    if (exists("M", inherits=FALSE) & length(cpars[["M"]])>0) message("M-at-age has been provided in OM. Overiding M from OM@cpars")
    # M is calculated as mean M of mature ages
    M <- rep(NA, nsim)
    for (sim in 1:nsim) M[sim] <- mean(Mage[sim,round(ageM[sim],0):maxage])
  }
  
  # == Mean Natural mortality by simulation and year ====
  if (!exists("Marray", inherits=FALSE)) {
    Marray <- gettempvar(M, Msd, Mgrad, nyears + proyears, nsim, Mrand)  # M by sim and year according to gradient and inter annual variability
  } else {
    if (any(dim(Marray) != c(nsim, nyears + proyears))) stop("'Marray' must be array with dimensions: nsim, nyears + proyears") 
  }
  
  # == Natural mortality by simulation, age and year ====
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
  
  # == Sample Discard Mortality ====
  if(!exists("Fdisc", inherits = FALSE)) Fdisc <- runif(nsim, min(Stock@Fdisc), max(Stock@Fdisc))
  StockOut$Fdisc <- Fdisc 
  
  
  # Check if M-at-age is constant that Maxage makes sense
  if (all(M_ageArray[1,,1] == mean(M_ageArray[1,,1]))) { # constant M at age
    calcMax <- ceiling(-log(0.01)/(min(M)))        # Age at which 1% of cohort survives
    if (maxage < 0.8*calcMax) {
      message("Note: Maximum age (", maxage, ") is lower than assuming 1% of cohort survives to maximum age (", calcMax, ")")
    }  
  }
  
  StockOut$ageM <- ageM
  StockOut$age95 <- age95
  StockOut$Linfarray <- Linfarray
  StockOut$Karray <- Karray
  StockOut$Agearray <- Agearray
  StockOut$Marray <- Marray
  StockOut$M_ageArray <- M_ageArray
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
  
  return(StockOut)
}
