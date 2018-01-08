utils::globalVariables(c("R0", "Mdb", "mod", "i", "nareas", "nsim", "dFfinal"))
tiny <- 1e-15  # define tiny variable
proportionMat <- TL <- Wa <- SurvWeiMat <- r <- lx <- logNormDensity <- sumlogNormDen <- NULL
proportionMat = vector()


#' Check that a DLM object is valid 
#' 
#' Check that all slots in Object are valid and contain values
#' 
#' @param OM An object of class OM, Stock, Fleet, Obs, or Imp
ChkObj <- function(OM) {
  if (!class(OM) %in% c("OM", "Stock", "Fleet", "Obs", "Imp"))
    stop("Argument must be of class: OM, Stock, Fleet, Obs, or Imp", call.=FALSE)
  
  # Add missing slots with default values 
  OM <- updateMSE(OM)
  slots <- slotNames(OM)
  Ok <- rep(TRUE, length(slots))
  for (sl in seq_along(slots)) {
    slotVal <- slot(OM, slots[sl])
    if (length(slotVal) == 0) Ok[sl] <- FALSE
    if (length(slotVal) > 0) {
      Ok[sl] <- class(slotVal) == class(slot(OM, slots[sl]))
      if (class(slotVal) != "character" & class(slotVal) != "list") Ok[sl] <- all(is.finite(slotVal)) & length(slotVal) > 0
    } 
  }
  optslots <- OptionalSlots()
  SelSlots <- optslots$SelSlots
  RecSlots <-  optslots$RecSlots
  
  # Slots ok to not contain values
  Ignore <- optslots$Ignore
  
  # if values present for one they need to be there for all! 
  if (any(SelSlots %in% slots[Ok])) Ignore <- Ignore[!Ignore %in% SelSlots] 
  if (any(RecSlots %in% slots[Ok])) Ignore <- Ignore[!Ignore %in% RecSlots] 
  
  probSlots <- slots[!Ok][!slots[!Ok] %in% Ignore]
  if (length(probSlots) > 0) 
    stop("Slots in Object have missing values:\n ", paste(probSlots, " "), call.=FALSE)
  return(OM)
  
}

OptionalSlots <- function() {
  SelSlots <- c("SelYears", "L5Lower", "L5Upper", "LFSLower",
                "LFSUpper", "VmaxLower", "VmaxUpper")
  RecSlots <-  c("Period", "Amplitude")
  
  OptPars <- c("M2", "AbsSelYears", 'MPA')
  
  
  # Slots ok to not contain values
  Ignore <- c("Name", "Source", "cpars", SelSlots, RecSlots, OptPars,
              "Agency", "Region", "Latitude", "Longitude", "Species", "Sponsor") 
  out <- list(SelSlots=SelSlots,
              RecSlots=RecSlots,
              OptPars=c(SelSlots, RecSlots, OptPars),
              Ignore=Ignore)
  out
}


# Demographic model
demofn <- function(log.r, M, amat, sigma, K, Linf, to, hR, maxage, a, b) {
  demographic2(log.r, M, amat, sigma, K, Linf, to, hR, maxage = maxage, a, b)$epsilon
}

demographic2 = function(log.r, M, amat, sigma, K, Linf, to, hR, maxage, a, b) {
  # switch on and off to use either S or m in MC simulations
  r = exp(log.r)
  lx = exp(-M)^((1:maxage) - 1)  #survivorship
  logNormDensity = (dnorm(x = log((1:maxage)), mean = log(amat), sd = sigma))/(1:maxage)  #Maturity ogive calculation
  logNormDensity[1] = 0
  sumlogNormDen = sum(logNormDensity)
  NormalisedMaturity = logNormDensity/sumlogNormDen
  proportionMat[1] = NormalisedMaturity[1]
  for (i in 2:maxage) proportionMat[i] = proportionMat[i - 1] + NormalisedMaturity[i]
  TL = Linf * (1 - exp(-K * ((1:maxage) - to)))  #length at age
  Wa = a * TL^b  #wegith at age
  SurvWeiMat = lx * Wa * proportionMat  #survivorship X weight X maturity
  SBPR = sum(SurvWeiMat)  #Spawner biomass per recruit
  RPS = 1/(SBPR * (1 - hR)/(4 * hR))  # Beverton Holt
  # RPS=(5*hR)^(5/4)/SBPR # Ricker Recruitment per spawner biomass
  RPF = Wa * proportionMat * RPS  #Recruits per female
  Lotka = lx * RPF * exp(-(1:maxage) * r)
  sumLotka = sum(Lotka)
  epsilon = (1 - sumLotka)^2  #objective function
  return(list(epsilon = epsilon, r = r))
}

#' Calculate historical fishing mortality
#' 
#' 
#' @param Esd vector of standard deviation 
#' @param nyears number of years 
#' @param EffYears index of years
#' @param EffLower vector of lower bound
#' @param EffUpper vector of upper bound
#' @export getEffhist
#' @keywords internal
#'  
getEffhist <- function(Esd, nyears, EffYears, EffLower, EffUpper) {
  if (length(EffLower) == length(EffUpper) & length(EffUpper) == length(EffYears)) {
    nsim <- length(Esd)  # get nsim 
    refYear <- floor(range01(EffYears + 0.5) * nyears) + 1  # standardize years 
    refYear[length(refYear)] <- nyears  # first year is year 1 
    if (any(EffLower > EffUpper)) {
      ind <- which(EffLower > EffUpper)
      message("Some values in 'EffLower' are higher than 'EffUpper': Years ", paste(ind, ""),
              "\nSetting 'EffLower' to the lower of the two values.")
      tt <- cbind(EffLower, EffUpper)
      EffLower <- apply(tt, 1, min)
      EffUpper <- apply(tt, 1, max)
    }
    
    # sample Effort
    # fmat <- rbind(EffLower, EffUpper)
    
    # nyrs <- length(EffLower)
    # Effs <- matrix(0, nsim, nyrs)
    
    # ind <- which(diff(fmat) > 0)[1]
    # for (X in 1:ind) {
    # Effs[,X] <- runif(nsim, min(fmat[,X]), max(fmat[,X]))  
    # }
    
    # val <- (Effs[,ind] - min(fmat[,ind]))/ diff(fmat[,ind])
    # for (X in 2:nyrs) Effs[,X] <- min(fmat[,X]) + diff(fmat[,X])*val
    
    Effs <- mapply(runif, n = nsim, min = EffLower, max = EffUpper)  # sample Effort
    if (nsim > 1) {
      effort <- t(sapply(1:nsim, function(x) approx(x = refYear, 
                                                    y = Effs[x, ], method = "linear", n = nyears)$y))  # linear interpolation
    }
    
    
    if (nsim == 1) {
      # Effs <- Effs/max(Effs)
      effort <- matrix(approx(x = refYear, y = Effs, method = "linear", n = nyears)$y, nrow = 1)
    }
    
    effort <- range01(effort)
    effort[effort == 0] <- 0.01
    
    Emu <- -0.5 * Esd^2
    Eerr <- array(exp(rnorm(nyears * nsim, rep(Emu, nyears), rep(Esd, nyears))), c(nsim, nyears))  # calc error
    out <- NULL
    eff <- effort * Eerr  # add error 
    out[[1]] <- eff
    out[[2]] <- (effort[, nyears] - effort[, nyears - 4])/5
    return(out)
  } else {
    message("Input vectors of effort years and bounds not of same length")
    return(NULL)
  }
}

#' get object class
#' 
#' Internal function for determining if object is of classy
#' 
#' 
#' @param x Character string object name
#' @param classy A class of object (character string, e.g. 'Fleet')
#' @author T. Carruthers
#' @return TRUE or FALSE
getclass <- function(x, classy) any(class(get(x)) == classy) # inherits(get(x), classy) - this gives a problem since we now inherit Stock etc in OM

#' Optimization function to find a movement model that matches user specified
#' movement characteristics modified for Rcpp.
#' 
#' The user specifies the probability of staying in the same area and spatial
#' heterogeneity (both in the unfished state).
#' 
#' This is paired with movfit to find the correct movement model.
#' 
#' @param x A position in vectors Prob_staying and Frac_area_1
#' @param Prob_staying User specified probability that individuals in area 1
#' remain in that area (unfished conditions)
#' @param Frac_area_1 User specified fraction of individuals found in area 1
#' (unfished conditions)
#' @return A markov movement matrix
#' @author T. Carruthers
#' @export
#' @examples
#' 
#' Prob_staying<-0.8 # probability  that individuals remain in area 1 between time-steps
#' Frac_area_1<-0.35 # the fraction of the stock found in area 1 under equilibrium conditions
#' markovmat<-getmov2(1,Prob_staying, Frac_area_1)
#' vec<-c(0.5,0.5) # initial guess at equilibrium distribution (2 areas)
#' for(i in 1:300)vec<-apply(vec*markovmat,2,sum) # numerical approximation to stable distribution
#' c(markovmat[1,1],vec[1]) # pretty close right?
#' 
#' 
getmov2 <- function(x, Prob_staying, Frac_area_1) {
  test <- optim(par = c(0, 0, 0), movfit_Rcpp, method = "L-BFGS-B", 
                lower = rep(-6, 3), upper = rep(6, 3), prb = Prob_staying[x], 
                frac = Frac_area_1[x])	
  mov <- array(c(test$par[1], test$par[2], 0, test$par[3]), dim = c(2, 2))
  mov <- exp(mov)
  mov/array(apply(mov, 1, sum), dim = c(2, 2))
}

#' Creates a time series per simulation that has gradient grad and random normal walk with sigma
#' 
#' @param targ mean
#' @param targsd standard deviation
#' @param targgrad gradient 
#' @param nyears number of historical years
#' @param nsim number of simulations
#' @keywords internal
#'  
gettempvar <- function(targ, targsd, targgrad, nyears, nsim, rands=NULL) {
  mutemp <- -0.5 * targsd^2
  temp <- array(1, dim = c(nsim, nyears))
  if (is.null(rands)) {
    for (i in 2:nyears) {
      temp[, i] <- temp[, i] * exp(rnorm(nsim, mutemp, targsd))
    }
  }
  if (!is.null(rands)) {
    for (i in 2:nyears) {
      temp[, i] <- temp[, i] * rands[,i]
    }
  }
  yarray <- array(rep((1:nyears) - 1, each = nsim), dim = c(nsim, nyears))
  temp <- temp * (1 + targgrad/100)^yarray
  targ * temp/apply(temp, 1, mean)
}


#' Linear interpolation of a y value at level xlev based on a vector x and y
#'
#' @param x A vector of x values
#' @param y A vector of y values (identical length to x)
#' @param xlev A the target level of x from which to guess y
#' @param ascending Are the the x values supposed to be ordered before interpolation
#' @param zeroint is there a zero-zero x-y intercept?
#' @author T. Carruthers
#' @export LinInterp
LinInterp<-function(x,y,xlev,ascending=F,zeroint=F){
  
  if(zeroint){
    x<-c(0,x)
    y<-c(0,y)
  } 
  
  if(ascending){
    cond<-(1:length(x))<which.max(x)
  }else{
    cond<-rep(TRUE,length(x))
  }
  
  close<-which.min((x[cond]-xlev)^2)
  ind<-c(close,close+(x[close]<xlev)*2-1)
  ind <- ind[ind <= length(x)]
  if (length(ind)==1) ind <- c(ind, ind-1)
  ind<-ind[order(ind)]
  pos<-(xlev-x[ind[1]])/(x[ind[2]]-x[ind[1]])
  y[ind[1]]+pos*(y[ind[2]]-y[ind[1]])
  
}

#'  Condition met?
#' 
#' @param vec vector of logical values 
condmet <- function(vec) TRUE %in% vec

#'  Sample vector
#' 
#' @param x vector of values 
sampy <- function(x) sample(x, 1, prob = !is.na(x))


#' Standardize values
#' 
#' Function to standardize to value relative to minimum and maximum values
#' @param x vector of values 
#' @param Max Maximum value
#' @param Min Minimum value 
Range <- function(x, Max, Min) {
  (x - Min)/(Max - Min)  
}

range01 <- function(x) {
  (x - min(x))/(max(x) - min(x)) 
}

#' runMSE with no messages - for testing 
#'
#' For testing purposes only 
#' @param ... Arguments to runMSE function 
#'
#' @keywords internal
#' @importFrom utils capture.output
#'
runMSEnomsg <- function(...) {
  capture.output(out <- suppressMessages(runMSE(...)))
  out
}


## Selectivity functions ####

#' Double-normal selectivity curve
#'
#' @param lens Vector of lengths 
#' @param lfs Length at full selection
#' @param sl Sigma of ascending limb
#' @param sr Sigma of descending limb
#'
#'
dnormal<-function(lens,lfs,sl,sr){
  cond<-lens<=lfs
  sel<-rep(NA,length(lens))
  sel[cond]<-2.0^-((lens[cond]-lfs)/sl*(lens[cond]-lfs)/sl)
  sel[!cond]<-2.0^-((lens[!cond]-lfs)/sr*(lens[!cond]-lfs)/sr)
  sel
}


#' Calculate selectivity curve
#'
#' @param x Simulation number
#' @param lens Matrix of lengths (nsim by nlengths)
#' @param lfs Vector of length at full selection (nsim long)
#' @param sls Vector of sigmas of ascending limb (nsim long)
#' @param srs Vector of sigmas of descending limb (nsim long)
#'
#'
getsel <- function(x, lens, lfs, sls, srs) {
  if (is.null(ncol(lens))) return(dnormal(lens, lfs[x], sls[x], srs[x]))
  dnormal(lens[x,], lfs[x], sls[x], srs[x])
}











