utils::globalVariables(c("R0", "Mdb", "mod", "i", "nareas", "nsim", "dFfinal"))
tiny <- 1e-15  # define tiny variable
proportionMat <- TL <- Wa <- SurvWeiMat <- r <- lx <- logNormDensity <- sumlogNormDen <- NULL
proportionMat = vector()

MPCheck <- function(MPs, Data, timelimit, silent=FALSE) {
  if(!silent) message("Determining available methods") 
  PosMPs <- Can(Data, timelimit = timelimit)  # list all the methods that could be applied
  if (is.na(MPs[1])) {
    MPs <- PosMPs  # if the user does not supply an argument MPs run the MSE for all available methods
    if(!silent) message("No MPs specified: running all available")
  }
  cant <- MPs[!MPs %in% PosMPs]
  if (length(cant) > 0) {
    if(!silent) message("Cannot run some MPs: ")
    if(!silent) print(DLMdiag(Data, "not available", funcs1=cant, timelimit = timelimit))
  }
  MPs <- MPs[MPs %in% PosMPs]  # otherwise run the MSE for all methods that are deemed possible
  if (length(MPs) == 0) {
    if(!silent) message(Cant(Data, timelimit = timelimit))
    stop("MSE stopped: no viable methods \n\n", call. = FALSE)  # if none of the user specified methods are possible stop the run
  }
  
  ok <- rep(TRUE, length(MPs))
  for (mm in seq_along(MPs)) {
    test <- try(get(MPs[mm]), silent=TRUE)
    if (!class(test) == 'MP') {
      ok[mm] <- FALSE
      if (class(test) == 'try-error') {
        message('Object ', paste(MPs[mm], ""), " does not exist - Ignoring")
      } else message('Dropping MP: ', paste(MPs[mm], ""), " - Not class 'MP'")
    }
  }
  MPs <- MPs[ok]
  MPs
}

#' What Data objects can be used to run this MP?
#'
#' @param MP Name of the MP
#'
#' @export
#'
#' @keywords internal
CanMP <- function(MP) {
  avData <- avail("Data")
  log <- rep(FALSE, length(avData))
  for (x in seq_along(avData)) {
    log[x] <- MP %in% Can(get(avData[x]))
  }
  return(avData[log])
}



#' Check that a DLM object is valid 
#' 
#' Check that all slots in Object are valid and contain values
#' 
#' @param OM An object of class OM, Stock, Fleet, Obs, or Imp
#' @param error Logical. Stop on missing parameter values? FALSE = warning
ChkObj <- function(OM, error=TRUE) {
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
  Ignore <- c(Ignore, "Mgrad", "Kgrad", "Linfgrad")
  
  # if values present for one they need to be there for all! 
  if (any(SelSlots %in% slots[Ok])) Ignore <- Ignore[!Ignore %in% SelSlots] 
  if (any(RecSlots %in% slots[Ok])) Ignore <- Ignore[!Ignore %in% RecSlots] 
  
  probSlots <- slots[!Ok][!slots[!Ok] %in% Ignore]
  probSlots <- probSlots[!probSlots %in% names(OM@cpars)]
  if ('Len_age' %in% names(OM@cpars)) {
    probSlots <- probSlots[!probSlots %in% c("Linf", "K", "t0")]
  }
  
  if (length(probSlots) > 0) {
    if (error) stop("Slots in Object have missing values:\n ", paste(probSlots, " "), call.=FALSE)
    if (!error) warning("Slots in Object have missing values:\n ", paste(probSlots, " "), call.=FALSE)
  }

  
  ## Add variability - variability required in all values otherwise
  # OMs are not reproducible. 
  # for (sl in slotNames(OM)) {
  #   slt <- slot(OM, sl)
  #   if (class(slt) == "numeric") {
  #     if (length(slt) ==2) {
  #       if (!all(is.na(slt)) && slt[1] == slt[2]) {
  #         slt[1] <- slt[1]
  #         slt[2] <- slt[2]+tiny
  #       }
  #      
  #     }
  #     if (length(slt)==1) {
  #       if (!all(is.na(slt)) && slt == 0) slt <- tiny
  #     }
  #     slot(OM,sl) <- slt
  #   }
  # }?
  OM
  
  
}

OptionalSlots <- function() {
  SelSlots <- c("SelYears", "L5Lower", "L5Upper", "LFSLower",
                "LFSUpper", "VmaxLower", "VmaxUpper")
  RecSlots <-  c("Period", "Amplitude")
  
  OptPars <- c("M2", "Mexp", "AbsSelYears", 'MPA')
  
  
  # Slots ok to not contain values
  Ignore <- c("Name", "Source", "cpars", SelSlots, RecSlots, OptPars,
              "Agency", "Region", "Latitude", "Longitude", "Species", "Sponsor", "Common_Name") 
  out <- list(SelSlots=SelSlots,
              RecSlots=RecSlots,
              OptPars=c(SelSlots, RecSlots, OptPars),
              Ignore=Ignore)
  out
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
      if (ncol(Effs) == 1) {
        effort <- matrix(Effs, nrow=nsim, ncol=nyears)
      } else {
        effort <- t(sapply(1:nsim, function(x) approx(x = refYear, 
                                                      y = Effs[x, ], method = "linear", n = nyears)$y))  # linear interpolation
      }
      
    } 
    if (nsim == 1) {
      if (length(Effs) == 1) {
        effort <- matrix(Effs, nrow=nsim, ncol=nyears)
      } else {
        effort <- matrix(approx(x = refYear, y = Effs, method = "linear", n = nyears)$y, nrow = 1)
      }
    }
  
    if (!all(effort == mean(effort))) effort <- range01(effort)  
    
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
#' @author T. Carruthers with nasty hacks from A. Hordyk
#' @return TRUE or FALSE
getclass <- function(x, classy) {
  return(any(class(get(x)) == classy)) # inherits(get(x), classy) - this gives a problem since we now inherit Stock etc in OM
} 
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

#' Creates a time series per simulation that has a random normal walk with sigma
#' 
#' @param targ mean
#' @param targsd standard deviation
#' @param targgrad gradient - no longer used
#' @param nyears number of years to simulate
#' @param nsim number of simulations
#' @keywords internal
#'  
gettempvar <- function(targ, targsd, targgrad, nyears, nsim, rands=NULL) {
  mutemp <- -0.5 * targsd^2
  
  if (!is.null(rands)) {
    temp <- rands
  } else {
    temp <- array(exp(rnorm(nsim*nyears, mutemp, targsd)),dim = c(nsim, nyears))
  }
  
  # yarray <- array(rep((1:nyears) - 1, each = nsim), dim = c(nsim, nyears))
  # temp <- temp * (1 + targgrad/100)^yarray
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
#' @keywords internal
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
#' @keywords internal
sampy <- function(x) sample(x, 1, prob = !is.na(x))


#' Standardize values
#' 
#' Function to standardize to value relative to minimum and maximum values
#' @param x vector of values 
#' @param Max Maximum value
#' @param Min Minimum value 
#' @keywords internal
Range <- function(x, Max, Min) {
  (x - Min)/(Max - Min)  
}

range01 <- function(x) {
  (x - min(x))/(max(x) - min(x)) 
}


runMSEnomsg <- function(...) {
  capture.output(out <- suppressMessages(runMSE(...)))
  out
}


run_parallel <- function(i, itsim, OM, MPs, CheckMPs, timelimit, Hist, ntrials, fracD, CalcBlow, 
                         HZN, Bfrac, AnnualMSY, silent, PPD, control, parallel=FALSE) {
  OM@nsim <- itsim[i]
  
  OM@seed <- OM@seed + i 
  mse <- runMSE_int(OM, MPs, CheckMPs, timelimit, Hist, ntrials, fracD, CalcBlow, 
                    HZN, Bfrac, AnnualMSY, silent, PPD=PPD, control=control, parallel=parallel)
  return(mse)
}

assign_DLMenv <- function() {
  DLMenv_list <- snowfall::sfClusterEval(mget(ls(DLMenv), envir = DLMenv)) # Grab objects from cores' DLMenv
  clean_env <- snowfall::sfClusterEval(rm(list = ls(DLMenv), envir = DLMenv)) # Remove cores' DLMenv objects
  env_names <- unique(do.call(c, lapply(DLMenv_list, names)))
  
  if(length(env_names) > 0) {
    for(i in 1:length(env_names)) {
      temp <- lapply(DLMenv_list, getElement, env_names[i])
      assign(env_names[i], do.call(c, temp), envir = DLMenv) # Assign objects to home DLMenv
    }
  }
  
  return(invisible(env_names))
}


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



# Generate size comps
genSizeCompWrap <- function(i, vn, CAL_binsmid, retL,
                            CAL_ESS, CAL_nsamp,
                            Linfarray, Karray, t0array,
                            LenCV, truncSD=2) {
  VulnN <- as.matrix(vn[i,,]) 
  if (ncol(VulnN)>1) {
    VulnN <- VulnN/rowSums(VulnN) * CAL_nsamp[i] # get relative numbers at age   
  } else {
    VulnN <- VulnN/sum(VulnN) * CAL_nsamp[i] # get relative numbers at age 
  }
  
  VulnN <- round(VulnN,0) # convert to integers
  nyrs <- nrow(as.matrix(Linfarray[i,]))
  if (nyrs == 1) VulnN <- t(VulnN)
  retLa <- as.matrix(retL[i,,])
  
  lens <- genSizeComp(VulnN, CAL_binsmid, retLa,
              CAL_ESS=CAL_ESS[i], CAL_nsamp=CAL_nsamp[i],
              Linfs=Linfarray[i,], Ks=Karray[i,], t0s=t0array[i,],
              LenCV=LenCV[i], truncSD)
  
  lens[!is.finite(lens)] <- 0
  lens
  
}

getfifth <- function(lenvec, CAL_binsmid) {
  temp <- rep(CAL_binsmid, lenvec)
  if(sum(lenvec)==0) return(NA)
  dens <- try(density(temp), silent=TRUE)
  if(class(dens)!="density") return(NA)
  dens$x[min(which(cumsum(dens$y/sum(dens$y)) >0.05))]
}

# genSizeCompWrap2<- function(i, vn, CAL_binsmid, 
#                             CAL_ESS, CAL_nsamp,
#                             Linfarray, Karray, t0array,
#                             LenCV, truncSD=2) {
#   
#   VulnN <- as.matrix(vn[i,,]) 
#   VulnN <- round(VulnN,0)
#   nyrs <- nrow(as.matrix(Linfarray[i,]))
#   if (nyrs == 1) VulnN <- t(VulnN)
#   
#   
#   lens <- genSizeComp2(VulnN, CAL_binsmid, 
#                       CAL_ESS=CAL_ESS[i], CAL_nsamp=CAL_nsamp[i],
#                       Linfs=Linfarray[i,], Ks=Karray[i,], t0s=t0array[i,],
#                       LenCV=LenCV[i], truncSD)
#   
#   lens
#   
# }


userguide_link <- function(url, ref=NULL) {
  url <- paste0('https://dlmtool.github.io/DLMtool/userguide/', url, '.html')
  if (ref!="NULL") url <- paste0(url, "#", ref)
  paste0("See relevant section of the \\href{", url, "}{DLMtool User Guide} for more information.")
}

  
  
#' Simulate Catch-at-Age Data
#'
#' CAA generated with a multinomial observation model from retained catch-at-age
#' data
#'
#' @param nsim Number of simulations
#' @param yrs Number of years 
#' @param maxage Maximum age
#' @param Cret Retained Catch at age in numbers - array(sim, years, maxage)
#' @param CAA_ESS CAA effective sample size 
#' @param CAA_nsamp CAA sample size
#'
#' @return CAA array 
simCAA <- function(nsim, yrs, maxage, Cret, CAA_ESS, CAA_nsamp) {
  # generate CAA from retained catch-at-age 

  CAA <- array(NA, dim = c(nsim, yrs, maxage))  # Catch  at age array
  
  # a multinomial observation model for catch-at-age data
  for (i in 1:nsim) {
    for (j in 1:yrs) {
      if (!sum(Cret[i, j,])) {
        CAA[i, j, ] <- 0 
      } else {
        CAA[i, j, ] <- ceiling(-0.5 + rmultinom(1, CAA_ESS[i], Cret[i, j,]) * CAA_nsamp[i]/CAA_ESS[i])   
      }
    }
  }
  CAA
}

#' Simulate Catch-at-Length Data
#' 
#' Simulate CAL and calculate length-at-first capture (LFC),
#' mean length (ML), modal length (Lc), and mean length over modal length (Lbar)  
#' 
#' @param nsim Number of simulations
#' @param nyears Number of years 
#' @param maxage Maximum age
#' @param CAL_ESS CAA effective sample size 
#' @param CAL_nsamp CAA sample size
#' @param nCALbins number of CAL bins
#' @param CAL_binsmid mid-points of CAL bins
#' @param vn Vulnerable numbers-at-age
#' @param retL Retention at length curve
#' @param Linfarray Array of Linf values by simulation and year
#' @param Karray Array of K values by simulation and year
#' @param t0array Array of t0 values by simulation and year
#' @param LenCV CV of length-at-age#'
#' @return named list with CAL array and LFC, ML, & Lc vectors
simCAL <- function(nsim, nyears, maxage,  CAL_ESS, CAL_nsamp, nCALbins, CAL_binsmid,  
                   vn, retL, Linfarray, Karray, t0array, LenCV) {
  # a multinomial observation model for catch-at-length data
  # assumed normally-distributed length-at-age truncated at 2 standard deviations from the mean
  CAL <- array(NA, dim=c(nsim,  nyears, nCALbins))
  
  # Generate size comp data with variability in age
 
  tempSize <- lapply(1:nsim, genSizeCompWrap, vn, CAL_binsmid, retL, CAL_ESS, CAL_nsamp,
                     Linfarray, Karray, t0array, LenCV, truncSD=2)
  CAL <- aperm(array(as.numeric(unlist(tempSize, use.names=FALSE)), 
                     dim=c(nyears, length(CAL_binsmid), nsim)), c(3,1,2))
  
  # calculate LFC - length-at-first capture - 5th percentile
  LFC <- rep(NA, nsim)
  LFC <- unlist(lapply(tempSize, function(x) getfifth(x[nyears, ], CAL_binsmid)))
  LFC[is.na(LFC)] <- 1
  LFC[LFC<1] <- 1
  
  # Mean Length 
  temp <- CAL * rep(CAL_binsmid, each = nsim * nyears)
  ML <- apply(temp, 1:2, sum)/apply(CAL, 1:2, sum)
  ML[!is.finite(ML)] <- 0 
  
  # Lc - modal length 
  Lc <- array(CAL_binsmid[apply(CAL, 1:2, which.max)], dim = c(nsim, nyears))
  
  # Lbar - mean length above Lc
  nuCAL <- CAL
  for (i in 1:nsim) {
    for (j in 1:nyears) {
      # nuCAL[i, j, 1:match(max(1, Lc[i, j]), CAL_binsmid, nomatch=1)] <- NA
      lcbin <- max(1,match(max(1, Lc[i, j]), CAL_binsmid, nomatch=1)-1)
      nuCAL[i, j, 1:lcbin] <- NA
    }
  }
  
  temp <- nuCAL * rep(CAL_binsmid, each = nsim * nyears)
  Lbar <- apply(temp, 1:2, sum, na.rm=TRUE)/apply(nuCAL, 1:2, sum, na.rm=TRUE)
  Lbar[!is.finite(Lbar)] <- 0 
  out <- list()
  out$CAL <- CAL
  out$LFC <- LFC
  out$ML <- ML 
  out$Lc <- Lc
  out$Lbar <- Lbar
  
  out
}


# initialize development mode
dev.mode <- function() {
  devtools::load_all()
  DFargs <- formals(runMSE_int)
  argNames <- names(DFargs)
  ind <- which(unlist(lapply(argNames, exists)))
  DFargs[ind] <- NULL
  ind <- which(lapply(DFargs, class) == 'call')
  for (x in ind) {
    DFargs[[x]] <- eval((DFargs[[x]]))
  }
  DFargs
}


getbeta<-function(beta,x,y)sum((y-x^beta)^2)

indfitwrap <- function(x, type, sim.indices, ind.type, Data, nyears, plot=FALSE) {
  sim.index <- sim.indices[[match(type, ind.type)]][x,]
  dat <- Data@RInd[x,match(type, Data@Type),]
  dat <- dat[!is.na(dat)]
  nyears <- min(nyears, length(dat))
  obs.ind <- Data@RInd[x,match(type, Data@Type), 1:nyears]
  sim.index <- sim.index[1:nyears]
  Year <- Data@Year
  indfit(sim.index,obs.ind, Year, plot, lcex=0.8)
}

lcs<-function(x){
  if (class(x) == "matrix") {
    nsim <- nrow(x)
    nyr <- ncol(x)
    x1 <- x/matrix(apply(x, 1, mean), nrow=nsim, ncol=nyr) # rescale to mean 1
    x2<- log(x1) # log it
    x3 <- x2 -matrix(apply(x2, 1, mean), nrow=nsim, ncol=nyr) # mean 0
    x3
  } else {
    x1<-x/mean(x) # rescale to mean 1
    x2<-log(x1)     # log it
    x3<-x2-mean(x2) # mean 0
    x3
  }

}



makeVec <- function(obj, row=1, nsim) {
  tt <- obj[row,] %>% unlist()
  if (is.null(tt)) return(rep(NA, nsim))
  return(tt)
}

indfit <- function(sim.index,obs.ind, Year, plot=FALSE, lcex=0.8){
  
  sim.index <- lcs(sim.index) # log space conversion of standardized simulated index
  obs.ind <- lcs(obs.ind) # log space conversion of standardized observed ind
  
  if(plot){
    par(mfrow=c(1,2),mai=c(0.7,0.5,0.05,0.01),omi=c(0.01,0.2,0.01,0.01))
    plot(exp(sim.index),exp(obs.ind),xlab="",ylab="",pch=19,col=rgb(0,0,0,0.5))
    mtext("Model estimate",1,line=2.2)
    mtext("Index",2,outer=T,line=0)
  }
  
  opt<-optimize(getbeta,x=exp(sim.index),y=exp(obs.ind),interval=c(0.1,10))
  res<-exp(obs.ind)-(exp(sim.index)^opt$minimum)
  ac<-acf(res,plot=F)$acf[2,1,1] # lag-1 autocorrelation
  
  res2<-obs.ind-sim.index                  # linear, without hyperdepletion / hyperstability
  ac2<-acf(res2,plot=F)$acf[2,1,1] # linear AC
  
  if(plot){
    SSBseq<-seq(min(exp(sim.index)),max(exp(sim.index)),length.out=1000)
    lines(SSBseq,SSBseq^opt$minimum,col='#0000ff90',pch=19)
    legend('bottomright',legend=round(c(sum((obs.ind-sim.index)^2),opt$objective),3),text.col=c("black","blue"),bty='n',title="SSQ",cex=lcex)
    legend('topleft',legend=round(opt$minimum,3),text.col="blue",bty='n',title='Hyper-stability, beta',cex=lcex)
    legend('left',legend=round(stats::cor(sim.index,obs.ind),3),bty='n',title='Correlation',cex=lcex)
    
    plot(Year,sim.index,ylab="",xlab="",ylim=range(c(obs.ind,sim.index)),type="l")
    mtext("Year",1,line=2.2)
    points(Year,obs.ind,col='#ff000090',pch=19)
    legend('topleft',legend=round(ac,3),text.col="red",bty='n',title="Lag 1 autocorrelation",cex=lcex)
    legend('bottomleft',legend=round(sd(res),3),text.col="red",bty='n',title="Residual StDev",cex=lcex)
    legend('topright',legend=c("Model estimate","Index"),text.col=c("black","red"),bty='n',cex=lcex)
  }
  
  data.frame(beta=opt$minimum,AC=ac,sd=sd(exp(obs.ind)/(exp(sim.index)^opt$minimum)),
             cor=stats::cor(sim.index,obs.ind),AC2=ac2,sd2=sd(obs.ind-sim.index))
  
  # list(stats=data.frame(beta=opt$minimum,AC=ac,sd=sd(exp(obs.ind)/(exp(sim.index)^opt$minimum)),
  #                       cor=cor(sim.index,obs.ind),AC2=ac2,sd2=sd(obs.ind-sim.index)),
  #      mult=exp(obs.ind)/(exp(sim.index)^opt$minimum))
  
}

generateRes <- function(df, nsim, proyears, lst.err) {
  sd <- df$sd 
  ac <- df$ac
  if (all(is.na(sd))) return(rep(NA, nsim))
  mu <- 0.5 * (sd)^2
  Res <- matrix(rnorm(proyears*nsim, mu, sd), nrow=proyears, ncol=nsim, byrow=TRUE) 
  # apply a pseudo AR1 autocorrelation 
  Res <- sapply(1:nsim, applyAC, res=Res, ac=ac, max.years=proyears, lst.err=lst.err) # log-space
  exp(t(Res))
}

applyAC <- function(x, res, ac, max.years, lst.err) {
  for (y in 1:max.years) {
    if (y == 1) {
      res[y,x] <- ac[x] * lst.err[x] + lst.err[x] * (1-ac[x] * ac[x])^0.5 
    } else {
      res[y,x] <- ac[x] * res[y-1,x] + res[y,x] * (1-ac[x] * ac[x])^0.5  
    }
    
  }
  res[,x]
}

addRealInd <- function(Data, SampCpars, ErrList, Biomass, VBiomass, SSB, nsim, nyears,
                       proyears, silent=FALSE) {
  if (!is.null(SampCpars$Data)) {
    # real data has been provided 
    types <- SampCpars$Data@Type
    chk <- types %in% c("Biomass", "VBiomass", "SpBiomass")
    if (any(!chk)) 
      if (!is.na(types)) 
        stop("Invalid index type in cpars$Data@Type. \nValid types are: 'Biomass', 'VBiomass', 'SpBiomass'", call.=FALSE)
    if (all(is.na(types))) {
      if (!silent) message("Data object provided in cpars but Data@Type is empty. Ignoring Data")
    } else {
      indices <- SampCpars$Data@RInd[1,,]
      if (all(is.na(indices))) {
        if (!silent) message("Data object provided in cpars but Data@RInd is empty. Ignoring Data")
      } else {
        if (!silent) message("Conditioning simulated RInd on cpars$Data@RInd")
        # add Real indices to Data object 
        Data@Type <- SampCpars$Data@Type
        RInd <- SampCpars$Data@RInd[1,,]
        if (class(RInd) == "numeric") {
          Data@RInd <- array(RInd, dim=c(nsim, 1, length(RInd)))
        } else {
          temp <- array(RInd, dim=c(nrow(RInd), ncol(RInd), nsim))
          Data@RInd <- aperm(temp, c(3,1,2)) 
        }
        
        # Calculate implied observation error
        sim.indices <- list(Biomass=apply(Biomass, c(1, 3), sum),
                            VBiomass=apply(VBiomass, c(1, 3), sum),
                            SpBiomass=apply(SSB, c(1, 3), sum))
        
        ind.type <- c('Biomass', 'VBiomass', 'SpBiomass')
        out.list <- list()
        for (type in Data@Type) {
          out.list[[type]] <- sapply(1:nsim, indfitwrap, type=type, 
                                     sim.indices=sim.indices, 
                                     ind.type=ind.type, Data=Data, nyears=nyears)
        } 
        
        # Make data.frame 
        stats.df <- data.frame(Index=rep(c("Biomass", "VBiomass", 'SpBiomass'), each=nsim),
                               beta=c(makeVec(out.list$Biomass, 1, nsim),
                                      makeVec(out.list$VBiomass, 1, nsim),
                                      makeVec(out.list$SpBiomass, 1, nsim)),
                               ac = c(makeVec(out.list$Biomass, 2, nsim),
                                      makeVec(out.list$VBiomass, 2, nsim),
                                      makeVec(out.list$SpBiomass, 2, nsim)),
                               sd = c(makeVec(out.list$Biomass, 3, nsim),
                                      makeVec(out.list$VBiomass, 3, nsim),
                                      makeVec(out.list$SpBiomass, 3, nsim)),
                               cor =c(makeVec(out.list$Biomass, 4, nsim),
                                      makeVec(out.list$VBiomass, 4, nsim),
                                      makeVec(out.list$SpBiomass, 4, nsim)))
        
        # Generate index random error
        BIerr <- array(NA, dim=c(nsim, nyears+proyears))
        VBIerr <- array(NA, dim=c(nsim, nyears+proyears))
        SBIerr <- array(NA, dim=c(nsim, nyears+proyears))
        
        # Calculate error where data exists
        
        bbeta <- stats.df %>% filter(Index=="Biomass") %>% select(beta)
        yrs1 <- min(nyears, length(Data@RInd[1,match("Biomass", Data@Type),]))
        BErr <- exp(lcs(Data@RInd[,match("Biomass", Data@Type),1:yrs1]))/exp(lcs(sim.indices$Biomass[, 1:yrs1]))^bbeta$beta
        
        bbeta <- stats.df %>% filter(Index=="VBiomass") %>% select(beta)
        yrs2 <- min(nyears, length(Data@RInd[1,match("VBiomass", Data@Type),]))
        VBErr <- exp(lcs(Data@RInd[,match("VBiomass", Data@Type),1:yrs2]))/exp(lcs(sim.indices$VBiomass[, 1:yrs2]))^bbeta$beta
        
        bbeta <- stats.df %>% filter(Index=="SpBiomass") %>% select(beta)
        yrs3 <- min(nyears, length(Data@RInd[1,match("SpBiomass", Data@Type),]))
        SpBErr <- exp(lcs(Data@RInd[,match("SpBiomass", Data@Type),1:yrs3]))/exp(lcs(sim.indices$SpBiomass[, 1:yrs3]))^bbeta$beta
        
        BIerr[, 1:ncol(BErr)] <- BErr
        VBIerr[, 1:ncol(VBErr)] <- VBErr
        SBIerr[, 1:ncol(SpBErr)] <- SpBErr
        
        # Generate error for projection years 
        BIerr[,(ncol(BErr)+1):ncol(BIerr)] <- generateRes(stats.df %>% filter(Index=="Biomass"), 
                                                          nsim, length((ncol(BErr)+1):ncol(BIerr)), 
                                                          lst.err=log(BIerr[,yrs1]))
        VBIerr[,(ncol(BErr)+1):ncol(VBIerr)] <- generateRes(stats.df %>% filter(Index=="VBiomass"),
                                                            nsim, length((ncol(BErr)+1):ncol(BIerr)), 
                                                            lst.err=log(VBIerr[,yrs2]))
        SBIerr[,(ncol(BErr)+1):ncol(SBIerr)] <- generateRes(stats.df %>% filter(Index=="SpBiomass"),
                                                            nsim, length((ncol(BErr)+1):ncol(BIerr)), 
                                                            lst.err=log(SBIerr[,yrs3]))
        
        ErrList$stats.df <- stats.df
        ErrList$RIerr <- list(BIerr=BIerr, VBIerr=VBIerr, SBIerr=SBIerr)
      }
    }
  }
  return(list(Data=Data, ErrList=ErrList))
}


# calculate average unfished ref points over first A50 years
CalcUnfishedRefs <- function(x, ageM, N0_a, SSN0_a, SSB0_a, B0_a, VB0_a, SSBpRa, SSB0a_a) {
  avg.ind <- 1:ceiling(ageM[x,1]) # unfished eq ref points averaged over these years 
  nyears <- dim(N0_a)[2]
  if (length(avg.ind) > nyears) avg.ind <- 1:nyears
  
  N0 <- mean(N0_a[x, avg.ind])
  SSN0  <- mean(SSN0_a[x, avg.ind])
  SSB0 <- mean(SSB0_a[x, avg.ind])
  B0 <- mean(B0_a[x, avg.ind])
  VB0 <- mean(VB0_a[x, avg.ind])
  SSBpR <- mean(SSBpRa[x, avg.ind])
  if (length(avg.ind)>1) {
    SSB0a <- apply(SSB0a_a[x,avg.ind,], 2, mean)  
  } else {
    SSB0a <- SSB0a_a[x,avg.ind,]
  }
  
  list(N0=N0, SSN0=SSN0, SSB0=SSB0, B0=B0, VB0=VB0, SSBpR=SSBpR, SSB0a=SSB0a)
}

CalcMSYRefs <- function(x, MSY_y, FMSY_y, SSBMSY_y, BMSY_y, VBMSY_y, ageM, OM) {
  n.yrs <- ceiling(ageM[x,OM@nyears]) # MSY ref points averaged over these years
  nyears <- dim(ageM)[2]
  minY <- floor(n.yrs/2) 
  maxY <- n.yrs - minY - 1 
  avg.ind <- (OM@nyears - minY):(OM@nyears + maxY)
  if (max(avg.ind) > nyears) avg.ind <- avg.ind[avg.ind < nyears]
  
  MSY <- mean(MSY_y[x, avg.ind])
  FMSY <- mean(FMSY_y[x, avg.ind])
  SSBMSY <- mean(SSBMSY_y[x, avg.ind])
  BMSY <- mean(BMSY_y[x, avg.ind])
  VBMSY <- mean(VBMSY_y[x, avg.ind])
  data.frame(MSY=MSY, FMSY=FMSY, SSBMSY=SSBMSY, BMSY=BMSY, VBMSY=VBMSY)
}

# 
# makeSizeCompW <- function(i, maxage, Linfarray, Karray, t0array, LenCV,
#                           CAL_bins, CAL_binsmid, retL, CAL_ESS, CAL_nsamp, 
#                           vn, truncSD=2, scaleR0=1) {
#   
#   if(length(scaleR0)==1) scaleR0 <- rep(scaleR0, i)
#   AgeVec <- 1:maxage 
#   SubAgeVec <- seq(from=0, to=maxage+1, length.out=101) # create pseudo sub-year classes
#   
#   Linfarray_c <- as.matrix(Linfarray[i,])
#   nyrs <- nrow(Linfarray_c)
#   VulnN <- as.matrix(vn[i,,]) * scaleR0[i]
#   VulnN <- round(VulnN,0)
#   if (nyrs == 1) VulnN <- t(VulnN)
#   # 
#   # out <- list(AgeVec=AgeVec, SubAgeVec=SubAgeVec,
#   #             Linfarray_c=Linfarray_c,
#   #             Karray_c=as.matrix(Karray[i,]),
#   #             t0array_c=as.matrix(t0array[i,]),
#   #             LenCV_c=LenCV[i],
#   #             CAL_bins, CAL_binsmid, retLength=as.matrix(retL[i,, ]), CAL_ESS=CAL_ESS[i],
#   #             CAL_nsamp=CAL_nsamp[i], VulnN=VulnN, truncSD=truncSD)
#   # saveRDS(out, 'out.rdata')
#  
#  makeLenComp(AgeVec, SubAgeVec,
#               Linfarray_c=Linfarray_c,
#               Karray_c=as.matrix(Karray[i,]),
#               t0array_c=as.matrix(t0array[i,]),
#               LenCV_c=LenCV[i],
#               CAL_bins, CAL_binsmid, retLength=as.matrix(retL[i,, ]), CAL_ESS=CAL_ESS[i],
#               CAL_nsamp=CAL_nsamp[i], VulnN=VulnN, truncSD)
#   
# }
# 
# 
# makeSizeCompW2 <- function(i, nyrs, maxage, Linfarray, Karray, t0array, LenCV,
#                           CAL_bins, CAL_binsmid, retL, CAL_ESS, CAL_nsamp,
#                           vn, truncSD=2) {
# 
# 
#   makeSizeComp(nyrs, maxage, Linfarray_c=as.matrix(Linfarray[i,]),
#                Karray_c=as.matrix(Karray[i,]),
#                t0array_c=as.matrix(t0array[i,]),
#                LenCV_c=LenCV[i],
#                CAL_bins, CAL_binsmid, retLength=as.matrix(retL[i,, ]), CAL_ess=CAL_ESS[i],
#                CAL_Nsamp=CAL_nsamp[i], VulnN=as.matrix(vn[i,,]), truncSD)
# 
# }
# 
# 
# makeSizeComp <- function(nyrs, maxage, Linfarray_c, Karray_c, t0array_c, LenCV_c,
#                          CAL_bins, CAL_binsmid, retLength, CAL_ess, CAL_Nsamp,
#                          VulnN, truncSD=2) {
# 
#   AgeVec <- seq(from=0, to=maxage+1, length.out=101) # create pseudo sub-year classes
# 
#   ageby <- AgeVec[2] - AgeVec[1]
#   ages <- seq(AgeVec[2]-0.5*ageby, by=ageby, length.out=length(AgeVec)-1)
# 
#   VulnN2 <- matrix(NA, nrow=nyrs, ncol=length(ages))
#   LenAge <- matrix(NA, nrow=length(ages), ncol=nyrs)
# 
#   tempfun <- function(x, VulnN, maxage, AgeVec) {
#     tempval <- rep(1:maxage, VulnN[x,]) + runif(sum(VulnN[x,]), -0.5, 0.5) # add variability to ages
#     hist(tempval, breaks=AgeVec, plot=FALSE)$counts
#   }
# 
#   if (nyrs > 1) {
#     temp <- lapply(1:nyrs, tempfun, VulnN=round(VulnN,0), maxage=maxage, AgeVec=AgeVec)
#     VulnN2 <- do.call("rbind", temp)
#   } else {
#     VulnN <- round(VulnN,0)
#     tempval <- rep(1:maxage, VulnN[,1]) + runif(sum(VulnN[,1]), -0.5, 0.5) # add variability to ages
#     VulnN2 <- matrix(hist(tempval, breaks=AgeVec, plot=FALSE)$counts, nrow=nyrs, ncol=length(AgeVec)-1)
#   }
# 
#   ind <- as.matrix(expand.grid(1:length(ages), 1:nyrs))  # an index for calculating Length at age
#   LenAge[ind] <- Linfarray_c[ind[,2]] * (1 - exp(-Karray_c[ind[, 2]] * (ages - t0array_c[ind[, 2]])))
#   LenAge[LenAge<0] <- 0
#   LenSD <- LenAge * LenCV_c
# 
#   genLenComp(CAL_bins, CAL_binsmid, retLength, CAL_ess, CAL_Nsamp,
#              VulnN2, LenAge, LenSD, truncSD)
# }
# 
# 






