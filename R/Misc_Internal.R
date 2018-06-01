utils::globalVariables(c("R0", "Mdb", "mod", "i", "nareas", "nsim", "dFfinal"))
tiny <- 1e-15  # define tiny variable
proportionMat <- TL <- Wa <- SurvWeiMat <- r <- lx <- logNormDensity <- sumlogNormDen <- NULL
proportionMat = vector()


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
  for (sl in slotNames(OM)) {
    slt <- slot(OM, sl)
    if (class(slt) == "numeric") {
      if (length(slt) ==2) {
        if (!all(is.na(slt)) && slt[1] == slt[2]) {
          slt[1] <- slt[1]
          slt[2] <- slt[2]+tiny
        }
       
      }
      if (length(slt)==1) {
        if (!all(is.na(slt)) && slt == 0) slt <- tiny
      }
      slot(OM,sl) <- slt
    }
  }
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


runMSEnomsg <- function(...) {
  capture.output(out <- suppressMessages(runMSE(...)))
  out
}


run_parallel <- function(i, itsim, OM, MPs, CheckMPs, timelimit, Hist, ntrials, fracD, CalcBlow, 
                         HZN, Bfrac, AnnualMSY, silent, PPD) {
  OM@nsim <- itsim[i]
  
  OM@seed <- OM@seed + i 
  mse <- runMSE_int(OM, MPs, CheckMPs, timelimit, Hist, ntrials, fracD, CalcBlow, 
                    HZN, Bfrac, AnnualMSY, silent, PPD)
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



# Generate size comps
genSizeCompWrap <- function(i, vn, CAL_binsmid, retL,
                            CAL_ESS, CAL_nsamp,
                            Linfarray, Karray, t0array,
                            LenCV, truncSD=2) {
  
  VulnN <- as.matrix(vn[i,,]) 
  VulnN <- round(VulnN,0)
  nyrs <- nrow(as.matrix(Linfarray[i,]))
  if (nyrs == 1) VulnN <- t(VulnN)

  lens <- genSizeComp(VulnN, CAL_binsmid, retL[i,,],
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

genSizeCompWrap2<- function(i, vn, CAL_binsmid, 
                            CAL_ESS, CAL_nsamp,
                            Linfarray, Karray, t0array,
                            LenCV, truncSD=2) {
  
  VulnN <- as.matrix(vn[i,,]) 
  VulnN <- round(VulnN,0)
  nyrs <- nrow(as.matrix(Linfarray[i,]))
  if (nyrs == 1) VulnN <- t(VulnN)
  
  
  lens <- genSizeComp2(VulnN, CAL_binsmid, 
                      CAL_ESS=CAL_ESS[i], CAL_nsamp=CAL_nsamp[i],
                      Linfs=Linfarray[i,], Ks=Karray[i,], t0s=t0array[i,],
                      LenCV=LenCV[i], truncSD)
  
  lens
  
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





