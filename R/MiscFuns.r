
# Collection of miscellaneous functions.  All functions have
# accompanying help files.

#' Setup parallel processing
#'
#' Sets up parallel processing using the snowfall package
#' @importFrom snowfall sfExport sfIsRunning sfSapply
#' @importFrom snow setDefaultClusterOptions
#' @export 
setup <- function() {
  snowfall::sfInit(parallel=TRUE,cpus=parallel::detectCores())  
  if (length(ls()) > 0) snowfall::sfExportAll() 
}

#' What objects of this class are available
#' 
#' Generic class finder
#' 
#' Finds objects of the specified class in the global environment or the
#' package:DLMtool
#' 
#' @usage avail(classy)
#' @param classy A class of object (character string, e.g. 'Fleet')
#' @author T. Carruthers
#' @export avail
avail <- function(classy) {
  return(unique(c(ls("package:DLMtool")[unlist(lapply(ls("package:DLMtool"), 
    getclass, classy = classy))], ls(envir = .GlobalEnv)[unlist(lapply(ls(envir = .GlobalEnv), 
    getclass, classy = classy))])))
}

#' get object class
#' 
#' Internal function for determining if object is of classy
#' 
#' 
#' @usage getclass(x, classy)
#' @param x Character string object name
#' @param classy A class of object (character string, e.g. 'Fleet')
#' @author T. Carruthers
#' @return TRUE or FALSE
#' @export getclass
getclass <- function(x, classy) inherits(get(x), classy)


#' What methods need what data
#' 
#' A function that finds all methods in the environment and searches the
#' function text for slots in the DLM data object
#' 
#' 
#' @usage Required(funcs = NA)
#' @param funcs A character vector of possible methods of class DLM quota, DLM
#' space or DLM size
#' @author T. Carruthers
#' @export Required
Required <- function(funcs = NA) {
  if (is.na(funcs[1])) 
    funcs <- c(avail("DLM_output"), avail("DLM_input"))
  slots <- slotNames("DLM_data")
  slotnams <- paste("DLM_data@", slotNames("DLM_data"), sep = "")
  repp <- rep("", length(funcs))
  
  for (i in 1:length(funcs)) {
    temp <- format(match.fun(funcs[i]))
    temp <- paste(temp[1:(length(temp))], collapse = " ")
    rec <- ""
    for (j in 1:length(slotnams)) if (grepl(slotnams[j], temp)) 
      rec <- c(rec, slots[j])
    if (length(rec) > 1) 
      repp[i] <- paste(rec[2:length(rec)], collapse = ", ")
  }
  cbind(funcs, repp, deparse.level = 0)
}

# A way of locating where the package was installed so you can find
# example data files and code etc.


#' Directory of the installed package on your computer
#' 
#' A way of locating where the package was installed so you can find example
#' data files and code etc.
#' 
#' 
#' @usage DLMDataDir(stock=NA)
#' @param stock Character string representing the name of a .csv file e.g.
#' 'Snapper', 'Rockfish'
#' @author T. Carruthers
#' @export DLMDataDir
DLMDataDir <- function(stock = NA) {
  if (is.na(stock)) {
    return(paste(searchpaths()[match("package:DLMtool", search())], 
      "/", sep = ""))
  } else {
    return(paste(searchpaths()[match("package:DLMtool", search())], 
      "/", stock, ".csv", sep = ""))
  }
}


#' MP feasibility diagnostic
#' 
#' What MPs may be run (best case scenario) for various data-availability
#' scenarios?
#' 
#' 
#' @usage Fease(feaseobj,outy='table')
#' @param feaseobj An object of class 'DLM_fease'
#' @param outy Determines whether you would like a full table or some column of
#' the table for a specific case of the feasibility object. When set equal to
#' table, the full table is produced. When set equal to an integer number the
#' names of MPs that are feasible for that case are returned.
#' @author T. Carruthers
#' @export Fease
Fease <- function(feaseobj, outy = "table") {
  
  if (class(feaseobj) != "DLM_fease") 
    stop("Incorrect format: you need an object of class DLM_fease")
  
  sloty <- c("Cat", "Ind", "AvC", "Dt", "Rec", "CAA", "CAL", "Mort", 
    "L50", "L95", "vbK", "vbLinf", "vbt0", "wla", "wlb", "steep", "LFC", 
    "LFS", "Cref", "Bref", "Iref", "Dep", "Abun", "ML")
  
  type <- c("Catch", "Index", "Catch", "Index", "Recruitment_index", 
    "Catch_at_age", "Catch_at_length", "Natural_mortality_rate", "Maturity_at_length", 
    "Maturity_at_length", "Growth", "Growth", "Growth", "Length_weight_conversion", 
    "Length_weight_conversion", "Stock_recruitment_relationship", "Fleet_selectivity", 
    "Fleet_selectivity", "Target_catch", "Target_biomass", "Target_index", 
    "Index", "Abundance")
  
  ncases <- length(feaseobj@Case)
  slots <- slotNames(feaseobj)
  ns <- length(slots)
  ftab <- array(TRUE, c(ns - 2, ncases))
  for (j in 3:ns) ftab[j - 2, ] <- as.logical(as.numeric(slot(feaseobj, 
    slots[j])))
  
  req <- Required()
  nMPs <- nrow(req)
  gridy <- array("", c(nMPs, ncases))
  for (i in 1:ncases) {
    types <- slotNames(feaseobj)[3:17][ftab[, i]]
    slots <- sloty[type %in% types]
    for (m in 1:nMPs) {
      brec <- unlist(strsplit(req[m, 2], ", "))
      brec <- brec[grep("CV_", brec, invert = T)]  #remove CV dependencies (we think we can guess these...)
      brec <- brec[brec != "Year" & brec != "MaxAge" & brec != "FMSY_M" & 
        brec != "BMSY_B0" & brec != "t" & brec != "OM" & brec != 
        "MPrec" & brec != "CAL_bins" & brec != "MPeff" & brec != 
        "LHYear"]
      nr <- length(brec)
      if (nr == 0) {
        gridy[m, i] <- "Yes"
      } else {
        cc <- 0
        for (r in 1:nr) {
          # loop over requirements
          if (brec[r] %in% slots) 
          cc <- cc + 1
        }
        if (cc == nr) 
          gridy[m, i] <- "Yes"
      }
    }
  }
  gridy <- as.data.frame(gridy)
  row.names(gridy) = req[, 1]
  names(gridy) = feaseobj@Case
  if (outy == "table") 
    return(gridy)
  if (outy != "table" & class(outy) != "numeric") 
    return(req[, 1][gridy[, 1] == "Yes"])
  if (class(outy) == "numeric") {
    if (outy < (ncases + 1)) {
      return(req[, 1][gridy[, as.integer(outy)] == "Yes"])
    } else {
      return(req[, 1][gridy[, 1] == "Yes"])
    }
  }
  
}



#' Optimization function to find a movement model that matches user specified
#' movement characteristics.
#' 
#' The user specifies the probability of staying in the same area and spatial
#' heterogeneity (both in the unfished state).
#' 
#' This is paired with movfit to find the correct movement model.
#' 
#' @usage getmov(x,Prob_staying,Frac_area_1)
#' @param x A position in vectors Prob_staying and Frac_area_1
#' @param Prob_staying User specified probability that individuals in area 1
#' remain in that area (unfished conditions)
#' @param Frac_area_1 User specified fraction of individuals found in area 1
#' (unfished conditions)
#' @return A markov movement matrix
#' @author T. Carruthers
#' @examples
#' 
#' Prob_staying<-0.8 # probability  that individuals remain in area 1 between time-steps
#' Frac_area_1<-0.35 # the fraction of the stock found in area 1 under equilibrium conditions
#' markovmat<-getmov(1,Prob_staying, Frac_area_1)
#' vec<-c(0.5,0.5) # initial guess at equilibrium distribution (2 areas)
#' for(i in 1:300)vec<-apply(vec*markovmat,2,sum) # numerical approximation to stable distribution
#' c(markovmat[1,1],vec[1]) # pretty close right?
#' 
#' 
#' @export getmov
getmov <- function(x, Prob_staying, Frac_area_1) {
  test <- optim(par = c(0, 0, 0), movfit, method = "L-BFGS-B", 
    lower = rep(-6, 3), upper = rep(6, 3), prb = Prob_staying[x], 
	frac = Frac_area_1[x])
  mov <- array(c(test$par[1], test$par[2], 0, test$par[3]), dim = c(2, 2))
  mov <- exp(mov)
  mov/array(apply(mov, 1, sum), dim = c(2, 2))
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
#' @examples
#' 
#' Prob_staying<-0.8 # probability  that individuals remain in area 1 between time-steps
#' Frac_area_1<-0.35 # the fraction of the stock found in area 1 under equilibrium conditions
#' markovmat<-getmov(1,Prob_staying, Frac_area_1)
#' vec<-c(0.5,0.5) # initial guess at equilibrium distribution (2 areas)
#' for(i in 1:300)vec<-apply(vec*markovmat,2,sum) # numerical approximation to stable distribution
#' c(markovmat[1,1],vec[1]) # pretty close right?
#' 
#' 
#' @export getmov2
getmov2 <- function(x, Prob_staying, Frac_area_1) {
  test <- optim(par = c(0, 0, 0), movfit_Rcpp, method = "L-BFGS-B", 
    lower = rep(-6, 3), upper = rep(6, 3), prb = Prob_staying[x], 
	frac = Frac_area_1[x])	
  mov <- array(c(test$par[1], test$par[2], 0, test$par[3]), dim = c(2, 2))
  mov <- exp(mov)
  mov/array(apply(mov, 1, sum), dim = c(2, 2))
}

#' Optimization function that returns the squared difference between user
#' specified and calculated movement parameters. (deprecated: now in Rcpp)
#' 
#' The user specifies the probability of staying in the same area and spatial
#' heterogeneity (both in the unfished state). This function returns the
#' squared difference between these values and those produced by the three
#' logit movement model.
#' 
#' This is paired with getmov to find the correct movement model.
#' 
#' @usage movfit(par,prb,frac)
#' @param par Three parameters in the logit space that control the four
#' probabilities of moving between 2 areas
#' @param prb User specified probability that individuals in area 1 remain in
#' that area (unfished conditions)
#' @param frac User specified fraction of individuals found in area 1 (unfished
#' conditions)
#' @author T. Carruthers
#' @export movfit
movfit <- function(par, prb, frac) {
  # .Deprecated("movfit_Rcpp")
  mov <- array(c(par[1], par[2], 0, par[3]), dim = c(2, 2))
  mov <- exp(mov)
  mov <- mov/array(apply(mov, 1, sum), dim = c(2, 2))
  dis <- c(frac, 1 - frac)
  for (i in 1:100) dis <-  apply(array(dis, c(2, 2)) * mov, 2, sum)
  (log(mov[1, 1]) - log(prb))^2 + (log(frac) - log(dis[1]))^2
}


#' Optimization function that find the catchability (q where F=qE) value
#' required to get to user-specified stock depletion (current biomass /
#' unfished biomass)
#' 
#' The user specifies the level of stock depleiton. This funciton takes the
#' derived effort trajectories and finds the catchabiltiy to get the stock
#' there.
#'
#' @param x internal parameter
#' @param dep internal parameter
#' @param Find internal parameter
#' @param Perr internal parameter
#' @param Marray internal parameter
#' @param hs internal parameter
#' @param Mat_age internal parameter
#' @param Wt_age internal parameter
#' @param R0 internal parameter
#' @param V internal parameter
#' @param nyears internal parameter
#' @param maxage internal parameter
#' @param mov internal parameter
#' @param Spat_targ internal parameter
#' @param SRrel internal parameter
#' @param aR internal parameter
#' @param bR internal parameter
#' 
#' Paired with qopt
#' @keywords internal
#' @export getq 
#'
#' @author T. Carruthers
getq <- function(x, dep, Find, Perr, Marray, hs, Mat_age, Wt_age, R0, V, 
  nyears, maxage, mov, Spat_targ, SRrel, aR, bR) {
  opt <- optimize(qopt, log(c(0.0075, 15)), depc = dep[x], Fc = Find[x, ], 
    Perrc = Perr[x, ], Mc = Marray[x, ], hc = hs[x], Mac = Mat_age[x, ], 
    Wac = Wt_age[x, , ], R0c = R0[x], Vc = V[x, , ], nyears = nyears, 
	maxage = maxage, movc = mov[x, , ], Spat_targc = Spat_targ[x], 
    SRrelc = SRrel[x], aRc = aR[x, ], bRc = bR[x, ])	
  return(exp(opt$minimum))
}

#' Optimization function that find the catchability (q where F=qE) value
#' required to get to user-specified stock depletion (current biomass /
#' unfished biomass) modified for Rcpp 
#' 
#' The user specifies the level of stock depleiton. This funciton takes the
#' derived effort trajectories and finds the catchabiltiy to get the stock
#' there.
#'
#' @param x internal parameter
#' @param dep internal parameter
#' @param Find internal parameter
#' @param Perr internal parameter
#' @param Marray internal parameter
#' @param hs internal parameter
#' @param Mat_age internal parameter
#' @param Wt_age internal parameter
#' @param R0 internal parameter
#' @param V internal parameter
#' @param nyears internal parameter
#' @param maxage internal parameter
#' @param mov internal parameter
#' @param Spat_targ internal parameter
#' @param SRrel internal parameter
#' @param aR internal parameter
#' @param bR internal parameter
#' 
#' Paired with qopt
#' @keywords internal
#' @export getq2 
#'
#' @author T. Carruthers
getq2 <- function(x, dep, Find, Perr, Marray, hs, Mat_age, Wt_age, R0, V, 
  nyears, maxage, mov, Spat_targ, SRrel, aR, bR) {
  opt <- optimize(projOpt_cpp, log(c(0.0075, 15)), depc = dep[x], Fc = Find[x, ], 
    Perrc = Perr[x, ], Mc = Marray[x, ], hc = hs[x], Mac = Mat_age[x, ], 
    Wac = Wt_age[x, , ], R0c = R0[x], Vc = V[x, , ], nyears = nyears, 
	maxage = maxage, movc = mov[x, , ], Spat_targc = Spat_targ[x], 
    SRrelc = SRrel[x], aRc = aR[x, ], bRc = bR[x, ], proyears=0, FMSY=0, Control=1)	
  return(exp(opt$minimum))
}



#' Internal optimization function that find the catchability (q where F=qE)
#' value required to get to user-specified stock depletion (current biomass /
#' unfished biomass)
#' 
#' The user specifies the level of stock depleiton. This funciton takes the
#' derived effort trajectories and finds the catchabiltiy to get the stock
#' there.
#' 
#' @param lnq internal parameter
#' @param depc internal parameter
#' @param Fc internal parameter
#' @param Perrc internal parameter
#' @param Mc internal parameter
#' @param hc internal parameter
#' @param Mac internal parameter
#' @param Wac internal parameter
#' @param R0c internal parameter
#' @param Vc internal parameter
#' @param nyears internal parameter
#' @param maxage internal parameter
#' @param movc internal parameter
#' @param Spat_targc internal parameter
#' @param SRrelc internal parameter
#' @param aRc internal parameter
#' @param bRc internal parameter
#' @param opt internal parameter
#' 
#' @export qopt 
#' @keywords internal
#' @author T. Carruthers 
qopt <- function(lnq, depc, Fc, Perrc, Mc, hc, Mac, Wac, R0c, Vc, nyears, 
  maxage, movc, Spat_targc, SRrelc, aRc, bRc, opt = T) {
  qc <- exp(lnq)
  nareas <- nrow(movc)
  # areasize<-c(asizec,1-asizec)
  idist <- rep(1/nareas, nareas)
  for (i in 1:300) idist <- apply(array(idist, c(2, 2)) * movc, 2, sum)
 
  N <- array(exp(-Mc[1] * ((1:maxage) - 1)) * R0c, dim = c(maxage, nareas)) * 
    array(rep(idist, each = maxage), dim = c(maxage, nareas))
  SSN <- Mac * N  # Calculate initial spawning stock numbers
  Biomass <- N * Wac[, 1]
  SSB <- SSN * Wac[, 1]  # Calculate spawning stock biomass
  
  B0 <- sum(Biomass)
  R0a <- idist * R0c
  SSB0 <- apply(SSB, 2, sum)
  SSBpR <- SSB0/R0a  # Calculate spawning stock biomass per recruit
  
  for (y in 1:nyears) {
    # set up some indices for indexed calculation
    targ <- (apply(Vc[, y] * Biomass, 2, sum)^Spat_targc)/mean(apply(Vc[, y] * 
	  Biomass, 2, sum)^Spat_targc)
    FMc <- array(qc * Fc[y] * Vc[, y], dim = c(maxage, nareas)) * array(rep(targ, each = maxage), 
	  dim = c(maxage, nareas))  # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
    Zc <- FMc + Mc[y]

    N[2:maxage, ] <- N[1:(maxage - 1), ] * exp(-Zc[1:(maxage - 1), ])  # Total mortality
    if (SRrelc == 1) {
      N[1, ] <- Perrc[y] * (0.8 * R0a * hc * apply(SSB, 2, sum))/(0.2 * 
        SSBpR * R0a * (1 - hc) + (hc - 0.2) * apply(SSB, 2, sum))  # Recruitment assuming regional R0 and stock wide steepness
    } else {
      N[1, ] <- aRc * apply(SSB, 2, sum) * exp(-bRc * apply(SSB, 2, sum)) # No error?
    }
	  
    # print(N[1])
    indMov <- as.matrix(expand.grid(1:nareas, 1:nareas, 1:maxage)[3:1])
    indMov2 <- indMov[, 1:2]
    indMov3 <- indMov[, 2:3]
    temp <- array(N[indMov2] * movc[indMov3], dim = c(nareas, nareas, maxage))
		
    N <- apply(temp, c(3, 1), sum)
    SSN <- N * Mac
    SSB <- SSN * Wac[, y]
    Biomass <- N * Wac[, y]
    SBiomass <- SSN * Wac[, y]
    # print(sum(Biomass))
  }  # end of year
  
return((log(depc) - log(sum(SBiomass)/sum(SSB0)))^2)
}



	
## Operating Model Functions
## --------------------------------------------------- These functions
## are used to manually specify, choose, or estimate various parameters
## of the Operating Model.  The functions typically take OM object (or
## Stock or Fleet) and return the same object with the relevant
## parameters populated.

#' Depletion and F estimation from mean length of catches
#' 
#' A highly dubious means of getting very uncertain estimates of current stock
#' biomass and (equilibrium) fishing mortality rate from growth, natural
#' mortality rate, recruitment and fishing selectivity.
#' 
#' 
#' @usage ML2D(OM,ML,nsim=100,ploty=T,Dlim=c(0.05,0.6))
#' @param OM An object of class 'OM'
#' @param ML A estimate of current mean length of catches
#' @param nsim Number of simulations
#' @param ploty Produce a plot of depletion and F
#' @param Dlim Limits on the depletion that is returned as a fraction of
#' unfished biomass.
#' @return A table of nsim rows and 2 columns (depletion, fishing mortality
#' rate)
#' @author T. Carruthers
#' @export ML2D
ML2D <- function(OM, ML, nsim = 100, ploty = T, Dlim = c(0.05, 0.6)) {
  
  maxage <- OM@maxage
  M <- runif(nsim, OM@M[1], OM@M[2])  # Natural mortality rate
  h <- runif(nsim, OM@h[1], OM@h[2])  # Steepness
  Linf <- runif(nsim, OM@Linf[1], OM@Linf[2])  # Maximum length
  K <- runif(nsim, OM@K[1], OM@K[2])  # Maximum growth rate
  t0 <- runif(nsim, OM@t0[1], OM@t0[2])  # Theorectical length at age zero
  
  if (OM@isRel == "0" | OM@isRel == "FALSE" | OM@isRel == FALSE) {
    if (max(OM@LFS) > 0) {
      LFS <- runif(nsim, OM@LFS[1], OM@LFS[2])
    } else {
      LFS <- runif(nsim, mean(OM@LFSLower), mean(OM@LFSUpper))
    }
  } else {
    if (max(OM@LFS) > 0) {
      LFS <- runif(nsim, OM@LFS[1], OM@LFS[2]) * mean(OM@L50)
    } else {
      LFS <- runif(nsim, mean(OM@LFSLower), mean(OM@LFSUpper)) * 
        mean(OM@L50)
    }
  }
  AFS <- L2A(t0, Linf, K, LFS, maxage)
  
  L5 <- runif(nsim, OM@L5[1], OM@L5[2]) * mean(OM@L50)
  age05 <- L2A(t0, Linf, K, L5, maxage)
  
  Vmaxage <- runif(nsim, OM@Vmaxlen[1], OM@Vmaxlen[2])  #runif(BT_fleet@Vmaxage[1],BT_fleet@Vmaxage[2]) # selectivity of oldest age class
  
  LM <- runif(nsim, OM@L50[1], OM@L50[2])
  AM <- L2A(t0, Linf, K, LM, maxage)
  
  # age at maturity
  a <- OM@a  # length-weight parameter a
  b <- OM@b  # length-weight parameter b
  
  mod <- AFS  # the age at modal (or youngest max) selectivity
  deriv <- getDNvulnS(mod, age05, Vmaxage, maxage, nsim)  # The vulnerability schedule
  vuln <- deriv[[1]]
  
  Agearray <- array(rep(1:maxage, each = nsim), c(nsim, maxage))
  mat <- 1/(1 + exp((AM - (Agearray))/(AM * 0.1)))  # Maturity at age array
  
  nyears <- 100
  # bootfun<-function(dat,ind)mean(dat[ind]) MLo<-boot(MLt,bootfun,nsim)
  # ML<-MLo$t
  out <- CSRA(M, h, Linf, K, t0, AM, a, b, vuln, mat, ML = rep(ML, nsim), 
    NA, NA, maxage, nyears)
  cond <- out[, 1] > Dlim[1] & out[, 1] < Dlim[2] & out[, 2] < 2.5  # Stock levels are unlikely to be above 80% unfished, F is unlikely to be above 2.5
  
  if (ploty & sum(cond) > 0) {
    par(mfrow = c(1, 2))
    plot(density(out[cond, 1], from = 0, adj = 0.4), main = "Depletion")
    plot(density(out[cond, 2], from = 0, adj = 0.4), main = "Fishing mortality rate")
    OM@D <- quantile(out[cond, 1], c(0.05, 0.95))
  }
  if (sum(cond) == 0) {
    message("All estimates of Depletion outside bounds of Dlim")
    message("Operating Model object not updated")
  }
  if (!ploty) 
    message("Operating Model object not updated")
  OM
}

# Composition stock reduction analysis


#' Catch at size reduction analysis
#' 
#' What depletion level and corresponding equlibrium F arise from data
#' regarding mean length of current catches, natural mortality rate, steepness
#' of the stock recruitment curve, maximum length, maximum growth rate, age at
#' maturity, age based vulnerability, maturity at age, maximum age and number
#' of historical years of fishing.
#' 
#' 
#' @usage CSRA(M,h,Linf,K,t0,AM,a,b,vuln,mat,ML,CAL,CAA,maxage,nyears)
#' @param M A vector of natural mortality rate estimates
#' @param h A vector of sampled steepness (Beverton-Holt stock recruitment)
#' @param Linf A vector of maximum length (von Bertalanffy growth)
#' @param K A vector of maximum growth rate (von Bertalanffy growth)
#' @param t0 A vector of theoretical age at length zero (von Bertalanffy
#' growth)
#' @param AM A vector of age at maturity
#' @param a Length-weight conversion parameter a (W=aL^b)
#' @param b Length-weight conversion parameter b (W=aL^b)
#' @param vuln A matrix nsim x nage of the vulnerabilty at age (max 1) to
#' fishing.
#' @param mat A matrix nsim x nage of the maturity at age (max 1)
#' @param ML A vector of current mean length estimates
#' @param CAL A catch-at-length matrix nyears x (1 Linf unit) length bins
#' @param CAA A catch-at-age matrix nyears x maximum age
#' @param maxage Maximum age
#' @param nyears Number of historical years of fishing
#' @author T. Carruthers
#' @export CSRA
CSRA <- function(M, h, Linf, K, t0, AM, a, b, vuln, mat, ML, CAL, CAA, 
  maxage, nyears) {
  nsim <- length(M)
  Dep <- rep(NA, nsim)
  Fm <- rep(NA, nsim)
  for (i in 1:nsim) {
    fit <- optimize(CSRAfunc, log(c(1e-04, 5)), Mc = M[i], hc = h[i], 
      maxage, nyears, Linfc = Linf[i], Kc = K[i], t0c = t0[i], AMc = AM[i], 
      ac = a, bc = b, vulnc = vuln[i, ], matc = mat[i, ], MLc = ML[i], 
      CAL = NA, CAA = NA, opt = T)
    
    
    out <- CSRAfunc(fit$minimum, Mc = M[i], hc = h[i], maxage, nyears, 
      Linfc = Linf[i], Kc = K[i], t0c = t0[i], AMc = AM[i], ac = a, 
      bc = b, vulnc = vuln[i, ], matc = mat[i, ], MLc = ML[i], CAL = NA, 
      CAA = NA, opt = 3)
    
    Dep[i] <- out[1]
    Fm[i] <- out[2]
    
    
  }
  cbind(Dep, Fm)
}

# The function that CSRA operates on

#' Optimization function for CSRA
#' 
#' What depletion level and corresponding equlibrium F arise from data
#' regarding mean length of current catches, natural mortality rate, steepness
#' of the stock recruitment curve, maximum length, maximum growth rate, age at
#' maturity, age based vulnerability, maturity at age, maximum age and number
#' of historical years of fishing.
#' 
#' 
#' @usage CSRAfunc(lnF,Mc,hc,maxage,nyears,AFSc,AFCc,Linfc,
#' Kc,t0c,AMc,ac,bc,vulnc,matc,MLc,CAL,CAA, opt=T,meth='ML')
#' @param lnF A proposed value of current instantaneous fishing mortality rate
#' @param Mc Natural mortality rate estimates
#' @param hc Steepness (Beverton-Holt stock recruitment)
#' @param maxage Maximum age
#' @param nyears Number of historical years of fishing
#' @param AFSc Age at full selection
#' @param AFCc Age at first capture
#' @param Linfc Maximum length (von Bertalanffy growth)
#' @param Kc Maximum growth rate (von Bertalanffy growth)
#' @param t0c Theoretical age at length zero (von Bertalanffy growth)
#' @param AMc Age at maturity
#' @param ac Length-weight conversion parameter a (W=aL^b)
#' @param bc Length-weight conversion parameter b (W=aL^b)
#' @param vulnc A vector (nage long) of the vulnerabilty at age (max 1) to
#' fishing.
#' @param matc A vector (nage long) of the maturity at age (max 1)
#' @param MLc A current mean length estimates
#' @param CAL A catch-at-length matrix nyears x (1 Linf unit) length bins
#' @param CAA A catch-at-age matrix nyears x maximum age
#' @param opt Should the measure of fit be returned?
#' @param meth Are we fitting to mean length or catch composition?
#' @author T. Carruthers
#' @export CSRAfunc
CSRAfunc <- function(lnF, Mc, hc, maxage, nyears, AFSc, AFCc, Linfc, Kc, 
  t0c, AMc, ac, bc, vulnc, matc, MLc, CAL, CAA, opt = T, meth = "ML") {
  
  Fm <- exp(lnF)
  Fc <- vulnc * Fm
  Lac <- Linfc * (1 - exp(-Kc * ((1:maxage) - t0c)))
  Wac <- ac * Lac^bc
  N <- exp(-Mc * ((1:maxage) - 1))
  SSN <- matc * N  # Calculate initial spawning stock numbers
  Biomass <- N * Wac
  SSB <- SSN * Wac  # Calculate spawning stock biomass
  
  B0 <- sum(Biomass)
  SSB0 <- sum(SSB)
  SSN0 <- SSN
  SSBpR <- sum(SSB)  # Calculate spawning stock biomass per recruit
  SSNpR <- SSN
  Zc <- Fc + Mc
  CN <- array(NA, dim = c(nyears, maxage))
  HR <- rep(0, maxage)
  pen <- 0
  for (y in 1:nyears) {
    VB <- Biomass * vulnc * exp(-Mc)
    CN[y, ] <- N * (1 - exp(-Zc)) * (Fc/Zc)
    N[2:maxage] <- N[1:(maxage - 1)] * exp(-Zc[1:(maxage - 1)])  # Total mortality
    N[1] <- (0.8 * hc * sum(SSB))/(0.2 * SSBpR * (1 - hc) + (hc - 0.2) * 
      sum(SSB))  # Recruitment assuming regional R0 and stock wide steepness
    Biomass <- N * Wac
    SSN <- N * matc
    SSB <- SSN * Wac
  }  # end of year
  
  pred <- sum((CN[nyears, ] * Lac))/sum(CN[nyears, ])
  fobj <- (pred - MLc)^2  # Currently a least squares estimator. Probably not worth splitting hairs WRT likelihood functions!
  if (opt == 1) {
    return(fobj)
  } else {
    c(sum(SSB)/sum(SSB0), Fm)
  }
}

# Stochastic inverse growth curve used to back-calculate age at first
# capture from length at first capture
#' Calculate age at first capture from length at first capture and growth
#' 
#' As title.
#' 
#' 
#' @usage getAFC(t0c,Linfc,Kc,LFC,maxage)
#' @param t0c A vector of theoretical age at length zero (von Bertalanffy
#' growth)
#' @param Linfc A vector of maximum length (von Bertalanffy growth)
#' @param Kc A vector of maximum growth rate (von Bertalanffy growth)
#' @param LFC A vector of length at first capture
#' @param maxage Maximum age
#' @author T. Carruthers
#' @export getAFC
getAFC <- function(t0c, Linfc, Kc, LFC, maxage) {
  nsim <- length(t0c)
  agev <- c(1e-04, 1:maxage)
  agearray <- matrix(rep(agev, each = nsim), nrow = nsim)
  Larray <- Linfc * (1 - exp(-Kc * (agearray - t0c)))
  matplot(agev, t(Larray), type = "l")
  abline(h = LFC, col = "#ff000030", lwd = 2)
  AFC <- (log(1 - (LFC/Linfc))/-Kc) + t0c
  abline(v = AFC, col = "#0000ff30", lwd = 2)
  AFC
}



#' Length to age conversion
#' 
#' Simple deterministic length to age conversion given inverse von Bertalanffy
#' growth.
#' 
#' 
#' @usage L2A(t0c,Linfc,Kc,Len,maxage)
#' @param t0c Theoretical age at length zero
#' @param Linfc Maximum length
#' @param Kc Maximum growth rate
#' @param Len Length
#' @param maxage Maximum age
#' @return An age (vector of ages, matrix of ages) corresponding with Len
#' @author T. Carruthers
#' @export L2A
L2A <- function(t0c, Linfc, Kc, Len, maxage) {
  nsim <- length(t0c)
  agev <- c(1e-04, 1:maxage)
  agearray <- matrix(rep(agev, each = nsim), nrow = nsim)
  Larray <- Linfc * (1 - exp(-Kc * (agearray - t0c)))
  matplot(agev, t(Larray), type = "l")
  abline(h = Len, col = "#ff000030", lwd = 2)
  age <- (log(1 - (Len/Linfc))/-Kc) + t0c
  abline(v = age, col = "#0000ff30", lwd = 2)
  age
}


#' Function to calculate cyclic recruitment pattern given user-specified values
#' of period and amplitude.
#' 
#' Calculates cyclic pattern in recruitment deviations for a simulation. Ranges
#' for Period and Amplitude are specified by user, and function produces cyclic
#' pattern from within these ranges. Default is a sine wave.
#' 
#' 
#' @usage SetRecruitCycle(x=1, Period, Amplitude, TotYears, Shape=c('sin',
#' 'shift'))
#' @param x Simulation number.
#' @param Period A vector of length 2 specifying the minimum and maximum values
#' for the period of the recruitment cycles. e.g., if Period = c(10,10), then
#' recruitment cycle occurs every 10 years exactly.
#' @param Amplitude A vector of length 2 specifying the minimum and maximum
#' values for the amplitude of the recruitment cycles. e.g., if Amplitude =
#' c(0,0.5), the average recruitment will increase (or decrease) by a factor
#' between 0 and 0.5 each cycle.
#' @param TotYears A numeric value specifying the total number of years (should
#' be nyears + proyears).
#' @param Shape Specifies whether cyclic recruitment pattern is sine wave
#' (default) or a step-change (shift).
#' @author A. Hordyk
#' @export SetRecruitCycle
SetRecruitCycle <- function(x = 1, Period, Amplitude, TotYears, Shape = c("sin", 
  "shift")) {
  Shape <- match.arg(Shape)
  Npers <- ceiling(TotYears/min(Period))
  pers <- round(runif(Npers, min(Period), max(Period)), 0)
  amp <- runif(Npers, min(Amplitude), max(Amplitude))
  ct <- 1
  Dir <- sample(c(-1, 1), 1)
  Rm <- rep(NA, Npers)
  if (Shape == "sin") {
    for (X in 1:length(pers)) {
      yrper <- pers[X]
      Period <- pi/yrper  #round(runif(1,-1, 1),0) * pi/yrper 
      xs <- seq(from = 0, by = 1, to = pers[X])
      Dir <- ifelse(Dir >= 0, -1, 1)  # change direction each cycle
      Rm[ct:(ct + pers[X])] <- amp[X] * sin(xs * Period) * Dir
      ct <- ct + pers[X]
    }
  }
  if (Shape == "shift") {
    for (X in 1:length(pers)) {
      Dir <- ifelse(Dir >= 0, -1, 1)  # change direction each cycle
      if (X == 1) 
        Rm[ct:(ct + pers[X])] <- 0
      if (X > 1) 
        Rm[ct:(ct + pers[X])] <- amp[X] * Dir
      Rm[ct:(ct + pers[X])] <- amp[X] * Dir
      ct <- ct + pers[X]
    }
  }
  Rm <- Rm[1:TotYears] + 1
  return(Rm)
}






#' Convert a OM object to one without observation or process error
#' 
#' Takes an existing OM object and converts it to one without any observation
#' error, and very little process error.  Used for debugging and testing that
#' MPs perform as expected under perfect conditions.
#' 
#' 
#' @usage makePerf(OMin, except = NULL)
#' @param OMin An object of class \code{OM}
#' @param except An optional vector of slot names in the OM that will not be
#' changed (not tested perfectly so watch out!)
#' @return A new \code{OM} object
#' @author A. Hordyk
#' @export makePerf
makePerf <- function(OMin, except = NULL) {
  nms <- slotNames(OMin)
  # exceptions
  if (is.null(except)) 
    except <- "EVERYTHING"
  exclude <- unique(grep(paste(except, collapse = "|"), nms, value = FALSE))
  
  vars <- c("grad", "cv", "sd", "inc")
  ind <- unique(grep(paste(vars, collapse = "|"), nms, value = FALSE))
  ind <- ind[(!(nms[ind] %in% exclude))]
  for (X in seq_along(ind)) {
    n <- length(slot(OMin, nms[ind[X]]))
    slot(OMin, nms[ind[X]]) <- rep(0, n)
  }
  
  if (!("Cobs" %in% exclude)) 
    OMin@Cobs <- c(0, 0)
  if (!("Perr" %in% exclude)) 
    OMin@Perr <- c(0, 0)
  if (!("Iobs" %in% exclude)) 
    OMin@Iobs <- c(0, 0)
  if (!("AC" %in% exclude)) 
    OMin@AC <- c(0, 0)
  if (!("Btbias" %in% exclude)) 
    OMin@Btbias <- c(1, 1)
  if (!("CAA_ESS" %in% exclude)) 
    OMin@CAA_ESS <- c(1000, 1000)
  if (!("CAA_nsamp" %in% exclude)) 
    OMin@CAA_nsamp <- c(2000, 2000)
  if (!("CAL_ESS" %in% exclude)) 
    OMin@CAL_ESS <- c(1000, 1000)
  if (!("CAL_nsamp" %in% exclude)) 
    OMin@CAL_nsamp <- c(2000, 2000)
  if (!("beta" %in% exclude)) 
    OMin@beta <- c(1, 1)
  return(OMin)
}



#' Print out plotting functions
#' 
#' This function prints out the available plotting functions for objects of
#' class MSE or DLM_data
#' 
#' 
#' @usage plotFun(class = c('MSE', 'DLM_data'), msg=TRUE)
#' @param class Character string. Prints out the plotting functions for objects
#' of this class.
#' @param msg Logical. Should the functions be printed to screen?
#' @note Basically the function looks for any functions in the DLMtool that
#' have the word `plot` in them.  There is a chance that some plotting
#' functions are missed. Let us know if you find any and we will add them.
#' @author A. Hordyk
#' @export plotFun
plotFun <- function(class = c("MSE", "DLM_data"), msg = TRUE) {
  class <- match.arg(class)
  tt <- lsf.str("package:DLMtool")
  p <- p2 <- rep(FALSE, length(tt))
  for (X in seq_along(tt)) {
    temp <- grep("plot", tolower(tt[[X]]))
    if (length(temp) > 0) 
      p[X] <- TRUE
    temp2 <- grep(class, paste(format(match.fun(tt[[X]])), collapse = " "))
    if (length(temp2) > 0) 
      p2[X] <- TRUE
  }
  if (msg) 
    message("DLMtool functions for plotting objects of class ", class, 
      " are:")
  out <- sort(tt[which(p & p2)])
  out <- out[-grep("plotFun", out)]
  if (class == "MSE") {
    out <- c(out, "barplot", "boxplot", "VOI", "VOI2")
    out <- sort(out)
  }
  if (class == "DLM_data") {
    out <- c(out, "boxplot", "Sense")
    out <- sort(out)
  }
  if (msg) {
    if (length(out) > 5) {
      sq <- seq(from = 1, to = length(out), by = 5)
      for (x in seq_along(sq)) {
        cat(out[sq[x]:min(sq[x + 1] - 1, length(out), na.rm = TRUE)])
        cat("\n")
      }
    } else cat(out)
    cat("\n")
  }
  invisible(out)
}


### Various Small Helper Functions ###


#' Is a value NA or zero.
#' 
#' As title
#' 
#' 
#' @usage NAor0(x)
#' @param x A numeric value.
#' @return TRUE or FALSE 
#' @author T. Carruthers
#' @export NAor0
NAor0 <- function(x) {
    if (length(x) == 0) 
        return(TRUE)
    if (length(x) > 0) 
        return(is.na(x[1]))
}




#' Calculate CV from vector of values 
#' 
#' 
#' @usage cv(x)
#' @param x vector of numeric values 
#' @author T. Carruthers
#' @return numeric
#' @export cv
cv <- function(x) sd(x)/mean(x)


#' Get log normal standard deviation from transformed space mean and standard deviation 
#' 
#' @usage sdconv(m, sd)
#' @param m mean 
#' @param sd standard deviation
#' @author T. Carruthers
#' @return numeric
#' @export sdconv
sdconv <- function(m, sd) (log(1 + ((sd^2)/(m^2))))^0.5

#' Get log normal mean from transformed space mean and standard deviation
#' 
#' @usage mconv(m, sd)
#' @param m mean 
#' @param sd standard deviation
#' @author T. Carruthers
#' @return numeric
#' @export mconv
mconv <- function(m, sd) log(m) - 0.5 * log(1 + ((sd^2)/(m^2)))

#' Calculate alpha parameter for beta distribution from mean and standard deviation 
#' 
#' @usage alphaconv(m, sd)
#' @param m mean 
#' @param sd standard deviation
#' @author T. Carruthers
#' @return numeric
#' @export alphaconv
alphaconv <- function(m, sd) m * (((m * (1 - m))/(sd^2)) - 1)

#' Calculate beta parameter for beta distribution from mean and standard deviation 
#' 
#' @usage betaconv(m, sd)
#' @param m mean 
#' @param sd standard deviation
#' @author T. Carruthers
#' @return numeric
#' @export betaconv
betaconv <- function(m, sd) (1 - m) * (((m * (1 - m))/(sd^2)) - 1)

#'  Generate log-normally distributed random numbers 
#' 
#' @usage trlnorm(reps, mu, cv)
#' @param reps number of random numbers 
#' @param mu mean 
#' @param cv coefficient of variation
#' @author T. Carruthers
#' @return numeric
#' @export trlnorm
trlnorm <- function(reps, mu, cv) {
    if (all(is.na(mu))) return(rep(NA, reps))
    if (all(is.na(cv))) return(rep(NA, reps))
    if (reps == 1)  return(mu)
    return(rlnorm(reps, mconv(mu, mu * cv), sdconv(mu, mu * cv)))
}


#'  Calculate density of log-normally distributed random numbers 
#' 
#' @usage tdlnorm(x, mu, cv)
#' @param x vector 
#' @param mu mean 
#' @param cv coefficient of variation
#' @author T. Carruthers
#' @return numeric
#' @export tdlnorm
tdlnorm <- function(x, mu, cv) dlnorm(x, mconv(mu, mu * cv), sdconv(mu, mu * cv))

#'  Condition met?
#' 
#' @usage condmet(vec)
#' @param vec vector of logical values 
#' @export condmet
condmet <- function(vec) TRUE %in% vec

#'  Sample vector
#' 
#' @usage sampy(x)
#' @param x vector of values 
#' @export sampy
sampy <- function(x) sample(x, 1, prob = !is.na(x))

#' Standardize values
#' 
#' Function to standardize to value relative to minimum and maximum values
#' 
#' @usage Range(x, Max, Min)
#' @param x vector of values 
#' @param Max Maximum value
#' @param Min Minimum value 
#' @export Range
#'  
Range <- function(x, Max, Min) (x - Min)/(Max - Min)  


#' @rdname Range
#' @export range01
range01 <- function(x) (x - min(x))/(max(x) - min(x))  


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
        Effs <- mapply(runif, n = nsim, min = EffLower, max = EffUpper)  # sample Effort
        if (nsim > 1) {
            effort <- t(sapply(1:nsim, function(x) approx(x = refYear, 
                y = Effs[x, ], method = "linear", n = nyears)$y))  # linear interpolation
        }
        if (nsim == 1) {
            # Effs <- Effs/max(Effs)
            effort <- matrix(approx(x = refYear, y = Effs, method = "linear", 
                n = nyears)$y, nrow = 1)
        }
        
        effort <- range01(effort)
        effort[effort == 0] <- 0.01
        
        Emu <- -0.5 * Esd^2
        Eerr <- array(exp(rnorm(nyears * nsim, rep(Emu, nyears), rep(Esd, 
            nyears))), c(nsim, nyears))  # calc error
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


#' Creates a time series per simulation that has gradient grad and random normal walk with sigma
#' 
#' @param targ mean
#' @param targsd standard deviation
#' @param targgrad gradient 
#' @param nyears number of historical years
#' @param nsim number of simulations
#' @export gettempvar
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




