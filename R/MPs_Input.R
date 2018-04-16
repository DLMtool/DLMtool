## Size limit MPs ####

#' Size limit management procedures
#' 
#' Examples of the implementation of input controls in the DLM toolkit. See details below.
#' 
#' 
#' @param x A position in a data-limited methods object
#' @param Data A data-limited methods object
#' @param buffer Parameter controlling the fraction of Lopt to set the minimum
#' length of fish caught: minlen=Lopt*(0.7+buffer).
#' @param ... Optional additional arguments that are ignored. Note arguments
#' \code{reps} or \code{...} are required for all input controls
#' @return A Rec object with input control recommendations 
#' @author T. Carruthers & A. Hordyk
#' @references
#' Hordyk, A., Ono, K., Sainsbury, K., Loneragan, N., and J.
#' Prince. 2014. Some explorations of the life history ratios to describe
#' length composition, spawning-per-recruit, and the spawning potential ratio
#' ICES Journal of Marine Science, doi:10.1093/icesjms/fst235.
#' @describeIn matlenlim A data-limited method in which fishing retention-at-length
#' is set equivalent to the maturity curve.
#' @export 
matlenlim <- function(x, Data, ...) {
  # Knife-edge vulnerability at estimated length-at-maturity  
  dependencies = "Data@L50"
  
  rec <- new("Rec") # create recommendation object
  rec@LR5 <- Data@L50[x] * 0.95 # new length at 5% retention  
  rec@LFR <-  Data@L50[x] # new length at full retention   
  
  # other slots aren't specified so remain unchanged
  rec
}
class(matlenlim) <- "MP"


#' @describeIn matlenlim Selectivity-at-length is set slightly higher (110\%) 
#' than the maturity-at-length.
#' @export matlenlim2
matlenlim2 <- function(x, Data, ...) {
  # Knife-edge vulnerability slightly higher than length at maturity
  dependencies = "Data@L50"
  
  rec <- new("Rec") # create recommendation object
  rec@LFR <-  1.1 * Data@L50[x]  # new length at full retention   
  rec@LR5 <-  0.95 * rec@LFR # new length at 5% retention
  # other slots aren't specified so remain unchanged
  return(rec)
}
class(matlenlim2) <- "MP"



#' @describeIn matlenlim This input control sets the minimum length of fish 
#' caught to a fraction of the length that maximises the biomass, Lopt. The aim 
#' of this simple MP is restrict the catch of small fish to rebuild
#' the stock biomass towards the optimal length, Lopt, expressed in terms of
#' the growth parameters Lopt=b/(M/k+b) (Hordyk et al. 2014) (Author: HF Geromont)
#' @export minlenLopt1
minlenLopt1 <- function(x, Data, buffer = 0.1, ...) {
  
  # Minimum length MPs: Fix length-at-full-selectivity to 0.8*Lopt and
  # set length-at-first-capture 10% below LFs
  
  dependencies = "Data@MPeff, Data@vbLinf, Data@wlb, Data@Mort, Data@vbK"
  Lopt <- Data@vbLinf[x] * Data@wlb[x]/((Data@Mort[x]/Data@vbK[x]) +  Data@wlb[x])
  
  rec <- new("Rec") # create recommendation object
  rec@LFR <- Lopt * (0.7 + buffer) # Lopt too precautionary, so set it to % below
  rec@LR5 <- rec@LFR * 0.9
  rec
  
}
class(minlenLopt1) <- "MP"


#' @describeIn matlenlim Selectivity-at-length is set using a slot limit; that is, a minimum and
#' maximum legal length.  The maximum limit is set here, quite arbitrarily, as
#' the 75th percentile between the new minimum legal length and the estimated
#' asymptotic length Linf.
#' @export slotlim
slotlim <- function(x, Data, ...) {
  # Example of slot limit between 0.95 and 1.25 * L50
  dependencies = "Data@L50, Data@vbLinf"
  
  rec <- new("Rec") # create recommendation object
  rec@LFR <- 1.1 * Data@L50[x]
  rec@LR5 <- 0.95 * rec@LFR 
  rec@HS <- as.numeric(quantile(c(rec@LFR , Data@vbLinf[x]), 0.75))
  
  rec
}
class(slotlim) <- "MP"



# --- Spatial Closure MPs ----

#' Spatial closure and allocation management procedures
#' 
#' Management procedures which can close areas to fishing and reallocate 
#' fishing effort spatially.
#' 
#' @describeIn MRreal A spatial control that prevents fishing in area 1 and reallocates this
#' fishing effort to area 2 (or over other areas).
#' 
#' 
#' @param x A position in data / simulation object DLM
#' @param Data A data limited methods data object
#' @param ... Optional additional arguments that are ignored. Note arguments
#' \code{reps} or \code{...} are required for all input controls
#' @author T. Carruthers
#' @export 
MRreal <- function(x, Data, ...) {
  # A Marine reserve in area 1 with spatial reallocation of effort
  
  rec <- new("Rec") # create recommendation object
  rec@Allocate <- 1
  rec@Spatial <- c(0, rep(1, Data@nareas-1))
  
  # other slots aren't specified so remain unchanged
  return(rec)
}
class(MRreal) <- "MP"


#' @describeIn MRreal A spatial control that prevents fishing in area 1 
#' (e.g., An marine reserve) and does not reallocate this fishing effort to area 2.
#' @export MRnoreal
MRnoreal <- function(x, Data, ...) {
  # A Marine reserve in area 1 with no spatial reallocation of effort
  
  rec <- new("Rec") # create recommendation object
  rec@Allocate <- 0
  rec@Spatial <- c(0, rep(1, Data@nareas-1))
  
  # other slots aren't specified so remain unchanged
  return(rec)
}
class(MRnoreal) <- "MP"



# --- Effort Control MPs ----

#' Fishing at current effort levels
#' 
#' Constant fishing effort set at final year of historical simulations subject
#' to changes in catchability determined by OM@@qinc and interannual variability
#' in catchability determined by OM@@qcv. This MP is intended to represent a
#' 'status quo' management approach.
#' 
#' @param x A position in a data-limited methods data object.
#' @param Data A data-limited methods data object.
#' @param ... Optional additional arguments (currently not used).
#' @note Made up for this package.
#' @author T. Carruthers.
#' @describeIn curE Set effort to 100\% of that in final year of historical simulations.
#' @export 
curE <- function(x, Data, ...) {
  # current effort
  rec <- new("Rec") # create recommendation object
  rec@Effort <- 1
  rec
}
class(curE) <- "MP"


#' @describeIn curE Set effort to 75\% of that in final year.
#' @export 
curE75 <- function(x, Data, ...) {
  # 75% current effort
  rec <- new("Rec") # create recommendation object
  rec@Effort <- 0.75
  rec
}
class(curE75) <- "MP"


#' @describeIn DD Effort-control version. The recommended effort is EMSY.
#' @importFrom stats approx 
#' @export 
DDe <- function(x, Data, reps = 100) {
  dependencies = "Data@vbLinf, Data@vbK, Data@vbt0, Data@Mort, Data@wla, Data@wlb, Data@Cat, Data@Ind, Data@L50, Data@MaxAge"
  Winf = Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
  age <- 1:Data@MaxAge
  la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
  wa <- Data@wla[x] * la^Data@wlb[x]
  a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], 
              Data@L50[x])
  a50V <- max(a50V, 1)
  #yind <- (1:length(Data@Cat[x, ]))[!is.na(Data@Cat[x, ] + Data@Ind[x, 
  #                                                                  ])]
  yind <- !is.na(Data@Cat[x, ] + Data@Ind[x, ])[1] # First year with both catch and index
  yind <- yind:length(Data@Cat[x, ])
  Year <- Data@Year[yind]
  C_hist <- Data@Cat[x, yind]
  I_hist <- Data@Ind[x, yind]
  if(any(is.na(C_hist))) { # Linear interpolation of any missing catch
    C.xind <- 1:length(C_hist)
    C_hist <- approx(x = C.xind[!is.na(C_hist)], y = C_hist[!is.na(C_hist)], n = length(C.xind))$y
  }
  if(any(is.na(I_hist))) { # Linear interpolation of any missing index
    I.xind <- 1:length(I_hist)
    I_hist <- approx(x = I.xind[!is.na(I_hist)], y = I_hist[!is.na(I_hist)], n = length(I.xind))$y
  }
  E_hist <- C_hist/I_hist
  E_hist <- E_hist/mean(E_hist)
  ny_DD <- length(C_hist)
  k_DD <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)  --
  k_DD[k_DD > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho_DD <- (wa[k_DD + 2] - Winf)/(wa[k_DD + 1] - Winf)
  Alpha_DD <- Winf * (1 - Rho_DD)
  So_DD <- exp(-Data@Mort[x])  # get So survival rate
  wa_DD <- wa[k_DD]
  UMSYpriorpar <- c(1 - exp(-Data@Mort[x] * 0.5), 0.3) # Prior for UMSY is that corresponding to F = 0.5 M with CV = 0.3
  UMSYprior <- c(alphaconv(UMSYpriorpar[1], prod(UMSYpriorpar)), betaconv(UMSYpriorpar[1], prod(UMSYpriorpar))) # Convert to beta parameters
  params <- log(c(UMSYpriorpar[1]/(1 - UMSYpriorpar[1]), 3*mean(C_hist, na.rm = T), Data@Mort[x]))
  opt <- optim(params, DD_R, opty = 1, So_DD = So_DD, Alpha_DD = Alpha_DD, 
               Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, wa_DD = wa_DD, E_hist = E_hist, 
               C_hist = C_hist, UMSYprior = UMSYprior, method = "BFGS") # hessian = FALSE
  
  #U_hist <- 1 - exp(-exp(opt$par[3]) * E_hist)
  UMSY <- 1/(1 + exp(-opt$par[1]))
  EMSY <- -log(1 - UMSY)/exp(opt$par[3]) # q in denominator
  
  Allocate <- 1
  #eff <- exp(opt$par[1])/U_hist[Data@LHYear]
  eff <- EMSY/E_hist[Data@LHYear - Year[1] + 1] # Effort advice is ratio of EMSY and obs. Eff in LHYear
  eff[!is.finite(eff)] <- 0.01
  eff[eff > 1e+05] <- 0.01
  rec <- new("Rec")
  rec@Effort <- max(0.01, eff)
  rec 
}
class(DDe) <- "MP"

#' @describeIn DD Variant of \code{DDe} that limits the maximum change in effort to 10 percent.
#' @importFrom stats approx
#' @export DDes
DDes <- function(x, Data, reps = 100, LB = 0.9, UB = 1.1) {
  dependencies = "Data@vbLinf, Data@vbK, Data@vbt0, Data@Mort, Data@wla, Data@wlb, Data@Cat, Data@Ind, Data@L50, Data@MaxAge"
  Winf = Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
  age <- 1:Data@MaxAge
  la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
  wa <- Data@wla[x] * la^Data@wlb[x]
  a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], 
              Data@L50[x])
  a50V <- max(a50V, 1)
  #yind <- (1:length(Data@Cat[x, ]))[!is.na(Data@Cat[x, ] + Data@Ind[x, 
  #                                                                  ])]
  yind <- !is.na(Data@Cat[x, ] + Data@Ind[x, ])[1] # First year with both catch and index
  yind <- yind:length(Data@Cat[x, ])
  Year <- Data@Year[yind]
  C_hist <- Data@Cat[x, yind]
  I_hist <- Data@Ind[x, yind]
  if(any(is.na(C_hist))) { # Linear interpolation of any missing catch
    C.xind <- 1:length(C_hist)
    C_hist <- approx(x = C.xind[!is.na(C_hist)], y = C_hist[!is.na(C_hist)], n = length(C.xind))$y
  }
  if(any(is.na(I_hist))) { # Linear interpolation of any missing index
    I.xind <- 1:length(I_hist)
    I_hist <- approx(x = I.xind[!is.na(I_hist)], y = I_hist[!is.na(I_hist)], n = length(I.xind))$y
  }
  E_hist <- C_hist/I_hist
  E_hist <- E_hist/mean(E_hist)
  ny_DD <- length(C_hist)
  k_DD <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)  --
  k_DD[k_DD > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho_DD <- (wa[k_DD + 2] - Winf)/(wa[k_DD + 1] - Winf)
  Alpha_DD <- Winf * (1 - Rho_DD)
  So_DD <- exp(-Data@Mort[x])  # get So survival rate
  wa_DD <- wa[k_DD]
  UMSYpriorpar <- c(1 - exp(-Data@Mort[x] * 0.5), 0.3) # Prior for UMSY is that corresponding to F = 0.5 M with CV = 0.3
  UMSYprior <- c(alphaconv(UMSYpriorpar[1], prod(UMSYpriorpar)), betaconv(UMSYpriorpar[1], prod(UMSYpriorpar))) # Convert to beta parameters
  params <- log(c(UMSYpriorpar[1]/(1 - UMSYpriorpar[1]), 3*mean(C_hist, na.rm = T), Data@Mort[x]))
  opt <- optim(params, DD_R, opty = 1, So_DD = So_DD, Alpha_DD = Alpha_DD, 
               Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, wa_DD = wa_DD, E_hist = E_hist, 
               C_hist = C_hist, UMSYprior = UMSYprior, method = "BFGS") # hessian = FALSE

  #U_hist <- 1 - exp(-exp(opt$par[3]) * E_hist)
  UMSY <- 1/(1 + exp(-opt$par[1]))
  EMSY <- -log(1 - UMSY)/exp(opt$par[3]) # q in denominator
  #fac <- UMSY/U_hist[Data@LHYear]
  fac <- EMSY/E_hist[Data@LHYear - Year[1] + 1] # Effort advice is ratio of EMSY and obs. Eff in LHYear
  
  if (fac < LB) 
    fac <- LB
  if (fac > UB) 
    fac <- UB
  
  rec <- new("Rec")
  rec@Effort <- max(0.01, Data@MPeff[x] * fac)
  rec
  
}
class(DDes) <- "MP"



#' @describeIn DD Variant of \code{DDe} where the recommended effort is 75\% EMSY.
#' @importFrom stats approx
#' @export 
DDe75 <- function(x, Data, reps = 100) {
  dependencies = "Data@vbLinf, Data@vbK, Data@vbt0, Data@Mort, Data@wla, Data@wlb, Data@Cat, Data@Ind, Data@L50, Data@MaxAge"
  Winf = Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
  age <- 1:Data@MaxAge
  la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
  wa <- Data@wla[x] * la^Data@wlb[x]
  a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], 
              Data@L50[x])
  a50V <- max(a50V, 1)
  #yind <- (1:length(Data@Cat[x, ]))[!is.na(Data@Cat[x, ] + Data@Ind[x, 
  #                                                                  ])]
  yind <- !is.na(Data@Cat[x, ] + Data@Ind[x, ])[1] # First year with both catch and index
  yind <- yind:length(Data@Cat[x, ])
  Year <- Data@Year[yind]
  C_hist <- Data@Cat[x, yind]
  I_hist <- Data@Ind[x, yind]
  if(any(is.na(C_hist))) { # Linear interpolation of any missing catch
    C.xind <- 1:length(C_hist)
    C_hist <- approx(x = C.xind[!is.na(C_hist)], y = C_hist[!is.na(C_hist)], n = length(C.xind))$y
  }
  if(any(is.na(I_hist))) { # Linear interpolation of any missing index
    I.xind <- 1:length(I_hist)
    I_hist <- approx(x = I.xind[!is.na(I_hist)], y = I_hist[!is.na(I_hist)], n = length(I.xind))$y
  }
  E_hist <- C_hist/I_hist
  E_hist <- E_hist/mean(E_hist)
  ny_DD <- length(C_hist)
  params <- log(c(Data@Mort[x], mean(C_hist, na.rm = T), Data@Mort[x]))
  k_DD <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)  --
  k_DD[k_DD > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho_DD <- (wa[k_DD + 2] - Winf)/(wa[k_DD + 1] - Winf)
  Alpha_DD <- Winf * (1 - Rho_DD)
  So_DD <- exp(-Data@Mort[x])  # get So survival rate
  wa_DD <- wa[k_DD]
  UMSYpriorpar <- c(1 - exp(-Data@Mort[x] * 0.5), 0.3) # Prior for UMSY is that corresponding to F = 0.5 M with CV = 0.3
  UMSYprior <- c(alphaconv(UMSYpriorpar[1], prod(UMSYpriorpar)), betaconv(UMSYpriorpar[1], prod(UMSYpriorpar))) # Convert to beta parameters
  params <- log(c(UMSYpriorpar[1]/(1 - UMSYpriorpar[1]), 3*mean(C_hist, na.rm = T), Data@Mort[x]))
  opt <- optim(params, DD_R, opty = 1, So_DD = So_DD, Alpha_DD = Alpha_DD, 
               Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, wa_DD = wa_DD, E_hist = E_hist, 
               C_hist = C_hist, UMSYprior = UMSYprior, method = "BFGS") #hessian = FALSE
  
  #U_hist <- 1 - exp(-exp(opt$par[3]) * E_hist)
  UMSY <- 1/(1 + exp(-opt$par[1]))
  EMSY <- -log(1 - UMSY)/exp(opt$par[3]) # q in denominator
  
  Allocate <- 1
  #eff <- 0.75 *exp( opt$par[1])/U_hist[Data@LHYear]
  eff <- 0.75 *EMSY/E_hist[Data@LHYear - Year[1] + 1] # Effort advice is ratio of EMSY and obs. Eff in LHYear
  eff[!is.finite(eff)] <- 0.01
  eff[eff > 1e+05] <- 0.01
  rec <- new("Rec")
  rec@Effort <- max(0.01, eff)
  rec 
}
class(DDe75) <- "MP"




#' Effort searching MP aiming for 40 per cent stock depletion
#' 
#' A very simple MP that modifies effort to reach 40 percent stock depletion
#' 
#' 
#' @usage DTe40(x, Data, reps = 100, alpha=0.4, LB=0.9, UB=1.1)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @param alpha The target level of depletion
#' @param LB The lowest permitted factor of previous fishing effort
#' @param UB The highest permitted factor of previous fishing effort
#' @author T. Carruthers
#' @export DTe40
DTe40 <- function(x, Data, reps = 100, alpha = 0.4, LB = 0.9, UB = 1.1) {
  
  dependencies = "Data@Dep"
  
  fac <- Data@Dep[x]/alpha
  
  if (fac < LB) 
    fac <- LB
  if (fac > UB) 
    fac <- UB
  
  rec <- new("Rec")
  rec@Effort <- max(0.01, Data@MPeff[x] * fac)
  rec
  
}
class(DTe40) <- "MP"



#' Effort searching MP aiming for 50 per cent stock depletion
#' 
#' A very simple MP that modifies effort to reach 50 percent stock depletion
#' 
#' 
#' @usage DTe50(x, Data, reps = 100, alpha=0.5, LB=0.9, UB=1.1)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @param alpha The target level of depletion
#' @param LB The lowest permitted factor of previous fishing effort
#' @param UB The highest permitted factor of previous fishing effort
#' @author T. Carruthers
#' @export DTe50
DTe50 <- function(x, Data, reps = 100, alpha = 0.5, LB = 0.9, UB = 1.1) {
  
  dependencies = "Data@Dep"
  
  fac <- Data@Dep[x]/alpha
  if (fac < LB) 
    fac <- LB
  if (fac > UB) 
    fac <- UB
  
  rec <- new("Rec")
  rec@Effort <- max(0.01, Data@MPeff[x] * fac)
  rec
  
}
class(DTe50) <- "MP"


#' Effort MP: adjust effort up/down if mean length above/below Ltarget
#' 
#' 
#' @usage EtargetLopt(x, Data, reps = 100, yrsmth=3, buffer=0.1)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrsmth Number of years to calculate average length
#' @param buffer Parameter controlling the fraction of mean catch to set the
#' reference (or target) TAC level - acts as a precautionary buffer
#' @return An adjustment for fishing effort
#' @author HF Geromont
#' @export EtargetLopt
EtargetLopt <- function(x, Data, reps = 100, yrsmth = 3, buffer = 0.1) {
  
  # Effort MP: adjust effort up/down if mean length above/below Ltarget
  
  dependencies = "Data@Year, Data@ML, Data@L50, Data@MPeff, Data@vbLinf, Data@wlb, Data@Mort, Data@vbK"
  
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 3 years
  Lrecent <- mean(Data@ML[ind])
  Lopt <- Data@vbLinf[x] * Data@wlb[x]/((Data@Mort[x]/Data@vbK[x]) + 
                                          Data@wlb[x])
  ratio <- Lrecent/Lopt
  
  rec <- new("Rec")
  w <- 0.5
  rec@Effort <- (1 - buffer) * (w + (1 - w) * ratio)
  rec 
}
class(EtargetLopt) <- "MP"


#' Index Target Effort-Based 5
#' 
#' An index target MP where the Effort is modified according to current index
#' levels (mean index over last 5 years) relative to a target level. Maximum
#' annual changes are 5 per cent.
#' 
#' 
#' @usage ITe5(x, Data, reps = 100, yrsmth = 5, mc = 0.05)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the quota recommendation
#' @param yrsmth The number of historical years over which to average the index
#' @param mc The maximum fractional change in the effort among years.
#' @return A numeric vector of input controls
#' @author T. Carruthers
#' @export ITe5
ITe5 <- function(x, Data, reps = 100, yrsmth = 5, mc = 0.05) {
  
  dependencies = "Data@Ind, Data@MPeff, Data@CV_Ind, Data@Iref"
  ind <- max(1, (length(Data@Year) - yrsmth + 1)):length(Data@Year)
  deltaI <- mean(Data@Ind[x, ind])/Data@Iref[x]
  if (deltaI < (1 - mc)) deltaI <- 1 - mc
  if (deltaI > (1 + mc)) deltaI <- 1 + mc
  
  Effort <- Data@MPeff[x] * deltaI * trlnorm(reps, 1, Data@CV_Ind[x])
  if (reps == 1)  Effort <- Data@MPeff[x] * deltaI
  Effort[Effort < 0.01] <- 0.01  # for simulations in case Effort goes negative
  rec <- new("Rec")
  rec@Effort <- mean(Effort) 
  rec 
}
class(ITe5) <- "MP"



#' Index Target Effort-Based 10
#' 
#' An index target MP where the Effort is modified according to current index
#' levels (mean index over last 5 years) relative to a target level. Maximum
#' annual changes are 10 per cent.
#' 
#' 
#' @usage ITe10(x, Data, reps = 100, yrsmth = 5, mc = 0.1)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the quota recommendation
#' @param yrsmth The number of historical years over which to average the index
#' @param mc The maximum fractional change in the Effort among years.
#' @return A numeric vector of input controls
#' @author T. Carruthers
#' @export ITe10
ITe10 <- function(x, Data, reps = 100, yrsmth = 5, mc = 0.1) {
  
  dependencies = "Data@Ind, Data@MPeff, Data@CV_Ind, Data@Iref"
  ind <- max(1, (length(Data@Year) - yrsmth + 1)):length(Data@Year)
  
  deltaI <- mean(Data@Ind[x, ind])/Data@Iref[x]
  if (deltaI < (1 - mc)) deltaI <- 1 - mc
  if (deltaI > (1 + mc)) deltaI <- 1 + mc
  
  Effort <- Data@MPeff[x] * deltaI * trlnorm(reps, 1, Data@CV_Ind[x])
  if (reps == 1) 
    Effort <- Data@MPeff[x] * deltaI
  Allocate <- 1
  Effort[Effort < 0.01] <- 0.01  # for simulations in case Effort goes negative
  rec <- new("Rec")
  rec@Effort <- mean(Effort) 
  rec 
}
class(ITe10) <- "MP"




#' A management procedure that incrementally adjusts the effort to reach a
#' target CPUE / relative abundance index
#' 
#' An effort-based version of the least biologically precautionary of two
#' index/CPUE target MPs proposed by Geromont and Butterworth 2014. Tested by
#' Carruthers et al. 2015
#' 
#' Tested by Carruthers et al. 2015.
#' 
#' @usage ItargetE1(x, Data, reps = 100, yrsmth = 5, xx = 0, Imulti = 1.5)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param Imulti Parameter controlling how much larger target CPUE / index is
#' compared with recent levels.
#' @return A numeric vector of input controls
#' @author T. Carruthers
#' @references Carruthers et al. 2015. Performance evaluation of simple
#' management procedures. Fish and Fisheries. In press.
#' 
#' Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for
#' data-poor fisheries; forecasting with few data. ICES J. Mar. Sci.
#' doi:10.1093/icesjms/fst232
#' @export ItargetE1
ItargetE1 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0, Imulti = 1.5) {
  
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 5 years
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  ind2 <- ((ylast - (yrsmth - 1)):ylast)  # historical 5 pre-projection years
  ind3 <- ((ylast - (yrsmth * 2 - 1)):ylast)  # historical 10 pre-projection years
  
  Irecent <- mean(Data@Ind[x, ind])
  Iave <- mean(Data@Ind[x, ind3])
  Itarget <- Iave * Imulti
  I0 <- 0.8 * Iave
  if (Irecent > I0) {
    Effort <- 0.5 * Data@MPeff[x] * (1 + ((Irecent - I0)/(Itarget - 
                                                            I0)))
  } else {
    Effort <- 0.5 * Data@MPeff[x] * (Irecent/I0)^2
  }
  
  Step <- (Effort/Data@MPeff[x])  # step change in effort 
  Step[Step < 0.85] <- 0.85
  Step[Step > 1.15] <- 1.15
  Allocate <- 1
  Effort <- Step * Data@MPeff[x]
  Effort[Effort < 0.01] <- 0.01  # for simulations in case Effort goes negative
  rec <- new("Rec")
  rec@Effort <- mean(Effort)
  rec
}
class(ItargetE1) <- "MP"



#' A management procedure that incrementally adjusts the Effort to reach a
#' target CPUE / relative abundance index
#' 
#' An effort-based version of the most biologically precautionary of two
#' index/CPUE target MPs proposed by Geromont and Butterworth 2014.
#' 
#' Tested by Carruthers et al. 2015.
#' 
#' @usage ItargetE4(x, Data, reps = 100, yrsmth = 5, xx = 0, Imulti = 2.5)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param Imulti Parameter controlling how much larger target CPUE / index is
#' compared with recent levels.
#' @return A numeric vector of input controls
#' @author T. Carruthers
#' @references Carruthers et al. 2015. Performance evaluation of simple
#' management procedures. Fish and Fisheries. In press.
#' 
#' Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for
#' data-poor fisheries; forecasting with few data. ICES J. Mar. Sci.
#' doi:10.1093/icesjms/fst232
#' @export ItargetE4
ItargetE4 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0, Imulti = 2.5) {
  
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 5 years
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  ind2 <- ((ylast - (yrsmth - 1)):ylast)  # historical 5 pre-projection years
  ind3 <- ((ylast - (yrsmth * 2 - 1)):ylast)  # historical 10 pre-projection years
  
  Irecent <- mean(Data@Ind[x, ind])
  Iave <- mean(Data@Ind[x, ind3])
  Itarget <- Iave * Imulti
  I0 <- 0.8 * Iave
  if (Irecent > I0) {
    Effort <- 0.5 * Data@MPeff[x] * (1 + ((Irecent - I0)/(Itarget - 
                                                            I0)))
  } else {
    Effort <- 0.5 * Data@MPeff[x] * (Irecent/I0)^2
  }
  Step <- (Effort/Data@MPeff[x])  # step change in effort 
  Step[Step < 0.8] <- 0.8
  Step[Step > 1.2] <- 1.2
  
  Allocate <- 1
  Effort <- Step * Data@MPeff[x]
  Effort[Effort < 0.01] <- 0.01  # for simulations in case Effort goes negative
  rec <- new("Rec")
  rec@Effort <- mean(Effort)
  rec
  
}
class(ItargetE4) <- "MP"




#' A management procedure that incrementally adjusts the TAC according to the
#' mean length of recent catches.
#' 
#' A effort-based version of least biologically precautionary of four adaptive
#' length-based MPs proposed by Geromont and Butterworth 2014. Tested by
#' Carruthers et al. 2015
#' 
#' 
#' @usage LstepCE1(x, Data, reps = 100, yrsmth = 5, xx = 0, stepsz = 0.05,
#' llim = c(0.96, 0.98, 1.05))
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of effort samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param stepsz Parameter controlling the size of the effort update increment.
#' @param llim A vector of length reference points that determine the
#' conditions for increasing, maintaining or reducing the effort.
#' @return A numeric vector of input controls
#' @author T. Carruthers
#' @export LstepCE1
LstepCE1 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0, stepsz = 0.05, 
                     llim = c(0.96, 0.98, 1.05)) {
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 5 years
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  # ind2<-((ylast-(yrsmth-1)):ylast) # historical 5 pre-projection years
  ind3 <- ((ylast - (yrsmth * 2 - 1)):ylast)  # historical 10 pre-projection years
  
  Lrecent <- mean(Data@ML[ind])
  Lave <- mean(Data@ML[ind3])
  rat <- Lrecent/Lave
  
  step <- stepsz
  
  if (rat < llim[1]) {
    Effort <- Data@MPeff[x] - 2 * (step * Data@MPeff[x])
  } else if (rat < llim[2]) {
    Effort <- Data@MPeff[x] - (step * Data@MPeff[x])
  } else if (rat > llim[3]) {
    Effort <- Data@MPeff[x] + (step * Data@MPeff[x])
  } else {
    Effort <- Data@MPeff[x]
  }
  
  Allocate <- 1
  Effort[Effort < 0.01] <- 0.01  # for simulations in case Effort goes negative
  rec <- new("Rec")
  rec@Effort <- mean(Effort)
  rec
  
}
class(LstepCE1) <- "MP"



#' A management procedure that incrementally adjusts the Effort according to
#' the mean length of recent catches.
#' 
#' A effort-based version of one of the four adaptive length-based MPs proposed
#' by Geromont and Butterworth 2014.
#' 
#' 
#' @usage LstepCE2(x, Data, reps = 100, yrsmth = 5, xx = 0, stepsz = 0.1,
#' llim = c(0.96, 0.98, 1.05))
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param stepsz Parameter controlling the size of the effort update increment.
#' @param llim A vector of length reference points that determine the
#' conditions for increasing, maintaining or reducing the effort.
#' @return A numeric vector of input controls
#' @author T. Carruthers
#' @export LstepCE2
LstepCE2 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0, stepsz = 0.1, 
                     llim = c(0.96, 0.98, 1.05)) {
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 5 years
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  # ind2<-((ylast-(yrsmth-1)):ylast) # historical 5 pre-projection years
  ind3 <- ((ylast - (yrsmth * 2 - 1)):ylast)  # historical 10 pre-projection years
  
  Lrecent <- mean(Data@ML[ind])
  Lave <- mean(Data@ML[ind3])
  rat <- Lrecent/Lave
  step <- stepsz
  
  if (rat < llim[1]) {
    Effort <- Data@MPeff[x] - 2 * (step * Data@MPeff[x])
  } else if (rat < llim[2]) {
    Effort <- Data@MPeff[x] - (step * Data@MPeff[x])
  } else if (rat > llim[3]) {
    Effort <- Data@MPeff[x] + (step * Data@MPeff[x])
  } else {
    Effort <- Data@MPeff[x]
  }
  Effort[Effort < 0.01] <- 0.01  # for simulations in case Effort goes negative
  rec <- new("Rec")
  rec@Effort <- mean(Effort)
  rec
}
class(LstepCE2) <- "MP"



#' A management procedure that incrementally adjusts the Effort to reach a
#' target mean length in catches.
#' 
#' A effort based version of the least biologically precautionary of four
#' target length MPs proposed by Geromont and Butterworth 2014.
#' 
#' 
#' @usage LtargetE1(x, Data, reps = 100, yrsmth = 5, xx = 0, xL = 1.05)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param xL Parameter controlling the magnitude of the target mean length of
#' catches relative to average length in catches.
#' @return A numeric vector of input controls
#' @author T. Carruthers
#' @export LtargetE1
LtargetE1 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0, xL = 1.05) {
  
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 5 years
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  ind2 <- ((ylast - (yrsmth - 1)):ylast)  # historical 5 pre-projection years
  ind3 <- ((ylast - (yrsmth * 2 - 1)):ylast)  # historical 10 pre-projection years
  
  Lrecent <- mean(Data@ML[ind])
  Lave <- mean(Data@ML[ind3])
  L0 <- 0.9 * Lave
  Ltarget <- xL * Lave
  if (Lrecent > L0) {
    Effort <- 0.5 * Data@MPeff[x] * (1 + ((Lrecent - L0)/(Ltarget - 
                                                            L0)))
  } else {
    Effort <- 0.5 * Data@MPeff[x] * (Lrecent/L0)^2
  }
  Step <- (Effort/Data@MPeff[x])  # step change in effort 
  Step[Step < 0.85] <- 0.85
  Step[Step > 1.15] <- 1.15
  
  Allocate <- 1
  Effort <- Step * Data@MPeff[x]
  Effort[Effort < 0.01] <- 0.01  # for simulations in case Effort goes negative
  rec <- new("Rec")
  rec@Effort <- mean(Effort)
  rec
}
class(LtargetE1) <- "MP"



#' A management procedure that incrementally adjusts the Effort to reach a
#' target mean length in catches.
#' 
#' A effort based version of the most biologically precautionary of four target
#' length MPs proposed by Geromont and Butterworth 2014.
#' 
#' 
#' @usage LtargetE4(x, Data, reps = 100, yrsmth = 5, xx = 0, xL = 1.15)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param xL Parameter controlling the magnitude of the target mean length of
#' catches relative to average length in catches.
#' @return A numeric vector of input controls
#' @author T. Carruthers
#' @export LtargetE4
LtargetE4 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0, xL = 1.15) {
  
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 5 years
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  ind2 <- ((ylast - (yrsmth - 1)):ylast)  # historical 5 pre-projection years
  ind3 <- ((ylast - (yrsmth * 2 - 1)):ylast)  # historical 10 pre-projection years
  
  Lrecent <- mean(Data@ML[ind])
  Lave <- mean(Data@ML[ind3])
  L0 <- 0.9 * Lave
  Ltarget <- xL * Lave
  if (Lrecent > L0) {
    Effort <- 0.5 * Data@MPeff[x] * (1 + ((Lrecent - L0)/(Ltarget - 
                                                            L0)))
  } else {
    Effort <- 0.5 * Data@MPeff[x] * (Lrecent/L0)^2
  }
  
  Step <- (Effort/Data@MPeff[x])  # step change in effort 
  Step[Step < 0.8] <- 0.8
  Step[Step > 1.2] <- 1.2
  Allocate <- 1
  Effort <- Step * Data@MPeff[x]
  Effort[Effort < 0.01] <- 0.01  # for simulations in case Effort goes negative
  rec <- new("Rec")
  rec@Effort <- mean(Effort)
  rec
}
class(LtargetE4) <- "MP"