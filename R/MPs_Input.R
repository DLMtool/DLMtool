## Size limit MPs ####

#' Size limit management procedures
#' 
#' Examples of the implementation of input controls in the DLM toolkit. See details below.
#' 
#' @templateVar mp matlenlim
#' @template MPtemplate
#' @template MPuses
#' @author T. Carruthers & A. Hordyk
#' @references
#' Hordyk, A., Ono, K., Sainsbury, K., Loneragan, N., and J.
#' Prince. 2015. Some explorations of the life history ratios to describe
#' length composition, spawning-per-recruit, and the spawning potential ratio
#' ICES Journal of Marine Science, doi:10.1093/icesjms/fst235.
#' @describeIn matlenlim A data-limited method in which fishing retention-at-length
#' is set equivalent to the maturity curve.
#' @examples 
#' matlenlim(1, DLMtool::Atlantic_mackerel, plot=TRUE)
#' @export 
matlenlim <- function(x, Data, reps, plot=FALSE) {
  # Knife-edge vulnerability at estimated length-at-maturity  
  rec <- new("Rec") # create recommendation object
  rec@LFR <- Data@L50[x] # new length at full retention   
  rec@LR5  <- rec@LFR * 0.95 # new length at 5% retention   
  
  # other slots aren't specified so remain unchanged
  rec
}
class(matlenlim) <- "MP"

#' @describeIn matlenlim Selectivity-at-length is set slightly higher (110\%) 
#' than the maturity-at-length.
#' 
#' @templateVar mp matlenlim2
#' @template MPuses
#' @export 
#' @examples 
#' matlenlim2(1, DLMtool::Atlantic_mackerel, plot=TRUE)
matlenlim2 <- function(x, Data, reps, plot=FALSE) {
  # Knife-edge vulnerability slightly higher than length at maturity
  dependencies = "Data@L50"
  
  rec <- new("Rec") # create recommendation object
  rec@LFR <-  1.1 * Data@L50[x]  # new length at full retention   
  rec@LR5 <-  0.95 * rec@LFR # new length at 5% retention
  # other slots aren't specified so remain unchanged
  return(rec)
}
class(matlenlim2) <- "MP"


#' @templateVar mp minlenLopt1
#' @template MPuses
#' @param buffer Parameter controlling the fraction of Lopt to set the minimum
#' length of fish caught: minlen=Lopt*(0.7+buffer).
#' 
#' @describeIn matlenlim This input control sets the minimum length of fish 
#' caught to a fraction of the length that maximises the biomass, Lopt. The aim 
#' of this simple MP is restrict the catch of small fish to rebuild
#' the stock biomass towards the optimal length, Lopt, expressed in terms of
#' the growth parameters Lopt=b/(M/k+b) (Hordyk et al. 2015) (Author: HF Geromont)
#' @export 
#' @examples 
#' minlenLopt1(1, DLMtool::Atlantic_mackerel, plot=TRUE)
minlenLopt1 <- function(x, Data, reps, plot=FALSE, buffer = 0.1) {
  
  # Minimum length MPs: Fix length-at-full-selectivity to 0.8*Lopt and
  # set length-at-first-capture 10% below LFs
  
  dependencies = "Data@vbLinf, Data@wlb, Data@Mort, Data@vbK"
  Lopt <- Data@vbLinf[x] * Data@wlb[x]/((Data@Mort[x]/Data@vbK[x]) +  Data@wlb[x])
  
  rec <- new("Rec") # create recommendation object
  rec@LFR <- Lopt * (0.7 + buffer) # Lopt too precautionary, so set it to % below
  rec@LR5 <- rec@LFR * 0.9
  rec
  
}
class(minlenLopt1) <- "MP"


#' @templateVar mp slotlim
#' @template MPuses
#' @describeIn matlenlim Selectivity-at-length is set using a slot limit; that is, a minimum and
#' maximum legal length.  The maximum limit is set here, quite arbitrarily, as
#' the 75th percentile between the new minimum legal length and the estimated
#' asymptotic length Linf.
#' @export 
#' @examples 
#' slotlim(1, DLMtool::Atlantic_mackerel, plot=TRUE)
slotlim <- function(x, Data, reps, plot=FALSE) {
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
#' Management procedures which close areas to fishing and reallocate 
#' fishing effort spatially.
#' 
#' @describeIn MRreal A spatial control that prevents fishing in area 1 and reallocates this
#' fishing effort to area 2 (or over other areas).
#'
#' @templateVar mp MRreal 
#' @template MPtemplate
#' @template MPuses
#' 
#' @author T. Carruthers
#' @export 
#' @examples 
#' MRreal(1, DLMtool::Atlantic_mackerel, plot=TRUE)
MRreal <- function(x, Data, reps, plot=FALSE) {
  # A Marine reserve in area 1 with spatial reallocation of effort
  
  rec <- new("Rec") # create recommendation object
  rec@Allocate <- 1
  rec@Spatial <- c(0, rep(1, Data@nareas-1))
  
  # other slots aren't specified so remain unchanged
  return(rec)
}
class(MRreal) <- "MP"

#' @templateVar mp MRnoreal
#' @template MPuses
#' @describeIn MRreal A spatial control that prevents fishing in area 1 
#' (e.g., A marine reserve) and does not reallocate this fishing effort to area 2.
#' @export 
#' @examples 
#' MRnoreal(1, DLMtool::Atlantic_mackerel, plot=TRUE)
MRnoreal <- function(x, Data, reps, plot=FALSE) {
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
#' @templateVar mp curE 
#' @template MPtemplate
#' @template MPuses
#' 
#' @author T. Carruthers.
#' @describeIn curE Set effort to 100\% of that in final year of historical simulations.
#' @export 
#' @examples 
#' curE(1, DLMtool::Atlantic_mackerel, plot=TRUE)
curE <- function(x, Data, reps, plot=FALSE) {
  # current effort
  rec <- new("Rec") # create recommendation object
  rec@Effort <- 1 * Data@MPeff[x] 
  rec
}
class(curE) <- "MP"

#' @templateVar mp curE75 
#' @template MPuses
#' @describeIn curE Set effort to 75\% of that in final year.
#' @export 
#' @examples 
#' curE75(1, DLMtool::Atlantic_mackerel, plot=TRUE)
curE75 <- function(x, Data, reps, plot=FALSE) {
  # 75% current effort
  rec <- new("Rec") # create recommendation object
  rec@Effort <- 0.75 * Data@MPeff[x]
  rec
}
class(curE75) <- "MP"



#### Delay-Difference Effort MPs ####

#' Effort-based Delay - Difference Stock Assessment with UMSY and MSY as leading parameters
#' 
#' A simple delay-difference assessment that estimates \eqn{E_{\textrm{MSY}}} using a
#' time-series of catches and a relative abundance index. 
#' 
#' This DD model is observation error only and has does not estimate
#' process error (recruitment deviations). Assumption is that knife-edge 
#' selectivity occurs at the age of 50% maturity. Similar to many other assessment
#' models it depends on a whole host of dubious assumptions such as temporally
#' stationary productivity and proportionality between the abundance index and
#' real abundance. Unsurprisingly the extent to which these assumptions are
#' violated tends to be the biggest driver of performance for this method.
#' 
#' The method is conditioned on effort and estimates catch. The effort is calculated
#' as the ratio of catch and index. Thus, to get a complete effort time series, a full
#' time series of catch and index is also needed. Missing values are linearly interpolated.
#' 
#' A detailed description of the delay-difference model can be found in Chapter 9 of Hilborn 
#' and Walters (1992).
#'
#' @templateVar mp DDe 
#' @template MPtemplate
#' @template MPuses
#'
#' @describeIn DDe Effort-control version. The recommended effort is EMSY.
#' @family Delay-Difference MPs
#' @examples 
#' DDe(1, Data=DLMtool::Atlantic_mackerel, plot=TRUE)
#' @export 
DDe <- function(x, Data, reps = 100, plot=FALSE) {
  
  runDD <- DD_(x, Data, reps)
  Eff <- runDD$eff 
  rec <- new("Rec")
  rec@Effort <- max(0.01, Data@MPeff[x] * Eff)
  rec 
}
class(DDe) <- "MP"



#' @param LB The lowest permitted factor of previous fishing effort
#' @param UB The highest permitted factor of previous fishing effort
#' @describeIn DDe Variant of `DDe` that limits the maximum change in effort to 10 percent.
#' @export 
#' @examples
#' DDes(1, Data=DLMtool::Atlantic_mackerel, plot=TRUE)
DDes <- function(x, Data, reps = 100, plot=FALSE, LB = 0.9, UB = 1.1) {
  runDD <- DD_(x, Data, reps)
  Eff <- runDD$eff 
  if (Eff < LB) Eff <- LB
  if (Eff > UB) Eff <- UB
  
  rec <- new("Rec")
  rec@Effort <- max(0.01, Data@MPeff[x] * Eff)
  rec
  
}
class(DDes) <- "MP"



#' @describeIn DDe Variant of `DDe` where the recommended effort is 75\% EMSY.
#' @export 
#' @examples 
#' DDe75(1, Data=DLMtool::Atlantic_mackerel, plot=TRUE)
DDe75 <- function(x, Data, reps = 100) {
  runDD <- DD_(x, Data, reps)
  Eff <- runDD$eff * 0.75
  rec <- new("Rec")
  rec@Effort <- max(0.01, Data@MPeff[x] * Eff)
  rec 
}
class(DDe75) <- "MP"




#' Effort searching MP aiming for a fixed stock depletion
#' 
#' Effort is adjusted using a simple rule that aims for a specified level of depletion.
#' 
#' The TAE is calculated as:
#' \deqn{\textrm{TAE}_y = \frac{D}{\alpha} \textrm{TAE}_{y-1}}
#' where \eqn{D} is estimated current level of depletion and \eqn{\alpha} is argument
#' `alpha` specifying the target level of depletion.
#' 
#' The maximum fractional change in TAE is specified with arguments `LB` and `UB`
#' 
#' @templateVar mp DTe40
#' @template MPtemplate
#' @template MPuses
#' 
#' @param alpha The target level of depletion
#' @param LB The lowest permitted factor of previous fishing effort
#' @param UB The highest permitted factor of previous fishing effort
#' @author T. Carruthers
#' @export 
#' @examples 
#' DTe40(1, DLMtool::Atlantic_mackerel, plot=TRUE)
#' @describeIn DTe40 Effort is adjusted to reach 40 percent stock depletion 
DTe40 <- function(x, Data, reps = 100, plot=FALSE, alpha = 0.4, LB = 0.9, UB = 1.1) {
  
  fac <- Data@Dep[x]/alpha
  
  if (fac < LB) fac <- LB
  if (fac > UB) fac <- UB
  
  rec <- new("Rec")
  rec@Effort <- max(0.01, Data@MPeff[x] * fac)
  rec
  
}
class(DTe40) <- "MP"


#' @describeIn DTe40 Effort is adjusted to reach 50 percent stock depletion
#' @export 
DTe50 <- function(x, Data, reps = 100, plot=FALSE, alpha = 0.5, LB = 0.9, UB = 1.1) {
  
  fac <- Data@Dep[x]/alpha
  if (fac < LB) fac <- LB
  if (fac > UB) fac <- UB
  
  rec <- new("Rec")
  rec@Effort <- max(0.01, Data@MPeff[x] * fac)
  rec
  
}
class(DTe50) <- "MP"




#' Effort Target Optimum Length
#' 
#' This MP adjusts effort limit based on the ratio of recent mean length (over 
#' last `yrsmth` years) and a target length defined as \eqn{L_{\textrm{opt}}}.
#' Effort MP: adjust effort up/down if mean length above/below Ltarget
#' 
#' The TAE is calculated as:
#' \deqn{\textrm{TAE}_y = \textrm{TAE}_{y-1} \left((1-\textrm{buffer}) (w + (1-w)r) \right)}
#' where \eqn{\textrm{buffer}} is specified in argument `buffer`, \eqn{w} is fixed at 0.5, and:
#' \deqn{r = \frac{L_{\textrm{recent}}}{L_{\textrm{opt}}}}
#' where \eqn{L_{\textrm{recent}}}is mean
#' length over last `yrmsth` years, and: 
#' \deqn{L_{\textrm{opt}} = \frac{L_\infty W_b}{\frac{M}{K} + W_b }}
#' where \eqn{L_\infty} is von Bertalanffy asymptotic length, \eqn{W_b} is 
#' exponent of the length-weight relationship, \eqn{M} is natural mortality, and 
#' \eqn{K} is von Bertalanffy growth coefficient.#' 
#' 
#' @templateVar mp EtargetLopt 
#' @template MPtemplate
#' @template MPuses
#' 
#' @param yrsmth Number of years to calculate average length
#' @param buffer Parameter controlling the fraction of calculated effor - acts as a precautionary buffer
#' 
#' @author HF Geromont
#' @export 
#' @examples 
#' EtargetLopt(1, DLMtool::SimulatedData, plot=TRUE)
EtargetLopt <- function(x, Data, reps = 100, plot=FALSE, yrsmth = 3, buffer = 0.1) {
  
  # Effort MP: adjust effort up/down if mean length above/below Ltarget
 
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 3 years
  Lrecent <- mean(Data@ML[ind])
  Lopt <- Data@vbLinf[x] * Data@wlb[x]/((Data@Mort[x]/Data@vbK[x]) + Data@wlb[x])
  ratio <- Lrecent/Lopt
  
  rec <- new("Rec")
  w <- 0.5
  eff <- (1 - buffer) * (w + (1 - w) * ratio)
  rec@Effort <-  max(0.01, Data@MPeff[x] * eff)
  rec 
}
class(EtargetLopt) <- "MP"





#' Index Target Effort-Based 
#' 
#' An index target MP where the Effort is modified according to current index
#' levels (mean index over last 5 years) relative to a target level.
#' 
#' The TAE is calculated as:
#' \deqn{\textrm{TAE}_y = \textrm{TAE}_{y-1} \delta}
#' where \eqn{\delta} is \eqn{\frac{I} {I_{\textrm{ref}}}} averaged over last
#' `yrsmth` years. \eqn{I_{\textrm{ref}}} is the index target (`Data@Iref`).
#' 
#' The maximum fractional change in TAE is specified in `mc`.
#' 
#' @templateVar mp ITe5 
#' @template MPtemplate
#' @template MPuses
#' 
#' @param yrsmth The number of historical years over which to average the index
#' @param mc The maximum fractional change in the effort among years.
#' 
#' @author T. Carruthers
#' @describeIn ITe5  Maximum annual changes are 5 per cent.
#' @examples
#' ITe5(1, DLMtool::SimulatedData, plot=TRUE)
#' @export ITe5
ITe5 <- function(x, Data, reps = 100, plot=FALSE, yrsmth = 5, mc = 0.05) {
  
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



#' Index Target Effort-Based 
#' 
#' 
#' @templateVar mp ITe10 
#' @template MPuses
#' 
#' @describeIn ITe5  Maximum annual changes are 10 per cent.
#' @examples
#' ITe10(1, DLMtool::SimulatedData, plot=TRUE)
#' @export 
ITe10 <- function(x, Data, reps = 100, plot=FALSE, yrsmth = 5, mc = 0.1) {
  
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



#' Incremental Index Target MP Internal Function - Effort MP
#'
#' @param x Iteration number
#' @param Data An object of class `Data`
#' @param reps Number of replicates
#' @param plot Logical. Show the plot?
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param Imulti Parameter controlling how much larger target CPUE / index is
#' compared with recent levels.
#'
#' @return A named list with Effort recommendations 
#' @export
#'
#' @keywords internal
Itargeteff_ <- function(x, Data, reps = 100, plot=FALSE, yrsmth = 5, xx = 0, Imulti = 1.5) {
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 5 years
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  ind2 <- ((ylast - (yrsmth - 1)):ylast)  # historical 5 pre-projection years
  ind3 <- ((ylast - (yrsmth * 2 - 1)):ylast)  # historical 10 pre-projection years
  
  Irecent <- mean(Data@Ind[x, ind])
  Iave <- mean(Data@Ind[x, ind3])
  Itarget <- Iave * Imulti
  I0 <- 0.8 * Iave
  if (Irecent > I0) {
    Effort <- 0.5 * Data@MPeff[x] * (1 + ((Irecent - I0)/(Itarget - I0)))
  } else {
    Effort <- 0.5 * Data@MPeff[x] * (Irecent/I0)^2
  }
  
  Step <- (Effort/Data@MPeff[x])  # step change in effort 
  Step[Step < 0.85] <- 0.85
  Step[Step > 1.15] <- 1.15
  Allocate <- 1
  Effort <- Step * Data@MPeff[x]
  Effort[Effort < 0.01] <- 0.01  # for simulations in case Effort goes negative
  
  
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(op)
    par(mfrow=c(1,1))
    
    ylim <- range(c(Data@Ind[x, ], Itarget, I0))
    plot(Data@Year, Data@Ind[x, ], type="l", lwd=2, bty="l",
         xlab="Year", ylab="Index", ylim=ylim)
    points(max(Data@Year), mean(Data@Ind[x, ind]), cex=2, pch=16,col='blue')
    text(max(Data@Year), mean(Data@Ind[x, ind]), cex=1, 'Irecent', pos=3, col='blue', xpd=NA)
    
    lines(Data@Year[ind3], rep(mean(Data@Ind[x, ind3]), length(ind3)), lty=2, col="orange")
    text(mean(Data@Year[ind3]), mean(Data@Ind[x, ind3]), "Iave", col="orange", pos=1)
    
    points(max(Data@Year), Itarget, cex=2, pch=16,col='green')
    text(max(Data@Year), Itarget, cex=1, 'Itarget', pos=3, col='green', xpd=NA)
    
    points(max(Data@Year), I0, cex=2, pch=16,col='red')
    text(max(Data@Year), I0, cex=1, 'I0', pos=3, col='red', xpd=NA)
  
    
  }
  
  list(Effort=Effort )
}



#' Incremental Index Target MP - Effort-Based
#'  
#' A management procedure that incrementally adjusts the fishing effort 
#' to reach a target CPUE / relative abundance index
#' 
#' Four index/CPUE target MPs proposed by Geromont and Butterworth 2014. 
#' 
#' The TAE is calculated as:
#' If  \eqn{I_\textrm{recent} \geq I_0}:
#' \deqn{\textrm{TAE_y}= 0.5 \textrm{TAE}_{y-1} \left[1+\left(\frac{I_\textrm{recent} - I_0}{I_\textrm{target} - I_0}\right)\right]}
#' 
#' else:
#' \deqn{\textrm{TAE_y}= 0.5 \textrm{TAE}_{y-1} \left[\frac{I_\textrm{recent}}{I_0}^2\right]}
#' 
#' where \eqn{I_0} is \eqn{0.8 I_{\textrm{ave}}} (the average index over the 2 x `yrsmth` years prior to the projection period), 
#' \eqn{I_\textrm{recent}} is the average index over the past `yrsmth` years, and 
#' \eqn{I_\textrm{target}} is `Imulti` times \eqn{I_{\textrm{ave}}}.
#' 
#' @templateVar mp ItargetE1 
#' @template MPtemplate
#' @template MPuses
#' 
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param Imulti Parameter controlling how much larger target CPUE / index is
#' compared with recent levels.

#' @author T. Carruthers
#' @references Carruthers et al. 2015. Performance evaluation of simple
#' management procedures. ICES J. Mar Sci. 73, 464-482.
#' 
#' Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for
#' data-poor fisheries; forecasting with few data. ICES J. Mar. Sci. 72, 251-261.
#' doi:10.1093/icesjms/fst232
#' @describeIn ItargetE1 The less precautionary TAE-based MP
#' @family Index methods
#' @export 
#' @examples 
#' ItargetE1(1, DLMtool::Atlantic_mackerel, plot=TRUE)
ItargetE1 <- function(x, Data, reps = 100, plot=FALSE, yrsmth = 5, Imulti = 1.5) {
  
  runItargetE <- Itargeteff_(x, Data, reps, plot, yrsmth, Imulti)
  
  rec <- new("Rec")
  rec@Effort <- mean(runItargetE$Effort)
  rec
}
class(ItargetE1) <- "MP"


#' @describeIn ItargetE1 Increasing biologically precautionary TAE-based MP
#' @export 
#' @examples 
#' ItargetE2(1, DLMtool::Atlantic_mackerel, plot=TRUE)
ItargetE2 <- ItargetE1
formals(ItargetE2)$Imulti <- 2
class(ItargetE2) <- "MP"

#' @describeIn ItargetE1 Increasing biologically precautionary TAE-based MP
#' @export 
#' @examples 
#' ItargetE3(1, DLMtool::Atlantic_mackerel, plot=TRUE)
ItargetE3 <- ItargetE1
formals(ItargetE3)$Imulti <- 2.5
class(ItargetE3) <- "MP"

#' @describeIn ItargetE1 The most biologically precautionary TAE-based MP
#' @export 
#' @examples 
#' ItargetE4(1, DLMtool::Atlantic_mackerel, plot=TRUE)
ItargetE4 <- ItargetE1
formals(ItargetE4)$xx <- 0.3 
formals(ItargetE4)$Imulti <- 2.5
class(ItargetE4) <- "MP"


## TO DO ####



#' @describeIn LstepCC1 The least biologically precautionary effort-based MP.
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




#' @describeIn LstepCC1 The most precautionary effort-based MP.
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




#' @describeIn Ltarget1 The least biologically precautionary effort-based MP.
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




#' @describeIn Ltarget1 The most biologically precautionary effort-based MP.
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



#### Length-Based SPR ####


#' Internal Estimation Function for LBSPR MP 
#'
#' @param x Iteration number
#' @param Data An object of class `Data`
#' @param reps Number of samples. Not used in this method.
#' @param n Numeric. Number of historical years to run the model.
#' @param smoother Logical. Should estimates be smoothed over multiple years?
#'
#' @export
#' @references  
#' Hordyk, A., Ono, K., Valencia, S., loneragan, N., and Prince J; 
#' A novel length-based empirical estimation method of spawning potential ratio (SPR),
#' and tests of its performance, for small-scale, data-poor fisheries, 
#' ICES Journal of Marine Science, 72 (1) 2015, 217–231, 
#' https://doi.org/10.1093/icesjms/fsu004
#' 
#' @keywords internal
LBSPR_ <- function(x, Data, reps, n=5, smoother=TRUE) {
  if (NAor0(Data@L50[x])) stop("Data@L50 is NA")
  if (NAor0(Data@L95[x])) stop("Data@L95 is NA")
  if (NAor0(Data@wlb[x])) stop("Data@wlb is NA")
  
  LenBins <- Data@CAL_bins
  By <- LenBins[2] - LenBins[1]
  LenMids <- seq(from=By*0.5, by=By, length.out = length(LenBins)-1)
  CVLinf <- Data@LenCV[x]
  MK <- Data@Mort[x]/Data@vbK[x]
  Linf <- Data@vbLinf[x]
  L50 <- Data@L50[x]
  L95 <- Data@L95[x]
  Beta <- Data@wlb[x]
  
  n <- min(n, nrow(Data@CAL[x,,]))
  if (length(Data@Misc) ==0) Data@Misc[[x]] <- NA 
  if (length(Data@Misc[[x]]) ==0) Data@Misc[[x]] <- NA 
  if (is.null(Data@Misc[[x]]) || is.na(Data@Misc[[x]])) { # first time it's being run
    
    # run model for n most recent years 
    yind <- match(Data@LHYear, Data@Year)
    CALdata <- Data@CAL[x, (yind-n+1):length(Data@Year),]
    if (class(CALdata) == 'numeric')  CALdata <- matrix(CALdata, ncol=length(LenMids))
    Ests <- matrix(NA, nrow=nrow(CALdata), ncol=4)
    Fit <- list()
    for (y in 1:nrow(CALdata)) {
      CAL <- CALdata[y,]
      sl50start <- LenMids[which.max(CAL)]
      starts <- log(c(sl50start/Linf, sl50start/Linf*0.1, 1))
      runOpt <- optim(starts, LBSPRopt, CAL=CAL, nage=101, nlen=length(LenMids), CVLinf=CVLinf, 
                      LenBins=LenBins, LenMids=LenMids, x=seq(0, to=1, length.out = 101),
                      MK=MK, Linf=Linf, P=0.01, L50=L50, L95=L95, Beta=Beta)
      SL50 <- exp(runOpt$par[1]) * Linf 
      dSL50 <- exp(runOpt$par[2])
      SL95 <- SL50 + dSL50 * SL50
      FM <- exp(runOpt$par[3])
      
      runMod <- LBSPRgen(SL50, SL95, FM, nage=101, nlen=length(LenMids), CVLinf, 
                         LenBins, LenMids, x=seq(0, to=1, length.out = 101), 
                         MK, Linf, P=0.01, L50, L95, Beta)
      
      Ests[y,] <- c(SL50, SL95, FM, runMod[[2]])
      Fit[[y]] <- runMod[[1]] * sum(CALdata[y,])
    }
    
    if (smoother && nrow(Ests) > 1) Ests <- apply(Ests, 2, FilterSmooth)
      
    Ests <- as.data.frame(Ests)
    names(Ests) <- c("SL50", "SL95", "FM", "SPR")
    Ests$Year <- (yind-n+1):length(Data@Year)

  } else {
    lastYr <- max(Data@Misc[[x]]$Year)
    curYr <- max(Data@Year)
    yrs <- (lastYr+1):curYr
    
    CALdata <- Data@CAL[x, (length(Data@Year)-length(yrs)+1):length(Data@Year),]
    if (class(CALdata) == 'numeric')  CALdata <- matrix(CALdata, ncol=length(LenMids))
    Ests <- matrix(NA, nrow=nrow(CALdata), ncol=4)
    Fit <- list()
    for (y in 1:nrow(CALdata)) {
      CAL <- CALdata[y,]
      sl50start <- LenMids[which.max(CAL)]
      starts <- log(c(sl50start/Linf, sl50start/Linf*0.1, 1))
      runOpt <- optim(starts, LBSPRopt, CAL=CAL, nage=101, nlen=length(LenMids), CVLinf=CVLinf, 
                      LenBins=LenBins, LenMids=LenMids, x=seq(0, to=1, length.out = 101),
                      MK=MK, Linf=Linf, P=0.01, L50=L50, L95=L95, Beta=Beta)
      SL50 <- exp(runOpt$par[1]) * Linf 
      dSL50 <- exp(runOpt$par[2])
      SL95 <- SL50 + dSL50 * SL50
      FM <- exp(runOpt$par[3])
      
      runMod <- LBSPRgen(SL50, SL95, FM, nage=101, nlen=length(LenMids), CVLinf, 
                         LenBins, LenMids, x=seq(0, to=1, length.out = 101), 
                         MK, Linf, P=0.01, L50, L95, Beta)
      
      Ests[y,] <- c(SL50, SL95, FM, runMod[[2]])
      Fit[[y]] <- runMod[[1]] * sum(CALdata[y,])
    }
    Ests <- as.data.frame(Ests)
    names(Ests) <- c("SL50", "SL95", "FM", "SPR")
    Ests$Year <- (length(Data@Year)-length(yrs)+1):length(Data@Year)
    AllEsts <- rbind(Data@Misc[[x]], Ests)
    if (smoother) {
      SmoothEsts <- apply(AllEsts[,1:4], 2, FilterSmooth)
      AllEsts[,1:4] <- SmoothEsts
    }
    Ests <-AllEsts
  }
  
  
 return(list(Ests=Ests, Fit=Fit))
  
  
}


#' Length-Based SPR Effort Control
#' 
#' The spawning potential ratio (SPR) is estimated using the LBSPR method 
#' and compared to a target of 0.4.
#' 
#' Effort is increased by 10 per cent if the ratio of \eqn{\frac{\textm{SPR}}{\textrm{SPR}_{\textrm{targ}}}} is 
#' > 1.25, reduced by 10 per cent if the ratio is < 0.75, and remains unchanged 
#' otherwise.
#' 
#' The effort HCR has not been tuned. The increase/decrease in effort can
#' be adjusted using the `frac` argument.#' 
#' 
#' @templateVar mp LBSPR 
#' @template MPtemplate
#' @template MPuses 
#' 
#' @param n Last number of years to run the model on.
#' @param smoother Logical. Should the SPR estimates be smoothed?
#' 
#' @export
#' @references  
#' Hordyk, A., Ono, K., Valencia, S., loneragan, N., and Prince J; 
#' A novel length-based empirical estimation method of spawning potential ratio (SPR),
#' and tests of its performance, for small-scale, data-poor fisheries, 
#' ICES Journal of Marine Science, 72 (1) 2015, 217–231, 
#' https://doi.org/10.1093/icesjms/fsu004
#' @examples 
#' LBSPR(1, Data=DLMtool::SimulatedData, plot=TRUE)
LBSPR <- function(x, Data, reps=NA, plot=FALSE, n=5, smoother=TRUE, frac=0.1) {

  runLBSPR <- LBSPR_(x, Data, reps, n, smoother)
  
  Ests <- runLBSPR[[1]]
  estSPR <- Ests$SPR[length(Ests$SPR)]
  SPRtarg <- 0.4 
  ratio <- estSPR/SPRtarg
  if (ratio > 1.25) {
    Eff <- Data@MPeff[x] * (1 + frac)
  } else if (ratio < 0.75 ) {
    Eff <- Data@MPeff[x] * (1 - frac)
  } else {
    Eff <- Data@MPeff[x]
  }
  
  if (plot) {
  
    nyr <- length(runLBSPR$Fit)
    
    CAL <- Data@CAL[x,,]
    nyears <- dim(CAL)[1]
    CAL <- CAL[(nyears-nyr+1):nyears,]
    LenBins <- Data@CAL_bins
    By <- LenBins[2] - LenBins[1]
    LenMids <- seq(from=By*0.5, by=By, length.out = length(LenBins)-1)
    
    op <- par(no.readonly = TRUE)
    on.exit(op)
    nrow <- ceiling(sqrt(nyr))
    ncol <- ceiling((nyr + 1)/nrow)
    par(mfrow=c(nrow,ncol)) 
    
    if (nyr > 1) {
      ymin <- min(unlist(apply(CAL > 0, 1, which)))
      ymax <- max(unlist(apply(CAL > 0, 1, which)))
      ind <- (ymin-1):(ymax+1)
      ylim <- c(0, max(c(CAL, unlist(lapply(runLBSPR$Fit, max)))))
      for (p in 1:nyr) {
        tt <- barplot(CAL[p,ind], xlab="Length", ylab="Count", bty="l", names=LenMids[ind], ylim=ylim)
        lines(tt, runLBSPR$Fit[[p]][ind], lwd=2)
        title(paste0("Year ", runLBSPR$Ests$Year[p]))
      }
    } else {
        ymin <- min(which(CAL > 0))
        ymax <- max(which(CAL > 0))
        ind <- (ymin-1):(ymax+1)
        ylim <- c(0, max(c(CAL, runLBSPR$Fit[[1]])))
        tt <- barplot(CAL[ind], xlab="Length", ylab="Count", bty="l", names=LenMids[ind], ylim=ylim)
        lines(tt, runLBSPR$Fit[[1]][ind], lwd=2)
        title(paste0("Year ", runLBSPR$Ests$Year[p]))
    }
   
    plot(runLBSPR$Ests$Year, runLBSPR$Ests$SPR, ylim=c(0,1), xlab="Year", 
         ylab="SPR", type="b", las=1, bty="l")
    
  }
  Rec <- new("Rec")
  Rec@Effort <- Eff
  Rec@Misc <- Ests
  Rec
  
}
class(LBSPR) <- 'MP'

 



#' Kalman filter and Rauch-Tung-Striebel smoother
#'
#' A function that applies a filter and smoother to estimates
#'
#' @param RawEsts a vector of estimated values
#' @param R variance of sampling noise
#' @param Q variance of random walk increments
#' @param Int covariance of initial uncertainty
#' @return a vector of smoothed values
#'
#' @export
#' @keywords internal
FilterSmooth <- function(RawEsts, R=1, Q=0.1, Int=100) {
  # Modified from \url{"http://read.pudn.com/downloads88/ebook/336360/Kalman%20Filtering%20Theory%20and%20Practice,%20Using%20MATLAB/CHAPTER4/RTSvsKF.m__.htm"}
  Ppred <-  rep(Int, length(RawEsts))
  nNA <- sum(is.na(RawEsts))
  while(nNA > 0) { # NAs get replaced with last non-NA
    RawEsts[is.na(RawEsts)] <- RawEsts[which(is.na(RawEsts))-1]
    nNA <- sum(is.na(RawEsts))
  }
  Pcorr <- xcorr <- xpred <- rep(0, length(RawEsts))
  # Kalman Filter
  for (X in 1:length(Ppred)) {
    if (X !=1) {
      Ppred[X] <- Pcorr[X-1] + Q
      xpred[X] <- xcorr[X-1]
    }
    W <- Ppred[X]/(Ppred[X] + R)
    xcorr[X] <- xpred[X] + W * (RawEsts[X] - xpred[X]) # Kalman filter estimate
    Pcorr[X] <- Ppred[X] - W * Ppred[X]
  }
  # Smoother
  xsmooth <- xcorr
  for (X in (length(Pcorr)-1):1) {
    A <- Pcorr[X]/Ppred[X+1]
    xsmooth[X] <- xsmooth[X] + A*(xsmooth[X+1] - xpred[X+1])
  }
  return(xsmooth)
}
  
  
  
  
  
  
  
