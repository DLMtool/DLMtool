
## Output Control MPs ####



# ----- Average Catch ----

#' Average Catch
#'
#' A simple average catch MP that is included to demonstrate a 'status quo' management option
#'
#' The average catch method is very simple. The mean historical catch is calculated and used to set 
#' a constant catch limit (TAC). If `reps` > 1 then the `reps` samples are drawn from a log-normal
#' distribution with mean `TAC` and standard deviation (in log-space) of 0.2.    
#' 
#' For completeness, the TAC is calculated by:
#' 
#' \deqn{\textrm{TAC} =\frac{\sum_{y=1}^{\textrm{n}}{C_y}}{\textrm{n}}}
#' 
#' where \eqn{\textrm{TAC}} is the the mean catch recommendation, \eqn{n} is the number of historical years, and
#' \eqn{C_y}  is the catch in historical year \eqn{y}
#' 
#' @templateVar mp AvC
#' @template MPtemplate 
#' @template MPuses
#' @export 
#' @author T. Carruthers
#' 
#' @family Average Catch MPs
#' 
#' @examples 
#' Rec <- AvC(1, DLMtool::Cobia, reps=1000, plot=TRUE) # 1,000 log-normal samples with CV = 0.2
#'         
AvC <- function(x, Data, reps = 100, plot=FALSE) {
  dependencies = "Data@Cat Data@LHYear"
  yrs <- min(Data@Year):(Data@Year[Data@Year==Data@LHYear])
  yr.ind <- match(yrs, Data@Year)
  histCatch <- Data@Cat[x, yr.ind]
  meanC <- mean(histCatch, na.rm = T)
  if (reps >1) {
    TAC <- rlnorm(reps, log(meanC), 0.2)
  } else {
    TAC <- meanC
  }
  Rec <- new("Rec")
  Rec@TAC <- TAC
  
  if (plot) AvC_plot(x, Data, Rec, meanC, histCatch, yr.ind, lwd=3, cex.lab=1.25)
  Rec
}
class(AvC) <- "MP"

#### Beddington-Kirkwood Fmax estimation ####

#' Beddington and Kirkwood life-history MP
#' 
#' Family of management procedures that sets the TAC by approximation of Fmax
#' based on the length at first capture relative to asymptotic length and the 
#' von Bertalanffy growth parameter *K*.
#' 
#' @details The TAC is calculated as:
#' \deqn{\textrm{TAC} = A F_{\textrm{max}}}
#' where \eqn{A} is (vulnerable) stock abundance, and \eqn{F_{\textrm{max}}} is calculated as:
#' \deqn{F_{\textrm{max}} = \frac{0.6K}{0.67-L_c/L_\infty}}
#' where \eqn{K} is the von Bertalanffy growth coefficient, \eqn{L_c} is the
#' length at first capture, and \eqn{L_\infty} is the von Bertalanffy asymptotic length
#' 
#' 
#' Abundance (*A*) is either assumed known (`BK`) or estimated (`BK_CC` and `BK_ML`):
#' \deqn{A = \frac{\bar{C}}{\left(1-e^{-F}\right)}}
#' where \eqn{\bar{C}} is the mean catch, and *F* is estimated. 
#' See Functions section below for the estimation of *F*.
#' 
#' 
#' @note 
#' Note that the Beddington-Kirkwood method is designed to estimate \eqn{F_\textrm{max}},
#' that is, the fishing mortality that produces the maximum yield *assuming constant 
#' recruitment independent of spawning biomass*. 
#' 
#' Beddington and Kirkwood (2005) 
#' recommend estimating *F* using other methods (e.g a catch curve) and comparing the 
#' estimated *F* to the estimated \eqn{F_\textrm{max}} and adjusting exploitation accordingly. 
#' These MPs have not been implemented that way.
#' 
#' 
#' @describeIn BK Assumes that abundance is known, i.e. `Data@Abun`
#' and `Data@CV_abun` contain values
#' 
#' @examples 
#' BK(1, DLMtool::SimulatedData, reps=1000, plot=TRUE)
#' 
#' @author T. Carruthers.
#' @references Beddington, J.R., Kirkwood, G.P., 2005. The estimation of
#' potential yield and stock status using life history parameters. Philos.
#' Trans. R. Soc. Lond. B Biol. Sci. 360, 163-170.
#' @templateVar mp BK
#' @template MPtemplate
#' @template MPuses
#' @export 
#' 
#' 
BK <- function(x, Data, reps = 100, plot=FALSE) {
  # Beddington and Kirkwood life-history analysis
  dependencies = "Data@LFC, Data@vbLinf, Data@CV_vbLinf, Data@Abun, Data@CV_Abun, Data@vbK, Data@CV_vbK"
  Lc <- trlnorm(reps * 10, Data@LFC[x], 0.2)
  Linfc <- trlnorm(reps * 10, Data@vbLinf[x], Data@CV_vbLinf[x])
  Ac <- trlnorm(reps * 10, Data@Abun[x], Data@CV_Abun[x])
  Kc <- trlnorm(reps * 10, Data@vbK[x], Data@CV_vbK[x])
  Fmax <- (0.6 * Kc)/(0.67 - (Lc/Linfc)) 
  TAC <- Ac * Fmax 
  ind <- TAC > 0 & Fmax > 0
  TAC <- TACfilter(TAC[ind][1:reps])  # Interval censor only those positive catch recommendations
  Rec <- new("Rec")
  Rec@TAC <- TAC
  
  if (plot) {
    DF <- data.frame(vals=c(Lc[ind][1:reps]/Linfc[ind][1:reps], Kc[ind][1:reps], 
                            Fmax[ind][1:reps], Ac[ind][1:reps], TAC[ind][1:reps]), 
                     vars=rep(c('Lc/Linf', 'K', 'Fmax', 'Abundance', 'TAC'), each=reps),
                     stringsAsFactors = FALSE)
    DF <- DF[!is.na(DF$vals),]
    DF$vars <- factor(DF$vars, levels=c('Lc/Linf', 'K', 'Fmax', 'Abundance', 'TAC'), ordered=TRUE)
    BK_plot(DF)
  }
  
  Rec
  
} 
class(BK) <- "MP"

#' @param Fmin The minimum fishing mortality rate that is derived from the
#' catch-curve (interval censor).
#' @describeIn BK Abundance is estimated. An age-based catch curve is used 
#' to estimate *Z* and *F*. 
#' @templateVar mp BK_CC
#' @template MPuses
#' @examples 
#' BK_CC(1, DLMtool::SimulatedData, reps=1000, plot=TRUE)
#' 
#' @export 
BK_CC <- function(x, Data, reps = 100, plot=FALSE, Fmin = 0.005) {
  dependencies = "Data@LFC, Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@CAA, Data@Mort"
  Lc <- trlnorm(reps, Data@LFC[x], 0.2)
  Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  Mdb <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])
  MuC <- Data@Cat[x, length(Data@Cat[x, ])]
  Cc <- trlnorm(reps, MuC, Data@CV_Cat[x])
  Zdb <- CC(x, Data, reps = reps * 10)
  Fdb <- Zdb - Mdb
  Ac <- Cc/(1 - exp(-Fdb))
  Fmax <- (0.6 * Kc)/(0.67 - (Lc/Linfc))
  
  ind <- (1:(reps * 10))[Fdb > Fmin][(1:reps)*5]
  Fdb <- Fdb[ind]
  Mdb <- Mdb[ind]
  SM <- sum(is.na(ind))
  if (SM > 0) {
    Mdb[is.na(ind)] <- trlnorm(SM, Data@Mort[x], Data@CV_Mort[x])
    Fdb[is.na(ind)] <- Fmin
  }
  Ac <- Cc/(1 - exp(-Fdb))
  
  TAC <- Ac * Fmax
  ind <- TAC > 0 & Fmax > 0
  TAC <- TACfilter(TAC[ind][1:reps])
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  
  if (plot) {
    DF <- data.frame(vals=c(Lc[ind][1:reps]/Linfc[ind][1:reps], Kc[ind][1:reps], 
                            Fmax[ind][1:reps], Ac[ind][1:reps], TAC[ind][1:reps]), 
                     vars=rep(c('Lc/Linf', 'K', 'Fmax', 'Abundance', 'TAC'), each=reps),
                     stringsAsFactors = FALSE)
    DF <- DF[!is.na(DF$vals),]
    
    DF$vars <- factor(DF$vars, levels=c('Lc/Linf', 'K', 'Fmax', 'Abundance', 'TAC'), ordered=TRUE)
    BK_plot(DF)
    
  }
  Rec
}  
class(BK_CC) <- "MP"

#' @describeIn BK Abundance is estimated. *Z* and *F* are estimated from mean length.
#' @templateVar mp BK_ML
#' @template MPuses
#' @examples 
#' BK_ML(1, DLMtool::SimulatedData, reps=1000, plot=TRUE)
#' 
#' @export 
BK_ML <- function(x, Data, reps = 100, plot=FALSE) {
  dependencies = "Data@LFC, Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@CAL, Data@Mort"
  Lc <- trlnorm(reps * 10, Data@LFC[x], 0.2)
  Linfc <- trlnorm(reps * 10, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps * 10, Data@vbK[x], Data@CV_vbK[x])
  Mdb <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])
  Z <- MLne(x, Data, Linfc = Linfc, Kc = Kc, ML_reps = reps * 10, MLtype = "F")
  if (all(is.na(Z))) {
    Rec <- new("Rec")
    Rec@TAC <- TACfilter(rep(NA, reps))
    return(Rec)
  } 
  FM <- Z - Mdb
  MuC <- Data@Cat[x, length(Data@Cat[x, ])]
  Cc <- trlnorm(reps * 10, MuC, Data@CV_Cat[x])
  Ac <- Cc/(1 - exp(-FM))
  Fmax <- (0.6 * Kc)/(0.67 - (Lc/Linfc))  # robustifying for use in MSE 
  TAC <- Ac * Fmax
  ind <- TAC > 0 & Fmax > 0
  TAC <- TACfilter(TAC[ind][1:reps])
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  
  if (plot) {
    DF <- data.frame(vals=c(Lc[ind][1:reps]/Linfc[ind][1:reps], Kc[ind][1:reps], 
                            Fmax[ind][1:reps], Ac[ind][1:reps], TAC[ind][1:reps]), 
                     vars=rep(c('Lc/Linf', 'K', 'Fmax', 'Abundance', 'TAC'), each=reps),
                     stringsAsFactors = FALSE)
    DF <- DF[!is.na(DF$vals),]
    DF$vars <- factor(DF$vars, levels=c('Lc/Linf', 'K', 'Fmax', 'Abundance', 'TAC'), ordered=TRUE)
    BK_plot(DF)
  }
  
  Rec
}
class(BK_ML) <- "MP"








#' Geromont and Butterworth (2015) Constant Catch 
#' 
#' The TAC is the average historical catch over the last `yrsmth` (default 5) years, 
#' multiplied by (1-`xx`)
#' 
#' The TAC is calculated as:
#' \deqn{\textrm{TAC} = (1-x)C_{\textrm{ave}}}
#' where *x* lies between 0 and 1, and \eqn{C_{\textrm{ave}}} is average historical
#' catch over the previous `yrsmth` years. 
#' 
#' The TAC is constant for all future projections.
#' 
#' @param yrsmth Years over which to calculate mean catches
#' @param xx Parameter controlling the TAC. Mean catches are multiplied by
#' (1-`xx`)
#' 
#' @template MPtemplate
#' @templateVar mp CC1
#' @template MPuses
#' @export
#' 
#' @author T. Carruthers
#' @references Geromont, H. F., and D. S. Butterworth. 2015. 
#' Generic Management Procedures for Data-Poor Fisheries: Forecasting with Few Data.
#'  ICES Journal of Marine Science: Journal Du Conseil 72 (1). 251-61.
#'  
#' @describeIn CC1 TAC is average historical catch from recent `yrsmth` years
#' @examples 
#' CC1(1, DLMtool::Cobia, plot=TRUE)
#' 

#### Constant Catch MPs ####
CC1 <- function(x, Data, reps = 100, plot=FALSE, yrsmth = 5, xx = 0) {
  dependencies = "Data@Cat, Data@CV_Cat"
  
  yrlast <- match(Data@LHYear, Data@Year)
  yrfirst <- yrlast - yrsmth + 1
  # C_dat <- Data@Cat[x, (length(Data@Year) - (yrsmth - 1)):length(Data@Year)]
  C_dat <- Data@Cat[x, yrfirst:yrlast]
  TAC <- (1 - xx) * trlnorm(reps, mean(C_dat), Data@CV_Cat/(yrsmth^0.5))  # mean catches over the interval
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(mfrow=c(1,1))
    ylim <- c(0, max(c(Data@Cat[x,], TACfilter(TAC))))
    plot(c(Data@Year, max(Data@Year)+1), c(Data@Cat[x,],NA), type="l", lwd=2, las=1, bty="l",
         xlab="Year", ylab=paste0("Catch (", Data@Units, ")"),
         cex.lab=1.5, cex.axis=1.25, ylim=ylim)
    abline(v=Data@LHYear, lty=3, col="darkgray")
    lines(Data@Year[yrfirst:yrlast], rep(mean(C_dat), yrsmth), col="blue", lwd=3)
    boxplot(Rec@TAC, at=max(Data@Year)+1, add=TRUE, axes=FALSE)
  }
  
  Rec
}
class(CC1) <- "MP"

#' @describeIn CC1 TAC is average historical catch from recent `yrsmth` years reduced by 10\%. 
#' @examples 
#' CC2(1, DLMtool::Cobia, plot=TRUE)
#' 
#' @export 
CC2 <- CC1
formals(CC2)$xx <- 0.1
class(CC2) <- "MP"

#' @describeIn CC1 TAC is average historical catch from recent `yrsmth` years reduced by 20\%. 
#' @examples 
#' CC3(1, DLMtool::Cobia, plot=TRUE)
#' 
#' @export 
CC3 <- CC1
formals(CC3)$xx <- 0.2
class(CC3) <- "MP"


#' @describeIn CC1 TAC is average historical catch from recent `yrsmth` years reduced by 30\%. 
#' @examples 
#' CC4(1, DLMtool::Cobia, plot=TRUE)
#' 
#' @export 
CC4 <- CC1
formals(CC4)$xx <- 0.3
class(CC4) <- "MP"

#' @describeIn CC1 TAC is average historical catch from recent `yrsmth` years reduced by 40\%. 
#' @examples 
#' CC5(1, DLMtool::Cobia, plot=TRUE)
#' 
#' @export 
CC5 <- CC1
formals(CC5)$xx <- 0.4
class(CC5) <- "MP"

 

#' Geromont and Butterworth Constant Catch Harvest Control Rule
#' 
#' A simple MP that aims for average historical catches (as a proxy for MSY)
#' subject to imperfect information.
#' 
#' Note that this is my interpretation of their MP and is now stochastic.
#' Currently it is generalized and is not 'tuned' to more detailed assessment
#' data which might explain why in some cases it leads to stock declines.
#' 
#' @template MPtemplate
#' @templateVar mp GB_CC
#' @template MPuses
#' 
#' @author T. Carruthers
#' @references Geromont, H.F. and Butterworth, D.S. 2014. Complex assessment or
#' simple management procedures for efficient fisheries management: a
#' comparative study. ICES J. Mar. Sci. doi:10.1093/icesjms/fsu017
#' @export
#' @family Average Catch MPs
GB_CC <- function(x, Data, reps = 100, plot=FALSE) {
  dependencies = "Data@Cref,Data@Cat"
  Catrec <- Data@Cat[x, length(Data@Cat[x, ])]
  TAC <- trlnorm(reps, Data@Cref[x], Data@CV_Cref)
  TAC[TAC > (1.2 * Catrec)] <- 1.2 * Catrec
  TAC[TAC < (0.8 * Catrec)] <- 0.8 * Catrec
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(GB_CC) <- "MP"

## TO DO  - add plot and equations for GB_CC ####

#### Age-Comp SRA ####

#' Internal function for CompSRA MP
#'
#' @param x Simulation number
#' @param Data An object of class Data
#' @param reps Number of reps
#'
#' @return A list
#' @keywords internal
#' @export
#'
CompSRA_ <- function(x, Data, reps=100) {
  maxage <- Data@MaxAge
  TAC <- Bt_K <- FMSY <- Ac <- rep(NA, reps)
  predout <- list()
  for (i in 1:reps) {
    Mc <- trlnorm(1, Data@Mort[x], Data@CV_Mort)
    hc <- sample_steepness2(1, Data@steep[x], Data@CV_steep[x])
    Linfc <- trlnorm(1, Data@vbLinf[x], Data@CV_vbLinf[x])
    Kc <- trlnorm(1, Data@vbK[x], Data@CV_vbK[x])
    if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
      t0c <- -trlnorm(1, -Data@vbt0[x], Data@CV_vbt0[x])
    } else {
      t0c <- Data@vbt0[x]
    }
    t0c[!is.finite(t0c)] <- 0
    LFSc <- trlnorm(1, Data@LFS[x], Data@CV_LFS[x])
    LFCc <- trlnorm(1, Data@LFC[x], Data@CV_LFC[x])
    AMc <- trlnorm(1, iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], Data@L50[x]), Data@CV_L50[x])
    ac <- trlnorm(1, Data@wla[x], Data@CV_wla[x])
    bc <- trlnorm(1, Data@wlb[x], Data@CV_wlb[x])
    Catch <- Data@Cat[x, ]
    ny <- length(Catch)
    nyCAA <- dim(Data@CAA)[2]
    CAA <- Data@CAA[x, max(nyCAA - 2, 1):nyCAA, ]  # takes last three years as the sample (or last year if there is only one)
    
    Nac <- exp(-Mc * ((1:maxage) - 1))  # put a rough range on estimate of R0 assuming a mean harvest rate of 10%
    Lac <- Linfc * (1 - exp(-Kc * ((1:maxage) - t0c)))
    Wac <- ac * Lac^bc
    AFC <- log(1 - min(0.99, LFCc/Linfc))/-Kc + t0c
    AFS <- log(1 - min(0.99, LFSc/Linfc))/-Kc + t0c
    if (AFC >= 0.7 * maxage) AFC <- 0.7 * maxage
    if (AFS >= 0.9 * maxage) AFS <- 0.9 * maxage
    
    KES <- max(2, ceiling(mean(c(AFC, AFS))))
    pred <- Nac * Wac
    pred[1:(KES - 1)] <- 0
    pred <- pred/sum(pred)
    pred <- ((mean(Catch)/0.1) * pred/Wac)/exp(-(1:maxage) * Mc)
    pred <- pred[pred > 0]
    R0range <- c(mean(pred)/1000, mean(pred) * 1000)
    
    fit <- optimize(SRAfunc, log(R0range), Mc, hc, maxage, LFSc, LFCc, Linfc, Kc, t0c, AMc, ac, bc, Catch, CAA)
    getvals <- SRAfunc(fit$minimum, Mc, hc, maxage, LFSc, LFCc, Linfc, Kc, t0c, AMc, ac, bc, Catch, CAA, opt = 2)
    Ac[i] <- getvals$B
    Bt_K[i] <-  getvals$D
    predout[[i]] <- getvals$pred
    
    fit2 <- optimize(SRAFMSY, log(c(1e-04, 3)), Mc, hc, maxage, LFSc, LFCc, Linfc, Kc, t0c, AMc, ac, bc)
    FMSY[i] <- exp(fit2$minimum)
    if ((FMSY[i]/Mc) > 3) FMSY[i] <- 3 * Mc
    TAC[i] <- Ac[i] * FMSY[i]
  }
  
  return(list(TAC=TAC, Bt_K=Bt_K, FMSY=FMSY, Ac=Ac, pred=predout, CAA=CAA))
  
}

#' Age-Composition Stock-Reduction Analysis
#' 
#' A stock reduction analysis (SRA) model is fitted to the age-composition 
#' from the last 3 years (or less if fewer data are available)
#' 
#' A stock reduction analysis (SRA) model is fitted to the age-composition 
#' from the last 3 years (or less if fewer data are available) assuming a constant 
#' total mortality rate (*Z*) and used to estimate current stock depletion (*D*),
#'  *FMSY*, and stock abundance (*A*). 
#' 
#' Fishing mortality is estimated as:
#' \deqn{F = -\log{\left(1-\frac{C}{V}\right)}}
#' where \eqn{C} and \eqn{V} are catch and vulnerable biomass of first age of full
#' selection respectively. 
#' 
#' Abundance is estimated in the SRa. \eqn{F_{\textrm{MSY}}} is calculated assuming 
#' knife-edge vulnerability at the age of full selection. 
#'   
#' The TAC is calculated as \eqn{F_\textrm{MSY} A}. `CompSRA4010` uses a 40-10 
#' harvest control rule to reduce TAC at low biomass.
#' 
#' @templateVar mp CompSRA
#' @template MPtemplate
#' @template MPuses
#' 
#' @describeIn CompSRA TAC is FMSY x Abundance
#' @author T. Carruthers
#' 
#' @examples
#' CompSRA(1, DLMtool::SimulatedData, plot=TRUE)
#' 
#' @export 
CompSRA <- function(x, Data, reps = 100, plot=FALSE) {
  # optimize for fixed F to get you to current depletion C/Fcur =
  # abundance
  dependencies = "Data@Mort, Data@CV_Mort, Data@vbK, Data@CV_vbK, Data@vbLinf, Data@CV_vbLinf, Data@vbt0, Data@CV_vbt0, Data@MaxAge, Data@wla, Data@CV_wla, Data@wlb, Data@CV_wlb, Data@L50, Data@CV_L50, Data@CAA, Data@steep, Data@CV_steep, Data@LFS, Data@CV_LFS, Data@LFC, Data@CV_LFC, Data@Cat"
 
  runCompSRA <- CompSRA_(x, Data, reps)
  TAC <- TACfilter(runCompSRA$TAC)

  if (plot) CompSRA_plot(runCompSRA, TAC)
    
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}
class(CompSRA) <- "MP"


#' @templateVar mp CompSRA
#' @template MPuses
#' @describeIn CompSRA With a 40-10 control rule based on estimated depletion
#' @examples
#' CompSRA4010(1, DLMtool::SimulatedData, plot=TRUE)
#' 
#' @export CompSRA4010
CompSRA4010 <- function(x, Data, reps = 100, plot=FALSE) {
  dependencies = "Data@Mort, Data@CV_Mort, Data@vbK, Data@CV_vbK, Data@vbLinf, Data@CV_vbLinf, Data@vbt0, Data@CV_vbt0, Data@MaxAge, Data@wla, Data@CV_wla, Data@wlb, Data@CV_wlb, Data@L50, Data@CV_L50, Data@CAA, Data@steep, Data@CV_steep, Data@LFS, Data@CV_LFS, Data@LFC, Data@CV_LFC, Data@Cat"
  
  runCompSRA <- CompSRA_(x, Data, reps)
  
  TAC <- runCompSRA$TAC
  Bt_K <- runCompSRA$Bt_K

  # 40-10 rule
  cond1 <- Bt_K < 0.4 & Bt_K > 0.1
  cond2 <- Bt_K < 0.1
  TAC[cond1] <- TAC[cond1] * (Bt_K[cond1] - 0.1)/0.3
  TAC[cond2] <- TAC[cond2] * tiny  # this has to still be stochastic albeit very small
  
  if (plot) CompSRA_plot(runCompSRA, TAC)
  
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
  
}
class(CompSRA4010) <- "MP"


# ---- Depletion Corrected Average Catch MPs -----

#' Depletion Corrected Average Catch 
#' 
#' Internal code to calculate dcac
#'
#' @param x Simulation number
#' @param Data An object of class `Data`
#' @param reps The number of reps
#' @param Bt_K Optional value of depletion (otherwise taken from Data object)
#' @param updateD Logical. Should depletion be updated in projection years?
#'
#' @return A list with dcac, Bt_K, and BMSY_K
#' @keywords internal
#' @export
#'
DCAC_ <- function(x, Data, reps=100, Bt_K=NULL, updateD=FALSE) {
  
  yr.lst <- match(Data@LHYear, Data@Year)
  yrs <- 1:yr.lst
  C_tot <-  Data@AvC[x] * Data@t[x]
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
  FMSY_M <- trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])  # standard deviation of 0.2 - referred to as 'standard error' in MacCall 2009
  if (is.null(Bt_K))  Bt_K <- trlnorm(reps, Data@Dt[x], Data@CV_Dt[x])
  if (!updateD) {
    if (Data@LHYear != max(Data@Year)) {
      dcac <- rep(Data@MPrec[x], reps) # catch limit is static for future
      return(list(dcac=dcac, Bt_K=Bt_K))
    }
  }
  
  Bt_K[Bt_K>1] <-1
  Bt_K[Bt_K<0] <-0
  if (any(is.na(c(Data@BMSY_B0[x], Data@CV_BMSY_B0[x])))) {
    warning("Data@BMSY_B0 or Data@CV_BMSY_B0 do not contain values")
    return(list(dcac=rep(NA, reps), Bt_K=Bt_K))
  }
  BMSY_K <- rbeta(reps, alphaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * Data@CV_BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * Data@CV_BMSY_B0[x]))  
  
  dcac <- C_tot/(yr.lst + ((1 - Bt_K)/(BMSY_K * FMSY_M * Mdb)))
  
  return(list(dcac=dcac, Bt_K=Bt_K, BMSY_K=BMSY_K))
  
}


#' Depletion Corrected Average Catch
#' 
#' This group of MPs calculates a catch limit (*dcac*; intended as an MSY proxy) based on 
#' average historical catch while accounting for the windfall catch that got the stock down to
#' its current depletion level (*D*).
#'
#' The method calculates the depletion-corrected average catch (*dcac*) as:
#'
#' \deqn{\textrm{dcac} = \frac{\sum_{y=1}^{n}{C_y}}{n+(1-D)/Y_{\textrm{pot}}}}
#' 
#' where
#' \deqn{Y_{\textrm{pot}} = \frac{B_{\textrm{MSY}}}{B_0}\frac{F_{\textrm{MSY}}}{M}M} 
#' and \eqn{C} is the historical catches; i.e \eqn{C} does not change in the future projections in the MSE
#' 
#' The methods differ in the assumptions of current depletion (*D*). See the **Functions** section below for details.
#' 
#' @templateVar mp DCAC
#' @template MPtemplate
#' 
#' 
#' @note It's probably worth noting that DCAC TAC recommendations do not tend
#' to zero as depletion tends to zero. It adjusts for depletion only in
#' calculating historical average catch. It follows that at stock levels much
#' below BMSY, DCAC tends to chronically overfish.
#'
#' @references 
#' MacCall, A.D., 2009. Depletion-corrected average catch: a simple
#' formula for estimating sustainable yields in data-poor situations. ICES J.
#' Mar. Sci. 66, 2267-2271.
#' 
#' @author T. Carruthers
#' @family Average Catch MPs
#' @describeIn DCAC Depletion is not updated in the future projections. 
#' The TAC is static and not updated in the future years. This represents 
#' an application of the DCAC method where a catch limit is calculated based on 
#' current estimate of depletion and time-series of catch from the beginning
#' of the fishery, and the TAC is fixed at this level for all future projections.
#' 
#' @templateVar mp DCACs
#' @template MPuses
#' @examples 
#' DCACs(1, DLMtool::Atlantic_mackerel, plot=TRUE) 
#' 
#' @export 
DCACs <- function(x, Data, reps = 100, plot=FALSE) {
  rundcac <- DCAC_(x, Data, reps)
  TAC <- TACfilter(rundcac$dcac)
  if (plot)  DCAC_plot(x, Data, dcac=rundcac$dcac, TAC, Bt_K=rundcac$Bt_K, yrs=1:length(Data@Year))
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}  
class(DCACs) <- "MP"


#' @describeIn DCAC Depletion is estimated each management interval and used to 
#' update the catch limit recommendation based on the historical catch (which is 
#' not updated in the future projections).
#' @templateVar mp DCAC
#' @template MPuses
#' @examples 
#' DCAC(1, DLMtool::Atlantic_mackerel, plot=TRUE) 
#' 
#' @export 
DCAC <- function(x, Data, reps = 100, plot=FALSE) {
  "dependenciesData@AvC, Data@t, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@Dt, Data@CV_Dt, Data@BMSY_B0, Data@CV_BMSY_B0"
  
  rundcac <- DCAC_(x, Data, reps, updateD=TRUE)
  
  TAC <- TACfilter(rundcac$dcac)
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}  
class(DCAC) <- "MP"



#' @describeIn DCAC Current stock biomass is assumed to be exactly at 40 per cent 
#' of unfished levels. The 40 percent depletion assumption may not 
#' really affect DCAC that much as it already makes TAC recommendations that are 
#' quite MSY-like.
#' @templateVar mp DCAC_40
#' @template MPuses
#' @examples 
#' DCAC_40(1, DLMtool::Atlantic_mackerel, plot=TRUE) 
#' 
#' @export 
DCAC_40 <- function(x, Data, reps = 100, plot=FALSE) {
  dependencies = "Data@AvC, Data@t, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@BMSY_B0, Data@CV_BMSY_B0"
  
  rundcac <- DCAC_(x, Data, reps, Bt_K=0.4)
  TAC <- TACfilter(rundcac$dcac)
  if (plot)  DCAC_plot(x, Data, dcac=rundcac$dcac, TAC, Bt_K=rundcac$Bt_K, yrs=1:length(Data@Year))
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}  
class(DCAC_40) <- "MP"


#' @describeIn DCAC The dynamic DCAC (depletion is updated) is paired with the 40-10 rule that throttles 
#' back the OFL to zero at 10 percent of unfished stock size (the OFL is not subject to downward
#' adjustment above 40 percent unfished). DCAC can overfish below BMSY levels. The 40-10 
#' harvest control rule largely resolves this problem providing an MP with surprisingly good
#' performance even at low stock levels.
#' @templateVar mp DCAC4010
#' @template MPuses
#' @examples 
#' Data <- DLMtool::Atlantic_mackerel
#' Data@LHYear <- 2005
#' DCAC4010(1, Data, plot=TRUE) 
#' 
#' @export 
DCAC4010 <- function(x, Data, reps = 100, plot=FALSE) {
  dependencies = "Data@AvC, Data@t, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@Dt, Data@CV_Dt, Data@BMSY_B0, Data@CV_BMSY_B0"
  rundcac <- DCAC_(x, Data, reps, updateD=TRUE)
  dcac <- rundcac$dcac 
  Bt_K <- rundcac$Bt_K 
  TAC <- dcac
  
  # 40-10 rule
  cond1 <- Bt_K < 0.4 & Bt_K > 0.1
  cond2 <- Bt_K < 0.1
  if (length(cond1) > 0) TAC[cond1] <- TAC[cond1] * (Bt_K[cond1] - 0.1)/0.3
  if (length(cond2) > 0)  TAC[cond2] <- TAC[cond2] * tiny  # this has to still be stochastic albeit very small
  if (length(cond1) < 1 & length(cond2) < 1)  return(NA)
  TAC <- TACfilter(TAC)
  
  if (plot)  DCAC_plot(x, Data, dcac=rundcac$dcac, TAC, Bt_K=rundcac$Bt_K, yrs=1:length(Data@Year))
  
  Rec <- new("Rec")
  Rec@TAC <- TAC 
  Rec
}  
class(DCAC4010) <- "MP"



#' @describeIn DCAC This variant uses the mean length estimator to calculate current stock
#' depletion. The mean length extension was programmed by Gary Nelson as part of his
#' excellent R package 'fishmethods'.
#' @templateVar mp DCAC_ML
#' @template MPuses
#' @examples 
#' DCAC_ML(1, DLMtool::SimulatedData, plot=TRUE) 
#' 
#' @export 
DCAC_ML <- function(x, Data, reps = 100, plot=FALSE) {
  dependencies = "Data@AvC, Data@t, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@BMSY_B0, Data@CV_BMSY_B0, Data@Year, Data@CAL, Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK"
  if (is.na(Data@BMSY_B0[x]) | is.na(Data@CV_BMSY_B0[x])) return(NA)
  if (is.na(Data@FMSY_M[x]) | is.na(Data@CV_FMSY_M[x])) return(NA)
  
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  
  FMSY_M <- trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])
  Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  Z <- MLne(x, Data, Linfc = Linfc, Kc = Kc, ML_reps = reps, MLtype = "dep")
  if (all(is.na(Z))) {
    Rec <- new("Rec")
    Rec@TAC <- TACfilter(rep(NA, reps))
    return(Rec)
  }
  FM <- Z - Mdb
  nyears <- length(Data@Year)
  Ct1 <- mean(Data@Cat[x, 1:3])
  Ct2 <- mean(Data@Cat[x, (nyears - 2):nyears])
  dep <- rep(c(Ct1, Ct2), each = reps)/(1 - exp(-FM))
  if (reps == 1)Bt_K <- dep[2]/dep[1]
  if (reps > 1) Bt_K <- dep[, 2]/dep[, 1]
  
  rundcac <- DCAC_(x, Data, reps, Bt_K=Bt_K)
  TAC <- TACfilter(rundcac$dcac)
  
  if (plot)  DCAC_plot(x, Data, dcac=rundcac$dcac, TAC, Bt_K=rundcac$Bt_K, yrs=1:length(Data@Year))
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
} 
class(DCAC_ML) <- "MP"


#' @describeIn DCAC Depletion Adjusted Average Catch: essentially DCAC (with updated Depletion)
#'  divided by BMSY/B0 (Bpeak) (Harford and Carruthers, 2017).
#' @references
#' Harford W. and Carruthers, T. 2017. Interim and long-term performance of 
#' static and adaptive management procedures. Fish. Res. 190, 84-94.
#' @templateVar mp DAAC
#' @template MPuses
#' 
#' @examples
#' Data <- DLMtool::Atlantic_mackerel
#' Data@LHYear <- 2005
#' DAAC(1, Data, plot=TRUE) 
#' 
#' @export 
DAAC <- function(x, Data, reps = 100, plot=FALSE) {
  # extended depletion-corrected average catch (Harford and Carruthers 2017)
  dependencies = "Data@AvC, Data@t, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@Dt, Data@CV_Dt, Data@BMSY_B0, Data@CV_BMSY_B0"
  
  rundcac <- DCAC_(x, Data, reps, updateD = TRUE)
  TAC <- rundcac$dcac * rundcac$Bt_K/rundcac$BMSY_K
  
  TAC <- TACfilter(TAC)
  if (plot)  DCAC_plot(x, Data, dcac=rundcac$dcac, TAC, Bt_K=rundcac$Bt_K, yrs=1:length(Data@Year))
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}
class(DAAC) <- "MP"



#' @describeIn DCAC Hybrid Depletion Adjusted Average Catch: essentially DCAC  (with updated Depletion) 
#' divided by BMSY/B0 (Bpeak) when below BMSY, 
#' and DCAC above BMSY (Harford and Carruthers 2017).
#' @templateVar mp HDAAC
#' @template MPuses
#' @examples
#' Data <- DLMtool::Atlantic_mackerel
#' Data@LHYear <- 2005
#' HDAAC(1, Data, plot=TRUE) 
#' 
#' @export 
HDAAC <- function(x, Data, reps = 100, plot=FALSE) {
  dependencies = "Data@AvC, Data@t, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@Dt, Data@CV_Dt, Data@BMSY_B0, Data@CV_BMSY_B0"
  
  rundcac <- DCAC_(x, Data, reps, updateD = TRUE)
  TAC <- rundcac$dcac
  ddcac <- rundcac$dcac * rundcac$Bt_K/rundcac$BMSY_K
  TAC[rundcac$Bt_K < rundcac$BMSY_K] <- ddcac[rundcac$Bt_K < rundcac$BMSY_K]
  
  TAC <- TACfilter(TAC)
  
  if (plot)  DCAC_plot(x, Data, dcac=rundcac$dcac, TAC, Bt_K=rundcac$Bt_K, yrs=1:length(Data@Year))
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
  
}
class(HDAAC) <- "MP"


#### Depletion Based SRA ####

#' Internal optimization function
#'
#' @param lnK log unfished biomass
#' @param C_hist Historical catches
#' @param nys number of years (length of C_hist)
#' @param Mdb M 
#' @param FMSY_M  FMSY/M
#' @param BMSY_K BMSY/K
#' @param Bt_K Current depletion
#' @param adelay age of maturity
#' @param opt Numeric. 1 for optimization. 2 to return Biomass trend
#'
#' @keywords internal
#' @export
#'
DBSRAopt <- function(lnK, C_hist, nys, Mdb, FMSY_M, BMSY_K, Bt_K, adelay, opt=1) {
  # the optimization for B0 given DBSRA assumptions
  Kc <- exp(lnK)
  fn <- function(n, BMSY_K) {
    # optimizer to find parameter n according to sampled BMSY/B0 (theta)
    thetapred <- n^(-1/(n - 1))
    (BMSY_K - thetapred)^2
  }
  n <- optimize(fn, c(0.01, 6), BMSY_K = BMSY_K)$minimum  #get the optimum getn(BMSY_K)
  g <-  (n^(n/(n - 1)))/(n - 1)# gety(n) get the y parameter for n
  FMSY <- FMSY_M * Mdb
  UMSY <- (FMSY/(FMSY + Mdb)) * (1 - exp(-(FMSY + Mdb)))
  MSY <- Kc * BMSY_K * UMSY
  # Bjoin rules from Dick & MacCall 2011
  Bjoin_K <- 0.5
  if (BMSY_K < 0.3) Bjoin_K <- 0.5 * BMSY_K
  if (BMSY_K > 0.3 & BMSY_K < 0.5) Bjoin_K <- 0.75 * BMSY_K - 0.075
  Bjoin <- Bjoin_K * Kc
  prodPTF <- function(depletion, n, MSY) {
    # Pella-Tomlinson production function required for DB-SRA
    y <- (n^(n/(n - 1)))/(n - 1)
    MSY * y * depletion - MSY * y * depletion^n
  }
  PBjoin <- prodPTF(Bjoin_K, n, MSY)
  cp <- (1 - n) * g * MSY * (Bjoin^(n - 2)) * Kc^-n
  Bc <- rep(NA, nys)
  Bc[1] <- Kc
  obj <- 0
  for (yr in 2:nys) {
    yref <- max(1, yr - adelay)
    if (Bc[yref] > Bjoin | BMSY_K > 0.5) {
      Bc[yr] <- Bc[yr - 1] + g * MSY * (Bc[yref]/Kc) - g * MSY *
        (Bc[yref]/Kc)^n - C_hist[yr - 1]
    } else {
      Bc[yr] <- Bc[yr - 1] + Bc[yref] * ((PBjoin/Bjoin) + cp * (Bc[yref] - Bjoin)) - C_hist[yr - 1]
    }
    if (Bc[yr] < 0)
      obj <- obj + log(-Bc[yr])
    Bc[yr] <- max(1e-06, Bc[yr])
  }
  if (opt ==1) return(obj + ((Bc[nys]/Kc) - Bt_K)^2)
  if (opt ==2) return(Bc)
  
}  # end of DBSRA optimization function

#' Depletion-based SRA internal function
#' 
#' @param x iteration number
#' @param Data object of class Data 
#' @param reps number of reps
#' @param depo Optional fixed depletion (single value)
#' @param hcr Optional harvest control rule for throttling catch as a function of B/B0. 
#' Numeric vector of length 2 specifying HCR break points - e.g c(0.4, 0.1) for 40-10 HCR
#'
#' @export
#'
#' @keywords internal
#' 
DBSRA_ <- function(x, Data, reps = 100, depo=NULL, hcr=NULL) {
  # returns a vector of DBSRA estimates of the TAC for a particular
  # simulation x for(x in 1:nsim){
  C_hist <- Data@Cat[x, ]
  TAC <- rep(NA, reps)
  Btrend <- matrix(NA, nrow=reps, ncol=length(C_hist))
  Bt_Kstore <- FMSY_Mstore <- BMSY_K_Mstore <- rep(NA, reps)
  DBSRAcount <- 1
  if (is.null(depo)) {
    if (is.na(Data@Dep[x]) | is.na(Data@CV_Dep[x])) return(new("Rec"))
  } else {
    Data@CV_Dep[x] <- tiny
  }
  if (is.na(Data@BMSY_B0[x]))return(new("Rec"))
    
  while (DBSRAcount < (reps + 1)) {
    if (is.null(depo)) depo <- max(0.01, min(0.99, Data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
    Bt_K <- rbeta(100, alphaconv(depo, min(depo * Data@CV_Dep[x], (1 - depo) * Data@CV_Dep[x])), 
                  betaconv(depo, min(depo * Data@CV_Dep[x], (1 - depo) * Data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
    Bt_K <- Bt_K[Bt_K > 0.00999 & Bt_K < 0.99001][1]  # interval censor (0.01,0.99)  as in Dick and MacCall 2011
    
    Mdb <- trlnorm(100, Data@Mort[x], Data@CV_Mort[x])
    Mdb <- Mdb[Mdb < 0.9][1]  # !!!! maximum M is 0.9   interval censor
    if (is.na(Mdb)) Mdb <- 0.9  # !!!! maximum M is 0.9   absolute limit
    FMSY_M <- trlnorm(1, Data@FMSY_M[x], Data@CV_FMSY_M[x])
    BMSY_K <- rbeta(100, alphaconv(Data@BMSY_B0[x], Data@CV_BMSY_B0[x] *  Data@BMSY_B0[x]), 
                    betaconv(Data@BMSY_B0[x], Data@CV_BMSY_B0[x] * Data@BMSY_B0[x]))  #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
    tryBMSY_K <- BMSY_K[BMSY_K > 0.05 & BMSY_K < 0.95][1]  # interval censor (0.05,0.95) as in Dick and MacCall, 2011
    if (is.na(tryBMSY_K)) {
      Min <- min(BMSY_K, na.rm = TRUE)
      Max <- max(BMSY_K, na.rm = TRUE)
      if (Max <= 0.05) BMSY_K <- 0.05
      if (Min >= 0.95) BMSY_K <- 0.95
    }
    if (!is.na(tryBMSY_K))  BMSY_K <- tryBMSY_K
    
    Bt_Kstore[DBSRAcount] <- Bt_K
    FMSY_Mstore[DBSRAcount] <- FMSY_M
    BMSY_K_Mstore[DBSRAcount] <- BMSY_K
    
    adelay <- max(floor(iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],  Data@L50[x])), 1)
    
    # scale catches for optimization
    scaler <- 1000/mean(C_hist)
    C_hist2 <- scaler * C_hist
    opt <- optimize(DBSRAopt, log(c(0.01 * mean(C_hist2), 1000 * mean(C_hist2))), C_hist = C_hist2, 
                    nys = length(C_hist2), Mdb = Mdb, FMSY_M = FMSY_M, BMSY_K = BMSY_K, 
                    Bt_K = Bt_K, adelay = adelay, tol = 0.01)
    
    Bctemp <- DBSRAopt(opt$minimum,  C_hist = C_hist2, 
                   nys = length(C_hist2), Mdb = Mdb, FMSY_M = FMSY_M, BMSY_K = BMSY_K, 
                   Bt_K = Bt_K, adelay = adelay, opt=2)
   
    Btrend[DBSRAcount,] <- Bctemp/scaler
    Kc <- exp(opt$minimum) / scaler
    BMSYc <- Kc * BMSY_K
    FMSYc <- Mdb * FMSY_M
    UMSYc <- (FMSYc/(FMSYc + Mdb)) * (1 - exp(-(FMSYc + Mdb)))
    MSYc <- Kc * BMSY_K * UMSYc
    TAC[DBSRAcount] <- UMSYc * Kc * Bt_K
    
    if(!is.null(hcr)) {
      if (length(hcr)!=2) stop("hcr must be numeric vector of length 2")
      # 40-10 rule
      if (Bt_K < hcr[1] & Bt_K > hcr[2])  TAC[DBSRAcount] <- TAC[DBSRAcount] * (Bt_K - hcr[2])/(hcr[1]-hcr[2])
      if (Bt_K < hcr[2]) TAC[DBSRAcount] <- TAC[DBSRAcount] * tiny  # this has to still be a numeric value, albeit very small
    }
    DBSRAcount <- DBSRAcount + 1

  }  # end of reps
  list(TAC=TAC, Btrend=Btrend, C_hist=C_hist, Bt_Kstore=Bt_Kstore, FMSY_Mstore=FMSY_Mstore, 
       BMSY_K_Mstore=BMSY_K_Mstore, hcr=hcr)
  
}  


#' Depletion-Based Stock Reduction Analysis 
#' 
#' Depletion-Based Stock Reduction Analysis (DB-SRA) is a method designed for 
#' determining a catch limit and management reference points for data-limited 
#' fisheries where catches are known from the beginning of exploitation. 
#' User prescribed BMSY/B0, M, FMSY/M are used to find B0 and therefore the a catch limit
#' by back-constructing the stock to match a user specified level of stock
#' depletion.
#' 
#' 
#' DB-SRA assumes that a complete time-series of catch from the beginning of 
#' exploitation is available. Users prescribe estimates of current depletion \eqn{(D)}, 
#' biomass at MSY relative to unfished \eqn{\left(\frac{B_\textrm{MSY}}{B_0}\right)}, 
#' the natural mortality rate \eqn{(M)}, and the ratio fishing mortality at 
#' MSY to M \eqn{\left(\frac{F_{\textrm{MSY}}}{M}\right)}.
#' 
#' 
#' You may have noticed that you -the user- specify three of the factors that
#' make the quota recommendation. So this can be quite a subjective method.
#' In the MSE the MSY reference points (e.g., \eqn{\left(\frac{F_\textrm{MSY}}{M}\right)})
#' are taken as the true value calculate in the MSE with added uncertainty specified
#' in the Obs object (e.g `Obs@FMSY_Mbiascv`).
#' 
#' The catch limit, for the Base Version, is calculated as:
#' \deqn{\textrm{TAC} = M . \frac{F_{\textrm{MSY}}}{M} . D . B_0}
#' 
#' 
#' @templateVar mp DBSRA
#' @template MPtemplate
#' @template MPuses
#' 
#' @note The DB-SRA method of this package isn't exactly the same as the
#' original method of Dick and MacCall (2011) because it has to work for
#' simulated depletions above BMSY/B0 and even on occasion over B0. It also  
#' doesn't have the modification for flatfish life histories that has
#' previously been applied by Dick and MacCall (2011).
#' 
#' @author T. Carruthers
#' @references  
#' Dick, E.J., MacCall, A.D., 2010. Estimates of sustainable yield for 50 data-poor 
#' stocks in the Pacific Coast groundfish fishery management plan. Technical memorandum. 
#' Southwest fisheries Science Centre, Santa Cruz, CA. National Marine Fisheries Service, 
#' National Oceanic and Atmospheric Administration of the U.S. Department of Commerce. 
#' NOAA-TM-NMFS-SWFSC-460.
#' 
#' Dick, E.J., MacCall, A.D., 2011. Depletion-Based Stock Reduction
#' Analysis: A catch-based method for determining sustainable yields for
#' data-poor fish stocks. Fish. Res. 110, 331-341.
#'
#' @describeIn DBSRA Base Version. TAC is calculated assumed MSY harvest rate 
#' multiplied by the estimated current abundance (estimated B0 x Depletion)
#' @examples 
#' DBSRA(1, DLMtool::ourReefFish, plot=TRUE)
#' 
#' @export 
DBSRA <- function(x, Data, reps = 100, plot=FALSE) {
  runDBSRA <- DBSRA_(x,Data, reps)
  TAC <- TACfilter(runDBSRA[[1]])
  if (plot) DBSRA_plot(runDBSRA, Data, TAC)

  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
} 
class(DBSRA) <- "MP"


#' @templateVar mp DBSRA_40
#' @template MPuses
#' @describeIn DBSRA Same as the Base Version but assumes 40 percent current depletion (Bcurrent/B0 = 0.4), which is 
#' more or less the most optimistic state for a stock (ie very close to BMSY/B0 for many stocks).
#' @examples 
#' DBSRA_40(1, DLMtool::ourReefFish, plot=TRUE)
#' @export 
DBSRA_40 <- function(x, Data, reps = 100, plot=FALSE) {
  runDBSRA <- DBSRA_(x,Data, reps, depo=0.4)
  TAC <- TACfilter(runDBSRA[[1]])
  if (plot)  DBSRA_plot(runDBSRA, Data, TAC)
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}  # end of DBSRA_apply
class(DBSRA_40) <- "MP"


#' @templateVar mp DBSRA4010
#' @template MPuses
#' @describeIn DBSRA Base version paired with the 40-10 rule that linearly throttles
#' back the TAC when depletion is below 0.4 down to zero at 10 percent of unfished biomass.
#' @examples 
#' DBSRA4010(1, DLMtool::ourReefFish, plot=TRUE)
#' @export 
#' 
DBSRA4010 <- function(x, Data, reps = 100, plot=FALSE) {
  runDBSRA <- DBSRA_(x,Data, reps, hcr=c(0.4, 0.1))
  TAC <- TACfilter(runDBSRA[[1]])
  if (plot) DBSRA_plot(runDBSRA, Data, TAC)

  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}  
class(DBSRA4010) <- "MP"



# #' Depletion-Based Stock Reduction Analysis using mean length estimator of
# #' stock depletion
# #' 
# #' DBSRA using the mean length estimator to calculate current stock depletion.
# #' 
# #'
# #' @param x A position in a data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of stochastic samples of the quota recommendation
# #' @note The mean length extension was programmed by Gary Nelson as part of his
# #' excellent R package 'fishmethods'
# #' @author T. Carruthers
# #' @references Dick, E.J., MacCall, A.D., 2011. Depletion-Based Stock Reduction
# #' Analysis: A catch-based method for determining sustainable yields for
# #' data-poor fish stocks. Fish. Res. 110, 331-341.
# #' @export 
# DBSRA_ML <- function(x, Data, reps = 100) {
#   dependencies = "Data@Cat, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@BMSY_B0, Data@CV_BMSY_B0, Data@L50, Data@CAL, Data@Year, Data@Cat"
#   C_hist <- Data@Cat[x, ]
#   TAC <- rep(NA, reps)
#   DBSRAcount <- 1
#   maxIts <- 200
#   nIts <- 0
#   if (is.na(Data@Dep[x]) | is.na(Data@CV_Dep[x])) return(NA)
#   while (DBSRAcount < (reps + 1) & nIts < maxIts) {
#     Linfc <- trlnorm(1, Data@vbLinf[x], Data@CV_vbLinf[x])
#     Kc <- trlnorm(1, Data@vbK[x], Data@CV_vbK[x])
#     Mdb <- trlnorm(100, Data@Mort[x], Data@CV_Mort[x])
#     Mdb <- Mdb[Mdb < 0.9][1]  # !!!! maximum M is 0.9   interval censor
#     if (is.na(Mdb)) Mdb <- 0.9  # !!!! maximum M is 0.9   absolute limit
#     Z <- MLne(x, Data, Linfc = Linfc, Kc = Kc, ML_reps = 1, MLtype = "dep")
#     if (all(is.na(Z))) {
#       Rec <- new("Rec")
#       Rec@TAC <- TACfilter(rep(NA, reps))
#       return(Rec)
#     } 
#     FM <- Z - Mdb
#     FM[FM < 0] <- 0.01
#     nyears <- length(Data@Year)
#     Ct1 <- mean(Data@Cat[x, 1:3])
#     Ct2 <- mean(Data@Cat[x, (nyears - 2):nyears])
#     # dep<-c(Ct1,Ct2)/(1-exp(-FM[,c(1,2)]))
#     dep <- c(Ct1, Ct2)/(1 - exp(-FM))
#     Bt_K <- dep[2]/dep[1]
#     
#     if (Bt_K < 0.01) Bt_K <- 0.01  # interval censor / temporary hack to avoid doing multiple depletion estimates that would take far too long
#     if (Bt_K > 0.99) Bt_K <- 0.99  # interval censor / temporary hack to avoid doing multiple depletion estimates that would take far too long
#     
#     
#     FMSY_M <- trlnorm(1, Data@FMSY_M[x], Data@CV_FMSY_M[x])
#     BMSY_K <- rbeta(100, alphaconv(Data@BMSY_B0[x], Data@CV_BMSY_B0[x] * 
#                                      Data@BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@CV_BMSY_B0[x] * 
#                                                                   Data@BMSY_B0[x]))  #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
#     tryBMSY_K <- BMSY_K[BMSY_K > 0.05 & BMSY_K < 0.95][1]  # interval censor (0.05,0.95) as in Dick and MacCall, 2011
#     
#     if (is.na(tryBMSY_K)) {
#       Min <- min(BMSY_K, na.rm = TRUE)
#       Max <- max(BMSY_K, na.rm = TRUE)
#       if (Max <= 0.05) 
#         BMSY_K <- 0.05
#       if (Min >= 0.95) 
#         BMSY_K <- 0.95
#     }
#     if (!is.na(tryBMSY_K))  BMSY_K <- tryBMSY_K
#     if (all(is.na(BMSY_K))) return(rep(NA, reps))
#     adelay <- max(floor(iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],  Data@L50[x])), 1)
#     opt <- optimize(DBSRAopt, log(c(0.1 * mean(C_hist), 1000 * mean(C_hist))), 
#                     C_hist = C_hist, nys = length(C_hist), Mdb = Mdb, FMSY_M = FMSY_M, 
#                     BMSY_K = BMSY_K, Bt_K = Bt_K, adelay = adelay, tol = 0.01)
#     nIts <- nIts + 1
#     if (opt$objective < 0.1) {
#       Kc <- exp(opt$minimum)
#       BMSYc <- Kc * BMSY_K
#       FMSYc <- Mdb * FMSY_M
#       UMSYc <- (FMSYc/(FMSYc + Mdb)) * (1 - exp(-(FMSYc + Mdb)))
#       MSYc <- Kc * BMSY_K * UMSYc
#       TAC[DBSRAcount] <- UMSYc * Kc * Bt_K
#       DBSRAcount <- DBSRAcount + 1
#       nIts <- 0
#     }
#   }  # end of reps
#   Rec <- new("Rec")
#   Rec@TAC <- TACfilter(TAC)
#   Rec
# }
# class(DBSRA_ML) <- "MP"
# 




## Delay-Difference MP #### 

#' Delay-Difference Internal Function
#'
#' @param x iteration number
#' @param Data object of class Data 
#' @param reps number of reps
#' @param hcr Optional harvest control rule for throttling catch as a function of B/B0. 
#' Numeric vector of length 2 specifying HCR break points - e.g c(0.4, 0.1) for 40-10 HCR
#' @export
#' @keywords internal
DD_ <- function(x, Data, reps = 100, hcr=NULL) {
  Winf <- Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
  age <- 1:Data@MaxAge
  la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
  wa <- Data@wla[x] * la^Data@wlb[x]
  a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],  Data@L50[x]) # assume vulnerability = maturity
  a50V <- max(a50V, 1)
  
  yind <- which(!is.na(Data@Cat[x, ] + Data@Ind[x, ]))[1] # First year with both catch and index
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
  k_DD <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)  
  k_DD[k_DD > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho_DD <- (wa[k_DD + 2] - Winf)/(wa[k_DD + 1] - Winf) # Brody growth coefficient
  Alpha_DD <- Winf * (1 - Rho_DD)
  So_DD <- exp(-Data@Mort[x])  # get So survival rate
  wa_DD <- wa[k_DD]
  UMSYpriorpar <- c(1 - exp(-Data@Mort[x] * 0.5), 0.3) # Prior for UMSY is that corresponding to F = 0.5 M with CV = 0.3
  UMSYprior <- c(alphaconv(UMSYpriorpar[1], prod(UMSYpriorpar)), betaconv(UMSYpriorpar[1], prod(UMSYpriorpar))) # Convert to beta parameters
  params <- log(c(UMSYpriorpar[1]/(1 - UMSYpriorpar[1]), 3*mean(C_hist, na.rm = T), Data@Mort[x]))
  opt <- optim(params, DD_R, opty = 1, So_DD = So_DD, Alpha_DD = Alpha_DD, 
               Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, wa_DD = wa_DD, E_hist = E_hist, 
               C_hist = C_hist, UMSYprior = UMSYprior, method = "BFGS", hessian = TRUE)
  
  if (reps > 1) {
    samps <- mvtnorm::rmvnorm(reps,opt$par,solve(opt$hessian)) # assuming log  
  } else {
    samps <- matrix(c(opt$par[1], opt$par[2], opt$par[3]), nrow = 1)
  }
  
  # parameters are multivariate normal hessian approximation
  # samps <- cbind(rnorm(reps, opt$par[1], ((opt$par[1])^2)^0.5 * 0.1), 
  #                rnorm(reps, opt$par[2], ((opt$par[2])^2)^0.5 * 0.1), 
  #                rnorm(reps, opt$par[3], ((opt$par[3])^2)^0.5 * 0.1))
 
  getVals <- sapply(1:reps, function(i) DD_R(params=samps[i, ], opty = 2, So_DD = So_DD, 
                                             Alpha_DD = Alpha_DD, Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, 
                                             wa_DD = wa_DD, E_hist = E_hist, C_hist = C_hist, UMSYprior = UMSYprior))
  
  TAC <- unlist(getVals[1,])
  dep <- unlist(getVals[2,])
  Cpredict <- do.call("cbind", getVals[3,])
  B_DD <- do.call("cbind", getVals[4,])

  if (!is.null(hcr)) {
    cond1 <- !is.na(dep) & dep < hcr[1] & dep > hcr[2]
    cond2 <- !is.na(dep) & dep < hcr[2]
    TAC[cond1] <- TAC[cond1] * (dep[cond1] - hcr[2])/(hcr[1]-hcr[2])
    TAC[cond2] <- TAC[cond2] * tiny  # this has to still be stochastic albeit very small
  }
  
  
  return(list(TAC=TAC, C_hist=C_hist, I_hist=I_hist, E_hist=E_hist, dep=dep, 
              Cpredict=Cpredict, B_DD=B_DD, hcr=hcr, Year=Data@Year[yind]))
}


#' Delay-Difference Internal Function
#'
#' @param params internal parameter
#' @param opty  internal parameter
#' @param So_DD  internal parameter
#' @param Alpha_DD  internal parameter
#' @param Rho_DD  internal parameter
#' @param ny_DD  internal parameter 
#' @param k_DD internal parameter
#' @param wa_DD  internal parameter
#' @param E_hist  internal parameter
#' @param C_hist internal parameter 
#' @param UMSYprior internal parameter 
#'
#' @export
#'
#' @keywords internal
DD_R <- function(params, opty, So_DD, Alpha_DD, Rho_DD, ny_DD, k_DD, wa_DD, E_hist,
                 C_hist, UMSYprior) {
  UMSY_DD = 1/(1 + exp(-params[1])) # Logit transform to constrain u between 0-1
  MSY_DD = exp(params[2])
  q_DD = exp(params[3])
  SS_DD = So_DD * (1 - UMSY_DD)  # Initialise for UMSY, MSY and q leading.
  Spr_DD = (SS_DD * Alpha_DD/(1 - SS_DD) + wa_DD)/(1 - Rho_DD * SS_DD)
  DsprDu_DD = ((Alpha_DD + Spr_DD * (1 + Rho_DD - 2 * Rho_DD * SS_DD))/((1 - Rho_DD * SS_DD) * (1 - SS_DD)) + 
                 Alpha_DD * SS_DD/((1 - Rho_DD * SS_DD) * (1 - SS_DD)^2) - Spr_DD/(1 - SS_DD)) * -So_DD
  Arec_DD = 1/(((1 - UMSY_DD)^2) * (Spr_DD + UMSY_DD * DsprDu_DD))
  Brec_DD = UMSY_DD * (Arec_DD * Spr_DD - 1/(1 - UMSY_DD))/MSY_DD
  Spr0_DD = (So_DD * Alpha_DD/(1 - So_DD) + wa_DD)/(1 - Rho_DD * So_DD)
  Ro_DD = (Arec_DD * Spr0_DD - 1)/(Brec_DD * Spr0_DD)
  Bo_DD = Ro_DD * Spr0_DD
  No_DD = Ro_DD/(1 - So_DD)
  
  B_DD <- rep(NA, ny_DD + 1)
  N_DD <- rep(NA, ny_DD + 1)
  R_DD <- rep(NA, ny_DD + k_DD)
  Cpred_DD <- rep(NA, ny_DD)
  
  B_DD[1] = Bo_DD
  N_DD[1] = No_DD
  R_DD[1:k_DD] = Ro_DD
  
  for (tt in 1:ny_DD) {
    Surv_DD = So_DD * exp(-q_DD * E_hist[tt])
    Cpred_DD[tt] = B_DD[tt] * (1 - exp(-q_DD * E_hist[tt]))
    Sp_DD = B_DD[tt] - Cpred_DD[tt]
    R_DD[tt + k_DD] = Arec_DD * Sp_DD/(1 + Brec_DD * Sp_DD)
    B_DD[tt + 1] = Surv_DD * (Alpha_DD * N_DD[tt] + Rho_DD * B_DD[tt]) +
      wa_DD * R_DD[tt + 1]
    N_DD[tt + 1] = Surv_DD * N_DD[tt] + R_DD[tt + 1]
  }
  tiny <- 1e-15
  Cpred_DD[Cpred_DD < tiny] <- tiny
  
  if (opty == 1) {
    # The following conditions must be met for positive values
    # of Arec_DD and Brec_DD, respectively:
    # umsy * DsprDu + Spr_DD > 0 and Arec_DD * Spr_DD * (1 - UMSY_DD) - 1 > 0
    # Thus, create a likelihood penalty of 100 if either condition is not met
    umsy_penalty <- ifelse(Spr_DD + UMSY_DD * DsprDu_DD > 0, 0, UMSY_DD * 100)
    alpha_penalty <- ifelse(Arec_DD * Spr_DD * (1 - UMSY_DD) - 1 > 0, 0, UMSY_DD * 100)
    
    sigma <- sqrt(sum((log(C_hist) - log(Cpred_DD))^2)/ny_DD) # Analytical solution
    
    test <- dnorm(log(C_hist), log(Cpred_DD), sigma, log = T)
 
    test2 <- dbeta(UMSY_DD, UMSYprior[1], UMSYprior[2], log = T)
    test[is.na(test)] <- -1000
    test[test == (-Inf)] <- -1000
    if (is.na(test2) | test2 == -Inf | test2 == Inf)  test2 <- 1000
    return(-sum(test, test2) + umsy_penalty + alpha_penalty)  # return objective function
  } else {
    return(list(TAC= UMSY_DD * B_DD[ny_DD + 1], dep=B_DD[tt + 1]/Bo_DD, 
                Cpred_DD=Cpred_DD,B_DD=B_DD))
  } 
  # } if (opty == 2) {
  #   # return MLE TAC estimate
  #  
  # } else if (opty == 3) {
  #   
  # } else {
  #   cbind(C_hist, Cpred_DD)  # return observations vs predictions
  # }
}

#' Delay - Difference Stock Assessment with UMSY and MSY as leading parameters
#' 
#' A simple delay-difference assessment that estimates the TAC using a
#' time-series of catches and a relative abundance index. Conditioned on effort.
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
#' @templateVar mp DD
#' @template MPtemplate
#' @template MPuses
#' 
#' @author T. Carruthers
#' @references  
#' Carruthers, T, Walters, C.J,, and McAllister, M.K. 2012. Evaluating methods that classify
#' fisheries stock status using only fisheries catch data. Fisheries Research 119-120:66-79.
#' 
#' Hilborn, R., and Walters, C. 1992. Quantitative Fisheries Stock Assessment: Choice,
#' Dynamics and Uncertainty. Chapman and Hall, New York. 
#' @describeIn DD Base version where the TAC = UMSY * Current Biomass.
#' 
#' @examples 
#' DD(1, Data=DLMtool::Atlantic_mackerel, plot=TRUE)
#' @family Delay-Difference MPs
#' 
#' @export 
DD <- function(x, Data, reps = 100, plot=FALSE) {
  runDD <- DD_(x, Data, reps)
  TAC <- TACfilter(runDD$TAC)
  
  if (plot) DD_plot(x, runDD, Data, TAC)
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}
class(DD) <- "MP"



#' @templateVar mp DD4010
#' @template MPuses
#' @examples 
#' DD4010(1, Data=DLMtool::Atlantic_mackerel, plot=TRUE)
#' @describeIn DD A 40-10 rule is imposed over the TAC recommendation.
#' @export 
DD4010 <- function(x, Data, reps = 100, plot=FALSE) {
  
  runDD <- DD_(x, Data, reps, hcr=c(.4, .1))
  TAC <- TACfilter(runDD$TAC)
  
  if (plot) DD_plot(x, runDD, Data, TAC)
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}
class(DD4010) <- "MP"







## Dynamic Fratio ####






#' Dynamic Fratio MP
#' 
#' The Fratio MP with a controller that changes the level of F according to the
#' estimated relationship between surplus production and biomass. Ie lower F when dSP/dB
#' is positive and higher F when dSP/dB is negative.
#' 
#' The method smoothes historical catches and biomass and then infers the
#' relationship between surplus production and biomass (as suggested by Mark
#' Maunder and Carl Walters). The approach then regulates a F based policy
#' according to this gradient in which F may range between two different
#' fractions of natural mortality rate.
#' 
#' The core advantage is the TAC(t) is not strongly determined by TAC(t-1) and
#' therefore errors are not as readily propagated. The result is method that
#' tends to perform alarmingly well and therefore requires debunking ASAP.
#'
#' The catch limit (TAC) is calculated as:
#' \deqn{\textrm{TAC}=F * B}
#' where \eqn{F} is fishing mortality and \eqn{B} is the estimated current biomass.
#' 
#' \eqn{F} is calculated as:
#' \deqn{F = F_{\textrm{MSY}} \exp{-gG}}
#' where \eqn{F_{\textrm{MSY}}} is calculated from assumed values of \eqn{\frac{F_{\textrm{MSY}}}{M}} and
#' \eqn{M}, *g* is a gain parameter and *G* is the estimated gradient in surplus
#' production (*SP*) as a function of biomass (*B*). Surplus production for year *y* is calculated as:
#' \deqn{SP_y = B_{y+1} - B_y + C_y}
#' Trends in historical catch (*C*) and biomass (*B*) are both estimated using a loess smoother, over the last `yrsmth` years,
#' of available catch and a time-series of abundance, calculated from an index of abundance (`Data@Ind`) 
#' and an estimate of abundance (`Data@Abun`) for the current year.   
#' 
#' 
#' @templateVar mp DynF
#' @template MPtemplate
#' @template MPuses
#' 
#' @param gg A gain parameter that modifies F according to the gradient in
#' surplus production with biomass
#' @param yrsmth The number of historical recent years used for smoothing catch
#' and biomass data
#' 
#' @author T. Carruthers
#' @references Made-up for this package.
#' @family Fmsy/M methods
#' @examples 
#' DynF(1, Data=DLMtool::Atlantic_mackerel, plot=TRUE)
#' @export 
DynF <- function(x, Data, reps = 100, plot=FALSE, yrsmth = 10, gg = 2) {
  
  dependencies = "Data@Year, Data@Cat, Data@Ind, Data@Abun, Data@Mort, Data@FMSY_M"
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
  years <- Data@Year[ind]
  C_dat <- log(Data@Cat[x, ind])
  C_dat[C_dat == -Inf] <- 0
  B_dat <- log(Data@Ind[x, ind]/Data@Ind[x, ind[yrsmth]] * Data@Abun[x])
  B_dat[B_dat == -Inf] <- 0
  C_hist <- exp(predict(loess(C_dat ~ ind, degree = 1)))
  B_hist <- exp(predict(loess(B_dat ~ ind, degree = 1)))
  
  
  ind <- 2:yrsmth
  ind1 <- 1:(yrsmth - 1)
  SP_hist <- B_hist[ind] - B_hist[ind1] + C_hist[ind1]
  
  Frat <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x]) * 
    trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])
  Flim <- matrix(NA, nrow = 2, ncol = reps)
  Flim[1, ] <- Frat * 0.5
  Flim[2, ] <- Frat * 2
  
  # yind <- 1:length(SP_hist)
  # SP_mu <- predict(lm(SP_hist ~ yind), newdat = list(yind = length(SP_hist) + 1))
  # SP_se <- predict(lm(SP_hist ~ yind), newdat = list(yind = length(SP_hist) + 1), se = T)$se.fit
  # SP_new <- rnorm(reps, SP_mu, SP_se/2)
  Glm <- summary(lm(SP_hist ~ B_hist[ind1]))$coefficients[2, 1:2]  # gradient in SP as function of biomass
  G_new <- rnorm(reps, Glm[1], Glm[2]/2) # sample gradients
  # G_new[G_new>2*Frat]<-2*Frat[G_new<(2*Frat)]
  # G_new[G_new<(-2*Frat)]<--2*Frat[G_new<(-2*Frat)]
  G_new[G_new > 0] <- G_new[G_new > 0] * 3
  newF <- Frat * exp(-G_new * gg)
  newF[newF < Flim[1]] <- Flim[1]
  newF[newF > Flim[2]] <- Flim[2]
  
  TAC <- TACfilter(newF * B_hist[yrsmth])
  
  if (plot) DynF_plot(C_dat, C_hist, TAC, yrsmth, B_dat, B_hist, Data, SP_hist, 
                      ind, ind1, G_new, Frat,newF, years)
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
  
}
class(DynF) <- "MP"



#' An adaptive MP that uses trajectory in inferred suplus production and
#' fishing mortality rate to update a TAC
#' 
#' Fishing rate is modified each year according to the gradient of surplus
#' production with biomass (aims for zero).  F is bounded by FMSY/2 and 2FMSY
#' and walks in the logit space according to dSP/dB. This is derived from the
#' theory of Maunder 2014.
#' 
#' Tested in Carruthers et al. 2015.
#' 
#' @templateVar mp Fadapt
#' @template MPtemplate
#' @template MPuses
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param gg A gain parameter controlling the speed in update in TAC.
#' @return A numeric vector of quota recommendations
#' @author T. Carruthers
#' @references Carruthers et al. 2015. Performance evaluation of simple
#' management procedures. ICES J. Mar Sci. 73, 464-482.
#' 
#' Maunder, M. 2014.
#' http://www.iattc.org/Meetings/Meetings2014/MAYSAC/PDFs/SAC-05-10b-Management-Strategy-Evaluation.pdf
#' @family Fmsy/M methods
#' @examples 
#' Fadapt(1, Data=DLMtool::Atlantic_mackerel, plot=TRUE)
#' @export 
Fadapt <- function(x, Data, reps = 100, plot=FALSE, yrsmth = 7, gg = 1) {
  
  dependencies = "Data@Year, Data@Cat, Data@Ind, Data@Abun, Data@Mort, Data@FMSY_M"
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
  years <- Data@Year[ind]
  C_dat <- log(Data@Cat[x, ind])
  C_dat[C_dat == -Inf] <- 0
  B_dat <- log(Data@Ind[x, ind]/Data@Ind[x, ind[yrsmth]] * Data@Abun[x])
  B_dat[B_dat == -Inf] <- 0
  C_hist <- exp(predict(loess(C_dat ~ ind, degree = 1)))
  B_hist <- exp(predict(loess(B_dat ~ ind, degree = 1)))
  
  ind <- 2:yrsmth
  ind1 <- 1:(yrsmth - 1)
  SP_hist <- B_hist[ind] - B_hist[ind1] + C_hist[ind1]
  
  Frat <- Data@Mort[x] * Data@FMSY_M[x]
  Flim <- Frat * c(0.5, 2)
  Flimr <- Flim[2] - Flim[1]
  
  yind <- 1:length(SP_hist)
  SP_mu <- predict(lm(SP_hist ~ yind), newdat = list(yind = length(SP_hist) + 1))
  SP_se <- predict(lm(SP_hist ~ yind), newdat = list(yind = length(SP_hist) + 1), se = T)$se.fit
  SP_new <- rnorm(reps, SP_mu, SP_se/2)
  Glm <- summary(lm(SP_hist ~ B_hist[ind1]))$coefficients[2, 1:2]  # plot(B_hist[ind1],SP_hist) # points(B_hist[ind1],SP_hist,col='green')
  G_new <- rnorm(reps, Glm[1], Glm[2])
  
  Fold <- mean(C_hist/B_hist)
  
  if (Fold < Flim[1]) Fmod1 <- (-2)
  if (Fold > Flim[2]) Fmod1 <- 2
  if (Fold > Flim[1] & Fold < Flim[2]) {
    Ffrac <- (Fold - Flim[1])/Flimr
    Fmod1 <- log(Ffrac/(1 - Ffrac))
  }
  Fmod2 <- Fmod1 + gg * -G_new
  newF <- Flim[1] + (exp(Fmod2)/(1 + exp(Fmod2))) * Flimr
  TAC <- TACfilter(newF * B_hist[yrsmth])
  
  if (plot) Fadapt_plot(C_dat, C_hist, TAC, yrsmth, B_dat, B_hist, Data, SP_hist, 
                      ind, ind1, G_new, Frat,newF, years)
  
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}
class(Fadapt) <- "MP"



#' Internal function to estimate r
#'
#' @param x iteration number 
#' @param Data Object of class `Data` 
#' @param Mvec M
#' @param Kvec K
#' @param Linfvec Linf 
#' @param t0vec t0
#' @param hvec steepness
#' @param maxage maximum age
#' @param r_reps number of reps
#'
#' @export
#'
#' @keywords internal 
getr <- function(x, Data, Mvec, Kvec, Linfvec, t0vec, hvec, maxage,
                 r_reps = 100) {
  r <- rep(NA, r_reps)
  for (i in 1:r_reps) {
    log.r <- log(0.3)
    opt <- optimize(demofn, lower = log(1e-04), upper = log(1.4), M = Mvec[i],
                   amat = iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],
                              Data@L50[x]), sigma = 0.2, K = Kvec[i], Linf = Linfvec[i],
                   to = t0vec[i], hR = hvec[i], maxage = maxage, a = Data@wla[x],
                   b = Data@wlb[x])
    # demographic2(opt$minimum,M[x],ageM[x],0.2,K[x],Linf,t0,steepness[x],maxage,a,b)$r
    r[i] <- exp(opt$minimum)
  }
  r
}



#' Internal demographic function
#'
#' @param log.r internal
#' @param M internal
#' @param amat internal
#' @param sigma internal
#' @param K internal
#' @param Linf internal
#' @param to internal
#' @param hR internal
#' @param maxage internal
#' @param a internal
#' @param b internal
#'
#' @export
#'
#' @keywords internal
demographic2 <- function(log.r, M, amat, sigma, K, Linf, to, hR, maxage, a, b) {
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

#' @describeIn demographic2 Internal function
#' @export
demofn <- function(log.r, M, amat, sigma, K, Linf, to, hR, maxage, a, b) {
  demographic2(log.r, M, amat, sigma, K, Linf, to, hR, maxage = maxage, a, b)$epsilon
}


#' Internal function for Fdem MP
#'
#' @param x internal
#' @param Data internal
#' @param reps internal
#' @param Ac Optional abundance value 
#'
#' @export
#' @keywords internal
Fdem_ <- function(x, Data, reps, Ac=NULL) {
  # Demographic FMSY estimate (FMSY=r/2)
  Mvec <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  Linfc = trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
    if (Data@vbt0[x]<0) {
      t0c <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
    } else {
      t0c <- trlnorm(reps, Data@vbt0[x], Data@CV_vbt0[x])
    }
  } else {
    t0c <- rep(Data@vbt0[x], reps)
  }
  t0c[!is.finite(t0c)] <- 0
  # hvec <- trlnorm(reps, Data@steep[x], Data@CV_steep[x])
  hvec <- sample_steepness2(reps, Data@steep[x], Data@CV_steep[x])
  
  if (is.null(Ac)) Ac <- trlnorm(reps, Data@Abun[x], Data@CV_Abun[x])
  FMSY <- getr(x, Data, Mvec, Kc, Linfc, t0c, hvec, maxage = Data@MaxAge, 
               r_reps = reps)/2
  TAC <- FMSY * Ac
  return(list(TAC=TAC, Ac=Ac, FMSY=FMSY))
}


## ADD description and equations and plot ####

#' Demographic FMSY method
#' 
#' FMSY is calculated as r/2 where r is calculated from a demographic approach
#' (inc steepness). Coupled with an estimate of current abundance that gives
#' you the OFL.
#' 
#' 
#' @templateVar mp Fdem
#' @template MPtemplate
#' @template MPuses
#' 
#' @author T. Carruthers
#' @references McAllister, M.K., Pikitch, E.K., and Babcock, E.A. 2001. Using
#' demographic methods to construct Bayesian priors for the intrinsic rate of
#' increase in the Schaefer model and implications for stock rebuilding. Can.
#' J. Fish. Aquat. Sci. 58: 1871-1890.
#' @describeIn Fdem This uses Murdoch McAllister's demographic r
#' method to derive FMSY (r/2) and then makes the quota r*current biomass / 2.
#' Easy.
#' @examples 
#' Fdem(1, DLMtool::Atlantic_mackerel, plot=TRUE)
#' @export 
Fdem <- function(x, Data, reps = 100, plot=FALSE) {
  runFdem <- Fdem_(x, Data, reps)
  TAC <- TACfilter(runFdem$TAC)
  
  if (plot) Fdem_plot()
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}
class(Fdem) <- "MP"

## TO DO - add plot ####
Fdem_plot <- function() {
  
}

#' @templateVar mp Fdem_CC
#' @template MPuses
#' 
#' @param Fmin The minimum fishing mortality rate derived from the catch-curve
#' analysis
#' 
#' @describeIn Fdem FMSY is calculated as r/2 from a demographic r prior method, current
#' abundance is estimated from catch curve analysis.
#' @examples 
#' Fdem_CC(1, DLMtool::SimulatedData, plot=TRUE)
#' @export 
Fdem_CC <- function(x, Data, reps = 100, plot=FALSE, Fmin = 0.005) {
  
  Mvec <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])
  MuC <- Data@Cat[x, length(Data@Cat[x, ])]
  Cc <- trlnorm(reps, MuC, Data@CV_Cat[x])
  Zdb <- CC(x, Data, reps = reps * 10)
  Fdb <- Zdb - Mvec
  ind <- (1:(reps * 10))[Fdb > Fmin][1:reps]
  
  Fdb <- Fdb[ind]
  SM <- sum(is.na(ind))
  if (SM > 0) {
    Fdb[is.na(ind)] <- Fmin
  }
  
  Ac <- Cc/(1 - exp(-Fdb))
  
  runFdem <- Fdem_(x, Data, reps, Ac=Ac)
  TAC <- TACfilter(runFdem$TAC)
  
  if (plot) Fdem_plot()
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}
class(Fdem_CC) <- "MP"

#' @describeIn Fdem Demographic F (r/2) method using the mean length estimator to calculate
#' current abundance.
#' @examples 
#' Fdem_ML(1, DLMtool::SimulatedData, plot=TRUE)
#' @export 
Fdem_ML <- function(x, Data, reps = 100, plot=FALSE, Fmin = 0.005) {
  Mvec <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])
  Kc <- trlnorm(reps * 10, Data@vbK[x], Data@CV_vbK[x])
  Linfc = trlnorm(reps * 10, Data@vbLinf[x], Data@CV_vbLinf[x])
  MuC <- Data@Cat[x, length(Data@Cat[x, ])]
  Cc <- trlnorm(reps * 10, MuC, Data@CV_Cat[x])
  Z <- MLne(x, Data, Linfc = Linfc, Kc = Kc, ML_reps = reps * 10, MLtype = "F")
  
  if (all(is.na(Z))) {
    Rec <- new("Rec")
    Rec@TAC <- TACfilter(rep(NA, reps))
    return(Rec)
  } 
  
  ind <- !is.na(Z)
  FM <- Z[ind] - Mvec[ind]
  ind <- (1:(reps * 10))[FM > Fmin][1:reps]
  FM <- FM[ind]
  SM <- sum(is.na(ind))
  if (SM > 0) {
    FM[is.na(ind)] <- Fmin
  }
  
  Ac <- Cc[ind]/(1 - exp(-FM))
  
  runFdem <- Fdem_(x, Data, reps, Ac=Ac)
  TAC <- TACfilter(runFdem$TAC)
  
  if (plot) Fdem_plot()
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}
class(Fdem_ML) <- "MP"



## ADD description and equations ####

#' An FMSY/M ratio method
#' 
#' Calculates the OFL based on a fixed ratio of FMSY to M multiplied by a
#' current estimate of abundance.
#'  
#' @templateVar mp Fratio
#' @template MPtemplate
#' @template MPuses
#' 
#' @author T. Carruthers
#' @references Gulland, J.A., 1971. The fish resources of the ocean. Fishing
#' News Books, West Byfleet, UK.
#' 
#' Martell, S., Froese, R., 2012. A simple method for estimating MSY from catch
#' and resilience. Fish Fish. doi: 10.1111/j.1467-2979.2012.00485.x.
#' 
#' @describeIn Fratio A simple method that tends to outperform many other approaches alarmingly
#' often even when current biomass is relatively poorly known. The low stock
#' crash potential is largely due to the quite large difference between Fmax
#' and FMSY for most stocks.
#' @family Fmsy/M methods
#' @examples 
#' Fratio(1, DLMtool::Atlantic_mackerel, plot=TRUE)
#' @export 
Fratio <- function(x, Data, reps = 100, plot=FALSE) {
  # FMSY / M ratio method e.g. Gulland
  
  runFrat <- Fratio_(x, Data, reps)
  TAC <- TACfilter(runFrat$TAC)

  if (plot) Fratio_plot(x, Data, TAC, runFrat)
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}  
class(Fratio) <- "MP"

#' @templateVar mp Fratio4010
#' @template MPuses
#' @describeIn Fratio Paired with the 40-10 rule that throttles back the OFL to zero at 10 percent of unfished biomass. Requires an estimate of current depletion.
#' @examples 
#' Fratio4010(1, DLMtool::Atlantic_mackerel, plot=TRUE)
#' @export 
Fratio4010 <- function(x, Data, reps = 100, plot=FALSE) {
  runFrat <- Fratio_(x, Data, reps)
  TAC <- TACfilter(runFrat$TAC)

  # 40-10 rule
  # if (is.na(Data@Dt[x]) || is.na(Data@CV_Dt[x])) return(new("Rec"))
  # Bt_K <- trlnorm(reps, Data@Dt[x], Data@CV_Dt[x])
  if (is.na(Data@Dep[x]) | is.na(Data@CV_Dep[x])) return(new("Rec"))
  depo <- max(0.01, min(0.99, Data@Dep[x]))  
  Bt_K <- rbeta(reps * 100, alphaconv(depo, min(depo * Data@CV_Dep[x], 
                                                (1 - depo) * Data@CV_Dep[x])), 
                betaconv(depo, min(depo * Data@CV_Dep[x], (1 - depo) * Data@CV_Dep[x]))) 
  Bt_K <- Bt_K[Bt_K >= 0.01 & Bt_K <= 0.99][1:reps]  # interval censor (0.01,0.99)  as in Dick and MacCall 2011
  
  cond1 <- Bt_K < 0.4 & Bt_K > 0.1
  cond2 <- Bt_K < 0.1
  TAC[cond1] <- TAC[cond1] * (Bt_K[cond1] - 0.1)/0.3
  TAC[cond2] <- TAC[cond2] * tiny  # this has to still be stochastic albeit very small
  
  runFrat$Bt_K <- Bt_K 
  if (plot) Fratio_plot(x, Data, TAC, runFrat)
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}  
class(Fratio4010) <- "MP"

#' @templateVar mp DepF
#' @template MPuses
#' 
#' @describeIn Fratio Depletion Corrected Fratio: the Fratio MP with a harvest control 
#' rule that reduces F according to the production curve given an estimate of current 
#' stock depletion (made-up for this package).
#' @export 
DepF <- function(x, Data, reps = 100, plot=FALSE) {
  
  runFrat <- Fratio_(x, Data, reps)
  TAC <- TACfilter(runFrat$TAC)
  
  if (is.na(Data@Dep[x]) | is.na(Data@CV_Dep[x])) return(new("Rec"))
  depo <- max(0.01, min(0.99, Data@Dep[x]))  
  Bt_K <- rbeta(reps * 100, alphaconv(depo, min(depo * Data@CV_Dep[x], 
                                                (1 - depo) * Data@CV_Dep[x])), 
                betaconv(depo, min(depo * Data@CV_Dep[x], (1 - depo) * Data@CV_Dep[x]))) 
  Bt_K <- Bt_K[Bt_K >= 0.01 & Bt_K <= 0.99][1:reps]  # interval censor (0.01,0.99)  as in Dick and MacCall 2011
  

  adj <- Bt_K * (1 - Bt_K) * 4
  adj[Bt_K > 0.5] <- 1
  
  TAC <- TAC * adj
  runFrat$Bt_K <- Bt_K
  
  if (plot) Fratio_plot(x, Data, TAC, runFrat)
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}
class(DepF) <- "MP"


#' @templateVar mp Fratio_CC
#' @template MPuses
#' @param Fmin Minimum current fishing mortality rate for the catch-curve analysis
#' @describeIn Fratio Current abundance is estimated using average catch and estimate of F from an age-based catch curve 
#' @examples 
#' Fratio_CC(1, DLMtool::SimulatedData, plot=TRUE)
#' @export
Fratio_CC <- function(x, Data, reps = 100, plot=TRUE, Fmin = 0.005) {
  # estimate abundance from average catch and F
  MuC <- Data@Cat[x, length(Data@Cat[x, ])]
  Cc <- trlnorm(reps, MuC, Data@CV_Cat[x])
  Mdb <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
  Zdb <- CC(x, Data, reps = reps * 10)
  Fdb <- Zdb - Mdb
  ind <- (1:(reps * 10))[Fdb > 0.005][1:reps]
  
  Fdb <- Fdb[ind]
  Mdb <- Mdb[ind]
  SM <- sum(is.na(ind))
  if (SM > 0) {
    Mdb[is.na(ind)] <- trlnorm(SM, Data@Mort[x], Data@CV_Mort[x])
    Fdb[is.na(ind)] <- Fmin
  }
  
  Ac <- Cc/(1 - exp(-Fdb))
  
  runFrat <- Fratio_(x, Data, reps, Abun=Ac)

  TAC <- TACfilter(runFrat$TAC)
  
  if (plot) Fratio_plot(x, Data, TAC, runFrat)
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}  
class(Fratio_CC) <- "MP"

#' @templateVar mp Fratio_ML
#' @template MPuses
#' @describeIn Fratio Current abundance is estimated using average catch and estimate of F from mean lengths 
#' @examples 
#' Fratio_ML(1, DLMtool::SimulatedData, plot=TRUE)
#' 
#' @export 
Fratio_ML <- function(x, Data, reps = 100, plot=FALSE) {
  # estimate abundance from average catch and F
  MuC <- Data@Cat[x, length(Data@Cat[x, ])]
  Cc <- trlnorm(reps * 10, MuC, Data@CV_Cat[x])
  Mdb <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
  Linfc <- trlnorm(reps * 10, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps * 10, Data@vbK[x], Data@CV_vbK[x])
  Z <- MLne(x, Data, Linfc = Linfc, Kc = Kc, ML_reps = reps * 10, MLtype = "F")
  if (all(is.na(Z))) return(new("Rec"))
  
  FM <- Z - Mdb
  ind <- which(FM>0)[1:reps]
  FM <- FM[ind]
  Ac <- Cc[ind]/(1 - exp(-FM))
  
  runFrat <- Fratio_(x, Data, reps, Abun=Ac)
  TAC <- TACfilter(runFrat$TAC)
  
  if (plot) Fratio_plot(x, Data, TAC, runFrat)
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}
class(Fratio_ML) <- "MP"




### TO DO - add Fratio plot ####
Fratio_plot <- function(x, Data, TAc, runFrat) {
  boxplot(runFrat)
  
}


#' Fratio internal function 
#'
#' @param x Iteration number
#' @param Data Object of class Data
#' @param reps Number of reps
#' @param Abun Optional estimate of abundance
#' 
#' @export
#'
#' @keywords internal 
Fratio_ <- function(x, Data, reps=100, Abun=NULL) {
  Frat <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x]) * trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x]) # estimate of Fmsy
  if (is.null(Abun)) Abun <- trlnorm(reps, Data@Abun[x], Data@CV_Abun[x])
  
  TAC <- Frat * Abun 
  
  return(list(TAC=TAC, Abun=Abun, Frat=Frat))
  
}













#' Geromont and Butterworth index slope Harvest Control Rule
#' 
#' An MP similar to SBT1 that modifies a time-series of catch recommendations
#' and aims for a stable catch rates.
#' 
#' Note that this is my interpretation of their approach and is now stochastic.
#' Currently it is generalized and is not 'tuned' to more detailed assessment
#' data which might explain why in some cases it leads to stock declines.
#' 
#' @usage GB_slope(x, Data, reps = 100, yrsmth = 5, lambda = 1)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrsmth Number of years for evaluating slope in relative abundance
#' index
#' @param lambda A gain parameter
#' @author T. Carruthers
#' @references Geromont, H.F. and Butterworth, D.S. 2014. Complex assessment or
#' simple management procedures for efficient fisheries management: a
#' comparative study. ICES J. Mar. Sci. doi:10.1093/icesjms/fsu017
#' @export GB_slope
GB_slope <- function(x, Data, reps = 100, yrsmth = 5, lambda = 1) {
  dependencies = "Data@Year, Data@Cat, Data@CV_Cat, Data@Ind"
  Catrec <- Data@Cat[x, length(Data@Cat[x, ])]
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
  I_hist <- Data@Ind[x, ind]
  yind <- 1:yrsmth
  slppar <- summary(lm(I_hist ~ yind))$coefficients[2, 1:2]
  Islp <- rnorm(reps, slppar[1], slppar[2])
  MuC <- Data@Cat[x, length(Data@Cat[x, ])]
  Cc <- stats::rlnorm(reps, mconv(MuC, Data@CV_Cat[x] * MuC), sdconv(MuC, 
                                                                     Data@CV_Cat[x] * MuC))
  TAC <- Cc * (1 + lambda * Islp)
  TAC[TAC > (1.2 * Catrec)] <- 1.2 * Catrec
  TAC[TAC < (0.8 * Catrec)] <- 0.8 * Catrec
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(GB_slope) <- "MP"



#' Geromont and Butterworth target CPUE and catch MP
#' 
#' An MP similar to SBT2 that modifies a time-series of catch recommendations
#' and aims for target catch rate and catch level based on BMSY/B0 and MSY,
#' respectively.
#' 
#' Note that this is my interpretation of their MP and is now stochastic.
#' Currently it is generalized and is not 'tuned' to more detailed assessment
#' data which might explain why in some cases it leads to stock declines.
#' 
#' @usage GB_target(x, Data, reps = 100, w = 0.5)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of quota samples
#' @param w A gain parameter
#' @author T. Carruthers
#' @references Geromont, H.F. and Butterworth, D.S. 2014. Complex assessment or
#' simple management procedures for efficient fisheries management: a
#' comparative study. ICES J. Mar. Sci. doi:10.1093/icesjms/fsu017
#' @export GB_target
GB_target <- function(x, Data, reps = 100, w = 0.5) {
  dependencies = "Data@Cat, Data@Cref, Data@Iref, Data@Ind"
  Catrec <- Data@Cat[x, length(Data@Cat[x, ])]
  TACtarg <- trlnorm(reps, Data@Cref[x], Data@CV_Cref)
  Itarg <- trlnorm(reps, Data@Iref[x], Data@CV_Iref)
  Iav <- mean(Data@Ind[x, (length(Data@Ind[x, ]) - 4):length(Data@Ind[x, 
                                                                      ])], na.rm = T)
  Irec <- mean(Data@Ind[x, (length(Data@Ind[x, ]) - 3):length(Data@Ind[x, 
                                                                       ])], na.rm = T)
  I0 <- 0.2 * Iav
  TAC <- rep(NA, reps)
  if (Irec > I0) 
    TAC <- TACtarg * (w + (1 - w) * ((Irec - I0)/(Itarg - I0)))
  if (Irec < I0) 
    TAC <- TACtarg * w * (Irec/I0)^2
  TAC[TAC > (1.2 * Catrec)] <- 1.2 * Catrec
  TAC[TAC < (0.8 * Catrec)] <- 0.8 * Catrec
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(GB_target) <- "MP"

#' G-control MP
#' 
#' A harvest control rule proposed by Carl Walters that uses trajectory in
#' inferred surplus production to make upward/downward adjustments to TAC
#' recommendations
#' 
#' 
#' @usage Gcontrol(x, Data, reps = 100, yrsmth = 10, gg = 2, glim = c(0.5,
#' 2))
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of quota samples
#' @param yrsmth The number of years over which to smooth catch and biomass
#' data
#' @param gg A gain parameter
#' @param glim A constraint limiting the maximum level of change in quota
#' recommendations
#' @author C. Walters and T. Carruthers
#' @references 
#' Carruthers et al. 2015. Performance evaluation of simple
#' management procedures. ICES J. Mar Sci. 73, 464-482.
#' @export Gcontrol
Gcontrol <- function(x, Data, reps = 100, yrsmth = 10, gg = 2, glim = c(0.5, 
                                                                        2)) {
  dependencies = "Data@Year, Data@Cat, Data@Ind, Data@Abun"
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
  C_dat <- log(Data@Cat[x, ind])
  C_dat[C_dat == -Inf] <- 0
  B_dat <- log(Data@Ind[x, ind]/Data@Ind[x, ind[yrsmth]] * Data@Abun[x])
  B_dat[B_dat == -Inf] <- 0
  C_hist <- exp(predict(loess(C_dat ~ ind, degree = 1)))
  B_hist <- exp(predict(loess(B_dat ~ ind, degree = 1)))
  ind <- 2:yrsmth
  ind1 <- 1:(yrsmth - 1)
  SP_hist <- B_hist[ind] - B_hist[ind1] + C_hist[ind1]
  yind <- 1:length(SP_hist)
  SP_mu <- predict(lm(SP_hist ~ yind), newdat = list(yind = length(SP_hist) + 
                                                       1))
  SP_se <- predict(lm(SP_hist ~ yind), newdat = list(yind = length(SP_hist) + 
                                                       1), se = T)$se.fit
  SP_new <- rnorm(reps, SP_mu, SP_se/2)
  Glm <- summary(lm(SP_hist ~ B_hist[ind1]))$coefficients[2, 1:2]
  G_new <- rnorm(reps, Glm[1], Glm[2]/2)
  
  TAC <- SP_new * (1 - gg * G_new)
  TAC[TAC < glim[1] * C_hist[yrsmth]] <- glim[1] * C_hist[yrsmth]
  TAC[TAC > glim[2] * C_hist[yrsmth]] <- glim[2] * C_hist[yrsmth]
  
  # Carr<-cbind(array(rep(Data@Cat[x,],each=reps),c(reps,length(Data@Cat[x,]))),TAC)
  # Warr<-(Data@Mort[x]*exp(-Data@Mort[x]*(1:ncol(Carr))))[ncol(Carr):1]
  # Warr<-Warr/sum(Warr)
  # TAC<-apply(t(matrix(Warr,nrow=ncol(Carr),ncol=reps))*Carr,1,sum)
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(Gcontrol) <- "MP"







#' Index Confidence Interval (ICI) MP by Jardim et al. (2015)
#' 
#' The MP adjusts catch based on the value of the index in the current year relative to the 
#' time series mean and standard error.
#'  
#' @describeIn ICI The mean and standard error of the index time series is calculated. There are two thresholds 
#' which delineates whether catch is reduced, held constant, or increased. The catch is reduced by 0.75
#' if the Z-score of the current year's index is less than -0.44. The catch is increased by 1.05
#' if the Z-score of the current year's index is greater than 1.96. Otherwise, the catch is held constant.
#' 
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @author Coded by Q. Huynh. Developed by Jardim et al. (2015)
#' @references Ernesto Jardim, Manuela Azevedo, Nuno M. Brites, Harvest control rules for 
#' data limited stocks using length-based reference points and survey biomass indices, 
#' Fisheries Research, Volume 171, November 2015, Pages 12-19, ISSN 0165-7836, 
#' https://doi.org/10.1016/j.fishres.2014.11.013
#' @export 
ICI <- function(x, Data, reps) {
  dependencies = "Data@Ind, Data@CV_Ind, Data@Cat, Data@CV_Cat"
  
  Index <- Data@Ind[x, ]
  Index <- Index[!is.na(Index)]
  nI <- length(Index)
  Ind.samp <- trlnorm(reps * nI, Index, Data@CV_Ind[x])
  Ind.samp <- matrix(Ind.samp, ncol = reps)
  
  muI <- apply(Ind.samp, 2, mean, na.rm = TRUE)
  sigmaI <- apply(Ind.samp, 2, sd, na.rm = TRUE)
  
  Ind <- Ind.samp[nI, ]
  z.low <- -0.44 #qnorm(0.33)
  z.upp <- 1.96 #qnorm(0.974)
  
  ci.low <- muI + z.low * sigmaI / sqrt(nI)
  ci.high <- muI + z.upp * sigmaI / sqrt(nI)
  
  alpha <- rep(NA, reps)
  alpha[Ind < ci.low] <- 0.75
  alpha[Ind > ci.high] <- 1.05
  alpha[Ind >= ci.low & Ind <= ci.high] <- 1
  
  Cat <- Data@Cat[x, length(Data@Cat[x, ])]
  Cc <- trlnorm(reps, Cat, Data@CV_Cat[x])
  
  TAC <- alpha * Cc
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(ICI) <- "MP"



#' @describeIn ICI Compared to \code{ICI}, this method is less precautionary of the two ICI MPs by allowing for a 
#' larger increase in TAC and a lower threshold of the index to decrease the TAC. The catch is reduced by 0.75
#' if the Z-score of the current year's index is less than -1.96. The catch is increased by 1.25
#' if the Z-score of the current year's index is greater than 1.96. Otherwise, the catch is held constant.
#' @export 
ICI2 <- function(x, Data, reps) {
  dependencies = "Data@Ind, Data@CV_Ind, Data@Cat, Data@CV_Cat"
  
  Index <- Data@Ind[x, ]
  Index <- Index[!is.na(Index)]
  nI <- length(Index)
  Ind.samp <- trlnorm(reps * nI, Index, Data@CV_Ind[x])
  Ind.samp <- matrix(Ind.samp, ncol = reps)
  
  muI <- apply(Ind.samp, 2, mean, na.rm = TRUE)
  sigmaI <- apply(Ind.samp, 2, sd, na.rm = TRUE)
  
  Ind <- Ind.samp[nI, ]
  z.low <- -1.96 #qnorm(0.025)
  z.upp <- 1.96 #qnorm(0.975)
  
  ci.low <- muI + z.low * sigmaI / sqrt(nI)
  ci.high <- muI + z.upp * sigmaI / sqrt(nI)
  
  alpha <- rep(NA, reps)
  alpha[Ind < ci.low] <- 0.75
  alpha[Ind > ci.high] <- 1.25
  alpha[Ind >= ci.low & Ind <= ci.high] <- 1
  
  Cat <- Data@Cat[x, length(Data@Cat[x, ])]
  Cc <- trlnorm(reps, Cat, Data@CV_Cat[x])
  
  TAC <- alpha * Cc
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(ICI2) <- "MP"




#' Mean index ratio MP from Jardim et al. 2015
#' 
#' The TAC is adjusted by the ratio alpha, where the numerator 
#' being the mean index in the most recent two years of the time series and the denominator
#' being the mean index in the three years prior to those in the numerator. This MP is the 
#' stochastic version of Method 3.2 used by ICES for Data-Limited Stocks (ICES 2012).
#' 
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrs Vector of length 2 specifying the reference years
#' @author Coded by Q. Huynh. Developed by Jardim et al. (2015)
#' @references Ernesto Jardim, Manuela Azevedo, Nuno M. Brites, Harvest control rules for 
#' data limited stocks using length-based reference points and survey biomass indices, 
#' Fisheries Research, Volume 171, November 2015, Pages 12-19, ISSN 0165-7836, 
#' https://doi.org/10.1016/j.fishres.2014.11.013
#' 
#' ICES. 2012. ICES Implementation of Advice for Data-limited Stocks in 2012 in its 2012
#' Advice. ICES CM 2012/ACOM 68. 42 pp.
#' @export
Iratio <- function(x, Data, reps, yrs = c(2, 5)) {
  dependencies = "Data@Ind, Data@CV_Ind, Data@Cat, Data@CV_Cat"
  
  ind.num <- (length(Data@Year) - yrs[1]+1):length(Data@Year)
  ind.den <- (length(Data@Year) - yrs[2]+1):(length(Data@Year) - yrs[1])
  
  I.num <- trlnorm(reps * length(ind.num), Data@Ind[x, ind.num], Data@CV_Ind[x])
  I.num <- matrix(I.num, ncol = reps)
  
  I.den <- trlnorm(reps * length(ind.den), Data@Ind[x, ind.den], Data@CV_Ind[x])
  I.den <- matrix(I.den, ncol = reps)
  
  alpha <- apply(I.num, 2, mean, na.rm = TRUE)/apply(I.den, 2, mean, na.rm = TRUE)
  Cat <- Data@Cat[x, length(Data@Cat[x, ])]
  Cc <- trlnorm(reps, Cat, Data@CV_Cat[x])
  TAC <- alpha * Cc
  
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(Iratio) <- "MP"


#' A management procedure that incrementally adjusts the TAC to maintain a
#' constant CPUE or relative abundance index
#' 
#' The least biologically precautionary of two constant index / CPUE methods
#' proposed by Geromont and Butterworth 2014. Tested by Carruthers et al. 2015
#' 
#' Tested by Carruthers et al. 2015.
#' 
#' @usage Islope1(x, Data, reps = 100, yrsmth = 5, lambda=0.4,xx=0.2)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrsmth Years over which to smooth recent estimates of surplus production
#' @param lambda A gain parameter controlling the speed in update in TAC.
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @return A numeric vector of quota recommendations
#' @author T. Carruthers
#' @references Carruthers et al. 2015. Performance evaluation of simple
#' management procedures. ICES J. Mar Sci. 73, 464-482.
#' 
#' Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for
#' data-poor fisheries; forecasting with few data. ICES J. Mar. Sci.
#' doi:10.1093/icesjms/fst232
#' @describeIn Islope1 The least biologically precautionary of two constant index / CPUE methods
#' proposed by Geromont and Butterworth 2014. Tested by Carruthers et al. 2015
#' @export Islope1
Islope1 <- function(x, Data, reps = 100, yrsmth = 5, lambda = 0.4,xx = 0.2) {
  dependencies = "Data@Year, Data@Cat, Data@CV_Cat, Data@Ind"
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  C_dat <- Data@Cat[x, ind]
  if (is.na(Data@MPrec[x]) || length(Data@Year) == ylast + 1) {
    TACstar <- (1 - xx) * trlnorm(reps, mean(C_dat), Data@CV_Cat/(yrsmth^0.5))
  } else {
    TACstar <- rep(Data@MPrec[x], reps)
  }
  
  I_hist <- Data@Ind[x, ind]
  yind <- 1:yrsmth
  # slppar <- summary(lm(I_hist ~ yind))$coefficients[2, 1:2]
  slppar <- summary(lm(log(I_hist) ~ yind))$coefficients[2, 1:2]
  Islp <- rnorm(reps, slppar[1], slppar[2])
  TAC <- TACstar * (1 + lambda * Islp)
  
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
  
}
class(Islope1) <- "MP"



#' @describeIn Islope1 The most biologically precautionary of two constant index / CPUE methods
#' proposed by Geromont and Butterworth 2014. Tested by Carruthers et al. 2015
#' @export Islope4
Islope4 <- function(x, Data, reps = 100, yrsmth = 5, lambda = 0.2, 
                    xx = 0.4) {
  dependencies = "Data@Year, Data@Cat, Data@CV_Cat, Data@Ind"
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  C_dat <- Data@Cat[x, ind]
  if (is.na(Data@MPrec[x]) || length(Data@Year) == ylast + 1) {
    TACstar <- (1 - xx) * trlnorm(reps, mean(C_dat), Data@CV_Cat/(yrsmth^0.5))
  } else {
    TACstar <- rep(Data@MPrec[x], reps)
  }
  I_hist <- Data@Ind[x, ind]
  yind <- 1:yrsmth
  # slppar <- summary(lm(I_hist ~ yind))$coefficients[2, 1:2]
  slppar <- summary(lm(log(I_hist) ~ yind))$coefficients[2, 1:2]
  Islp <- rnorm(reps, slppar[1], slppar[2])
  TAC <- TACstar * (1 + lambda * Islp)
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(Islope4) <- "MP"


#' Index Target Management Procedure
#' 
#' An index target MP where the TAC is modified according to current index
#' levels (mean index over last 5 years) relative to a target level. Maximum
#' annual changes are 5 per cent.
#' 
#' 
#' @usage IT5(x, Data, reps = 100,yrsmth=5,mc=0.05)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the quota recommendation
#' @param yrsmth The number of historical years over which to average the index
#' @param mc The maximum fractional change in the TAC among years.
#' @return A numeric vector of TAC recommendations
#' @author T. Carruthers
#' @describeIn IT5 Maximum annual changes in TAC are 5 per cent.
#' @export IT5
IT5 <- function(x, Data, reps = 100, yrsmth = 5, mc = 0.05) {
  dependencies = "Data@Ind, Data@MPrec, Data@CV_Ind, Data@Iref"
  ind <- max(1, (length(Data@Year) - yrsmth + 1)):length(Data@Year)
  deltaI <- mean(Data@Ind[x, ind])/Data@Iref[x]
  if (deltaI < (1 - mc)) deltaI <- 1 - mc
  if (deltaI > (1 + mc)) deltaI <- 1 + mc
  TAC <- Data@MPrec[x] * deltaI * trlnorm(reps, 1, Data@CV_Ind[x])
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(IT5) <- "MP"

#' @describeIn IT5 Maximum annual changes are 10 per cent.
#' @export IT10
IT10 <- function(x, Data, reps = 100, yrsmth = 5, mc = 0.1) {
  dependencies = "Data@Ind, Data@MPrec, Data@CV_Ind, Data@Iref"
  ind <- max(1, (length(Data@Year) - yrsmth + 1)):length(Data@Year)
  deltaI <- mean(Data@Ind[x, ind])/Data@Iref[x]
  if (deltaI < (1 - mc)) deltaI <- 1 - mc
  if (deltaI > (1 + mc)) deltaI <- 1 + mc
  TAC <- Data@MPrec[x] * deltaI * trlnorm(reps, 1, Data@CV_Ind[x])
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(IT10) <- "MP"


#' A management procedure that incrementally adjusts the TAC (starting from
#' reference level that is a fraction of mean recent catches) or effort (from
#' a reference level) to reach a target CPUE / relative abundance index
#' 
#' Four index/CPUE target MPs proposed by Geromont and Butterworth 2014. 
#' Tested by Carruthers et al. 2015.
#' 
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param Imulti Parameter controlling how much larger target CPUE / index is
#' compared with recent levels.
#' @param plot Logical. Show the plot?
#' @param ... other arguments passed to plotting function
#' @return A Rec object
#' @author T. Carruthers
#' @references Carruthers et al. 2015. Performance evaluation of simple
#' management procedures. ICES J. Mar Sci. 73, 464-482.
#' 
#' Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for
#' data-poor fisheries; forecasting with few data. ICES J. Mar. Sci. 72, 251-261.
#' doi:10.1093/icesjms/fst232
#' @describeIn Itarget1 The less precautionary TAC-based MP
#' @export Itarget1
Itarget1 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0, Imulti = 1.5, plot=FALSE, ...) {
  dependencies = "Data@Cat, Data@CV_Cat"
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 5 years
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  ind2 <- ((ylast - (yrsmth - 1)):ylast)  # historical 5 pre-projection years
  ind3 <- ((ylast - (yrsmth * 2 - 1)):ylast)  # historical 10 pre-projection years
  C_dat <- Data@Cat[x, ind2]
  TACstar <- (1 - xx) * trlnorm(reps, mean(C_dat), Data@CV_Cat/(yrsmth^0.5))
  Irecent <- mean(Data@Ind[x, ind])
  Iave <- mean(Data@Ind[x, ind3])
  Itarget <- Iave * Imulti
  I0 <- 0.8 * Iave
  if (Irecent > I0) {
    TAC <- 0.5 * TACstar * (1 + ((Irecent - I0)/(Itarget - I0)))
  } else {
    TAC <- 0.5 * TACstar * (Irecent/I0)^2
  }
  TAC <- TACfilter(TAC)
  
  if (plot) {
    environment(Itarget_p) <- environment()
    Itarget_p(...)
  }
    
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}
class(Itarget1) <- "MP"




#' @describeIn Itarget1 The most biologically precautionary TAC-based MP
#' @export Itarget4
Itarget4 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0.3, Imulti = 2.5, plot=FALSE, ...) {
  dependencies = "Data@Cat, Data@CV_Cat"
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 5 years
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  ind2 <- ((ylast - (yrsmth - 1)):ylast)  # historical 5 pre-projection years
  ind3 <- ((ylast - (yrsmth * 2 - 1)):ylast)  # historical 10 pre-projection years
  C_dat <- Data@Cat[x, ind2]
  TACstar <- (1 - xx) * trlnorm(reps, mean(C_dat), Data@CV_Cat/(yrsmth^0.5))
  Irecent <- mean(Data@Ind[x, ind])
  Iave <- mean(Data@Ind[x, ind3])
  Itarget <- Iave * Imulti
  I0 <- 0.8 * Iave
  if (Irecent > I0) {
    TAC <- 0.5 * TACstar * (1 + ((Irecent - I0)/(Itarget - I0)))
  } else {
    TAC <- 0.5 * TACstar * (Irecent/I0)^2
  }
  
  if (plot) {
    environment(Itarget_p) <- environment()
    Itarget_p(...)
  }
    
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(Itarget4) <- "MP"


#' Index Target based on natural mortality rate
#' 
#' An index target MP where the TAC is modified according to current index
#' levels (mean index over last yrsmth years) relative to a target level.
#' Maximum fractional annual changes are mc.  mc=(5+M*25)/100
#' yrsmth=4*(1/M)^(0.25)
#' 
#' 
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the quota recommendation
#' @return A numeric vector of TAC recommendations
#' @author T. Carruthers
#' @export 
ITM <- function(x, Data, reps = 100) {
  
  dependencies = "Data@Ind, Data@Cat, Data@CV_Ind, Data@Iref, Data@Mort"
  mc <- (5 + Data@Mort[x] * 25)/100
  if (mc > 0.2) 
    mc <- 0.2
  yrsmth <- floor(4 * (1/Data@Mort[x])^(1/4))
  ind <- max(1, (length(Data@Year) - yrsmth + 1)):length(Data@Year)
  
  deltaI <- mean(Data@Ind[x, ind])/Data@Iref[x]
  if (deltaI < (1 - mc)) deltaI <- 1 - mc
  if (deltaI > (1 + mc)) deltaI <- 1 + mc
  
  TAC <- Data@MPrec[x] * deltaI * trlnorm(reps, 1, Data@CV_Ind[x])
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(ITM) <- "MP"

#' A management procedure that adjusts the TAC up/down from reference (target)
#' level (that is a fraction of mean recent premanagement catches) to reach a
#' target mean length of fish caught.
#' 
#' This MP is based on Ltarget1 proposed by Geromont and Butterworth 2014, but
#' here the target and limit mean lengths are based on the length at maturity
#' distribution rather than an arbitrary multiplicative of the mean length.
#' 
#' 
#' @usage L95target(x, Data, reps = 100, yrsmth = 5, buffer=0)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrsmth Years over which to calculate the mean historical catch
#' @param buffer Parameter controlling the fraction of mean catch to set the
#' reference (or target) TAC level - acts as a precautionary buffer
#' @return A numeric vector of TAC recommendations
#' @author HF Geromont
#' @references Geromont, H.F., Butterworth, D.S. 2014. Generic management
#' procedures for data-poor fisheries; forecasting with few data. ICES J. Mar.
#' Sci. doi:10.1093/icesjms/fst232
#' @seealso \link{Ltarget1}
#' @export L95target
L95target <- function(x, Data, reps = 100, yrsmth = 5, buffer = 0) {
  
  # Set target length to L95 and limit to 0.9L50 Alternative targets
  # TACstar depending on buffer value
  
  dependencies = "Data@Cat, Data@Year, Data@LHYear, Data@ML, Data@L50, Data@L95"
  
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 5 years
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  ind2 <- ((ylast - (yrsmth - 1)):ylast)  # historical 5 pre-projection years
  C_dat <- Data@Cat[x, ind2]
  TACstar <- (1 - buffer) * trlnorm(reps, mean(C_dat), Data@CV_Cat/(yrsmth^0.5))
  Lrecent <- mean(Data@ML[ind])
  Ltarget <- Data@L95[x]
  L0 <- 0.9 * Data@L50[x]
  if (Lrecent > L0) {
    TAC <- 0.5 * TACstar * (1 + ((Lrecent - L0)/(Ltarget - L0)))
  } else {
    TAC <- 0.5 * TACstar * (Lrecent/L0)^2
  }
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(L95target) <- "MP"


#' Mean length-based indicator MP of Jardim et al. 2015
#' 
#' The TAC is adjusted by the ratio alpha, where the numerator 
#' is the mean length of the catch (of lengths larger than Lc) and 
#' the denominator is the mean length expected at MSY. Here, Lc is the length at 
#' full selection (LFS).
#' 
#' Argument yrsmth currently takes the mean length of the most recent 3 years of data 
#' as a smoother.
#' 
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrsmth The most recent years of data to smooth the calculation of the mean length
#' @author Coded by Q. Huynh. Developed by Jardim et al. (2015)
#' @references Ernesto Jardim, Manuela Azevedo, Nuno M. Brites, Harvest control rules for 
#' data limited stocks using length-based reference points and survey biomass indices, 
#' Fisheries Research, Volume 171, November 2015, Pages 12-19, ISSN 0165-7836, 
#' https://doi.org/10.1016/j.fishres.2014.11.013
#' @describeIn Lratio_BHI Assumes M/K = 1.5 and FMSY/M = 1. Natural mortality M and von Bertalanffy 
#' K are not used in this MP (see Appendix A of Jardim et al. 2015). 
#' @seealso \link{Ltarget1} \link{L95target}
#' @export
Lratio_BHI <- function(x, Data, reps, yrsmth = 3) {
  dependencies = "Data@vb_Linf, Data@CV_vbLinf, Data@Cat, Data@CV_Cat, Data@CAL, Data@CAL_bins,
  Data@LFS, Data@CV_LFS"
  Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Lc <- trlnorm(reps, Data@LFS[x], Data@CV_LFS[x])
  Lref <- 0.75 * Lc + 0.25 * Linfc
  
  Cat <- Data@Cat[x, length(Data@Cat[x, ])]
  Cc <- trlnorm(reps, Cat, Data@CV_Cat[x])
  
  LYear <- dim(Data@CAL)[2]
  nlbin <- ncol(Data@CAL[x,,])
  mlbin <- 0.5 * (Data@CAL_bins[1:nlbin] + Data@CAL_bins[2:(nlbin + 1)])
  
  ind.year <- (LYear - yrsmth + 1):LYear
  CAL <- colSums(Data@CAL[x, ind.year, ])
  
  LSQ <- rep(NA, reps)
  for(i in 1:reps) {
    lensamp <- sample(mlbin, 0.5*sum(CAL), replace = T, prob = CAL)
    LSQ[i] <- mean(lensamp)
  }
  
  TAC <- (LSQ/Lref) * Cc
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(Lratio_BHI) <- "MP"


#' @describeIn Lratio_BHI More general version that calculates the reference mean length as a
#' function of M, K, and presumed FMSY/M.
#' @export
Lratio_BHI2 <- function(x, Data, reps, yrsmth = 3) {
  
  dependencies = "Data@vb_Linf, Data@CV_vbLinf, Data@Cat, Data@CV_Cat, Data@Mort, Data@CV_Mort,
  Data@vb_K, Data@CV_vbK, Data@FMSY_M, Data@CV_FMSY_M, Data@CAL, Data@CAL_bins,
  Data@LFS, Data@CV_LFS"
  Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Lc <- trlnorm(reps, Data@LFS[x], Data@CV_LFS[x])
  
  Mvec <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  Kvec <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  gamma <- trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])
  theta <- Kvec/Mvec
  
  Lref <- (theta * Linfc + Lc * (gamma + 1)) / (gamma + theta + 1)
  
  Cat <- Data@Cat[x, length(Data@Cat[x, ])]
  Cc <- trlnorm(reps, Cat, Data@CV_Cat[x])
  
  LYear <- dim(Data@CAL)[2]
  nlbin <- ncol(Data@CAL[x,,])
  mlbin <- 0.5 * (Data@CAL_bins[1:nlbin] + Data@CAL_bins[2:(nlbin + 1)])
  
  ind.year <- (LYear - yrsmth + 1):LYear
  CAL <- colSums(Data@CAL[x, ind.year, ])
  
  LSQ <- rep(NA, reps)
  for(i in 1:reps) {
    lensamp <- sample(mlbin, 0.5*sum(CAL), replace = T, prob = CAL)
    LSQ[i] <- mean(lensamp)
  }
  
  TAC <- (LSQ/Lref) * Cc
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(Lratio_BHI2) <- "MP"


#' A management procedure that incrementally adjusts the TAC or effort according to the
#' mean length of recent catches.
#' 
#' Four adaptive length-based MPs proposed by Geromont and Butterworth 2014. 
#' Tested by Carruthers et al. 2015.
#' 
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param stepsz Parameter controlling the size of update increment in TAC or effort.
#' @param llim A vector of length reference points that determine the
#' conditions for increasing, maintaining or reducing the TAC or effort.
#' @return A \linkS4class{Rec} object
#' @author T. Carruthers
#' @references Carruthers et al. 2015. Performance evaluation of simple
#' management procedures. ICES J. Mar Sci. 73, 464-482.
#' 
#' Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for
#' data-poor fisheries; forecasting with few data. ICES J. Mar. Sci.
#' doi:10.1093/icesjms/fst232
#' @describeIn LstepCC1 The least biologically precautionary TAC-based MP.
#' @export LstepCC1
LstepCC1 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0, stepsz = 0.05, 
                     llim = c(0.96, 0.98, 1.05)) {
  dependencies = "Data@Cat, Data@CV_Cat, Data@ML"
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 5 years
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  ind2 <- ((ylast - (yrsmth - 1)):ylast)  # historical 5 pre-projection years
  ind3 <- ((ylast - (yrsmth * 2 - 1)):ylast)  # historical 10 pre-projection years
  C_dat <- Data@Cat[x, ind2]
  if (is.na(Data@MPrec[x]) || length(Data@Year) == ylast + 1) {
    TACstar <- (1 - xx) * trlnorm(reps, mean(C_dat), Data@CV_Cat/(yrsmth^0.5))
  } else {
    TACstar <- rep(Data@MPrec[x], reps)
  }
  step <- stepsz * TACstar
  Lrecent <- mean(Data@ML[ind])
  Lave <- mean(Data@ML[ind3])
  rat <- Lrecent/Lave
  if (rat < llim[1]) {
    TAC <- TACstar - 2 * step
  } else if (rat < llim[2]) {
    TAC <- TACstar - step
  } else if (rat > llim[3]) {
    TAC <- TACstar + step
  } else {
    TAC <- TACstar
  }
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(LstepCC1) <- "MP"




#' @describeIn LstepCC1 The most precautionary TAC-based MP.
#' @export LstepCC4
LstepCC4 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0.3, stepsz = 0.05, 
                     llim = c(0.96, 0.98, 1.05)) {
  dependencies = "Data@Cat, Data@CV_Cat, Data@ML"
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 5 years
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  ind2 <- ((ylast - (yrsmth - 1)):ylast)  # historical 5 pre-projection years
  ind3 <- ((ylast - (yrsmth * 2 - 1)):ylast)  # historical 10 pre-projection years
  C_dat <- Data@Cat[x, ind2]
  if (is.na(Data@MPrec[x]) || length(Data@Year) == ylast + 1) {
    TACstar <- (1 - xx) * trlnorm(reps, mean(C_dat), Data@CV_Cat/(yrsmth^0.5))
  } else {
    TACstar <- rep(Data@MPrec[x], reps)
  }
  step <- stepsz * TACstar
  Lrecent <- mean(Data@ML[ind])
  Lave <- mean(Data@ML[ind3])
  rat <- Lrecent/Lave
  if (rat < llim[1]) {
    TAC <- TACstar - 2 * step
  } else if (rat < llim[2]) {
    TAC <- TACstar - step
  } else if (rat > llim[3]) {
    TAC <- TACstar + step
  } else {
    TAC <- TACstar
  }
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(LstepCC4) <- "MP"




#' A management procedure that incrementally adjusts the TAC or effort to reach 
#' a target mean length in catches.
#' 
#' Four target length MPs proposed by Geromont and Butterworth 2014.
#' Tested by Carruthers et al. 2015.
#' 
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param xL Parameter controlling the magnitude of the target mean length of
#' catches relative to average length in catches.
#' @return A \linkS4class{Rec} object
#' @author T. Carruthers
#' @references Carruthers et al. 2015. Performance evaluation of simple
#' management procedures. ICES J. Mar Sci. 73, 464-482.
#' 
#' Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for
#' data-poor fisheries; forecasting with few data. ICES J. Mar. Sci.
#' doi:10.1093/icesjms/fst232
#' @describeIn Ltarget1 The least biologically precautionary TAC-based MP.
#' @seealso \link{L95target}
#' @export Ltarget1
Ltarget1 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0, xL = 1.05) {
  dependencies = "Data@Cat, Data@CV_Cat, Data@ML"
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 5 years
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  ind2 <- ((ylast - (yrsmth - 1)):ylast)  # historical 5 pre-projection years
  ind3 <- ((ylast - (yrsmth * 2 - 1)):ylast)  # historical 10 pre-projection years
  C_dat <- Data@Cat[x, ind2]
  TACstar <- (1 - xx) * trlnorm(reps, mean(C_dat), Data@CV_Cat/(yrsmth^0.5))
  Lrecent <- mean(Data@ML[ind])
  Lave <- mean(Data@ML[ind3])
  L0 <- 0.9 * Lave
  Ltarget <- xL * Lave
  if (Lrecent > L0) {
    TAC <- 0.5 * TACstar * (1 + ((Lrecent - L0)/(Ltarget - L0)))
  } else {
    TAC <- 0.5 * TACstar * (Lrecent/L0)^2
  }
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(Ltarget1) <- "MP"



#' @describeIn Ltarget1 The most biologically precautionary TAC-based MP.
#' @export Ltarget4
Ltarget4 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0.2, xL = 1.15) {
  dependencies = "Data@Cat, Data@CV_Cat, Data@ML"
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 5 years
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  ind2 <- ((ylast - (yrsmth - 1)):ylast)  # historical 5 pre-projection years
  ind3 <- ((ylast - (yrsmth * 2 - 1)):ylast)  # historical 10 pre-projection years
  C_dat <- Data@Cat[x, ind2]
  TACstar <- (1 - xx) * trlnorm(reps, mean(C_dat), Data@CV_Cat/(yrsmth^0.5))
  Lrecent <- mean(Data@ML[ind])
  Lave <- mean(Data@ML[ind3])
  L0 <- 0.9 * Lave
  Ltarget <- xL * Lave
  if (Lrecent > L0) {
    TAC <- 0.5 * TACstar * (1 + ((Lrecent - L0)/(Ltarget - L0)))
  } else {
    TAC <- 0.5 * TACstar * (Lrecent/L0)^2
  }
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(Ltarget4) <- "MP"


#' Mean Catch Depletion
#' 
#' A simple average catch-depletion MP that was included to demonstrate just
#' how informative an estimate of current stock depletion can be. 
#'
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the quota recommendation
#' @return A Rec object with TAC recommendations
#' @describeIn MCD The calculated TAC = 2 \* depletion \* AvC
#' @author T. Carruthers
#' @export MCD
MCD <- function(x, Data, reps = 100) {
  # Daft method to demonstrate the relative value of information of
  # current depletion
  dependencies = "Data@Dep, Data@CV_Dep, Data@Cat"
  if (is.na(Data@Dep[x])) 
    return(NA)
  depo <- max(0.01, min(0.99, Data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
  Bt_K <- rbeta(reps * 100, alphaconv(depo, min(depo * Data@CV_Dep[x], 
                                                (1 - depo) * Data@CV_Dep[x])), betaconv(depo, min(depo * Data@CV_Dep[x], 
                                                                                                  (1 - depo) * Data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
  Bt_K <- Bt_K[Bt_K > 0.00999 & Bt_K < 0.99001][1:reps]  # interval censor (0.01,0.99)  as in Dick and MacCall 2011
  AvC <- stats::rlnorm(reps, log(mean(Data@Cat[x, ], na.rm = T)), Data@CV_Cat[x])
  TAC <- AvC * 2 * Bt_K
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(MCD) <- "MP"



#' @describeIn MCD Linked to a 40-10 harvest control rule
#' @author T. Carruthers
#' @export MCD4010
MCD4010 <- function(x, Data, reps = 100) {
  # Daft method to demonstrate the relative value of information of
  # current depletion
  dependencies = "Data@Dep, Data@CV_Dep, Data@Cat"
  if (is.na(Data@Dep[x])) 
    return(NA)
  depo <- max(0.01, min(0.99, Data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
  Bt_K <- rbeta(reps * 100, alphaconv(depo, min(depo * Data@CV_Dep[x], 
                                                (1 - depo) * Data@CV_Dep[x])), betaconv(depo, min(depo * Data@CV_Dep[x], 
                                                                                                  (1 - depo) * Data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
  Bt_K <- Bt_K[Bt_K > 0.00999 & Bt_K < 0.99001][1:reps]  # interval censor (0.01,0.99)  as in Dick and MacCall 2011
  AvC <- stats::rlnorm(reps, log(mean(Data@Cat[x, ], na.rm = T)), Data@CV_Cat[x])
  TAC <- AvC * 2 * Bt_K
  
  # 40-10 HCR
  cond1 <- Bt_K < 0.4 & Bt_K > 0.1
  cond2 <- Bt_K < 0.1
  TAC[cond1] <- TAC[cond1] * (Bt_K[cond1] - 0.1)/0.3
  TAC[cond2] <- TAC[cond2] * tiny  # this has to still be stochastic albeit very small
  
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(MCD4010) <- "MP"



#' Harvest Control Rule using prior for intrinsic rate of increase
#' 
#' An MP proposed by Carl Walters that modifies TACs according to trends in
#' apparent surplus production that includes information from a demographically
#' derived prior for intrinsic rate of increase
#' 
#' 
#' @usage Rcontrol(x, Data, reps = 100, yrsmth = 10, gg = 2, glim = c(0.5,
#' 2))
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of quota samples
#' @param yrsmth The number of years for smoothing catch and biomass data
#' @param gg A gain parameters
#' @param glim Limits for the change in TAC among years
#' @author C. Walters and T. Carruthers
#' @references Made-up for this package.
#' @export Rcontrol
Rcontrol <- function(x, Data, reps = 100, yrsmth = 10, gg = 2, glim = c(0.5, 
                                                                        2)) {
  dependencies = "Data@Mort, Data@CV_Mort, Data@vbK, Data@CV_vbK, Data@vbLinf, Data@CV_vbLinf, Data@vbt0, Data@CV_vbt0 Data@steep, Data@CV_steep, Data@MaxAge, Data@Dep, Data@CV_Dep, Data@Cat, Data@Ind"
  Mvec <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  Kvec <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  Linfvec <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
    t0vec <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
  } else {
    t0vec <- rep(Data@vbt0[x], reps)
  }
  t0vec[!is.finite(t0vec)] <- 0
  
  # hvec <- trlnorm(reps, Data@steep[x], Data@CV_steep[x])
  hvec <- sample_steepness2(reps, Data@steep[x], Data@CV_steep[x])
  rsamp <- getr(x, Data, Mvec, Kvec, Linfvec, t0vec, hvec, maxage = Data@MaxAge, 
                r_reps = reps)
  
  depo <- max(0.01, min(0.99, Data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
  if (any(is.na(c(Data@Dep[x], Data@CV_Dep[x])))) 
    return(NA)
  
  Bt_K <- rbeta(100, alphaconv(depo, min(depo * Data@CV_Dep[x], 
                                         (1 - depo) * Data@CV_Dep[x])), betaconv(depo, min(depo * Data@CV_Dep[x], 
                                                                                           (1 - depo) * Data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
  Bt_K <- Bt_K[Bt_K > 0.01 & Bt_K < 0.99][1]  # interval censor (0.01,0.99)  as in Dick and MacCall 2011
  
  G_new <- rsamp * (1 - 2 * Bt_K)  # here is a big difference from SPHCR
  
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
  C_dat <- log(Data@Cat[x, ind])
  C_dat[C_dat == -Inf] <- 0
  B_dat <- log(Data@Ind[x, ind]/Data@Ind[x, ind[yrsmth]] * Data@Abun[x])
  B_dat[B_dat == -Inf] <- 0
  C_hist <- exp(predict(loess(C_dat ~ ind, degree = 1)))
  B_hist <- exp(predict(loess(B_dat ~ ind, degree = 1)))
  ind <- 2:yrsmth
  ind1 <- 1:(yrsmth - 1)
  SP_hist <- B_hist[ind] - B_hist[ind1] + C_hist[ind1]
  yind <- 1:length(SP_hist)
  SP_mu <- predict(lm(SP_hist ~ yind), newdat = list(yind = length(SP_hist) + 
                                                       1))
  SP_se <- predict(lm(SP_hist ~ yind), newdat = list(yind = length(SP_hist) + 
                                                       1), se = T)$se.fit
  SP_new <- rnorm(reps, SP_mu, SP_se/2)
  
  TAC <- SP_new * (1 - gg * G_new)
  TAC[TAC < glim[1] * C_hist[yrsmth]] <- glim[1] * C_hist[yrsmth]
  TAC[TAC > glim[2] * C_hist[yrsmth]] <- glim[2] * C_hist[yrsmth]
  
  # Carr<-cbind(array(rep(Data@Cat[x,],each=reps),c(reps,length(Data@Cat[x,]))),TAC)
  # Warr<-(Data@Mort[x]*exp(-Data@Mort[x]*(1:ncol(Carr))))[ncol(Carr):1]
  # Warr<-Warr/sum(Warr)
  # TAC<-apply(t(matrix(Warr,nrow=ncol(Carr),ncol=reps))*Carr,1,sum)
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(Rcontrol) <- "MP"



#' @describeIn Rcontrol This is different from
#' Rcontrol because it includes a quadratic approximation of recent trend in
#' surplus production given biomass
#' @export Rcontrol2
Rcontrol2 <- function(x, Data, reps = 100, yrsmth = 10, gg = 2, glim = c(0.5, 
                                                                         2)) {
  dependencies = "Data@Mort, Data@CV_Mort, Data@vbK, Data@CV_vbK, Data@vbLinf, Data@CV_vbLinf, Data@vbt0, Data@CV_vbt0, Data@steep, Data@CV_steep, Data@MaxAge, Data@Dep, Data@CV_Dep, Data@Cat, Data@Ind"
  Mvec <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  Kvec <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  Linfvec <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
    t0vec <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
  } else {
    t0vec <- rep(Data@vbt0[x], reps)
  }
  t0vec[!is.finite(t0vec)] <- 0
  # hvec <- trlnorm(reps, Data@steep[x], Data@CV_steep[x])
  hvec <- sample_steepness2(reps, Data@steep[x], Data@CV_steep[x])
  rsamp <- getr(x, Data, Mvec, Kvec, Linfvec, t0vec, hvec, maxage = Data@MaxAge, 
                r_reps = reps)
  
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
  C_dat <- log(Data@Cat[x, ind])
  C_dat[C_dat == -Inf] <- 0
  B_dat <- log(Data@Ind[x, ind]/Data@Ind[x, ind[yrsmth]] * Data@Abun[x])
  B_dat[B_dat == -Inf] <- 0
  C_hist <- exp(predict(loess(C_dat ~ ind, degree = 1)))
  B_hist <- exp(predict(loess(B_dat ~ ind, degree = 1)))
  ind <- 2:yrsmth
  ind1 <- 1:(yrsmth - 1)
  SP_hist <- B_hist[ind] - B_hist[ind1] + C_hist[ind1]
  yind <- 1:length(SP_hist)
  SP_mu <- predict(lm(SP_hist ~ yind), newdat = list(yind = length(SP_hist) + 
                                                       1))
  SP_se <- predict(lm(SP_hist ~ yind), newdat = list(yind = length(SP_hist) + 
                                                       1), se = T)$se.fit
  SP_new <- rnorm(reps, SP_mu, SP_se/2)
  SParr <- array(rep(SP_hist, each = reps), dim = c(reps, yrsmth - 1))
  Barr <- array(rep(B_hist[ind], each = reps), dim = c(reps, yrsmth - 
                                                         1))
  rarr <- array(rep(rsamp, yrsmth - 1), dim = c(reps, yrsmth - 1))
  b2 <- apply(SParr/Barr - rarr, 1, sum) * apply(Barr, 1, sum)/apply(Barr^2, 
                                                                     1, sum)
  G_new <- rsamp - 2 * b2 * B_hist[yrsmth]
  
  TAC <- SP_new * (1 - gg * G_new)
  TAC[TAC < glim[1] * C_hist[yrsmth]] <- glim[1] * C_hist[yrsmth]
  TAC[TAC > glim[2] * C_hist[yrsmth]] <- glim[2] * C_hist[yrsmth]
  # Carr<-cbind(array(rep(Data@Cat[x,],each=reps),c(reps,length(Data@Cat[x,]))),TAC)
  # Warr<-(Data@Mort[x]*exp(-Data@Mort[x]*(1:ncol(Carr))))[ncol(Carr):1]
  # Warr<-Warr/sum(Warr)
  # TAC<-apply(t(matrix(Warr,nrow=ncol(Carr),ncol=reps))*Carr,1,sum)
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(Rcontrol2) <- "MP"


#' SBT simple MP
#' 
#' An MP that makes incremental adjustments to TAC recommendations based on the
#' apparent trend in CPUE
#' 
#' This isn't exactly the same as the proposed methods and is stochastic in
#' this implementation. The method doesn't tend to work too well under many
#' circumstances possibly due to the lack of 'tuning' that occurs in the real
#' SBT assessment environment. You could try asking Rich Hillary at CSIRO about
#' this approach.
#' 
#' @usage SBT1(x, Data, reps = 100, yrsmth=10, k1=1.5, k2=3, gamma=1)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of samples of the TAC
#' @param yrsmth The number of years for evaluating trend in relative abundance
#' indices
#' @param k1 Control parameter
#' @param k2 Control parameter
#' @param gamma Control parameter
#' @author T. Carruthers
#' @references http://www.ccsbt.org/site/recent_assessment.php
#' @seealso \link{SBT2}
#' @export SBT1
SBT1 <- function(x, Data, reps = 100, yrsmth = 10, k1 = 1.5, k2 = 3, gamma = 1) {
  dependencies = "Data@Cat, Data@Year, Data@Ind"
  Cr <- length(Data@Cat[x, ])
  cct <- trlnorm(reps, Data@Cat[x, Cr], Data@CV_Cat)
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
  I_hist <- Data@Ind[x, ind]
  test <- summary(lm(I_hist ~ ind))$coefficients[2, 1:2]
  lambda <- rnorm(reps, test[1], test[2])
  # TAC <- cct * 1 + k2 * lambda
  # see https://github.com/DLMtool/DLMtool/issues/17
  TAC <- cct * (1 + k2 * lambda)
  cond <- lambda < 0
  # TAC[cond] <- cct[cond] * 1 - k1 * -lambda[cond]^gamma
  TAC[cond] <- cct[cond] * (1 - k1 * -lambda[cond]^gamma)
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(SBT1) <- "MP"



#' SBT complex MP
#' 
#' An MP that makes incremental adjustments to TAC recommendations based on
#' index levels relative to target levels (BMSY/B0) and catch levels relative
#' to target levels (MSY)
#' 
#' This isn't exactly the same as the proposed methods and is stochastic in
#' this implementation. The method doesn't tend to work too well under many
#' circumstances possibly due to the lack of 'tuning' that occurs in the real
#' SBT assessment environment. You could try asking Rich Hillary at CSIRO about
#' this approach.
#' 
#' @usage SBT2(x, Data, reps = 100,
#' epsB=0.25,epsR=0.75,tauR=5,tauB=7,gamma=1)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of samples of the TAC
#' @param epsB Control parameter
#' @param epsR Control parameter
#' @param tauR Control parameter
#' @param tauB Control parameter
#' @param gamma Control parameter
#' @author T. Carruthers
#' @references http://www.ccsbt.org/site/recent_assessment.php
#' @seealso \link{SBT1}
#' @export SBT2
SBT2 <- function(x, Data, reps = 100, epsB = 0.25, epsR = 0.75, tauR = 5, 
                 tauB = 7, gamma = 1) {
  dependencies = "Data@Cref, Data@Rec, Data@Cat"
  # Bnow<-trlnorm(reps,Data@Abun[x],Data@CV_Abun)
  # testrat<-Bnow/Data@Bref Ctarg<-rep(NA,reps)
  # Ctarg[testrat>1]<-delta*testrat[testrat>1]^(1-epsB)
  # Ctarg[testrat<1]<-detla*testrat[testrat<1]^(1+epsB)
  Ctarg <- trlnorm(reps, Data@Cref[x], Data@CV_Cref)
  muR <- mean(Data@Rec[x, (length(Data@Rec[x, ]) - tauR + 1):length(Data@Rec[x, ])])
  phi <- mean(Data@Rec[x, (length(Data@Rec[x, ]) - 9):length(Data@Rec[x,])])
  Rrat <- muR/phi
  deltaR <- rep(NA, reps)
  deltaR[Rrat > 1] <- Rrat[Rrat > 1]^(1 - epsR)
  deltaR[Rrat < 1] <- Rrat[Rrat < 1]^(1 + epsR)
  TAC <- 0.5 * (Data@Cat[x, length(Data@Cat[x, ])] + Ctarg *  deltaR)
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(SBT2) <- "MP"

#' Surplus production based catch-limit modifier
#' 
#' An MP that makes incremental adjustments to TAC recommendations based on the
#' apparent trend in surplus production. Based on the theory of Mark Maunder
#' (IATTC)
#' 
#' Note that this isn't exactly what Mark has previously suggested and is
#' stochastic in this implementation.
#' 
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of quota samples
#' @param alp Condition for modifying the TAC (bounds on change in abundance)
#' @param bet Limits for how much the TAC can change among years
#' @return A numeric vector of TAC recommendations
#' @author T. Carruthers
#' @references
#' http://www.iattc.org/Meetings/Meetings2014/MAYSAC/PDFs/SAC-05-10b-Management-Strategy-Evaluation.pdf
#' @seealso \link{SPslope}
#' @export 
SPmod <- function(x, Data, reps = 100, alp = c(0.8, 1.2), bet = c(0.8, 1.2)) {
  dependencies = "Data@Cat, Data@Ind, Data@Abun, Data@CV_Ind, Data@CV_Cat,  Data@CV_Abun"
  Ir <- length(Data@Ind[x, ])
  Cr <- length(Data@Cat[x, ])
  rat <- trlnorm(reps, Data@Ind[x, Ir], Data@CV_Ind[x])/trlnorm(reps, Data@Ind[x, Ir - 1], Data@CV_Ind[x])
  cct <- trlnorm(reps, Data@Cat[x, Cr], Data@CV_Cat[x])
  Abun <- trlnorm(reps, Data@Abun[x], Data@CV_Abun[x])
  TAC <- rep(NA, reps)
  TAC[rat < alp[1]] <- cct[rat < alp[1]] * bet[1]
  TAC[rat > alp[1] & rat < alp[2]] <- cct[rat > alp[1] & rat < alp[2]]
  
  cond <- rat > alp[2]
  reps2 <- sum(cond)
  if (reps2 > 0) {
    qq1 <- trlnorm(reps2, Data@Ind[x, Ir]/Abun[cond], Data@CV_Ind[x])
    bio1 <- Data@Ind[x, Ir - 1]/qq1
    bio2 <- Data@Ind[x, Ir]/qq1
    cct1 <- trlnorm(reps2, Data@Cat[x, Cr - 1], Data@CV_Cat[x])
    PP <- bio2 - bio1 + cct1
    TAC[cond] <- bet[2] * PP
  }
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(SPmod) <- "MP"


#' Catch trend Surplus Production MSY MP
#' 
#' An MP that uses Martell and Froese (2012) method for estimating MSY to
#' determine the OFL. Since their approach estimates stock trajectories based
#' on catches and a rule for intrinsic rate of increase it also returns
#' depletion. Given their surplus production model predicts K, r and depletion
#' it is straighforward to calculate the OFL based on the Schaefer productivity
#' curve. OFL = dep \* (1-dep) \* r \* K \* 2
#' 
#' Requires the assumption that catch is proportional to abundance.
#' Occasionally the rule that limits r and K ranges does not allow r-K pairs to
#' be found that lead to the depletion inferred by the catch trajectories. In
#' this case this method widens the search.
#' 
#' @usage SPMSY(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of samples of the TAC
#' @author T. Carruthers
#' @references Martell, S. and Froese, R. 2012. A simple method for estimating
#' MSY from catch and resilience. Fish and Fisheries. DOI:
#' 10.1111/j.1467-2979.2012.00485.x
#' @export SPMSY
SPMSY <- function(x, Data, reps = 100) {
  # Martell and Froese 2012 Schaefer SP estimate of MSY given priors on
  # r, k and depletion for(x in 1:100){
  dependencies = "Data@MaxAge, Data@vbK, Data@L50, Data@Cat"
  maxage <- Data@MaxAge
  nsamp <- reps * 200
  
  # Froese 2012 http://www.fishbase.de/rfroese/Catch-MSY_SummaryFinal.doc
  rule <- rep(4, 3)
  
  if (Data@vbK[x] > 0.3) {
    # K rules
    rule[1] <- 1
  } else if (Data@vbK[x] < 0.3 & Data@vbK[x] > 0.16) {
    rule[1] <- 2
  } else if (Data@vbK[x] < 0.16 & Data@vbK[x] > 0.05) {
    rule[1] <- 3
  }
  AM <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], Data@L50[x])
  if (AM < 1.5) {
    # Age at maturity rules
    rule[2] <- 1
  } else if (AM < 4.5 & AM > 1.5) {
    rule[2] <- 2
  } else if (AM < 10 & AM > 4.5) {
    rule[2] <- 3
  }
  
  if (Data@MaxAge < 4) {
    # Maximum age rules
    rule[3] <- 1
  } else if (Data@MaxAge < 11 & Data@MaxAge > 3) {
    rule[3] <- 2
  } else if (Data@MaxAge < 31 & Data@MaxAge > 10) {
    rule[3] <- 3
  }
  
  if (mean(rule) < 1.5)   rsamp <- runif(nsamp, 0.6, 1.5)
  if (mean(rule) > 1.5 & mean(rule) < 2.5)   rsamp <- runif(nsamp, 0.2, 1)
  if (mean(rule) > 2.5 & mean(rule) < 3.5)  rsamp <- runif(nsamp, 0.05, 0.5)
  if (mean(rule) > 3.5) rsamp <- runif(nsamp, 0.015, 0.1)
  
  Ksamp <- runif(nsamp, mean(Data@Cat[x, ])/rsamp, (10 * mean(Data@Cat[x, ]))/rsamp)
  nyears <- length(Data@Cat[x, ])
  B <- array(NA, dim = c(nsamp, nyears))
  
  if (Data@Cat[x, 1] < (0.5 * max(Data@Cat[x, ]))) {
    # Martell and Froese decision rules (makes absolutely no sense to me!)
    B[, 1] <- Ksamp * runif(nsamp, 0.5, 0.9)
  } else {
    B[, 1] <- Ksamp * runif(nsamp, 0.3, 0.6)
  }
  
  if (Data@Cat[x, nyears] < (0.5 * max(Data@Cat[x, ]))) {
    # Martell and Froese decision rules (makes absolutely no sense to me!)
    LB <- 0.01
    UB <- 0.4
  } else {
    LB <- 0.3
    UB <- 0.7
  }
  
  for (i in 2:nyears) {
    B[, i] <- B[, i - 1] - Data@Cat[x, i - 1]
    B[, i] <- B[, i] + rsamp * B[, i] * (1 - B[, i]/Ksamp)
  }
  B <- B/rep(Ksamp, nyears)
  cond <- (B[, nyears] >= LB) & (B[, nyears] <= UB)
  if (sum(cond) < 1) {
    B[B[, nyears] >= UB, nyears] <- UB
    cond <- (B[, nyears] >= LB) & (B[, nyears] <= UB)
  }
  dep <- B[cond, nyears][1:reps]
  MSY <- rsamp[cond][1:reps] * Ksamp[cond][1:reps]/4
  Kc <- Ksamp[cond][1:reps]
  rc <- rsamp[cond][1:reps]
  TAC <- Kc * dep * rc/2
  
  if (sum(!is.na(TAC)) < ceiling(reps/10)) {
    # a fudge of the original method that widens current depletion to the
    # lowest and higest bounds to get an TAC sample
    cond <- (B[, nyears] >= 0.01) & (B[, nyears] <= 0.7)
    dep <- B[cond, nyears][1:reps]
    MSY <- rsamp[cond][1:reps] * Ksamp[cond][1:reps]/4
    Kc <- Ksamp[cond][1:reps]
    rc <- rsamp[cond][1:reps]
    TAC <- Kc * dep * rc/2
    
  }
  # }
  
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}  # end of SPMSY
class(SPMSY) <- "MP"

#' Slope in surplus production MP
#' 
#' A management procedure that makes incremental adjustments to TAC
#' recommendations based on the apparent trend in recent surplus production.
#' Based on the theory of Mark Maunder (IATTC)
#' 
#' Note that this isn't exactly what Mark has previously suggested and is
#' stochastic in this implementation.
#' 
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of quota samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param alp Condition for modifying the Data (bounds on change in
#' abundance)
#' @param bet Limits for how much the Data can change among years
#' @return A numeric vector of Data recommendations
#' @author T. Carruthers
#' @seealso \link{SPmod}
#' @references
#' http://www.iattc.org/Meetings/Meetings2014/MAYSAC/PDFs/SAC-05-10b-Management-Strategy-Evaluation.pdf
#' @export 
SPslope <- function(x, Data, reps = 100, yrsmth = 4, alp = c(0.9, 1.1), 
                    bet = c(1.5, 0.9)) {
  
  dependencies = "Data@Year, Data@Cat, Data@Ind, Data@Abun"
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
  yind <- 1:yrsmth
  C_dat <- Data@Cat[x, ind]
  B_dat <- Data@Ind[x, ind]/Data@Ind[x, ind[yrsmth]] * Data@Abun[x]
  Pt_mu <- max(B_dat[yrsmth] - B_dat[yrsmth - 1] + C_dat[yrsmth - 1], 
               tiny)
  Pt_1 <- trlnorm(reps, Pt_mu, Data@CV_Cat[x])
  It <- exp(predict(lm(log(B_dat) ~ yind), newdat = list(yind = yrsmth + 
                                                           1)))
  Ilast <- B_dat[yrsmth]
  MC <- max(mean(C_dat), tiny)
  Ct_1 <- trlnorm(reps, MC, Data@CV_Cat[x]/(yrsmth^0.5))  # mean catches over the interval
  
  rat <- It/Ilast
  
  mult <- max((1 - bet[1] * (Ilast - It)/Ilast), tiny)
  if (rat < alp[1]) TAC <- mult * Ct_1
  if (rat > alp[1] & rat < alp[2]) TAC <- Ct_1
  if (rat > alp[2]) TAC <- bet[2] * Pt_1
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(SPslope) <- "MP"

#' Surplus Production Stock Reduction Analysis
#' 
#' A surplus production equivalent of DB-SRA that uses a demographically
#' derived prior for intrinsic rate of increase (McAllister method, below)
#' 
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object (class DLM)
#' @param reps The number of samples of the TAC taken for the calculation of
#' the quota
#' @author T. Carruthers
#' @references McAllister, M.K., Pikitch, E.K., and Babcock, E.A. 2001. Using
#' demographic methods to construct Bayesian priors for the intrinsic rate of
#' increase in the Schaefer model and implications for stock rebuilding. Can.
#' J. Fish. Aquat. Sci. 58: 1871-1890.
#' @describeIn SPSRA Base version
#' @export 
SPSRA <- function(x, Data, reps = 100) {
  # Surplus productin stock reduction analysis T.Carruthers - basically
  # an SP version of DBSRA
  dependencies = "Data@Mort, Data@CV_Mort, Data@vbK, Data@CV_vbK, Data@vbLinf, Data@CV_vbLinf, Data@vbt0, Data@CV_vbt0, Data@Dep, Data@CV_Dep, Data@Cat, Data@steep"
  Mvec <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  Kvec <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  Linfvec <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
    t0vec <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
  } else {
    t0vec <- rep(Data@vbt0[x], reps)
  }
  t0vec[!is.finite(t0vec)] <- 0
  if (all(is.nan(t0vec))) 
    t0vec <- rep(0, reps)
  # hvec <- trlnorm(reps, Data@steep[x], Data@CV_steep[x])
  hvec <- sample_steepness2(reps, Data@steep[x], Data@CV_steep[x])
  if (all(!is.finite(hvec))) 
    return(NA)
  rsamp <- getr(x, Data, Mvec, Kvec, Linfvec, t0vec, hvec, maxage = Data@MaxAge, 
                r_reps = reps)
  dep <- trlnorm(reps, Data@Dep[x], Data@CV_Dep[x])
  Ct <- Data@Cat[x, ]
  Csamp <- array(rep(Ct, each = reps) * trlnorm(length(Ct) * reps, 1, 
                                                Data@CV_Cat[x]), dim = c(reps, length(Ct)))
  Psamp <- array(trlnorm(length(Ct) * reps, 1, 0.1), dim = c(reps, length(Ct)))
  Ksamp <- rep(NA, reps)
  for (i in 1:reps) Ksamp[i] <- exp(optimize(SPSRAopt, log(c(mean(Csamp[i, 
                                                                        ]), 1000 * mean(Csamp[i, ]))), dep = dep[i], r = rsamp[i], Ct = Csamp[i, 
                                                                                                                                              ], PE = Psamp[i, ])$minimum)
  MSY <- Ksamp * rsamp/4
  TAC <- Ksamp * dep * rsamp/2
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(SPSRA) <- "MP"

#' @describeIn SPSRA Variant that uses a mean-length mortality estimator to obtain
#' a prior for current stock depletion. The mean length extension was programmed by 
#' Gary Nelson as part of his excellent R package 'fishmethods'.
#' @export SPSRA_ML
SPSRA_ML <- function(x, Data, reps = 100) {
  dependencies = "Data@Mort, Data@CV_Mort, Data@vbK, Data@CV_vbK, Data@vbLinf, Data@CV_vbLinf, Data@vbt0, Data@CV_vbt0, Data@CAL, Data@Cat, Data@steep"
  Mvec <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  Kvec <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  Linfvec = trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
    t0vec <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
  } else {
    t0vec <- rep(Data@vbt0[x], reps)
  }
  t0vec[!is.finite(t0vec)] <- 0
  # hvec <- trlnorm(reps, Data@steep[x], Data@CV_steep[x])
  hvec <- sample_steepness2(reps, Data@steep[x], Data@CV_steep[x])
  rsamp <- getr(x, Data, Mvec, Kvec, Linfvec, t0vec, hvec, maxage = Data@MaxAge, 
                r_reps = reps)
  Z <- MLne(x, Data, Linfc = Linfvec, Kc = Kvec, ML_reps = reps,  MLtype = "dep")
  if (all(is.na(Z))) {
    Rec <- new("Rec")
    Rec@TAC <- TACfilter(rep(NA, reps))
    return(Rec)
  } 
  FM <- Z - Mvec
  nyears <- length(Data@Year)
  Ct1 <- mean(Data@Cat[x, 1:3])
  Ct2 <- mean(Data@Cat[x, (nyears - 2):nyears])
  dep <- rep(c(Ct1, Ct2), each = reps)/(1 - exp(-FM))
  if (reps == 1) 
    dep <- dep[2]/dep[1]
  if (reps > 1) 
    dep <- dep[, 2]/dep[, 1]
  Ksamp <- rep(NA, reps)
  Ct <- Data@Cat[x, ]
  Csamp <- array(rep(Ct, each = reps) * trlnorm(length(Ct) * reps, 1, 
                                                Data@CV_Cat[x]), dim = c(reps, length(Ct)))
  Psamp <- array(trlnorm(length(Ct) * reps, 1, 0.1), dim = c(reps, length(Ct)))
  for (i in 1:reps) Ksamp[i] <- exp(optimize(SPSRAopt, log(c(mean(Csamp[i, 
                                                                        ]), 1000 * mean(Csamp[i, ]))), dep = dep[i], r = rsamp[i], Ct = Csamp[i, 
                                                                                                                                              ], PE = Psamp[i, ])$minimum)
  MSY <- Ksamp * rsamp/4
  TAC <- Ksamp * dep * rsamp/2
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(SPSRA_ML) <- "MP"



# A generic VPA (Walters and Licandeo UBC)
# VPA <- function(x, Data, reps = reps) {
#   
#   # now do optimization for FMSY
#   dependencies = "Data@Mort, Data@CV_Mort, Data@vbK, Data@CV_vbK, Data@vbLinf, Data@CV_vbLinf, Data@vbt0, Data@CV_vbt0, Data@MaxAge, Data@wla, Data@CV_wla, Data@wlb, Data@CV_wlb, Data@L50, Data@CV_L50, Data@CAA, Data@steep, Data@CV_steep, Data@LFS, Data@CV_LFS, Data@LFC, Data@CV_LFC, Data@Cat"
#   CAAind <- (Data@CAA[x, , ] == 0) * array(rep(1:Data@MaxAge, 
#                                                each = length(Data@CAA[x, , 1])), dim(Data@CAA[x, , ]))
#   maxage <- min(CAAind[CAAind != 0])
#   maxage <- which.min(abs(cumsum(apply(Data@CAA[x, , ], 2, sum))/sum(Data@CAA[x, 
#                                                                               , ]) - 0.75))
#   CAAv <- Data@CAA[x, , 1:maxage]
#   CAAv[, maxage] <- CAAv[, maxage] + apply(Data@CAA[x, , (maxage + 
#                                                             1):length(Data@CAA[x, 1, ])], 1, sum)
#   
#   TAC <- Bt_K <- rep(NA, reps)
#   
#   for (i in 1:reps) {
#     
#     Mc <- trlnorm(1, Data@Mort[x], Data@CV_Mort[x])
#     # hc <- trlnorm(1, Data@steep[x], Data@CV_steep[x])
#     hc <- sample_steepness2(1, Data@steep[x], Data@CV_steep[x])
#     Linfc <- trlnorm(1, Data@vbLinf[x], Data@CV_vbLinf[x])
#     Kc <- trlnorm(1, Data@vbK[x], Data@CV_vbK[x])
#     t0c <- -trlnorm(1, -Data@vbt0[x], Data@CV_vbt0[x])
#     LFSc <- trlnorm(1, Data@LFS[x], Data@CV_LFS[x])
#     LFCc <- trlnorm(1, Data@LFC[x], Data@CV_LFC[x])
#     AMc <- trlnorm(1, iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], 
#                           Data@L50[x]), Data@CV_L50[x])
#     ac <- trlnorm(1, Data@wla[x], Data@CV_wla[x])
#     bc <- trlnorm(1, Data@wlb[x], Data@CV_wlb[x])
#     
#     pmat <- rep(1, maxage)
#     pmat[1:ceiling(AMc)] <- 0
#     age <- 1:maxage
#     la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - 
#                                                        Data@vbt0[x]))))
#     wa <- ac * la^bc
#     
#     Cat <- Data@Cat[x, ]
#     Cat[1] <- Cat[2]  # temporary fix until effort simulation function gets sorted
#     
#     CAAv[, maxage][CAAv[, maxage] == 0] <- 1
#     
#     opt = optim(c(-3, -2), VPAopt, Cat = CAAv, yt = Data@Ind[x, ],
#                 S = exp(-Mc), maxage = maxage, wa = wa, pmat = pmat, method = "L-BFGS-B", 
#                 lower = c(-5, -5), upper = c(5, 5))
#     out = VPAopt(opt$par, Cat = CAAv, yt = Data@Ind[x, ], S = exp(-Data@Mort[x]), 
#                  maxage = maxage, wa = wa, pmat = pmat, opt = F)
#     
#     fit2 <- optimize(VPAFMSY, log(c(1e-04, 3)), Mc = Mc, hc = hc, maxage = maxage, 
#                      vul = out$va, Linfc = Linfc, Kc = Kc, t0c = t0c, AMc = AMc, 
#                      ac = ac, bc = bc)
#     FMSY <- VPAFMSY(fit2$minimum, Mc = Mc, hc = hc, maxage = maxage, 
#                     vul = out$va, Linfc = Linfc, Kc = Kc, t0c = t0c, AMc = AMc, 
#                     ac = ac, bc = bc, opt = F)
#     if ((FMSY/Mc) > 3) FMSY <- 3 * Mc
#     TAC[i] <- out$bt[length(out$bt)] * FMSY
#   }
#   
#   Rec <- new("Rec")
#   Rec@TAC <- TACfilter(TAC)
#   Rec
#   
# }
# class(VPA) <- "MP"
# 

#### Yield-per-recruit MPs ####

#' Internal function for YPR MPs
#'
#' @param x Iteration number
#' @param Data Object of class Data
#' @param reps Number of reps
#' @param Abun Optional numeric of length `reps` of abundance
#'
#' @export
#'
#' @keywords internal
YPR_ <- function(x, Data, reps = 100, Abun=NULL) {
  # Yield per recruit analysis F01 - Meaghan Bryan
  Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
    t0c <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
  } else {
    t0c <- rep(Data@vbt0[x], reps)
  }
  t0c[!is.finite(t0c)] <- 0
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  LFS <- trlnorm(reps, Data@LFS[x], Data@CV_LFS[x])
  a <- Data@wla[x]
  b <- Data@wlb[x]
  FMSY <- YPRopt(Linfc, Kc, t0c, Mdb, a, b, LFS, maxage=Data@MaxAge, reps)
  
  if (is.null(Abun)){
    Ac <- trlnorm(reps, Data@Abun[x], Data@CV_Abun[x])
  }else{
    Ac<-Abun
  }
  
  TAC <- Ac * FMSY
  
  return(list(TAC=TAC))
  
}



#' @describeIn YPR_ YPR calculation function
#' @note Yield per recruit estimate of FMSY Meaghan Bryan 2013
#' @param Linfc von Bertalanffy  Linf 
#' @param Kc von Bertalanffy K
#' @param t0c von Bertalanffy t0 
#' @param a alpha parameter of L-W relationship
#' @param b beta parameter of L-W relationship
#' @param LFS Length at full selection
#' @param maxage Maximum age 
YPRopt <- function(Linfc, Kc, t0c, Mdb, a, b, LFS, maxage, reps = 100) {
  
  nf <- 200
  frates <- seq(0, 3, length.out = nf)
  Winf = a * Linfc^b
  rat <- LFS/Linfc
  rat[rat > 0.8] <- 0.8  # need to robustify this for occasionally very high samples of LFS
  tc = log(1 - rat)/-Kc + t0c
  tc = round(tc, 0)
  tc[tc < 1] <- 1
  tc[tc > maxage] <- maxage
  
  vul <- array(0, dim = c(reps, maxage))
  mat <- array(0, dim = c(reps, maxage))
  lx <- array(NA, dim = c(reps, maxage))
  lxo <- array(NA, dim = c(reps, maxage))
  
  ypr <- array(NA, dim = c(reps, nf))
  sbpr <- array(NA, dim = c(reps, nf))
  sbpr.ratio <- array(NA, dim = c(reps, nf))
  sbpr.dif <- array(NA, dim = c(reps, nf))
  
  f.max <- array(NA, dim = c(reps, maxage))
  
  # average weight at age - follow von Bertalanffy growth
  age <- array(rep(1:maxage, each = reps), dim = c(reps, maxage))
  la <- Linfc * (1 - exp(-Kc * ((age - t0c))))
  wa <- a * la^b
  
  # vulnerability schedule - assumes knife-edge vulnerability, where all
  # individuals age tc to maxage are fully vulnerbale all individulas
  # less than age tc are not vulnerable
  for (i in 1:reps) {
    if (tc[i] > 0) vul[i, tc[i]:maxage] <- 1
    if (tc[i] > 1) mat[i, max(1, tc[i] - 1):maxage] <- 1
  }
  
  lx[, 1] <- 1
  lxo[, 1] <- 1
  for (k in 1:nf) {
    for (i in 2:maxage) {
      lx[, i] = lx[, i - 1] * exp(-(Mdb + vul[, i - 1] * frates[k]))
      lxo[, i] = lx[, i] * exp(-Mdb)
    }
    phi_vb = apply(lx * wa * vul, 1, sum)
    sbpro = apply(lxo * wa * mat, 1, sum)
    
    ypr[, k] = (1 - exp(-frates[k])) * phi_vb
    sbpr[, k] = apply(lx * wa * mat, 1, sum)
    sbpr.ratio[, k] = sbpr[, k]/sbpro
    sbpr.dif[, k] = abs(sbpr.ratio[, k] - 0.3)  #hard code comparison ratio
  }
  
  # frates[apply(ypr,1,which.max)] Fmaxypr
  
  # More code that derived F0.1 in 'per recruit analysis.R' (Meaghan
  # Bryan)
  slope.origin = (ypr[, 2] - ypr[, 1])/(frates[2] - frates[1])
  slope.10 = round(0.1 * slope.origin, 2)
  
  slope = array(NA, dim = dim(ypr))  #vector(length=length(ypr))
  slope[, 1] = slope.origin
  for (i in 3:ncol(ypr)) {
    slope[, i - 1] = round((ypr[, i] - ypr[, i - 1])/(frates[i] - frates[i - 1]), 2)
  }
  dif = abs(slope - slope.10)
  dif[is.na(dif)] <- 1e+11
  frates[apply(dif, 1, which.min)]  #frates[which.min(dif)]
}



## Describe, equations, and add plot ####


#' Yield Per Recruit analysis to get FMSY proxy F01
#' 
#' A simple yield per recruit approximation to FMSY (F01) which is the position
#' of the ascending YPR curve for which dYPR/dF = 0.1(dYPR/d0)
#' 
#' @templateVar mp YPR
#' @template MPtemplate
#' @template MPuses
#' 
#' @note Based on the code of Meaghan Bryan
#' @author Meaghan Bryan and Tom Carruthers
#' @references Beverton and Holt. 1954.
#' @describeIn YPR Requires an external estimate of abundance.
#' @examples 
#' YPR(1, DLMtool::Atlantic_mackerel, plot=TRUE)
#' @export 
YPR <- function(x, Data, reps = 100, plot=FALSE) {
  runYPR <- YPR_(x, Data, reps = reps, Abun=NULL)
  TAC <- TACfilter(runYPR$TAC)
  
  if (plot) YPR_plot()
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec

} 
class(YPR) <- "MP"


YPR_plot <- function() {
  
}


#' @templateVar mp YPR_CC
#' @template MPuses
#' @param Fmin The minimum fishing mortality rate inferred from the catch-curve
#' analysis
#' @describeIn YPR A catch-curve analysis is used to determine recent Z which given M (Mort)
#' gives F and thus abundance = Ct/(1-exp(-F))
#' @examples 
#' YPR_CC(1, DLMtool::SimulatedData, plot=TRUE)
#' @export 
YPR_CC <- function(x, Data, reps = 100, plot=FALSE, Fmin = 0.005) {

  MuC <- Data@Cat[x, length(Data@Cat[x, ])]
  Cc <- trlnorm(reps, MuC, Data@CV_Cat[x])
  
  Mdb <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])
  Zdb <- CC(x, Data, reps = reps * 10)
  Fdb <- Zdb - Mdb
  ind <- (1:(reps * 10))[Fdb > Fmin][1:reps]
  
  Fdb <- Fdb[ind]
  Mdb <- Mdb[ind]
  SM <- sum(is.na(ind))
  if (SM > 0) {
    Mdb[is.na(ind)] <- trlnorm(SM, Data@Mort[x], Data@CV_Mort[x])
    Fdb[is.na(ind)] <- Fmin
  }
  
  Ac <- Cc/(1 - exp(-Fdb))
  
  runYPR <- YPR_(x, Data, reps = reps, Abun=Ac)
  TAC <- TACfilter(runYPR$TAC)
  
  if (plot) YPR_plot()
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}
class(YPR_CC) <- "MP"


#' @templateVar mp YPR_ML
#' @template MPuses
#' @describeIn YPR A mean-length estimate of recent Z is used to infer current 
#' abundance.
#' @examples 
#' YPR_ML(1, DLMtool::SimulatedData, plot=TRUE) 
#' @export 
YPR_ML <- function(x, Data, reps = 100) {

  MuC <- Data@Cat[x, length(Data@Cat[x, ])]
  Cc <- trlnorm(reps * 10, MuC, Data@CV_Cat[x])
  Z <- MLne(x, Data, Linfc = Linfc, Kc = Kc, ML_reps = reps * 10, MLtype = "F")
  if (all(is.na(Z)))     return(new("Rec"))
  
  
  FM <- Z - Mdb
  ind <- which(FM>0)[1:reps]
  Ac <- Cc[ind]/(1 - exp(-FM[ind]))

  runYPR <- YPR_(x, Data, reps = reps, Abun=Ac)
  TAC <- TACfilter(runYPR$TAC)
  
  if (plot) YPR_plot()
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}
class(YPR_ML) <- "MP"








