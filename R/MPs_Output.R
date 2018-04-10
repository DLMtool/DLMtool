
# #' Example name
# #' 
# #' A brief description of the management procedure.
# #' 
# #' @details more details in paragraphs
# #' 
# #' Second paragraph
# #' 
# #' italics text: \emph{italics}
# #' bold text:  \strong{bold}
# #' r code: \code{r_function_call(with = "arguments")}, \code{NULL}, \code{TRUE}
# #'
# #' @param x A position in data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of TAC samples
# #' @author MP author
# #' @references references (can be links) to papers, etc
# #' @family other similar functions in same family
# #' @seealso links to other functions or urls
# #'   \code{\link{prod}} for products, \code{\link{cumsum}} for cumulative
# #'   sums, and \code{\link{colSums}}/\code{\link{rowSums}} marginal sums over
# #'   high-dimensional arrays.
# #' @examples 
# #' @export 
# testFunction <- function(x, Data, reps = 100) {
#   rec <- new("Rec") # create recommendation object
#   rec@TAC <- trlnorm(reps, Data@OM$A[x] * (1 - exp(-Data@OM$FMSY[x])), 0.01)
#   rec
# }




## Output Control MPs ####

#' Average Catch
#'
#' A simple average catch MP that is included to demonstrate a 'status quo' management option
#'
#' 
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @author T. Carruthers
#' 
#' @examples 
#' Data <- DLMtool::Cobia
#' # Plot the historical catches 
#' plot(Data@Year, Data@Cat[1,], type="l", 
#'      xlab="Year", ylab=paste0("Catch (", Data@Units, ")"), lwd=2)
#' abline(h=mean(Data@Cat[1,]), lty=2) # plot mean catches
#' 
#' # Apply the AvC MP to the Data
#' Rec <- AvC(1, Data, reps=1000) # 1,000 log-normal samples with CV = 0.2
#' 
#' # Distribution of TACs
#' boxplot(Rec@TAC, add=TRUE, at=max(Data@Year), col="grey", 
#'         width=1, outline=FALSE)
#'         
#' @export 
#'
#' @importFrom abind abind
#' @importFrom graphics abline axis barplot boxplot hist identify layout legend
#' lines matplot mtext par plot plot.new points polygon segments text title text
#' @importFrom grDevices col2rgb colorRampPalette rainbow rgb xy.coords
#' @importFrom methods getClassDef .hasSlot new slot slot<- slotNames
#' @importFrom stats approx coef dbeta density dnorm dlnorm lm loess loess.smooth
#' median nlm optim optimise optimize plogis pnorm predict qlnorm quantile rbeta
#' rlnorm rmultinom rnorm runif sd
#' @importFrom utils packageVersion lsf.str read.csv read.csv2
AvC <- function(x, Data, reps = 100) {
  dependencies = "Data@Cat"
  Rec <- new("Rec")
  Rec@TAC <- rlnorm(reps, log(mean(Data@Cat[x, ], na.rm = T)), 0.2)
  Rec
}
class(AvC) <- "MP"


#' Beddington and Kirkwood life-history MP
#' 
#' Family of management procedures that sets the TAC by approximation of FMSY
#' based on the length at first capture.
#' 
#' @param x A position in a data-limited methods data object.
#' @param Data A data-limited methods data object.
#' @param reps The number of stochastic samples of the TAC recommendation.
#' @param Fmin The minimum fishing mortality rate that is derived from the
#' catch-curve (interval censor).
#' @note The mean length extension was programmed by Gary Nelson as part of his
#' excellent R package 'fishmethods'
#' @describeIn BK This is the simple version of the BK MP which requires an estimate
#' of abundance. The paper has a more complex approach that might work better.
#' @author T. Carruthers.
#' @references Beddington, J.R., Kirkwood, G.P., 2005. The estimation of
#' potential yield and stock status using life history parameters. Philos.
#' Trans. R. Soc. Lond. B Biol. Sci. 360, 163-170.
#' @export BK
BK <- function(x, Data, reps = 100) {
  # Beddington and Kirkwood life-history analysis
  dependencies = "Data@LFC, Data@vbLinf, Data@CV_vbLinf, Data@Abun, Data@CV_Abun, Data@vbK, Data@CV_vbK"
  Lc <- trlnorm(reps * 10, Data@LFC[x], 0.2)
  Linfc <- trlnorm(reps * 10, Data@vbLinf[x], Data@CV_vbLinf[x])
  Ac <- trlnorm(reps * 10, Data@Abun[x], Data@CV_Abun[x])
  Kc <- trlnorm(reps * 10, Data@vbK[x], Data@CV_vbK[x])
  TAC <- Ac * (0.6 * Kc)/(0.67 - (Lc/Linfc))  # robustifying for use in MSE
  TAC <- (TAC[TAC > 0][1:reps])  # Interval censor only those positive catch recommendations
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec

}  # end of BK
class(BK) <- "MP"

#' @describeIn BK Calculates an OFL using an approximation of FMSY based 
#' on length at first capture and a catch curve estimate of current F.
#' @export BK_CC
BK_CC <- function(x, Data, reps = 100, Fmin = 0.005) {
  dependencies = "Data@LFC, Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@CAA, Data@Mort"
  Lc <- trlnorm(reps, Data@LFC[x], 0.2)
  Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  Mdb <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])
  MuC <- Data@Cat[x, length(Data@Cat[x, ])]
  Cc <- trlnorm(reps, MuC, Data@CV_Cat[x])
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
  TAC <- Ac * (0.6 * Kc)/(0.67 - (Lc/Linfc))
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}  # end of BK_CC
class(BK_CC) <- "MP"

#' @describeIn BK Uses an approximation to FMSY based on length at first capture and an
#' estimate of current abundance based on a mean-length estimator.
#' @export BK_ML
BK_ML <- function(x, Data, reps = 100) {
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
  FMSY <- (0.6 * Kc)/(0.67 - (Lc/Linfc))  # robustifying for use in MSETAC<-Ac*FMSY
  TAC <- Ac * FMSY
  TAC <- TAC[TAC > 0 & TAC < (mean(TAC, na.rm = T) + 3 * stats::sd(TAC, na.rm = T))][1:reps]
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(BK_ML) <- "MP"


#' Constant catch management procedure of Geromont and Butterworth (2014)
#' 
#' The TAC is the average catch over the last yrsmth (by default, 5) years. 
#' This is one of four constant catch rules of Geromont and Butterworth 2014.
#' 
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrsmth Years over which to calculate mean catches
#' @param xx Parameter controlling the TAC. Mean catches are multiplied by
#' (1-xx)
#' @return A Rec object containing a numeric vector of TAC recommendations
#' @author T. Carruthers
#' @references Geromont, H.F., Butterworth, D.S. 2014. Generic management
#' procedures for data-poor fisheries; forecasting with few data. ICES J. Mar.
#' Sci. doi:10.1093/icesjms/fst232
#' @describeIn CC1 The TAC is the average catch over last \code{yrsmth} years.
#' @export
CC1 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0) {
  dependencies = "Data@Cat, Data@CV_Cat"
  C_dat <- Data@Cat[x, (length(Data@Year) - (yrsmth - 1)):length(Data@Year)]
  TAC <- (1 - xx) * trlnorm(reps, mean(C_dat), Data@CV_Cat/(yrsmth^0.5))  # mean catches over the interval
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(CC1) <- "MP"



#' @describeIn CC1 An additional 30\% reduction in catch is taken compared to \code{CC1}.
#' @export 
CC4 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0.3) {
  dependencies = "Data@Cat, Data@CV_Cat"
  C_dat <- Data@Cat[x, (length(Data@Year) - (yrsmth - 1)):length(Data@Year)]
  TAC <- (1 - xx) * trlnorm(reps, mean(C_dat), Data@CV_Cat/(yrsmth^0.5))  # mean catches over the interval
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(CC4) <- "MP"


#' Age-composition-based estimate of current stock depletion given constant Z
#' linked to an FMSY estimate to provide OFL
#' 
#' Estimates an OFL based on a Stock Reduction analysis fitted to current
#' age-composition data. Knife-edge vulnerability at age at maturity allows for
#' an FMSY estimate. 
#' 
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC.
#' @note Given a fixed historical F, What level of depletion gives you this
#' length composition?
#' @return A Rec object with a vector of TAC recommendation
#' @describeIn CompSRA Basic control rule where OFL=FMSY*F/C
#' @author T. Carruthers
#' @export CompSRA
CompSRA <- function(x, Data, reps = 100) {
  # optimize for fixed F to get you to current depletion C/Fcur =
  # abundance
  dependencies = "Data@Mort, Data@CV_Mort, Data@vbK, Data@CV_vbK, Data@vbLinf, Data@CV_vbLinf, Data@vbt0, Data@CV_vbt0, Data@MaxAge, Data@wla, Data@CV_wla, Data@wlb, Data@CV_wlb, Data@L50, Data@CV_L50, Data@CAA, Data@steep, Data@CV_steep, Data@LFS, Data@CV_LFS, Data@LFC, Data@CV_LFC, Data@Cat"
  maxage <- Data@MaxAge
  TAC <- rep(NA, reps)
  for (i in 1:reps) {
    Mc <- trlnorm(1, Data@Mort[x], Data@CV_Mort)
    # hc <- trlnorm(1, Data@steep[x], Data@CV_steep[x])
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
    AMc <- trlnorm(1, iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], 
                          Data@L50[x]), Data@CV_L50[x])
    ac <- trlnorm(1, Data@wla[x], Data@CV_wla[x])
    bc <- trlnorm(1, Data@wlb[x], Data@CV_wlb[x])
    Catch <- Data@Cat[x, ]
    ny <- length(Catch)
    nyCAA <- dim(Data@CAA)[2]
    CAA <- Data@CAA[x, max(nyCAA - 2, 1):nyCAA, ]  # takes last two years as the sample (or last year if there is only one
    
    Nac <- exp(-Mc * ((1:maxage) - 1))  # put a rough range on estimate of R0 assuming a mean harvest rate of 10%
    Lac <- Linfc * (1 - exp(-Kc * ((1:maxage) - t0c)))
    Wac <- ac * Lac^bc
    AFC <- log(1 - min(0.99, LFCc/Linfc))/-Kc + t0c
    AFS <- log(1 - min(0.99, LFSc/Linfc))/-Kc + t0c
    if (AFC >= 0.7 * maxage) 
      AFC <- 0.7 * maxage
    if (AFS >= 0.9 * maxage) 
      AFS <- 0.9 * maxage
    KES <- max(2, ceiling(mean(c(AFC, AFS))))
    pred <- Nac * Wac
    pred[1:(KES - 1)] <- 0
    pred <- pred/sum(pred)
    pred <- ((mean(Catch)/0.1) * pred/Wac)/exp(-(1:maxage) * Mc)
    pred <- pred[pred > 0]
    R0range <- c(mean(pred)/1000, mean(pred) * 1000)
    
    fit <- optimize(SRAfunc, log(R0range), Mc, hc, maxage, LFSc, LFCc, 
                    Linfc, Kc, t0c, AMc, ac, bc, Catch, CAA)
    Ac <- SRAfunc(fit$minimum, Mc, hc, maxage, LFSc, LFCc, Linfc, Kc, 
                  t0c, AMc, ac, bc, Catch, CAA, opt = 2)
    fit2 <- optimize(SRAFMSY, log(c(1e-04, 3)), Mc, hc, maxage, LFSc, 
                     LFCc, Linfc, Kc, t0c, AMc, ac, bc)
    # FMSY<-SRAFMSY(fit2$minimum,Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc,opt=F)
    FMSY <- exp(fit2$minimum)
    if ((FMSY/Mc) > 3) 
      FMSY <- 3 * Mc
    TAC[i] <- Ac * FMSY
    # message(i, ' of ', reps) flush.console()
  }
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(CompSRA) <- "MP"




#' @describeIn CompSRA With a 40-10 control rule
#' @export CompSRA4010
CompSRA4010 <- function(x, Data, reps = 100) {
  # optimize for fixed F to get you to current depletion C/Fcur =
  # abundance for (x in 1:nsim) {
  dependencies = "Data@Mort, Data@CV_Mort, Data@vbK, Data@CV_vbK, Data@vbLinf, Data@CV_vbLinf, Data@vbt0, Data@CV_vbt0, Data@MaxAge, Data@wla, Data@CV_wla, Data@wlb, Data@CV_wlb, Data@L50, Data@CV_L50, Data@CAA, Data@steep, Data@CV_steep, Data@LFS, Data@CV_LFS, Data@LFC, Data@CV_LFC, Data@Cat"
  maxage <- Data@MaxAge
  TAC <- Bt_K <- rep(NA, reps)
  for (i in 1:reps) {
    Mc <- trlnorm(1, Data@Mort[x], Data@CV_Mort)
    # hc <- trlnorm(1, Data@steep[x], Data@CV_steep[x])
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
    AMc <- trlnorm(1, iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], 
                          Data@L50[x]), Data@CV_L50[x])
    ac <- trlnorm(1, Data@wla[x], Data@CV_wla[x])
    bc <- trlnorm(1, Data@wlb[x], Data@CV_wlb[x])
    Catch <- Data@Cat[x, ]
    ny <- length(Catch)
    nyCAA <- dim(Data@CAA)[2]
    CAA <- Data@CAA[x, max(nyCAA - 2, 1):nyCAA, ]  # takes last two years as the sample (or last year if there is only one
    
    Nac <- exp(-Mc * ((1:maxage) - 1))  # put a rough range on estimate of R0 assuming a mean harvest rate of 10%
    Lac <- Linfc * (1 - exp(-Kc * ((1:maxage) - t0c)))
    Wac <- ac * Lac^bc
    
    AFC <- log(1 - min(0.99, LFCc/Linfc))/-Kc + t0c
    AFS <- log(1 - min(0.99, LFSc/Linfc))/-Kc + t0c
    if (AFC >= 0.7 * maxage) 
      AFC <- 0.7 * maxage
    if (AFS >= 0.9 * maxage) 
      AFS <- 0.9 * maxage
    
    KES <- max(2, ceiling(mean(c(AFC, AFS))))
    pred <- Nac * Wac
    pred[1:(KES - 1)] <- 0
    pred <- pred/sum(pred)
    pred <- ((mean(Catch)/0.1) * pred/Wac)/exp(-(1:maxage) * Mc)
    pred <- pred[pred > 0]
    R0range <- c(mean(pred)/1000, mean(pred) * 1000)
    
    fit <- optimize(SRAfunc, log(R0range), Mc, hc, maxage, LFSc, LFCc, 
                    Linfc, Kc, t0c, AMc, ac, bc, Catch, CAA)
    Ac <- SRAfunc(fit$minimum, Mc, hc, maxage, LFSc, LFCc, Linfc, Kc, 
                  t0c, AMc, ac, bc, Catch, CAA, opt = 2)
    Bt_K[i] <- SRAfunc(fit$minimum, Mc, hc, maxage, LFSc, LFCc, Linfc, 
                       Kc, t0c, AMc, ac, bc, Catch, CAA, opt = 3)
    fit2 <- optimize(SRAFMSY, log(c(1e-04, 3)), Mc, hc, maxage, LFSc, 
                     LFCc, Linfc, Kc, t0c, AMc, ac, bc)
    FMSY <- SRAFMSY(fit2$minimum, Mc, hc, maxage, LFSc, LFCc, Linfc, 
                    Kc, t0c, AMc, ac, bc, opt = F)
    if ((FMSY/Mc) > 3) 
      FMSY <- 3 * Mc
    TAC[i] <- Ac * FMSY
  }
  # 40-10 rule
  cond1 <- Bt_K < 0.4 & Bt_K > 0.1
  cond2 <- Bt_K < 0.1
  TAC[cond1] <- TAC[cond1] * (Bt_K[cond1] - 0.1)/0.3
  TAC[cond2] <- TAC[cond2] * tiny  # this has to still be stochastic albeit very small
  
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
  
  # message(x, ' of ', nsim) flush.console() }
}
class(CompSRA4010) <- "MP"


#' @describeIn DCAC Depletion Adjusted Average Catch: essentially DCAC multiplied 
#' by 2*depletion and divided by BMSY/B0 (Bpeak) (Harford and Carruthers, 2017).
#' @export DAAC
DAAC <- function(x, Data, reps = 100) {
  # extended depletion-corrected average catch (Harford and Carruthers
  # 2017)
  dependencies = "Data@AvC, Data@t, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@Dt, Data@CV_Dt, Data@BMSY_B0, Data@CV_BMSY_B0"
  C_tot <- Data@AvC[x] * Data@t[x]
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  FMSY_M <- trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])
  Bt_K <- trlnorm(reps, Data@Dt[x], Data@CV_Dt[x])
  if (any(is.na(c(Data@BMSY_B0[x], Data@CV_BMSY_B0[x])))) 
    return(NA)
  BMSY_K <- rbeta(reps, alphaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
                                    Data@CV_BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
                                                                    Data@CV_BMSY_B0[x]))
  dcac <- C_tot/(Data@t[x] + ((1 - Bt_K)/(BMSY_K * FMSY_M * Mdb)))
  TAC <- dcac * Bt_K/BMSY_K
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(DAAC) <- "MP"

#' Depletion-Based Stock Reduction Analysis
#' 
#' User prescribed BMSY/B0, M, FMSY/M are used to find B0 and therefore the OFL
#' by back-constructing the stock to match a user specified level of stock
#' depletion (OFL = M * FMSY/M * depletion* B0).
#' 
#' 
#' @param x A position in a data-limited methods object.
#' @param Data A data-limited methods object.
#' @param reps The number of samples of the TAC (OFL) recommendation.
#' @return A Rec object with a vector of TAC (OFL) values.
#' @note This is set up to return the OFL (FMSY * current biomass).
#' 
#' You may have noticed that you -the user- specify three of the factors that
#' make the quota recommendation. So this can be quite a subjective method.
#' 
#' Also the DB-SRA method of this package isn't exactly the same as the
#' original method of Dick and MacCall (2011) because it has to work for
#' simulated depletions above BMSY/B0 and even on occasion over B0. Also it
#' doesn't have the modification for flatfish life histories that has
#' previously been applied by Dick and MacCall.
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
#' @describeIn DBSRA Base version. You specify a range of stock depletion and, 
#' given historical catches DB-SRA calculates what unfished biomass must have 
#' been to get you here given samples for M, FMSY relative to M and also BMSY 
#' relative to Bunfished.
#' @export 
DBSRA <- function(x, Data, reps = 100) {
  # returns a vector of DBSRA estimates of the TAC for a particular
  # simulation x for(x in 1:nsim){
  dependencies = "Data@Cat, Data@Dep, Data@CV_Dep, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M,Data@BMSY_B0, Data@CV_BMSY_B0, Data@L50"
  C_hist <- Data@Cat[x, ]
  TAC <- rep(NA, reps)
  DBSRAcount <- 1
  if (is.na(Data@Dep[x]) | is.na(Data@CV_Dep[x])) return(NA)
  while (DBSRAcount < (reps + 1)) {
    depo <- max(0.01, min(0.99, Data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
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
    
    adelay <- max(floor(iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],  Data@L50[x])), 1)
 
    # opt <- optimize(DBSRAopt, log(c(0.01 * mean(C_hist), 1000 * mean(C_hist))), C_hist = C_hist, 
                    # nys = length(C_hist), Mdb = Mdb, FMSY_M = FMSY_M, BMSY_K = BMSY_K, 
                    # Bt_K = Bt_K, adelay = adelay, tol = 0.01)
    # scale catches for optimization
    scaler <- 1000/mean(C_hist)
    C_hist2 <- scaler * C_hist
    opt <- optimize(DBSRAopt, log(c(0.01 * mean(C_hist2), 1000 * mean(C_hist2))), C_hist = C_hist2, 
                    nys = length(C_hist2), Mdb = Mdb, FMSY_M = FMSY_M, BMSY_K = BMSY_K, 
                    Bt_K = Bt_K, adelay = adelay, tol = 0.01)
    
    # if(opt$objective<0.1){
    Kc <- exp(opt$minimum) / scaler
    BMSYc <- Kc * BMSY_K
    FMSYc <- Mdb * FMSY_M
    UMSYc <- (FMSYc/(FMSYc + Mdb)) * (1 - exp(-(FMSYc + Mdb)))
    MSYc <- Kc * BMSY_K * UMSYc
    TAC[DBSRAcount] <- UMSYc * Kc * Bt_K
    DBSRAcount <- DBSRAcount + 1
    # }
    
    # print(DBSRAcount)
  }  # end of reps
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
  
}  # end of DBSRA_apply
class(DBSRA) <- "MP"


#' @describeIn DBSRA Assumes 40 percent current depletion (Bcurrent/B0 = 0.4), which is 
#' more or less the most optimistic state for a stock (ie very close to BMSY/B0 for many stocks).
#' @export 
DBSRA_40 <- function(x, Data, reps = 100) {
  # returns a vector of DBSRA estimates of the TAC for a particular
  # simulation x
  dependencies = "Data@Cat, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@BMSY_B0, Data@CV_BMSY_B0, Data@L50"
  C_hist <- Data@Cat[x, ]
  TAC <- rep(NA, reps)
  DBSRAcount <- 1
  if (is.na(Data@Dep[x]) | is.na(Data@CV_Dep[x]))   return(NA)
  while (DBSRAcount < (reps + 1)) {
    depo <- 0.4
    Bt_K <- rbeta(100, alphaconv(depo, min(depo * Data@CV_Dep[x], 
                                           (1 - depo) * Data@CV_Dep[x])), betaconv(depo, min(depo * 
                                                                                               Data@CV_Dep[x], (1 - depo) * Data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
    Bt_K <- Bt_K[Bt_K > 0.00999 & Bt_K < 0.99001][1]  # interval censor (0.01,0.99)  as in Dick and MacCall 2011
    Mdb <- stats::rlnorm(100, mconv(Data@Mort[x], Data@CV_Mort[x] * 
                                      Data@Mort[x]), sdconv(Data@Mort[x], Data@CV_Mort[x] * 
                                                              Data@Mort[x]))  # log space stdev 0.4 as in Dick and MacCall 2011
    Mdb <- Mdb[Mdb < 0.9][1]  # !!!! maximum M is 0.9   interval censor
    if (is.na(Mdb)) 
      Mdb <- 0.9  # !!!! maximum M is 0.9   absolute limit
    FMSY_M <- trlnorm(1, Data@FMSY_M[x], Data@CV_FMSY_M[x])
    BMSY_K <- rbeta(100, alphaconv(Data@BMSY_B0[x], Data@CV_BMSY_B0[x] * 
                                     Data@BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@CV_BMSY_B0[x] * 
                                                                  Data@BMSY_B0[x]))  #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
    tryBMSY_K <- BMSY_K[BMSY_K > 0.05 & BMSY_K < 0.95][1]  # interval censor (0.05,0.95) as in Dick and MacCall, 2011
    if (is.na(tryBMSY_K)) {
      Min <- min(BMSY_K, na.rm = TRUE)
      Max <- max(BMSY_K, na.rm = TRUE)
      if (Max <= 0.05) 
        BMSY_K <- 0.05
      if (Min >= 0.95) 
        BMSY_K <- 0.95
    }
    if (!is.na(tryBMSY_K))  BMSY_K <- tryBMSY_K
    adelay <- max(floor(iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],  Data@L50[x])), 1)
    # opt <- optimize(DBSRAopt, log(c(0.1 * mean(C_hist), 1000 * mean(C_hist))), 
    #                 C_hist = C_hist, nys = length(C_hist), Mdb = Mdb, FMSY_M = FMSY_M, 
    #                 BMSY_K = BMSY_K, Bt_K = Bt_K, adelay = adelay, tol = 0.01)
    # scale catches for optimization
    scaler <- 1000/mean(C_hist)
    C_hist2 <- scaler * C_hist
    opt <- optimize(DBSRAopt, log(c(0.01 * mean(C_hist2), 1000 * mean(C_hist2))), C_hist = C_hist2, 
                    nys = length(C_hist2), Mdb = Mdb, FMSY_M = FMSY_M, BMSY_K = BMSY_K, 
                    Bt_K = Bt_K, adelay = adelay, tol = 0.01)
    
    # if(opt$objective<0.1){
    Kc <- exp(opt$minimum)
    BMSYc <- Kc * BMSY_K
    FMSYc <- Mdb * FMSY_M
    UMSYc <- (FMSYc/(FMSYc + Mdb)) * (1 - exp(-(FMSYc + Mdb)))
    MSYc <- Kc * BMSY_K * UMSYc
    TAC[DBSRAcount] <- UMSYc * Kc * Bt_K
    DBSRAcount <- DBSRAcount + 1
    # }
  }  # end of reps
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
  
}  # end of DBSRA_apply
class(DBSRA_40) <- "MP"



#' @describeIn DBSRA Base version paired with the 40-10 rule that throttles
#' back the OFL to zero at 10 percent of unfished biomass.
#' @export DBSRA4010
DBSRA4010 <- function(x, Data, reps = 100) {
  # returns a vector of DBSRA estimates of the TAC for a particular
  # simulation x for(x in 1:nsim){
  dependencies = "Data@Cat, Data@Dep, Data@CV_Dep, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M,Data@BMSY_B0, Data@CV_BMSY_B0, Data@L50"
  C_hist <- Data@Cat[x, ]
  TAC <- rep(NA, reps)
  DBSRAcount <- 1
  if (is.na(Data@Dep[x]) | is.na(Data@CV_Dep[x])) 
    return(NA)
  while (DBSRAcount < (reps + 1)) {
    depo <- max(0.01, min(0.99, Data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
    Bt_K <- rbeta(100, alphaconv(depo, min(depo * Data@CV_Dep[x], 
                                           (1 - depo) * Data@CV_Dep[x])), betaconv(depo, min(depo * 
                                                                                               Data@CV_Dep[x], (1 - depo) * Data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
    Bt_K <- Bt_K[Bt_K > 0.00999 & Bt_K < 0.99001][1]  # interval censor (0.01,0.99)  as in Dick and MacCall 2011
    Mdb <- trlnorm(100, Data@Mort[x], Data@CV_Mort[x])
    Mdb <- Mdb[Mdb < 0.9][1]  # !!!! maximum M is 0.9   interval censor
    if (is.na(Mdb)) 
      Mdb <- 0.9  # !!!! maximum M is 0.9   absolute limit
    FMSY_M <- trlnorm(1, Data@FMSY_M[x], Data@CV_FMSY_M[x])
    BMSY_K <- rbeta(100, alphaconv(Data@BMSY_B0[x], Data@CV_BMSY_B0[x] * 
                                     Data@BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@CV_BMSY_B0[x] * 
                                                                  Data@BMSY_B0[x]))  #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
    tryBMSY_K <- BMSY_K[BMSY_K > 0.05 & BMSY_K < 0.95][1]  # interval censor (0.05,0.95) as in Dick and MacCall, 2011
    if (is.na(tryBMSY_K)) {
      Min <- min(BMSY_K, na.rm = TRUE)
      Max <- max(BMSY_K, na.rm = TRUE)
      if (Max <= 0.05) 
        BMSY_K <- 0.05
      if (Min >= 0.95) 
        BMSY_K <- 0.95
    }
    if (!is.na(tryBMSY_K))  BMSY_K <- tryBMSY_K
    adelay <- max(floor(iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], 
                            Data@L50[x])), 1)
    # opt <- optimize(DBSRAopt, log(c(0.01 * mean(C_hist), 1000 * mean(C_hist))), 
    #                 C_hist = C_hist, nys = length(C_hist), Mdb = Mdb, FMSY_M = FMSY_M, 
    #                 BMSY_K = BMSY_K, Bt_K = Bt_K, adelay = adelay, tol = 0.01)
    # scale catches for optimization
    scaler <- 1000/mean(C_hist)
    C_hist2 <- scaler * C_hist
    opt <- optimize(DBSRAopt, log(c(0.01 * mean(C_hist2), 1000 * mean(C_hist2))), C_hist = C_hist2, 
                    nys = length(C_hist2), Mdb = Mdb, FMSY_M = FMSY_M, BMSY_K = BMSY_K, 
                    Bt_K = Bt_K, adelay = adelay, tol = 0.01)
    
    # if(opt$objective<0.1){
    Kc <- exp(opt$minimum)
    BMSYc <- Kc * BMSY_K
    FMSYc <- Mdb * FMSY_M
    UMSYc <- (FMSYc/(FMSYc + Mdb)) * (1 - exp(-(FMSYc + Mdb)))
    MSYc <- Kc * BMSY_K * UMSYc
    TAC[DBSRAcount] <- UMSYc * Kc * Bt_K
    # 40-10 rule
    if (Bt_K < 0.4 & Bt_K > 0.1) 
      TAC[DBSRAcount] <- TAC[DBSRAcount] * (Bt_K - 0.1)/0.3
    if (Bt_K < 0.1) 
      TAC[DBSRAcount] <- TAC[DBSRAcount] * tiny  # this has to still be stochastic albeit very small
    DBSRAcount <- DBSRAcount + 1
    # }
  }  # end of reps
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
  
  # }
}  # end of DBSRA_apply
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

#' Depletion Corrected Average Catch
#' 
#' A method of calculating an MSY proxy (FMSY * BMSY and therefore the OFL at
#' most productive stock size) based on average catches accounting for the
#' windfall catch that got the stock down to BMSY levels.
#' 
#' 
#' @usage DCAC(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @note It's probably worth noting that DCAC TAC recommendations do not tend
#' to zero as depletion tends to zero. It adjusts for depletion only in
#' calculating historical average catch. It follows that at stock levels much
#' below BMSY, DCAC tends to chronically overfish.
#' @author T. Carruthers
#' @references 
#' MacCall, A.D., 2009. Depletion-corrected average catch: a simple
#' formula for estimating sustainable yields in data-poor situations. ICES J.
#' Mar. Sci. 66, 2267-2271.
#' 
#' Harford W. and Carruthers, T. 2017. Interim and long-term performance of 
#' static and adaptive management procedures. Fish. Res. 190, 84-94.
#' @export DCAC
DCAC <- function(x, Data, reps = 100) {
  dependencies = "Data@AvC, Data@t, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@Dt, Data@CV_Dt, Data@BMSY_B0, Data@CV_BMSY_B0"
  C_tot <- Data@AvC[x] * Data@t[x]
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
  FMSY_M <- trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])  # standard deviation of 0.2 - referred to as 'standard error' in MacCall 2009
  Bt_K <- trlnorm(reps, Data@Dt[x], Data@CV_Dt[x])
  Bt_K[Bt_K>1] <-1
  if (any(is.na(c(Data@BMSY_B0[x], Data@CV_BMSY_B0[x])))) return(NA)
  BMSY_K <- rbeta(reps, alphaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
                                    Data@CV_BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
                                                                    Data@CV_BMSY_B0[x]))  #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(C_tot/(Data@t[x] + ((1 - Bt_K)/(BMSY_K * FMSY_M * Mdb))))
  Rec
}  # end of DCAC
class(DCAC) <- "MP"


#' @describeIn DCAC In this method, DCAC is paired with the 40-10 rule that throttles 
#' back the OFL to zero at 10 percent of unfished stock size (the OFL is not subject to downward
#' adjustment above 40 percent unfished). DCAC can overfish below BMSY levels. The 40-10 
#' harvest control rule largely resolves this problem providing an MP with surprisingly good
#' performance even at low stock levels.
#' @export DCAC4010
DCAC4010 <- function(x, Data, reps = 100) {
  dependencies = "Data@AvC, Data@t, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@Dt, Data@CV_Dt, Data@BMSY_B0, Data@CV_BMSY_B0"
  C_tot <- Data@AvC[x] * Data@t[x]
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
  FMSY_M <- trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])  # standard deviation of 0.2 - referred to as 'standard error' in MacCall 2009
  Bt_K <- trlnorm(reps, Data@Dt[x], Data@CV_Dt[x])
  Bt_K[Bt_K>1] <-1
  if (any(is.na(c(Data@BMSY_B0[x], Data@CV_BMSY_B0[x])))) 
    return(NA)
  BMSY_K <- rbeta(reps, alphaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
                                    Data@CV_BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
                                                                    Data@CV_BMSY_B0[x]))  #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
  TAC <- C_tot/(Data@t[x] + ((1 - Bt_K)/(BMSY_K * FMSY_M * Mdb)))
  # 40-10 rule
  cond1 <- Bt_K < 0.4 & Bt_K > 0.1
  cond2 <- Bt_K < 0.1
  if (length(cond1) > 0) 
    TAC[cond1] <- TAC[cond1] * (Bt_K[cond1] - 0.1)/0.3
  if (length(cond2) > 0) 
    TAC[cond2] <- TAC[cond2] * tiny  # this has to still be stochastic albeit very small
  if (length(cond1) < 1 & length(cond2) < 1) 
    return(NA)
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
  
}  # end of DCAC
class(DCAC4010) <- "MP"



#' @describeIn DCAC This variant assumes that current stock biomass is exactly 
#' 40 per cent of unfished levels. The 40 percent depletion assumption may not 
#' really affect DCAC that much as it already makes TAC recommendations that are 
#' quite MSY-like.
#' @export DCAC_40
DCAC_40 <- function(x, Data, reps = 100) {
  dependencies = "Data@AvC, Data@t, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@BMSY_B0, Data@CV_BMSY_B0"
  C_tot <- Data@AvC[x] * Data@t[x]
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  FMSY_M <- trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])
  Bt_K <- 0.4
  if (any(is.na(c(Data@BMSY_B0[x], Data@CV_BMSY_B0[x])))) 
    return(NA)
  BMSY_K <- rbeta(reps, alphaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
                                    Data@CV_BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
                                                                    Data@CV_BMSY_B0[x]))  #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
  TAC <- (C_tot/(Data@t[x] + ((1 - Bt_K)/(BMSY_K * FMSY_M * Mdb))))
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}  # end of DCAC40
class(DCAC_40) <- "MP"



#' @describeIn DCAC This variant uses the mean length estimator to calculate current stock
#' depletion. The mean length extension was programmed by Gary Nelson as part of his
#' excellent R package 'fishmethods'.
#' @export DCAC_ML
DCAC_ML <- function(x, Data, reps = 100) {
  dependencies = "Data@AvC, Data@t, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@BMSY_B0, Data@CV_BMSY_B0, Data@Year, Data@CAL, Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK"
  if (is.na(Data@BMSY_B0[x]) | is.na(Data@CV_BMSY_B0[x])) return(NA)
  if (is.na(Data@FMSY_M[x]) | is.na(Data@CV_FMSY_M[x])) return(NA)
  C_tot <- Data@AvC[x] * Data@t[x]
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # default CV of 0.5 as in MacCall 2009
  FMSY_M <- trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])  # standard deviation of 0.2 - referred to as 'standard error' in MacCall 2009
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
  Bt_K[Bt_K>1] <-1
  if (any(is.na(c(Data@BMSY_B0[x], Data@CV_BMSY_B0[x])))) return(NA)
  BMSY_K <- rbeta(reps, alphaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
                                    Data@CV_BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
                                                                    Data@CV_BMSY_B0[x]))  #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
  TAC <- C_tot/(Data@t[x] + ((1 - Bt_K)/(BMSY_K * FMSY_M * Mdb)))
  
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}  # end of DCAC_ML
class(DCAC_ML) <- "MP"


#' Delay - Difference Stock Assessment with UMSY and MSY as leading parameters
#' 
#' A simple delay-difference assessment that estimates the TAC using a
#' time-series of catches and a relative abundance index. Conditioned on effort.
#' 
#' 
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @param LB The lowest permitted factor of previous fishing effort
#' @param UB The highest permitted factor of previous fishing effort
#' @return A Rec object of either TAC or effort recommendations
#' @note This DD model is observation error only and has does not estimate
#' process error (recruitment deviations). Similar to many other assessment
#' models it depends on a whole host of dubious assumptions such as temporally
#' stationary productivity and proportionality between the abundance index and
#' real abundance. Unsurprisingly the extent to which these assumptions are
#' violated tends to be the biggest driver of performance for this method.
#' @author T. Carruthers
#' @references  
#' Carruthers, T, Walters, C.J,, and McAllister, M.K. 2012. Evaluating methods that classify
#' fisheries stock status using only fisheries catch data. Fisheries Research 119-120:66-79.
#' 
#' Hilborn, R., and Walters, C. 1992. Quantitative Fisheries Stock Assessment: Choice,
#' Dynamics and Uncertainty. Chapman and Hall, New York. 
#' @describeIn DD Base version where the TAC = UMSY * Current Biomass.
#' @export 
DD <- function(x, Data, reps = 100) {
  dependencies = "Data@vbLinf, Data@vbK, Data@vbt0, Data@Mort, Data@wla, Data@wlb, Data@Cat, Data@Ind, Data@L50, Data@MaxAge"
  Winf = Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
  age <- 1:Data@MaxAge
  la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
  wa <- Data@wla[x] * la^Data@wlb[x]
  a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],  Data@L50[x])
  a50V <- max(a50V, 1)
  yind <- (1:length(Data@Cat[x, ]))[!is.na(Data@Cat[x, ] + Data@Ind[x,   ])]
  C_hist <- Data@Cat[x, yind]
  E_hist <- C_hist/Data@Ind[x, yind]
  E_hist <- E_hist/mean(E_hist)
  ny_DD <- length(C_hist)
  k_DD <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)  
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
               C_hist = C_hist, UMSYprior = UMSYprior, method = "BFGS", hessian = TRUE)
  
  TAC <- rep(NA, reps)
  # samps<-rmvnorm(reps,opt$par,solve(opt$hessian)) # assuming log
  # parameters are multivariate normal hessian approximation
  samps <- cbind(rnorm(reps, opt$par[1], ((opt$par[1])^2)^0.5 * 0.1), 
                 rnorm(reps, opt$par[2], ((opt$par[2])^2)^0.5 * 0.1), rnorm(reps, 
                                                                            opt$par[3], ((opt$par[3])^2)^0.5 * 0.1))
  if (reps == 1) 
    samps <- matrix(c(opt$par[1], opt$par[2], opt$par[3]), nrow = 1)
  for (i in 1:reps) TAC[i] <- DD_R(samps[i, ], opty = 2, So_DD = So_DD, 
                                   Alpha_DD = Alpha_DD, Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, 
                                   wa_DD = wa_DD, E_hist = E_hist, C_hist = C_hist, UMSYprior = UMSYprior)
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(DD) <- "MP"


#' @describeIn DD In this version, a 40-10 rule is imposed over the TAC recommendation.
#' @export DD4010
DD4010 <- function(x, Data, reps = 100) {
  dependencies = "Data@vbLinf, Data@vbK, Data@vbt0, Data@Mort, Data@wla, Data@wlb, Data@Cat, Data@Ind, Data@L50, Data@MaxAge"
  Winf = Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
  age <- 1:Data@MaxAge
  la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
  wa <- Data@wla[x] * la^Data@wlb[x]
  a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], 
              Data@L50[x])
  a50V <- max(a50V, 1)
  yind <- (1:length(Data@Cat[x, ]))[!is.na(Data@Cat[x, ] + Data@Ind[x, 
                                                                    ])]
  C_hist <- Data@Cat[x, yind]
  E_hist <- Data@Ind[x, yind]
  E_hist <- C_hist/E_hist
  E_hist <- E_hist/mean(E_hist)
  ny_DD <- length(C_hist)
  k_DD <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)
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
               C_hist = C_hist, UMSYprior = UMSYprior, method = "BFGS", hessian = TRUE)
  
  TAC <- rep(NA, reps)
  dep <- rep(NA, reps)
  # samps<-rmvnorm(reps,opt$par,solve(opt$hessian)) # assuming log
  # parameters are multivariate normal hessian approximation
  samps <- cbind(rnorm(reps, opt$par[1], ((opt$par[1])^2)^0.5 * 0.1), 
                 rnorm(reps, opt$par[2], ((opt$par[2])^2)^0.5 * 0.1), rnorm(reps, 
                                                                            opt$par[3], ((opt$par[3])^2)^0.5 * 0.1))
  if (reps == 1) 
    samps <- matrix(c(opt$par[1], opt$par[2], opt$par[3]), nrow = 1)
  for (i in 1:reps) TAC[i] <- DD_R(samps[i, ], opty = 2, So_DD = So_DD, 
                                   Alpha_DD = Alpha_DD, Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, 
                                   wa_DD = wa_DD, E_hist = E_hist, C_hist = C_hist, UMSYprior = UMSYprior)
  for (i in 1:reps) dep[i] <- DD_R(samps[i, ], opty = 3, So_DD = So_DD, 
                                   Alpha_DD = Alpha_DD, Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, 
                                   wa_DD = wa_DD, E_hist = E_hist, C_hist = C_hist, UMSYprior = UMSYprior)
  cond1 <- !is.na(dep) & dep < 0.4 & dep > 0.1
  cond2 <- !is.na(dep) & dep < 0.1
  TAC[cond1] <- TAC[cond1] * (dep[cond1] - 0.1)/0.3
  TAC[cond2] <- TAC[cond2] * tiny  # this has to still be stochastic albeit very small
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(DD4010) <- "MP"


#' @describeIn Fratio Depletion Corrected Fratio: the Fratio MP with a harvest control 
#' rule that reduces F according to the production curve given an estimate of current 
#' stock depletion (made-up for this package).
#' @export DepF
DepF <- function(x, Data, reps = 100) {
  dependencies = "Data@Year, Data@Dep, Data@Mort, Data@FMSY_M, Data@BMSY_B0"
  Frat <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x]) * trlnorm(reps, 
                                                                 Data@FMSY_M[x], Data@CV_FMSY_M[x])
  if (is.na(Data@Dep[x]) | is.na(Data@CV_Dep[x])) 
    return(NA)
  depo <- max(0.01, min(0.99, Data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
  Bt_K <- rbeta(reps * 100, alphaconv(depo, min(depo * Data@CV_Dep[x], 
                                                (1 - depo) * Data@CV_Dep[x])), betaconv(depo, min(depo * Data@CV_Dep[x], 
                                                                                                  (1 - depo) * Data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
  Bt_K <- Bt_K[Bt_K >= 0.01 & Bt_K <= 0.99][1:reps]  # interval censor (0.01,0.99)  as in Dick and MacCall 2011
  adj <- Bt_K * (1 - Bt_K) * 4
  adj[Bt_K > 0.5] <- 1
  TAC <- Frat * Data@Abun[x] * adj
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(DepF) <- "MP"

#' Dynamic Fratio MP
#' 
#' The Fratio MP with a controller that changes the level of F according to the
#' relationship between Surplus production and biomass. Ie lower F when dSP/dB
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
#' @usage DynF(x, Data, yrsmth=10, gg=2, reps = 100)
#' @param x A position in a data-limited methods object
#' @param Data A data-limited methods object
#' @param yrsmth The number of historical recent years used for smoothing catch
#' and biomass data
#' @param gg A gain parameter that modifies F according to the gradient in
#' surplus production with biomass
#' @param reps The number samples of the TAC
#' @return A numeric vector of TAC recommendations
#' @author T. Carruthers
#' @references Made-up for this package.
#' @seealso \link{Fratio}
#' @export DynF
DynF <- function(x, Data, yrsmth = 10, gg = 2, reps = 100) {
  
  dependencies = "Data@Year, Data@Cat, Data@Ind, Data@Abun, Data@Mort, Data@FMSY_M"
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
  
  Frat <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x]) * trlnorm(reps, 
                                                                 Data@FMSY_M[x], Data@CV_FMSY_M[x])
  Flim <- matrix(NA, nrow = 2, ncol = reps)
  Flim[1, ] <- Frat * 0.5
  Flim[2, ] <- Frat * 2
  
  yind <- 1:length(SP_hist)
  SP_mu <- predict(lm(SP_hist ~ yind), newdat = list(yind = length(SP_hist) + 
                                                       1))
  SP_se <- predict(lm(SP_hist ~ yind), newdat = list(yind = length(SP_hist) + 
                                                       1), se = T)$se.fit
  SP_new <- rnorm(reps, SP_mu, SP_se/2)
  Glm <- summary(lm(SP_hist ~ B_hist[ind1]))$coefficients[2, 1:2]  # plot(B_hist[ind1],SP_hist) # points(B_hist[ind1],SP_hist,col='green')
  G_new <- rnorm(reps, Glm[1], Glm[2]/2)
  # G_new[G_new>2*Frat]<-2*Frat[G_new<(2*Frat)]
  # G_new[G_new<(-2*Frat)]<--2*Frat[G_new<(-2*Frat)]
  G_new[G_new > 0] <- G_new[G_new > 0] * 3
  newF <- Frat * exp(-G_new * gg)
  newF[newF < Flim[1]] <- Flim[1]
  newF[newF > Flim[2]] <- Flim[2]
  
  TAC <- newF * B_hist[yrsmth]
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
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
#' @usage Fadapt(x, Data, reps = 100, yrsmth = 7, gg=1)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
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
#' @export Fadapt
Fadapt <- function(x, Data, reps = 100, yrsmth = 7, gg = 1) {
  
  dependencies = "Data@Year, Data@Cat, Data@Ind, Data@Abun, Data@Mort, Data@FMSY_M"
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
  
  Frat <- Data@Mort[x] * Data@FMSY_M[x]
  Flim <- Frat * c(0.5, 2)
  Flimr <- Flim[2] - Flim[1]
  
  yind <- 1:length(SP_hist)
  SP_mu <- predict(lm(SP_hist ~ yind), newdat = list(yind = length(SP_hist) + 
                                                       1))
  SP_se <- predict(lm(SP_hist ~ yind), newdat = list(yind = length(SP_hist) + 
                                                       1), se = T)$se.fit
  SP_new <- rnorm(reps, SP_mu, SP_se/2)
  Glm <- summary(lm(SP_hist ~ B_hist[ind1]))$coefficients[2, 1:2]  # plot(B_hist[ind1],SP_hist) # points(B_hist[ind1],SP_hist,col='green')
  G_new <- rnorm(reps, Glm[1], Glm[2])
  
  Fold <- mean(C_hist/B_hist)
  
  if (Fold < Flim[1]) 
    Fmod1 <- (-2)
  if (Fold > Flim[2]) 
    Fmod1 <- 2
  if (Fold > Flim[1] & Fold < Flim[2]) {
    Ffrac <- (Fold - Flim[1])/Flimr
    Fmod1 <- log(Ffrac/(1 - Ffrac))
  }
  Fmod2 <- Fmod1 + gg * -G_new
  newF <- Flim[1] + (exp(Fmod2)/(1 + exp(Fmod2))) * Flimr
  TAC <- newF * B_hist[yrsmth]
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(Fadapt) <- "MP"


#' Demographic FMSY method
#' 
#' FMSY is calculated as r/2 where r is calculated from a demographic approach
#' (inc steepness). Coupled with an estimate of current abundance that gives
#' you the OFL.
#' 
#' Made up for this package. This uses Murdoch McAllister's demographic r
#' method to derive FMSY (r/2) and then makes the quota r*current biomass / 2.
#' Easy.
#' 
#' @usage Fdem(x, Data, reps = 100)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @author T. Carruthers
#' @references McAllister, M.K., Pikitch, E.K., and Babcock, E.A. 2001. Using
#' demographic methods to construct Bayesian priors for the intrinsic rate of
#' increase in the Schaefer model and implications for stock rebuilding. Can.
#' J. Fish. Aquat. Sci. 58: 1871-1890.
#' @export Fdem
Fdem <- function(x, Data, reps = 100) {
  # Demographic FMSY estimate (FMSY=r/2)
  dependencies = "Data@Mort, Data@CV_Mort, Data@vbK, Data@CV_vbK, Data@vbLinf, Data@CV_vbLinf, Data@vbt0, Data@CV_vbt0, Data@MaxAge, Data@wla, Data@wlb, Data@Abun, Data@CV_Abun, Data@steep, Data@CV_steep"
  Mvec <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  Linfc = trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
    t0c <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
  } else {
    t0c <- rep(Data@vbt0[x], reps)
  }
  t0c[!is.finite(t0c)] <- 0
  # hvec <- trlnorm(reps, Data@steep[x], Data@CV_steep[x])
  hvec <- sample_steepness2(reps, Data@steep[x], Data@CV_steep[x])
  
  Ac <- trlnorm(reps, Data@Abun[x], Data@CV_Abun[x])
  FMSY <- getr(x, Data, Mvec, Kc, Linfc, t0c, hvec, maxage = Data@MaxAge, 
               r_reps = reps)/2
  TAC <- FMSY * Ac
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(Fdem) <- "MP"



#' Demographic FMSY method using catch-curve analysis to estimate recent Z
#' 
#' FMSY is calculated as r/2 from a demographic r prior method, current
#' abudnance is estimated from naive catch curve analysis.
#' 
#' 
#' @usage Fdem_CC(x, Data, reps = 100, Fmin=0.005)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param Fmin The minimum fishing mortality rate derived from the catch-curve
#' analysis
#' @author T. Carruthers
#' @references McAllister, M.K., Pikitch, E.K., and Babcock, E.A. 2001. Using
#' demographic methods to construct Bayesian priors for the intrinsic rate of
#' increase in the Schaefer model and implications for stock rebuilding. Can.
#' J. Fish. Aquat. Sci. 58: 1871-1890.
#' @export Fdem_CC
Fdem_CC <- function(x, Data, reps = 100, Fmin = 0.005) {
  dependencies = "Data@Mort, Data@CV_Mort, Data@vbK, Data@CV_vbK, Data@vbLinf, Data@CV_vbLinf, Data@vbt0, Data@CV_vbt0, Data@MaxAge, Data@wla, Data@wlb, Data@CAA, Data@steep, Data@CV_steep"
  Mvec <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])
  Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  Linfc = trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
    t0c <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
  } else {
    t0c <- rep(Data@vbt0[x], reps)
  }
  t0c[!is.finite(t0c)] <- 0
  # hvec <- trlnorm(reps, Data@steep[x], Data@CV_steep[x])
  hvec <- sample_steepness2(reps, Data@steep[x], Data@CV_steep[x])
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
  FMSY <- getr(x, Data, Mvec, Kc, Linfc, t0c, hvec, maxage = Data@MaxAge, r_reps = reps)/2
  TAC <- FMSY * Ac
  
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(Fdem_CC) <- "MP"

#' Demographic FMSY method that uses mean length data to estimate recent Z
#' 
#' Demographic F (r/2) method using the mean length estimator to calculate
#' current abundance.
#' 
#' 
#' @usage Fdem_ML(x, Data, reps = 100)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @note The mean length extension was programmed by Gary Nelson as part of his
#' excellent R package 'fishmethods'
#' @author T. Carruthers
#' @references McAllister, M.K., Pikitch, E.K., and Babcock, E.A. 2001. Using
#' demographic methods to construct Bayesian priors for the intrinsic rate of
#' increase in the Schaefer model and implications for stock rebuilding. Can.
#' J. Fish. Aquat. Sci. 58: 1871-1890.
#' @export 
Fdem_ML <- function(x, Data, reps = 100) {
  dependencies = "Data@Mort, Data@CV_Mort, Data@vbK, Data@CV_vbK, Data@vbLinf, Data@CV_vbLinf, Data@vbt0, Data@CV_vbt0, Data@MaxAge, Data@wla, Data@wlb, Data@CAL, Data@steep, Data@CV_steep"
  Mvec <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])
  Kc <- trlnorm(reps * 10, Data@vbK[x], Data@CV_vbK[x])
  Linfc = trlnorm(reps * 10, Data@vbLinf[x], Data@CV_vbLinf[x])
  t0c <- -trlnorm(reps * 10, -Data@vbt0[x], Data@CV_vbt0[x])
  t0c[!is.finite(t0c)] <- 0
  # hvec <- trlnorm(reps * 10, Data@steep[x], Data@CV_steep[x])
  hvec <- sample_steepness2(reps*10, Data@steep[x], Data@CV_steep[x])
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
  Ac <- Cc[ind]/(1 - exp(-FM))
  FMSY <- getr(x, Data, Mvec[ind], Kc[ind], Linfc[ind], t0c[ind], hvec[ind], 
               maxage = Data@MaxAge, r_reps = sum(ind))/2
  TAC <- FMSY * Ac
  TAC <- TAC[TAC > 0][1:reps]
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(Fdem_ML) <- "MP"


#' An FMSY/M ratio method
#' 
#' Calculates the OFL based on a fixed ratio of FMSY to M multiplied by a
#' current estimate of abundance.
#'  
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of samples of the TAC recommendation
#' @param Fmin Minimum current fishing mortality rate for the catch-curve
#' analysis
#' @return A Rec object with TAC recommendations
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
#' @seealso \link{DynF}
#' @export Fratio
Fratio <- function(x, Data, reps = 100) {
  # FMSY / M ratio method e.g. Gulland
  depends = "Data@Abun,Data@CV_Abun,Data@FMSY_M, Data@CV_FMSY_M,Data@Mort,Data@CV_Mort"
  Ac <- trlnorm(reps, Data@Abun[x], Data@CV_Abun[x])
  TAC <- (Ac * trlnorm(reps, Data@Mort[x], Data@CV_Mort[x]) * 
            trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x]))
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}  # end of Fratio
class(Fratio) <- "MP"

#' @describeIn Fratio Paired with the 40-10
#' rule that throttles back the OFL to zero at 10 percent of unfished biomass.
#' @export Fratio4010
Fratio4010 <- function(x, Data, reps = 100) {
  # FMSY / M ratio method e.g. Gulland
  dependencies = "Data@Abun, Data@CV_Abun, Data@FMSY_M, Data@CV_FMSY_M, Data@Mort, Data@CV_Mort, Data@Dep"
  Ac <- trlnorm(reps, Data@Abun[x], Data@CV_Abun[x])
  TAC <- Ac * trlnorm(reps, Data@Mort[x], Data@CV_Mort[x]) * 
    trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])
  Bt_K <- trlnorm(reps, Data@Dt[x], Data@CV_Dt[x])
  # 40-10 rule
  cond1 <- Bt_K < 0.4 & Bt_K > 0.1
  cond2 <- Bt_K < 0.1
  TAC[cond1] <- TAC[cond1] * (Bt_K[cond1] - 0.1)/0.3
  TAC[cond2] <- TAC[cond2] * tiny  # this has to still be stochastic albeit very small
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}  # end of Fratio
class(Fratio4010) <- "MP"

#' @describeIn Fratio Pairs with a catch curve estimate of current stock size.
#' @export Fratio_CC
Fratio_CC <- function(x, Data, reps = 100, Fmin = 0.005) {
  # FMSY / M ratio method using catch curve analysis to determine current
  # abundance ================================== for (x in 1:nsim) {
  dependencies = " Data@FMSY_M, Data@CV_FMSY_M, Data@Mort, Data@CV_Mort, Data@Cat, Data@CV_Cat, Data@CAA"
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
  TAC <- Ac * Mdb * trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])
  
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
  
}  # end of Fratio_CC
class(Fratio_CC) <- "MP"

#' @describeIn Fratio Pairs with an estimate of
#' current stock size from a mean-length estimator.
#' @export Fratio_ML
Fratio_ML <- function(x, Data, reps = 100) {
  dependencies = " Data@FMSY_M, Data@CV_FMSY_M, Data@Mort, Data@CV_Mort, Data@Cat, Data@CV_Cat, Data@CAL"
  MuC <- Data@Cat[x, length(Data@Cat[x, ])]
  Cc <- trlnorm(reps * 10, MuC, Data@CV_Cat[x])
  Mdb <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
  Linfc <- trlnorm(reps * 10, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps * 10, Data@vbK[x], Data@CV_vbK[x])
  Z <- MLne(x, Data, Linfc = Linfc, Kc = Kc, ML_reps = reps * 10, MLtype = "F")
  if (all(is.na(Z))) {
    Rec <- new("Rec")
    Rec@TAC <- TACfilter(rep(NA, reps))
    return(Rec)
  } 
  FM <- Z - Mdb
  Ac <- Cc/(1 - exp(-FM))
  TAC <- Ac * trlnorm(reps * 10, Data@FMSY_M[x], Data@CV_FMSY_M[x]) * 
    trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  TAC <- TAC[TAC > 0][1:reps]
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(Fratio_ML) <- "MP"



#' Geromont and Butterworth Constant Catch Harvest Control Rule
#' 
#' A simple MP that aims for average historical catches (as a proxy for MSY)
#' subject to imperfect information.
#' 
#' Note that this is my interpretation of their MP and is now stochastic.
#' Currently it is generalized and is not 'tuned' to more detailed assessment
#' data which might explain why in some cases it leads to stock declines.
#' 
#' @usage GB_CC(x, Data, reps = 100)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @author T. Carruthers
#' @references Geromont, H.F. and Butterworth, D.S. 2014. Complex assessment or
#' simple management procedures for efficient fisheries management: a
#' comparative study. ICES J. Mar. Sci. doi:10.1093/icesjms/fsu017
#' @export GB_CC
GB_CC <- function(x, Data, reps = 100) {
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
#' @references Made-up for this package. Carruthers et al. 2015. Performance of
#' Simple Management Procedures.
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



#' @describeIn DCAC Hybrid Depletion Adjusted Average Catch: essentially DCAC 
#' multiplied by 2*depletion and divided by BMSY/B0 (Bpeak) when below BMSY, 
#' and DCAC above BMSY (Harford and Carruthers 2017).
#' @export 
HDAAC <- function(x, Data, reps = 100) {
  dependencies = "Data@AvC, Data@t, Data@Mort, Data@CV_Mort, Data@Dt, Data@CV_Dt, Data@BMSY_B0, Data@CV_BMSY_B0"
  C_tot <- Data@AvC[x] * Data@t[x]
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  FMSY_M <- trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])
  Bt_K <- trlnorm(reps, Data@Dt[x], Data@CV_Dt[x])
  if (any(is.na(c(Data@BMSY_B0[x], Data@CV_BMSY_B0[x])))) 
    return(NA)
  BMSY_K <- rbeta(reps, alphaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
                                    Data@CV_BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * Data@CV_BMSY_B0[x]))
  dcac <- C_tot/(Data@t[x] + ((1 - Bt_K)/(BMSY_K * FMSY_M * Mdb)))
  ddcac <- dcac * Bt_K/BMSY_K
  TAC <- dcac
  TAC[Bt_K < BMSY_K] <- ddcac[Bt_K < BMSY_K]
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(HDAAC) <- "MP"



#' Index Confidence Interval (ICI) MP by Jardim et al. (2015)
#' 
#' The MP adjusts catch based on the value of the index in the current year relative to the 
#' time series mean and standard error.
#'  
#' The mean and standard error of the index time series is calculated. There are two thresholds 
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



#' Less Precautionary Index Confidence Interval (ICI) MP by Jardim et al. (2015)
#' 
#' The MP adjusts catch based on the value of the index in the current year relative to the 
#' time series mean and standard error. This method is less precautionary of the two ICI MPs by allowing for a larger increase in TAC
#' and a lower threshold of the index to decrease the TAC (see Jardim et al. 2015).
#' 
#' The mean and standard error of the index time series is calculated. There are two thresholds 
#' which delineates whether catch is reduced, held constant, or increased. The catch is reduced by 0.75
#' if the Z-score of the current year's index is less than -1.96. The catch is increased by 1.25
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
#' being the mean index in the three years prior to those in the numerator.
#' 
#' This MP is the stochastic version of Method 3.2 used by ICES for Data-Limited Stocks (ICES 2012).
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
#' @references Carruthers et al. 2015. Performance review of simple management
#' procedures. Fish and Fisheries. In press.
#' 
#' Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for
#' data-poor fisheries; forecasting with few data. ICES J. Mar. Sci.
#' doi:10.1093/icesjms/fst232
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



#' A management procedure that incrementally adjusts the TAC to maintain a
#' constant CPUE or relative abundance index
#' 
#' The most biologically precautionary of two constant index / CPUE methods
#' proposed by Geromont and Butterworth 2014. Tested by Carruthers et al. 2015
#' 
#' Tested by Carruthers et al. 2015.
#' 
#' @usage Islope4(x, Data, reps = 100, yrsmth = 5, lambda=0.2,xx=0.4)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param lambda A gain parameter controlling the speed in update in TAC.
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @return A numeric vector of quota recommendations
#' @author T. Carruthers
#' @references Carruthers et al. 2015. Performance evaluation of simple
#' management procedures. Fish and Fisheries. In press.
#' 
#' Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for
#' data-poor fisheries; forecasting with few data. ICES J. Mar. Sci.
#' doi:10.1093/icesjms/fst232
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


#' Index Target 5
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

#' Index Target 10
#' 
#' An index target MP where the TAC is modified according to current index
#' levels (mean index over last 5 years) relative to a target level. Maximum
#' annual changes are 10 per cent.
#' 
#' 
#' @usage IT10(x, Data, reps = 100,yrsmth=5,mc=0.1)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the quota recommendation
#' @param yrsmth The number of historical years over which to average the index
#' @param mc The maximum fractional change in the TAC among years.
#' @return A numeric vector of TAC recommendations
#' @author T. Carruthers
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
#' reference level that is a fraction of mean recent catches) to reach a target
#' CPUE / relative abundance index
#' 
#' The least biologically precautionary of two index/CPUE target MPs proposed
#' by Geromont and Butterworth 2014. Tested by Carruthers et al. 2015
#' 
#' Tested by Carruthers et al. 2015.
#' 
#' @usage Itarget1(x, Data, reps = 100, yrsmth = 5, xx=0, Imulti=1.5)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param Imulti Parameter controlling how much larger target CPUE / index is
#' compared with recent levels.
#' @return A numeric vector of TAC recommendations
#' @author T. Carruthers
#' @references Carruthers et al. 2015. Performance evaluation of simple
#' management procedures. Fish and Fisheries. In press.
#' 
#' Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for
#' data-poor fisheries; forecasting with few data. ICES J. Mar. Sci.
#' doi:10.1093/icesjms/fst232
#' @export Itarget1
Itarget1 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0, Imulti = 1.5) {
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
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(Itarget1) <- "MP"



#' A management procedure that incrementally adjusts the TAC (starting from
#' reference level that is a fraction of mean recent catches) to reach a target
#' CPUE / relative abundance index
#' 
#' The most biologically precautionary of two index/CPUE target MPs proposed by
#' Geromont and Butterworth 2014. Tested by Carruthers et al. 2015
#' 
#' Tested by Carruthers et al. 2015.
#' 
#' @usage Itarget4(x, Data, reps = 100, yrsmth = 5, xx=0.3, Imulti=2.5)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param Imulti Parameter controlling how much larger target CPUE / index is
#' compared with recent levels.
#' @return A numeric vector of TAC recommendations
#' @author T. Carruthers
#' @references Carruthers et al. 2015. Performance evaluation of simple
#' management procedures. Fish and Fisheries. In press.
#' 
#' Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for
#' data-poor fisheries; forecasting with few data. ICES J. Mar. Sci.
#' doi:10.1093/icesjms/fst232
#' @export Itarget4
Itarget4 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0.3, Imulti = 2.5) {
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


#' Mean length-based indicator MP of Jardim et al. 2015 using Beverton-Holt invariant 
#' M/K ratio = 1.5 and assumes FMSY = M.
#' 
#' The TAC is adjusted by the ratio alpha, where the numerator 
#' is the mean length of the catch (of lengths larger than Lc) and 
#' the denominator is the mean length expected when FMSY = M and M/K = 1.5. 
#' Natural mortality M and von Bertalanffy K are not used in this MP 
#' (see Appendix A of Jardim et al. 2015). Here, Lc is the length at 
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


#' The more general version of the mean length-based indicator MP of Jardim et al. 2015.
#' 
#' The TAC is adjusted by the ratio alpha, where the numerator 
#' is the mean length of the catch (of lengths larger than Lc) 
#' and the denominator is the mean length as a function of Linf,
#' FMSY/M, and M/K (see Appendix A of Jardim et al. 2015). Here, 
#' Lc is the length at full selection (LFS).
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


#' A management procedure that incrementally adjusts the TAC according to the
#' mean length of recent catches.
#' 
#' The least biologically precautionary of four adaptive length-based MPs
#' proposed by Geromont and Butterworth 2014. Tested by Carruthers et al. 2015
#' 
#' Tested by Carruthers et al. 2015.
#' 
#' @usage LstepCC1(x, Data, reps = 100, yrsmth = 5, xx=0, stepsz=0.05,
#' llim=c(0.96,0.98,1.05))
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param stepsz Parameter controlling the size of the TAC update increment.
#' @param llim A vector of length reference points that determine the
#' conditions for increasing, maintaining or reducing the TAC.
#' @return A numeric vector of TAC recommendations
#' @author T. Carruthers
#' @references Carruthers et al. 2015. Performance evaluation of simple
#' management procedures. Fish and Fisheries. In press.
#' 
#' Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for
#' data-poor fisheries; forecasting with few data. ICES J. Mar. Sci.
#' doi:10.1093/icesjms/fst232
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




#' A management procedure that incrementally adjusts the TAC according to the
#' mean length of recent catches.
#' 
#' The most biologically precautionary of four adaptive length-based MPs
#' proposed by Geromont and Butterworth 2014. Tested by Carruthers et al. 2015
#' 
#' Tested by Carruthers et al. 2015.
#' 
#' @usage LstepCC4(x, Data, reps = 100, yrsmth = 5, xx=0.3, stepsz=0.05,
#' llim=c(0.96,0.98,1.05))
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param stepsz Parameter controlling the size of the TAC update increment.
#' @param llim A vector of length reference points that determine the
#' conditions for increasing, maintaining or reducing the TAC.
#' @return A numeric vector of TAC recommendations
#' @author T. Carruthers
#' @references Carruthers et al. 2015. Performance evaluation of simple
#' management procedures. Fish and Fisheries. In press.
#' 
#' Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for
#' data-poor fisheries; forecasting with few data. ICES J. Mar. Sci.
#' doi:10.1093/icesjms/fst232
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




#' A management procedure that incrementally adjusts the TAC to reach a target
#' mean length in catches.
#' 
#' The least biologically precautionary of four target length MPs proposed by
#' Geromont and Butterworth 2014. Tested by Carruthers et al. 2015
#' 
#' Tested by Carruthers et al. 2015.
#' 
#' @usage Ltarget1(x, Data, reps = 100, yrsmth = 5, xx=0, xL=1.05)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param xL Parameter controlling the magnitude of the target mean length of
#' catches relative to average length in catches.
#' @return A numeric vector of TAC recommendations
#' @author T. Carruthers
#' @references Carruthers et al. 2015. Performance evaluation of simple
#' management procedures. Fish and Fisheries. In press.
#' 
#' Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for
#' data-poor fisheries; forecasting with few data. ICES J. Mar. Sci.
#' doi:10.1093/icesjms/fst232
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



#' A management procedure that incrementally adjusts the TAC to reach a target
#' mean length in catches.
#' 
#' The most biologically precautionary of four target length MPs proposed by
#' Geromont and Butterworth 2014. Tested by Carruthers et al. 2015
#' 
#' Tested by Carruthers et al. 2015.
#' 
#' @usage Ltarget4(x, Data, reps = 100, yrsmth = 5, xx=0.2, xL=1.15)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param xL Parameter controlling the magnitude of the target mean length of
#' catches relative to average length in catches.
#' @return A numeric vector of TAC recommendations
#' @author T. Carruthers
#' @references Carruthers et al. 2015. Performance evaluation of simple
#' management procedures. Fish and Fisheries. In press.
#' 
#' Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for
#' data-poor fisheries; forecasting with few data. ICES J. Mar. Sci.
#' doi:10.1093/icesjms/fst232
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
#' how informative an estimate of current stock depletion can be. TAC=2*D*AvC
#' 
#' 
#' @usage MCD(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the quota recommendation
#' @return A numeric vector of TAC recommendations
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



#' Mean Catch Depletion
#' 
#' A simple average catch-depletion MP linked to a 40-10 harvest controle rule
#' that was included to demonstrate just how informative an estimate of current
#' stock depletion can be. TAC=d(1-d)AvC
#' 
#' 
#' @usage MCD4010(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the quota recommendation
#' @return A numeric vector of TAC recommendations
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



#' MP using prior for intrinsic rate of increase with a quadratic approximation
#' to surplus production
#' 
#' An MP proposed by Carl Walters that modifies quotas according to trends in
#' apparent surplus production that includes information from a demographically
#' derived prior for intrinsic rate of increase. This is different from
#' Rcontrol because it includes a quadratic approximation of recent trend in
#' surplus production given biomass
#' 
#' 
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrsmth The number of years for smoothing catch and biomass data
#' @param gg A gain parameters
#' @param glim Limits for the change in TAC among years
#' @author C. Walters and T. Carruthers
#' @references Made-up for this package.
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
#' curve. OFL = dep x (1-dep) x r x K x 2
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

#' Surplus Production Stock Reduction Analysis using a mean-length estimate of
#' current stock depletion
#' 
#' A surplus production equivalent of DB-SRA that uses a demographically
#' derived prior for intrinsic rate of increase. A prior for depletion is
#' calculated from a mean-length estimator
#' 
#' 
#' @usage SPSRA_ML(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object (class DLM)
#' @param reps The number of samples of the TAC taken
#' @note The mean length extension was programmed by Gary Nelson as part of his
#' excellent R package 'fishmethods'
#' @author T. Carruthers
#' @references McAllister, M.K., Pikitch, E.K., and Babcock, E.A. 2001. Using
#' demographic methods to construct Bayesian priors for the intrinsic rate of
#' increase in the Schaefer model and implications for stock rebuilding. Can.
#' J. Fish. Aquat. Sci. 58: 1871-1890.
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
#' Yield Per Recruit analysis to get FMSY proxy F01
#' 
#' A simple yield per recruit approximation to FMSY (F01) which is the position
#' of the ascending YPR curve for which dYPR/dF = 0.1(dYPR/d0)
#' 
#' 
#' @usage YPR(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of samples of the TAC
#' @return A numeric vector of TAC samples
#' @note Based on the code of Meaghan Bryan
#' @author Meaghan Bryan and Tom Carruthers
#' @references Beverton and Holt. 1954.
#' @export YPR
YPR <- function(x, Data, reps = 100) {
  # Yield per recruit analysis F01 - Meaghan Bryan for(x in 1:10){
  dependencies = "Data@Mort, Data@CV_Mort, Data@vbK, Data@CV_vbK, Data@vbLinf, Data@CV_vbLinf, Data@vbt0, Data@CV_vbt0, Data@MaxAge, Data@Abun, Data@CV_Abun, Data@wla, Data@wlb"
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
  Ac <- trlnorm(reps, Data@Abun[x], Data@CV_Abun[x])
  FMSY <- YPRopt(Linfc, Kc, t0c, Mdb, a, b, LFS, Data@MaxAge, reps)
  TAC <- Ac * FMSY
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
  # }
}  # end of YPR
class(YPR) <- "MP"



#' Yield Per Recruit analysis to get FMSY proxy F01 paired to a naive catch
#' curve estimate of recent Z
#' 
#' A simple yield per recruit approximation to FMSY (F01) which is the position
#' of the ascending YPR curve for which dYPR/dF = 0.1(dYPR/d0) A naive
#' catch-curve analysis is used to determine recent Z which given M (Mort)
#' gives F and thus abundance = Ct/(1-exp(-F))
#' 
#' 
#' @usage YPR_CC(x, Data, reps = 100, Fmin=0.005)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object (class DLM)
#' @param reps The number of samples of the TAC
#' @param Fmin The minimum fishing mortality rate inferred from the catch-curve
#' analysis
#' @author Meaghan Bryan and T. Carruthers
#' @export YPR_CC
YPR_CC <- function(x, Data, reps = 100, Fmin = 0.005) {
  # for(x in 1:16){
  dependencies = "Data@Mort, Data@CV_Mort, Data@vbK, Data@CV_vbK, Data@vbLinf, Data@CV_vbLinf, Data@vbt0, Data@CV_vbt0, Data@MaxAge, Data@wla, Data@wlb, Data@CAA, Data@Cat"
  Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
    t0c <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
  } else {
    t0c <- rep(Data@vbt0[x], reps)
  }
  t0c[!is.finite(t0c)] <- 0
  LFS <- trlnorm(reps, Data@LFS[x], Data@CV_LFS[x])
  a <- Data@wla[x]
  b <- Data@wlb[x]
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
  FMSY <- YPRopt(Linfc, Kc, t0c, Mdb, a, b, LFS, Data@MaxAge, reps)
  TAC <- Ac * FMSY
  # }
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(YPR_CC) <- "MP"



#' Yield Per Recruit analysis to get FMSY proxy F01 paired with a mean-length
#' estimate of current stock size
#' 
#' A simple yield per recruit approximation to FMSY (F01) which is the position
#' of the ascending YPR curve for which dYPR/dF = 0.1(dYPR/d0) A mean-length
#' estimate of recent Z is used to infer current abundance
#' 
#' 
#' @usage YPR_ML(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of samples of the TAC
#' @note The mean length extension was programmed by Gary Nelson as part of his
#' excellent R package 'fishmethods'
#' @author Meaghan Bryan and T. Carruthers
#' @export YPR_ML
YPR_ML <- function(x, Data, reps = 100) {
  dependencies = "Data@Mort, Data@CV_Mort, Data@vbK, Data@CV_vbK, Data@vbLinf, Data@CV_vbLinf, Data@vbt0, Data@CV_vbt0, Data@MaxAge, Data@wla, Data@wlb, Data@CAL, Data@Cat"
  Mdb <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])
  Linfc <- trlnorm(reps * 10, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps * 10, Data@vbK[x], Data@CV_vbK[x])
  Mdb <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])
  t0c <- -trlnorm(reps * 10, -Data@vbt0[x], Data@CV_vbt0[x])
  t0c[!is.finite(t0c)] <- 0
  LFS <- trlnorm(reps * 10, Data@LFS[x], Data@CV_LFS[x])
  a <- Data@wla[x]
  b <- Data@wlb[x]
  MuC <- Data@Cat[x, length(Data@Cat[x, ])]
  Cc <- trlnorm(reps * 10, MuC, Data@CV_Cat[x])
  Z <- MLne(x, Data, Linfc = Linfc, Kc = Kc, ML_reps = reps * 10, MLtype = "F")
  if (all(is.na(Z))) {
    Rec <- new("Rec")
    Rec@TAC <- TACfilter(rep(NA, reps))
    return(Rec)
  } 
  FM <- Z - Mdb
  Ac <- Cc/(1 - exp(-FM))
  FMSY <- YPRopt(Linfc, Kc, t0c, Mdb, a, b, LFS, Data@MaxAge, reps * 
                   10)
  TAC <- Ac * FMSY
  TAC <- TAC[TAC > 0][1:reps]
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
  
}
class(YPR_ML) <- "MP"