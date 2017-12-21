
#' Example name
#' 
#' A brief description of the management procedure.
#' 
#' @details more details in paragraphs
#' 
#' Second paragraph
#' 
#' italics text: \emph{italics}
#' bold text:  \strong{bold}
#' r code: \code{r_function_call(with = "arguments")}, \code{NULL}, \code{TRUE}
#'
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @author T. Carruthers
#' @family other similar functions in same famile
#' @seealso links to other functions or urls
#'   \code{\link{prod}} for products, \code{\link{cumsum}} for cumulative
#'   sums, and \code{\link{colSums}}/\code{\link{rowSums}} marginal sums over
#'   high-dimensional arrays.
#' @examples 
#' @export 
testFunction <- function(x, Data, reps = 100) {
  rec <- new("Rec") # create recommendation object
  rec@TAC <- trlnorm(reps, Data@OM$A[x] * (1 - exp(-Data@OM$FMSY[x])), 0.01)
  rec
}

## Reference MPs ####

#' A reference FMSY method (uses perfect information about FMSY)
#' 
#' FMSY is taken from the operating model stored at DLM@OM$FMSY
#' 
#' Note that you can out-perform this MP even though it has perfect
#' information of FMSY and current abundance. The requirement for fixed F is
#' actually quite strict and is by no means the upper limit in terms of yield.
#' Don't panic if your method beats this one for yield, especially for
#' short-lived species of high temporal variability in productivity!
#' 
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @author T. Carruthers
#' @export 
FMSYref <- function(x, Data, reps = 100) {
  rec <- new("Rec") # create recommendation object
  rec@TAC <- trlnorm(reps, Data@OM$A[x] * (1 - exp(-Data@OM$FMSY[x])), 0.01)
  rec
  
}
class(FMSYref) <- "MP"




## Output Control MPs ####

#' Average Catch
#'
#' A simple average catch MP that is included to demonstrate a 'status quo' management option
#'
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @author T. Carruthers
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
#' @importFrom utils packageVersion lsf.str read.csv
AvC <- function(x, Data, reps = 100) {
  dependencies = "Data@Cat"
  Rec <- new("Rec")
  Rec@TAC <- rlnorm(reps, log(mean(Data@Cat[x, ], na.rm = T)), 0.2)
  Rec
}
class(AvC) <- "MP"


#' Beddington and Kirkwood life-history MP (simple version)
#' 
#' Sets an OFL according to current abundance and an approximation of Fmax
#' based on length at first capture.
#' 
#' 
#' @usage BK(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object.
#' @param Data A data-limited methods data object.
#' @param reps The number of stochastic samples of the TAC recommendation
#' @note This is the simple version of the BK MP. The paper has a more complex
#' approach that might work better.
#' @author T. Carruthers.
#' @references Beddington, J.R., Kirkwood, G.P., 2005. The estimation of
#' potential yield and stock status using life history parameters. Philos.
#' Trans. R. Soc. Lond. B Biol. Sci. 360, 163-170.
#' @export BK
BK <- function(x, Data, reps = 100) {
  # Beddington and Kirkwood life-history analysis
  # ==============================================
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

#' Beddington and Kirkwood life-history method combined with catch curve
#' analysis
#' 
#' Calculates an OFL using a catch curve estimate of current F and an
#' approximation of FMSY based on length at first capture.
#' 
#' 
#' @usage BK_CC(x, Data, reps = 100, Fmin=0.005)
#' @param x Position in a data-limited methods data object
#' @param Data A data-limited methods data object (class Data)
#' @param reps The number of samples of the TAC recommendation
#' @param Fmin The minimum fishing mortality rate that is derived from the
#' catch-curve (interval censor)
#' @author T. Carruthers
#' @references Beddington, J.R., Kirkwood, G.P., 2005. The estimation of
#' potential yield and stock status using life history parameters. Philos.
#' Trans. R. Soc. Lond. B Biol. Sci. 360, 163-170.
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

#' Beddington and Kirkwood life-history analysis with mean-length estimator of
#' current abundance
#' 
#' Uses an approximation to FMSY based on length at first capture and an
#' estimate of current abundance based on a mean-length estimator.
#' 
#' 
#' @usage BK_ML(x, Data, reps = 100)
#' @param x Position in a data-limited methods data object
#' @param Data A data-limited methods data object (class Data)
#' @param reps The number of samples of the TAC recommendation
#' @note The mean length extension was programmed by Gary Nelson as part of his
#' excellent R package 'fishmethods'
#' @author T. Carruthers
#' @references Beddington, J.R., Kirkwood, G.P., 2005. The estimation of
#' potential yield and stock status using life history parameters. Philos.
#' Trans. R. Soc. Lond. B Biol. Sci. 360, 163-170.
#' @export BK_ML
BK_ML <- function(x, Data, reps = 100) {
  dependencies = "Data@LFC, Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@CAL, Data@Mort"
  Lc <- trlnorm(reps * 10, Data@LFC[x], 0.2)
  Linfc <- trlnorm(reps * 10, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps * 10, Data@vbK[x], Data@CV_vbK[x])
  Mdb <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])
  Z <- MLne(x, Data, Linfc = Linfc, Kc = Kc, ML_reps = reps * 10, MLtype = "F")
  if (all(is.na(Z))) return(rep(NA, reps))
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
#' The TAC is the average catch over last yrsmth years.
#' 
#' This is one of four constant catch rules of Geromont and Butterworth 2014.
#' 
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrsmth Years over which to calculate mean catches
#' @param xx Parameter controlling the TAC. Mean catches are multiplied by
#' (1-xx)
#' @return A numeric vector of TAC recommendations
#' @author T. Carruthers
#' @references Geromont, H.F., Butterworth, D.S. 2014. Generic management
#' procedures for data-poor fisheries; forecasting with few data. ICES J. Mar.
#' Sci. doi:10.1093/icesjms/fst232
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


#' Constant catch management procedure of Geromont and Butterworth (2014)
#' 
#' The TAC is the average catch over last yrsmth years reduced by 30%.
#' 
#' This is one of four constant catch MPs of Geromont and Butterworth 2014.
#' 
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrsmth Years over which to average catches
#' @param xx Parameter controlling the TAC. Mean catches are multiplied by
#' (1-xx)
#' @return A numeric vector of TAC recommendations
#' @author T. Carruthers
#' @references Geromont, H.F., Butterworth, D.S. 2014. Generic management
#' procedures for data-poor fisheries; forecasting with few data. ICES J. Mar.
#' Sci. doi:10.1093/icesjms/fst232
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


#' Depletion-Based Stock Reduction Analysis
#' 
#' User prescribed BMSY/B0, M, FMSY/M are used to find B0 and therefore the OFL
#' by back-constructing the stock to match a user specified level of stock
#' depletion (OFL = M * FMSY/M * depletion* B0).
#' 
#' You specify a range of stock depletion and, given historical catches DB-SRA
#' calculates what unfished biomass must have been to get you here given
#' samples for M, FMSY relative to M and also BMSY relative to Bunfished.
#' 
#' @param x A position in a data-limited methods object.
#' @param Data A data-limited methods object.
#' @param reps The number of samples of the TAC (OFL) recommendation.
#' @return A vector of TAC (OFL) values.
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
#' @references Dick, E.J., MacCall, A.D., 2011. Depletion-Based Stock Reduction
#' Analysis: A catch-based method for determining sustainable yields for
#' data-poor fish stocks. Fish. Res. 110, 331-341.
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
      if (Max <= 0.05) 
        BMSY_K <- 0.05
      if (Min >= 0.95) 
        BMSY_K <- 0.95
    }
    if (!is.na(tryBMSY_K))  BMSY_K <- tryBMSY_K
    
    adelay <- max(floor(iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],  Data@L50[x])), 1)
    opt <- optimize(DBSRAopt, log(c(0.01 * mean(C_hist), 1000 * mean(C_hist))), C_hist = C_hist, 
                    nys = length(C_hist), Mdb = Mdb, FMSY_M = FMSY_M, BMSY_K = BMSY_K, 
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
    
    # print(DBSRAcount)
  }  # end of reps
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
 
}  # end of DBSRA_apply
class(DBSRA) <- "MP"


#' Depletion-Based Stock Reduction Analysis assuming 40 per cent stock
#' depletion
#' 
#' DBSRA assuming that current stock depletion is exactly 40 per cent of
#' unfished stock levels.
#' 
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @note A 40 percent assumption for current depletion is more or less the most
#' optimistic state for a stock (ie very close to BMSY/B0 for many stocks).
#' @author T. Carruthers.
#' @references Dick, E.J., MacCall, A.D., 2010. Estimates of sustainable yield
#' for 50 data-poor stocks in the Pacific Coast groundfish fishery management
#' plan. Technical memorandum. Southwest fisheries Science Centre, Santa Cruz,
#' CA. National Marine Fisheries Service, National Oceanic and Atmospheric
#' Administration of the U.S. Department of Commerce. NOAA-TM-NMFS-SWFSC-460.
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
    opt <- optimize(DBSRAopt, log(c(0.1 * mean(C_hist), 1000 * mean(C_hist))), 
                    C_hist = C_hist, nys = length(C_hist), Mdb = Mdb, FMSY_M = FMSY_M, 
                    BMSY_K = BMSY_K, Bt_K = Bt_K, adelay = adelay, tol = 0.01)
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



#' Depletion-Based Stock Reduction Analysis paired with 40-10 harvest control
#' rule
#' 
#' User prescribed BMSY/B0, M, FMSY/M are used to find B0 and therefore the OFL
#' by back-constructing the stock to match a user specified level of stock
#' depletion (OFL = M * FMSY/M * depletion* B0). In this method DBSRA is paried
#' with the 40-10 rule that throttles back the OFL to zero at 10 percent of
#' unfished biomass.
#' 
#' 
#' @usage DBSRA4010(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @author T. Carruthers
#' @references Dick, E.J., MacCall, A.D., 2011. Depletion-Based Stock Reduction
#' Analysis: A catch-based method for determining sustainable yields for
#' data-poor fish stocks. Fish. Res. 110, 331-341.
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
    opt <- optimize(DBSRAopt, log(c(0.01 * mean(C_hist), 1000 * mean(C_hist))), 
                    C_hist = C_hist, nys = length(C_hist), Mdb = Mdb, FMSY_M = FMSY_M, 
                    BMSY_K = BMSY_K, Bt_K = Bt_K, adelay = adelay, tol = 0.01)
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



#' Depletion-Based Stock Reduction Analysis using mean length estimator of
#' stock depletion
#' 
#' DBSRA using the mean length estimator to calculate current stock depletion.
#' 
#'
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the quota recommendation
#' @note The mean length extension was programmed by Gary Nelson as part of his
#' excellent R package 'fishmethods'
#' @author T. Carruthers
#' @references Dick, E.J., MacCall, A.D., 2011. Depletion-Based Stock Reduction
#' Analysis: A catch-based method for determining sustainable yields for
#' data-poor fish stocks. Fish. Res. 110, 331-341.
#' @export 
DBSRA_ML <- function(x, Data, reps = 100) {
  dependencies = "Data@Cat, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@BMSY_B0, Data@CV_BMSY_B0, Data@L50, Data@CAL, Data@Year, Data@Cat"
  C_hist <- Data@Cat[x, ]
  TAC <- rep(NA, reps)
  DBSRAcount <- 1
  maxIts <- 200
  nIts <- 0
  if (is.na(Data@Dep[x]) | is.na(Data@CV_Dep[x])) return(NA)
  while (DBSRAcount < (reps + 1) & nIts < maxIts) {
    Linfc <- trlnorm(1, Data@vbLinf[x], Data@CV_vbLinf[x])
    Kc <- trlnorm(1, Data@vbK[x], Data@CV_vbK[x])
    Mdb <- trlnorm(100, Data@Mort[x], Data@CV_Mort[x])
    Mdb <- Mdb[Mdb < 0.9][1]  # !!!! maximum M is 0.9   interval censor
    if (is.na(Mdb)) Mdb <- 0.9  # !!!! maximum M is 0.9   absolute limit
    Z <- MLne(x, Data, Linfc = Linfc, Kc = Kc, ML_reps = 1, MLtype = "dep")
    if (all(is.na(Z))) return(rep(NA, reps))
    FM <- Z - Mdb
    FM[FM < 0] <- 0.01
    nyears <- length(Data@Year)
    Ct1 <- mean(Data@Cat[x, 1:3])
    Ct2 <- mean(Data@Cat[x, (nyears - 2):nyears])
    # dep<-c(Ct1,Ct2)/(1-exp(-FM[,c(1,2)]))
    dep <- c(Ct1, Ct2)/(1 - exp(-FM))
    Bt_K <- dep[2]/dep[1]
    
    if (Bt_K < 0.01) Bt_K <- 0.01  # interval censor / temporary hack to avoid doing multiple depletion estimates that would take far too long
    if (Bt_K > 0.99) Bt_K <- 0.99  # interval censor / temporary hack to avoid doing multiple depletion estimates that would take far too long
    
    
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
    if (all(is.na(BMSY_K))) return(rep(NA, reps))
    adelay <- max(floor(iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],  Data@L50[x])), 1)
    opt <- optimize(DBSRAopt, log(c(0.1 * mean(C_hist), 1000 * mean(C_hist))), 
                    C_hist = C_hist, nys = length(C_hist), Mdb = Mdb, FMSY_M = FMSY_M, 
                    BMSY_K = BMSY_K, Bt_K = Bt_K, adelay = adelay, tol = 0.01)
    nIts <- nIts + 1
    if (opt$objective < 0.1) {
      Kc <- exp(opt$minimum)
      BMSYc <- Kc * BMSY_K
      FMSYc <- Mdb * FMSY_M
      UMSYc <- (FMSYc/(FMSYc + Mdb)) * (1 - exp(-(FMSYc + Mdb)))
      MSYc <- Kc * BMSY_K * UMSYc
      TAC[DBSRAcount] <- UMSYc * Kc * Bt_K
      DBSRAcount <- DBSRAcount + 1
      nIts <- 0
    }
  }  # end of reps
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}
class(DBSRA_ML) <- "MP"


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
#' @references MacCall, A.D., 2009. Depletion-corrected average catch: a simple
#' formula for estimating sustainable yields in data-poor situations. ICES J.
#' Mar. Sci. 66, 2267-2271.
#' @export DCAC
DCAC <- function(x, Data, reps = 100) {
  dependencies = "Data@AvC, Data@t, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@Dt, Data@CV_Dt, Data@BMSY_B0, Data@CV_BMSY_B0"
  C_tot <- Data@AvC[x] * Data@t[x]
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
  FMSY_M <- trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])  # standard deviation of 0.2 - referred to as 'standard error' in MacCall 2009
  Bt_K <- trlnorm(reps, Data@Dt[x], Data@CV_Dt[x])
  if (any(is.na(c(Data@BMSY_B0[x], Data@CV_BMSY_B0[x])))) 
    return(NA)
  BMSY_K <- rbeta(reps, alphaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
                                    Data@CV_BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
                                                                    Data@CV_BMSY_B0[x]))  #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(C_tot/(Data@t[x] + ((1 - Bt_K)/(BMSY_K * FMSY_M * Mdb))))
  Rec
}  # end of DCAC
class(DCAC) <- "MP"


#' Depletion Corrected Average Catch paired with the 40-10 rule
#' 
#' A method of calculating an MSY proxy (FMSY * BMSY and therefore the OFL at
#' most productive stock size) based on average catches accounting for the
#' windfall catch that got the stock down to BMSY levels. In this method DCAC
#' is paired with the 40-10 rule that throttles back the OFL to zero at 10
#' percent of unfished stock size (the OFL is not subject to downward
#' adjustment above 40 percent unfished)
#' 
#' 
#' @usage DCAC4010(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @note DCAC can overfish below BMSY levels. The 40-10 harvest control rule
#' largely resolves this problem providing an MP with surprisingly good
#' performance even at low stock levels.
#' @author T. Carruthers
#' @references MacCall, A.D., 2009. Depletion-corrected average catch: a simple
#' formula for estimating sustainable yields in data-poor situations. ICES J.
#' Mar. Sci. 66, 2267-2271.
#' @export DCAC4010
DCAC4010 <- function(x, Data, reps = 100) {
  dependencies = "Data@AvC, Data@t, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@Dt, Data@CV_Dt, Data@BMSY_B0, Data@CV_BMSY_B0"
  C_tot <- Data@AvC[x] * Data@t[x]
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
  FMSY_M <- trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])  # standard deviation of 0.2 - referred to as 'standard error' in MacCall 2009
  Bt_K <- trlnorm(reps, Data@Dt[x], Data@CV_Dt[x])
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



#' Depletion Corrected Average Catch assuming 40 per cent stock depletion
#' 
#' DCAC assuming that current stock biomass is exactly 40 per cent of unfished
#' levels.
#' 
#' 
#' @usage DCAC_40(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @note The 40 percent depletion assumption doesn't really affect DCAC that
#' much as it already makes TAC recommendations that are quite MSY-like.
#' @author T. Carruthers
#' @references MacCall, A.D., 2009. Depletion-corrected average catch: a simple
#' formula for estimating sustainable yields in data-poor situations. ICES J.
#' Mar. Sci. 66, 2267-2271.
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



#' Depletion-Based Stock Reduction Analysis using mean-length estimator of
#' current depletion
#' 
#' DCAC that uses the mean length estimator to calculate current stock
#' depletion.
#' 
#' 
#' @usage DCAC_ML(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @note The mean length extension was programmed by Gary Nelson as part of his
#' excellent R package 'fishmethods'
#' @author T. Carruthers
#' @references MacCall, A.D., 2009. Depletion-corrected average catch: a simple
#' formula for estimating sustainable yields in data-poor situations. ICES J.
#' Mar. Sci. 66, 2267-2271.
#' @export DCAC_ML
DCAC_ML <- function(x, Data, reps = 100) {
  
  dependencies = "Data@AvC, Data@t, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@BMSY_B0, Data@CV_BMSY_B0, Data@Year, Data@CAL, Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK"
  if (is.na(Data@BMSY_B0[x]) | is.na(Data@CV_BMSY_B0[x])) 
    return(NA)
  if (is.na(Data@FMSY_M[x]) | is.na(Data@CV_FMSY_M[x])) 
    return(NA)
  C_tot <- Data@AvC[x] * Data@t[x]
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # default CV of 0.5 as in MacCall 2009
  FMSY_M <- trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])  # standard deviation of 0.2 - referred to as 'standard error' in MacCall 2009
  Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  Z <- MLne(x, Data, Linfc = Linfc, Kc = Kc, ML_reps = reps, MLtype = "dep")
  if (all(is.na(Z))) return(rep(NA, reps))
  FM <- Z - Mdb
  nyears <- length(Data@Year)
  Ct1 <- mean(Data@Cat[x, 1:3])
  Ct2 <- mean(Data@Cat[x, (nyears - 2):nyears])
  dep <- rep(c(Ct1, Ct2), each = reps)/(1 - exp(-FM))
  if (reps == 1) 
    Bt_K <- dep[2]/dep[1]
  if (reps > 1) 
    Bt_K <- dep[, 2]/dep[, 1]
  if (any(is.na(c(Data@BMSY_B0[x], Data@CV_BMSY_B0[x])))) 
    return(NA)
  BMSY_K <- rbeta(reps, alphaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
                                    Data@CV_BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
                                                                    Data@CV_BMSY_B0[x]))  #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
  TAC <- C_tot/(Data@t[x] + ((1 - Bt_K)/(BMSY_K * FMSY_M * Mdb)))
  
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  Rec
}  # end of DCAC_ML
class(DCAC_ML) <- "MP"


#' Delay - Difference Stock Assessment with UMSY and MSY leading
#' 
#' A simple delay-difference assessment that estimates the TAC using a
#' time-series of catches and a relative abundance index.
#' 
#' 
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @return A numeric vector of TAC recommendations
#' @note This DD model is observation error only and has does not estimate
#' process error (recruitment deviations). Similar to many other assessment
#' models it depends on a whole host of dubious assumptions such as temporally
#' stationary productivity and proportionality between the abundance index and
#' real abundance. Unsurprisingly the extent to which these assumptions are
#' violated tends to be the biggest driver of performance for this method.
#' @author T. Carruthers
#' @references Method based on equations of Carl Walters (bug him with
#' questions and expect colourful responses)
#' @export 
DD <- function(x, Data, reps = 100) {
  # for(x in 1:nsim){
  dependencies = "Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@vbt0, Data@CV_vbt0, Data@Mort, Data@CV_Mort, Data@wla, Data@wlb"
  Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
    t0c <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
  } else {
    t0c <- rep(Data@vbt0[x], reps)
  }
  t0c[!is.finite(t0c)] <- 0
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
  a <- Data@wla[x]
  b <- Data@wlb[x]
  
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
  params <- log(c(Data@Mort[x], mean(C_hist, na.rm = T), Data@Mort[x]))
  k_DD <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)  
  k_DD[k_DD > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho_DD <- (wa[k_DD + 2] - Winf)/(wa[k_DD + 1] - Winf)
  Alpha_DD <- Winf * (1 - Rho_DD)
  So_DD <- exp(-Data@Mort[x])  # get So survival rate
  wa_DD <- wa[k_DD]
  UMSYprior <- c(1 - exp(-Data@Mort[x] * 0.5), 0.3)
  opt <- optim(params, DD_R, opty = 1, So_DD = So_DD, Alpha_DD = Alpha_DD, 
               Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, wa_DD = wa_DD, E_hist = E_hist, 
               C_hist = C_hist, UMSYprior = UMSYprior, method = "L-BFGS-B", lower = log(exp(params)/20), 
               upper = log(exp(params) * 20), hessian = TRUE)
  
  # Catfit<-DD_R(opt$par,opty=3,So_DD=So_DD,Alpha_DD=Alpha_DD,Rho_DD=Rho_DD,ny_DD=ny_DD,k_DD=k_DD,wa_DD=wa_DD,E_hist=E_hist,C_hist=C_hist,UMSYprior=UMSYprior)
  # plot(Catfit[,1],ylim=c(0,max(Catfit))) lines(Catfit[,2],col='red')
  
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


#' Delay - Difference Stock Assessment with UMSY and MSY leading coupled with a
#' 40-10 harvest control rule
#' 
#' A simple delay-difference assessment that estimates the OFL using a
#' time-series of catches and a relative abundance index. In this version of
#' the DD MP a 40-10 rule is imposed over the OFL recommendation.
#' 
#' 
#' @usage DD4010(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @return A numeric vector of TAC recommendations
#' @author T. Carruthers
#' @references Method based on equations of Carl Walters
#' @export DD4010
DD4010 <- function(x, Data, reps = 100) {
  dependencies = "Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@vbt0, Data@CV_vbt0, Data@Mort, Data@CV_Mort. Data@wla, Data@ wlb"
  Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
    t0c <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
  } else {
    t0c <- rep(Data@vbt0[x], reps)
  }
  t0c[!is.finite(t0c)] <- 0
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
  a <- Data@wla[x]
  b <- Data@wlb[x]
  
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
  params <- log(c(Data@Mort[x], mean(C_hist, na.rm = T), Data@Mort[x]))
  k_DD <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)
  k_DD[k_DD > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho_DD <- (wa[k_DD + 2] - Winf)/(wa[k_DD + 1] - Winf)
  Alpha_DD <- Winf * (1 - Rho_DD)
  So_DD <- exp(-Data@Mort[x])  # get So survival rate
  wa_DD <- wa[k_DD]
  UMSYprior <- c(1 - exp(-Data@Mort * 0.5), 0.3)
  opt <- optim(params, DD_R, opty = 1, So_DD = So_DD, Alpha_DD = Alpha_DD, 
               Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, wa_DD = wa_DD, E_hist = E_hist, 
               C_hist = C_hist, UMSYprior = UMSYprior, method = "L-BFGS-B", lower = log(exp(params)/20), 
               upper = log(exp(params) * 20), hessian = TRUE)
  
  # Catfit<-DD_R(opt$par,opty=3,So_DD=So_DD,Alpha_DD=Alpha_DD,Rho_DD=Rho_DD,ny_DD=ny_DD,k_DD=k_DD,wa_DD=wa_DD,E_hist=E_hist,C_hist=C_hist,UMSYprior=UMSYprior)
  # plot(Catfit[,1],ylim=c(0,max(Catfit))) lines(Catfit[,2],col='red')
  
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
  
  if (all(is.na(Z))) return(rep(NA, reps))
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
#' A simple method that tends to outperform many other approaches alarmingly
#' often even when current biomass is relatively poorly known. The low stock
#' crash potential is largely due to the quite large difference between Fmax
#' and FMSY for most stocks.
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of samples of the TAC recommendation
#' @return A numeric vector of TAC recommendations
#' @author T. Carruthers
#' @references Gulland, J.A., 1971. The fish resources of the ocean. Fishing
#' News Books, West Byfleet, UK.
#' 
#' Martell, S., Froese, R., 2012. A simple method for estimating MSY from catch
#' and resilience. Fish Fish. doi: 10.1111/j.1467-2979.2012.00485.x.
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

#' An FMSY/M ratio method paired with the 40-10 rule
#' 
#' Calculates the OFL based on a fixed ratio of FMSY to M multiplied by a
#' current estimate of abundance. In this method DBSRA is paired with the 40-10
#' rule that throttles back the OFL to zero at 10 percent of unfished biomass.
#' 
#' 
#' @usage Fratio4010(x, Data, reps = 100)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @author T. Carruthers
#' @references Gulland, J.A., 1971. The fish resources of the ocean. Fishing
#' News Books, West Byfleet, UK.
#' 
#' Martell, S., Froese, R., 2012. A simple method for estimating MSY from catch
#' and resilience. Fish Fish. doi: 10.1111/j.1467-2979.2012.00485.x.
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

#' A data-limited method that uses FMSY/M ratio and a naive catch-curve
#' estimate of recent Z
#' 
#' Calculates the OFL based on a fixed ratio of FMSY to M and a catch curve
#' estimate of current stock size.
#' 
#' 
#' @usage Fratio_CC(x, Data, reps = 100, Fmin = 0.005)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param Fmin Minimum current fishing mortality rate for the catch-curve
#' analysis
#' @author T. Carruthers
#' @references Gulland, J.A., 1971. The fish resources of the ocean. Fishing
#' News Books, West Byfleet, UK.
#' 
#' Martell, S., Froese, R., 2012. A simple method for estimating MSY from catch
#' and resilience. Fish Fish. doi: 10.1111/j.1467-2979.2012.00485.x.
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

#' An FMSY/M ratio MP that uses a mean length estimator of recent Z
#' 
#' Calculates the OFL based on a fixed ratio of FMSY/M and an estimate of
#' current stock size from a mean-length estimator.
#' 
#' 
#' @usage Fratio_ML(x, Data, reps = 100)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @note The mean length extension was programmed by Gary Nelson as part of his
#' excellent R package 'fishmethods'
#' @author T. Carruthers
#' @references Gulland, J.A., 1971. The fish resources of the ocean. Fishing
#' News Books, West Byfleet, UK.
#' 
#' Martell, S., Froese, R., 2012. A simple method for estimating MSY from catch
#' and resilience. Fish Fish. doi: 10.1111/j.1467-2979.2012.00485.x.
#' @export Fratio_ML
Fratio_ML <- function(x, Data, reps = 100) {
  dependencies = " Data@FMSY_M, Data@CV_FMSY_M, Data@Mort, Data@CV_Mort, Data@Cat, Data@CV_Cat, Data@CAL"
  MuC <- Data@Cat[x, length(Data@Cat[x, ])]
  Cc <- trlnorm(reps * 10, MuC, Data@CV_Cat[x])
  Mdb <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
  Linfc <- trlnorm(reps * 10, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps * 10, Data@vbK[x], Data@CV_vbK[x])
  Z <- MLne(x, Data, Linfc = Linfc, Kc = Kc, ML_reps = reps * 10, MLtype = "F")
  if (all(is.na(Z))) return(rep(NA, reps))
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




## Size limit MPs ####

#' A data-limited method in which fishing retention is set according to the
#' maturity curve NEW
#' 
#' An example of the implementation of input controls in the DLM toolkit, where
#' retention-at-length is set equivalent to maturity-at-length
#' 
#' @param x A position in a data-limited methods object
#' @param Data A data-limited methods object
#' @param ... Optional additional arguments that are ignored. Note arguments
#' \code{reps} or \code{...} are required for all input controls
#' @return A input control recommendation object 
#' @author T. Carruthers
#' @references Made-up for this package
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


#' A data-limited method in which fishing vulnerability is set slightly higher
#' than the maturity curve
#' 
#' An example of the implementation of input controls in the DLM toolkit, where
#' selectivity-at-length is set slightly higher than the maturity-at-length
#' 
#' 
#' @usage matlenlim2(x, Data, ...)
#' @param x A position in a data-limited methods object
#' @param Data A data-limited methods object
#' @param ... Optional additional arguments that are ignored. Note arguments
#' \code{reps} or \code{...} are required for all input controls
#' @return A vector of input control recommendations, with values for length at
#' first capture and full selection
#' @author A. Hordyk
#' @references Made-up for this package
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



#' This input control sets the minimum length of fish caught to a fraction of
#' the length that maximises the biomass, Lopt.
#' 
#' This aim of this simple MP is restrict the catch of small fish to rebuild
#' the stock biomass towards the optimal length, Lopt, expressed in terms of
#' the growth parameters Lopt=b/(M/k+b) (Hordyk et al. (2014)
#' 
#' 
#' @usage minlenLopt1(x, Data, reps = 100, buffer = 0.1)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param buffer Parameter controlling the fraction of Lopt to set the minimum
#' length of fish caught: minlen=Lopt*(0.7+buffer).
#' @return The length at first caprture, LFC, and length at full selectivity
#' @author HF Geromont
#' @references Hordyk, A., Ono, K., Sainsbury, K., Loneragan, N., and J.
#' Prince. 2014. Some explorations of the life history ratios to describe
#' length composition, spawning-per-recruit, and the spawning potential ratio
#' ICES Journal of Marine Science, doi:10.1093/icesjms/fst235.
#' @export minlenLopt1
minlenLopt1 <- function(x, Data, reps = 100, buffer = 0.1) {
  
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


#' An data-limited method which sets a slot limit
#' 
#' An example of the implementation of input controls in the DLM toolkit, where
#' selectivity-at-length is set using a slot limit; that is, a minimum and
#' maximum legal length.  The maximum limit is set here, quite arbitrarily, as
#' the 75th percentile between the new minimum legal length and the estimated
#' asymptotic length.
#' 
#' 
#' @usage slotlim(x, Data, ...)
#' @param x A position in a data-limited methods object
#' @param Data A data-limited methods object
#' @param ... Optional additional arguments that are ignored. Note arguments
#' \code{reps} or \code{...} are required for all input controls
#' @return An object of class 'Rec'
#' @author A. Hordyk
#' @references Made-up for this package
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

#' An marine reserve in area 1 with full reallocation of fishing effort
#' 
#' A spatial control that prevents fishing in area 1 and reallocates this
#' fishing effort to area 2.
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
  rec@Spatial <- c(0,1)
  
  # other slots aren't specified so remain unchanged
  return(rec)
}
class(MRreal) <- "MP"


#' An marine reserve in area 1 with no spatial reallocation of fishing effort
#' 
#' A spatial control that prevents fishing in area 1 and does not reallocate
#' this fishing effort to area 2.
#' 
#' 
#' @usage MRnoreal(x, Data, ...)
#' @param x A position in data / simulation object DLM
#' @param Data A data limited methods data object
#' @param ... Optional additional arguments that are ignored. Note arguments
#' \code{reps} or \code{...} are required for all input controls
#' @author T. Carruthers
#' @export MRnoreal
MRnoreal <- function(x, Data, ...) {
  # A Marine reserve in area 1 with no spatial reallocation of effort
  
  rec <- new("Rec") # create recommendation object
  rec@Allocate <- 0
  rec@Spatial <- c(0,1)
  
  # other slots aren't specified so remain unchanged
  return(rec)
}
class(MRnoreal) <- "MP"



# --- Effort Control MPs ----

#' Fishing at current effort levels
#' 
#' Constant fishing effort set at final year of historical simulations subject
#' to changes in catchability determined by OM@qinc and interannual variability
#' in catchability determined by OM@qcv. This MP is intended to represent a
#' 'status quo' management approach.
#' 
#' @param x A position in a data-limited methods data object.
#' @param Data A data-limited methods data object.
#' @param ... Optional additional arguments that are ignored. Note arguments
#' \code{reps} or \code{...} are required for all input controls
#' @note Made up for this package.
#' @author T. Carruthers.
#' @export 
curE <- function(x, Data, ...) {
  # current effort
  rec <- new("Rec") # create recommendation object
  rec@Effort <- 1
  rec
}
class(curE) <- "MP"


#' Fishing at 75 per cent of current effort levels
#' 
#' Constant fishing effort set at 75 per cent of final year of historical
#' simulations subject to changes in catchability determined by OM@qinc and
#' interannual variability in catchability determined by OM@qcv. This MP is
#' intended to represent a 'status quo' management approach.
#' 
#' @param x A position in a data-limited methods data object.
#' @param Data A data-limited methods data object.
#' @param ... Optional additional arguments that are ignored. Note arguments
#' \code{reps} or \code{...} are required for all input controls
#' @note Made up for this package.
#' @author T. Carruthers.
#' @export 
curE75 <- function(x, Data, ...) {
  # 75% current effort
  rec <- new("Rec") # create recommendation object
  rec@Effort <- 0.75
  rec
}
class(curE75) <- "MP"


#' Effort control version of DD - Delay - Difference Stock Assessment with UMSY
#' and MSY leading
#' 
#' A simple delay-difference assessment that estimates and recommends FMSY
#' using a time-series of catches and a relative abundance index.
#' 
#' 
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @note This DD model is observation error only and has does not estimate
#' process error (recruitment deviations). Similar to many other assessment
#' models it depends on a whole host of dubious assumptions such as temporally
#' stationary productivity and proportionality between the abundance index and
#' real abundance. Unsurprisingly the extent to which these assumptions are
#' violated tends to be the biggest driver of performance for this method.
#' @author T. Carruthers
#' @references Method based on equations of Carl Walters (bug him with
#' questions and expect colourful responses)
#' @export 
DDe <- function(x, Data, reps = 100) {
  # for(x in 1:nsim){
  dependencies = "Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@vbt0, Data@CV_vbt0, Data@Mort, Data@CV_Mort. Data@wla, Data@ wlb"
  Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
    t0c <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
  } else {
    t0c <- rep(Data@vbt0[x], reps)
  }
  t0c[!is.finite(t0c)] <- 0
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
  a <- Data@wla[x]
  b <- Data@wlb[x]
  
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
  E_hist <- C_hist/Data@Ind[x, yind]
  E_hist <- E_hist/mean(E_hist)
  ny_DD <- length(C_hist)
  params <- log(c(Data@Mort[x], mean(C_hist, na.rm = T), Data@Mort[x]))
  k_DD <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)  --
  k_DD[k_DD > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho_DD <- (wa[k_DD + 2] - Winf)/(wa[k_DD + 1] - Winf)
  Alpha_DD <- Winf * (1 - Rho_DD)
  So_DD <- exp(-Data@Mort[x])  # get So survival rate
  wa_DD <- wa[k_DD]
  UMSYprior <- c(1 - exp(-Data@Mort[x] * 0.5), 0.3)
  opt <- optim(params, DD_R, opty = 1, So_DD = So_DD, Alpha_DD = Alpha_DD, 
               Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, wa_DD = wa_DD, E_hist = E_hist, 
               C_hist = C_hist, UMSYprior = UMSYprior, method = "L-BFGS-B", lower = log(exp(params)/20), 
               upper = log(exp(params) * 20), hessian = TRUE)
  
  U_hist <- 1 - exp(-exp(opt$par[3]) * E_hist)
  
  Allocate <- 1
  eff <- exp(opt$par[1])/U_hist[Data@LHYear]
  eff[!is.finite(eff)] <- 0.01
  eff[eff > 1e+05] <- 0.01
  rec <- new("Rec")
  rec@Effort <- max(0.01, eff)
  rec 
}
class(DDe) <- "MP"

#' Effort searching version of DD - Delay - Difference Stock Assessment with
#' UMSY and MSY leading that fishes at 75 per cent of FMSY
#' 
#' A simple delay-difference assessment that estimates FMSY using a time-series
#' of catches and a relative abundance index. The MP provides a change in
#' effort in the direction of FMSY up to a maximum change of 10 percent.
#' 
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @param LB The lowest permitted factor of previous fishing effort
#' @param UB The highest permitted factor of previous fishing effort
#' @note This DD model is observation error only and has does not estimate
#' process error (recruitment deviations). Similar to many other assessment
#' models it depends on a whole host of dubious assumptions such as temporally
#' stationary productivity and proportionality between the abundance index and
#' real abundance. Unsurprisingly the extent to which these assumptions are
#' violated tends to be the biggest driver of performance for this method.
#' @author T. Carruthers
#' @references Method based on equations of Carl Walters (bug him with
#' questions and expect colourful responses)
#' @export DDes
DDes <- function(x, Data, reps = 100, LB = 0.9, UB = 1.1) {
  # for(x in 1:nsim){
  dependencies = "Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@vbt0, Data@CV_vbt0, Data@Mort, Data@CV_Mort. Data@wla, Data@ wlb"
  Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
    t0c <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
  } else {
    t0c <- rep(Data@vbt0[x], reps)
  }
  t0c[!is.finite(t0c)] <- 0
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
  a <- Data@wla[x]
  b <- Data@wlb[x]
  
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
  E_hist <- C_hist/Data@Ind[x, yind]
  E_hist <- E_hist/mean(E_hist)
  ny_DD <- length(C_hist)
  params <- log(c(Data@Mort[x], mean(C_hist, na.rm = T), Data@Mort[x]))
  k_DD <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)  --
  k_DD[k_DD > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho_DD <- (wa[k_DD + 2] - Winf)/(wa[k_DD + 1] - Winf)
  Alpha_DD <- Winf * (1 - Rho_DD)
  So_DD <- exp(-Data@Mort[x])  # get So survival rate
  wa_DD <- wa[k_DD]
  UMSYprior <- c(1 - exp(-Data@Mort[x] * 0.5), 0.3)
  opt <- optim(params, DD_R, opty = 1, So_DD = So_DD, Alpha_DD = Alpha_DD, 
               Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, wa_DD = wa_DD, E_hist = E_hist, 
               C_hist = C_hist, UMSYprior = UMSYprior, method = "L-BFGS-B", lower = log(exp(params)/20), 
               upper = log(exp(params) * 20), hessian = TRUE)
  
  U_hist <- 1 - exp(-exp(opt$par[3]) * E_hist)
  fac <- exp(opt$par[1])/U_hist[Data@LHYear]  # ratio of UMSY to reference U
  fac <- fac * (U_hist[Data@LHYear]/U_hist[length(U_hist)])  # ratio of last U to reference U
  
  if (fac < LB) 
    fac <- LB
  if (fac > UB) 
    fac <- UB
  
  rec <- new("Rec")
  rec@Effort <- max(0.01, Data@MPeff[x] * fac)
  rec
  
}
class(DDes) <- "MP"



#' Effort control version of DD - Delay - Difference Stock Assessment with UMSY
#' and MSY leading that fishes at 75 per cent of FMSY
#' 
#' A simple delay-difference assessment that estimates and recommends 75 per
#' cent FMSY using a time-series of catches and a relative abundance index.
#' 
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @note This DD model is observation error only and has does not estimate
#' process error (recruitment deviations). Similar to many other assessment
#' models it depends on a whole host of dubious assumptions such as temporally
#' stationary productivity and proportionality between the abundance index and
#' real abundance. Unsurprisingly the extent to which these assumptions are
#' violated tends to be the biggest driver of performance for this method.
#' @author T. Carruthers
#' @references Method based on equations of Carl Walters (bug him with
#' questions and expect colourful responses)
#' @export 
DDe75 <- function(x, Data, reps = 100) {
  # for(x in 1:nsim){
  dependencies = "Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@vbt0, Data@CV_vbt0, Data@Mort, Data@CV_Mort. Data@wla, Data@ wlb"
  Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
    t0c <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
  } else {
    t0c <- rep(Data@vbt0[x], reps)
  }
  t0c[!is.finite(t0c)] <- 0
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
  a <- Data@wla[x]
  b <- Data@wlb[x]
  
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
  E_hist <- C_hist/Data@Ind[x, yind]
  E_hist <- E_hist/mean(E_hist)
  ny_DD <- length(C_hist)
  params <- log(c(Data@Mort[x], mean(C_hist, na.rm = T), Data@Mort[x]))
  k_DD <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)  --
  k_DD[k_DD > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho_DD <- (wa[k_DD + 2] - Winf)/(wa[k_DD + 1] - Winf)
  Alpha_DD <- Winf * (1 - Rho_DD)
  So_DD <- exp(-Data@Mort[x])  # get So survival rate
  wa_DD <- wa[k_DD]
  UMSYprior <- c(1 - exp(-Data@Mort[x] * 0.5), 0.3)
  opt <- optim(params, DD_R, opty = 1, So_DD = So_DD, Alpha_DD = Alpha_DD, 
               Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, wa_DD = wa_DD, E_hist = E_hist, 
               C_hist = C_hist, UMSYprior = UMSYprior, method = "L-BFGS-B", lower = log(exp(params)/20), 
               upper = log(exp(params) * 20), hessian = TRUE)
  
  U_hist <- 1 - exp(-exp(opt$par[3]) * E_hist)
  
  Allocate <- 1
  eff <- 0.75 *exp( opt$par[1])/U_hist[Data@LHYear]
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



## Combined Management MPs ####


