# Output methods

#' Check Data object is valid for a MP
#' 
#' Checks that all slots in Data object required by the MP contain 
#' finite values
#' 
#' @param Data An object of class Data
#' @param dependencies A string of slots in the Data object required for the MP
#' @author A. Hordyk
#' @export
ChkDatNA <- function(Data, dependencies) {
  slots <- trimws(gsub(",", "", unlist(strsplit(dependencies, "Data@"))))
  slots <- slots[nchar(slots) > 0]
  chk <- rep(FALSE, length(slots))
  for (x in seq_along(slots)) {
    chk[x] <- any(!is.finite(slot(Data, slots[x])))
  }
  if (any(chk)) stop("Missing data in slots -  ", paste(slots[chk], " "), call.=FALSE)
}


#' TAC Filter
#' 
#' Filters vector of TAC recommendations by replacing negatives with NA and
#' and values beyond five standard deviations from the mean as NA
#' 
#' @usage TACfilter(TAC)
#' @param TAC A numeric vector of TAC recommendations
#' @author T. Carruthers
#' @export TACfilter
TACfilter <- function(TAC) {
  TAC[TAC < 0] <- NA  # Have to robustify due to R optmization problems.. work in progress.
  TAC[TAC > (mean(TAC, na.rm = T) + 5 * stats::sd(TAC, na.rm = T))] <- NA  # remove very large TAC samples
  return(TAC)
}

prodPTF <- function(depletion, n, MSY) {
  # Pella-Tomlinson production function required for DB-SRA
  y <- (n^(n/(n - 1)))/(n - 1)
  MSY * y * depletion - MSY * y * depletion^n
}

fn <- function(n, BMSY_K) {
  # optimizer to find parameter n according to sampled BMSY/B0 (theta)
  thetapred <- n^(-1/(n - 1))
  (BMSY_K - thetapred)^2
}

getn <- function(BMSY_K) {
  # wrapper for n finder
  optimize(fn, c(0.01, 6), BMSY_K = BMSY_K)$minimum  #get the optimum
}

gety <- function(n) (n^(n/(n - 1)))/(n - 1)  # More DBSRA code: get the y parameter for n


# #' Average Catch
# #'
# #' A simple average catch MP that is included to demonstrate a 'status quo' management option
# #'
# #' @usage AvC(x, Data, reps = 100)
# #' @param x A position in a data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of stochastic samples of the TAC recommendation
# #' @author T. Carruthers
# #' @export AvC
# #'
# #' @importFrom abind abind
# #' @importFrom graphics abline axis barplot boxplot hist identify layout legend
# #' lines matplot mtext par plot plot.new points polygon segments text title text
# #' @importFrom grDevices col2rgb colorRampPalette rainbow rgb xy.coords
# #' @importFrom methods getClassDef .hasSlot new slot slot<- slotNames
# #' @importFrom stats approx coef dbeta density dnorm dlnorm lm loess loess.smooth
# #' median nlm optim optimise optimize plogis pnorm predict qlnorm quantile rbeta
# #' rlnorm rmultinom rnorm runif sd
# #' @importFrom utils packageVersion lsf.str read.csv
# AvC <- function(x, Data, reps = 100) {
#   dependencies = "Data@Cat"
#   rlnorm(reps, log(mean(Data@Cat[x, ], na.rm = T)), 0.2)
# }
# class(AvC) <- "Output"



# #' A reference FMSY method (uses perfect information about FMSY)
# #' 
# #' FMSY is taken from the operating model stored at DLM@OM$FMSY
# #' 
# #' Note that you can out-perform this MP even though it has perfect
# #' information of FMSY and current abundance. The requirment for fixed F is
# #' actually quite strict and is by no means the upper limit in terms of yield.
# #' Don't panic if your method beats this one for yield, especially for
# #' short-lived species of high temporal variability in productivity!
# #' 
# #' @usage FMSYref(x, Data, reps = 100)
# #' @param x A position in data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of TAC samples
# #' @author T. Carruthers
# #' @export FMSYref
# FMSYref <- function(x, Data, reps = 100) {
#   trlnorm(reps, Data@OM$A[x] * (1 - exp(-Data@OM$FMSY[x])), 0.01)
# }
# class(FMSYref) <- "Output"



#' A reference FMSY method that fishes at half of FMSY (uses perfect
#' information about FMSY)
#' 
#' FMSY is taken from the operating model stored at DLM@OM$FMSY
#' 
#' Note that you can out-performm this method easily. The requirement for fixed
#' F is actually quite strict and is by no means the upper limit in terms of
#' yield. Don't panic if your method beats this one for yield!
#' 
#' Interesting that the reduction in yield is no way near commensurate with the
#' reduction in F - as predicted by a yield curve and expressed in the pretty
#' good yield theory.
#' 
#' @usage FMSYref50(x, Data, reps = 100)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC (OFL) samples
#' @author T. Carruthers
#' @export FMSYref50
FMSYref50 <- function(x, Data, reps = 100) 
  trlnorm(reps, Data@OM$A[x] * (1 - exp(-Data@OM$FMSY[x]*0.5)) , 0.01)
class(FMSYref50) <- "Output"



#' A reference FMSY method that fishes at three quarters of FMSY (uses perfect
#' information about FMSY)
#' 
#' FMSY is taken from the operating model stored at DLM@OM$FMSY
#' 
#' Note that you can out-performm this method easily. The requirement for fixed
#' F is actually quite strict and is by no means the upper limit in terms of
#' yield. Don't panic if your method beats this one for yield!
#' 
#' Interesting that the reduction in yield is no way near commensurate with the
#' reduction in F as predicted by a yield curve and expressed in the pretty
#' good yield theory.
#' 
#' @usage FMSYref75(x, Data, reps = 100)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @author T. Carruthers
#' @export FMSYref75
FMSYref75 <- function(x, Data, reps = 100) 
  trlnorm(reps, Data@OM$A[x] * (1 - exp(-Data@OM$FMSY[x]*0.75)) , 0.01)
class(FMSYref75) <- "Output"



# #' Dynamic Fratio MP
# #' 
# #' The Fratio MP with a controller that changes the level of F according to the
# #' relationship between Surplus production and biomass. Ie lower F when dSP/dB
# #' is positive and higher F when dSP/dB is negative.
# #' 
# #' The method smoothes historical catches and biomass and then infers the
# #' relationship between surplus production and biomass (as suggested by Mark
# #' Maunder and Carl Walters). The approach then regulates a F based policy
# #' according to this gradient in which F may range between two different
# #' fractions of natural mortality rate.
# #' 
# #' The core advantage is the TAC(t) is not strongly determined by TAC(t-1) and
# #' therefore errors are not as readily propagated. The result is method that
# #' tends to perform alarmingly well and therefore requires debunking ASAP.
# #' 
# #' @usage DynF(x, Data, yrsmth=10, gg=2, reps = 100)
# #' @param x A position in a data-limited methods object
# #' @param Data A data-limited methods object
# #' @param yrsmth The number of historical recent years used for smoothing catch
# #' and biomass data
# #' @param gg A gain parameter that modifies F according to the gradient in
# #' surplus production with biomass
# #' @param reps The number samples of the TAC
# #' @return A numeric vector of TAC recommendations
# #' @author T. Carruthers
# #' @references Made-up for this package.
# #' @export DynF
# DynF <- function(x, Data, yrsmth = 10, gg = 2, reps = 100) {
#   
#   dependencies = "Data@Year, Data@Cat, Data@Ind, Data@Abun, Data@Mort, Data@FMSY_M"
#   ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
#   
#   C_dat <- log(Data@Cat[x, ind])
#   C_dat[C_dat == -Inf] <- 0
#   B_dat <- log(Data@Ind[x, ind]/Data@Ind[x, ind[yrsmth]] * Data@Abun[x])
#   B_dat[B_dat == -Inf] <- 0
#   C_hist <- exp(predict(loess(C_dat ~ ind, degree = 1)))
#   B_hist <- exp(predict(loess(B_dat ~ ind, degree = 1)))
#   
#   ind <- 2:yrsmth
#   ind1 <- 1:(yrsmth - 1)
#   SP_hist <- B_hist[ind] - B_hist[ind1] + C_hist[ind1]
#   
#   Frat <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x]) * trlnorm(reps, 
#     Data@FMSY_M[x], Data@CV_FMSY_M[x])
#   Flim <- matrix(NA, nrow = 2, ncol = reps)
#   Flim[1, ] <- Frat * 0.5
#   Flim[2, ] <- Frat * 2
#   
#   yind <- 1:length(SP_hist)
#   SP_mu <- predict(lm(SP_hist ~ yind), newdat = list(yind = length(SP_hist) + 
#     1))
#   SP_se <- predict(lm(SP_hist ~ yind), newdat = list(yind = length(SP_hist) + 
#     1), se = T)$se.fit
#   SP_new <- rnorm(reps, SP_mu, SP_se/2)
#   Glm <- summary(lm(SP_hist ~ B_hist[ind1]))$coefficients[2, 1:2]  # plot(B_hist[ind1],SP_hist) # points(B_hist[ind1],SP_hist,col='green')
#   G_new <- rnorm(reps, Glm[1], Glm[2]/2)
#   # G_new[G_new>2*Frat]<-2*Frat[G_new<(2*Frat)]
#   # G_new[G_new<(-2*Frat)]<--2*Frat[G_new<(-2*Frat)]
#   G_new[G_new > 0] <- G_new[G_new > 0] * 3
#   newF <- Frat * exp(-G_new * gg)
#   newF[newF < Flim[1]] <- Flim[1]
#   newF[newF > Flim[2]] <- Flim[2]
#   
#   TAC <- newF * B_hist[yrsmth]
#   TACfilter(TAC)
#   
# }
# class(DynF) <- "Output"



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
#' management procedures. Fish and Fisheries. In press.  Maunder. 2014.
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
  TACfilter(TAC)
}
class(Fadapt) <- "Output"



#' Depletion Corrected Fratio
#' 
#' The Fratio MP with a harvest control rule that reduces F according to the
#' production curve given an estimate of current stock depletion.
#' 
#' 
#' @usage DepF(x, Data, reps = 100)
#' @param x A position in data-limited methods data object DLM
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @return A numeric vector of TAC recommendations
#' @author T. Carruthers
#' @references Made-up for this package.
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
  TACfilter(TAC)
}
class(DepF) <- "Output"



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
  TACfilter(TAC)
}
class(Gcontrol) <- "Output"




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
  
  Bt_K <- rbeta(100, alphaconv(depo, min(depo * Data@CV_Dep[x], (1 - 
    depo) * Data@CV_Dep[x])), betaconv(depo, min(depo * Data@CV_Dep[x], 
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
  TACfilter(TAC)
}
class(Rcontrol) <- "Output"



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
#' @usage Rcontrol2(x, Data, reps = 100, yrsmth = 10, gg = 2, glim = c(0.5,
#' 2))
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
  TACfilter(TAC)
}
class(Rcontrol2) <- "Output"



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
  TACfilter(TAC)
}
class(GB_CC) <- "Output"



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
  TACfilter(TAC)
}
class(GB_slope) <- "Output"



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
  TACfilter(TAC)
}
class(GB_target) <- "Output"



# #' Constant catch management procedure of Geromont and Butterworth (2014)
# #' 
# #' The TAC is the average catch over last yrsmth years.
# #' 
# #' This is one of four constant catch rules of Geromont and Butterworth 2014.
# #' 
# #' @usage CC1(x, Data, reps = 100, yrsmth = 5, xx=0)
# #' @param x A position in data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of TAC samples
# #' @param yrsmth Years over which to calculate mean catches
# #' @param xx Parameter controlling the TAC. Mean catches are multiplied by
# #' (1-xx)
# #' @return A numeric vector of TAC recommendations
# #' @author T. Carruthers
# #' @references Geromont, H.F., Butterworth, D.S. 2014. Generic management
# #' procedures for data-poor fisheries; forecasting with few data. ICES J. Mar.
# #' Sci. doi:10.1093/icesjms/fst232
# #' @export CC1
# CC1 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0) {
#   dependencies = "Data@Cat, Data@CV_Cat"
#   C_dat <- Data@Cat[x, (length(Data@Year) - (yrsmth - 1)):length(Data@Year)]
#   TAC <- (1 - xx) * trlnorm(reps, mean(C_dat), Data@CV_Cat/(yrsmth^0.5))  # mean catches over the interval
#   TACfilter(TAC)
# }
# class(CC1) <- "Output"



# #' Constant catch management procedure of Geromont and Butterworth (2014)
# #' 
# #' The TAC is the average catch over last yrsmth years reduced by 30%.
# #' 
# #' This is one of four constant catch MPs of Geromont and Butterworth 2014.
# #' 
# #' @usage CC4(x, Data, reps = 100, yrsmth = 5, xx=0.3)
# #' @param x A position in data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of TAC samples
# #' @param yrsmth Years over which to average catches
# #' @param xx Parameter controlling the TAC. Mean catches are multiplied by
# #' (1-xx)
# #' @return A numeric vector of TAC recommendations
# #' @author T. Carruthers
# #' @references Geromont, H.F., Butterworth, D.S. 2014. Generic management
# #' procedures for data-poor fisheries; forecasting with few data. ICES J. Mar.
# #' Sci. doi:10.1093/icesjms/fst232
# #' @export CC4
# CC4 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0.3) {
#   dependencies = "Data@Cat, Data@CV_Cat"
#   C_dat <- Data@Cat[x, (length(Data@Year) - (yrsmth - 1)):length(Data@Year)]
#   TAC <- (1 - xx) * trlnorm(reps, mean(C_dat), Data@CV_Cat/(yrsmth^0.5))  # mean catches over the interval
#   TACfilter(TAC)
# }
# class(CC4) <- "Output"




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
  TACfilter(TAC)
}
class(LstepCC1) <- "Output"




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
  TACfilter(TAC)
}
class(LstepCC4) <- "Output"




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
  TACfilter(TAC)
}
class(Ltarget1) <- "Output"



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
  TACfilter(TAC)
}
class(Ltarget4) <- "Output"




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
  TACfilter(TAC)
}
class(Itarget1) <- "Output"



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
  TACfilter(TAC)
}
class(Itarget4) <- "Output"




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
 
  TACfilter(TAC)
  
}
class(Islope1) <- "Output"



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
  TACfilter(TAC)
}
class(Islope4) <- "Output"





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
  
  dependencies = "Data@Ind, Data@Cat, Data@CV_Ind, Data@Iref"
  ind <- max(1, (length(Data@Year) - yrsmth + 1)):length(Data@Year)
  
  
  
  deltaI <- mean(Data@Ind[x, ind])/Data@Iref[x]
  if (deltaI < (1 - mc)) 
    deltaI <- 1 - mc
  if (deltaI > (1 + mc)) 
    deltaI <- 1 + mc
  
  TAC <- Data@MPrec[x] * deltaI * trlnorm(reps, 1, Data@CV_Ind[x])
  TAC
}
class(IT10) <- "Output"



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
  
  dependencies = "Data@Ind, Data@Cat, Data@CV_Ind, Data@Iref"
  ind <- max(1, (length(Data@Year) - yrsmth + 1)):length(Data@Year)
  deltaI <- mean(Data@Ind[x, ind])/Data@Iref[x]
  if (deltaI < (1 - mc)) 
    deltaI <- 1 - mc
  if (deltaI > (1 + mc)) 
    deltaI <- 1 + mc
  
  TAC <- Data@MPrec[x] * deltaI * trlnorm(reps, 1, Data@CV_Ind[x])
  TAC
}
class(IT5) <- "Output"



#' Index Target based on natural mortality rate
#' 
#' An index target MP where the TAC is modified according to current index
#' levels (mean index over last yrsmth years) relative to a target level.
#' Maximum fractional annual changes are mc.  mc=(5+M*25)/100
#' yrsmth=4*(1/M)^(0.25)
#' 
#' 
#' @usage ITM(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the quota recommendation
#' @return A numeric vector of TAC recommendations
#' @author T. Carruthers
#' @export ITM
ITM <- function(x, Data, reps = 100) {
  
  dependencies = "Data@Ind, Data@Cat, Data@CV_Ind, Data@Iref, Data@Mort"
  mc <- (5 + Data@Mort[x] * 25)/100
  if (mc > 0.2) 
    mc <- 0.2
  yrsmth <- floor(4 * (1/Data@Mort[x])^(1/4))
  ind <- max(1, (length(Data@Year) - yrsmth + 1)):length(Data@Year)
  
  deltaI <- mean(Data@Ind[x, ind])/Data@Iref[x]
  if (deltaI < (1 - mc)) 
    deltaI <- 1 - mc
  if (deltaI > (1 + mc)) 
    deltaI <- 1 + mc
  
  TAC <- Data@MPrec[x] * deltaI * trlnorm(reps, 1, Data@CV_Ind[x])
  TAC
}
class(ITM) <- "Output"



#' Surplus production based catch-limit modifier
#' 
#' An MP that makes incremental adjustments to TAC recommendations based on the
#' apparent trend in surplus production. Based on the theory of Mark Maunder
#' (IATTC)
#' 
#' Note that this isn't exactly what Mark has previously suggested and is
#' stochastic in this implementation.
#' 
#' @usage SPmod(x, Data, reps = 100, alp = c(0.8, 1.2), bet = c(0.8, 1.2))
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of quota samples
#' @param alp Condition for modifying the TAC (bounds on change in abundance)
#' @param bet Limits for how much the TAC can change among years
#' @return A numeric vector of TAC recommendations
#' @author T. Carruthers
#' @references
#' http://www.iattc.org/Meetings/Meetings2014/MAYSAC/PDFs/SAC-05-10b-Management-Strategy-Evaluation.pdf
#' @export SPmod
SPmod <- function(x, Data, reps = 100, alp = c(0.8, 1.2), bet = c(0.8, 
  1.2)) {
  dependencies = "Data@Cat, Data@Ind, Data@Abun"
  Ir <- length(Data@Ind[x, ])
  Cr <- length(Data@Cat[x, ])
  rat <- trlnorm(reps, Data@Ind[x, Ir], Data@CV_Ind)/trlnorm(reps, 
    Data@Ind[x, Ir - 1], Data@CV_Ind)
  cct <- trlnorm(reps, Data@Cat[x, Cr], Data@CV_Cat)
  Abun <- trlnorm(reps, Data@Abun[x], Data@CV_Abun)
  TAC <- rep(NA, reps)
  TAC[rat < alp[1]] <- cct[rat < alp[1]] * bet[1]
  TAC[rat > alp[1] & rat < alp[2]] <- cct[rat > alp[1] & rat < alp[2]]
  
  cond <- rat > alp[2]
  reps2 <- sum(cond)
  if (reps2 > 0) {
    qq1 <- trlnorm(reps2, Data@Ind[x, Ir]/Abun, Data@CV_Ind)
    bio1 <- Data@Ind[x, Ir - 1]/qq1
    bio2 <- Data@Ind[x, Ir]/qq1
    cct1 <- trlnorm(reps2, Data@Cat[x, Cr - 1], Data@CV_Cat)
    PP <- bio2 - bio1 + cct1
    TAC[cond] <- bet[2] * PP
  }
  TACfilter(TAC)
}
class(SPmod) <- "Output"



#' Slope in surplus production MP
#' 
#' A management procedure that makes incremental adjustments to TAC
#' recommendations based on the apparent trend in recent surplus production.
#' Based on the theory of Mark Maunder (IATTC)
#' 
#' Note that this isn't exactly what Mark has previously suggested and is
#' stochastic in this implementation.
#' 
#' @usage SPslope(x, Data, reps = 100, yrsmth = 4, alp = c(0.9, 1.1), bet =
#' c(1.5, 0.9))
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
#' @export SPslope
SPslope <- function(x, Data, reps = 100, yrsmth = 4, alp = c(0.9, 1.1), 
  bet = c(1.5, 0.9)) {
  
  dependencies = "Data@Year, Data@Cat, Data@Ind, Data@Abun"
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
  yind <- 1:yrsmth
  C_dat <- Data@Cat[x, ind]
  B_dat <- Data@Ind[x, ind]/Data@Ind[x, ind[yrsmth]] * Data@Abun[x]
  Pt_mu <- max(B_dat[yrsmth] - B_dat[yrsmth - 1] + C_dat[yrsmth - 1], 
    tiny)
  Pt_1 <- trlnorm(reps, Pt_mu, Data@CV_Cat)
  It <- exp(predict(lm(log(B_dat) ~ yind), newdat = list(yind = yrsmth + 
    1)))
  Ilast <- B_dat[yrsmth]
  MC <- max(mean(C_dat), tiny)
  Ct_1 <- trlnorm(reps, MC, Data@CV_Cat/(yrsmth^0.5))  # mean catches over the interval
  
  rat <- It/Ilast
  
  mult <- max((1 - bet[1] * (Ilast - It)/Ilast), tiny)
  if (rat < alp[1]) 
    TAC <- mult * Ct_1
  if (rat > alp[1] & rat < alp[2]) 
    TAC <- Ct_1
  if (rat > alp[2]) 
    TAC <- bet[2] * Pt_1
  TACfilter(TAC)
}
class(SPslope) <- "Output"



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
  TACfilter(TAC)
}
class(SBT1) <- "Output"



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
  TACfilter(TAC)
}
class(SBT2) <- "Output"



# #' Delay - Difference Stock Assessment with UMSY and MSY leading
# #' 
# #' A simple delay-difference assessment that estimates the TAC using a
# #' time-series of catches and a relative abundance index.
# #' 
# #' 
# #' @usage DD(x, Data, reps = 100)
# #' @param x A position in a data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of stochastic samples of the TAC recommendation
# #' @return A numeric vector of TAC recommendations
# #' @note This DD model is observation error only and has does not estimate
# #' process error (recruitment deviations). Similar to many other assessment
# #' models it depends on a whole host of dubious assumptions such as temporally
# #' stationary productivity and proportionality between the abundance index and
# #' real abundance. Unsurprisingly the extent to which these assumptions are
# #' violated tends to be the biggest driver of performance for this method.
# #' @author T. Carruthers
# #' @references Method based on equations of Carl Walters (bug him with
# #' questions and expect colourful responses)
# #' @export DD
# DD <- function(x, Data, reps = 100) {
#   # for(x in 1:nsim){
#   dependencies = "Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@vbt0, Data@CV_vbt0, Data@Mort, Data@CV_Mort, Data@wla, Data@wlb"
#   Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
#   Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
#   if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
#     t0c <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
#   } else {
#     t0c <- rep(Data@vbt0[x], reps)
#   }
#   t0c[!is.finite(t0c)] <- 0
#   Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
#   a <- Data@wla[x]
#   b <- Data@wlb[x]
#   
#   Winf = Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
#   age <- 1:Data@MaxAge
#   la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
#   wa <- Data@wla[x] * la^Data@wlb[x]
#   a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],  Data@L50[x])
#   a50V <- max(a50V, 1)
#   yind <- (1:length(Data@Cat[x, ]))[!is.na(Data@Cat[x, ] + Data@Ind[x,   ])]
#   C_hist <- Data@Cat[x, yind]
#   E_hist <- C_hist/Data@Ind[x, yind]
#   E_hist <- E_hist/mean(E_hist)
#   ny_DD <- length(C_hist)
#   params <- log(c(Data@Mort[x], mean(C_hist, na.rm = T), Data@Mort[x]))
#   k_DD <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)  
#   k_DD[k_DD > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
#   Rho_DD <- (wa[k_DD + 2] - Winf)/(wa[k_DD + 1] - Winf)
#   Alpha_DD <- Winf * (1 - Rho_DD)
#   So_DD <- exp(-Data@Mort[x])  # get So survival rate
#   wa_DD <- wa[k_DD]
#   UMSYprior <- c(1 - exp(-Data@Mort[x] * 0.5), 0.3)
#   opt <- optim(params, DD_R, opty = 1, So_DD = So_DD, Alpha_DD = Alpha_DD, 
#     Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, wa_DD = wa_DD, E_hist = E_hist, 
#     C_hist = C_hist, UMSYprior = UMSYprior, method = "L-BFGS-B", lower = log(exp(params)/20), 
#     upper = log(exp(params) * 20), hessian = TRUE)
#   
#   # Catfit<-DD_R(opt$par,opty=3,So_DD=So_DD,Alpha_DD=Alpha_DD,Rho_DD=Rho_DD,ny_DD=ny_DD,k_DD=k_DD,wa_DD=wa_DD,E_hist=E_hist,C_hist=C_hist,UMSYprior=UMSYprior)
#   # plot(Catfit[,1],ylim=c(0,max(Catfit))) lines(Catfit[,2],col='red')
#   
#   TAC <- rep(NA, reps)
#   # samps<-rmvnorm(reps,opt$par,solve(opt$hessian)) # assuming log
#   # parameters are multivariate normal hessian approximation
#   samps <- cbind(rnorm(reps, opt$par[1], ((opt$par[1])^2)^0.5 * 0.1), 
#     rnorm(reps, opt$par[2], ((opt$par[2])^2)^0.5 * 0.1), rnorm(reps, 
#       opt$par[3], ((opt$par[3])^2)^0.5 * 0.1))
#   if (reps == 1) 
#     samps <- matrix(c(opt$par[1], opt$par[2], opt$par[3]), nrow = 1)
#   for (i in 1:reps) TAC[i] <- DD_R(samps[i, ], opty = 2, So_DD = So_DD, 
#     Alpha_DD = Alpha_DD, Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, 
#     wa_DD = wa_DD, E_hist = E_hist, C_hist = C_hist, UMSYprior = UMSYprior)
#   TACfilter(TAC)
# }
# class(DD) <- "Output"




# #' Delay - Difference Stock Assessment with UMSY and MSY leading coupled with a
# #' 40-10 harvest control rule
# #' 
# #' A simple delay-difference assessment that estimates the OFL using a
# #' time-series of catches and a relative abundance index. In this version of
# #' the DD MP a 40-10 rule is imposed over the OFL recommendation.
# #' 
# #' 
# #' @usage DD4010(x, Data, reps = 100)
# #' @param x A position in a data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of stochastic samples of the TAC recommendation
# #' @return A numeric vector of TAC recommendations
# #' @author T. Carruthers
# #' @references Method based on equations of Carl Walters
# #' @export DD4010
# DD4010 <- function(x, Data, reps = 100) {
#   dependencies = "Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@vbt0, Data@CV_vbt0, Data@Mort, Data@CV_Mort. Data@wla, Data@ wlb"
#   Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
#   Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
#   if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
#     t0c <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
#   } else {
#     t0c <- rep(Data@vbt0[x], reps)
#   }
#   t0c[!is.finite(t0c)] <- 0
#   Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
#   a <- Data@wla[x]
#   b <- Data@wlb[x]
#   
#   Winf = Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
#   age <- 1:Data@MaxAge
#   la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
#   wa <- Data@wla[x] * la^Data@wlb[x]
#   a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], 
#     Data@L50[x])
#   a50V <- max(a50V, 1)
#   yind <- (1:length(Data@Cat[x, ]))[!is.na(Data@Cat[x, ] + Data@Ind[x, 
#     ])]
#   C_hist <- Data@Cat[x, yind]
#   E_hist <- Data@Ind[x, yind]
#   E_hist <- C_hist/E_hist
#   E_hist <- E_hist/mean(E_hist)
#   ny_DD <- length(C_hist)
#   params <- log(c(Data@Mort[x], mean(C_hist, na.rm = T), Data@Mort[x]))
#   k_DD <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)  
#   k_DD[k_DD > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
#   Rho_DD <- (wa[k_DD + 2] - Winf)/(wa[k_DD + 1] - Winf)
#   Alpha_DD <- Winf * (1 - Rho_DD)
#   So_DD <- exp(-Data@Mort[x])  # get So survival rate
#   wa_DD <- wa[k_DD]
#   UMSYprior <- c(1 - exp(-Data@Mort * 0.5), 0.3)
#   opt <- optim(params, DD_R, opty = 1, So_DD = So_DD, Alpha_DD = Alpha_DD, 
#     Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, wa_DD = wa_DD, E_hist = E_hist, 
#     C_hist = C_hist, UMSYprior = UMSYprior, method = "L-BFGS-B", lower = log(exp(params)/20), 
#     upper = log(exp(params) * 20), hessian = TRUE)
#   
#   # Catfit<-DD_R(opt$par,opty=3,So_DD=So_DD,Alpha_DD=Alpha_DD,Rho_DD=Rho_DD,ny_DD=ny_DD,k_DD=k_DD,wa_DD=wa_DD,E_hist=E_hist,C_hist=C_hist,UMSYprior=UMSYprior)
#   # plot(Catfit[,1],ylim=c(0,max(Catfit))) lines(Catfit[,2],col='red')
#   
#   TAC <- rep(NA, reps)
#   dep <- rep(NA, reps)
#   # samps<-rmvnorm(reps,opt$par,solve(opt$hessian)) # assuming log
#   # parameters are multivariate normal hessian approximation
#   samps <- cbind(rnorm(reps, opt$par[1], ((opt$par[1])^2)^0.5 * 0.1), 
#     rnorm(reps, opt$par[2], ((opt$par[2])^2)^0.5 * 0.1), rnorm(reps, 
#       opt$par[3], ((opt$par[3])^2)^0.5 * 0.1))
#   if (reps == 1) 
#     samps <- matrix(c(opt$par[1], opt$par[2], opt$par[3]), nrow = 1)
#   for (i in 1:reps) TAC[i] <- DD_R(samps[i, ], opty = 2, So_DD = So_DD, 
#     Alpha_DD = Alpha_DD, Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, 
#     wa_DD = wa_DD, E_hist = E_hist, C_hist = C_hist, UMSYprior = UMSYprior)
#   for (i in 1:reps) dep[i] <- DD_R(samps[i, ], opty = 3, So_DD = So_DD, 
#     Alpha_DD = Alpha_DD, Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, 
#     wa_DD = wa_DD, E_hist = E_hist, C_hist = C_hist, UMSYprior = UMSYprior)
#   cond1 <- !is.na(dep) & dep < 0.4 & dep > 0.1
#   cond2 <- !is.na(dep) & dep < 0.1
#   TAC[cond1] <- TAC[cond1] * (dep[cond1] - 0.1)/0.3
#   TAC[cond2] <- TAC[cond2] * tiny  # this has to still be stochastic albeit very small
#   TACfilter(TAC)
# }
# class(DD4010) <- "Output"




#' Internal optimization function for delay-difference MPs
#'
#' @param params vector of length 3 with log parameters
#' @param opty optimization option
#' @param So_DD internal parameter
#' @param Alpha_DD  internal parameter
#' @param Rho_DD  internal parameter
#' @param ny_DD  internal parameter
#' @param k_DD  internal parameter
#' @param wa_DD  internal parameter
#' @param E_hist  internal parameter
#' @param C_hist  internal parameter
#' @param UMSYprior  internal parameter
#'
#' @author T. Carruthers
#' @keywords internal
#' @export
#'
DD_R <- function(params, opty, So_DD, Alpha_DD, Rho_DD, ny_DD, k_DD, wa_DD, E_hist, 
                 C_hist, UMSYprior) {
  UMSY_DD = exp(params[1])
  MSY_DD = exp(params[2])
  q_DD = exp(params[3])
  SS_DD = So_DD * (1 - UMSY_DD)  # Initialise for UMSY, MSY and q leading.
  Spr_DD = (SS_DD * Alpha_DD/(1 - SS_DD) + wa_DD)/(1 - Rho_DD * SS_DD)
  DsprDu_DD = -So_DD * (Rho_DD/(1 - Rho_DD * SS_DD) * Spr_DD + 1/(1 - 
    Rho_DD * SS_DD) * (Alpha_DD/(1 - SS_DD) + SS_DD * Alpha_DD/(1 - 
    SS_DD)^2))
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
  Cpred_DD[Cpred_DD < tiny] <- tiny
  
  if (opty == 1) {
    test <- dnorm(log(Cpred_DD), log(C_hist), 0.25, log = T)
    test2 <- dlnorm(UMSY_DD, log(UMSYprior[1]), UMSYprior[2], log = T)
    test[is.na(test)] <- -1000
    test[test == (-Inf)] <- -1000
    if (is.na(test2) | test2 == -Inf | test2 == Inf) 
      test2 <- 1000
    return(-sum(test, test2))  # return objective function
  } else if (opty == 2) {
    # return MLE TAC estimate
    UMSY_DD * B_DD[ny_DD]
  } else if (opty == 3) {
    B_DD[tt + 1]/Bo_DD
  } else {
    cbind(C_hist, Cpred_DD)  # return observations vs predictions
  }
}


# 
# #' Depletion-Based Stock Reduction Analysis
# #' 
# #' User prescribed BMSY/B0, M, FMSY/M are used to find B0 and therefore the OFL
# #' by back-constructing the stock to match a user specified level of stock
# #' depletion (OFL = M * FMSY/M * depletion* B0).
# #' 
# #' You specify a range of stock depletion and, given historical catches DB-SRA
# #' calculates what unfished biomass must have been to get you here given
# #' samples for M, FMSY relative to M and also BMSY relative to Bunfished.
# #' 
# #' @usage DBSRA(x, Data, reps = 100)
# #' @param x A position in a data-limited methods object.
# #' @param Data A data-limited methods object.
# #' @param reps The number of samples of the TAC (OFL) recommendation.
# #' @return A vector of TAC (OFL) values.
# #' @note This is set up to return the OFL (FMSY * current biomass).
# #' 
# #' You may have noticed that you -the user- specify three of the factors that
# #' make the quota recommendation. So this can be quite a subjective method.
# #' 
# #' Also the DB-SRA method of this package isn't exactly the same as the
# #' original method of Dick and MacCall (2011) because it has to work for
# #' simulated depletions above BMSY/B0 and even on occasion over B0. Also it
# #' doesn't have the modification for flatfish life histories that has
# #' previously been applied by Dick and MacCall.
# #' @author T. Carruthers
# #' @references Dick, E.J., MacCall, A.D., 2011. Depletion-Based Stock Reduction
# #' Analysis: A catch-based method for determining sustainable yields for
# #' data-poor fish stocks. Fish. Res. 110, 331-341.
# #' @export DBSRA
# DBSRA <- function(x, Data, reps = 100) {
#   # returns a vector of DBSRA estimates of the TAC for a particular
#   # simulation x for(x in 1:nsim){
#   dependencies = "Data@Cat, Data@Dep, Data@CV_Dep, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M,Data@BMSY_B0, Data@CV_BMSY_B0, Data@L50"
#   C_hist <- Data@Cat[x, ]
#   TAC <- rep(NA, reps)
#   DBSRAcount <- 1
#   if (is.na(Data@Dep[x]) | is.na(Data@CV_Dep[x])) return(NA)
#   while (DBSRAcount < (reps + 1)) {
#     depo <- max(0.01, min(0.99, Data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
#     Bt_K <- rbeta(100, alphaconv(depo, min(depo * Data@CV_Dep[x], (1 - depo) * Data@CV_Dep[x])), 
#                   betaconv(depo, min(depo * Data@CV_Dep[x], (1 - depo) * Data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
#     Bt_K <- Bt_K[Bt_K > 0.00999 & Bt_K < 0.99001][1]  # interval censor (0.01,0.99)  as in Dick and MacCall 2011
#     Mdb <- trlnorm(100, Data@Mort[x], Data@CV_Mort[x])
#     Mdb <- Mdb[Mdb < 0.9][1]  # !!!! maximum M is 0.9   interval censor
#     if (is.na(Mdb)) Mdb <- 0.9  # !!!! maximum M is 0.9   absolute limit
#     FMSY_M <- trlnorm(1, Data@FMSY_M[x], Data@CV_FMSY_M[x])
#     BMSY_K <- rbeta(100, alphaconv(Data@BMSY_B0[x], Data@CV_BMSY_B0[x] *  Data@BMSY_B0[x]), 
#                     betaconv(Data@BMSY_B0[x], Data@CV_BMSY_B0[x] * Data@BMSY_B0[x]))  #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
#     tryBMSY_K <- BMSY_K[BMSY_K > 0.05 & BMSY_K < 0.95][1]  # interval censor (0.05,0.95) as in Dick and MacCall, 2011
#     if (is.na(tryBMSY_K)) {
#       Min <- min(BMSY_K, na.rm = TRUE)
#       Max <- max(BMSY_K, na.rm = TRUE)
#       if (Max <= 0.05) 
#         BMSY_K <- 0.05
#       if (Min >= 0.95) 
#         BMSY_K <- 0.95
#     }
#     if (!is.na(tryBMSY_K))  BMSY_K <- tryBMSY_K
#     
#     adelay <- max(floor(iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],  Data@L50[x])), 1)
#     opt <- optimize(DBSRAopt, log(c(0.01 * mean(C_hist), 1000 * mean(C_hist))), C_hist = C_hist, 
#                     nys = length(C_hist), Mdb = Mdb, FMSY_M = FMSY_M, BMSY_K = BMSY_K, 
#                     Bt_K = Bt_K, adelay = adelay, tol = 0.01)
#     # if(opt$objective<0.1){
#     Kc <- exp(opt$minimum)
#     BMSYc <- Kc * BMSY_K
#     FMSYc <- Mdb * FMSY_M
#     UMSYc <- (FMSYc/(FMSYc + Mdb)) * (1 - exp(-(FMSYc + Mdb)))
#     MSYc <- Kc * BMSY_K * UMSYc
#     TAC[DBSRAcount] <- UMSYc * Kc * Bt_K
#     DBSRAcount <- DBSRAcount + 1
#     # }
#     
#     # print(DBSRAcount)
#   }  # end of reps
#   TACfilter(TAC)
#   # print(x) }
# }  # end of DBSRA_apply
# class(DBSRA) <- "Output"



# #' Depletion-Based Stock Reduction Analysis assuming 40 per cent stock
# #' depletion
# #' 
# #' DBSRA assuming that current stock depletion is exactly 40 per cent of
# #' unfished stock levels.
# #' 
# #' 
# #' @usage DBSRA_40(x, Data, reps = 100)
# #' @param x A position in a data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of stochastic samples of the TAC recommendation
# #' @note A 40 percent assumption for current depletion is more or less the most
# #' optimistic state for a stock (ie very close to BMSY/B0 for many stocks).
# #' @author T. Carruthers.
# #' @references Dick, E.J., MacCall, A.D., 2010. Estimates of sustainable yield
# #' for 50 data-poor stocks in the Pacific Coast groundfish fishery management
# #' plan. Technical memorandum. Southwest fisheries Science Centre, Santa Cruz,
# #' CA. National Marine Fisheries Service, National Oceanic and Atmospheric
# #' Administration of the U.S. Department of Commerce. NOAA-TM-NMFS-SWFSC-460.
# #' @export DBSRA_40
# DBSRA_40 <- function(x, Data, reps = 100) {
#   # returns a vector of DBSRA estimates of the TAC for a particular
#   # simulation x
#   dependencies = "Data@Cat, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@BMSY_B0, Data@CV_BMSY_B0, Data@L50"
#   C_hist <- Data@Cat[x, ]
#   TAC <- rep(NA, reps)
#   DBSRAcount <- 1
#   if (is.na(Data@Dep[x]) | is.na(Data@CV_Dep[x]))   return(NA)
#   while (DBSRAcount < (reps + 1)) {
#     depo <- 0.4
#     Bt_K <- rbeta(100, alphaconv(depo, min(depo * Data@CV_Dep[x], 
#       (1 - depo) * Data@CV_Dep[x])), betaconv(depo, min(depo * 
#       Data@CV_Dep[x], (1 - depo) * Data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
#     Bt_K <- Bt_K[Bt_K > 0.00999 & Bt_K < 0.99001][1]  # interval censor (0.01,0.99)  as in Dick and MacCall 2011
#     Mdb <- stats::rlnorm(100, mconv(Data@Mort[x], Data@CV_Mort[x] * 
#       Data@Mort[x]), sdconv(Data@Mort[x], Data@CV_Mort[x] * 
#       Data@Mort[x]))  # log space stdev 0.4 as in Dick and MacCall 2011
#     Mdb <- Mdb[Mdb < 0.9][1]  # !!!! maximum M is 0.9   interval censor
#     if (is.na(Mdb)) 
#       Mdb <- 0.9  # !!!! maximum M is 0.9   absolute limit
#     FMSY_M <- trlnorm(1, Data@FMSY_M[x], Data@CV_FMSY_M[x])
#     BMSY_K <- rbeta(100, alphaconv(Data@BMSY_B0[x], Data@CV_BMSY_B0[x] * 
#       Data@BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@CV_BMSY_B0[x] * 
#       Data@BMSY_B0[x]))  #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
#     tryBMSY_K <- BMSY_K[BMSY_K > 0.05 & BMSY_K < 0.95][1]  # interval censor (0.05,0.95) as in Dick and MacCall, 2011
#     if (is.na(tryBMSY_K)) {
#       Min <- min(BMSY_K, na.rm = TRUE)
#       Max <- max(BMSY_K, na.rm = TRUE)
#       if (Max <= 0.05) 
#         BMSY_K <- 0.05
#       if (Min >= 0.95) 
#         BMSY_K <- 0.95
#     }
#     if (!is.na(tryBMSY_K))  BMSY_K <- tryBMSY_K
#     adelay <- max(floor(iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],  Data@L50[x])), 1)
#     opt <- optimize(DBSRAopt, log(c(0.1 * mean(C_hist), 1000 * mean(C_hist))), 
#       C_hist = C_hist, nys = length(C_hist), Mdb = Mdb, FMSY_M = FMSY_M, 
#       BMSY_K = BMSY_K, Bt_K = Bt_K, adelay = adelay, tol = 0.01)
#     # if(opt$objective<0.1){
#     Kc <- exp(opt$minimum)
#     BMSYc <- Kc * BMSY_K
#     FMSYc <- Mdb * FMSY_M
#     UMSYc <- (FMSYc/(FMSYc + Mdb)) * (1 - exp(-(FMSYc + Mdb)))
#     MSYc <- Kc * BMSY_K * UMSYc
#     TAC[DBSRAcount] <- UMSYc * Kc * Bt_K
#     DBSRAcount <- DBSRAcount + 1
#     # }
#   }  # end of reps
#   TACfilter(TAC)
# }  # end of DBSRA_apply
# class(DBSRA_40) <- "Output"



# #' Depletion-Based Stock Reduction Analysis using mean length estimator of
# #' stock depletion
# #' 
# #' DBSRA using the mean length estimator to calculate current stock depletion.
# #' 
# #' 
# #' @usage DBSRA_ML(x, Data, reps = 100)
# #' @param x A position in a data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of stochastic samples of the quota recommendation
# #' @note The mean length extension was programmed by Gary Nelson as part of his
# #' excellent R package 'fishmethods'
# #' @author T. Carruthers
# #' @references Dick, E.J., MacCall, A.D., 2011. Depletion-Based Stock Reduction
# #' Analysis: A catch-based method for determining sustainable yields for
# #' data-poor fish stocks. Fish. Res. 110, 331-341.
# #' @export DBSRA_ML
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
#     if (all(is.na(Z))) return(rep(NA, reps))
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
#       Data@BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@CV_BMSY_B0[x] * 
#       Data@BMSY_B0[x]))  #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
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
#       C_hist = C_hist, nys = length(C_hist), Mdb = Mdb, FMSY_M = FMSY_M, 
#       BMSY_K = BMSY_K, Bt_K = Bt_K, adelay = adelay, tol = 0.01)
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
#   TACfilter(TAC)
# }
# class(DBSRA_ML) <- "Output"



# #' Depletion-Based Stock Reduction Analysis paired with 40-10 harvest control
# #' rule
# #' 
# #' User prescribed BMSY/B0, M, FMSY/M are used to find B0 and therefore the OFL
# #' by back-constructing the stock to match a user specified level of stock
# #' depletion (OFL = M * FMSY/M * depletion* B0). In this method DBSRA is paried
# #' with the 40-10 rule that throttles back the OFL to zero at 10 percent of
# #' unfished biomass.
# #' 
# #' 
# #' @usage DBSRA4010(x, Data, reps = 100)
# #' @param x A position in a data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of stochastic samples of the TAC recommendation
# #' @author T. Carruthers
# #' @references Dick, E.J., MacCall, A.D., 2011. Depletion-Based Stock Reduction
# #' Analysis: A catch-based method for determining sustainable yields for
# #' data-poor fish stocks. Fish. Res. 110, 331-341.
# #' @export DBSRA4010
# DBSRA4010 <- function(x, Data, reps = 100) {
#   # returns a vector of DBSRA estimates of the TAC for a particular
#   # simulation x for(x in 1:nsim){
#   dependencies = "Data@Cat, Data@Dep, Data@CV_Dep, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M,Data@BMSY_B0, Data@CV_BMSY_B0, Data@L50"
#   C_hist <- Data@Cat[x, ]
#   TAC <- rep(NA, reps)
#   DBSRAcount <- 1
#   if (is.na(Data@Dep[x]) | is.na(Data@CV_Dep[x])) 
#     return(NA)
#   while (DBSRAcount < (reps + 1)) {
#     depo <- max(0.01, min(0.99, Data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
#     Bt_K <- rbeta(100, alphaconv(depo, min(depo * Data@CV_Dep[x], 
#       (1 - depo) * Data@CV_Dep[x])), betaconv(depo, min(depo * 
#       Data@CV_Dep[x], (1 - depo) * Data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
#     Bt_K <- Bt_K[Bt_K > 0.00999 & Bt_K < 0.99001][1]  # interval censor (0.01,0.99)  as in Dick and MacCall 2011
#     Mdb <- trlnorm(100, Data@Mort[x], Data@CV_Mort[x])
#     Mdb <- Mdb[Mdb < 0.9][1]  # !!!! maximum M is 0.9   interval censor
#     if (is.na(Mdb)) 
#       Mdb <- 0.9  # !!!! maximum M is 0.9   absolute limit
#     FMSY_M <- trlnorm(1, Data@FMSY_M[x], Data@CV_FMSY_M[x])
#     BMSY_K <- rbeta(100, alphaconv(Data@BMSY_B0[x], Data@CV_BMSY_B0[x] * 
#       Data@BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@CV_BMSY_B0[x] * 
#       Data@BMSY_B0[x]))  #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
#     tryBMSY_K <- BMSY_K[BMSY_K > 0.05 & BMSY_K < 0.95][1]  # interval censor (0.05,0.95) as in Dick and MacCall, 2011
#     if (is.na(tryBMSY_K)) {
#       Min <- min(BMSY_K, na.rm = TRUE)
#       Max <- max(BMSY_K, na.rm = TRUE)
#       if (Max <= 0.05) 
#         BMSY_K <- 0.05
#       if (Min >= 0.95) 
#         BMSY_K <- 0.95
#     }
#     if (!is.na(tryBMSY_K))  BMSY_K <- tryBMSY_K
#     adelay <- max(floor(iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], 
#       Data@L50[x])), 1)
#     opt <- optimize(DBSRAopt, log(c(0.01 * mean(C_hist), 1000 * mean(C_hist))), 
#       C_hist = C_hist, nys = length(C_hist), Mdb = Mdb, FMSY_M = FMSY_M, 
#       BMSY_K = BMSY_K, Bt_K = Bt_K, adelay = adelay, tol = 0.01)
#     # if(opt$objective<0.1){
#     Kc <- exp(opt$minimum)
#     BMSYc <- Kc * BMSY_K
#     FMSYc <- Mdb * FMSY_M
#     UMSYc <- (FMSYc/(FMSYc + Mdb)) * (1 - exp(-(FMSYc + Mdb)))
#     MSYc <- Kc * BMSY_K * UMSYc
#     TAC[DBSRAcount] <- UMSYc * Kc * Bt_K
#     # 40-10 rule
#     if (Bt_K < 0.4 & Bt_K > 0.1) 
#       TAC[DBSRAcount] <- TAC[DBSRAcount] * (Bt_K - 0.1)/0.3
#     if (Bt_K < 0.1) 
#       TAC[DBSRAcount] <- TAC[DBSRAcount] * tiny  # this has to still be stochastic albeit very small
#     DBSRAcount <- DBSRAcount + 1
#     # }
#   }  # end of reps
#   TACfilter(TAC)
#   # }
# }  # end of DBSRA_apply
# class(DBSRA4010) <- "Output"

DBSRAopt <- function(lnK, C_hist, nys, Mdb, FMSY_M, BMSY_K, Bt_K, adelay) {
  # the optimization for B0 given DBSRA assumptions
  Kc <- exp(lnK)
  n <- getn(BMSY_K)
  g <- gety(n)
  FMSY <- FMSY_M * Mdb
  UMSY <- (FMSY/(FMSY + Mdb)) * (1 - exp(-(FMSY + Mdb)))
  MSY <- Kc * BMSY_K * UMSY
  # Bjoin rules from Dick & MacCall 2011
  Bjoin_K <- 0.5
  if (BMSY_K < 0.3) 
    Bjoin_K <- 0.5 * BMSY_K
  if (BMSY_K > 0.3 & BMSY_K < 0.5) 
    Bjoin_K <- 0.75 * BMSY_K - 0.075
  Bjoin <- Bjoin_K * Kc
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
      Bc[yr] <- Bc[yr - 1] + Bc[yref] * ((PBjoin/Bjoin) + cp * (Bc[yref] - 
        Bjoin)) - C_hist[yr - 1]
    }
    if (Bc[yr] < 0) 
      obj <- obj + log(-Bc[yr])
    Bc[yr] <- max(1e-06, Bc[yr])
  }
  obj + ((Bc[nys]/Kc) - Bt_K)^2
}  # end of DBSRA optimization function

# DCAC

# Variables for parallel processing of DCAC_apply function using
# sfSapply()
C_tot <- nyearsDCAC <- NULL



# #' Depletion Corrected Average Catch
# #' 
# #' A method of calculating an MSY proxy (FMSY * BMSY and therefore the OFL at
# #' most productive stock size) based on average catches accounting for the
# #' windfall catch that got the stock down to BMSY levels.
# #' 
# #' 
# #' @usage DCAC(x, Data, reps = 100)
# #' @param x A position in a data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of stochastic samples of the TAC recommendation
# #' @note It's probably worth noting that DCAC TAC recommendations do not tend
# #' to zero as depletion tends to zero. It adjusts for depletion only in
# #' calculating historical average catch. It follows that at stock levels much
# #' below BMSY, DCAC tends to chronically overfish.
# #' @author T. Carruthers
# #' @references MacCall, A.D., 2009. Depletion-corrected average catch: a simple
# #' formula for estimating sustainable yields in data-poor situations. ICES J.
# #' Mar. Sci. 66, 2267-2271.
# #' @export DCAC
# DCAC <- function(x, Data, reps = 100) {
#   dependencies = "Data@AvC, Data@t, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@Dt, Data@CV_Dt, Data@BMSY_B0, Data@CV_BMSY_B0"
#   C_tot <- Data@AvC[x] * Data@t[x]
#   Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
#   FMSY_M <- trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])  # standard deviation of 0.2 - referred to as 'standard error' in MacCall 2009
#   Bt_K <- trlnorm(reps, Data@Dt[x], Data@CV_Dt[x])
#   if (any(is.na(c(Data@BMSY_B0[x], Data@CV_BMSY_B0[x])))) 
#     return(NA)
#   BMSY_K <- rbeta(reps, alphaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
#     Data@CV_BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
#     Data@CV_BMSY_B0[x]))  #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
#   TACfilter(C_tot/(Data@t[x] + ((1 - Bt_K)/(BMSY_K * FMSY_M * Mdb))))
# }  # end of DCAC
# class(DCAC) <- "Output"



# #' Depletion Corrected Average Catch paired with the 40-10 rule
# #' 
# #' A method of calculating an MSY proxy (FMSY * BMSY and therefore the OFL at
# #' most productive stock size) based on average catches accounting for the
# #' windfall catch that got the stock down to BMSY levels. In this method DCAC
# #' is paired with the 40-10 rule that throttles back the OFL to zero at 10
# #' percent of unfished stock size (the OFL is not subject to downward
# #' adjustment above 40 percent unfished)
# #' 
# #' 
# #' @usage DCAC4010(x, Data, reps = 100)
# #' @param x A position in a data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of stochastic samples of the TAC recommendation
# #' @note DCAC can overfish below BMSY levels. The 40-10 harvest control rule
# #' largely resolves this problem providing an MP with surprisingly good
# #' performance even at low stock levels.
# #' @author T. Carruthers
# #' @references MacCall, A.D., 2009. Depletion-corrected average catch: a simple
# #' formula for estimating sustainable yields in data-poor situations. ICES J.
# #' Mar. Sci. 66, 2267-2271.
# #' @export DCAC4010
# DCAC4010 <- function(x, Data, reps = 100) {
#   dependencies = "Data@AvC, Data@t, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@Dt, Data@CV_Dt, Data@BMSY_B0, Data@CV_BMSY_B0"
#   C_tot <- Data@AvC[x] * Data@t[x]
#   Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
#   FMSY_M <- trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])  # standard deviation of 0.2 - referred to as 'standard error' in MacCall 2009
#   Bt_K <- trlnorm(reps, Data@Dt[x], Data@CV_Dt[x])
#   if (any(is.na(c(Data@BMSY_B0[x], Data@CV_BMSY_B0[x])))) 
#     return(NA)
#   BMSY_K <- rbeta(reps, alphaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
#     Data@CV_BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
#     Data@CV_BMSY_B0[x]))  #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
#   TAC <- C_tot/(Data@t[x] + ((1 - Bt_K)/(BMSY_K * FMSY_M * Mdb)))
#   # 40-10 rule
#   cond1 <- Bt_K < 0.4 & Bt_K > 0.1
#   cond2 <- Bt_K < 0.1
#   if (length(cond1) > 0) 
#     TAC[cond1] <- TAC[cond1] * (Bt_K[cond1] - 0.1)/0.3
#   if (length(cond2) > 0) 
#     TAC[cond2] <- TAC[cond2] * tiny  # this has to still be stochastic albeit very small
#   if (length(cond1) < 1 & length(cond2) < 1) 
#     return(NA)
#   TACfilter(TAC)
#   
# }  # end of DCAC
# class(DCAC4010) <- "Output"
# 
# 
# 
# #' Depletion Corrected Average Catch assuming 40 per cent stock depletion
# #' 
# #' DCAC assuming that current stock biomass is exactly 40 per cent of unfished
# #' levels.
# #' 
# #' 
# #' @usage DCAC_40(x, Data, reps = 100)
# #' @param x A position in a data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of stochastic samples of the TAC recommendation
# #' @note The 40 percent depletion assumption doesn't really affect DCAC that
# #' much as it already makes TAC recommendations that are quite MSY-like.
# #' @author T. Carruthers
# #' @references MacCall, A.D., 2009. Depletion-corrected average catch: a simple
# #' formula for estimating sustainable yields in data-poor situations. ICES J.
# #' Mar. Sci. 66, 2267-2271.
# #' @export DCAC_40
# DCAC_40 <- function(x, Data, reps = 100) {
#   dependencies = "Data@AvC, Data@t, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@BMSY_B0, Data@CV_BMSY_B0"
#   C_tot <- Data@AvC[x] * Data@t[x]
#   Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
#   FMSY_M <- trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])
#   Bt_K <- 0.4
#   if (any(is.na(c(Data@BMSY_B0[x], Data@CV_BMSY_B0[x])))) 
#     return(NA)
#   BMSY_K <- rbeta(reps, alphaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
#     Data@CV_BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
#     Data@CV_BMSY_B0[x]))  #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
#   TACfilter(C_tot/(Data@t[x] + ((1 - Bt_K)/(BMSY_K * FMSY_M * Mdb))))
# }  # end of DCAC40
# class(DCAC_40) <- "Output"
# 
# 
# 
# #' Depletion-Based Stock Reduction Analysis using mean-length estimator of
# #' current depletion
# #' 
# #' DCAC that uses the mean length estimator to calculate current stock
# #' depletion.
# #' 
# #' 
# #' @usage DCAC_ML(x, Data, reps = 100)
# #' @param x A position in a data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of stochastic samples of the TAC recommendation
# #' @note The mean length extension was programmed by Gary Nelson as part of his
# #' excellent R package 'fishmethods'
# #' @author T. Carruthers
# #' @references MacCall, A.D., 2009. Depletion-corrected average catch: a simple
# #' formula for estimating sustainable yields in data-poor situations. ICES J.
# #' Mar. Sci. 66, 2267-2271.
# #' @export DCAC_ML
# DCAC_ML <- function(x, Data, reps = 100) {
#   
#   dependencies = "Data@AvC, Data@t, Data@Mort, Data@CV_Mort, Data@FMSY_M, Data@CV_FMSY_M, Data@BMSY_B0, Data@CV_BMSY_B0, Data@Year, Data@CAL, Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK"
#   if (is.na(Data@BMSY_B0[x]) | is.na(Data@CV_BMSY_B0[x])) 
#     return(NA)
#   if (is.na(Data@FMSY_M[x]) | is.na(Data@CV_FMSY_M[x])) 
#     return(NA)
#   C_tot <- Data@AvC[x] * Data@t[x]
#   Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # default CV of 0.5 as in MacCall 2009
#   FMSY_M <- trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])  # standard deviation of 0.2 - referred to as 'standard error' in MacCall 2009
#   Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
#   Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
#   Z <- MLne(x, Data, Linfc = Linfc, Kc = Kc, ML_reps = reps, MLtype = "dep")
#   if (all(is.na(Z))) return(rep(NA, reps))
#   FM <- Z - Mdb
#   nyears <- length(Data@Year)
#   Ct1 <- mean(Data@Cat[x, 1:3])
#   Ct2 <- mean(Data@Cat[x, (nyears - 2):nyears])
#   dep <- rep(c(Ct1, Ct2), each = reps)/(1 - exp(-FM))
#   if (reps == 1) 
#     Bt_K <- dep[2]/dep[1]
#   if (reps > 1) 
#     Bt_K <- dep[, 2]/dep[, 1]
#   if (any(is.na(c(Data@BMSY_B0[x], Data@CV_BMSY_B0[x])))) 
#     return(NA)
#   BMSY_K <- rbeta(reps, alphaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
#     Data@CV_BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
#     Data@CV_BMSY_B0[x]))  #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
#   TAC <- C_tot/(Data@t[x] + ((1 - Bt_K)/(BMSY_K * FMSY_M * Mdb)))
#   
#   TACfilter(TAC)
# }  # end of DCAC_ML
# class(DCAC_ML) <- "Output"
# 


# #' Beddington and Kirkwood life-history MP (simple version)
# #' 
# #' Sets an OFL according to current abundance and an approximation of Fmax
# #' based on length at first capture.
# #' 
# #' 
# #' @usage BK(x, Data, reps = 100)
# #' @param x A position in a data-limited methods data object.
# #' @param Data A data-limited methods data object.
# #' @param reps The number of stochastic samples of the TAC recommendation
# #' @note This is the simple version of the BK MP. The paper has a more complex
# #' approach that might work better.
# #' @author T. Carruthers.
# #' @references Beddington, J.R., Kirkwood, G.P., 2005. The estimation of
# #' potential yield and stock status using life history parameters. Philos.
# #' Trans. R. Soc. Lond. B Biol. Sci. 360, 163-170.
# #' @export BK
# BK <- function(x, Data, reps = 100) {
#   # Beddington and Kirkwood life-history analysis
#   # ==============================================
#   dependencies = "Data@LFC, Data@vbLinf, Data@CV_vbLinf, Data@Abun, Data@CV_Abun, Data@vbK, Data@CV_vbK"
#   Lc <- trlnorm(reps * 10, Data@LFC[x], 0.2)
#   Linfc <- trlnorm(reps * 10, Data@vbLinf[x], Data@CV_vbLinf[x])
#   Ac <- trlnorm(reps * 10, Data@Abun[x], Data@CV_Abun[x])
#   Kc <- trlnorm(reps * 10, Data@vbK[x], Data@CV_vbK[x])
#   TAC <- Ac * (0.6 * Kc)/(0.67 - (Lc/Linfc))  # robustifying for use in MSE
#   TACfilter(TAC[TAC > 0][1:reps])  # Interval censor only those positive catch recommendations
#   
# }  # end of BK
# class(BK) <- "Output"
# 
# #' Beddington and Kirkwood life-history method combined with catch curve
# #' analysis
# #' 
# #' Calculates an OFL using a catch curve estimate of current F and an
# #' approximation of FMSY based on length at first capture.
# #' 
# #' 
# #' @usage BK_CC(x, Data, reps = 100, Fmin=0.005)
# #' @param x Position in a data-limited methods data object
# #' @param Data A data-limited methods data object (class Data)
# #' @param reps The number of samples of the TAC recommendation
# #' @param Fmin The minimum fishing mortality rate that is derived from the
# #' catch-curve (interval censor)
# #' @author T. Carruthers
# #' @references Beddington, J.R., Kirkwood, G.P., 2005. The estimation of
# #' potential yield and stock status using life history parameters. Philos.
# #' Trans. R. Soc. Lond. B Biol. Sci. 360, 163-170.
# #' @export BK_CC
# BK_CC <- function(x, Data, reps = 100, Fmin = 0.005) {
#   dependencies = "Data@LFC, Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@CAA, Data@Mort"
#   Lc <- trlnorm(reps, Data@LFC[x], 0.2)
#   Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
#   Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
#   Mdb <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])
#   MuC <- Data@Cat[x, length(Data@Cat[x, ])]
#   Cc <- trlnorm(reps, MuC, Data@CV_Cat[x])
#   Zdb <- CC(x, Data, reps = reps * 10)
#   Fdb <- Zdb - Mdb
#   ind <- (1:(reps * 10))[Fdb > Fmin][1:reps]
#   Fdb <- Fdb[ind]
#   Mdb <- Mdb[ind]
#   SM <- sum(is.na(ind))
#   if (SM > 0) {
#     Mdb[is.na(ind)] <- trlnorm(SM, Data@Mort[x], Data@CV_Mort[x])
#     Fdb[is.na(ind)] <- Fmin
#   }
#   
#   Ac <- Cc/(1 - exp(-Fdb))
#   TAC <- Ac * (0.6 * Kc)/(0.67 - (Lc/Linfc))
#   TACfilter(TAC)
#   
# }  # end of BK_CC
# class(BK_CC) <- "Output"
# 


# #' Beddington and Kirkwood life-history analysis with mean-length estimator of
# #' current abundance
# #' 
# #' Uses an approximation to FMSY based on length at first capture and an
# #' estimate of current abundance based on a mean-length estimator.
# #' 
# #' 
# #' @usage BK_ML(x, Data, reps = 100)
# #' @param x Position in a data-limited methods data object
# #' @param Data A data-limited methods data object (class Data)
# #' @param reps The number of samples of the TAC recommendation
# #' @note The mean length extension was programmed by Gary Nelson as part of his
# #' excellent R package 'fishmethods'
# #' @author T. Carruthers
# #' @references Beddington, J.R., Kirkwood, G.P., 2005. The estimation of
# #' potential yield and stock status using life history parameters. Philos.
# #' Trans. R. Soc. Lond. B Biol. Sci. 360, 163-170.
# #' @export BK_ML
# BK_ML <- function(x, Data, reps = 100) {
#   dependencies = "Data@LFC, Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@CAL, Data@Mort"
#   Lc <- trlnorm(reps * 10, Data@LFC[x], 0.2)
#   Linfc <- trlnorm(reps * 10, Data@vbLinf[x], Data@CV_vbLinf[x])
#   Kc <- trlnorm(reps * 10, Data@vbK[x], Data@CV_vbK[x])
#   Mdb <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])
#   Z <- MLne(x, Data, Linfc = Linfc, Kc = Kc, ML_reps = reps * 10, MLtype = "F")
#   if (all(is.na(Z))) return(rep(NA, reps))
#   FM <- Z - Mdb
#   MuC <- Data@Cat[x, length(Data@Cat[x, ])]
#   Cc <- trlnorm(reps * 10, MuC, Data@CV_Cat[x])
#   Ac <- Cc/(1 - exp(-FM))
#   FMSY <- (0.6 * Kc)/(0.67 - (Lc/Linfc))  # robustifying for use in MSETAC<-Ac*FMSY
#   TAC <- Ac * FMSY
#   TAC <- TAC[TAC > 0 & TAC < (mean(TAC, na.rm = T) + 3 * stats::sd(TAC, na.rm = T))][1:reps]
#   TACfilter(TAC)
# }
# class(BK_ML) <- "Output"



# #' An FMSY/M ratio method
# #' 
# #' Calculates the OFL based on a fixed ratio of FMSY to M multiplied by a
# #' current estimate of abundance.
# #' 
# #' A simple method that tends to outperform many other approaches alarmingly
# #' often even when current biomass is relatively poorly known. The low stock
# #' crash potential is largely due to the quite large difference between Fmax
# #' and FMSY for most stocks.
# #' 
# #' @usage Fratio(x, Data, reps = 100)
# #' @param x A position in a data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of samples of the TAC recommendation
# #' @return A numeric vector of TAC recommendations
# #' @author T. Carruthers
# #' @references Gulland, J.A., 1971. The fish resources of the ocean. Fishing
# #' News Books, West Byfleet, UK.
# #' 
# #' Martell, S., Froese, R., 2012. A simple method for estimating MSY from catch
# #' and resilience. Fish Fish. doi: 10.1111/j.1467-2979.2012.00485.x.
# #' @export Fratio
# Fratio <- function(x, Data, reps = 100) {
#   # FMSY / M ratio method e.g. Gulland
#   # ===============================================================================
#   depends = "Data@Abun,Data@CV_Abun,Data@FMSY_M, Data@CV_FMSY_M,Data@Mort,Data@CV_Mort"
#   Ac <- trlnorm(reps, Data@Abun[x], Data@CV_Abun[x])
#   TACfilter(Ac * trlnorm(reps, Data@Mort[x], Data@CV_Mort[x]) * 
#     trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x]))
# }  # end of Fratio
# class(Fratio) <- "Output"
# 


# #' An FMSY/M ratio method paired with the 40-10 rule
# #' 
# #' Calculates the OFL based on a fixed ratio of FMSY to M multiplied by a
# #' current estimate of abundance. In this method DBSRA is paired with the 40-10
# #' rule that throttles back the OFL to zero at 10 percent of unfished biomass.
# #' 
# #' 
# #' @usage Fratio4010(x, Data, reps = 100)
# #' @param x A position in data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of TAC samples
# #' @author T. Carruthers
# #' @references Gulland, J.A., 1971. The fish resources of the ocean. Fishing
# #' News Books, West Byfleet, UK.
# #' 
# #' Martell, S., Froese, R., 2012. A simple method for estimating MSY from catch
# #' and resilience. Fish Fish. doi: 10.1111/j.1467-2979.2012.00485.x.
# #' @export Fratio4010
# Fratio4010 <- function(x, Data, reps = 100) {
#   # FMSY / M ratio method e.g. Gulland
#   # ===============================================================================
#   dependencies = "Data@Abun, Data@CV_Abun, Data@FMSY_M, Data@CV_FMSY_M, Data@Mort, Data@CV_Mort, Data@Dep"
#   Ac <- trlnorm(reps, Data@Abun[x], Data@CV_Abun[x])
#   TAC <- Ac * trlnorm(reps, Data@Mort[x], Data@CV_Mort[x]) * 
#     trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])
#   Bt_K <- trlnorm(reps, Data@Dt[x], Data@CV_Dt[x])
#   # 40-10 rule
#   cond1 <- Bt_K < 0.4 & Bt_K > 0.1
#   cond2 <- Bt_K < 0.1
#   TAC[cond1] <- TAC[cond1] * (Bt_K[cond1] - 0.1)/0.3
#   TAC[cond2] <- TAC[cond2] * tiny  # this has to still be stochastic albeit very small
#   TACfilter(TAC)
# }  # end of Fratio
# class(Fratio4010) <- "Output"
# 


# #' A data-limited method that uses FMSY/M ratio and a naive catch-curve
# #' estimate of recent Z
# #' 
# #' Calculates the OFL based on a fixed ratio of FMSY to M and a catch curve
# #' estimate of current stock size.
# #' 
# #' 
# #' @usage Fratio_CC(x, Data, reps = 100, Fmin = 0.005)
# #' @param x A position in data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of TAC samples
# #' @param Fmin Minimum current fishing mortality rate for the catch-curve
# #' analysis
# #' @author T. Carruthers
# #' @references Gulland, J.A., 1971. The fish resources of the ocean. Fishing
# #' News Books, West Byfleet, UK.
# #' 
# #' Martell, S., Froese, R., 2012. A simple method for estimating MSY from catch
# #' and resilience. Fish Fish. doi: 10.1111/j.1467-2979.2012.00485.x.
# #' @export Fratio_CC
# Fratio_CC <- function(x, Data, reps = 100, Fmin = 0.005) {
#   # FMSY / M ratio method using catch curve analysis to determine current
#   # abundance ================================== for (x in 1:nsim) {
#   dependencies = " Data@FMSY_M, Data@CV_FMSY_M, Data@Mort, Data@CV_Mort, Data@Cat, Data@CV_Cat, Data@CAA"
#   MuC <- Data@Cat[x, length(Data@Cat[x, ])]
#   Cc <- trlnorm(reps, MuC, Data@CV_Cat[x])
#   Mdb <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
#   Zdb <- CC(x, Data, reps = reps * 10)
#   Fdb <- Zdb - Mdb
#   ind <- (1:(reps * 10))[Fdb > 0.005][1:reps]
#   
#   Fdb <- Fdb[ind]
#   Mdb <- Mdb[ind]
#   SM <- sum(is.na(ind))
#   if (SM > 0) {
#     Mdb[is.na(ind)] <- trlnorm(SM, Data@Mort[x], Data@CV_Mort[x])
#     Fdb[is.na(ind)] <- Fmin
#   }
#   
#   Ac <- Cc/(1 - exp(-Fdb))
#   TAC <- Ac * Mdb * trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])
#   
#   TACfilter(TAC)
#   # }
# }  # end of Fratio_CC
# class(Fratio_CC) <- "Output"
# 



# #' An FMSY/M ratio MP that uses a mean length estimator of recent Z
# #' 
# #' Calculates the OFL based on a fixed ratio of FMSY/M and an estimate of
# #' current stock size from a mean-length estimator.
# #' 
# #' 
# #' @usage Fratio_ML(x, Data, reps = 100)
# #' @param x A position in data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of TAC samples
# #' @note The mean length extension was programmed by Gary Nelson as part of his
# #' excellent R package 'fishmethods'
# #' @author T. Carruthers
# #' @references Gulland, J.A., 1971. The fish resources of the ocean. Fishing
# #' News Books, West Byfleet, UK.
# #' 
# #' Martell, S., Froese, R., 2012. A simple method for estimating MSY from catch
# #' and resilience. Fish Fish. doi: 10.1111/j.1467-2979.2012.00485.x.
# #' @export Fratio_ML
# Fratio_ML <- function(x, Data, reps = 100) {
#   dependencies = " Data@FMSY_M, Data@CV_FMSY_M, Data@Mort, Data@CV_Mort, Data@Cat, Data@CV_Cat, Data@CAL"
#   MuC <- Data@Cat[x, length(Data@Cat[x, ])]
#   Cc <- trlnorm(reps * 10, MuC, Data@CV_Cat[x])
#   Mdb <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
#   Linfc <- trlnorm(reps * 10, Data@vbLinf[x], Data@CV_vbLinf[x])
#   Kc <- trlnorm(reps * 10, Data@vbK[x], Data@CV_vbK[x])
#   Z <- MLne(x, Data, Linfc = Linfc, Kc = Kc, ML_reps = reps * 10, MLtype = "F")
#   if (all(is.na(Z))) return(rep(NA, reps))
#   FM <- Z - Mdb
#   Ac <- Cc/(1 - exp(-FM))
#   TAC <- Ac * trlnorm(reps * 10, Data@FMSY_M[x], Data@CV_FMSY_M[x]) * 
#     trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
#   TAC <- TAC[TAC > 0][1:reps]
#   TACfilter(TAC)
# }
# class(Fratio_ML) <- "Output"



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
  
  TACfilter(TAC)
}  # end of SPMSY
class(SPMSY) <- "Output"



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
  TACfilter(TAC)
}
class(MCD) <- "Output"



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
  
  TACfilter(TAC)
}
class(MCD4010) <- "Output"



#' Surplus Production Stock Reduction Analysis
#' 
#' A surplus production equivalent of DB-SRA that uses a demographically
#' derived prior for intrinsic rate of increase (McAllister method, below)
#' 
#' 
#' @usage SPSRA(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object (class DLM)
#' @param reps The number of samples of the TAC taken for the calculation of
#' the quota
#' @author T. Carruthers
#' @references McAllister, M.K., Pikitch, E.K., and Babcock, E.A. 2001. Using
#' demographic methods to construct Bayesian priors for the intrinsic rate of
#' increase in the Schaefer model and implications for stock rebuilding. Can.
#' J. Fish. Aquat. Sci. 58: 1871-1890.
#' @export SPSRA
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
  TACfilter(TAC)
}
class(SPSRA) <- "Output"



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
  if (all(is.na(Z))) return(rep(NA, reps))
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
  TACfilter(TAC)
}
class(SPSRA_ML) <- "Output"

SPSRAopt <- function(lnK, dep, r, Ct, PE) {
  nyears <- length(Ct)
  B <- rep(NA, nyears)
  B[1] <- exp(lnK)
  OBJ <- 0
  for (y in 2:nyears) {
    if ((B[y - 1] - Ct[y - 1]) < 0) 
      OBJ <- OBJ + (B[y - 1] - Ct[y - 1])^2
    B[y] <- max(0.01, B[y - 1] - Ct[y - 1])
    B[y] <- B[y] + r * B[y] * (1 - B[y]/B[1]) * PE[y]
  }
  return(OBJ + ((B[nyears]/B[1]) - dep)^2)
}




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
  TACfilter(TAC)
  # }
}  # end of YPR
class(YPR) <- "Output"



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
  TACfilter(TAC)
}
class(YPR_CC) <- "Output"



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
  if (all(is.na(Z))) return(rep(NA, reps))
  FM <- Z - Mdb
  Ac <- Cc/(1 - exp(-FM))
  FMSY <- YPRopt(Linfc, Kc, t0c, Mdb, a, b, LFS, Data@MaxAge, reps * 
    10)
  TAC <- Ac * FMSY
  TAC <- TAC[TAC > 0][1:reps]
  TACfilter(TAC)
  
}
class(YPR_ML) <- "Output"



# #' Demographic FMSY method
# #' 
# #' FMSY is calculated as r/2 where r is calculated from a demographic approach
# #' (inc steepness). Coupled with an estimate of current abundance that gives
# #' you the OFL.
# #' 
# #' Made up for this package. This uses Murdoch McAllister's demographic r
# #' method to derive FMSY (r/2) and then makes the quota r*current biomass / 2.
# #' Easy.
# #' 
# #' @usage Fdem(x, Data, reps = 100)
# #' @param x A position in data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of TAC samples
# #' @author T. Carruthers
# #' @references McAllister, M.K., Pikitch, E.K., and Babcock, E.A. 2001. Using
# #' demographic methods to construct Bayesian priors for the intrinsic rate of
# #' increase in the Schaefer model and implications for stock rebuilding. Can.
# #' J. Fish. Aquat. Sci. 58: 1871-1890.
# #' @export Fdem
# Fdem <- function(x, Data, reps = 100) {
#   # Demographic FMSY estimate (FMSY=r/2)
#   dependencies = "Data@Mort, Data@CV_Mort, Data@vbK, Data@CV_vbK, Data@vbLinf, Data@CV_vbLinf, Data@vbt0, Data@CV_vbt0, Data@MaxAge, Data@wla, Data@wlb, Data@Abun, Data@CV_Abun, Data@steep, Data@CV_steep"
#   Mvec <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
#   Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
#   Linfc = trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
#   if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
#     t0c <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
#   } else {
#     t0c <- rep(Data@vbt0[x], reps)
#   }
#   t0c[!is.finite(t0c)] <- 0
#   # hvec <- trlnorm(reps, Data@steep[x], Data@CV_steep[x])
#   hvec <- sample_steepness2(reps, Data@steep[x], Data@CV_steep[x])
#   
#   Ac <- trlnorm(reps, Data@Abun[x], Data@CV_Abun[x])
#   FMSY <- getr(x, Data, Mvec, Kc, Linfc, t0c, hvec, maxage = Data@MaxAge, 
#     r_reps = reps)/2
#   TAC <- FMSY * Ac
#   TACfilter(TAC)
# }
# class(Fdem) <- "Output"
# 


# #' Demographic FMSY method using catch-curve analysis to estimate recent Z
# #' 
# #' FMSY is calculated as r/2 from a demographic r prior method, current
# #' abudnance is estimated from naive catch curve analysis.
# #' 
# #' 
# #' @usage Fdem_CC(x, Data, reps = 100, Fmin=0.005)
# #' @param x A position in data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of TAC samples
# #' @param Fmin The minimum fishing mortality rate derived from the catch-curve
# #' analysis
# #' @author T. Carruthers
# #' @references McAllister, M.K., Pikitch, E.K., and Babcock, E.A. 2001. Using
# #' demographic methods to construct Bayesian priors for the intrinsic rate of
# #' increase in the Schaefer model and implications for stock rebuilding. Can.
# #' J. Fish. Aquat. Sci. 58: 1871-1890.
# #' @export Fdem_CC
# Fdem_CC <- function(x, Data, reps = 100, Fmin = 0.005) {
#   dependencies = "Data@Mort, Data@CV_Mort, Data@vbK, Data@CV_vbK, Data@vbLinf, Data@CV_vbLinf, Data@vbt0, Data@CV_vbt0, Data@MaxAge, Data@wla, Data@wlb, Data@CAA, Data@steep, Data@CV_steep"
#   Mvec <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])
#   Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
#   Linfc = trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
#   if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
#     t0c <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
#   } else {
#     t0c <- rep(Data@vbt0[x], reps)
#   }
#   t0c[!is.finite(t0c)] <- 0
#   # hvec <- trlnorm(reps, Data@steep[x], Data@CV_steep[x])
#   hvec <- sample_steepness2(reps, Data@steep[x], Data@CV_steep[x])
#   MuC <- Data@Cat[x, length(Data@Cat[x, ])]
#   Cc <- trlnorm(reps, MuC, Data@CV_Cat[x])
#   Zdb <- CC(x, Data, reps = reps * 10)
#   Fdb <- Zdb - Mvec
#   ind <- (1:(reps * 10))[Fdb > Fmin][1:reps]
#   
#   Fdb <- Fdb[ind]
#   SM <- sum(is.na(ind))
#   if (SM > 0) {
#     Fdb[is.na(ind)] <- Fmin
#   }
#   
#   Ac <- Cc/(1 - exp(-Fdb))
#   FMSY <- getr(x, Data, Mvec, Kc, Linfc, t0c, hvec, maxage = Data@MaxAge, r_reps = reps)/2
#   TAC <- FMSY * Ac
#   
#   TACfilter(TAC)
# }
# class(Fdem_CC) <- "Output"
# 
# Catch curve estimate of recent F (naive)
# ========================================================================================
CC <- function(x, Data, reps = 100) {
  ny <- dim(Data@CAA)[2]
  CAA <- apply(Data@CAA[x, max(ny - 2, 1):ny, ], 2, sum)  # takes last two years as the sample (or last year if there is only one)
  maxageobs <- length(CAA)
  AFS <- which.max(CAA)
  AFS[AFS > (maxageobs - 3)] <- maxageobs - 3  # provides at least three datapoints
  
  nS <- ceiling(sum(CAA)/2)
  y <- log(CAA[AFS:maxageobs]/sum(CAA[AFS:maxageobs], na.rm = T))
  xc <- 1:length(y)
  y[y == "-Inf"] <- NA
  mod <- lm(y ~ xc)
  chk <- sum(is.na(coef(mod)))  # check if model failed
  if (chk) {
    return(NA)
  } else {
    coefs <- summary(mod, weights = CAA[AFS:maxageobs])$coefficients[2, 1:2]
    coefs[is.nan(coefs)] <- tiny
    return(-rnorm(reps, coefs[1], coefs[2]))
  }
}
# class(CC)<-'Output'


# 
# #' Demographic FMSY method that uses mean length data to estimate recent Z
# #' 
# #' Demographic F (r/2) method using the mean length estimator to calculate
# #' current abundance.
# #' 
# #' 
# #' @usage Fdem_ML(x, Data, reps = 100)
# #' @param x A position in data-limited methods data object
# #' @param Data A data-limited methods data object
# #' @param reps The number of TAC samples
# #' @note The mean length extension was programmed by Gary Nelson as part of his
# #' excellent R package 'fishmethods'
# #' @author T. Carruthers
# #' @references McAllister, M.K., Pikitch, E.K., and Babcock, E.A. 2001. Using
# #' demographic methods to construct Bayesian priors for the intrinsic rate of
# #' increase in the Schaefer model and implications for stock rebuilding. Can.
# #' J. Fish. Aquat. Sci. 58: 1871-1890.
# #' @export Fdem_ML
# Fdem_ML <- function(x, Data, reps = 100) {
#   dependencies = "Data@Mort, Data@CV_Mort, Data@vbK, Data@CV_vbK, Data@vbLinf, Data@CV_vbLinf, Data@vbt0, Data@CV_vbt0, Data@MaxAge, Data@wla, Data@wlb, Data@CAL, Data@steep, Data@CV_steep"
#   Mvec <- trlnorm(reps * 10, Data@Mort[x], Data@CV_Mort[x])
#   Kc <- trlnorm(reps * 10, Data@vbK[x], Data@CV_vbK[x])
#   Linfc = trlnorm(reps * 10, Data@vbLinf[x], Data@CV_vbLinf[x])
#   t0c <- -trlnorm(reps * 10, -Data@vbt0[x], Data@CV_vbt0[x])
#   t0c[!is.finite(t0c)] <- 0
#   # hvec <- trlnorm(reps * 10, Data@steep[x], Data@CV_steep[x])
#   hvec <- sample_steepness2(reps*10, Data@steep[x], Data@CV_steep[x])
#   MuC <- Data@Cat[x, length(Data@Cat[x, ])]
#   Cc <- trlnorm(reps * 10, MuC, Data@CV_Cat[x])
#   Z <- MLne(x, Data, Linfc = Linfc, Kc = Kc, ML_reps = reps * 10, MLtype = "F")
# 
#   if (all(is.na(Z))) return(rep(NA, reps))
#   ind <- !is.na(Z)
#   FM <- Z[ind] - Mvec[ind]
#   Ac <- Cc[ind]/(1 - exp(-FM))
#   FMSY <- getr(x, Data, Mvec[ind], Kc[ind], Linfc[ind], t0c[ind], hvec[ind], 
#                maxage = Data@MaxAge, r_reps = sum(ind))/2
#   TAC <- FMSY * Ac
#   TAC <- TAC[TAC > 0][1:reps]
#   TACfilter(TAC)
# }
# class(Fdem_ML) <- "Output"



#' Age-composition-based estimate of current stock depletion given constant Z
#' linked to an FMSY estimate to provide OFL
#' 
#' Estimates an OFL based on a Stock Reduction analysis fitted to current
#' age-composition data. Knife-edge vulnerability at age at maturity allows for
#' an FMSY estimate. OFL=FMSY*F/C
#' 
#' 
#' @usage CompSRA(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC.
#' @note Given a fixed historical F, What level of depletion gives you this
#' length composition?
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
  TACfilter(TAC)
}
class(CompSRA) <- "Output"




#' Age-composition-based estimate of current stock depletion given constant Z
#' linked to an FMSY estimate to provide OFL (with a 40-10 rule)
#' 
#' Estimates an OFL based on a Stock Reduction analysis fitted to current
#' age-composition data. Knife-edge vulnerability at age at maturity allows for
#' an FMSY estimate. OFL=FMSY*F/C
#' 
#' 
#' @usage CompSRA4010(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC.
#' @note Given a fixed historical F, What level of depletion gives you this
#' length composition?
#' @author T. Carruthers
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
  
  TACfilter(TAC)
  
  # message(x, ' of ', nsim) flush.console() }
}
class(CompSRA4010) <- "Output"

# options(warn=2) options(warn=1) options(warn=1)

SRAfunc <- function(lnR0c, Mc, hc, maxage, LFSc, LFCc, Linfc, Kc, t0c, 
  AMc, ac, bc, Catch, CAA, opt = 1) {
  
  ny <- length(Catch)
  AFC <- log(1 - min(0.99, LFCc/Linfc))/-Kc + t0c
  AFS <- log(1 - min(0.99, LFSc/Linfc))/-Kc + t0c
  if (AFC >= 0.7 * maxage) 
    AFC <- 0.7 * maxage
  if (AFS >= 0.9 * maxage) 
    AFS <- 0.9 * maxage
  KES <- max(2, ceiling(mean(c(AFC, AFS))))
  vul <- rep(1, maxage)
  vul[1:(KES - 1)] <- 0
  Mac <- rep(1, maxage)
  Mac[1:max(1, floor(AMc))] <- 0
  Lac <- Linfc * (1 - exp(-Kc * ((1:maxage) - t0c)))
  Wac <- ac * Lac^bc
  R0c <- exp(lnR0c)
  N <- exp(-Mc * ((1:maxage) - 1)) * R0c
  SSN <- Mac * N  # Calculate initial spawning stock numbers
  Biomass <- N * Wac
  SSB <- SSN * Wac  # Calculate spawning stock biomass
  
  B0 <- sum(Biomass)
  SSB0 <- sum(SSB)
  SSN0 <- SSN
  SSBpR <- sum(SSB)/R0c  # Calculate spawning stock biomass per recruit
  SSNpR <- SSN/R0c
  
  CN <- array(NA, dim = c(ny, maxage))
  HR <- rep(0, maxage)
  pen <- 0
  for (y in 1:ny) {
    # set up some indices for indexed calculation
    VB <- Biomass[KES:maxage] * exp(-Mc)
    CB <- Catch[y] * VB/sum(VB)
    testHR <- CB[1]/VB[1]
    if (testHR > 0.8) 
      pen <- pen + (testHR - 0.8)^2
    HR[KES:maxage] <- min(testHR, 0.8)
    FMc <- -log(1 - HR)  # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
    Zc <- FMc + Mc
    
    CN[y, ] <- N * (1 - exp(-Zc)) * (FMc/Zc)
    N[2:maxage] <- N[1:(maxage - 1)] * exp(-Zc[1:(maxage - 1)])  # Total mortality
    N[1] <- (0.8 * R0c * hc * sum(SSB))/(0.2 * SSBpR * R0c * (1 - hc) + 
      (hc - 0.2) * sum(SSB))  # Recruitment assuming regional R0 and stock wide steepness
    # print(N[1])
    Biomass <- N * Wac
    SSN <- N * Mac
    SSB <- SSN * Wac
    
  }  # end of year
  
  CN[CN < 0] <- 0  # stop any negative catches
  syear <- ny - dim(CAA)[1] + 1
  pred <- CN[syear:ny, ]
  pred <- pred/array(apply(pred, 1, sum), dim = c(dim(CAA)[1], maxage))
  
  fobj <- pen - sum(log(pred + tiny) * CAA, na.rm = T)
  if (opt == 1) {
    return(fobj)
  } else if (opt == 2) {
    return(sum(Biomass))
  } else if (opt == 3) {
    sum(SSB)/sum(SSB0)
  }
  # CBc<-sum(CB)
}

SRAFMSY <- function(lnFMc, Mc, hc, maxage, LFSc, LFCc, Linfc, Kc, t0c, 
  AMc, ac, bc, opt = T) {
  
  FMc <- exp(lnFMc)
  ny <- 100
  AFC <- log(1 - min(0.99, LFCc/Linfc))/-Kc + t0c
  AFS <- log(1 - min(0.99, LFSc/Linfc))/-Kc + t0c
  if (AFC >= 0.7 * maxage) 
    AFC <- 0.7 * maxage
  if (AFS >= 0.9 * maxage) 
    AFS <- 0.9 * maxage
  
  KES <- max(2, ceiling(mean(c(AFC, AFS))))
  vul <- rep(1, maxage)
  vul[1:(KES - 1)] <- 0
  Mac <- rep(1, maxage)
  Mac[1:max(1, floor(AMc))] <- 0
  Lac <- Linfc * (1 - exp(-Kc * ((1:maxage) - t0c)))
  Wac <- ac * Lac^bc
  R0c <- 1
  N <- exp(-Mc * ((1:maxage) - 1)) * R0c
  SSN <- Mac * N  # Calculate initial spawning stock numbers
  Biomass <- N * Wac
  SSB <- SSN * Wac  # Calculate spawning stock biomass
  
  B0 <- sum(Biomass)
  SSB0 <- sum(SSB)
  SSN0 <- SSN
  SSBpR <- sum(SSB)/R0c  # Calculate spawning stock biomass per recruit
  SSNpR <- SSN/R0c
  
  N <- N/2
  SSN <- Mac * N  # Calculate initial spawning stock numbers
  Biomass <- N * Wac
  SSB <- SSN * Wac
  
  for (y in 1:ny) {
    # set up some indices for indexed calculation Fishing mortality rate
    # determined by effort, catchability, vulnerability and spatial
    # preference according to biomass
    Zc <- FMc * vul + Mc
    CN <- N * (1 - exp(-Zc)) * (FMc/Zc)
    CB <- CN * Wac
    Biomass <- N * Wac
    N[2:maxage] <- N[1:(maxage - 1)] * exp(-Zc[1:(maxage - 1)])  # Total mortality
    N[1] <- (0.8 * R0c * hc * sum(SSB))/(0.2 * SSBpR * R0c * (1 - hc) + 
      (hc - 0.2) * sum(SSB))  # Recruitment assuming regional R0 and stock wide steepness
    # print(N[1])
    SSN <- N * Mac
    SSB <- SSN * Wac
    
  }  # end of year
  
  if (opt) {
    return(-sum(CB))
  } else {
    return(FMc)
  }
}


# Yield per recruit estimate of FMSY Meaghan Bryan 2013
# ==================================
YPRopt = function(Linfc, Kc, t0c, Mdb, a, b, LFS, maxage, reps = 100) {
  
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
    if (tc[i] > 0) 
      vul[i, tc[i]:maxage] <- 1
    if (tc[i] > 1) 
      mat[i, max(1, tc[i] - 1):maxage] <- 1
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
    slope[, i - 1] = round((ypr[, i] - ypr[, i - 1])/(frates[i] - frates[i - 
      1]), 2)
  }
  dif = abs(slope - slope.10)
  dif[is.na(dif)] <- 1e+11
  frates[apply(dif, 1, which.min)]  #frates[which.min(dif)]
}


MLne <- function(x, Data, Linfc, Kc, ML_reps = 100, MLtype = "dep") {
  year <- 1:dim(Data@CAL)[2]
  nlbin <- ncol(Data@CAL[x, , ])
  nlyr <- nrow(Data@CAL[x, , ])
  mlbin <- (Data@CAL_bins[1:nlbin] + Data@CAL_bins[2:(nlbin + 1)])/2
  nbreaks <- 1
  Z <- matrix(NA, nrow = ML_reps, ncol = nbreaks + 1)
  Z2 <- rep(NA, ML_reps)
  # temp<-apply(Data@CAL[x,,],2,sum) Lc<-mlbin[which.max(temp)] #
  # modal length
  
  # dd <- dim(Data@CAL[x,,]) curLen <- Data@CAL[x,dd[1],] Lc <-
  # mlbin[which.max(curLen)] Lc <- Data@LFS[,x] Lc <- Lc[length(Lc)]
  Lc <- Data@LFS[x]
  
  for (i in 1:ML_reps) {
    mlen <- rep(NA, length(year))
    ss <- ceiling(apply(Data@CAL[x, , ], 1, sum)/2)
    if (MLtype == "dep") {
      for (y in 1:length(year)) {
        if (sum(Data@CAL[x, y, ] > 0) > 0.25 * length(Data@CAL[x, y, ])) {
          temp2 <- sample(mlbin, ceiling(sum(Data@CAL[x, y, ])/2), replace = T, prob = Data@CAL[x, y, ])
          mlen[y] <- mean(temp2[temp2 >= Lc], na.rm = TRUE)
        }
      }
      
      fitmod <- bhnoneq(year = year, mlen = mlen, ss = ss, K = Kc[i], Linf = Linfc[i], 
                        Lc = Lc, nbreaks = nbreaks, styrs = ceiling(length(year) * ((1:nbreaks)/(nbreaks + 1))), 
                        stZ = rep(Data@Mort[x], nbreaks + 1))
      if (all(fitmod == FALSE)) {
        Z[i, ] <- NA
      } else Z[i, ] <- fitmod
    } else {
      
      # ind<-(which.min(((Data@CAL_bins-Data@LFS[x])^2)^0.5)-1):(length(Data@CAL_bins)-1)
      for (y in 1:length(year)) {
        if (sum(Data@CAL[x, y, ] > 0) > 0.25 * length(Data@CAL[x, y, ])) {
          temp2 <- sample(mlbin, ceiling(sum(Data@CAL[x, y, ])/2), replace = T, prob = Data@CAL[x, y, ])
          mlen[y] <- mean(temp2[temp2 >= Lc], na.rm = TRUE)
        }
      }
      mlen <- mean(mlen[(length(mlen) - 2):length(mlen)], na.rm = TRUE)
      Z2[i] <- bheq(K = Kc[i], Linf = Linfc[i], Lc = Lc, Lbar = mlen)
      
    }
  }
  # Z <- Z[,ncol(Z)] # last estimate of Z? Z needs to be vector reps long
  if (MLtype == "F") {
    Z2[Z2<0] <- NA
    return(Z2)
  }
  if (MLtype == "dep") return(Z)
}

bheq <- function(K, Linf, Lc, Lbar) {
  K * (Linf - Lbar)/(Lbar - Lc)
}

bhnoneq <- function(year, mlen, ss, K, Linf, Lc, nbreaks, styrs, stZ) {
 
  mlen[mlen <= 0 | is.na(mlen)] <- -99
  ss[ss <= 0 | is.na(ss) | mlen == -99] <- 0
  stpar <- c(stZ, styrs)
  
  # results <-
  # optim(stpar,bhnoneq_LL,method='BFGS',year=year,Lbar=mlen,ss=ss,
  # nbreaks=nbreaks,K=K,Linf=Linf,Lc=Lc,control=list(maxit=1e6))
  results <- try(optim(stpar, bhnoneq_LL, method = "Nelder-Mead", year = year, 
    Lbar = mlen, ss = ss, nbreaks = nbreaks, K = K, Linf = Linf, Lc = Lc, 
    control = list(maxit = 1e+06), hessian = FALSE), silent=TRUE)
  if (class(results) == "try-error") {
    return(FALSE)
  } else return(results$par[1:(nbreaks + 1)])
}

getdep <- function(lnFF, targ, Md, Linfd, Kd, t0d, AFSd, ad, bd, maxage, 
  opt) {
  
  FF <- exp(lnFF)
  Z <- rep(Md, maxage)
  Z[1] <- 0
  Z[AFSd:maxage] <- Z[AFSd:maxage] + FF
  for (a in 2:maxage) Z[a] <- Z[a - 1] + Z[a]
  Nd <- exp(-Z)
  Nobs <- Nd
  Nobs[1:max(1, (AFSd - 1))] <- 0
  Ld <- Linfd * (1 - exp(-Kd * (((1:maxage) - 0.5) - t0d)))
  
  if (opt) 
    return(((sum(Nobs * Ld)/sum(Nobs)) - targ)^2)
  if (!opt) {
    Nd0 <- exp(-Md * ((1:maxage) - 0.5))
    Wd <- ad * Ld^bd
    return((sum(Nd * Wd)/sum(Nd))/(sum(Nd0 * Wd)/sum(Nd0)))
  }
}

getr <- function(x, Data, Mvec, Kvec, Linfvec, t0vec, hvec, maxage, 
  r_reps = 100) {
  r <- rep(NA, r_reps)
  for (i in 1:r_reps) {
    log.r = log(0.3)
    
    opt = optimize(demofn, lower = log(1e-04), upper = log(1.4), M = Mvec[i], 
      amat = iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], 
        Data@L50[x]), sigma = 0.2, K = Kvec[i], Linf = Linfvec[i], 
      to = t0vec[i], hR = hvec[i], maxage = maxage, a = Data@wla[x], 
      b = Data@wlb[x])
    # demographic2(opt$minimum,M[x],ageM[x],0.2,K[x],Linf,t0,steepness[x],maxage,a,b)$r
    r[i] <- exp(opt$minimum)
  }
  r
}

iVB <- function(t0, K, Linf, L) max(1, ((-log(1 - L/Linf))/K + t0))  # Inverse Von-B



#' Depletion Adjusted Average Catch
#' 
#' Essentially DCAC multiplied by 2*depletion and divided by BMSY/B0 (Bpeak)
#' 
#' 
#' @usage DAAC(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @author W. Harford and T. Carruthers
#' @references MacCall, A.D., 2009. Depletion-corrected average catch: a simple
#' formula for estimating sustainable yields in data-poor situations. ICES J.
#' Mar. Sci. 66, 2267-2271. Harford W. and Carruthers, T. 2016. Simulation
#' testing novel catch-based fisheries management. In draft, intended for Fish.
#' Bull.
#' @export DAAC
DAAC <- function(x, Data, reps = 100) {
  # extended depletion-corrected average catch (Harford and Carruthers
  # 2015)
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
  TACfilter(TAC)
}
class(DAAC) <- "Output"



#' Hybrid Depletion Adjusted Average Catch
#' 
#' Essentially DCAC multiplied by 2*depletion and divided by BMSY/B0 (Bpeak)
#' when below BMSY, and DCAC above BMSY
#' 
#' 
#' @usage HDAAC(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @author W. Harford and T. Carruthers
#' @references MacCall, A.D., 2009. Depletion-corrected average catch: a simple
#' formula for estimating sustainable yields in data-poor situations. ICES J.
#' Mar. Sci. 66, 2267-2271. Harford W. and Carruthers, T. 2016. Testing novel
#' catch-based fisheries management procedures.
#' @export HDAAC
HDAAC <- function(x, Data, reps = 100) {
  dependencies = "Data@AvC, Data@t, Data@Mort, Data@CV_Mort, Data@Dt, Data@CV_Dt, Data@BMSY_B0, Data@CV_BMSY_B0"
  C_tot <- Data@AvC[x] * Data@t[x]
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  FMSY_M <- trlnorm(reps, Data@FMSY_M[x], Data@CV_FMSY_M[x])
  Bt_K <- trlnorm(reps, Data@Dt[x], Data@CV_Dt[x])
  if (any(is.na(c(Data@BMSY_B0[x], Data@CV_BMSY_B0[x])))) 
    return(NA)
  BMSY_K <- rbeta(reps, alphaconv(Data@BMSY_B0[x], Data@BMSY_B0[x] * 
    Data@CV_BMSY_B0[x]), betaconv(Data@BMSY_B0[x], Data@CV_BMSY_B0[x]))
  dcac <- C_tot/(Data@t[x] + ((1 - Bt_K)/(BMSY_K * FMSY_M * Mdb)))
  ddcac <- dcac * Bt_K/BMSY_K
  TAC <- dcac
  TAC[Bt_K < BMSY_K] <- ddcac[Bt_K < BMSY_K]
  TACfilter(TAC)
}
class(HDAAC) <- "Output"



# A generic VPA (Walters and Licandeo UBC)
VPA <- function(x, Data, reps = reps) {
  
  # now do optimization for FMSY
  dependencies = "Data@Mort, Data@CV_Mort, Data@vbK, Data@CV_vbK, Data@vbLinf, Data@CV_vbLinf, Data@vbt0, Data@CV_vbt0, Data@MaxAge, Data@wla, Data@CV_wla, Data@wlb, Data@CV_wlb, Data@L50, Data@CV_L50, Data@CAA, Data@steep, Data@CV_steep, Data@LFS, Data@CV_LFS, Data@LFC, Data@CV_LFC, Data@Cat"
  CAAind <- (Data@CAA[x, , ] == 0) * array(rep(1:Data@MaxAge, 
    each = length(Data@CAA[x, , 1])), dim(Data@CAA[x, , ]))
  maxage <- min(CAAind[CAAind != 0])
  maxage <- which.min(abs(cumsum(apply(Data@CAA[x, , ], 2, sum))/sum(Data@CAA[x, 
    , ]) - 0.75))
  CAAv <- Data@CAA[x, , 1:maxage]
  CAAv[, maxage] <- CAAv[, maxage] + apply(Data@CAA[x, , (maxage + 
    1):length(Data@CAA[x, 1, ])], 1, sum)
  
  TAC <- Bt_K <- rep(NA, reps)
  
  for (i in 1:reps) {
    
    Mc <- trlnorm(1, Data@Mort[x], Data@CV_Mort[x])
    # hc <- trlnorm(1, Data@steep[x], Data@CV_steep[x])
    hc <- sample_steepness2(1, Data@steep[x], Data@CV_steep[x])
    Linfc <- trlnorm(1, Data@vbLinf[x], Data@CV_vbLinf[x])
    Kc <- trlnorm(1, Data@vbK[x], Data@CV_vbK[x])
    t0c <- -trlnorm(1, -Data@vbt0[x], Data@CV_vbt0[x])
    LFSc <- trlnorm(1, Data@LFS[x], Data@CV_LFS[x])
    LFCc <- trlnorm(1, Data@LFC[x], Data@CV_LFC[x])
    AMc <- trlnorm(1, iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], 
      Data@L50[x]), Data@CV_L50[x])
    ac <- trlnorm(1, Data@wla[x], Data@CV_wla[x])
    bc <- trlnorm(1, Data@wlb[x], Data@CV_wlb[x])
    
    pmat <- rep(1, maxage)
    pmat[1:ceiling(AMc)] <- 0
    age <- 1:maxage
    la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - 
      Data@vbt0[x]))))
    wa <- ac * la^bc
    
    Cat <- Data@Cat[x, ]
    Cat[1] <- Cat[2]  # temporary fix until effort simulation function gets sorted
    
    
    CAAv[, maxage][CAAv[, maxage] == 0] <- 1
    
    
    opt = optim(c(-3, -2), VPAopt, Cat = CAAv, yt = Data@Ind[x, 
      ], S = exp(-Mc), maxage = maxage, wa = wa, pmat = pmat, method = "L-BFGS-B", 
      lower = c(-5, -5), upper = c(5, 5))
    out = VPAopt(opt$par, Cat = CAAv, yt = Data@Ind[x, ], S = exp(-Data@Mort[x]), 
      maxage = maxage, wa = wa, pmat = pmat, opt = F)
    
    fit2 <- optimize(VPAFMSY, log(c(1e-04, 3)), Mc = Mc, hc = hc, maxage = maxage, 
      vul = out$va, Linfc = Linfc, Kc = Kc, t0c = t0c, AMc = AMc, 
      ac = ac, bc = bc)
    FMSY <- VPAFMSY(fit2$minimum, Mc = Mc, hc = hc, , maxage = maxage, 
      vul = out$va, Linfc = Linfc, Kc = Kc, t0c = t0c, AMc = AMc, 
      ac = ac, bc = bc, opt = F)
    if ((FMSY/Mc) > 3) 
      FMSY <- 3 * Mc
    TAC[i] <- out$bt[length(out$bt)] * FMSY
  }
  
  TACfilter(TAC)
  
}
class(VPA) <- "Output"

# VPAFMSY<-function(lnFMc,Mc,hc,maxage,vul,Linfc,Kc,t0c,AMc,ac,bc,opt=T,ny=50){

VPAopt = function(theta, Cat, yt, S, maxage, wa, pmat, opt = T, usewat = F) {
  
  Uterm <- exp(theta[1])/(1 + exp(theta[1]))
  minagecom = 1
  minagesel = ceiling(maxage * 0.66)
  avg_yrsel <- max(2, (min(5, floor(dim(Cat)[1]/2))))
  
  sig = exp(theta[2])
  tiny = 1e-10
  n = dim(Cat)[1]
  A = dim(Cat)[2]
  Nat = matrix(NA, n, A)  # Numbers-at-age matrix
  Ut = rep(NA, length = n)
  ai = 1:(A - 2)
  va = c(plogis(1:(A - 4), 2, 0.2), rep(1, 4))  # Initial values at the terminal selectivy
  Ut[n] = Uterm
  
  for (j in 1:15) {
    # Numerical convergence to terminal F print(Ut)
    Nat[n, ] = Cat[n, ]/(Uterm * va)  # Initialize the terminal year
    
    for (i in (n - 1):1) {
      Nat[i, ai] = Nat[i + 1, ai + 1]/S + Cat[i, ai]
      Nat[i, A - 1] = (Nat[i + 1, A]/S + Cat[i, A - 1] + Cat[i, A]) * 
        (Cat[i, A - 1]/(Cat[i, A - 1] + Cat[i, A] + tiny))
      Nat[i, A] = (Nat[i + 1, A]/S + Cat[i, A - 1] + Cat[i, A]) * 
        (Cat[i, A]/(Cat[i, A - 1] + Cat[i, A] + tiny))
    }
    ############################################## modify this parameters if need it #####
    
    # minagesel = 8
    minagecom = 1
    
    Ut = rowSums(Cat[, minagesel:(A - 1)])/rowSums(Nat[, minagesel:(A - 
      1)])  # Exploitation rate for fully recruited fish
    # Ut[n] = 0.4
    vat = Cat/Nat/Ut  # relative vulnerablility at age
    va = colMeans(vat[(n - avg_yrsel):(n - minagecom), ])  # update terminal vul
    va[minagesel:A] = 1
    
  }
  
  vat[is.na(vat)] <- 1
  Ut[n] = Uterm
  if (usewat == T) 
    vbt = rowSums((Nat * vat) * wa) else vbt = (Nat * vat) %*% wa
  if (usewat == T) 
    bt = as.vector(rowSums(Nat * wa)) else bt = as.vector(Nat %*% wa)
  fec = pmat * wa
  zt = log(yt/vbt)
  epsilon = zt - mean(zt)
  if (usewat == T) 
    ssb = as.vector(rowSums(Nat * fec)) else ssb = as.vector(Nat %*% fec)
  predcpue = exp(mean(zt)) * vbt  ### check again if bt or vbt
  cpue_q = yt/exp(mean(zt))
  qhat = exp(mean(epsilon))
  
  lnl = sum(dnorm(epsilon, mean = 0, sd = sig, log = T))
  
  if (opt) {
    return(-lnl)
  } else {
    # ss = sum(epsilon^2) lnl = 0.5*n*log(ss)
    return(list(Uterm = Uterm, va = va, rt = Nat[, 1], ssb = ssb, yt = yt, 
      vbt = vbt, cpue_q = cpue_q, Nat = Nat, vat = vat, Ut = Ut, 
      bt = bt, predcpue = predcpue, epsilon = epsilon/sig, lnl = lnl, 
      qhat = qhat, minagesel = minagesel, minagecom = minagecom, 
      avg_yrsel = avg_yrsel))
  }
  
  
}

VPAFMSY <- function(lnFMc, Mc, hc, maxage, vul, Linfc, Kc, t0c, AMc, ac, 
  bc, opt = T, ny = 50) {
  
  FMc <- exp(lnFMc)
  
  Mac <- rep(1, maxage)
  Mac[1:max(1, floor(AMc))] <- 0
  Lac <- Linfc * (1 - exp(-Kc * ((1:maxage) - t0c)))
  Wac <- ac * Lac^bc
  R0c <- 1
  N <- exp(-Mc * ((1:maxage) - 1)) * R0c
  SSN <- Mac * N  # Calculate initial spawning stock numbers
  Biomass <- N * Wac
  SSB <- SSN * Wac  # Calculate spawning stock biomass
  
  B0 <- sum(Biomass)
  SSB0 <- sum(SSB)
  SSN0 <- SSN
  SSBpR <- sum(SSB)/R0c  # Calculate spawning stock biomass per recruit
  SSNpR <- SSN/R0c
  
  N <- N/2
  SSN <- Mac * N  # Calculate initial spawning stock numbers
  Biomass <- N * Wac
  SSB <- SSN * Wac
  
  for (y in 1:ny) {
    # set up some indices for indexed calculation Fishing mortality rate
    # determined by effort, catchability, vulnerability and spatial
    # preference according to biomass
    Zc <- FMc * vul + Mc
    CN <- N * (1 - exp(-Zc)) * (FMc/Zc)
    CB <- CN * Wac
    Biomass <- N * Wac
    N[2:maxage] <- N[1:(maxage - 1)] * exp(-Zc[1:(maxage - 1)])  # Total mortality
    N[1] <- (0.8 * R0c * hc * sum(SSB))/(0.2 * SSBpR * R0c * (1 - hc) + 
      (hc - 0.2) * sum(SSB))  # Recruitment assuming regional R0 and stock wide steepness
    # print(N[1])
    SSN <- N * Mac
    SSB <- SSN * Wac
    
  }  # end of year
  
  if (opt) {
    return(-sum(CB))
  } else {
    return(FMc)
  }
}

# A reference MP with zero(ish) catch


#' No Fishing Reference MP
#' 
#' A reference MP that sets annual catch to zero (or very close to it). Used
#' for looking at variability in stock with no fishing.
#' 
#' 
#' @usage NFref(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the quota recommendation
#' @return A TAC of 0.01
#' @author A. Hordyk
#' @export NFref
NFref <- function(x, Data, reps = 100) {
  rep(0.01, reps)
}
class(NFref) <- "Output"



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
  TACfilter(TAC)
}
class(L95target) <- "Output"


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
  
  TACfilter(TAC)
}
class(Iratio) <- "Output"

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
  TACfilter(TAC)
}
class(ICI) <- "Output"



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
  TACfilter(TAC)
}
class(ICI2) <- "Output"


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
  TACfilter(TAC)
}
class(Lratio_BHI) <- "Output"


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
  TACfilter(TAC)
}
class(Lratio_BHI2) <- "Output"





# SCA<-function(x,Data,reps=100){ # Requires a character string
# DLMexe (e.g. 'C:/DLMexe') that represents the

# dependencies=''

# ny<-dim(Data@CAA)[2] na<-dim(Data@CAA)[3]

# write starter file
# ---------------------------------------------------------------------------------------------------------------------------------

# starterfile=paste(DLMexe,'SCA/Starter.ss',sep='/')
# write('SCA.dat',starterfile,1,append=F)
# write('SCA.ctl',starterfile,1,append=T)
# write(0,starterfile,1,append=T) # use init values
# write(1,starterfile,1,append=T) # run display detail (0,1,2)
# write(1,starterfile,1,append=T) # detailed age-structured reports in
# REPORT.SSO (0,1)

# write(0,starterfile,1,append=T) # write detailed checkup.sso file
# (0,1) write(1,starterfile,1,append=T) # write parm values to
# ParmTrace.sso write(0,starterfile,1,append=T) # write to
# cumreport.sso (0=no, 1=like$tiemseries; 2=add survey fits)

# write(1,starterfile,1,append=T) # include prior like for
# non-estimated parameters (0,1) write(0,starterfile,1,append=T) # use
# soft boundaries to aid convergence (0,1) (recommended)
# write(3,starterfile,1,append=T) # Number of data files to produce:
# 1st is input, 2nd is estimates, 3rd and higher are bootstrap

# write(10,starterfile,1,append=T) # Turn off estimation for parameters
# entering after this phase write(10,starterfile,1,append=T) # MCeval
# burn interval write(2,starterfile,1,append=T) # MCeval thin interval

# write(0,starterfile,1,append=T) # jitter initial parm value by this
# fraction write(-1,starterfile,1,append=T) # min yr for sdreport
# outputs (-1 for styr) write(-2,starterfile,1,append=T) # max yr for
# sdreport outputs

# write(0,starterfile,1,append=T) # N individual STD years
# write(0.001,starterfile,1,append=T) # final convergence criteria
# (e.g. 1.0e-04) write(0,starterfile,1,append=T) # retrospective year
# relative to end year (e.g. -4)

# write(1,starterfile,1,append=T) # min age for calc of summary biomass
# write(2,starterfile,1,append=T) # Depletion basis: denom is: 0=skip;
# 1=rel X*B0; 2=rel X*Bmsy; 3=rel X*B_styr !!!!!!!!!!!!!!!!
# write('1.0',starterfile,1,append=T) # Fraction (X) for Depletion
# denominator (e.g. 0.4)

# write(2,starterfile,1,append=T) # SPR_report_basis: 0=skip;
# 1=(1-SPR)/(1-SPR_tgt); 2=(1-SPR)/(1-SPR_MSY);
# 3=(1-SPR)/(1-SPR_Btarget); 4=rawSPR write(4,starterfile,1,append=T) #
# F_report_units: 0=skip; 1=exploitation(Bio); 2=exploitation(Num);
# 3=sum(Frates); 4=true F for range of ages
# write(c(floor(max(1,Data@MaxAge/4)),
# ceiling(Data@MaxAge/3)),starterfile,2,append=T) #_min and max age
# over which average F will be calculated

# write(2,starterfile,1,append=T) # F_report_basis: 0=raw; 1=F/Fspr;
# 2=F/Fmsy ; 3=F/Fbtgt write(999,starterfile,1,append=T) # check value
# for end of file



# write control file
# ---------------------------------------------------------------------------------------------------------------------------------

# ctlfile=paste(DLMexe,'SCA/SCA.ctl',sep='/')

# write(1,ctlfile,1,append=F) #_N_Growth_Patterns
# write(1,ctlfile,1,append=T) #_N_Morphs_Within_GrowthPattern
# write(0,ctlfile,1,append=T) #_Nblock_Patterns
# write(0.5,ctlfile,1,append=T) #_fracfemale
# write(0,ctlfile,1,append=T) #_natM_type:_0=1Parm;
# 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate

# write(1,ctlfile,1,append=T) # GrowthModel: 1=vonBert with L1&L2;
# 2=Richards with L1&L2; 3=age_speciific_K; 4=not implemented
# write(0,ctlfile,1,append=T) #_Growth_Age_for_L1 !!!!!!!
# write(999,ctlfile,1,append=T) #_Growth_Age_for_L2 (999 to use as
# Linf)

# write(0,ctlfile,1,append=T) #_SD_add_to_LAA (set to 0.1 for SS2 V1.x
# compatibility) write(0,ctlfile,1,append=T) #_CV_Growth_Pattern: 0
# CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
# write(1,ctlfile,1,append=T) #_maturity_option: 1=length logistic;
# 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read
# age-fecundity; 5=read fec and wt from wtatage.ss

# write(1,ctlfile,1,append=T) #_First_Mature_Age !!!!!!!!!!!!!!!!!!!!
# write(1,ctlfile,1,append=T) #_fecundity
# option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b;
# (4)eggs=a+b*L; (5)eggs=a+b*W write(0,ctlfile,1,append=T)
# #_hermaphroditism option: 0=none; 1=age-specific fxn

# write(1,ctlfile,1,append=T) #_parameter_offset_approach (1=none, 2=
# M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
# !!!!!!!!!!!!!!!!!!!!  write(2,ctlfile,1,append=T)
# #_env/block/dev_adjust_method (1=standard; 2=logistic transform keeps
# in base parm bounds; 3=standard w/ no bound check)

# female write(paste(round(Data@Mort[x]*0.9,3),
# round(Data@Mort[x]*1.1,3),round(Data@Mort[x],3),round(Data@Mort[x],3),-1,0.1,-3,0,0,0,0,0,0,0,sep='
# '),ctlfile,1,append=T) # NatM_p_1_Fem_GP_1
# write(paste(0,0,0,0,-1,10,-3,0,0,0,0,0,0,0,sep='
# '),ctlfile,1,append=T) # L_at_Amin_Fem_GP_1
# write(paste(round(Data@vbLinf[x]*0.9,3),round(Data@vbLinf[x]*1.1,3),
# round(Data@vbLinf[x],3),round(Data@vbLinf[x],3),-1,0.1,-3,0,0,0,0,0,0,0,sep='
# '),ctlfile,1,append=T) # L_at_Amax_Fem_GP_1

# write(paste(round(Data@vbK[x]*0.9,3),round(Data@vbK[x]*1.1,3),
# round(Data@vbK[x],3),round(Data@vbK[x],3),-1,0.1,-3,0,0,0,0,0,0,0,sep='
# '),ctlfile,1,append=T) # VonBert_K_Fem_GP_1
# write(paste(0.05,0.25,0.1,0.1,-1,0.8,-3,0,0,0,0,0,0,0,sep='
# '),ctlfile,1,append=T) # CV_young_Fem_GP_1
# write(paste(0.05,0.25,0.1,0.1,-1,0.8,-3,0,0,0,0,0,0,0,sep='
# '),ctlfile,1,append=T) # CV_old_Fem_GP_1

# write(paste(Data@wla[x]*0.9,
# Data@wla[x]*1.1,Data@wla[x],Data@wla[x],-1,0.1,-3,0,0,0,0,0,0,0,sep='
# '),ctlfile,1,append=T) # Wtlen_1_Fem write(paste(Data@wlb[x]*0.9,
# Data@wlb[x]*1.1,Data@wlb[x],Data@wlb[x],-1,0.1,-3,0,0,0,0,0,0,0,sep='
# '),ctlfile,1,append=T) # Wtlen_2_Fem

# write(paste(round(Data@L50[x]*0.9,3),
# round(Data@L50[x]*1.1,3),round(Data@L50[x],3),round(Data@L50[x],3),-1,0.1,-3,0,0,0,0,0,0,0,sep='
# '),ctlfile,1,append=T) # Mat50%_Fem
# Lsd=-2.94439/(Data@L95[x]-Data@L50[x]) #
# slope=log(1/0.95-1)/(L95-L50)
# write(paste(-3,3,round(Lsd,3),round(Lsd,3),-1,0.8,-3,0,0,0,0,0,0,0,sep='
# '),ctlfile,1,append=T) # Mat_slope_Fem

# write(paste(Data@wla[x]*0.9,
# Data@wla[x]*1.1,Data@wla[x],Data@wla[x],-1,0.1,-3,0,0,0,0,0,0,0,sep='
# '),ctlfile,1,append=T) # Wtlen_1_Mal write(paste(Data@wlb[x]*0.9,
# Data@wlb[x]*1.1,Data@wlb[x],Data@wlb[x],-1,0.1,-3,0,0,0,0,0,0,0,sep='
# '),ctlfile,1,append=T) # Wtlen_2_Mal

# write(paste(0,3,1,1,-1,0.8,-3,0,0,0,0,0,0,0,sep='
# '),ctlfile,1,append=T) # RecrDist_Area_1
# write(paste(0,3,0,0,-1,0.8,-3,0,0,0,0,0,0,0,sep='
# '),ctlfile,1,append=T) # RecrDist_Seas_1

# write(paste(0,0,0,0,-1,0,-4,0,0,0,0,0,0,0,sep='
# '),ctlfile,1,append=T) # RecrDist_Area_1
# write(paste(0,0,0,0,-1,0,-4,0,0,0,0,0,0,0,sep='
# '),ctlfile,1,append=T) # RecrDist_Seas_1
# write(paste(0,0,0,0,-1,0,-4,0,0,0,0,0,0,0,sep='
# '),ctlfile,1,append=T) # CohortGrowDev

# write(paste(0.5,1.5,1,1,-1,0,-4,0,0,0,0,0,0,0,sep='
# '),ctlfile,1,append=T) # CohortGrowDev


# write(paste(0,0,0,0,0,0,0,0,0,0,sep=' '),ctlfile,1,append=T)
# #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K


# Reconstruct a possible catch-at-age matrix to get to R0 given M
# la<-Data@vbLinf[x]*(1-exp(-Data@vbK[x]*(((1:na)-Data@vbt0[x]))))
# wa<-Data@wla[x]*la^Data@wlb[x]
# Cw<-t(array(wa,c(na,ny)))*Data@CAA[x,,] # weight of the observed
# CAA Cwtot<-apply(Cw,1,sum) # summation by year
# CAAup<-Data@CAA[x,,]*Data@Cat[x,]/Cwtot # uprate CAA to total
# catch weight c1<-apply(CAAup,2,mean)
# plusgroup<-which.min((cumsum(c1)/sum(c1)-0.95)^2)
# c1[plusgroup]<-c1[plusgroup]+sum(c1[(plusgroup+1):na])
# c1<-c1[plusgroup:1] aa<-rep(NA,plusgroup) aa[1]<-c1[1] for(i in
# 2:plusgroup)aa[i]=aa[i-1]/exp(-Data@Mort[x])+c1[i]
# R0est<-aa[plusgroup] _Spawner-Recruitment write(3,ctlfile,1,append=T)
# #_SR_function: 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop;
# 7=survival_3Parm _LO HI INIT PRIOR PR_type SD PHASE

# write(paste(round(log(R0est/10),3),round(log(R0est*10),3),round(log(R0est),3),round(log(R0est),3),-1,10,1,sep='
# '),ctlfile,1,append=T) # SR_LN(R0)
# write(paste(round(Data@steep[x]*0.9,3),round(Data@steep[x],3),
# round(Data@steep[x],3),round(Data@steep[x],3),1,0.04,-3,sep='
# '),ctlfile,1,append=T) # SR_BH_steep
# write(paste(0,2,0.6,0.8,-1,0.8,-4,sep=' '),ctlfile,1,append=T) #
# SR_sigmaR

# write(paste(-5,5,0,0,-1,1,-3,sep=' '),ctlfile,1,append=T) #
# SR_envlink write(paste(-5,5,0,0,-1,1,-2,sep=' '),ctlfile,1,append=T)
# # SR_R1_offset write(paste(0,0,0,0,-1,0,-99,sep='
# '),ctlfile,1,append=T) # SR_autocorr

# write(0,ctlfile,1,append=T) #_SR_env_link write(0,ctlfile,1,append=T)
# #_SR_env_target_0=none;1=devs;_2=R0;_3=steepness
# write(1,ctlfile,1,append=T) # do_recdev write(1,ctlfile,1,append=T) #
# first year of main recr_devs; early devs can preceed this era
# write(ny,ctlfile,1,append=T) # last year of main recr_devs; forecast
# devs start in following year write(2,ctlfile,1,append=T) #_recdev
# phase

# write(1,ctlfile,1,append=T) # (0/1) to read 13 advanced options

# write(-na+2,ctlfile,1,append=T) # write(4,ctlfile,1,append=T) #
# write(0,ctlfile,1,append=T) # write(1,ctlfile,1,append=T) #
# write(1,ctlfile,1,append=T) write(10,ctlfile,1,append=T) #
# write(ny-8,ctlfile,1,append=T) # write(ny-1,ctlfile,1,append=T) #
# write(0.8,ctlfile,1,append=T) # write(0,ctlfile,1,append=T) #
# write(-5,ctlfile,1,append=T) # write(5,ctlfile,1,append=T) #
# write(0,ctlfile,1,append=T) #


# write(Data@Mort[x],ctlfile,1,append=T) # F ballpark for tuning
# early phases write(-ny,ctlfile,1,append=T) # F ballpark year (neg
# value to disable) write(3,ctlfile,1,append=T) # F_Method: 1=Pope;
# 2=instan. F; 3=hybrid (hybrid is recommended)
# write(2.9,ctlfile,1,append=T) # max F or harvest rate, depends on
# F_Method write(4,ctlfile,1,append=T) # N iterations for tuning F in
# hybrid method (recommend 3 to 7)

# _LO HI INIT PRIOR PR_type SD PHASE
# write(paste(0,round(Data@Mort[x]*3,3),0,0.01,0,99,-1,sep='
# '),ctlfile,1,append=T) #InitF_1FISHERY1

# _Den-dep env-var extra_se Q_type write(paste(0,0,0,0,sep='
# '),ctlfile,1,append=T) # FISHERY write(paste(0,0,0,2,sep='
# '),ctlfile,1,append=T) # SURVEY

# LOq<-1/(mean(Data@Cat[x,])/(Data@Mort[x]/10))
# HIq<-1/(mean(Data@Cat[x,])/(Data@Mort[x]*5))
# muq<-1/(mean(Data@Cat[x,])/Data@Mort[x]) LO HI INIT PRIOR
# PR_type SD PHASE
# write(paste(round(log(LOq),3),round(log(HIq),3),round(log(muq),3),round(log(muq),3),0,1,-1,sep='
# '),ctlfile,1,append=T) # Q_base_SURVEY

# _Pattern Discard Male Special write(paste(0,0,0,0,sep='
# '),ctlfile,1,append=T) # FISHERY write(paste(0,0,0,0,sep='
# '),ctlfile,1,append=T) # SURVEY

# _age_selex_types write(paste(12,0,0,0,sep=' '),ctlfile,1,append=T) #
# FISHERY write(paste(12,0,0,0,sep=' '),ctlfile,1,append=T) # SURVEY

# LO HI INIT PRIOR PR_type SD PHASE env-var use_dev dev_minyr dev_maxyr
# dev_stddev Block Block_Fxn
# write(paste(1,floor(na*0.666),5,5,-1,10,1,0,0,0,0,0,0,0,sep='
# '),ctlfile,1,append=T) # AgeSel 1 Fishery write(paste(0,0.5,
# 0.25,0.25,-1,10,1,0,0,0,0,0,0,0,sep=' '),ctlfile,1,append=T) # AgeSel
# 2 Fishery
# write(paste(1,floor(na*0.666),5,5,-1,10,1,0,0,0,0,0,0,0,sep='
# '),ctlfile,1,append=T) # AgeSel 1 Fishery write(paste(0,0.5,
# 0.25,0.25,-1,10,1,0,0,0,0,0,0,0,sep=' '),ctlfile,1,append=T) # AgeSel
# 2 Fishery

# write(paste(0,floor(na*0.666),floor(na/3),floor(na/3),0,2,2,0,0,0,0,0.5,0,0,sep='
# '),ctlfile,1,append=T) # AgeSel 1 SURVEY
# write(paste(0.01,floor(na*0.8),floor(na/2),floor(na/2),0,2,2,0,0,0,0,0.5,0,0,sep='
# '),ctlfile,1,append=T) # AgeSel 1 SURVEY

# Tag loss and Tag reporting parameters go next
# write(0,ctlfile,1,append=T) # TG_custom: 0=no read; 1=read if tags
# exist

# write(0,ctlfile,1,append=T)#_Variance_adjustments_to_input_values
# _fleet: 1 2 3 write(0,ctlfile,1,append=T) #_add_to_survey_CV
# write(0,ctlfile,1,append=T) #_add_to_discard_stddev
# write(0,ctlfile,1,append=T) #_add_to_bodywt_CV
# write(1,ctlfile,1,append=T) #_mult_by_lencomp_N
# write(1,ctlfile,1,append=T) #_mult_by_agecomp_N
# write(1,ctlfile,1,append=T) #_mult_by_size-at-age_N

# write(2,ctlfile,1,append=T) #_maxlambdaphase
# write(1,ctlfile,1,append=T) #_sd_offset write(1,ctlfile,1,append=T) #
# number of changes to make to default Lambdas (default value is 1.0)
# write(paste(1,2,1,0,1,sep=' '),ctlfile,1,append=T) # Like_comp SURVEY

# 1 2 2 1 1 4 2 2 1 1 4 2 3 1 1

# write(0,ctlfile,1,append=T) # (0/1) read specs for more stddev
# reporting write(999,ctlfile,1,append=T) # end-of-file



# write data file
# ---------------------------------------------------------------------------------------------------------------------------------

# datfile=paste(DLMexe,'SCA/SCA.dat',sep='/')

# write(1,datfile,1,append=F) #_styr write(ny,datfile,1,append=T)
# #_endyr write(1,datfile,1,append=T) #_nseas
# write(12,datfile,1,append=T) # months/season
# write(1,datfile,1,append=T) #_spawn_seas write(1,datfile,1,append=T)
# #_Nfleet write(1,datfile,1,append=T) #_Nsurveys
# write(1,datfile,1,append=T) #_N_areas
# write('FISHERY%SURVEY',datfile,1,append=T)
# write(paste(c(-1,0.5),collapse=' '),datfile,1,append=T)
# #_surveytiming_in_season write(paste(c(1,1),collapse='
# '),datfile,1,append=T) #_area_assignments_for_each_fishery_and_survey
# write(1,datfile,1,append=T) #_units of catch: 1=bio; 2=num
# write(0.1,datfile,1,append=T) #_se of log(catch) only used for
# init_eq_catch and for Fmethod 2 and 3; use -1 for discard only fleets
# write(1,datfile,1,append=T) #_Ngenders write(na,datfile,1,append=T)
# #_Nages write(0,datfile,1,append=T)
# #_init_equil_catch_for_each_fishery write(ny,datfile,1,append=T)
# #_N_lines_of_catch_to_read
# _catch_biomass(mtons):_columns_are_fisheries,year,season for(i in
# 1:ny) write(paste(round(Data@Cat[x,i],2),i,1,sep='
# '),datfile,1,append=T)

# write(ny,datfile,1,append=T) #_N_cpue_and_surveyabundance_observation
# _Units: 0=numbers; 1=biomass; 2=F _Errtype: -1=normal; 0=lognormal;
# >0=T _Fleet Units Errtype write(paste(1,1,0,sep='
# '),datfile,1,append=T) # SURVEY write(paste(2,1,0,sep='
# '),datfile,1,append=T) # FISHERY

# _year seas index obs err for(i in 1:ny)
# write(paste(i,1,2,round(Data@Ind[x,i],2),0.2,sep='
# '),datfile,1,append=T) # index 2 is SURVEY

# write(0,datfile,1,append=T) #_N_fleets_with_discard
# write(0,datfile,1,append=T) # N discard obs
# write(0,datfile,1,append=T) #_N_meanbodywt_obs
# write(30,datfile,1,append=T) #_DF_for_meanbodywt_T-distribution_like

# write(2,datfile,1,append=T) # length bin method: 1=use databins;
# 2=generate from binwidth,min,max below; 3=read vector
# write(2,datfile,1,append=T) # binwidth for population size comp
# write(Data@CAL_bins[2],datfile,1,append=T) # minimum size in the
# population (lower edge of first bin and size at age 0.00)
# write(Data@CAL_bins[length(Data@CAL_bins)],datfile,1,append=T)
# # maximum size in the population (lower edge of last bin)

# write(0.0001,datfile,1,append=T) #_comp_tail_compression
# write(1e-007,datfile,1,append=T) #_add_to_comp
# write(0,datfile,1,append=T) #_combine males into females at or below
# this bin number write(length(Data@CAL_bins),datfile,1,append=T)
# #_N_LengthBins write(paste(Data@CAL_bins,collapse='
# '),datfile,1,append=T) # Length bins write(0,datfile,1,append=T)
# #_N_Length_obs


# maxage<-which.min(abs(cumsum(apply(Data@CAA[x,,],2,sum))/sum(Data@CAA[x,,])-0.95))
# # plus group at cumulative 75th percentile
# CAA<-Data@CAA[x,,1:maxage]
# CAA[,maxage]<-CAA[,maxage]+apply(Data@CAA[x,,maxage:na],1,sum)

# write(maxage,datfile,1,append=T) #_N_age_bins
# write(paste(1:maxage,collapse=' '),datfile,1,append=T) #age_bins

# write(1,datfile,1,append=T) #_N_ageerror_definitions
# write(paste(rep(-1,na+1),collapse=' '),datfile,1,append=T)
# write(paste(rep(0.001,na+1),collapse=' '),datfile,1,append=T)

# write(ny*2,datfile,1,append=T) #_N_Agecomp_obs
# write(1,datfile,1,append=T) #_Lbin_method: 1=poplenbins;
# 2=datalenbins; 3=lengths write(1,datfile,1,append=T) #_combine males
# into females at or below this bin number

# Yr Seas Flt/Svy Gender Part Ageerr Lbin_lo Lbin_hi Nsamp
# datavector(female-male)

# for(i in 1:ny) write(paste('
# ',i,1,1,0,0,1,-1,-1,sum(CAA[i,]),paste(CAA[i,],collapse=' '),sep='
# '),datfile,1,append=T) for(i in 1:ny) write(paste('
# ',i,1,2,0,0,1,-1,-1,sum(CAA[i,]),paste(CAA[i,],collapse=' '),sep='
# '),datfile,1,append=T)

# write(0,datfile,1,append=T) #_N_MeanSize-at-Age_obs
# write(0,datfile,1,append=T) #_N_environ_variables
# write(0,datfile,1,append=T) #_N_environ_obs
# write(0,datfile,1,append=T) # N sizefreq methods to read
# write(0,datfile,1,append=T) # no tag data write(0,datfile,1,append=T)
# # no morphcomp data write(999,datfile,1,append=T) # end of file



# write forecast file
# ---------------------------------------------------------------------------------------------------------------------------------

# forefile=paste(DLMexe,'SCA/Forecast.ss',sep='/')

# write(1,forefile,1,append=F) # Benchmarks: 0=skip; 1=calc
# F_spr,F_btgt,F_msy write(2,forefile,1,append=T) # MSY: 1= set to
# F(SPR); 2=calc F(MSY); 3=set to F(Btgt); 4=set to F(endyr)
# write(0.40,forefile,1,append=T) # SPR target (e.g. 0.40)
# write(0.35,forefile,1,append=T) # Biomass target (e.g. 0.40)
# _Bmark_years: beg_bio, end_bio, beg_selex, end_selex, beg_relF,
# end_relF (enter actual year, or values of 0 or -integer to be rel.
# endyr) write(paste(0,0,0,0,0,0,sep=' '),forefile,1,append=T) #
# FISHERY write(1,forefile,1,append=T) #Bmark_relF_Basis: 1 = use year
# range; 2 = set relF same as forecast below
# write(2,forefile,1,append=T) # Forecast: 0=none; 1=F(SPR); 2=F(MSY)
# 3=F(Btgt); 4=Ave F (uses first-last relF yrs); 5=input annual F
# scalar write(3,forefile,1,append=T) # N Forecast years
# write(0.2,forefile,1,append=T) # F scalar (only used for
# Do_Forecast==5) _Fcast_years: beg_selex, end_selex, beg_relF,
# end_relF (enter actual year, or values of 0 or -integer to be rel.
# endyr) write(paste(0,0,-10,0,sep=' '),forefile,1,append=T) # FISHERY
# write(1,forefile,1,append=T) # Control rule method (1=catch=f(SSB)
# west coast; 2=F=f(SSB) ) write(0.4,forefile,1,append=T) # Control
# rule Biomass level for constant F (as frac of Bzero, e.g. 0.40)
# write(0.1,forefile,1,append=T) # Control rule Biomass level for no F
# (as frac of Bzero, e.g. 0.10) write(0.75,forefile,1,append=T) #
# Control rule target as fraction of Flimit (e.g. 0.75)
# write(3,forefile,1,append=T) #_N forecast loops (1-3) (fixed at 3 for
# now) write(3,forefile,1,append=T) #_First forecast loop with
# stochastic recruitment write(0,forefile,1,append=T) #_Forecast loop
# control #3 (reserved for future bells&whistles)
# write(0,forefile,1,append=T) #_Forecast loop control #4 (reserved for
# future bells&whistles) write(0,forefile,1,append=T) #_Forecast loop
# control #5 (reserved for future bells&whistles)
# write(ny,forefile,1,append=T) #FirstYear for caps and allocations
# (should be after years with fixed inputs)
# write(0.0,forefile,1,append=T) # stddev of log(realized catch/target
# catch) in forecast (set value>0.0 to cause active impl_error)
# write(0,forefile,1,append=T) # Do West Coast gfish rebuilder output
# (0/1) write(ny-10,forefile,1,append=T) # Rebuilder: first year catch
# could have been set to zero (Ydecl)(-1 to set to 1999)
# write(ny-5,forefile,1,append=T) # Rebuilder: year for current age
# structure (Yinit) (-1 to set to endyear+1)
# write(1,forefile,1,append=T) # fleet relative F: 1=use first-last
# alloc year; 2=read seas(row) x fleet(col) below
# write(2,forefile,1,append=T) # basis for fcast catch tuning and for
# fcast catch caps and allocation (2=deadbio; 3=retainbio; 5=deadnum;
# 6=retainnum) write(-1,forefile,1,append=T) # max totalcatch by fleet
# (-1 to have no max) write(-1,forefile,1,append=T) # max totalcatch by
# area (-1 to have no max) write(0,forefile,1,append=T) # fleet
# assignment to allocation group (enter group ID# for each fleet, 0 for
# not included in an alloc group) write(0,forefile,1,append=T) # Number
# of forecast catch levels to input (else calc catch from forecast F)
# write(2,forefile,1,append=T) # basis for input Fcast catch: 2=dead
# catch; 3=retained catch; 99=input Hrate(F) (units are from
# fleetunits; note new codes in SSV3.20) write(999,forefile,1,append=T)
# # End of input



# Run SS3 and read outputs
# ----------------------------------------------------------------------------------------------------------------------------

# system(paste(DLMexe,'/SCA/SS3.exe',sep=''),wait=T,show.output.on.console=T)
# rl <- SS_output(dir=paste(DLMexe,'/SCA/',sep=''))
# F_FMSY<-rl$Kobe$F.Fmsy[1:ny] B_BMSY<-rl$Kobe$B.Bmsy[1:ny]
# ind<-match(paste('F_',1:ny,sep=''),row.names(rl$derived_quants))
# Fs<-rl$derived_quants[ind,2] Fs/F_FMSY

# SS_plots(replist=myreplist)



# } class(SCA)<-'Output'




