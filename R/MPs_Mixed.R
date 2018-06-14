
## Mixed Management MPs ####

#' Itarget1 with an MPA
#' 
#' A example mixed control MP that uses the Itarget1 output control MP together with a 
#' spatial closure. 
#' 
#' This MP has been included for demonstration purposes of a mixed control MP
#' 
#' @templateVar mp Itarget1_MPA
#' @template MPtemplate
#' @template MPuses
#' 
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param Imulti Parameter controlling how much larger target CPUE / index is
#' compared with recent levels.
#' 
#' @export 
#' @examples 
#' Itarget1_MPA(1, DLMtool::Atlantic_mackerel, plot=TRUE)
#' @family Index methods
Itarget1_MPA <- function(x, Data, reps = 100, plot=FALSE, yrsmth = 5, xx = 0, Imulti = 1.5) {
  runItarget <- Itarget_(x, Data, reps, plot, yrsmth, xx, Imulti)
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(runItarget$TAC)
  Rec@Spatial <- c(0, rep(1, Data@nareas-1))
  Rec
}
class(Itarget1_MPA) <- "MP"


#' Average Catch with a size limit
#' 
#' A example mixed control MP that uses the average catch output control MP together with a 
#' minimul size limit set at the size of maturity. 
#' 
#' This MP has been included for demonstration purposes of a mixed control MP
#' Included demonstration purposes of a mixed control MP
#'
#' @templateVar mp AvC_MLL
#' @template MPtemplate
#' @template MPuses
#' 
#' @export 
#' @family Average Catch MPs
#' 
#' @examples 
#' Rec <- AvC_MLL(1, DLMtool::Cobia, reps=1000, plot=TRUE) # 1,000 log-normal samples with CV = 0.2
AvC_MLL <- function(x, Data, reps = 100, plot=FALSE) {
  if (length(Data@Year)<1 | is.na(Data@LHYear)) {
    Rec <- new("Rec")
    Rec@TAC <- rep(as.numeric(NA), reps)
    return(Rec)
  }
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
  Rec@LR5 <- Data@L50[x] * 0.95 # new length at 5% retention  
  Rec@LFR <-  Data@L50[x] # new length at full retention 
  Rec
}
class(AvC_MLL) <- "MP"






#' Create an MP that averages the results of multiple MPs
#'
#' @param MPs A vector of MPs names 
#'
#' @return A function of class MP 
#' @export
#'
#' @examples
#' \dontrun{
#' MeanMP <- makeMeanMP(c("AvC", "DCAC"))
#' MSE <- runMSE(DLMtool::testOM, MPs=c("AvC", "DCAC", "MeanMP"))
#' Tplot2(MSE)
#' 
#' MeanMP <- makeMeanMP(c("matlenlim", "matlenlim2")) 
#' Data <- DLMtool::SimulatedData
#' matlenlim(1, Data)
#' matlenlim2(1, Data)
#' MeanMP(1, Data)
#' }
makeMeanMP <- function(MPs) {
  if (length(MPs)<2) stop("Must provide more than one MP")
  if (class(MPs) != 'character') stop("MPs must be a character vector")
  for (x in MPs) {
    if(class(get(x)) != "MP") stop(x, " is not class MP")
  }
  
  MP <- function(x, Data, reps=100) {
    nareas <- Data@nareas
    nMPs <- length(MPs)
    TAC <- matrix(NA, nrow=nMPs, ncol=reps)
    Effort <- matrix(NA, nrow=nMPs, ncol=reps)
    Spatial <- matrix(NA, nrow=nMPs, ncol=nareas)
    Allocate <- matrix(NA, nrow=nMPs, ncol=1)
    LR5 <- matrix(NA, nrow=nMPs, ncol=1)
    LFR <- matrix(NA, nrow=nMPs, ncol=1)
    HS <- matrix(NA, nrow=nMPs, ncol=1)
    Rmaxlen <- matrix(NA, nrow=nMPs, ncol=1) 
    L5  <- matrix(NA, nrow=nMPs, ncol=1)
    LFS <- matrix(NA, nrow=nMPs, ncol=1)
    Vmaxlen <- matrix(NA, nrow=nMPs, ncol=1)
    Fdisc  <- matrix(NA, nrow=nMPs, ncol=1)

    for (mm in seq_along(MPs)) {
      mod <- get(MPs[mm])
      rec <- mod(x, Data, reps)
      if (length(rec@TAC)>0) TAC[mm,] <- rec@TAC 
      if (length(rec@Effort)>0) Effort[mm,] <- rec@Effort 
      if (length(rec@Spatial)>0) Spatial[mm,] <- rec@Spatial 
      if (length(rec@Allocate)>0) Allocate[mm,] <- rec@Allocate 
      if (length(rec@LR5)>0) LR5[mm,] <- rec@LR5 
      if (length(rec@LFR)>0) LFR[mm,] <- rec@LFR 
      if (length(rec@HS)>0) HS[mm,] <- rec@HS 
      if (length(rec@Rmaxlen)>0) Rmaxlen[mm,] <- rec@Rmaxlen 
      if (length(rec@L5)>0) L5[mm,] <- rec@L5
      if (length(rec@LFS)>0) LFS[mm,] <- rec@LFS
      if (length(rec@Vmaxlen)>0) Vmaxlen[mm,] <- rec@Vmaxlen
      if (length(rec@Fdisc)>0) Fdisc[mm,] <- rec@Fdisc
    }
    
    rec <- new("Rec")
    rec@TAC <- apply(TAC, 2, mean, na.rm=TRUE)
    rec@Effort <- apply(Effort, 2, mean, na.rm=TRUE)
    rec@Spatial <- apply(Spatial, 2, mean, na.rm=TRUE)
    rec@Allocate <- apply(Allocate, 2, mean, na.rm=TRUE)
    rec@LR5 <- apply(LR5, 2, mean, na.rm=TRUE)
    rec@LFR <- apply(LFR, 2, mean, na.rm=TRUE)
    rec@HS <- apply(HS, 2, mean, na.rm=TRUE)
    rec@Rmaxlen <- apply(Rmaxlen, 2, mean, na.rm=TRUE)
    rec@L5 <- apply(L5, 2, mean, na.rm=TRUE)
    rec@LFS <- apply(LFS, 2, mean, na.rm=TRUE)
    rec@Vmaxlen <- apply(Vmaxlen, 2, mean, na.rm=TRUE)
    rec@Fdisc <- apply(Fdisc, 2, mean, na.rm=TRUE)
    for (sl in slotNames(rec)) {
      if (sl !="Misc")  if (all(!is.finite(slot(rec, sl)))) slot(rec, sl) <- numeric(0)
    }
    rec@Misc <- list(MPs)
    rec
  }
  class(MP) <- "MP"
  MP
}






