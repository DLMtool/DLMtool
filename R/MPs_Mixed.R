
## Mixed Management MPs ####

#' A example mixed control MP that uses the Itarget1 output control MP together with a 
#' spatial closure.
#' 
#' Included demonstration purposes of a mixed control MP
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
#' 
#' @export 
Itarget1_MPA <- function(x, Data, reps = 100, yrsmth = 5, xx = 0, Imulti = 1.5) {
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
  Rec@Spatial <- c(0, rep(1, Data@nareas-1))
  Rec
}
class(Itarget1_MPA) <- "MP"


#' Average Catch with a size limit
#'
#' Included demonstration purposes of a mixed control MP
#'
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @export 
#'
AvC_MLL <- function(x, Data, reps = 100) {
  dependencies = "Data@Cat, L50"
  Rec <- new("Rec")
  Rec@TAC <- rlnorm(reps, log(mean(Data@Cat[x, ], na.rm = T)), 0.2)
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






