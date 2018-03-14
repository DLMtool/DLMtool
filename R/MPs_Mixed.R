
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
