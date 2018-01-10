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
FMSYref50 <- function(x, Data, reps = 100) {
  rec <- new("Rec") # create recommendation object
  rec@TAC <- trlnorm(reps, Data@OM$A[x] * (1 - exp(-Data@OM$FMSY[x]*0.5)) , 0.01)
  rec
}
class(FMSYref50) <- "MP"



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
FMSYref75 <- function(x, Data, reps = 100) {
  rec <- new("Rec") # create recommendation object
  rec@TAC <- trlnorm(reps, Data@OM$A[x] * (1 - exp(-Data@OM$FMSY[x]*0.75)) , 0.01)
  rec
}
class(FMSYref75) <- "MP"

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
  rec <- new("Rec") # create recommendation object
  rec@TAC <- rep(0.01, reps)
  rec
}
class(NFref) <- "MP"