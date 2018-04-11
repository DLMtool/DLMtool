## Reference MPs ####

#' Reference management procedures
#' 
#' Several reference MPs for your operating model to use in the management strategy
#' evaluation. FMSYref (and related) assume perfect information about FMSY (FMSY 
#' is taken from the operating model stored at Data@@OM$FMSY). NFref sets annual catch 
#' to zero (or close to it). Used for looking at variability in stock with no fishing.
#'
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples 
#' @details 
#' Note that you can out-perform \code{FMSYref} easily. The requirement for fixed
#' F is actually quite strict and is by no means the upper limit in terms of
#' yield. Don't panic if your method beats this one for yield, especially for
#' short-lived species of high temporal variability in productivity!
#' 
#' \code{NFref} sets annual catch to 0.01.
#' @author T. Carruthers, A. Hordyk
#' @describeIn FMSYref A reference FMSY method that fishes at FMSY
#' @export 
FMSYref <- function(x, Data, reps = 100) {
  rec <- new("Rec") # create recommendation object
  rec@TAC <- trlnorm(reps, Data@OM$A[x] * (1 - exp(-Data@OM$FMSY[x])), 0.01)
  rec
  
}
class(FMSYref) <- "MP"

#' @describeIn FMSYref A reference FMSY method that fishes at half of FMSY
#' @export FMSYref50 
FMSYref50 <- function(x, Data, reps = 100) {
  rec <- new("Rec") # create recommendation object
  rec@TAC <- trlnorm(reps, Data@OM$A[x] * (1 - exp(-Data@OM$FMSY[x]*0.5)) , 0.01)
  rec
}
class(FMSYref50) <- "MP"

#' @describeIn FMSYref A reference FMSY method that fishes at 75\% of FMSY
#' @export FMSYref75
FMSYref75 <- function(x, Data, reps = 100) {
  rec <- new("Rec") # create recommendation object
  rec@TAC <- trlnorm(reps, Data@OM$A[x] * (1 - exp(-Data@OM$FMSY[x]*0.75)) , 0.01)
  rec
}
class(FMSYref75) <- "MP"


#' @describeIn FMSYref A reference MP that sets annual catch to zero 
#' (or very close to it).
#' @export NFref
NFref <- function(x, Data, reps = 100) {
  rec <- new("Rec") # create recommendation object
  rec@TAC <- rep(0.01, reps)
  rec
}
class(NFref) <- "MP"
