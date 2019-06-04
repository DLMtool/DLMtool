## Reference MPs ####

#' Reference management procedures
#' 
#' Several reference MPs for your operating model to use in the management strategy
#' evaluation. FMSYref (and related) assume perfect information about FMSY (FMSY 
#' is taken from the operating model stored at Data@@OM$FMSY). NFref sets annual catch 
#' to zero (or close to it) and is used for looking at variability in stock with no fishing.
#'
#' @templateVar mp FMSYref 
#' @template MPtemplate
#' 
#' @details 
#' Note that you can out-perform \code{FMSYref} easily. The requirement for fixed
#' F is actually quite strict and is by no means the upper limit in terms of
#' yield. Don't panic if your method beats this one for yield, especially for
#' short-lived species of high temporal variability in productivity!
#' 
#' @author T. Carruthers, A. Hordyk
#' @describeIn FMSYref A reference FMSY method that fishes at FMSY
#' @examples 
#' FMSYref(1, DLMtool::SimulatedData, plot=TRUE)
#' @export 
FMSYref <- function(x, Data, reps = 100, plot=FALSE) {
  rec <- new("Rec") # create recommendation object
  rec@TAC <- trlnorm(reps, Data@OM$A[x] * (1 - exp(-Data@OM$FMSY[x])), 0.01)
  if (plot) boxplot(rec@TAC, ylab=paste0("TAC (", Data@Units, ")"))
  rec
  
}
class(FMSYref) <- "MP"

#' @describeIn FMSYref A reference FMSY method that fishes at 50\% of FMSY
#' @examples 
#' FMSYref50(1, DLMtool::SimulatedData, plot=TRUE)
#' @export  
FMSYref50 <- function(x, Data, reps = 100, plot=FALSE) {
  rec <- new("Rec") # create recommendation object
  rec@TAC <- trlnorm(reps, Data@OM$A[x] * (1 - exp(-Data@OM$FMSY[x]*0.5)) , 0.01)
  if (plot) boxplot(rec@TAC, ylab=paste0("TAC (", Data@Units, ")"))
  rec
}
class(FMSYref50) <- "MP"

#' @describeIn FMSYref A reference FMSY method that fishes at 75\% of FMSY
#' @examples 
#' FMSYref75(1, DLMtool::SimulatedData, plot=TRUE)
#' @export 
FMSYref75 <- function(x, Data, reps = 100, plot=FALSE) {
  rec <- new("Rec") # create recommendation object
  rec@TAC <- trlnorm(reps, Data@OM$A[x] * (1 - exp(-Data@OM$FMSY[x]*0.75)) , 0.01)
  if (plot) boxplot(rec@TAC, ylab=paste0("TAC (", Data@Units, ")"))
  rec
}
class(FMSYref75) <- "MP"


#' @describeIn FMSYref A reference MP that sets annual catch to almost zero (0.01)
#' @examples 
#' NFref(1, DLMtool::SimulatedData, plot=TRUE)
#' @export 
NFref <- function(x, Data, reps = 100, plot=FALSE) {
  rec <- new("Rec") # create recommendation object
  rec@TAC <- rep(0.01, reps)
  if (plot) boxplot(rec@TAC, ylab=paste0("TAC (", Data@Units, ")"))
  rec
}
class(NFref) <- "MP"

