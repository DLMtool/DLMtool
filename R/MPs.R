

## Output Control MPs ####
#' Average Catch NEW
#' 
#' A simple average catch MP that is included to demonstrate a 'status quo' management option
#' 
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @author T. Carruthers
#' @export 
#' 
AvCNEW <- function(x, Data, reps = 100) {
  dependencies = "Data@Cat"
  Rec <- new("Rec")
  Rec@TAC <- rlnorm(reps, log(mean(Data@Cat[x, ], na.rm = T)), 0.2)
  Rec
}
class(AvCNEW) <- "MP"


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
matlenlimNEW <- function(x, Data, ...) {
  # Knife-edge vulnerability at estimated length-at-maturity  
  dependencies = "Data@L50"
  
  rec <- new("Rec") # create recommendation object
  rec@LR5 <- Data@L50[x] * 0.95 # new length at 5% retention  
  rec@LFR <-  Data@L50[x] # new length at full retention   
  
  # other slots aren't specified so remain unchanged
  rec
}
class(matlenlimNEW) <- "MP"


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
MRrealNEW <- function(x, Data, ...) {
  # A Marine reserve in area 1 with spatial reallocation of effort
  
  rec <- new("Rec") # create recommendation object
  rec@Allocate <- 1
  rec@Spatial <- c(0,1)
  
  # other slots aren't specified so remain unchanged
  return(rec)
}
class(MRrealNEW) <- "MP"

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
curENEW <- function(x, Data, ...) {
  # current effort
  rec <- new("Rec") # create recommendation object
  rec@Effort <- 1
  rec
}
class(curENEW) <- "MP"

## Combined Management MPs ####

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
FMSYrefNEW <- function(x, Data, reps = 100) {
  rec <- new("Rec") # create recommendation object
  rec@Spatial <- c(0,1)
  rec@TAC <- trlnorm(reps, Data@OM$A[x] * (1 - exp(-Data@OM$FMSY[x])), 0.01)
  rec
  
}
class(FMSYrefNEW) <- "MP"

