
#' Performance Metric: Probability B > 0.1 BMSY
#'
#' @param MSEobj An object of class MSE
#'
#' @return An object of class PMobj
#' @export
#'
#' @examples 
#' \dontrun{
#' P10(myMSE)
#' }
P10 <- function(MSEobj=NULL) {
  if (class(MSEobj)!='MSE') stop('Require object of class MSE')
  PMobj <- new("PMobj")
  PMobj@name <- "Probability Biomass above 10% BMSY" 
  PMobj@caption <- 'Prob. B > 0.1 BMSY'
  
  y.st <- 1 
  y.end <- MSEobj@proyears
  
  PMobj@Stat <- MSEobj@B_BMSY[,,y.st:y.end] # Performance Metric statistic of interest - here BMSY/FMSY 
  PMobj@Prob <- calcProb(PMobj@Stat > 0.1, MSEobj) # calculate probability Stat > 0.1 nsim by nMP
  
  PMobj@Mean <- calcMean(PMobj@Prob, MSEobj) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
  
}
class(P10) <- "PM"


#' Performance Metric: Probability B > 0.5 BMSY
#'
#' @param MSEobj An object of class MSE
#'
#' @return An object of class PMobj
#' @export
#'
#' @examples 
#' \dontrun{
#' P50(myMSE)
#' }
P50 <- function(MSEobj=NULL) {
  if (class(MSEobj)!='MSE') stop('Require object of class MSE')
  PMobj <- new("PMobj")
  PMobj@name <- "Probability Not Overfished" 
  PMobj@caption <- 'Prob. SB > 0.5 SBMSY'
  
  y.st <- 1 
  y.end <- MSEobj@proyears
  
  PMobj@Stat <- MSEobj@B_BMSY[,,y.st:y.end] # Performance Metric statistic of interest - here B/BMSY 
  PMobj@Prob <- calcProb(PMobj@Stat > 0.5, MSEobj) # calculate probability Stat > 0.5 nsim by nMP
  
  PMobj@Mean <- calcMean(PMobj@Prob, MSEobj) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
}
class(P50) <- "PM"

#' Performance Metric: Probability B > BMSY
#'
#' @param MSEobj An object of class MSE
#'
#' @return An object of class PMobj
#' @export
#'
#' @examples 
#' \dontrun{
#' P100(myMSE)
#' }
P100 <- function(MSEobj=NULL) {
  if (class(MSEobj)!='MSE') stop('Require object of class MSE')
  PMobj <- new("PMobj")
  PMobj@name <- "Probability Biomass > BMSY" 
  PMobj@caption <- 'Prob. SB > SBMSY'
  
  y.st <- 1 
  y.end <- MSEobj@proyears
  
  PMobj@Stat <- MSEobj@B_BMSY[,,y.st:y.end] # Performance Metric statistic of interest - here B/BMSY 
  PMobj@Prob <- calcProb(PMobj@Stat > 1, MSEobj) # calculate probability Stat > 1 nsim by nMP
  
  PMobj@Mean <- calcMean(PMobj@Prob, MSEobj) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
}
class(P100) <- "PM"

#' Performance Metric: Probability F < FMSY
#'
#' @param MSEobj An object of class MSE
#'
#' @return An object of class PMobj
#' @export
#'
#' @examples 
#' \dontrun{
#' POF(myMSE)
#' }
POF <- function(MSEobj=NULL) {
  if (class(MSEobj)!='MSE') stop('Require object of class MSE')
  PMobj <- new("PMobj")
  PMobj@name <- "Probability Not Overfishing" 
  PMobj@caption <- 'Prob. F < FMSY'
  
  y.st <- 1 
  y.end <- MSEobj@proyears
  
  PMobj@Stat <- MSEobj@F_FMSY[,,y.st:y.end] # Performance Metric statistic of interest - here F/FMSY 
  PMobj@Prob <- calcProb(PMobj@Stat < 1, MSEobj) # calculate probability Stat < 1 nsim by nMP
  
  
  PMobj@Mean <- calcMean(PMobj@Prob, MSEobj) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
  
}
class(POF) <- "PM"

#' Performance Metric: Probability Long-Term Yield > 0.5 Relative Yield
#'
#' @param MSEobj An object of class MSE
#'
#' @return An object of class PMobj
#' @export
#'
#' @examples 
#' \dontrun{
#' LTY(myMSE)
#' }
LTY <- function(MSEobj=NULL) {
  if (class(MSEobj)!='MSE') stop('Require object of class MSE')
  PMobj <- new("PMobj")
  PMobj@name <- "Average Long-Term Yield relative to Reference Yield" 
  PMobj@caption <- 'Prob. LTY > 0.5 Ref. Yield'
  
  y.st <- max(MSEobj@proyears-9,1) 
  y.end <- MSEobj@proyears # last 10 years
  
  RefYd <- array(MSEobj@OM$RefY, dim=dim(MSEobj@C[,,y.st:y.end]))
  
  PMobj@Stat <- MSEobj@C[,,y.st:y.end]/RefYd
  PMobj@Prob <- calcProb(PMobj@Stat > 0.5, MSEobj)  # probability LTY > 0.5 Ref Yield
  
  
  PMobj@Mean <- calcMean(PMobj@Prob, MSEobj) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
  
}
class(LTY) <- "PM"

#' Performance Metric: Probability Short-Term Yield > 0.5 Relative Yield
#'
#' @param MSEobj An object of class MSE
#'
#' @return An object of class PMobj
#' @export
#'
#' @examples 
#' \dontrun{
#' STY(myMSE)
#' }
STY <- function(MSEobj=NULL) {
  if (class(MSEobj)!='MSE') stop('Require object of class MSE')
  PMobj <- new("PMobj")
  PMobj@name <- "Average Short-Term Yield relative to Reference Yield" 
  PMobj@caption <- 'Prob. STY > 0.5 Ref. Yield'
  
  y.st <- 1
  y.end <- 10 # first 10 years
  
  RefYd <- array(MSEobj@OM$RefY, dim=dim(MSEobj@C[,,y.st:y.end]))
  
  PMobj@Stat <- MSEobj@C[,,y.st:y.end]/RefYd
  PMobj@Prob <- calcProb(PMobj@Stat > 0.5, MSEobj)  # probability STY > 0.5 Ref Yield
  
  PMobj@Mean <- calcMean(PMobj@Prob, MSEobj) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
  
}
class(STY) <- "PM"

#' Performance Metric: Probability AAVY < 0.2 
#'
#' @param MSEobj An object of class MSE
#'
#' @return An object of class PMobj
#' @export
#'
#' @examples 
#' \dontrun{
#' AAVY(myMSE)
#' }
AAVY <- function(MSEobj=NULL) {
  if (class(MSEobj)!='MSE') stop('Require object of class MSE')
  PMobj <- new("PMobj")
  PMobj@name <- "Average Annual Variability in Yield" 
  PMobj@caption <- 'Prob. AAVY < 0.2'
  
  y1<-1:(MSEobj@proyears-1) # year index
  y2<-2:MSEobj@proyears # 
  
  if (MSEobj@nMPs > 1) {
    AAVY <- apply(((MSEobj@C[,,y1]-MSEobj@C[,,y2])^2)^0.5,c(1,2),mean)/apply(MSEobj@C[,,y2],c(1,2),mean) 
  } else {
    AAVY <- array(apply(((MSEobj@C[,1,y1]-MSEobj@C[,1,y2])^2)^0.5,c(1),mean)/apply(MSEobj@C[,1,y2],c(1),mean))
  }
  
  PMobj@Stat <- AAVY
  PMobj@Prob <- calcProb(PMobj@Stat < 0.2, MSEobj)  # probability AAVY < 0.2 
  
  PMobj@Mean <- calcMean(PMobj@Prob, MSEobj) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
  
}
class(AAVY) <- "PM"

#' Performance Metric: Average Yield
#'
#' @param MSEobj An object of class MSE
#'
#' @return An object of class PMobj
#' @export
#'
#' @examples 
#' \dontrun{
#' Yield(myMSE)
#' }
Yield <- function(MSEobj=NULL) {
  if (class(MSEobj)!='MSE') stop('Require object of class MSE')
  PMobj <- new("PMobj")
  PMobj@name <- "Yield relative to Reference Yield" 
  PMobj@caption <- 'Mean Yield'
  
  y.st <- 1 
  y.end <- MSEobj@proyears 
  
  RefYd <- array(MSEobj@OM$RefY, dim=dim(MSEobj@C[,,y.st:y.end]))
  
  PMobj@Stat <- MSEobj@C[,,y.st:y.end]/RefYd
  PMobj@Prob <- calcProb(PMobj@Stat, MSEobj) # no probability to calculate
  
  PMobj@Mean <- calcMean(PMobj@Prob, MSEobj) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
  
}
class(Yield) <- "PM"


