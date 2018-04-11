
#' Performance Metric: Probability B > 0.1 BMSY
#'
#' @param MSEobj An object of class MSE
#' @param Ref Reference Point for biomass (fraction of BMSY)
#' @return An object of class PMobj
#' @export
#'
#' @examples 
#' \dontrun{
#' P10(myMSE)
#' }
P10 <- function(MSEobj=NULL, Ref=0.1) {
  if (class(MSEobj)!='MSE') stop('Require object of class MSE')
  PMobj <- new("PMobj")
  PMobj@Name <- "Spawning Biomass relative to SBMSY"
  if (Ref !=1) {
    PMobj@Caption <- paste0('Prob. SB > ', Ref, ' SBMSY')
  } else {
    PMobj@Caption <- 'Prob. SB > SBMSY'
  }
  
  y.st <- 1 
  y.end <- MSEobj@proyears
  
  PMobj@Stat <- MSEobj@B_BMSY[,,y.st:y.end] # Performance Metric statistic of interest - here SB/SBMSY 
  PMobj@Ref <- Ref
  PMobj@Prob <- calcProb(PMobj@Stat > PMobj@Ref, MSEobj) # calculate probability Stat > 0.1 nsim by nMP
  
  PMobj@Mean <- calcMean(PMobj@Prob, MSEobj) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
  
}
class(P10) <- "PM"



#' Performance Metric: Probability B > 0.5 BMSY
#'
#' @param MSEobj An object of class MSE
#' @param Ref Reference Point for biomass (fraction of BMSY)
#' @return An object of class PMobj
#' @export
#'
#' @examples 
#' \dontrun{
#' P50(myMSE)
#' }
P50 <- P10
formals(P50)$Ref <- 0.5
class(P50) <- "PM"


#' Performance Metric: Probability B > BMSY
#'
#' @param MSEobj An object of class MSE
#' @param Ref Reference Point for biomass (fraction of BMSY)
#' @return An object of class PMobj
#' @export
#'
#' @examples 
#' \dontrun{
#' P100(myMSE)
#' }
P100 <- P10 
formals(P100)$Ref <- 1
class(P100) <- "PM"




#' Performance Metric: Probability F < FMSY
#'
#' @param MSEobj An object of class MSE
#' @param Ref Reference Point for fishing mortality (fraction of FMSY)
#' @return An object of class PMobj
#' @export
#'
#' @examples 
#' \dontrun{
#' POF(myMSE)
#' }
POF <- function(MSEobj=NULL, Ref=1) {
  if (class(MSEobj)!='MSE') stop('Require object of class MSE')
  PMobj <- new("PMobj")
  PMobj@Name <- "Fishing Mortality relative to FMSY"
  if (Ref !=1) {
    PMobj@Caption <- paste0('Prob. F < ', Ref, ' FMSY')
  } else {
    PMobj@Caption <- 'Prob. F < FMSY'
  }

  y.st <- 1 
  y.end <- MSEobj@proyears
  
  PMobj@Stat <- MSEobj@F_FMSY[,,y.st:y.end] # Performance Metric statistic of interest - here F/FMSY
  PMobj@Ref <- Ref
  PMobj@Prob <- calcProb(PMobj@Stat < PMobj@Ref, MSEobj) # calculate probability Stat < 1 nsim by nMP
  
  
  PMobj@Mean <- calcMean(PMobj@Prob, MSEobj) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
  
}
class(POF) <- "PM"

#' Performance Metric: Probability Long-Term Yield > 0.5 Relative Yield
#'
#' @param MSEobj An object of class MSE
#' @param Ref Reference Point for long-term yield (fraction of reference yield)
#' @return An object of class PMobj
#' @export
#'
#' @examples 
#' \dontrun{
#' LTY(myMSE)
#' }
LTY <- function(MSEobj=NULL, Ref=0.5) {
  if (class(MSEobj)!='MSE') stop('Require object of class MSE')
  PMobj <- new("PMobj")
  PMobj@Name <- "Average Long-Term Yield relative to Reference Yield" 
  if (Ref != 1) {
    PMobj@Caption <- paste0('Prob. LTY > ', Ref, ' Ref. Yield')
  } else {
    PMobj@Caption <- 'Prob. LTY > Ref. Yield'
  }
 
  
  y.st <- max(MSEobj@proyears-9,1) 
  y.end <- MSEobj@proyears # last 10 years
  
  RefYd <- array(MSEobj@OM$RefY, dim=dim(MSEobj@C[,,y.st:y.end]))
  
  PMobj@Stat <- MSEobj@C[,,y.st:y.end]/RefYd
  PMobj@Ref <- 0.5
  PMobj@Prob <- calcProb(PMobj@Stat > PMobj@Ref, MSEobj)  # probability LTY > 0.5 Ref Yield
  
  
  PMobj@Mean <- calcMean(PMobj@Prob, MSEobj) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
  
}
class(LTY) <- "PM"

#' Performance Metric: Probability Short-Term Yield > 0.5 Relative Yield
#'
#' @param MSEobj An object of class MSE
#' @param Ref Reference Point for short-term yield (fraction of reference yield)
#' @return An object of class PMobj
#' @export
#'
#' @examples 
#' \dontrun{
#' STY(myMSE)
#' }
STY <- function(MSEobj=NULL, Ref=0.5) {
  if (class(MSEobj)!='MSE') stop('Require object of class MSE')
  PMobj <- new("PMobj")
  PMobj@Name <- "Average Short-Term Yield relative to Reference Yield" 
  if (Ref != 1) {
    PMobj@Caption <- paste0('Prob. STY > ', Ref, ' Ref. Yield')
  } else {
    PMobj@Caption <- 'Prob. STY > Ref. Yield'
  }
  
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
#' @param Ref Reference Point for AAVY (fraction)
#' @return An object of class PMobj
#' @export
#'
#' @examples 
#' \dontrun{
#' AAVY(myMSE)
#' }
AAVY <- function(MSEobj=NULL, Ref=0.2) {
  if (class(MSEobj)!='MSE') stop('Require object of class MSE')
  PMobj <- new("PMobj")
  PMobj@Name <- "Average Annual Variability in Yield" 
  PMobj@Caption <- paste0('Prob. AAVY < ', Ref*100, "%")
  
  y1<-1:(MSEobj@proyears-1) # year index
  y2<-2:MSEobj@proyears # 
  
  if (MSEobj@nMPs > 1) {
    AAVY <- apply(((MSEobj@C[,,y1]-MSEobj@C[,,y2])^2)^0.5,c(1,2),mean)/apply(MSEobj@C[,,y2],c(1,2),mean) 
  } else {
    AAVY <- array(apply(((MSEobj@C[,1,y1]-MSEobj@C[,1,y2])^2)^0.5,c(1),mean)/apply(MSEobj@C[,1,y2],c(1),mean))
  }
  
  PMobj@Stat <- AAVY
  PMobj@Ref <- Ref
  PMobj@Prob <- calcProb(PMobj@Stat < Ref, MSEobj)  # probability AAVY < 0.2 
  
  PMobj@Mean <- calcMean(PMobj@Prob, MSEobj) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
  
}
class(AAVY) <- "PM"

#' Performance Metric: Average Yield
#'
#' @param MSEobj An object of class MSE
#' @param Ref Reference Point for yield (fraction of reference yield)
#' @return An object of class PMobj
#' @export
#'
#' @examples 
#' \dontrun{
#' Yield(myMSE)
#' }
Yield <- function(MSEobj=NULL, Ref=1) {
  if (class(MSEobj)!='MSE') stop('Require object of class MSE')
  PMobj <- new("PMobj")
  PMobj@Name <- "Yield relative to Reference Yield" 
  PMobj@Caption <- 'Mean Relative Yield'
  
  y.st <- 1 
  y.end <- MSEobj@proyears 
  
  RefYd <- array(MSEobj@OM$RefY, dim=dim(MSEobj@C[,,y.st:y.end]))
  
  PMobj@Stat <- MSEobj@C[,,y.st:y.end]/RefYd
  PMobj@Ref <- Ref
  PMobj@Prob <- calcProb(PMobj@Stat, MSEobj) # no probability to calculate
  
  PMobj@Mean <- calcMean(PMobj@Prob, MSEobj) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
  
}
class(Yield) <- "PM"


