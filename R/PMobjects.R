# All performance metric functions share a single help file: PerformanceMetric.Rd

#' Performance Metrics
#' 
#' Performance metrics (PMs) for your management strategy evaluation.
#' 
#' @name PerformanceMetric
#' @param MSEobj An object of class MSE
#' @param Ref Reference point for calculating the performance metric. See details.
#' @param Yrs Numeric vector of length 2 with year indices to summarize performance. 
#' If NULL, the performance is summarized over all projection years.
#'  
#' @details Performance Metric definitions:
#' 
#' \tabular{ll}{
#' \code{P10} \tab Probability B > 0.1 BMSY \cr
#' \code{P50} \tab Probability B > 0.5 BMSY \cr
#' \code{P100} \tab Probability B > BMSY \cr
#' \code{POF} \tab Probability F < FMSY \cr
#' \code{LTY} \tab Probability Long-Term Yield > 0.5 Relative Yield \cr
#' \code{STY} \tab Probability Short-Term Yield > 0.5 Relative Yield \cr
#' \code{AAVY} \tab Probability AAVY < 0.2 (Average Annual Variability in Yield) \cr
#' \code{Yield} \tab Average Yield (relative to Reference Yield) \cr
#' }
#' 
#' Argument \code{Ref} provides the ratio relative to the reference point for calculating
#' the performance metric. For biomass-based PMs (P10, P50, P100), this is the fraction of 
#' BMSY. For POF, the fraction of FMSY. For Yield (and LTY/STY), the fraction of the 
#' Reference Yield.
#' 
#' Long-Term Yield is the Yield in the last ten years of the projection period in the MSE.
#' Short-Term Yield is that in the first 10 years of the projection period
#' @return An object of class PMobj
#' @examples 
#' \dontrun{
#' P10(myMSE)
#' P50(myMSE)
#' P100(myMSE)
#' POF(myMSE)
#' LTY(myMSE)
#' STY(myMSE)
#' AAVY(myMSE)
#' Yield(myMSE)
#' }
NULL

#' Check the years to summarize performance
#'
#' @param Yrs Numeric vector of length 2 with year indices to summarize performance. 
#' If NULL, the performance is summarized over all projection years.
#' @param MSEobj An object of class `MSE` 
#'
#' @return A numeric vector of length 2 with year indices to summarize performance
#' @examples 
#' \dontrun{
#' MSE <- runMSE()
#' ChkYrs(NULL, MSE) # returns c(1, MSE@proyears)
#' ChkYrs(c(2,5), MSE) # returns c(2,5)
#' ChkYrs(c(70,80), MSE) # returns c(MSE@proyears-10,MSE@proyears)
#' }
#' 
#' @keywords internal
#' @export
ChkYrs <- function(Yrs, MSEobj) {
  if (class(MSEobj) !='MSE') stop('Require object of class MSE', call.=FALSE)
  if (is.null(Yrs)) {
    y.st <- 1 
    y.end <- MSEobj@proyears
  } else {
    y.st <- Yrs[1]
    y.end <- Yrs[2]
    if (length(Yrs)!=2) stop("Yrs must be numeric vector of length 2", call.=FALSE)
    if (Yrs[1] > Yrs[2]) stop("Yrs[1] is > Yrs[2]", call.=FALSE)
    if (any(Yrs < 1)) stop("Yrs must be positive", call.=FALSE)
    
    if (Yrs[2] > MSEobj@proyears) {
      message('Yrs[2] is greater than MSEobj@proyears. Setting Yrs[2] = MSEobj@proyears')
      y.end <- MSEobj@proyears
    }  
    if (Yrs[1] > MSEobj@proyears) {
      message('Yrs[1] is greater than MSEobj@proyears. Setting Yrs[1] = Yrs[2] - Yrs[1]')
      y.st <- max(1,y.end - (Yrs[2] - Yrs[1]))
    }
  }
  return(c(y.st, y.end))
}



#' @rdname PerformanceMetric 
#' @export
P10 <- function(MSEobj=NULL, Ref=0.1, Yrs=NULL) {
  Yrs <- ChkYrs(Yrs, MSEobj)
  
  PMobj <- new("PMobj")
  PMobj@Name <- "Spawning Biomass relative to SBMSY"
 
  if (Ref !=1) {
    PMobj@Caption <- paste0('Prob. SB > ', Ref, ' SBMSY (Years ', Yrs[1], ' - ', Yrs[2], ')')
  } else {
    PMobj@Caption <- paste0('Prob. SB > SBMSY (Years ', Yrs[1], ' - ', Yrs[2], ')')
  }
  
  PMobj@Ref <- Ref
  PMobj@Stat <- MSEobj@B_BMSY[,,Yrs[1]:Yrs[2]] # Performance Metric statistic of interest - here SB/SBMSY 
  PMobj@Prob <- calcProb(PMobj@Stat > PMobj@Ref, MSEobj) # calculate probability Stat > 0.1 nsim by nMP
  
  PMobj@Mean <- calcMean(PMobj@Prob, MSEobj) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
  
}
class(P10) <- "PM"


#' @rdname PerformanceMetric 
#' @export
P50 <- P10
formals(P50)$Ref <- 0.5
class(P50) <- "PM"

#' @rdname PerformanceMetric 
#' @export
P100 <- P10 
formals(P100)$Ref <- 1
class(P100) <- "PM"


#' @rdname PerformanceMetric 
#' @export
POF <- function(MSEobj=NULL, Ref=1, Yrs=NULL) {
  Yrs <- ChkYrs(Yrs, MSEobj)
  PMobj <- new("PMobj")
  PMobj@Name <- "Fishing Mortality relative to FMSY"
  if (Ref !=1) {
    PMobj@Caption <- paste0('Prob. F < ', Ref, ' FMSY (Years ', Yrs[1], ' - ', Yrs[2], ')')
  } else {
    PMobj@Caption <- paste0('Prob. F < FMSY (Years ', Yrs[1], ' - ', Yrs[2], ')')
  }

  PMobj@Stat <- MSEobj@F_FMSY[,,Yrs[1]:Yrs[2]] # Performance Metric statistic of interest - here F/FMSY
  PMobj@Ref <- Ref
  PMobj@Prob <- calcProb(PMobj@Stat < PMobj@Ref, MSEobj) # calculate probability Stat < 1 nsim by nMP
  
  
  PMobj@Mean <- calcMean(PMobj@Prob, MSEobj) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
  
}
class(POF) <- "PM"

#' @rdname PerformanceMetric 
#' @export
LTY <- function(MSEobj=NULL, Ref=0.5, Yrs=c(max(MSEobj@proyears-9,1), MSEobj@proyears)) {
  Yrs <- ChkYrs(Yrs, MSEobj)
  PMobj <- new("PMobj")
  PMobj@Name <- paste0("Average Yield relative to Reference Yield (Years ", Yrs[1], "-", Yrs[2], ")") 
  if (Ref != 1) {
    PMobj@Caption <- paste0('Prob. Yield > ', Ref, ' Ref. Yield (Years ', Yrs[1], "-", Yrs[2], ")") 
  } else {
    PMobj@Caption <- paste0('Prob. Yield > Ref. Yield (Years ', Yrs[1], "-", Yrs[2], ")") 
  }

  RefYd <- array(MSEobj@OM$RefY, dim=dim(MSEobj@C[,,Yrs[1]:Yrs[2]]))
  
  PMobj@Stat <- MSEobj@C[,,Yrs[1]:Yrs[2]]/RefYd
  PMobj@Ref <- 0.5
  PMobj@Prob <- calcProb(PMobj@Stat > PMobj@Ref, MSEobj)  
  
  PMobj@Mean <- calcMean(PMobj@Prob, MSEobj) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
  
}
class(LTY) <- "PM"

#' @rdname PerformanceMetric 
#' @export
STY <- LTY 
formals(STY)$Yrs <- c(1,10)
class(STY) <- "PM"



#' @rdname PerformanceMetric 
#' @export
Yield <- function(MSEobj=NULL, Ref=1, Yrs=NULL) {
  Yrs <- ChkYrs(Yrs, MSEobj)
  PMobj <- new("PMobj")
  PMobj@Name <- paste0("Yield relative to Reference Yield (Years ", Yrs[1], "-", Yrs[2], ")") 
  PMobj@Caption <- paste0("Mean Relative Yield (Years ", Yrs[1], "-", Yrs[2], ")")
  
  RefYd <- array(MSEobj@OM$RefY, dim=dim(MSEobj@C[,,Yrs[1]:Yrs[2]]))
  
  PMobj@Stat <- MSEobj@C[,,Yrs[1]:Yrs[2]]/RefYd
  PMobj@Ref <- Ref
  PMobj@Prob <- calcProb(PMobj@Stat, MSEobj) # no probability to calculate
  
  PMobj@Mean <- calcMean(PMobj@Prob, MSEobj) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
  
}
class(Yield) <- "PM"

#' @rdname PerformanceMetric 
#' @export
AAVY <- function(MSEobj=NULL, Ref=0.2, Yrs=NULL) {
  Yrs <- ChkYrs(Yrs, MSEobj)
  PMobj <- new("PMobj")
  PMobj@Name <- paste0("Average Annual Variability in Yield (Years ", Yrs[1], "-", Yrs[2], ")") 
  PMobj@Caption <- paste0('Prob. AAVY < ', Ref*100, "% (Years ", Yrs[1], "-", Yrs[2], ")")
  
  y1<- Yrs[1]:(Yrs[2]-1) # year index
  y2<-(Yrs[1]+1):Yrs[2] 
  
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



