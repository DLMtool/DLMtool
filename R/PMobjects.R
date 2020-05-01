# All performance metric functions share a single help file: PerformanceMetric.Rd

#' Performance Metrics Methods
#' 
#' Performance metric (PMs) methods for your management strategy evaluation.
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
#' \code{PNOF} \tab Probability F < FMSY \cr
#' \code{LTY} \tab Probability Long-Term Yield > 0.5 Relative Yield \cr
#' \code{STY} \tab Probability Short-Term Yield > 0.5 Relative Yield \cr
#' \code{AAVY} \tab Probability AAVY < 0.2 (Average Annual Variability in Yield) \cr
#' \code{AAVE} \tab Probability AAVE < 0.2 (Average Annual Variability in Effort) \cr
#' \code{Yield} \tab Average Yield (relative to Reference Yield) \cr
#' }
#' 
#' Argument `Ref` provides the ratio relative to the reference point for calculating
#' the performance metric. For biomass-based PMs (`P10`, `P50`, `P100`), this is the fraction of 
#' BMSY. For `PNOF`, the fraction of FMSY. For `Yield` (and `LTY`/`STY`), the fraction of the 
#' Reference Yield. For `AAVY` is it the maximum acceptable variability in yield (i.e, default 
#' for `AAVY` is `Ref=0.2`)
#' 
#' The `Yrs` argument defines the number of years to calculate the performance statistic over. 
#' A value of `NULL`, the default for `AAVY`, `AAVE`, `P10`, `P50`, `P100`, and `PNOF`, means that the 
#' performance metric is calculated over all projection years. A numeric vector of length two is used 
#' to specify the first and last year, e.g, if `Yrs=c(1,10)` the performance statistic is calculated 
#' over the first 10 projection years. A numeric vector of length one with positive or negative value 
#' respectively can be used to specify the first *x* or last *x* years, e.g, `Yrs=10` is first 10 years,
#' and `Yrs=-10` is the last 10 years. See \code{\link{ChkYrs}} for more details.
#' 
#' By default Long-Term Yield (`LTY`) is the Yield in the last ten years of the projection period in the MSE, 
#' and Short-Term Yield (`STY`) is that in the first 10 years of the projection period.
#'
#' @templateVar url performance-metrics
#' @templateVar ref NULL
#' @template userguide_link
#' 
#' @return An object of class `PMobj`
#' @examples 
#' \dontrun{
#' myMSE <- runMSE()
#' P10(myMSE)
#' P50(myMSE)
#' P100(myMSE)
#' PNOF(myMSE)
#' LTY(myMSE)
#' STY(myMSE)
#' AAVY(myMSE)
#' AAVE(myMSE)
#' Yield(myMSE)
#' }
NULL

#' Check the years to summarize performance
#'
#' @param Yrs Numeric vector of length 2 with year indices to summarize performance. 
#' If NULL, the performance is summarized over all projection years. 
#' `Yrs` can also be length one, in which case if it is positive it is the first `Yrs` and 
#' if negative the last `Yrs` of the projection period. 
#' @param MSEobj An object of class `MSE` 
#'
#' @return A numeric vector of length 2 with year indices to summarize performance
#' @examples 
#' \dontrun{
#' MSE <- runMSE()
#' ChkYrs(NULL, MSE) # returns c(1, MSE@proyears)
#' ChkYrs(c(2,5), MSE) # returns c(2,5)
#' ChkYrs(c(70,80), MSE) # returns c(MSE@proyears-10,MSE@proyears)
#' ChkYrs(5, MSE) # returns c(1,5)
#' ChkYrs(-5, MSE) # returns c(46,50)
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
    if (length(Yrs) == 1) {
      if (Yrs == 0) stop("Yrs must be postive or negative", call.=FALSE)
      if (Yrs < 0) {
        y.st <- MSEobj@proyears + Yrs[1] + 1 
        y.end <- MSEobj@proyears
      } else {
        y.st <- 1 
        y.end <- y.st + Yrs[1] - 1 
      }
    } else {
      if (length(Yrs)>2) stop("Yrs must be numeric vector of length 1 or 2", call.=FALSE)
      
      y.st <- Yrs[1]
      y.end <- Yrs[2]
      
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
  
  PMobj@Mean <- calcMean(PMobj@Prob) # calculate mean probability by MP
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
PNOF <- function(MSEobj=NULL, Ref=1, Yrs=NULL) {
  Yrs <- ChkYrs(Yrs, MSEobj)
  PMobj <- new("PMobj")
  PMobj@Name <- "Probability of not overfishing (F<FMSY)"
  if (Ref !=1) {
    PMobj@Caption <- paste0('Prob. F < ', Ref, ' FMSY (Years ', Yrs[1], ' - ', Yrs[2], ')')
  } else {
    PMobj@Caption <- paste0('Prob. F < FMSY (Years ', Yrs[1], ' - ', Yrs[2], ')')
  }

  PMobj@Stat <- MSEobj@F_FMSY[,,Yrs[1]:Yrs[2]] # Performance Metric statistic of interest - here F/FMSY
  PMobj@Ref <- Ref
  PMobj@Prob <- calcProb(PMobj@Stat < PMobj@Ref, MSEobj) # calculate probability Stat < 1 nsim by nMP
  
  
  PMobj@Mean <- calcMean(PMobj@Prob) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
  
}
class(PNOF) <- "PM"

#' @rdname PerformanceMetric 
#' @export
LTY <- function(MSEobj=NULL, Ref=0.5, Yrs=-10) {
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
  PMobj@Ref <- Ref
  PMobj@Prob <- calcProb(PMobj@Stat > PMobj@Ref, MSEobj)  
  
  PMobj@Mean <- calcMean(PMobj@Prob) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
  
}
class(LTY) <- "PM"

#' @rdname PerformanceMetric 
#' @export
STY <- LTY 
formals(STY)$Yrs <- 10
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
  
  PMobj@Mean <- calcMean(PMobj@Prob) # calculate mean probability by MP
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
    AAVY <- apply(((((MSEobj@C[, , y1] - MSEobj@C[, , y2])/MSEobj@C[, , y2])^2)^0.5), c(1, 2), mean)
    # AAVY <- apply(((MSEobj@C[,,y1]-MSEobj@C[,,y2])^2)^0.5,c(1,2),mean)/apply(MSEobj@C[,,y2],c(1,2),mean) 
  } else {
    AAVY <- array(apply(((((MSEobj@C[,1,y1]-MSEobj@C[,1,y2])/MSEobj@C[,1,y2])^2)^0.5),c(1),mean))
  }
  
  PMobj@Stat <- AAVY
  PMobj@Ref <- Ref
  PMobj@Prob <- calcProb(PMobj@Stat < Ref, MSEobj)  # probability AAVY < 0.2 
  
  PMobj@Mean <- calcMean(PMobj@Prob) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
  
}
class(AAVY) <- "PM"

#' @rdname PerformanceMetric 
#' @export
AAVE <- function(MSEobj=NULL, Ref=0.2, Yrs=NULL) {
  Yrs <- ChkYrs(Yrs, MSEobj)
  PMobj <- new("PMobj")
  PMobj@Name <- paste0("Average Annual Variability in Effort (Years ", Yrs[1], "-", Yrs[2], ")") 
  PMobj@Caption <- paste0('Prob. AAVE < ', Ref*100, "% (Years ", Yrs[1], "-", Yrs[2], ")")
  
  y1<- Yrs[1]:(Yrs[2]-1) # year index
  y2<-(Yrs[1]+1):Yrs[2] 
  
  if (MSEobj@nMPs > 1) {
    AAVE <- apply(((((MSEobj@Effort[, , y1] - MSEobj@Effort[, , y2])/MSEobj@Effort[, , y2])^2)^0.5), c(1, 2), mean)
  } else {
    AAVE <- array(apply(((((MSEobj@Effort[,1,y1]-MSEobj@Effort[,1,y2])/MSEobj@Effort[,1,y2])^2)^0.5),c(1),mean))
  }
  
  PMobj@Stat <- AAVE
  PMobj@Ref <- Ref
  PMobj@Prob <- calcProb(PMobj@Stat < Ref, MSEobj)  # probability AAVE < 0.2 
  
  PMobj@Mean <- calcMean(PMobj@Prob) # calculate mean probability by MP
  PMobj@MPs <- MSEobj@MPs
  PMobj
  
}
class(AAVE) <- "PM"




#' Calculate Probability
#' 
#' @param PM A PM method 
#'
#' @export
#' @keywords internal
#'
calcProb <- function(PM, MSEobj) {
  if (MSEobj@nMPs > 1) {
    mar <- 2 
  } else mar <- 1
  mar <- 1:mar
  apply(PM, mar, mean, na.rm=TRUE)
}


#' Calculate Mean Probability
#' 
#' @param Prob Prob slot from an object of class PMobj 
#'
#' @export
#' @keywords internal
#'
calcMean <- function(Prob) {
  if ('matrix' %in% class(Prob)) return(apply(Prob , 2, mean, na.rm=TRUE))
  if ('numeric' %in% class(Prob)) return(mean(Prob, na.rm=TRUE))
}

