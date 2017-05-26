#'  Long-Term Yield Peformance Metric Function
#'
#'  An object of class PM
#'
"LTY"

#'  Short-Term Yield Peformance Metric Function
#'
#'  An object of class PM
#'
"STY"

#'  Spawning Biomass Relative to SBMSY
#'
#'  An object of class PM
#'
"SB_BMSY"

#'  Spawning Biomass Relative to SB0
#'
#'  An object of class PM
#'
"SB_B0"

#'  Fishing Mortality Relative to FMSY
#'
#'  An object of class PM
#'
"F_FMSY"




setClassUnion("label.class", c("call", "character", "function"))
setClassUnion("yr.class", c("numeric", "logical"))
setClassUnion("stat.class", c("character", "function"))

#' Class \code{'PM'}
#' 
#' A Performance Metrics Function Object
#' 
#' @name PM-class
#' @docType class
#' @section Create Object: Objects of class PM can be created by calls of the form \code{new('PM')}
#' @slot Name The name of the Performance Metric object 
#' @details ADD DETAILS HERE 

#' @keywords classes
#' @examples
#' 
#' showClass('PM')
#' 
setClass("PM", representation(Name = "character",  Description="character",  
                              Func = "function",
                              Ref="numeric",  Y1="yr.class", Y2="yr.class", Stat = "stat.class", 
                              LRP="numeric", TRP="numeric",
                              Label="label.class", Var="character"))




#' Internal function to check or correct years in performance metrics
#'
#' @param MSEobj An object of class 'MSE'
#' @param y1 A numeric value indicating the first year to summarize performance
#' @param y2 A numeric value indicating the last year to summarize performance
#'
#' @return A data.frame containing y1 and y2
#'
ChkYrs <- function(MSEobj, y1=NULL, y2=NULL) {
  if (!is.numeric(y2) || is.null(y2) || is.na(y2)) y2 <- MSEobj@proyears
  if (y2 > MSEobj@proyears) {
    y2 <- MSEobj@proyears
    warning("y2 is greater than proyears. Defaulting to y2=proyears", call.=FALSE)
  }
  if (!is.numeric(y1) || is.null(y2) || is.na(y1)) {
    y1 <- 1 
  } else {
    if (y1 < 0) y1 <- y2 - abs(y1) + 1 
    if (y1 == 0) y1 <- 1
  }
  if (y1 <= 0) {
    y1 <- 1
    warning("y1 is negative. Defaulting to y1=1", call.=FALSE)
  }
  if (y1 >= y2) {
    y1 <- y2 - 9 
    warning("y1 is greater or equal to y2. Defaulting to y1=y2-9", call.=FALSE)
    
  }
  data.frame(y1=y1, y2=y2)
}





