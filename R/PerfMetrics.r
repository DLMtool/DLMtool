#' Label class union for performance metric objects
#' 
#' @description Used internally. Nothing to see here!
#'  
#' @export
#' 
setClassUnion(name="label.class", members=c("call", "character", "function"))


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
                              Func = "function", Ref="numeric", Stat="character",
                              RP="numeric",  Label="label.class"))


setMethod("initialize", "PM", function(.Object) {
  .Object@Stat <- "mean"
  .Object
  
})



#' Calculate Trade-Off Performance Metrics
#' 
#' Function takes two performance metrics objects and returns a data.frame of calculated
#' performance metrics - used internally by the trade-off plot function
#'
#' @param PM1 An object of class 'PM'
#' @param PM2 An object of class 'PM'
#' @param MSEobj An object of class 'MSE' or a list of objects of class 'MSE'
#' @param MSEname An optional character string of names for the MSE objects if 'MSEobj' is 
#' a list
#' @param Can An optional vector with the output of the 'Can' function, or a logical vector
#' indicating which MPs can be run 
#' @param Class An optional vector indicating the class of MPs. 
#'
#' @return An object of class 'PMtrade' used in the trade-off plot function
#' @export
PMtrade <- function(PM1, PM2, MSEobj, MSEname=NULL, Can=NULL, Class=NULL) {
  if (class(PM1) != "PM" | class(PM2) != "PM") stop("First two arguments must be class 'PM'", call.=FALSE)
  lst <- list()
  if (is.null(MSEname)) MSEname <- 1:length(MSEobj)
  if (class(MSEname) == "character") MSEname <- factor(MSEname, ordered=TRUE, levels=MSEname)
  # Create a x,y data frame
  for (X in 1:length(MSEobj)) {
    if (class(MSEobj) == "MSE") MSE <- MSEobj
    if (class(MSEobj) == "list") MSE <- MSEobj[[X]]
    if (is.null(Can)) Can <- rep(NA, MSE@nMPs)
    if (is.null(Class)) Class <- MPclass(MSE@MPs)
    if (class(Can) == "character") {  # character string of MP names that can be used
      Can <- MSE@MPs %in% Can 
    }
    x <- apply(PM1@Func(MSE, PM1), 2, get(PM1@Stat))
    y <- apply(PM2@Func(MSE, PM2), 2, get(PM2@Stat)) 
    x.rp <- PM1@RP
    x.pass <- x >=x.rp
    if (length(x.rp) == 0) x.rp <- NA
    if (length(x.pass) == 0) x.pass <- NA
    y.rp <- PM2@RP
    y.pass <- y >= y.rp
    if (length(y.rp) == 0) y.rp <- NA
    if (length(y.pass) == 0) y.pass <- NA
    pass <- as.logical(x.pass * y.pass)
    lst[[X]] <- data.frame(MP=MSE@MPs, x=x, y=y, x.rp=x.rp, x.pass=x.pass, 
                           y.pass=y.pass, y.rp=y.rp, pass=pass,
                           MSE=MSEname[X], Can=Can, Class=Class, 
                           stringsAsFactors = FALSE)
  }
  DF <- do.call("rbind", lst)
  out <- list()
  out$DF <- DF 
  if (class(PM1@Label) == "function") out$xlab <- PM1@Label(MSE, PM1)
  if (class(PM1@Label) != "function") out$xlab <- PM1@Label
  if (class(PM2@Label) == "function") out$ylab <- PM2@Label(MSE, PM2)
  if (class(PM2@Label) != "function") out$ylab <- PM2@Label
  if (all(!grepl("Prob.", out$xlab))) {
    out$xlab <- bquote(.(simpleCap(PM1@Stat)) ~ .(out$xlab))
  } 
  if (all(!grepl("Prob.", out$ylab))) {
    out$ylab <- bquote(.(simpleCap(PM2@Stat)) ~ .(out$ylab))
  }
  class(out) <- "PMtrade"
  out
}

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

