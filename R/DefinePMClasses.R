#' Label class union for performance metric objects
#' 
#' @description Used internally. Nothing to see here!
#'  
#' @export
setClassUnion(name="label.class", members=c("call", "character", "function"))

#' Prob class union for performance metric objects
#' 
#' @description Used internally. Nothing to see here!
#'  
#' @export
setClassUnion(name="prob.class", members=c("matrix", "numeric", "data.frame"))



#' An object for storing data for analysis using data-limited methods
#' 
#' Used interally
#' 
#' @name PMobj-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('PMobj')} 
#' @slot name Name of the Performance Metric. Character 
#' @slot caption A caption to be used in plots. Character, call, or function.
#' @slot Stat Statistic of interest for the PM. Dimensions: nsim, nMP, yrs. Array 
#' @slot Prob Probability (mean over years) Dimensions: nsim by MP.  Matrix, numeric or data.frame  
#' @slot Mean Mean probability (mean over years and simulations). Numeric. Length nMPs 
#' @slot MPs Name of MPs. Single value. Character string  
#' @author  A. Hordyk
#' @keywords classes

setClass("PMobj", representation(name = "character",  caption='label.class', 
                                 Stat='array', Prob='prob.class', Mean='numeric',
                                 MPs="character"))


calcProb <- function(PM,  MSEobj) {
  mar <- ifelse(MSEobj@nMPs>1, 2, 1) # set margins for apply
  mar <- 1:mar
  apply(PM, mar, mean)
}

calcMean <- function(Prob, MSEobj) {
  if (class(Prob) == 'matrix') return( apply(Prob , 2, mean))
  if (class(Prob) == 'numeric') return(mean(Prob))
}


#' Show the output of a PM
#'
#' @param object object of class MSE
#' @rdname show-MSE
#' @export
setMethod("show", signature = (object="PMobj"), function(object) {
  cat(object@name)
  cat("\n", object@caption)
  cat("\n")
  
  nMP <- length(object@MPs)
  if (nMP > 1) nsim <- dim(object@Prob)[1]
  if (nMP == 1) nsim <- length(object@Prob)
  
  nprint <- min(nsim, 10)
  if (nMP > 1) df <- data.frame(object@Prob[1:nprint,])
  if (nMP == 1) df <- data.frame(object@Prob[1:nprint])
  if (nMP > 1) lst <- object@Prob[nprint+1,]
  if (nMP == 1) lst <- object@Prob[nprint+1]
  df <- signif(df,2)
  lst <- signif(lst,2)
  colnames(df) <- object@MPs
  names(lst) <- object@MPs
  if (nsim > (nprint+1)) {
    df <- rbind(df,
                rep(".", nMP),
                rep(".", nMP),
                rep(".", nMP),
                lst)
    rownames(df) <- c(1:(nprint+3), nsim)
  }
  print(df)
  
  cat("\nMean\n")
  print(signif(object@Mean,2))
})


#' Summary of MSE object
#'
#' @param object object of class MSE
#' @param ... a list of PM methods
#' @rdname summary-MSE
#' @export
setMethod('summary', signature="MSE", function(object, ...) {
  cl <- try(class(...), silent=TRUE)
  if (class(cl) == "try-error") PMlist <- list(...) 
  if (cl == "list") PMlist <- unlist(list(...))
  if (cl == "PM") PMlist <- list(...)
  # check
  if (any(unlist(lapply(PMlist, class)) != "PM")) {
    stop("Arguments must be PM methods")
  }
  
  if(length(PMlist) == 0) PMlist <- list(DLMtool::Yield, DLMtool::POF, DLMtool::P10, 
                                         DLMtool::P50, DLMtool::P100, DLMtool::LTY, 
                                         DLMtool::STY, DLMtool::AAVY)
  
  message("Calculating Performance Metrics")
  storeMean <- vector('list', length(PMlist))
  storeName <- vector('list', length(PMlist))
  storeCaption <- vector('list', length(PMlist))
  storeMP <- vector('list', length(PMlist))
  for (X in 1:length(PMlist)) {
    flush.console()
    runPM <- PMlist[[X]](object)
    storeMean[[X]] <- runPM@Mean
    storeName[[X]] <- runPM@name
    storeCaption[[X]] <- runPM@caption
    storeMP[[X]] <- runPM@MPs
  }
  
  df <- data.frame('MP'=storeMP[[1]], signif(do.call('cbind', storeMean),2))
  caps <- do.call('rbind', storeCaption)
  colnames(df)[2:(length(PMlist)+1)] <- 1:length(PMlist) #caps # gsub(" ", "", caps)
  print(data.frame('Performance Metrics' = do.call('rbind', storeName)))
  cat("\n")
  print(df)
  invisible(df)
  
})

