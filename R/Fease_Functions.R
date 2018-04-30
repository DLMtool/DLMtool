#' MP feasibility diagnostic
#' 
#' What MPs may be run (best case scenario) for various data-availability
#' scenarios and management constraints?
#' 
#' 
#' @param Data An object of class 'Data'. Optional. If Data object is included, the returned MPs are both feasible (in terms of management)
#' and possible (sufficient data to run MP)
#' @param TAC Logical. Are catch limits feasible for this fishery?
#' @param TAE Logical. Are effort controls feasible for this fishery?
#' @param SL Logical. Are size-selectivity regulations (either gear changes or size-retention regulations) feasible for this fishery?
#' @param Spatial Logical. Are spatial closures feasible for this fishery?
#' @param names.only Logical. Should only the names of the feasible MPs be returned (default)? If FALSE, a data frame with MP name, and two columns
#' of logical values: Can (possible given data) and Fease (feasible given management constraints) is returned 
#' @param msg Logical. Should messages be printed to the console?
#' @param include.ref Logical. Should reference MPs (e.g. FMSYref) be included as feasible methods? Default is FALSE
#' 
#' @return Either a vector of MP names that are feasible for the fishery (default) or a 3 column data frame (`names.only=FALSE`). 

#' @author T. Carruthers & A. Hordyk
#' @examples 
#' Fease(TAC=FALSE)
#' Fease(SL=FALSE, Spatial=FALSE)
#' Data <- get(avail("Data")[1])
#' Fease(Data, TAE=FALSE, names.only=FALSE)
#' @export 
Fease <- function(Data=NULL, TAC=TRUE, TAE=TRUE, SL=TRUE, Spatial=TRUE, names.only=TRUE, msg=TRUE, include.ref=FALSE) {
  if (msg) {
    message("Feasible management: ")
    if (TAC) message("TAC - total allowable catch")
    if (TAE) message("TAE - total allowable effort")
    if (SL) message("SL - size selectivity")
    if (Spatial) message("Spatial - spatial closures")
  }
  if (!(TAC | TAE | SL | Spatial)) stop("No feasible management options!", call.=FALSE)
  MPs <- avail('MP')
  if (class(Data) == "Data") {
    if (msg) message("Data object provided. Returning feasible and available MPs")
    canMPs <- Can(Data)  
  } else {
    if (msg) message("No Data object provided. Returning feasible MPs")
    canMPs <- MPs
  }
  mptypes <- MPtype(MPs)
  mprecs <- mptypes[,3]
  isfease <- rep(TRUE, length(MPs))
  isfease[17]
  cbind(MPs, mprecs)
  
  if (!TAC) isfease[grepl("TAC", mprecs)] <- FALSE
  if (!TAE) isfease[grepl("TAE", mprecs)] <- FALSE
  if (!SL) isfease[grepl("SL", mprecs)] <- FALSE
  if (!Spatial) isfease[grepl("Spatial", mprecs)] <- FALSE
  
  
  df <- data.frame(MP=mptypes[,1], Can=mptypes[,1]%in%canMPs, Fease=isfease, stringsAsFactors = FALSE)
  df <- df[order(df$MP),]
  if (!include.ref)df <- df[mptypes[,2] != "Reference",]
  
  if (names.only) {
    return(df$MP[df$Can & df$Fease])
  } else {
    return(df)
  }
}

