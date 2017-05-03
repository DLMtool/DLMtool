#' Replace an existing Stock, Fleet, Obs, or Imp object 
#' 
#' A function that replaces a Stock, Fleet, Obs, or Imp object from an 
#' OM with one from another OM. Mainly used for internal functions.
#' 
#' @param OM An operating model object (class OM) which will be updated with a sub-model from another OM
#' @param from The OM object from which the sub-model is being taken
#' @param Sub A character string specifying what object type to replace
#' "Stock", "Fleet", "Obs" or "Imp" (default is all four which is probably not what you want to do)
#' @return An object of class OM
#' @author A. Hordyk
#' @export 
Replace <- function(OM, from, Sub=c("Stock", "Fleet", "Obs", "Imp")) {
  if (class(OM) =="character") OM <- get(OM)
  if (class(from) !="OM") fromOM <- get(from)
  if (class(OM) !="OM") stop("OM must be of class OM ", call.=FALSE)
  if (class(from) !="OM") stop("''from' must be of class OM ", call.=FALSE)
  Sub <- match.arg(Sub, several.ok=TRUE)
  
  Stock <- SubOM(OM, "Stock")
  Fleet <- SubOM(OM, "Fleet")
  Obs <- SubOM(OM, "Obs")
  Imp <- SubOM(OM, "Imp")
  
  message("Replacing sub-models:", paste0(" ", Sub))
  for (x in 1:length(Sub)) {
    assign(Sub[x], SubOM(from, Sub[x]))
  }
  
  outOM <- new("OM", Stock, Fleet, Obs, Imp) 
  outOM@nsim <- OM@nsim 
  outOM@cpars <- OM@cpars 
  outOM@seed <- OM@seed 
  outOM 
} 
