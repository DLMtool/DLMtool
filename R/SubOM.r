#' Subset a Stock, Fleet, or Obs object from an OM object
#' 
#' A function that strips out a Stock, Fleet, or Obs object from a 
#' complete OM object. Mainly used for internal functions.
#' 
#' @param OM An operating model object (class OM)
#' @param Sub A character string specifying what object type to strip out
#' "Stock", "Fleet", or "Obs"
#' @return An object of class Stock, Fleet, or Obs
#' @author A. Hordyk
#' @export SubOM
SubOM <- function(OM, Sub=c("Stock", "Fleet", "Obs")) {
  if (class(OM) !="OM") stop("OM must be of class OM ", call.=FALSE)
  Sub <- match.arg(Sub)
  temp <- new(Sub)
  
  slots <- slotNames(temp)
  for (X in seq_along(slots)) 
    slot(temp, slots[X]) <- slot(OM, slots[X]) 

  colon <- gregexpr(":", temp@Name)
  space <- gregexpr("  ", temp@Name)
  ind <- switch(Sub, Stock=1, Fleet=2, Obs=3)
  
  if (ind < 3) temp@Name <- substr(temp@Name, colon[[1]][ind]+1, space[[1]][ind]-1)
  if (ind == 3) temp@Name <- substr(temp@Name, colon[[1]][ind]+1, nchar(temp@Name))
 
  temp 
}
