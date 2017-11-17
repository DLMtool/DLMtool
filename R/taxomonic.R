
# check for internet 
# https://stackoverflow.com/a/5078002/2885462
havingIP <- function() {
  if (.Platform$OS.type == "windows") {
    ipmessage <- system("ipconfig", intern = TRUE)
  } else {
    ipmessage <- system("ifconfig", intern = TRUE)
  }
  validIP <- "((25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)[.]){3}(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)"
  any(grep(validIP, ipmessage))
}


#' Retrieve taxonomic information from online database
#'
#' @param Species 
#'
#' @return A data.frame with the taxonomic information 
#' @author A. Hordyk
#' @export
getTaxo <- function(Species) {
  if (havingIP()) { # check for internet connection
    return(taxize::classification(Species, "worms", rows=1))
  } else message("No internet connection")
  
}

