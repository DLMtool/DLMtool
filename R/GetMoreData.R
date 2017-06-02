
#' Load more data from DLMdata package
#'
#' Downloads the DLMdata package from GitHub 
#' @export
#'
#' @importFrom devtools install_github
GetMoreData <- function() {
  message("Downloading 'DLMdata' from GitHub")
  tt <- devtools::install_github("DLMtool/DLMdata", quiet=TRUE)
  if (tt) {
    message("Use 'library(DLMdata)' to load additional data into workspace")  
  } else {
    message("Package 'DLMdata' already up to date\n Use 'library(DLMdata)' to load additional data into workspace")
  }
}