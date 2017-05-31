
#' Load more data from DLMdata package
#'
#' Downloads the DLMdata package from GitHub 
#' @export
#'
#' @importFrom devtools install_github
GetMoreData <- function() {
  message("Downloading 'DLMdata' from GitHub")
  devtools::install_github("DLMtool/DLMdata")
  message("Use 'library(DLMdata)' to load additional data into workspace")
}