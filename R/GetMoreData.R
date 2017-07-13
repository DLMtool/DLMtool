
#' Load more data from DLMdata package
#'
#' Downloads the DLMdata package from GitHub 
#' @param silent Logical. Should messages to printed?
#' @export
#'
#' @importFrom devtools install_github
GetMoreData <- function(silent=FALSE) {
  if (!silent) message("\nDownloading 'DLMdata' from GitHub")
  tt <- devtools::install_github("DLMtool/DLMdata", quiet=TRUE)
  if (tt) {
    if (!silent) message("Use 'library(DLMdata)' to load additional data into workspace")
  } else {
    if (!silent) message("Package 'DLMdata' already up to date\n Use 'library(DLMdata)' to load additional data into workspace")
  }
  # d <- data(package = "DLMdata")
  # DataObjs <- d$results[,3]
  # for (X in 1:length(DataObjs)) {
  #   dat <- eval(parse(text=paste0("DLMdata::",(DataObjs[X]))))
  #   AAname <- DataObjs[X]
  #   assign(AAname, dat)
  #   # environment(AAname) <- asNamespace('DLMtool')
  #   # environment(AAname) <- as.environment("package:DLMtool")
  # }

}