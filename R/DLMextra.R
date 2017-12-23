
#' Load more data from DLMextra package
#'
#' Downloads the DLMextra package from GitHub 
#' @param silent Logical. Should messages to printed?
#' @export
#'
#' @importFrom devtools install_github
DLMextra <- function(silent=FALSE) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if (!silent) message("\nDownloading 'DLMextra' from GitHub")
  tt <- devtools::install_github("DLMtool/DLMextra", quiet=TRUE)
  if (tt) {
    if (!silent) message("Use 'library(DLMextra)' to load additional data into workspace")
  } else {
    if (!silent) message("Package 'DLMextra' already up to date\n Use 'library(DLMextra)' to load additional data into workspace")
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