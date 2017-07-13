
#' runMSE with no messages - for testing 
#'
#' For testing purposes only 
#' @param ... Arguments to runMSE function 
#'
#' @export
#' @keywords internal
#' @importFrom utils capture.output
#'
runMSEnomsg <- function(...) {
  capture.output(out <- suppressMessages(runMSE(...)))
  out
}