#' Open the DLMtool User Guide
#'
#' Opens the DLMtool User Guide website (requires internet connection)
#' 
#' @export
#' @importFrom utils browseURL
#' @examples
#' \dontrun{
#' userguide()
#' }
userguide <- function() {
  utils::browseURL("https://dlmtool.github.io/DLMtool/userguide/index.html")
}