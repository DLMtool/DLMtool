# ---- Plot Data Object -----
#' Plot Data object
#'
#' @param x object of class Data
#' @param upq Upper quantile of TACs for max ylim
#' @param lwq Lower quantile of TACs for min ylim
#' @param outline Logical. Include outliers in plot?
#' @param ...  Optional additional arguments passed to \code{boxplot}
#' @export
plot.Data <- function(x, upq=0.9, lwq=0.1, outline = FALSE, ...) {
  boxplot.Data(x, upq, lwq, outline = FALSE, ...)
}

#' Boxplot of TAC recommendations
#' 
#' @param x An object of class MSE
#' @param upq Upper quantile of TACs for max ylim
#' @param lwq Lower quantile of TACs for min ylim
#' @param ylim Optional numeric vector of length 2 to specify limits of y-axis.
#' @param outline Logical. Include outliers in plot?
#' @param col Optional colours to pass to \code{boxplot}
#' @param ...  Optional additional arguments passed to \code{boxplot}
#' @return Returns a data frame containing the information shown in the plot
#' @author A. Hordyk
#' @export
boxplot.Data <- function(x, upq=0.9, lwq=0.1, ylim=NULL, outline = FALSE, col = NULL, ...) {
  Data <- updateMSE(x)
  if (class(Data) != "Data")  stop("Object must be of class 'Data'")
  tacs <- t(Data@TAC[, , 1])
  
  if (all(is.na(tacs))) {
    message("Nothing found in TAC slot")
    return(invisible(NULL))
  }
  units <- TRUE
  if (length(nchar(x@Units)) < 1) units <- FALSE
  MPs <- Data@MPs
  ind <- grep("ref", MPs)
  if (length(ind) > 0) {
    tacs <- tacs[, -ind, drop=FALSE]
    MPs <- MPs[-ind]
  }
  
  # exclude NAs 
  nMPs <- dim(Data@TAC)[1]
  
  if (nMPs>1){
    allNAs <- colSums(apply(tacs, 2, is.na)) == nrow(tacs)
    tacs <- tacs[,!allNAs, drop=FALSE]
    MPs <- MPs[!allNAs]
    nMPs<-length(MPs)
  }
  
  if (nMPs>1) {
    if (is.null(col)) col <- rainbow(30)
    ord <- order(apply(tacs, 2, median, na.rm = TRUE))
    MPs <- MPs[ord]
    tacs <- tacs[, ord]
    ymax <- max(apply(tacs, 2, quantile, upq, na.rm = TRUE))
    ymin <- min(apply(tacs, 2, quantile, lwq, na.rm = TRUE))
    if (is.null(ylim)) ylim <- c(ymin, ymax)
    Median <- round(apply(tacs, 2, median, na.rm = TRUE), 2)
    SD <- round(apply(tacs, 2, sd, na.rm = TRUE), 2)
  } else {
    if (is.null(ylim)) ylim <- c(min(tacs), max(tacs))
    Median <- median(tacs)
    SD <- sd(tacs)
    tacs <- as.numeric(tacs)
    if (is.null(col)) col <- "darkgray"
  }
  
  par(mfrow = c(1, 1), oma = c(2, 4, 1, 0), mar = c(3, 3, 0, 0))
  if (nMPs>1) {
    boxplot(tacs, names = MPs, las = 1, col = col, outline = outline, 
          frame = FALSE, ylim = ylim, horizontal = TRUE, ...)
    if (units) mtext(paste("TAC (", Data@Units, ")", sep = ""), side = 1, outer = T, 
                     line = 0.5, cex = 1.25)
    if (!units) mtext("TAC (no units supplied)", side = 1, outer = T, 
                      line = 0.5, cex = 1.25)
    mtext(side = 2, "Management Procedures", outer = TRUE, line = 3, cex = 1.25)
  } else {
    boxplot(tacs, names = MPs, las = 1, col = col, outline = outline, 
            frame = FALSE, ylim = ylim, horizontal = FALSE, ...)
    if (units) mtext(paste("TAC (", Data@Units, ")", sep = ""), side = 2, outer = T, 
                     line = 0.5, cex = 1.25)
    if (!units) mtext("TAC (no units supplied)", side = 2, outer = T, 
                      line = 0.5, cex = 1.25)
    mtext(side = 3, MPs, outer = TRUE, line=-1, cex = 1.25, xpd=NA)
  }
 
  if (units) {
      data.frame(MP = MPs, Median = Median, SD = SD, Units = Data@Units)
  } else {
      data.frame(MP = MPs, Median = Median, SD = SD)
  }
}