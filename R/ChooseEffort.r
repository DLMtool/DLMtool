
#' Manually map the historical relative fishing effort trajectory.
#' 
#' Interactive plot which allows users to specify the relative trajectory and
#' variability in the historical fishing effort and populates Fleet object.
#' 
#' 
#' @usage ChooseEffort(FleetObj, Years=NULL)
#' @param FleetObj A fleet object.
#' @param Years An optional vector of years. Should be nyears long.
#' @author A. Hordyk
#' @export ChooseEffort
ChooseEffort <- function(FleetObj, Years = NULL) {
  nyears <- FleetObj@nyears
  runSketch <- SketchFun(nyears, Years)
  if (!is.null(Years)) FleetObj@nyears <- length(Years)
  FleetObj@EffYears <- runSketch[, 1]
  FleetObj@EffLower <- runSketch[, 2]
  FleetObj@EffUpper <- runSketch[, 3]
  return(FleetObj)
}

#' Manually map the historical relative fishing effort trajectory.
#' 
#' Internal function for innteractive plot which allows users to specify the relative trajectory and
#' variability in the historical fishing effort.
#' 
#' 
#' @usage SketchFun(nyears, Years=NULL)
#' @param nyears Number of years
#' @param Years An optional vector of years. Should be nyears long.
#' @author A. Hordyk
#' @export SketchFun
SketchFun <- function(nyears = NULL, Years = NULL) {
    
    if (length(Years) == 0 & length(nyears) == 0) 
        stop()
    if (length(Years) > 0) {
        nyears <- length(Years)
        years <- Years
    }
    if (length(Years) == 0) years <- 1:nyears
    par(mfrow = c(1, 1), mar = c(5, 4, 5, 2))
    ys <- seq(from = 0, to = 1, by = 0.05)
    years1 <- seq(from = years[1], by = 2, to = max(years))
    years2 <- seq(from = years[2], by = 2, to = max(years))
    grd1 <- expand.grid(years1, ys)
    grd2 <- expand.grid(years2, ys)
    grd3 <- expand.grid(years, ys)
    Xs <- grd3[, 1]
    Ys <- grd3[, 2]
    
    plot(grd1[, 1], grd1[, 2], col = "black", pch = 19, cex = 0.2, xlab = "Years", 
        ylab = "Relative Fishing Effort", cex.lab = 1.5)
    points(grd2[, 1], grd2[, 2], col = "darkgrey", pch = 19, cex = 0.2)
    
    mtext(side = 3, "Click 'Finish' or press 'Escape' to Finish.", xpd = NA, 
        cex = 1.25)
    line1 <- "Use mouse to select points on the grid"
    line2 <- "First and last year must be selected."
    line3 <- "Select two points in a single year to represent range of uncertainty"
    par(xpd = TRUE)
    text(years[1], par("usr")[4] + 0.15 * (par("usr")[4] - par("usr")[3]), 
        line1, cex = 1, pos = 4)
    text(years[1], par("usr")[4] + 0.1125 * (par("usr")[4] - par("usr")[3]), 
        line2, cex = 1, pos = 4)
    text(years[1], par("usr")[4] + 0.075 * (par("usr")[4] - par("usr")[3]), 
        line3, cex = 1, pos = 4)
    par(xpd = FALSE)
    message(line1, "\n", line2, "\n", line3, "\n")
    flush.console()
    
    
    out <- NULL
    out <- identifyPch(x = Xs, y = Ys, tolerance = 0.1)
    while (is.null(dim(out))) {
        message("Must choose more than one point")
        flush.console()
        out <- identifyPch(x = Xs, y = Ys, tolerance = 0.1)
    }
    while (min(out[, 1]) != years[1]) {
        message("Choose point(s) for first year (usually 0)")
        flush.console()
        dat <- rbind(out, identifyPch(x = Xs, y = Ys, tolerance = 0.1))
        out <- dat[order(dat[, 1]), ]
    }
    while (max(out[, 1]) != years[length(years)]) {
        message("Choose point(s) for last year (nyear)")
        flush.console()
        dat <- rbind(out, identifyPch(x = Xs, y = Ys, tolerance = 0.1))
        out <- dat[order(dat[, 1]), ]
    }
    ord <- order(out[, 1], out[, 2])
    out <- out[ord, ]
    yrs <- unique(out[, 1])
    mat <- matrix(NA, nrow = length(yrs), ncol = 3)
    mat[, 1] <- yrs
    ind <- which(!duplicated(out[, 1]))
    mat[, 2:3] <- out[ind, 2]
    for (X in seq_along(yrs)) {
        chk <- out[, 1] %in% yrs[X]
        ind <- range(which(chk))
        if (sum(chk) > 1) {
            mat[X, 2:3] <- out[ind, 2]
        }
    }
    
    lines(mat[, 1], mat[, 2])
    lines(mat[, 1], mat[, 3])
    
    colnames(mat) <- c("Years", "Lower", "Upper")
    return(mat)
}

identifyPch <- function(x, y = NULL, n = length(x), pch = 19, ...) {
  xy <- xy.coords(x, y)
  x <- xy$x
  y <- xy$y
  sel <- rep(FALSE, length(x))
  res <- integer(0)
  while (sum(sel) < n) {
    ans <- identify(x[!sel], y[!sel], n = 1, plot = FALSE, ...)
    if (!length(ans)) 
      break
    ans <- which(!sel)[ans]
    points(x[ans], y[ans], pch = pch)
    sel[ans] <- TRUE
    res <- c(res, ans)
  }
  out <- cbind(x[res], y[res])
  out <- out[order(out[, 1]), ]
  return(out)
}





