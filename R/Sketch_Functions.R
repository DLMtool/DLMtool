##### ChooseEffort, ChooseM, and ChooseSelect use one help file: Choose.Rd

#' Manually map parameters for the historical period of operating model
#' 
#' Interactive plots to specify trends and variability in 
#' fishing effort, fleet selectivity, and natural mortality for the
#' operating model.
#' 
#' @name Choose 
#' @param Fleet A fleet object.
#' @param Years An optional vector of years. Should be nyears long.
#' @param OM An object of class 'OM'
#' @param type A character string - is M to be mapped by 'age' or 'length'?
#' @param x Optional vector for x-axis
#' @param y Optional vector for y-axis
#' @param Stock Optional Stock object. If provided, average length-at-maturity
#' is included on plot for reference.
#' @param FstYr Optional value for first historical year. If empty, user must
#' specify the year in console.
#' @param SelYears Optional vector of values for each year where selectivity
#' pattern changed. If empty, user must specify the years in console (comma
#' separated).
#' @details 
#' \tabular{ll}{
#' \code{ChooseEffort} \tab Interactive plot which allows users to specify the 
#' relative trajectory and variability in the historical fishing effort and 
#' populates Fleet object. \cr
#' \code{ChooseM} \tab Interactive plot which allows users to specify M by age 
#' or size class \cr
#' \code{ChooseSelect} \tab Input the first historical year, and all years where 
#' selectivity pattern
#' changed (separated by comma).  Interactive plot which allows users to
#' specify a range for the length at 5\% and full selection (LFS), as well as
#' selectivity at maximum length for each year.  Produces a simple plot which
#' shows the range in selectivity pattern for each break-point year.
#' Selectivity-at-length is fixed in between break-point years.  Note that this
#' function replaces 'nyears' in the Fleet object with the value defined here
#' (FstYr:current year). \cr
#' }
#' @return \code{ChooseEffort} and \code{ChooseSelect} return a Fleet object while
#' \code{ChooseM} returns an OM object.
#' @author A. Hordyk
NULL


#' @rdname Choose
#' @export ChooseEffort
ChooseEffort <- function(Fleet, Years = NULL) {
  nyears <- Fleet@nyears
  runSketch <- SketchFun(nyears, Years)
  if (!is.null(Years)) Fleet@nyears <- length(Years)
  Fleet@EffYears <- runSketch[, 1]
  Fleet@EffLower <- runSketch[, 2]
  Fleet@EffUpper <- runSketch[, 3]
  return(Fleet)
}

#' Manually map the historical relative fishing effort trajectory.
#' 
#' Internal function for interactive plot which allows users to specify the relative trajectory and
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
  try(dev.off(), silent=TRUE)

  op <- par(mfrow = c(1, 1), mar = c(5, 4, 5, 2),no.readonly = TRUE)
  on.exit(par(op))
  
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
  
  out <- NULL
  out <- identifyPch(x = Xs, y = Ys, tolerance = 0.1)
  while (is.null(dim(out))) {
    message("Must choose more than one point")
    out <- identifyPch(x = Xs, y = Ys, tolerance = 0.1)
  }
  while (min(out[, 1]) != years[1]) {
    message("Choose point(s) for first year (usually 0)")
    dat <- rbind(out, identifyPch(x = Xs, y = Ys, tolerance = 0.1))
    out <- dat[order(dat[, 1]), ]
  }
  while (max(out[, 1]) != years[length(years)]) {
    message("Choose point(s) for last year (nyear)")
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
    if (!length(ans)) {
      break
    } 
    ans <- which(!sel)[ans]
    points(x[ans], y[ans], pch = pch)
    sel[ans] <- TRUE
    res <- c(res, ans)
  }
  out <- cbind(x[res], y[res])
  out <- out[order(out[, 1]), ]
  return(out)
}



#' @rdname Choose
#' @export ChooseM
ChooseM <- function(OM, type=c("age", "length"), x=NULL, y=NULL) {
  
  type <- match.arg(type)
  if (class(OM) != "OM") stop("First argument must be object of class 'OM'", call.=FALSE)
  
  ylab <- "Natural mortality (M)"
  if (is.null(y)) y <- seq(0, 1, 0.05)
  if (type == "age") {
    if (is.null(x)) x <- 1:OM@maxage
    xlab <- "Age Class"
    runSketch <- SketchM(x, y, xlab, ylab)  
    M1 <- approx(x = runSketch[,1], y = runSketch[,2], xout=x, method = "linear")$y  # linear interpolation
    M2 <- approx(x = runSketch[,1], y = runSketch[,3], xout=x, method = "linear")$y  # linear interpolation
    OM@M <- M1 
    OM@M2 <- M2 
    return(OM)
  }
  if (type == "length") {
    if (is.null(x)) x <- seq(from=0, to=max(OM@Linf), 5)
    xlab <- "Length Class"
    runSketch <- SketchM(x, y, xlab, ylab)  
    M1 <- approx(x = runSketch[,1], y = runSketch[,2], xout=x, method = "linear")$y  # linear interpolation
    M2 <- approx(x = runSketch[,1], y = runSketch[,3], xout=x, method = "linear")$y  # linear interpolation
    M_at_Length <- data.frame(Lens=x, M1, M2)
    OM@cpars[["M_at_Length"]] <- M_at_Length
    return(OM)
  }
}


SketchM <- function(x, y, xlab, ylab) {
  x1 <- x[seq(1, by=2, to=length(x))]
  x2 <- x[seq(2, by=2, to=length(x))]
  grd1 <- expand.grid(x1, y)
  grd2 <- expand.grid(x2, y)
  grd3 <- expand.grid(x, y)
  Xs <- grd3[, 1]
  Ys <- grd3[, 2]
  
  op <- par(mfrow=c(1,1), las=1, oma=c(1,2,1,2))
  plot(range(x), range(y), type="n", xlab = xlab, ylab = ylab, cex.lab = 1.5)
  points(grd1[, 1], grd1[, 2], col = "black", pch = 19, cex = 0.2)
  points(grd2[, 1], grd2[, 2], col = "darkgrey", pch = 19, cex = 0.2)
  axis(side=3)
  axis(side=4)  
  mtext(side = 3, "Click 'Finish' or press 'Escape' to Finish.", xpd = NA, 
        cex = 1.25, line=2.5)
  line1 <- "Use mouse to select points on the grid"
  line2 <- "First and last age/size classes must be selected."
  line3 <- "Select two points for each age/size class to represent range of uncertainty"
  message(line1, "\n", line2, "\n", line3, "\n")
  
  
  out <- NULL
  out <- identifyPch(x = Xs, y = Ys, tolerance = 0.1)
  while (is.null(dim(out))) {
    message("Must choose more than one point")
    out <- identifyPch(x = Xs, y = Ys, tolerance = 0.1)
  }
  while (min(out[, 1]) != x[1]) {
    message("Choose point(s) for first age/size class")
    dat <- rbind(out, identifyPch(x = Xs, y = Ys, tolerance = 0.1))
    out <- dat[order(dat[, 1]), ]
  }
  while (max(out[, 1]) != max(x)) {
    message("Choose point(s) for last age/size class")
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
  
  par(op)
  return(mat)
}


#' @rdname Choose
#' @export ChooseSelect
ChooseSelect <- function(Fleet, Stock, FstYr = NULL, SelYears = NULL) {
  # chk <- class(Fleet@isRel)
  # if (chk == "character") {
  # chkrel <- tolower(Fleet@isRel)
  # if (chkrel == "true" | Fleet@isRel == "1") isRel <- TRUE
  # if (chkrel == "false" | Fleet@isRel == "0") isRel <- FALSE
  # }
  # if (chk == "numeric") {
  # if (Fleet@isRel == 1) isRel <- TRUE
  # if (Fleet@isRel == 0) isRel <- FALSE
  # }
  if (class(Fleet) != "Fleet") stop("First argument must be object of class 'Fleet'", call.=FALSE)
  if (class(Stock) != "Stock") stop("First argument must be object of class 'Stock'", call.=FALSE)
  # if ((!isRel) & is.null(Stock))  stop("Require Stock object")
  Fleet@isRel <- 'FALSE'
  isRel <- FALSE
  if (is.null(FstYr)) {
    message("*****************")
    message("Enter first historical year")
    FstYr <- as.numeric(readline("First Historical Year: "))
    message("\n")
  }
  
  message("Enter last historical year")
  Fleet@CurrentYr <- as.numeric(readline("Last Historical Year: "))
  message("\n")
  
  LastYr <- Fleet@CurrentYr 
  if (LastYr < FstYr) stop("Last year before first year", call.=FALSE)
  
  if (is.null(SelYears)) {
    message("Enter each selectivity break point year, separated by a comma")
    message("Note: break points are the years where selectivity pattern changed")
    message("Note: break points must be within year range")
    inString <- readline("Enter each selectivity break point year, separated by a comma: ")
    options(warn = -1)
    SelYears <- as.numeric(unlist(strsplit(inString, ",")))
    if (is.na(SelYears)) SelYears <- as.numeric(unlist(strsplit(inString, " ")))
    options(warn = 0)
  }
  
  SelYears <- sort(SelYears)
  if (SelYears[1] != FstYr)  SelYears <- c(FstYr, SelYears)
  message("Break Points Years are: ")
  print(SelYears)
  if (length(SelYears) < 2) stop("Must be more than one year")
  if (max(SelYears) > LastYr) stop("Must specify historical year")
  if (min(SelYears) < FstYr) stop("Year before first year")
  flush.console()
  Selnyears <- length(SelYears)
  
  Years <- FstYr:LastYr  #SelYears[1]:SelYears[length(SelYears)]
  Fleet@nyears <- length(Years)
  ind <- round((Range(SelYears, Max = LastYr, Min = FstYr)) * Fleet@nyears, 0) + 1
  ind[length(ind)] <- max(ind) - 1
  Fleet@AbsSelYears <- SelYears
  Fleet@SelYears <- ind
  
  tempL5 <- matrix(0, nrow = Selnyears, ncol = 2)
  tempLFS <- matrix(0, nrow = Selnyears, ncol = 2)
  tempmaxlen <- matrix(0, nrow = Selnyears, ncol = 2)
  
  # if(is.null(Stock)) Stock <- NA
  set.par <- par(no.readonly = TRUE)
  message("Select selectivity points on plot")
  flush.console()
  for (N in 1:Selnyears) {
    BlankSelPlot(Stock = Stock, Yr = SelYears[N], N = N, isRel = isRel)
    L5Out <- ChooseL5(Fleet, Stock, isRel)
    tempL5[N, ] <- sort(L5Out[, 1])
    LFSout <- ChooseLFS(L5Out, Fleet, Stock, isRel)
    tempLFS[N, ] <- sort(LFSout[, 1])
    Vmaxout <- ChooseVmaxlen(Fleet, Stock, isRel)
    tempmaxlen[N, ] <- sort(Vmaxout[, 2])
    if (isRel) {
      polygon(x = c(0, max(tempL5[N, ]), max(tempLFS[N, ]), 3, rev(c(0, 
                                                                     min(tempL5[N, ]), min(tempLFS[N, ]), 3))), y = c(0, 0.05, 
                                                                                                                      1, min(tempmaxlen[N, ]), rev(c(0, 0.05, 1, max(tempmaxlen[N, 
                                                                                                                                                                                ])))), col = "grey")
      par(ask = TRUE)
    } else {
      polygon(x = c(0, max(tempL5[N, ]), max(tempLFS[N, ]), mean(Stock@Linf), 
                    rev(c(0, min(tempL5[N, ]), min(tempLFS[N, ]), mean(Stock@Linf)))), 
              y = c(0, 0.05, 1, min(tempmaxlen[N, ]), rev(c(0, 0.05, 
                                                            1, max(tempmaxlen[N, ])))), col = "grey")
      par(ask = TRUE)
    }
  }
  par(set.par)
  # CheckSelect(Fleet, Stock)
  Fleet@LFS <- rep(-9999, 2) # code for clearly non-sensical value 
  Fleet@L5 <- rep(-9999, 2) 
  Fleet@Vmaxlen <- rep(-9999, 2)
  Fleet@L5Lower <- tempL5[, 1]
  Fleet@L5Upper <- tempL5[, 2]
  Fleet@LFSLower <- tempLFS[, 1]
  Fleet@LFSUpper <- tempLFS[, 2]
  Fleet@VmaxLower <- tempmaxlen[, 1]
  Fleet@VmaxUpper <- tempmaxlen[, 2]
  Fleet
}

#' Internal function to create a blank plot for mapping selectivity at length
#' 
#' @usage BlankSelPlot(Stock = NULL, Yr = NULL, N = NULL, isRel)
#' @param Stock Stock object. If provided, average length-at-maturity
#' is included on plot for reference.
#' @param Yr Year number
#' @param N Year index
#' @param isRel Logical. Are selectivity parameters relative to length of maturity?
#' @author A. Hordyk
#' @keywords internal
#' @export BlankSelPlot
BlankSelPlot <- function(Stock = NULL, Yr = NULL, N = NULL, isRel) {
  if (isRel) {
    Max <- 3
    AxCex <- 1.3
    By <- 0.05
    par(mfrow = c(1, 1), mai = c(2, 1, 0.5, 0.3), oma = c(1, 1, 1, 1))
    plot(c(0, 3), c(0, 1), type = "n", xlab = "", ylab = "", axes = FALSE)
    mtext(side = 2, line = 3, "Selectivity", cex = AxCex)
    axis(side = 2)
    Xax <- seq(from = 0, to = Max - By, by = 2 * By)
    axis(side = 1, at = Xax)
    axis(side = 1, at = c(2, Max), labels = c("", "Lmax"), xpd = NA)
    mtext(side = 1, line = 3.5, "Relative Length", cex = AxCex)
    axis(side = 1, at = 1, line = 1.5, labels = "L50")
    if (!is.null(Stock) & class(Stock) == "Stock") {
      L50 <- mean(Stock@L50)  # mean length at maturity
      MatAx <- L50 * Xax
      axis(side = 1, line = 5.5, at = Xax, labels = MatAx)
      axis(side = 1, line = 5.5, at = c(2, Max), labels = c("", "Lmax"), 
           xpd = NA)
      mtext(side = 1, line = 8.5, "Approx. Length", cex = AxCex)
      axis(side = 1, at = 1, line = 6.5, labels = "Mean L50")
    }
  } else {
    Max <- max(Stock@Linf)
    AxCex <- 1.3
    By <- 2.5
    par(mfrow = c(1, 1), mai = c(2, 1, 0.5, 0.3), oma = c(1, 1, 1, 1))
    plot(c(0, Max), c(0, 1), type = "n", xlab = "", ylab = "", axes = FALSE)
    mtext(side = 2, line = 3, "Selectivity", cex = AxCex)
    axis(side = 2)
    Xax <- seq(from = 0, to = Max, by = 2 * By)
    axis(side = 1, at = Xax, labels = Xax)
    mtext(side = 1, line = 3.5, "Length", cex = AxCex)
    abline(v=min(Stock@L50), lty=2)
    abline(v=max(Stock@L50), lty=2)
    text(mean(Stock@L50), 0.1, "L50")
  }
  if (N == 1) {
    title(paste("Choose selectivity points for Year", Yr, "(First Year)"))
  } else {
    title(paste("Choose selectivity points for Year", Yr))
  }
}

# #' Internal function for mapping selectivity
# #'
# #' @author A. Hordyk
# #' @keywords internal
ChooseL5 <- function(Fleet, Stock, isRel) {
  if (isRel) {
    By <- 0.05
    Xs <- seq(from = 0, to = 1.5, by = By)
    Ys <- rep(0.05, length(Xs))
    points(Xs, Ys, col = "gray", cex = 0.5)
    text(0.5, 0.2, "Choose two points for L5")
    L5out <- identifyPch(x = Xs, y = Ys, tolerance = 0.1, n = 2)
  } else {
    By <- 2.5
    Xs <- seq(from = 0, to = max(Stock@L50) * 1.5, by = By)
    Ys <- rep(0.05, length(Xs))
    points(Xs, Ys, col = "gray", cex = 0.5)
    text(Xs[5], 0.2, "Choose two points for L5", pos = 4)
    L5out <- identifyPch(x = Xs, y = Ys, tolerance = 0.1, n = 2)
  }
  L5out
}

# #' Internal function for mapping selectivity
# #'
# #' @author A. Hordyk
# #' @keywords internal
ChooseLFS <- function(L5out, Fleet, Stock, isRel) {
  if (isRel) {
    Max <- 3
    By <- 0.05
    Xs <- seq(from = max(L5out[, 1]), to = Max, by = By)
    Ys <- rep(1, length(Xs))
    points(Xs, Ys, col = "gray", cex = 0.5)
    text(1.5, 0.5, "Choose two points for LFS", pos = 4)
    LFSout <- identifyPch(x = Xs, y = Ys, tolerance = 0.1, n = 2)
  } else {
    Max <- 3
    By <- 2.5
    Xs <- seq(from = max(L5out[, 1]), to = max(Stock@L50) * 2, by = By)
    Ys <- rep(1, length(Xs))
    points(Xs, Ys, col = "gray", cex = 0.5)
    text(Xs[5], 0.5, "Choose two points for LFS", pos = 4)
    LFSout <- identifyPch(x = Xs, y = Ys, tolerance = 0.1, n = 2)
  }
  LFSout
}

# #' Internal function for mapping selectivity
# #'
# #' @author A. Hordyk
# #' @keywords internal
ChooseVmaxlen <- function(Fleet, Stock, isRel) {
  if (isRel) {
    Max <- 3
    Ys <- seq(from = 0, to = 1, by = 0.05)
    Xs <- rep(Max, length(Ys))
    points(Xs, Ys, col = "gray", cex = 0.5)
    text(2, 0.8, "Choose two points for selectivity\n at maximum length")
    Vmaxout <- identifyPch(x = Xs, y = Ys, tolerance = 0.1, n = 2)
  } else {
    Max <- mean(Stock@Linf)
    By <- 2.5
    Ys <- seq(from = 0, to = 1, by = 0.05)
    Xs <- rep(Max, length(Ys))
    points(Xs, Ys, col = "gray", cex = 0.5)
    text(Xs[1], 0.8, "Choose two points for selectivity\n at maximum length", 
         pos = 2)
    Vmaxout <- identifyPch(x = Xs, y = Ys, tolerance = 0.1, n = 2)
  }
  Vmaxout
}

# #' Rough Plot of Historical Selectivity Patterns
# #'
# #' @author A. Hordyk
# #' @keywords internal
CheckSelect <- function(Fleet, Stock = NULL) {
  # NEEDS TO BE FIXED
  if (length(Fleet@SelYears) < 1) 
    stop("No break points in selectivity pattern")
  n <- length(Fleet@SelYears)
  if (n < 4) 
    par(mfrow = c(n, 1), mar = c(4, 4, 1, 1), oma = c(2, 3, 1, 1), 
        bty = "l")
  if (n >= 4) {
    N <- ceiling(n/2)
    par(mfrow = c(N, 2), mar = c(4, 4, 1, 1), oma = c(2, 3, 1, 1), 
        bty = "l")
  }
  for (X in 1:n) {
    plot(c(0, 3), c(0, 1), type = "n", xlab = "", ylab = "")
    if (length(Fleet@AbsSelYears) > 0) 
      title(Fleet@AbsSelYears[X])
    if (length(Fleet@AbsSelYears) == 0) 
      title(Fleet@SelYears[X])
    polygon(x = c(0, max(Fleet@L5[X, ]), max(Fleet@LFS[X, ]), 3, rev(c(0, 
                                                                       min(Fleet@L5[X, ]), min(Fleet@LFS[X, ]), 3))), y = c(0, 0.05, 
                                                                                                                            1, min(Fleet@Vmaxlen[X, ]), rev(c(0, 0.05, 1, max(Fleet@Vmaxlen[X, 
                                                                                                                                                                                            ])))), col = "grey")
    lines(c(1, 1), c(0, 1), lty = 3)
    text(1.1, 0.2, "L50", cex = 1.25)
  }
  mtext(side = 2, outer = TRUE, "Selectivity", cex = 1.25, xpd = NA)
  mtext(side = 1, outer = TRUE, "Relative Length", cex = 1.25, xpd = NA)
}



