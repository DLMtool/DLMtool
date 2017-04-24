
#' Manually choose the historical selectivity pattern
#' 
#' Input the first historical year, and all years where selectivity pattern
#' changed (separated by comma).  Interactive plot which allows users to
#' specify a range for the length at 5\% and full selection (LFS), as well as
#' selectivity at maximum length for each year.  Produces a simple plot which
#' shows the range in selectivity pattern for each break-point year.
#' Selectivity-at-length is fixed in between break-point years.  Note that this
#' function replaces 'nyears' in the Fleet object with the value defined here
#' (FstYr:current year).
#' 
#' @param Fleet A fleet object.
#' @param Stock Optional Stock object. If provided, average length-at-maturity
#' is included on plot for reference.
#' @param FstYr Optional value for first historical year. If empty, user must
#' specify the year in console.
#' @param SelYears Optional vector of values for each year where selectivity
#' pattern changed. If empty, user must specify the years in console (comma
#' separated).
#' @author A. Hordyk
#' @importFrom utils flush.console
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
        par(mfrow = c(N, 2), , mar = c(4, 4, 1, 1), oma = c(2, 3, 1, 1), 
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

