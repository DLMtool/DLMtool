# ---- Plot MSE object ----
#' Plot MSE object
#' @param x object of class MSE
#' @param ... other parameters passed to plot (currently ignored)
#' @export
plot.MSE <- function(x, ...) {
  Pplot(x, nam = deparse(substitute(x)))
  Kplot(x)
  Tplot(x, nam = deparse(substitute(x)))
}
     



#' Plot several plots with a shared legend
#'
#' @param plots list of plot objects of class `gg` or `ggplot`  
#' @param ncol Optional number of columns
#' @param nrow Optional number of rows
#' @param position position of the legend ("bottom" or "right")
#' @param legend Logical. Use a legend?
#' @export
#'
#' @note modified from https://github.com/tidyverse/ggplot2/wiki/share-a-legend-between-two-ggplot2-graphs
join_plots <- function(plots, ncol = length(plots), nrow = 1, position = c("right", "bottom"),
                       legend=TRUE) {
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package \"gridExtra\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  position <- match.arg(position)
  if (legend) {
    g <- ggplot2::ggplotGrob(plots[[1]] + ggplot2::theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + ggplot2::theme(legend.position="none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    combined <- switch(position,
                       "bottom" = gridExtra::arrangeGrob(do.call(gridExtra::arrangeGrob, gl),
                                                         legend,
                                                         ncol = 1,
                                                         heights = grid::unit.c(grid::unit(1, "npc") - lheight, lheight)),
                       "right" = gridExtra::arrangeGrob(do.call(gridExtra::arrangeGrob, gl),
                                                        legend,
                                                        ncol = 2,
                                                        widths = grid::unit.c(grid::unit(1, "npc") - lwidth, lwidth)))
  } else {
    gl <- lapply(plots, function(x) x + ggplot2::theme(legend.position="none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    combined <- gridExtra::arrangeGrob(do.call(gridExtra::arrangeGrob, gl))
  }

  grid::grid.newpage()
  grid::grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}


#' Plot a barplot of MSE results
#' 
#' @param height An object of class MSE. Generic function must have argument
#' height. But note that this must be an MSE object.
#' @param MSEobj Optional. An object of class MSE. Overides \code{height}
#' @param PMs List of performance metrics. Options are \code{c('SSB_SSB0',
#' 'B_BMSY', 'F_FMSY', 'AAVE', 'AAVY')}
#' @param PLim Probability threshold
#' @param lastYrs Last number of years in projection to calculate statistics
#' @param maxMP Maximum number of MPs to include in each plot
#' @param MPs Optional subset MSE object by MP
#' @param Title Optional title for plot
#' @param sims Optional subset MSE object by simulation
#' @param msg Logical. Print out messages?
#' @param incRef Logical. Include the reference methods?
#' @param cex.names Size of names
#' @param ...  Optional additional arguments passed to \code{barplot}
#' @author A. Hordyk
#' @export
barplot.MSE <- function(height, MSEobj = NULL, PMs = list(B_BMSY = 0.5, 
                                                          SSB_SSB0 = 0.2), PLim = 0.8, lastYrs = 10, maxMP = 14, MPs = NA, Title = NULL, 
                        sims = NULL, msg = TRUE, cex.names = 1.3, incRef = FALSE, ...) {
  
  MSEobj <- match.arg(MSEobj)
  if (is.null(MSEobj)) 
    MSEobj <- height
  if (!is.null(sims) & all(is.na(MPs))) 
    MSEobj <- Sub(MSEobj, sims = sims)
  if (!is.null(sims) & all(!is.na(MPs))) 
    MSEobj <- Sub(MSEobj, sims = sims, MPs = MPs)
  if (is.null(sims) & !all(is.na(MPs))) 
    MSEobj <- Sub(MSEobj, MPs = MPs)
  
  if (!incRef) {
    mps <- MSEobj@MPs[!grepl("ref", MSEobj@MPs)]
    MSEobj <- Sub(MSEobj, MPs = mps)
  }
  DF <- list()
  if (length(lastYrs) > 1) {
    for (xx in seq_along(lastYrs)) DF[[xx]] <- MPStats(MSEobj, PMRefs = PMs, 
                                                       lastYrs = lastYrs[xx], UseMean = TRUE, msg = msg)$Perf
  } else {
    DF[[1]] <- MPStats(MSEobj, PMRefs = PMs, lastYrs = lastYrs, UseMean = TRUE, 
                       msg = msg)$Perf
  }
  lastYrs[lastYrs >= MSEobj@proyears] <- 10
  PosPMs <- c("SSB_SSB0", "B_BMSY", "F_FMSY", "AAVE", "AAVY")
  PMNames <- names(PMs)
  ind <- match(PMNames, PosPMs)
  
  pms <- paste0(names(PMs), "p")
  
  vars <- which(!is.na(match(names(DF[[1]]), pms)))
  
  temp <- lapply(DF, "[", vars)
  DF2 <- cbind(DF[[1]][, 1], do.call(cbind.data.frame, temp))
  names(DF2)[1] <- "MP"
  
  B0Ref <- unique(DF[[1]]$SSB_SSB0Ref)
  B_BMSYRef <- unique(DF[[1]]$B_BMSYRef)
  F_FMSYRef <- unique(DF[[1]]$F_FMSYRef)
  AAVERef <- unique(DF[[1]]$AAVERef)
  AAVYRef <- unique(DF[[1]]$AAVYRef)
  
  Years <- paste("Years", (MSEobj@proyears - lastYrs) + 1, "-", MSEobj@proyears)
  
  B0Leg <- bquote(italic(B) > ~.(B0Ref) ~ italic(B[0]))
  BMSYLeg <- bquote(italic(B) > ~.(B_BMSYRef) ~ italic(B[MSY]))
  FMSYLeg <- bquote(italic(F) < ~.(F_FMSYRef) ~ italic(F[MSY]))
  AAVELeg <- bquote(italic(AAVE) > ~.(AAVERef) ~ "%")
  AAVYRef <- bquote(italic(AAVY) > ~.(AAVYRef) ~ "%")
  LegList <- list(B0Leg, BMSYLeg, FMSYLeg, AAVELeg, AAVYRef)
  Legend <- NULL
  Legend <- append(Legend, as.expression(LegList[ind]))
  
  if (length(Years) > 1) {
    temp <- NULL
    count <- 1
    for (xx in seq_along(Years)) {
      for (yy in seq_along(LegList[ind])) {
        temp[[count]] <- bquote(.(LegList[ind][[yy]]) ~ (Years ~ 
                                                           .((MSEobj@proyears - lastYrs[xx]) + 1) ~ "-" ~ .(MSEobj@proyears)))
        count <- count + 1
      }
    }
    Legend <- append(NULL, as.expression(temp))
  }
  
  MPs <- as.character(DF2$MP)
  nMPs <- length(MPs)
  if (is.null(MPs)) 
    stop("Dataframe must a column named `MP`")
  ndat <- ncol(DF2) - 1
  # if (length(ProbLims) != ndat) stop('Must be a probablility limit for
  # each variable')
  Probs <- t(DF2[, 2:(ndat + 1)])
  colnames(Probs) <- MPs
  if (max(Probs) <= 1) 
    Probs <- Probs * 100
  ProbLims <- VLine <- PLim
  if (max(ProbLims) <= 1) 
    ProbLims <- ProbLims * 100
  if (nrow(Probs) > 1) {
    Ord <- order(apply(Probs, 2, mean))
    Probs <- Probs[, Ord]
    MPnames <- colnames(Probs)
    Pass <- as.logical(apply(Probs >= ProbLims, 2, prod))
  } else {
    MPnames <- colnames(Probs)
    Ord <- order(Probs)
    Probs <- t(as.matrix(Probs[Ord, drop = FALSE]))
    MPnames <- MPnames[Ord]
    colnames(Probs) <- MPnames
    Pass <- as.logical(Probs >= ProbLims)
  }
  
  # if (length(ProbLims) == 1) Pass <- Probs > ProbLims
  
  # bcols <- brewer.pal(8, 'Dark2')[1:length(vars)]
  cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", 
            "#A6761D", "#666666")
  bcols <- cols[1:ndat]
  
  nplots <- ceiling(nMPs/maxMP)
  Ncol <- ceiling(sqrt(nplots))
  Nrow <- ceiling(nplots/Ncol)
  if (!(Ncol * Nrow >= nplots)) 
    stop("Error in number of plots (you've found a bug in the function!)")
  
  tempmat <- matrix(1:(Ncol * Nrow), nrow = Nrow, byrow = TRUE)
  bspace <- max(nchar(MPs)) - 1
  
  op <- par(mfrow = c(Nrow, Ncol), oma = c(2.5, 1, 5, 2), mar = c(3, bspace, 
                                                                  0, 0))
  if (nplots == 1) {
    barplot(Probs, beside = TRUE, horiz = TRUE, names = MPnames, las = 1, 
            xlim = c(0, 100), cex.axis = 1.5, cex.lab = 2, col = bcols, 
            xlab = "Probability", xpd = NA, cex.names = cex.names, ...)
    if (!is.null(VLine)) {
      if (VLine < 1) 
        VLine <- VLine * 100
      if (VLine > 0) 
        abline(v = VLine, lty = 2, lwd = 2, col = "gray")
    }
    rng <- par("usr")
    lg <- legend(rng[1], rng[2], legend = Legend, bty = "n", cex = 1.25, 
                 horiz = TRUE, xpd = NA, title = "", plot = FALSE)
    legend(rng[1], rng[4] + lg$rect$h, legend = Legend, bty = "n", 
           cex = 1.25, horiz = TRUE, xpd = NA, fill = bcols, title = "")
  }
  # Split the MPs over multiple plots
  if (length(vars) == 1) 
    colnames(Probs) <- MPnames
  if (nplots > 1) {
    npplot <- min(ceiling(nMPs/nplots), maxMP)  # number per plot 
    xx <- 1
    xx2 <- npplot
    for (X in 1:nplots) {
      # if (length(Legend) > 1) splitdat <- Probs[,xx:xx2] if (length(Legend)
      # == 1) splitdat <- Probs[,xx:xx2]
      splitdat <- Probs[, xx:xx2]
      if (X == 1) {
        tt <- barplot(splitdat, beside = TRUE, ylab = "", xlim = c(0, 
                                                                   100), col = bcols, names.arg = NULL, cex.axis = 1.5, 
                      cex.lab = 2, las = 1, xpd = NA, horiz = TRUE, cex.names = cex.names, 
                      ...)
        if (!is.null(VLine)) {
          if (VLine < 1) 
            VLine <- VLine * 100
          if (VLine > 0) 
            abline(v = VLine, lty = 2, lwd = 2, col = "gray")
        }
        rng <- par("usr")
        lg <- legend(rng[1], rng[2], legend = Legend, bty = "n", 
                     cex = 1.25, horiz = TRUE, xpd = NA, fill = bcols, title = "", 
                     plot = FALSE)
        legend(rng[1], rng[4] + lg$rect$h, legend = Legend, bty = "n", 
               cex = 1.5, horiz = TRUE, xpd = NA, fill = bcols, title = "", 
               ...)
      } else {
        barplot(splitdat, beside = TRUE, ylab = "", xlim = c(0, 
                                                             100), names.arg = NULL, cex.axis = 1.5, cex.lab = 2, 
                las = 1, axes = TRUE, col = bcols, horiz = TRUE, cex.names = cex.names, 
                ...)
        # abline(v=PLim, lty=2, col='gray', xpd=FALSE)
      }
      xx <- xx2 + 1
      xx2 <- min(xx + npplot - 1, nMPs)
      mtext(side = 1, outer = TRUE, "Probability", line = 1, cex = 1.4)
      if (!is.null(VLine)) {
        if (VLine < 1) 
          VLine <- VLine * 100
        if (VLine > 0) 
          abline(v = VLine, lty = 2, lwd = 2, col = "gray")
      }
    }
  }
  if (!is.null(Title)) 
    mtext(side = 3, outer = TRUE, title, line = 3.5, cex = 1.25)
  
  if (length(Years) == 1) 
    mtext(side = 3, outer = TRUE, Years, line = 1.5)
  Pout <- t(Probs)
  if (is.null(rownames(Pout))) {
    MP <- colnames(Pout)
    colnames(Pout) <- NULL
  } else {
    MP <- rownames(Pout)
    rownames(Pout) <- NULL
  }
  MPclass <- MPtype(MPs)[,2]
  par(op)
  invisible(data.frame(MP = MP, Pout, Pass = Pass, MPClass = MPclass, 
                       stringsAsFactors = FALSE))
}


# #' Boxplot of MP performance from MSE object
# #' 
# #' @param x An object of class MSE
# #' @param MPs Optional subset MSE object by MP
# #' @param maxMP Maximum number of MPs to plot
# #' @param PMRefs List containing the Performance Metrics reference points.
# #' Options are \code{'SSB_SSB0', 'B_BMSY', 'F_FMSY', 'AAVE', 'AAVY'}
# #' @param lastYrs Last number of years in projection to calculate statistics
# #' @param cex.lab Size of axis label text
# #' @param cex.PM Size of performacne metric text
# #' @param canMPs Optional character vector of MPs that can be applied (plotted
# #' in different colour)
# #' @param cols Optional vector of colours
# #' @param outline Logical. Include outliers in boxplot?
# #' @param CexName Size of the names
# #' @param incLine Logical. Include vertical line?
# #' @param incref Logical. Include reference methods?
# #' @param Names Logical. Include MP names in plot?
# #' @param ...  Additional arguments to be passed to plotting functions
# #' @author A. Hordyk
# #' @export
# boxplot.MSE <- function(x, MPs = NA, maxMP = 8, 
#                         PMRefs = list(B_BMSY = 1, 
#                                       SSB_SSB0 = 0.2, F_FMSY = 1, AAVY = 30, AAVE = 30), 
#                         lastYrs = 10, cex.lab = 1.2, 
#                         cex.PM = 0.75, canMPs = NULL, cols = TRUE, outline = FALSE, CexName = 1.25, 
#                         incLine = TRUE, incref = FALSE, Names = TRUE, ...) {
#   
#   MSEobj <- x
#   if (!all(is.na(MPs))) {
#     if (!incref) {
#       ind <- grep("FMSYref", MPs)  # FMSYref methods 
#       ind <- c(ind, grep("NFref", MPs))  # No fishing reference
#       if (length(ind) > 0) 
#         MPs <- MPs[-ind]
#     }
#     MSEobj <- Sub(MSEobj, MPs = MPs)
#   }
#   nmp <- MSEobj@nMPs
#   Nyears <- MSEobj@nyears
#   Pyears <- MSEobj@proyears
#   nMPs <- MSEobj@nMPs
#   MPs <- MSEobj@MPs
#   nsim <- MSEobj@nsim
#   RefYd <- MSEobj@OM$RefY
#   nyrs <- MSEobj@proyears
#   
#   perf <- MPStats(MSEobj, lastYrs = lastYrs)
#   if (lastYrs >= MSEobj@proyears) 
#     lastYrs <- 10
#   yrs <- (nyrs - lastYrs + 1):nyrs
#   perfdat <- perf[[1]]
#   Years <- paste("Years", (MSEobj@proyears - lastYrs) + 1, "-", MSEobj@proyears, 
#                  "(last", lastYrs, "years)")
#   
#   # MPtype <- perfdat$MPtype
#   MPtypes <- MPtype(MPs)[,2]
#   outmps <- perfdat[MPtypes == "Output", ]$MP
#   inmps <- perfdat[MPtypes == "Input", ]$MP
#   nOut <- ceiling(nmp/2)
#   nIn <- floor(nmp - nOut)
#   if (nMPs > nmp) {
#     outDist <- perfdat[MPtypes == "Output", ]$Dist
#     OutMPs <- perfdat[MPtypes == "Output", ]$MP[order(outDist)[1:nOut]]
#     inDist <- perfdat[MPtypes == "Input", ]$Dist
#     InMPs <- perfdat[MPtypes == "Input", ]$MP[order(inDist)[1:nIn]]
#   } else {
#     OutMPs <- perfdat[MPtypes == "Output", ]$MP
#     InMPs <- perfdat[MPtypes == "Input", ]$MP
#   }
#   
#   mseobj <- Sub(MSEobj, MPs = c(OutMPs, InMPs))
#   nMPs <- mseobj@nMPs
#   perf <- MPStats(mseobj, msg = FALSE)
#   perfdat <- perf$Perf
#   MPtypes <- MPtype(mseobj@MPs)[,2]
#   # MPtype <- perfdat$MPtype
#   rawVals <- perf$BySim
#   MPs <- perfdat$MP
#   
#   # Reporting Distribution in last years
#   BBMSY <- rawVals$B_BMSY[, , yrs]
#   BBMSY <- aperm(BBMSY, c(1, 3, 2))
#   dim(BBMSY) <- c(nsim * lastYrs, nMPs)
#   
#   BB0 <- rawVals$SSB_SSB0[, , yrs]
#   BB0 <- aperm(BB0, c(1, 3, 2))
#   dim(BB0) <- c(nsim * lastYrs, nMPs)
#   
#   FFMSY <- rawVals$F_FMSY[, , yrs]
#   FFMSY <- aperm(FFMSY, c(1, 3, 2))
#   dim(FFMSY) <- c(nsim * lastYrs, nMPs)
#   
#   LTY <- rawVals$LTY
#   LTY <- aperm(LTY, c(1, 3, 2))
#   dim(LTY) <- c(nsim * 10, nMPs)
#   
#   AAVY <- rawVals$AAVY
#   AAVY <- aperm(AAVY, c(1, 3, 2))
#   dim(AAVY) <- c(nsim * 10, nMPs)
#   
#   AAVE <- rawVals$AAVE
#   # apply(AAVE[,,10], 2, mean, na.rm=TRUE)
#   AAVE <- aperm(AAVE, c(1, 3, 2))
#   dim(AAVE) <- c(nsim * 10, nMPs)
#   # apply(AAVE, 2, mean, na.rm=TRUE)
#   
#   # Colors - to make the plot a bit more cheerful
#   inputs <- which(MPtypes == "Input")
#   outputs <- which(MPtypes == "Output")
#   fonts <- rep(2, nmp)
#   fonts[MPtypes == "Input"] <- 4
#   NameCol <- rep("black", nmp)
#   NameCol[MPs %in% canMPs] <- "green"
#   index <- grep("ref", MPs)
#   NameCol[index] <- "darkgray"
#   
#   if (class(cols) == "character") 
#     Col <- cols
#   if (class(cols) == "logical" && !(cols)) 
#     Col <- "lightblue"
#   
#   if (class(cols) == "logical" && cols) {
#     cols1 <- rep(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
#                    "#FFFF33", "#A65628", "#F781BF"), 4)  # brewer.pal(8, 'Set1')
#     cols2 <- rep(c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", 
#                    "#FFD92F", "#E5C494", "#B3B3B3"), 4)  # brewer.pal(8, 'Set2')
#     Col <- rep(cols1, length.out = nMPs)[1:nMPs]
#     if (length(inputs) > 0) 
#       Col[inputs] <- cols2[1:length(inputs)]
#   }
#   
#   Line <- 3.5
#   CexName <- 1.25
#   lncol <- "slategray"
#   topcex <- cex.PM
#   ncharmp <- max(nchar(MPs)) - 1.5
#   multi <- 1.15
#   op <- par(mfrow = c(2, 3), bty = "n", oma = c(4, ncharmp, 2, 0), mar = c(3, 
#                                                                            2, 4, 3))
#   # if (class(Title) == 'character') par(mfrow=c(2,3), bty='n',
#   # oma=c(4,ncharmp,2,0), mar=c(3,2,4,3)) if (class(Title) !=
#   # 'character') par(mfrow=c(2,3), bty='n', oma=c(4,ncharmp,0,0),
#   # mar=c(3,2,4,3)) B/B0
#   YLim <- c(0, 1)
#   tt <- boxplot(BB0, names = FALSE, las = 1, outline = outline, horizontal = TRUE, 
#                 xlab = "", bty = "n", col = Col, xpd = NA, cex.names = CexName, 
#                 axes = FALSE, ylim = YLim)
#   axis(side = 1, cex.axis = CexName)
#   axis(side = 2, at = outputs, labels = FALSE, las = 2)
#   axis(side = 2, at = inputs, labels = FALSE, las = 2)
#   if (Names) 
#     text(x = rep(-0.07, nmp), y = 1:nmp, labels = MPs, col = NameCol, 
#          cex = CexName, font = fonts, pos = 2, xpd = NA, srt = 0)
#   mtext(side = 1, expression(italic(B/B[0])), cex = cex.lab, line = Line)
#   mtext(side = 3, outer = TRUE, Years, line = -0.5)
#   PosPMs <- c("SSB_SSB0", "B_BMSY", "F_FMSY", "AAVE", "AAVY")
#   ind <- match(names(PMRefs), PosPMs)
#   refs <- unlist(PMRefs[ind])
#   
#   if (incLine && refs[1] > 0) {
#     ps <- round(apply(BB0 > refs[1], 2, sum, na.rm = TRUE)/(nrow(BB0)), 
#                 2)
#     font <- rep(1, length(ps))
#     ps <- sprintf("%3.2f", ps)
#     if (refs[1] != 1) 
#       txt <- bquote(Prob. ~ italic(B) > ~.(refs[1]) ~ italic(B[0]))
#     if (refs[1] == 1) 
#       txt <- bquote(Prob. ~ italic(B) > ~italic(B[0]))
#     
#     text(x = YLim[2] * 1.1, y = 1:nMPs, labels = ps, xpd = NA, font = font)
#     text(x = YLim[2] * 1.2, y = nMPs * multi, txt, xpd = NA, pos = 2)
#     abline(v = refs[1], lty = 2, lwd = 2, col = lncol)
#   }
#   
#   text(-0.2, nmp * 1.25, "Input Control", font = 4, xpd = NA, cex = 1.2, 
#        pos = 4)
#   text(-0.2, nmp * 1.35, "Output Control", font = 2, xpd = NA, cex = 1.2, 
#        pos = 4)
#   if (!is.null(canMPs)) {
#     text(1, nmp * 1.25, "Available", font = 2, xpd = NA, cex = 1.2, 
#          col = "green", pos = 2)
#     text(1, nmp * 1.35, "Not Available", font = 2, xpd = NA, cex = 1.2, 
#          pos = 2)
#   }
#   
#   # B/BMSY
#   YLim <- c(0, 3)
#   tt <- boxplot(BBMSY, names = FALSE, las = 1, outline = outline, horizontal = TRUE, 
#                 xlab = "", bty = "n", col = Col, xpd = NA, axes = FALSE, ylim = YLim)
#   
#   mtext(side = 1, expression(italic(B/B[MSY])), cex = cex.lab, line = Line)
#   axis(side = 1, cex.axis = CexName)
#   axis(side = 2, at = outputs, labels = FALSE)
#   axis(side = 2, at = inputs, labels = FALSE)
#   if (incLine && refs[2] > 0) {
#     ps <- round(apply(BBMSY > refs[2], 2, sum, na.rm = TRUE)/(nrow(BBMSY)), 
#                 2)
#     ps <- sprintf("%3.2f", ps)
#     axis(side = 4, at = 1:nMPs, labels = ps, las = 1, tick = FALSE, 
#          line = -1, xpd = NA)
#     if (refs[2] != 1) 
#       txt <- bquote(Prob. ~ italic(B) > ~.(refs[2]) ~ italic(B[MSY]))
#     if (refs[2] == 1) 
#       txt <- bquote(Prob. ~ italic(B) > ~italic(B[MSY]))
#     # text(x=YLim[2]*1.1, y=1:nMPs, labels=ps, xpd=NA, font=font)
#     text(x = YLim[2] * 1.2, y = nMPs * multi, txt, xpd = NA, pos = 2)
#     abline(v = refs[2], lty = 2, lwd = 2, col = lncol)
#   }
#   
#   
#   # F/FMSY
#   YLim <- c(0, 2)
#   tt <- boxplot(FFMSY, names = FALSE, las = 1, outline = outline, horizontal = TRUE, 
#                 xlab = "", bty = "n", col = Col, xpd = NA, axes = FALSE, ylim = YLim)
#   
#   mtext(side = 1, expression(italic(F/F[MSY])), cex = cex.lab, line = Line)
#   axis(side = 1, cex.axis = CexName)
#   axis(side = 2, at = outputs, labels = FALSE)
#   axis(side = 2, at = inputs, labels = FALSE)
#   if (incLine && refs[3] > 0) {
#     ps <- round(apply(FFMSY < refs[3], 2, sum, na.rm = TRUE)/(nrow(FFMSY)), 
#                 2)
#     ps <- sprintf("%3.2f", ps)
#     axis(side = 4, at = 1:nMPs, labels = ps, las = 1, tick = FALSE, 
#          line = -1, xpd = NA)
#     abline(v = refs[3], lty = 2, lwd = 2, col = lncol)
#     if (refs[3] != 1) 
#       txt <- bquote(Prob. ~ italic(F) < ~.(refs[3]) ~ italic(F[MSY]))
#     if (refs[3] == 1) 
#       txt <- bquote(Prob. ~ italic(F) < ~italic(F[MSY]))
#     text(x = YLim[2] * 1.2, y = nMPs * multi, txt, xpd = NA, pos = 2)
#   }
#   
#   
#   # AAVY
#   if (refs[4] < 1) 
#     refs[4] <- refs[4] * 100
#   AAVY <- AAVY * 100
#   YLim <- c(0, 60)
#   tt <- boxplot(AAVY, names = FALSE, las = 1, outline = outline, horizontal = TRUE, 
#                 xlab = "", bty = "n", col = Col, xpd = NA, axes = FALSE, ylim = YLim)
#   mtext(side = 1, "Average Annual\n Variation Yield (%)", cex = cex.lab, 
#         line = Line + 1)
#   axis(side = 1, cex.axis = CexName)
#   axis(side = 2, at = outputs, labels = FALSE, las = 2)
#   axis(side = 2, at = inputs, labels = FALSE, las = 2)
#   text(x = rep(-5, nmp), y = 1:nmp, labels = MPs, col = NameCol, cex = CexName, 
#        font = fonts, pos = 2, xpd = NA, srt = 0)
#   if (incLine && refs[4] > 0) {
#     ps <- round(apply(AAVY < refs[4], 2, sum, na.rm = TRUE)/(nrow(AAVY)), 
#                 2)
#     ps <- sprintf("%3.2f", ps)
#     axis(side = 4, at = 1:nMPs, labels = ps, las = 1, tick = FALSE, 
#          line = -1, xpd = NA)
#     abline(v = refs[4], lty = 2, lwd = 2, col = lncol)
#     if (refs[4] != 1) 
#       txt <- bquote(Prob. ~ italic(AAVY) < ~.(refs[4]) ~ "%")
#     if (refs[4] == 1) 
#       txt <- bquote(Prob. ~ italic(AAVY) < ~"%")
#     text(x = YLim[2] * 1.2, y = nMPs * multi, txt, xpd = NA, pos = 2)
#   }
#   abline(v = refs[4], lty = 2, lwd = 2, col = lncol)
#   
#   # AAVE
#   if (refs[5] < 1) 
#     refs[5] <- refs[5] * 100
#   AAVE <- AAVE * 100
#   tt <- boxplot(AAVE, names = FALSE, las = 1, outline = outline, horizontal = TRUE, 
#                 xlab = "", bty = "n", col = Col, xpd = NA, axes = FALSE, ylim = c(0, 
#                                                                                   60))
#   mtext(side = 1, "Average Annual\n Variation Effort (%)", cex = cex.lab, 
#         line = Line + 1)
#   axis(side = 1, cex.axis = CexName)
#   axis(side = 2, at = outputs, labels = FALSE)
#   axis(side = 2, at = inputs, labels = FALSE)
#   if (incLine && refs[4] > 0) {
#     ps <- round(apply(AAVE < refs[4], 2, sum, na.rm = TRUE)/(nrow(AAVE)), 
#                 2)
#     ps <- sprintf("%3.2f", ps)
#     axis(side = 4, at = 1:nMPs, labels = ps, las = 1, tick = FALSE, 
#          line = -1, xpd = NA)
#     abline(v = refs[4], lty = 2, lwd = 2, col = lncol)
#     if (refs[5] != 1) 
#       txt <- bquote(Prob. ~ italic(AAVE) < ~.(refs[5]) ~ "%")
#     if (refs[5] == 1) 
#       txt <- bquote(Prob. ~ italic(AAVE) < ~"%")
#     text(x = YLim[2] * 1.2, y = nMPs * multi, txt, xpd = NA, pos = 2)
#     
#   }
#   abline(v = refs[4], lty = 2, lwd = 2, col = lncol)
#   
#   
#   # LTY
#   tt <- boxplot(LTY, names = FALSE, las = 1, outline = outline, horizontal = TRUE, 
#                 xlab = "", bty = "n", col = Col, xpd = NA, axes = FALSE)
#   mtext(side = 1, "Relative Long-Term\n Yield", cex = cex.lab, line = Line + 
#           1)
#   axis(side = 1, cex.axis = CexName)
#   axis(side = 2, at = outputs, labels = FALSE)
#   axis(side = 2, at = inputs, labels = FALSE)
#   if (incLine && refs[5] > 0) 
#     abline(v = refs[4], lty = 2, lwd = 2, col = "darkgray")
#   par(op)
#   
# }
# 
# 
# #' Plot the median biomass and yield relative to last historical year
# #' 
# #' Compare median biomass and yield in first year and last 5 years of
# #' projection
# #' 
# #' @param MSEobj An object of class MSE
# #' @param MPs Optional subset by MP
# #' @param lastYrs Last number of years of projection to calculate median
# #' @param XMin Optional minimum for the x-axis
# #' @param YMin Optional minimum for the y-axis
# #' @param ShowLabs Logical. Show the MP labels? Otherwise only plot points
# #' @return Invisibly returns a data frame containing information shown in the
# #' plot
# #' @author A. Hordyk
# #' @export Cplot
# Cplot <- function(MSEobj, MPs = NA, lastYrs = 5, XMin = NULL, YMin = NULL, 
#                   ShowLabs = FALSE) {
#   if (!all(is.na(MPs))) 
#     MSEobj <- Sub(MSEobj, MPs = MPs)
#   nsim <- MSEobj@nsim
#   Alpha <- 60
#   if (nsim < 10) 
#     Alpha <- 180
#   nMPs <- MSEobj@nMPs
#   MPs <- MSEobj@MPs
#   nyears <- MSEobj@nyears
#   proyears <- MSEobj@proyears
#   
#   Stat <- MPStats(MSEobj, lastYrs = lastYrs)$BySim
#   ny <- dim(Stat$Yield)[3]
#   Stat$Yield <- Stat$Yield[, , , drop = FALSE]/Stat$Yield[, , rep(1, 
#                                                                   ny), drop = FALSE]
#   
#   RelYield <- apply(Stat$Yield, 2, median, na.rm = TRUE)
#   
#   Bcurr <- Stat$B_BMSY[, , 1]  # Biomass at start of projections
#   Bend <- apply((Stat$B_BMSY[, , (proyears - lastYrs + 1):proyears]), 
#                 c(1, 2), median, na.rm = TRUE)  # median biomass in last years
#   RelBio <- apply(Bend/Bcurr, 2, median, na.rm = TRUE)
#   
#   XMin <- ifelse(is.null(XMin), 0, XMin)
#   YMin <- ifelse(is.null(YMin), 0, YMin)
#   XLim <- c(YMin, ceiling(max(RelBio)/0.5) * 0.5) * c(0.95, 1.05)
#   YLim <- c(XMin, ceiling(max(RelYield)/0.5) * 0.5) * c(0.95, 1.05)
#   op <- par(mfrow = c(1, 1), oma = c(3, 5, 1, 1), mar = c(2, 2, 0, 0))
#   plot(RelBio, RelYield, xlim = XLim, ylim = YLim, type = "n", bty = "l", 
#        xlab = "", ylab = "", xaxs = "i", yaxs = "i", las = 1)
#   if (ShowLabs) 
#     text(RelBio, RelYield, MSEobj@MPs)
#   if (!ShowLabs) 
#     points(RelBio, RelYield, pch = 21, cex = 2, bg = "lightgray")
#   abline(h = 1, lty = 3, col = "lightgray")
#   abline(v = 1, lty = 3, col = "lightgray")
#   mtext(side = 1, line = 3.5, paste("Median Biomass (last", lastYrs, 
#                                     "years)\n relative to current"), cex = 1.25)
#   mtext(side = 2, line = 3, paste("Median Yield (last", lastYrs, "years)\n relative to current"), 
#         cex = 1.25)
#   par(op)
#   DF <- data.frame(MP = MSEobj@MPs, Biomass = RelBio, Catch = RelYield, 
#                    stringsAsFactors = FALSE)
#   invisible(DF)
#   
# }
# 
# 

# #' Joint probability plot
# #' 
# #' Calculates and plots the joint probability of meeting all performance
# #' metrics simultaneously
# #' 
# #' 
# #' @param MSEobj An object of class MSE
# #' @param PLim Probability limit (acceptable risk threshold; e.g., 0.8 for 80
# #' percent)
# #' @param YVar What to plot of the y-axis: choose from \code{c('LTY', 'STY',
# #' 'avgSSB_SSB0', 'avgB_BMSY')}
# #' @param PMRefs List containing the reference limits for each metric
# #' @param UseMean Logical. Calculate mean (TRUE) or median (FALSE)
# #' @param lastYrs Last number of years in projection period to calculate
# #' summary statistics
# #' @param AvailMPs Optional character vector of available MPs (plotted in a
# #' different colour)
# #' @param XLim Optional limits for the x-axis
# #' @param ShowCols Logical. Show the background colours?
# #' @param ShowLabs Logical. Show the MP labels?
# #' @param All Logical. Plot all MPs (TRUE) or only those above the probability
# #' limit (\code{PLim}))?
# #' @return Invisibly returns data frame containing statistics shown in the plot
# #' @author A. Hordyk
# #' @export Jplot
# Jplot <- function(MSEobj, PLim = 0.8, YVar = c("LTY", "STY", "avgSSB_SSB0", 
#                                                "avgB_BMSY"), PMRefs = list(B_BMSY = 0.5, SSB_SSB0 = 0.2), UseMean = TRUE, 
#                   lastYrs = 10, AvailMPs = NULL, XLim = NULL, ShowCols = TRUE, ShowLabs = FALSE, 
#                   All = TRUE) {
#   
#   nsim <- MSEobj@nsim
#   nMPs <- MSEobj@nMPs
#   MPs <- MSEobj@MPs
#   
#   PMs <- names(PMRefs)
#   YVar <- match.arg(YVar, several.ok = FALSE)
#   mYVar <- YVar
#   mYVar[mYVar == "avgSSB_SSB0"] <- "SSB_SSB0m"
#   mYVar[mYVar == "avgB_BMSY"] <- "B_BMSYm"
#   
#   perf <- MPStats(MSEobj, PMRefs = PMRefs, lastYrs = lastYrs, UseMean = UseMean)
#   Probs <- perf$Probs
#   
#   if (lastYrs >= MSEobj@proyears) 
#     lastYrs <- 10
#   index <- paste0(PMs, "ref")
#   # Joint Prob above B refs and below F
#   jointP <- Probs[[index[1]]]
#   if (length(PMs) > 1) {
#     for (X in 2:length(PMs)) {
#       jointP <- jointP * Probs[[index[X]]]
#     }
#   }
#   JP <- apply(jointP, 2, sum, na.rm = TRUE)/(lastYrs * nsim)
#   
#   # Probability Thresholds
#   vl <- PLim
#   hl <- 0
#   if (vl < 1) 
#     vl <- vl * 100
#   
#   x <- JP
#   if (max(x) <= 1) 
#     x <- x * 100
#   y <- perf[[1]][, mYVar]
#   if (max(y) <= 10) 
#     y <- y * 100
#   adjj <- c(0.9, 1.1)
#   if (is.null(XLim)) 
#     XLim <- c(min(c(-10, min(x, na.rm = T) * adjj)), max(c(max(x, na.rm = T) * 
#                                                              adjj, 110)))
#   if (!All) 
#     XLim <- c(vl * 0.95, 105)
#   YLim <- c(min(c(-10, min(y, na.rm = T) * adjj)), max(c(max(y, na.rm = T) * 
#                                                            adjj, 110)))
#   
#   ## Legend
#   BmsyRef <- unique(perf[[1]]$B_BMSYRef)
#   B0Ref <- unique(perf[[1]]$SSB_SSB0Ref)
#   FRef <- unique(perf[[1]]$F_FMSYRef)
#   # legend text
#   if (FRef == 1) 
#     leg1 <- bquote(italic(F) < ~italic(F[MSY]))
#   if (FRef != 1) 
#     leg1 <- bquote(italic(F) < ~.(FRef) ~ italic(F[MSY]))
#   if (B0Ref == 1) 
#     leg2 <- bquote(italic(B) > ~italic(B[0]))
#   if (B0Ref != 1) 
#     leg2 <- bquote(italic(B) > ~.(B0Ref) ~ italic(B[0]))
#   if (BmsyRef == 1) 
#     leg3 <- bquote(italic(B) > ~~italic(B[MSY]))
#   if (BmsyRef != 1) 
#     leg3 <- bquote(italic(B) > ~.(BmsyRef) ~ italic(B[MSY]))
#   
#   legtex <- list(F_FMSY = leg1, SSB_SSB0 = leg2, B_BMSY = leg3)
#   ind <- match(PMs, names(legtex))
#   Legend <- NULL
#   Legend <- append(Legend, as.expression(legtex[ind]))
#   
#   # Plot
#   xlab <- "Probability > All Performance Limits"
#   if (UseMean) 
#     ylab <- switch(YVar, LTY = paste0("Long-Term Yield (last ", lastYrs, 
#                                       " years)"), STY = paste0("Short-Term Yield (first ", lastYrs, 
#                                                                " years)"), avgB_BMSY = paste0("Mean B/BMSY (%) (last ", lastYrs, 
#                                                                                               " years)"), avgSSB_SSB0 = paste0("Mean B/B0 (%) (last ", lastYrs, 
#                                                                                                                                " years)"))
#   if (!UseMean) 
#     ylab <- switch(YVar, LTY = paste0("Long-Term Yield (last ", lastYrs, 
#                                       " years)"), STY = paste0("Short-Term Yield (first ", lastYrs, 
#                                                                " years)"), avgB_BMSY = paste0("Median B/BMSY (%) (last ", 
#                                                                                               lastYrs, " years)"), avgSSB_SSB0 = paste0("Median B/B0 (%) (last ", 
#                                                                                                                                         lastYrs, " years)"))
#   op <- par(mfrow = c(1, 1), mar = c(4, 4, 1, 1), oma = c(1, 1, 0, 6))
#   plot(NA, xlim = XLim, ylim = YLim, xlab = xlab, ylab = ylab, bty = "l", 
#        las = 1, cex.lab = 1.25)
#   
#   abline(v = vl, col = "#99999940", lwd = 2)
#   Alpha <- 15
#   # polygons
#   LeftCol <- rgb(red = 255, green = 0, blue = 0, alpha = Alpha, names = NULL, 
#                  maxColorValue = 255)
#   RightCol <- rgb(red = 0, green = 255, blue = 0, alpha = Alpha, names = NULL, 
#                   maxColorValue = 255)
#   if (ShowCols) {
#     polygon(x = c(0, vl, vl, 0), y = c(0, 0, hl, hl), col = LeftCol, 
#             border = NA)
#     polygon(x = c(0, vl, vl, 0), y = c(0, 0, max(YLim), max(YLim)), 
#             col = LeftCol, border = NA)
#     polygon(x = c(vl, max(XLim), max(XLim), vl), y = c(0, 0, max(YLim), 
#                                                        max(YLim)), col = RightCol, border = NA)
#     polygon(x = c(vl, max(XLim), max(XLim), vl), y = c(hl, hl, max(YLim), 
#                                                        max(YLim)), col = RightCol, border = NA)
#   }
#   
#   MPs <- MSEobj@MPs
#   Pch <- rep(21, length(MPs))
#   Pch[grep("FMSY", MPs)] <- 24
#   Pch[grep("NFref", MPs)] <- 24
#   
#   # Which MPs meet minimum PMs
#   ind <- which(x >= vl & y >= hl)
#   coly <- rep("darkgray", length(MPs))
#   coly[ind] <- "black"
#   if (!is.null(AvailMPs)) 
#     coly[MPs %in% AvailMPs & !(x >= vl & y >= hl)] <- "Medium Sea Green"
#   if (!is.null(AvailMPs)) 
#     coly[MPs %in% AvailMPs & (x >= vl & y >= hl)] <- "green"
#   coly[grep("FMSY", MPs)] <- "lightgray"
#   coly[grep("NFref", MPs)] <- "lightgray"
#   
#   MPtype <- sapply(1:nMPs, function(X) class(get(MPs[X])))
#   fonts <- rep(2, nMPs)
#   fonts[MPtype == "Input"] <- 4
#   Cex <- 1.5
#   if (!ShowLabs) 
#     points(x, y, bg = coly, pch = Pch, cex = Cex, col = "black")
#   if (ShowLabs) 
#     text(x, y, MPs, font = fonts, col = coly, cex = 1)
#   
#   # Add Legend
#   text(min(XLim), max(YLim), "Performance Limits", cex = 1.25, pos = 4)
#   temp <- 0.95
#   for (XX in 1:length(Legend)) {
#     text(min(XLim), max(YLim) * temp, Legend[XX], pos = 4)
#     temp <- temp - 0.05
#   }
#   Cvec <- c("darkgray", "black", "Medium Sea Green", "green", "lightgray")
#   if (!ShowLabs) 
#     legend(100, max(YLim), pch = c(21, 21, 21, 21, 24), legend = c("Not Acceptable", 
#                                                                    "Acceptable", "Available", "Acceptable & Available", "Reference"), 
#            bty = "n", pt.bg = Cvec, pt.cex = Cex, xpd = NA)
#   
#   if (ShowLabs) 
#     legend(100, max(YLim), text.font = c(2, 4), legend = c("Output Control", 
#                                                            "Input Control", "Not Acceptable", "Acceptable", "Available", 
#                                                            "Acceptable & Available", "Reference"), bty = "n", text.col = c("black", 
#                                                                                                                            "black", Cvec), pt.cex = Cex, xpd = NA)
#   
#   Years <- paste("Years", (MSEobj@proyears - lastYrs) + 1, "-", MSEobj@proyears, 
#                  "(last", lastYrs, "years)")
#   mtext(side = 3, Years, cex = 1.25)
#   
#   if (YVar == "B_BMSYm") 
#     abline(h = 100, col = "#99999940", lwd = 2)
#   if (YVar == "SSB_SSB0m") 
#     abline(h = 50, col = "#99999940", lwd = 2)
#   par(op)
#   DF <- data.frame(MPs = MPs, Yield = y, Prob = x, Pass = x >= vl, stringsAsFactors = FALSE)
#   invisible(DF[order(DF$Prob, decreasing = TRUE), ])
# }
# 
#' KOBE plot: a projection by projection plot of F/FMSY and B/BMSY
#' 
#' A standard KOBE plot by each method that also shows the percentage of
#' methods that ended up in each quadrant.
#' 
#' 
#' @param MSEobj An object of class MSE
#' @param maxsim Maximum number of simulations (lines) to plot on each panel.
#' @param MPs Optional subset MSE object by MP
#' @param sims Optional subset MSE object by simulation
#' @param maxMP Maximum number of MPs to include in plot
#' @param nam The name of the plot
#' @param cex.leg Size of legend
#' @note Apologies for the nauseating shading.
#' @author T. Carruthers with some additions from A. Hordyk
#' @export Kplot
Kplot <- function(MSEobj, maxsim = 60, MPs = NA, sims = NULL, maxMP = 9, 
                  nam = NA, cex.leg = 1.5) {
  
  # png('Kplot.png')
  
  if (!is.null(sims) & all(is.na(MPs))) 
    MSEobj <- Sub(MSEobj, sims = sims)
  if (!is.null(sims) & all(!is.na(MPs))) 
    MSEobj <- Sub(MSEobj, sims = sims, MPs = MPs)
  if (is.null(sims) & !all(is.na(MPs))) 
    MSEobj <- Sub(MSEobj, MPs = MPs)
  
  nMPs <- MSEobj@nMPs
  nsim <- MSEobj@nsim
  if (is.null(sims) & nsim > maxsim) 
    MSEobj <- Sub(MSEobj, sims = 1:maxsim)
  
  if (nMPs > maxMP) {
    message("MSE object has more than ", maxMP, " MPs. Plotting the first ", 
            maxMP)
    MSEobj <- Sub(MSEobj, MPs = 1:maxMP)
    nMPs <- MSEobj@nMPs
  }
  
  nr <- floor((MSEobj@nMPs)^0.5)
  nc <- ceiling((MSEobj@nMPs)/nr)
  Cex <- 1.5
  TitleCex <- 1.5
  
  FMSYr <- quantile(MSEobj@F_FMSY, c(0.001, 0.9), na.rm = T)
  BMSYr <- quantile(MSEobj@B_BMSY, c(0.001, 0.975), na.rm = T)
  
  # dev.new2(width=nc*3,height=nr*3.6)
  # par(mfrow=c(nr,nc),mai=c(0.45,0.45,0.45,0.01),omi=c(0.45,0.3,0.35,0.01))
  # par(mfcol=c(nr,nc),mai=c(0.2,0.35,0.3,0.01),omi=c(0.5,0.4,0.4,0.05))
  if (is.na(nam)) 
    op <- par(mfrow = c(nr, nc), mar = c(2, 2, 3, 1), oma = c(3, 3.5, 1.2, 
                                                              0))
  if (!is.na(nam)) 
    op <- par(mfrow = c(nr, nc), mar = c(2, 2, 3, 1), oma = c(3, 3.5, 3, 0))
  colsse <- rainbow(MSEobj@proyears, start = 0.63, end = 0.95)[1:MSEobj@proyears]
  colsse <- makeTransparent(colsse, 95)
  
  XLim <- c(0, 3)
  YLim <- c(0, 2.5)
  pmat <- matrix(NA, nrow = nc, ncol = nr, byrow = FALSE)
  pmat[1:nMPs] <- 1:nMPs
  pmat <- t(pmat)
  
  for (mm in 1:MSEobj@nMPs) {
    plot(MSEobj@B_BMSY[1, mm, 1], MSEobj@F_FMSY[1, 
                                                mm, 1], xlim = XLim, ylim = YLim, 
         col = colsse[1],  bty = "n", axes = FALSE)
    
    if (nrow(pmat) > 1) {
      if (mm %in% pmat[, 1]) {
        axis(side = 2, labels = TRUE, las = 1)
      } else axis(side = 2, labels = FALSE)
      if (mm %in% pmat[nr, ]) {
        axis(side = 1, labels = TRUE)
      } else axis(side = 1, labels = FALSE)
      nas <- apply(pmat, 1, sum)
      rr <- which.max(nas)
      nas2 <- apply(pmat, 2, sum)
      cc <- which(is.na(nas2))
      if (mm %in% pmat[rr, cc]) 
        axis(side = 1, labels = TRUE)
    } else {
      if (mm == 1) 
        axis(side = 2, labels = TRUE, las = 1)
      if (mm != 1) 
        axis(side = 2, labels = FALSE)
      axis(side = 1, labels = TRUE)
    }
    
    
    OO <- round(sum(MSEobj@B_BMSY[, mm, MSEobj@proyears] < 1 & MSEobj@F_FMSY[, 
                                                                             mm, MSEobj@proyears] > 1, na.rm = T)/MSEobj@nsim * 100, 1)
    OU <- round(sum(MSEobj@B_BMSY[, mm, MSEobj@proyears] > 1 & MSEobj@F_FMSY[, 
                                                                             mm, MSEobj@proyears] > 1, na.rm = T)/MSEobj@nsim * 100, 1)
    UO <- round(sum(MSEobj@B_BMSY[, mm, MSEobj@proyears] < 1 & MSEobj@F_FMSY[, 
                                                                             mm, MSEobj@proyears] < 1, na.rm = T)/MSEobj@nsim * 100, 1)
    UU <- round(sum(MSEobj@B_BMSY[, mm, MSEobj@proyears] > 1 & MSEobj@F_FMSY[, 
                                                                             mm, MSEobj@proyears] < 1, na.rm = T)/MSEobj@nsim * 100, 1)
    
    # alp<-80
    # polygon(c(1,-1000,-1000,1),c(1,1,1000,1000),col=makeTransparent('orange',alp),border=makeTransparent('orange',alp))
    # polygon(c(1,1000,1000,1),c(1,1,1000,1000),col=makeTransparent('yellow',alp),border=makeTransparent('yellow',alp))
    # polygon(c(1,-1000,-1000,1),c(1,1,-1000,-1000),col=makeTransparent('yellow',alp),border=makeTransparent('yellow',alp))
    # polygon(c(1,1000,1000,1),c(1,1,-1000,-1000),col=makeTransparent('green',alp),border=makeTransparent('yellow',alp))
    
    
    abline(h = 1, col = "grey", lwd = 3)
    abline(v = 1, col = "grey", lwd = 3)
    # abline(v=c(0.1,0.5),col='grey',lwd=2)
    y <- 1:(MSEobj@proyears - 1)
    #y1 <- y + 1
    #x0 <- as.vector(MSEobj@B_BMSY[, mm, y])
    #x1 <- as.vector(MSEobj@B_BMSY[, mm, y1])
    #  y0 <- as.vector(MSEobj@F_FMSY[, mm, y])
    # y1 <- as.vector(MSEobj@F_FMSY[, mm, y1])
    #segments(x0, y0, x1, y1, col = rep(colsse,each=nsim))
    
    rng <- 1:min(maxsim, MSEobj@nsim)
    points(MSEobj@B_BMSY[rng, mm, 1], MSEobj@F_FMSY[rng, mm, 1], pch = 19, 
           cex = 0.8, col = colsse[1])
    points(MSEobj@B_BMSY[rng, mm, MSEobj@proyears], MSEobj@F_FMSY[rng, 
                                                                  mm, MSEobj@proyears], pch = 19, cex = 0.8, col = colsse[MSEobj@proyears])
    
    if (mm == 1) 
      legend("right", c("Start", "End"), bty = "n", text.col = c(colsse[1], 
                                                                 colsse[MSEobj@proyears]), pch = 19, col = c(colsse[1], 
                                                                                                             colsse[MSEobj@proyears]))
    legend("topleft", paste(OO, "%", sep = ""), bty = "n", text.font = 2, 
           cex = cex.leg)
    legend("topright", paste(OU, "%", sep = ""), bty = "n", text.font = 2, 
           cex = cex.leg)
    legend("bottomleft", paste(UO, "%", sep = ""), bty = "n", text.font = 2, 
           cex = cex.leg)
    legend("bottomright", paste(UU, "%", sep = ""), bty = "n", text.font = 2, 
           cex = cex.leg)
    
    mtext(MSEobj@MPs[mm], 3, line = 0.6, cex = TitleCex)
  }
  
  mtext(expression(B/B[MSY]), 1, outer = T, line = 2, cex = Cex)
  mtext(expression(F/F[MSY]), 2, outer = T, line = 1.2, cex = Cex)
  # if(is.na(nam))mtext(deparse(substitute(MSEobj)),3,outer=T,line=0.25,font=2,
  # cex=TitleCex)
  # if(!is.na(nam))mtext(MSEobj@Name,3,outer=T,line=0.25,font=2,
  # cex=TitleCex)
  if (!is.na(nam)) 
    mtext(nam, 3, outer = TRUE, line = 0.25, font = 2, cex = TitleCex)
  # dev.off()
  par(op)
}



#' National Oceanographic and Atmospheric Administration default plot 1
#' 
#' A preliminary plot for returning trade-offs plots and performance table for
#' total yield, variability in yield, probability of overfishing and likelihood
#' of biomass dropping below 50 per cent BMSY
#' 
#' 
#' @param MSEobj An object of class MSE
#' @param nam Title of plot
#' @param type Plots full range of data if NA. Plots a subset that meet
#' thresholds if not NA.
#' @param panel Should a two panel plot be made or should plots be made in
#' sequence.
#' @return A table of performance metrics.
#' @author T. Carruthers
#' @export NOAA_plot
NOAA_plot <- function(MSEobj, nam = NA, type = NA, panel = T) {
  Yd <- rep(NA, MSEobj@nMPs)
  B50 <- rep(NA, MSEobj@nMPs)
  PNOF <- rep(NA, MSEobj@nMPs)
  LTY <- rep(NA, MSEobj@nMPs)
  STY <- rep(NA, MSEobj@nMPs)
  VY <- rep(NA, MSEobj@nMPs)
  
  y1 <- 1:(MSEobj@proyears - 1)
  y2 <- 2:MSEobj@proyears
  
  yend <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
  
  RefYd <- MSEobj@OM$RefY
  
  for (mm in 1:MSEobj@nMPs) {
    
    PNOF[mm] <- round(mean(MSEobj@F_FMSY[, mm, ] <= 1, na.rm = T) * 100, 1)
    B50[mm] <- round(mean(MSEobj@B_BMSY[, mm, ] >= 0.5, na.rm = T) * 100, 1)
    LTY[mm] <- round(mean(MSEobj@C[, mm, yend]/RefYd >= 0.5, na.rm = T), 3) * 100
    AAVY <- apply((((MSEobj@C[, mm, y1] - MSEobj@C[, mm, y2])/MSEobj@C[, mm, y2])^2)^0.5, 1, mean, na.rm = T)
    VY[mm] <- round(mean(AAVY <= 0.15, na.rm = T), 3) * 100
    
  }
  
  # dev.new2(width=7,height=7)
  if (panel) 
    op <- par(mfrow = c(1, 2), mai = c(1.5, 1.5, 0.1, 0.1), omi = c(0.1, 0.1, 0.4, 0))
  
  if (is.na(type)) {
    tradeoffplot(PNOF, LTY, "Prob. of not overfishing (%)", "Long-term yield ", 
                 MSEobj@MPs[1:MSEobj@nMPs], vl = 50, hl = 100)
    tradeoffplot(B50, VY, "Prob. biomass above half BMSY (%)", "Prob. AAVY less than 15%", 
                 MSEobj@MPs[1:MSEobj@nMPs], vl = 80, hl = 50)
  } else {
    tradeoffplot3(PNOF, LTY, "Prob. of not overfishing (%)", "Long-term yield", 
                  MSEobj@MPs[1:MSEobj@nMPs], vl = 50, hl = 100, xlim = c(45, 
                                                                         105), ylim = c(0, 105))
    tradeoffplot3(B50, VY, "Prob. biomass above half BMSY (%)", "Prob. AAVY less than 15%", 
                  MSEobj@MPs[1:MSEobj@nMPs], vl = 80, hl = 50, xlim = c(75, 105), 
                  ylim = c(45, 105))
    
  }
  
  # if(is.na(nam))mtext(deparse(substitute(MSEobj)),3,outer=T,line=0.3,font=2)
  # if(!is.na(nam) &
  # !is.character(nam))mtext(MSEobj@Name,3,outer=T,line=0.3,font=2)
  # if(!is.na(nam) &
  # is.character(nam))mtext(nam,3,outer=T,line=0.3,font=2)
  
  if (panel) par(op)
  temp <- data.frame(PNOF, B50, LTY, VY)
  row.names(temp) <- MSEobj@MPs[1:MSEobj@nMPs]
  temp
  
}



#' A projection by projection plot of F/FMSY and B/BMSY
#' 
#' A shorter version of the plot method for MSEs that just shows the projected
#' trends in stock status and over exploitation
#' 
#' 
#' @param MSEobj An object of class MSE
#' @param nam Title of plot
#' @param maxMP The maximum number of MPs to plot (defaults to the first 10)
#' @param MPs A character vector of MPs to plot
#' @param maxsims Integer, the maximum number of simulations to plot
#' @author T. Carruthers 
#' @export Pplot
Pplot <- function(MSEobj, nam = NA, maxMP = 10,MPs=NA,maxsims=20) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  if(!is.na(MPs)){
    maxMP<-length(MPs)
    MSEobj<-Sub(MSEobj,MPs=MPs)
  } else{
    
    if(MSEobj@nMPs>maxMP)MSEobj<-Sub(MSEobj,MPs=MSEobj@MPs[1:maxMP])
    
  }
  maxsims <- min(maxsims, MSEobj@nsim)
  MSEobj<-Sub(MSEobj,sims=1:maxsims)
  
  FMSYr <- quantile(MSEobj@F_FMSY, c(0.001, 0.9), na.rm = T)
  BMSYr <- quantile(MSEobj@B_BMSY, c(0.001, 0.975), na.rm = T)
  
  colsse <- rainbow(100, start = 0, end = 0.36)[1:100]
  colB <- rep(colsse[100], ceiling(BMSYr[2] * 100))
  colB[1:100] <- colsse
  colB <- makeTransparent(colB, 60)
  colsse <- rainbow(200, start = 0, end = 0.36)[200:1]
  colF <- rep(colsse[200], ceiling(FMSYr[2] * 100))
  colF[1:200] <- colsse
  colF <- makeTransparent(colF, 60)
  
  Yd <- rep(NA, MSEobj@nMPs)
  P10 <- rep(NA, MSEobj@nMPs)
  P50 <- rep(NA, MSEobj@nMPs)
  P100 <- rep(NA, MSEobj@nMPs)
  POF <- rep(NA, MSEobj@nMPs)
  yind <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
  RefYd <- MSEobj@OM$RefY
  
  for (mm in 1:MSEobj@nMPs) {
    Yd[mm] <- round(mean(apply(MSEobj@C[, mm, yind], 1, mean, na.rm = T)/RefYd, 
                         na.rm = T) * 100, 1)
    # cbind(MSEobj@C[,mm,yind],unlist(MSEobj@OM$MSY))
    POF[mm] <- round(sum(MSEobj@F_FMSY[, mm, ] > 1, na.rm = T)/prod(dim(MSEobj@F_FMSY[, 
                                                                                      mm, ]), na.rm = T) * 100, 1)
    P10[mm] <- round(sum(MSEobj@B_BMSY[, mm, ] <= 0.1, na.rm = T)/prod(dim(MSEobj@B_BMSY[, 
                                                                                         mm, ])) * 100, 1)
    P50[mm] <- round(sum(MSEobj@B_BMSY[, mm, ] <= 0.5, na.rm = T)/prod(dim(MSEobj@B_BMSY[, 
                                                                                         mm, ])) * 100, 1)
    P100[mm] <- round(sum(MSEobj@B_BMSY[, mm, ] <= 1, na.rm = T)/prod(dim(MSEobj@B_BMSY[, 
                                                                                        mm, ])) * 100, 1)
  }
  
  nr <- ceiling(MSEobj@nMPs/8)
  nc <- ceiling(MSEobj@nMPs/nr)
  nr <- nr * 2
  MSEcols <- c("red", "green", "blue", "orange", "brown", "purple", "dark grey", 
               "violet", "dark red", "pink", "dark blue", "grey")
  temp <- array(0, c(nr * 2 + (nr/2 - 1), nc * 2))
  i <- 0
  for (c in 1:nc) {
    for (r in 1:nr) {
      i <- i + 1
      temp[(ceiling(r/2) - 1) + (1:2) + (r - 1) * 2, (1:2) + (c - 1) * 2] <- ((c - 1) * nr) + r
    }
  }
  
  par(mfcol = c(nr, nc), mar = c(2, 2, 2, 1), oma = c(3, 2, 2, 0))
  layout(temp)
  # dev.new2(width=nc*3,height=nr*3)
  lwdy <- 2.5
  
  for (mm in 1:MSEobj@nMPs) {
    plot(MSEobj@F_FMSY[1, mm, ], ylim = FMSYr, col = colF[ceiling(mean(MSEobj@F_FMSY[1, 
                                                                                     mm, ], na.rm = T) * 100)], type = "l", lwd = lwdy)
    for (i in 1:MSEobj@nsim) lines(MSEobj@F_FMSY[i, mm, ], col = colF[ceiling(mean(MSEobj@F_FMSY[i, 
                                                                                                 mm, ], na.rm = T) * 100)], lwd = lwdy)
    abline(h = 100, col = "grey", lwd = 3)
    mtext(MSEobj@MPs[mm], 3, outer = F, line = 0.6)
    legend("topright", c(paste(POF[mm], "% POF", sep = ""), paste(Yd[mm], 
                                                                  "% FMSY yield", sep = "")), bty = "n", cex = 0.8)
    if (mm %in% (1:(nr/2))) 
      mtext("F/FMSY", 2, line = 2.5, outer = F)
    abline(h = 1, col = makeTransparent("grey", 30), lwd = 2.5)
    
    plot(MSEobj@B_BMSY[1, mm, ], ylim = BMSYr, col = colB[ceiling(MSEobj@B_BMSY[1, 
                                                                                mm, MSEobj@proyears] * 100)], type = "l", lwd = lwdy)
    for (i in 1:MSEobj@nsim) lines(MSEobj@B_BMSY[i, mm, ], col = colB[ceiling(MSEobj@B_BMSY[i, 
                                                                                            mm, MSEobj@proyears] * 100)], lwd = lwdy)
    abline(h = 100, col = "grey", lwd = 3)
    legend("topright", c(paste(P100[mm], "% < BMSY", sep = ""), paste(P50[mm], 
                                                                      "% < 0.5BMSY", sep = ""), paste(P10[mm], "% < 0.1BMSY", sep = "")), 
           bty = "n", cex = 0.8)
    if (mm %in% (1:(nr/2))) 
      mtext("B/BMSY", 2, line = 2.5, outer = F)
    abline(h = 1, col = makeTransparent("grey", 30), lwd = 2.5)
    
  }
  mtext("Projection year", 1, outer = T, line = 1.2)
  if (is.na(nam)) 
    mtext(deparse(quote(MSEobj)), 3, outer = T, line = 0.3, font = 2)
  if (!is.na(nam) & !is.character(nam)) 
    mtext(MSEobj@Name, 3, outer = T, line = 0.3, font = 2)
  if (!is.na(nam) & is.character(nam)) 
    mtext(nam, 3, outer = T, line = 0.3, font = 2)
  return(invisible())
}


#' A projection by projection plot of F/FMSY, B/BMSY, B/B0, and yield
#' 
#' 
#' @param MSEobj An object of class MSE
#' @param YVar What to plot on the y-axis? Options are: \code{c('SSB_SSB0',
#' 'SSB_SSBMSY', 'F_FMSY', 'Yield')}
#' @param MPs Optional subset by MP
#' @param sims Optional subset by simulation
#' @param traj Plot all projections (\code{all}), only quantiles
#' (\code{quant}), or both projections and median (\code{both})
#' @param quants Numeric vector of length 2 specifying the quantiles (e.g.,
#' 10th and 90th. Median is always included)
#' @param incquant Logical. Include the quantiles or only plot median?
#' @param quantcol Colour of the quantile polygon
#' @param ref.lines Numeric vector of y-values for horizontal reference lines. Set to NULL to remove lines.
#' @param RefYield Should yield be relative to long-term optimum (\code{lto})
#' or last historical year (\code{curr})
#' @param LastYr Logical. Include the last historical year in the yield
#' projections?
#' @param maxMP Maximum number of MPs to plot
#' @param alpha Alpha for transparency of lines
#' @param cex.axis Size of axis text
#' @param cex.lab Size of axis label
#' @param YLab Optional label for y-axis
#' @param incMP Logical. Include name of MP?
#' @param MPcex Size of MP label
#' @param MPcol Optional character vector of colors for MP labels
#' @param incLeg Logical. Include a legend?
#' @param cex.leg Size of legend text
#' @param legPos Legend position
#' @param yline Optional horizontal lines
#' @param xline Optional vertical lines
#' @param parOR Logical to over-ride the par parameters
#' @param xaxis Logical. Should x-axis labels be displayed?
#' @param yaxis Logical. Should y-axis labels be displayed?
#' @param oneIt Logical. Should one iteration be plotted on the quantile plot?
#' @param ...  Additional arguments to be passed to plotting functions
#' @author T. Carruthers & A.Hordyk
#' @export Pplot2
Pplot2 <- function(MSEobj, YVar = c("F_FMSY", "SSB_SSBMSY"), MPs = NA, sims = NULL, 
                   traj = c("all", "quant", "both"), quants = c(0.1, 0.9), incquant = TRUE, 
                   quantcol = "lightgray", 
                   RefYield = c("lto", "curr"), LastYr = TRUE, 
                   ref.lines=c(0.5, 1, 1.5), maxMP = 6, alpha = 60, 
                   cex.axis = 1, cex.lab = 1, YLab = NULL, incMP = TRUE, MPcex = 1, 
                   MPcol='black',
                   incLeg = TRUE, cex.leg = 1.5, legPos = "topleft", yline = NULL, xline=NULL, parOR = FALSE, 
                   xaxis = TRUE, yaxis = TRUE, oneIt=TRUE, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  YVars <- c("SSB_SSB0", "SSB_SSBMSY", "F_FMSY", "Yield")
  YVar <- match.arg(YVar, choices = YVars, several.ok = TRUE)
  op <- par(no.readonly=TRUE)
  if (!is.null(YLab) & length(YLab) != length(YVar)) 
    stop("Length of YLab must equal length of YVar")
  if (!is.null(sims) & all(is.na(MPs))) 
    MSEobj <- Sub(MSEobj, sims = sims)
  if (!is.null(sims) & all(!is.na(MPs))) 
    MSEobj <- Sub(MSEobj, sims = sims, MPs = MPs)
  if (is.null(sims) & !all(is.na(MPs))) 
    MSEobj <- Sub(MSEobj, MPs = MPs)
  
  nsim <- MSEobj@nsim
  if (nsim < 10) 
    alpha <- 180
  nMPs <- MSEobj@nMPs
  if (nMPs > maxMP) {
    message("MSE object has more than ", maxMP, " MPs. Plotting the first ", 
            maxMP)
    MSEobj <- Sub(MSEobj, MPs = 1:maxMP)
    nMPs <- MSEobj@nMPs
  }
  MPs <- MSEobj@MPs
  proyears <- MSEobj@proyears
  RefYd <- MSEobj@OM$RefY
  RefYield <- match.arg(RefYield)
  traj <- match.arg(traj)
  # Calculate Statistics Biomass/B0
  temp <- as.matrix(expand.grid(1:nsim, 1:nMPs, 1:proyears))
  Deplet <- array(NA, dim = dim(MSEobj@B_BMSY))
  Deplet[temp] <- (MSEobj@B_BMSY[temp] * MSEobj@OM$SSBMSY_SSB0[temp[, 1]])
  
  # Yield - including last historical year (current year)
  pastC <- apply(MSEobj@CB_hist[, , , , drop = FALSE], c(1, 3), sum, 
                 na.rm = TRUE)/RefYd
  temp <- aperm(replicate(nMPs, pastC), c(1, 3, 2))
  lastYr <- temp[, , MSEobj@nyears, drop = FALSE]
  Yield <- abind::abind(lastYr, MSEobj@C[, , , drop = FALSE]/RefYd, along = 3)
  
  Dat <- list(SSB_SSB0 = Deplet, SSB_SSBMSY = MSEobj@B_BMSY, F_FMSY = MSEobj@F_FMSY, 
              Yield = Yield)
  Dat <- Dat[YVar]
  
  if ("Yield" %in% YVar & RefYield == "curr") {
    ny <- dim(Dat$Yield)[3]
    Dat$Yield <- Dat$Yield[, , , drop = FALSE]/Dat$Yield[, , rep(1, 
                                                                 ny), drop = FALSE]
  }
  if ("Yield" %in% YVar & !LastYr) {
    Dat$Yield <- Dat$Yield[, , 2:proyears, drop = FALSE]
  }
  
  nr <- length(Dat)
  nc <- nMPs
  
  dots <- list(...)
  
  ylims <- cbind(0, unlist(lapply(Dat, quantile, 0.9, na.rm = TRUE)))
  if ("SSB_SSB0" %in% YVar) {
    index <- which(YVar == "SSB_SSB0")
    ylims[index, ] <- c(0, max(1, max(ylims[index, ])))
  }
  if (length(dots$ylim) != 0) ylims <- matrix(rep((dots$ylim), length(Dat)), nrow = nr, byrow = TRUE)
  colrange <- matrix(unlist(lapply(Dat, quantile, c(0.001, 0.975), na.rm = TRUE)), 
                     nrow = nr, byrow = TRUE)
  
  colsse <- rainbow(100, start = 0, end = 0.36)[1:100]
  # Col<-rep(colsse[100],ceiling(colrange[2]*100)) Col[1:100]<-colsse
  Col <- makeTransparent(colsse, alpha)
  
  if (length(dots$lwd) == 0) lwd <- 3
  if (length(dots$lwd) != 0) lwd <- dots$lwd
  
  YLabs <- list(expression(SSB/SSB[0]), expression(SSB/SSB[MSY]), expression(F/F[MSY]), 
                "Yield relative\n to Long-Term\n Optimum")
  if ("Yield" %in% YVar & RefYield == "curr") 
    YLabs[[4]] <- expression(Yield/Yield[current])
  YLabs <- YLabs[match(YVar, YVars)]
  if (!is.null(YLab)) 
    YLabs <- YLab
  if (!parOR) {
    if ("Yield" %in% YVar & RefYield != "curr") {
      op <- par(mfrow = c(nr, nc), bty = "n", mar = c(2, 2, 0, 0), oma = c(4, 8, 2, 1))
    } else op <- par(mfrow = c(nr, nc), bty = "n", mar = c(2, 2, 0, 0), oma = c(4, 4, 2, 1))
  }
  if (parOR) {
    nr <- par()$mfrow[1]
    nc <- par()$mfrow[2]
  }
  
  for (X in 1:length(Dat)) {
    Col2 <- Col
    dat <- Dat[[X]]
    ylim <- ylims[X, ]
    ylab <- YLabs[[X]]
    if (grepl("F_FMSY", YVar[X])) Col2 <- rev(Col)
    for (mm in 1:nMPs) {
      plot(1:length(dat[1, mm, ]), dat[1, mm, ], ylim = ylim, type = "n", 
           axes = FALSE, xlab = "", ylab = "")
      
      # add reference lines 
      if (!is.null(ref.lines)) {
        for (h in ref.lines) abline(h=h, lty=3, lwd=1, col="darkgray")
      }
      

      if (traj == "all" | traj=="both") 
        for (i in 1:MSEobj@nsim) lines(dat[i, mm, ], 
                                       col = Col2[min(100, ceiling(dat[i, mm, length(dat[i, mm, ])] * 100))], lwd = lwd)
      if (traj == "quant" | traj=="both") {
        stats <- apply(dat[, mm, , drop = FALSE], 3, quantile, 
                       c(quants[1], 0.5, quants[2]), na.rm = TRUE)
        if (length(quants) == 4) 
          stats2 <- apply(dat[, mm, , drop = FALSE], 3, quantile, 
                          c(quants[3], quants[4]), na.rm = TRUE)
        if (length(quants) != 4) stats2 <- NULL
        if (traj=="both") {
          stats2 <- NULL
          incquant <- FALSE
        }
          
        if (!incquant) lines(1:length(stats[2, ]), stats[2, ], lwd = 3)
        if (incquant) {
          if (!is.null(stats2)) {
            if (length(quantcol) < 2) 
              quantcol <- c(quantcol, "darkgray")
            polygon(x = c(1:length(stats2[2, ]), length(stats2[2, ]):1), 
                    y = c(stats2[1, ], rev(stats2[2, ])), col = quantcol[2], 
                    border = FALSE)
          }
          polygon(x = c(1:length(stats[2, ]), length(stats[2, ]):1), 
                  y = c(stats[1, ], rev(stats[3, ])), col = quantcol[1], 
                  border = FALSE)
          lines(1:length(stats[2, ]), stats[2, ], lwd = 3)
          if (oneIt) lines(dat[2, mm, , drop = FALSE] , lwd=1)
        }
      }
      if (X == nr & xaxis) 
        axis(side = 1, labels = TRUE, cex.axis = cex.axis)
      if (X != nr | !xaxis) 
        axis(side = 1, labels = FALSE)
      if (mm == 1 & yaxis) {
        axis(side = 2, labels = TRUE, cex.axis = cex.axis, las = 1)
        mtext(side = 2, ylab, cex = cex.lab, line = 3)
      }
      if (parOR) {
        if (xaxis) 
          axis(side = 1, labels = TRUE, cex.axis = cex.axis, las = 1)
        if (yaxis & mm == 1) 
          axis(side = 2, labels = TRUE, cex.axis = cex.axis, las = 1)
      }
      axis(side = 2, labels = FALSE)
      
      MPcol <- rep(MPcol, MSEobj@nMPs)[1:MSEobj@nMPs]
      if (incMP & X == 1 & !parOR) 
        mtext(side = 3, MSEobj@MPs[mm], cex = MPcex, col=MPcol[mm])
      if (incMP & parOR) 
        mtext(side = 3, MSEobj@MPs[mm], cex = MPcex, col=MPcol[mm])
      
      # Legend #
      if (mm == 1 & incLeg & (traj == "quant"||traj=="both") & X == 1) {
        if (incquant) {
          if (is.null(stats2)) {
            pquants <- quants
            if (max(quants) < 1) pquants <- quants * 100
            legend(legPos, legend = c(expression("50"^"th" * " (median)"), 
                                      bquote(.(pquants[1])^th ~ and ~ .(pquants[2])^th)), 
                   pch = 22, title = "Percentile", pt.bg = c("black", quantcol[1],  quantcol[2]), 
                   bty = "n", cex = cex.leg, xpd = NA, col = "black")
            
          }
          if (!is.null(stats2)) {
            pquants <- quants
            if (max(quants) < 1) 
              pquants <- quants * 100
            legend(legPos, legend = c(expression("50"^"th" * " (median)"), 
                                      bquote(.(pquants[1])^th ~ and ~ .(pquants[2])^th), 
                                      bquote(.(pquants[3])^th ~ and ~ .(pquants[4])^th)), 
                   pch = 22, title = "Percentile", pt.bg = c("black", quantcol[1], quantcol[2]), 
                   bty = "n", cex = cex.leg, xpd = NA, col = "black")
          }
        } else {
          legend(legPos, legend = "Median", pch = 15, col = "black", bty = "n", cex = 1.25)
        }
      }
      if (!is.null(yline)) 
        abline(h = yline[X], lwd = 2, col = "darkgray", lty = 2)
      if (!is.null(xline)) {
        for (xx in 1:length(xline)) {
          abline(v = xline[xx], lwd = 2, col = "darkgray", lty = xx+1)
        }
      }
    }
  }
  mtext(side = 1, "Projection Years", line = 2, cex = cex.lab, outer = TRUE)
  if (!parOR) par(op)
  invisible(Dat)
}



#' Performance Whisker Plot
#' 
#' A NAFO / ICCAT / SSB style MSE performance whisker plot
#' 
#' 
#' @param MSEobj An object of class MSE
#' @return A box plot of performance 
#' @author T. Carruthers
#' @export PWhisker
PWhisker<-function(MSEobj){#},Pnames=c("POF","C30","D30","LD","DNC","LDNC","PGK","AAVC")){
  
  P10 <- P50 <- PNOF <- LTY <- STY <- AAVY <- VY <- array(NA,c(MSEobj@nsim, MSEobj@nMPs))
  Pnames<-c("B10","B50","PNOF","LTY","AAVY")
  
  nsim<-MSEobj@nsim
  nMPs<-MSEobj@nMPs
  nperf<-length(Pnames)
  
  store<-array(NA,c(MSEobj@nMPs,nperf,5))
  
  y1 <- 1:(MSEobj@proyears - 1)
  y2 <- 2:MSEobj@proyears
  
  yend <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
  
  RefYd <- MSEobj@OM$RefY
  
  for (mm in 1:MSEobj@nMPs) {
    
    P10[,mm] <- round(apply(MSEobj@B_BMSY[, mm, ] <= 0.1,1,mean, na.rm = T) * 100, 1)
    P50[,mm] <- round(apply(MSEobj@B_BMSY[, mm, ] <= 0.5,1,mean, na.rm = T) * 100, 1)
    PNOF[,mm] <- round(apply(MSEobj@F_FMSY[, mm, ] <= 1,1,mean, na.rm = T) * 100, 1)
    LTY[,mm] <- round(apply(MSEobj@C[, mm, yend]/RefYd >= 0.5,1,mean, na.rm = T), 3) * 100
    AAVY[,mm] <- apply((((MSEobj@C[, mm, y1] - MSEobj@C[, mm, y2])/MSEobj@C[, mm, y2])^2)^0.5, 1, mean, na.rm = T)
    #VY[,mm] <- round(apply(AAVY <= 0.15, na.rm = T), 3) * 100
    
    store[mm,1,]<-quantile(P10[,mm],c(0.05,0.25,0.5,0.72,0.95))
    store[mm,2,]<-quantile(P50[,mm],c(0.05,0.25,0.5,0.72,0.95))
    store[mm,3,]<-quantile(PNOF[,mm],c(0.05,0.25,0.5,0.72,0.95))
    store[mm,4,]<-quantile(LTY[,mm],c(0.05,0.25,0.5,0.72,0.95))
    store[mm,5,]<-quantile(AAVY[,mm],c(0.05,0.25,0.5,0.72,0.95))
    #store[mm,6,]<-quantile(VY[,mm],c(0.05,0.25,0.5,0.72,0.95))
    
  }
  
  par(mfrow=c(nperf,1),mai=c(0.1,0.3,0.01,0.05),omi=c(1.2,0.4,0.05,0.01))
  
  for(i in 1:nperf){
    
    xlab=F
    if(i==nperf)xlab=T
    custombar(dat=store[,i,],MPnams=MSEobj@MPs,xlab=xlab)
    mtext(Pnames[i],2,line=3,font=2)
    
  }
  
  mtext("Management Procedure",1,line=6.8,font=2,outer=T)
  
}


# #' Scatter plot of B/BMSY or B/B0 and F/FMSY for lastYrs
# #' 
# #' Scatter plot of B/BMSY or B/B0 and F/FMSY for lastYrs
# #' 
# #' 
# #' @param MSEobj An object of class MSE
# #' @param MPs Optional subset by MP
# #' @param All Logical. Plot all points or just the mean?
# #' @param Var What to plot on the y-axis: \code{B_BMSY} or \code{SSB_SSB0}
# #' @param lastYrs Last number of years in projection to calculate statistics
# #' @param Fref Location of F statistic reference line
# #' @param BMSYref Location of B_MSY statistic reference line
# #' @param B0ref Location of B_0 statistic reference line
# #' @param cex.MP size of MP label
# #' @param Fbg Logical. Include background colours for F-statistic?
# #' @param Bbg Logical. Include background colours for B-statistic?
# #' @param Props Logical. Display the proportion of points in each quadrant?
# #' @param TP Logical. Use transparent colours?
# #' @author A. Hordyk
# #' @export Splot
# Splot <- function(MSEobj = NULL, MPs = NA, All = TRUE, Var = c("B_BMSY", 
#                                                                "SSB_SSB0"), lastYrs = 10, Fref = 1, BMSYref = 1, B0ref = 0.4, cex.MP = 1, 
#                   Fbg = FALSE, Bbg = FALSE, Props = FALSE, TP = FALSE) {
#   
#   Var <- match.arg(Var)
#   if (!any(is.na(MPs))) 
#     MSEobj <- Sub(MSEobj, MPs = MPs)
#   
#   nMPs <- MSEobj@nMPs
#   nyrs <- MSEobj@proyears
#   nsim <- MSEobj@nsim
#   
#   perf <- MPStats(MSEobj)
#   if (lastYrs > MSEobj@proyears) 
#     lastYrs <- 10
#   yrs <- (nyrs - lastYrs + 1):nyrs
#   
#   ord <- perf[[1]]
#   MPord <- match(ord$MP, MSEobj@MPs)
#   sims <- perf$BySim
#   Nrow <- ceiling(sqrt(nMPs))
#   Ncol <- ceiling(nMPs/Nrow)
#   Yvar <- sims[["F_FMSY"]]
#   Xvar <- sims[[Var]]
#   YLim <- c(0, 2)
#   XMax <- ceiling(min(max(apply(Xvar[, , yrs], c(1, 2), mean), na.rm = TRUE), 
#                       3)/0.5) * 0.5
#   XLim <- c(0, XMax)
#   YLab <- expression(F/F[MSY])
#   XLab <- switch(Var, SSB_SSB0 = expression(B/B[0]), B_BMSY = expression(B/B[MSY]))
#   
#   Bref <- switch(Var, SSB_SSB0 = B0ref, B_BMSY = BMSYref)
#   
#   ColVec <- colorRampPalette(c("green", "orange", "red"))(nsim)
#   op <- par(mfrow = c(Nrow, Ncol), mar = c(1, 1, 2, 0), oma = c(4, 5, 2, 1))
#   
#   pmat <- matrix(NA, nrow = Ncol, ncol = Nrow)
#   pmat[1:nMPs] <- 1:nMPs
#   pmat <- t(pmat)
#   
#   
#   for (mm in seq_along(MPord)) {
#     xx <- MPord[mm]
#     if (!All) 
#       Xs <- apply(Xvar[, xx, yrs], 1, mean)
#     if (!All) 
#       Ys <- apply(Yvar[, xx, yrs], 1, mean)
#     if (All) 
#       Xs <- as.numeric(Xvar[, xx, yrs])
#     if (All) 
#       Ys <- as.numeric(Yvar[, xx, yrs])
#     colN <- ceiling((Xs/max(XLim)) * nsim) + 1
#     colN[colN > nsim] <- nsim
#     Cols <- ColVec[colN]
#     if (TP) 
#       Cols <- makeTransparent(Cols, alpha = 95)
#     plot(Xs, Ys, xlab = "", ylab = "", bty = "n", xlim = XLim, ylim = YLim, 
#          axes = FALSE, pch = 18, col = Cols)
#     
#     if (nrow(pmat) > 1) {
#       if (mm %in% pmat[, 1]) {
#         axis(side = 2, labels = TRUE, las = 1)
#       } else axis(side = 2, labels = FALSE)
#       if (mm %in% pmat[Nrow, ]) {
#         axis(side = 1, labels = TRUE)
#       } else axis(side = 1, labels = FALSE)
#       nas <- apply(pmat, 1, sum)
#       rr <- which.max(nas)
#       nas2 <- apply(pmat, 2, sum)
#       cc <- which(is.na(nas2))
#       if (mm %in% pmat[rr, cc]) 
#         axis(side = 1, labels = TRUE)
#     } else {
#       if (mm == 1) 
#         axis(side = 2, labels = TRUE, las = 1)
#       if (mm != 1) 
#         axis(side = 2, labels = FALSE)
#       axis(side = 1, labels = TRUE)
#     }
#     
#     abline(v = Fref, lty = 1, col = makeTransparent("grey", 150))
#     abline(h = Bref, lty = 1, col = makeTransparent("grey", 150))
#     mtext(side = 3, MSEobj@MPs[xx], cex = cex.MP)
#     
#     # Background colours
#     Alpha <- 30
#     # polygons
#     OutCol <- rgb(red = 255, green = 0, blue = 0, alpha = Alpha, names = NULL, 
#                   maxColorValue = 255)
#     vl <- Fref
#     hl <- Bref
#     if (Bbg) 
#       polygon(x = c(0, max(max(Xs), XLim[2]), max(max(Xs), XLim[2]), 
#                     0), y = c(0, 0, hl, hl), col = OutCol, border = NA)
#     if (Fbg) 
#       polygon(x = c(vl, max(max(Xs), XLim[2]), max(max(Xs), XLim[2]), 
#                     vl), y = c(0, 0, max(max(Ys), YLim[2]), max(max(Ys), YLim[2])), 
#               col = OutCol, border = NA)
#     
#     if (Props) {
#       P1 <- paste0(round(sum(Xs <= Fref & Ys >= Bref)/length(Xs) * 
#                            100, 1), "%")
#       P2 <- paste0(round(sum(Xs > Fref & Ys >= Bref)/length(Xs) * 
#                            100, 1), "%")
#       P3 <- paste0(round(sum(Xs > Fref & Ys < Bref)/length(Xs) * 
#                            100, 1), "%")
#       P4 <- paste0(round(sum(Xs <= Fref & Ys < Bref)/length(Xs) * 
#                            100, 1), "%")
#       legend("topleft", P1, bty = "n", cex = 1.25)
#       legend("topright", P2, bty = "n", cex = 1.25)
#       legend("bottomright", P3, bty = "n", cex = 1.25)
#       legend("bottomleft", P4, bty = "n", cex = 1.25)
#     }
#   }
#   mtext(side = 1, outer = TRUE, XLab, cex = 1.25, line = 2.5)
#   mtext(side = 2, outer = TRUE, YLab, cex = 1.25, line = 2)
#   Years <- paste("Years", (MSEobj@proyears - lastYrs) + 1, "-", MSEobj@proyears, 
#                  "(last", lastYrs, "years)")
#   mtext(side = 3, Years, cex = 1.25, outer = TRUE)
#   par(op)
# }
# 

#' Trade-off plots for an MSE object
#' 
#' Three figures showing trade-offs between fishing mortality, biomass, and yield.
#' 
#' @param MSEobj An object of class 'MSE'
#' @param nam Name of the plot
#' @author T. Carruthers & A. Hordyk
#' @seealso \link{TradePlot} \link{PerformanceMetric}
#' @describeIn Tplot Used in the plot method for MSE objects that shows trade-off between
#' yield versus probability of overfishing and biomass levels (relative to BMSY).
#' @export
Tplot_old <- function(MSEobj, nam = NA) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  FMSYr <- quantile(MSEobj@F_FMSY, c(0.001, 0.9), na.rm = T)
  BMSYr <- quantile(MSEobj@B_BMSY, c(0.001, 0.975), na.rm = T)
  
  colsse <- rainbow(100, start = 0, end = 0.36)[1:100]
  colB <- rep(colsse[100], ceiling(BMSYr[2] * 100))
  colB[1:100] <- colsse
  colB <- makeTransparent(colB, 60)
  colsse <- rainbow(200, start = 0, end = 0.36)[200:1]
  colF <- rep(colsse[200], ceiling(FMSYr[2] * 100))
  colF[1:200] <- colsse
  colF <- makeTransparent(colF, 60)
  
  Yd <- rep(NA, MSEobj@nMPs)
  P10 <- rep(NA, MSEobj@nMPs)
  P50 <- rep(NA, MSEobj@nMPs)
  P100 <- rep(NA, MSEobj@nMPs)
  POF <- rep(NA, MSEobj@nMPs)
  yind <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
  RefYd <- MSEobj@OM$RefY
  
  for (mm in 1:MSEobj@nMPs) {
    Yd[mm] <- round(mean(apply(MSEobj@C[, mm, yind], 1, mean, na.rm = T)/RefYd, 
                         na.rm = T) * 100, 1)
    # cbind(MSEobj@C[,mm,yind],unlist(MSEobj@OM$MSY))
    POF[mm] <- round(sum(MSEobj@F_FMSY[, mm, ] >= 1, na.rm = T)/prod(dim(MSEobj@F_FMSY[, 
                                                                                       mm, ]), na.rm = T) * 100, 1)
    P10[mm] <- round(sum(MSEobj@B_BMSY[, mm, ] <= 0.1, na.rm = T)/prod(dim(MSEobj@B_BMSY[, 
                                                                                         mm, ])) * 100, 1)
    P50[mm] <- round(sum(MSEobj@B_BMSY[, mm, ] <= 0.5, na.rm = T)/prod(dim(MSEobj@B_BMSY[, 
                                                                                         mm, ])) * 100, 1)
    P100[mm] <- round(sum(MSEobj@B_BMSY[, mm, ] <= 1, na.rm = T)/prod(dim(MSEobj@B_BMSY[, 
                                                                                        mm, ])) * 100, 1)
  }
  
  # dev.new2(width=7,height=7)
  par(mfrow = c(2, 2), mar = c(5, 4, 1, 1), oma = c(0, 0, 2, 0))
  
  tradeoffplot(POF, Yd, "Prob. of overfishing (%)", "Relative yield", 
               MSEobj@MPs[1:MSEobj@nMPs], vl = 50, hl = 100)
  tradeoffplot(P100, Yd, "Prob. biomass < BMSY (%)", "Relative yield", 
               MSEobj@MPs[1:MSEobj@nMPs], vl = 50, hl = 100)
  tradeoffplot(P50, Yd, "Prob. biomass < 0.5BMSY (%)", "Relative yield", 
               MSEobj@MPs[1:MSEobj@nMPs], vl = 50, hl = 100)
  tradeoffplot(P10, Yd, "Prob. biomass < 0.1BMSY (%)", "Relative yield", 
               MSEobj@MPs[1:MSEobj@nMPs], vl = 50, hl = 100)
  
  if (is.na(nam)) 
    mtext(deparse(substitute(MSEobj)), 3, outer = T, line = 0.3, font = 2)
  if (!is.na(nam) & !is.character(nam)) 
    mtext(MSEobj@Name, 3, outer = T, line = 0.3, font = 2)
  if (!is.na(nam) & is.character(nam)) 
    mtext(nam, 3, outer = T, line = 0.3, font = 2)
  return(invisible())
}



#' @describeIn Tplot Simpler plot that compares long-term yield (LTY:
#' fraction of simulations getting over half FMSY yield in the last ten years
#' of the projection), short-term yield (STY: fraction of simulations getting
#' over half FMSY yield in the first ten years of the projection), variability
#' in yield (VY: fraction of simulations where average annual variability in
#' yield is less than 10 per cent) and biomass level (B10: the fraction of
#' simulations in which biomass stays above 10 percent of BMSY).
#' @export
Tplot2_old  <- function(MSEobj, nam = NA) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  LTY <- rep(NA, MSEobj@nMPs)
  STY <- rep(NA, MSEobj@nMPs)
  VY <- rep(NA, MSEobj@nMPs)
  B10 <- rep(NA, MSEobj@nMPs)
  yend <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
  ystart <- 1:5
  RefYd <- MSEobj@OM$RefY
  y1 <- 1:(MSEobj@proyears - 1)
  y2 <- 2:MSEobj@proyears
  for (mm in 1:MSEobj@nMPs) {
    LTY[mm] <- round(sum(MSEobj@C[, mm, yend]/RefYd > 0.5, na.rm = T)/(MSEobj@nsim * 
                                                                         length(yend)), 3) * 100
    STY[mm] <- round(sum(MSEobj@C[, mm, ystart]/RefYd > 0.5, na.rm = T)/(MSEobj@nsim * 
                                                                           length(ystart)), 3) * 100
    AAVY <- apply(((MSEobj@C[, mm, y1] - MSEobj@C[, mm, y2])^2)^0.5, 
                  1, mean, na.rm = T)/apply(MSEobj@C[, mm, y2], 1, mean, na.rm = T)
    VY[mm] <- round(sum(AAVY < 0.1, na.rm = T)/MSEobj@nsim, 3) * 100
    B10[mm] <- round(sum(MSEobj@B_BMSY[, mm, ] > 0.1, na.rm = T)/prod(dim(MSEobj@B_BMSY[, 
                                                                                        mm, ])), 3) * 100
  }
  par(mfrow = c(1, 2), mar = c(5, 4, 1, 1), oma = c(0, 0, 2, 0))

  tradeoffplot(STY, LTY, "P(Short term yield > 0.5 FMSY)", "P(Long term yield > 0.5 FMSY)", 
               MSEobj@MPs[1:MSEobj@nMPs], vl = 1, hl = 1)
  tradeoffplot(B10, VY, "P(Biomass > 0.1 BMSY)", "P(CV in yield < 0.1)", 
               MSEobj@MPs[1:MSEobj@nMPs], vl = 1, hl = 1)
  if (is.na(nam)) 
    mtext(deparse(substitute(MSEobj)), 3, outer = T, line = 0.3, font = 2)
  if (!is.na(nam)) 
    mtext(MSEobj@Name, 3, outer = T, line = 0.3, font = 2)
  return(invisible())
}


# #' @describeIn Tplot By default, trade-off plots among LTY, STY, and biomass level B50 
# #' (fraction of simulations in which biomass stays above 50 percent of BMSY), and 
# #' Average Annual Variability in Yield (AAVY).
# 
# #' @export
# Tplot3_old  <- function(MSEobj, ..., lims=c(0.2, 0.2, 0.8, 0.8)) {
#   PMlist <- unlist(list(...))
#   if(length(PMlist) == 0) PMlist <- c("LTY", "STY", "P50", "AAVY")
#   if (class(PMlist) != 'character') stop("Must provide names of PM methods")
#   # check
#   
#   for (X in seq_along(PMlist))
#     if (!PMlist[X] %in% avail("PM")) stop(PMlist[X], " is not a valid PM method")
#   if (length(PMlist)<2) stop("Must provided more than 1 PM method")
#   
#   runPM <- vector("list", length(PMlist))
#   for (X in 1:length(PMlist)) runPM[[X]] <- eval(call(PMlist[X], MSEobj))
#   
#   PlotList <- combn(unique(PMlist), 2)
#   lims <- rep(lims, 100)[1:length(PMlist)]
#   
#   n.col <- ceiling(sqrt(ncol(PlotList)))
#   n.row <- ceiling(ncol(PlotList)/n.col)
#   
#   m <- matrix(1:(n.col*n.row), ncol=n.col, nrow=n.row, byrow=FALSE)
#   xmin <- xmax <- ymin <- ymax <- x <- y <- Class <- label <- fontface <- NULL
#   plots <- listout <- list()
#   for (pp in 1:ncol(PlotList)) {
#     yPM <- PlotList[1,pp]
#     yvals <- runPM[[match(yPM, PMlist)]]@Mean
#     ycap <-  runPM[[match(yPM, PMlist)]]@Caption
#     yname <-  runPM[[match(yPM, PMlist)]]@Name
#     yline <- lims[match(yPM, PMlist)]
#     
#     xPM <- PlotList[2,pp]
#     xvals <- runPM[[match(xPM, PMlist)]]@Mean
#     xcap <-  runPM[[match(xPM, PMlist)]]@Caption
#     xname <-  runPM[[match(xPM, PMlist)]]@Name
#     xline <- lims[match(xPM, PMlist)]
#     
#     xlim <- c(0, max(max(xvals, 1)))
#     ylim <- c(0, max(max(yvals, 1)))
#     
#     xrect <- data.frame(xmin=0, xmax=xline, ymin=0, ymax=max(ylim))
#     yrect <- data.frame(xmin=0, xmax=max(xlim), ymin=0, ymax=yline)
#     
#     MPType <- MPtype(MSEobj@MPs)
#     Class <- MPType[match(MSEobj@MPs, MPType[,1]),2]
#     
#     df <- data.frame(x=xvals, y=yvals, label=MSEobj@MPs, Class=Class,
#                      pass=xvals>xline & yvals>yline, fontface="plain", xPM=xPM, yPM=yPM)
#     df$fontface <- as.character(df$fontface)
#     df$fontface[!df$pass] <- "italic"
#     df$fontface <- factor(df$fontface)
#     listout[[pp]] <- df
#     plots[[pp]] <- ggplot2::ggplot() + 
#       ggplot2::geom_rect(data=xrect, ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.4) +
#       ggplot2::geom_rect(data=yrect, ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.4)
#     
#     # plots[[pp]] <- ggplot(df, aes(x, y, shape=Class, color=Class, label=label)) +
#     plots[[pp]] <-   plots[[pp]] + 
#       ggplot2::geom_point(data=df, ggplot2::aes(x, y, shape=Class, color=Class), size=2) +
#       ggrepel::geom_text_repel(data=df, ggplot2::aes(x, y, color=Class, label=label, fontface = fontface), show.legend=FALSE) + 
#       ggplot2::xlab(xcap) + ggplot2::ylab(ycap) +
#       ggplot2::xlim(xlim) + ggplot2::ylim(ylim) +
#       ggplot2::theme_classic() +
#       ggplot2::theme(axis.title.x = ggplot2::element_text(size=11),
#                      axis.title.y = ggplot2::element_text(size=11),
#                      legend.text=ggplot2::element_text(size=12)) + 
#       ggplot2::labs(shape= "MP Type", color="MP Type")
#     
#   }
#   out <- do.call("rbind", listout)
#   tab <- table(out$label, out$pass)
#   passall <- rownames(tab)[tab[,ncol(tab)] == ncol(PlotList)]
#   tt <- summary(MSEobj, PMlist, silent=TRUE)
#   tt$Satisificed <- FALSE
#   tt$Satisificed[match(passall, tt$MP)] <- TRUE
#   
#   join_plots(plots, n.col, n.row)
#   tt
#   
# }
# 
# 
# 
#' Generic Trade-off Plot
#' 
#' Creates a trade-off plot (up to four panels) of built-in performance
#' metrics.
#' 
#' Returns a list containing the names of performance metrics that meet the
#' minimum performance metrics for each trade-off, and ranks the MPs by
#' increasing distance from the top-right corner.
#' 
#' @param MSEobj Object of class MSE, output of the runMSE function
#' @param XAxis Character string describing the performance metrics for the
#' x-axis (or x-axes if vector; max 4). Must be chosen for list of existing PMs
#' and same length as YAxis. See \code{PMs}
#' @param YAxis Character string describing the performance metrics for the
#' y-axis (or y-axes if vector; max 4). Must be chosen for list of existing PMs
#' and same length as XAxis. See \code{PMs}
#' @param XThresh Minimum threshold values in percent (i.e., 50 = 50\%) for the
#' x-axes (must be same length as XAxis)
#' @param YThresh Minimum threshold values in percent (i.e., 50 = 50\%) for the
#' y-axes (must be same length as YAxis)
#' @param maxVar Reference for average annual variability in yield in percent
#' @param BmsyRef Reference level of BMSY, in proportion, i.e., 0.5 = 0.5BMSY
#' @param B0Ref Reference level of B0, in proportion, i.e., 0.2 = 0.2B0
#' @param AvailMPs vector of MPs that *could* be applied to the fishery, i.e.,
#' sufficient data exists. These a plotted with different symbol
#' @param ShowLabs Logical to specify if MP labels are shown
#' @param ShowCols Logical to specify if background colors are shown
#' @author A. Hordyk
#' @export 
TradePlot_old <- function(MSEobj, XAxis = c("Overfishing", "Biomass:BMSY"), 
                      YAxis = c("Long-term Yield", "AnnualVar"), XThresh = c(30, 80), 
                      YThresh = c(0, 50), maxVar = 15, BmsyRef = 0.5, B0Ref = 0.2, 
                      AvailMPs = NULL, ShowLabs = FALSE, ShowCols = TRUE) {
  PMs <- c("Long-term Yield", "Short-term Yield", "Overfishing", "Biomass:BMSY", 
           "Biomass:B0", "AnnualVar")
  op <- par(no.readonly=TRUE)
  on.exit(par(op))
  # Error Checks
  if (prod(XAxis %in% PMs) != 1) {
    message("Available Performance Metrics")
    print(PMs)
    stop("Invalid XAxis Performance Metrics")
  }
  if (prod(YAxis %in% PMs) != 1) {
    message("Available Performance Metrics")
    print(PMs)
    stop("Invalid YAxis Performance Metrics")
  }
  if (length(XAxis) > 4) 
    stop("Too many Performance Metrics (max 4)")
  if (length(YAxis) > 4) 
    stop("Too many Performance Metrics (max 4)")
  if (length(XAxis) != length(YAxis)) 
    stop("XAxis must be of same length as YAxis")
  if (length(XThresh) != length(XAxis) | length(YThresh) != length(XAxis)) 
    warning("Risk Threshold not same length as number of PMs")
  
  Yd <- rep(NA, MSEobj@nMPs)
  BMSYref <- rep(NA, MSEobj@nMPs)
  B0ref <- rep(NA, MSEobj@nMPs)
  PNOF <- rep(NA, MSEobj@nMPs)
  LTY <- rep(NA, MSEobj@nMPs)
  STY <- rep(NA, MSEobj@nMPs)
  VY <- rep(NA, MSEobj@nMPs)
  
  y1 <- 1:(MSEobj@proyears - 1)
  y2 <- 2:MSEobj@proyears
  
  ystart <- 1:5
  yend <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
  
  RefYd <- MSEobj@OM$RefY
  if (maxVar < 1) 
    maxVar <- maxVar * 100
  
  for (mm in 1:MSEobj@nMPs) {
    PNOF[mm] <- round(sum(MSEobj@F_FMSY[, mm, ] <= 1, na.rm = T)/prod(dim(MSEobj@F_FMSY[, 
                                                                                        mm, ]), na.rm = T) * 100, 1)
    BMSYref[mm] <- round(sum(MSEobj@B_BMSY[, mm, ] > BmsyRef, na.rm = T)/prod(dim(MSEobj@B_BMSY[, mm, ])) * 100, 1)
    B0ref[mm] <- round(sum((MSEobj@B_BMSY[, mm, ] * MSEobj@OM$SSBMSY_SSB0) > 
                             B0Ref, na.rm = T)/prod(dim(MSEobj@B_BMSY[, mm, ])) * 100, 1)
    # LTY[mm]<-round(sum(MSEobj@C[,mm,yend]/RefYd>0.5,na.rm=T)/(MSEobj@nsim*length(yend)),3)*100
    # STY[mm]<-round(sum(MSEobj@C[,mm,ystart]/RefYd>0.5,na.rm=T)/(MSEobj@nsim*length(ystart)),3)*100
    LTY[mm] <- round(mean(apply(MSEobj@C[, mm, yend], 1, mean, na.rm = T)/RefYd, 
                          na.rm = T) * 100, 1)
    STY[mm] <- round(mean(apply(MSEobj@C[, mm, ystart], 1, mean, na.rm = T)/RefYd, 
                          na.rm = T) * 100, 1)
    AAVY <- apply((((MSEobj@C[, mm, y1] - MSEobj@C[, mm, y2])/MSEobj@C[, 
                                                                       mm, y2])^2)^0.5, 1, mean, na.rm = T)
    VY[mm] <- round(sum(AAVY < (maxVar/100), na.rm = T)/MSEobj@nsim, 
                    3) * 100
  }
  
  for (xx in seq_along(XAxis)) {
    name <- paste0("X", xx)
    name1 <- paste0("XLab", xx)
    assign(name, GetStat(XAxis[xx], LTY, STY, PNOF, BMSYref, B0ref, 
                         VY))
    assign(name1, StatLab(XAxis[xx], maxVar, BmsyRef, B0Ref))
    name <- paste0("Y", xx)
    name1 <- paste0("YLab", xx)
    assign(name, GetStat(YAxis[xx], LTY, STY, PNOF, BMSYref, B0ref, 
                         VY))
    assign(name1, StatLab(YAxis[xx], maxVar, BmsyRef, B0Ref))
  }
  
  Nplot <- length(XAxis)
  if (Nplot == 1) 
    par(mfrow = c(1, 1), mar = c(4, 4.5, 1, 1), oma = c(1, 1, 0, 0))
  if (Nplot == 2) 
    par(mfrow = c(1, 2), mar = c(4, 4.5, 1, 1), oma = c(1, 1, 0, 0))
  if (Nplot == 3) 
    par(mfrow = c(1, 3), mar = c(4, 4.5, 1, 1), oma = c(1, 1, 0, 0))
  if (Nplot == 4) 
    par(mfrow = c(2, 2), mar = c(4, 4.5, 1, 1), oma = c(1, 1, 0, 0))
  
  OutList <- list()
  for (xx in seq_along(XAxis)) {
    Xname <- paste0("X", xx)
    XLab <- paste0("XLab", xx)
    Yname <- paste0("Y", xx)
    YLab <- paste0("YLab", xx)
    rr <- tradeoffplot4(x = get(Xname), y = get(Yname), get(XLab), 
                        get(YLab), labs = MSEobj@MPs[1:MSEobj@nMPs], vl = XThresh[xx], 
                        hl = YThresh[xx], ShowLabs = ShowLabs, ShowCols = ShowCols, 
                        AvailMPs = AvailMPs)
    
    labs <- MSEobj@MPs[1:MSEobj@nMPs]
    ind <- which(labs %in% rr)
    tempDF <- data.frame(MP = rr, X = get(Xname)[ind], Y = get(Yname)[ind])
    Dist <- NULL  # calculate distance from corner
    for (X in 1:length(tempDF[, 2])) Dist[X] <- euc.dist(c(tempDF[X, 
                                                                  2], tempDF[X, 3]), c(100, 100))
    tempDF <- tempDF[order(Dist), ]
    rownames(tempDF) <- 1:nrow(tempDF)
    OutList[[xx]] <- tempDF
  }
  
  
  OutList
  
}

#' Calculate Value Of Information
#' 
#' A function that relates operating model parameters and parameters of the
#' observation model to yield (by default). A user can also specific their own
#' utility values (Ut) which is arranged in a matrix of nsim rows and nMP
#' columns.
#' 
#' 
#' @param MSEobj An object of class MSE
#' @param ncomp Maximum number of variables to examine per MP
#' @param nbins Number of percentile bins for sampled parameters of the
#' operating model or observation model, which is used for calculating
#' variability in utility across the sampled range of each parameter
#' @param maxrow maximum number of MPs per plot
#' @param Ut A matrix of user-specified utility values of nsim rows and nMPs
#' columns
#' @param Utnam The name of the utility measure for plotting
#' @param plot Logical. Show the plot?
#' @author T. Carruthers
#' @export VOI
VOI <- function(MSEobj, ncomp = 6, nbins = 8, maxrow = 8, Ut = NA, Utnam = "Utility", plot=TRUE) {
  objnam <- deparse(substitute(MSEobj))
  nsim <- MSEobj@nsim
  
  if (is.na(Ut[1])) {
    Ut <- array(NA, c(nsim, MSEobj@nMPs))
    yind <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
    RefYd <- MSEobj@OM$RefY
    
    for (mm in 1:MSEobj@nMPs) {
      Ut[, mm] <- apply(MSEobj@C[, mm, yind], 1, mean, na.rm = T)/RefYd * 100
      # POF[,mm]<-apply(MSEobj@F_FMSY[,mm,]>1,1,sum)/MSEobj@proyears
      # P10[,mm]<-apply(MSEobj@B_BMSY[,mm,]<0.1,1,sum)/MSEobj@proyears
    }
    Utnam <- "Long-term yield relative to MSY (%)"
  }
  
  MPs <- MSEobj@MPs
  nMPs <- MSEobj@nMPs
  
  onlycor <- c("RefY", "A", "MSY", "Linf", "t0", "OFLreal", "Spat_targ")
  vargood <- (apply(MSEobj@OM, 2, sd)/(apply(MSEobj@OM, 2, mean)^2)^0.5) >  0.005
  # MSEobj@OM<- MSEobj@OM[,(!names(MSEobj@OM)%in%onlycor)&vargood]
  vargood[grep("qvar", names(MSEobj@OM))] <- FALSE
  MSEobj@OM <- MSEobj@OM[, which((!names(MSEobj@OM) %in% onlycor) & vargood)]
  OMp <- apply(MSEobj@OM, 2, quantile, p = seq(0, 1, length.out = nbins + 1), na.rm = TRUE)
  Obsp <- apply(MSEobj@Obs, 2, quantile, p = seq(0, 1, length.out = nbins + 1), na.rm = TRUE)
  OMv <- array(NA, c(nMPs, ncol(MSEobj@OM), nbins))
  Obsv <- array(NA, c(nMPs, ncol(MSEobj@Obs), nbins))
  
  for (mm in 1:nMPs) {
    for (j in 1:nbins) {
      for (i in 1:ncol(MSEobj@OM)) {
        cond <- MSEobj@OM[, i] > OMp[j, i] & MSEobj@OM[, i] < OMp[j + 1, i]
        OMv[mm, i, j] <- mean(Ut[cond, mm], na.rm = T)
      }
      for (i in 1:ncol(MSEobj@Obs)) {
        cond <- MSEobj@Obs[, i] > Obsp[j, i] & MSEobj@Obs[, i] < Obsp[j + 1, i]
        Obsv[mm, i, j] <- mean(Ut[cond, mm], na.rm = T)
      }
    }
  }
  
  # cbind(names(MSEobj@OM), OMs[mm,])
  
  # -- Operating model variables
  OMs <- apply(OMv, 1:2, sd, na.rm = T)
  OMstr <- array("", c(nMPs * 2, ncomp + 1))
  
  for (mm in 1:nMPs) {
    ind <- order(OMs[mm, ], decreasing = T)[1:ncomp]
    OMstr[1 + (mm - 1) * 2, 1] <- MPs[mm]
    OMstr[1 + (mm - 1) * 2, 2:(1 + ncomp)] <- names(MSEobj@OM[ind])
    OMstr[2 + (mm - 1) * 2, 2:(1 + ncomp)] <- round(OMs[mm, ind], 2)
  }
  OMstr <- data.frame(OMstr)
  names(OMstr) <- c("MP", 1:ncomp)
  
  
  # -- Observation model variables
  slots <- c("Cat", "Cat", "AvC", "AvC", "CAA", "CAA", "CAL", "CAL", 
             "Ind", "Dep", "Dep", "Dt", "Dt", "Mort", "FMSY_M", "SSBMSY_SSB0", "L50", 
             "L95", "LFC", "LFS", "Abun", "Abun", "vbK", "vbt0", "vbLinf", "Steep", 
             "Iref", "Cref", "Bref", "ML")
  Obsnam <- c("Cbias", "Csd", "Cbias", "Csd", "CAA_nsamp", "CAA_ESS", 
              "CAL_nsamp", "CAL_ESS", "Isd", "Dbias", "Derr", "Dbias", "Derr", 
              "Mbias", "FMSY_Mbias", "BMSY_B0bias", "lenMbias", "lenMbias", "LFCbias", 
              "LFSbias", "Abias", "Aerr", "Kbias", "t0bias", "Linfbias", "hbias", 
              "Irefbias", "Crefbias", "Brefbias", "")
  Obss <- apply(Obsv, 1:2, sd, na.rm = T)
  Obsstr <- array("", c(nMPs * 2, ncomp + 1))
  for (mm in 1:nMPs) {
    relobs <- Obsnam[slots %in% unlist(strsplit(Required(MPs[mm])[, 
                                                                  2], split = ", "))]
    ind <- (1:ncol(MSEobj@Obs))[match(relobs, names(MSEobj@Obs))]
    pos <- names(MSEobj@Obs)[ind]  # possible observation processes
    maxy <- min(max(1, length(pos)), ncomp, na.rm = T)
    ind2 <- order(Obss[mm, ind], decreasing = T)[1:maxy]
    Obsstr[1 + (mm - 1) * 2, 1] <- MPs[mm]
    Obsstr[1 + (mm - 1) * 2, 2:(1 + maxy)] <- pos[ind2]
    Obsstr[2 + (mm - 1) * 2, 2:(1 + maxy)] <- round(Obss[mm, ind][ind2], 
                                                    2)
  }
  Obsstr <- data.frame(Obsstr)
  names(Obsstr) <- c("MP", 1:ncomp)
  
  ncols <- 40
  # colsse<-makeTransparent(rainbow(ncols,start=0,end=0.36),95)[ncols:1]
  colsse <- makeTransparent(rainbow(ncols, start = 0, end = 0.36), 90)[ncols:1]
  minsd <- 0
  maxsd <- max(OMs, na.rm = T)
  coly <- ceiling(OMs/maxsd * ncols)
  
  # Operating model variables
  mbyp <- split(1:nMPs, ceiling(1:nMPs/maxrow))
  ylimy = c(0, max(OMv, na.rm = T) * 1.2)
  
  if (plot) {
    for (pp in 1:length(mbyp)) {
      
      op <- par(mfrow = c(length(mbyp[[pp]]), ncomp), mai = c(0.15, 0.1, 0.15, 
                                                              0.05), omi = c(0.1, 0.9, 0.3, 0.05))
      
      for (mm in mbyp[[pp]]) {
        for (cc in 1:ncomp) {
          rind <- (mm - 1) * 2 + 1
          y <- Ut[, mm]
          cind <- match(OMstr[rind, 1 + cc], names(MSEobj@OM))
          x <- MSEobj@OM[, cind]
          plot(x, y, col = "white", axes = F, ylim = ylimy)
          axis(1, pretty(OMp[, cind]), pretty(OMp[, cind]), cex.axis = 0.8, 
               padj = -1.5)
          abline(v = OMp[, cind], col = "#99999960")
          points(x, y, col = colsse[coly[mm, cind]], pch = 19, cex = 0.8)
          x2 <- (OMp[1:nbins, cind] + OMp[2:(nbins + 1), cind])/2
          y2 <- OMv[mm, cind, ]
          lines(x2, y2)
          legend("bottomright", legend = round(OMs[mm, cind], 2), bty = "n", cex = 0.8)
          legend("topleft", legend = OMstr[rind, 1 + cc], bty = "n", cex = 0.85)
          if (cc == 1)  {
            mtext(MPs[mm], 2, font = 2, outer = F, cex = 0.8, line = 2)
            ytick <- pretty(seq(ylimy[1], ylimy[2] * 1.3, length.out = 10))
            axis(2, ytick, ytick, cex.axis = 0.8)
          }  # only first column
        }  # parameters (columns)
      }  # MPs (rows)
      
      mtext(Utnam, 2, outer = T, cex = 0.9, line = 3.5)
      mtext(paste("Operating model parameters: ", objnam, "@OM", sep = ""), 3, outer = T, font = 2, cex = 0.9)
      
    }  # Plots
    
    # Observation model values
    
    ylimy = c(0, max(Obsv, na.rm = T) * 1.2)
    minsd <- 0
    maxsd <- max(Obss)
    coly <- ceiling(Obss/maxsd * ncols)
    
    if (sum(is.na(Obsstr) | Obsstr == "") < (ncomp + 1) * nMPs * 2 - nMPs) 
    {
      # only if there is data to plot
      
      for (pp in 1:length(mbyp)) {
        
        op <- par(mfrow = c(length(mbyp[[pp]]), ncomp), mai = c(0.15, 
                                                                0.1, 0.15, 0.05), omi = c(0.1, 0.9, 0.3, 0.05),no.readonly=TRUE)
        on.exit(par(op))
        
        for (mm in mbyp[[pp]]) {
          rind <- (mm - 1) * 2 + 1
          npres <- sum(Obsstr[rind + 1, ] != "")
          for (cc in 1:ncomp) {
            if (!is.na(npres) & cc < (npres + 1)) {
              y <- Ut[, mm]
              cind <- match(Obsstr[rind, 1 + cc], names(MSEobj@Obs))
              x <- MSEobj@Obs[, cind]
              plot(x, y, col = "white", axes = F, ylim = ylimy)
              axis(1, pretty(Obsp[, cind]), pretty(Obsp[, cind]), cex.axis = 0.8, padj = -2)
              abline(v = Obsp[, cind], col = "#99999960")
              points(x, y, col = colsse[coly[mm, cind]], pch = 19, cex = 0.8)
              x2 <- (Obsp[1:nbins, cind] + Obsp[2:(nbins + 1), cind])/2
              y2 <- Obsv[mm, cind, ]
              lines(x2, y2)
              legend("bottomright", legend = round(Obss[mm, cind], 2), bty = "n", cex = 0.8)
              legend("topleft", legend = Obsstr[rind, 1 + cc], bty = "n", cex = 0.75)
              if (cc == 1)    {
                mtext(MPs[mm], 2, font = 2, outer = F, cex = 0.6, line = 2)
                ytick <- pretty(seq(ylimy[1], ylimy[2] * 1.3, length.out = 10))
                axis(2, ytick, ytick, cex.axis = 0.8)
              }  # only first column
            } else {
              plot(0, type = "n", axes = FALSE, ann = FALSE)
              if (cc == 1)  {
                mtext(MPs[mm], 2, font = 2, outer = F, cex = 0.6, line = 2)
              }  # only first column
            }
          }  # parameters (columns)
        }  # MPs (rows)
        
        mtext(Utnam, 2, outer = T, cex = 0.9, line = 3.5)
        mtext(paste("Observation model parameters: ", objnam, "@Obs", sep = ""), 3, outer = T, font = 2, cex = 0.9)
        
      }  # Plots
    }  # if there is data to plot
    
  }

  list(OMstr, Obsstr)
  
}  # VOI



#' Calculate Value Of Information 2
#' 
#' A function that relates operating model parameters and parameters of the
#' observation model to relative yield (yield over last 5 years of projection
#' relative to a 'best F' scenario that maximizes yield).
#' 
#' @param MSEobj An object of class MSE
#' @param ncomp Maximum number of observation variables to examine per MP
#' @param nbins Number of bins for sampled observation variables used for
#' calculating variability in utility across the sampled range of each
#' parameter
#' @param Ut A matrix of user-specified utility values of nsim rows and nMPs
#' columns
#' @param Utnam The name of the utility measure for plotting
#' @param lay Controls whether labels are in lay terms or not
#' @note VOI2 assumes that relative cost for each type of improvement in data
#' is linearly related to the number of samples (e.g. nCAAobs) or square
#' function of improved precision and bias e.g.: relative cost=
#' 1/(newCV/oldCV)^2
#' @author T. Carruthers
#' @export VOI2
VOI2 <- function(MSEobj, ncomp = 6, nbins = 4, Ut = NA, Utnam = "yield", 
                 lay = F) {
  op <- par(no.readonly=TRUE)
  on.exit(par(op))
  objnam <- deparse(substitute(MSEobj))
  nsim <- MSEobj@nsim
  
  if (is.na(Ut[1])) {
    Ut <- array(NA, c(nsim, MSEobj@nMPs))
    yind <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
    RefYd <- MSEobj@OM$RefY
    
    for (mm in 1:MSEobj@nMPs) {
      Ut[, mm] <- apply(MSEobj@C[, mm, yind], 1, mean, na.rm = T)/RefYd * 
        100
      # POF[,mm]<-apply(MSEobj@F_FMSY[,mm,]>1,1,sum)/MSEobj@proyears
      # P10[,mm]<-apply(MSEobj@B_BMSY[,mm,]<0.1,1,sum)/MSEobj@proyears
    }
    
  }
  
  MPs <- MSEobj@MPs
  nMPs <- MSEobj@nMPs
  
  # -- Observation model variables
  slots <- c("Cat", "Cat", "AvC", "AvC", "CAA", "CAA", "CAL", "CAL", 
             "Ind", "Ind", "Dep", "Dep", "Dt", "Dt", "Mort", "FMSY_M", "SSBMSY_SSB0", 
             "L50", "L95", "LFC", "LFS", "Abun", "Abun", "vbK", "vbt0", "vbLinf", 
             "Steep", "Iref", "Cref", "Bref")
  Obsnam <- c("Cbias", "Csd", "Cbias", "Csd", "CAA_nsamp", "CAA_ESS", 
              "CAL_nsamp", "CAL_ESS", "Isd", "betas", "Dbias", "Derr", "Dbias", 
              "Derr", "Mbias", "FMSY_Mbias", "BMSY_B0bias", "lenMbias", "lenMbias", 
              "LFCbias", "LFSbias", "Abias", "Aerr", "Kbias", "t0bias", "Linfbias", 
              "hbias", "Irefbias", "Crefbias", "Brefbias")
  
  
  Obsnam2 <- c("Cbias", "Csd", "CAA_nsamp", "CAA_ESS", "CAL_nsamp", "CAL_ESS", 
               "Isd", "betas", "Dbias", "Derr", "Mbias", "FMSY_Mbias", "BMSY_B0bias", 
               "lenMbias", "LFCbias", "LFSbias", "Abias", "Aerr", "Kbias", "t0bias", 
               "Linfbias", "hbias", "Irefbias", "Crefbias", "Brefbias")
  Obsnam3 <- c("Catch bias", "Catch error", "n CAA samples", "CAA ESS", 
               "n CAL samples", "CAL ESS", "Abun. Ind. error", "Hyperstability", 
               "Depln. bias", "Depln. error", "M bias", "FMSY/M bias", "BMSY/B0 bias", 
               "lenMbias", "Len 1st Cap bias", "Len full sel bias", "Cur Abund bias", 
               "Cur Abun err", "vB K bias", "vB t0 bias", "vB Linf bias", "Steepness bias", 
               "Ref index bias", "Ref catch bias", "Ref biomass bias")
  # Types of observation error model 1:lognorm 2:percentile 3:replicates
  # (higher is better) ##4:uniform on log 5:logit space
  oem <- c(1, 2, 3, 3, 3, 3, 2, 4, 1, 2, 1, 1, 1, 1, 1, 1, 4, 2, 1, 1, 
           1, 2, 1, 1, 1)
  # oem<-c( 2, 2, 3, 3, 3, 3, 2, 4, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2, 2, 2,
  # 2, 2, 1, 1, 1)
  
  Obsd <- apply(MSEobj@Obs, 2, sd, na.rm = TRUE)
  Obm <- apply(MSEobj@Obs, 2, mean, na.rm = TRUE)
  Obmd <- apply(MSEobj@Obs, 2, quantile, p = 0.5, na.rm = TRUE)
  
  maxcomp <- length(Obsnam2)
  Obsv <- array(NA, c(nMPs, maxcomp, nbins))
  Obsval <- array(NA, c(nMPs, maxcomp, nbins))
  Obscost <- array(NA, c(nMPs, maxcomp, nbins))
  Obsname <- list()
  
  div <- seq(1, 2, length.out = nbins + 1)[2:(nbins + 1)]  # for distributions
  percs <- seq(0.5, 1, length.out = nbins + 1)[1:nbins]  # for samples
  percsCAA <- seq(0, 1, length.out = nbins + 2)[2:(nbins + 1)]
  percUL <- seq(0, 0.25, length.out = nbins + 1)[2:(nbins + 1)]
  percUU <- 1 - percUL
  for (mm in 1:nMPs) {
    Y1 <- Ut[, mm]
    relobs <- Obsnam[slots %in% unlist(strsplit(Required(MPs[mm])[, 
                                                                  2], split = ", "))]
    Obsname[[mm]] <- relobs
    nr <- length(relobs)
    if (length(relobs) > 0) 
    {
      
      for (r in 1:nr) {
        oemi <- match(relobs[r], Obsnam2)
        obsi <- match(relobs[r], names(MSEobj@Obs))
        for (cc in 1:nbins) {
          if (oem[oemi] == 1) {
            # Redundant SIR code for log-normal biases
            T1 <- tdlnorm(MSEobj@Obs[, obsi], Obm[obsi], Obsd[obsi]/Obm[obsi])
            # plot(MSEobj@Obs[,obsi],T1) # check
            T2 <- tdlnorm(MSEobj@Obs[, obsi], Obm[obsi], Obsd[obsi]/(div[cc] * 
                                                                       Obm[obsi]))
            W <- T2/T1
            nrep2 <- nsim * 20
            Y2 <- sample(Y1, nrep2 * 5, replace = T, prob = W)
            Obsv[mm, r, cc] <- (mean(Y2) - mean(Y1))/mean(Y1) * 
              100
            Obsval[mm, r, cc] <- Obsd[obsi]/(div[cc] * Obm[obsi])
            Obscost[mm, r, cc] <- div[cc]^2
          } else if (oem[oemi] == 2) {
            refval <- quantile(MSEobj@Obs[, obsi], percs[nbins:1][cc], 
                               na.rm = TRUE)
            ind <- MSEobj@Obs[, obsi] < refval
            Obsv[mm, r, cc] <- (mean(Y1[ind]) - mean(Y1))/mean(Y1) * 
              100
            Obsval[mm, r, cc] <- mean(MSEobj@Obs[ind, obsi])
            Obscost[mm, r, cc] <- 1/(Obsval[mm, r, cc]/mean(MSEobj@Obs[, 
                                                                       obsi]))^2
          } else if (oem[oemi] == 3) {
            refval <- quantile(MSEobj@Obs[, obsi], percsCAA[cc], 
                               na.rm = TRUE)
            ind <- MSEobj@Obs[, obsi] > refval
            Obsv[mm, r, cc] <- (mean(Y1[ind]) - mean(Y1))/mean(Y1) * 
              100
            Obsval[mm, r, cc] <- mean(MSEobj@Obs[ind, obsi])
            Obscost[mm, r, cc] <- Obsval[mm, r, cc]
          } else if (oem[oemi] == 4) {
            refval <- quantile(MSEobj@Obs[, obsi], percUL[cc], 
                               na.rm = TRUE)
            refval2 <- quantile(MSEobj@Obs[, obsi], percUU[cc], 
                                na.rm = TRUE)
            ind <- (MSEobj@Obs[, obsi] > refval) & (MSEobj@Obs[, 
                                                               obsi] < refval2)
            Obsv[mm, r, cc] <- (mean(Y1[ind]) - mean(Y1))/mean(Y1) * 
              100
            Obsval[mm, r, cc] <- sd(MSEobj@Obs[ind, obsi])
            Obscost[mm, r, cc] <- 1/(Obsval[mm, r, cc]/sd(MSEobj@Obs[, 
                                                                     obsi]))^2
          }
          # observation model type
        }  # loop over bins
      }  # loop over r
    }  # observation variables?
  }  # loop over MPs
  
  cb <- array(NA, c(MSEobj@nMPs, maxcomp))
  for (mm in 1:MSEobj@nMPs) {
    if (sum(!is.na(Obscost[mm, , ])) > 0) {
      for (r in 1:length(Obsname[[mm]])) {
        dat <- data.frame(x = Obscost[mm, r, ], y = Obsv[mm, r, 
                                                         ])
        if (prod(apply(dat, 2, is.finite)) > 0) {
          # plot(dat$x,dat$y)
          cb[mm, r] <- lm(y ~ x - 1, data = dat)$coefficients[1]
        }
      }
    }
  }
  
  ncols <- 100
  # colsse<-makeTransparent(rainbow(ncols,start=0,end=0.36),95)[ncols:1]
  colt <- rainbow(ncols, start = 0, end = 0.36)[1:ncols]
  colsse <- makeTransparent(colt, 98)
  
  cb[cb < 0 | is.na(cb)] <- 0
  coly <- ceiling((cb/max(cb, na.rm = T))^0.5 * ncols)
  coly[coly == 0] <- 1
  
  ncol <- ceiling(MSEobj@nMPs^0.5)
  nrow <- ceiling(MSEobj@nMPs/ncol)
  
  op <- par(mfrow = c(nrow, ncol), mar = c(2.4, 2.4, 0.1, 0.1), omi = c(0.4, 
                                                                        0.35, 0.3, 0))
  
  gcol1 <- "#99999960"
  gcol2 <- "#99999940"
  gcol3 <- "#99999920"
  
  for (mm in 1:MSEobj@nMPs) {
    if (sum(!is.na(Obscost[mm, , ])) > 0) 
    {
      
      plot(c(1, 5), range(Obsv, na.rm = T), col = "white", main = "")
      legend("topleft", legend = MSEobj@MPs[mm], bty = "n", text.font = 2, 
             cex = 1.4)
      
      abline(h = (-20:50) * 4, col = gcol2, lwd = 1.5)
      abline(h = (-20:50) * 4 + 2, col = gcol3, lwd = 1)
      
      abline(v = 1:4, col = gcol2, lwd = 1.5)
      abline(v = (1:4) + 0.5, col = gcol3, lwd = 1)
      abline(h = 0, col = gcol1, lwd = 3)
      
      
      no <- length(Obsname[[mm]])
      ind <- order(cb[mm, 1:no], decreasing = T)[1:ncomp]
      ind <- ind[!is.na(ind)]
      
      ind2 <- order(Obsv[mm, 1:no, nbins])[1:ncomp]
      ind2 <- ind2[!is.na(ind2)]
      
      
      lpos <- Obsv[mm, ind2, nbins]
      ppos <- seq(min(Obsv, na.rm = T), max(Obsv, na.rm = T), 
                  length.out = length(ind2))
      wt <- (max(Obsv[mm, 1:no, nbins], na.rm = T) - min(Obsv[mm, 
                                                              1:no, nbins], na.rm = T))/(max(Obsv, na.rm = T) - min(Obsv, 
                                                                                                                    na.rm = T))/(no/ncomp)
      wt <- wt^0.66
      nupos <- wt * lpos + (1 - wt) * ppos
      
      for (r2 in 1:length(ind)) {
        r <- ind2[r2]
        lines(c(1, Obscost[mm, r, ]), c(0, Obsv[mm, r, ]), col = colsse[coly[mm, 
                                                                             r]], lwd = 3)
        if (!lay) {
          text(4.5, nupos[r2], Obsname[[mm]][r], col = colt[coly[mm, 
                                                                 r]], font = 2, cex = 1.2)
        } else {
          text(4.5, nupos[r2], Obsnam3[match(Obsname[[mm]][r], 
                                             Obsnam2)], col = colt[coly[mm, r]], font = 2, cex = 1.2)
        }
      }  # observation quantities (lines)
      
      
      # legend('topleft',legend=Obsname[[mm]][ind],text.col=colt[coly[mm,ind]],text.font=2,cex=1.2,bty='n')
      
    }  # if there is data to plot
    
    
  }  # MPs (plots)
  
  
  mtext("Cost relative to today", 1, outer = T, cex = 0.9, line = 1, 
        font = 2)
  # mtext(paste('Operating model parameters:
  # ',objnam,'@OM',sep=''),3,outer=T,font=2,cex=0.9)
  mtext(paste("% Change in ", Utnam, " relative to today", sep = ""), 
        2, outer = T, line = 0.6, font = 2, cex = 0.9)
  par(op)
  invisible(list(Obscost, Obsv, Obsval, cb, Obsname, MSEobj@MPs))
  
}  # VOI2

#' Yet another Value of Information Plot
#' 
#' A function that relates parameters of the observation model and the
#' operating model parameters to yield.
#' 
#' 
#' @param MSEobj An object of class MSE
#' @param MPs The MPs to plot. If NA it will plot the first nMP from MSEobj
#' @param nvars The number of observation or operating model parameters to plot
#' (number of columns)
#' @param nMP The maximum number of MPs to plot (number of rows)
#' @param Par Plot Operating Model (OM) or Observation (Obs) parameters?
#' @param YVar Variable for Y-Axis: Yield (Y) or Biomass (B) (relative to BMSY)
#' @param doPlot Output the plot?
#' @param incStat Include a print out of statistic describing the curviness of
#' the line?
#' @param availMP Optional character string of MPs that are available. These
#' names are colored black
#' @param acceptMP Optional character string of MPs that are acceptable. These
#' names are colored green if they are also in availMP
#' @param incNames Include the names?
#' @param labcex Character size of the label
#' @param quants Quantiles to calculate
#' @return A list of all the information included in the plot
#' @author A. Hordyk
#' @export VOIplot
VOIplot <- function(MSEobj, MPs = NA, nvars = 5, nMP = 4, 
                    Par = c("Obs", "OM"), YVar = c("Y", "B"), doPlot = TRUE, incStat = FALSE, 
                    availMP = NULL, acceptMP = NULL, incNames = TRUE, labcex = 0.8, 
                    quants = c(0.05, 0.95)) {
  
  YVar <- match.arg(YVar)
  nvars <- max(nvars, 2)  # maximum number of variables 
  Par <- match.arg(Par)  # Operating Model or Observation 
  nMPs <- MSEobj@nMPs  # Number of MPs  
  # Subset to specified MPs
  if (any(is.na(MPs))) MPs <- MSEobj@MPs
  if (class(MPs) == "numeric" | class(MPs) == "integer")   MPs <- MSEobj@MPs[MPs]
  if (length(MPs) < 1)   stop("No MPs found")
  nMPss <- length(MPs) 
  if (nMP > nMPs)     nMP <- nMPs
  if (!all(MSEobj@MPs %in% MPs)) {
    mse <- Sub(MSEobj, MPs = MPs)
    nMPs <- mse@nMPs
  } else {
    mse <- MSEobj
  }
  
  # Calculate MSE sensitivities per MP
  if (length(MPs) > 1) senseDat <- sapply(1:nMPs, calcMSESense, MSEobj = mse, 
                                          YVar = YVar, Par = Par, 
                                          simplify = FALSE, quants = quants)
  if (length(MPs) == 1) senseDat <- calcMSESense(MP = MPs, MSEobj = mse, 
                                                 YVar = YVar, Par = Par,
                                                 quants = quants)
  
  # Operating Model or Observation Statistics
  
  if (nMPs == 1) {
    varNames <- senseDat$OMNames
  } else varNames <- senseDat[[1]]$OMNames
  
  used <- matrix(FALSE, nrow = length(varNames), ncol = nMPs)
  if (Par == "OM") {
    used <- matrix(TRUE, nrow = length(varNames), ncol = nMPs)  # all OM parameters used 
    Obsnam <- varNames
    LnName <- varNames
    # LnName <- c("Reference yield", "Natural mortality", "Depletion", "Abundance",
    #             "BMSY/B0", "FMSY/M", "M gradient", "Inter-annual variability M", 
    #             "Recruitment variability", "Inter-annual variability effort", 
    #             "Final effort",  "MSY", "Average change in catchability", "Inter-annual variabilility in catchability", 
    #             "FMSY", "von Bert. Linf", "von Bert. K", "von Bert. t0", "Steepness", 
    #             "Linf gradient", "K gradient", "Inter-annual variability in Linf", 
    #             "Recruitment gradient", "Inter-annual variability in K", "Age at maturity", 
    #             "Length at 5% selection", "Length at full selection", "Dome-shaped selectivity",
    #             "Length at first capture", "Auto-correlation recruitment", 
    #             "Length 50% maturity", "Length 95% maturity", "B0", "N0", "SSB0", "BMSY_B0",
    #             "TACSD", "TACFrac", "TAESD", "TAEFrac", "SizeLimSD", "SizeLimFrac", "Blow", "BMSY",
    #             "SSBMSY", "Mexp", "Discard mortality", "LR5", "LFR", "DR", "Lm/SL")
    #  # cbind(Obsnam, LnName)
    
  }
  if (Par == "Obs") {
    slots <- c("Cat", "Cat", "AvC", "AvC", "CAA", "CAA", "CAL", "CAL", 
               "Ind", "Ind", "Dep", "Dep", "Dt", "Dt", "Mort", "FMSY_M", "SSBMSY_SSB0", 
               "L50", "L95", "LFC", "LFS", "Abun", "Abun", "vbK", "vbt0", 
               "vbLinf", "Steep", "Iref", "Cref", "Bref", "ML", "ML")
    Obsnam <- c("Cbias", "Csd", "Cbias", "Csd", "CAA_nsamp", "CAA_ESS", 
                "CAL_nsamp", "CAL_ESS", "Isd", "betas", "Dbias", "Derr", "Dbias", 
                "Derr", "Mbias", "FMSY_Mbias", "BMSY_B0Bias", "lenMbias", "lenMbias", 
                "LFCbias", "LFSbias", "Abias", "Aerr", "Kbias", "t0bias", "Linfbias", 
                "hbias", "Irefbias", "Crefbias", "Brefbias", "CAL_nsamp", "CAL_ESS")
    LnName <- Obsnam
    # LnName <- c("Catch bias", "Catch error", "Catch bias", "Catch error", 
    #             "n CAA samples", "CAA effective sample size", "n CAL samples", 
    #             "CAL effective sample size", "Index Abundance error", "Hyperstability/hyperdepletion", 
    #             "Depletion bias", "Depletion error", "Depletion bias", "Depletion error", 
    #             "M bias", "FMSY/M bias", "BMSY/B0 bias", "Length maturity bias", 
    #             "Length maturity bias", "Length first capture bias", "Length full capture bias", 
    #             "Current abundance bias", "Current abundance error", "vB K bias", 
    #             "vB t0 bias", "vB Linf bias", "Steepness bias", "Reference index bias", 
    #             "Reference catch bias", "Reference biomass bias", "Mean length", 
    #             "Mean length")
    # print(cbind(slots, Obsnam, LnName))
    for (mm in 1:nMPs) {
      ids <- Obsnam[slots %in% unlist(strsplit(Required(MPs[mm])[, 2], split = ", "))]
      used[match(ids, varNames), mm] <- TRUE
    }
  }
  colnames(used) <- MPs
  rownames(used) <- varNames
  
  
  # Find the highest Stat for each variable
  if (nMPs > 1) {
    stat <- lapply(senseDat, "[[", "OMStat")
    Stat <- matrix(unlist(stat), ncol = nMPs) * used
    if (max(Stat, na.rm = TRUE) > 100) 
      Stat <- Stat/100
    rownames(Stat) <- varNames
    statord <- apply(Stat, 2, order, decreasing = TRUE)
    topStat <- statord[1:nvars, ]  # highest nvars for each MP
  } else {
    stat <- senseDat$OMStat
    Stat <- matrix(stat, ncol = nMPs)
    rownames(Stat) <- varNames
    colnames(Stat) <- MPs
    statord <- order(Stat, decreasing = TRUE)
    topStat <- statord[1:nvars]
  }
  
  
  if (doPlot) {
    ## Create Plotting Space ##
    Ncol <- nvars
    if (nMPs < 2) 
      Ncol <- min(sum(used[, MPs]), nvars)
    if (nMPs > 1) {
      temp <- apply(used[, MPs], 2, sum)
      if (all(temp < nvars)) 
        Ncol <- max(temp)
    }
    Nrow <- min(nMP, sum(apply(used, 2, sum) > 0))
    if (sum(apply(used, 2, sum) > 0) == 0) 
      print(paste("No", Par, "used for these MPs"))
    
    mat <- matrix(1:(Nrow * Ncol), nrow = Nrow, byrow = TRUE)
    op <- par(mfrow = c(Nrow, Ncol), oma = c(3, 6, 2, 0), mar = c(3, 2, 2, 1))
    if (Par == "OM") Title <- "Operating Model Parameters"
    if (Par == "Obs") Title <- "Observation Parameters"
    
    # Colors and Controls
    ncols <- length(topStat)  # nrow(Stat) * ncol(Stat)
    Cols <- colorRampPalette(c("green", "red"))(ncols)
    # rev(rainbow(ncols,start=0,end=0.36))
    highest <- max(Stat, na.rm = TRUE)
    pch <- 18
    LWD <- 3
    LCol <- "black"
    count <- 1
    mm <- 1
    AxCex <- 1.15
    doneMP <- 1
    
    # Make MP colors for Available, Acceptable, and Not-Acceptable
    availCol <- "green"
    acceptCol <- "black"
    nonAAcol <- "darkgray"
    mpcol <- "black"  # default
    mpCols <- data.frame(MPs = MPs, col = rep(mpcol, nMPs), stringsAsFactors = FALSE)
    if (is.null(acceptMP)) 
      acceptMP <- MPs
    if (!is.null(availMP) & (!is.null(acceptMP))) {
      mpCols[, 2] <- nonAAcol
      mpCols[MPs %in% acceptMP, 2] <- acceptCol
      mpCols[MPs %in% acceptMP & MPs %in% availMP, 2] <- availCol
    }
    
    AllMPs <- 1:nMPs
    AllMPs <- AllMPs[apply(used, 2, sum) > 0]  # only include MPs which use the parameter 
    AllMPs <- AllMPs[1:Nrow]  # first nMPs 
    
    # find max y 
    maxY <- 1
    for (mm in AllMPs) {
      for (vr in 1:Ncol) {
        if (nMPs > 1)  varind <- topStat[vr, mm]  # Variable index 
        if (nMPs == 1) varind <- topStat[vr]
        if (nMPs > 1) {
          ys <- senseDat[[mm]]$OMPoints[[varind]][, 2]
        } else {
          ys <- senseDat$OMPoints[[varind]][, 2]
        }
      }
      maxY <- max(maxY, max(ys))
    }
    YLim <- c(0, maxY)
    for (mm in AllMPs) {
      # Loop along MPs Loop along variables
      for (vr in 1:Ncol) {
        if (nMPs > 1)  varind <- topStat[vr, mm]  # Variable index 
        if (nMPs == 1) varind <- topStat[vr]
        varSN <- varNames[varind]
        varLN <- LnName[match(varSN, Obsnam)]
        if (nMPs > 1) {
          xs <- senseDat[[mm]]$OMPoints[[varind]][, 1]
          ys <- senseDat[[mm]]$OMPoints[[varind]][, 2]
        } else {
          xs <- senseDat$OMPoints[[varind]][, 1]
          ys <- senseDat$OMPoints[[varind]][, 2]
        }
        if (used[varSN, MPs[mm]]) {
          # variable is used
          Col <- makeTransparent(Cols[ceiling(Stat[varind, MPs[mm]]/highest * ncols)])
          # ylim <- c(0, quantile(ys, 0.95, na.rm = TRUE))
          plot(xs, ys, col = Col, pch = pch, bty = "n", axes = FALSE, 
               xlab = "", ylab = "", ylim = YLim)
          if (vr == 1) {
            MyCol <- mpCols[match(MPs[mm], mpCols[, 1]), 2]
            axis(side = 2, las = 1, cex.axis = AxCex)
            mtext(side = 2, MPs[mm], line = 2.75, cex = 1.2, col = MyCol)
          }
          if (vr != 1) 
            axis(side = 2, labels = FALSE)
          axis(side = 1, cex.axis = AxCex)
          if (incStat) 
            text(max(xs), 0.05 * max(ys), round(Stat[varind, MPs[mm]], 
                                                2), pos = 2)
          # Smoother line
          if (nMPs > 1) {
            smX <- senseDat[[mm]]$OMSmooth[[varind]]$x
            smY <- senseDat[[mm]]$OMSmooth[[varind]]$y
          } else {
            smX <- senseDat$OMSmooth[[varind]]$x
            smY <- senseDat$OMSmooth[[varind]]$y
          }
          lines(smX, smY, lwd = LWD, col = LCol)
          # Variable Name
          if (!incNames) 
            mtext(side = 1, varSN, cex = 1, line = 2.5)
          if (incNames) 
            mtext(side = 1, varLN, cex = labcex, line = 2.5)
          
        } else {
          plot(c(0, 1), axes = FALSE, type = "n", xlab = "", ylab = "")
        }
      }
    }
    mtext(side = 3, outer = TRUE, Title, cex = 1.5)
    if (YVar == "Y") 
      mtext(side = 2, outer = TRUE, "Long-term yield relative to MSY (%)", 
            cex = 1.25, line = 3)
    if (YVar == "B") 
      mtext(side = 2, outer = TRUE, "B/BMSY in last 5 years", cex = 1.25, 
            line = 3)
  }
  par(op)
  # invisible(Out)
  
}




#' Biomass wormplot
#' 
#' A worm plot for plotting the likelihood of meeting biomass targets in future
#' years.
#' 
#' Returns a matrix of nMPs rows and proyears columns which is the fraction of
#' simulations for which biomass was above Bref.
#' 
#' @param MSEobj Object of class MSE, output of the runMSE function
#' @param Bref The reference fraction of BMSY (to evaluate the probability of
#' exceeding this level)
#' @param LB The lower bound probability that seperates red (bad) and yellow
#' (O.K.) colored segments
#' @param UB The upper bound probability that seperates yellow (O.K.) and green
#' (good) colored segments
#' @author T. Carruthers
#' @export wormplot
wormplot <- function(MSEobj, Bref = 0.5, LB = 0.25, UB = 0.75) {
  
  if (UB < LB) 
    stop("LB parameter must be lower than UB parameter")
  if (LB < 0 | LB > 1) 
    stop("LB parameter must be in the range of 0 to 1")
  if (UB < 0 | UB > 1) 
    stop("UB parameter must be in the range of 0 to 1")
  
  ncol <- ceiling(MSEobj@nMPs^0.3)
  nrow <- ceiling(MSEobj@nMPs/ncol)
  
  op <- par(mfcol = c(nrow, ncol), mar = c(0.1, 0.1, 0.1, 0.1), omi = c(0.6, 
                                                                        0.25, 0.3, 0))
  
  Bprob <- apply(MSEobj@B_BMSY > Bref, 2:3, sum)/MSEobj@nsim
  
  ind <- order(apply(Bprob, 1, sum), decreasing = T)
  
  BLB <- Bprob > LB
  BUB <- Bprob > UB
  
  col <- array("red", dim(Bprob))
  col[BLB & !BUB] = "yellow"
  col[BUB] = "green"
  
  for (i in 1:(nrow * ncol)) {
    if (i < (MSEobj@nMPs + 1)) {
      MP <- ind[i]
      plot(c(1, MSEobj@proyears + 2), c(-1, 1), col = "white", axes = F)
      # abline(h=0)
      
      for (ys in 1:MSEobj@proyears) {
        x <- c(ys - 1, ys, ys, ys - 1)
        y <- c(rep(Bprob[MP, ys], 2), rep(-Bprob[MP, ys], 2))
        pol <- data.frame(x, y)
        polygon(pol, col = col[MP, ys], border = NA)
      }
      
      legend("top", legend = MSEobj@MPs[MP], bty = "n")
      if ((i/nrow) == round(i/nrow, 0)) 
        axis(1, pretty(1:MSEobj@proyears), pretty(1:MSEobj@proyears))
      
      
    } else {
      
      plot.new()
      
    }
    
    if (i == (nrow * ncol)) {
      legend("topright", fill = c("green", "red"), legend = c(paste(">", 
                                                                    round(UB * 100, 0), "% prob.", sep = ""), paste("<", round(LB * 
                                                                                                                                 100, 0), "% prob.", sep = "")), bty = "n")
      
    }
    
  }
  
  mtext(paste("Probability of biomass above ", round(Bref * 100, 0), 
              "% BMSY for ", deparse(substitute(MSE)), sep = ""), 3, outer = T, 
        line = 0.5)
  mtext("Projection year", 1, outer = T, line = 2.5)
  mtext(paste("Fraction of simulations above ", round(Bref * 100, 0), 
              "% BMSY", sep = ""), 2, outer = T, line = 0.25)
  par(op)  
  invisible(Bprob)
  
}

## old functions - may not be used anymore ####
#' Calculate Statistics for MP Performance
#' 
#' Function calculates probabilities and other statistics for a range of
#' performance metrics
#' 
#' @param MSEobj An object of class MSE
#' @param PMRefs A list of reference points for the performance metrics (must
#' be named)
#' @param lastYrs The last number of years in the projection to calculate the
#' statistics
#' @param UseMean Logical. Calculate mean (TRUE) or median (FALSE)?
#' @param msg Logical. Print out messages?
#' @author A. Hordyk
#' @export MPStats
MPStats <- function(MSEobj, PMRefs = list(B_BMSY = 0.5, SSB_SSB0 = 0.2, F_FMSY = 1, 
                                          AAVY = 30, AAVE = 30), lastYrs = 10, UseMean = TRUE, msg = TRUE) {
  
  if (msg) 
    message("Calculating MP Performance for last ", lastYrs, " years")
  flush.console()
  
  sumFun <- ifelse(UseMean, mean, median)
  
  # Assign Variables
  Nyears <- MSEobj@nyears
  Pyears <- MSEobj@proyears
  nMPs <- MSEobj@nMPs
  MPs <- MSEobj@MPs
  nsim <- MSEobj@nsim
  RefYd <- MSEobj@OM$RefY
  
  trefs <- PMRefs
  if (length(names(trefs)) != 5) {
    PMnames <- c("B_BMSY", "SSB_SSB0", "F_FMSY", "AAVY", "AAVE")
    DF <- c(0.5, 0.2, 1, 30, 30)
    ind <- which(!PMnames %in% names(trefs))
    if (length(ind) > 0) {
      for (xx in ind) {
        trefs[[PMnames[xx]]] <- DF[xx]
      }
    }
  }
  
  # Error Checks
  if (!class(lastYrs) == "character" & lastYrs >= Pyears) {
    message("lastYrs set too high. Defaulting to last 10 years")
    lastYrs <- 10
  }
  if (lastYrs <= 0 | lastYrs == FALSE | lastYrs == "all" | lastYrs == 
      "All") 
    lastYrs <- Pyears
  
  yrs <- (Pyears - lastYrs + 1):Pyears  # Years to summarize performance
  yrs <- yrs[!yrs == 1]
  ystart <- 1:lastYrs  # First 10 years
  y1 <- (yrs[1] - 1):(yrs[length(yrs)] - 1)  # for calculating interannual variability
  y2 <- y1 + 1
  
  # Biomass/BMSY
  B_BMSYm <- apply(MSEobj@B_BMSY[, , yrs, drop = FALSE], 2, sumFun, na.rm = TRUE)  # median/mean in last yrs
  B_BMSYsd <- apply(MSEobj@B_BMSY[, , yrs, drop = FALSE], 2, sd, na.rm = TRUE)  # sd in last yrs 
  B_BMSYref <- MSEobj@B_BMSY[, , yrs, drop = FALSE] > trefs$B_BMSY  #  above reference?
  B_BMSYp <- round(apply(B_BMSYref, 2, sum, na.rm = TRUE)/(lastYrs * 
                                                             nsim), 2)  # prob above ref
  
  # Biomass/B0
  temp <- as.matrix(expand.grid(1:nsim, 1:nMPs, 1:Pyears))
  Deplet <- array(NA, dim = dim(MSEobj@B_BMSY))
  Deplet[temp] <- (MSEobj@B_BMSY[temp] * MSEobj@OM$SSBMSY_SSB0[temp[, 1]])
  
  SSB_SSB0ref <- Deplet[, , yrs, drop = FALSE] > trefs$SSB_SSB0  #  above reference?
  SSB_SSB0m <- apply(Deplet[, , yrs, drop = FALSE], 2, sumFun, na.rm = TRUE)  # median/mean in last yrs
  SSB_SSB0sd <- apply(Deplet[, , yrs, drop = FALSE], 2, sd, na.rm = TRUE)  # sd in last yrs 
  SSB_SSB0p <- round(apply(SSB_SSB0ref, 2, sum, na.rm = TRUE)/(lastYrs * nsim), 
                     2)  # prob above ref
  
  # F/FMSY
  F_FMSYm <- apply(MSEobj@F_FMSY[, , yrs, drop = FALSE], 2, sumFun, na.rm = TRUE)  # median/mean in last yrs
  F_FMSYsd <- apply(MSEobj@F_FMSY[, , yrs, drop = FALSE], 2, sd, na.rm = TRUE)  # sd in last yrs 
  F_FMSYref <- MSEobj@F_FMSY[, , yrs, drop = FALSE] < trefs$F_FMSY  #  below reference?
  F_FMSYp <- round(apply(F_FMSYref, 2, sum, na.rm = TRUE)/(lastYrs * nsim), 2)  # prob below ref
  
  # AAVY - Interannual variability in yield
  maxVar <- ifelse(trefs$AAVY > 1, trefs$AAVY/100, trefs$AAVY)
  MSEobj@C[(!is.finite(MSEobj@C[, , , drop = FALSE]))] <- 0  # if catch is NAN or NA, make it 0
  aavy <- (((MSEobj@C[, , y1, drop = FALSE] - MSEobj@C[, , y2, drop = FALSE])/MSEobj@C[, 
                                                                                       , y2, drop = FALSE])^2)^0.5
  
  aavy[is.nan(aavy)] <- NA
  if (sum(aavy == Inf, na.rm=TRUE) > 0) aavy[aavy == Inf] <- NA
  
  # aavy[aavy > 10 * median(aavy, na.rm = TRUE)] <- NA  # remove huge spikes from F=0
  AAVYm <- apply(aavy, 2, sumFun, na.rm = TRUE)  # median/mean in last yrs
  
  AAVYsd <- apply(aavy, 2, sd, na.rm = TRUE)  # sd in last yrs
  AAVYref <- aavy < maxVar
  AAVYp <- round(apply(AAVYref, 2, sum, na.rm = TRUE)/(lastYrs * nsim),  2)
  
  # AAVE - Interannual varialility in effort
  maxVar <- ifelse(trefs$AAVE > 1, trefs$AAVE/100, trefs$AAVE)
  eff <- MSEobj@Effort
  if (all(is.na(eff))) {
    message("Variability in effort unable to be calculated from this version of MSE object")
    aave <- NA
    AAVEm <- rep(NA, length(AAVYm))
    AAVEsd <- rep(NA, length(AAVYsd))
    AAVEref <- array(NA, dim = dim(AAVYref))
    AAVEp <- rep(NA, length(AAVYp))
    flush.console()
  } else {
    aave <- (((eff[, , y1 - 1, drop = FALSE] - eff[, , y2 - 1, drop = FALSE])/eff[, 
                                                                                  , y2 - 1, drop = FALSE])^2)^0.5
    aave[is.nan(aave)] <- NA
    aave[aave == Inf] <- NA
    # aave[aave > 10 * apply(aave, 2, median, na.rm = TRUE)] <- NA  # remove huge spikes from F=0
    AAVEm <- apply(aave, 2, sumFun, na.rm = TRUE)  # median/mean in last yrs
    AAVEsd <- apply(aave, 2, sd, na.rm = TRUE)  # sd in last yrs
    AAVEref <- aave < maxVar
    AAVEp <- round(apply(AAVEref, 2, sum, na.rm = TRUE)/(lastYrs * 
                                                           nsim), 2)
  }
  
  # Yield
  LTY <- round(apply(MSEobj@C[, , yrs, drop = FALSE]/RefYd, 2, sumFun, 
                     na.rm = TRUE) * 100, 1)  # long-term yield
  STY <- round(apply(MSEobj@C[, , ystart, drop = FALSE]/RefYd, 2, sumFun, 
                     na.rm = TRUE) * 100, 1)  # short-term yield
  
  MPtype <- sapply(1:nMPs, function(X) class(get(MPs[X])))
  
  DF <- data.frame(MP = MPs, B_BMSYm = B_BMSYm, B_BMSYsd = B_BMSYsd, 
                   B_BMSYp = B_BMSYp, SSB_SSB0m = SSB_SSB0m, SSB_SSB0sd = SSB_SSB0sd, SSB_SSB0p = SSB_SSB0p, 
                   F_FMSYm = F_FMSYm, F_FMSYsd = F_FMSYsd, F_FMSYp = F_FMSYp, AAVYm = AAVYm, 
                   AAVYsd = AAVYsd, AAVYp = AAVYp, AAVEm = AAVEm, AAVEsd = AAVEsd, 
                   AAVEp = AAVEp, LTY = LTY, STY = STY, lastYrs = lastYrs, B_BMSYRef = trefs$B_BMSY, 
                   SSB_SSB0Ref = trefs$SSB_SSB0, F_FMSYRef = trefs$F_FMSY, AAVYRef = trefs$AAVY, 
                   AAVERef = trefs$AAVE, MPtype = MPtype, stringsAsFactors = FALSE)
  
  Dist <- NULL  # calculate distance from corner
  for (X in 1:nrow(DF)) Dist[X] <- euc.dist(c(DF[X, "B_BMSYp"] * 100, 
                                              DF[X, "LTY"]), c(100, max(DF[, "LTY"], na.rm = TRUE)))
  
  DF$Dist <- Dist
  
  Probs <- list(F_FMSYref = F_FMSYref, SSB_SSB0ref = SSB_SSB0ref, B_BMSYref = B_BMSYref, 
                AAVYref = AAVYref, AAVEref = AAVEref, Effort = eff)
  
  pastC <- apply(MSEobj@CB_hist[, , , , drop = FALSE], c(1, 3), sum, 
                 na.rm = TRUE)/RefYd
  temp <- replicate(nMPs, pastC)
  temp <- aperm(temp, c(1, 3, 2))
  
  lastYr <- temp[, , Nyears, drop = FALSE]
  lastYr2 <- replicate(Pyears, lastYr)
  Yield <- abind::abind(lastYr, MSEobj@C[, , , drop = FALSE]/RefYd, along = 3)
  # dim(Yield)
  
  # totC <- abind(temp, MSEobj@C[,,, drop=FALSE]/RefYd, along=3)
  
  bySim <- list(SSB_SSB0 = Deplet, B_BMSY = MSEobj@B_BMSY, F_FMSY = MSEobj@F_FMSY, 
                AAVY = aavy, AAVE = aave, LTY = MSEobj@C[, , yrs, drop = FALSE]/RefYd, 
                STY = MSEobj@C[, , ystart, drop = FALSE]/RefYd, Yield = Yield)
  
  Out <- list(Perf = DF, BySim = bySim, Probs = Probs)
  Out
}


## Supporting functions - not exported ####


calcStat <- function(rr, evalbreaks) {
  # supporting function
  ind <- as.integer(evalbreaks/2)
  ind2 <- as.integer(0.1 * evalbreaks)
  ind3 <- as.integer(0.9 * evalbreaks)
  if (all(is.na(rr))) return(0)
  if (all(rr$x == 0)) return(0)
  sum((rr$y - mean(rr$y, na.rm = TRUE))^2)
}


calcMSESense <- function(MP = 1, MSEobj, YVar = c("Y", "B"), Par = c("Obs","OM"), 
                         evalbreaks = NULL, quants = c(0.05, 0.95)) {
  # supporting function
  
  YVar <- match.arg(YVar)
  # Calculate for a single MP
  if (length(MP) > 1) stop("Only one MP")
  nMPs <- MSEobj@nMPs
  MPs <- MSEobj@MPs
  if (class(MP) == "character") mm <- which(MPs %in% MP)
  if (class(MP) == "numeric" | class(MP) == "integer") {
    mm <- MP
    MP <- MPs[mm]
  }
  nsims <- MSEobj@nsim
  RefYd <- MSEobj@OM$RefY
  yind <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
  if (is.null(evalbreaks)) evalbreaks <- as.integer(nsims/4)  # number of breaks for loess smoother
  
  if (YVar == "Y") {
    if (length(dim(MSEobj@C)) > 2) {
      yout <- apply(MSEobj@C[, mm, yind], 1, mean, na.rm = T)/RefYd * 100
    } else {
      yout <- apply(MSEobj@C[, yind], 1, mean, na.rm = T)/RefYd * 100
    }
  }
  if (YVar == "B") {
    if (length(dim(MSEobj@B_BMSY)) > 2) {
      yout <- apply(MSEobj@B_BMSY[, mm, yind], 1, mean, na.rm = T)
    } else {
      yout <- apply(MSEobj@B_BMSY[, yind], 1, mean, na.rm = T)
    }
  }
  
  # Operating Model names to include
  MSEobj@OM$Lm_SL <- MSEobj@OM$L50/MSEobj@OM$LFS
  varnames <- names(MSEobj@OM)
  vars <- MSEobj@OM
  vargood <- (apply(vars, 2, sd, na.rm = TRUE)/(apply(vars, 2, mean, na.rm = TRUE)^2)^0.5) > 0.005
  # vargood[grep("qvar", varnames)] <- FALSE
  vargood[is.na(vargood)] <- TRUE
  varnames <- varnames[vargood]
  
  omvals <- MSEobj@OM[, varnames]
  
  # Ignore these parameters from VOI plot
  omvals$RefY <- 0
  omvals$A <- 0
  omvals$Asp <- 0
  omvals$OFLreal <- 0
  omvals$FMSY <- 0
  omvals$MSY <- 0
  omvals$B0 <- 0
  omvals$N0 <- 0
  omvals$SSB0 <- 0
  omvals$BMSY <- 0
  omvals$SSBMSY <- 0
  omvals$LFC <- 0 
  omvals$FMSY_M <- 0 
  omvals$BMSY_B0 <- 0
  omvals$SSBMSY_SSB0 <- 0 
  omvals$MGT <- 0 
  omvals$Blow <- 0
  omvals$FMSY_M <- 0 
  omvals$maxlen <- 0 
  omvals$SRrel <- 0 
  omvals$FinF <- 0 
  omvals[is.na(omvals)] <- 0
  
  OMSmooth <- OMStat <- obvals <- OBSmooth <- OMPoints <- OMNames <- NULL
  if (Par == "OM") {
    # Apply loess smoother to Operating Model parameters
    rr <- as.list(omvals)
    OMNames <- names(rr)
    OMSmooth <- list()
    for (xx in 1:length(rr)) {
      m <- matrix(c(rr[[xx]], yout), ncol = 2)
      m <- m[order(m[, 1]), ]
      qnts <- quantile(m[, 1], quants, na.rm = TRUE)
      if (all(m[, 1] == 0)) {
        OMSmooth[[xx]] <- NA
      } else {
        nas <- which(apply(!apply(m, 2, is.na), 1, prod) == 0)
        if (length(nas) > 0) m <- m[-nas, ]
        ind <- findInterval(qnts, m[, 1])
        m <- m[ind[1]:ind[2], ]
        rr[[xx]] <- m
        OMSmooth[[xx]] <- suppressWarnings(loess.smooth(rr[[xx]][, 1], rr[[xx]][, 2]))
        OMPoints[[xx]] <- m
      }
    }
    # Calculate stat for OM curve
    OMStat <- unlist(lapply(OMSmooth, calcStat, evalbreaks = evalbreaks))
  }
  
  if (Par == "Obs") {
    # Observation Parameters
    varnames <- names(MSEobj@Obs)
    vars <- MSEobj@Obs
    vargood <- (apply(vars, 2, sd, na.rm = TRUE)/(apply(vars, 2, mean, na.rm = TRUE)^2)^0.5) > 0.005
    
    varnames <- varnames[vargood]
    varnames <- varnames[!is.na(varnames)]
    obvals <- MSEobj@Obs[, varnames]
    
    rr <- as.list(obvals)
    OMNames <- names(rr)
    OMSmooth <- list()
    for (xx in 1:length(rr)) {
      m <- matrix(c(rr[[xx]], yout), ncol = 2)
      m <- m[order(m[, 1]), ]
      qnts <- quantile(m[, 1], quants, na.rm = TRUE)
      if (all(m[, 1] == 0)) {
        OMSmooth[[xx]] <- NULL
      } else {
        nas <- which(apply(!apply(m, 2, is.na), 1, prod) == 0)
        if (length(nas) > 0) 
          m <- m[-nas, ]
        
        ind <- findInterval(qnts, m[, 1])
        m <- m[ind[1]:ind[2], ]
        rr[[xx]] <- m
        OMSmooth[[xx]] <- suppressWarnings(loess.smooth(rr[[xx]][, 1], rr[[xx]][, 2]))
        OMPoints[[xx]] <- m
      }
    }
    # Calculate stat for OM curve
    OMStat <- unlist(lapply(OMSmooth, calcStat, evalbreaks = evalbreaks))
  }
  Out <- list(OMPoints = OMPoints, OMSmooth = OMSmooth, OMStat = OMStat, OMNames = OMNames)
  Out
}


custombar<-function(dat,MPnams,tickwd1=0.1,tickwd2=0.05,lwd1=2,lwd2=2,xlab=T,bcol="azure2",barcol='dark grey'){
  
  nMPer<-nrow(dat)
  incr<-(max(dat)-min(dat))*0.05
  plot(dat[,5],ylim=c(min(dat)-incr,max(dat)+incr),xlim=c(0.25,nMPer+0.25),col='white',axes=F,ylab="",xlab="")
  
  if(xlab){
    axis(1,-1:(nMPer+1),c("","",MPnams,""),las=2,font=2,cex.axis=1.2)
  }else{
    axis(1,-1:(nMPer+1),rep("",nMPer+3))
  }
  
  incr<-(max(dat)-min(dat))*0.2
  yp<-pretty(seq(min(dat)-incr,max(dat)+incr,length.out=12))
  axis(2,yp,yp,las=2)
  big<-1E20
  polygon(c(-big,big,big,-big),c(-big,-big,big,big),col=bcol)# grey92
  abline(h=yp,col='white')
  
  for(i in 1:nMPer){
    
    lines(c(i,i),c(dat[i,1],dat[i,5]),lwd=lwd2) # 80%
    lines(c(i,i),c(dat[i,2],dat[i,4]),lwd=lwd1) # 80%
    polygon(c(i-tickwd1/2,i+tickwd1/2,i+tickwd1/2,i-tickwd1/2),c(dat[i,2],dat[i,2],dat[i,4],dat[i,4]),lwd=lwd1,col=barcol) # lower interquartile
    #lines(c(i-tickwd1/2,i+tickwd1/2),c(dat[i,4],dat[i,4]),lwd=lwd1) # upper interquartile
    
    lines(c(i-tickwd2/2,i+tickwd2/2),c(dat[i,1],dat[i,1]),lwd=lwd1) # lower 80%
    lines(c(i-tickwd2/2,i+tickwd2/2),c(dat[i,5],dat[i,5]),lwd=lwd2) # upper 80%
    
    
  }
  
  points(dat[,3],pch=19,cex=1.1)
  
  
}

# A simulation by simulation approach to plotting results
comp <- function(MSEobj, MPs = NA) {
  if (is.na(MPs)) 
    MPs <- MSEobj@MPs
  notm <- MPs[!(MPs %in% MSEobj@MPs)]
  canm <- MPs[MPs %in% MSEobj@MPs]
  if (length(notm) > 0) 
    print(paste("Methods", paste(notm, collapse = ", "), "were not carried out in MSE", 
                deparse(substitute(MSEobj)), sep = " "))
  
  if (length(canm) == 0) 
    stop(paste("None of the methods you specified were carried out in the MSE", 
               deparse(substitute(MSEobj)), sep = ""))
  
  if (length(canm) > 4) {
    print(paste("A maximum of four methods can be compared at once. Plotting first four:", 
                paste(canm[1:4], collapse = ", "), sep = " "))
    canm <- canm[1:4]
  }
  
  mind <- match(canm, MSEobj@MPs)
  nm <- length(mind)
  nsim <- MSEobj@nsim
  proyears <- MSEobj@proyears
  
  Yd <- array(NA, c(nm, nsim))
  P10 <- array(NA, c(nm, nsim))
  P50 <- array(NA, c(nm, nsim))
  P100 <- array(NA, c(nm, nsim))
  POF <- array(NA, c(nm, nsim))
  yind <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
  RefYd <- MSEobj@OM$RefY
  
  for (m in 1:nm) {
    mm <- mind[m]
    Yd[m, ] <- round(apply(MSEobj@C[, mm, yind], 1, mean, na.rm = T)/RefYd * 
                       100, 1)
    POF[m, ] <- round(apply(MSEobj@F_FMSY[, mm, ] > 1, 1, sum, na.rm = T)/proyears * 
                        100, 1)
    P10[m, ] <- round(apply(MSEobj@B_BMSY[, mm, ] < 0.1, 1, sum, na.rm = T)/proyears * 
                        100, 1)
    P50[m, ] <- round(apply(MSEobj@B_BMSY[, mm, ] < 0.5, 1, sum, na.rm = T)/proyears * 
                        100, 1)
    P100[m, ] <- round(apply(MSEobj@B_BMSY[, mm, ] < 1, 1, sum, na.rm = T)/proyears * 
                         100, 1)
  }
  
  MSEcols <- c("red", "green", "blue", "orange")
  
  # dev.new2(width=7,height=7)
  op <- par(mfrow = c(2, 2), mai = c(0.85, 0.7, 0.1, 0.1), omi = rep(0.01, 4))
  
  tradeoffplot2(POF, Yd, "Prob. of overfishing (%)", "Relative yield", 
                vl = 50, hl = 100, coly = MSEcols, leg = NA)
  tradeoffplot2(P100, Yd, "Prob. biomass < BMSY (%)", "Relative yield", 
                vl = 50, hl = 100, coly = MSEcols, leg = canm)
  tradeoffplot2(P50, Yd, "Prob. biomass < 0.5BMSY (%)", "Relative yield", 
                vl = 50, hl = 100, coly = MSEcols, leg = NA)
  tradeoffplot2(P10, Yd, "Prob. biomass < 0.1BMSY (%)", "Relative yield", 
                vl = 50, hl = 100, coly = MSEcols, leg = NA)
  par(op)
}


# Supporting Functions for Plots
tradeoffplot <- function(x, y, xlab, ylab, labs, cex, vl, hl) {
  adjj <- c(0.7, 1.3)
  XLim <- c(min(c(-10, min(x, na.rm = T) * adjj)), max(c(max(x, na.rm = T) * 
                                                           adjj, 110)))
  YLim <- c(min(c(-10, min(y, na.rm = T) * adjj)), max(c(max(y, na.rm = T) * 
                                                           adjj, 110)))
  coly <- rep(c("#0000ff95", "#ff000095", "#20ff1095"), 50)[1:length(labs)]
  coly[labs %in% c("AvC", "curE", "FMSYref")] <- "black"
  # plot(NA,xlim=range(x,na.rm=T)*adjj,ylim=range(y,na.rm=T)*adjj,xlab=xlab,ylab=ylab)
  plot(NA, xlim = XLim, ylim = YLim, xlab = xlab, ylab = ylab)
  abline(v = vl, col = "#99999940", lwd = 2)
  abline(h = hl, col = "#99999940", lwd = 2)
  text(x, y, labs, font = 2, col = coly, cex = 0.9)
}

tradeoffplot2 <- function(x, y, xlab, ylab, cex = 1, vl, hl, coly, leg) {
  adjj <- c(0.7, 1.3)
  plot(NA, xlim = range(x, na.rm = T) * adjj, ylim = range(y, na.rm = T) * 
         adjj, xlab = xlab, ylab = ylab)
  abline(v = vl, col = "grey", lwd = 2)
  abline(h = hl, col = "grey", lwd = 2)
  for (m in 1:nrow(x)) points(x[m, ], y[m, ], col = makeTransparent(coly[m], 
                                                                    50), pch = 19, cex = cex)
  if (!is.na(leg[1])) 
    legend("topright", legend = leg, text.col = coly, bty = "n")
}

tradeoffplot3 <- function(x, y, xlab, ylab, labs, cex, vl, hl, xlim, ylim) {
  coly <- rep(c("#0000ff95", "#ff000095", "#20ff1095"), 10)
  plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
  abline(v = vl, col = "#99999940", lwd = 2)
  abline(h = hl, col = "#99999940", lwd = 2)
  text(x, y, labs, font = 2, col = coly, cex = 1)
}

tradeoffplot4 <- function(x, y, xlab, ylab, labs, cex, vl, hl, ShowLabs = FALSE, 
                          ShowCols = FALSE, AvailMPs = NULL) {
  adjj <- c(0.9, 1.1)
  XLim <- c(min(c(-10, min(x, na.rm = T) * adjj)), max(c(max(x, na.rm = T) * 
                                                           adjj, 110)))
  YLim <- c(min(c(-10, min(y, na.rm = T) * adjj)), max(c(max(y, na.rm = T) * 
                                                           adjj, 110)))
  
  # Which MPs meet minimum PMs
  ind <- which(x >= vl & y >= hl)
  coly <- rep("darkgray", length(labs))
  coly[ind] <- "black"
  # coly[labs%in%c('AvC','curE','FMSYref')]<-'black'
  
  Pch <- rep(21, length(labs))
  # Pch[labs%in%c('AvC','curE','FMSYref')] <- 17
  coly[grep("FMSY", labs)] <- "black"
  Pch[grep("FMSY", labs)] <- 24
  if (!is.null(AvailMPs)) 
    Pch[labs %in% AvailMPs] <- 21
  if (!is.null(AvailMPs)) 
    coly[labs %in% AvailMPs & (x >= vl & y >= hl)] <- "green"
  # coly<-rep(c('#0000ff95','#ff000095','#20ff1095'),50)[1:length(labs)]
  
  plot(NA, xlim = XLim, ylim = YLim, xlab = xlab, ylab = ylab, bty = "l", 
       las = 1)
  abline(v = vl, col = "#99999940", lwd = 2)
  abline(h = hl, col = "#99999940", lwd = 2)
  
  Alpha <- 30
  # polygons
  LeftCol <- rgb(red = 255, green = 0, blue = 0, alpha = Alpha, names = NULL, 
                 maxColorValue = 255)
  RightCol <- rgb(red = 0, green = 255, blue = 0, alpha = Alpha, names = NULL, 
                  maxColorValue = 255)
  
  if (ShowCols) {
    polygon(x = c(0, vl, vl, 0), y = c(0, 0, hl, hl), col = LeftCol, 
            border = NA)
    polygon(x = c(0, vl, vl, 0), y = c(0, 0, 100, 100), col = LeftCol, 
            border = NA)
    polygon(x = c(vl, 100, 100, vl), y = c(0, 0, 100, 100), col = RightCol, 
            border = NA)
    polygon(x = c(vl, 100, 100, vl), y = c(hl, hl, 100, 100), col = RightCol, 
            border = NA)
  }
  
  Cex <- 1.5
  if (!ShowLabs) 
    points(x, y, bg = coly, pch = Pch, cex = Cex, col = "black")
  if (ShowLabs) 
    text(x, y, labs, font = 2, col = coly, cex = 1)
  # if(IdPoints) { message('Click points on plot to display MP name')
  # message('Click Stop to finish') flush.console() identify(x,y,
  # labels=labs) }
  
  labs[ind]
  
}


#' Make colors transparent
#' 
#' 
#' @param someColor Character string describing color
#' @param alpha transparency 
#' @author T. Carruthers
#' @export makeTransparent
makeTransparent <- function(someColor, alpha = 100) {
  newColor <- col2rgb(someColor)
  apply(newColor, 2, function(curcoldata) {
    rgb(red = curcoldata[1], green = curcoldata[2], blue = curcoldata[3], 
        alpha = alpha, maxColorValue = 255)
  })
}


euc.dist <- function(x1, x2) sqrt(sum((x1 - x2)^2))

GetStat <- function(PM, LTY, STY, PNOF, BMSYref, B0ref, VY) {
  switch(PM, `Long-term Yield` = LTY, `Short-term Yield` = STY, Overfishing = PNOF, 
         `Biomass:BMSY` = BMSYref, `Biomass:B0` = B0ref, AnnualVar = VY)
}

StatLab <- function(PM, maxVar, BmsyRef, B0Ref) {
  switch(PM, `Long-term Yield` = "Long-term Yield", `Short-term Yield` = "Short-term Yield", 
         Overfishing = "Prob. of Not Overfishing (%)", `Biomass:BMSY` = paste0("Prob. Biomass >", 
                                                                               BmsyRef, "BMSY (%)"), `Biomass:B0` = paste0("Prob. Biomass >", 
                                                                                                                           B0Ref, "B0 (%)"), AnnualVar = paste0("Prob. AAVY <", maxVar, 
                                                                                                                                                                "%"))
}