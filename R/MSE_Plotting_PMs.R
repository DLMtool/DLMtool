


#' Generic Trade-Plot Function - IN DEVELOPMENT
#'
#' @param MSEobj An object of class `MSE`
#' @param ... Names of Performance Metrics (PMs). First PM is recycled if number of PMs is not even
#' @param lims A numeric vector of acceptable risk/minimum probability thresholds. Recycled if not equal to number of PMs.
#' @param point.size Numeric. Size of the MP points
#' @param lab.size Numeric. Size of MP label
#' @param axis.title.size Numeric. Size of axis titles
#' @param axis.text.size Numeric. Size of axis text
#' @param legend.title.size Numeric. Size of legend title text
#' @param position Character. Position of legend - 'right' or 'bottom'
#' @param fill Character. Color of the fill
#' @param alpha Numeric. Transparency of fill
#' @param PMlist Optional list of PM names. Overrides any supplied in ... above
#' @param Refs An optional named list (matching the PM names) with numeric values to override the default `Ref` values. See examples.
#'
#' @author A. Hordyk
#' @return A summary table of MP performance
#' @examples
#' @export
#'
TradePlot_n <- function(MSEobj, ..., lims=c(0.2, 0.2, 0.8, 0.8), 
                      point.size=2,
                      lab.size=4,
                      axis.title.size=12,
                      axis.text.size=10,
                      legend.title.size=12,
                      position = c("right", "bottom"),
                      fill="gray80",
                      alpha=0.4,
                      PMlist=NULL,
                      Refs=NULL
                      ) {
  if (is.null(PMlist)) {
    PMlist <- unlist(list(...))
  } else {
    PMlist <- unlist(PMlist)
  }
  position <- match.arg(position)
  if(length(PMlist) == 0) PMlist <- c("STY", "LTY", "P10", "AAVY")
  if (class(PMlist) != 'character') stop("Must provide names of PM methods")
  # check
  
  for (X in seq_along(PMlist))
    if (!PMlist[X] %in% avail("PM")) stop(PMlist[X], " is not a valid PM method")
  if (length(PMlist)<2) stop("Must provided more than 1 PM method")
  
  nPMs <- length(PMlist)
  if (nPMs %% 2 != 0) {
    message("Odd number of PMs. Recycling first PM")
    PMlist <- c(PMlist, PMlist[1])
    nPMs <- length(PMlist)
  }
  if (length(lims) < nPMs) {
    message("Recycling limits")
    lims <- rep(lims,10)[1:nPMs]
  }
  if (length(lims) > nPMs) {
    lims <- lims[1:nPMs]
  }
  
  runPM <- vector("list", length(PMlist))
  for (X in 1:length(PMlist)) {
    ref <- Refs[[PMlist[X]]]
    if (is.null(ref)) {
      runPM[[X]] <- eval(call(PMlist[X], MSEobj))  
    } else {
      runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref=ref))
    }
    
  }
  nplots <- nPMs/2
  n.col <- ceiling(sqrt(nplots))
  n.row <- ceiling(nplots/n.col)
  
  m <- matrix(1:(n.col*n.row), ncol=n.col, nrow=n.row, byrow=FALSE)
  xmin <- xmax <- ymin <- ymax <- x <- y <- Class <- label <- fontface <- NULL
  plots <- listout <- list()
  
  xInd <- seq(1, by=2, length.out=nplots)
  yInd <- xInd + 1
  
  for (pp in 1:nplots) {
    yPM <- PMlist[yInd[pp]]
    yvals <- runPM[[match(yPM, PMlist)]]@Mean
    ycap <-  runPM[[match(yPM, PMlist)]]@Caption
    yname <-  runPM[[match(yPM, PMlist)]]@Name
    yline <- lims[match(yPM, PMlist)]
    
    xPM <- PMlist[xInd[pp]]
    xvals <- runPM[[match(xPM, PMlist)]]@Mean
    xcap <-  runPM[[match(xPM, PMlist)]]@Caption
    xname <-  runPM[[match(xPM, PMlist)]]@Name
    xline <- lims[match(xPM, PMlist)]
    
    xlim <- c(0, max(max(xvals, 1)))
    ylim <- c(0, max(max(yvals, 1)))
    
    xrect <- data.frame(xmin=0, xmax=xline, ymin=0, ymax=max(ylim))
    yrect <- data.frame(xmin=0, xmax=max(xlim), ymin=0, ymax=yline)
    
    MPType <- MPtype(MSEobj@MPs)
    Class <- MPType[match(MSEobj@MPs, MPType[,1]),2]
    
    df <- data.frame(x=xvals, y=yvals, label=MSEobj@MPs, Class=Class,
                     pass=xvals>xline & yvals>yline, fontface="plain", xPM=xPM, yPM=yPM)
    df$fontface <- as.character(df$fontface)
    df$fontface[!df$pass] <- "italic"
    df$fontface <- factor(df$fontface)
    listout[[pp]] <- df
    plots[[pp]] <- ggplot2::ggplot() + 
      ggplot2::geom_rect(data=xrect, ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=fill, alpha=alpha) +
      ggplot2::geom_rect(data=yrect, ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=fill, alpha=alpha)
    
    plots[[pp]] <-   plots[[pp]] + 
      ggplot2::geom_point(data=df, ggplot2::aes(x, y, shape=Class, color=Class), size=point.size) +
      ggrepel::geom_text_repel(data=df, ggplot2::aes(x, y, color=Class, label=label, fontface = fontface), show.legend=FALSE, size=lab.size) + 
      ggplot2::xlab(xcap) + ggplot2::ylab(ycap) +
      ggplot2::xlim(xlim) + ggplot2::ylim(ylim) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.title = ggplot2::element_text(size=axis.title.size),
                     axis.text = ggplot2::element_text(size=axis.text.size),
                     legend.text=ggplot2::element_text(size=legend.title.size),
                     legend.title = ggplot2::element_text(size=legend.title.size)) + 
      ggplot2::labs(shape= "MP Type", color="MP Type")
    
  }
  out <- do.call("rbind", listout)
  tab <- table(out$label, out$pass)
  passall <- rownames(tab)[tab[,ncol(tab)] == nplots]
  tt <- summary(MSEobj, PMlist, silent=TRUE, Refs=Refs)
  tt$Satisificed <- FALSE
  tt$Satisificed[match(passall, tt$MP)] <- TRUE
  
  DLMtool:::grid_arrange_shared_legend(plots, n.col, n.row,  position = position)
  tt
  
}


#' @describeIn TradePlot_n A trade-off plot showing probabilities that:
#' \itemize{
#' \item overfishing (POF) against long-term yield is > 50\% of reference yield (LTY)
#' \item spawning biomass is below BMSY (P100) against LTY
#' \item spawning biomass is below 0.5BMSY (P50) against LTY
#' \item spawning biomass is below 0.1BMSY (P10) against LTY
#' }
#' 
#' @export 
Tplot_n <- function(MSEobj, lims=c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5), ...) {
  if (class(lims)!="numeric") stop("Second argument must be numeric")
  TradePlot_n(MSEobj, lims=lims, PMlist=list("POF", "LTY", "P100", "LTY", "P50", "LTY", "P10", "LTY"),  ...)
}

#' @describeIn TradePlot_n A trade-off plot showing probabilities that:
#' \itemize{
#' \item short-term yield is > 50% of reference yield(STY) against long-term yield is > 50\% of reference yield (LTY)
#' \item spawning biomass is below 0.1BMSY (P10) against average annual variability in yield is < 20\% (AAVY)
#' }
#' 
#' @export
Tplot2_n <- function(MSEobj, lims=c(0.2, 0.2, 0.8, 0.8), ...) {
  if (class(lims)!="numeric") stop("Second argument must be numeric")
  TradePlot_n(MSEobj, lims=lims, PMlist=list("STY", "LTY", "P10", "AAVY"), ...)
}

#' @describeIn TradePlot_n A trade-off plot showing probabilities that:
#' \itemize{
#' \item overfishing (POF) against long-term yield is > 50% of reference yield (LTY)
#' \item spawning biomass is below 0.1BMSY (P10) against average annual variability in yield is < 20\% (AAVY)
#' }
#' 
#' @export
Tplot3_n <- function(MSEobj, lims=c(0.5, 0.5, 0.8, 0.5), ...) {
  if (class(lims)!="numeric") stop("Second argument must be numeric")
  TradePlot_n(MSEobj, lims=lims, PMlist=list("POF", "LTY", "P50", "AAVY"), ...)
}

