
# modified from
# https://github.com/tidyverse/ggplot2/wiki/share-a-legend-between-two-ggplot2-graphs
grid_arrange_shared_legend <- function(plots, ncol = length(plots), nrow = 1, position = c("bottom", "right")) {
  

  position <- match.arg(position)
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
  
  grid::grid.newpage()
  grid::grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

#' Test Trade-Off Plot
#'
#' @param MSEobj An object of class MSE
#' @param ... Names of PM methods to plot
#' @param lims Numeric vector of satisficing limits. Recycled to number of PM methods
#' @return produces a plot 
#' @author A. Hordyk
#' @importFrom ggplot2 ggplot aes geom_rect geom_point xlim ylim xlab ylab theme theme_classic labs ggplotGrob
#' @importFrom ggrepel geom_text_repel
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid unit.c unit grid.newpage grid.draw
#' @importFrom utils combn
#' @export
#'
#' @examples 
#' \dontrun{
#'  Tplot3{myMSE}
#' }
Tplot3 <- function(MSEobj, ..., lims=c(0.2, 0.2, 0.8, 0.8)) {
  PMlist <- unlist(list(...))
  if(length(PMlist) == 0) PMlist <- c("LTY", "STY", "P50", "AAVY")
  if (class(PMlist) != 'character') stop("Must provide names of PM methods")
  # check
  
  for (X in seq_along(PMlist))
    if (!PMlist[X] %in% avail("PM")) stop(PMlist[X], " is not a valid PM method")
  if (length(PMlist)<2) stop("Must provided more than 1 PM method")
  
  runPM <- vector("list", length(PMlist))
  for (X in 1:length(PMlist)) runPM[[X]] <- eval(call(PMlist[X], MSEobj))

  PlotList <- combn(unique(PMlist), 2)
  lims <- rep(lims, 100)[1:length(PMlist)]

  n.col <- ceiling(sqrt(ncol(PlotList)))
  n.row <- ceiling(ncol(PlotList)/n.col)

  m <- matrix(1:(n.col*n.row), ncol=n.col, nrow=n.row, byrow=FALSE)
  xmin <- xmax <- ymin <- ymax <- x <- y <- Class <- label <- fontface <- NULL
  plots <- listout <- list()
  for (pp in 1:ncol(PlotList)) {
    yPM <- PlotList[1,pp]
    yvals <- runPM[[match(yPM, PMlist)]]@Mean
    ycap <-  runPM[[match(yPM, PMlist)]]@caption
    yname <-  runPM[[match(yPM, PMlist)]]@name
    yline <- lims[match(yPM, PMlist)]
    
    xPM <- PlotList[2,pp]
    xvals <- runPM[[match(xPM, PMlist)]]@Mean
    xcap <-  runPM[[match(xPM, PMlist)]]@caption
    xname <-  runPM[[match(xPM, PMlist)]]@name
    xline <- lims[match(xPM, PMlist)]
    
    xlim <- c(0, max(max(xvals, 1)))
    ylim <- c(0, max(max(yvals, 1)))
    
    xrect <- data.frame(xmin=0, xmax=xline, ymin=0, ymax=max(ylim))
    yrect <- data.frame(xmin=0, xmax=max(xlim), ymin=0, ymax=yline)
    
    df <- data.frame(x=xvals, y=yvals, label=MSEobj@MPs, Class=MPtype(MSEobj@MPs)[,2],
                    pass=xvals>xline & yvals>yline, fontface="plain", xPM=xPM, yPM=yPM)
    df$fontface <- as.character(df$fontface)
    df$fontface[!df$pass] <- "italic"
    df$fontface <- factor(df$fontface)
    listout[[pp]] <- df
    plots[[pp]] <- ggplot2::ggplot() + 
      ggplot2::geom_rect(data=xrect, ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.4) +
      ggplot2::geom_rect(data=yrect, ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.4)
    
  
    
    
    # plots[[pp]] <- ggplot(df, aes(x, y, shape=Class, color=Class, label=label)) +
    plots[[pp]] <-   plots[[pp]] + 
      ggplot2::geom_point(data=df, ggplot2::aes(x, y, shape=Class, color=Class), size=2) +
      ggrepel::geom_text_repel(data=df, ggplot2::aes(x, y, color=Class, label=label, fontface = fontface), show.legend=FALSE) + 
      ggplot2::xlab(xcap) + ggplot2::ylab(ycap) +
      ggplot2::xlim(xlim) + ggplot2::ylim(ylim) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.title.x = ggplot2::element_text(size=11),
            axis.title.y = ggplot2::element_text(size=11),
            legend.text=ggplot2::element_text(size=12)) + 
      ggplot2::labs(shape= "MP Class", color="MP Class")
    
  }
  out <- do.call("rbind", listout)
  tab <- table(out$label, out$pass)
  passall <- rownames(tab)[tab[,ncol(tab)] == ncol(PlotList)]
  tt <- summary(MSEobj, PMlist, silent=TRUE)
  tt$Satisificed <- FALSE
  tt$Satisificed[match(passall, tt$MP)] <- TRUE
  
  grid_arrange_shared_legend(plots, n.col, n.row)
  tt
  
}
