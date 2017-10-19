

# modified from
# https://github.com/tidyverse/ggplot2/wiki/share-a-legend-between-two-ggplot2-graphs
grid_arrange_shared_legend <- function(plots, ncol = length(plots), nrow = 1, position = c("bottom", "right")) {
  

  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
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
#' @export
#'
#' @examples 
#' \dontrun{
#'  testplot{myMSE}
#' }
testplot <- function(MSEobj, ..., lims=c(0.2, 0.2, 0.8, 0.8)) {
  PMlist <- unlist(list(...))
  if(length(PMlist) == 0) PMlist <- c("LTY", "STY", "P50", "AAVY")
  if (class(PMlist) != 'character') stop("Must provide names of PM methods")
  # check
  
  for (X in seq_along(PMlist))
    if (!PMlist[X] %in% avail("PM")) stop(PMlist[X], " is not a valid PM method")
  if (length(PMlist)<2) stop("Must provided more than 1 PM method")
  
  # PMlist <<- PMlist
  runPM <- vector("list", length(PMlist))
  for (X in 1:length(PMlist)) runPM[[X]] <- eval(call(PMlist[X], MSEobj))

  PlotList <- combn(unique(PMlist), 2 )
  lims <- rep(lims, 100)[1:length(PMlist)]

  n.col <- ceiling(sqrt(ncol(PlotList)))
  n.row <- ceiling(ncol(PlotList)/n.col)

  m <- matrix(1:(n.col*n.row), ncol=n.col, nrow=n.row, byrow=FALSE)

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
    
    df <- data.frame(x=xvals, y=yvals, label=MSEobj@MPs, Class=MPclass(MSEobj@MPs),
                    pass=xvals>xline & yvals>yline, fontface="plain", xPM=xPM, yPM=yPM)
    df$fontface <- as.character(df$fontface)
    df$fontface[!df$pass] <- "italic"
    df$fontface <- factor(df$fontface)
    listout[[pp]] <- df
    plots[[pp]] <- ggplot() + 
      geom_rect(data=xrect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.4) +
      geom_rect(data=yrect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.4)
    
  
    
    
    # plots[[pp]] <- ggplot(df, aes(x, y, shape=Class, color=Class, label=label)) +
    plots[[pp]] <-   plots[[pp]] + 
      geom_point(data=df, aes(x, y, shape=Class, color=Class), size=2) +
      geom_text_repel(data=df, aes(x, y, color=Class, label=label, fontface = fontface), show.legend=FALSE) + 
      xlab(xcap) + ylab(ycap) +
      xlim(xlim) + ylim(ylim) +
      theme_classic() +
      theme(axis.title.x = element_text(size=11),
            axis.title.y = element_text(size=11),
            legend.text=element_text(size=12)) + 
      labs(shape= "MP Class", color="MP Class")
    
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
