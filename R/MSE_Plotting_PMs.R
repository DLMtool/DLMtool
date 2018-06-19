


#' Generic Trade-Plot Function - IN DEVELOPMENT
#'
#' @param MSEobj An object of class `MSE`
#' @param ... Names of Performance Metrics (PMs), or other arguments to `TradePlot`. First PM is recycled if number of PMs is not even
#' @param Lims A numeric vector of acceptable risk/minimum probability thresholds. Recycled if not equal to number of PMs.
#' @param Title Optional title for each plot. Character vector of `length(PMs)`/2. Recycled.
#' @param Satisficed Logical. Show only the MPs that meet minimum acceptable thresholds (specified in `Lims`)
#' @param Show Logical. Show the Results table and plots?
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
#' @param Yrs An optional named list (matching the PM names) with numeric values to override the default `Yrs` values. See examples.


#' @author A. Hordyk
#' @return Invisibily returns a list with summary table of MP performance and 
#' the ggplot objects for the plots
#' @export
#'
TradePlot <- function(MSEobj, ..., Lims=c(0.2, 0.2, 0.8, 0.8), 
                      Title=NULL,
                      Satisficed=FALSE,
                      Show=TRUE,
                      point.size=2,
                      lab.size=4,
                      axis.title.size=12,
                      axis.text.size=10,
                      legend.title.size=12,
                      position = c("right", "bottom"),
                      fill="gray80",
                      alpha=0.4,
                      PMlist=NULL,
                      Refs=NULL,
                      Yrs=NULL
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
  if (length(Lims) < nPMs) {
    message("Recycling limits")
    Lims <- rep(Lims,10)[1:nPMs]
  }
  if (length(Lims) > nPMs) {
    Lims <- Lims[1:nPMs]
  }
  
  runPM <- vector("list", length(PMlist))
  for (X in 1:length(PMlist)) {
    ref <- Refs[[PMlist[X]]]
    yrs <- Yrs[[PMlist[X]]]
    if (is.null(ref)) {
      if (is.null(yrs)) {
        runPM[[X]] <- eval(call(PMlist[X], MSEobj))    
      } else {
        runPM[[X]] <- eval(call(PMlist[X], MSEobj, Yrs=yrs))  
      }
      
    } else {
      if (is.null(yrs)) {
        runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref=ref))    
      } else {
        runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref=ref, Yrs=yrs))  
      }
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
  
  if (!(is.null(Title))) Title <- rep(Title, nplots)[1:nplots]
    
  for (pp in 1:nplots) {
    yPM <- PMlist[yInd[pp]]
    yvals <- runPM[[match(yPM, PMlist)]]@Mean
    ycap <-  runPM[[match(yPM, PMlist)]]@Caption
    yname <-  runPM[[match(yPM, PMlist)]]@Name
    yline <- Lims[match(yPM, PMlist)]
    
    xPM <- PMlist[xInd[pp]]
    xvals <- runPM[[match(xPM, PMlist)]]@Mean
    xcap <-  runPM[[match(xPM, PMlist)]]@Caption
    xname <-  runPM[[match(xPM, PMlist)]]@Name
    xline <- Lims[match(xPM, PMlist)]
    
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
    
    if (Satisficed) {
     xlim <- c(xline, 1)
     ylim <- c(yline, 1)
     plots[[pp]] <- ggplot2::ggplot() 
    } else {
      plots[[pp]] <- ggplot2::ggplot() + 
        ggplot2::geom_rect(data=xrect, ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=fill, alpha=alpha) +
        ggplot2::geom_rect(data=yrect, ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=fill, alpha=alpha)
    }
      
    plots[[pp]] <-  plots[[pp]] + 
      ggplot2::geom_point(data=df, ggplot2::aes(x, y, shape=Class, color=Class), size=point.size, na.rm=TRUE) +
      ggrepel::geom_text_repel(data=df, ggplot2::aes(x, y, color=Class, label=label, 
                                                     fontface = fontface), 
                               show.legend=FALSE, size=lab.size, na.rm=TRUE) + 
      ggplot2::xlab(xcap) + ggplot2::ylab(ycap) +
      ggplot2::xlim(xlim) + ggplot2::ylim(ylim) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.title = ggplot2::element_text(size=axis.title.size),
                     axis.text = ggplot2::element_text(size=axis.text.size),
                     legend.text=ggplot2::element_text(size=legend.title.size),
                     legend.title = ggplot2::element_text(size=legend.title.size)) + 
      ggplot2::labs(shape= "MP Type", color="MP Type")
    
    if (!is.null(Title)) 
      plots[[pp]] <-  plots[[pp]] + ggplot2::labs(title=Title[pp])
  }

  out <- do.call("rbind", listout)
  tab <- table(out$label, out$pass)
  passall <- rownames(tab)[tab[,ncol(tab)] == nplots]
  Results <- summary(MSEobj, PMlist, silent=TRUE, Refs=Refs)
  Results$Satisificed <- FALSE
  Results$Satisificed[match(passall, Results$MP)] <- TRUE

  
  if (Show) {
    join_plots(plots, n.col, n.row,  position = position)
    print(Results)  
  } 
  
  out <- list(Results=Results, Plots=plots)
  
  invisible(out)
  
}


#' @describeIn TradePlot A trade-off plot showing probabilities that:
#' \itemize{
#' \item overfishing (POF) against long-term yield is > 50\% of reference yield (LTY)
#' \item spawning biomass is below BMSY (P100) against LTY
#' \item spawning biomass is below 0.5BMSY (P50) against LTY
#' \item spawning biomass is below 0.1BMSY (P10) against LTY
#' }
#' 
#' @export 
Tplot <- function(MSEobj, Lims=c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5), ...) {
  if (class(Lims)!="numeric") stop("Second argument must be numeric")
  TradePlot(MSEobj, Lims=Lims, PMlist=list("POF", "LTY", "P100", "LTY", "P50", "LTY", "P10", "LTY"),  ...)
}

#' @describeIn TradePlot A trade-off plot showing probabilities that:
#' \itemize{
#' \item short-term yield is > 50% of reference yield(STY) against long-term yield is > 50\% of reference yield (LTY)
#' \item spawning biomass is below 0.1BMSY (P10) against average annual variability in yield is < 20\% (AAVY)
#' }
#' 
#' @export
Tplot2 <- function(MSEobj, Lims=c(0.2, 0.2, 0.8, 0.8), ...) {
  if (class(Lims)!="numeric") stop("Second argument must be numeric")
  TradePlot(MSEobj, Lims=Lims, PMlist=list("STY", "LTY", "P10", "AAVY"), ...)
}

#' @describeIn TradePlot A trade-off plot showing probabilities that:
#' \itemize{
#' \item overfishing (POF) against long-term yield is > 50% of reference yield (LTY)
#' \item spawning biomass is below 0.1BMSY (P10) against average annual variability in yield is < 20\% (AAVY)
#' }
#' 
#' @export
Tplot3 <- function(MSEobj, Lims=c(0.5, 0.5, 0.8, 0.5), ...) {
  if (class(Lims)!="numeric") stop("Second argument must be numeric")
  TradePlot(MSEobj, Lims=Lims, PMlist=list("POF", "LTY", "P50", "AAVY"), ...)
}







#' Plot the median biomass and yield relative to last historical year
#' 
#' Compare median biomass and yield in first year and last 5 years of
#' projection
#' 
#' @param MSEobj An object of class MSE 
#' @param MPs Optional vector of MPs to plot
#' @param lastYrs Numeric. Last number of years to summarize results.
#' @param point.size Size of the points
#' @param lab.size Size of labels
#' @param axis.title.size Axis title size
#' @param axis.text.size  Axis text size
#' @param legend.title.size Legend title size
#'
#' @export
#' 
#' @examples 
#' \dontrun{
#' MSE <- runMSE()
#' Cplot(MSE)
#' }
Cplot <- function(MSEobj, MPs = NA, lastYrs = 5,
                  point.size=2,
                  lab.size=4,
                  axis.title.size=12,
                  axis.text.size=10,
                  legend.title.size=12) {
  if (!all(is.na(MPs))) MSEobj <- Sub(MSEobj, MPs = MPs)
  mp <- Catch <- Biomass <- mB <- mC <- NULL # cran check hacks 
  nsim <- MSEobj@nsim
  nMPs <- MSEobj@nMPs
  MPs <- MSEobj@MPs
  nyears <- MSEobj@nyears
  proyears <- MSEobj@proyears
  RefYd <- MSEobj@OM$RefY
  
  MPType <- MPtype(MSEobj@MPs)
  Class <- MPType[match(MSEobj@MPs, MPType[,1]),2]
  
  pastC <- apply(MSEobj@CB_hist[, , , , drop = FALSE], c(1, 3), sum, na.rm = TRUE)/RefYd # relative catch in last historical year
  temp <- aperm( replicate(nMPs, pastC), c(1, 3, 2))
  
  lastYr <- temp[, , nyears, drop = FALSE]
  Yield <- abind::abind(lastYr, MSEobj@C[, , , drop = FALSE]/RefYd, along = 3) # 
  
  ny <- MSEobj@nyears + 1
  relYield <- Yield[, , , drop = FALSE]/Yield[, , rep(1, ny), drop = FALSE] # catch relative to last historical year
  relYield <- relYield[,,(proyears - lastYrs + 1):proyears]
  
  # Biomass 
  bio <- MSEobj@B_BMSY[,,(proyears - lastYrs + 1):proyears] # biomass in lastyrs
  histSSB <- apply(MSEobj@SSB_hist[, , , , drop = FALSE], c(1, 3), sum, na.rm = TRUE)
  relSSB <- histSSB[,nyears]/ MSEobj@OM$SSBMSY # SSB/SSBmsy in last historical year
  temp <- array(replicate(nMPs, relSSB), dim=dim(bio)) # array with biomass in last projection year
  
  relbio <- bio/temp # biomass relative to biomass at start of projections
  
  dimnames(relbio) <- list(1:nsim, MPs, 1:lastYrs)
  Bdf <- as.data.frame.table(relbio)
  names(Bdf) <- c("sim", "mp", "yr", 'Biomass')

  dimnames(relYield) <- list(1:nsim, MPs, 1:lastYrs)
  Cdf <- as.data.frame.table(relYield)
  names(Cdf) <- c("sim", "mp", "yr", 'Catch')
  
  DF <- dplyr::left_join(Cdf, Bdf,by = c("sim", "mp", "yr"))
  DF <- DF %>% dplyr::group_by(mp) %>% dplyr::summarize(mC=median(Catch),
                                    upC=quantile(Catch, 0.95),
                                    lowC=quantile(Catch, 0.05),
                                    mB=median(Biomass),
                                    upB=quantile(Biomass, 0.95),
                                    lowB=quantile(Biomass, 0.05))
  DF$class <- Class
  
  
  p1 <- ggplot2::ggplot(DF, ggplot2::aes(x=mB, y=mC, color=class, shape=class)) + 
    ggplot2::geom_point() +
    ggplot2::expand_limits(x=0, y=0) + 
    ggplot2::geom_vline(xintercept = 1, color="gray") +
    ggplot2::geom_hline(yintercept = 1, color="gray") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title = ggplot2::element_text(size=axis.title.size),
                   axis.text = ggplot2::element_text(size=axis.text.size),
                   legend.text=ggplot2::element_text(size=legend.title.size),
                   legend.title = ggplot2::element_text(size=legend.title.size)) +
    ggrepel::geom_text_repel(ggplot2::aes(label=mp), show.legend=FALSE) +
    ggplot2::labs(x=paste("Median Biomass (last", lastYrs, 
                 "years)\n relative to current"),
         y=paste("Median Yield (last", lastYrs, "years)\n relative to current"),
         shape= "MP Type", color="MP Type")

  print(p1)
  
}



