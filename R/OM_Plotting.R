
render_plot <- function(Object, Class, Stock=NULL, RMD=NULL, nsamp=3, nsim=200, nyears=50, 
                        proyears=28, output_file=NULL, output_dir=getwd(), 
                        quiet=TRUE, tabs=TRUE, title=NULL, date=NULL,
                        plotPars =NULL, open=TRUE, dev=FALSE, parallel=TRUE) {
  
  SampCpars <- list() # empty list
  
  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("Package \"knitr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package \"rmarkdown\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if (is.null(plotPars)) plotPars <- list(breaks=10, col="darkgray", axes=FALSE, 
                                          cex.main=1, lwd=2)
 
  if (class(Object) == "OM") {
    nsim <- Object@nsim 
    nyears <- Object@nyears
    proyears <- Object@proyears
    SampCpars <- if(length(Object@cpars)>0) SampCpars <- SampleCpars(Object@cpars, nsim, msg=FALSE)
    set.seed(Object@seed)
    Stock <- SubOM(Object, "Stock")
    # Class <- "OM"
  }
  
 
  if (Class == "Stock") {
    if (is.null(title)) title <- "Stock Object Plots"
    Pars <- SampleStockPars(Object, nsim, nyears, proyears, SampCpars, 
                            msg=FALSE)
    Pars$Name <- gsub(" ", "_", Object@Name)  
  } else if (Class == "Fleet") {
    if (is.null(title)) title <- "Fleet Object Plots"
    if (class(Stock)!="Stock") 
      stop("Must provide object of class 'Stock'", call. = FALSE)
    StockPars <- SampleStockPars(Stock, nsim, nyears, proyears, SampCpars, 
                                 msg=FALSE)
    FleetPars <- SampleFleetPars(Object, StockPars, nsim, nyears, proyears,
                                 SampCpars, msg=FALSE)
  
    Pars <- c(StockPars, FleetPars)
    Pars$Name <- gsub(" ", "_", Object@Name)
    Pars$CurrentYr <- Object@CurrentYr
    Pars$MPA <- Object@MPA
    
  } else if (Class == "Obs") {
    if (is.null(title)) title <- "Obs Object Plots"
    ObsPars <- SampleObsPars(Object, nsim, cpars=SampCpars)
    BMSY_B0bias <- array(rlnorm(nsim, 
                                mconv(1, Object@BMSY_B0biascv), sdconv(1, Object@BMSY_B0biascv)), 
                         dim = c(nsim))  # trial samples of BMSY relative to unfished  

    ObsPars$BMSY_B0bias <- BMSY_B0bias
    
    Pars <- c(ObsPars)
    
  } else if (Class == "Imp") {
    if (is.null(title)) title <- "Imp Object Plots"
    ImpPars <- SampleImpPars(Object, nsim, cpars=SampCpars)
    Pars <- c(ImpPars)
  } else if (Class == "OM") {
    if (is.null(title)) title <- "OM Object Plots"
    message("Sampling Stock, Fleet, Obs, and Imp parameters")
    StockPars <- SampleStockPars(SubOM(Object, "Stock"), nsim, nyears, proyears, SampCpars, msg=FALSE)
    FleetPars <- SampleFleetPars(SubOM(Object, "Fleet"), StockPars, nsim, nyears, proyears, SampCpars, msg=FALSE)
    ObsPars <- SampleObsPars(Object, nsim, cpars=SampCpars)
    BMSY_B0bias <- array(rlnorm(nsim, 
                                mconv(1, Object@BMSY_B0biascv), sdconv(1, Object@BMSY_B0biascv)), 
                         dim = c(nsim))  # trial samples of BMSY relative to unfished  
    
    ObsPars$BMSY_B0bias <- BMSY_B0bias
    ImpPars <- SampleImpPars(SubOM(Object, "Imp"), nsim, cpars=SampCpars)

    Pars <- c(StockPars, FleetPars, ObsPars, ImpPars)
    Pars$CurrentYr <- Object@CurrentYr
    
    if (!parallel) dopar <- FALSE
    if (nsim>=48 & parallel) dopar <- TRUE
    if (nsim<48& parallel) dopar <- FALSE
    message("Running Historical Simulations")
    Hist <- runMSE(Object, Hist=TRUE, silent=TRUE, parallel = dopar)
    Pars$Hist <- Hist
    Pars$Name <- "OM"
    Pars$MPA <- Object@MPA
    Pars$CurrentYr <- Object@CurrentYr
  } else if (Class == "Hist") {
    Pars <- list()
    Pars$Hist <- Object
    Pars$CurrentYr <- Object@Misc$CurrentYr
    nyears <- length(Object@Data@Year)
    if (is.null(title)) title <- "Historical Simulations"

  } else {
    stop("Object must be class 'Stock', 'Fleet', 'Obs', or 'Imp'", call.=FALSE)  
  }
  

  if (Class !="Hist" & Class !="OM") {
    Pars$Name <- gsub(" ", "_", Object@Name)  
  } 
  its <- sample(1:nsim, nsamp)
  # Pars <<- Pars
  Params <- list(
    title = title,
    Pars = Pars,
    plotPars=plotPars,
    tabs = tabs,
    its = its,
    nyears=nyears,
    proyears=proyears,
    date=NULL
  )

  outname <- paste0("_", RMD, ".html")
  if (Class !="Hist" & Class !="OM") {
    if (is.null(output_file)) output_file <- paste0(Pars$Name, outname)
  } else {
    if (is.null(output_file)) output_file <-  paste0(RMD, ".html")
  }
  message("Rendering HTML file")
  
  RMD <- paste0(RMD, ".Rmd")
  if (dev) {
    input <- file.path('inst/Rmd', Class, RMD) 
  } else {
    input <- file.path(system.file(package = 'DLMtool'),'Rmd', Class, RMD)  
  }
 
  knitr::knit_meta(class=NULL, clean = TRUE)
  rend <- try(rmarkdown::render(input, params=Params,
                                output_file=output_file,
                                output_dir=output_dir,
                                quiet=quiet), silent=TRUE)
  if (class(rend) == "try-error") {
    print(rend)
  } else {
    message("Rendered ", output_file, " in ", output_dir)
    if (open) utils::browseURL(file.path(output_dir, output_file))
  }
}


#' @method plot character
#' @export
#' @keywords internal
plot.character <- function(x, Object, ...) {
  plot.pars(x, Object, ...)  
}

#' @param Object An object of class `Stock` or `Fleet`
#' @param Stock An object of class `Stock` required for `Fleet` parameters
#' @rdname plot.Stock
#' @export
plot.pars <- function(x, Object, Stock=NULL, nsamp=3, nsim=200, nyears=50, 
                      proyears=28, output_file=NULL, output_dir=getwd(), 
                      quiet=TRUE, tabs=TRUE, title=NULL, date=NULL,
                      plotPars =NULL, open=TRUE, dev=FALSE, ...) {
  
  StockDF <- data.frame(chr=c("M",
                              "Growth",
                              "Maturity",
                              "Recruitment",
                              "Spatial",
                              "Depletion"),
                        Class="Stock",
                        stringsAsFactors = FALSE)
  
 FleetDF <- data.frame(chr=c("Effort",
                              "Catchability",
                              "MPA",
                              "Selectivity"),
                        Class="Fleet",
                        stringsAsFactors = FALSE)
  
 DF <- dplyr::bind_rows(StockDF, FleetDF)
  
  if (!x %in% DF[,1]) 
    stop("Invalid argument. Valid arguments are: ", paste0(DF[,1], sep=" "), call.=FALSE)
  
  Class <- DF$Class[match(x, DF[,1])]
  if (class(Object) !="OM" & class(Object) != Class) 
    stop("Incorrect class object for this parameter", call.=FALSE)
  
  
  if (x == "M") x <- "NaturalMortality"
  
  render_plot(Object=Object, Class=Class, Stock=Stock, RMD=x, nsamp=nsamp, nsim=nsim, 
              nyears=nyears, proyears=proyears,
              output_file=output_file, output_dir=output_dir, quiet=quiet,
              tabs=tabs, title=title, date=date,
              plotPars=plotPars, open=open, dev=dev)
}  



#' @title Plot Operating Model Object
#' 
#' @description Generate HTML reports with plots of operating model components ("Stock",
#' "Fleet", "Obs", and "Imp"), the historical simulations ("Hist"), or the complete OM ("OM").
#' 
#' The individual component plots of objects of class `Stock` and `Fleet` can also be generated by 
#' using the generic `plot.pars` function. See Examples below.
#' 
#' @param x An object of class `Stock`, `Fleet`, `Obs`, `Imp`, `Hist`, or `OM`, OR one 
#' of the following character strings for `Object` of class `Stock`: "M", "Growth", "Maturity", "Recruitment", "Spatial",
#' or "Depletion" and for `Object` of class `Fleet`: "Effort", "Catchability", "MPA",
#' and "Selectivity".
#' @param nsamp The number of random samples to show in the plot
#' @param nsim The number of simulations (only used for objects not of class `OM`)
#' @param nyears The number of historical years (only used for objects not of class `OM`)
#' @param proyears The number of projection years (only used for objects not of class `OM`)
#' @param output_file Name of the output html file (without file extension)
#' @param output_dir Output directory. Defaults to `getwd()`
#' @param quiet An option to suppress printing of the pandoc command line 
#' @param tabs Include tabs in the HTML file?
#' @param title Optional title for the markdown report
#' @param date Optional date for the markdown report
#' @param plotPars A named list with options for plots:
#' \itemize{
#'   \item breaks - numeric. Number of breaks in histograms.
#'   \item col - character. Color of histograms.
#'   \item axes - logical. Include axes in histogram?
#'   \item cex.main - numeric. Size of main title in plots.
#'   \item lwd - numeric. Line width for time-series plots.
#' }
#' @param open Logical. Open the html file?
#' @param dev Logical. For development use only.
#' @param ... Not used
#'
#' @method plot Stock
#' @export
#' @examples
#' \dontrun{
#' # Plot Stock Object:
#' Stock <- DLMtool::Albacore
#' plot(Stock)
#' 
#' # Individual plots:
#' plot("M", Stock)
#' plot("Growth", Stock)
#' plot("Maturity", Stock)
#' plot("Recruitment", Stock)
#' plot("Spatial", Stock)
#' plot("Depletion", Stock)
#' 
#' # Plot Fleet Object
#' Fleet <- DLMtool::Generic_DecE
#' plot(Fleet, Stock)
#' 
#' # Individual plots:
#' plot("Effort", Fleet, Stock)
#' plot("Catchability", Fleet, Stock)
#' plot("MPA", Fleet, Stock)
#' plot("Selectivity", Fleet, Stock)
#' 
#' 
#' # Plot Obs Object
#' Obs <- DLMtool::Imprecise_Unbiased
#' plot(Obs)
#' 
#' # Plot Imp Object
#' Imp <- DLMtool::Overages
#' plot(Imp)
#' 
#' 
#' # Plot Hist Object
#' OM <- DLMtool::testOM 
#' Hist <- runMSE(OM, Hist=TRUE)
#' plot(Hist)
#' 
#' # Plot OM Object
#' plot(OM)
#' }
plot.Stock <- function(x, nsamp=3, nsim=200, nyears=50, 
                       proyears=28, output_file=NULL, output_dir=getwd(), 
                       quiet=TRUE, tabs=TRUE, title=NULL, date=NULL,
                       plotPars =NULL, open=TRUE, dev=FALSE, ...){
  
  render_plot(Object=x, Class="Stock", RMD='Stock', nsamp=nsamp, nsim=nsim, 
              nyears=nyears, proyears=proyears,
              output_file=output_file, output_dir=output_dir, quiet=quiet,
              tabs=tabs, title=title, date=date,
              plotPars=plotPars, open=open, dev=dev)
}


#' @rdname plot.Stock
#' @method plot Fleet
#' @export
plot.Fleet <- function(x, Stock=NULL, nsamp=3, nsim=200, nyears=50, 
                       proyears=28, output_file=NULL, output_dir=getwd(), 
                       quiet=TRUE, tabs=TRUE, title=NULL, date=NULL,
                       plotPars =NULL, open=TRUE, dev=FALSE, ...){
  if (class(Stock) !="Stock" & class(x) !="OM")
    stop("Must provide object of class 'Stock'")
  
  render_plot(Object=x, Class="Fleet", Stock=Stock, RMD='Fleet', nsamp=nsamp, nsim=nsim, 
              nyears=nyears, proyears=proyears,
              output_file=output_file, output_dir=output_dir, quiet=quiet,
              tabs=tabs, title=title, date=date,
              plotPars=plotPars, open=open, dev=dev)
}

#' @rdname plot.Stock
#' @method plot Obs
#' @export
plot.Obs <- function(x, nsamp=3, nsim=200, nyears=50, 
                       proyears=28, output_file=NULL, output_dir=getwd(), 
                       quiet=TRUE, tabs=TRUE, title=NULL, date=NULL,
                       plotPars =NULL, open=TRUE, dev=FALSE, ...){
  
  render_plot(Object=x, Class="Obs", Stock=NULL, RMD='Obs', nsamp=nsamp, nsim=nsim, 
              nyears=nyears, proyears=proyears,
              output_file=output_file, output_dir=output_dir, quiet=quiet,
              tabs=tabs, title=title, date=date,
              plotPars=plotPars, open=open, dev=dev)
}

#' @rdname plot.Stock
#' @method plot Imp
#' @export
plot.Imp <- function(x, nsamp=3, nsim=200, nyears=50, 
                     proyears=28, output_file=NULL, output_dir=getwd(), 
                     quiet=TRUE, tabs=TRUE, title=NULL, date=NULL,
                     plotPars =NULL, open=TRUE, dev=FALSE, ...){
  
  render_plot(Object=x, Class="Imp", Stock=NULL, RMD='Imp', nsamp=nsamp, nsim=nsim, 
              nyears=nyears, proyears=proyears,
              output_file=output_file, output_dir=output_dir, quiet=quiet,
              tabs=tabs, title=title, date=date,
              plotPars=plotPars, open=open, dev=dev)
}

#' @rdname plot.Stock
#' @method plot Hist
#' @export
plot.Hist <- function(x, nsamp=3, nsim=200, nyears=50, 
                      proyears=28, output_file=NULL, output_dir=getwd(), 
                      quiet=TRUE, tabs=TRUE, title=NULL, date=NULL,
                      plotPars =NULL, open=TRUE, dev=FALSE, ...) {
  render_plot(Object=x, Class="Hist", Stock=NULL, RMD='Hist', nsamp=nsamp, nsim=nsim, 
              nyears=nyears, proyears=proyears,
              output_file=output_file, output_dir=output_dir, quiet=quiet,
              tabs=tabs, title=title, date=date,
              plotPars=plotPars, open=open, dev=dev)
}

#' @rdname plot.Stock
#' @method plot OM
#' @export
plot.OM <- function(x, nsamp=3, nsim=200, nyears=50, 
                    proyears=28, output_file=NULL, output_dir=getwd(), 
                    quiet=TRUE, tabs=TRUE, title=NULL, date=NULL,
                    plotPars =NULL, open=TRUE, dev=FALSE, ...) {
  render_plot(Object=x, Class="OM", Stock=NULL, RMD='OM', nsamp=nsamp, nsim=nsim, 
              nyears=nyears, proyears=proyears,
              output_file=output_file, output_dir=output_dir, quiet=quiet,
              tabs=tabs, title=title, date=date,
              plotPars=plotPars, open=open, dev=dev)
  
  
}

#### --- Old Code ----------------------------------------------------------####
#' Plot the Historical Spatial Closures
#'
#' @param OM An object of class OM
#' @param sim Optional. Simulation number to plot
#'
#' @export
#' @author A. Hordyk
#'
#' @examples
#' \dontrun{
#' OM <- new("OM", Albacore, Generic_Fleet, Perfect_Info, Perfect_Imp)
#' 
#' ## 50% of Area 1 was closed 30 years ago 
#' cl1 <- c(OM@nyears-30, 0.5, 1)
#' ## 80% of Area 1 was closed 15 years ago
#' cl2 <- c(OM@nyears-15, 0.2, 1)
#' ## 100% of Area 1 was closed last year
#' cl3 <- c(OM@nyears-1, 0, 1)
#' 
#' OM@MPA <- matrix(c(cl1, cl2, cl3), ncol=3, byrow=TRUE)
#' plotMPA(OM)
#' }
#' 
plotMPA <- function(OM, sim=NA) {
  .Deprecated('plot("MPA", Fleet, Stock')
  if (class(OM)!="OM") stop("Object must be class 'OM'")
  if (class(OM)=="OM") {
    proyears <- OM@proyears
  } else {
    proyears <- 0
  }
  set.seed(OM@seed)
  if (is.na(sim)) sim <- ceiling(runif(1, 1,OM@nsim))
  if (sim > OM@nsim) sim <- OM@nsim
  
  nyears <- OM@nyears
  if (all(!is.na(OM@MPA)) && sum(OM@MPA) != 0) { # historical spatial closures have been specified
    nareas <- ncol(OM@MPA)-1
    MPA <- matrix(1, nyears+proyears, ncol=nareas)
    yrindex <- OM@MPA[,1]
    if (max(yrindex)>nyears) stop("Invalid year index for spatial closures: must be <= nyears")
    if (min(yrindex)<1) stop("Invalid year index for spatial closures: must be > 1")
    for (xx in seq_along(yrindex)) {
      MPA[yrindex[xx]:nrow(MPA),] <- matrix(OM@MPA[xx, 2:ncol(OM@MPA)], nrow=length(yrindex[xx]:nrow(MPA)),ncol=nareas, byrow = TRUE)
    }
  } else {
    stop("No historical MPAs. MPA slot is empty")
  }
  
  x <- 1:(nyears+proyears)
  nyrs <- length(x)
  op <- par(mfrow=c(1,1), mar=c(3,3,0,0), oma=c(0,0,1,0), no.readonly = TRUE)
  on.exit(par(op))
  
  OM@Prob_staying <- c(0.5,0.5)
  Stock <- SampleStockPars(OM, cpars=OM@cpars, msg = FALSE)
 
  area_sizes <- Stock$Asize[sim,]

  plot(c(1, nyears+proyears), c(0,sum(area_sizes)), type="n", bty="n", xlab="", ylab="", axes=FALSE)

  origin <- cumsum(c(0, area_sizes)) #  seq(0, by=1, length.out = nareas)
  for (aa in 1:nareas) {
    polygon(x=c(x, rev(x)), y=c(rep(origin[aa], nyrs), origin[aa+1]*rev(MPA[,aa])), 
            col='lightgray', border = TRUE)
  }
  
  if (OM@CurrentYr < 1000) years <- (OM@CurrentYr - nyears+1) : (OM@CurrentYr+proyears) -OM@CurrentYr
  if (OM@CurrentYr > 1000) years <- (OM@CurrentYr - nyears+1) : (OM@CurrentYr+proyears) 
  
  xp <- seq(from=min(x), to=max(x), by=5)
  ind <- match(xp, x)
  axis(side=1, at=x[ind], labels=years[ind])
  
  mtext(side=1, "Years", line=2, xpd=NA, cex=1.25)
  
  axis(side=2, at=origin[1:(length(origin)-1)] + 0.5 * area_sizes, labels=1:nareas, las=1, col = "white", tcl = 0)
  
  mtext(side=2, "Areas", line=2, xpd=NA, cex=1.25, las=3)
  abline(v=nyears, lty=2, col="darkgray")
  
  title(paste0('Fraction open to fishing (grey) (sim = ', sim, ")"), outer=TRUE)
  
}

#' Plot the vulnerability and retention curves 
#'
#' @param OM An object of class 'OM' 
#' @param Pars Named list of sampled parameters
#' @param pyears number of years to plot
#' @param sim the simulation to plot. default is NA to plot a random simulation.
#' Set to 1 for reproducible plot
#' @param type plot type - line "l", point "p", or both "b"
#' @author A. Hordyk
#' @export
#'
plotSelect <- function(OM, Pars=NULL, pyears=4, sim=NA, type="l") {
  .Deprecated('plot("Selectivity", Fleet, Stock')
  if (class(OM) != "OM") stop("Object must be class 'OM' ")
  nsim <- OM@nsim
  if (!is.na(sim)) OM@nsim <- nsim <- max(10, sim) # OM@nsim 
  years <- OM@nyears + OM@proyears
  yr.vert <- round(seq(1, years, length.out=pyears),0)
  if (is.null(Pars)) {
    stckPars <- SampleStockPars(OM, msg=FALSE)
    Pars <- c(stckPars,SampleFleetPars(OM, Stock=stckPars, msg=FALSE))
  }
  
  set.seed(OM@seed)
  if (is.na(sim)) sim <- sample(1:nsim, 1)
  
  if (pyears > 1) gr <- c(2, pyears)
  if (pyears == 1) gr <- c(1, 2)
  op <- par(mfcol=gr, bty="l", las=1, mar=c(3,3,2,0), oma=c(3,3,2,1), xpd=NA)
  on.exit(par(op))
  
  for (yr in yr.vert) {
    # plot vulnerability & selection at age
    plot(1:Pars$maxage, Pars$V2[sim,, yr], type=type, ylim=c(0,1), lwd=2, 
         axes=FALSE, ylab="", xlab="")
    axis(side=1)
    mtext(side=1, "Age", line=2.5)
    if (yr == yr.vert[1]) {
      axis(side=2)
      mtext(side=2, "Vulnerabilty/Retention", las=3, line=3)
    }
    if (yr != yr.vert[1]) axis(side=2, labels=FALSE)
    if (pyears > 1) title(paste("Year", yr))
    if (pyears == 1) title(paste("Year", yr), outer=TRUE)
    
    polygon(x=c(1:Pars$maxage, rev(1:Pars$maxage)), 
            y=c(Pars$V[sim,, yr], rev(Pars$retA[sim,, yr])), col="gray", border=FALSE)
    lines(1:Pars$maxage, Pars$V[sim,, yr], col=2, lwd=2, lty=2, type=type)
    lines(1:Pars$maxage, Pars$retA[sim,, yr], col=4, lwd=2, lty=3, type=type)
    
    if (yr == yr.vert[pyears]) {
      minval <- min(c(Pars$V[sim,Pars$maxage, yr],  Pars$retA[sim,Pars$maxage, yr]))
      if (minval >= 0.5) loc <- "bottomright"
      if (minval < 0.5) loc <- "topright"
      legend(loc, legend = c("Vulnerability", "Realized Selection", "Retention"),
             lwd=2, col=c(1, 2, 4), bty="n", lty=c(1,2,3))
    }
    
    # plot vulnerability & selection at length
    plot(Pars$CAL_binsmid, Pars$SLarray2[sim,, yr], type=type, ylim=c(0,1), lwd=2, 
         axes=FALSE, ylab="", xlab="")
    axis(side=1)
    mtext(side=1, "Length", line=2.5)
    if (yr == yr.vert[1]) {
      axis(side=2)
      mtext(side=2, "Vulnerabilty/Retention", las=3, line=3)
    }
    if (yr != yr.vert[1]) axis(side=2, labels=FALSE)
    
    polygon(x=c(Pars$CAL_binsmid, rev(Pars$CAL_binsmid)), 
            y=c(Pars$SLarray[sim,, yr], rev(Pars$retL[sim,, yr])), col="gray", border=FALSE)
    lines(Pars$CAL_binsmid, Pars$SLarray[sim,, yr], col=2, lwd=2, lty=2,type=type)
    lines(Pars$CAL_binsmid, Pars$retL[sim,, yr], col=4, lwd=2, lty=3, type=type)
    
  }
  mtext(side=3, paste0("Selection and Retention curves for simulation: ", sim),  outer=TRUE)
  
  invisible(list(V=Pars$V, V2=Pars$V2 ,retA=Pars$retA, DR=Pars$DR, 
                 Fdisc=Pars$Fdisc, sim=sim))
  
  
}

#' Plot M-at-Age and Size
#'
#' @param Stock An object of class 'Stock' or 'OM' 
#' @param nsim The number of simulations to plot
#'
#' @author A. Hordyk
#' @export
#'
#' @examples 
#' \dontrun{
#' plotM(Albacore)
#' }
plotM <- function(Stock, nsim=5) {
  .Deprecated('plot("M", Stock')
  if (class(Stock) != "Stock" && class(Stock) != "OM") stop("Must supply object of class 'Stock' or 'OM'")
  
  nyears <- 30
  proyears <- 30
  SampCpars <- list() # empty list 
  if (class(Stock) == "OM") {
    # custom parameters exist - sample and write to list
    if(length(Stock@cpars)>0){
      # ncparsim<-cparscheck(Stock@cpars)   # check each list object has the same length and if not stop and error report
      SampCpars <- SampleCpars(Stock@cpars, nsim) 
    }
    nyears <- Stock@nyears 
    proyears <- Stock@proyears
    Stock@nsim <- nsim
  }
  
  StockPars <- SampleStockPars(Stock, nsim, nyears, proyears, SampCpars, msg=FALSE)
  # Assign Stock pars to function environment
  for (X in 1:length(StockPars)) assign(names(StockPars)[X], StockPars[[X]])
  
  M_at_age <- StockPars$M_ageArray
  Len_at_age <- StockPars$Len_age
  Wt_at_age <- StockPars$Wt_age
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  par(mfrow=c(3,3), bty="l", las=1, mar=c(3,3,2,1), oma=c(2,2,0,0))
  ylim <- c(0, max(M_at_age))
  lwd <- 2
  matplot(t(M_at_age[,,1]), type="l", xlab="", ylab="", ylim=ylim, lwd=lwd)
  mtext(side=2, "M", line=3)
  # mtext(side=1, "Age", line=2.5)
  matplot(t(Len_at_age[,,1]), t(M_at_age[,,1]), type="l", xlab="", ylab="", ylim=ylim, lwd=lwd)
  # mtext(side=1, "Length", line=2.5)
  title("First Historical Year")
  matplot(t(Wt_at_age[,,1]), t(M_at_age[,,1]), type="l", xlab="", ylab="", ylim=ylim, lwd=lwd)
  # mtext(side=1, "Weight", line=2.5)
  
  matplot(t(M_at_age[,,nyears]), type="l", xlab="", ylab="", ylim=ylim, lwd=lwd)
  mtext(side=2, "M", line=3)
  # mtext(side=1, "Age", line=2.5)
  matplot(t(Len_at_age[,,nyears]), t(M_at_age[,,nyears]), type="l", xlab="", ylab="", ylim=ylim, lwd=lwd)
  # mtext(side=1, "Length", line=2.5)
  title("Last Historical Year")
  matplot(t(Wt_at_age[,,nyears]), t(M_at_age[,,nyears]), type="l", xlab="", ylab="", ylim=ylim, lwd=lwd)
  # mtext(side=1, "Weight", line=2.5)
  
  matplot(t(M_at_age[,,proyears]), type="l", xlab="", ylab="", ylim=ylim, lwd=lwd)
  mtext(side=2, "M", line=3)
  mtext(side=1, "Age", line=2.5)
  matplot(t(Len_at_age[,,proyears]), t(M_at_age[,,proyears]), type="l", xlab="", ylab="", ylim=ylim, lwd=lwd)
  mtext(side=1, "Length", line=2.5)
  title("Last Projection Year")
  matplot(t(Wt_at_age[,,proyears]), t(M_at_age[,,proyears]), type="l", xlab="", ylab="", ylim=ylim, lwd=lwd)
  mtext(side=1, "Weight", line=2.5)
  
}




#' Wrapper for histogram function 
#' 
#' Produces a blank plot if all values in x are equal 
#'
#' @param x A vector of values
#' @param col Colour of the histogram
#' @param axes Logical - should axes be included?
#' @param main Character - main title
#' @param breaks Number of breaks. See ?hist for more details
#' @param cex.main Text size of the main title
#' 
#' @export
#'
hist2 <- function(x, col, axes=FALSE, main="", breaks=10,cex.main=1) {
  if (mean(x) == x[1]) {
    
    plot(mean(x)*c(0.9,1.1),c(0,1000),col="white",axes=F,main=main,xlab="",ylab="")
    axis(1)
    #hist(x, border="white", xlim=range(x),xlab="", ...)
    abline(v=mean(x),col=col,lwd=2)
    
  } else {
    col="dark grey"
    hist(x, border='white',xlab="",col=col,axes=axes,main=main,breaks=breaks, ylab="")
  }
}







# #' @method plot Stock
# #' @export
# plot.Stock <- function(x, ...)  plotStock(x, ...)

#' Plot the Stock object parameters 
#' 
#' A function that plots histograms of samples from the Stock object parameters,
#' and time-series plots of `nsamp` samples of time-varying parameters. Used to 
#' visually examine the parameter values and ranges entered into the Stock object.
#' 
#' @param x An object of class Stock (or of class OM) 
#' @param nsamp Number of random samples for time-series plots
#' @param nsim Number of iterations for histograms. Ignored if x is class 'OM'
#' @param nyears Number of historical years. Ignored if x is class 'OM'
#' @param proyears Number of projection years. Ignored if x is class 'OM'
#' @param col Color of histograms 
#' @param breaks Number of breaks for histograms 
#' @param lwd line width 
#' @param ask Ask before displaying next page?
#' @param incVB Show the sampled von Bertalanffy growth curves on second page?
#' @param ...  Optional additional arguments passed to \code{plot}
#' @author A. Hordyk
#' @export 
plotStock <- function(x, nsamp=3, nsim=500, nyears=50, proyears=28, 
                      col="darkgray", breaks=10, lwd=2, ask=FALSE, incVB=TRUE, ...) {
  Stock <- x 
  .Deprecated('plot.Stock')
  SampCpars <- list() # empty list 
  if (class(Stock) == "OM") {
    Stock <- updateMSE(x) 
    if (is.finite(Stock@nyears)) nyears <- Stock@nyears
    if (is.finite(Stock@proyears)) proyears <- Stock@proyears
    if (is.finite(Stock@nsim)) nsim <- Stock@nsim	
    
    if(length(Stock@cpars)>0){ # custom parameters exist - sample and write to list
      #ncparsim<-cparscheck(Stock@cpars)   # check each list object has the same length and if not stop and error report
      SampCpars <- SampleCpars(Stock@cpars, nsim, msg=FALSE) 
    }
    Stock <- SubOM(Stock)
  }
  its <- sample(1:nsim, nsamp)
 
  
  # --- Sample Stock Parameters ----

  StockPars <- SampleStockPars(Stock, nsim, nyears, proyears, SampCpars, msg=FALSE)
 
  # Assign Stock pars to function environment
  for (X in 1:length(StockPars)) assign(names(StockPars)[X], StockPars[[X]])
  
  
  ncol <- 8
  
  m <- layout(matrix(c(c(1, 2, 3, 0, 4, 4, 4, 5),
                       c(6, 7, 8, 0, 9, 10, 11, 31),
                       c(12, 13, 14, 0, 15, 15, 15, 16),			  
                       c(17, 18, 19, 0, 20, 20, 20, 0),
                       c(21, 22, 23 ,0, 24, 24, 24, 25),
                       c(26, 27, 28, 0, 29, 29, 30, 30)), 
                     ncol=ncol, byrow = TRUE), 
              widths=c(1, 1, 1, 0.5, 1, 1, 1, 1))
  
  # layout.show(m)
  # stop()
  op <- par(mar = c(2, 1, 3, 1), oma=c(1,2,4,1), ask=FALSE, las=1)
  on.exit(par(op))
  
  # Row 1 -- Natural Mortality ---- 
  hist2(M, col=col, axes=FALSE, main="M", breaks=breaks)
  abline(v=M[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  hist2(Msd, col=col, axes=FALSE, main="Msd", breaks=breaks)
  abline(v=Msd[its], col=1:nsamp, lwd=lwd)
  axis(side=1) 
  plot(c(0,1), c(0,1), type="n", axes=FALSE, xlab="", ylab="")
  # hist2(Mgrad, col=col, axes=FALSE, main="Mgrad", breaks=breaks)
  # abline(v=Mgrad[its], col=1:nsamp, lwd=lwd)
  # axis(side=1)
  
  # M traj
  matplot(t(Marray[its,]), type="l", lty=1, bty="l", main="M by Year", lwd=lwd)
  
  # M/K
  hist2(M/K, col=col, axes=FALSE, main="M/K", breaks=breaks)
  abline(v=(M/K)[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  #  Row 2 -- M-at-age/length ---- 
  # M-at-length
  lims <- range(M_ageArray[its,, ])
  xlims <- range(Len_age[its,,])
  matplot(t(Len_age[its,,1]), t(M_ageArray[its,,1]), type="l", lty=1, bty="l", lwd=lwd, ylim=lims, xlim=xlims)
  mtext(side=3, "First historical year", cex=0.8, line=-1)
  mtext(side=1, "Length", line=2, cex=0.7)
  
  matplot(t(Len_age[its,,nyears]), t(M_ageArray[its,,nyears]), type="l", lty=1, bty="l", main="M-at-length", lwd=lwd, ylim=lims, xlim=xlims, axes=FALSE)
  mtext(side=3, "Last historical year", cex=0.8, line=-1)
  axis(side=1)
  mtext(side=1, "Length", line=2, cex=0.7)
  matplot(t(Len_age[its,,nyears+proyears]), t(M_ageArray[its,,nyears+proyears]), type="l", lty=1, bty="l", lwd=lwd, ylim=lims, axes=FALSE, xlim=xlims)
  mtext(side=3, "Last projected year", cex=0.8, line=-1)
  axis(side=1)
  mtext(side=1, "Length", line=2, cex=0.7)
  
  # M-at-age
  lims <- range(M_ageArray[its,, ])
  matplot(t(M_ageArray[its,,1]), type="l", lty=1, bty="l", lwd=lwd, ylim=lims)
  mtext(side=3, "First historical year", cex=0.8, line=-1)
  mtext(side=1, "Age", line=2, cex=0.7)
  
  matplot(t(M_ageArray[its,,nyears]), type="l", lty=1, bty="l", main="M-at-age", lwd=lwd, ylim=lims, axes=FALSE)
  mtext(side=3, "Last historical year", cex=0.8, line=-1)
  axis(side=1)
  mtext(side=1, "Age", line=2, cex=0.7)
  matplot(t(M_ageArray[its,,nyears+proyears]), type="l", lty=1, bty="l", lwd=lwd, ylim=lims, axes=FALSE)
  mtext(side=3, "Last projected year", cex=0.8, line=-1)
  axis(side=1)
  mtext(side=1, "Age", line=2, cex=0.7)
  
  
  
  # Row 3 -- Linf ---- 
  hist2(Linf, col=col, axes=FALSE, main="Linf", breaks=breaks)
  abline(v=Linf[its], col=1:nsamp, lwd=lwd)
  axis(side=1) 
  hist2(Linfsd, col=col, axes=FALSE, main="Linfsd", breaks=breaks)
  abline(v=Linfsd[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  plot(c(0,1), c(0,1), type="n", axes=FALSE, xlab="", ylab="")
  # hist2(Linfgrad, col=col, axes=FALSE, main="Linfgrad", breaks=breaks)
  # abline(v=Linfgrad[its], col=1:nsamp, lwd=lwd)
  # axis(side=1)
  
  # Linf traj 
  matplot(t(Linfarray[its,]), type="l", bty="l", main="Linf by Year", lwd=lwd, lty=1)
  
  hist2(t0, col=col, axes=FALSE, main="t0", breaks=breaks)
  abline(v=t0[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  
  # Row 4 -- K ----
  hist2(K, col=col, axes=FALSE, main="K", breaks=breaks)
  abline(v=K[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  hist2(Ksd, col=col, axes=FALSE, main="Ksd", breaks=breaks)
  abline(v=Ksd[its], col=1:nsamp, lwd=lwd)
  axis(side=1) 
  plot(c(0,1), c(0,1), type="n", axes=FALSE, xlab="", ylab="")
  # hist2(Kgrad, col=col, axes=FALSE, main="Kgrad", breaks=breaks)
  # abline(v=Kgrad[its], col=1:nsamp, lwd=lwd)
  # axis(side=1)  
  
  # K traj 
  matplot(t(Karray[its,]), type="l", bty="l", main="K by Year", lwd=lwd, lty=1)
  
  
  # Row 5 -- Recruitment ----
  hist2(hs, col=col, axes=FALSE, main="Steepness (h)", breaks=breaks)
  abline(v=hs[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  hist2(procsd, col=col, axes=FALSE, main="procsd", breaks=breaks)
  abline(v=procsd[its], col=1:nsamp, lwd=lwd)
  axis(side=1) 
  hist2(AC, col=col, axes=FALSE, main="AC", breaks=breaks)
  abline(v=AC[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  # matplot(t(Perr[its,]), type="l", bty="l", main="Rec Devs by Year", lwd=lwd, lty=1)
  matplot(t(Perr_y[its,]), type="l", bty="l", main="Rec Devs by Year", lwd=lwd, lty=1)
  
  biomass <- seq(0, 1.5, by=0.05)
  
  SSB0a <- 1
  SSBpR <- 1 
  bR <- matrix(log(5 * hs)/(0.8 * SSB0a), nrow=nsim)  # Ricker SR params
  aR <- matrix(exp(bR * SSB0a)/SSBpR, nrow=nsim)  # Ricker SR params
  
  if (SRrel[1] == 1) {
    recs <- sapply(its, function(X) (0.8  * hs[X] * biomass)/(0.2 * (1 - hs[X]) + (hs[X] - 0.2) * biomass)) # BH SRR 
  } else {
    # most transparent form of the Ricker uses alpha and beta params
    recs <- sapply(its, function(X) aR[X] * biomass * exp(-bR[X] * biomass))
  }
  matplot(biomass, recs, type="l", bty="l", main="Stock-Recruit", 
          ylim=c(0,max(1, max(recs))), xlim=c(0,max(biomass)), ylab="", xlab="", lwd=lwd, lty=1)
  
  
  
  # Row 6 -- Depletion and Maturity ----
  hist2(D, col=col, axes=FALSE, main="Depletion", breaks=breaks)
  abline(v=D[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  # Maturity 
  hist2(L50, col=col, axes=FALSE, main="L50", breaks=breaks)
  abline(v=L50[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  hist2(L95, col=col, axes=FALSE, main="L95", breaks=breaks)
  abline(v=L95[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  slope <- log(19)/(L95-L50)
  Ls <- seq(0, to=max(Linf), length.out=200)
  
  Mat_len <- sapply(its, function(X)  plogis(Ls, L50[X], 1/slope[X]))
  matplot(Ls, Mat_len, type="l", bty="l", main="Maturity-at-length (current year)", lwd=lwd, lty=1)
  
  if (length(dim(Mat_age)) == 2) matplot(t(Mat_age[its,]), type="l", bty="l", main="Maturity-at-age", lwd=lwd, lty=1, axes=FALSE, xlim=c(0, maxage))
  if (length(dim(Mat_age)) == 3) matplot(t(Mat_age[its,, nyears]), type="l", bty="l", main="Maturity-at-age (current year)", lwd=lwd, lty=1, axes=FALSE, xlim=c(0, maxage))
  axis(side=1)
  axis(side=2, labels=FALSE)
  
  # Add Mexp 
  hist2(Mexp, col=col, axes=FALSE, main="Mexp", breaks=breaks)
  abline(v=Mexp[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  
  title(Stock@Name, outer=TRUE)
  title(paste("nyears =", nyears, "  proyears =", proyears, "  ", nsamp, "sampled iterations"), outer=TRUE, line=0)
  
  if (incVB) {
    par(mfrow=c(1,3), ask=ask, mar=c(5,4,1,1))
    # vB   
    cex.lab <- 2 
    fstYr <- Len_age[its,,1]
    curYr <- Len_age[its,,nyears]
    lstYr <- Len_age[its,,proyears+nyears]
    MaxL <- max(Len_age)
    matplot(t(fstYr), type="l", bty="l", main="First historical year", ylim=c(0, MaxL), xlab="Age", ylab="Length", cex.lab=cex.lab, lwd=lwd, lty=1)
    matplot(t(curYr), type="l", bty="l", main="Last historical year", ylim=c(0, MaxL),  axes=FALSE, xlab="Age", ylab="", cex.lab=cex.lab, lwd=lwd, lty=1)
    axis(side=1)
    axis(side=2, labels=FALSE)  
    matplot(t(lstYr), type="l", bty="l", main="Last projected year", ylim=c(0, MaxL), axes=FALSE, xlab="Age", ylab="", cex.lab=cex.lab, lwd=lwd, lty=1)	
    axis(side=1)
    axis(side=2, labels=FALSE)  
    title("Sampled length-at-age curves", outer=TRUE, cex.main=2)
  }
  
  invisible(StockPars)
}




# #' @method plot Fleet
# #' @export
# plot.Fleet <- function(x, ...)  plotFleet(x, ...)

#' Plot the Fleet object parameters 
#' 
#' @param x An object of class Fleet (or of class OM)
#' @param Stock  An object of class Stock 
#' @param nsamp Number of random samples for time-series plots
#' @param nsim Number of iterations for histograms
#' @param proyears Number of projection years 
#' @param col Color of histograms 
#' @param breaks Number of breaks for histograms 
#' @param lwd line width 
#' @param ...  Optional additional arguments passed to \code{plot}
#' @rdname plot-Fleet 
#' @author A. Hordyk
#' @export 
plotFleet <- function(x, Stock=NULL, nsamp=3, nsim=500, proyears=28, col="darkgray", 
                      breaks=10, lwd=2, ...) { 
  .Deprecated('plot.Fleet')
  Fleet <- updateMSE(x) # add missing slots with default values 
  SampCpars <- list() # empty list 
  nyears <- Fleet@nyears
  if (class(Fleet) == "OM") {
    if (is.finite(Fleet@proyears)) proyears <- Fleet@proyears
    if (is.finite(Fleet@nsim)) nsim <- Fleet@nsim	
    if (length(Fleet@cpars) > 0) {
      # ncparsim <-cparscheck(Fleet@cpars)   # check each list object has the same length and if not stop and error report
      SampCpars <- SampleCpars(Fleet@cpars, nsim, msg=FALSE) 
    }
    Stock <- SubOM(Fleet, "Stock")
    Fleet <- SubOM(Fleet, "Fleet")
  }
  if (class(Stock) != "Stock") stop("Must include a Stock object", call.=FALSE)
  
  its <- sample(1:nsim, nsamp)  
  
  StockPars <- SampleStockPars(Stock, nsim, nyears, proyears, SampCpars, msg=FALSE)
  # Assign Stock pars to function environment
  for (X in 1:length(StockPars)) assign(names(StockPars)[X], StockPars[[X]])
  
  FleetPars <- SampleFleetPars(Fleet, StockPars, nsim, nyears, proyears, SampCpars)
  # Assign Stock pars to function environment
  for (X in 1:length(FleetPars)) assign(names(FleetPars)[X], FleetPars[[X]])
  
  # Start plotting shenanigans 
  ncol <- 12
  
  m <- layout(matrix(c(c(rep(1,3), rep(2,3), rep(3,3), rep(4,3)),
                       c(rep(0, ncol)),
                       c(rep(5,6), rep(6,6)),
                       c(rep(0, ncol)),
                       c(rep(7,4), rep(8, 4), rep(9,4)),
                       c(rep(0, ncol)),
                       c(0, 10, 11, 12, 0, 13, 14, 15, 0, 16, 17, 18),
                       c(rep(0, ncol)),
                       c(rep(19,4), rep(20,4), rep(21,4))
  ), ncol=ncol, byrow = TRUE),
  heights=c(1,0.2, 1,0.3, 1, 0.3, 1, 0.3, 1))
  # layout.show(m)				   
  
  op <- par(mar = c(2,2, 2, 1), oma=c(2,2,4,1), las=1)
  on.exit(par(op))	 
  hist2(Esd, col=col, axes=FALSE, main="Esd", breaks=breaks)
  abline(v=Esd[its], col=1:nsamp, lwd=lwd)
  axis(side=1) 
  hist2(qinc, col=col, axes=FALSE, main="qinc", breaks=breaks)
  abline(v=qinc[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  hist2(qcv, col=col, axes=FALSE, main="qcv", breaks=breaks)
  abline(v=qcv[its], col=1:nsamp, lwd=lwd)
  axis(side=1) 
  hist2(Spat_targ, col=col, axes=FALSE, main="Spat_targ", breaks=breaks)
  abline(v=Spat_targ[its], col=1:nsamp, lwd=lwd)
  axis(side=1) 
  # title("not currently used", line=0)  
  
  
  # Effort 
  matplot(t(Find[its,]), type="l", lwd=lwd, bty="l", main="Fishing effort\ (historical)", lty=1)
  
  # Future catchability
  ind <- as.matrix(expand.grid(its, 1:proyears, 1:nsamp))
  Qfuture <- matrix(NA, nrow=proyears, ncol=nsamp)
  X <- 0 
  for (sim in its) {
    X <- X + 1 
    Qfuture[,X] <- qvar[sim,1] * (1 + qinc[sim]/100)^(1:proyears)
    Qfuture[,X] <- Qfuture[,X]/Qfuture[1,X]
  }
  matplot(Qfuture, type="l", lwd=lwd, bty="l", main="Change in future\n catchability", lty=1)  
  
  # Selectivity at length
  sampV <- V[its,,]
  sampVL <- SLarray[its,,]
  matplot(CAL_binsmid, t(sampVL[,,1]), type="l", lwd=lwd, bty="l", main="First historical\n year", 
          xlab="Length", ylim=c(0,1), lty=1)  
  matplot(CAL_binsmid, t(sampVL[,,nyears]), type="l", lwd=lwd, bty="l", main="Last historical\n year", 
          xlab="Length", ylim=c(0,1), lty=1)  
  title(line=3, cex.main=1.5, "Selectivity-at-length", xpd=NA)
  matplot(CAL_binsmid, t(sampVL[,,nyears+proyears]), type="l", lwd=lwd, bty="l", 
          main="Last projected\n year", xlab="Length", ylim=c(0,1), lty=1)  
  
  # Selectivity Parameters #
  hist2(L5[1, ], col=col, axes=FALSE, main="L5", breaks=breaks)
  axis(side=1) 
  abline(v=L5[1, its], col=1:nsamp, lwd=lwd)
  hist2(LFS[1,], col=col, axes=FALSE, main="LFS", breaks=breaks)
  axis(side=1) 
  abline(v=LFS[1, its], col=1:nsamp, lwd=lwd)
  hist2(Vmaxlen[1,], col=col, axes=FALSE, main="Vmaxlen", breaks=breaks)
  axis(side=1) 
  abline(v=Vmaxlen[1, its], col=1:nsamp, lwd=lwd)
  
  
  hist2(L5[nyears, ], col=col, axes=FALSE, main="L5", breaks=breaks)
  axis(side=1) 
  abline(v=L5[nyears, its], col=1:nsamp, lwd=lwd)
  hist2(LFS[nyears,], col=col, axes=FALSE, main="LFS", breaks=breaks)
  axis(side=1) 
  abline(v=LFS[nyears, its], col=1:nsamp, lwd=lwd)
  hist2(Vmaxlen[nyears,], col=col, axes=FALSE, main="Vmaxlen", breaks=breaks)
  axis(side=1) 
  abline(v=Vmaxlen[nyears, its], col=1:nsamp, lwd=lwd)
  
  
  hist2(L5[nyears+proyears, ], col=col, axes=FALSE, main="L5", breaks=breaks)
  axis(side=1) 
  abline(v=L5[nyears+proyears, its], col=1:nsamp, lwd=lwd)
  hist2(LFS[nyears+proyears,], col=col, axes=FALSE, main="LFS", breaks=breaks)
  axis(side=1) 
  abline(v=LFS[nyears+proyears, its], col=1:nsamp, lwd=lwd)
  hist2(Vmaxlen[nyears+proyears,], col=col, axes=FALSE, main="Vmaxlen", breaks=breaks)
  axis(side=1) 
  abline(v=Vmaxlen[nyears+proyears, its], col=1:nsamp, lwd=lwd)
  
  
  # Selecitivty at age 
  matplot(t(sampV[,,1]), type="l", lwd=lwd, bty="l", main="First historical\n year", xlab="Age", ylim=c(0,1), lty=1)  
  matplot(t(sampV[,,nyears]), type="l", lwd=lwd, bty="l", main="Last historical\n year", xlab="Age", ylim=c(0,1), lty=1)  
  title(line=3, cex.main=1.5, "Selectivity-at-age", xpd=NA)
  matplot(t(sampV[,,nyears+proyears]), type="l", lwd=lwd, bty="l", main="Last projected\n year", xlab="Age", 
          ylim=c(0,1), lty=1)  
  
  title(paste("Fleet:", Fleet@Name, "Stock:", Stock@Name), outer=TRUE, line=1.5)
  title(paste("nyears =", nyears, "  proyears =", proyears, "  ", nsamp, "sampled iterations"), outer=TRUE, line=0)
  
  # om <- new("OM", Stock, Fleet, Perfect_Info, Perfect_Imp)
  # plotSelect(om)
  
  invisible(FleetPars)
}


# #' @method plot Obs
# #' @export
# plot.Obs <- function(x, ...)  plotObs(x, ...)

#' Plot the Observation object parameters 
#' 
#' A function that plots histograms of samples from the observation object parameters,
#' and time-series plots of `nsamp` samples of time-series examples. Used to 
#' visually examine the parameter values and ranges entered into the Obs object.
#' 
#' @param x An object of class Obs (or of class OM) 
#' @param nsim Number of iterations for histograms
#' @param nyears Number of historical years
#' @param col Color of histograms 
#' @param breaks Number of breaks for histograms 
#' @param ...  Optional additional arguments passed to \code{plot}
#' @rdname plot-Obs 
#' @author T. Carruthers and A. Hordyk
#' @export 
plotObs <- function(x, nsim=500, nyears=50, 
                    col="darkgray", breaks=10, ...) {
  .Deprecated('plot.Obs')
  Obs <- x
  SampCpars <- list() # empty list 
  if (class(Obs) == "OM") {
    if (is.finite(Obs@nyears)) nyears <- Obs@nyears
    if (is.finite(Obs@nsim)) nsim <- Obs@nsim	
    if (length(Obs@cpars) > 0) {
      # ncparsim <-cparscheck(Obs@cpars)   # check each list object has the same length and if not stop and error report
      SampCpars <- SampleCpars(Obs@cpars, nsim, msg=FALSE) 
    }
    Obs <- SubOM(Obs,"Obs")
  }

  nsamp <- 3
  its <- sample(1:nsim, nsamp)
  
  # === Sample Observation Model Parameters ====
  ObsPars <- SampleObsPars(Obs, nsim, cpars = SampCpars)
  # Assign Obs pars to function environment
  for (X in 1:length(ObsPars)) assign(names(ObsPars)[X], ObsPars[[X]])
  
  
  # === Non time series ==================================================================== 
  cex.main <- 0.5
  op <- par(mfrow=c(4,4),mai=c(0.6,0.6,0.25,0.01),omi=c(0.01,0.01,0.4,0.01))
  on.exit(par(op))	 
  
  hist2(CAA_nsamp,col=col,axes=FALSE, main="No. annual catch-at-age obs (CAA_samp)", breaks=breaks,cex.main=cex.main)
  axis(side=1) 
  
  hist2(CAA_ESS,col=col, axes=FALSE, main="Effective sample size CAA obs (CAA_ESS)", breaks=breaks,cex.main=cex.main)
  axis(side=1) 
  
  hist2(CAL_nsamp,col=col, axes=FALSE, main="No. annual catch-at-length obs (CAL_samp)", breaks=breaks,cex.main=cex.main)
  axis(side=1) 
  
  hist2(CAL_ESS,col=col, axes=FALSE, main="Effective sample size CAL obs (CAL_ESS)", breaks=breaks,cex.main=cex.main)
  axis(side=1) 
  
  hist2(Mbias,col=col, axes=FALSE, main="Natural mortality rate bias (Mbias)", breaks=breaks,cex.main=cex.main)
  axis(side=1) 
  
  hist2(FMSY_Mbias,col=col, axes=FALSE, main="FMSY/M bias (FMSY_Mbias)", breaks=breaks,cex.main=cex.main)
  axis(side=1)
  
  hist2(lenMbias,col=col, axes=FALSE, main="Bias in length at maturity (lenMbias)", breaks=breaks,cex.main=cex.main)
  axis(side=1)
  
  hist2(LFCbias,col=col, axes=FALSE, main="Bias in length at first capture (LFCbias)", breaks=breaks,cex.main=cex.main)
  axis(side=1)
  
  hist2(LFSbias,col=col, axes=FALSE, main="Bias in length at full selection (LFSbias)", breaks=breaks,cex.main=cex.main)
  axis(side=1)
  
  hist2(Kbias,col=col, axes=FALSE, main="Bias in von B. K (Kbias)", breaks=breaks,cex.main=cex.main)
  axis(side=1)
  
  hist2(t0bias,col=col, axes=FALSE, main="Bias in von B. t0 (t0bias)", breaks=breaks,cex.main=cex.main)
  axis(side=1)
  
  hist2(Linfbias,col=col, axes=FALSE, main="Bias in von B. Linf (Linfbias)", breaks=breaks,cex.main=cex.main)
  axis(side=1)
  
  hist2(Irefbias,col=col, axes=FALSE, main="Bias in index at MSY (Irefbias)", breaks=breaks,cex.main=cex.main)
  axis(side=1)
  
  hist2(Crefbias,col=col, axes=FALSE, main="Bias in MSY catch (Crefbias)", breaks=breaks,cex.main=cex.main)
  axis(side=1)
  
  hist2(Brefbias,col=col, axes=FALSE, main="Bias in MSY biomass (Brefbias)", breaks=breaks,cex.main=cex.main)
  axis(side=1)
  
  hist2(Recsd,col=col, axes=FALSE, main="Bias in recent recruitment strength (Recsd)", breaks=breaks,cex.main=cex.main)
  axis(side=1)
  
  mtext(paste0("Observation biases and sample sizes for observation object ",Obs@Name),3,outer=T,line= 0.7,font=2)
  
  
  # ============= Time series ==============================================================
  
  op <- par(mfrow=c(3,3),mai=c(0.6,0.6,0.2,0.01),omi=c(0.01,0.01,0.4,0.01))
  
  m <- layout(matrix(c(c(1, 2, 3, 3),
                       c(4, 5, 6, 6),
                       c(7,8,9,9),
                       c(10,11,12,12)
  ), ncol=4, byrow = TRUE),
  widths=c(1, 1, 1, 1))
  
  
  # ---- Catches --------------
  
  
  ObsTSplot(Cbias,Csd,nyears,labs=c("Catch bias (Cbias)",
                                    "Catch error (Csd)","Catch discrepancy for three samples",
                                    "Cbias","Csd"), breaks=breaks, its=its, nsamp=nsamp, col=col)
  
  # --- Depletion -------------
  
  ObsTSplot(Dbias,Derr,nyears,labs=c("Depletion bias (Dbias)",
                                     "Depletion error (Derr)","Depletion discrepancy for three samples",
                                     "Depletion bias","Depletion error"), breaks=breaks, its=its, nsamp=nsamp, col=col)
  
  # --- Abundance -------------
  
  ObsTSplot(Abias,Aerr,nyears,labs=c("Current abundance bias (Abias)",
                                     "Current abundance error (Derr)","Abundance discrepancy for three samples",
                                     "Abias","Aerr"), breaks=breaks, its=its, nsamp=nsamp, col=col)
  
  # --- Indices --------------
  
  
  hist2(betas,col=col, axes=FALSE, main="Index hyper stability (betas)", breaks=breaks,cex.main=0.95)
  axis(side=1) 
  abline(v=1)
  abline(v=betas[its],col=makeTransparent(c("Black","Red","Green"),80),lwd=2)
  
  hist2(Isd,col=col, axes=FALSE, main="Index error (Isd)", breaks=breaks,cex.main=0.95)
  abline(v=Isd[its],col=makeTransparent(c("Black","Red","Green"),80),lwd=2)
  axis(side=1)  
  
  ind<-seq(1,0.1,length.out=nyears)
  Ierr <- array(rlnorm(nyears * nsamp, mconv(1, rep(Isd[its], nyears)), sdconv(1, rep(Isd[its], nyears))),  c(nsamp, nyears))  # composite of bias and observation error
  Imu<-array(rep(ind,each=nsamp)^rep(betas[its],nyears),c(nsamp,nyears))*Ierr
  Imu<-Imu/apply(Imu,1,mean)
  
  matplot(t(Imu),type='l',main="Three example indices",xlab="Year",ylab="Index (mean 1)",cex.main=0.95)
  lines(1:nyears,ind,col=makeTransparent("grey",70),lwd=3)
  legend('topleft',legend=round(betas[its],2),text.col=c("Black","Red","Green"),title="betas",bty='n')
  legend('topright',legend=round(Isd[its],2),text.col=c("Black","Red","Green"),title="Isd",bty='n')
  
  if (!is.na(Obs@Name)) mtext(paste0("Observation time series plots for observation object ",Obs@Name),3,outer=T,line= 0.7,font=2)
  if (is.na(Obs@Name)) mtext(paste0("Observation time series plots for observation object "),3,outer=T,line= 0.7,font=2)
  
  invisible(ObsPars)
}


# #' @method plot Imp
# #' @export
# plot.Imp <- function(x, ...)  plotImp(x, ...)

#' Plot the Implementation object parameters 
#' 
#' A function that plots histograms of samples from the implementation object parameters,
#' and time-series plots of `nsamp` samples of time-series examples. Used to 
#' visually examine the parameter values and ranges entered into the Obs object.
#' 
#' @param x An object of class Imp (or of class OM) 
#' @param nsim Number of iterations for histograms
#' @param nyears Number of historical years
#' @param col Color of histograms 
#' @param breaks Number of breaks for histograms 
#' @param ...  Optional additional arguments passed to \code{plot}
#' @rdname plot-Imp 
#' @author T. Carruthers and A. Hordyk
#' @export 
plotImp<-function(x,nsim=500, nyears=50, 
                  col="darkgray", breaks=10, ...){
  
  .Deprecated('plot.Imp')
  Imp <- x
  if (class(Imp) == "OM") {
    if (is.finite(Imp@nyears)) nyears <- Imp@nyears
    if (is.finite(Imp@nsim)) nsim <- Imp@nsim	
    Imp <- SubOM(Imp,"Imp")
  }
  
  # === Sample Imp Model Parameters ====
  ImpPars <- SampleImpPars(Imp, nsim)
  # Assign Imp pars to function environment
  for (X in 1:length(ImpPars)) assign(names(ImpPars)[X], ImpPars[[X]])
  
  nsamp=3
  its <- sample(1:nsim, nsamp)
  
  op <- par(mfrow=c(4,3),mai=c(0.6,0.6,0.2,0.01),omi=c(0.01,0.01,0.4,0.01))
  on.exit(par(op))	 
  
  ObsTSplot(TACFrac,TACSD,nyears,labs=c("Fraction of TAC (TACFrac)",
                                        "TAC error (TACSD)","TAC discrepancy for three samples",
                                        "TACFrac","TACSD"), breaks=breaks, its=its, nsamp=nsamp, col=col)
  
  ObsTSplot(TAEFrac,TAESD,nyears,labs=c("Fraction of effort (TAEFrac)",
                                        "Effort error (TAESD)","Effort discrepancy for three samples",
                                        "TAEFrac","TAESD"), breaks=breaks, its=its, nsamp=nsamp, col=col)
  
  ObsTSplot(SizeLimFrac,SizeLimSD,nyears,labs=c("Fraction of Size Limit (SizeLimFrac)",
                                                "Size Limit error (SizeLimSD)","Size limit discrepancy for three samples",
                                                "SizeLimFrac","SizeLimSD"), breaks=breaks, its=its, nsamp=nsamp, col=col)
  
  mtext(paste0("Implementation error time series plots for implementation object ",Imp@Name),3,outer=T,line= 0.7,font=2)
  
  invisible(ImpPars)
}



ObsTSplot<-function(Cbias,Csd,nyears,labs, breaks, its, nsamp, col){
  
  hist2(Cbias,col=col, axes=FALSE, main=labs[1], breaks=breaks,cex.main=0.95)
  if(sd(Cbias)>0.01)axis(side=1) 
  abline(v=1)
  abline(v=Cbias[its],col=makeTransparent(c("Black","Red","Green"),80),lwd=2)
  hist2(Csd,col=col, axes=FALSE, main=labs[2], breaks=breaks,cex.main=0.95)
  abline(v=Csd[its],col=makeTransparent(c("Black","Red","Green"),80),lwd=2)
  if(sd(Csd)>0.01)axis(side=1)  
  
  Cbiasa <- array(Cbias[its], c(nsamp, nyears))  # Bias array
  Cerr <- array(rlnorm(nyears * nsamp, mconv(1, rep(Csd[its], nyears)), sdconv(1, rep(Csd[its], nyears))),  c(nsamp, nyears))  # composite of bias and observation error
  matplot(t(Cbiasa*Cerr),type='l',main=labs[3],xlab="Year",ylab="Observed/real",cex.main=0.95)
  abline(h=1,col=makeTransparent("grey",70),lwd=3)
  legend('topleft',legend=round(Cbias[its],2),text.col=c("Black","Red","Green"),title=labs[4],bty='n')
  legend('topright',legend=round(Csd[its],2),text.col=c("Black","Red","Green"),title=labs[5],bty='n')
  
  abline(h=Cbias[its],col=makeTransparent(c("Black","Red","Green"),50),lwd=2)
  
}








#' Plot the operating model (OM) object parameters 
#' 
#' A function that plots the parameters and resulting time series of an operating model.
#' 
#' @param x An object of class OM or an object of class Hist (ie runMSE(OM, Hist=TRUE))
#' @param rmd Logical. Used in a rmd file?
#' @param head Character. Heading for rmd file. Default is '##' (second level heading)
#' @param ...  Optional additional arguments passed to \code{plot}
#' @rdname plot-OM
#' @author T. Carruthers
#' @export 
plotOM <-function(x, rmd=FALSE, head="##", ...){
  .Deprecated('plot.OM')
 op <- par(no.readonly = TRUE)
 on.exit(par(op))
 if (class(x) == "OM") {
   OM <- updateMSE(x) # update and add missing slots with default values
   out<-runMSE(OM,Hist=T, ...)
   nsim<-OM@nsim
   nyears<-OM@nyears
   if (rmd) {
     cat('\n')
     cat('\n')
     cat(head, 'Stock Object')
     cat('\n')
   }
   plotStock(OM)
   if (rmd) {
     cat('\n')
     cat('\n')
     cat(head, 'Fleet Object')
     cat('\n')
   }
   plotFleet(OM)
   if (rmd) {
     cat('\n')
     cat('\n')
     cat(head, 'Obs Object')
     cat('\n')
   }
   plotObs(OM)
   if (rmd) {
     cat('\n')
     cat('\n')
     cat(head, 'Imp Object')
     cat('\n')
   }
   plotImp(OM)
   yrlab<-OM@CurrentYr-((nyears-1):0)
 } else if (class(x) == "Hist") {
   out <- x 
   nyears <- dim(out@TSdata[[1]])[2]
   nsim <- dim(out@TSdata[[1]])[1]
   yrlab<-nyears-((nyears-1):0)
 } else stop("argument must be class 'OM' or 'Hist' ")
 
 if (rmd) {
   cat('\n')
   cat('\n')
   cat(head, 'OM Simulations')
   cat('\n')
 }
 
 # Time series
 par(mfrow=c(4,2),mai=c(0.7,0.7,0.05,0.05),omi=c(0.01,0.01,0.3,0.01))
 
 # SSB
 TSplot(yrlab,out@TSdata$SSB,xlab="Historical year",ylab="Spawning biomass")
 
 # Depletion
 TSplot(yrlab,out@TSdata$SSB/rep(out@Ref$SSB0,each=nyears),xlab="Historical year",ylab="Stock depletion (SSB)")
 
 # Apical F
 FM<-t(out@TSdata$Find*out@OM$qs)
 FM[FM > out@OM$maxF[1]] <- out@OM$maxF[1] # add maxF constraint
 TSplot(yrlab,t(FM),xlab="Historical year",ylab="Fishing mortality rate (apical)")
 
 # Catches
 TSplot(yrlab,out@TSdata$Catch,xlab="Historical year",ylab="Annual catches")
 
 # Recruitment
 TSplot(yrlab,out@TSdata$Rec,xlab="Historical year",ylab="Recruitment")
 
 # SSB-Rec
 TSplot(x=out@TSdata$SSB[,2:nyears],y=out@TSdata$Rec[,2:nyears],
        xlab="Spawning biomass",ylab="Recruitment",mat=F,type='p')
 
 F_FMSY<-FM/matrix(out@Ref$FMSY, nrow=nyears, ncol=nsim, byrow=TRUE)
 B_BMSY<-out@TSdata$SSB/matrix(out@Ref$SSBMSY, nrow=nsim, ncol=nyears)
 
 TSKplot(B_BMSY,t(F_FMSY),yrlab)
 
 # Age vulnerability
 maxage<-dim(out@AtAge$Select)[2]
 colors <- c("green","blue","grey45")
 for (x in 1:3) {
   Zvals <- t(out@AtAge$Select[x,,1:nyears])
   if(sd(Zvals, na.rm=TRUE) != 0) {
     if (x==1)contour(x=yrlab,y=1:maxage,z=Zvals,levels=c(0.25,0.75),col=colors[x],drawlabels=F,lwd=c(1,2))
     if (x!=1)contour(x=yrlab,y=1:maxage,z=Zvals,levels=c(0.25,0.75),col=colors[x],drawlabels=F, add=T,lwd=c(1,2))
   }
 }
 
 legend('topright',legend=c(paste("Simulation",1:3)),text.col=c("green","blue","grey45"),bty='n')
 legend('topleft',legend="Age vulnerability (0.25, 0.75)",bty='n')
 
 mtext("Historical year", 1, line = 2.5, cex = 1)
 mtext("Age", 2, line = 2.3, cex = 1)
 
 if (class(x) == 'OM') mtext(paste0("Time series plots for operating model ",OM@Name),3,outer=T,line= 0.2,font=2)
 
 return(invisible(out))
}

TSplot<-function(x,y,xlab=NA,ylab=NA,zeroy=T,incx=T,incy=T,type='l',mat=T){
  
  cols<-rep(makeTransparent(c("grey35","blue","orange","green")),100)
  nsim<-ncol(y)
  vlarg<-1e20
  
  rx<-range(x)
  ry<-range(y)
  if(zeroy)ry<-range(0,ry)
  
  rx1<-rx+c(-vlarg,vlarg)
  ry1<-ry+c(0,vlarg)
  
  plot(rx,ry,axes=F,col="white",xlab="",ylab="")
  polygon(rx1[c(1,1,2,2)],ry1[c(1,2,2,1)],col='grey94',border="grey94")
  xl<-pretty(x)
  abline(v=xl,col='white')
  yl<-pretty(seq(ry[1],ry[2],length.out=12))
  abline(h=yl,col='white')
  
  if(mat){
    matplot(x,t(y),type=type,col=cols,xlab="",ylab="",add=T)
  }else{
    if(type=='p')for(i in 1:nsim)points(x[,i],y[,i],col=cols[i],pch=19)
    if(type=='l')for(i in 1:nsim)lines(x[,i],y[,i],col=cols[i])
  }
  
  axis(1,c(-vlarg,vlarg),c(-vlarg,vlarg))
  axis(2,c(-vlarg,vlarg),c(-vlarg,vlarg))
  
  if(incx)axis(1,xl,xl)
  if(incy)axis(2,yl,yl)
  
  if(!is.na(xlab))mtext(xlab,1,line=2.5)
  if(!is.na(ylab))mtext(ylab,2,line=2.5)
  
}




TSKplot<-function(B_BMSY,F_FMSY,yrlab,maxsim=10){
  
  nyears<-ncol(B_BMSY)
  nsim<-nrow(B_BMSY)
  cex.leg<-0.9
  vlarg<-1e20
  
  FMSYr <- quantile(F_FMSY, c(0.001, 0.9), na.rm = T)
  BMSYr <- quantile(B_BMSY, c(0.001, 0.975), na.rm = T)
  
  colsse <- cols<-rainbow(nyears, start = 0.63, end = 0.95)[1:nyears]
  colsse <- makeTransparent(cols, 95)
  
  XLim <- c(0, 3)
  YLim <- c(0, 2.5)
  YLim[2]<-max(YLim[2],F_FMSY*1.1)
  
  plot(c(B_BMSY[1, 1], B_BMSY[1, 2]), c(F_FMSY[1,1], F_FMSY[1, 2]), xlim = XLim, ylim = YLim, col = colsse[1], type = "l", axes = FALSE,xlab="",ylab="")
  
  polygon(c(-100,100,100,-100),c(-100,-100,100,100),col='grey94',border="grey94")
  
  axis(1,c(-vlarg,vlarg),c(-vlarg,vlarg))
  axis(2,c(-vlarg,vlarg),c(-vlarg,vlarg))
  
  axis(side = 2, labels = TRUE, las = 1)
  axis(side = 1, labels = TRUE)
  
  OO <- round(sum(B_BMSY[,nyears] < 1 & F_FMSY[,nyears] > 1, na.rm = T)/nsim * 100, 1)
  OU <- round(sum(B_BMSY[,nyears] > 1 & F_FMSY[,nyears] > 1, na.rm = T)/nsim * 100, 1)
  UO <- round(sum(B_BMSY[,nyears] < 1 & F_FMSY[,nyears] < 1, na.rm = T)/nsim * 100, 1)
  UU <- round(sum(B_BMSY[,nyears] > 1 & F_FMSY[,nyears] < 1, na.rm = T)/nsim * 100, 1)
  
  abline(h = c(0,1), col = "white", lwd = 2)
  abline(v = c(0,1), col = "white", lwd = 2)
  
  y <- 1:(nyears - 1)
  y1 <- y + 1
  x0 <- as.vector(B_BMSY[, y])
  x1 <- as.vector(B_BMSY[, y1])
  y0 <- as.vector(F_FMSY[, y])
  y1 <- as.vector(F_FMSY[, y1])
  segments(x0, y0, x1, y1, col = rep(colsse,each=nsim))
  
  rng <- 1:min(maxsim, nsim)
  points(B_BMSY[rng, 1], F_FMSY[rng, 1], pch = 19, cex = 0.8, col = colsse[1])
  points(B_BMSY[rng, nyears],F_FMSY[rng,nyears], pch = 19, cex = 0.8, col = colsse[nyears])
  
  text(B_BMSY[1, ],F_FMSY[1,],yrlab,cex=0.7,font=2)  
  
  legend("right", legend=c(yrlab[1], yrlab[nyears]), bty = "n", text.col = c(cols[1],  cols[nyears]), pch = 19, col = c(cols[1], cols[nyears]))
  
  legend("topleft", paste(OO, "%", sep = ""), bty = "n", text.font = 2, cex = cex.leg,text.col='grey39')
  legend("topright", paste(OU, "%", sep = ""), bty = "n", text.font = 2,  cex = cex.leg,text.col='grey39')
  legend("bottomleft", paste(UO, "%", sep = ""), bty = "n", text.font = 2, cex = cex.leg,text.col='grey39')
  legend("bottomright", paste(UU, "%", sep = ""), bty = "n", text.font = 2, cex = cex.leg,text.col='grey39')
  
  mtext(expression(B/B[MSY]), 1, line = 2.5, cex = 1)
  mtext(expression(F/F[MSY]), 2, line = 2, cex = 1)
  
}
