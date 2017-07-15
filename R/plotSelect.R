#' Plot the vulnerability and retention curves 
#'
#' @param OM An object of class 'OM' 
#' @param pyears number of years to plot
#' @param sim the simulation to plot. default is NA to plot a random simulation
#' @param type plot type - line "l", point "p", or both "b"
#' @author A. Hordyk
#' @export
#'
plotSelect <- function(OM, pyears=4, sim=NA, type="l") {
  if (class(OM) != "OM") stop("Object must be class 'OM' ")
  nsim <- OM@nsim 
  years <- OM@nyears + OM@proyears
  yr.vert <- round(seq(1, years, length.out=pyears),0)
  StockPars <- SampleStockPars(OM)
  FleetPars <- SampleFleetPars(OM, Stock=StockPars)
  set.seed(OM@seed)
  if (is.na(sim)) sim <- sample(1:nsim, 1)
  message("Plotting selection and retention curve for simulation ", sim)
  
  if (pyears > 1) gr <- c(2, pyears)
  if (pyears == 1) gr <- c(1, 2)
  op <- par(mfcol=gr, bty="l", las=1, mar=c(3,3,2,0), oma=c(3,3,1,1), xpd=NA)
  
  for (yr in yr.vert) {
    # plot vulnerability & selection at age
    plot(1:StockPars$maxage, FleetPars$V2[sim,, yr], type=type, ylim=c(0,1), lwd=2, 
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
    
    polygon(x=c(1:StockPars$maxage, rev(1:StockPars$maxage)), 
            y=c(FleetPars$V[sim,, yr], rev(FleetPars$retA[sim,, yr])), col="gray", border=FALSE)
    lines(1:StockPars$maxage, FleetPars$V[sim,, yr], col=2, lwd=2, lty=2, type=type)
    lines(1:StockPars$maxage, FleetPars$retA[sim,, yr], col=4, lwd=2, lty=3, type=type)
    
    if (yr == yr.vert[pyears]) {
      minval <- min(c(FleetPars$V[sim,StockPars$maxage, yr],  FleetPars$retA[sim,StockPars$maxage, yr]))
      if (minval >= 0.5) loc <- "bottomright"
      if (minval < 0.5) loc <- "topright"
      legend(loc, legend = c("Vulnerability", "Realized Selection", "Retention"),
             lwd=2, col=c(1, 2, 4), bty="n", lty=c(1,2,3))
    }
    
    # plot vulnerability & selection at length
    plot(StockPars$CAL_binsmid, FleetPars$SLarray2[sim,, yr], type=type, ylim=c(0,1), lwd=2, 
         axes=FALSE, ylab="", xlab="")
    axis(side=1)
    mtext(side=1, "Length", line=2.5)
    if (yr == yr.vert[1]) {
      axis(side=2)
      mtext(side=2, "Vulnerabilty/Retention", las=3, line=3)
    }
    if (yr != yr.vert[1]) axis(side=2, labels=FALSE)
    
    polygon(x=c(StockPars$CAL_binsmid, rev(StockPars$CAL_binsmid)), 
            y=c(FleetPars$SLarray[sim,, yr], rev(FleetPars$retL[sim,, yr])), col="gray", border=FALSE)
    lines(StockPars$CAL_binsmid, FleetPars$SLarray[sim,, yr], col=2, lwd=2, lty=2,type=type)
    lines(StockPars$CAL_binsmid, FleetPars$retL[sim,, yr], col=4, lwd=2, lty=3, type=type)
    
  }
  
  on.exit(par(op))
  
  invisible(list(V=FleetPars$V, V2=FleetPars$V2 ,retA=FleetPars$retA, DR=FleetPars$DR, 
                 Fdisc=FleetPars$Fdisc, sim=sim))
  
  
}