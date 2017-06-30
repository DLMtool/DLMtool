#' Plot the vulnerability and retention curves 
#'
#' @param OM An object of class 'OM' 
#' 
#' @export
#'
plotSelect <- function(OM) {
  if (class(OM) != "OM") stop("Object must be class 'OM' ")
  nsim <- OM@nsim 
  years <- OM@nyears + OM@proyears
  yr.vert <- seq(1, years, length.out=4)
  StockPars <- SampleStockPars(OM)
  FleetPars <- SampleFleetPars(OM, Stock=StockPars)
  
  sim <- sample(1:nsim, 1)
  message("Plotting selection and retention curve for a random simulation")
  op <- par(mfcol=c(2,4), bty="l", las=1, mar=c(3,3,2,0), oma=c(3,3,1,1), xpd=NA)
  
  for (yr in yr.vert) {
    # plot vulnerability & selection at age
    plot(1:StockPars$maxage, FleetPars$V2[sim,, yr], type="l", ylim=c(0,1), lwd=2, 
         axes=FALSE, ylab="", xlab="")
    axis(side=1)
    mtext(side=1, "Age", line=2.5)
    if (yr == yr.vert[1]) {
      axis(side=2)
      mtext(side=2, "Vulnerabilty/Retention", las=3, line=3)
    }
    if (yr != yr.vert[1]) axis(side=2, labels=FALSE)
    title(paste("Year", yr))
    
    polygon(x=c(1:StockPars$maxage, rev(1:StockPars$maxage)), 
            y=c(FleetPars$V[sim,, yr], rev(FleetPars$retA[sim,, yr])), col="gray", border=FALSE)
    lines(1:StockPars$maxage, FleetPars$V[sim,, yr], col=2, lwd=2, lty=2)
    lines(1:StockPars$maxage, FleetPars$retA[sim,, yr], col=4, lwd=2, lty=3)
    
    if (yr == yr.vert[1]) {
      legend("bottomleft", legend = c("Vulnerability", "Realized Selection", "Retention"),
             lwd=2, col=c(1, 2, 4), bty="n", lty=c(1,2,3))
    }
    
    # plot vulnerability & selection at length
    plot(StockPars$CAL_binsmid, FleetPars$SLarray2[sim,, yr], type="l", ylim=c(0,1), lwd=2, 
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
    lines(StockPars$CAL_binsmid, FleetPars$SLarray[sim,, yr], col=2, lwd=2, lty=2)
    lines(StockPars$CAL_binsmid, FleetPars$retL[sim,, yr], col=4, lwd=2, lty=3)
    
  }
  
  on.exit(par(op))
  
}