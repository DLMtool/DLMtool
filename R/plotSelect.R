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
  if (class(OM) != "OM") stop("Object must be class 'OM' ")
  nsim <- OM@nsim
  if (!is.na(sim)) OM@nsim <- nsim <- max(10, sim) # OM@nsim 
  years <- OM@nyears + OM@proyears
  yr.vert <- round(seq(1, years, length.out=pyears),0)
  if (is.null(Pars)) {
    stckPars <- SampleStockPars(OM)
    Pars <- c(stckPars, SampleFleetPars(OM, Stock=stckPars))
  }
  
  set.seed(OM@seed)
  if (is.na(sim)) sim <- sample(1:nsim, 1)
  
  if (pyears > 1) gr <- c(2, pyears)
  if (pyears == 1) gr <- c(1, 2)
  op <- par(mfcol=gr, bty="l", las=1, mar=c(3,3,2,0), oma=c(3,3,2,1), xpd=NA)
  
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
  on.exit(par(op))
  
  invisible(list(V=Pars$V, V2=Pars$V2 ,retA=Pars$retA, DR=Pars$DR, 
                 Fdisc=Pars$Fdisc, sim=sim))
  
  
}