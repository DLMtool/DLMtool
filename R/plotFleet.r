#' @method plot Fleet
#' @export
plot.Fleet <- function(x, ...)  plotFleet(x, ...)

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

  Fleet <- updateMSE(x) # add missing slots with default values 
  SampCpars <- list() # empty list 
  nyears <- Fleet@nyears
  if (class(Fleet) == "OM") {
	  if (is.finite(Fleet@proyears)) proyears <- Fleet@proyears
	  if (is.finite(Fleet@nsim)) nsim <- Fleet@nsim	
	  if (length(Fleet@cpars) > 0) {
	    ncparsim <-cparscheck(Fleet@cpars)   # check each list object has the same length and if not stop and error report
	    SampCpars <- SampleCpars(Fleet@cpars, nsim) 
	  }
	  Stock <- SubOM(Fleet, "Stock")
	  Fleet <- SubOM(Fleet, "Fleet")
  }
  if (class(Stock) != "Stock") stop("Must include a Stock object", call.=FALSE)
  
  its <- sample(1:nsim, nsamp)  
  
  # --- Sample Stock Parameters ----
  StockPars <- SampleStockPars(Stock, nsim, nyears, proyears, SampCpars)
  # Assign Stock pars to function environment
  for (X in 1:length(StockPars)) assign(names(StockPars)[X], StockPars[[X]])
  
  # --- Sample Fleet Parameters ----
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
  
  title(paste("Fleet:", Fleet@Name, "Stock:", Stock@Name), outer=TRUE)
  title(paste("nyears =", nyears, "  proyears =", proyears, "  ", nsamp, "sampled iterations"), outer=TRUE, line=0)
    
  # om <- new("OM", Stock, Fleet, Perfect_Info, Perfect_Imp)
  # plotSelect(om)
  
  on.exit(par(op))	  
}
