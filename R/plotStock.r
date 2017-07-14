
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
    hist(x, border='white',xlab="",col=col,axes=axes,main=main,breaks=breaks)
  }
}


#' @method plot Stock
#' @export
plot.Stock <- function(x, ...)  plotStock(x, ...)

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
  SampCpars <- list() # empty list 
  if (class(Stock) == "OM") {
    if (is.finite(Stock@nyears)) nyears <- Stock@nyears
	  if (is.finite(Stock@proyears)) proyears <- Stock@proyears
	  if (is.finite(Stock@nsim)) nsim <- Stock@nsim	
	  
	  if(length(Stock@cpars)>0){ # custom parameters exist - sample and write to list
	    ncparsim<-cparscheck(Stock@cpars)   # check each list object has the same length and if not stop and error report
	    SampCpars <- SampleCpars(Stock@cpars, nsim) 
	  }
	  Stock <- SubOM(Stock)
  }
  its <- sample(1:nsim, nsamp)
  
  
  # --- Sample Stock Parameters ----
  StockPars <- SampleStockPars(Stock, nsim, nyears, proyears, SampCpars)
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
  
  # Row 1 -- Natural Mortality ---- 
  hist2(M, col=col, axes=FALSE, main="M", breaks=breaks)
  abline(v=M[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  hist2(Msd, col=col, axes=FALSE, main="Msd", breaks=breaks)
  abline(v=Msd[its], col=1:nsamp, lwd=lwd)
  axis(side=1) 
  hist2(Mgrad, col=col, axes=FALSE, main="Mgrad", breaks=breaks)
  abline(v=Mgrad[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
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
  hist2(Linfgrad, col=col, axes=FALSE, main="Linfgrad", breaks=breaks)
  abline(v=Linfgrad[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
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
  hist2(Kgrad, col=col, axes=FALSE, main="Kgrad", breaks=breaks)
  abline(v=Kgrad[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
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
  
  matplot(t(Perr[its,]), type="l", bty="l", main="Rec Devs by Year", lwd=lwd, lty=1)
	
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
  hist2(dep, col=col, axes=FALSE, main="Depletion", breaks=breaks)
  abline(v=dep[its], col=1:nsamp, lwd=lwd)
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
  matplot(Ls, Mat_len, type="l", bty="l", main="Maturity-at-length", lwd=lwd, lty=1)

  matplot(t(Mat_age[its,]), type="l", bty="l", main="Maturity-at-age", lwd=lwd, lty=1, axes=FALSE, xlim=c(0, maxage))
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

  on.exit(par(op))
  
}

