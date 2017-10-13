plotText <- function(OM, slots, RMDfile) {
  if (any(c("M", "h", "Linf", "L50", "D", "EffUpper", "qcv", "Vmaxlen", "DR") %in% slots)) {
    # slotstext <- paste("c(", paste(slots, sep=",", collapse = ","), ")")
    slotstext <- slots[slots %in% c("M", "h", "Linf", "L50", "D", "EffUpper", "qcv", "Vmaxlen", 
                                    "DR")]
    fig.asp <- switch(slotstext,
                      "M" = 1.5,
                      "h" = 1,
                      "Linf" =1.5,
                      "L50" = 1,
                      "D" = 0.5,
                      "EffUpper" = 1/2,
                      "qcv" = 1,
                      "Vmaxlen"=1,
                      "DR" = 0.75)
    cat("```{r plot.", slotstext, ", echo=FALSE, fig.asp=", fig.asp, "}\n", append=TRUE, file=RMDfile, sep="")
    cat("plotSlot(OM, Pars, slot='", slotstext, "')\n", append=TRUE, file=RMDfile, sep="")
    cat("```\n\n\n", append=TRUE, file=RMDfile, sep="")   
    
  } else if ('Obs' %in% slots) {
    cat("\n### Obs Plots\n", append=TRUE, file=RMDfile, sep="")
    cat("```{r plot.Obs, echo=FALSE, fig.asp=1}\n", append=TRUE, file=RMDfile, sep="")
    cat("plotObs(OM)\n", append=TRUE, file=RMDfile, sep="")
    cat("```\n\n", append=TRUE, file=RMDfile, sep="")   
    
  } else if ("Imp" %in% slots) {
    cat("\n### Imp Plots\n", append=TRUE, file=RMDfile, sep="")
    cat("```{r plot.Imp, echo=FALSE, fig.asp=1}\n", append=TRUE, file=RMDfile, sep="")
    cat("plotImp(OM)\n", append=TRUE, file=RMDfile, sep="")
    cat("```\n\n", append=TRUE, file=RMDfile, sep="") 
  } 
}


plotSlot <- function(OM, Pars, slot) {
  
  if (slot == 'M') plotM2(OM, Pars) 
  if (slot == "h") plotRec(OM, Pars) 
  if (slot == "Linf") plotGrowth(OM, Pars) 
  if (slot == "L50") plotMat(OM, Pars) 
  if (slot == "D") plotDep(OM, Pars) 
  if (slot == "EffUpper") plotEff(OM, Pars) 
  if (slot == "qcv") plotQcv(OM, Pars)
  if (slot == "Vmaxlen") plotSelHists(OM, Pars)
  if (slot == "DR") plotSelect(OM, Pars)

}

plotSelHists <- function(OM, Pars, nsamp=3, col="darkgray", 
                         breaks=10, lwd=2) {
 
  
  ncol <- 3
  m <- layout(matrix(c(1,2,3,
                       4,5,6,
                       7,0,0), ncol=ncol, byrow = TRUE))
  
  op <- par(mar = c(3, 2, 2, 1), oma=c(3,3,2,1), las=1, no.readonly = TRUE)
  on.exit(par(op))
  
  its <- sample(1:OM@nsim, nsamp)
  
  hist2(Pars$L5, col=col, axes=FALSE, main="L5", breaks=breaks)
  abline(v=Pars$L5[,its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  hist2(Pars$LFS, col=col, axes=FALSE, main="LFS", breaks=breaks)
  abline(v=Pars$LFS[,its], col=1:nsamp, lwd=lwd)
  axis(side=1)
 
  hist2(Pars$Vmaxlen, col=col, axes=FALSE, main="Vmaxlen", breaks=breaks)
  abline(v=Pars$Vmaxlen[,its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  
  hist2(Pars$LR5, col=col, axes=FALSE, main="LR5", breaks=breaks)
  abline(v=Pars$LR5[,its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  hist2(Pars$LFR, col=col, axes=FALSE, main="LFR", breaks=breaks)
  abline(v=Pars$LFR[,its], col=1:nsamp, lwd=lwd)
  axis(side=1)
   
  hist2(Pars$Rmaxlen, col=col, axes=FALSE, main="Rmaxlen", breaks=breaks)
  abline(v=Pars$Rmaxlen[,its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  hist2(Pars$DR, col=col, axes=FALSE, main="DR", breaks=breaks)
  abline(v=Pars$DR[,its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
}
plotQcv <- function(OM, Pars, nsamp=3, col="darkgray", 
                    breaks=10, lwd=2) {
  
  ncol <- 2
  m <- layout(matrix(c(1,2,
                       3,3), ncol=ncol, byrow = TRUE))
  
  op <- par(mar = c(3, 2, 2, 1), oma=c(3,3,2,1), las=1, no.readonly = TRUE)
  on.exit(par(op))
  
  its <- sample(1:OM@nsim, nsamp)

  hist2(Pars$qinc, col=col, axes=FALSE, main="qinc", breaks=breaks)
  abline(v=Pars$qinc[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  hist2(Pars$qcv, col=col, axes=FALSE, main="qcv", breaks=breaks)
  abline(v=Pars$qcv[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  matplot(t(Pars$qvar), type="l", bty="l", xlab="Projection Years",
          ylab="Catchability", las=1, ylim=c(0, max(Pars$qvar)), xpd=NA, cex.lab=1.5)
  matplot(t(Pars$qvar)[,its], type="l", bty="l", xlab="", add=TRUE,
          ylab="", las=1, ylim=c(0, max(Pars$qvar)), lwd=3, lty=1)
  
}

plotEff <- function(OM, Pars, nsamp=3, col="darkgray", 
                    breaks=10, lwd=2) {
  
  ncol <- 3
  m <- layout(matrix(c(1,1,2), ncol=ncol, byrow = TRUE))
  
  op <- par(mar = c(3, 2, 2, 1), oma=c(3,3,2,1), las=1, no.readonly = TRUE)
  on.exit(par(op))
  
  its <- sample(1:OM@nsim, nsamp)
  
  yrs <- (OM@CurrentYr -  OM@nyears + 1) : OM@CurrentYr
  matplot(yrs, t(Pars$Find), type="l", bty="l", xlab="Historical Years",
          ylab="Fishing Effort", las=1, ylim=c(0, max(Pars$Find)), cex.lab=1.5, xpd=NA)
  matplot(t(Pars$Find)[,its], type="l", bty="l", xlab="", add=TRUE,
          ylab="", las=1, ylim=c(0, max(Pars$Find)), lwd=3, lty=1)
  
  hist2(Pars$Esd, col=col, axes=FALSE, main="Esd", breaks=breaks)
  abline(v=Pars$Esd[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
}





plotDep <- function(OM, Pars=NULL, nsim=48, nyears=50, proyears=50, nsamp=3, col="darkgray", 
                    breaks=10, lwd=2) {
  if (class(OM) != "OM") stop("Must supply object of class 'OM'")
  
  if (is.finite(OM@nyears)) nyears <- OM@nyears
  if (is.finite(OM@proyears)) proyears <- OM@proyears
  if (is.finite(OM@nsim)) nsim <- OM@nsim	
  
  if (is.null(Pars)) {
    OM <- updateMSE(OM) # update and add missing slots with default values
    out<- runMSE(OM,Hist=T)
    Pars <- c(out$SampPars, out$TSdata, out$MSYs)
  }
  
  its <- sample(1:nsim, nsamp)
  
  ncol <- 2
  m <- layout(matrix(c(1,2), ncol=ncol, byrow = TRUE))
  
  op <- par(mar = c(3, 2, 2, 1), oma=c(2,3,2,1), las=1, no.readonly = TRUE)
  on.exit(par(op))
  
  ssb0 <- matrix(rep(Pars$SSB0, nyears), nrow=nyears, byrow=TRUE)
  dep <- Pars$SSB/ssb0
  ylim <- c(0, max(dep))
  matplot(dep,  type="l", bty="l", ylab="SB/SB0", xlab="Historical Years", xpd=NA, ylim=ylim)
  matplot(dep[, its],  type="l", bty="l", ylab="", xlab="", add=TRUE, lwd=4, col=1:nsamp, 
          lty=1, ylim=ylim)
  
  hist2(Pars$dep, col=col, axes=FALSE, main="Depletion (SB/SB0)", breaks=breaks)
  abline(v=Pars$dep[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  

  
  
}



plotMat <- function(OM, Pars=NULL, nsim=48, nyears=50, proyears=50, nsamp=3, col="darkgray", breaks=10, lwd=2) {
  if (class(OM) != "OM") stop("Must supply object of class 'OM'")
  
  if (is.finite(OM@nyears)) nyears <- OM@nyears
  if (is.finite(OM@proyears)) proyears <- OM@proyears
  if (is.finite(OM@nsim)) nsim <- OM@nsim	
  
  if (is.null(Pars)) {
    OM <- updateMSE(OM) # update and add missing slots with default values
    out<- runMSE(OM,Hist=T)
    Pars <- c(out$SampPars, out$TSdata, out$MSYs)
  }
  
  its <- sample(1:nsim, nsamp)
  
  ncol <- 2
  m <- layout(matrix(c(1,2,3,4), ncol=ncol, byrow = TRUE))
  
  op <- par(mar = c(3, 2, 2, 1), oma=c(2,3,2,1), las=1, no.readonly = TRUE)
  on.exit(par(op))
  
  
  # Maturity 
  hist2(Pars$L50, col=col, axes=FALSE, main="L50", breaks=breaks)
  abline(v=Pars$L50[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  hist2(Pars$L95, col=col, axes=FALSE, main="L95", breaks=breaks)
  abline(v=Pars$L95[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  
  slope <- log(19)/(Pars$L95-Pars$L50)
  Ls <- seq(0, to=max(Pars$Linf), length.out=200)
  
  Mat_len <- sapply(its, function(X)  plogis(Ls, Pars$L50[X], 1/slope[X]))
  matplot(Ls, Mat_len, type="l", bty="l", main="Maturity-at-length", lwd=lwd, lty=1, 
          ylab="Probability", xlab="Length", ylim=c(0,1), xpd=NA)
  
  matplot(t(Pars$Mat_age[its,,nyears]), type="l", bty="l", main="Maturity-at-age", lwd=lwd, 
          lty=1, axes=FALSE, xlim=c(0, Pars$maxage), ylab="", xlab="Age", ylim=c(0,1), xpd=NA)
  axis(side=1)
  axis(side=2, labels=FALSE)
  
}



plotGrowth <- function(OM, Pars=NULL, nsim=48, nyears=50, proyears=50, nsamp=3, col="darkgray", 
                       breaks=10, lwd=2) {
  if (class(OM) != "OM") stop("Must supply object of class 'OM'")
  
  if (is.finite(OM@nyears)) nyears <- OM@nyears
  if (is.finite(OM@proyears)) proyears <- OM@proyears
  if (is.finite(OM@nsim)) nsim <- OM@nsim	
  
  if (is.null(Pars)) {
    OM <- updateMSE(OM) # update and add missing slots with default values
    out<- runMSE(OM,Hist=T)
    Pars <- c(out$SampPars, out$TSdata, out$MSYs)
  }
  
  its <- sample(1:nsim, nsamp)
  
  ncol <- 9
  m <- layout(matrix(c(1,1, 2,2, 0, 3,3, 4,4,
                       5,5, 6,6, 0, 7,7, 8,8,
                       9,9,9,9,  0, 10,10,10,10,
                       11,11,11, 12,12,12, 13,13,13,
                       14,14,14, 15,15,15, 16,16,16), ncol=ncol, byrow = TRUE))
  
  op <- par(mar = c(5, 3, 3, 1), oma=c(4,6,2,1), las=0, no.readonly = TRUE)
  on.exit(par(op))
  
  # Histograms #### 
  # Linf
  hist2(Pars$Linf, col=col, axes=FALSE, main="Linf", breaks=breaks)
  abline(v=Pars$Linf[its], col=1:nsamp, lwd=lwd)
  axis(side=1)   
  
  # K 
  hist2(Pars$K, col=col, axes=FALSE, main="K", breaks=breaks)
  abline(v=Pars$K[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
  # t0 
  hist2(Pars$t0, col=col, axes=FALSE, main="t0", breaks=breaks)
  abline(v=Pars$t0[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
  # LenCV 
  hist2(Pars$LenCV, col=col, axes=FALSE, main="LenCV", breaks=breaks)
  abline(v=Pars$LenCV[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
  # Linfsd 
  hist2(Pars$Linfsd, col=col, axes=FALSE, main="Linfsd", breaks=breaks)
  abline(v=Pars$Linfsd[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
  # Linfgrad 
  hist2(Pars$Linfgrad, col=col, axes=FALSE, main="Linfgrad", breaks=breaks)
  abline(v=Pars$Linfgrad[its], col=1:nsamp, lwd=lwd)
  axis(side=1) 
  
  # Ksd 
  hist2(Pars$Ksd, col=col, axes=FALSE, main="Ksd", breaks=breaks)
  abline(v=Pars$Ksd[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
  # Kgrad 
  hist2(Pars$Kgrad, col=col, axes=FALSE, main="Kgrad", breaks=breaks)
  abline(v=Pars$Kgrad[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
 
  
  # By year 
  # Linf 
  matplot(t(Pars$Linfarray[its,]), type="l", bty="l", main="Linf by Year", lwd=lwd, lty=1, ylab="", xpd=NA) 
  # K 
  matplot(t(Pars$Karray[its,]), type="l", bty="l", main="K by Year", lwd=lwd, lty=1, ylab="", xpd=NA) 
  
  
  # Growth curves
  Len_age <- Pars$Len_age
  Wt_age <- Pars$Wt_age
  cex.lab <- 1.25
  fstYr <- Len_age[its,,1]
  curYr <- Len_age[its,,nyears]
  lstYr <- Len_age[its,,proyears+nyears]
  MaxL <- max(Len_age)
  matplot(t(fstYr), type="l", bty="l", main="First historical year", ylim=c(0, MaxL), 
          xlab="Age", ylab="Length", cex.lab=cex.lab, lwd=lwd, lty=1, xpd=NA)
  matplot(t(curYr), type="l", bty="l", main="Last historical year", ylim=c(0, MaxL),  
          axes=FALSE, xlab="Age", ylab="", cex.lab=cex.lab, lwd=lwd, lty=1, xpd=NA)
  axis(side=1)
  axis(side=2, labels=FALSE)  
  matplot(t(lstYr), type="l", bty="l", main="Last projected year", ylim=c(0, MaxL), axes=FALSE, 
          xlab="Age", ylab="", cex.lab=cex.lab, lwd=lwd, lty=1, xpd=NA)	
  axis(side=1)
  axis(side=2, labels=FALSE)  
  title("Sampled length-at-age curves", outer=TRUE, cex.main=2)
  
  fstYr <- Wt_age[its,,1]
  curYr <- Wt_age[its,,nyears]
  lstYr <- Wt_age[its,,proyears+nyears]
  MaxL <- max(Wt_age)
  matplot(t(fstYr), type="l", bty="l", main="First historical year", ylim=c(0, MaxL), 
          xlab="Age", ylab="Weight", cex.lab=cex.lab, lwd=lwd, lty=1, xpd=NA)
  matplot(t(curYr), type="l", bty="l", main="Last historical year", ylim=c(0, MaxL), 
          axes=FALSE, xlab="Age", ylab="", cex.lab=cex.lab, lwd=lwd, lty=1, xpd=NA)
  axis(side=1)
  axis(side=2, labels=FALSE)  
  matplot(t(lstYr), type="l", bty="l", main="Last projected year", ylim=c(0, MaxL), 
          axes=FALSE, xlab="Age", ylab="", cex.lab=cex.lab, lwd=lwd, lty=1, xpd=NA)	
  axis(side=1)
  axis(side=2, labels=FALSE)  
  title("Sampled length-at-age curves", outer=TRUE, cex.main=2)
  
}


plotRec <- function(OM, Pars=NULL, nsim=48, nyears=50, proyears=50, nsamp=3, col="darkgray", breaks=10, lwd=2) {
  if (class(OM) != "OM") stop("Must supply object of class 'OM'")
  
  if (is.finite(OM@nyears)) nyears <- OM@nyears
  if (is.finite(OM@proyears)) proyears <- OM@proyears
  if (is.finite(OM@nsim)) nsim <- OM@nsim	
    
  if (is.null(Pars)) {
    OM <- updateMSE(OM) # update and add missing slots with default values
    out<- runMSE(OM,Hist=T)
    Pars <- c(out$SampPars, out$TSdata, out$MSYs)
  }
  
  its <- sample(1:nsim, nsamp)
  
  ncol <- 3
  m <- layout(matrix(c(1, 2, 3, 
                       4, 5, 0,
                       6, 6, 6), ncol=ncol, byrow = TRUE))
  
  op <- par(mar = c(3, 1, 3, 1), oma=c(3,4,2,1), las=0, no.readonly = TRUE)
  on.exit(par(op))
  
  ## Histograms ####
  
  # h 
  hist2(Pars$hs, col=col, axes=FALSE, main="h", breaks=breaks)
  abline(v=Pars$hs[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
  # Perr 
  hist2(Pars$procsd, col=col, axes=FALSE, main="Perr", breaks=breaks)
  abline(v=Pars$procsd[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
  # AC 
  hist2(Pars$AC, col=col, axes=FALSE, main="AC", breaks=breaks)
  abline(v=Pars$AC[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
  # recgrad - not used 
  # hist2(Pars$recgrad, col=col, axes=FALSE, main="recgrad- currently not used", breaks=breaks)
  # abline(v=Pars$recgrad[its], col=1:nsamp, lwd=lwd)
  # axis(side=1)
  
  # Period - if used
  if (!is.null(Pars$Period)) {
    hist2(Pars$Period, col=col, axes=FALSE, main="Period", breaks=breaks)
    abline(v=Pars$Period[its], col=1:nsamp, lwd=lwd)
    axis(side=1)  
  } else {
    hist2(0, col=col, axes=FALSE, main="Period", breaks=breaks)
  }

  # Amplitude
  if (!is.null(Pars$Amplitude)) {
    hist2(Pars$Amplitude, col=col, axes=FALSE, main="Amplitude", breaks=breaks)
    abline(v=Pars$Amplitude[its], col=1:nsamp, lwd=lwd)
    axis(side=1)  
  } else {
    hist2(0, col=col, axes=FALSE, main="Amplitude", breaks=breaks)
  }

  
  # Recruitment
  matplot(t(Pars$Perr[its,]), type="l", bty="l", main="Rec Devs by Year", lwd=lwd, lty=1, ylab="")
  
 
}





# Plot Natural Mortality 
plotM2 <- function(OM, Pars=NULL, nsim=48, nyears=50, proyears=50, nsamp=3, col="darkgray", breaks=10, lwd=2) {
  if (class(OM) != "OM") stop("Must supply object of class 'OM'")
  
  if (is.finite(OM@nyears)) nyears <- OM@nyears
  if (is.finite(OM@proyears)) proyears <- OM@proyears
  if (is.finite(OM@nsim)) nsim <- OM@nsim	
  
  if (is.null(Pars)) {
    OM <- updateMSE(OM) # update and add missing slots with default values
    out<- runMSE(OM,Hist=T)
    Pars <- c(out$SampPars, out$TSdata, out$MSYs)
  }
  
  its <- sample(1:nsim, nsamp)
  
  
  ncol <- 4
  m <- layout(matrix(c(1, 2, 3, 4,
                       rep(5, 4),
                       6, 7, 8, 0,
                       9, 10, 11, 0,
                       12, 13, 14, 0), ncol=ncol, byrow = TRUE))
  
  op <- par(mar = c(3, 1, 3, 1), oma=c(3,4,2,1), las=1, no.readonly = TRUE)
  on.exit(par(op))
  
  ## Plot histograms of M parameters ####
  hist2(Pars$M, col=col, axes=FALSE, main="M", breaks=breaks)
  abline(v=Pars$M[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
  hist2(Pars$Mexp, col=col, axes=FALSE, main="Mexp", breaks=breaks)
  abline(v=Pars$Mexp[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  hist2(Pars$Msd, col=col, axes=FALSE, main="Msd", breaks=breaks)
  abline(v=Pars$Msd[its], col=1:nsamp, lwd=lwd)
  axis(side=1) 
  
  hist2(Pars$Mgrad, col=col, axes=FALSE, main="Mgrad", breaks=breaks)
  abline(v=Pars$Mgrad[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  # M by year 
  ylims <- range(Pars$M_ageArray[its,, ]) * c(0.95, 1.05)
  matplot(t(Pars$Marray[its,]), type="l", lty=1, bty="l", main="M by Year", lwd=lwd, ylab="M", ylim=ylims)
  abline(v=nyears, col="gray", lty=2)
  text(nyears, min(Pars$Marray[its,]), "Last historical year", pos=4, col="gray")
  
  # M at age 
  M_ageArray <- Pars$M_ageArray
  Len_age <- Pars$Len_age
  Wt_at_age <- Pars$Wt_age
  
  
  matplot(t(M_ageArray[its,,1]), type="l", lty=1, bty="l", lwd=lwd, ylim=ylims, ylab="M")
  mtext(side=3, "First historical year", cex=0.8, line=-1)
  mtext(side=1, "Age", line=2, cex=0.7)
  
  matplot(t(M_ageArray[its,,nyears]), type="l", lty=1, bty="l", main="M-at-age", lwd=lwd, ylim=ylims, axes=FALSE, ylab="")
  mtext(side=3, "Last historical year", cex=0.8, line=-1)
  axis(side=1)
  mtext(side=1, "Age", line=2, cex=0.7)
  
  matplot(t(M_ageArray[its,,nyears+proyears]), type="l", lty=1, bty="l", lwd=lwd, ylim=ylims, axes=FALSE, ylab="")
  mtext(side=3, "Last projected year", cex=0.8, line=-1)
  axis(side=1)
  mtext(side=1, "Age", line=2, cex=0.7)
  
  
  # M at length 
  xlims <- range(Len_age[its,, c(1, nyears, nyears+proyears)]) * c(0.95, 1.05)
  
  matplot(t(Len_age[its,,1]), t(M_ageArray[its,,1]), type="l", lty=1, bty="l", lwd=lwd, 
          ylim=ylims, xlim=xlims, ylab="M", xlab="")
  mtext(side=3, "First historical year", cex=0.8, line=-1)
  mtext(side=1, "Length", line=2, cex=0.7)
  
  matplot(t(Len_age[its,,nyears]), t(M_ageArray[its,,nyears]), type="l", lty=1, bty="l", 
          main="M-at-length", lwd=lwd, ylim=ylims, xlim=xlims, axes=FALSE, ylab="", xlab="")
  mtext(side=3, "Last historical year", cex=0.8, line=-1)
  axis(side=1)
  mtext(side=1, "Length", line=2, cex=0.7)
  
  matplot(t(Len_age[its,,nyears+proyears]), t(M_ageArray[its,,nyears+proyears]), type="l", 
          lty=1, bty="l", lwd=lwd, ylim=ylims, axes=FALSE, xlim=xlims, ylab="", xlab="")
  mtext(side=3, "Last projected year", cex=0.8, line=-1)
  axis(side=1)
  mtext(side=1, "Length", line=2, cex=0.7)
  
  
  # M at weight
  xlims <- range(Wt_at_age[its,, c(1, nyears, nyears+proyears)]) * c(0.95, 1.05)
  
  matplot(t(Wt_at_age[its,,1]), t(M_ageArray[its,,1]), type="l", lty=1, bty="l", lwd=lwd, 
          ylim=ylims, xlim=xlims, ylab="M", xlab="")
  mtext(side=3, "First historical year", cex=0.8, line=-1)
  mtext(side=1, "Weight", line=2, cex=0.7)
  
  matplot(t(Wt_at_age[its,,nyears]), t(M_ageArray[its,,nyears]), type="l", lty=1, bty="l", 
          main="M-at-weight", lwd=lwd, ylim=ylims, xlim=xlims, axes=FALSE, ylab="", xlab="")
  mtext(side=3, "Last historical year", cex=0.8, line=-1)
  axis(side=1)
  mtext(side=1, "Weight", line=2, cex=0.7)
  
  matplot(t(Wt_at_age[its,,nyears+proyears]), t(M_ageArray[its,,nyears+proyears]), type="l", 
          lty=1, bty="l", lwd=lwd, ylim=ylims, axes=FALSE, xlim=xlims, ylab="", xlab="")
  mtext(side=3, "Last projected year", cex=0.8, line=-1)
  axis(side=1)
  mtext(side=1, "Weight", line=2, cex=0.7)
  
  
}

