#' @method plot Obs
#' @export
plot.Obs <- function(x, ...)  plotObs(x, ...)

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
  
  Obs <- x
  if (class(Obs) == "OM") {
    if (is.finite(Obs@nyears)) nyears <- Obs@nyears
    if (is.finite(Obs@nsim)) nsim <- Obs@nsim	
    Obs <- SubOM(Obs,"Obs")
  }
  
  nsamp <- 3
  its <- sample(1:nsim, nsamp)
  
  # === Sample Observation Model Parameters ====
  ObsPars <- SampleObsPars(Obs, nsim)
  # Assign Obs pars to function environment
  for (X in 1:length(ObsPars)) assign(names(ObsPars)[X], ObsPars[[X]])
  

  # === Non time series ==================================================================== 
  cex.main <- 0.5
  op <- par(mfrow=c(4,4),mai=c(0.6,0.6,0.25,0.01),omi=c(0.01,0.01,0.4,0.01))
  
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

  
  # -----
  
  if (!is.na(Obs@Name)) mtext(paste0("Observation time series plots for observation object ",Obs@Name),3,outer=T,line= 0.7,font=2)
  if (is.na(Obs@Name)) mtext(paste0("Observation time series plots for observation object "),3,outer=T,line= 0.7,font=2)
  
  on.exit(par(op))
  
}


#' @method plot Imp
#' @export
plot.Imp <- function(x, ...)  plotImp(x, ...)

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
  
  on.exit(par(op))

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


