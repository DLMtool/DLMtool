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
#' @method plot Obs 
#' @author T. Carruthers and A. Hordyk
#' @export 
plot.Obs <- function(x, nsim=500, nyears=50, 
                       col="darkgray", breaks=10, ...) {
  
  Obs <- x
  if (class(Obs) == "OM") {
    if (is.finite(Obs@nyears)) nyears <- Obs@nyears
    if (is.finite(Obs@nsim)) nyears <- Obs@nsim	
    Obs <- SubOM(Obs,"Obs")
  }
  
  nsamp=3
  its <- sample(1:nsim, nsamp)
  
 
  # Observation model Parameters 

  Csd <- runif(nsim, Obs@Cobs[1], Obs@Cobs[2])  # Sampled catch observation error (lognormal sd)
  Cbias <- rlnorm(nsim, mconv(1, Obs@Cbiascv), sdconv(1, Obs@Cbiascv))  # Sampled catch bias (log normal sd)
  CAA_nsamp <- ceiling(runif(nsim, Obs@CAA_nsamp[1], Obs@CAA_nsamp[2]))  # Number of catch-at-age observations
  CAA_ESS <- ceiling(runif(nsim, Obs@CAA_ESS[1], Obs@CAA_ESS[2]))  # Effective sample size
  CAL_nsamp <- runif(nsim, Obs@CAL_nsamp[1], Obs@CAL_nsamp[2])  # Observation error standard deviation for single catch at age by area
  CAL_ESS <- ceiling(runif(nsim, Obs@CAL_ESS[1], Obs@CAL_ESS[2]))  # Effective sample size
  CALcv <- runif(nsim, Obs@CALcv[1], Obs@CALcv[2])  # Observation error standard deviation for single catch at age by area
  betas <- exp(runif(nsim, log(Obs@beta[1]), log(Obs@beta[2])))  # the sampled hyperstability / hyperdepletion parameter beta>1 (hyperdepletion) beta<1 (hyperstability)
  Isd <- runif(nsim, Obs@Iobs[1], Obs@Iobs[2])  # Abundance index observation error (log normal sd)
  Derr <- runif(nsim, Obs@Dcv[1], Obs@Dcv[2])
  Dbias <- rlnorm(nsim, mconv(1, Obs@Dbiascv), sdconv(1, Obs@Dbiascv))  # sample of depletion bias
  Mbias <- rlnorm(nsim, mconv(1, Obs@Mcv), sdconv(1, Obs@Mcv))  # sample of M bias
  FMSY_Mbias <- rlnorm(nsim, mconv(1, Obs@FMSY_Mcv), sdconv(1, Obs@FMSY_Mcv))  # sample of FMSY/M bias
  
  lenMbias <- rlnorm(nsim, mconv(1, Obs@LenMcv), sdconv(1, Obs@LenMcv))  # sample of length at maturity bias - assume same error as age based maturity
  LFCbias <- rlnorm(nsim, mconv(1, Obs@LFCcv), sdconv(1, Obs@LFCcv))  # sample of length at first capture bias
  LFSbias <- rlnorm(nsim, mconv(1, Obs@LFScv), sdconv(1, Obs@LFScv))  # sample of length at full selection bias
  Aerr <- runif(nsim, Obs@Btcv[1], Obs@Btcv[2])
  Abias <- exp(runif(nsim, log(Obs@Btbias[1]), log(Obs@Btbias[2])))  #rlnorm(nsim,mconv(1,Obs@Btbiascv),sdconv(1,Obs@Btbiascv))    # smaple of current abundance bias
  Kbias <- rlnorm(nsim, mconv(1, Obs@Kcv), sdconv(1, Obs@Kcv))  # sample of von B. K parameter bias
  t0bias <- rlnorm(nsim, mconv(1, Obs@t0cv), sdconv(1, Obs@t0cv))  # sample of von B. t0 parameter bias
  Linfbias <- rlnorm(nsim, mconv(1, Obs@Linfcv), sdconv(1, Obs@Linfcv))  # sample of von B. maximum length bias
  Irefbias <- rlnorm(nsim, mconv(1, Obs@Irefcv), sdconv(1, Obs@Irefcv))  # sample of bias in reference (target) abundance index
  Crefbias <- rlnorm(nsim, mconv(1, Obs@Crefcv), sdconv(1, Obs@Crefcv))  # sample of bias in reference (target) catch index
  Brefbias <- rlnorm(nsim, mconv(1, Obs@Brefcv), sdconv(1, Obs@Brefcv))  # sample of bias in reference (target) biObsass index
  Recsd <- runif(nsim, Obs@Reccv[1], Obs@Reccv[2])  # Recruitment deviation 
  
  
  # === Non time series ==================================================================== 
  
  op <- par(mfrow=c(4,4),mai=c(0.6,0.6,0.25,0.01),omi=c(0.01,0.01,0.4,0.01))
  
  hist(CAA_nsamp,col=col,border="white", axes=FALSE, main="No. annual catch-at-age obs (CAA_samp)", breaks=breaks,cex.main=0.9)
  axis(side=1) 
  
  hist(CAA_ESS,col=col,border="white", axes=FALSE, main="Effective sample size CAA obs (CAA_ESS)", breaks=breaks,cex.main=0.9)
  axis(side=1) 
  
  hist(CAL_nsamp,col=col,border="white", axes=FALSE, main="No. annual catch-at-length obs (CAL_samp)", breaks=breaks,cex.main=0.9)
  axis(side=1) 
  
  hist(CAL_ESS,col=col,border="white", axes=FALSE, main="Effective sample size CAL obs (CAL_ESS)", breaks=breaks,cex.main=0.9)
  axis(side=1) 
  
  hist(Mbias,col=col,border="white", axes=FALSE, main="Natural mortality rate bias (Mbias)", breaks=breaks,cex.main=0.9)
  axis(side=1) 
  
  hist(FMSY_Mbias,col=col,border="white", axes=FALSE, main="FMSY/M bias (FMSY_Mbias)", breaks=breaks,cex.main=0.9)
  axis(side=1)
  
  hist(lenMbias,col=col,border="white", axes=FALSE, main="Bias in length at maturity (lenMbias)", breaks=breaks,cex.main=0.9)
  axis(side=1)
  
  hist(LFCbias,col=col,border="white", axes=FALSE, main="Bias in length at first capture (LFCbias)", breaks=breaks,cex.main=0.9)
  axis(side=1)
  
  hist(LFSbias,col=col,border="white", axes=FALSE, main="Bias in length at full selection (LFSbias)", breaks=breaks,cex.main=0.9)
  axis(side=1)
  
  hist(Kbias,col=col,border="white", axes=FALSE, main="Bias in von B. K (Kbias)", breaks=breaks,cex.main=0.9)
  axis(side=1)
  
  hist(t0bias,col=col,border="white", axes=FALSE, main="Bias in von B. t0 (t0bias)", breaks=breaks,cex.main=0.9)
  axis(side=1)
  
  hist(Linfbias,col=col,border="white", axes=FALSE, main="Bias in von B. Linf (Linfbias)", breaks=breaks,cex.main=0.9)
  axis(side=1)
  
  hist(Irefbias,col=col,border="white", axes=FALSE, main="Bias in index at MSY (Irefbias)", breaks=breaks,cex.main=0.9)
  axis(side=1)
  
  hist(Crefbias,col=col,border="white", axes=FALSE, main="Bias in MSY catch (Crefbias)", breaks=breaks,cex.main=0.9)
  axis(side=1)
  
  hist(Brefbias,col=col,border="white", axes=FALSE, main="Bias in MSY biomass (Brefbias)", breaks=breaks,cex.main=0.9)
  axis(side=1)
  
  hist(Recsd,col=col,border="white", axes=FALSE, main="Bias in recent recruitment strength (Recsd)", breaks=breaks,cex.main=0.9)
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
  
  hist(betas,col=col,border="white", axes=FALSE, main="Index hyper stability (betas)", breaks=breaks,cex.main=0.95)
  axis(side=1) 
  abline(v=1)
  abline(v=betas[its],col=makeTransparent(c("Black","Red","Green"),80),lwd=2)
  
  hist(Isd,col=col,border="white", axes=FALSE, main="Index error (Isd)", breaks=breaks,cex.main=0.95)
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
  
  mtext(paste0("Observation time series plots for observation object ",Obs@Name),3,outer=T,line= 0.7,font=2)
  
  on.exit(par(op))
  
}

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
#' @method plot Imp 
#' @author T. Carruthers and A. Hordyk
#' @export 
plot.Imp<-function(x,nsim=500, nyears=50, 
                   col="darkgray", breaks=10, ...){
  Imp <- x
  if (class(Imp) == "OM") {
    if (is.finite(Imp@nyears)) nyears <- Imp@nyears
    if (is.finite(Imp@nsim)) nyears <- Imp@nsim	
    Imp <- SubOM(Imp,"Imp")
  }

  nsamp=3
  its <- sample(1:nsim, nsamp)

  op <- par(mfrow=c(4,3),mai=c(0.6,0.6,0.2,0.01),omi=c(0.01,0.01,0.4,0.01))
   
  
  TACSD <- runif(nsim, Imp@TACSD[1], Imp@TACSD[2])  # Sampled TAC error (lognormal sd)
  TACFrac <- runif(nsim, Imp@TACFrac[1], Imp@TACFrac[2])  # Sampled TAC fraction (log normal sd)
  
  ESD <- runif(nsim, Imp@ESD[1], Imp@ESD[2])  # Sampled Effort error (lognormal sd)
  EFrac <- runif(nsim, Imp@EFrac[1], Imp@EFrac[2])  # Sampled Effort fraction (log normal sd)
  
  SizeLimSD<-runif(nsim,Imp@SizeLimSD[1],Imp@SizeLimSD[2])
  SizeLimFrac<-runif(nsim,Imp@SizeLimFrac[1],Imp@SizeLimFrac[2])
  
  DiscMort<-runif(nsim,Imp@DiscMort[1],Imp@DiscMort[2])
  
  ObsTSplot(TACFrac,TACSD,nyears,labs=c("Fraction of TAC (TACFrac)",
                                     "TAC error (TACSD)","TAC discrepancy for three samples",
                                     "TACFrac","TACSD"), breaks=breaks, its=its, nsamp=nsamp, col=col)
  
  ObsTSplot(EFrac,ESD,nyears,labs=c("Fraction of effort (EFrac)",
                                        "Effort error (ESD)","Effort discrepancy for three samples",
                                        "EFrac","ESD"), breaks=breaks, its=its, nsamp=nsamp, col=col)
  
  ObsTSplot(SizeLimFrac,SizeLimSD,nyears,labs=c("Fraction of Size Limit (SizeLimFrac)",
                                        "Size Limit error (SizeLimSD)","Size limit discrepancy for three samples",
                                        "SizeLimFrac","SizeLimSD"), breaks=breaks, its=its, nsamp=nsamp, col=col)
 
  hist(DiscMort,col=col,border="white", axes=FALSE, main="Discard Mortality rate (DiscMort)", breaks=breaks,cex.main=0.9)
  if(sd(DiscMort)>0.01)axis(side=1)
  
   
  mtext(paste0("Implementation error time series plots for implementation object ",Imp@Name),3,outer=T,line= 0.7,font=2)
  
  on.exit(par(op))

}



ObsTSplot<-function(Cbias,Csd,nyears,labs, breaks, its, nsamp, col){
  
  hist(Cbias,col=col,border="white", axes=FALSE, main=labs[1], breaks=breaks,cex.main=0.95)
  if(sd(Cbias)>0.01)axis(side=1) 
  abline(v=1)
  abline(v=Cbias[its],col=makeTransparent(c("Black","Red","Green"),80),lwd=2)
  hist(Csd,col=col,border="white", axes=FALSE, main=labs[2], breaks=breaks,cex.main=0.95)
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


