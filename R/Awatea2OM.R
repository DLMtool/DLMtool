# === OM specification using SS3 (Methot 2012) stock assessments ====================

#' Reads MCMC estimates from Awatea (Paul Starr) processed r file structure into an operating model
#'
#' @description A function that uses the file location of a fitted Awatea model post-processed into a set of rmd files
#' @param AwateaDir A folder with Awatea files
#' @param nsim The number of simulations to take for parameters with uncertainty (for OM@cpars custom parameters)
#' @param proyears The number of projection years for MSE
#' @param Name The name of the operating model
#' @param Source Reference to assessment documentation e.g. a url
#' @param Author Who did the assessment
#' @param verbose Should the r4ss function SS_ouput return detailed messages?
#' @author T. Carruthers
#' @export Awatea2OM
Awatea2OM<-function(AwateaDir,nsim=48,proyears=50,Name=NULL,Source="No source provided",
                 Author="No author provided",verbose=T){

  bmcmc_file<-list.files(AwateaDir)[grep("Bmcmc",list.files(AwateaDir))]
  Awateatxt<-paste0(AwateaDir,"/",list.files(AwateaDir)[grep(".txt",list.files(AwateaDir))])

  load(paste0(AwateaDir,"/",bmcmc_file))
  load(paste0(AwateaDir,"/currentMCMC.rda"))
  load(paste0(AwateaDir,"/currentMSY.rda"))
  load(paste0(AwateaDir,"/currentProj.rda"))
  load(paste0(AwateaDir,"/currentRes.rda"))

  bmcmc<-Bmcmc[[1]][[1]] # Second tier list is the outputs


  #  ---------- Create new OM and set ranges ------------------------------
  OM<-new('OM',DLMtool::Albacore, DLMtool::Generic_Fleet, DLMtool::Generic_Obs, DLMtool::Perfect_Imp)
  OM@nsim<-nsim
  OM@proyears<-proyears
  OM@nyears<-nyears<-ncol(currentMCMC$U)
  OM@CurrentYr<-as.integer(substr(names(currentMCMC$U)[nyears],1,4))
  OM@maxage<-maxage<-max(currentRes$Sel$Age)
  OM@Fdisc<-rep(0.5,2)
  OM@DR<-rep(0.1,2)
  OM@Mgrad<-OM@Kgrad<-OM@Linfgrad<-rep(0,2)


  # Select posterior samples
  nmcmc<-nrow(currentMCMC$P)
  if(nsim<nmcmc){
    samp<-sample(1:nmcmc,nsim)
  }else{
    samp<-sample(1:nmcmc,nsim,replace=T)
    message(paste0("nsim (",nsim,") > nmcmc (",nmcmc,"), sampling posterior with replacement"))
  }


  # Awatea inputs to growth parameters -----------------------------------
  out<-AwateaRead(Awateatxt)
  OM@a<-out$a[1]
  OM@b<-out$b[1]
  OM@Linf<-rep(out$Linf[1],2)
  OM@K<-rep(out$K[1],2)
  OM@t0<-rep(out$t0[1],2)
  OM@LenCV<-rep(out$Lsd[1]/out$Linf[1],2)
  OM@isRel<-"FALSE"
  OM@L50<-2/3*OM@Linf

  # ------------ Custom Parameters --------------------------------------

  # ----- R0 ------------------------------------------
  R0<-currentMCMC$P$R_0[samp]
  OM@R0<-range(R0)

  # ----- Steepness -----------------------------------
  hs<-currentMCMC$P$h[samp]
  OM@h<-range(hs)

  # ----- Depletion -----------------------------------
  B0<-bmcmc$B0.MCMC[samp]
  D<-currentMCMC$B[samp,nyears]/B0
  OM@D<-range(D)

  # ----- F index -------------------------------------
  Find<-(-log(1-currentMCMC$U[samp,]))
  Find<-Find/apply(Find,1,mean)
  Find<-as.matrix(Find)
  OM@EffLower<-apply(Find,2,quantile,0.1)
  OM@EffUpper<-apply(Find,2,quantile,0.9)
  OM@EffYears<-1:nyears

  # ----- M -------------------------------------------
  M<-currentMCMC$P[samp,grep("M",names(currentMCMC$P))[1]]
  OM@M<-range(M)

  # ----- Process error recruitment calcs -------------
  Recs<-bmcmc$R.MCMC[samp,]
  B<-currentMCMC$B[samp,1:ncol(Recs)]
  SSBpR<-B0/R0
  predR<-(4*R0 * hs * B)/(SSBpR * R0 * (1-hs) + (5*hs-1) * B)

  inInd<-1:(maxage-1)
  histInd<-maxage+(1:nyears)-1
  projInd<-maxage+nyears+(0:(proyears-1))

  recdevs<-log(Recs/predR)
  ACfunc<-function(x)acf(x,plot=F)$acf[2,1,1]
  AC<-apply(recdevs,1,ACfunc)
  OM@AC<-range(AC)
  procsd<-apply(recdevs,1,sd,na.rm=T)
  OM@Perr<-range(procsd)
  OM@Perr<-quantile(procsd,c(0.025,0.975)) # uniform range is a point estimate from assessment MLE
  procmu <- -0.5 * (procsd)^2  # adjusted log normal mean
  Perr<-matrix(rnorm(nsim*(proyears+nyears+maxage-1),procmu,procsd),nrow=nsim)
  Perr[,histInd]<-as.matrix(recdevs)
  for (y in c(inInd[2:(maxage-1)],projInd)) Perr[, y] <- AC * Perr[, y - 1] +   Perr[, y] * (1 - AC * AC)^0.5
  Perr<-exp(Perr)


  # ----- Selectivity --------------------------------
  Sel<-currentRes$Sel

  if(length(unique(Sel$Sex))>1){
    Sel<-Sel[Sel$Sex=="Female",]
    message("Two sex assessment: using selectivity for females only")
  }

  GSel<-Sel[grep("Gear",Sel$Series),]
  Gnames<-unique(GSel$Series)

  if(length(GSel)>1){
    GSel<-GSel[GSel$Series==GSel[1],]
    message(paste("More than one gear selected, using selectivity for",Gnames[1]))
  }

  Selectivity<-GSel[,4]
  V <- array(rep(Selectivity,each=nsim), dim = c(nsim, maxage, nyears + proyears))

  # ----- Maturity ------------------------------------
  Mat<-Sel[Sel$Series=="Maturity",4]
  #plot(Mat,type='l',xlab="Age",ylab="Female fraction mature")
  Mat_age <- array(rep(Mat,each=nsim), dim=c(nsim, maxage, nyears+proyears))
  OM@L50<-2/3*OM@Linf
  OM@L50_95<-0.05*OM@Linf

  # ----- populate cpars ------------------------------
  OM@cpars<-list(M=M, R0=R0, D=D, Find=Find, hs=hs, Perr=Perr, Mat_age=Mat_age, V=V)

  OM@seed <- 1

  message("Awatea .txt input file (growth parameters) and process R output files read successfully")
  OM

}


AwateaRead<-function(Awateatxt,quiet=T){

  out<-list()
  uniquetext<-c("Bi-scalar","bii exponent","L-infinity","k of the von","t0 of the von","S.d. of length at age of oldest")
  param<-c("a","b","Linf","K","t0","Lsd")
  alltext<-readLines(Awateatxt)

  for(i in 1:length(param)){
    loc<-grep(uniquetext[i],alltext)
    out[[i]]<-scan(Awateatxt,skip=loc,nlines=1,quiet=quiet)
  }

  names(out)<-param
  out

}

