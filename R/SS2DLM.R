# === OM specification using SS3 (Methot 2012) stock assessments ====================

#' Reads MLE estimates from Stock Synthesis file structure into an operating model using package r4ss
#'
#' @description A function that uses the file location of a fitted SS3 model including input files to population the various slots of an operating model with MLE parameter estimates 
#' @param SSdir A folder with Stock Synthesis input and output files in it
#' @param nsim The number of simulations to take for parameters with uncertainty (for OM@cpars custom parameters)
#' @param proyears The number of projection years for MSE
#' @param length_timestep The duration (in years) of each timestep in the model (if an quarterly model is used this is 0.25)
#' @param Name The name of the operating model
#' @param Source Reference to assessment documentation e.g. a url
#' @param Author Who did the assessment
#' @param printstats Should the r4ss function SS_output return info on data that was read in?
#' @param verbose Should the r4ss function SS_ouput return detailed messages?
#' @author T. Carruthers
#' @export SS2DLM
#' @importFrom r4ss SS_output
SS2DLM<-function(SSdir,nsim=48,proyears=50,length_timestep=NA,Name=NULL,Source="No source provided",
                 Author="No author provided",printstats=F,verbose=T){
  
  message("-- Using function SS_output of package r4ss to extract data from SS file structure --")
  
  replist <- r4ss::SS_output(dir=SSdir,covar=F,ncols=1000,printstats=printstats,verbose=verbose)
  #replistB <- SS_output(dir="F:/Base",covar=F,ncols=1000,printstats=printstats,verbose=verbose)
  #replistY <- SS_output(dir="F:/Base3",covar=F,ncols=1000,printstats=printstats,verbose=verbose)
  
  message("-- End of r4ss operations --")
  
  OM<-new('OM',DLMtool::Albacore, DLMtool::Generic_Fleet, DLMtool::Generic_Obs, DLMtool::Perfect_Imp)
  OM@nsim<-nsim
  OM@proyears<-proyears
  
  
  # The trouble is that DLMtool is an annual model so if the SS model is seasonal we need to aggregate overseasons
  # The reason for the code immediately below (and length_timestep) is that some SS models just assume quarterly timesteps with no way of knowing this (other than maxage and M possibly)!
  if(is.na(length_timestep)){
    
    if((replist$endyr-replist$startyr+1)>100){
      nseas<-1/replist$seasduration[1] # too many  years for industrial fishing so assumes a quarterly model
    }else{
      nseas=1
    }
    
  }else{
    nseas<-1/length_timestep
  }
  
  if(is.null(Name)){
    OM@Name=SSdir
  }else{
    OM@Name=Name
  }
  
  OM@nyears<-nyears<-(replist$endyr-replist$startyr+1)/nseas
  if(nseas==1){
    OM@CurrentYr<-replist$endyr
  }else{
    OM@CurrentYr<-nyears
  }
  
  yind<-(1:nyears)*nseas
  yind2<-rep(1:nyears,each=nseas)
  
  # === Stock parameters =========================================================================================================
  
  growdat<-getGpars(replist)
  OM@maxage<-maxage<-ceiling(length(growdat$Age)/nseas)
  aind<-rep(1:OM@maxage,each=nseas)[1:length(growdat$Age)]
  
  GP <- replist$Growth_Parameters
  
  if(nrow(GP)>1){
    print(paste("Warning:", nrow(GP),"different sets of growth parameters were reported they look like this:"))
    print(GP)
    print("Only the first row of values was used here")
    print("")
  }
  
  muLinf=GP$Linf[1]
  cvLinf=GP$CVmax[1]
  muK=GP$K[1]*nseas
  out<-negcorlogspace(muLinf,muK,cvLinf,nsim) # K and Linf negatively correlated 90% with identifcal CV to Linf
  Linf<-out[,1]
  K<-out[,2]
  OM@K<-quantile(K,c(0.025,0.975))
  OM@Linf<-quantile(Linf,c(0.025,0.975))
  OM@t0=rep(GP$A_a_L0[1],2) # t0 is not 
  OM@Msd<-OM@Ksd<-OM@Linfsd<-OM@Mgrad<-OM@Kgrad<-OM@Linfgrad<-c(0,0)
  L50=GP$Mat1[1]
  OM@a=GP$WtLen1[1]
  OM@b=GP$WtLen2[1]
  
  M<-aggregate(growdat$M,by=list(aind),mean)$x
  Wt<-aggregate(growdat$Wt_Mid,by=list(aind),mean)$x
  Wt_age<-array(rep(Wt,each=nsim)*trlnorm(nsim,1,cvLinf), dim = c(nsim, maxage, nyears)) 
  Mat<-aggregate(growdat$Age_Mat,by=list(aind),mean)$x
  
  # SO FAR: Wt_age K Linf  
  
  surv<-c(1,exp(cumsum(-c(M[2:maxage]))))# currently M and survival ignore age 0 fish
  OM@M<-rep(sum(M[2:maxage]*surv[2:maxage])/sum(surv[2:maxage]),2)
  SSB0<-as.numeric(replist$Dynamic_Bzero$SPB[replist$Dynamic_Bzero$Era=="VIRG"])
  SpR<-sum(Wt*Mat*surv)
  OM@R0<-rep(SSB0/SpR,2)
  
  SSB<-replist$recruit$spawn_bio
  rec<-replist$recruit$pred_recr
  SSBpR<-SSB[1]/rec[1]
  hs<-SRopt(nsim,SSB,rec,SSBpR,plot=F,type="BH")
  OM@h<-quantile(hs,c(0.025,0.975))
  OM@SRrel<-1 # This is the default 
  
  # SO FAR OM:     Name, nsim, proyears, nyears, maxage, R0, M, Msd, Mgrad, h, SRrel, Linf, K, t0, Ksd, Kgrad, Linfsd, Linfgrad, recgrad, a, b,
  # SO FAR cpars:  Wt_age, K Linf hs
  
  #SSBcur<-as.numeric(replist$Dynamic_Bzero[nrow(replist$Dynamic_Bzero),3])
  #othdep<-SSBcur/SSB0
  OM@D<-rep(replist$current_depletion,2)
  
  # Movement modelling ----------------------------
  
  if(nrow(replist$movement)>0){
    mov<-movdistil(replist$movement)
    vec<-rep(1/2,2)
    for(i in 1:200)vec<-vec%*%mov
    OM@Frac_area_1<-OM@Size_area_1<-rep(vec[1],2)
    OM@Prob_staying<-rep(mean(mov[cbind(1:2,1:2)]),2)
  }else{
    OM@Frac_area_1<-OM@Size_area_1<-rep(0.5,2)
    OM@Prob_staying<-rep(0.5,2)
  }
  
  OM@Source <-paste0(Source,". Author: ",Author,".")
  
  # Maturity --------------------------------------
  
  Mat_age<-growdat$Age_Mat
  if(max(Mat_age)>1)Mat_age<-as.numeric(Mat_age>1)
  
  Len_age<-growdat$Len_Mid
  
  # Currently using linear interpolation of mat vs len, is robust and very close to true logistic model predictions
  
  L50<-LinInterp(Mat_age,Len_age,0.5+1E-6)
  OM@L50<-rep(L50,2)
  
  L95<-LinInterp(Mat_age,Len_age,0.95)
  OM@L50_95<-rep(L95-L50,2)
  
  
  # Fleet parameters ============================================================================================================
  
  # Vulnerability --------------------------------------------
  
  ages<-growdat$Age
  cols<-match(ages,names(replist$Z_at_age))
  F_at_age=t(replist$Z_at_age[,cols]-replist$M_at_age[,cols])
  F_at_age[nrow(F_at_age),]<-F_at_age[nrow(F_at_age)-1,]# ad-hoc mirroring to deal with SS missing predicitons of F in terminal age
  Ftab<-cbind(expand.grid(1:dim(F_at_age)[1],1:dim(F_at_age)[2]),as.vector(F_at_age))
  
  if(nseas>1){
    sumF<-aggregate(Ftab[,3],by=list(aind[Ftab[,1]],Ftab[,2]),mean,na.rm=T)
    sumF<-aggregate(sumF[,3],by=list(sumF[,1],yind2[sumF[,2]]),sum,na.rm=T)
  }else{
    sumF<-Ftab
  }
  
  sumF<-sumF[sumF[,2]<nyears,] # generic solution: delete partial observation of final F predictions in seasonal model (last season of last year is NA)
  V <- array(NA, dim = c(nsim, maxage, nyears + proyears)) 
  V[,,1:(nyears-1)]<-rep(sumF[,3],each=nsim) # for some reason SS doesn't predict F in final year
  V[,,nyears:(nyears+proyears)]<-V[,,nyears-1]
  
  Find<-apply(V,c(1,3),max,na.rm=T) # get apical F
  
  ind<-as.matrix(expand.grid(1:nsim,1:maxage,1:(nyears+proyears)))
  V[ind]<-V[ind]/Find[ind[,c(1,3)]]
  
  # guess at length parameters # this is over ridden anyway
  
  muFage<-as.vector(apply(F_at_age,1,mean))
  Vuln<-muFage/max(muFage,na.rm=T)
  
  OM@L5<-rep(LinInterp(Vuln,Len_age,0.05,ascending=T,zeroint=T),2)                            # not used if V is in cpars
  OM@LFS<-rep(Len_age[which.min((exp(Vuln)-exp(1.05))^2 * 1:length(Vuln))],2)  # not used if V is in cpars
  OM@Vmaxlen<-rep(mean(Vuln[(length(Vuln)-(nseas+1)):length(Vuln)],na.rm=T),2)  # not used if V is in cpars
  
  OM@isRel="FALSE" # these are real lengths not relative to length at 50% maturity
  
  # -- Recruitment -----------------------------------------------
  
  rind<-rep(1:1000,each=nseas)[1:length(replist$recruit$dev)]
  recs<-aggregate(replist$recruit$dev,list(rind),mean,na.rm=T)$x
  recs<-recs[!is.na(recs)&recs!=0]
  nrecs<-length(recs)
  recdevs<-rep(0,nyears+maxage-1)
  recdevs[(nyears+maxage)-(nrecs:1)]<-recs# last year is mean recruitment
  recdevs<-recdevs[1:(nyears+maxage-2)]
  #recdevs<-replist[length(replist$recruit$dev)-nyears)]
  #recdevs[is.na(recdevs)]<-0
  OM@AC<-rep(acf(recdevs[!is.na(recdevs)],plot=F)$acf[2,1,1],2)
  
  Perr<-array(NA,c(nsim,nyears+proyears+maxage-1))
  Perr[,1:(maxage+nyears-2)]<-matrix(rnorm(nsim*(maxage+nyears-2),rep(recdevs,each=nsim),0.2),nrow=nsim) # generate a bunch of simulations with uncertainty
  procsd<-apply(Perr,1,sd,na.rm=T)
  OM@Perr<-quantile(procsd,c(0.025,0.975)) # uniform range is a point estimate from assessment MLE
  procmu <- -0.5 * (procsd)^2  # adjusted log normal mean
  Perr[,(maxage+nyears-2)+(1:(proyears+1))]<-matrix(rnorm(nsim*(proyears+1),rep(procmu,proyears+1),rep(procsd,proyears+1)),nrow=nsim)
  AC<-mean(OM@AC)
  for (y in (nyears-1):(nyears + proyears)) Perr[, y] <- AC * Perr[, y - 1] +   Perr[, y] * (1 - AC * AC)^0.5  
  Perr<-exp(Perr)
  
  
  # -- Fishing mortality rate index ----------------------------
  
  Find<-Find[,1:nyears] # is only historical years
  Find<-Find/apply(Find,1,mean,na.rm=T)
  
  # SO FAR OM:     Name, nsim, proyears, nyears, maxage, R0, M, Msd, Mgrad, h, SRrel, Linf, K, t0, Ksd, Kgrad, Linfsd, Linfgrad, recgrad, a, b,
  #                Size_area_1, Frac_area_1, Prob_Staying, Source, L50, L50_95, SelYears, AbsSelYears, L5, LFS, Vmaxlen, L5Lower, L5Upper,
  #                LFSLower, LFSUpper, VmaxLower, VmaxUpper, isRel,
  # SO FAR cpars:  Wt_age, K Linf hs, Perr, Find
  
  
  #plot(replist$cpue$Obs,replist$cpue$Exp)
  
  OM@Spat_targ<-rep(1,2)
  OM@Esd<-quantile(apply((Find[,1:(nyears-1)]-Find[,2:nyears])/Find[,2:nyears],1,sd),c(0.05,0.95))
  
  OM@Period<-rep(NaN,2)
  OM@Amplitude<-rep(NaN,2)
  OM@EffYears<-1:nyears
  OM@EffLower<-Find[1,]
  OM@EffUpper<-Find[1,]
  OM@qinc<-c(0,0)
  OM@qcv<-OM@Esd
  
  
  # Observation model parameters ==============================================================================
  
  # Index observations -------------------------------------------------------
  
  OM@Iobs<-rep(sd(log(replist$cpue$Obs)-log(replist$cpue$Exp)),2)
  getbeta<-function(beta,x,y)sum((y-x^beta)^2)
  OM@beta<-rep(optimize(getbeta,x=replist$cpue$Obs,y=replist$cpue$Exp,interval=c(0.1,10))$minimum,2)
  
  #F_FMSY<-replist$derived_quants[grep("F_",replist$derived_quants[,1]),2]
  #Fref<-F_FMSY[length(F_FMSY)]
  #Bref<-replist$Kobe[nrow(replist$Kobe),2]
  #OBJ<-obj$likelihoods_used[1,1]
  #gmax<-obj$maximum_gradient_component
  
  Wt_age2<-array(NA, dim = c(nsim, maxage, nyears+proyears))
  Wt_age2[,,1:nyears]<-Wt_age
  Wt_age2[,,nyears+1:proyears]<-rep(Wt_age[,,nyears],proyears)
  
  OM@cpars<-list(V=V,Perr=Perr,Wt_age=Wt_age2,K=K,Linf=Linf,hs=hs,Find=Find)
  # attr(OM, "build") <- "SS2DLM"
  OM@seed <- 1
  OM
  
}


#' Linear interpolation of a y value at level xlev based on a vector x and y
#'
#' @param x A vector of x values
#' @param y A vector of y values (identical length to x)
#' @param xlev A the target level of x from which to guess y
#' @param ascending Are the the x values supposed to be ordered before interpolation
#' @param zeroint is there a zero-zero x-y intercept?
#' @author T. Carruthers
#' @export LinInterp
LinInterp<-function(x,y,xlev,ascending=F,zeroint=F){
  
  if(zeroint){
    x<-c(0,x)
    y<-c(0,y)
  } 
  
  if(ascending){
    cond<-(1:length(x))<which.max(x)
  }else{
    cond<-rep(TRUE,length(x))
  }
  
  close<-which.min((x[cond]-xlev)^2)
  
  ind<-c(close,close+(x[close]<xlev)*2-1)
  ind <- ind[ind <= length(x)]
  if (length(ind)==1) ind <- c(ind, ind-1)
  if (min(ind)==0) ind <- 1:2
  ind<-ind[order(ind)]

  pos<-(xlev-x[ind[1]])/(x[ind[2]]-x[ind[1]])
  max(1,y[ind[1]]+pos*(y[ind[2]]-y[ind[1]]))
  
}



#' A function that samples multivariate normal (logspace) variables 
#'
#' @param xmu The mean (normal space) of the first (x) variable 
#' @param ymu The mean (normal space) of the second (y) variable
#' @param xcv The coefficient of variation (normal space, log normal sd) of the x variable
#' @param nsim The number of random draws
#' @param cor The off-diagonal (symmetrical) correlation among x and y
#' @param ploty Whether a plot of the sampled variables should be produced
#' @author T. Carruthers
#' @export negcorlogspace
#' @importFrom mvtnorm rmvnorm
negcorlogspace<-function(xmu,ymu,xcv=0.1,nsim,cor=-0.9,ploty=F){
  
  varcov=matrix(c(1,cor,cor,1),nrow=2)
  out<-mvtnorm::rmvnorm(nsim,c(0,0),sigma=varcov)
  out<-out/rep(apply(out,2,sd)/xcv,each=nsim)
  out<-exp(out)
  out<-out/rep(apply(out,2,mean),each=nsim)
  out[,1]<-out[,1]*xmu
  out[,2]<-out[,2]*ymu
  if(ploty)plot(out[,1],out[,2])
  out  
  
}

#' Simplified a multi-area transition matrix into the best 2 x 2 representation
#'
#' @description A Function that takes a larger movement matrix, identifies the most parsimonious representation of 2 non-mixed areas and returns the final unfished movement matrix
#' @param movtab a table of estimated movements 
#' @author T. Carruthers
#' @export getGpars
movdistil<-function(movtab){

  nareas<-max(movtab$Source_area,movtab$Dest_area)
  mov<-array(0,c(nareas,nareas))
  mov[cbind(movtab$Source_area,movtab$Dest_area)]<-movtab[,ncol(movtab)]
  
  vec<-rep(1/nareas,nareas)
  for(i in 1:100)vec<-vec%*%mov
  endmov<-array(vec,c(nareas,nareas))*mov
  
  listy<-new('list')
  for(i in 1:nareas)listy[[i]]<-c(1,2)
  
  combins<-expand.grid(listy)
  combins<-combins[!apply(combins,1,sum)%in%c(nareas*1,nareas*2),]
  nc<-nrow(combins)/2
  combins<-combins[(nc+1):nrow(combins),]
  
  base<-cbind(expand.grid(1:nareas,1:nareas),as.vector(endmov))
  
  emt<-NULL
  out<-rep(NA,nc)
 
  for(i in 1:nc){
    
    vec<-combins[i,]
    vec<-c(-1,1)[as.numeric(vec)]
    out[i]<-sum((vec-(vec%*%mov))^2) 
 
  } 
  
  best<-as.numeric(combins[which.min(out),])
  
  aggdat<-cbind(expand.grid(best,best),as.vector(endmov))
  agg<-aggregate(aggdat[,3],by=list(aggdat[,1],aggdat[,2]),sum)
  newmov<-array(NA,c(2,2))
  newmov[as.matrix(agg[,1:2])]<-agg[,3]
  newmov/apply(newmov,1,sum)
  
}

#' Predict Beverton-Holt recruitment and return fit to S-R observations
#'
#' @description Internal function to optBH 
#' @param pars an initial guess at model parameters steepness and R0 
#' @param SSB 'observations' of spawning biomass
#' @param rec 'observations' (model predictions) of recruitment
#' @param SSBpR spawning stock biomass per recruit at unfished conditions
#' @param mode should fit or recruitment deviations be returned
#' @param plot should a plot of the model fit be produced? 
#' @author T. Carruthers
#' @export getBH
getBH<-function(pars,SSB,rec,SSBpR,mode=1,plot=F){
  
  h<-0.2+1/(1+exp(-pars[1]))*0.8
  R0<-exp(pars[2])
  
  recpred<-((0.8*R0*h*SSB)/(0.2*SSBpR*R0*(1-h)+(h-0.2)*SSB))
  
  if(plot){
    ord<-order(SSB)
    plot(SSB[ord],rec[ord],ylim=c(0,max(rec,R0)),xlim=c(0,max(SSB,R0*SSBpR)),xlab="",ylab="")
    SSB2<-seq(0,R0*SSBpR,length.out=500)
    recpred2<-((0.8*R0*h*SSB2)/(0.2*SSBpR*R0*(1-h)+(h-0.2)*SSB2))
    lines(SSB2,recpred2,col='blue')
    abline(v=c(0.2*R0*SSBpR,R0*SSBpR),lty=2,col='red')
    abline(h=c(R0,R0*h),lty=2,col='red')
    legend('topright',legend=c(paste0("h = ",round(h,3)),paste0("lnR0 = ",round(log(R0),3))),bty='n')
  }
  
  if(mode==1){
    #return(sum(((recpred-rec)/10000)^2))
    return(-sum(dnorm(log(recpred)-log(rec),0,0.5,log=T))-dnorm(pars[1],0,6,log=T)) # add a vague prior on h = 0.6
    #return(-sum(dnorm(recpred,rec,rec*0.5,log=T)))
  }else{
    return(rec-recpred)
  }
  
}

#' Wrapper for estimating stock recruitment parameters from resampled stock-recruitment data
#'
#' @param x position (currently redundant)
#' @param SSB 'observations' of spawning biomass
#' @param rec 'observations' (model predictions) of recruitment
#' @param SSBpR spawning stock biomass per recruit at unfished conditions
#' @param R0temp an initial guess at the level of unfished recruitment
#' @param pars an initial guess at model parameters steepness and R0 
#' @param frac the fraction of observations for resampling
#' @param plot should a plot of model fit be produced?
#' @author T. Carruthers
#' @export optBH
optBH<-function(x,SSB,rec,SSBpR,R0temp,pars,frac=0.5,plot=F){

  samp<-sample(1:length(SSB),size=ceiling(length(SSB)*frac),replace=F)
  opt<-optim(pars,getBH,method="L-BFGS-B",lower=c(-6.,log(R0temp/50)),upper=c(6.,log(R0temp*50)), hessian=T,
           SSB=SSB[samp],rec=rec[samp],SSBpR=SSBpR,mode=1, plot=F)
  if(plot)getBH(opt$par,SSB,rec,SSBpR,mode=2,plot=plot)
  0.2+1/(1+exp(-opt$par[1]))*0.8

}

#' Function that returns a stochastic estimate of steepness given observed stock recruitment data
#'
#' @param nsim number of samples of steepness to generate
#' @param SSB 'observations' of spawning biomass
#' @param rec 'observations' (model predictions) of recruitment
#' @param SSBpR spawning stock biomass per recruit at unfished conditions
#' @param plot should plots of model fit be produced?
#' @param type what type of stock recruitment curve is being fitted BH = Beverton-Holt
#' @author T. Carruthers
#' @export SRopt
SRopt<-function(nsim,SSB,rec,SSBpR,plot=F,type="BH"){
  
    R0temp<-rec[1] # have a guess at R0 for initializing nlm
    pars<-c(0,log(R0temp))
    SSBpR=SSB[1]/rec[1]
    if(sfIsRunning()){
      sfSapply(1:nsim,optBH,SSB=SSB,rec=rec,SSBpR=SSBpR,R0temp=R0temp,pars=pars,frac=0.8,plot=plot)
    }else{
      sapply(1:nsim,optBH,SSB=SSB,rec=rec,SSBpR=SSBpR,R0temp=R0temp,pars=pars,frac=0.8,plot=plot)
    }

}


#' Extracts growth parameters from a SS3 r4ss replist
#'
#' @param replist the list output of the r4ss SS_output function (a list of assessment inputs / outputs)
#' @param seas The reference season for the growth (not actually sure what this does yet) 
#' @author T. Carruthers
#' @export getGpars
getGpars<-function(replist, seas = 1) { # This is a rip-off of SSPlotBiology
  
  nseasons <- replist$nseasons
  growdat <- replist$endgrowth[replist$endgrowth$Seas == seas, ]
  growdat$CV_Beg <- growdat$SD_Beg/growdat$Len_Beg
  growthCVtype <- replist$growthCVtype
  biology <- replist$biology
  startyr <- replist$startyr
  FecType <- replist$FecType
  FecPar1name <- replist$FecPar1name
  FecPar2name <- replist$FecPar2name
  FecPar1 <- replist$FecPar1
  FecPar2 <- replist$FecPar2
  parameters <- replist$parameters
  nsexes <- replist$nsexes
  mainmorphs <- replist$mainmorphs
  accuage <- replist$accuage
  startyr <- replist$startyr
  endyr <- replist$endyr
  growthvaries <- replist$growthvaries
  growthseries <- replist$growthseries
  ageselex <- replist$ageselex
  MGparmAdj <- replist$MGparmAdj
  wtatage <- replist$wtatage
  Growth_Parameters <- replist$Growth_Parameters
  Grow_std <- replist$derived_quants[grep("Grow_std_", replist$derived_quants$LABEL), ]
  if (nrow(Grow_std) == 0) {
    Grow_std <- NULL
  }  else {
    Grow_std$pattern <- NA
    Grow_std$sex_char <- NA
    Grow_std$sex <- NA
    Grow_std$age <- NA
    for (irow in 1:nrow(Grow_std)) {
      tmp <- strsplit(Grow_std$LABEL[irow], split = "_")[[1]]
      Grow_std$pattern[irow] <- as.numeric(tmp[3])
      Grow_std$sex_char[irow] <- tmp[4]
      Grow_std$age[irow] <- as.numeric(tmp[6])
    }
    Grow_std$sex[Grow_std$sex_char == "Fem"] <- 1
    Grow_std$sex[Grow_std$sex_char == "Mal"] <- 2
  }
  if (!is.null(replist$wtatage_switch)) {
    wtatage_switch <- replist$wtatage_switch
  } else{ stop("SSplotBiology function doesn't match SS_output function. Update one or both functions.")
  }
  if (wtatage_switch) 
    cat("Note: this model uses the empirical weight-at-age input.\n", 
        "     Therefore many of the parametric biology quantities which are plotted\n", 
        "     are not used in the model.\n")
  if (!seas %in% 1:nseasons) 
    stop("'seas' input should be within 1:nseasons")
  
  if (length(mainmorphs) > nsexes) {
    cat("!Error with morph indexing in SSplotBiology function.\n", 
        " Code is not set up to handle multiple growth patterns or birth seasons.\n")
  }
  if (FecType == 1) {
    fec_ylab <- "Eggs per kg"
    FecX <- biology$Wt_len_F
    FecY <- FecPar1 + FecPar2 * FecX
  }
  
  growdatF <- growdat[growdat$Gender == 1 & growdat$Morph == 
                        mainmorphs[1], ]
  growdatF$Sd_Size <- growdatF$SD_Beg
  
  if (growthCVtype == "logSD=f(A)") {
    growdatF$high <- qlnorm(0.975, meanlog = log(growdatF$Len_Beg), 
                            sdlog = growdatF$Sd_Size)
    growdatF$low <- qlnorm(0.025, meanlog = log(growdatF$Len_Beg), 
                           sdlog = growdatF$Sd_Size)
  }  else {
    growdatF$high <- qnorm(0.975, mean = growdatF$Len_Beg, 
                           sd = growdatF$Sd_Size)
    growdatF$low <- qnorm(0.025, mean = growdatF$Len_Beg, 
                          sd = growdatF$Sd_Size)
  }
  if (nsexes > 1) {
    growdatM <- growdat[growdat$Gender == 2 & growdat$Morph == 
                          mainmorphs[2], ]
    xm <- growdatM$Age_Beg
    growdatM$Sd_Size <- growdatM$SD_Beg
    if (growthCVtype == "logSD=f(A)") {
      growdatM$high <- qlnorm(0.975, meanlog = log(growdatM$Len_Beg), 
                              sdlog = growdatM$Sd_Size)
      growdatM$low <- qlnorm(0.025, meanlog = log(growdatM$Len_Beg), 
                             sdlog = growdatM$Sd_Size)
    }    else {
      growdatM$high <- qnorm(0.975, mean = growdatM$Len_Beg, 
                             sd = growdatM$Sd_Size)
      growdatM$low <- qnorm(0.025, mean = growdatM$Len_Beg, 
                            sd = growdatM$Sd_Size)
    }
  }
  
  growdatF
  
}  



someplot<-function (replist, yrs = "all", Ftgt = NA, ylab = "Summary Fishing Mortality", 
          plot = TRUE, print = FALSE, plotdir = "default", verbose = TRUE, 
          uncertainty = TRUE, pwidth = 6.5, pheight = 5, punits = "in", 
          res = 300, ptsize = 10) 
{
  pngfun <- function(file, caption = NA) {
    png(filename = file, width = pwidth, height = pheight, 
        units = punits, res = res, pointsize = ptsize)
    plotinfo <- rbind(plotinfo, data.frame(file = file, caption = caption))
    return(plotinfo)
  }
  plotinfo <- NULL
  if (plotdir == "default") 
    plotdir <- replist$inputs$dir
  if (yrs[1] == "all") {
    yrs <- replist$startyr:replist$endyr
  }
  Ftot <- replist$derived_quants[match(paste("F_", yrs, sep = ""), 
                                       replist$derived_quants$LABEL), ]
  if (all(is.na(Ftot$Value))) {
    warning("Skipping SSplotSummaryF because no real values found in DERIVED_QUANTITIES\n", 
            "    Values with labels like F_2012 may not be real.\n")
    return()
  }
  Fmax <- max(c(Ftot$Value, Ftgt + 0.01), na.rm = TRUE)
  if (uncertainty) {
    uppFtot <- Ftot$Value + 1.96 * Ftot$StdDev
    lowFtot <- Ftot$Value - 1.96 * Ftot$StdDev
    Fmax <- max(c(uppFtot, Ftgt + 0.01), na.rm = TRUE)
  }
  plotfun <- function() {
    plot(0, type = "n", , xlab = "Year", ylab = ylab, xlim = range(yrs), 
         ylim = c(0, Fmax), cex.lab = 1, cex.axis = 1, cex = 0.7)
    abline(h = 0, col = "grey")
    if (uncertainty) 
      segments(as.numeric(substring(Ftot$LABEL, 3, 6)), 
               uppFtot, as.numeric(substring(Ftot$LABEL, 3, 
                                             6)), lowFtot, col = gray(0.5))
    points(as.numeric(substring(Ftot$LABEL, 3, 6)), Ftot$Value, 
           pch = 16, type = "p")
    abline(h = Ftgt, col = "red")
  }
  if (plot) 
    plotfun()
  if (print) {
    file <- file.path(plotdir, "ts_summaryF.png")
    caption <- "Summary F (definition of F depends on setting in starter.ss)"
    plotinfo <- pngfun(file = file, caption = caption)
    plotfun()
    dev.off()
    if (!is.null(plotinfo)) 
      plotinfo$category <- "Timeseries"
  }
  if (verbose) 
    cat("Plotting Summary F\n")
  return(invisible(plotinfo))
}



#' Reads data Stock Synthesis file structure into an data object using package r4ss
#'
#' @description A function that uses the file location of a fitted SS3 model including input files to population the various slots of an data object 
#' @param SSdir A folder with Stock Synthesis input and output files in it
#' @param Source Reference to assessment documentation e.g. a url
#' @param length_timestep The duration (in years) of each timestep in the model (if an quarterly model is used this is 0.25)
#' @param Name The name of the operating model
#' @param Author Who did the assessment
#' @param printstats Should the r4ss function SS_output return info on data that was read in?
#' @param verbose Should the r4ss function SS_ouput return detailed messages?
#' @author T. Carruthers
#' @export SS2Data
#' @importFrom r4ss SS_output
SS2Data<-function(SSdir,Source="No source provided",length_timestep=NA,Name="",
                  Author="No author provided",printstats=F,verbose=T){
  
  message("-- Using function SS_output of package r4ss to extract data from SS file structure --")
  
  replist <- r4ss::SS_output(dir=SSdir, covar=F, ncols=1000, printstats=printstats, verbose=verbose)
  #replistB <- SS_output(dir="F:/Base",covar=F,ncols=1000,printstats=printstats,verbose=verbose)
  #replistY <- SS_output(dir="F:/Base3",covar=F,ncols=1000,printstats=printstats,verbose=verbose)
  
  message("-- End of r4ss operations --")
  
  GP <- replist$Growth_Parameters
  
  if(nrow(GP)>1){
    print(paste("Warning:", nrow(GP),"different sets of growth parameters were reported they look like this:"))
    print(GP)
    print("Only the first row of values was used here")
    print("")
  }
  
  Data <- new("Data") 
  
  # The trouble is that DLMtool is an annual model so if the SS model is seasonal we need to aggregate overseasons
  # The reason for the code immediately below (and length_timestep) is that some SS models just assume quarterly timesteps with no way of knowing this (other than maxage and M possibly)!
  if(is.na(length_timestep)){
    
    if((replist$endyr-replist$startyr+1)>100){
      nseas<-1/replist$seasduration[1] # too many  years for industrial fishing so assumes a quarterly model
    }else{
      nseas=1
    }
    
  }else{
    nseas<-1/length_timestep
  }
  
  nyears<-(replist$endyr-replist$startyr+1)/nseas
  Data@Year=(1:nyears)
  
  if(is.null(Name)){
    Data@Name=SSdir
  }else{
    Data@Name=Name
  }
  
  yind<-(1:nyears)*nseas
  yind2<-rep(1:nyears,each=nseas)
  
  
  
  cat<-rep(0,length(replist$timeseries$"obs_cat:_1"))
  for(i in 1:length(replist$timeseries)){
    if(grepl("obs_cat",names(replist$timeseries)[i]))cat<-cat+replist$timeseries[[i]]
    
  }
  if(nseas>1){
    
    ind<-rep(1:(length(cat)/nseas),nseas)
    cat2<-aggregate(cat,by=list(ind),sum)
    cat3<-cat2[1:length(yind2),2]
    cat4<-aggregate(cat3,by=list(yind2),sum)
  }
  
  Data@Cat <- matrix(cat4[,2],nrow=1)
  SSB<-replist$recruit$spawn_bio[1:(nyears*nseas)]
  SSB<-aggregate(SSB,by=list(yind2),mean)[,2]
  Ind<-SSB
  Ind<-Ind/mean(Ind)
  Data@Ind <- matrix(Ind,nrow=1) 
  rec<-replist$recruit$pred_recr
  rec<-rec[1:(nyears*nseas)]
  Recind<- aggregate(rec,by=list(yind2),mean)[,2]
  Recind<-Recind/mean(Recind)
  Data@Rec<-matrix(Recind,nrow=1)
  Data@t<-nyears
  Data@AvC<-mean(Data@Cat)
  Data@Dt<-Data@Dep<-Ind[nyears]/Ind[1]
  
  growdat<-getGpars(replist)
  growdat<-getGpars(replist)
  Data@MaxAge<-maxage<-ceiling(length(growdat$Age)/nseas)
  aind<-rep(1:maxage,each=nseas)[1:length(growdat$Age)]
  M<-aggregate(growdat$M,by=list(aind),mean)$x
  Wt<-aggregate(growdat$Wt_Mid,by=list(aind),mean)$x
  Mat<-aggregate(growdat$Age_Mat,by=list(aind),mean)$x
  
  surv<-c(1,exp(cumsum(-c(M[2:maxage]))))# currently M and survival ignore age 0 fish
  Data@Mort<-sum(M[2:maxage]*surv[2:maxage])/sum(surv[2:maxage])
  
  SSB0<-replist$derived_quants[replist$derived_quants$LABEL=="SSB_Unfished",2]
  
  Data@Bref<-replist$derived_quants[replist$derived_quants$LABEL=="SSB_MSY",2]
  Data@Cref<-replist$derived_quants[replist$derived_quants$LABEL=="TotYield_MSY",2]
  FMSY<-replist$derived_quants[replist$derived_quants$LABEL=="Fstd_MSY",2]*nseas
  Data@Iref<-Data@Ind[1]*Data@Bref/SSB[1]
  
  Data@FMSY_M=FMSY/Data@Mort
  Data@BMSY_B0<-Data@Bref/SSB0
  
  Data@Bref<-replist$derived_quants[replist$derived_quants$LABEL=="SSB_MSY",2]
  
  Len_age<-growdat$Len_Mid
  Mat_age<-growdat$Age_Mat
  
  Data@L50<-LinInterp(Mat_age,Len_age,0.5+1E-6)
  Data@L95<-LinInterp(Mat_age,Len_age,0.95)
  
  
  ages<-growdat$Age
  cols<-match(ages,names(replist$Z_at_age))
  F_at_age=t(replist$Z_at_age[,cols]-replist$M_at_age[,cols])
  F_at_age[nrow(F_at_age),]<-F_at_age[nrow(F_at_age)-1,]# ad-hoc mirroring to deal with SS missing predicitons of F in terminal age
  Ftab<-cbind(expand.grid(1:dim(F_at_age)[1],1:dim(F_at_age)[2]),as.vector(F_at_age))
  
  if(nseas>1){
    sumF<-aggregate(Ftab[,3],by=list(aind[Ftab[,1]],Ftab[,2]),mean,na.rm=T)
    sumF<-aggregate(sumF[,3],by=list(sumF[,1],yind2[sumF[,2]]),sum,na.rm=T)
  }else{
    sumF<-Ftab
  }
  
  sumF<-sumF[sumF[,2]<nyears,] # generic solution: delete partial observation of final F predictions in seasonal model (last season of last year is NA)
  V <- array(0, dim = c(maxage, nyears)) 
  V[,1:(nyears-1)]<-sumF[,3] # for some reason SS doesn't predict F in final year
  
  
  Find<-apply(V,2,max,na.rm=T) # get apical F
  
  ind<-as.matrix(expand.grid(1:maxage,1:nyears))
  V[ind]<-V[ind]/Find[ind[,2]]
  
  # guess at length parameters # this is over ridden anyway
  ind<-((nyears-3)*nseas)+(0:((nseas*3)-1))
  muFage<-as.vector(apply(F_at_age[,ind],1,mean))
  Vuln<-muFage/max(muFage,na.rm=T)
  
  Data@LFC<-LinInterp(Vuln,Len_age,0.05,ascending=T,zeroint=T)                            # not used if V is in cpars
  Data@LFS<-Len_age[which.min((exp(Vuln)-exp(1.05))^2 * 1:length(Vuln))]  # not used if V is in cpars
  
  ages<-growdat$Age
  cols<-match(ages,names(replist$Z_at_age))
  
  
  cols<-match(ages,names(replist$catage))
  CAA=as.matrix(replist$catage[,cols])
  yr<-rep(replist$catage$Yr,dim(CAA)[2])
  age<-rep(ages,each=dim(CAA)[1])
  CAAagg<-aggregate(as.vector(CAA),by=list(yr,age),sum)
  CAAagg<-CAAagg[CAAagg[,2]!=0,]
  CAAagg<-CAAagg[CAAagg[,1]<=(nyears*nseas),]
  CAAagg<-CAAagg[CAAagg[,2]<=length(aind),]
  
  yind3<-yind2[CAAagg[,1]]
  aind3<-aind[CAAagg[,2]]
  
  CAA2<-aggregate(CAAagg[,3],by=list(yind3,aind3),sum)
  Data@CAA<-array(0,c(1,nyears,maxage))
  Data@CAA[as.matrix(cbind(rep(1,nrow(CAA2)),CAA2[,1:2]))]<-CAA2[,3]
  
  
  
  # CAL data
  Data@CAL_bins<-c(0,replist$lbins)
  
  Binno<-match(replist$lendbase$Bin,Data@CAL_bins)
  Yrno<-yind2[replist$lendbase$Yr]
  CALagg<-aggregate(replist$lendbase$Obs*replist$lendbase$N,by=list(Yrno,Binno),sum)
  Data@CAL<-array(0,c(1,nyears,length(Data@CAL_bins)))
  Data@CAL[as.matrix(cbind(rep(1,nrow(CALagg)),CALagg[,1:2]))]<-CALagg[,3]
  Data@CAL<-ceiling(Data@CAL)
  
  Data@ML<-matrix(apply(Data@CAL[1,,]*rep(Data@CAL_bins,each=nyears),1,mean),nrow=1)
  Data@ML[Data@ML==0]<-NA
  lcpos<-apply(Data@CAL[1,,],1,which.max)
  Data@Lc<-matrix(Data@CAL_bins[lcpos],nrow=1)
  Data@Lc[Data@Lc==0]<-NA
  
  Data@Lbar<-matrix(rep(NA,nyears),nrow=1)
  
  for(i in 1:nyears){
    if(!is.na(Data@Lc[i])){
      ind<-Data@CAL_bins>Data@Lc[i]
      Data@Lbar[1,i]<-mean(Data@CAL[1,i,ind]*Data@CAL_bins[ind])
    } 
  }
  
  Data@Abun<-Data@SpAbun<-SSB[nyears]
  
  Data@vbLinf=GP$Linf[1]
  Data@CV_vbLinf=GP$CVmax[1]
  Data@vbK=GP$K[1]*nseas
  Data@vbt0=GP$A_a_L0
  Data@wla=GP$WtLen1[1]
  Data@wlb=GP$WtLen2[1]
  SSBpR<-SSB[1]/rec[1]  
  hs<-mean(SRopt(100,SSB,rec,SSBpR,plot=F,type="BH"))
  Data@steep<-hs
  
  Data@Units<-""
  Data@Ref<-mean(replist$derived_quants[grepl("OFLCatch",replist$derived_quants$LABEL),2])
  Data@Ref_type = "OFL"
  Data@LHYear<-nyears
  Data@MPeff<-1
  Data@Log<-list(warning="This data was created by an alpha version of SS2Data and should not be trusted!")
  Data@Misc<-list(warning="This data was created by an alpha version of SS2Data and should not be trusted!")
  Data@Name<-"This data was created by an alpha version of SS2Data and should not be trusted!"
  
  Data@MPrec<-Data@Ref<-Data@Cat[1,nyears]
  #Data@Ref<-mean(replist$derived_quants[grepl("OFLCatch",replist$derived_quants$LABEL),2])
  Data@Ref_type = "Most recent annual catch"
  Data
  
}



library(r4ss)



