# === OM specification using iSCAM stock assessment ====================

#' Reads MLE estimates from iSCAM file structure into an operating model 
#'
#' @description A function that uses the file location of a fitted SS3 model including input files to population the various slots of an operating model with MLE parameter estimates 
#' @param iscamdir A folder with Stock Synthesis input and output files in it
#' @param nsim The number of simulations to take for parameters with uncertainty (for OM@cpars custom parameters)
#' @param Name The name of the operating model
#' @param Source Reference to assessment documentation e.g. a url
#' @param Author Who did the assessment
#' @author T. Carruthers
#' @export iSCAM2DLM
iSCAM2DLM<-function(iSCAMdir,nsim=48,proyears=50,Name=NULL,Source="No source provided",
                 Author="No author provided"){
  
  message("-- Using function SS_output of package r4ss to extract data from SS file structure --")

 # replist <- SS_output(dir=SSdir,covar=F,ncols=1000,printstats=printstats,verbose=verbose)
  #replistB <- SS_output(dir="F:/Base",covar=F,ncols=1000,printstats=printstats,verbose=verbose)
  #replistY <- SS_output(dir="F:/Base3",covar=F,ncols=1000,printstats=printstats,verbose=verbose)

  message("-- End of r4ss operations --")
  
  OM<-new('OM',testOM,Generic_fleet,Generic_obs)
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
  OM@Msd<-OM@Ksd<-OM@Linfsd<-OM@Mgrad<-OM@Kgrad<-OM@Linfgrad<-OM@recgrad<-c(0,0)
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
    
  OM@Source<-paste0(Source,". Author: ",Author,".")
  
  # Maturity --------------------------------------
  
  Mat_age<-growdat$Age_Mat
  Len_age<-growdat$Len_Mid
  
  # Currently using linear interpolation of mat vs len, is robust and very close to true logistic model predictions
  
  L50<-LinInterp(Mat_age,Len_age,0.5)
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
  
  muFage<-as.vector(apply(F_at_age[,ceiling(ncol(F_at_age)*0.75):ncol(F_at_age)],1,mean))
  Vuln<-muFage/max(muFage,na.rm=T)
  
  OM@L5<-rep(LinInterp(Vuln,Len_age,0.05,ascending=T,zeroint=T),2)                            # not used if V is in cpars
  OM@LFS<-rep(Len_age[which.min((exp(Vuln)-exp(1.05))^2 * 1:length(Vuln))],2)  # not used if V is in cpars
  OM@Vmaxlen<-rep(mean(Vuln[(length(Vuln)-(nseas+1)):length(Vuln)],na.rm=T),2)  # not used if V is in cpars
  
  OM@isRel="FALSE" # these are real lengths not relative to length at 50% maturity
  
  # -- Recruitment -----------------------------------------------
  
  nrecs<-length(replist$recruit$dev)
  recdevs<-replist$recruit$dev[(nrecs-nyears+1):(nrecs-1)]# last year is mean recruitment
  #recdevs<-replist[length(replist$recruit$dev)-nyears)]
  #recdevs[is.na(recdevs)]<-0
  OM@AC<-rep(acf(recdevs)$acf[2,1,1],2)
  
  Perr<-array(NA,c(nsim,nyears+proyears))
  Perr[,1:(nyears-1)]<-matrix(rnorm(nsim*(nyears-1),rep(recdevs,each=nsim),0.2),nrow=nsim) # generate a bunch of simulations with uncertainty
  procsd<-apply(Perr,1,sd,na.rm=T)
  OM@Perr<-quantile(procsd,c(0.025,0.975)) # uniform range is a point estimate from assessment MLE
  procmu <- -0.5 * (procsd)^2  # adjusted log normal mean
  Perr[,nyears:(nyears+proyears)]<-matrix(rnorm(nsim*(proyears+1),rep(procmu,proyears+1),rep(procsd,proyears+1)),nrow=nsim)
  AC<-mean(OM@AC)
  for (y in nyears:(nyears + proyears)) Perr[, y] <- AC * Perr[, y - 1] +   Perr[, y] * (1 - AC * AC)^0.5  
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
  OM@Fsd<-quantile(apply((Find[,1:(nyears-1)]-Find[,2:nyears])/Find[,2:nyears],1,sd),c(0.05,0.95))
  
  OM@Period<-rep(NaN,2)
  OM@Amplitude<-rep(NaN,2)
  OM@EffYears<-1:nyears
  OM@EffLower<-Find[1,]
  OM@EffUpper<-Find[1,]
  OM@qinc<-c(0,0)
  OM@qcv<-OM@Fsd
  
  
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
  OM
 
}




