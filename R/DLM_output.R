# DLM_output methods
TACfilter<-function(TAC) {
  TAC[TAC<0]<-NA    # Have to robustify due to R optmization problems.. work in progress.
  TAC[TAC>(mean(TAC,na.rm=T)+5*sd(TAC,na.rm=T))]<-NA  # remove very large TAC samples
  return(TAC)
}

prodPTF<-function(depletion,n,MSY){   # Pella-Tomlinson production function required for DB-SRA
   y<-(n^(n/(n-1)))/(n-1)
   MSY*y*depletion-MSY*y*depletion^n
}

fn<-function(n,BMSY_K){               # optimizer to find parameter n according to sampled BMSY/B0 (theta)
   thetapred<-n^(-1/(n-1))
   (BMSY_K-thetapred)^2
}

getn<-function(BMSY_K){               # wrapper for n finder
   optimize(fn,c(0.01,6),BMSY_K=BMSY_K)$minimum #get the optimum
}

gety<-function(n)  (n^(n/(n-1)))/(n-1) # More DBSRA code: get the y parameter for n

FMSYref<-function(x,DLM_data,reps=100)trlnorm(reps,DLM_data@OM$A[x]*(1-exp(-DLM_data@OM$FMSY[x])),0.01)
class(FMSYref)<-"DLM_output"

FMSYref50<-function(x,DLM_data,reps=100)trlnorm(reps,DLM_data@OM$A[x]*(1-exp(-DLM_data@OM$FMSY[x]))*0.5,0.01)
class(FMSYref50)<-"DLM_output"

FMSYref75<-function(x,DLM_data,reps=100)trlnorm(reps,DLM_data@OM$A[x]*(1-exp(-DLM_data@OM$FMSY[x]))*0.75,0.01)
class(FMSYref75)<-"DLM_output"

DynF<-function(x,DLM_data,yrsmth=10,gg=2,reps=100){
    
  dependencies="DLM_data@Year, DLM_data@Cat, DLM_data@Ind, DLM_data@Abun, DLM_data@Mort, DLM_data@FMSY_M"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  
  C_dat<-log(DLM_data@Cat[x,ind])
  B_dat<-log(DLM_data@Ind[x,ind]/DLM_data@Ind[x,ind[yrsmth]]*DLM_data@Abun[x])
  C_hist<-exp(predict(loess(C_dat~ind,degree=1)))
  B_hist<-exp(predict(loess(B_dat~ind,degree=1)))
  
  ind<-2:yrsmth
  ind1<-1:(yrsmth-1)
  SP_hist<-B_hist[ind]-B_hist[ind1]+C_hist[ind1]
  
  Frat<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])*trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])
  Flim <- matrix(NA, nrow=2, ncol=reps)
  Flim[1,] <- Frat * 0.5 
  Flim[2,] <- Frat * 2

  yind<-1:length(SP_hist)
  SP_mu<-predict(lm(SP_hist~yind),newdat=list(yind=length(SP_hist)+1))
  SP_se<-predict(lm(SP_hist~yind),newdat=list(yind=length(SP_hist)+1),se=T)$se.fit
  SP_new<-rnorm(reps,SP_mu,SP_se/2)
  Glm<-summary(lm(SP_hist~B_hist[ind1]))$coefficients[2,1:2] # plot(B_hist[ind1],SP_hist) # points(B_hist[ind1],SP_hist,col='green')
  G_new<-rnorm(reps,Glm[1],Glm[2]/2)
  #G_new[G_new>2*Frat]<-2*Frat[G_new<(2*Frat)]
  #G_new[G_new<(-2*Frat)]<--2*Frat[G_new<(-2*Frat)]
  G_new[G_new>0]<-G_new[G_new>0]*3
  newF<-Frat*exp(-G_new*gg)
  newF[newF<Flim[1]]<-Flim[1]
  newF[newF>Flim[2]]<-Flim[2]
  
  TAC<-newF*B_hist[yrsmth]
  TACfilter(TAC)
  
}
class(DynF)<-"DLM_output"

Fadapt<-function(x,DLM_data,reps=100,yrsmth=7,gg=1){
  
  dependencies="DLM_data@Year, DLM_data@Cat, DLM_data@Ind, DLM_data@Abun, DLM_data@Mort, DLM_data@FMSY_M"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  
  C_dat<-log(DLM_data@Cat[x,ind])
  B_dat<-log(DLM_data@Ind[x,ind]/DLM_data@Ind[x,ind[yrsmth]]*DLM_data@Abun[x])
  C_hist<-exp(predict(loess(C_dat~ind,degree=1)))
  B_hist<-exp(predict(loess(B_dat~ind,degree=1)))
  
  ind<-2:yrsmth
  ind1<-1:(yrsmth-1)
  SP_hist<-B_hist[ind]-B_hist[ind1]+C_hist[ind1]
  
  Frat<-DLM_data@Mort[x]*DLM_data@FMSY_M[x]
  Flim<-Frat*c(0.5,2)
  Flimr<-Flim[2]-Flim[1]
  
  yind<-1:length(SP_hist)
  SP_mu<-predict(lm(SP_hist~yind),newdat=list(yind=length(SP_hist)+1))
  SP_se<-predict(lm(SP_hist~yind),newdat=list(yind=length(SP_hist)+1),se=T)$se.fit
  SP_new<-rnorm(reps,SP_mu,SP_se/2)
  Glm<-summary(lm(SP_hist~B_hist[ind1]))$coefficients[2,1:2] # plot(B_hist[ind1],SP_hist) # points(B_hist[ind1],SP_hist,col='green')
  G_new<-rnorm(reps,Glm[1],Glm[2])
  
  Fold<-mean(C_hist/B_hist)
  
  if(Fold<Flim[1])Fmod1<-(-2)
  if(Fold>Flim[2])Fmod1<-2
  if(Fold>Flim[1]&Fold<Flim[2]){
    Ffrac<-(Fold-Flim[1])/Flimr
    Fmod1<-log(Ffrac/(1-Ffrac))
  }
  Fmod2<-Fmod1+gg*-G_new
  newF<-Flim[1]+(exp(Fmod2)/(1+exp(Fmod2)))*Flimr
  TAC<-newF*B_hist[yrsmth]
  TACfilter(TAC)
}
class(Fadapt)<-"DLM_output"

DepF<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@Year, DLM_data@Dep, DLM_data@Mort, DLM_data@FMSY_M, DLM_data@BMSY_B0"
  Frat<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])*trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])
  if (is.na(DLM_data@Dep[x]) | is.na(DLM_data@CV_Dep[x])) return(NA)
  depo<-max(0.01,min(0.99,DLM_data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
  Bt_K<-rbeta(reps*100,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
  Bt_K<-Bt_K[Bt_K>=0.01&Bt_K<=0.99][1:reps] # interval censor (0.01,0.99)  as in Dick and MacCall 2011
  adj<-Bt_K*(1-Bt_K)*4
  adj[Bt_K>0.5]<-1
  TAC<-Frat*DLM_data@Abun[x]*adj
  TACfilter(TAC)
}
class(DepF)<-"DLM_output"

Gcontrol<-function(x,DLM_data,reps=100,yrsmth=10,gg=2,glim=c(0.5,2)){
  dependencies="DLM_data@Year, DLM_data@Cat, DLM_data@Ind, DLM_data@Abun"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  C_dat<-log(DLM_data@Cat[x,ind])
  B_dat<-log(DLM_data@Ind[x,ind]/DLM_data@Ind[x,ind[yrsmth]]*DLM_data@Abun[x])
  C_hist<-exp(predict(loess(C_dat~ind,degree=1)))
  B_hist<-exp(predict(loess(B_dat~ind,degree=1)))
  ind<-2:yrsmth
  ind1<-1:(yrsmth-1)
  SP_hist<-B_hist[ind]-B_hist[ind1]+C_hist[ind1]
  yind<-1:length(SP_hist)
  SP_mu<-predict(lm(SP_hist~yind),newdat=list(yind=length(SP_hist)+1))
  SP_se<-predict(lm(SP_hist~yind),newdat=list(yind=length(SP_hist)+1),se=T)$se.fit
  SP_new<-rnorm(reps,SP_mu,SP_se/2)
  Glm<-summary(lm(SP_hist~B_hist[ind1]))$coefficients[2,1:2]
  G_new<-rnorm(reps,Glm[1],Glm[2]/2)

  TAC<-SP_new*(1-gg*G_new)
  TAC[TAC<glim[1]*C_hist[yrsmth]]<-glim[1]*C_hist[yrsmth]
  TAC[TAC>glim[2]*C_hist[yrsmth]]<-glim[2]*C_hist[yrsmth]
  
  #Carr<-cbind(array(rep(DLM_data@Cat[x,],each=reps),c(reps,length(DLM_data@Cat[x,]))),TAC)
  #Warr<-(DLM_data@Mort[x]*exp(-DLM_data@Mort[x]*(1:ncol(Carr))))[ncol(Carr):1]
  #Warr<-Warr/sum(Warr)
  #TAC<-apply(t(matrix(Warr,nrow=ncol(Carr),ncol=reps))*Carr,1,sum)
  TACfilter(TAC)
}
class(Gcontrol)<-"DLM_output"

Rcontrol<-function(x,DLM_data,reps=100,yrsmth=10,gg=2,glim=c(0.5,2)){
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0DLM_data@steep, DLM_data@CV_steep, DLM_data@MaxAge, DLM_data@Dep, DLM_data@CV_Dep, DLM_data@Cat, DLM_data@Ind"
  Mvec<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  Kvec<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Linfvec<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  t0vec<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  hvec<-trlnorm(reps,DLM_data@steep[x],DLM_data@CV_steep[x])
  rsamp<-getr(x,DLM_data,Mvec,Kvec,Linfvec,t0vec,hvec,maxage=DLM_data@MaxAge,r_reps=reps)

  depo<-max(0.01,min(0.99,DLM_data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
  Bt_K<-rbeta(100,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
  Bt_K<-Bt_K[Bt_K>0.01&Bt_K<0.99][1] # interval censor (0.01,0.99)  as in Dick and MacCall 2011

  G_new<-rsamp*(1-2*Bt_K)          # here is a big difference from SPHCR

  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  C_dat<-log(DLM_data@Cat[x,ind])
  B_dat<-log(DLM_data@Ind[x,ind]/DLM_data@Ind[x,ind[yrsmth]]*DLM_data@Abun[x])
  C_hist<-exp(predict(loess(C_dat~ind,degree=1)))
  B_hist<-exp(predict(loess(B_dat~ind,degree=1)))
  ind<-2:yrsmth
  ind1<-1:(yrsmth-1)
  SP_hist<-B_hist[ind]-B_hist[ind1]+C_hist[ind1]
  yind<-1:length(SP_hist)
  SP_mu<-predict(lm(SP_hist~yind),newdat=list(yind=length(SP_hist)+1))
  SP_se<-predict(lm(SP_hist~yind),newdat=list(yind=length(SP_hist)+1),se=T)$se.fit
  SP_new<-rnorm(reps,SP_mu,SP_se/2)

  TAC<-SP_new*(1-gg*G_new)
  TAC[TAC<glim[1]*C_hist[yrsmth]]<-glim[1]*C_hist[yrsmth]
  TAC[TAC>glim[2]*C_hist[yrsmth]]<-glim[2]*C_hist[yrsmth]
  
  #Carr<-cbind(array(rep(DLM_data@Cat[x,],each=reps),c(reps,length(DLM_data@Cat[x,]))),TAC)
  #Warr<-(DLM_data@Mort[x]*exp(-DLM_data@Mort[x]*(1:ncol(Carr))))[ncol(Carr):1]
  #Warr<-Warr/sum(Warr)
  #TAC<-apply(t(matrix(Warr,nrow=ncol(Carr),ncol=reps))*Carr,1,sum)
  TACfilter(TAC)

}
class(Rcontrol)<-"DLM_output"

Rcontrol2<-function(x,DLM_data,reps=100,yrsmth=10,gg=2,glim=c(0.5,2)){
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@steep, DLM_data@CV_steep, DLM_data@MaxAge, DLM_data@Dep, DLM_data@CV_Dep, DLM_data@Cat, DLM_data@Ind"
  Mvec<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  Kvec<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Linfvec<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  t0vec<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  hvec<-trlnorm(reps,DLM_data@steep[x],DLM_data@CV_steep[x])
  rsamp<-getr(x,DLM_data,Mvec,Kvec,Linfvec,t0vec,hvec,maxage=DLM_data@MaxAge,r_reps=reps)

  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  C_dat<-log(DLM_data@Cat[x,ind])
  B_dat<-log(DLM_data@Ind[x,ind]/DLM_data@Ind[x,ind[yrsmth]]*DLM_data@Abun[x])
  C_hist<-exp(predict(loess(C_dat~ind,degree=1)))
  B_hist<-exp(predict(loess(B_dat~ind,degree=1)))
  ind<-2:yrsmth
  ind1<-1:(yrsmth-1)
  SP_hist<-B_hist[ind]-B_hist[ind1]+C_hist[ind1]
  yind<-1:length(SP_hist)
  SP_mu<-predict(lm(SP_hist~yind),newdat=list(yind=length(SP_hist)+1))
  SP_se<-predict(lm(SP_hist~yind),newdat=list(yind=length(SP_hist)+1),se=T)$se.fit
  SP_new<-rnorm(reps,SP_mu,SP_se/2)
  SParr<-array(rep(SP_hist,each=reps),dim=c(reps,yrsmth-1))
  Barr<-array(rep(B_hist[ind],each=reps),dim=c(reps,yrsmth-1))
  rarr<-array(rep(rsamp,yrsmth-1),dim=c(reps,yrsmth-1))
  b2<-apply(SParr/Barr-rarr,1,sum)*apply(Barr,1,sum)/apply(Barr^2,1,sum)
  G_new<-rsamp-2*b2*B_hist[yrsmth]

  TAC<-SP_new*(1-gg*G_new)
  TAC[TAC<glim[1]*C_hist[yrsmth]]<-glim[1]*C_hist[yrsmth]
  TAC[TAC>glim[2]*C_hist[yrsmth]]<-glim[2]*C_hist[yrsmth]
  #Carr<-cbind(array(rep(DLM_data@Cat[x,],each=reps),c(reps,length(DLM_data@Cat[x,]))),TAC)
  #Warr<-(DLM_data@Mort[x]*exp(-DLM_data@Mort[x]*(1:ncol(Carr))))[ncol(Carr):1]
  #Warr<-Warr/sum(Warr)
  #TAC<-apply(t(matrix(Warr,nrow=ncol(Carr),ncol=reps))*Carr,1,sum)
  TACfilter(TAC)
}
class(Rcontrol2)<-"DLM_output"

GB_CC<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@Cref,DLM_data@Cat"
  Catrec<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  TAC<-trlnorm(reps,DLM_data@Cref[x],DLM_data@CV_Cref)
  TAC[TAC>(1.2*Catrec)]<-1.2*Catrec
  TAC[TAC<(0.8*Catrec)]<-0.8*Catrec
  TACfilter(TAC)
}
class(GB_CC)<-"DLM_output"

GB_slope<-function(x,DLM_data,reps=100,yrsmth=5,lambda=1){
  dependencies="DLM_data@Year, DLM_data@Cat, DLM_data@CV_Cat, DLM_data@Ind"
  Catrec<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  I_hist<-DLM_data@Ind[x,ind]
  yind<-1:yrsmth
  slppar<-summary(lm(I_hist~yind))$coefficients[2,1:2]
  Islp <-rnorm(reps,slppar[1],slppar[2])
  MuC<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  Cc<-rlnorm(reps,mconv(MuC,DLM_data@CV_Cat[x]*MuC),sdconv(MuC,DLM_data@CV_Cat[x]*MuC))
  TAC<-Cc*(1+lambda*Islp)
  TAC[TAC>(1.2*Catrec)]<-1.2*Catrec
  TAC[TAC<(0.8*Catrec)]<-0.8*Catrec
  TACfilter(TAC)
}
class(GB_slope)<-"DLM_output"

GB_target<-function(x,DLM_data,reps=100,w=0.5){
  dependencies="DLM_data@Cat, DLM_data@Cref, DLM_data@Iref, DLM_data@Ind"
  Catrec<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  TACtarg<-trlnorm(reps,DLM_data@Cref[x],DLM_data@CV_Cref)
  Itarg<-trlnorm(reps,DLM_data@Iref[x],DLM_data@CV_Iref)
  Iav<-mean(DLM_data@Ind[x,(length(DLM_data@Ind[x,])-4):length(DLM_data@Ind[x,])],na.rm=T)
  Irec<-mean(DLM_data@Ind[x,(length(DLM_data@Ind[x,])-3):length(DLM_data@Ind[x,])],na.rm=T)
  I0<-0.2*Iav
  TAC<-rep(NA,reps)
  if(Irec>I0)TAC<-TACtarg*(w+(1-w)*((Irec-I0)/(Itarg-I0)))
  if(Irec<I0)TAC<-TACtarg*w*(Irec/I0)^2
  TAC[TAC>(1.2*Catrec)]<-1.2*Catrec
  TAC[TAC<(0.8*Catrec)]<-0.8*Catrec
  TACfilter(TAC)
}
class(GB_target)<-"DLM_output"

CC1<-function(x,DLM_data,reps=100,yrsmth=5,xx=0){
  dependencies="DLM_data@Cat, DLM_data@CV_Cat"
  C_dat<-DLM_data@Cat[x,(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)]
  TAC<-(1-xx)*trlnorm(reps,mean(C_dat),DLM_data@CV_Cat/(yrsmth^0.5)) # mean catches over the interval
  TACfilter(TAC)
}  
class(CC1)<-"DLM_output"

CC4<-function(x,DLM_data,reps=100,yrsmth=5,xx=0.3){
  dependencies="DLM_data@Cat, DLM_data@CV_Cat"
  C_dat<-DLM_data@Cat[x,(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)]
  TAC<-(1-xx)*trlnorm(reps,mean(C_dat),DLM_data@CV_Cat/(yrsmth^0.5)) # mean catches over the interval
  TACfilter(TAC)
}  
class(CC4)<-"DLM_output"


LstepCC1<-function(x,DLM_data,reps=100,yrsmth=5,xx=0,stepsz=0.05,llim=c(0.96,0.98,1.05)){
  dependencies="DLM_data@Cat, DLM_data@CV_Cat, DLM_data@ML"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year) # recent 5 years
  ylast<-(DLM_data@LHYear-DLM_data@Year[1])+1 #last historical year
  ind2<-((ylast-(yrsmth-1)):ylast) # historical 5 pre-projection years
  ind3<-((ylast-(yrsmth*2-1)):ylast) # historical 10 pre-projection years
  C_dat<-DLM_data@Cat[x,ind2]
  if(length(DLM_data@Year)==ylast+1) {TACstar<-(1-xx)*trlnorm(reps,mean(C_dat),DLM_data@CV_Cat/(yrsmth^0.5))
  }else{TACstar<-rep(DLM_data@MPrec[x],reps)}
  step<-stepsz*TACstar
  Lrecent<-mean(DLM_data@ML[ind])
  Lave<-mean(DLM_data@ML[ind3])
  rat<-Lrecent/Lave
  if(rat<llim[1]){TAC<-TACstar-2*step
  }else if(rat<llim[2]){TAC<-TACstar-step
  }else if(rat>llim[3]){TAC<-TACstar+step
  }else{TAC<-TACstar
  }
  TACfilter(TAC)
}  
class(LstepCC1)<-"DLM_output"

LstepCC4<-function(x,DLM_data,reps=100,yrsmth=5,xx=0.3,stepsz=0.05,llim=c(0.96,0.98,1.05)){
  dependencies="DLM_data@Cat, DLM_data@CV_Cat, DLM_data@ML"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year) # recent 5 years
  ylast<-(DLM_data@LHYear-DLM_data@Year[1])+1 #last historical year
  ind2<-((ylast-(yrsmth-1)):ylast) # historical 5 pre-projection years
  ind3<-((ylast-(yrsmth*2-1)):ylast) # historical 10 pre-projection years
  C_dat<-DLM_data@Cat[x,ind2]
  if(length(DLM_data@Year)==ylast+1) {TACstar<-(1-xx)*trlnorm(reps,mean(C_dat),DLM_data@CV_Cat/(yrsmth^0.5))
  }else{TACstar<-rep(DLM_data@MPrec[x],reps)}
  step<-stepsz*TACstar
  Lrecent<-mean(DLM_data@ML[ind])
  Lave<-mean(DLM_data@ML[ind3])
  rat<-Lrecent/Lave
  if(rat<llim[1]){TAC<-TACstar-2*step
  }else if(rat<llim[2]){TAC<-TACstar-step
  }else if(rat>llim[3]){TAC<-TACstar+step
  }else{TAC<-TACstar
  }
  TACfilter(TAC)
}  
class(LstepCC4)<-"DLM_output"

Ltarget1<-function(x,DLM_data,reps=100,yrsmth=5,xx=0,xL=1.05){
  dependencies="DLM_data@Cat, DLM_data@CV_Cat, DLM_data@ML"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year) # recent 5 years
  ylast<-(DLM_data@LHYear-DLM_data@Year[1])+1 #last historical year
  ind2<-((ylast-(yrsmth-1)):ylast) # historical 5 pre-projection years
  ind3<-((ylast-(yrsmth*2-1)):ylast) # historical 10 pre-projection years
  C_dat<-DLM_data@Cat[x,ind2]
  TACstar<-(1-xx)*trlnorm(reps,mean(C_dat),DLM_data@CV_Cat/(yrsmth^0.5))
  Lrecent<-mean(DLM_data@ML[ind])
  Lave<-mean(DLM_data@ML[ind3])
  L0<-0.9*Lave
  Ltarget<-xL*Lave
  if(Lrecent>L0){TAC<-0.5*TACstar*(1+((Lrecent-L0)/(Ltarget-L0)))
  }else{TAC<-0.5*TACstar*(Lrecent/L0)^2                  
  }
  TACfilter(TAC)
}  
class(Ltarget1)<-"DLM_output"

Ltarget4<-function(x,DLM_data,reps=100,yrsmth=5,xx=0.2,xL=1.15){
  dependencies="DLM_data@Cat, DLM_data@CV_Cat, DLM_data@ML"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year) # recent 5 years
  ylast<-(DLM_data@LHYear-DLM_data@Year[1])+1 #last historical year
  ind2<-((ylast-(yrsmth-1)):ylast) # historical 5 pre-projection years
  ind3<-((ylast-(yrsmth*2-1)):ylast) # historical 10 pre-projection years
  C_dat<-DLM_data@Cat[x,ind2]
  TACstar<-(1-xx)*trlnorm(reps,mean(C_dat),DLM_data@CV_Cat/(yrsmth^0.5))
  Lrecent<-mean(DLM_data@ML[ind])
  Lave<-mean(DLM_data@ML[ind3])
  L0<-0.9*Lave
  Ltarget<-xL*Lave
  if(Lrecent>L0){TAC<-0.5*TACstar*(1+((Lrecent-L0)/(Ltarget-L0)))
  }else{TAC<-0.5*TACstar*(Lrecent/L0)^2                  
  }
  TACfilter(TAC)
}  
class(Ltarget4)<-"DLM_output"


Itarget1<-function(x,DLM_data,reps=100,yrsmth=5,xx=0,Imulti=1.5){
  dependencies="DLM_data@Cat, DLM_data@CV_Cat"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year) # recent 5 years
  ylast<-(DLM_data@LHYear-DLM_data@Year[1])+1 #last historical year
  ind2<-((ylast-(yrsmth-1)):ylast) # historical 5 pre-projection years
  ind3<-((ylast-(yrsmth*2-1)):ylast) # historical 10 pre-projection years
  C_dat<-DLM_data@Cat[x,ind2] 
  TACstar<-(1-xx)*trlnorm(reps,mean(C_dat),DLM_data@CV_Cat/(yrsmth^0.5))
  Irecent<-mean(DLM_data@Ind[x,ind])
  Iave<-mean(DLM_data@Ind[x,ind3])
  Itarget<-Iave*Imulti
  I0<-0.8*Iave
  if(Irecent>I0){TAC<-0.5*TACstar*(1+((Irecent-I0)/(Itarget-I0)))
  }else{TAC<-0.5*TACstar*(Irecent/I0)^2}
  TACfilter(TAC)  
}  
class(Itarget1)<-"DLM_output"

Itarget4<-function(x,DLM_data,reps=100,yrsmth=5,xx=0.3,Imulti=2.5){
  dependencies="DLM_data@Cat, DLM_data@CV_Cat"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year) # recent 5 years
  ylast<-(DLM_data@LHYear-DLM_data@Year[1])+1 #last historical year
  ind2<-((ylast-(yrsmth-1)):ylast) # historical 5 pre-projection years
  ind3<-((ylast-(yrsmth*2-1)):ylast) # historical 10 pre-projection years
  C_dat<-DLM_data@Cat[x,ind2] 
  TACstar<-(1-xx)*trlnorm(reps,mean(C_dat),DLM_data@CV_Cat/(yrsmth^0.5))
  Irecent<-mean(DLM_data@Ind[x,ind])
  Iave<-mean(DLM_data@Ind[x,ind3])
  Itarget<-Iave*Imulti
  I0<-0.8*Iave
  if(Irecent>I0){TAC<-0.5*TACstar*(1+((Irecent-I0)/(Itarget-I0)))
  }else{TAC<-0.5*TACstar*(Irecent/I0)^2}
  TACfilter(TAC)  
}  
class(Itarget4)<-"DLM_output"


Islope1<-function(x,DLM_data,reps=100,yrsmth=5,lambda=0.4,xx=0.2){
  dependencies="DLM_data@Year, DLM_data@Cat, DLM_data@CV_Cat, DLM_data@Ind"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  ylast<-(DLM_data@LHYear-DLM_data@Year[1])+1 #last historical year
  C_dat<-DLM_data@Cat[x,ind]
  if(length(DLM_data@Year)==ylast+1) {TACstar<-(1-xx)*trlnorm(reps,mean(C_dat),DLM_data@CV_Cat/(yrsmth^0.5))
  }else{TACstar<-rep(DLM_data@MPrec[x],reps)}
    
  I_hist<-DLM_data@Ind[x,ind]
  yind<-1:yrsmth
  slppar<-summary(lm(I_hist~yind))$coefficients[2,1:2]
  Islp <-rnorm(reps,slppar[1],slppar[2])
  TAC<-TACstar*(1+lambda*Islp)
  TACfilter(TAC)
}
class(Islope1)<-"DLM_output"

Islope4<-function(x,DLM_data,reps=100,yrsmth=5,lambda=0.2,xx=0.4){
  dependencies="DLM_data@Year, DLM_data@Cat, DLM_data@CV_Cat, DLM_data@Ind"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  ylast<-(DLM_data@LHYear-DLM_data@Year[1])+1 #last historical year
  C_dat<-DLM_data@Cat[x,ind]
  if(length(DLM_data@Year)==ylast+1) {TACstar<-(1-xx)*trlnorm(reps,mean(C_dat),DLM_data@CV_Cat/(yrsmth^0.5))
  }else{TACstar<-rep(DLM_data@MPrec[x],reps)}
  I_hist<-DLM_data@Ind[x,ind]
  yind<-1:yrsmth
  slppar<-summary(lm(I_hist~yind))$coefficients[2,1:2]
  Islp <-rnorm(reps,slppar[1],slppar[2])
  TAC<-TACstar*(1+lambda*Islp)
  TACfilter(TAC)
}
class(Islope4)<-"DLM_output"




IT10<-function(x,DLM_data,reps=100,yrsmth=5,mc=0.1){
  
  dependencies="DLM_data@Ind, DLM_data@Cat, DLMdata@CV_Ind, DLMdata@Iref"
  ind<-max(1,(length(DLM_data@Year)-yrsmth+1)):length(DLM_data@Year)
  
  
  
  deltaI<-mean(DLM_data@Ind[x,ind])/DLM_data@Iref[x]
  if(deltaI<(1-mc))deltaI<-1-mc
  if(deltaI>(1+mc))deltaI<-1+mc
  
  TAC<-DLM_data@MPrec[x]*deltaI*trlnorm(reps,1,DLM_data@CV_Ind[x])
  TAC
}
class(IT10)<-"DLM_output"

IT5<-function(x,DLM_data,reps=100,yrsmth=5,mc=0.05){
  
  dependencies="DLM_data@Ind, DLM_data@Cat, DLMdata@CV_Ind, DLMdata@Iref"
  ind<-max(1,(length(DLM_data@Year)-yrsmth+1)):length(DLM_data@Year)
  deltaI<-mean(DLM_data@Ind[x,ind])/DLM_data@Iref[x]
  if(deltaI<(1-mc))deltaI<-1-mc
  if(deltaI>(1+mc))deltaI<-1+mc
  
  TAC<-DLM_data@MPrec[x]*deltaI*trlnorm(reps,1,DLM_data@CV_Ind[x])
  TAC
}
class(IT5)<-"DLM_output"

ITM<-function(x,DLM_data,reps=100){
  
  dependencies="DLM_data@Ind, DLM_data@Cat, DLMdata@CV_Ind, DLMdata@Iref, DLMdata@Mort"
  mc<-(5+DLM_data@Mort[x]*25)/100
  if(mc>0.2)mc<-0.2
  yrsmth<-floor(4*(1/DLM_data@Mort[x])^(1/4))
  ind<-max(1,(length(DLM_data@Year)-yrsmth+1)):length(DLM_data@Year)
  
  deltaI<-mean(DLM_data@Ind[x,ind])/DLM_data@Iref[x]
  if(deltaI<(1-mc))deltaI<-1-mc
  if(deltaI>(1+mc))deltaI<-1+mc
  
  TAC<-DLM_data@MPrec[x]*deltaI*trlnorm(reps,1,DLM_data@CV_Ind[x])
  TAC
}
class(ITM)<-"DLM_output"

SPmod<-function(x,DLM_data,reps=100,alp=c(0.8,1.2),bet=c(0.8,1.2)){
  dependencies="DLM_data@Cat, DLM_data@Ind, DLM_data@Abun"
  Ir<-length(DLM_data@Ind[x,])
  Cr<-length(DLM_data@Cat[x,])
  rat<-trlnorm(reps,DLM_data@Ind[x,Ir],DLM_data@CV_Ind)/trlnorm(reps,DLM_data@Ind[x,Ir-1],DLM_data@CV_Ind)
  cct<-trlnorm(reps,DLM_data@Cat[x,Cr],DLM_data@CV_Cat)
  Abun<-trlnorm(reps,DLM_data@Abun[x],DLM_data@CV_Abun)
  TAC<-rep(NA,reps)
  TAC[rat<alp[1]]<-cct[rat<alp[1]]*bet[1]
  TAC[rat>alp[1]& rat<alp[2]]<-cct[rat>alp[1]& rat<alp[2]]
 
  cond<-rat>alp[2]
  reps2<-sum(cond)
  if(reps2>0){
    qq1<-trlnorm(reps2,DLM_data@Ind[x,Ir]/Abun,DLM_data@CV_Ind)
    bio1<-DLM_data@Ind[x,Ir-1]/qq1
    bio2<-DLM_data@Ind[x,Ir]/qq1
    cct1<-trlnorm(reps2,DLM_data@Cat[x,Cr-1],DLM_data@CV_Cat)
    PP<-bio2-bio1+cct1
    TAC[cond]<-bet[2]*PP
  }
  TACfilter(TAC)
}
class(SPmod)<-"DLM_output"

SPslope<-function(x,DLM_data,reps=100,yrsmth=4,alp=c(0.9,1.1),bet=c(1.5,0.9)){
  
  dependencies="DLM_data@Year, DLM_data@Cat, DLM_data@Ind, DLM_data@Abun"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  yind<-1:yrsmth
  C_dat<-DLM_data@Cat[x,ind]
  B_dat<-DLM_data@Ind[x,ind]/DLM_data@Ind[x,ind[yrsmth]]*DLM_data@Abun[x]
  Pt_mu<-max(B_dat[yrsmth]-B_dat[yrsmth-1]+C_dat[yrsmth-1],tiny)
  Pt_1<-trlnorm(reps,Pt_mu,DLM_data@CV_Cat)
  It<-exp(predict(lm(log(B_dat)~yind),newdat=list(yind=yrsmth+1)))
  Ilast<-B_dat[yrsmth]
  MC <- max(mean(C_dat), tiny)
  Ct_1<-trlnorm(reps,MC,DLM_data@CV_Cat/(yrsmth^0.5)) # mean catches over the interval
 
  rat<-It/Ilast
  
  mult <- max((1-bet[1]*(Ilast-It)/Ilast), tiny)
  if(rat<alp[1])TAC<-mult*Ct_1
  if(rat>alp[1]&rat<alp[2])TAC<-Ct_1
  if(rat>alp[2])TAC<-bet[2]*Pt_1
  TACfilter(TAC)
}
class(SPslope)<-"DLM_output"

SBT1<-function(x,DLM_data,reps=100,yrsmth=10,k1=1.5,k2=3,gamma=1){
  dependencies="DLM_data@Cat, DLM_data@Year, DLM_data@Ind"
  Cr<-length(DLM_data@Cat[x,])
  cct<-trlnorm(reps,DLM_data@Cat[x,Cr],DLM_data@CV_Cat)
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  I_hist<-DLM_data@Ind[x,ind]
  test<-summary(lm(I_hist~ind))$coefficients[2,1:2]
  lambda<-rnorm(reps,test[1],test[2])
  TAC<-cct*1+k2*lambda
  cond<-lambda<0
  TAC[cond]<-cct[cond]*1-k1*-lambda[cond]^gamma
  TACfilter(TAC) 
}
class(SBT1)<-"DLM_output"

SBT2<-function(x,DLM_data,reps=100,epsB=0.25,epsR=0.75,tauR=5,tauB=7,gamma=1){
  dependencies="DLM_data@Cref, DLM_data@Rec, DLM_data@Cat"
  #Bnow<-trlnorm(reps,DLM_data@Abun[x],DLM_data@CV_Abun)
  #testrat<-Bnow/DLM_data@Bref
  #Ctarg<-rep(NA,reps)
  #Ctarg[testrat>1]<-delta*testrat[testrat>1]^(1-epsB)
  #Ctarg[testrat<1]<-detla*testrat[testrat<1]^(1+epsB)
  Ctarg<-trlnorm(reps,DLM_data@Cref[x],DLM_data@CV_Cref)
  muR<-mean(DLM_data@Rec[x,(length(DLM_data@Rec[x,])-tauR+1):length(DLM_data@Rec[x,])])
  phi<-mean(DLM_data@Rec[x,(length(DLM_data@Rec[x,])-9):length(DLM_data@Rec[x,])])
  Rrat<-muR/phi
  deltaR<-rep(NA,reps)
  deltaR[Rrat>1]<-Rrat[Rrat>1]^(1-epsR)
  deltaR[Rrat<1]<-Rrat[Rrat<1]^(1+epsR)
  TAC<-0.5*(DLM_data@Cat[x,length(DLM_data@Cat[x,])]+Ctarg*deltaR)
  TACfilter(TAC) 
}
class(SBT2)<-"DLM_output"

DD<-function(x,DLM_data,reps=100){
#for(x in 1:nsim){
  dependencies="DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@Mort, DLM_data@CV_Mort. DLM_data@wla, DLM_data@ wlb"
  Linfc<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  t0c <- -trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  t0c[!is.finite(t0c)] <- 0 
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])   # CV of 0.5 as in MacCall 2009
  a<-DLM_data@wla[x]
  b<-DLM_data@wlb[x]

  Winf=DLM_data@wla[x]*DLM_data@vbLinf[x]^DLM_data@wlb[x]
	age<-1:DLM_data@MaxAge
  la<-DLM_data@vbLinf[x]*(1-exp(-DLM_data@vbK[x]*((age-DLM_data@vbt0[x]))))
  wa<-DLM_data@wla[x]*la^DLM_data@wlb[x]
  a50V<-iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x])
  a50V <- max(a50V, 1)
  yind<-(1:length(DLM_data@Cat[x,]))[!is.na(DLM_data@Cat[x,]+DLM_data@Ind[x,])]
  C_hist<-DLM_data@Cat[x,yind]
  E_hist<-C_hist/DLM_data@Ind[x,yind]
  E_hist<-E_hist/mean(E_hist)
  ny_DD<-length(C_hist)
  params<-log(c(DLM_data@Mort[x],mean(C_hist,na.rm=T),DLM_data@Mort[x]))
  k_DD<-ceiling(a50V)   # get age nearest to 50% vulnerability (ascending limb)  -------------
  k_DD[k_DD>DLM_data@MaxAge/2]<-ceiling(DLM_data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho_DD<-(wa[k_DD+2]-Winf)/(wa[k_DD+1]-Winf)
  Alpha_DD<-Winf*(1-Rho_DD)
  So_DD<-exp(-DLM_data@Mort[x]) # get So survival rate
  wa_DD<-wa[k_DD]
  UMSYprior<-c(1-exp(-DLM_data@Mort[x]*0.5),0.3)
  opt<-optim(params,DD_R,opty=1,So_DD=So_DD,Alpha_DD=Alpha_DD,Rho_DD=Rho_DD,
                    ny_DD=ny_DD,k_DD=k_DD,wa_DD=wa_DD,E_hist=E_hist,
                    C_hist=C_hist,UMSYprior=UMSYprior,
                    method="L-BFGS-B",
                    lower=log(exp(params)/20),upper=log(exp(params)*20),
                    hessian=TRUE)

  #Catfit<-DD_R(opt$par,opty=3,So_DD=So_DD,Alpha_DD=Alpha_DD,Rho_DD=Rho_DD,ny_DD=ny_DD,k_DD=k_DD,wa_DD=wa_DD,E_hist=E_hist,C_hist=C_hist,UMSYprior=UMSYprior)
  #plot(Catfit[,1],ylim=c(0,max(Catfit)))
  #lines(Catfit[,2],col="red")

  TAC<-rep(NA,reps)
  #samps<-rmvnorm(reps,opt$par,solve(opt$hessian)) # assuming log parameters are multivariate normal hessian approximation
  samps<-cbind(rnorm(reps,opt$par[1],((opt$par[1])^2)^0.5*0.1),rnorm(reps,opt$par[2],((opt$par[2])^2)^0.5*0.1),rnorm(reps,opt$par[3],((opt$par[3])^2)^0.5*0.1))
  if(reps==1)samps<-matrix(c(opt$par[1],opt$par[2],opt$par[3]),nrow=1)
  for(i in 1:reps)TAC[i]<-DD_R(samps[i,],opty=2,So_DD=So_DD,Alpha_DD=Alpha_DD,Rho_DD=Rho_DD,ny_DD=ny_DD,k_DD=k_DD,wa_DD=wa_DD,E_hist=E_hist,C_hist=C_hist,UMSYprior=UMSYprior)
  TACfilter(TAC)
}
class(DD)<-"DLM_output"


DD4010<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@Mort, DLM_data@CV_Mort. DLM_data@wla, DLM_data@ wlb"
  Linfc<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  t0c<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  t0c[!is.finite(t0c)] <- 0 
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])   # CV of 0.5 as in MacCall 2009
  a<-DLM_data@wla[x]
  b<-DLM_data@wlb[x]
  
  Winf=DLM_data@wla[x]*DLM_data@vbLinf[x]^DLM_data@wlb[x]
  age<-1:DLM_data@MaxAge
  la<-DLM_data@vbLinf[x]*(1-exp(-DLM_data@vbK[x]*((age-DLM_data@vbt0[x]))))
  wa<-DLM_data@wla[x]*la^DLM_data@wlb[x]
  a50V<-iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x])
  a50V <- max(a50V, 1)
  yind<-(1:length(DLM_data@Cat[x,]))[!is.na(DLM_data@Cat[x,]+DLM_data@Ind[x,])]
  C_hist<-DLM_data@Cat[x,yind]
  E_hist<-DLM_data@Ind[x,yind]
  E_hist<-C_hist/E_hist
  E_hist<-E_hist/mean(E_hist)
  ny_DD<-length(C_hist)
  params<-log(c(DLM_data@Mort[x],mean(C_hist,na.rm=T),DLM_data@Mort[x]))
  k_DD<-ceiling(a50V)   # get age nearest to 50% vulnerability (ascending limb)  -------------
  k_DD[k_DD>DLM_data@MaxAge/2]<-ceiling(DLM_data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho_DD<-(wa[k_DD+2]-Winf)/(wa[k_DD+1]-Winf)
  Alpha_DD<-Winf*(1-Rho_DD)
  So_DD<-exp(-DLM_data@Mort[x]) # get So survival rate
  wa_DD<-wa[k_DD]
  UMSYprior<-c(1-exp(-DLM_data@Mort*0.5),0.3)
  opt<-optim(params,DD_R,opty=1,So_DD=So_DD,Alpha_DD=Alpha_DD,Rho_DD=Rho_DD,
             ny_DD=ny_DD,k_DD=k_DD,wa_DD=wa_DD,E_hist=E_hist,
             C_hist=C_hist,UMSYprior=UMSYprior,
             method="L-BFGS-B",
             lower=log(exp(params)/20),upper=log(exp(params)*20),
             hessian=TRUE)
  
  #Catfit<-DD_R(opt$par,opty=3,So_DD=So_DD,Alpha_DD=Alpha_DD,Rho_DD=Rho_DD,ny_DD=ny_DD,k_DD=k_DD,wa_DD=wa_DD,E_hist=E_hist,C_hist=C_hist,UMSYprior=UMSYprior)
  #plot(Catfit[,1],ylim=c(0,max(Catfit)))
  #lines(Catfit[,2],col="red")
  
  TAC<-rep(NA,reps)
  dep<-rep(NA,reps)
  #samps<-rmvnorm(reps,opt$par,solve(opt$hessian)) # assuming log parameters are multivariate normal hessian approximation
  samps<-cbind(rnorm(reps,opt$par[1],((opt$par[1])^2)^0.5*0.1),rnorm(reps,opt$par[2],((opt$par[2])^2)^0.5*0.1),rnorm(reps,opt$par[3],((opt$par[3])^2)^0.5*0.1))
  if(reps==1)samps<-matrix(c(opt$par[1],opt$par[2],opt$par[3]),nrow=1)
  for(i in 1:reps)TAC[i]<-DD_R(samps[i,],opty=2,So_DD=So_DD,Alpha_DD=Alpha_DD,Rho_DD=Rho_DD,ny_DD=ny_DD,k_DD=k_DD,wa_DD=wa_DD,E_hist=E_hist,C_hist=C_hist,UMSYprior=UMSYprior)
  for(i in 1:reps)dep[i]<-DD_R(samps[i,],opty=3,So_DD=So_DD,Alpha_DD=Alpha_DD,Rho_DD=Rho_DD,ny_DD=ny_DD,k_DD=k_DD,wa_DD=wa_DD,E_hist=E_hist,C_hist=C_hist,UMSYprior=UMSYprior)
  cond1<-!is.na(dep) & dep<0.4 & dep>0.1
  cond2<-!is.na(dep) & dep<0.1
  TAC[cond1]<-TAC[cond1]*(dep[cond1]-0.1)/0.3
  TAC[cond2]<-TAC[cond2]*tiny # this has to still be stochastic albeit very small
  TACfilter(TAC)
}
class(DD4010)<-"DLM_output"

DD_R<-function(params,opty,So_DD,Alpha_DD,Rho_DD,ny_DD,k_DD,wa_DD,E_hist,C_hist,UMSYprior){
  UMSY_DD=exp(params[1])
  MSY_DD=exp(params[2])
  q_DD=exp(params[3])
  SS_DD=So_DD*(1-UMSY_DD)    # Initialise for UMSY, MSY and q leading.
  Spr_DD=(SS_DD*Alpha_DD/(1-SS_DD)+wa_DD)/(1-Rho_DD*SS_DD)
  DsprDu_DD=-So_DD*(Rho_DD/(1-Rho_DD*SS_DD)*Spr_DD+1/(1-Rho_DD*SS_DD)*(Alpha_DD/(1-SS_DD)+SS_DD*Alpha_DD/(1-SS_DD)^2))
  Arec_DD=1/(((1-UMSY_DD)^2)*(Spr_DD+UMSY_DD*DsprDu_DD))
  Brec_DD=UMSY_DD*(Arec_DD*Spr_DD-1/(1-UMSY_DD))/MSY_DD
  Spr0_DD=(So_DD*Alpha_DD/(1-So_DD)+wa_DD)/(1-Rho_DD*So_DD)
  Ro_DD=(Arec_DD*Spr0_DD-1)/(Brec_DD*Spr0_DD)
  Bo_DD=Ro_DD*Spr0_DD
  No_DD=Ro_DD/(1-So_DD)

  B_DD<-rep(NA,ny_DD+1)
  N_DD<-rep(NA,ny_DD+1)
  R_DD<-rep(NA,ny_DD+k_DD)
  Cpred_DD<-rep(NA,ny_DD)

  B_DD[1]=Bo_DD
  N_DD[1]=No_DD
  R_DD[1:k_DD]=Ro_DD

  for(tt in 1:ny_DD){

    Surv_DD=So_DD*exp(-q_DD*E_hist[tt])
    Cpred_DD[tt]=B_DD[tt]*(1-exp(-q_DD*E_hist[tt]))
    Sp_DD=B_DD[tt]-Cpred_DD[tt]
    R_DD[tt+k_DD]=Arec_DD*Sp_DD/(1+Brec_DD*Sp_DD);
    B_DD[tt+1]=Surv_DD*(Alpha_DD*N_DD[tt]+Rho_DD*B_DD[tt])+wa_DD*R_DD[tt+1]
    N_DD[tt+1]=Surv_DD*N_DD[tt]+R_DD[tt+1]

  }
  Cpred_DD[Cpred_DD<tiny]<-tiny

  if(opty==1){
    test<-dnorm(log(Cpred_DD),log(C_hist),0.25,log=T)
    test2<-dlnorm(UMSY_DD,log(UMSYprior[1]),UMSYprior[2],log=T)
    test[is.na(test)]<--1000
    test[test==(-Inf)]<--1000
    if(is.na(test2)|test2==-Inf|test2==Inf)test2<-1000
    return(-sum(test,test2))      # return objective function
  }else if(opty==2){                                  # return MLE TAC estimate
    UMSY_DD*B_DD[ny_DD]
  }else if(opty==3){ 
    B_DD[tt+1]/Bo_DD
  }else{
    cbind(C_hist,Cpred_DD)                           # return observations vs predictions
  }
}

DBSRA<-function(x,DLM_data,reps=100){  # returns a vector of DBSRA estimates of the TAC for a particular simulation x
  # for(x in 1:nsim){
  dependencies="DLM_data@Cat, DLM_data@Dep, DLM_data@CV_Dep, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@FMSY_M, DLM_data@CV_FMSY_M,DLM_data@BMSY_B0, DLM_data@CV_BMSY_B0, DLM_data@L50"
  C_hist<-DLM_data@Cat[x,]
  TAC<-rep(NA,reps)
  DBSRAcount<-1
  if (is.na(DLM_data@Dep[x]) | is.na(DLM_data@CV_Dep[x])) return(NA)
  while(DBSRAcount<(reps+1)){
    depo<-max(0.01,min(0.99,DLM_data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
    Bt_K<-rbeta(100,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
    Bt_K<-Bt_K[Bt_K>0.00999&Bt_K<0.99001][1] # interval censor (0.01,0.99)  as in Dick and MacCall 2011
    Mdb<-trlnorm(100,DLM_data@Mort[x],DLM_data@CV_Mort[x])
    Mdb<-Mdb[Mdb<0.9][1]    # !!!! maximum M is 0.9   interval censor
    if(is.na(Mdb))Mdb<-0.9  # !!!! maximum M is 0.9   absolute limit
    FMSY_M<-trlnorm(1,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])
    BMSY_K<-rbeta(100,alphaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x]*DLM_data@BMSY_B0[x]),betaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x]*DLM_data@BMSY_B0[x])) #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
    tryBMSY_K <- BMSY_K[BMSY_K > 0.05 & BMSY_K < 0.95][1] # interval censor (0.05,0.95) as in Dick and MacCall, 2011
	if (is.na(tryBMSY_K)) {
	  Min <- min(BMSY_K, na.rm=TRUE)
	  Max <- max(BMSY_K, na.rm=TRUE)
	  if (Max <= 0.05) BMSY_K <- 0.05
	  if (Min >= 0.95) BMSY_K <- 0.95
	}
	if (!is.na(tryBMSY_K)) BMSY_K <- tryBMSY_K
	
    adelay<-max(floor(iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x])),1)
    opt<-optimize(DBSRAopt,log(c(0.01*mean(C_hist),1000*mean(C_hist))),C_hist=C_hist,nys=length(C_hist),Mdb=Mdb,
              FMSY_M=FMSY_M,BMSY_K=BMSY_K,Bt_K=Bt_K,adelay=adelay,tol=0.01)
    # if(opt$objective<0.1){
    Kc<-exp(opt$minimum)
    BMSYc<-Kc*BMSY_K
    FMSYc<-Mdb*FMSY_M
    UMSYc<-(FMSYc/(FMSYc+Mdb))*(1-exp(-(FMSYc+Mdb)))
    MSYc<-Kc*BMSY_K*UMSYc
    TAC[DBSRAcount]<-UMSYc*Kc*Bt_K
    DBSRAcount<-DBSRAcount+1
    # }
	 
	# print(DBSRAcount)
  } # end of reps
  TACfilter(TAC)
   # print(x)
  # }
}  # end of DBSRA_apply
class(DBSRA)<-"DLM_output"

DBSRA_40<-function(x,DLM_data,reps=100){  # returns a vector of DBSRA estimates of the TAC for a particular simulation x
  dependencies="DLM_data@Cat, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@FMSY_M, DLM_data@CV_FMSY_M, DLM_data@BMSY_B0, DLM_data@CV_BMSY_B0, DLM_data@L50"
  C_hist<-DLM_data@Cat[x,]
  TAC<-rep(NA,reps)
  DBSRAcount<-1
  if (is.na(DLM_data@Dep[x]) | is.na(DLM_data@CV_Dep[x])) return(NA)
  while(DBSRAcount<(reps+1)){
    depo<-0.4
    Bt_K<-rbeta(100,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
    Bt_K<-Bt_K[Bt_K>0.00999&Bt_K<0.99001][1] # interval censor (0.01,0.99)  as in Dick and MacCall 2011
    Mdb<-rlnorm(100,mconv(DLM_data@Mort[x],DLM_data@CV_Mort[x]*DLM_data@Mort[x]),sdconv(DLM_data@Mort[x],DLM_data@CV_Mort[x]*DLM_data@Mort[x]))   # log space stdev 0.4 as in Dick and MacCall 2011
    Mdb<-Mdb[Mdb<0.9][1]    # !!!! maximum M is 0.9   interval censor
    if(is.na(Mdb))Mdb<-0.9  # !!!! maximum M is 0.9   absolute limit
    FMSY_M<-trlnorm(1,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])
    BMSY_K<-rbeta(100,alphaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x]*DLM_data@BMSY_B0[x]),betaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x]*DLM_data@BMSY_B0[x])) #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
    BMSY_K<-BMSY_K[BMSY_K>0.05&BMSY_K<0.95][1] # interval censor (0.05,0.95) as in Dick and MacCall, 2011
    adelay<-max(floor(iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x])),1)
    opt<-optimize(DBSRAopt,log(c(0.1*mean(C_hist),1000*mean(C_hist))),C_hist=C_hist,nys=length(C_hist),Mdb=Mdb,
              FMSY_M=FMSY_M,BMSY_K=BMSY_K,Bt_K=Bt_K,adelay=adelay,tol=0.01)
    #if(opt$objective<0.1){
      Kc<-exp(opt$minimum)
      BMSYc<-Kc*BMSY_K
      FMSYc<-Mdb*FMSY_M
      UMSYc<-(FMSYc/(FMSYc+Mdb))*(1-exp(-(FMSYc+Mdb)))
      MSYc<-Kc*BMSY_K*UMSYc
      TAC[DBSRAcount]<-UMSYc*Kc*Bt_K
      DBSRAcount<-DBSRAcount+1
    #}
  } # end of reps
  TACfilter(TAC)
}  # end of DBSRA_apply
class(DBSRA_40)<-"DLM_output"

DBSRA_ML<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@Cat, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@FMSY_M, DLM_data@CV_FMSY_M, DLM_data@BMSY_B0, DLM_data@CV_BMSY_B0, DLM_data@L50, DLM_data@CAL, DLM_data@Year, DLM_data@Cat"
  C_hist<-DLM_data@Cat[x,]
  TAC<-rep(NA,reps)
  DBSRAcount<-1
  if (is.na(DLM_data@Dep[x]) | is.na(DLM_data@CV_Dep[x])) return(NA)
  while(DBSRAcount<(reps+1)){
    Linfc<-trlnorm(1,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
    Kc<-trlnorm(1,DLM_data@vbK[x],DLM_data@CV_vbK[x])
    Mdb<-trlnorm(100,DLM_data@Mort[x],DLM_data@CV_Mort[x])
    Mdb<-Mdb[Mdb<0.9][1]    # !!!! maximum M is 0.9   interval censor
    if(is.na(Mdb))Mdb<-0.9  # !!!! maximum M is 0.9   absolute limit
    Z<-MLne(x,DLM_data,Linfc=Linfc,Kc=Kc,ML_reps=1,MLtype="dep")
    FM<-Z-Mdb
    FM[FM<0]<-0.01
    nyears<-length(DLM_data@Year)
    Ct1<-mean(DLM_data@Cat[x,1:3])
    Ct2<-mean(DLM_data@Cat[x,(nyears-2):nyears])
    dep<-c(Ct1,Ct2)/(1-exp(-FM[,c(1,2)]))
    Bt_K<-dep[2]/dep[1]
    if(Bt_K<0.01)Bt_K<-0.01       # interval censor / temporary hack to avoid doing multiple depletion estimates that would take far too long
    if(Bt_K>0.99)Bt_K<-0.99       # interval censor / temporary hack to avoid doing multiple depletion estimates that would take far too long

    FMSY_M<-trlnorm(1,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])
    BMSY_K<-rbeta(100,alphaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x]*DLM_data@BMSY_B0[x]),betaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x]*DLM_data@BMSY_B0[x])) #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
    BMSY_K<-BMSY_K[BMSY_K>0.05&BMSY_K<0.95][1] # interval censor (0.05,0.95) as in Dick and MacCall, 2011
    adelay<-max(floor(iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x])),1)
    opt<-optimize(DBSRAopt,log(c(0.1*mean(C_hist),1000*mean(C_hist))),C_hist=C_hist,nys=length(C_hist),Mdb=Mdb,
              FMSY_M=FMSY_M,BMSY_K=BMSY_K,Bt_K=Bt_K,adelay=adelay,tol=0.01)
    # if(opt$objective<0.1){
      Kc<-exp(opt$minimum)
      BMSYc<-Kc*BMSY_K
      FMSYc<-Mdb*FMSY_M
      UMSYc<-(FMSYc/(FMSYc+Mdb))*(1-exp(-(FMSYc+Mdb)))
      MSYc<-Kc*BMSY_K*UMSYc
      TAC[DBSRAcount]<-UMSYc*Kc*Bt_K
      DBSRAcount<-DBSRAcount+1
    # }
  } # end of reps
  TACfilter(TAC)
}
class(DBSRA_ML)<-"DLM_output"

DBSRA4010<-function(x,DLM_data,reps=100){  # returns a vector of DBSRA estimates of the TAC for a particular simulation x
  #for(x in 1:nsim){
  dependencies="DLM_data@Cat, DLM_data@Dep, DLM_data@CV_Dep, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@FMSY_M, DLM_data@CV_FMSY_M,DLM_data@BMSY_B0, DLM_data@CV_BMSY_B0, DLM_data@L50"
  C_hist<-DLM_data@Cat[x,]
  TAC<-rep(NA,reps)
  DBSRAcount<-1
  if (is.na(DLM_data@Dep[x]) | is.na(DLM_data@CV_Dep[x])) return(NA)
  while(DBSRAcount<(reps+1)){
    depo<-max(0.01,min(0.99,DLM_data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
    Bt_K<-rbeta(100,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
    Bt_K<-Bt_K[Bt_K>0.00999&Bt_K<0.99001][1] # interval censor (0.01,0.99)  as in Dick and MacCall 2011
    Mdb<-trlnorm(100,DLM_data@Mort[x],DLM_data@CV_Mort[x])
    Mdb<-Mdb[Mdb<0.9][1]    # !!!! maximum M is 0.9   interval censor
    if(is.na(Mdb))Mdb<-0.9  # !!!! maximum M is 0.9   absolute limit
    FMSY_M<-trlnorm(1,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])
    BMSY_K<-rbeta(100,alphaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x]*DLM_data@BMSY_B0[x]),betaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x]*DLM_data@BMSY_B0[x])) #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
    BMSY_K<-BMSY_K[BMSY_K>0.05&BMSY_K<0.95][1] # interval censor (0.05,0.95) as in Dick and MacCall, 2011
    adelay<-max(floor(iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x])),1)
    opt<-optimize(DBSRAopt,log(c(0.01*mean(C_hist),1000*mean(C_hist))),C_hist=C_hist,nys=length(C_hist),Mdb=Mdb,
                  FMSY_M=FMSY_M,BMSY_K=BMSY_K,Bt_K=Bt_K,adelay=adelay,tol=0.01)
    #if(opt$objective<0.1){
    Kc<-exp(opt$minimum)
    BMSYc<-Kc*BMSY_K
    FMSYc<-Mdb*FMSY_M
    UMSYc<-(FMSYc/(FMSYc+Mdb))*(1-exp(-(FMSYc+Mdb)))
    MSYc<-Kc*BMSY_K*UMSYc
    TAC[DBSRAcount]<-UMSYc*Kc*Bt_K
    # 40-10 rule
    if(Bt_K<0.4 & Bt_K>0.1)TAC[DBSRAcount]<-TAC[DBSRAcount]*(Bt_K-0.1)/0.3
    if(Bt_K<0.1)TAC[DBSRAcount]<-TAC[DBSRAcount]*tiny # this has to still be stochastic albeit very small
    DBSRAcount<-DBSRAcount+1
    #}
  } # end of reps
  TACfilter(TAC)
  #}
}  # end of DBSRA_apply
class(DBSRA4010)<-"DLM_output"

DBSRAopt<-function(lnK,C_hist,nys,Mdb,FMSY_M,BMSY_K,Bt_K,adelay){         # the optimization for B0 given DBSRA assumptions
  Kc<-exp(lnK)
  n<-getn(BMSY_K)
  g<-gety(n)
  FMSY<-FMSY_M*Mdb
  UMSY<-(FMSY/(FMSY+Mdb))*(1-exp(-(FMSY+Mdb)))
  MSY<-Kc*BMSY_K*UMSY
  # Bjoin rules from Dick & MacCall 2011  ---------
  Bjoin_K<-0.5
  if(BMSY_K<0.3)Bjoin_K<-0.5*BMSY_K
  if(BMSY_K>0.3&BMSY_K<0.5)Bjoin_K<-0.75*BMSY_K-0.075
  Bjoin<-Bjoin_K*Kc
  PBjoin<-prodPTF(Bjoin_K,n,MSY)
  cp<-(1-n)*g*MSY*(Bjoin^(n-2))*Kc^-n
  Bc<-rep(NA,nys)
  Bc[1]<-Kc
  obj<-0
  for(yr in 2:nys){
    yref<-max(1,yr-adelay)
    if(Bc[yref]>Bjoin|BMSY_K>0.5){
      Bc[yr]<-Bc[yr-1]+g*MSY*(Bc[yref]/Kc)-g*MSY*(Bc[yref]/Kc)^n-C_hist[yr-1]
    }else{
      Bc[yr]<-Bc[yr-1]+Bc[yref]*((PBjoin/Bjoin)+cp*(Bc[yref]-Bjoin))-C_hist[yr-1]
    }
    if(Bc[yr]<0)obj<-obj+log(-Bc[yr])
    Bc[yr]<-max(0.000001,Bc[yr])
  }
  obj+((Bc[nys]/Kc)-Bt_K)^2
}  # end of DBSRA optimization function

# DCAC =====================================================================================================================
# Variables for parallel processing of DCAC_apply function using sfSapply()
C_tot<-nyearsDCAC<-NULL

DCAC<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@AvC, DLM_data@t, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@FMSY_M, DLM_data@CV_FMSY_M, DLM_data@Dt, DLM_data@CV_Dt, DLM_data@BMSY_B0, DLM_data@CV_BMSY_B0"
  C_tot<-DLM_data@AvC[x]*DLM_data@t[x]
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])   # CV of 0.5 as in MacCall 2009
  FMSY_M<-trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x]) # standard deviation of 0.2 - referred to as 'standard error' in MacCall 2009
  Bt_K<-trlnorm(reps,DLM_data@Dt[x],DLM_data@CV_Dt[x])
  BMSY_K<-rbeta(reps,alphaconv(DLM_data@BMSY_B0[x],DLM_data@BMSY_B0[x]*DLM_data@CV_BMSY_B0[x]),betaconv(DLM_data@BMSY_B0[x],DLM_data@BMSY_B0[x]*DLM_data@CV_BMSY_B0[x])) #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
  TACfilter(C_tot/(DLM_data@t[x]+((1-Bt_K)/(BMSY_K*FMSY_M*Mdb))))
} # end of DCAC
class(DCAC)<-"DLM_output"

DCAC4010<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@AvC, DLM_data@t, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@FMSY_M, DLM_data@CV_FMSY_M, DLM_data@Dt, DLM_data@CV_Dt, DLM_data@BMSY_B0, DLM_data@CV_BMSY_B0"
  C_tot<-DLM_data@AvC[x]*DLM_data@t[x]
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])   # CV of 0.5 as in MacCall 2009
  FMSY_M<-trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x]) # standard deviation of 0.2 - referred to as 'standard error' in MacCall 2009
  Bt_K<-trlnorm(reps,DLM_data@Dt[x],DLM_data@CV_Dt[x])
  BMSY_K<-rbeta(reps,alphaconv(DLM_data@BMSY_B0[x],DLM_data@BMSY_B0[x]*DLM_data@CV_BMSY_B0[x]),betaconv(DLM_data@BMSY_B0[x],DLM_data@BMSY_B0[x]*DLM_data@CV_BMSY_B0[x])) #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
  TAC<-C_tot/(DLM_data@t[x]+((1-Bt_K)/(BMSY_K*FMSY_M*Mdb)))
  # 40-10 rule
  cond1<-Bt_K<0.4 & Bt_K>0.1
  cond2<-Bt_K<0.1
  if (length(cond1)>0) TAC[cond1]<-TAC[cond1]*(Bt_K[cond1]-0.1)/0.3
  if (length(cond2)>0) TAC[cond2]<-TAC[cond2]*tiny # this has to still be stochastic albeit very small
  if (length(cond1) <1 & length(cond2) < 1) return(NA)
  TACfilter(TAC)
  
} # end of DCAC
class(DCAC4010)<-"DLM_output"

DCAC_40<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@AvC, DLM_data@t, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@FMSY_M, DLM_data@CV_FMSY_M, DLM_data@BMSY_B0, DLM_data@CV_BMSY_B0"
  C_tot<-DLM_data@AvC[x]*DLM_data@t[x]
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  FMSY_M<-trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])
  Bt_K<-0.4
  BMSY_K<-rbeta(reps,alphaconv(DLM_data@BMSY_B0[x],DLM_data@BMSY_B0[x]*DLM_data@CV_BMSY_B0[x]),betaconv(DLM_data@BMSY_B0[x],DLM_data@BMSY_B0[x]*DLM_data@CV_BMSY_B0[x])) #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
  TACfilter(C_tot/(DLM_data@t[x]+((1-Bt_K)/(BMSY_K*FMSY_M*Mdb))))
} # end of DCAC40
class(DCAC_40)<-"DLM_output"

DCAC_ML<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@AvC, DLM_data@t, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@FMSY_M, DLM_data@CV_FMSY_M, DLM_data@BMSY_B0, DLM_data@CV_BMSY_B0, DLM_data@Year, DLM_data@CAL, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbK, DLM_data@CV_vbK"
  if (is.na(DLM_data@BMSY_B0[x]) | is.na(DLM_data@CV_BMSY_B0[x])) return(NA)
  if (is.na(DLM_data@FMSY_M[x]) | is.na(DLM_data@CV_FMSY_M[x])) return(NA)
  C_tot<-DLM_data@AvC[x]*DLM_data@t[x]
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])   # default CV of 0.5 as in MacCall 2009
  FMSY_M<-trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x]) # standard deviation of 0.2 - referred to as 'standard error' in MacCall 2009
  Linfc<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Z<-MLne(x,DLM_data,Linfc=Linfc,Kc=Kc,ML_reps=reps,MLtype="dep")
  FM<-Z-Mdb
  FM[FM<0]<-0.01
  nyears<-length(DLM_data@Year)
  Ct1<-mean(DLM_data@Cat[x,1:3])
  Ct2<-mean(DLM_data@Cat[x,(nyears-2):nyears])
  dep<-rep(c(Ct1,Ct2),each=reps)/(1-exp(-FM[,c(1,2)]))
  if (reps == 1) Bt_K<-dep[2]/dep[1]
  if (reps > 1)  Bt_K<-dep[,2]/dep[,1]
  BMSY_K<-rbeta(reps,alphaconv(DLM_data@BMSY_B0[x],DLM_data@BMSY_B0[x]*DLM_data@CV_BMSY_B0[x]),betaconv(DLM_data@BMSY_B0[x],DLM_data@BMSY_B0[x]*DLM_data@CV_BMSY_B0[x])) #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
  TAC<-C_tot/(DLM_data@t[x]+((1-Bt_K)/(BMSY_K*FMSY_M*Mdb)))
  TACfilter(TAC)
} # end of DCAC_ML
class(DCAC_ML)<-"DLM_output"

BK<-function(x,DLM_data,reps=100){   # Beddington and Kirkwood life-history analysis ==============================================
  dependencies="DLM_data@LFC, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@Abun, DLM_data@CV_Abun, DLM_data@vbK, DLM_data@CV_vbK"
  Lc<-trlnorm(reps*10,DLM_data@LFC[x],0.2)
  Linfc<-trlnorm(reps*10,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Ac<-trlnorm(reps*10,DLM_data@Abun[x],DLM_data@CV_Abun[x])
  Kc<-trlnorm(reps*10,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  TAC<-Ac*(0.6*Kc)/(0.67-(Lc/Linfc))         # robustifying for use in MSE
  TACfilter(TAC[TAC>0][1:reps])              # Interval censor only those positive catch recommendations

}  # end of BK
class(BK)<-"DLM_output"


BK_CC<-function(x,DLM_data,reps=100,Fmin=0.005){
  dependencies="DLM_data@LFC, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@CAA, DLM_data@Mort"
  Lc<-trlnorm(reps,DLM_data@LFC[x],0.2)
  Linfc<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Mdb<-trlnorm(reps*10,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  MuC<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  Cc<-trlnorm(reps,MuC,DLM_data@CV_Cat[x])
  Zdb<-CC(x,DLM_data,reps=reps*10)
  Fdb<-Zdb-Mdb
  ind<-(1:(reps*10))[Fdb>Fmin][1:reps]
  Fdb<-Fdb[ind]
  Mdb<-Mdb[ind]
  SM <- sum(is.na(ind))
  if (SM > 0 ) {
    Mdb[is.na(ind)] <- trlnorm(SM,DLM_data@Mort[x],DLM_data@CV_Mort[x])
    Fdb[is.na(ind)] <- Fmin
  }	
  
  Ac<-Cc/(1-exp(-Fdb))
  TAC<-Ac*(0.6*Kc)/(0.67-(Lc/Linfc))  
  TACfilter(TAC) 
    
}  # end of BK_CC
class(BK_CC)<-"DLM_output"

BK_ML<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@LFC, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@CAL, DLM_data@Mort"
  Lc<-trlnorm(reps*10,DLM_data@LFC[x],0.2)
  Linfc<-trlnorm(reps*10,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Kc<-trlnorm(reps*10,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Mdb<-trlnorm(reps*10,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  Z<-MLne(x,DLM_data,Linfc=Linfc,Kc=Kc,ML_reps=reps*2,MLtype="F")
  FM<-Z-Mdb
  MuC<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  Cc<-trlnorm(reps,MuC,DLM_data@CV_Cat[x])
  Ac<-Cc/(1-exp(-FM))
  FMSY<-(0.6*Kc)/(0.67-(Lc/Linfc))  # robustifying for use in MSETAC<-Ac*FMSY
  TAC<-Ac*FMSY
  print(TAC[TAC>0&TAC<(mean(TAC,na.rm=T)+3*sd(TAC,na.rm=T))][1:reps])
}
class(BK_ML)<-"DLM_output"

Fratio<-function(x,DLM_data,reps=100){  # FMSY / M ratio method e.g. Gulland ===============================================================================
  depends="DLM_data@Abun,DLM_data@CV_Abun,DLM_data@FMSY_M, DLM_data@CV_FMSY_M,DLM_data@Mort,DLM_data@CV_Mort"
  Ac<-trlnorm(reps,DLM_data@Abun[x],DLM_data@CV_Abun[x])
  TACfilter(Ac*trlnorm(reps,DLM_data@Mort[x],
            DLM_data@CV_Mort[x])*trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x]))
} # end of Fratio
class(Fratio)<-"DLM_output"

Fratio4010<-function(x,DLM_data,reps=100){  # FMSY / M ratio method e.g. Gulland ===============================================================================
  dependencies="DLM_data@Abun, DLM_data@CV_Abun, DLM_data@FMSY_M, DLM_data@CV_FMSY_M, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@Dep"
  Ac<-trlnorm(reps,DLM_data@Abun[x],DLM_data@CV_Abun[x])
  TAC<-Ac*trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])*trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])
  Bt_K<-trlnorm(reps,DLM_data@Dt[x],DLM_data@CV_Dt[x])
  # 40-10 rule
  cond1<-Bt_K<0.4 & Bt_K>0.1
  cond2<-Bt_K<0.1
  TAC[cond1]<-TAC[cond1]*(Bt_K[cond1]-0.1)/0.3
  TAC[cond2]<-TAC[cond2]*tiny # this has to still be stochastic albeit very small
  TACfilter(TAC)  
} # end of Fratio
class(Fratio4010)<-"DLM_output"

Fratio_CC<-function(x,DLM_data,reps=100,Fmin=0.005){ # FMSY / M ratio method using catch curve analysis to determine current abundance ==================================
# for (x in 1:nsim) {
  dependencies=" DLM_data@FMSY_M, DLM_data@CV_FMSY_M, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@Cat, DLM_data@CV_Cat, DLM_data@CAA"
  MuC<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  Cc<-trlnorm(reps,MuC,DLM_data@CV_Cat[x])
  Mdb<-trlnorm(reps*10,DLM_data@Mort[x],DLM_data@CV_Mort[x])   # CV of 0.5 as in MacCall 2009
  Zdb<-CC(x,DLM_data,reps=reps*10)
  Fdb<-Zdb-Mdb
  ind <-(1:(reps*10))[Fdb>0.005][1:reps]
  
  Fdb<-Fdb[ind]
  Mdb<-Mdb[ind]
  SM <- sum(is.na(ind))
  if (SM > 0 ) {
    Mdb[is.na(ind)] <- trlnorm(SM,DLM_data@Mort[x],DLM_data@CV_Mort[x])
    Fdb[is.na(ind)] <- Fmin
  }	
  
  Ac<-Cc/(1-exp(-Fdb))
  TAC<-Ac*Mdb*trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])

  TACfilter(TAC)
# }
} # end of Fratio_CC
class(Fratio_CC)<-"DLM_output"


Fratio_ML<-function(x,DLM_data,reps=100){
  dependencies=" DLM_data@FMSY_M, DLM_data@CV_FMSY_M, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@Cat, DLM_data@CV_Cat, DLM_data@CAL"
  MuC<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  Cc<-trlnorm(reps,MuC,DLM_data@CV_Cat[x])
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])   # CV of 0.5 as in MacCall 2009
  Linfc<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Z<-MLne(x,DLM_data,Linfc=Linfc,Kc=Kc,ML_reps=reps,MLtype="F")
  FM<-Z-Mdb
  Ac<-Cc/(1-exp(-FM))
  TAC<-Ac*trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])*trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  TACfilter(TAC)
}
class(Fratio_ML)<-"DLM_output"

SPMSY<-function(x,DLM_data,reps=100){  # Martell and Froese 2012 Schaefer SP estimate of MSY given priors on r, k and depletion
  #for(x in 1:100){
  dependencies="DLM_data@MaxAge, DLM_data@vbK, DLM_data@L50, DLM_data@Cat"
  maxage<-DLM_data@MaxAge
  nsamp<-reps*200

  # Froese 2012 http://www.fishbase.de/rfroese/Catch-MSY_SummaryFinal.doc
  rule<-rep(4,3)

  if(DLM_data@vbK[x]>0.3){   # K rules
   rule[1]<-1
  }else if(DLM_data@vbK[x]<0.3 & DLM_data@vbK[x]>0.16){
   rule[1]<-2
  }else if(DLM_data@vbK[x]<0.16 & DLM_data@vbK[x]>0.05){
   rule[1]<-3
  }
  AM<-iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x])
  if(AM<1.5){  # Age at maturity rules
   rule[2]<-1
  }else if(AM<4.5 & AM>1.5){
   rule[2]<-2
  }else if(AM<10 & AM>4.5){
   rule[2]<-3
  }

  if(DLM_data@MaxAge<4){   # Maximum age rules
   rule[3]<-1
  }else if(DLM_data@MaxAge<11 & DLM_data@MaxAge>3){
   rule[3]<-2
  }else if(DLM_data@MaxAge<31 & DLM_data@MaxAge>10){
   rule[3]<-3
  }

  if(mean(rule)<1.5) rsamp<-runif(nsamp,0.6,1.5)
  if(mean(rule)>1.5&mean(rule)<2.5)rsamp<-runif(nsamp,0.2,1)
  if(mean(rule)>2.5&mean(rule)<3.5)rsamp<-runif(nsamp,0.05,0.5)
  if(mean(rule)>3.5) rsamp<-runif(nsamp,0.015,0.1)

  Ksamp<-runif(nsamp,mean(DLM_data@Cat[x,])/rsamp,(10*mean(DLM_data@Cat[x,]))/rsamp)
  nyears<-length(DLM_data@Cat[x,])
  B<-array(NA,dim=c(nsamp,nyears))

  if(DLM_data@Cat[x,1]<(0.5*max(DLM_data@Cat[x,]))){    # Martell and Froese decision rules (makes absolutely no sense to me!)
    B[,1]<-Ksamp*runif(nsamp,0.5,0.9)
  }else{
    B[,1]<-Ksamp*runif(nsamp,0.3,0.6)
  }

  if(DLM_data@Cat[x,nyears]<(0.5*max(DLM_data@Cat[x,]))){    # Martell and Froese decision rules (makes absolutely no sense to me!)
    LB<-0.01
    UB<-0.4
  }else{
    LB<-0.3
    UB<-0.7
  }

  for(i in 2:nyears){
   B[,i]<-B[,i-1]-DLM_data@Cat[x,i-1]
   B[,i]<-B[,i]+rsamp*B[,i]*(1-B[,i]/Ksamp)
  }
  B<-B/rep(Ksamp,nyears)
  cond<-(B[,nyears]>=LB)&(B[,nyears]<=UB)
  if (sum(cond) < 1) {
	B[B[,nyears]>=UB ,nyears]<- UB
	cond<-(B[,nyears]>=LB)&(B[,nyears]<=UB)
  }	
  dep<-B[cond,nyears][1:reps]
  MSY<-rsamp[cond][1:reps]*Ksamp[cond][1:reps]/4
  Kc<-Ksamp[cond][1:reps]
  rc<-rsamp[cond][1:reps]
  TAC<-Kc*dep*rc/2

  if(sum(!is.na(TAC))<ceiling(reps/10)){ # a fudge of the original method that widens current depletion to the lowest and higest bounds to get an TAC sample
    cond<-(B[,nyears]>=0.01)&(B[,nyears]<=0.7)
    dep<-B[cond,nyears][1:reps]
    MSY<-rsamp[cond][1:reps]*Ksamp[cond][1:reps]/4
    Kc<-Ksamp[cond][1:reps]
    rc<-rsamp[cond][1:reps]
    TAC<-Kc*dep*rc/2
	
  }
  #}
  
  TACfilter(TAC)
} # end of SPMSY
class(SPMSY)<-"DLM_output"

MCD<-function(x,DLM_data,reps=100){  # Daft method to demonstrate the relative value of information of current depletion
  dependencies="DLM_data@Dep, DLM_data@CV_Dep, DLM_data@Cat"
  depo<-max(0.01,min(0.99,DLM_data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
  Bt_K<-rbeta(reps*100,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
  Bt_K<-Bt_K[Bt_K>0.00999&Bt_K<0.99001][1:reps] # interval censor (0.01,0.99)  as in Dick and MacCall 2011
  AvC<-rlnorm(reps,log(mean(DLM_data@Cat[x,],na.rm=T)),DLM_data@CV_Cat[x])
  TAC<-AvC*2*Bt_K
  TACfilter(TAC)
}
class(MCD)<-"DLM_output"

MCD4010<-function(x,DLM_data,reps=100){  # Daft method to demonstrate the relative value of information of current depletion
  dependencies="DLM_data@Dep, DLM_data@CV_Dep, DLM_data@Cat"
  depo<-max(0.01,min(0.99,DLM_data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
  Bt_K<-rbeta(reps*100,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
  Bt_K<-Bt_K[Bt_K>0.00999&Bt_K<0.99001][1:reps] # interval censor (0.01,0.99)  as in Dick and MacCall 2011
  AvC<-rlnorm(reps,log(mean(DLM_data@Cat[x,],na.rm=T)),DLM_data@CV_Cat[x])
  TAC<-AvC*2*Bt_K
  
  #40-10 HCR
  cond1<-Bt_K<0.4 & Bt_K>0.1
  cond2<-Bt_K<0.1
  TAC[cond1]<-TAC[cond1]*(Bt_K[cond1]-0.1)/0.3
  TAC[cond2]<-TAC[cond2]*tiny # this has to still be stochastic albeit very small
  
  TACfilter(TAC)
}
class(MCD4010)<-"DLM_output"

SPSRA<-function(x,DLM_data,reps=100){  # Surplus productin stock reduction analysis T.Carruthers - basically an SP version of DBSRA
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@Dep, DLM_data@CV_Dep, DLM_data@Cat, DLM_data@steep"
  Mvec<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  Kvec<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Linfvec <- trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  t0vec<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  if(all(is.nan(t0vec))) t0vec <- rep(0,reps) 
  hvec<-trlnorm(reps,DLM_data@steep[x],DLM_data@CV_steep[x])
  if (all(!is.finite(hvec))) return(NA)
  rsamp<-getr(x,DLM_data,Mvec,Kvec,Linfvec,t0vec,hvec,maxage=DLM_data@MaxAge,r_reps=reps)
  dep<-trlnorm(reps,DLM_data@Dep[x],DLM_data@CV_Dep[x])
  Ct<-DLM_data@Cat[x,]
  Csamp<-array(rep(Ct,each=reps)*trlnorm(length(Ct)*reps,1,DLM_data@CV_Cat[x]),dim=c(reps,length(Ct)))
  Psamp<-array(trlnorm(length(Ct)*reps,1,0.1),dim=c(reps,length(Ct)))
  Ksamp<-rep(NA,reps)
  for(i in 1:reps)Ksamp[i]<-exp(optimize(SPSRAopt,log(c(mean(Csamp[i,]),1000*mean(Csamp[i,]))),dep=dep[i],r=rsamp[i],Ct=Csamp[i,],PE=Psamp[i,])$minimum)
  MSY<-Ksamp*rsamp/4
  TAC<-Ksamp*dep*rsamp/2
  TACfilter(TAC)
}
class(SPSRA)<-"DLM_output"

SPSRA_ML<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@CAL, DLM_data@Cat, DLM_data@steep"
  Mvec<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  Kvec<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Linfvec=trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  t0vec<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  hvec<-trlnorm(reps,DLM_data@steep[x],DLM_data@CV_steep[x])
  rsamp<-getr(x,DLM_data,Mvec,Kvec,Linfvec,t0vec,hvec,maxage=DLM_data@MaxAge,r_reps=reps)
  Z<-MLne(x,DLM_data,Linfc=Linfvec,Kc=Kvec,ML_reps=reps,MLtype="dep")
  FM<-Z-Mvec
  FM[FM<0]<-0.01
  nyears<-length(DLM_data@Year)
  Ct1<-mean(DLM_data@Cat[x,1:3])
  Ct2<-mean(DLM_data@Cat[x,(nyears-2):nyears])
  dep<-rep(c(Ct1,Ct2),each=reps)/(1-exp(-FM[,c(1,2)]))
  if (reps == 1) dep<-dep[2]/dep[1]
  if (reps > 1)  dep<-dep[,2]/dep[,1]
  Ksamp<-rep(NA,reps)
  Ct<-DLM_data@Cat[x,]
  Csamp<-array(rep(Ct,each=reps)*trlnorm(length(Ct)*reps,1,DLM_data@CV_Cat[x]),dim=c(reps,length(Ct)))
  Psamp<-array(trlnorm(length(Ct)*reps,1,0.1),dim=c(reps,length(Ct)))
  for(i in 1:reps)Ksamp[i]<-exp(optimize(SPSRAopt,log(c(mean(Csamp[i,]),1000*mean(Csamp[i,]))),dep=dep[i],r=rsamp[i],Ct=Csamp[i,],PE=Psamp[i,])$minimum)
  MSY<-Ksamp*rsamp/4
  TAC<-Ksamp*dep*rsamp/2
  TACfilter(TAC)
}
class(SPSRA_ML)<-"DLM_output"

SPSRAopt<-function(lnK,dep,r,Ct,PE){
  nyears<-length(Ct)
  B<-rep(NA,nyears)
  B[1]<-exp(lnK)
  OBJ<-0
  for(y in 2:nyears){
    if((B[y-1]-Ct[y-1])<0)OBJ<-OBJ+(B[y-1]-Ct[y-1])^2
    B[y]<-max(0.01,B[y-1]-Ct[y-1])
    B[y]<-B[y]+r*B[y]*(1-B[y]/B[1])*PE[y]
  }
  return(OBJ+((B[nyears]/B[1])-dep)^2)
}


YPR<-function(x,DLM_data,reps=100){   # Yield per recruit analysis F01 - Meaghan Bryan
  #for(x in 1:10){
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@MaxAge, DLM_data@Abun, DLM_data@CV_Abun, DLM_data@wla, DLM_data@wlb"
  Linfc<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  t0c<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  t0c[!is.finite(t0c)] <- 0 
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  LFS<-trlnorm(reps,DLM_data@LFS[x],DLM_data@CV_LFS[x])
  a<-DLM_data@wla[x]
  b<-DLM_data@wlb[x]
  Ac<-trlnorm(reps,DLM_data@Abun[x],DLM_data@CV_Abun[x])
  FMSY<-YPRopt(Linfc,Kc,t0c,Mdb,a,b,LFS,DLM_data@MaxAge,reps)
  TAC<-Ac*FMSY
  TACfilter(TAC)
  #}
}   # end of YPR
class(YPR)<-"DLM_output"

YPR_CC<-function(x,DLM_data,reps=100,Fmin=0.005){
  #for(x in 1:16){
    dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@MaxAge, DLM_data@wla, DLM_data@wlb, DLM_data@CAA, DLM_data@Cat"
  Linfc<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  t0c<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  t0c[!is.finite(t0c)] <- 0 
  LFS<-trlnorm(reps,DLM_data@LFS[x],DLM_data@CV_LFS[x])
  a<-DLM_data@wla[x]
  b<-DLM_data@wlb[x]
  MuC<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  Cc<-trlnorm(reps,MuC,DLM_data@CV_Cat[x])
  
  Mdb<-trlnorm(reps*10,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  Zdb<-CC(x,DLM_data,reps=reps*10)
  Fdb<-Zdb-Mdb
  ind<-(1:(reps*10))[Fdb>Fmin][1:reps]
  
  Fdb<-Fdb[ind]
  Mdb<-Mdb[ind]
  SM <- sum(is.na(ind))
  if (SM > 0 ) {
    Mdb[is.na(ind)] <- trlnorm(SM,DLM_data@Mort[x],DLM_data@CV_Mort[x])
    Fdb[is.na(ind)] <- Fmin
  }	
  
  Ac<-Cc/(1-exp(-Fdb))
  FMSY<-YPRopt(Linfc,Kc,t0c,Mdb,a,b,LFS,DLM_data@MaxAge,reps)
  TAC<-Ac*FMSY
 # }
  TACfilter(TAC)
}
class(YPR_CC)<-"DLM_output"

YPR_ML<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@MaxAge, DLM_data@wla, DLM_data@wlb, DLM_data@CAL, DLM_data@Cat"
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  Linfc<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  t0c<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  t0c[!is.finite(t0c)] <- 0 
  LFS<-trlnorm(reps,DLM_data@LFS[x],DLM_data@CV_LFS[x])
  a<-DLM_data@wla[x]
  b<-DLM_data@wlb[x]
  MuC<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  Cc<-trlnorm(reps,MuC,DLM_data@CV_Cat[x])
  Z<-MLne(x,DLM_data,Linfc=Linfc,Kc=Kc,ML_reps=reps,MLtype="F")
  FM<-Z-Mdb
  Ac<-Cc/(1-exp(-FM))
  FMSY<-YPRopt(Linfc,Kc,t0c,Mdb,a,b,LFS,DLM_data@MaxAge,reps)
  TAC<-Ac*FMSY
  TACfilter(TAC)
}
class(YPR_ML)<-"DLM_output"

Fdem<-function(x,DLM_data,reps=100){   # Demographic FMSY estimate (FMSY=r/2)
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@MaxAge, DLM_data@wla, DLM_data@wlb, DLM_data@Abun, DLM_data@CV_Abun, DLM_data@steep, DLM_data@CV_steep"
  Mvec<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Linfc=trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  t0c<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  t0c[!is.finite(t0c)] <- 0 
  hvec<-trlnorm(reps,DLM_data@steep[x],DLM_data@CV_steep[x])
  Ac<-trlnorm(reps,DLM_data@Abun[x],DLM_data@CV_Abun[x])
  FMSY<-getr(x,DLM_data,Mvec,Kc,Linfc,t0c,hvec,maxage=DLM_data@MaxAge,r_reps=reps)/2
  TAC<-FMSY*Ac
  TACfilter(TAC)
}
class(Fdem)<-"DLM_output"

Fdem_CC<-function(x,DLM_data,reps=100,Fmin=0.005){
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@MaxAge, DLM_data@wla, DLM_data@wlb, DLM_data@CAA, DLM_data@steep, DLM_data@CV_steep"
  Mvec<-trlnorm(reps*10,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Linfc=trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  t0c<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  t0c[!is.finite(t0c)] <- 0 
  hvec<-trlnorm(reps,DLM_data@steep[x],DLM_data@CV_steep[x])
  MuC<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  Cc<-trlnorm(reps,MuC,DLM_data@CV_Cat[x])
  Zdb<-CC(x,DLM_data,reps=reps*10)
  Fdb<-Zdb-Mvec
  ind<-(1:(reps*10))[Fdb>Fmin][1:reps]
  
  Fdb<-Fdb[ind]
  SM <- sum(is.na(ind))
  if (SM > 0 ) {
    Fdb[is.na(ind)] <- Fmin
  }	

  Ac<-Cc/(1-exp(-Fdb))
  FMSY<-getr(x,DLM_data,Mvec,Kc,Linfc,t0c,hvec,maxage=DLM_data@MaxAge,r_reps=reps)/2
  TAC<-FMSY*Ac
  
  TACfilter(TAC)
}
class(Fdem_CC)<-"DLM_output"

# Catch curve estimate of recent F (naive) ========================================================================================
CC<-function(x,DLM_data,reps=100){
  ny<-dim(DLM_data@CAA)[2]
  CAA<-apply(DLM_data@CAA[x,max(ny-2,1):ny,],2,sum) # takes last two years as the sample (or last year if there is only one
  maxageobs<-length(CAA)
  AFS<-which.max(CAA)
  AFS[AFS>(maxageobs-3)]<-maxageobs-3   # provides at least three datapoints
  
  nS<-ceiling(sum(CAA)/2)
  y <-log(CAA[AFS:maxageobs]/sum(CAA[AFS:maxageobs],na.rm=T))
  xc<-1:length(y)
  y[y=='-Inf']<-NA
  mod <- lm(y~xc)
  chk <- sum(is.na(coef(mod))) # check if model failed
  if (chk) {
    return(NA)
  } else {
    coefs <-summary(mod,weights=CAA[AFS:maxageobs])$coefficients[2,1:2]
	coefs[is.nan(coefs)] <- tiny
   return(-rnorm(reps,coefs[1],coefs[2]))
  } 
}
# class(CC)<-"DLM_output"

Fdem_ML<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@MaxAge, DLM_data@wla, DLM_data@wlb, DLM_data@CAL, DLM_data@steep, DLM_data@CV_steep"
  Mvec<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Linfc=trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  t0c<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  t0c[!is.finite(t0c)] <- 0 
  hvec<-trlnorm(reps,DLM_data@steep[x],DLM_data@CV_steep[x])
  MuC<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  Cc<-trlnorm(reps,MuC,DLM_data@CV_Cat[x])
  Z<-MLne(x,DLM_data,Linfc=Linfc,Kc=Kc,ML_reps=reps,MLtype="F")
  FM<-Z-Mvec
  Ac<-Cc/(1-exp(-FM))
  FMSY<-getr(x,DLM_data,Mvec,Kc,Linfc,t0c,hvec,maxage=DLM_data@MaxAge,r_reps=reps)/2
  TAC<-FMSY*Ac
  TACfilter(TAC)
}
class(Fdem_ML)<-"DLM_output"

CompSRA<-function(x,DLM_data,reps=100){    # optimize for fixed F to get you to current depletion C/Fcur = abundance
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@MaxAge, DLM_data@wla, DLM_data@CV_wla, DLM_data@wlb, DLM_data@CV_wlb, DLM_data@L50, DLM_data@CV_L50, DLM_data@CAA, DLM_data@steep, DLM_data@CV_steep, DLM_data@LFS, DLM_data@CV_LFS, DLM_data@LFC, DLM_data@CV_LFC, DLM_data@Cat"
  maxage<-DLM_data@MaxAge
  TAC<-rep(NA,reps)
  for(i in 1:reps){
    Mc<-trlnorm(1,DLM_data@Mort[x],DLM_data@CV_Mort)
    hc<-trlnorm(1,DLM_data@steep[x],DLM_data@CV_steep[x])
    Linfc<-trlnorm(1,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
    Kc<-trlnorm(1,DLM_data@vbK[x],DLM_data@CV_vbK[x])
    t0c<--trlnorm(1,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
	t0c[!is.finite(t0c)] <- 0 
    LFSc<-trlnorm(1,DLM_data@LFS[x],DLM_data@CV_LFS[x])
    LFCc<-trlnorm(1,DLM_data@LFC[x],DLM_data@CV_LFC[x])
    AMc<-trlnorm(1,iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x]),DLM_data@CV_L50[x])
    ac<-trlnorm(1,DLM_data@wla[x],DLM_data@CV_wla[x])
    bc<-trlnorm(1,DLM_data@wlb[x],DLM_data@CV_wlb[x])
    Catch<-DLM_data@Cat[x,]
    ny<-length(Catch)
    nyCAA<-dim(DLM_data@CAA)[2]
    CAA<-DLM_data@CAA[x,max(nyCAA-2,1):nyCAA,] # takes last two years as the sample (or last year if there is only one

    Nac<-exp(-Mc*((1:maxage)-1)) # put a rough range on estimate of R0 assuming a mean harvest rate of 10%
    Lac<-Linfc*(1-exp(-Kc*((1:maxage)-t0c)))
    Wac<-ac*Lac^bc
    AFC<-log(1-min(0.99,LFCc/Linfc))/-Kc+t0c
    AFS<-log(1-min(0.99,LFSc/Linfc))/-Kc+t0c
    if (AFC >= 0.7 * maxage) AFC <- 0.7 * maxage
    if (AFS >= 0.9 * maxage) AFS <- 0.9 * maxage
    KES<-max(2,ceiling(mean(c(AFC,AFS))))
    pred<-Nac*Wac
    pred[1:(KES-1)]<-0
    pred<-pred/sum(pred)
    pred<-((mean(Catch)/0.1)*pred/Wac)/exp(-(1:maxage)*Mc)
    pred<-pred[pred>0]
    R0range<-c(mean(pred)/1000,mean(pred)*1000)

    fit<-optimize(SRAfunc,log(R0range),Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc,Catch,CAA)
    Ac<-SRAfunc(fit$minimum,Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc,Catch,CAA,opt=2)
    fit2<-optimize(SRAFMSY,log(c(0.0001,3)),Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc)
    # FMSY<-SRAFMSY(fit2$minimum,Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc,opt=F)
	FMSY <- exp(fit2$minimum)
    if ((FMSY / Mc) > 3) FMSY <- 3 * Mc
    TAC[i]<-Ac*FMSY
	# message(i, " of ", reps)
    # flush.console()
  }
  TACfilter(TAC)
}
class(CompSRA)<-"DLM_output"

CompSRA4010<-function(x,DLM_data,reps=100){    # optimize for fixed F to get you to current depletion C/Fcur = abundance
 # for (x in 1:nsim) {
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@MaxAge, DLM_data@wla, DLM_data@CV_wla, DLM_data@wlb, DLM_data@CV_wlb, DLM_data@L50, DLM_data@CV_L50, DLM_data@CAA, DLM_data@steep, DLM_data@CV_steep, DLM_data@LFS, DLM_data@CV_LFS, DLM_data@LFC, DLM_data@CV_LFC, DLM_data@Cat"
  maxage<-DLM_data@MaxAge
  TAC<-Bt_K<-rep(NA,reps)
  for(i in 1:reps){
    Mc<-trlnorm(1,DLM_data@Mort[x],DLM_data@CV_Mort)
    hc<-trlnorm(1,DLM_data@steep[x],DLM_data@CV_steep[x])
    Linfc<-trlnorm(1,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
    Kc<-trlnorm(1,DLM_data@vbK[x],DLM_data@CV_vbK[x])
    t0c<--trlnorm(1,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
	t0c[!is.finite(t0c)] <- 0 
    LFSc<-trlnorm(1,DLM_data@LFS[x],DLM_data@CV_LFS[x])
    LFCc<-trlnorm(1,DLM_data@LFC[x],DLM_data@CV_LFC[x])
    AMc<-trlnorm(1,iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x]),DLM_data@CV_L50[x])
    ac<-trlnorm(1,DLM_data@wla[x],DLM_data@CV_wla[x])
    bc<-trlnorm(1,DLM_data@wlb[x],DLM_data@CV_wlb[x])
    Catch<-DLM_data@Cat[x,]
    ny<-length(Catch)
    nyCAA<-dim(DLM_data@CAA)[2]
    CAA<-DLM_data@CAA[x,max(nyCAA-2,1):nyCAA,] # takes last two years as the sample (or last year if there is only one
    
    Nac<-exp(-Mc*((1:maxage)-1)) # put a rough range on estimate of R0 assuming a mean harvest rate of 10%
    Lac<-Linfc*(1-exp(-Kc*((1:maxage)-t0c)))
    Wac<-ac*Lac^bc
	
	AFC<-log(1-min(0.99,LFCc/Linfc))/-Kc+t0c
	AFS<-log(1-min(0.99,LFSc/Linfc))/-Kc+t0c
	if (AFC >= 0.7 * maxage) AFC <- 0.7 * maxage
	if (AFS >= 0.9 * maxage) AFS <- 0.9 * maxage
	
    KES<-max(2,ceiling(mean(c(AFC,AFS))))
    pred<-Nac*Wac
    pred[1:(KES-1)]<-0
    pred<-pred/sum(pred)
    pred<-((mean(Catch)/0.1)*pred/Wac)/exp(-(1:maxage)*Mc)
    pred<-pred[pred>0]
    R0range<-c(mean(pred)/1000,mean(pred)*1000)
    
    fit<-optimize(SRAfunc,log(R0range),Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc,Catch,CAA)
    Ac<-SRAfunc(fit$minimum,Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc,Catch,CAA,opt=2)
    Bt_K[i]<-SRAfunc(fit$minimum,Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc,Catch,CAA,opt=3)
    fit2<-optimize(SRAFMSY,log(c(0.0001,3)),Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc)
    FMSY<-SRAFMSY(fit2$minimum,Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc,opt=F)
    if ((FMSY / Mc) > 3) FMSY <- 3 * Mc
    TAC[i]<-Ac*FMSY
  }   
  # 40-10 rule
  cond1<-Bt_K<0.4 & Bt_K>0.1
  cond2<-Bt_K<0.1
  TAC[cond1]<-TAC[cond1]*(Bt_K[cond1]-0.1)/0.3
  TAC[cond2]<-TAC[cond2]*tiny # this has to still be stochastic albeit very small

  TACfilter(TAC)
  
  # message(x, " of ", nsim)
  # flush.console()
  # }
}
class(CompSRA4010)<-"DLM_output"

# options(warn=2)
# options(warn=1)

SRAfunc<-function(lnR0c,Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc,Catch,CAA,opt=1){

  ny<-length(Catch)
  AFC<-log(1-min(0.99,LFCc/Linfc))/-Kc+t0c
  AFS<-log(1-min(0.99,LFSc/Linfc))/-Kc+t0c
  if (AFC >= 0.7 * maxage) AFC <- 0.7 * maxage
  if (AFS >= 0.9 * maxage) AFS <- 0.9 * maxage
  KES<-max(2,ceiling(mean(c(AFC,AFS))))
  vul<-rep(1,maxage)
  vul[1:(KES-1)]<-0
  Mac<-rep(1,maxage)
  Mac[1:max(1,floor(AMc))]<-0
  Lac<-Linfc*(1-exp(-Kc*((1:maxage)-t0c)))
  Wac<-ac*Lac^bc
  R0c<-exp(lnR0c)
  N<-exp(-Mc*((1:maxage)-1))*R0c
  SSN<-Mac*N                                 # Calculate initial spawning stock numbers
  Biomass<-N*Wac
  SSB<-SSN*Wac                               # Calculate spawning stock biomass

  B0<-sum(Biomass)
  SSB0<-sum(SSB)
  SSN0<-SSN
  SSBpR<-sum(SSB)/R0c                              # Calculate spawning stock biomass per recruit
  SSNpR<-SSN/R0c

  CN<-array(NA,dim=c(ny,maxage))
  HR<-rep(0,maxage)
  pen<-0
  for(y in 1:ny){
    # set up some indices for indexed calculation
    VB<-Biomass[KES:maxage]*exp(-Mc)
    CB<-Catch[y]*VB/sum(VB)
    testHR<-CB[1]/VB[1]
    if(testHR>0.8)pen<-pen+(testHR-0.8)^2
    HR[KES:maxage]<-min(testHR,0.8)
    FMc<--log(1-HR)                                        # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
    Zc<-FMc+Mc

    CN[y,]<-N * (1-exp(-Zc))*(FMc/Zc)
    N[2:maxage]<-N[1:(maxage-1)]*exp(-Zc[1:(maxage-1)])         # Total mortality
    N[1]<-(0.8*R0c*hc*sum(SSB))/(0.2*SSBpR*R0c*(1-hc)+(hc-0.2)*sum(SSB))  # Recruitment assuming regional R0 and stock wide steepness
    #print(N[1])
    Biomass<-N*Wac
    SSN<-N*Mac
    SSB<-SSN*Wac

  } # end of year
  
  CN[CN<0] <- 0 # stop any negative catches
  syear<-ny-dim(CAA)[1]+1
  pred<-CN[syear:ny,]
  pred<-pred/array(apply(pred,1,sum),dim=c(dim(CAA)[1],maxage))
 
  fobj<-pen-sum(log(pred+tiny)*CAA,na.rm=T)
  if(opt==1){return(fobj)
  }else if(opt==2){return(sum(Biomass))
  }else if(opt==3){sum(SSB)/sum(SSB0)
  }
  #CBc<-sum(CB)
}

SRAFMSY<-function(lnFMc,Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc,opt=T){

  FMc<-exp(lnFMc)
  ny<-100
  AFC<-log(1-min(0.99,LFCc/Linfc))/-Kc+t0c
  AFS<-log(1-min(0.99,LFSc/Linfc))/-Kc+t0c
  if (AFC >= 0.7 * maxage) AFC <- 0.7 * maxage
  if (AFS >= 0.9 * maxage) AFS <- 0.9 * maxage
  
  KES<-max(2,ceiling(mean(c(AFC,AFS))))
  vul<-rep(1,maxage)
  vul[1:(KES-1)]<-0
  Mac<-rep(1,maxage)
  Mac[1:max(1,floor(AMc))]<-0
  Lac<-Linfc*(1-exp(-Kc*((1:maxage)-t0c)))
  Wac<-ac*Lac^bc
  R0c<-1
  N<-exp(-Mc*((1:maxage)-1))*R0c
  SSN<-Mac*N   # Calculate initial spawning stock numbers
  Biomass<-N*Wac
  SSB<-SSN*Wac                               # Calculate spawning stock biomass

  B0<-sum(Biomass)
  SSB0<-sum(SSB)
  SSN0<-SSN
  SSBpR<-sum(SSB)/R0c                              # Calculate spawning stock biomass per recruit
  SSNpR<-SSN/R0c

  N<-N/2
  SSN<-Mac*N   # Calculate initial spawning stock numbers
  Biomass<-N*Wac
  SSB<-SSN*Wac

  for(y in 1:ny){
    # set up some indices for indexed calculation
                                           # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
    Zc<-FMc*vul+Mc
    CN<-N*(1-exp(-Zc))*(FMc/Zc)
    CB<-CN*Wac
    Biomass<-N*Wac
    N[2:maxage]<-N[1:(maxage-1)]*exp(-Zc[1:(maxage-1)])         # Total mortality
    N[1]<-(0.8*R0c*hc*sum(SSB))/(0.2*SSBpR*R0c*(1-hc)+(hc-0.2)*sum(SSB))  # Recruitment assuming regional R0 and stock wide steepness
    #print(N[1])
    SSN<-N*Mac
    SSB<-SSN*Wac

  } # end of year
  
  if(opt){return(-sum(CB))
  }else{return(FMc)
  }
}


# Yield per recruit estimate of FMSY Meaghan Bryan 2013 ==================================
YPRopt=function(Linfc,Kc,t0c,Mdb,a,b,LFS,maxage,reps=100) {

  nf<-200
  frates<-seq(0,3,length.out=nf)

  Winf=a*Linfc^b
	rat<-LFS/Linfc
	rat[rat>0.8]<-0.8     # need to robustify this for occasionally very high samples of LFS
  tc=log(1-rat)/-Kc+t0c
	tc=round(tc,0)
  tc[tc<1]<-1
  tc[tc>maxage]<-maxage

	vul<-array(0,dim=c(reps,maxage))
	mat<-array(0,dim=c(reps,maxage))
	lx<-array(NA,dim=c(reps,maxage))
	lxo<-array(NA,dim=c(reps,maxage))

  ypr<-array(NA,dim=c(reps,nf))
	sbpr<-array(NA,dim=c(reps,nf))
	sbpr.ratio<-array(NA,dim=c(reps,nf))
	sbpr.dif<-array(NA,dim=c(reps,nf))

  f.max<-array(NA,dim=c(reps,maxage))

	#average weight at age - follow von Bertalanffy growth
	age<-array(rep(1:maxage,each=reps),dim=c(reps,maxage))
  la<-Linfc*(1-exp(-Kc*((age-t0c))))
  wa<-a*la^b

	#vulnerability schedule - assumes knife-edge vulnerability, where all individuals age tc to maxage are fully vulnerbale
	#all individulas less than age tc are not vulnerable
	for(i in 1:reps){
    if(tc[i]>0)vul[i,tc[i]:maxage]<-1
  	if(tc[i]>1)mat[i,max(1,tc[i]-1):maxage]<-1
  }

  lx[,1]<-1
  lxo[,1]<-1
  for(k in 1:nf){
    for(i in 2:maxage){
      lx[,i]=lx[,i-1]*exp(-(Mdb+vul[,i-1]*frates[k]))
      lxo[,i]=lx[,i]*exp(-Mdb)
    }
    phi_vb=apply(lx*wa*vul,1,sum)
		sbpro=apply(lxo*wa*mat,1,sum)

		ypr[,k]=(1-exp(-frates[k]))*phi_vb
		sbpr[,k]=apply(lx*wa*mat,1,sum)
		sbpr.ratio[,k]=sbpr[,k]/sbpro
		sbpr.dif[,k]=abs(sbpr.ratio[,k]-0.3)	#hard code comparison ratio
  }

  #frates[apply(ypr,1,which.max)]    Fmaxypr

  # More code that derived F0.1 in 'per recruit analysis.R' (Meaghan Bryan)
  slope.origin=(ypr[,2]-ypr[,1])/(frates[2]-frates[1])
	slope.10=round(0.1*slope.origin,2)

	slope=array(NA,dim=dim(ypr))#vector(length=length(ypr))
	slope[,1]=slope.origin
	for(i in 3:ncol(ypr))
	{
		slope[,i-1]=round((ypr[,i]-ypr[,i-1])/(frates[i]-frates[i-1]),2)
	}
	dif=abs(slope-slope.10)
  dif[is.na(dif)]<-10e10
	frates[apply(dif,1,which.min)]#frates[which.min(dif)]
}


MLne<-function(x,DLM_data,Linfc,Kc,ML_reps=100,MLtype="F"){
  year<-1:dim(DLM_data@CAL)[2]
  nlbin<-ncol(DLM_data@CAL[x,,])
  nlyr<-nrow(DLM_data@CAL[x,,])
  mlbin<-(DLM_data@CAL_bins[1:nlbin]+DLM_data@CAL_bins[2:(nlbin+1)])/2
  nbreaks<-1
  Z<-matrix(NA,nrow=ML_reps,ncol=nbreaks+1)
  Z2<-rep(NA,ML_reps)
  temp<-apply(DLM_data@CAL[x,,],2,sum)
  Lc<-mlbin[which.max(temp)] # modal length  
  
  for(i in 1:ML_reps){
    mlen<-rep(NA,length(year))
    ss<-ceiling(apply(DLM_data@CAL[x,,],1,sum)/2)
    if(MLtype=="dep"){
      for(y in 1:length(year)) {
	    if (sum(DLM_data@CAL[x,y,] > 0) > 0.25 * length(DLM_data@CAL[x,y,])) {
	      temp2<-sample(mlbin,ceiling(sum(DLM_data@CAL[x,y,])/2),replace=T,prob=DLM_data@CAL[x,y,])
              mlen[y]<-mean(temp2[temp2>=Lc],na.rm=TRUE)
	    }  
      }
      Z[i,]<-bhnoneq(year=year,mlen=mlen,ss=ss,K=Kc[i],Linf=Linfc[i],Lc=Lc,nbreaks=nbreaks,
           styrs=ceiling(length(year)*((1:nbreaks)/(nbreaks+1))),stZ=rep(0.5,nbreaks+1))
    }else{
      
      #ind<-(which.min(((DLM_data@CAL_bins-DLM_data@LFS[x])^2)^0.5)-1):(length(DLM_data@CAL_bins)-1)
      for(y in 1:length(year)) {
	    if (sum(DLM_data@CAL[x,y,] > 0) > 0.25 * length(DLM_data@CAL[x,y,])) {
	      temp2<-sample(mlbin,ceiling(sum(DLM_data@CAL[x,y,])/2),replace=T,prob=DLM_data@CAL[x,y,])
              mlen[y]<-mean(temp2[temp2>=Lc],na.rm=TRUE)
	    }  
      }		
      mlen<-mean(mlen[(length(mlen)-2):length(mlen)], na.rm=TRUE)
      Z2<-bheq(K=Kc[i],Linf=Linfc[i],Lc=Lc,Lbar=mlen)
    }
  }
  if(MLtype=="F")return(Z2)
  if(MLtype=="dep")return(Z)
}

bheq<-function(K,Linf,Lc,Lbar){
  K*(Linf-Lbar)/(Lbar-Lc)
}

bhnoneq<-function(year,mlen,ss,K,Linf,Lc,nbreaks,styrs,stZ) {
  mlen[mlen<=0|is.na(mlen)]<--99
  ss[ss<=0|is.na(ss)|mlen==-99]<-0
  stpar<-c(stZ,styrs)
  # results <- optim(stpar,bhnoneq_LL,method="BFGS",year=year,Lbar=mlen,ss=ss,
                   # nbreaks=nbreaks,K=K,Linf=Linf,Lc=Lc,control=list(maxit=1e6))
  results <- optim(stpar,bhnoneq_LL,method="Nelder-Mead", year=year,Lbar=mlen,ss=ss,
                   nbreaks=nbreaks,K=K,Linf=Linf,Lc=Lc,control=list(maxit=1e6),
				   hessian=FALSE)					   
  return(results$par[1:(nbreaks+1)])
}

getdep<-function(lnFF,targ,Md,Linfd,Kd,t0d,AFSd,ad,bd,maxage,opt){

   FF<-exp(lnFF)
   Z<-rep(Md,maxage)
   Z[1]<-0
   Z[AFSd:maxage]<-Z[AFSd:maxage]+FF
   for(a in 2:maxage)Z[a]<-Z[a-1]+Z[a]
   Nd<-exp(-Z)
   Nobs<-Nd
   Nobs[1:max(1,(AFSd-1))]<-0
   Ld<-Linfd*(1-exp(-Kd*(((1:maxage)-0.5)-t0d)))

   if(opt)return(((sum(Nobs*Ld)/sum(Nobs))-targ)^2)
   if(!opt){
     Nd0<-exp(-Md*((1:maxage)-0.5))
     Wd<-ad*Ld^bd
     return((sum(Nd*Wd)/sum(Nd))/(sum(Nd0*Wd)/sum(Nd0)))
   }
}

getr <- function(x,DLM_data,Mvec,Kvec,Linfvec,t0vec,hvec,maxage,r_reps=100){
  r<-rep(NA,r_reps)
  for(i in 1:r_reps){
    log.r=log(0.3)

    opt=optimize(demofn,lower=log(0.0001),upper=log(1.4),
    M=Mvec[i],
    amat=iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x]),
    sigma=0.2,
    K=Kvec[i],
    Linf=Kvec[i],
    to=t0vec[i],
    hR=hvec[i],
    maxage=maxage,
    a=DLM_data@wla[x],
    b=DLM_data@wlb[x])
    #demographic2(opt$minimum,M[x],ageM[x],0.2,K[x],Linf,t0,steepness[x],maxage,a,b)$r
    r[i]<-exp(opt$minimum)
  }
  r
}

iVB<-function(t0,K,Linf,L)((-log(1-L/Linf))/K+t0) # Inverse Von-B

EDCAC<-function (x, DLM_data, reps = 100) # extended depletion-corrected average catch (Harford and Carruthers 2015)
{
  dependencies = "DLM_data@AvC, DLM_data@t, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@FMSY_M, DLM_data@CV_FMSY_M, DLM_data@Dt, DLM_data@CV_Dt, DLM_data@BMSY_B0, DLM_data@CV_BMSY_B0"
  C_tot <- DLM_data@AvC[x] * DLM_data@t[x]
  Mdb <- trlnorm(reps, DLM_data@Mort[x], DLM_data@CV_Mort[x])
  FMSY_M <- trlnorm(reps, DLM_data@FMSY_M[x], DLM_data@CV_FMSY_M[x])
  Bt_K <- trlnorm(reps, DLM_data@Dt[x], DLM_data@CV_Dt[x])
  BMSY_K <- rbeta(reps, alphaconv(DLM_data@BMSY_B0[x], DLM_data@BMSY_B0[x]*DLM_data@CV_BMSY_B0[x]), 
                  betaconv(DLM_data@BMSY_B0[x], DLM_data@BMSY_B0[x]*DLM_data@CV_BMSY_B0[x]))
  dcac<-C_tot/(DLM_data@t[x] + ((1 - Bt_K)/(BMSY_K * FMSY_M * Mdb)))
  TAC<-dcac*Bt_K/BMSY_K
  TACfilter(TAC)
}
class(EDCAC)<-"DLM_output"

AvC<-function(x,DLM_data,reps=100)rlnorm(reps,log(mean(DLM_data@Cat[x,],na.rm=T)),0.2)
class(AvC)<-"DLM_output"


LBSPR_ItTAC <- function(x, DLM_data, yrsmth=1,reps=reps) {
 dependencies="DLM_data@CAL, DLM_data@CAL_bins, DLM_data@vbLinf, DLM_data@vbK, DLM_data@Mort, LM_data@vbK, DLM_data@L50, DLM_data@L95, DLM_data@wlb" 
  
  MiscList <- LBSPR(x, DLM_data, yrsmth=yrsmth,reps=reps)
  if(all(is.na(MiscList[[1]]))) return(rep(NA, 6))
  if(all(is.na(MiscList[[1]][,2]))) return(rep(NA, 6))
  XX <- 1:4 
  YY <- MiscList[[1]][,2][(length(MiscList[[1]][,2]) - (max(XX)-1)):length(MiscList[[1]][,2])]
  
  EstSPR <- YY[length(YY)]
  
  TgSPR <- 0.4
  h <- DLM_data@steep[x]
  SPRLim <- -(2*(h-1))/(3*h+1) # SPR that results in 0.5 R0
  
  phi1 <- 6
  phi2 <- 1
  
  MaxDw <- -0.3
  MaxUp <- 0.3
  
  minSlope <- 0.01
  
  Slope <- coef(lm(YY~XX))[2]  
  # if (abs(Slope) < minSlope) Slope <- 0 
  Dist <- EstSPR - TgSPR 
  
  # Control Rule #
  Mod <- 0 
  Buff <- 0.1
  Buffer <- c(TgSPR - Buff,  TgSPR + Buff)
  inBuff <- FALSE
  belowTG <- FALSE 
  aboveTG <- FALSE
  slopeUp <- FALSE
  slopeDw <- FALSE 
  belowLim <- FALSE
  if (Dist < 0) belowTG <- TRUE 
  if (Dist > 0) aboveTG <- TRUE 
  if (EstSPR > min(Buffer) & EstSPR < max(Buffer)) inBuff <- TRUE
  if (Slope <= 0) slopeDw <- TRUE
  if (Slope > 0) slopeUp <- TRUE
  if (EstSPR < SPRLim) belowLim <- TRUE
   
  # If within buffer zone - only slope
  if (inBuff) Mod <- phi1 * Slope
  if (slopeUp & aboveTG) Mod <- phi1 * Slope +  phi2 * Dist
  if (slopeUp & belowTG) Mod <- phi1 * Slope 
  
  if (slopeDw & aboveTG) Mod <- phi1 * Slope 
  if (slopeDw & belowTG) Mod <- phi1 * Slope +  phi2 * Dist
  
  if (belowLim) Mod <- MaxDw
  
  Mod[Mod > MaxUp] <- MaxUp
  Mod[Mod < MaxDw] <- MaxDw
  Mod <- Mod + 1 
 
  TAC <- DLM_data@MPrec[x] * Mod
  TAC <- TACfilter(TAC)
 
  Out <- list()
  Out[[1]] <- TAC 
  Out[[2]] <- MiscList
 
  return(Out) 
}
class(LBSPR_ItTAC)<-"DLM_output"

# A generic VPA (Walters and Licandeo UBC)
VPA<-function(x, DLM_data, reps=reps) {
   
  # now do optimization for FMSY

  
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@MaxAge, DLM_data@wla, DLM_data@CV_wla, DLM_data@wlb, DLM_data@CV_wlb, DLM_data@L50, DLM_data@CV_L50, DLM_data@CAA, DLM_data@steep, DLM_data@CV_steep, DLM_data@LFS, DLM_data@CV_LFS, DLM_data@LFC, DLM_data@CV_LFC, DLM_data@Cat"
  CAAind<-(DLM_data@CAA[x,,]==0)*array(rep(1:DLM_data@MaxAge,each=length(DLM_data@CAA[x,,1])),dim(DLM_data@CAA[x,,]))
  maxage<-min(CAAind[CAAind!=0]) 
  maxage<-which.min(abs(cumsum(apply(DLM_data@CAA[x,,],2,sum))/sum(DLM_data@CAA[x,,])-0.75))
  CAAv<-DLM_data@CAA[x,,1:maxage]
  CAAv[,maxage]<-CAAv[,maxage]+apply(DLM_data@CAA[x,,(maxage+1):length(DLM_data@CAA[x,1,])],1,sum)
  
  
  
  
  
  TAC<-Bt_K<-rep(NA,reps)
  
  for(i in 1:reps){
  
    Mc<-trlnorm(1,DLM_data@Mort[x],DLM_data@CV_Mort[x])
    hc<-trlnorm(1,DLM_data@steep[x],DLM_data@CV_steep[x])
    Linfc<-trlnorm(1,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
    Kc<-trlnorm(1,DLM_data@vbK[x],DLM_data@CV_vbK[x])
    t0c<--trlnorm(1,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
    LFSc<-trlnorm(1,DLM_data@LFS[x],DLM_data@CV_LFS[x])
    LFCc<-trlnorm(1,DLM_data@LFC[x],DLM_data@CV_LFC[x])
    AMc<-trlnorm(1,iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x]),DLM_data@CV_L50[x])
    ac<-trlnorm(1,DLM_data@wla[x],DLM_data@CV_wla[x])
    bc<-trlnorm(1,DLM_data@wlb[x],DLM_data@CV_wlb[x])
    
    pmat<-rep(1,maxage)
    pmat[1:ceiling(AMc)]<-0
    age<-1:maxage
    la<-DLM_data@vbLinf[x]*(1-exp(-DLM_data@vbK[x]*((age-DLM_data@vbt0[x]))))
    wa<-ac*la^bc
    
    Cat<-DLM_data@Cat[x,]
    Cat[1]<-Cat[2] # temporary fix until effort simulation function gets sorted
    
    
    CAAv[,maxage][CAAv[,maxage]==0]<-1
    
    
    opt=optim(c(-3,-2),VPAopt,Cat=CAAv,yt=DLM_data@Ind[x,],S=exp(-Mc),maxage=maxage,wa=wa,pmat=pmat,method="L-BFGS-B",lower=c(-5,-5),upper=c(5,5))
    out=VPAopt(opt$par,Cat=CAAv,yt=DLM_data@Ind[x,],S=exp(-DLM_data@Mort[x]),maxage=maxage,wa=wa,pmat=pmat,opt=F)
    
    fit2<-optimize(VPAFMSY,log(c(0.0001,3)),Mc=Mc,hc=hc,maxage=maxage,vul=out$va,Linfc=Linfc,Kc=Kc,t0c=t0c,AMc=AMc,ac=ac,bc=bc)
    FMSY<-VPAFMSY(fit2$minimum,Mc=Mc,hc=hc,,maxage=maxage,vul=out$va,Linfc=Linfc,Kc=Kc,t0c=t0c,AMc=AMc,ac=ac,bc=bc,opt=F)
    if ((FMSY / Mc) > 3) FMSY <- 3 * Mc
    TAC[i]<-out$bt[length(out$bt)]*FMSY
  }   

  TACfilter(TAC)
  
}

class(VPA)<-"DLM_output"

#VPAFMSY<-function(lnFMc,Mc,hc,maxage,vul,Linfc,Kc,t0c,AMc,ac,bc,opt=T,ny=50){
  


VPAopt=function(theta,Cat,yt,S,maxage,wa,pmat,opt=T,usewat = F) {  
 
  Uterm<-exp(theta[1])/(1+exp(theta[1]))
  minagecom=1
  minagesel=ceiling(maxage*0.66)
  avg_yrsel<-max(2,(min(5,floor(dim(Cat)[1]/2))))
  
  sig = exp(theta[2]); 
  tiny = 1e-10
  n = dim(Cat)[1]; A=dim(Cat)[2]
  Nat = matrix(NA,n,A)	 # Numbers-at-age matrix
  Ut = rep(NA,length=n)
  ai = 1:(A-2)
  va = c(plogis(1:(A-4),2,0.2),rep(1,4));	# Initial values at the terminal selectivy
  Ut[n] = Uterm
  
  for(j in 1:15) {                  # Numerical convergence to terminal F
    #print(Ut)
    Nat[n,]=Cat[n,]/(Uterm*va)  		# Initialize the terminal year
    
    for(i in (n-1):1) {
      Nat[i,ai]  = Nat[i+1,ai+1]/S + Cat[i,ai]
      Nat[i,A-1] = (Nat[i+1,A]/S + Cat[i,A-1] + Cat[i,A])*( Cat[i,A-1]/( Cat[i,A-1] + Cat[i,A] + tiny) )
      Nat[i,A]   = (Nat[i+1,A]/S + Cat[i,A-1] + Cat[i,A])*( Cat[i,A]  /( Cat[i,A-1] + Cat[i,A] + tiny) ) 
    }
    ##############################################
    ### modify this parameters if need it   #####
    ##############################################
    
    # minagesel = 8
    minagecom = 1
    
    Ut = rowSums(Cat[,minagesel:(A-1)])/rowSums(Nat[,minagesel:(A-1)])  	# Exploitation rate for fully recruited fish
    #    Ut[n] = 0.4
    vat=Cat/Nat/Ut					                     	                        # relative vulnerablility at age
    va=colMeans(vat[(n-avg_yrsel):(n-minagecom),])				                # update terminal vul
    va[minagesel:A] = 1
  
  }
  
  vat[is.na(vat)]<-1
  Ut[n]=Uterm
  if(usewat==T) vbt = rowSums((Nat*vat)*wa) else vbt = (Nat*vat)%*%wa
  if(usewat==T) bt = as.vector(rowSums(Nat * wa)) else bt = as.vector(Nat %*% wa)
  fec = pmat*wa
  zt = log(yt/vbt)
  epsilon = zt - mean(zt)
  if(usewat==T) ssb = as.vector(rowSums(Nat * fec)) else ssb = as.vector(Nat %*% fec)
  predcpue = exp(mean(zt))*vbt ### check again if bt or vbt
  cpue_q = yt/exp(mean(zt))
  qhat = exp(mean(epsilon))
  
  lnl=sum(dnorm(epsilon,mean=0,sd=sig,log=T))
  
  if(opt){
    return(-lnl)
  }else{
  #   ss = sum(epsilon^2)
  #   lnl = 0.5*n*log(ss) 
   return(list(Uterm=Uterm,va=va,rt=Nat[,1],ssb=ssb,yt=yt,vbt=vbt,cpue_q=cpue_q,Nat=Nat,vat=vat
              ,Ut=Ut,bt=bt,predcpue=predcpue,epsilon=epsilon/sig,lnl=lnl,qhat=qhat,
              #               ,Ut=Ut,bt=bt,predcpue=predcpue,lnl=lnl,
              minagesel=minagesel,minagecom=minagecom,avg_yrsel=avg_yrsel))
  }
  
  
}

VPAFMSY<-function(lnFMc,Mc,hc,maxage,vul,Linfc,Kc,t0c,AMc,ac,bc,opt=T,ny=50){
  
  FMc<-exp(lnFMc)
  
  Mac<-rep(1,maxage)
  Mac[1:max(1,floor(AMc))]<-0
  Lac<-Linfc*(1-exp(-Kc*((1:maxage)-t0c)))
  Wac<-ac*Lac^bc
  R0c<-1
  N<-exp(-Mc*((1:maxage)-1))*R0c
  SSN<-Mac*N   # Calculate initial spawning stock numbers
  Biomass<-N*Wac
  SSB<-SSN*Wac                               # Calculate spawning stock biomass
  
  B0<-sum(Biomass)
  SSB0<-sum(SSB)
  SSN0<-SSN
  SSBpR<-sum(SSB)/R0c                              # Calculate spawning stock biomass per recruit
  SSNpR<-SSN/R0c
  
  N<-N/2
  SSN<-Mac*N   # Calculate initial spawning stock numbers
  Biomass<-N*Wac
  SSB<-SSN*Wac
  
  for(y in 1:ny){
    # set up some indices for indexed calculation
    # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
    Zc<-FMc*vul+Mc
    CN<-N*(1-exp(-Zc))*(FMc/Zc)
    CB<-CN*Wac
    Biomass<-N*Wac
    N[2:maxage]<-N[1:(maxage-1)]*exp(-Zc[1:(maxage-1)])         # Total mortality
    N[1]<-(0.8*R0c*hc*sum(SSB))/(0.2*SSBpR*R0c*(1-hc)+(hc-0.2)*sum(SSB))  # Recruitment assuming regional R0 and stock wide steepness
    #print(N[1])
    SSN<-N*Mac
    SSB<-SSN*Wac
    
  } # end of year
  
  if(opt){return(-sum(CB))
  }else{return(FMc)
  }
}

SCA<-function(x,DLM_data,reps=100){ # Requires a character string DLMexe (e.g. "C:/DLMexe") that represents the 
  
  dependencies=""
  
  ny<-dim(DLM_data@CAA)[2]
  na<-dim(DLM_data@CAA)[3]
  
  # write starter file ---------------------------------------------------------------------------------------------------------------------------------
  
  starterfile=paste(DLMexe,"SCA/Starter.ss",sep="/")
  write("SCA.dat",starterfile,1,append=F)
  write("SCA.ctl",starterfile,1,append=T)
  write(0,starterfile,1,append=T) # use init values
  write(1,starterfile,1,append=T) # run display detail (0,1,2)
  write(1,starterfile,1,append=T) # detailed age-structured reports in REPORT.SSO (0,1) 
  
  write(0,starterfile,1,append=T) # write detailed checkup.sso file (0,1)
  write(1,starterfile,1,append=T) # write parm values to ParmTrace.sso
  write(0,starterfile,1,append=T) # write to cumreport.sso (0=no, 1=like$tiemseries; 2=add survey fits)
  
  write(1,starterfile,1,append=T) # include prior like for non-estimated parameters (0,1)
  write(0,starterfile,1,append=T) # use soft boundaries to aid convergence (0,1) (recommended)
  write(3,starterfile,1,append=T) # Number of data files to produce: 1st is input, 2nd is estimates, 3rd and higher are bootstrap
  
  write(10,starterfile,1,append=T) # Turn off estimation for parameters entering after this phase
  write(10,starterfile,1,append=T) # MCeval burn interval
  write(2,starterfile,1,append=T) # MCeval thin interval
  
  write(0,starterfile,1,append=T) # jitter initial parm value by this fraction
  write(-1,starterfile,1,append=T) # min yr for sdreport outputs (-1 for styr)
  write(-2,starterfile,1,append=T) # max yr for sdreport outputs
  
  write(0,starterfile,1,append=T) # N individual STD years 
  write(0.001,starterfile,1,append=T) # final convergence criteria (e.g. 1.0e-04)
  write(0,starterfile,1,append=T) # retrospective year relative to end year (e.g. -4)
  
  write(1,starterfile,1,append=T) # min age for calc of summary biomass
  write(2,starterfile,1,append=T) # Depletion basis:  denom is: 0=skip; 1=rel X*B0; 2=rel X*Bmsy; 3=rel X*B_styr
  # !!!!!!!!!!!!!!!! 
  write("1.0",starterfile,1,append=T) # Fraction (X) for Depletion denominator (e.g. 0.4)
  
  write(2,starterfile,1,append=T) # SPR_report_basis:  0=skip; 1=(1-SPR)/(1-SPR_tgt); 2=(1-SPR)/(1-SPR_MSY); 3=(1-SPR)/(1-SPR_Btarget); 4=rawSPR
  write(4,starterfile,1,append=T) # F_report_units: 0=skip; 1=exploitation(Bio); 2=exploitation(Num); 3=sum(Frates); 4=true F for range of ages
  write(c(floor(max(1,DLM_data@MaxAge/4)), ceiling(DLM_data@MaxAge/3)),starterfile,2,append=T) #_min and max age over which average F will be calculated
  
  write(2,starterfile,1,append=T) # F_report_basis: 0=raw; 1=F/Fspr; 2=F/Fmsy ; 3=F/Fbtgt
  write(999,starterfile,1,append=T) # check value for end of file
  
  
  
  # write control file ---------------------------------------------------------------------------------------------------------------------------------
  
  ctlfile=paste(DLMexe,"SCA/SCA.ctl",sep="/")
  
  write(1,ctlfile,1,append=F) #_N_Growth_Patterns
  write(1,ctlfile,1,append=T) #_N_Morphs_Within_GrowthPattern
  write(0,ctlfile,1,append=T) #_Nblock_Patterns
  write(0.5,ctlfile,1,append=T) #_fracfemale
  write(0,ctlfile,1,append=T) #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate
  
  write(1,ctlfile,1,append=T) # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_speciific_K; 4=not implemented
  write(0,ctlfile,1,append=T) #_Growth_Age_for_L1  !!!!!!!
  write(999,ctlfile,1,append=T) #_Growth_Age_for_L2 (999 to use as Linf)
  
  write(0,ctlfile,1,append=T) #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
  write(0,ctlfile,1,append=T) #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
  write(1,ctlfile,1,append=T) #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=read fec and wt from wtatage.ss
  
  write(1,ctlfile,1,append=T) #_First_Mature_Age
  # !!!!!!!!!!!!!!!!!!!!
  write(1,ctlfile,1,append=T) #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
  write(0,ctlfile,1,append=T) #_hermaphroditism option:  0=none; 1=age-specific fxn
  
  write(1,ctlfile,1,append=T) #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
  # !!!!!!!!!!!!!!!!!!!!
  write(2,ctlfile,1,append=T) #_env/block/dev_adjust_method (1=standard; 2=logistic transform keeps in base parm bounds; 3=standard w/ no bound check)
  
  #female
  write(paste(round(DLM_data@Mort[x]*0.9,3), round(DLM_data@Mort[x]*1.1,3),round(DLM_data@Mort[x],3),round(DLM_data@Mort[x],3),-1,0.1,-3,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)                 # NatM_p_1_Fem_GP_1
  write(paste(0,0,0,0,-1,10,-3,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)                                                                                       # L_at_Amin_Fem_GP_1
  write(paste(round(DLM_data@vbLinf[x]*0.9,3),round(DLM_data@vbLinf[x]*1.1,3), round(DLM_data@vbLinf[x],3),round(DLM_data@vbLinf[x],3),-1,0.1,-3,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)         # L_at_Amax_Fem_GP_1 
  
  write(paste(round(DLM_data@vbK[x]*0.9,3),round(DLM_data@vbK[x]*1.1,3), round(DLM_data@vbK[x],3),round(DLM_data@vbK[x],3),-1,0.1,-3,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)                     # VonBert_K_Fem_GP_1
  write(paste(0.05,0.25,0.1,0.1,-1,0.8,-3,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)                                                                            # CV_young_Fem_GP_1
  write(paste(0.05,0.25,0.1,0.1,-1,0.8,-3,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)                                                                            # CV_old_Fem_GP_1
  
  write(paste(DLM_data@wla[x]*0.9, DLM_data@wla[x]*1.1,DLM_data@wla[x],DLM_data@wla[x],-1,0.1,-3,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)                     # Wtlen_1_Fem
  write(paste(DLM_data@wlb[x]*0.9, DLM_data@wlb[x]*1.1,DLM_data@wlb[x],DLM_data@wlb[x],-1,0.1,-3,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)                     # Wtlen_2_Fem
  
  write(paste(round(DLM_data@L50[x]*0.9,3), round(DLM_data@L50[x]*1.1,3),round(DLM_data@L50[x],3),round(DLM_data@L50[x],3),-1,0.1,-3,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)                     # Mat50%_Fem
  Lsd=-2.94439/(DLM_data@L95[x]-DLM_data@L50[x]) # slope=log(1/0.95-1)/(L95-L50)
  write(paste(-3,3,round(Lsd,3),round(Lsd,3),-1,0.8,-3,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)                                                                                 # Mat_slope_Fem
  
  #write(paste(DLM_data@wla[x]*0.9, DLM_data@wla[x]*1.1,DLM_data@wla[x],DLM_data@wla[x],-1,0.1,-3,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)                     # Wtlen_1_Mal
  #write(paste(DLM_data@wlb[x]*0.9, DLM_data@wlb[x]*1.1,DLM_data@wlb[x],DLM_data@wlb[x],-1,0.1,-3,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)                     # Wtlen_2_Mal
  
  write(paste(0,3,1,1,-1,0.8,-3,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)                                                                                        # RecrDist_Area_1
  write(paste(0,3,0,0,-1,0.8,-3,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)                                                                                        # RecrDist_Seas_1
  
  write(paste(0,0,0,0,-1,0,-4,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)                                                                                        # RecrDist_Area_1
  write(paste(0,0,0,0,-1,0,-4,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)                                                                                        # RecrDist_Seas_1
  write(paste(0,0,0,0,-1,0,-4,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)                                                                                        # CohortGrowDev
  
  write(paste(0.5,1.5,1,1,-1,0,-4,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)                                                                                        # CohortGrowDev
  
  
  write(paste(0,0,0,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)    #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
  
  
  # Reconstruct a possible catch-at-age matrix to get to R0 given M
  la<-DLM_data@vbLinf[x]*(1-exp(-DLM_data@vbK[x]*(((1:na)-DLM_data@vbt0[x]))))
  wa<-DLM_data@wla[x]*la^DLM_data@wlb[x]
  Cw<-t(array(wa,c(na,ny)))*DLM_data@CAA[x,,] # weight of the observed CAA
  Cwtot<-apply(Cw,1,sum) # summation by year
  CAAup<-DLM_data@CAA[x,,]*DLM_data@Cat[x,]/Cwtot # uprate CAA to total catch weight
  c1<-apply(CAAup,2,mean)
  plusgroup<-which.min((cumsum(c1)/sum(c1)-0.95)^2)
  c1[plusgroup]<-c1[plusgroup]+sum(c1[(plusgroup+1):na])
  c1<-c1[plusgroup:1]
  aa<-rep(NA,plusgroup)
  aa[1]<-c1[1]
  for(i in 2:plusgroup)aa[i]=aa[i-1]/exp(-DLM_data@Mort[x])+c1[i]
  R0est<-aa[plusgroup]
  #_Spawner-Recruitment
  write(3,ctlfile,1,append=T)                                     #_SR_function: 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm
  #_LO HI INIT PRIOR PR_type SD PHASE
  
  write(paste(round(log(R0est/10),3),round(log(R0est*10),3),round(log(R0est),3),round(log(R0est),3),-1,10,1,sep=" "),ctlfile,1,append=T)                                      # SR_LN(R0)
  write(paste(round(DLM_data@steep[x]*0.9,3),round(DLM_data@steep[x],3), round(DLM_data@steep[x],3),round(DLM_data@steep[x],3),1,0.04,-3,sep=" "),ctlfile,1,append=T)         # SR_BH_steep
  write(paste(0,2,0.6,0.8,-1,0.8,-4,sep=" "),ctlfile,1,append=T)                                                                          # SR_sigmaR
  
  write(paste(-5,5,0,0,-1,1,-3,sep=" "),ctlfile,1,append=T)                                                                             # SR_envlink
  write(paste(-5,5,0,0,-1,1,-2,sep=" "),ctlfile,1,append=T)                                                                               # SR_R1_offset
  write(paste(0,0,0,0,-1,0,-99,sep=" "),ctlfile,1,append=T)                                                                               # SR_autocorr
  
  write(0,ctlfile,1,append=T)                                     #_SR_env_link
  write(0,ctlfile,1,append=T)                                     #_SR_env_target_0=none;1=devs;_2=R0;_3=steepness
  write(1,ctlfile,1,append=T)                                     # do_recdev
  write(1,ctlfile,1,append=T)                                     # first year of main recr_devs; early devs can preceed this era
  write(ny,ctlfile,1,append=T)                                    # last year of main recr_devs; forecast devs start in following year
  write(2,ctlfile,1,append=T)                                     #_recdev phase 
  
  write(1,ctlfile,1,append=T)                                     # (0/1) to read 13 advanced options
  
  write(-na+2,ctlfile,1,append=T)                                     #
  write(4,ctlfile,1,append=T)                                    #
  write(0,ctlfile,1,append=T)                                     #
  write(1,ctlfile,1,append=T)                                     #
  write(1,ctlfile,1,append=T) 
  write(10,ctlfile,1,append=T)                                     #
  write(ny-8,ctlfile,1,append=T)                                  #
  write(ny-1,ctlfile,1,append=T)                                    #
  write(0.8,ctlfile,1,append=T)                                     #
  write(0,ctlfile,1,append=T)                                     #
  write(-5,ctlfile,1,append=T)                                    #
  write(5,ctlfile,1,append=T)                                     #
  write(0,ctlfile,1,append=T)                                     #
  
  
  write(DLM_data@Mort[x],ctlfile,1,append=T)                      # F ballpark for tuning early phases
  write(-ny,ctlfile,1,append=T)                                   # F ballpark year (neg value to disable)
  write(3,ctlfile,1,append=T)                                     # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
  write(2.9,ctlfile,1,append=T)                                   # max F or harvest rate, depends on F_Method
  write(4,ctlfile,1,append=T)                                     # N iterations for tuning F in hybrid method (recommend 3 to 7)
  
  #_LO HI INIT PRIOR PR_type SD PHASE
  write(paste(0,round(DLM_data@Mort[x]*3,3),0,0.01,0,99,-1,sep=" "),ctlfile,1,append=T)     #InitF_1FISHERY1
  
  #_Den-dep  env-var  extra_se  Q_type
  write(paste(0,0,0,0,sep=" "),ctlfile,1,append=T)          # FISHERY                   
  write(paste(0,0,0,2,sep=" "),ctlfile,1,append=T)          # SURVEY                    
  
  LOq<-1/(mean(DLM_data@Cat[x,])/(DLM_data@Mort[x]/10))
  HIq<-1/(mean(DLM_data@Cat[x,])/(DLM_data@Mort[x]*5))
  muq<-1/(mean(DLM_data@Cat[x,])/DLM_data@Mort[x])
  # LO HI INIT PRIOR PR_type SD PHASE
  write(paste(round(log(LOq),3),round(log(HIq),3),round(log(muq),3),round(log(muq),3),0,1,-1,sep=" "),ctlfile,1,append=T)      # Q_base_SURVEY
  
  #_Pattern Discard Male Special
  write(paste(0,0,0,0,sep=" "),ctlfile,1,append=T)      # FISHERY                                                                         
  write(paste(0,0,0,0,sep=" "),ctlfile,1,append=T)      # SURVEY                                                                         
  
  #_age_selex_types
  write(paste(12,0,0,0,sep=" "),ctlfile,1,append=T)     # FISHERY                                                                         
  write(paste(12,0,0,0,sep=" "),ctlfile,1,append=T)     # SURVEY                                                                         
  
  #LO HI INIT PRIOR PR_type SD PHASE env-var use_dev dev_minyr dev_maxyr dev_stddev Block Block_Fxn
  write(paste(1,floor(na*0.666),5,5,-1,10,1,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)      # AgeSel 1 Fishery     
  write(paste(0,0.5,  0.25,0.25,-1,10,1,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)        # AgeSel 2 Fishery     
  write(paste(1,floor(na*0.666),5,5,-1,10,1,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)      # AgeSel 1 Fishery     
  write(paste(0,0.5,  0.25,0.25,-1,10,1,0,0,0,0,0,0,0,sep=" "),ctlfile,1,append=T)        # AgeSel 2 Fishery     
  
  # write(paste(0,floor(na*0.666),floor(na/3),floor(na/3),0,2,2,0,0,0,0,0.5,0,0,sep=" "),ctlfile,1,append=T)      # AgeSel 1 SURVEY     
  #  write(paste(0.01,floor(na*0.8),floor(na/2),floor(na/2),0,2,2,0,0,0,0,0.5,0,0,sep=" "),ctlfile,1,append=T)      # AgeSel 1 SURVEY     
  
  # Tag loss and Tag reporting parameters go next
  write(0,ctlfile,1,append=T)  # TG_custom:  0=no read; 1=read if tags exist    
  
  write(0,ctlfile,1,append=T)#_Variance_adjustments_to_input_values
  #_fleet: 1 2 3 
  #write(0,ctlfile,1,append=T)  #_add_to_survey_CV
  #write(0,ctlfile,1,append=T)  #_add_to_discard_stddev
  #write(0,ctlfile,1,append=T)  #_add_to_bodywt_CV
  #write(1,ctlfile,1,append=T)  #_mult_by_lencomp_N 
  #write(1,ctlfile,1,append=T)  #_mult_by_agecomp_N
  #write(1,ctlfile,1,append=T)  #_mult_by_size-at-age_N
  
  write(2,ctlfile,1,append=T)  #_maxlambdaphase
  write(1,ctlfile,1,append=T)  #_sd_offset
  write(1,ctlfile,1,append=T)  # number of changes to make to default Lambdas (default value is 1.0)
  write(paste(1,2,1,0,1,sep=" "),ctlfile,1,append=T)     # Like_comp SURVEY                                                                         
  
  # 1 2 2 1 1
  # 4 2 2 1 1
  # 4 2 3 1 1
  
  write(0,ctlfile,1,append=T)  # (0/1) read specs for more stddev reporting 
  write(999,ctlfile,1,append=T)  # end-of-file
  
  
  
  # write data file ---------------------------------------------------------------------------------------------------------------------------------
  
  datfile=paste(DLMexe,"SCA/SCA.dat",sep="/")
  
  write(1,datfile,1,append=F)     #_styr
  write(ny,datfile,1,append=T)    #_endyr
  write(1,datfile,1,append=T)     #_nseas
  write(12,datfile,1,append=T)    # months/season
  write(1,datfile,1,append=T)     #_spawn_seas
  write(1,datfile,1,append=T)     #_Nfleet
  write(1,datfile,1,append=T)     #_Nsurveys
  write(1,datfile,1,append=T)     #_N_areas
  write("FISHERY%SURVEY",datfile,1,append=T)
  write(paste(c(-1,0.5),collapse=" "),datfile,1,append=T)   #_surveytiming_in_season
  write(paste(c(1,1),collapse=" "),datfile,1,append=T)     #_area_assignments_for_each_fishery_and_survey
  write(1,datfile,1,append=T)     #_units of catch:  1=bio; 2=num
  write(0.1,datfile,1,append=T)   #_se of log(catch) only used for init_eq_catch and for Fmethod 2 and 3; use -1 for discard only fleets
  write(1,datfile,1,append=T)     #_Ngenders
  write(na,datfile,1,append=T)    #_Nages
  write(0,datfile,1,append=T)     #_init_equil_catch_for_each_fishery
  write(ny,datfile,1,append=T)    #_N_lines_of_catch_to_read
  #_catch_biomass(mtons):_columns_are_fisheries,year,season
  for(i in 1:ny) write(paste(round(DLM_data@Cat[x,i],2),i,1,sep=" "),datfile,1,append=T)   
  
  write(ny,datfile,1,append=T)    #_N_cpue_and_surveyabundance_observation
  #_Units:  0=numbers; 1=biomass; 2=F
  #_Errtype:  -1=normal; 0=lognormal; >0=T
  #_Fleet Units Errtype
  write(paste(1,1,0,sep=" "),datfile,1,append=T)     # SURVEY                                                                         
  write(paste(2,1,0,sep=" "),datfile,1,append=T)     # FISHERY                                                                         
  
  #_year seas index obs err
  for(i in 1:ny) write(paste(i,1,2,round(DLM_data@Ind[x,i],2),0.2,sep=" "),datfile,1,append=T)    # index 2 is SURVEY
  
  write(0,datfile,1,append=T)     #_N_fleets_with_discard
  write(0,datfile,1,append=T)     # N discard obs
  write(0,datfile,1,append=T)     #_N_meanbodywt_obs
  write(30,datfile,1,append=T)    #_DF_for_meanbodywt_T-distribution_like
  
  write(2,datfile,1,append=T)     # length bin method: 1=use databins; 2=generate from binwidth,min,max below; 3=read vector
  write(2,datfile,1,append=T)     # binwidth for population size comp 
  write(DLM_data@CAL_bins[2],datfile,1,append=T)    # minimum size in the population (lower edge of first bin and size at age 0.00) 
  write(DLM_data@CAL_bins[length(DLM_data@CAL_bins)],datfile,1,append=T)    # maximum size in the population (lower edge of last bin) 
  
  write(0.0001,datfile,1,append=T)    #_comp_tail_compression
  write(1e-007,datfile,1,append=T)    #_add_to_comp
  write(0,datfile,1,append=T)    #_combine males into females at or below this bin number
  write(length(DLM_data@CAL_bins),datfile,1,append=T)    #_N_LengthBins
  write(paste(DLM_data@CAL_bins,collapse=" "),datfile,1,append=T)    # Length bins
  write(0,datfile,1,append=T)    #_N_Length_obs
  
  
  maxage<-which.min(abs(cumsum(apply(DLM_data@CAA[x,,],2,sum))/sum(DLM_data@CAA[x,,])-0.95)) # plus group at cumulative 75th percentile
  CAA<-DLM_data@CAA[x,,1:maxage]
  CAA[,maxage]<-CAA[,maxage]+apply(DLM_data@CAA[x,,maxage:na],1,sum)
  
  write(maxage,datfile,1,append=T)    #_N_age_bins
  write(paste(1:maxage,collapse=" "),datfile,1,append=T)    #age_bins
  
  write(1,datfile,1,append=T)    #_N_ageerror_definitions
  write(paste(rep(-1,na+1),collapse=" "),datfile,1,append=T)
  write(paste(rep(0.001,na+1),collapse=" "),datfile,1,append=T)
  
  write(ny*2,datfile,1,append=T)   #_N_Agecomp_obs
  write(1,datfile,1,append=T)    #_Lbin_method: 1=poplenbins; 2=datalenbins; 3=lengths
  write(1,datfile,1,append=T)    #_combine males into females at or below this bin number
  
  #Yr Seas Flt/Svy Gender Part Ageerr Lbin_lo Lbin_hi Nsamp datavector(female-male)
  
  for(i in 1:ny) write(paste(" ",i,1,1,0,0,1,-1,-1,sum(CAA[i,]),paste(CAA[i,],collapse=" "),sep=" "),datfile,1,append=T)   
  for(i in 1:ny) write(paste(" ",i,1,2,0,0,1,-1,-1,sum(CAA[i,]),paste(CAA[i,],collapse=" "),sep=" "),datfile,1,append=T)   
  
  write(0,datfile,1,append=T)      #_N_MeanSize-at-Age_obs
  write(0,datfile,1,append=T)      #_N_environ_variables
  write(0,datfile,1,append=T)      #_N_environ_obs
  write(0,datfile,1,append=T)      # N sizefreq methods to read 
  write(0,datfile,1,append=T)      # no tag data 
  write(0,datfile,1,append=T)      # no morphcomp data
  write(999,datfile,1,append=T)    # end of file
  
  
  
  # write forecast file ---------------------------------------------------------------------------------------------------------------------------------
  
  forefile=paste(DLMexe,"SCA/Forecast.ss",sep="/")
  
  write(1,forefile,1,append=F)       # Benchmarks: 0=skip; 1=calc F_spr,F_btgt,F_msy
  write(2,forefile,1,append=T)       # MSY: 1= set to F(SPR); 2=calc F(MSY); 3=set to F(Btgt); 4=set to F(endyr)
  write(0.40,forefile,1,append=T)    # SPR target (e.g. 0.40)
  write(0.35,forefile,1,append=T)    # Biomass target (e.g. 0.40)
  #_Bmark_years: beg_bio, end_bio, beg_selex, end_selex, beg_relF, end_relF (enter actual year, or values of 0 or -integer to be rel. endyr)
  write(paste(0,0,0,0,0,0,sep=" "),forefile,1,append=T)     # FISHERY                                                                         
  write(1,forefile,1,append=T)       #Bmark_relF_Basis: 1 = use year range; 2 = set relF same as forecast below
  write(2,forefile,1,append=T)       # Forecast: 0=none; 1=F(SPR); 2=F(MSY) 3=F(Btgt); 4=Ave F (uses first-last relF yrs); 5=input annual F scalar
  write(3,forefile,1,append=T)       # N Forecast years
  write(0.2,forefile,1,append=T)     # F scalar (only used for Do_Forecast==5) 
  #_Fcast_years:  beg_selex, end_selex, beg_relF, end_relF  (enter actual year, or values of 0 or -integer to be rel. endyr)
  write(paste(0,0,-10,0,sep=" "),forefile,1,append=T)     # FISHERY                                                                         
  write(1,forefile,1,append=T)       # Control rule method (1=catch=f(SSB) west coast; 2=F=f(SSB) )                                                                       
  write(0.4,forefile,1,append=T)     # Control rule Biomass level for constant F (as frac of Bzero, e.g. 0.40)                                                                      
  write(0.1,forefile,1,append=T)     # Control rule Biomass level for no F (as frac of Bzero, e.g. 0.10)
  write(0.75,forefile,1,append=T)    # Control rule target as fraction of Flimit (e.g. 0.75)
  write(3,forefile,1,append=T)       #_N forecast loops (1-3) (fixed at 3 for now)
  write(3,forefile,1,append=T)       #_First forecast loop with stochastic recruitment
  write(0,forefile,1,append=T)       #_Forecast loop control #3 (reserved for future bells&whistles)
  write(0,forefile,1,append=T)       #_Forecast loop control #4 (reserved for future bells&whistles)
  write(0,forefile,1,append=T)       #_Forecast loop control #5 (reserved for future bells&whistles)
  write(ny,forefile,1,append=T)      #FirstYear for caps and allocations (should be after years with fixed inputs)
  write(0.0,forefile,1,append=T)     # stddev of log(realized catch/target catch) in forecast (set value>0.0 to cause active impl_error)
  write(0,forefile,1,append=T)       # Do West Coast gfish rebuilder output (0/1)
  write(ny-10,forefile,1,append=T)   # Rebuilder:  first year catch could have been set to zero (Ydecl)(-1 to set to 1999)
  write(ny-5,forefile,1,append=T)    # Rebuilder:  year for current age structure (Yinit) (-1 to set to endyear+1)
  write(1,forefile,1,append=T)       # fleet relative F:  1=use first-last alloc year; 2=read seas(row) x fleet(col) below
  write(2,forefile,1,append=T)       # basis for fcast catch tuning and for fcast catch caps and allocation  (2=deadbio; 3=retainbio; 5=deadnum; 6=retainnum)
  write(-1,forefile,1,append=T)      # max totalcatch by fleet (-1 to have no max)
  write(-1,forefile,1,append=T)      # max totalcatch by area (-1 to have no max)
  write(0,forefile,1,append=T)       # fleet assignment to allocation group (enter group ID# for each fleet, 0 for not included in an alloc group)
  write(0,forefile,1,append=T)       # Number of forecast catch levels to input (else calc catch from forecast F)
  write(2,forefile,1,append=T)       # basis for input Fcast catch:  2=dead catch; 3=retained catch; 99=input Hrate(F) (units are from fleetunits; note new codes in SSV3.20)
  write(999,forefile,1,append=T)     # End of input
  
  
  
  # Run SS3 and read outputs ----------------------------------------------------------------------------------------------------------------------------
  
  system(paste(DLMexe,"/SCA/SS3.exe",sep=""),wait=T,show.output.on.console=T)
  rl <- SS_output(dir=paste(DLMexe,"/SCA/",sep=""))
  F_FMSY<-rl$Kobe$F.Fmsy[1:ny]
  B_BMSY<-rl$Kobe$B.Bmsy[1:ny]
  ind<-match(paste("F_",1:ny,sep=""),row.names(rl$derived_quants))
  Fs<-rl$derived_quants[ind,2]
  Fs/F_FMSY
  
  #SS_plots(replist=myreplist)
  
  
  
}
class(SCA)<-"DLM_output"




