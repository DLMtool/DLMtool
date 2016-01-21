
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
