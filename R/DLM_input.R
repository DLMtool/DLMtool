# DLM_input MPs

matlenlim<-function(x,DLM_data, ...){ # Length at maturity is knife-edge vulnerability
  dependencies="DLM_data@LFC, DLM_data@LFS"
  Allocate<-1
  Effort<-1
  Spatial<-c(1,1)
  
  newLFC <- DLM_data@L50[x] * 0.95
  newLFS <- DLM_data@L50[x] 
  Vuln <-c(newLFC, newLFS)
  c(Allocate, Effort, Spatial, Vuln)
}
class(matlenlim)<-"DLM_input"

# matagelim<-function(x,DLM_data, ...){ # Age at maturity is knife-edge vulnerability
  # dependencies="DLM_data@AM, DLM_data@MaxAge"
  # Allocate<-1
  # Effort<-1
  # Spatial<-c(1,1)
  # Vuln<-1/(1+exp((DLM_data@AM[x]-(1:DLM_data@MaxAge))/(DLM_data@AM[x]*DLM_data@AM[x]*0.05)))
  # c(Allocate, Effort, Spatial, Vuln)
# }
# class(matagelim)<-"DLM_input"

MRreal<-function(x,DLM_data, ...){ # A Marine reserve in area 1 with spatial reallocation of effort
  dependencies="DLM_data@MaxAge"
  Allocate<-1
  Effort<-1
  Spatial<-c(0,1)
  # Vuln<-rep(NA,DLM_data@MaxAge)
  Vuln<-rep(NA,2)
  c(Allocate, Effort, Spatial, Vuln)
}
class(MRreal)<-"DLM_input"

MRnoreal<-function(x,DLM_data, ...){ # A Marine reserve in area 1 with no spatial reallocation of effort
  dependencies="DLM_data@MaxAge"
  Allocate<-0
  Effort<-1
  Spatial<-c(0,1)
  # Vuln<-rep(NA,DLM_data@MaxAge)
  Vuln<-rep(NA,2)
  c(Allocate, Effort, Spatial, Vuln)
}
class(MRnoreal)<-"DLM_input"

curE<-function(x,DLM_data, ...){ # current effort
  dependencies="DLM_data@MaxAge"
  Allocate<-1
  Effort<-1
  Spatial<-c(1,1)
  # Vuln<-rep(NA,DLM_data@MaxAge)
  Vuln<-rep(NA,2)  
  c(Allocate, Effort, Spatial, Vuln)
}
class(curE)<-"DLM_input"

curE75<-function(x,DLM_data, ...){ #75% current effort
  dependencies="DLM_data@MaxAge"
  Allocate<-1
  Effort<-0.75
  Spatial<-c(1,1)
  # Vuln<-rep(NA,DLM_data@MaxAge)
  Vuln<-rep(NA,2)  
  c(Allocate, Effort, Spatial, Vuln)
}
class(curE75)<-"DLM_input"


## - Test Input controls

YPR_CC_Input <-function(x,DLM_data,reps=100, Fmin=0.005){
  #for(x in 1:16){
    dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@MaxAge, DLM_data@wla, DLM_data@wlb, DLM_data@CAA, DLM_data@Cat"
  Linfc<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  t0c<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
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
  
  Fdb # Estimated F
 
  Allocate <- 1
  Effort<-0.75
  Spatial<-c(1,1)
  Vuln<-rep(NA,2)
 
  if (median(Fdb) > 0.5) Effort <- 0.5 
  
  c(Allocate, Effort, Spatial, Vuln)
  
}
class(YPR_CC_Input)<-"DLM_input"


