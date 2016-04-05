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

matlenlim2 <-function(x,DLM_data, ...){ # Knife-edge vulnerability slightly higher than length at maturity 
  dependencies="DLM_data@LFC, DLM_data@LFS"
  Allocate<-1
  Effort<-1
  Spatial<-c(1,1)
  newLFS <- 14 + 1.1 * DLM_data@L50[x] 
  newLFC <- 0.95 * newLFS
  Vuln <-c(newLFC, newLFS)
  c(Allocate, Effort, Spatial, Vuln)
}
class(matlenlim2)<-"DLM_input"

slotlim <-function(x,DLM_data, ...){ # Example of slot limit between 0.95 and 1.25 * L50 
  dependencies="DLM_data@LFC, DLM_data@LFS"
  Allocate<-1
  Effort<-1
  Spatial<-c(1,1)
  
  newLFS <- 14 + 1.1 * DLM_data@L50[x]
  newLFC <- 0.95 * newLFS  
  UppLim <- as.numeric(quantile(c(newLFS, DLM_data@vbLinf[x]), 0.75))
  Vuln <-c(newLFC, newLFS, UppLim)
  c(Allocate, Effort, Spatial, Vuln)
}
class(slotlim)<-"DLM_input"

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

LBSPR_ItEff <- function(x, DLM_data, yrsmth=1, reps=reps) {
 dependencies="DLM_data@CAL, DLM_data@CAL_bins, DLM_data@vbLinf, 
	DLM_data@vbK, DLM_data@Mort, LM_data@vbK, DLM_data@L50, DLM_data@L95, 
	DLM_data@wlb"
  MiscList <- LBSPR(x, DLM_data, yrsmth=yrsmth, reps=reps)
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
  
  Allocate <- 1
  if (is.na(DLM_data@MPrec[x])) DLM_data@MPrec[x] <- 1 
  Effort <- DLM_data@MPrec[x] * Mod
  MiscList[[2]] <- append(MiscList[[2]], Effort)
  Spatial <- c(1,1)
  Vuln <- rep(NA,2)
  out <- c(Allocate, Effort, Spatial, Vuln)
   
  Out <- list()
  Out[[1]] <- out 
  Out[[2]] <- MiscList
 
  return(Out) 
}
class(LBSPR_ItEff)<-"DLM_input"

LBSPR_ItSel <- function(x, DLM_data, yrsmth=1, reps=reps) {
 dependencies="DLM_data@CAL, DLM_data@CAL_bins, DLM_data@vbLinf, 
	DLM_data@vbK, DLM_data@Mort, LM_data@vbK, DLM_data@L50, DLM_data@L95, 
	DLM_data@wlb"
  MiscList <- LBSPR(x, DLM_data, yrsmth=yrsmth,reps=reps)
  if(all(is.na(MiscList[[1]]))) return(rep(NA, 6))
  if(all(is.na(MiscList[[1]][,2]))) return(rep(NA, 6))
  XX <- 1:4 
  YY <- MiscList[[1]][,2][(length(MiscList[[1]][,2]) - (max(XX)-1)):length(MiscList[[1]][,2])]
  
  EstSPR <- YY[length(YY)]
  
  TgSPR <- 0.4
  h <- DLM_data@steep[x]
  SPRLim <- -(2*(h-1))/(3*h+1) # SPR that results in 0.5 R0
 
  Allocate <- 1
  Effort <- 1
  Spatial <- c(1,1)

  if (EstSPR < TgSPR) {
    newLFC <- DLM_data@L50[x] * 1.05
    newLFS <- DLM_data@L50[x] * 1.1
    Vuln <-c(newLFC, newLFS)
  }
  if (EstSPR < SPRLim) {
    newLFC <- DLM_data@L50[x] * 1.2
    newLFS <- DLM_data@L50[x] * 1.25
    Vuln <-c(newLFC, newLFS)  
  }
  if (EstSPR >= TgSPR) {
    newLFC <- DLM_data@L50[x] * 0.85
    newLFS <- DLM_data@L50[x] * 0.9
    Vuln <-c(newLFC, newLFS)  
  }
   
 
  out <- c(Allocate, Effort, Spatial, Vuln)
   
  Out <- list()
  Out[[1]] <- out 
  Out[[2]] <- MiscList
 
  return(Out) 
}
class(LBSPR_ItSel)<-"DLM_input"

DDe<-function(x,DLM_data,reps=100){
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
  
  U_hist<-1-exp(-exp(opt$par[3])*E_hist)
 
  Allocate <- 1
  Effort<-exp(opt$par[1])/U_hist[DLM_data@LHYear]
  Spatial <- c(1,1)
  Vuln<-rep(NA,3)
  out <- c(Allocate, Effort, Spatial, Vuln)
  
  #Out <- list()
  #Out[[1]] <- out 
  #Out[[2]] <- MiscList
  
  return(out) 

}
class(DDe)<-"DLM_input"

DDe75<-function(x,DLM_data,reps=100){
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
  
  U_hist<-1-exp(-exp(opt$par[3])*E_hist)
  
  Allocate <- 1
  Effort<-0.75*exp(opt$par[1])/U_hist[DLM_data@LHYear]
  Spatial <- c(1,1)
  Vuln<-rep(NA,3)
  out <- c(Allocate, Effort, Spatial, Vuln)
  
  #Out <- list()
  #Out[[1]] <- out 
  #Out[[2]] <- MiscList
  
  return(out) 
  
}
class(DDe75)<-"DLM_input"

DDes<-function(x,DLM_data,reps=100,LB=0.9,UB=1.1){
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
  
  U_hist<-1-exp(-exp(opt$par[3])*E_hist)
  fac<-exp(opt$par[1])/U_hist[DLM_data@LHYear] # ratio of UMSY to reference U
  fac<-fac*(U_hist[DLM_data@LHYear]/U_hist[length(U_hist)]) # ratio of last U to reference U
 
  if(fac<LB)fac<-LB
  if(fac>UB)fac<-UB
  
  Allocate <- 1
  Effort<-DLM_data@MPrec[x]*fac
  Spatial <- c(1,1)
  Vuln<-rep(NA,3)
  out <- c(Allocate, Effort, Spatial, Vuln)
  
  #Out <- list()
  #Out[[1]] <- out 
  #Out[[2]] <- MiscList
  
  return(out) 
  
}
class(DDes)<-"DLM_input"


DTe40<-function(x,DLM_data,reps=100,alpha=0.4,LB=0.9,UB=1.1){
  
  dependencies="DLM_data@Dep"
  
  fac<-DLM_data@Dep[x]/alpha
  if(fac<LB)fac<-LB
  if(fac>UB)fac<-UB
  
  Allocate <- 1
  Effort<-DLM_data@MPrec[x]*fac 
  Spatial <- c(1,1)
  Vuln<-rep(NA,3)
  out <- c(Allocate, Effort, Spatial, Vuln)
  
  #MiscList<-Effort
  
  #Out <- list()
  #Out[[1]] <- out 
  #Out[[2]] <- MiscList
  
  return(out) 
  
}
class(DTe40)<-"DLM_input"


DTe50<-function(x,DLM_data,reps=100,alpha=0.5,LB=0.9,UB=1.1){
  
  dependencies="DLM_data@Dep"
  
  fac<-DLM_data@Dep[x]/alpha
  if(fac<LB)fac<-LB
  if(fac>UB)fac<-UB
  
  Allocate <- 1
  Effort<-DLM_data@MPrec[x]*fac 
  Spatial <- c(1,1)
  Vuln<-rep(NA,3)
  out <- c(Allocate, Effort, Spatial, Vuln)
  
  #MiscList<-Effort
  
  #Out <- list()
  #Out[[1]] <- out 
  #Out[[2]] <- MiscList
  
  return(out) 
  
}
class(DTe50)<-"DLM_input"


