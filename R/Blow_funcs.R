# Blow optimization code

#' Blow parallel optimization function
#'
#' Find the current biomass at which it would take HZN mean generation times
#' to reach Bfrac x SSBMSY biomass level given zero catches
#'
#' @param x position in a vector
#' @param SSBMSY vector nsim long of spawning biomass at MSY
#' @param MGThorizon vector nsim long of MGT x HZN
#' @param Find matrix of fishing mortality rate nsim x nyears 
#' @param Perr matrix of recruitment devitions nsim x nyears + maxage -1
#' @param M_ageArray array of natural mortality rate nsim x maxage x nyears + proyears
#' @param hs vector nsim long of steepness values
#' @param Mat_age matrix nsim x nages of maturity at age
#' @param Wt_age matrix nsim x nages of weight at age
#' @param R0 vector nsim long of unfished recruitment
#' @param V array of vulnerability nsim x maxage x nyears 
#' @param nyears integer: number of historical years
#' @param maxage integer: maximum age
#' @param mov array of movement nsim x 2 x 2
#' @param Spat_targ vector of spatial targetting parameters
#' @param SRrel integer representing recruitmetn dynamics type 1: Bev Holt 2: Ricker
#' @param aR vector of recruitment parameters
#' @param bR vector of recruitment parameters
#' @param Bfrac fraction of SSBMSY that is the target 
#' @param ploty logical: should a plot be produced 
#' @author T. Carruthers
#' @export getBlow
getBlow<-function(x,SSBMSY,MGThorizon,Find,Perr,M_ageArray,hs,Mat_age,Wt_age,R0,V,nyears,maxage,mov,Spat_targ,SRrel,aR,bR,Bfrac=0.5,ploty=F){
  
  opt<-optimize(Blow_opt,log(c(0.0075,15)),SSBMSYc=SSBMSY[x],MGThorizon=MGThorizon[x],Fc=Find[x,],Perrc=Perr[x,],
                Mc=M_ageArray[x,,],hc=hs[x],Mac=Mat_age[x,],Wac=Wt_age[x,,],
                R0c=R0[x],Vc=V[x,,],nyears=nyears,maxage=maxage,movc=mov[x,,],
                Spat_targc=Spat_targ[x],SRrelc=SRrel[x],aRc=aR[x,],bRc=bR[x,],Bfrac,mode=1)
  
  if(ploty){Blow_opt(opt$minimum,SSBMSYc=SSBMSY[x],MGThorizon=MGThorizon[x],Fc=Find[x,],Perrc=Perr[x,],
                     Mc=M_ageArray[x,,],hc=hs[x],Mac=Mat_age[x,],Wac=Wt_age[x,,],
                     R0c=R0[x],Vc=V[x,,],nyears=nyears,maxage=maxage,movc=mov[x,,],
                     Spat_targc=Spat_targ[x],SRrelc=SRrel[x],aRc=aR[x,],bRc=bR[x,],Bfrac,mode=3)
  }
  
  Blow_opt(opt$minimum,SSBMSYc=SSBMSY[x],MGThorizon=MGThorizon[x],Fc=Find[x,],Perrc=Perr[x,],
           Mc=M_ageArray[x,,],hc=hs[x],Mac=Mat_age[x,],Wac=Wt_age[x,,],
           R0c=R0[x],Vc=V[x,,],nyears=nyears,maxage=maxage,movc=mov[x,,],
           Spat_targc=Spat_targ[x],SRrelc=SRrel[x],aRc=aR[x,],bRc=bR[x,],mode=2)
  
  
}

#' Blow internal parallel optimization function
#'
#' Find the current biomass at which it would take HZN mean generation times
#' to reach Bfrac x SSBMSY biomass level given zero catches
#'
#' @param lnq number: estimate of log catchability
#' @param SSBMSYc number: spawning biomass at MSY
#' @param MGThorizon number: MGT x HZN
#' @param Fc vector nyears long of fishing mortality rate
#' @param Perrc vector nyears+maxage-1 long of recruitment devitions 
#' @param Mc matrix maxage by nyears+proyears of natural mortality rate
#' @param hc number: steepness values
#' @param Mac vector nages long of maturity at age
#' @param Wac vector nages long  of weight at age
#' @param R0c number: unfished recruitment
#' @param Vc matrix of vulnerability maxage x nyears 
#' @param nyears integer: number of historical years
#' @param maxage integer: maximum age
#' @param movc matrix of movement 2 x 2
#' @param Spat_targc number: spatial targetting parameters
#' @param SRrelc integer representing recruitmetn dynamics type 1: Bev Holt 2: Ricker
#' @param aRc number: recruitment parameter
#' @param bRc number: recruitment parameter
#' @param Bfrac fraction of SSBMSY that is the target 
#' @param mode 1: find Blow 2:report blow  3:plot results 
#' @author T. Carruthers
#' @export Blow_opt
Blow_opt<-function(lnq,SSBMSYc,MGThorizon,Fc,Perrc,Mc,hc,Mac,Wac,R0c,Vc,nyears,maxage,movc,Spat_targc,SRrelc,aRc,bRc,Bfrac,mode=1){
  
  qc<-exp(lnq)
  nareas<-nrow(movc)
  #areasize<-c(asizec,1-asizec)
  idist<-rep(1/nareas,nareas)
  for(i in 1:300)idist<-apply(array(idist,c(2,2))*movc,2,sum)
  
  # N<-array(exp(-Mc[1]*((1:maxage)-1))*R0c,dim=c(maxage,nareas))*array(rep(idist,each=maxage),dim=c(maxage,nareas))\
  surv <- rep(1, maxage)
  surv[2:maxage] <-exp(-cumsum(Mc[,1]))[1:(maxage-1)]
  N <-array(surv * R0c,dim=c(maxage,nareas))*array(rep(idist,each=maxage),dim=c(maxage,nareas))

  SSN<-Mac*N   # Calculate initial spawning stock numbers
  Biomass<-N*Wac[,1]
  SSB<-SSN*Wac[,1]                               # Calculate spawning stock biomass
  
  B0<-sum(Biomass)
  R0a<-idist*R0c
  SSB0<-apply(SSB,2,sum)
  SSBpR<-SSB0/R0a                              # Calculate spawning stock biomass per recruit
  
  # N<-Perrc[maxage:1]*array(exp(-Mc[1]*((1:maxage)-1))*R0c,dim=c(maxage,nareas))*array(rep(idist,each=maxage),dim=c(maxage,nareas))
  N<-Perrc[maxage:1]*array(surv * R0c,dim=c(maxage,nareas))*array(rep(idist,each=maxage),dim=c(maxage,nareas))
  SSN<-Mac*N   # Calculate initial spawning stock numbers
  Biomass<-N*Wac[,1]
  SSB<-SSN*Wac[,1]                               # Calculate spawning stock biomass
  
  SSBstore<-rep(NA,nyears+MGThorizon)
  
  for(y in 1:(nyears+MGThorizon)){
   
    if(y<=nyears){ # historical catches
      targ<-(apply(Vc[,y]*Biomass,2,sum)^Spat_targc)/mean(apply(Vc[,y]*Biomass,2,sum)^Spat_targc)
      FMc<-array(qc*Fc[y]*Vc[,y],dim=c(maxage,nareas))*array(rep(targ,each=maxage),dim=c(maxage,nareas))
      Zc<-FMc+Mc[,y]# Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
    }else{ # zero catches
      targ<-(apply(Vc[,nyears]*Biomass,2,sum)^Spat_targc)/mean(apply(Vc[,nyears]*Biomass,2,sum)^Spat_targc)
      FMc<-array(0,dim=c(maxage,nareas)) # zero catch
      Zc<-FMc+Mc[,nyears]
    }
    
    N[2:maxage,]<-N[1:(maxage-1),]*exp(-Zc[1:(maxage-1),])         # Total mortality
    
    if(SRrelc==1){
      if(y<=nyears){
        N[1,]<-Perrc[y+maxage-1]*(0.8*R0a*hc*apply(SSB,2,sum))/(0.2*SSBpR*R0a*(1-hc)+(hc-0.2)*apply(SSB,2,sum))  # Recruitment assuming regional R0 and stock wide steepness
      }else{
        N[1,]<-(0.8*R0a*hc*apply(SSB,2,sum))/(0.2*SSBpR*R0a*(1-hc)+(hc-0.2)*apply(SSB,2,sum))
      }
    }else{
      N[1,]<- aRc*apply(SSB,2,sum)*exp(-bRc*apply(SSB,2,sum))
    }
    
    #print(N[1])
    indMov<-as.matrix(expand.grid(1:nareas,1:nareas,1:maxage)[3:1])
    indMov2<-indMov[,1:2]
    indMov3<-indMov[,2:3]
    temp<-array(N[indMov2]*movc[indMov3],dim=c(nareas,nareas,maxage))
    N<-apply(temp,c(3,1),sum)
    SSN<-N*Mac
    
    if(y<=nyears){
      SSB<-SSN*Wac[,y]
      Biomass<-N*Wac[,y]
      SBiomass<-SSN*Wac[,y]
    }else{
      SSB<-SSN*Wac[,nyears]
      Biomass<-N*Wac[,nyears]
      SBiomass<-SSN*Wac[,nyears]
      
    }
    
    SSBstore[y]<-sum(SSB)
    #print(sum(Biomass))
  } # end of year
  
  if(mode==1){
    pen<-0
    if(SSBstore[nyears]>(0.8*SSBMSYc))pen<-(SSBstore[nyears]-(0.8*SSBMSYc))^2 # penalty to keep inital depletion under SSB for rebuilding
    return(pen+(log(sum(SBiomass))-log(SSBMSYc*Bfrac))^2)
  }else if(mode==2){
    return(SSBstore[nyears])
  }else{
    plot(SSBstore,ylim=c(0,max(SSBstore)),col='red',xlab="Year",ylab="SSB")
    abline(h=0)
    abline(h=SSBMSYc*Bfrac,lty=2)
  }
  
}


