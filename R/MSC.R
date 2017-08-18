# ==================================================================
# == Code for Marine Stewardship Council - DLMtool application =====
# ==================================================================


MSCplot<-function(MSEobj,LTL=FALSE,zoom=NA,plot=T) {
  par(mfrow=c(1,2),mai = c(0.8, 0.8, 0.4, 0.02))
  
  yend <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
  MPcex<-0.8
  # - Status plots -------------------------------------------------
  
  Status=rep("None",MSEobj@nMPs)
  
  if(!LTL){
    
    B50<-100*apply(MSEobj@B_BMSY[, , yend]>0.5, 2, mean, na.rm = T)
    B100<-100*apply(MSEobj@B_BMSY[, , yend]>1, 2, mean, na.rm = T)
    
    if(plot){
      plot(B50, B100, col = "white", xlab = "", ylab = "", axes = F,ylim=c(0,100),xlim=c(60,100), main="Status")
      MSC_status_zones()
      xs <- pretty(seq(min(B50), max(B50), length.out = 8))
      ys <- pretty(seq(min(B100), max(B100), length.out = 8))
      axis(1, xs, xs)
      axis(2, ys, ys)
      text(B50, B100, MSEobj@MPs, font = 2, cex = MPcex)
      mtext("Probability Biomass > 50% BMSY", 1, line = 2.5)
      mtext("Probability Biomass > BMSY", 2, line = 2.5)
    }
    
    Status[B50>50]<-"SG60"
    Status[B50>50&B100>50]<-"SG80"
    Status[B50>95&B100>50]<-"SG100"
    
  }else{
    
    SSB20<-100*apply((MSEobj@SSB[, , yend]/MSEobj@OM$SSB0)>0.2, 2, mean, na.rm = T)
    SSB75<-100*apply((MSEobj@SSB[, , yend]/MSEobj@OM$SSB0)>0.75, 2, mean, na.rm = T)
    
    if(plot){
      plot(SSB20, SSB75, col = "white", xlab = "", ylab = "",  axes = F,xlim=c(60,100),ylim=c(0,100), main="Status")
      MSC_status_zones_LTL()
      xs <- pretty(seq(60, 100, length.out = 8))
      ys <- pretty(seq(0,100, length.out = 8))
      axis(1, xs, xs)
      axis(2, ys, ys)
      text(SSB20, SSB75, MSEobj@MPs, font = 2, cex = MPcex)
      mtext("Probability Spawning biomass > 20% SSB0", 1, line = 2.5)
      mtext("Probability Spawning biomass > 75% SSB0", 2, line = 2.5)
    }
    
    Status[SSB20>70]<-"SG60"
    Status[SSB20>80 & SSB75>50]<-"SG80"
    Status[SSB20>95 & SSB75>95]<-"SG100"
    
  }
  
  # - Rebuilding plots -------------------------------------------------
  
  Rebuilding=rep("SG60",MSEobj@nMPs)
  
  na<-dim(MSEobj@CAA)[3]
  Z<-array(MSEobj@OM$M,c(MSEobj@nsim,na))
  MGTsurv<-t(exp(-apply(Z,1,cumsum)))
  agearray<-array(rep(1:na,each=MSEobj@nsim),c(MSEobj@nsim,na))
  ageMsd<-0.15
  Mat_age <- 1/(1 + exp((MSEobj@OM$ageM - agearray)/(MSEobj@OM$ageM * ageMsd)))
  MGT<-apply(agearray*(Mat_age*MGTsurv),1,sum)/apply(Mat_age*MGTsurv,1,sum)
  HZN<-ceiling(MGT*2)
  HZN[HZN>20]<-20
  HZN[HZN<5]<-5
  
  if(sum(HZN>MSEobj@proyears)>0){
    warning("At least one simulation needs a greater number of projected years to reach the rebuilding time horizon")
    print(paste("Projected number of years:",MSEobj@proyears))
    print("Required time horizon for projection for each simulation:")
    print(HZN)
    cond<-HZN>MSEobj@proyears
    print(paste(sum(cond),"simulations had rebuilding time artificially set to the last projection year"))
    HZN[cond]<-MSEobj@proyears
  }
  
  HZN1<-ceiling(MGT)
  HZN1<-max(HZN1,2)
  HZN1<-min(HZN1,10)
  
  RB1<-100*apply(MSEobj@B_BMSY[, , HZN1]>1, 2, mean, na.rm = T)
  RB2<-100*apply(MSEobj@B_BMSY[, , HZN]>1, 2, mean, na.rm = T)
  
  if(plot){
    plot(RB1, RB2, col = "white", xlab = "", ylab = "", axes = F,xlim=c(0,100),ylim=c(20,100), main="Rebuilding")
    MSC_rebuilding_zones()
    
    xs <- pretty(seq(0,100, length.out = 8))
    ys <- pretty(seq(0,100, length.out = 8))
    axis(1, xs, xs)
    axis(2, ys, ys)
    text(RB1, RB2, MSEobj@MPs, font = 2, cex = MPcex)
    mtext("Probability Biommass > BMSY (2 MGT)", 1, line = 2.5)
    mtext("Probability Biomass > BMSY (1 MGT)", 2, line = 2.5)
  }
  Rebuilding[RB2>70]<-"SG80"
  Rebuilding[RB1>95]<-"SG100"
 
  as.data.frame(cbind(MSEobj@MPs,Status,Rebuilding))
  
}

MSC_status_zones<-function(){
 
  SG60<-'chartreuse'
  SG80<-'green2'
  SG100<-'green3'
  
  textcol='white'
  textcex<-0.9
  
  polygon(c(70,80,80,70),c(0,0,100,100),col=SG60,border=SG60)
  polygon(c(80,100,100,80),c(0,0,50,50),col=SG60,border=SG60)
  polygon(c(80,95,95,80),c(50,50,100,100),col=SG80,border=SG80)
  polygon(c(95,100,100,95),c(50,50,100,100),col=SG100,border=SG100)
  
  
  text(97.5,75,"SG100",col=textcol,font=2,srt=90,cex=textcex)
  text(87.5,75,"SG80",col=textcol,font=2,srt=90,cex=textcex)
  text(85,35,"SG60",col=textcol,font=2,cex=textcex)
 
}  
 
MSC_status_zones_LTL<-function(){
  
  SG60<-'chartreuse'
  SG80<-'green2'
  SG100<-'green3'
  
  textcol='white'
  textcex<-0.9
  
  polygon(c(70,80,80,70),c(0,0,100,100),col=SG60,border=SG60)
  polygon(c(80,100,100,80),c(0,0,50,50),col=SG60,border=SG60)
  polygon(c(80,95,95,80),c(50,50,100,100),col=SG80,border=SG80)
  polygon(c(95,100,100,95),c(50,50,95,95),col=SG80,border=SG80)
  polygon(c(95,100,100,95),c(95,95,100,100),col=SG100,border=SG100)
  
  
  text(97.5,97.5,"SG100",col=textcol,font=2,cex=textcex*0.6)
  text(90,72.5,"SG80",col=textcol,font=2,cex=textcex)
  text(85,35,"SG60",col=textcol,font=2,cex=textcex)
  
}   

  
MSC_rebuilding_zones<-function(){
  
  SG60<-'chartreuse'
  SG80<-'green2'
  SG100<-'green3'
  
  textcol='white'
  textcex<-0.9
  
  polygon(c(0,100,100,0),c(0,0,100,100),col=SG60,border=SG60)
  polygon(c(70,100,100,70),c(0,0,100,100),col=SG80,border=SG80)
  polygon(c(0,100,100,0),c(95,95,100,100),col=SG100,border=SG100)
  
  
  text(50,97.5,"SG100",col=textcol,font=2,cex=textcex*0.6)
  text(85,55,"SG80",col=textcol,font=2,cex=textcex)
  text(35,55,"SG60",col=textcol,font=2,cex=textcex)
  
} 