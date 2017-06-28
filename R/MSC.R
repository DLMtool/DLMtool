# ==================================================================
# == Code for Marine Stewardship Council - DLMtool application =====
# ==================================================================


MSCplot<-function(MSEobj,LTL=FALSE,zoom=NA) {
  par(mfrow=c(1,2),mai = c(1, 1, 0.02, 0.02))
  
  yend <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
 
  # - Status plots -------------------------------------------------
  
  if(!LTL){
    
    B50<-100*apply(MSEobj@B_BMSY[, , yend]>0.5, 2, mean, na.rm = T)
    B100<-100*apply(MSEobj@B_BMSY[, , yend]>1, 2, mean, na.rm = T)
    
    plot(B50, B100, col = "white", xlab = "", ylab = "", main = "", axes = F,xlim=c(60,100))
    MSC_status_zones()
    xs <- pretty(seq(min(B50), max(B50), length.out = 8))
    ys <- pretty(seq(min(B100), max(B100), length.out = 8))
    axis(1, xs, xs)
    axis(2, ys, ys)
    text(B50, B100, MSEobj@MPs, font = 2, cex = 0.9)
    mtext("Probability Biommass > 50% BMSY", 1, line = 2.5)
    mtext("Probability Biomass > BMSY", 2, line = 2.5)
  
  }else{
    
    SSB20<-100*apply((MSEobj@SSB[, , yend]/MSEobj@OM$SSB0)>0.2, 2, mean, na.rm = T)
    SSB75<-100*apply((MSEobj@SSB[, , yend]/MSEobj@OM$SSB0)>0.75, 2, mean, na.rm = T)
    
    plot(SSB20, SSB75, col = "white", xlab = "", ylab = "", main = "", axes = F,xlim=c(60,100),ylim=c(0,100))
    MSC_status_zones_LTL()
    xs <- pretty(seq(60, 100, length.out = 8))
    ys <- pretty(seq(0,100, length.out = 8))
    axis(1, xs, xs)
    axis(2, ys, ys)
    text(SSB20, SSB75, MSEobj@MPs, font = 2, cex = 0.9)
    mtext("Probability Spawning biomass > 20% SSB0", 1, line = 2.5)
    mtext("Probability Spawning biomass > 75% SSB0", 2, line = 2.5)
    
  }
  
  
  # - Rebuilding plots -------------------------------------------------
  
  
  
  
  
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
  
