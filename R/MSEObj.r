
# A collection of functions which operate on MSE object (class MSE) which is 
# produced by the DLMtool MSE.
# These functions are used to summarize, plot, and manipulate the results of the
# MSE. 
# Most of these functions have accompanying help files (?function.name), 
# although, for purpose of clarity, some supporting internal functions have also
# been included here.

# January 2016
# Tom Carruthers UBC (t.carruthers@fisheries.ubc.ca)
# Adrian Hordyk (a.hordyk@murdoch.edu.au)

# Have I undertaken enough simulations (nsim)? 
# Has my MSE converged on stable (reliable) peformance metrics?
CheckConverg <- function(MSEobj, thresh=2, Plot=TRUE) {
  nm<-MSEobj@nMPs
  nsim<-MSEobj@nsim
  proyears<-MSEobj@proyears
  
  Yd <- CumlYd <- array(NA,c(nm,nsim))
  P10 <- CumlP10 <- array(NA,c(nm,nsim))
  P50 <- CumlP50 <- array(NA,c(nm,nsim))
  P100 <- CumlP100 <- array(NA,c(nm,nsim))
  POF <- CumlPOF <- array(NA,c(nm,nsim))
  yind <-max(MSEobj@proyears-4,1):MSEobj@proyears
  RefYd<-MSEobj@OM$RefY
  
  for(m in 1:nm){
    Yd[m,] <-round(apply(MSEobj@C[,m,yind],1,mean,na.rm=T)/RefYd*100,1)
    POF[m,] <-round(apply(MSEobj@F_FMSY[,m,]>1,1,sum,na.rm=T)/proyears*100,1)
    P10[m,] <-round(apply(MSEobj@B_BMSY[,m,]<0.1,1,sum,na.rm=T)/proyears*100,1)
    P50[m,] <-round(apply(MSEobj@B_BMSY[,m,]<0.5,1,sum,na.rm=T)/proyears*100,1)
    P100[m,] <-round(apply(MSEobj@B_BMSY[,m,]<1,1,sum,na.rm=T)/proyears*100,1)
    CumlYd[m,] <- cumsum(Yd[m,]) / seq_along(Yd[m,])#/ mean(Yd[m,], na.rm=TRUE) 
    CumlPOF[m,] <- cumsum(POF[m,]) / seq_along(POF[m,])# / mean(POF[m,], na.rm=TRUE)
    CumlP10[m,] <- cumsum(P10[m,]) / seq_along(P10[m,])# / mean(P10[m,], na.rm=TRUE)
    CumlP50[m,] <- cumsum(P50[m,]) / seq_along(P50[m,])# / mean(P50[m,], na.rm=TRUE)
    CumlP100[m,] <- cumsum(P100[m,]) / seq_along(P100[m,])# / mean(P100[m,], na.rm=TRUE)
  }
  
  # CumlYd[is.nan(CumlYd)] <- 1
  # CumlPOF[is.nan(CumlPOF)] <- 1
  # CumlP10[is.nan(CumlP10)] <- 1
  # CumlP50[is.nan(CumlP50)] <- 1
  # CumlP100[is.nan(CumlP100)] <- 1
  if (Plot) {
    par(mfrow=c(2,3), cex.axis=1.5, cex.lab=1.7, oma=c(1,1,0,0), mar=c(5,5,1,1), bty="l")
    matplot(t(CumlYd), type="l", xlab="Iteration", ylab="Rel. Yield")
    matplot(t(CumlPOF), type="l", xlab="Iteration", ylab="Prob. F/FMSY > 1")
    matplot(t(CumlP10), type="l", xlab="Iteration", ylab="Prob. B/BMSY < 0.1")
    matplot(t(CumlP50), type="l", xlab="Iteration", ylab="Prob. B/BMSY < 0.5")
    matplot(t(CumlP100), type="l", xlab="Iteration", ylab="Prob. B/BMSY < 1")
  }
  
  Chk <- function(X) {
    # checks if difference in
    # last 10 iterations is greater than thresh
    L <- length(X)
    Y <- 1: min(nsim, 10)
    
    # return(mean(abs((X[L-1:10] - X[L]))/X[L], na.rm=TRUE) > thresh)
    return(mean(abs((X[L-Y] - X[L])), na.rm=TRUE) > thresh)
  }
  
  NonCon <- sort(unique(c(which(apply(CumlYd, 1, Chk)), which(apply(CumlPOF, 1, Chk)),
                          which(apply(CumlP10, 1, Chk)), which(apply(CumlP50, 1, Chk)), 	
                          which(apply(CumlP100, 1, Chk)))))
  
  if (length(NonCon) == 1) NonCon <- rep(NonCon, 2)
  if (length(NonCon) > 0) {
    if (Plot) {
      plot(c(0,1), c(0,1), type="n", axes=FALSE, xlab="", ylab="")
      text(0.5, 0.5, "Some MPs have not converged", cex=1)
      #ev.new()
      par(mfrow=c(2,3), cex.axis=1.5, cex.lab=1.7, oma=c(1,1,0,0), mar=c(5,5,1,1), bty="l")
      matplot(t(CumlYd[NonCon,]), type="b", xlab="Iteration", ylab="Rel. Yield", lwd=2)
      matplot(t(CumlPOF[NonCon,]), type="b", xlab="Iteration", ylab="Prob. F/FMSY > 1", lwd=2)
      matplot(t(CumlP10[NonCon,]), type="b", xlab="Iteration", ylab="Prob. B/BMSY < 0.1", lwd=2)
      matplot(t(CumlP50[NonCon,]), type="b", xlab="Iteration", ylab="Prob. B/BMSY < 0.5", lwd=2)
      matplot(t(CumlP100[NonCon,]), type="b", xlab="Iteration", ylab="Prob. B/BMSY < 1", lwd=2)
      legend(nsim*1.25, 50, legend=MSEobj@MPs[NonCon], col=1:length(NonCon), bty="n", 
             xpd=NA, lty=1, lwd=2, cex=1.25)
    }  
    
    message("Some MPs may not have converged in ", nsim, " iterations (threshold = ", 
            thresh, "%)")
    message("MPs are: ", paste(MSEobj@MPs[NonCon], " "))
    message("MPs #: ", paste(NonCon, " "))
    return(data.frame(Num=NonCon, MP=MSEobj@MPs[NonCon]))
  }
  if (length(NonCon) == 0) {
    if (Plot) {
      plot(c(0,1), c(0,1), type="n", axes=FALSE, xlab="", ylab="")
      text(0.5, 0.5, "All MPs converged", cex=1)
    }	
    message("All MPs appear to have converged in ", nsim, " iterations (threshold = ", 
            thresh, "%)")
  }

}

# Trade-Off Plots --------------------------------------------------------------
Tplot<-function(MSEobj,nam=NA){
  FMSYr<-quantile(MSEobj@F_FMSY,c(0.001,0.90),na.rm=T)
  BMSYr<-quantile(MSEobj@B_BMSY,c(0.001,0.975),na.rm=T)

  colsse<-rainbow(100,start=0,end=0.36)[1:100]
  colB<-rep(colsse[100],ceiling(BMSYr[2]*100))
  colB[1:100]<-colsse
  colB<-makeTransparent(colB,60)
  colsse<-rainbow(200,start=0,end=0.36)[200:1]
  colF<-rep(colsse[200],ceiling(FMSYr[2]*100))
  colF[1:200]<-colsse
  colF<-makeTransparent(colF,60)

  Yd<-rep(NA,MSEobj@nMPs)
  P10<-rep(NA,MSEobj@nMPs)
  P50<-rep(NA,MSEobj@nMPs)
  P100<-rep(NA,MSEobj@nMPs)
  POF<-rep(NA,MSEobj@nMPs)
  yind<-max(MSEobj@proyears-4,1):MSEobj@proyears
  RefYd<-MSEobj@OM$RefY

  for(mm in 1:MSEobj@nMPs){
    Yd[mm]<-round(mean(apply(MSEobj@C[,mm,yind],1,mean,na.rm=T)/RefYd,na.rm=T)*100,1)
  #cbind(MSEobj@C[,mm,yind],unlist(MSEobj@OM$MSY))
    POF[mm]<-round(sum(MSEobj@F_FMSY[,mm,]>1,na.rm=T)/prod(dim(MSEobj@F_FMSY[,mm,]),na.rm=T)*100,1)
    P10[mm]<-round(sum(MSEobj@B_BMSY[,mm,]<0.1,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
    P50[mm]<-round(sum(MSEobj@B_BMSY[,mm,]<0.5,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
    P100[mm]<-round(sum(MSEobj@B_BMSY[,mm,]<1,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
  }
  
  #dev.new2(width=7,height=7)
  par(mfrow=c(2,2),mai=c(0.9,1,0.1,0.1),omi=c(0.1,0.1,0.4,0))

  tradeoffplot(POF,Yd,"Prob. of overfishing (%)", "Relative yield",MSEobj@MPs[1:MSEobj@nMPs],vl=50,hl=100)
  tradeoffplot(P100,Yd,"Prob. biomass < BMSY (%)", "Relative yield",MSEobj@MPs[1:MSEobj@nMPs],vl=50,hl=100)
  tradeoffplot(P50,Yd,"Prob. biomass < 0.5BMSY (%)", "Relative yield",MSEobj@MPs[1:MSEobj@nMPs],vl=50,hl=100)
  tradeoffplot(P10,Yd,"Prob. biomass < 0.1BMSY (%)", "Relative yield",MSEobj@MPs[1:MSEobj@nMPs],vl=50,hl=100)

  if(is.na(nam))mtext(deparse(substitute(MSEobj)),3,outer=T,line=0.3,font=2)
  if(!is.na(nam) & !is.character(nam))mtext(MSEobj@Name,3,outer=T,line=0.3,font=2)
  if(!is.na(nam) & is.character(nam))mtext(nam,3,outer=T,line=0.3,font=2)
}

Tplot2<-function(MSEobj,nam=NA){
  LTY<-rep(NA,MSEobj@nMPs)
  STY<-rep(NA,MSEobj@nMPs)
  VY<-rep(NA,MSEobj@nMPs)
  B10<-rep(NA,MSEobj@nMPs)
  yend<-max(MSEobj@proyears-4,1):MSEobj@proyears
  ystart<-1:5
  RefYd<-MSEobj@OM$RefY
  y1<-1:(MSEobj@proyears-1)
  y2<-2:MSEobj@proyears
  for(mm in 1:MSEobj@nMPs){
    LTY[mm]<-round(sum(MSEobj@C[,mm,yend]/RefYd>0.5,na.rm=T)/(MSEobj@nsim*length(yend)),3)*100
    STY[mm]<-round(sum(MSEobj@C[,mm,ystart]/RefYd>0.5,na.rm=T)/(MSEobj@nsim*length(ystart)),3)*100
    AAVY<-apply(((MSEobj@C[,mm,y1]-MSEobj@C[,mm,y2])^2)^0.5,1,mean,na.rm=T)/apply(MSEobj@C[,mm,y2],1,mean,na.rm=T)
    VY[mm]<-round(sum(AAVY<0.1,na.rm=T)/MSEobj@nsim,3)*100
    B10[mm]<-round(sum(MSEobj@B_BMSY[,mm,]>0.1,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,])),3)*100
  }
  par(mfrow=c(1,2),mai=c(1.5,1.5,0.1,0.1),omi=c(0.1,0.1,0.4,0))
  tradeoffplot(STY,LTY,"P(Short term yield over half FMSY)","P(Long term yield over half FMSY)",MSEobj@MPs[1:MSEobj@nMPs],vl=1,hl=1)
  tradeoffplot(B10,VY,"P(Biomass > 0.1 BMSY)", "P(CV in yield less than 0.1)",MSEobj@MPs[1:MSEobj@nMPs],vl=1,hl=1)
  if(is.na(nam))mtext(deparse(substitute(MSEobj)),3,outer=T,line=0.3,font=2)
  if(!is.na(nam))mtext(MSEobj@Name,3,outer=T,line=0.3,font=2)
}

NOAA_plot<-function(MSEobj,nam=NA,type=NA,panel=T){
  
  Yd<-rep(NA,MSEobj@nMPs)
  B50<-rep(NA,MSEobj@nMPs)
  PNOF<-rep(NA,MSEobj@nMPs)
  LTY<-rep(NA,MSEobj@nMPs)
  STY<-rep(NA,MSEobj@nMPs)
  VY<-rep(NA,MSEobj@nMPs)
  
  y1<-1:(MSEobj@proyears-1)
  y2<-2:MSEobj@proyears
  
  yend<-max(MSEobj@proyears-4,1):MSEobj@proyears
  
  RefYd<-MSEobj@OM$RefY
  
  for(mm in 1:MSEobj@nMPs){
    
    PNOF[mm]<-round(sum(MSEobj@F_FMSY[,mm,]<1,na.rm=T)/prod(dim(MSEobj@F_FMSY[,mm,]),na.rm=T)*100,1)
    B50[mm]<-round(sum(MSEobj@B_BMSY[,mm,]>0.5,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
    LTY[mm]<-round(sum(MSEobj@C[,mm,yend]/RefYd>0.5,na.rm=T)/(MSEobj@nsim*length(yend)),3)*100
    AAVY<-apply((((MSEobj@C[,mm,y1]-MSEobj@C[,mm,y2])/MSEobj@C[,mm,y2])^2)^0.5,1,mean,na.rm=T) 
    VY[mm]<-round(sum(AAVY<0.15,na.rm=T)/MSEobj@nsim,3)*100

  }
  
  #dev.new2(width=7,height=7)
  if(panel)par(mfrow=c(1,2),mai=c(1.5,1.5,0.1,0.1),omi=c(0.1,0.1,0.4,0))
  
  if(is.na(type)){
    tradeoffplot(PNOF,LTY,"Prob. of not overfishing (%)", "Long-term yield ",MSEobj@MPs[1:MSEobj@nMPs],vl=50,hl=100)
    tradeoffplot(B50,VY,"Prob. biomass above half BMSY (%)", "Prob. AAVY less than 15%",MSEobj@MPs[1:MSEobj@nMPs],vl=80,hl=50)
  }else{
    tradeoffplot3(PNOF,LTY,"Prob. of not overfishing (%)", "Long-term yield",MSEobj@MPs[1:MSEobj@nMPs],vl=50,hl=100,xlim=c(45,105),ylim=c(0,105))
    tradeoffplot3(B50,VY,"Prob. biomass above half BMSY (%)", "Prob. AAVY less than 15%",MSEobj@MPs[1:MSEobj@nMPs],vl=80,hl=50,xlim=c(75,105),ylim=c(45,105))
  }
  
  #if(is.na(nam))mtext(deparse(substitute(MSEobj)),3,outer=T,line=0.3,font=2)
  #if(!is.na(nam) & !is.character(nam))mtext(MSEobj@Name,3,outer=T,line=0.3,font=2)
  #if(!is.na(nam) & is.character(nam))mtext(nam,3,outer=T,line=0.3,font=2)
  
  
  temp<-data.frame(PNOF,B50,LTY,VY)
  row.names(temp)<-MSEobj@MPs[1:MSEobj@nMPs]
  temp
  
}

# Projection plot - plots F/FMSY and B/BMSY for projection years and all MPs in 
# MSEobj.  
Pplot<-function(MSEobj,nam=NA){
  
  FMSYr<-quantile(MSEobj@F_FMSY,c(0.001,0.90),na.rm=T)
  BMSYr<-quantile(MSEobj@B_BMSY,c(0.001,0.975),na.rm=T)
  
  colsse<-rainbow(100,start=0,end=0.36)[1:100]
  colB<-rep(colsse[100],ceiling(BMSYr[2]*100))
  colB[1:100]<-colsse
  colB<-makeTransparent(colB,60)
  colsse<-rainbow(200,start=0,end=0.36)[200:1]
  colF<-rep(colsse[200],ceiling(FMSYr[2]*100))
  colF[1:200]<-colsse
  colF<-makeTransparent(colF,60)
  
  Yd<-rep(NA,MSEobj@nMPs)
  P10<-rep(NA,MSEobj@nMPs)
  P50<-rep(NA,MSEobj@nMPs)
  P100<-rep(NA,MSEobj@nMPs)
  POF<-rep(NA,MSEobj@nMPs)
  yind<-max(MSEobj@proyears-4,1):MSEobj@proyears
  RefYd<-MSEobj@OM$RefY
  
  for(mm in 1:MSEobj@nMPs){
    Yd[mm]<-round(mean(apply(MSEobj@C[,mm,yind],1,mean,na.rm=T)/RefYd,na.rm=T)*100,1)
    #cbind(MSEobj@C[,mm,yind],unlist(MSEobj@OM$MSY))
    POF[mm]<-round(sum(MSEobj@F_FMSY[,mm,]>1,na.rm=T)/prod(dim(MSEobj@F_FMSY[,mm,]),na.rm=T)*100,1)
    P10[mm]<-round(sum(MSEobj@B_BMSY[,mm,]<0.1,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
    P50[mm]<-round(sum(MSEobj@B_BMSY[,mm,]<0.5,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
    P100[mm]<-round(sum(MSEobj@B_BMSY[,mm,]<1,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
  }
    
  nr<-ceiling(MSEobj@nMPs/8)
  nc<-ceiling(MSEobj@nMPs/nr)
  nr<-nr*2
  MSEcols<-c('red','green','blue','orange','brown','purple','dark grey','violet','dark red','pink','dark blue','grey')
  temp<-array(0,c(nr*2+(nr/2-1),nc*2))
  i<-0
  for(c in 1:nc){
    for(r in 1:nr){
      i<-i+1
      temp[(ceiling(r/2)-1)+(1:2)+(r-1)*2,(1:2)+(c-1)*2]<-((c-1)*nr)+r
    }
  }
  par(mfcol=c(nr,nc),mai=c(0.2,0.35,0.3,0.01),omi=c(0.5,0.4,0.4,0.05))
  layout(temp)
  #dev.new2(width=nc*3,height=nr*3)
  #
  lwdy<-2.5

  for(mm in 1:MSEobj@nMPs){
    plot(MSEobj@F_FMSY[1,mm,],ylim=FMSYr,col=colF[ceiling(mean(MSEobj@F_FMSY[1,mm,],na.rm=T)*100)],type='l',lwd=lwdy)
    for(i in 1:MSEobj@nsim)lines(MSEobj@F_FMSY[i,mm,],col=colF[ceiling(mean(MSEobj@F_FMSY[i,mm,],na.rm=T)*100)],lwd=lwdy)
    abline(h=100,col="grey",lwd=3)
    mtext(MSEobj@MPs[mm],3,outer=F,line=0.6)
    legend('topright',c(paste(POF[mm],"% POF",sep=""),
                 paste(Yd[mm],"% FMSY yield",sep="")),bty='n',cex=0.8)
    if(mm%in%(1:(nr/2)))mtext("F/FMSY",2,line=2.5,outer=F)
    abline(h=1,col=makeTransparent("grey",30),lwd=2.5)
  
    plot(MSEobj@B_BMSY[1,mm,],ylim=BMSYr,col=colB[ceiling(MSEobj@B_BMSY[1,mm,MSEobj@proyears]*100)],type='l',lwd=lwdy)
    for(i in 1:MSEobj@nsim)lines(MSEobj@B_BMSY[i,mm,],col=colB[ceiling(MSEobj@B_BMSY[i,mm,MSEobj@proyears]*100)],lwd=lwdy)
    abline(h=100,col="grey",lwd=3)
    legend('topright',c(paste(P100[mm],"% < BMSY",sep=""),
                 paste(P50[mm],"% < 0.5BMSY",sep=""),
                 paste(P10[mm],"% < 0.1BMSY",sep="")),bty='n',cex=0.8)
    if(mm%in%(1:(nr/2)))mtext("B/BMSY",2,line=2.5,outer=F)
    abline(h=1,col=makeTransparent("grey",30),lwd=2.5)
  
  }
  mtext("Projection year",1,outer=T,line=1.2)
  if(is.na(nam))mtext(deparse(substitute(MSEobj)),3,outer=T,line=0.3,font=2)
  if(!is.na(nam))mtext(MSEobj@Name,3,outer=T,line=0.3,font=2)
}

# Kobe plot
Kplot<-function(MSEobj,maxsim=60,nam=NA){
  nr<-floor((MSEobj@nMPs)^0.5)
  nc<-ceiling((MSEobj@nMPs)/nr)
    
  FMSYr<-quantile(MSEobj@F_FMSY,c(0.001,0.90),na.rm=T)
  BMSYr<-quantile(MSEobj@B_BMSY,c(0.001,0.975),na.rm=T)
    
  #dev.new2(width=nc*3,height=nr*3.6)
  par(mfrow=c(nr,nc),mai=c(0.45,0.45,0.45,0.01),omi=c(0.45,0.3,0.35,0.01))
  
  colsse<-rainbow(MSEobj@proyears,start=0.63,end=0.95)[1:MSEobj@proyears]
  colsse<-makeTransparent(colsse,95)
  
  for(mm in 1:MSEobj@nMPs){
    plot(c(MSEobj@B_BMSY[1,mm,1],MSEobj@B_BMSY[1,mm,2]),
         c(MSEobj@F_FMSY[1,mm,1],MSEobj@F_FMSY[1,mm,2]),xlim=BMSYr,ylim=FMSYr,
         col=colsse[1],type='l')
    
    OO<-round(sum(MSEobj@B_BMSY[,mm,MSEobj@proyears]<1&MSEobj@F_FMSY[,mm,MSEobj@proyears]>1,na.rm=T)/MSEobj@nsim*100,1)
    OU<-round(sum(MSEobj@B_BMSY[,mm,MSEobj@proyears]>1&MSEobj@F_FMSY[,mm,MSEobj@proyears]>1,na.rm=T)/MSEobj@nsim*100,1)
    UO<-round(sum(MSEobj@B_BMSY[,mm,MSEobj@proyears]<1&MSEobj@F_FMSY[,mm,MSEobj@proyears]<1,na.rm=T)/MSEobj@nsim*100,1)
    UU<-round(sum(MSEobj@B_BMSY[,mm,MSEobj@proyears]>1&MSEobj@F_FMSY[,mm,MSEobj@proyears]<1,na.rm=T)/MSEobj@nsim*100,1)
    
    #alp<-80
    #polygon(c(1,-1000,-1000,1),c(1,1,1000,1000),col=makeTransparent("orange",alp),border=makeTransparent("orange",alp))
    #polygon(c(1,1000,1000,1),c(1,1,1000,1000),col=makeTransparent("yellow",alp),border=makeTransparent("yellow",alp))
    #polygon(c(1,-1000,-1000,1),c(1,1,-1000,-1000),col=makeTransparent("yellow",alp),border=makeTransparent("yellow",alp))
    #polygon(c(1,1000,1000,1),c(1,1,-1000,-1000),col=makeTransparent("green",alp),border=makeTransparent("yellow",alp))
    
    
    abline(h=1,col="grey",lwd=3)
    abline(v=1,col="grey",lwd=3)
    #abline(v=c(0.1,0.5),col="grey",lwd=2)
    rng<-1:min(maxsim,MSEobj@nsim)
    for(i in rng){
      for(y in 1:(MSEobj@proyears-1)){
        lines(c(MSEobj@B_BMSY[i,mm,y],MSEobj@B_BMSY[i,mm,y+1]),
              c(MSEobj@F_FMSY[i,mm,y],MSEobj@F_FMSY[i,mm,y+1]),
              col=colsse[y],lwd=1.6)
      }
    }
    
    points(MSEobj@B_BMSY[rng,mm,1],MSEobj@F_FMSY[rng,mm,1],pch=19,cex=0.8,col=colsse[1])
    points(MSEobj@B_BMSY[rng,mm,MSEobj@proyears],MSEobj@F_FMSY[rng,mm,MSEobj@proyears],pch=19,cex=0.8,col=colsse[MSEobj@proyears])
    
    if(mm==1)legend('right',c("Start","End"),bty='n',text.col=c(colsse[1],colsse[MSEobj@proyears]),pch=19,col=c(colsse[1],colsse[MSEobj@proyears]))
    legend('topleft',paste(OO,"%",sep=""),bty='n',text.font=2)
    legend('topright',paste(OU,"%",sep=""),bty='n',text.font=2)
    legend('bottomleft',paste(UO,"%",sep=""),bty='n',text.font=2)
    legend('bottomright',paste(UU,"%",sep=""),bty='n',text.font=2)
    
    mtext(MSEobj@MPs[mm],3,line=0.45)
  }
  mtext("B/BMSY",1,outer=T,line=1.4)
  mtext("F/FMSY",2,outer=T,line=0.2)
  if(is.na(nam))mtext(deparse(substitute(MSEobj)),3,outer=T,line=0.25,font=2)
  if(!is.na(nam))mtext(MSEobj@Name,3,outer=T,line=0.25,font=2)
}

# A simulation by simulation approach to plotting results
comp<-function(MSEobj,MPs=NA){
  
  if(is.na(MPs))MPs<-MSEobj@MPs
  notm<-MPs[!(MPs%in%MSEobj@MPs)]
  canm<-MPs[MPs%in%MSEobj@MPs]
  if(length(notm)>0)print(paste("Methods",paste(notm,collapse=", "),"were not carried out in MSE",deparse(substitute(MSEobj)),sep=" "))
  
  if(length(canm)==0)stop(paste('None of the methods you specified were carried out in the MSE', deparse(substitute(MSEobj)),sep=""))
  
  if(length(canm)>4){
    print(paste('A maximum of four methods can be compared at once. Plotting first four:',paste(canm[1:4],collapse=", "),sep=" "))
    canm<-canm[1:4]
  } 
  
  mind<-match(canm,MSEobj@MPs)
  nm<-length(mind)
  nsim<-MSEobj@nsim
  proyears<-MSEobj@proyears
  
  Yd<-array(NA,c(nm,nsim))
  P10<-array(NA,c(nm,nsim))
  P50<-array(NA,c(nm,nsim))
  P100<-array(NA,c(nm,nsim))
  POF<-array(NA,c(nm,nsim))
  yind<-max(MSEobj@proyears-4,1):MSEobj@proyears
  RefYd<-MSEobj@OM$RefY
  
  for(m in 1:nm){
    mm<-mind[m]
    Yd[m,]<-round(apply(MSEobj@C[,mm,yind],1,mean,na.rm=T)/RefYd*100,1)
    POF[m,]<-round(apply(MSEobj@F_FMSY[,mm,]>1,1,sum,na.rm=T)/proyears*100,1)
    P10[m,]<-round(apply(MSEobj@B_BMSY[,mm,]<0.1,1,sum,na.rm=T)/proyears*100,1)
    P50[m,]<-round(apply(MSEobj@B_BMSY[,mm,]<0.5,1,sum,na.rm=T)/proyears*100,1)
    P100[m,]<-round(apply(MSEobj@B_BMSY[,mm,]<1,1,sum,na.rm=T)/proyears*100,1)
  }
  
  MSEcols<-c('red','green','blue','orange')
  
  #dev.new2(width=7,height=7)
  par(mfrow=c(2,2),mai=c(0.85,0.7,0.1,0.1),omi=rep(0.01,4))
  
  tradeoffplot2(POF,Yd,"Prob. of overfishing (%)", "Relative yield",vl=50,hl=100,coly=MSEcols,leg=NA)
  tradeoffplot2(P100,Yd,"Prob. biomass < BMSY (%)", "Relative yield",vl=50,hl=100,coly=MSEcols,leg=canm)
  tradeoffplot2(P50,Yd,"Prob. biomass < 0.5BMSY (%)", "Relative yield",vl=50,hl=100,coly=MSEcols,leg=NA)
  tradeoffplot2(P10,Yd,"Prob. biomass < 0.1BMSY (%)", "Relative yield",vl=50,hl=100,coly=MSEcols,leg=NA)
  
}


# Supporting Functions for Plots
tradeoffplot<-function(x,y,xlab,ylab,labs,cex,vl,hl){
   adjj<-c(0.7,1.3)
   XLim <- c(min(c(-10, min(x,na.rm=T)*adjj)), max(c(max(x,na.rm=T)*adjj, 110)))
   YLim <- c(min(c(-10, min(y,na.rm=T)*adjj)), max(c(max(y,na.rm=T)*adjj, 110)))
   coly<-rep(c('#0000ff95','#ff000095','#20ff1095'),50)[1:length(labs)]
   coly[labs%in%c("AvC","curE","FMSYref")]<-'black'
   # plot(NA,xlim=range(x,na.rm=T)*adjj,ylim=range(y,na.rm=T)*adjj,xlab=xlab,ylab=ylab)
   plot(NA,xlim=XLim,ylim=YLim,xlab=xlab,ylab=ylab)
   abline(v=vl,col="#99999940",lwd=2)
   abline(h=hl,col="#99999940",lwd=2)
   text(x,y,labs,font=2,col=coly,cex=0.9)
}

tradeoffplot2<-function(x,y,xlab,ylab,cex=1,vl,hl,coly,leg){
  adjj<-c(0.7,1.3)
  plot(NA,xlim=range(x,na.rm=T)*adjj,ylim=range(y,na.rm=T)*adjj,xlab=xlab,ylab=ylab)
  abline(v=vl,col="grey",lwd=2)
  abline(h=hl,col="grey",lwd=2)
  for(m in 1:nrow(x))points(x[m,],y[m,],col=makeTransparent(coly[m],50),pch=19,cex=cex)
  if(!is.na(leg[1]))legend('topright',legend=leg,text.col=coly,bty='n')
}

tradeoffplot3<-function(x,y,xlab,ylab,labs,cex,vl,hl,xlim,ylim){
  coly<-rep(c('#0000ff95','#ff000095','#20ff1095'),10)
  plot(NA,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
  abline(v=vl,col="#99999940",lwd=2)
  abline(h=hl,col="#99999940",lwd=2)
  text(x,y,labs,font=2,col=coly,cex=1)
}

makeTransparent<-function(someColor, alpha=100){
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

# Value of information analysis ------------------------------------------------
# Value of information
VOI<-function(MSEobj,ncomp=6,nbins=8,maxrow=8,Ut=NA,Utnam="Utility"){
  objnam<-deparse(substitute(MSEobj))
  nsim<-MSEobj@nsim
 
  if(is.na(Ut[1])){
    Ut<-array(NA,c(nsim,MSEobj@nMPs))
    yind<-max(MSEobj@proyears-4,1):MSEobj@proyears
    RefYd<-MSEobj@OM$RefY
  
    for(mm in 1:MSEobj@nMPs){
      Ut[,mm]<-apply(MSEobj@C[,mm,yind],1,mean,na.rm=T)/RefYd*100
    #POF[,mm]<-apply(MSEobj@F_FMSY[,mm,]>1,1,sum)/MSEobj@proyears
    #P10[,mm]<-apply(MSEobj@B_BMSY[,mm,]<0.1,1,sum)/MSEobj@proyears
    }
    Utnam<-"Long-term yield relative to MSY (%)"
  }

  MPs<-MSEobj@MPs
  nMPs<-MSEobj@nMPs
 
  onlycor<-c("RefY","A","MSY","Linf","t0","OFLreal","Spat_targ")
  vargood<-(apply(MSEobj@OM,2,sd)/(apply(MSEobj@OM,2,mean)^2)^0.5)>0.005
  # MSEobj@OM<- MSEobj@OM[,(!names(MSEobj@OM)%in%onlycor)&vargood]
  MSEobj@OM[,which((!names(MSEobj@OM)%in%onlycor)&vargood)]  
  OMp<-apply(MSEobj@OM,2,quantile,p=seq(0,1,length.out=nbins+1))
  Obsp<-apply(MSEobj@Obs,2,quantile,p=seq(0,1,length.out=nbins+1))
  OMv<-array(NA,c(nMPs,ncol(MSEobj@OM),nbins))
  Obsv<-array(NA,c(nMPs,ncol(MSEobj@Obs),nbins))
  
  for(mm in 1:nMPs){
    for(j in 1:nbins){
      for(i in 1:ncol(MSEobj@OM)){
        cond<-MSEobj@OM[,i]>OMp[j,i]&MSEobj@OM[,i]<OMp[j+1,i]
        OMv[mm,i,j]<-mean(Ut[cond,mm],na.rm=T)
      }
      for(i in 1:ncol(MSEobj@Obs)){
        cond<-MSEobj@Obs[,i]>Obsp[j,i]&MSEobj@Obs[,i]<Obsp[j+1,i]
        Obsv[mm,i,j]<-mean(Ut[cond,mm],na.rm=T)
      }
    }
  }
  
  # -- Operating model variables
  OMs<-apply(OMv,1:2,sd,na.rm=T)
  OMstr<-array("",c(nMPs*2,ncomp+1))
   
  for(mm in 1:nMPs){
    ind<-order(OMs[mm,],decreasing=T)[1:ncomp]
    OMstr[1+(mm-1)*2,1]<-MPs[mm]
    OMstr[1+(mm-1)*2,2:(1+ncomp)]<-names(MSEobj@OM[ind])
    OMstr[2+(mm-1)*2,2:(1+ncomp)]<-round(OMs[mm,ind],2)
  }
  OMstr<-data.frame(OMstr)
  names(OMstr)<-c("MP",1:ncomp)
  
  
  # -- Observation model variables
  slots<-c( "Cat",  "Cat","AvC",  "AvC","CAA",      "CAA",    "CAL",      "CAL",    "Ind","Dep",  "Dep", "Dt",   "Dt", "Mort", "FMSY_M",    "BMSY_B0",     "L50",      "L95",    "LFC",    "LFS",    "Abun",  "Abun","vbK",  "vbt0",  "vbLinf",  "Steep","Iref",    "Cref",    "Bref")
  Obsnam<-c("Cbias","Csd","Cbias","Csd","CAA_nsamp","CAA_ESS","CAL_nsamp","CAL_ESS","Isd","Dbias","Derr","Dbias","Derr","Mbias","FMSY_Mbias","BMSY_B0bias", "lenMbias","lenMbias","LFCbias","LFSbias","Abias","Aerr","Kbias","t0bias","Linfbias","hbias","Irefbias","Crefbias","Brefbias")
  Obss<-apply(Obsv,1:2,sd,na.rm=T)
  Obsstr<-array("",c(nMPs*2,ncomp+1))
  for(mm in 1:nMPs){
    relobs<-Obsnam[slots%in%unlist(strsplit(Required(MPs[mm])[,2],split=", "))]
    ind<-(1:ncol(MSEobj@Obs))[match(relobs,names(MSEobj@Obs))]
    pos<-names(MSEobj@Obs)[ind]# possible observation processes
    maxy<-min(max(1,length(pos)),ncomp,na.rm=T)
    ind2<-order(Obss[mm,ind],decreasing=T)[1:maxy]
    Obsstr[1+(mm-1)*2,1]<-MPs[mm]
    Obsstr[1+(mm-1)*2,2:(1+maxy)]<-pos[ind2]
    Obsstr[2+(mm-1)*2,2:(1+maxy)]<-round(Obss[mm,ind][ind2],2)
  }
  Obsstr<-data.frame(Obsstr)
  names(Obsstr)<-c("MP",1:ncomp)

  ncols<-40
  #colsse<-makeTransparent(rainbow(ncols,start=0,end=0.36),95)[ncols:1]
  colsse<-makeTransparent(rainbow(ncols,start=0,end=0.36),90)[ncols:1]
  minsd<-0
  maxsd<-max(OMs,na.rm=T)
  coly<-ceiling(OMs/maxsd*ncols)
 
  # Operating model variables
  mbyp <- split(1:nMPs, ceiling(1:nMPs/maxrow))
  ylimy=c(0,max(OMv,na.rm=T)*1.2)              
  
  
  for(pp in 1:length(mbyp)){
    
    par(mfrow=c(length(mbyp[[pp]]),ncomp),mai=c(0.15,0.1,0.15,0.05),omi=c(0.1,0.9,0.3,0.05))
   
    for(mm in mbyp[[pp]]){
      for(cc in 1:ncomp){
        rind<-(mm-1)*2+1
        y<-Ut[,mm]
        cind<-match(OMstr[rind,1+cc],names(MSEobj@OM))
        x<-MSEobj@OM[,cind]
        plot(x,y,col="white",axes=F,ylim=ylimy)
        axis(1,pretty(OMp[,cind]),pretty(OMp[,cind]),cex.axis=0.8,padj=-1.5)
        abline(v=OMp[,cind],col="#99999960")
        points(x,y,col=colsse[coly[mm,cind]],pch=19,cex=0.8)
        x2<-(OMp[1:nbins,cind]+OMp[2:(nbins+1),cind])/2
        y2<-OMv[mm,cind,]
        lines(x2,y2)
        legend('bottomright',legend=round(OMs[mm,cind],2),bty='n',cex=0.8)
        legend('topleft',legend=OMstr[rind,1+cc],bty='n',cex=0.85)
        if(cc==1){ 
          mtext(MPs[mm],2,font=2,outer=F,cex=0.8,line=2)
          ytick<-pretty(seq(ylimy[1],ylimy[2]*1.3,length.out=10))
          axis(2,ytick,ytick,cex.axis=0.8)
        } # only first column
      } # parameters (columns)
    } # MPs (rows)
    
    mtext(Utnam,2,outer=T,cex=0.9,line=3.5)
    mtext(paste("Operating model parameters: ",objnam,"@OM",sep=""),3,outer=T,font=2,cex=0.9)
    
  } # Plots
  
  # Observation model values
  
  ylimy=c(0,max(Obsv,na.rm=T)*1.2)              
  minsd<-0
  maxsd<-max(Obss)
  coly<-ceiling(Obss/maxsd*ncols)
  
  if(sum(is.na(Obsstr)|Obsstr=="")<(ncomp+1)*nMPs*2-nMPs){ # only if there is data to plot
  
  for(pp in 1:length(mbyp)){
    
    par(mfrow=c(length(mbyp[[pp]]),ncomp),mai=c(0.15,0.1,0.15,0.05),omi=c(0.1,0.9,0.3,0.05))
    
    for(mm in mbyp[[pp]]){
      rind<-(mm-1)*2+1
      npres<-sum(Obsstr[rind+1,]!="")
      for(cc in 1:ncomp){
        if(!is.na(npres)&cc<(npres+1)){
          y<-Ut[,mm]
          cind<-match(Obsstr[rind,1+cc],names(MSEobj@Obs))
          x<-MSEobj@Obs[,cind]
          plot(x,y,col="white",axes=F,ylim=ylimy)
          axis(1,pretty(Obsp[,cind]),pretty(Obsp[,cind]),cex.axis=0.8,padj=-2)
          abline(v=Obsp[,cind],col="#99999960")
          points(x,y,col=colsse[coly[mm,cind]],pch=19,cex=0.8)
          x2<-(Obsp[1:nbins,cind]+Obsp[2:(nbins+1),cind])/2
          y2<-Obsv[mm,cind,]
          lines(x2,y2)
          legend('bottomright',legend=round(Obss[mm,cind],2),bty='n',cex=0.8)
          legend('topleft',legend=Obsstr[rind,1+cc],bty='n',cex=0.75)
          if(cc==1){ 
            mtext(MPs[mm],2,font=2,outer=F,cex=0.6,line=2)
            ytick<-pretty(seq(ylimy[1],ylimy[2]*1.3,length.out=10))
            axis(2,ytick,ytick,cex.axis=0.8)
          } # only first column
        }else{
          plot(0,type='n',axes=FALSE,ann=FALSE)
          if(cc==1){ 
            mtext(MPs[mm],2,font=2,outer=F,cex=0.6,line=2)
          } # only first column
        }  
      } # parameters (columns)
    } # MPs (rows)
    
    mtext(Utnam,2,outer=T,cex=0.9,line=3.5)
    mtext(paste("Observation model parameters: ",objnam,"@Obs",sep=""),3,outer=T,font=2,cex=0.9)
   
  } # Plots
  } # if there is data to plot
  
  list(OMstr,Obsstr)
  
} # VOI

# Manipulation of MSE Object ---------------------------------------------------
# Subset the MSE object by particular MPs (either MP number or name), 
#  or particular simulations
Sub <- function(MSEobj, MPs=NULL, sims=NULL, years=NULL) {
  Class <- class(MPs)
  if(Class == "NULL") subMPs <- MSEobj@MPs
  if(Class == "integer" | Class == "numeric") subMPs <- MSEobj@MPs[as.integer(MPs)]
  if(Class == "character") subMPs <- MPs
  SubMPs <- which(MSEobj@MPs %in% subMPs)
  not <- (subMPs %in% MSEobj@MPs) # Check for MPs misspelled
  ind <- which(not == FALSE)
  newMPs <- MSEobj@MPs[SubMPs]
  if (length(SubMPs) < 1) stop("MPs not found")
  if (length(ind > 0)) {
    message("Warning: MPs not found - ", paste0(subMPs[ind], " "))
	message("Subsetting by MPs: ", paste0(newMPs, " "))
  }
  
  ClassSims <- class(sims)
  if (ClassSims == "NULL") SubIts <- 1:MSEobj@nsim
  if (ClassSims == "integer" | ClassSims == "numeric") {
    # sims <- 1:min(MSEobj@nsim, max(sims))
	SubIts <- as.integer(sims)
  }	
  if (ClassSims == "logical") SubIts <- which(sims)

  ClassYrs <- class(years)
  AllNYears <- MSEobj@proyears
  if (ClassYrs == "NULL")  Years <- 1:AllNYears
  if (ClassYrs == "integer" | ClassYrs == "numeric") Years <- years
  if (max(Years) > AllNYears) stop("years exceeds number of years in MSE")
  if (min(Years) <= 0) stop("years must be positive")
  if (min(Years) != 1) {
    message("Not starting from first year. Are you sure you want to do this?")
    message("Probably a bad idea!")
  }
  if (!all(diff(Years) == 1)) stop("years are not consecutive")
  if (length(Years) <= 1) stop("You are going to want more than 1 projection year")
  MSEobj@proyears <- max(Years)
  
  SubF <- MSEobj@F_FMSY[SubIts,SubMPs,Years]
  SubB <- MSEobj@B_BMSY[SubIts,SubMPs,Years]
  OutOM <- MSEobj@OM[SubIts,]
  
  SubResults <- new('MSE',Name=MSEobj@Name, nyears=MSEobj@nyears, 
    proyears=MSEobj@proyears, nMPs=length(SubMPs), MPs=newMPs, 
	nsim=length(SubIts), OMtable=OutOM, Obs=MSEobj@Obs[SubIts,], 
	B_BMSYa=SubB, F_FMSYa=SubF, Ba=MSEobj@B[SubIts,SubMPs,Years], 
	FMa=MSEobj@FM[SubIts,SubMPs,Years], Ca=MSEobj@C[SubIts,SubMPs,Years], 
	TACa=MSEobj@TAC[SubIts,SubMPs,Years], SSB_hist=MSEobj@SSB_hist[SubIts,,,],
	CB_hist=MSEobj@CB_hist[SubIts,,,], FM_hist=MSEobj@FM_hist[SubIts,,,])
  
 return(SubResults)
}

# Join two or more MSE objects together 
joinMSE <- function(MSEobj1=NULL, MSEobj2=NULL, MSEobjs=NULL) {
  # join two MSE objects 
  if (class(MSEobj1) == "MSE" & class(MSEobj2) == "MSE") {
    obj1 <- MSEobj1
	obj2 <- MSEobj2
	if(!all(slotNames(MSEobj1) == slotNames(MSEobj2))) stop("The MSE objects don't have the same slots")
	if (any(MSEobj1@MPs != MSEobj2@MPs)) stop("MPs must be the same for all MSE objects")
	Nobjs <- 2 
  }
  if (class(MSEobjs) == "character") {
    Nobjs <- length(MSEobjs)
	for (X in 1:Nobjs) {
	  tt <- get(MSEobjs[X])
	  assign(paste0("obj", X), tt)
	}
	# Check 
	if (X > 1) {
	  tt <- get(MSEobjs[X])
	  tt2 <- get(MSEobjs[X-1])
	  if(!all(slotNames(tt) == slotNames(tt2))) stop("The MSE objects don't have the same slots")
	  if (any(tt@MPs != tt2@MPs)) stop("MPs must be the same for all MSE objects")
	}
  }
  if (class(MSEobjs) == "list") {
    Nobjs <- length(MSEobjs)
	for (X in 1:Nobjs) {
	  tt <- MSEobjs[[X]]
	  assign(paste0("obj", X), tt)
	}
	# Check 
	if (X > 1) {
	  tt <- MSEobjs[[X]]
	  tt2 <- MSEobjs[[X-1]]
	  if(!all(slotNames(tt) == slotNames(tt2))) stop("The MSE objects don't have the same slots")
	  if (any(tt@MPs != tt2@MPs)) stop("MPs must be the same for all MSE objects")
	}
  }
  
  # Check that nyears and proyears are the same for all 
  chkmat <- matrix(NA, nrow=Nobjs, ncol=2)
  nms <- NULL
  for (X in 1:Nobjs) {
    tt <- get(paste0("obj", X))
	chkmat[X, ] <- c(tt@nyears, tt@proyears)
	if (X > 1) if (!any(grepl(tt@Name, nms))) stop("MSE objects have different names")
	nms <- append(nms, tt@Name)
  }
  chk <- all(colSums(chkmat) == chkmat[1,] * Nobjs)
  if (!chk) stop("The MSE objects have different number of nyears or proyears")
  
  # Join them together
  Allobjs <- mget(paste0("obj", 1:Nobjs))
  sns <- slotNames(obj1)
  outlist <- vector("list", length(sns))
  for (sn in 1:length(sns)) {
    templs <- lapply(Allobjs, slot, name=sns[sn])
    if (class(templs[[1]]) == "character") {
	  outlist[[sn]] <- templs[[1]]
	}
	if (class(templs[[1]]) == "numeric" | class(templs[[1]]) == "integer") {
	  outlist[[sn]] <- do.call(c, templs)
	}
	if (class(templs[[1]]) == "matrix" | class(templs[[1]]) == "data.frame") {
	  outlist[[sn]] <- do.call(rbind, templs)
	}
	if (class(templs[[1]]) == "array") {
	  outlist[[sn]] <- abind(templs, along=1)
	}
  }
  names(outlist) <- sns
 
  newMSE <- new('MSE', Name=outlist$Name,nyears=unique(outlist$nyears),
    proyears=unique(outlist$proyears), nMP=unique(outlist$nMP), 
	MPs=unique(outlist$MPs), nsim=sum(outlist$nsim),OM=outlist$OM,
	Obs=outlist$Obs, B_BMSY=outlist$B_BMSY, F_FMSY=outlist$F_FMSY,
    outlist$B, outlist$FM, outlist$C, outlist$TAC, 
	outlist$SSB_hist, outlist$CB_hist, outlist$FM_hist)
 
  newMSE
}	
  
# Evaluate Peformance of MPs ---------------------------------------------------
# Function examines how consistently an MP outperforms another. 
DOM <- function(MSEobj, MPtg=NA) {
  if (any(is.na(MPtg))) MPtg <- MSEobj@MPs 
  proyears<-MSEobj@proyears
  nMP <- MSEobj@nMPs
  ind <- which(MSEobj@MPs %in%  MPtg)
  MPr <- which(!(MSEobj@MPs %in%  MPtg))
  yind<-max(MSEobj@proyears-4,1):MSEobj@proyears
  y1<-1:(MSEobj@proyears-1)
  y2<-2:MSEobj@proyears
  Mat <- matrix(0, nrow=length(MPtg), ncol=nMP)
  rownames(Mat) <- MPtg
  colnames(Mat) <- MSEobj@MPs
  POF <- P100 <- YieldMat <- IAVmat <- Mat
  for (X in 1:length(MPtg)) {
    # Overfishing (F > FMSY)
    ind1 <- as.matrix(expand.grid(1:nsim, ind[X], 1:proyears))
    ind2 <- as.matrix(expand.grid(1:nsim,  1:nMP, 1:proyears))
    t1 <- apply(array(MSEobj@F_FMSY[ind1] > 1, dim=c(nsim, 1, proyears)), c(1,2), sum, na.rm=TRUE) 
    t2 <- apply(array(MSEobj@F_FMSY[ind2] > 1, dim=c(nsim, nMP, proyears)), c(1,2), sum, na.rm=TRUE)
    POF[X,] <- round(apply(matrix(rep(t1, nMP), nrow=nsim) < t2, 2, sum) / nsim * 100, 0)
	# B < BMSY
    t1 <- apply(array(MSEobj@B_BMSY[ind1] < 1, dim=c(nsim, 1, proyears)), c(1,2), sum, na.rm=TRUE) 
    t2 <- apply(array(MSEobj@B_BMSY[ind2] < 1, dim=c(nsim, nMP, proyears)), c(1,2), sum, na.rm=TRUE) 
    P100[X,] <- round(apply(matrix(rep(t1, nMP), nrow=nsim) < t2, 2, sum, na.rm=TRUE) / nsim * 100, 0)
    # Relative yield in last 5 years
    ind1 <- as.matrix(expand.grid(1:nsim, ind[X], yind))
    ind2 <- as.matrix(expand.grid(1:nsim,  1:nMP, yind))
    t1 <- apply(array(MSEobj@C[ind1], dim=c(nsim, 1, length(yind))), c(1,2), sum, na.rm=TRUE)
    t2 <- apply(array(MSEobj@C[ind2], dim=c(nsim, nMP, length(yind))), c(1,2), sum, na.rm=TRUE)
    YieldMat[X,] <- round(apply(matrix(rep(t1, nMP), nrow=nsim) > t2, 2, sum, na.rm=TRUE) / nsim * 100, 0)
    # interannual variation in catch 
    ind1 <- as.matrix(expand.grid(1:nsim, ind[X], y1))
    ind2 <- as.matrix(expand.grid(1:nsim, ind[X], y2)) 
    AAVY1 <- apply(array(((MSEobj@C[ind1]-MSEobj@C[ind2])^2)^0.5, dim=c(nsim, 1, length(y1))),1,mean,na.rm=T)/apply(array(MSEobj@C[ind2], dim=c(nsim, 1, length(y1))),1,mean,na.rm=T) 
    ind1 <- as.matrix(expand.grid(1:nsim, 1:nMP, y1))
    ind2 <- as.matrix(expand.grid(1:nsim, 1:nMP, y2)) 
    AAVY2 <- apply(array(((MSEobj@C[ind1]-MSEobj@C[ind2])^2)^0.5, dim=c(nsim, nMP, length(y1))),c(1,2),mean,na.rm=T)/apply(array(MSEobj@C[ind2], dim=c(nsim, nMP, length(y1))),c(1,2),mean,na.rm=T)
    IAVmat[X,] <- round(apply(matrix(rep(AAVY1, nMP), nrow=nsim) < AAVY2, 2, sum, na.rm=TRUE) / nsim * 100, 0)
  }
  out <- list() 
  out$POF <- POF
  out$P100 <- P100
  out$Yd <- YieldMat
  out$AAVY <- IAVmat
  return(out)
}

# Trade-Off Plot Function ------------------------------------------------------
TradePlot <- function(MSEobj, XAxis=c("Overfishing", "Biomass:BMSY"), 
	YAxis=c("Long-term Yield", "AnnualVar"), XThresh=c(30, 80), YThresh=c(0,50),
	maxVar=15, BmsyRef=0.5, B0Ref=0.2, AvailMPs=NULL, ShowLabs=FALSE, 
	ShowCols=TRUE) {
  PMs <- c("Long-term Yield", "Short-term Yield", "Overfishing", "Biomass:BMSY",
	"Biomass:B0", "AnnualVar")
  # Error Checks 	
  if (prod(XAxis %in% PMs)!=1) {
    message("Available Performance Metrics")
    print(PMs)
    stop("Invalid XAxis Performance Metrics")
  }	
  if (prod(YAxis %in% PMs)!=1) {
    message("Available Performance Metrics")
    print(PMs)
    stop("Invalid YAxis Performance Metrics")
  }	
  if (length(XAxis) > 4) stop("Too many Performance Metrics (max 4)")
  if (length(YAxis) > 4) stop("Too many Performance Metrics (max 4)")
  if (length(XAxis) != length(YAxis)) stop("XAxis must be of same length as YAxis")
  if (length(XThresh) != length(XAxis) | length(YThresh) != length(XAxis)) 
	warning("Risk Threshold not same length as number of PMs")
  
  Yd<-rep(NA,MSEobj@nMPs)
  BMSYref<-rep(NA,MSEobj@nMPs)
  B0ref<-rep(NA,MSEobj@nMPs)
  PNOF<-rep(NA,MSEobj@nMPs)
  LTY<-rep(NA,MSEobj@nMPs)
  STY<-rep(NA,MSEobj@nMPs)
  VY<-rep(NA,MSEobj@nMPs)

  y1 <- 1:(MSEobj@proyears-1)
  y2 <- 2:MSEobj@proyears
  
  ystart<-1:5
  yend<-max(MSEobj@proyears-4,1):MSEobj@proyears
  
  RefYd<-MSEobj@OM$RefY
  if (maxVar < 1) maxVar <- maxVar * 100
  
  for(mm in 1:MSEobj@nMPs){  
    PNOF[mm]<-round(sum(MSEobj@F_FMSY[,mm,]<1,na.rm=T)/prod(dim(MSEobj@F_FMSY[,mm,]),na.rm=T)*100,1)
    BMSYref[mm]<-round(sum(MSEobj@B_BMSY[,mm,]>BmsyRef,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
	B0ref[mm]<-round(sum((MSEobj@B_BMSY[,mm,] * MSEobj@OM$BMSY_B0) >B0Ref,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
    # LTY[mm]<-round(sum(MSEobj@C[,mm,yend]/RefYd>0.5,na.rm=T)/(MSEobj@nsim*length(yend)),3)*100
	# STY[mm]<-round(sum(MSEobj@C[,mm,ystart]/RefYd>0.5,na.rm=T)/(MSEobj@nsim*length(ystart)),3)*100
	LTY[mm]<-round(mean(apply(MSEobj@C[,mm,yend],1,mean,na.rm=T)/RefYd,na.rm=T)*100,1)
	STY[mm]<-round(mean(apply(MSEobj@C[,mm,ystart],1,mean,na.rm=T)/RefYd,na.rm=T)*100,1)
    AAVY<-apply((((MSEobj@C[,mm,y1]-MSEobj@C[,mm,y2])/MSEobj@C[,mm,y2])^2)^0.5,1,mean,na.rm=T) 
    VY[mm]<-round(sum(AAVY<(maxVar/100),na.rm=T)/MSEobj@nsim,3)*100
  }
  
  for (xx in seq_along(XAxis)) {
    name <- paste0("X", xx)
	name1 <- paste0("XLab", xx)
    assign(name, GetStat(XAxis[xx], LTY, STY, PNOF, BMSYref, B0ref, VY))
	assign(name1, StatLab(XAxis[xx], maxVar, BmsyRef, B0Ref))
	name <- paste0("Y", xx)
	name1 <- paste0("YLab", xx)
	assign(name, GetStat(YAxis[xx], LTY, STY, PNOF, BMSYref, B0ref, VY))
	assign(name1, StatLab(YAxis[xx], maxVar, BmsyRef, B0Ref))
  }
  
  Nplot <- length(XAxis)
  if (Nplot == 1) par(mfrow=c(1,1), mar=c(4,4.5,1,1), oma=c(1,1,0,0))
  if (Nplot == 2) par(mfrow=c(1,2), mar=c(4,4.5,1,1), oma=c(1,1,0,0))
  if (Nplot == 3) par(mfrow=c(1,3), mar=c(4,4.5,1,1), oma=c(1,1,0,0))
  if (Nplot == 4) par(mfrow=c(2,2), mar=c(4,4.5,1,1), oma=c(1,1,0,0))
  
  OutList <- list()
  for (xx in seq_along(XAxis)) {
    Xname <- paste0("X", xx)
	XLab <- paste0("XLab", xx)
	Yname <- paste0("Y", xx)
	YLab <- paste0("YLab", xx)
    rr <- tradeoffplot4(x=get(Xname), y=get(Yname), get(XLab), get(YLab), 
		labs=MSEobj@MPs[1:MSEobj@nMPs],vl=XThresh[xx],hl=YThresh[xx], 
		ShowLabs=ShowLabs,  ShowCols=ShowCols, AvailMPs=AvailMPs)
	
	labs <- MSEobj@MPs[1:MSEobj@nMPs]
	ind <- which(labs %in% rr)
    tempDF <- data.frame(MP=rr, X=get(Xname)[ind], Y=get(Yname)[ind])
	Dist <- NULL # calculate distance from corner
    for (X in 1:length(tempDF[,2])) Dist[X] <- euc.dist(c(tempDF[X,2], tempDF[X,3]), c(100, 100))
	tempDF <- tempDF[order(Dist),]
	rownames(tempDF) <- 1:nrow(tempDF)
	OutList[[xx]] <- tempDF
  }
 
  print(OutList)
  invisible(OutList)
  
}

# Supporting functions 
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

GetStat <- function(PM, LTY, STY, PNOF, BMSYref, B0ref, VY) {
  switch(PM,
    "Long-term Yield" = LTY,
	"Short-term Yield" = STY,
	"Overfishing" = PNOF,
	"Biomass:BMSY" = BMSYref,
	"Biomass:B0" = B0ref,
    "AnnualVar" = VY)
}

StatLab <- function(PM, maxVar, BmsyRef, B0Ref) {
  switch(PM,
    "Long-term Yield" = "Long-term Yield",
	"Short-term Yield" = "Short-term Yield",
	"Overfishing" = "Prob. of Not Overfishing (%)",
	"Biomass:BMSY" = paste0("Prob. Biomass >",BmsyRef, "BMSY (%)"),
	"Biomass:B0" = paste0("Prob. Biomass >",B0Ref, "B0 (%)"),,
    "AnnualVar" = paste0("Prob. AAVY <", maxVar, "%")
	)
}

tradeoffplot4<-function(x,y,xlab,ylab,labs,cex,vl,hl, 
	ShowLabs=FALSE,  ShowCols=FALSE, AvailMPs=NULL){
   adjj<-c(0.9,1.1)
   XLim <- c(min(c(-10, min(x,na.rm=T)*adjj)), max(c(max(x,na.rm=T)*adjj, 110)))
   YLim <- c(min(c(-10, min(y,na.rm=T)*adjj)), max(c(max(y,na.rm=T)*adjj, 110)))
   
   # Which MPs meet minimum PMs 
   ind <- which(x >= vl & y >=hl)
   coly <- rep("darkgray", length(labs)) 
   coly[ind] <- "black" 
   # coly[labs%in%c("AvC","curE","FMSYref")]<-'black'
   
   Pch <- rep(21, length(labs))
   # Pch[labs%in%c("AvC","curE","FMSYref")] <- 17
   coly[grep("FMSY", labs)]<-'black'
   Pch[grep("FMSY", labs)] <- 24
   if (!is.null(AvailMPs)) Pch[labs%in%AvailMPs] <- 21
   if (!is.null(AvailMPs)) coly[labs%in%AvailMPs & (x >= vl & y >=hl)] <- "green"
   # coly<-rep(c('#0000ff95','#ff000095','#20ff1095'),50)[1:length(labs)]

   plot(NA,xlim=XLim,ylim=YLim,xlab=xlab,ylab=ylab, bty="l", las=1)
   abline(v=vl,col="#99999940",lwd=2)
   abline(h=hl,col="#99999940",lwd=2)
   
   Alpha <- 30
   # polygons 
   LeftCol <- rgb(red=255, green=0, blue=0, alpha=Alpha, names = NULL, 
	maxColorValue = 255)
   RightCol <- rgb(red=0, green=255, blue=0, alpha=Alpha, names = NULL, 
	maxColorValue = 255)   

   if(ShowCols) {
     polygon(x=c(0, vl,  vl, 0), y=c(0, 0, hl, hl), col=LeftCol, border=NA)
     polygon(x=c(0, vl,  vl, 0), y=c(0, 0, 100, 100), col=LeftCol, border=NA)
     polygon(x=c(vl,  100, 100, vl), y=c(0, 0, 100, 100), col=RightCol, border=NA)
     polygon(x=c(vl, 100,  100, vl), y=c(hl, hl, 100, 100), col=RightCol, border=NA)
    }
   
    Cex <- 1.5
   if(!ShowLabs) points(x,y, bg=coly, pch=Pch, cex=Cex, col="black" )
   if(ShowLabs) text(x,y,labs,font=2,col=coly,cex=1)
   # if(IdPoints) {
    # message("Click points on plot to display MP name")
	# message("Click Stop to finish")
	# flush.console()
	# identify(x,y, labels=labs)
   # }	
   
   labs[ind]
   
}

wormplot<-function(MSEobj,Bref=0.5,LB=0.25,UB=0.75){
  
  if(UB<LB)stop("LB parameter must be lower than UB parameter")
  if(LB<0|LB>1)stop("LB parameter must be in the range of 0 to 1")
  if(UB<0|UB>1)stop("UB parameter must be in the range of 0 to 1")
  
  ncol<-ceiling(MSEobj@nMPs^0.3)
  nrow<-ceiling(MSEobj@nMPs/ncol)
  
  par(mfcol=c(nrow,ncol),mar=c(0.1,0.1,0.1,0.1),omi=c(0.6,0.25,0.3,0))
  
  Bprob<-apply(MSEobj@B_BMSY>Bref,2:3,sum)/MSEobj@nsim
  
  ind<-order(apply(Bprob,1,sum),decreasing=T)
  
  BLB<-Bprob>LB
  BUB<-Bprob>UB
  
  col<-array('red',dim(Bprob))
  col[BLB&!BUB]="yellow"
  col[BUB]="green"
  
  for(i in 1:(nrow*ncol)){
    if(i<(MSEobj@nMPs+1)){
      MP<-ind[i]
      plot(c(1,MSEobj@proyears+2),c(-1,1),col='white',axes=F)
      # abline(h=0)
    
      for(ys in 1:MSEobj@proyears){
        x<-c(ys-1,ys,ys,ys-1)
        y<-c(rep(Bprob[MP,ys],2),rep(-Bprob[MP,ys],2))
        pol<-data.frame(x,y)
        polygon(pol,col=col[MP,ys],border=NA)
      }
    
      legend('top',legend=MSEobj@MPs[MP],bty='n')
      if((i/nrow)==round(i/nrow,0)) axis(1,pretty(1:MSEobj@proyears),pretty(1:MSEobj@proyears))
      
      
    }else{
      
      plot.new()
 
    }
    
    if(i==(nrow*ncol)){
      legend('topright',fill=c("green","red"),
             legend=c(paste(">",round(UB*100,0),"% prob.",sep=""),
                      paste("<",round(LB*100,0),"% prob.",sep="")
             ),
             bty="n")
      
    }
    
  }
  
  mtext(paste("Probability of biomass above ",round(Bref*100,0),"% BMSY for ",deparse(substitute(MSE)),sep=""),3,outer=T,line=0.5)
  mtext("Projection year",1,outer=T,line=2.5)
  mtext(paste("Fraction of simulations above ",round(Bref*100,0),"% BMSY",sep=""),2,outer=T,line=0.25)
  Bprob
  
}

VOI2<-function(MSEobj,ncomp=6,nbins=4,Ut=NA,Utnam="yield",lay=F){
  
  objnam<-deparse(substitute(MSEobj))
  nsim<-MSEobj@nsim
  
  if(is.na(Ut[1])){
    Ut<-array(NA,c(nsim,MSEobj@nMPs))
    yind<-max(MSEobj@proyears-4,1):MSEobj@proyears
    RefYd<-MSEobj@OM$RefY
    
    for(mm in 1:MSEobj@nMPs){
      Ut[,mm]<-apply(MSEobj@C[,mm,yind],1,mean,na.rm=T)/RefYd*100
      #POF[,mm]<-apply(MSEobj@F_FMSY[,mm,]>1,1,sum)/MSEobj@proyears
      #P10[,mm]<-apply(MSEobj@B_BMSY[,mm,]<0.1,1,sum)/MSEobj@proyears
    }
    
  }
  
  MPs<-MSEobj@MPs
  nMPs<-MSEobj@nMPs
  
  # -- Observation model variables
  slots<-c( "Cat",  "Cat","AvC",  "AvC","CAA",      "CAA",    "CAL",      "CAL",    "Ind","Ind",  "Dep",  "Dep", "Dt",   "Dt", "Mort", "FMSY_M",    "BMSY_B0",     "L50",      "L95",    "LFC",    "LFS",    "Abun",  "Abun","vbK",  "vbt0",  "vbLinf",  "Steep","Iref",    "Cref",    "Bref")
  Obsnam<-c("Cbias","Csd","Cbias","Csd","CAA_nsamp","CAA_ESS","CAL_nsamp","CAL_ESS","Isd","betas","Dbias","Derr","Dbias","Derr","Mbias","FMSY_Mbias","BMSY_B0bias", "lenMbias","lenMbias","LFCbias","LFSbias","Abias","Aerr","Kbias","t0bias","Linfbias","hbias","Irefbias","Crefbias","Brefbias")
  
  
  Obsnam2<-c("Cbias","Csd","CAA_nsamp","CAA_ESS","CAL_nsamp","CAL_ESS","Isd","betas","Dbias","Derr","Mbias","FMSY_Mbias","BMSY_B0bias", "lenMbias","LFCbias","LFSbias","Abias","Aerr","Kbias","t0bias","Linfbias","hbias","Irefbias","Crefbias","Brefbias")
  Obsnam3<-c("Catch bias","Catch error","n CAA samples", "CAA ESS","n CAL samples","CAL ESS","Abun. Ind. error","Hyperstability","Depln. bias",
             "Depln. error","M bias","FMSY/M bias","BMSY/B0 bias","lenMbias","Len 1st Cap bias","Len full sel bias","Cur Abund bias","Cur Abun err","vB K bias","vB t0 bias","vB Linf bias","Steepness bias","Ref index bias","Ref catch bias", "Ref biomass bias")
  #Types of observation error model   1:lognorm   2:percentile  3:replicates (higher is better) ##4:uniform on log  5:logit space
  oem<-c(     1,      2,    3,          3,        3,          3,        2,    4,      1,      2,     1,      1,            1,             1,        1,        1,        4,      2,     1,      1,       1,         2,      1,         1,         1)
  #oem<-c(     2,      2,    3,          3,        3,          3,        2,    4,      2,      2,     2,      2,            2,             2,        2,        2,        4,      2,     2,      2,       2,         2,      1,         1,         1)
  
  Obsd<-apply(MSEobj@Obs,2,sd)
  Obm<-apply(MSEobj@Obs,2,mean)
  Obmd<-apply(MSEobj@Obs,2,quantile,p=0.5)
  
  maxcomp<-length(Obsnam2) 
  Obsv<-array(NA,c(nMPs,maxcomp,nbins))
  Obsval<-array(NA,c(nMPs,maxcomp,nbins))
  Obscost<-array(NA,c(nMPs,maxcomp,nbins))
  Obsname<-list()
  
  div<-seq(1,2,length.out=nbins+1)[2:(nbins+1)] # for distributions
  percs<-seq(0.5,1,length.out=nbins+1)[1:nbins] # for samples
  percsCAA<-seq(0,1,length.out=nbins+2)[2:(nbins+1)]
  percUL<-seq(0,0.25,length.out=nbins+1)[2:(nbins+1)]
  percUU<-1-percUL
  for(mm in 1:nMPs){
    Y1<-Ut[,mm]
    relobs<-Obsnam[slots%in%unlist(strsplit(Required(MPs[mm])[,2],split=", "))]
    Obsname[[mm]]<-relobs
    nr<-length(relobs)
    if(length(relobs)>0){
      
      for(r in 1:nr){
        oemi<-match(relobs[r],Obsnam2)
        obsi<-match(relobs[r],names(MSEobj@Obs))
        for(cc in 1:nbins){
          if(oem[oemi]==1){ # Redundant SIR code for log-normal biases
            T1<-tdlnorm(MSEobj@Obs[,obsi],Obm[obsi],Obsd[obsi]/Obm[obsi]) 
            #plot(MSEobj@Obs[,obsi],T1) # check
            T2<-tdlnorm(MSEobj@Obs[,obsi],Obm[obsi],Obsd[obsi]/(div[cc]*Obm[obsi])) 
            W<-T2/T1
            nrep2<-nsim*20
            Y2<-sample(Y1,nrep2*5,replace=T,prob=W)
            Obsv[mm,r,cc]<-(mean(Y2)-mean(Y1))/mean(Y1)*100
            Obsval[mm,r,cc]<-Obsd[obsi]/(div[cc]*Obm[obsi])
            Obscost[mm,r,cc]<-div[cc]^2
          }else if(oem[oemi]==2){
            refval<-quantile(MSEobj@Obs[,obsi],percs[nbins:1][cc])
            ind<-MSEobj@Obs[,obsi]<refval
            Obsv[mm,r,cc]<-(mean(Y1[ind])-mean(Y1))/mean(Y1)*100
            Obsval[mm,r,cc]<-mean(MSEobj@Obs[ind,obsi])
            Obscost[mm,r,cc]<-1/(Obsval[mm,r,cc]/mean(MSEobj@Obs[,obsi]))^2
          }else if(oem[oemi]==3){  
            refval<-quantile(MSEobj@Obs[,obsi],percsCAA[cc])
            ind<-MSEobj@Obs[,obsi]>refval
            Obsv[mm,r,cc]<-(mean(Y1[ind])-mean(Y1))/mean(Y1)*100
            Obsval[mm,r,cc]<-mean(MSEobj@Obs[ind,obsi])
            Obscost[mm,r,cc]<-Obsval[mm,r,cc]
          }else if(oem[oemi]==4){
            refval<-quantile(MSEobj@Obs[,obsi],percUL[cc])
            refval2<-quantile(MSEobj@Obs[,obsi],percUU[cc])
            ind<-(MSEobj@Obs[,obsi]>refval)&(MSEobj@Obs[,obsi]<refval2)
            Obsv[mm,r,cc]<-(mean(Y1[ind])-mean(Y1))/mean(Y1)*100
            Obsval[mm,r,cc]<-sd(MSEobj@Obs[ind,obsi])
            Obscost[mm,r,cc]<-1/(Obsval[mm,r,cc]/sd(MSEobj@Obs[,obsi]))^2
          }
          # observation model type
        } # loop over bins
      } # loop over r
    } # observation variables?
  } # loop over MPs
  
  cb<-array(NA,c(MSEobj@nMPs,maxcomp))
  for(mm in 1:MSEobj@nMPs){
    if(sum(!is.na(Obscost[mm,,]))>0){
      for(r in 1:length(Obsname[[mm]])){ 
        dat<-data.frame(x=Obscost[mm,r,],y=Obsv[mm,r,])
		if (prod(apply(dat, 2, is.finite)) > 0) {
        #plot(dat$x,dat$y)
          cb[mm,r]<-lm(y~x-1,data=dat)$coefficients[1]
		}  
      }
    }
  } 
  
  ncols<-100
  #colsse<-makeTransparent(rainbow(ncols,start=0,end=0.36),95)[ncols:1]
  colt<-rainbow(ncols,start=0,end=0.36)[1:ncols]
  colsse<-makeTransparent(colt,98)
  
  cb[cb<0|is.na(cb)]<-0
  coly<-ceiling((cb/max(cb,na.rm=T))^0.5*ncols)
  coly[coly==0]<-1
  
  ncol<-ceiling(MSEobj@nMPs^0.5)
  nrow<-ceiling(MSEobj@nMPs/ncol)
  
  par(mfrow=c(nrow,ncol),mar=c(2.4,2.4,0.1,0.1),omi=c(0.4,0.35,0.3,0))
  
  gcol1<-"#99999960"
  gcol2<-"#99999940"
  gcol3<-"#99999920"
  
  for(mm in 1:MSEobj@nMPs){
    if(sum(!is.na(Obscost[mm,,]))>0){
      
      plot(c(1,5),range(Obsv,na.rm=T),col='white',main="")
      legend('topleft',legend=MSEobj@MPs[mm],bty='n',text.font=2,cex=1.4)
      
      abline(h=(-20:50)*4,col=gcol2,lwd=1.5)
      abline(h=(-20:50)*4+2,col=gcol3,lwd=1)
      
      abline(v=1:4,col=gcol2,lwd=1.5)
      abline(v=(1:4)+0.5,col=gcol3,lwd=1)
      abline(h=0,col=gcol1,lwd=3)
      
      
      no<-length(Obsname[[mm]])
      ind<-order(cb[mm,1:no],decreasing=T)[1:ncomp]
      ind<-ind[!is.na(ind)]
      
      ind2<-order(Obsv[mm,1:no,nbins])[1:ncomp]
      ind2<-ind2[!is.na(ind2)]
      
      
      lpos<-Obsv[mm,ind2,nbins]
      ppos<-seq(min(Obsv,na.rm=T),max(Obsv,na.rm=T),length.out=length(ind2))
      wt<-(max(Obsv[mm,1:no,nbins],na.rm=T)-min(Obsv[mm,1:no,nbins],na.rm=T))/(max(Obsv,na.rm=T)-min(Obsv,na.rm=T))/(no/ncomp)
      wt<-wt^0.66
      nupos<-wt*lpos+(1-wt)*ppos
      
      for(r2 in 1:length(ind)){
        r<-ind2[r2]
        lines(c(1,Obscost[mm,r,]),c(0,Obsv[mm,r,]),col=colsse[coly[mm,r]],lwd=3)
        if(!lay){
          text(4.5,nupos[r2],Obsname[[mm]][r],col=colt[coly[mm,r]],font=2,cex=1.2)
        }else{
          text(4.5,nupos[r2],Obsnam3[match(Obsname[[mm]][r],Obsnam2)],col=colt[coly[mm,r]],font=2,cex=1.2)
        }
      } # observation quantities (lines)
      
      
      #legend('topleft',legend=Obsname[[mm]][ind],text.col=colt[coly[mm,ind]],text.font=2,cex=1.2,bty='n')
      
    } # if there is data to plot
    
    
  } # MPs (plots)
  
  
  mtext("Cost relative to today",1,outer=T,cex=0.9,line=1,font=2)
  #mtext(paste("Operating model parameters: ",objnam,"@OM",sep=""),3,outer=T,font=2,cex=0.9)
  mtext(paste("% Change in ",Utnam," relative to today",sep=""),2,outer=T,line=0.6,font=2,cex=0.9)
  
  list(Obscost,Obsv,Obsval,cb,Obsname,MSEobj@MPs)
  
} # VOI2

##############
# These are temporary plotting functions
# Need to be cleaned up!


calcperf <- function(MSEobj, 
  Lims=c("B_BMSY", "B_B0", "F_FMSY", "AAVY", "AAVE"), 
  Refs=c(0.5, 0.2, 1.5, 30, 30), lastYrs=10, allyrs=FALSE) {
 
 # Assign Variables
  Nyears <- MSEobj@nyears
  Pyears <- MSEobj@proyears
  nMPs <- MSEobj@nMPs
  MPs <- MSEobj@MPs 
  nsim <- MSEobj@nsim
  RefYd <- MSEobj@OM$RefY
 
  # Error Checks 
  if (lastYrs >= Pyears) {
    message("lastYrs set too high. Defaulting to last 10 years")
    lastYrs <- 10
  }
  if (lastYrs <= 0 | lastYrs == FALSE | lastYrs == "all" | lastYrs == "All") 
    lastYrs <- Pyears   
  if (length(Refs) != length(Lims)) stop("Length of Lims and Refs must be equal")
  # Performance metrics 
  allMets <- c("B_BMSY", "B_B0", "F_FMSY", "AAVY", "AAVE")
  dfRefs <- c(0.5, 0.2, 1.25, 20, 20) # Default reference points
  PMdf <- data.frame(allMets, dfRefs) # Default performance metrics 
  
  Lims <- match.arg(Lims, allMets, several.ok=TRUE)
  myPMdf <- PMdf[match(Lims, PMdf[,1]),]
  myPMdf[,2] <- Refs
  
  BmsyRef <- myPMdf[grep("B_BMSY", myPMdf[,1]),2]
  B0Ref <- myPMdf[grep("B_B0", myPMdf[,1]),2]
  FRef <- myPMdf[grep("F_FMSY", myPMdf[,1]),2]
  maxVar <- myPMdf[grep("AAVY", myPMdf[,1]),2]
  maxVarE <- myPMdf[grep("AAVE", myPMdf[,1]),2]
  if (length(BmsyRef)> 1 | length(B0Ref)> 1 | length(B0Ref)> 1
    | length(maxVar) > 1 | length(maxVarE) > 1) 
	stop("Repeated performance metrics")
  
  if(length(BmsyRef) < 1) BmsyRef <- PMdf[grep("B_BMSY", PMdf[,1]),2]
  if(length(B0Ref) < 1) B0Ref <- PMdf[grep("B_B0", PMdf[,1]),2]
  if(length(FRef) < 1) FRef <- PMdf[grep("F_FMSY", PMdf[,1]),2]
  if(length(maxVar) < 1) maxVar <- PMdf[grep("AAVY", PMdf[,1]),2]
  if(length(maxVarE) < 1) maxVarE <- PMdf[grep("AAVE", PMdf[,1]),2]
  
  Names <- list(5)
  Names[[1]] <- bquote(italic(B) > ~ .(BmsyRef) ~ italic(B[MSY]))
  Names[[2]] <- bquote(italic(B) > ~ .(B0Ref) ~ italic(B[0]))
  Names[[3]] <- bquote(italic(F) < ~ .(FRef) ~ italic(F[MSY]))
  Names[[4]] <- bquote(AAVY < ~ .(maxVar)* '%')
  Names[[5]] <- bquote(AAVE < ~ .(maxVarE)* '%')
  
  myRefs <- list(BmsyRef=BmsyRef, B0Ref=B0Ref, FRef=FRef, maxVar=maxVar, maxVarE=maxVarE)
  myRefsind <- as.logical(unlist(lapply(myRefs, length)))
    
  if (length(maxVar) > 0 && maxVar < 1) maxVar <- maxVar * 100 # maximum variability in % 
  if (length(maxVarE) > 0 && maxVarE < 1) maxVarE <- maxVarE * 100 # maximum variability in %
  
  yrs <- (Pyears-lastYrs+1):Pyears 	# Years to summarize performance
  ystart <- 1:lastYrs			 	# First 10 years
  y1 <- yrs[1]:yrs[length(yrs)-1] 	# for calculating interannual variability
  y2 <- y1+1
 
  # Satisficing matrices 
  B_BMSY  <- array(NA, dim=c(nsim, nMPs))
  B_B0    <- array(NA, dim=c(nsim, nMPs))
  F_FMSY  <- array(NA, dim=c(nsim, nMPs))
  AAVY    <- array(NA, dim=c(nsim, nMPs))
  AAVE    <- array(NA, dim=c(nsim, nMPs)) # Need to add this to MSE obj 

  # Maximizing & Reporting 
  LTY <- rep(NA, nMPs)
  STY <- rep(NA, nMPs)
 
  # Calculate Joint Probability
  # Probability that MP exceeds metric for all of lastYrs  
  print("Calculating probability of meeting performance metrics")
  flush.console()
  
  if (!allyrs) {
    B_BMSY2 <- MSEobj@B_BMSY[,,yrs] > BmsyRef
    B_B02 <- (MSEobj@B_BMSY[,,yrs] * MSEobj@OM$BMSY_B0) > B0Ref
    F_FMSY2 <- MSEobj@F_FMSY[,,yrs] < FRef
    
    # Interannual variability in yield
    MSEobj@C[(!is.finite(MSEobj@C[,,]))] <- 0 # if catch is NAN or NA, make it 0
    aavy <- (((MSEobj@C[,,y1]-MSEobj@C[,,y2])/MSEobj@C[,,y2])^2)^0.5
    aavy[is.nan(aavy)] <- 1
    AAVY2	<- aavy < (maxVar/100)
    
    # Interannual variability in effort (actually fishing mortality here) - need to change MSE object
    aave <- (((MSEobj@FM[,,y1]-MSEobj@C[,,y2])/MSEobj@FM[,,y2])^2)^0.5
    aave[is.nan(aave)] <- 1
    AAVE2	<- aave < (maxVar/100) 
    
    Stats2 <- list(B_BMSY=B_BMSY2, B_B0=B_B02, F_FMSY=F_FMSY2, AAVY=AAVY2, AAVE=AAVE2)  
    
    Ind2 <- match(myPMdf[,1], names(Stats2)) # Specified Joint Performance Metrics
    
    Legend <- NULL
    for (X in 1:length(Ind2)) Legend <- append(Legend, as.expression(Names[[Ind2[X]]]))
    
    LTY <-sapply(1:nMPs, function(mm) 
      round(mean(apply(MSEobj@C[,mm,yrs],1,mean,na.rm=TRUE)/RefYd)*100,1))
    STY <-sapply(1:nMPs, function(mm) 
     round(mean(apply(MSEobj@C[,mm,ystart],1,mean,na.rm=TRUE)/RefYd)*100,1))
     
    # individual probabilities
    bmsy <- apply(B_BMSY2, 2, sum, na.rm=TRUE)/(lastYrs * nsim)
    b0 <- apply(B_B02, 2, sum, na.rm=TRUE)/(lastYrs * nsim)
    fmsy <- apply(F_FMSY2, 2, sum, na.rm=TRUE)/(lastYrs * nsim)
    aavy <- apply(AAVY2, 2, sum, na.rm=TRUE)/(lastYrs * nsim)
    aave <- apply(AAVE2, 2, sum, na.rm=TRUE)/(lastYrs * nsim)
    
    indivProbs2 <- rbind(bmsy, b0, fmsy, aavy, aave)
    Ind2 <- match(myPMdf[,1], names(Stats2)) # Specified Joint Performance Metrics
    indivProbs2 <- indivProbs2[Ind2,]
    colnames(indivProbs2) <- MPs 
    indivProbs2 <- indivProbs2[,order(indivProbs2[1,])] * 100
    
    # Joint probabilities
	total2 <- Reduce("*", Stats2[Ind2]) 
    prob2 <- apply(total2, 2, sum, na.rm=TRUE)/(lastYrs * nsim) * 100
    jointDF <- data.frame(MPs=MPs, Prob=prob2, LTY=LTY, STY=STY)
  } 
  if (allyrs) { 
    B_BMSY <-sapply(1:nMPs, function(mm) 
      apply(MSEobj@B_BMSY[,mm,yrs] > BmsyRef, 1, prod, na.rm=TRUE))
    B_B0 <-sapply(1:nMPs, function(mm) 
      apply(MSEobj@B_BMSY[,mm,yrs] * MSEobj@OM$BMSY_B0 > B0Ref, 1, prod,na.rm=T))
    F_FMSY <-sapply(1:nMPs, function(mm) 	
      apply(MSEobj@F_FMSY[,mm,yrs] < FRef, 1, prod, na.rm=TRUE))
    # Interannual variability in yield
    MSEobj@C[(!is.finite(MSEobj@C[,,]))] <- 0 # if catch is NAN or NA, make it 0
    aavy <- sapply(1:nMPs, function(mm)
      apply((((MSEobj@C[,mm,y1]-MSEobj@C[,mm,y2])/MSEobj@C[,mm,y2])^2)^0.5,1,mean,na.rm=T)) 
    aavy[is.nan(aavy)] <- 1
    AAVY	<- aavy < (maxVar/100)
    # Interannual variability in effort (actually fishing mortality here) - need to change MSE object
    aave <- sapply(1:nMPs, function(mm)
      apply((((MSEobj@FM[,mm,y1]-MSEobj@FM[,mm,y2])/MSEobj@FM[,mm,y2])^2)^0.5,1,mean,na.rm=T)) 
    aave[is.nan(aave)] <- 1
    AAVE	<- aave < (maxVarE/100)
    
    # Maximizing - average over simulations and years 
    LTY <-sapply(1:nMPs, function(mm) 
      round(mean(apply(MSEobj@C[,mm,yrs],1,mean,na.rm=TRUE)/RefYd)*100,1))
    STY <-sapply(1:nMPs, function(mm) 
     round(mean(apply(MSEobj@C[,mm,ystart],1,mean,na.rm=TRUE)/RefYd)*100,1))
    
    
    Stats <- list(B_BMSY=B_BMSY, B_B0=B_B0, F_FMSY=F_FMSY, AAVY=AAVY, AAVE=AAVE)  
    
    Ind <- match(myPMdf[,1], names(Stats)) # Specified Joint Performance Metrics
    
    Legend <- NULL
    for (X in 1:length(Ind)) Legend <- append(Legend, as.expression(Names[[Ind[X]]]))
    
    # individual probabilities
    inProbs <- lapply(lapply(Stats[Ind], colSums), "/" , nsim) 
    indivProbs <- matrix(unlist(inProbs), ncol = nMPs, byrow = TRUE) * 100 
    colnames(indivProbs) <- MPs 
    indivProbs <- indivProbs[,order(indivProbs[1,])]
    
    # Joint probabilities
    # does each sim meet the joint performance criteria?
    total <- Reduce("*", Stats[Ind]) 
    prob <- round(colSums(total)/nrow(total) * 100, 2) # probability of meeting PMs
    
    jointDF <- data.frame(MPs=MPs, Prob=prob, LTY=LTY, STY=STY)
  }
  
  out <- NULL 
  out$Leg <- Legend 
  if (allyrs) out$iprobs <- indivProbs
  if (!allyrs) out$iprobs <- indivProbs2
  out$jprobs <- jointDF
  out$MPs <- MPs 
  out$LTY <- LTY 
  out$STY <- STY
  out$lastYrs <- lastYrs
  return(out) 
}

bplot2 <- function(perfobj, PLim=0.5, bcol=NULL, maxmp=20, lastYrs=NULL) {
  iprobs <- perfobj$iprobs 
  Legend <- perfobj$Leg
  MPs <- perfobj$MPs
  nMPs <- length(MPs)
  if (is.null(lastYrs)) lastYrs <- perfobj$lastYrs 
  if (PLim < 1) PLim <- PLim * 100
  linet <- 2 
  Col <- "gray"
  if (is.null(bcol)) bcol <- gray.colors(length(Legend))
  if (length(bcol) < length(Legend)) warning("Colors are not same length as PMs")
  if (length(bcol) > length(Legend)) bcol <- bcol[1:length(Legend)]
  
  
  nplots <- ceiling(nMPs / maxmp)
  Ncol <- ceiling(sqrt(nplots))
  Nrow <- ceiling(nplots/Ncol)
  if (!(Ncol * Nrow >= nplots)) stop()
  
  tempmat <- matrix(1:(Ncol*Nrow), nrow=Nrow, byrow=TRUE)
  bspace <- max(nchar(MPs)) - 1 
  
  par(mfrow=c(Nrow, Ncol), oma=c(2.5,1,5,2), mar=c(3,bspace,0,0))
 
  if (nplots == 1) {
    tt <- barplot(iprobs, beside=TRUE, ylab="Probability",
  	  las=3, xlim=c(0,100),cex.axis=1.5, cex.lab=2, las=2, 
	  cex.names=1.5,xpd=NA, col=bcol, names=names(iprobs), horiz=TRUE)
	
	mtext(side=3, outer=TRUE, paste("Last", lastYrs, 
	  "years of projection period"), line=1, cex=1.25)	
    abline(h=PLim, lty=linet, col=Col)		
	rng <- par("usr")
    lg <- legend(rng[1], rng[2], legend=Legend,bty="n", 
	  cex=1.25, horiz=TRUE, xpd=NA, fill=bcol, title="Limit Reference Points", 
	  plot=FALSE)
    legend(rng[1], rng[4]+lg$rect$h, legend=Legend,bty="n", 
	  cex=1.25, horiz=TRUE, xpd=NA, fill=bcol, title="Limit Reference Points")
  } else {
	npplot <- min(ceiling(nMPs/nplots), maxmp)
    xx <- 1 
    xx2 <- npplot 
    for (X in 1:nplots) {
	  if (length(Legend) > 1) splitdat <- iprobs[,xx:xx2] 
	  if (length(Legend) == 1) splitdat <- iprobs[xx:xx2]
	  if (X == 1) {
    	tt <- barplot(splitdat, beside=TRUE, ylab="",
    	xlim=c(0,100), col=bcol, names=names(splitdat),
    	cex.axis=1.5, cex.lab=2, las=1, cex.names=1.5, xpd=NA, horiz=TRUE)
		abline(v=PLim, lty=linet, col=Col, xpd=FALSE)
		rng <- par("usr")
        lg <- legend(rng[1], rng[2], legend=Legend,bty="n", 
	      cex=1.25, horiz=TRUE, xpd=NA, fill=bcol, title="Limit Reference Points", 
	      plot=FALSE)
        legend(rng[1], rng[4]+lg$rect$h, legend=Legend,bty="n", 
	    cex=1.25, horiz=TRUE, xpd=NA, fill=bcol, title="Limit Reference Points")
	  } else {
    	barplot(splitdat, beside=TRUE, ylab="", xlim=c(0,100), 
		names=names(splitdat), cex.axis=1.5, cex.lab=2, las=1, cex.names=1.5, 
		axes=TRUE, col=bcol, horiz=TRUE)
		# if (X %in% tempmat[,1]) 
		  # axis(side=2, label=TRUE, cex.axis=1.5, cex.lab=2, las=2)
		# if (!X %in% tempmat[,1]) axis(side=2, label=FALSE)
		abline(v=PLim, lty=linet, col=Col, xpd=FALSE)
      }		
      xx <- xx2 + 1 
      xx2 <- min(xx + npplot -1, nMPs)
	  mtext(side=1, outer=TRUE, "Probability", line=1, cex=1.4)
	  mtext(side=3, outer=TRUE, paste("Last", lastYrs, "years of projection period"),
	  line=3.5, cex=1.25)
    }
  }

  
  if (!is.null(nrow(iprobs))) {
    rownames(iprobs) <- Legend
  }	
  perfobj$iprobs <- iprobs 
  perfobj$ipass <- iprobs > PLim
  perfobj
}

bplot <- function(perfobj, PLim=0.5, bcol=NULL, maxmp=10, lastYrs=NULL) {
  iprobs <- perfobj$iprobs 
  Legend <- perfobj$Leg
  MPs <- perfobj$MPs
  nMPs <- length(MPs)
  if (is.null(lastYrs)) lastYrs <- perfobj$lastYrs 
  if (PLim < 1) PLim <- PLim * 100
  linet <- 2 
  Col <- "gray"
  if (is.null(bcol)) bcol <- gray.colors(length(Legend))
  if (length(bcol) < length(Legend)) warning("Colors are not same length as PMs")
  if (length(bcol) > length(Legend)) bcol <- bcol[1:length(Legend)]
  
  nplots <- ceiling(nMPs / maxmp)
  Ncol <- ceiling(sqrt(nplots))
  Nrow <- ceiling(nplots/Ncol)
  if (!(Ncol * Nrow >= nplots)) stop()
  
  tempmat <- matrix(1:(Ncol*Nrow), nrow=Nrow, byrow=TRUE)
  bspace <- max(nchar(MPs)) - 1 
  
  par(mfrow=c(Nrow, Ncol), oma=c(0,5,5,0), mar=c(bspace,1,0,1))
 
  if (nplots == 1) {
    tt <- barplot(iprobs, beside=TRUE, ylab="Probability",
  	  las=3, ylim=c(0,100),cex.axis=1.5, cex.lab=2, las=2, 
	  cex.names=1.5,xpd=NA, col=bcol, names=names(iprobs))
	
	mtext(side=3, outer=TRUE, paste("Last", lastYrs, 
	  "years of projection period"), line=1, cex=1.25)	
    abline(h=PLim, lty=linet, col=Col)		
	rng <- par("usr")
    lg <- legend(rng[1], rng[2], legend=Legend,bty="n", 
	  cex=1.25, horiz=TRUE, xpd=NA, fill=bcol, title="Limit Reference Points", 
	  plot=FALSE)
    legend(rng[1], rng[4]+lg$rect$h, legend=Legend,bty="n", 
	  cex=1.25, horiz=TRUE, xpd=NA, fill=bcol, title="Limit Reference Points")
  } else {
	npplot <- min(ceiling(nMPs/nplots), maxmp)
    xx <- 1 
    xx2 <- npplot 
    for (X in 1:nplots) {
	  if (length(Legend) > 1) splitdat <- iprobs[,xx:xx2] 
	  if (length(Legend) == 1) splitdat <- iprobs[xx:xx2]
	  if (X == 1) {
    	tt <- barplot(splitdat, beside=TRUE, ylab="",
    	las=3, ylim=c(0,100), col=bcol, names=names(splitdat),
    	cex.axis=1.5, cex.lab=2, las=2, cex.names=1.5, xpd=NA)
		abline(h=PLim, lty=linet, col=Col, xpd=FALSE)
		rng <- par("usr")
        lg <- legend(rng[1], rng[2], legend=Legend,bty="n", 
	      cex=1.25, horiz=TRUE, xpd=NA, fill=bcol, title="Limit Reference Points", 
	      plot=FALSE)
        legend(rng[1], rng[4]+lg$rect$h, legend=Legend,bty="n", 
	    cex=1.25, horiz=TRUE, xpd=NA, fill=bcol, title="Limit Reference Points")
	  } else {
    	barplot(splitdat, beside=TRUE, ylab="", las=3, ylim=c(0,100), 
		names=names(splitdat), cex.axis=1.5, cex.lab=2, las=2, cex.names=1.5, 
		axes=FALSE, col=bcol)
		if (X %in% tempmat[,1]) 
		  axis(side=2, label=TRUE, cex.axis=1.5, cex.lab=2, las=2)
		if (!X %in% tempmat[,1]) axis(side=2, label=FALSE)
		  abline(h=PLim, lty=linet, col=Col, xpd=FALSE)
      }		
      xx <- xx2 + 1 
      xx2 <- min(xx + npplot -1, nMPs)
	  mtext(side=2, outer=TRUE, "Probability", line=3, cex=1.4)
	  mtext(side=3, outer=TRUE, paste("Last", lastYrs, "years of projection period"),
	  line=3, cex=1.25)
    }
  }

  
  if (!is.null(nrow(iprobs))) {
    rownames(iprobs) <- Legend
  }	
  perfobj$iprobs <- iprobs 
  perfobj$ipass <- iprobs > PLim
  perfobj
}

jplot <- function(perfobj, Yield=c("LTY","STY"), PLim=0.5, lastYrs=10,
  ShowCols=TRUE, ShowLabs=TRUE, ShowPMs=TRUE, AvailMPs=NULL, xlim=NULL, 
  labcol=NULL, lcex=0.9, onlyAvail=FALSE, doPlot=TRUE) {
 
  jointDF <- perfobj$jprobs
  MPs <- as.character(jointDF[,1])
  nMPs <- length(MPs)
  prob <- jointDF[,2]
  LTY <- jointDF[,3]
  STY <- jointDF[,4]
  Leg <- perfobj$Leg
  
  if (PLim < 1) PLim <- PLim * 100 
  if (!is.null(xlim)) if (max(xlim) <= 1) xlim <- xlim * 100 
  
  Yield <- match.arg(Yield)
  yield <- switch(Yield, LTY=LTY, STY=STY)
  ylb <- switch(Yield, LTY="Long-Term", STY="Short-Term")
  ylb2 <- switch(Yield, LTY=paste(" (last", lastYrs, "years)"), 
    STY=paste(" (first", lastYrs, "years)"))
  # trade-off plot of yield and prob of meeting PMs 
  x <- prob
  y <- yield
  labs <- MPs 
  xlab <- "Joint probability of exceeding Limit Reference Points"
  ylab <- paste0("Average ", ylb, " Yield", ylb2)
  par(mfrow=c(1,1), mar=c(5,5,2,1), oma=c(1,1,1,3))
  adjj<-c(0.9,1.1)
  XLim <- c(min(c(-10, min(x,na.rm=T)*adjj)), max(c(max(x,na.rm=T)*adjj, 110)))
  if (!is.null(xlim)) XLim <- c(min(xlim)*0.975, max(xlim*1.025))
  YLim <- c(min(c(-10, min(y,na.rm=T)*adjj)), max(c(max(y,na.rm=T)*adjj, 110)))
  # Which MPs meet minimum PMs 
  ind <- which(x >= PLim)
  coly <- rep("darkgray", length(labs)) 
  coly[ind] <- "black"   
  Pch <- rep(21, length(labs))
  coly[grep("FMSY", labs)]<-'darkgray'
  Pch[grep("FMSY", labs)] <- 24
  if (!is.null(AvailMPs)) Pch[labs%in%AvailMPs] <- 21
  if (!is.null(AvailMPs)) coly[labs%in%AvailMPs & (x >= PLim)] <- "green"
  pmtype <- sapply(1:length(coly), function(X) class(get(labs[X])))
  if (doPlot) {
    plot(NA,xlim=XLim,ylim=YLim,xlab=xlab,ylab=ylab, bty="l", las=1, 
    	xaxs="i", yaxs="i", cex.lab=1.5, cex.axis=1.25)
    if (PLim > 0) abline(v=PLim,col="#99999940",lwd=2)
    Alpha <- 30
    # polygons 
    LeftCol <- rgb(red=255, green=0, blue=0, alpha=Alpha, names = NULL, 
      maxColorValue = 255)
    RightCol <- rgb(red=0, green=255, blue=0, alpha=Alpha, names = NULL, 
      maxColorValue = 255) 
    if(ShowCols) {
      polygon(x=c(max(0, XLim[1]), PLim,  PLim, max(0, XLim[1])),
    	  y=c(0, 0, max(YLim), max(YLim)), col=LeftCol, border=NA)
      polygon(x=c(PLim,  100, 100, PLim), y=c(0, 0, max(YLim), max(YLim)), 
	    col=RightCol, border=NA)
    }
    if(ShowPMs){
      temp <- 1 
      text(quantile(XLim,0.05), max(YLim)*1.05, "Limit Reference Points", 
	    xpd=NA, pos=4, cex=1.2)
      for (xx in 1:length(Leg)) {
        temp <- temp - 0.05
        text(quantile(XLim,0.05), max(YLim*temp)*1.05, Leg[xx], xpd=NA, pos=4)  
      }
    }
    
    Cex <- 1.5
    if(!ShowLabs) {
	  points(x,y, bg=coly, pch=Pch, cex=Cex, col="black")
	  if (PLim <= 0) {
	    legend(-5, 0.92*max(YLim), legend=c("MPs", "Reference"), 
	      col="black", bty="n", pch=c(21, 24), pt.bg=c("black", "darkgray"), 
		  pt.cex=Cex, xpd=NA)
	  } else {
	    legend(-5, 0.92*max(YLim), legend=c("Accept", "Reject", "Reference"), 
	      col="black", bty="n", pch=c(21, 21, 24), pt.bg=c("black", "darkgray", "darkgray"), 
		  pt.cex=Cex, xpd=NA)
      }	  
	}  
    if (is.null(labcol)) coly2 <- coly 
    if (!is.null(labcol)) coly2 <- rep(labcol, length.out=length(coly))[1:length(coly)]
    
    labcex <- rep(lcex, nMPs)
    if (!is.null(AvailMPs)) {
      labcex[labs%in%AvailMPs & (x >= PLim)] <- lcex * 1.5
      coly2[!(labs%in%AvailMPs) | (x <= PLim)] <- "#d3d3d3" 
    }	
    
    coly2[(x <= PLim)] <- "#d3d3d3"
    coly2[grep("FMSY", labs)] <- '#000000'
    labcex[grep("FMSY", labs)] <- lcex * 0.9 
    font <- rep(2, length(coly)) # make all bold 
    font[pmtype == "DLM_input"] <- 4
    font[grep("FMSY", labs)] <- 1 
    if(ShowLabs) {
      if (onlyAvail & !is.null(AvailMPs)) {
	    text(x[labs%in%AvailMPs & (x >= PLim)],
	    y[labs%in%AvailMPs & (x >= PLim)],
	    labs[labs%in%AvailMPs & (x >= PLim)],
	    font=font[labs%in%AvailMPs & (x >= PLim)],
	    col=coly2[labs%in%AvailMPs & (x >= PLim)],
	    cex=labcex[labs%in%AvailMPs & (x >= PLim)])
		temp <- length(grep("FMSY", labs[labs%in%AvailMPs & (x >= PLim)]))
		if (!is.finite(temp)) temp <- 0 
		Num <- table(pmtype[labs%in%AvailMPs & (x >= PLim)])
		intext <- paste0("Input (n=", Num[1],")")
		outtext <- paste0("Output (n=", Num[2]-temp,")")
		legend(99, max(YLim), legend=c(intext, outtext), text.font=c(4,2), xpd=NA, 
	    bty="n", title="Method Type")
	  } else {
        text(x,y,labs,font=font,col=coly2,cex=labcex)
		temp <- length(grep("FMSY", labs[x >= PLim]))
		if (!is.finite(temp)) temp <- 0 
		Num <- table(pmtype[x > PLim])
		intext <- paste0("Input (n=", Num[1],")")
		outtext <- paste0("Output (n=", Num[2]-temp,")")
		legend(99, max(YLim), legend=c(intext, outtext), text.font=c(4,2), xpd=NA, 
	    bty="n", title="Method Type")
	  }
	  
    }	
  } 
  mat <- cbind(prob, switch(Yield, LTY=LTY, STY=STY))
  rownames(mat) <- MPs
  colnames(mat) <- c("Prob", "Yield")
  
  Avail <- rep(NA, nMPs) 
  if(!is.null(AvailMPs)) { 
    Avail <- rep(FALSE, nMPs)
    Avail[labs%in%AvailMPs] <- TRUE
  }	
  Dist <- NULL # calculate distance from corner
  for (X in 1:nrow(mat)) Dist[X] <- euc.dist(c(mat[X,1], mat[X,2]), 
    c(100, max(mat[,2])))
  out <- data.frame(Prob=mat[,1], Yield=mat[,2], Name=rownames(mat), Dist=Dist, 
    Type=pmtype, Avail=Avail, stringsAsFactors=FALSE)
  refind <- grep("FMSY", out$Name)
  if (length(refind) > 0) {
    out <- out[-refind,] # remove the reference methods from the list of MPs 
  }
  out <- out[order(out$Dist),]
  # out <- data.frame(Prob=mat[,1], Yield=mat[,2], Name=rownames(mat), 
    # stringsAsFactors=FALSE)
	
  perfobj$pass <- out[out[,1] >= PLim,]
  perfobj$perf <- out
  perfobj
  
  
}

dplot <- function(MSEobj, perfobj=NULL, MPs=NULL, outline=FALSE,
   maxN=3, bytype=TRUE, cols=TRUE, onlyAvail=TRUE, 
   refs=c(0.3, 1, 1, 0, 0), incline=TRUE, PLim=0.8) {
  
  if (is.null(perfobj) | is.null(perfobj$pass)) stop("Run jplot first")
  if (onlyAvail) Avail <- perfobj$pass$Avail
  if (!onlyAvail) Avail <- rep(TRUE, nrow(perfobj$pass))
  if(!is.null(perfobj) && bytype) { 
    output <- which(perfobj$pass$Type == "DLM_output" & Avail)[1:maxN]
	input <- which(perfobj$pass$Type == "DLM_input" & Avail)[1:maxN]
	mps <- c(perfobj$pass$Name[input], perfobj$pass$Name[output])
	mps <- mps[!is.na(mps)]
    mseobj <- Sub(MSEobj, MPs=mps)	
	lastYrs <- perfobj$lastYrs
	
  }
  if (is.null(MPs) & (!bytype)) {
    mps <- perfobj$pass$Name[perfobj$pass$Avail][1:maxN]
    mseobj <- Sub(MSEobj, MPs=mps)	
  }	
  if (!exists("lastYrs")) lastYrs <- 10 
  if (!is.null(MPs)) mseobj <- Sub(MSEobj, MPs=MPs)	
  if (is.null(MPs) & is.null(perfobj)) mseobj <- Sub(MSEobj, MPs=1:maxN)	
  
  Nyears <- mseobj@nyears
  Pyears <- mseobj@proyears
  nMPs <- mseobj@nMPs
  MPs <- mseobj@MPs 
  nsim <- mseobj@nsim
  RefYd <- mseobj@OM$RefY
  
  yrs <- (Pyears-lastYrs+1):Pyears 	# Years to summarize performance
  ystart <- 1:lastYrs			 	# First 10 years
  y1 <- yrs[1]:yrs[length(yrs)-1] 	# for calculating interannual variability
  y2 <- y1+1
  
  if (lastYrs >= Pyears) {
    message("lastYrs set too high. Defaulting to last 10 years")
    lastYrs <- 10
  }
  if (lastYrs <= 0 | lastYrs == FALSE | lastYrs == "all" | lastYrs == "All") 
    lastYrs <- Pyears 
 
  # Reporting Distribution in last years 
  BBMSY <- mseobj@B_BMSY[,,yrs]
  BBMSY <- aperm(BBMSY, c(1,3,2))
  dim(BBMSY) <- c(nsim*lastYrs,nMPs)
  
  BB0 <- mseobj@B_BMSY[,,yrs] * mseobj@OM$BMSY_B0
  BB0 <- aperm(BB0, c(1,3,2))
  dim(BB0) <- c(nsim*lastYrs,nMPs) 

  FFMSY <- mseobj@F_FMSY[,,yrs] 
  FFMSY <- aperm(FFMSY, c(1,3,2))
  dim(FFMSY) <- c(nsim*lastYrs,nMPs) 
  
  LTY <- mseobj@C[,,yrs]/RefYd
  LTY <- aperm(LTY, c(1,3,2))
  dim(LTY) <- c(nsim*lastYrs,nMPs)

  aavy <- sapply(1:nMPs, function(mm)
    apply((((mseobj@C[,mm,y1]-mseobj@C[,mm,y2])/mseobj@C[,mm,y2])^2)^0.5,1,mean,na.rm=T)) 
  # aavy[is.nan(aavy)] <- 1
  # aavy[aavy > 5] <- 5
  #Interannual variability in effort (actually fishing mortality here) - need to change MSE object  
  aave <- sapply(1:nMPs, function(mm)
    apply((((mseobj@FM[,mm,y1]-mseobj@FM[,mm,y2])/mseobj@FM[,mm,y2])^2)^0.5,1,mean,na.rm=T)) 
  # aave[is.nan(aave)] <- 1
  # aave[aave > 5] <- 5 
 
  # Colors - to make the plot a bit more cheerful 
  pmtype <- sapply(1:nMPs, function(X) class(get(MPs[X])))
  font1 <- 2 
  font2 <- 4
  inputs <- which(pmtype=="DLM_input")
  outputs <- which(pmtype=="DLM_output")
  if (!(cols)) Col <- "lightblue"
  if (cols) {
    cols1 <- brewer.pal(8, "Set1")
	cols2 <-  brewer.pal(8, "Set2")
    Col <- rep(cols1, length.out=nMPs)[1:nMPs]
	if (length(inputs) > 0) Col[inputs] <- cols2[1:length(inputs)]
  }

  Line <- 3.5
  Cex <- 1.05 
  CexName <- 1.25
  lncol <- "slategray"
  topcex <- 0.75
  ncharmp <- max(nchar(MPs)) - 2
  par(mfrow=c(2,3), bty="n", oma=c(4,ncharmp,0,0), mar=c(3,2,4,3))
  # B/B0 
  boxplot(BB0, names=FALSE, las=1, outline=outline, horizontal=TRUE,
      xlab="", bty="n", col=Col, xpd=NA, cex.names=CexName, axes=FALSE)
  axis(side=1, cex.axis=CexName)
  axis(side=2, at=outputs, labels=MPs[outputs], cex.axis=CexName, font=font1, las=2)	  
  axis(side=2, at=inputs, labels=MPs[inputs], cex.axis=CexName, font=font2, las=2)
  
  mtext(side=1, expression(italic(B/B[0])), cex=Cex, line=Line)
  if(incline && refs[1] > 0) {
    ps <- apply(BB0 > refs[1], 2, sum, na.rm=TRUE)/(nrow(BB0))
	ps <- sprintf("%3.2f", ps)
	ps <- c(ps, "Prob")
    axis(side=4, at=1:(nMPs+1), labels=ps, las=1, tick=FALSE, line=-1, xpd=NA)
    abline(v=refs[1], lty=2, lwd=2, col=lncol)
	mtext(side=3, bquote(italic(B) > ~ .(refs[1]) ~ italic(B[0])), cex=topcex)
  }
  # B/BMSY
  boxplot(BBMSY, names=FALSE, las=1, outline=outline, horizontal=TRUE,
      xlab="", bty="n", col=Col, xpd=NA, axes=FALSE)
  mtext(side=1, expression(italic(B/B[MSY])), cex=Cex, line=Line)  
  axis(side=1, cex.axis=CexName)
  axis(side=2, at=outputs, labels=FALSE, cex.axis=CexName, font=font1, las=2)	  
  axis(side=2, at=inputs, labels=FALSE, cex.axis=CexName, font=font2, las=2)
  if(incline && refs[2] > 0) { 
    ps <- apply(BBMSY > refs[2], 2, sum, na.rm=TRUE)/(nrow(BBMSY))
	ps <- sprintf("%3.2f", ps)
	ps <- c(ps, "Prob")
    axis(side=4, at=1:(nMPs+1), labels=ps, las=1, tick=FALSE, line=-1, xpd=NA)
    abline(v=refs[2], lty=2, lwd=2, col=lncol)
	mtext(side=3, bquote(italic(B) > ~ .(refs[2]) ~ italic(B[MSY])), cex=topcex)
  }
  # F/FMSY 
  boxplot(FFMSY, names=FALSE, las=1, outline=outline, horizontal=TRUE,
      xlab="", bty="n", col=Col, xpd=NA, axes=FALSE)
  mtext(side=1, expression(italic(F/F[MSY])), cex=Cex, line=Line)
  axis(side=1, cex.axis=CexName)
  axis(side=2, at=outputs, labels=FALSE, cex.axis=CexName, font=font1, las=2)	  
  axis(side=2, at=inputs, labels=FALSE, cex.axis=CexName, font=font2, las=2)
  if(incline && refs[3] > 0) {
    ps <- apply(FFMSY < refs[3], 2, sum, na.rm=TRUE)/(nrow(FFMSY))
	ps <- sprintf("%3.2f", ps)
	ps <- c(ps, "Prob")
    axis(side=4, at=1:(nMPs+1), labels=ps, las=1, tick=FALSE, line=-1, xpd=NA)
    abline(v=refs[3], lty=2, lwd=2, col=lncol)
	mtext(side=3, bquote(italic(F) < ~ .(refs[3]) ~ italic(F[MSY])), cex=topcex)
  }
  # AAVY 
  if (refs[4] < 1) refs[4] <- refs[4] * 100 
  aavy <- aavy * 100 
  boxplot(aavy, names=FALSE, las=1, outline=outline, horizontal=TRUE,
      xlab="", bty="n", col=Col, xpd=NA, axes=FALSE)
  mtext(side=1, "Average Annual Variation Yield (%)", cex=Cex, line=Line)
  axis(side=1, cex.axis=CexName)
  axis(side=2, at=outputs, labels=MPs[outputs], cex.axis=CexName, font=font1, las=2)	  
  axis(side=2, at=inputs, labels=MPs[inputs], cex.axis=CexName, font=font2, las=2)  
  if(incline && refs[4] > 0) {
    ps <- apply(aavy < refs[4], 2, sum, na.rm=TRUE)/(nrow(aavy))
	ps <- sprintf("%3.2f", ps)
	ps <- c(ps, "Prob")
    axis(side=4, at=1:(nMPs+1), labels=ps, las=1, tick=FALSE, line=-1, xpd=NA)
    abline(v=refs[4], lty=2, lwd=2, col=lncol)
	mtext(side=3, bquote(italic(AAVY) < ~ .(refs[4]) ~ "%"), cex=topcex) 
  }
  abline(v=refs[4], lty=2, lwd=2, col=lncol)
 # # AAVE 
  # boxplot(aave, names=FALSE, las=1, outline=outline, horizontal=TRUE,
      # xlab="", bty="n", col=Col, xpd=NA)
  # mtext(side=1, "Average Annual Variation Effort", cex=1.25, line=Line)  
  
  # LTY 
  boxplot(LTY, names=FALSE, las=1, outline=outline, horizontal=TRUE,
      xlab="", bty="n", col=Col, xpd=NA, axes=FALSE)
  mtext(side=1, "Relative Long-Term Yield", cex=Cex, line=Line)  
  axis(side=1, cex.axis=CexName)
  axis(side=2, at=outputs, labels=FALSE, cex.axis=CexName, font=font1, las=2)	  
  axis(side=2, at=inputs, labels=FALSE, cex.axis=CexName, font=font2, las=2)
  if(incline && refs[5] > 0) abline(v=refs[4], lty=2, lwd=2, col="darkgray")
}

PerfPlot <- function(MSEobj, 
  myPMs=c("B_BMSY", "B_B0", "F_FMSY", "AAVY", "AAVE"), 
  myrefs=c(0.5, 0.2, 1.5, 30, 30), Yield=c("LTY","STY"), PLim=0.5, lastYrs=10,
  ShowCols=TRUE, ShowLabs=TRUE, ShowPMs=TRUE, bcol=NULL,
  AvailMPs=NULL, bplot=TRUE, tplot=TRUE, maxmp=12, txlim=NULL) {
  
  # Assign Variables
  Nyears <- MSEobj@nyears
  Pyears <- MSEobj@proyears
  nMPs <- MSEobj@nMPs
  MPs <- MSEobj@MPs 
  nsim <- MSEobj@nsim
  RefYd <- MSEobj@OM$RefY
  Yield <- match.arg(Yield)
  
  # Error Checks 
  if (lastYrs >= Pyears) {
    message("lastYrs set too high. Defaulting to all years")
    lastYrs <- 0 
  }
  if (lastYrs <= 0 | lastYrs == FALSE | lastYrs == "all" | lastYrs == "All") 
    lastYrs <- Pyears - 1  
  if (length(myrefs) != length(myPMs)) stop("Length of myPMs and myrefs must be equal")
  if (PLim < 1) PLim <- PLim * 100
  # Performance metrics 
  allPMs <- c("B_BMSY", "B_B0", "F_FMSY", "AAVY", "AAVE", "AAVB")
  dfRefs <- c(0.5, 0.2, 1.25, 20, 20, 20) # Default reference points
  PMdf <- data.frame(allPMs, dfRefs) # Default performance metrics 
  
  myPMs <- match.arg(myPMs, allPMs, several.ok=TRUE)
  myPMdf <- PMdf[match(myPMs, PMdf[,1]),]
  myPMdf[,2] <- myrefs
  
  BmsyRef <- myPMdf[grep("B_BMSY", myPMdf[,1]),2]
  B0Ref <- myPMdf[grep("B_B0", myPMdf[,1]),2]
  FRef <- myPMdf[grep("F_FMSY", myPMdf[,1]),2]
  maxVar <- myPMdf[grep("AAVY", myPMdf[,1]),2]
  maxVarE <- myPMdf[grep("AAVE", myPMdf[,1]),2]
  maxVarB <- myPMdf[grep("AAVB", myPMdf[,1]),2]
  
  if(length(BmsyRef) < 1) BmsyRef <- PMdf[grep("B_BMSY", PMdf[,1]),2]
  if(length(B0Ref) < 1) B0Ref <- PMdf[grep("B_B0", PMdf[,1]),2]
  if(length(FRef) < 1) FRef <- PMdf[grep("F_FMSY", PMdf[,1]),2]
  if(length(maxVar) < 1) maxVar <- PMdf[grep("AAVY", PMdf[,1]),2]
  if(length(maxVarE) < 1) maxVarE <- PMdf[grep("AAVE", PMdf[,1]),2]
  if(length(maxVarB) < 1) maxVarB <- PMdf[grep("AAVB", PMdf[,1]),2]
  
  Names <- list(6)
  Names[[1]] <- bquote(italic(B) > ~ .(BmsyRef) ~ italic(B[MSY]))
  Names[[2]] <- bquote(italic(B) > ~ .(B0Ref) ~ italic(B[0]))
  Names[[3]] <- bquote(italic(F) < ~ .(FRef) ~ italic(F[MSY]))
  Names[[4]] <- bquote(AAVY < ~ .(maxVar)* '%')
  Names[[5]] <- bquote(AAVE < ~ .(maxVarE)* '%')
  Names[[6]] <- bquote(AAVB < ~ .(maxVarB)* '%')
    
  if (length(BmsyRef)> 1 | length(B0Ref)> 1 | length(B0Ref)> 1
    | length(maxVar) > 1 | length(maxVarE) > 1) 
	stop("Repeated performance metrics")
  if (maxVar < 1) maxVar <- maxVar * 100 # maximum variability in % 
  if (maxVarE < 1) maxVarE <- maxVarE * 100 # maximum variability in %
  
  yrs <- (Pyears-lastYrs+1):Pyears 	# Years to summarize performance
  ystart <- 1:lastYrs			 	# First 10 years
  y1 <- yrs[1]:yrs[length(yrs)-1] 	# for calculating interannual variability
  y2 <- y1+1
 
  # Satisficing matrices 
  B_BMSY  <- array(NA, dim=c(nsim, nMPs))
  B_B0    <- array(NA, dim=c(nsim, nMPs))
  F_FMSY  <- array(NA, dim=c(nsim, nMPs))
  AAVY    <- array(NA, dim=c(nsim, nMPs))
  AAVE    <- array(NA, dim=c(nsim, nMPs))
  VarYE   <- array(NA, dim=c(nsim, nMPs)) # variability in catch or effort
    
  # Maximizing & Reporting 
  LTY <- rep(NA, nMPs)
  STY <- rep(NA, nMPs)
  BBMSY <- array(NA, dim=c(nsim, nMPs))
  BB0 <- array(NA, dim=c(nsim, nMPs))
  FFMSY <- array(NA, dim=c(nsim, nMPs))
 
  # Calculate Stats  
  print("Calculating probability of meeting performance metrics")
  for (mm in 1:nMPs) {
    cat(".")
	flush.console()
    # Statisicing - proportion of simulations above threshold 
	B_BMSY[,mm] 	<- apply(MSEobj@B_BMSY[,mm,yrs] > BmsyRef, 1, prod, na.rm=TRUE) 
    B_B0[,mm] 		<- apply((MSEobj@B_BMSY[,mm,yrs] * MSEobj@OM$BMSY_B0 > B0Ref), 1, prod,na.rm=T)
    F_FMSY[,mm] 	<- apply(MSEobj@F_FMSY[,mm,yrs] < FRef, 1, prod, na.rm=TRUE)
	
	# Interannual variability in yield
	MSEobj@C[(!is.finite(MSEobj@C[,,]))] <- 0 # if catch is NAN or NA, make it 0
	
    aavy 			<- apply((((MSEobj@C[,mm,y1]-MSEobj@C[,mm,y2])/MSEobj@C[,mm,y2])^2)^0.5,1,mean,na.rm=T) 
    aavy[is.nan(aavy)] <- 1
    AAVY[,mm]		<- as.numeric(aavy<(maxVar/100))
	# Interannual variability in effort (actually fishing mortality here) - need to change MSE object
	aave 			<- apply((((MSEobj@FM[,mm,y1]-MSEobj@FM[,mm,y2])/MSEobj@FM[,mm,y2])^2)^0.5,1,mean,na.rm=T) 
    aave[is.nan(aave)] <- 1
    AAVE[,mm]		<- as.numeric(aave<(maxVarE/100))
	
	if (class(get(MSEobj@MPs[mm])) == "DLM_input") {
	  aave 	<- apply((((MSEobj@FM[,mm,y1]-MSEobj@FM[,mm,y2])/MSEobj@FM[,mm,y2])^2)^0.5,1,mean,na.rm=T)
	  VarYE[,mm]  <- as.numeric(aave<(maxVarB/100))
	} else {
	  aavy 	<- apply((((MSEobj@C[,mm,y1]-MSEobj@C[,mm,y2])/MSEobj@C[,mm,y2])^2)^0.5,1,mean,na.rm=T)
	  aave[is.nan(aave)] <- 1
	  aavy[is.nan(aavy)] <- 1
	  VarYE[,mm]  <- as.numeric(aavy<(maxVarB/100))
	}

	# Maximizing - average over simulations and years 
	LTY[mm] <- round(mean(apply(MSEobj@C[,mm,yrs],1,mean,na.rm=TRUE)/RefYd)*100,1)
	STY[mm] <- round(mean(apply(MSEobj@C[,mm,ystart],1,mean,na.rm=TRUE)/RefYd)*100,1)
	# Reporting Distribution in last years 
	BBMSY[,mm] <- apply(MSEobj@B_BMSY[,mm,yrs], 1, mean, na.rm=TRUE)
	BB0[,mm] <- apply(MSEobj@B_BMSY[,mm,yrs] * MSEobj@OM$BMSY_B0, 1, mean, na.rm=TRUE)
	FFMSY[,mm] <- apply(MSEobj@F_FMSY[,mm,yrs], 1, mean, na.rm=TRUE)	
  }
  cat("\n")
  Stats <- list(B_BMSY=B_BMSY, B_B0=B_B0, F_FMSY=F_FMSY, AAVY=AAVY, AAVE=AAVE, AAVB=VarYE)
  
  Names <- list(6)
  Names[[1]] <- bquote(italic(B) > ~ .(BmsyRef) ~ italic(B[MSY]))
  Names[[2]] <- bquote(italic(B) > ~ .(B0Ref) ~ italic(B[0]))
  Names[[3]] <- bquote(italic(F) < ~ .(FRef) ~ italic(F[MSY]))
  Names[[4]] <- bquote(AAVY < ~ .(maxVar)* '%')
  Names[[5]] <- bquote(AAVE < ~ .(maxVarE)* '%')
  Names[[6]] <- bquote(AAVB < ~ .(maxVarB)* '%')
  
  Ind <- match(myPMdf[,1], names(Stats))
  Legend <- NULL
  for (X in 1:length(Ind)) Legend <- append(Legend, as.expression(Names[[Ind[X]]]))
	
  # individual probabilities
  inProbs <- lapply(lapply(Stats[Ind], colSums), "/" , nsim) 
  output <- matrix(unlist(inProbs), ncol = nMPs, byrow = TRUE) * 100 
  colnames(output) <- MPs 
  rownames(output) <- Names[Ind]
  plotout <- output[,order(output[1,])]
  linet <- 2 
  Col <- "gray"
  if (is.null(bcol)) bcol <- gray.colors(nrow(output))
  if (length(bcol) != length(Ind)) warning("Colors are not same length as PMs")
  if (bplot) { # show barplots of individual probabilities
	if (nMPs > maxmp) maxmp <- max(which(nMPs %% 1:maxmp==0))
	nplots <- ceiling(nMPs/maxmp)
    Ncol <- ceiling(sqrt(nplots))
    Nrow <- ceiling(nplots/Ncol)
    tempmat <- matrix(1:(Ncol*Nrow), nrow=Nrow, byrow=TRUE)	
    plotlist <- list()

    par(mfrow=c(Nrow, Ncol), oma=c(0,5,5,0), mar=c(10,1,1,1))
    if (nplots == 1) {
      tt <- barplot(plotout, beside=TRUE, ylab="Probability",
    	las=3, ylim=c(0,100),cex.axis=1.5, cex.lab=2, las=2, 
		cex.names=1.5,xpd=NA, col=bcol)
	  legend(x=tt[1], y=105, legend=Legend,bty="n", 
		  cex=1.25, horiz=TRUE, xpd=NA, fill=bcol)
	  mtext(side=3, outer=TRUE, paste("Last", lastYrs, 
	    "years of projection period"), line=1, cex=1.25)	
      abline(h=PLim, lty=linet, col=Col)		
    } else {
	  npplot <- min(ceiling(nMPs/nplots), maxmp)
      xx <- 1 
      xx2 <- npplot 
      for (X in 1:nplots) {
		if (length(Ind) > 1) splitdat <- plotout[,xx:xx2] 
		if (length(Ind) == 1) splitdat <- plotout[xx:xx2]
		if (X == 1) {
    	  tt <- barplot(splitdat, beside=TRUE, ylab="",
    	  las=3, ylim=c(0,100), col=bcol,
    	  cex.axis=1.5, cex.lab=2, las=2, cex.names=1.5, xpd=NA)
		  legend(x=tt[1], y=130, legend=Legend,bty="n", horiz=TRUE,
		  cex=1.25, xpd=NA, fill=bcol)
		  abline(h=PLim, lty=linet, col=Col)
		} else {
    	  barplot(splitdat, beside=TRUE, ylab="", las=3, ylim=c(0,100), 
		  cex.axis=1.5, cex.lab=2, las=2, cex.names=1.5, axes=FALSE, col=bcol)
		  if (X %in% tempmat[,1]) axis(side=2, label=TRUE, 
		    cex.axis=1.5, cex.lab=2, las=2)
		  if (!X %in% tempmat[,1]) axis(side=2, label=FALSE)
		  abline(h=PLim, lty=linet, col=Col)
        }		
        xx <- xx2 + 1 
    	xx2 <- min(xx + npplot -1, nMPs)
		mtext(side=2, outer=TRUE, "Probability", line=3, cex=1.4)
		mtext(side=3, outer=TRUE, paste("Last", lastYrs, "years of projection period"),
		line=3, cex=1.25)
      }
    }
  }
  # Overall performance relative to chosen PMs 
  
  total <- Reduce("*", Stats[Ind]) # does each sim meet the performance criteria?

  prob <- round(colSums(total)/nrow(total) * 100, 2) # probability of meeting PMs
  
  yield <- switch(Yield, LTY=LTY, STY=STY)
  ylb <- switch(Yield, LTY="Long-Term", STY="Short-Term")
  ylb2 <- switch(Yield, LTY=paste(" (last", lastYrs, "years)"), 
    STY=paste(" (first", lastYrs, "years)"))
	
  if (tplot) {
    # trade-off plot of yield and prob of meeting PMs 
    x <- prob
    y <- yield
    labs <- MPs 
    xlab <- "Joint probability of meeting Performance Metrics"
    ylab <- paste0("Average ", ylb, " Yield", ylb2)
    par(mfrow=c(1,1), mar=c(5,5,2,1), oma=c(1,1,1,0))
    adjj<-c(0.9,1.1)
    XLim <- c(min(c(-10, min(x,na.rm=T)*adjj)), max(c(max(x,na.rm=T)*adjj, 110)))
	if (!is.null(txlim)) XLim <- txlim
    YLim <- c(min(c(-10, min(y,na.rm=T)*adjj)), max(c(max(y,na.rm=T)*adjj, 110)))
    # Which MPs meet minimum PMs 
    ind <- which(x >= PLim)
    coly <- rep("darkgray", length(labs)) 
    coly[ind] <- "black"   
    Pch <- rep(21, length(labs))
    coly[grep("FMSY", labs)]<-'black'
    Pch[grep("FMSY", labs)] <- 24
    if (!is.null(AvailMPs)) Pch[labs%in%AvailMPs] <- 21
    if (!is.null(AvailMPs)) coly[labs%in%AvailMPs & (x >= PLim)] <- "green"
    
    plot(NA,xlim=XLim,ylim=YLim,xlab=xlab,ylab=ylab, bty="l", las=1, 
    	xaxs="i", yaxs="i", cex.lab=1.5, cex.axis=1.25)
    if (PLim > 0) abline(v=PLim,col="#99999940",lwd=2)
    Alpha <- 30
    # polygons 
    LeftCol <- rgb(red=255, green=0, blue=0, alpha=Alpha, names = NULL, 
      maxColorValue = 255)
    RightCol <- rgb(red=0, green=255, blue=0, alpha=Alpha, names = NULL, 
      maxColorValue = 255) 
    if(ShowCols) {
      polygon(x=c(max(0, XLim[1]), PLim,  PLim, max(0, XLim[1])), y=c(0, 0, 100, 100), col=LeftCol, border=NA)
      polygon(x=c(PLim,  100, 100, PLim), y=c(0, 0, 100, 100), col=RightCol, border=NA)
    }
 
    if(ShowPMs){
      temp <- 1 
      text(quantile(XLim,0.05), max(YLim)*1.05, "Performance Metrics", xpd=NA, pos=4, cex=1.2)
      for (xx in 1:length(Ind)) {
        temp <- temp - 0.05
        text(quantile(XLim,0.05), max(YLim*temp)*1.05, Legend[xx], xpd=NA, pos=4)  
      }
    }

    Cex <- 1.5
    if(!ShowLabs) points(x,y, bg=coly, pch=Pch, cex=Cex, col="black" )
    if(ShowLabs) text(x,y,labs,font=2,col=coly,cex=1)
  }
  
  mat <- cbind(prob, switch(Yield, LTY=LTY, STY=STY))
  rownames(mat) <- MPs
  colnames(mat) <- c("Prob", "Yield")
  
  Dist <- NULL # calculate distance from corner
  for (X in 1:nrow(mat)) Dist[X] <- euc.dist(c(mat[X,1], mat[X,2]), c(100, 100))
  mat <- mat[order(Dist),]
  out <- data.frame(Prob=mat[,1], Yield=mat[,2], Name=rownames(mat), 
    stringsAsFactors=FALSE)
  OutList <- list() 
  OutList$IndP <- t(plotout)
  OutList$JointP <- out 
  OutList$Pass <- out[out[,1] >= PLim,]
  OutList$MPs <- MPs 
  OutList$Means <- array(c(BBMSY=BBMSY, BB0=BB0, FFMSY=FFMSY), 
    dim=c(nsim, nMPs, 3))

  OutList 	
}
