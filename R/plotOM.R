

#' Plot the operating model (OM) object parameters 
#' 
#' A function that plots the parameters and resulting time series of an operating model.
#' 
#' @param x An object of class OM 
#' @param ...  Optional additional arguments passed to \code{plot}
#' @rdname plot-OM
#' @method plot OM 
#' @author T. Carruthers
#' @export 
plot.OM <-function(x, ...){  
    OM <- updateMSE(x) # update and add missing slots with default values
    out<-runMSE(OM,Hist=T)

    nsim<-OM@nsim
    nyears<-OM@nyears
    
    plotStock(OM)
    plotFleet(OM)
    plotObs(OM)
    plotImp(OM)
    
    
    # Time series
    yrlab<-OM@CurrentYr-((nyears-1):0)
    op <-par(mfrow=c(4,2),mai=c(0.7,0.7,0.05,0.05),omi=c(0.01,0.01,0.3,0.01))
    
    # SSB
    TSplot(yrlab,out$TSdata$SSB,xlab="Historical year",ylab="Spawning biomass")
    
    # Depletion
    TSplot(yrlab,out$TSdata$SSB/rep(out$MSYs$SSB0,each=nyears),xlab="Historical year",ylab="Stock depletion (SSB)")
    
    # Apical F
    FM<-out$SampPars$Find*out$SampPars$qs
    TSplot(yrlab,t(FM),xlab="Historical year",ylab="Fishing mortality rate (apical)")
    
    # Catches
    TSplot(yrlab,out$TSdata$Catch,xlab="Historical year",ylab="Annual catches")
    
    # Recruitment
    TSplot(yrlab,out$TSdata$Rec,xlab="Historical year",ylab="Recruitment")
    
    # SSB-Rec
    TSplot(x=out$TSdata$SSB[2:nyears,],y=out$TSdata$Rec[2:nyears,],xlab="Spawning biomass",ylab="Recruitment",mat=F,type='p')
    
    F_FMSY<-FM/out$MSYs$FMSY
    B_BMSY<-t(out$TSdata$SSB)/out$MSYs$SSBMSY
    
    TSKplot(B_BMSY,F_FMSY,yrlab)
    
    # Age vulnerability
    maxage<-dim(out$SampPars$V)[2]
    contour(x=yrlab,y=1:maxage,z=t(out$SampPars$V[1,,1:nyears]),levels=c(0.25,0.75),col='green',drawlabels=F,lwd=c(1,2))
    contour(x=yrlab,y=1:maxage,z=t(out$SampPars$V[2,,1:nyears]),levels=c(0.25,0.75),col='blue',drawlabels=F,add=T,lwd=c(1,2))
    contour(x=yrlab,y=1:maxage,z=t(out$SampPars$V[3,,1:nyears]),levels=c(0.25,0.75),col='grey45',drawlabels=F,add=T,lwd=c(1,2))
    legend('topright',legend=c(paste("Simulation",1:3)),text.col=c("green","blue","grey45"),bty='n')
    legend('topleft',legend="Age vulnerability (0.25, 0.75)",bty='n')
    
    mtext("Historical year", 1, line = 2.5, cex = 1)
    mtext("Age", 2, line = 2.3, cex = 1)
    
    mtext(paste0("Time series plots for operating model ",OM@Name),3,outer=T,line= 0.2,font=2)
    
    on.exit(par(op))       
}

TSplot<-function(x,y,xlab=NA,ylab=NA,zeroy=T,incx=T,incy=T,type='l',mat=T){
  
  cols<-rep(makeTransparent(c("grey35","blue","orange","green")),100)
  nsim<-ncol(y)
  vlarg<-1e20
  
  rx<-range(x)
  ry<-range(y)
  if(zeroy)ry<-range(0,ry)
  
  rx1<-rx+c(-vlarg,vlarg)
  ry1<-ry+c(0,vlarg)
  
  plot(rx,ry,axes=F,col="white",xlab="",ylab="")
  polygon(rx1[c(1,1,2,2)],ry1[c(1,2,2,1)],col='grey94',border="grey94")
  xl<-pretty(x)
  abline(v=xl,col='white')
  yl<-pretty(seq(ry[1],ry[2],length.out=12))
  abline(h=yl,col='white')
  
  if(mat){
    matplot(x,y,type=type,col=cols,xlab="",ylab="",add=T)
  }else{
    if(type=='p')for(i in 1:nsim)points(x[,i],y[,i],col=cols[i],pch=19)
    if(type=='l')for(i in 1:nsim)lines(x[,i],y[,i],col=cols[i])
  }
  
  axis(1,c(-vlarg,vlarg),c(-vlarg,vlarg))
  axis(2,c(-vlarg,vlarg),c(-vlarg,vlarg))
  
  if(incx)axis(1,xl,xl)
  if(incy)axis(2,yl,yl)
  
  if(!is.na(xlab))mtext(xlab,1,line=2.5)
  if(!is.na(ylab))mtext(ylab,2,line=2.5)
  
}




TSKplot<-function(B_BMSY,F_FMSY,yrlab,maxsim=10){
  
  nyears<-ncol(B_BMSY)
  nsim<-nrow(B_BMSY)
  cex.leg<-0.9
  vlarg<-1e20
  
  FMSYr <- quantile(F_FMSY, c(0.001, 0.9), na.rm = T)
  BMSYr <- quantile(B_BMSY, c(0.001, 0.975), na.rm = T)
  
  colsse <- cols<-rainbow(nyears, start = 0.63, end = 0.95)[1:nyears]
  colsse <- makeTransparent(cols, 95)
  
  XLim <- c(0, 3)
  YLim <- c(0, 2.5)
  YLim[2]<-max(YLim[2],F_FMSY*1.1)
  
  plot(c(B_BMSY[1, 1], B_BMSY[1, 2]), c(F_FMSY[1,1], F_FMSY[1, 2]), xlim = XLim, ylim = YLim, col = colsse[1], type = "l", axes = FALSE,xlab="",ylab="")
  
  polygon(c(-100,100,100,-100),c(-100,-100,100,100),col='grey94',border="grey94")
  
  axis(1,c(-vlarg,vlarg),c(-vlarg,vlarg))
  axis(2,c(-vlarg,vlarg),c(-vlarg,vlarg))
  
  axis(side = 2, labels = TRUE, las = 1)
  axis(side = 1, labels = TRUE)
  
  OO <- round(sum(B_BMSY[,nyears] < 1 & F_FMSY[,nyears] > 1, na.rm = T)/nsim * 100, 1)
  OU <- round(sum(B_BMSY[,nyears] > 1 & F_FMSY[,nyears] > 1, na.rm = T)/nsim * 100, 1)
  UO <- round(sum(B_BMSY[,nyears] < 1 & F_FMSY[,nyears] < 1, na.rm = T)/nsim * 100, 1)
  UU <- round(sum(B_BMSY[,nyears] > 1 & F_FMSY[,nyears] < 1, na.rm = T)/nsim * 100, 1)
  
  abline(h = c(0,1), col = "white", lwd = 2)
  abline(v = c(0,1), col = "white", lwd = 2)
  
  y <- 1:(nyears - 1)
  y1 <- y + 1
  x0 <- as.vector(B_BMSY[, y])
  x1 <- as.vector(B_BMSY[, y1])
  y0 <- as.vector(F_FMSY[, y])
  y1 <- as.vector(F_FMSY[, y1])
  segments(x0, y0, x1, y1, col = rep(colsse,each=nsim))
  
  rng <- 1:min(maxsim, nsim)
  points(B_BMSY[rng, 1], F_FMSY[rng, 1], pch = 19, cex = 0.8, col = colsse[1])
  points(B_BMSY[rng, nyears],F_FMSY[rng,nyears], pch = 19, cex = 0.8, col = colsse[nyears])
  
  text(B_BMSY[1, ],F_FMSY[1,],yrlab,cex=0.7,font=2)  
  
  legend("right", legend=c(yrlab[1], yrlab[nyears]), bty = "n", text.col = c(cols[1],  cols[nyears]), pch = 19, col = c(cols[1], cols[nyears]))
  
  legend("topleft", paste(OO, "%", sep = ""), bty = "n", text.font = 2, cex = cex.leg,text.col='grey39')
  legend("topright", paste(OU, "%", sep = ""), bty = "n", text.font = 2,  cex = cex.leg,text.col='grey39')
  legend("bottomleft", paste(UO, "%", sep = ""), bty = "n", text.font = 2, cex = cex.leg,text.col='grey39')
  legend("bottomright", paste(UU, "%", sep = ""), bty = "n", text.font = 2, cex = cex.leg,text.col='grey39')
  
  mtext(expression(B/B[MSY]), 1, line = 2.5, cex = 1)
  mtext(expression(F/F[MSY]), 2, line = 2, cex = 1)
  
}
