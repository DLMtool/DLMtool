
# ==========================================================================
# === Tools for COSEWIC designation ========================================
# ==========================================================================


#' Deparment of Fisheries and Oceans historical plot
#'
#' A plot of current and historical stock status by simulation according to the
#' stock status zones and reference points of DFO. 
#' http://www.dfo-mpo.gc.ca/reports-rapports/regs/sff-cpd/precaution-eng.htm
#'
#' @param OM An operating model object of class OM
#' @param panel should the plots be seperate or in two panels?
#' @param nsim how many simulations should be plotted (over-ridden 
#' by OM@nsim where cpars is specified) 
#' @author T. Carruthers
#' @export DFO_hist
#' @importFrom MASS kde2d
DFO_hist <- function(OM, panel= T,nsim=48) {
  
  if(length(OM@cpars)>0){
    message("cpars slot is specified, OM@nsim used for historical simulations")
    nsim<-OM@nsim
  }else{
    OM@nsim=nsim
  }
  
  out<-runMSE(OM,Hist=T)
  Brel<-t(out$TSdata$SSB)/out$MSYs$SSBMSY
  Frel<-t(-log(1-out$TSdata$Catch/(out$TSdata$VB+out$TSdata$Catch)))/out$MSYs$FMSY
 
  if(panel)op<-par(mfrow=c(1,2),mai=c(0.7,0.8,0.5,0.1),omi=rep(0.01,4))
  if(!panel)opt<-par(mai=c(0.7,0.8,0.5,0.1),omi=rep(0.01,4))
  
  now<-ncol(Brel)
  DFO_Kobe(Br=Brel[,now],Fr=Frel[,now])
  mtext("Current status",3,line=0.7,font=2)
  DFO_Kobe_TS(Brel,Frel)
  mtext("Historical time series",3,line=0.7,font=2)
  
  on.exit(par(op))

}

#' Deparment of Fisheries and Oceans projection plot
#'
#' A projection plot of MP performance by simulation according to the
#' stock status zones and reference points of DFO. 
#' http://www.dfo-mpo.gc.ca/reports-rapports/regs/sff-cpd/precaution-eng.htm
#'
#' @param MSEobj An operating model object of class MSE
#' @param maxplot The maximum number of MPs to be plotted per figure 
#' @author T. Carruthers
#' @export DFO_proj
#' @importFrom MASS kde2d
#' @importFrom graphics arrows contour
#' @importFrom grDevices dev.off jpeg
#' @importFrom stats acf
DFO_proj <- function(MSEobj,maxplot=6) {
  
  maxplot<-min(maxplot,MSEobj@nMPs)
  nsim<-MSEobj@nsim
  nMPs<-MSEobj@nMPs
  proyears<-MSEobj@proyears
  plotorg<-split(1:nMPs, ceiling(seq_along(1:nMPs)/maxplot))
  ncol<-ceiling(maxplot*0.33)
  nrow<-ceiling(maxplot/ncol)
  
  for(j in 1:length(plotorg)){
    op<-par(mfrow=c(nrow,ncol),mai=c(0.3,0.3,0.05,0.1),omi=c(0.3,0.4,0.4,0.01))
    
    for(i in plotorg[[1]]){
    DFO_Kobe(Br=MSEobj@B_BMSY[,i,proyears],Fr=MSEobj@F_FMSY[,i,proyears],xlab="",ylab="")
    legend('top',MSEobj@MPs[i],text.font=2,bty='n')
    #DFO_Kobe_TS(Brel=MSEobj@B_BMSY[,i,],Frel=MSEobj@F_FMSY[,i,],labs=c("Current","Projection"))
    #legend('top',MSEobj@MPs[i],text.font=2,bty='n')
    
     
    
    }
    mtext("Conditions at the end of the MSE projection", 3,line=0.8,outer=T,font=2)
    mtext("B/BMSY",1,line=0.7,outer=T,font=2)
    mtext("F/FMSY",2,line=0.7,outer=T,font=2)
    
  }

  on.exit(par(op))
  
}

#' Deparment of Fisheries and Oceans trade-off plot
#'
#' A plot of mean biomass relative to BMSY and fishing mortality rate relative to FMSY 
#' over the final 5 years of the projection
#' http://www.dfo-mpo.gc.ca/reports-rapports/regs/sff-cpd/precaution-eng.htm
#'
#' @param MSEobj An MSE object of class MSE produced by DLMtool function runMSE
#' @param zero_origin Logical: should plots have a zero-zero origin?
#' @author T. Carruthers
#' @export DFO_plot
DFO_plot<-function(MSEobj,zero_origin=T){
 
  op<-par(mai=c(1,1,0.02,0.02))
  yend <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
  POF<-apply(MSEobj@F_FMSY[,,yend],2,mean,na.rm=T)
  
  POFed<-apply(MSEobj@B_BMSY[,,yend] ,2,mean, na.rm = T)
  
  ylim<-range(POF)
  xlim<-range(POFed)
  if(zero_origin)ylim[1]<-xlim[1]<-0
  col<-makeTransparent(c("red","dark green","blue","orange","black"),99)
  plot(POFed,POF,xlim=xlim,ylim=ylim,col="white",xlab="",ylab="",main="",axes=F)

  add_zones(textpos=quantile(POF,0.95))
  xs<-pretty(seq(xlim[1],xlim[2],length.out=8))
  ys<-pretty(seq(ylim[1],ylim[2],length.out=8))
  axis(1,xs,xs)
  axis(2,ys,ys)
  text(POFed,POF,MSEobj@MPs,col=col,font=2,cex=0.9)
  
  mtext("B/BMSY",1,line=2.5)
  mtext("F/FMSY",2,line=2.5)
  on.exit(par(op))
  
}

#' Deparment of Fisheries and Oceans stock status bar plot
#'
#' A plot of biomass relative to BMSY over projected years
#'
#' @param MSEobj An MSE object of class MSE produced by DLMtool function runMSE
#' @param yres Integer: the year interval over which to calculate B/BMSY in future years
#' @author T. Carruthers
#' @export DFO_bar
DFO_bar<-function(MSEobj,yres=10){
  
  sections<-(0:floor(MSEobj@proyears/yres))*yres
  nsec<-length(sections)-1
  op<-par(mfrow=c(nsec,1),mai=c(0.05,0.7,0.02,0.02),omi=c(0.6,0.35,0.01,0.01))
  nsim<-MSEobj@nsim
  nMPs<-MSEobj@nMPs
  POF<-array(NA,c(nMPs,5,nsec))
  
  for(i in 1:nsec){
    yind <- (sections[i]+1):sections[i+1]
    POF[,,i]<-t(apply(MSEobj@B_BMSY[,,yind] ,2,quantile,p=c(0.05,0.25,0.5,0.75,0.95), na.rm = T))
  }
 
  xlim<-c(0,min(2.5,max(POF)))
  xs<-pretty(seq(xlim[1],xlim[2],length.out=8))
  
  barpos<-seq(0,0.98,length.out=nMPs+2)[2:(nMPs+1)]
  
  for(i in 1:nsec){
    
    plot(xlim,c(0,1),col='white',axes=F,main="",xlab="",ylab="")
    add_zones_bar()
    legend('top',paste0("Projection years ",sections[i]+1,"-",sections[i+1]),bty='n')
    points(POF[,3,i],barpos,pch=19,cex=1.5)
    
    for(MP in 1:nMPs){
      lines(POF[MP,c(2,4),i],rep(barpos[MP],2),lwd=3)
      lines(POF[MP,c(1,5),i],rep(barpos[MP],2))
    }
    
    axis(2,MSEobj@MPs,at=barpos,las=2,tick=F)
    if(i==nsec)axis(1,xs,xs)
    if(i<nsec)axis(1,xs,rep("",length(xs)))
    
  }
  
  mtext("B/BMSY",1,line=2.5,outer=T)
  mtext("Management Procedure",2,line=0.2,outer=T)
  on.exit(par(op))
}

DFO_Kobe_TS<-function(Brel,Frel,labs=c("Unfished","Current")){
  
  nsim<-nrow(Brel)
  ny<-ncol(Brel)
  Brange<-c(0,quantile(Brel,0.99))
  Frange<-c(0,quantile(Frel,0.99))
  
  sampcols<-c("red","blue","green")
  plot(Brange,Frange,axes=F,col="white",xlab="",ylab="")
  
  textpos<-Frange[1]+0.85*(Frange[2]-Frange[1])
  add_zones(textpos)
  
  xp<-pretty(seq(0,max(Brange),length.out=10))
  yp<-pretty(seq(0,max(Frange),length.out=10))
  axis(1,xp,xp)
  axis(2,yp,yp)
  
  mtext("B/BMSY",1,line=2.5)
  mtext("F/FMSY",2,line=2.5)
  
  for(i in 1:2)lines(Brel[i,],Frel[i,],col=sampcols[i])  
  for(i in 1:2){
    points(Brel[i,1],Frel[i,1],pch=3,lwd=2,col=sampcols[i],cex=1.5)
    points(Brel[i,ny],Frel[i,ny],pch=19,lwd=2,col=sampcols[i],cex=1.5)
  }
  
  pointcol<-makeTransparent('black',80)
  linecol<-"black"#makeTransparent('black',95)
 
  encircle(Brel[,1],Frel[,1],col=linecol,xrange=Brange+c(0,0.1),yrange=Frange+c(0,0.1),perc=0.1,lwd=1.2,labels=labs[1],lty=2)
  
  encircle(Brel[,ny],Frel[,ny],col=linecol,xrange=Brange+c(0,0.1),yrange=Frange+c(0,0.1),perc=0.1,lwd=2,labels=labs[2],lty=2)
  
  Bmu<-apply(Brel,2,quantile,0.5)
  Fmu<-apply(Frel,2,quantile,0.5)
  lines(Bmu,Fmu,col="grey40",lwd=4)
  
  afun<-function(x1,x2)x1-(x1-x2)*c(0.25,0.75)
  arrows(afun(Bmu[1],Bmu[ny])[1],
         afun(Fmu[1],Fmu[ny])[1],
         afun(Bmu[1],Bmu[ny])[2],
         afun(Fmu[1],Fmu[ny])[2],col='black',lwd=2)
  
  points(Bmu[1],Fmu[1],pch=3,lwd=3,col="black",cex=2)
  points(Bmu[ny],Fmu[ny],pch=19,lwd=3,col="black",cex=2)
  
  legend("topright",legend=c("Mean trend","Sim 1","Sim 2"),bty='n',text.col=c("black","red","blue"),text.font=c(2,1,1),cex=0.9)
  legend("right",legend=labs,pch=c(3,19),cex=0.9)
  
  
}


DFO_Kobe<-function(Br,Fr,xlab=NA,ylab=NA){
  
  if(is.na(xlab))xlab="B/BMSY"
  if(is.na(ylab))ylab="F/FMSY"
  
  Brange<-c(0,quantile(Br,0.99))
  Frange<-c(0,quantile(Fr,0.99))
  
  nsim<-length(Br)
  plot(Brange,Frange,axes=F,col="white",xlab="",ylab="")
  textpos<-Frange[1]+0.85*(Frange[2]-Frange[1])
  
  add_zones(textpos)
  
  xp<-pretty(seq(0,max(Brange),length.out=10))
  yp<-pretty(seq(0,max(Frange),length.out=10))
  axis(1,xp,xp)
  axis(2,yp,yp)
  
  mtext(xlab,1,line=2.5)
  mtext(ylab,2,line=2.5)
  
  pointcol<-makeTransparent('black',80)
  linecol<-"black"#makeTransparent('black',95)
  
  points(Br,Fr,pch=19,cex=0.9,col=pointcol)
  
  encircle(Br,Fr,col=linecol,xrange=Brange+c(0,0.1),yrange=Frange+c(0,0.1),perc=0.1,lty=2,lwd=1.2)
  encircle(Br,Fr,col=linecol,xrange=Brange+c(0,0.1),yrange=Frange+c(0,0.1),perc=0.50,lwd=1.2)
  points(quantile(Br,0.5),quantile(Fr,0.5),col='black',pch=3,cex=1.2,lwd=2)
 
  fracs<-c(sum(Br<0.4&Fr>1),sum(0.4<Br&Br<0.8&Fr>1),sum(0.8<Br&Fr>1),
           sum(Br<0.4&Fr<1),sum(0.4<Br&Br<0.8&Fr<1),sum(0.8<Br&Fr<1))
  fracs<-round(fracs/nsim*100,1)
  
  fposH<-Frange[1]+0.99*(Frange[2]-Frange[1])
  fposL<-Frange[1]+0.02*(Frange[2]-Frange[1])
  
  text(0.2,fposH,paste(fracs[1],"%"),col="white",cex=0.9,font=2)
  text(0.6,fposH,paste(fracs[2],"%"),col="grey73",cex=0.9,font=2)
  text(1.1,fposH,paste(fracs[3],"%"),col="grey73",cex=0.9,font=2)
  text(0.2,fposL,paste(fracs[4],"%"),col="white",cex=0.9,font=2)
  text(0.6,fposL,paste(fracs[5],"%"),col="grey73",cex=0.9,font=2)
  text(1.1,fposL,paste(fracs[6],"%"),col="grey73",cex=0.9,font=2)
  legend('right',legend=c("A simulation","Median"),pch=c(19,3),col=pointcol,cex=0.9)
  legend('topright',legend=c("50%","90%"),lty=c(1,2),col=linecol,bty='n',cex=0.9)
  
}


add_zones_bar<-function(){
  
  cols<-c("grey86","grey94","white",
          "grey84","grey92","grey97")
  
  polygon(c(-0,0.4,0.4,0),c(0,0,1,1),col=cols[1],border=cols[1])
  polygon(c(0.4,0.8,0.8,0.4),c(0,0,1,1),col=cols[2],border=cols[2])
  polygon(c(0.8,1000,1000,0.8),c(0,0,1,1),col=cols[3],border=cols[3])
  
  text(0.2,0.5,"Critical",col="white",font=2,srt=270,cex=1.1)
  text(0.6,0.5,"Cautious",col="grey73",font=2,srt=270,cex=1.1)
  text(1.1,0.5,"Healthy",col="grey73",font=2,srt=270,cex=1.1)
  
}

add_zones_bar_horiz<-function(textpos=0.2){
  
  cols<-c("grey86","grey94","white",
          "grey84","grey92","grey97")
  
  polygon(c(0,0,10E10,10E10),c(0,0.4,0.4,0),,col=cols[1],border=cols[1])
  polygon(c(0,0,10E10,10E10),c(0.4,0.8,0.8,0.4),col=cols[2],border=cols[2])
  polygon(c(0,0,10E10,10E10),c(0.8,1000,1000,0.8),,col=cols[3],border=cols[3])
  
  text(textpos,0.1,"Critical",col="white",font=2,cex=1.1)
  text(textpos,0.6,"Cautious",col="grey73",font=2,cex=1.1)
  text(textpos,1.1,"Healthy",col="grey73",font=2,cex=1.1)
  
}

add_zones<-function(textpos){

  cols<-c("grey86","grey94","white",
          "grey84","grey92","grey97")
  
  polygon(c(-0,0.4,0.4,0),c(0,0,1,1),col=cols[1],border=cols[1])
  polygon(c(0.4,0.8,0.8,0.4),c(0,0,1,1),col=cols[2],border=cols[2])
  polygon(c(0.8,1000,1000,0.8),c(0,0,1,1),col=cols[3],border=cols[3])
  
  polygon(c(0,0.4,0.4,0),c(1,1,1000,1000),col=cols[4],border=cols[4])
  polygon(c(0.4,0.8,0.8,0.4),c(1,1,1000,1000),col=cols[5],border=cols[5])
  polygon(c(0.8,1000,1000,0.8),c(1,1,1000,1000),col=cols[6],border=cols[6])
  
  text(0.2,textpos,"Critical",col="white",font=2,srt=270,cex=1.1)
  text(0.6,textpos,"Cautious",col="grey73",font=2,srt=270,cex=1.1)
  text(1.1,textpos,"Healthy",col="grey73",font=2,srt=270,cex=1.1)
  
}

encircle<-function(x,y,col="red",perc=0.05,xrange=NA,yrange=NA,log=F,lty=1,lwd=1,labels=NA,drawlabels=F){
  nsim<-length(x)
  
  if(log){
    x<-log(x)
    y<-log(y)
  }
  
  if(is.na(xrange[1]))xrange<-range(x)
  if(is.na(yrange[1]))yrange<-range(y)
  
  if(log){
    xrange[xrange==0]<-0.01
    yrange[yrange==0]<-0.01
    xrange<-log(xrange)
    yrange<-log(yrange)
  }

  kerneld <- MASS::kde2d(x, y, n = 100, lims = c(xrange, yrange))

  pp <- array()
  for (i in 1:nsim){
    z.x <- max(which(kerneld$x < x[i]))
    z.y <- max(which(kerneld$y < y[i]))
    pp[i] <- kerneld$z[z.x, z.y]
  }
  
  confidencebound <- quantile(pp, perc, na.rm = TRUE)
  
  if(log){
    kerneld$x<-exp(kerneld$x)
    kerneld$y<-exp(kerneld$y)
  }
  
 if(is.na(labels)){
   
   contour(kerneld, levels = confidencebound, col = col, add = TRUE,drawlabels=drawlabels,lty=lty,lwd=lwd)
 
 }else{
   
   contour(kerneld, levels = confidencebound, col = col, add = TRUE,drawlabels=T,lty=lty,lwd=lwd,labels=labels)
 
 }

}

#' Deparment of Fisheries and Oceans default plot 2
#'
#' A preliminary plot for returning trade-offs plots and performance table for
#' probability of obtaining half reference (FMSY) yield and probability of biomass 
#' dropping below 50 per cent BMSY
#'
#' @param MSEobj An object of class MSE
#' @param nam Title of plot
#' @param panel Should the plots be organized in many panels in a single figure
#' @param Bcut The cutoff biomass for satisficing (relative to BMSY)
#' @param Ycut the cutoff yield for satisficing (relative to reference yield)
#' @return A table of performance metrics.
#' @author T. Carruthers
#' @export DFO_plot2
DFO_plot2 <- function(MSEobj, nam = NA,panel = T,Bcut=50, Ycut=50) {

  Yd <- rep(NA, MSEobj@nMPs)
  B50 <- rep(NA, MSEobj@nMPs)
  LTY <- rep(NA, MSEobj@nMPs)
  yend <- max(MSEobj@proyears - 4, 1):MSEobj@proyears

  #refMP<-match("FMSYref",MSEobj@MPs)
  #RefYd<-apply(MSEobj@C[, refMP, yend],1,mean)

  RefYd <- MSEobj@OM$RefY

  for (mm in 1:MSEobj@nMPs) {

    B50[mm] <- round(sum(MSEobj@B_BMSY[, mm, ] >= 0.5, na.rm = T)/prod(dim(MSEobj@B_BMSY[, mm, ])) * 100, 1)
    LTY[mm] <- round(sum(MSEobj@C[, mm, yend]/RefYd >= 0.5, na.rm = T)/(MSEobj@nsim * length(yend)), 3) * 100

  }

  cond<-B50>Bcut&LTY>Ycut

  op<-par(mfrow=c(1,2),mai = c(1, 0.95, 0.1, 0.1), omi = c(0.1, 0.1, 0.4, 0))

  xlab="B50. Prob. biomass above half BMSY (%)"
  ylab="LTY. Prob. yield greater than half FMSY (%)"
  labs=MSEobj@MPs[1:MSEobj@nMPs]
  vl = Bcut
  hl = Ycut
  x<-B50
  y<-LTY
  adjj <- c(0.9, 1.1)
  XLim <- c(min(c(-10, min(x, na.rm = T) * adjj)), max(c(max(x, na.rm = T) *
                                                           adjj, 110)))
  YLim <- c(min(c(-10, min(y, na.rm = T) * adjj)), max(c(max(y, na.rm = T) *
                                                           adjj, 110)))
  coly <- rep(c("#0000ff95", "#ff000095", "#20ff1095"), 50)[1:length(labs)]
  #coly[labs %in% c("AvC", "curE", "FMSYref")] <- "black"
  coly[grepl("AvC",labs)]<-"black"
  coly[grepl("curE",labs)]<-"black"
  coly[grepl("FMSYref",labs)]<-"black"
  coly[!cond]<-makeTransparent('dark grey',60)
  # plot(NA,xlim=range(x,na.rm=T)*adjj,ylim=range(y,na.rm=T)*adjj,xlab=xlab,ylab=ylab)
  plot(NA, xlim = XLim, ylim = YLim, xlab = xlab, ylab = ylab)
  abline(v = vl, col = "#99999940", lwd = 2,lty=2)
  abline(h = hl, col = "#99999940", lwd = 2,lty=2)
  text(x, y, labs, font = 2, col = coly, cex = 0.9)

  mtext("All MPs",3,line=0.6)


  x<-B50[cond]
  y<-LTY[cond]

  adjj <- c(0.9, 1.1)
  XLim <- c(min(min(x, na.rm = T) * adjj), max(c(max(x, na.rm = T) *
                                                           adjj, 105)))
  YLim <- c(min(min(y, na.rm = T) * adjj), max(c(max(y, na.rm = T) *
                                                           adjj)))
  plot(NA, xlim = XLim, ylim = YLim, xlab = xlab, ylab = ylab)
  abline(v = vl, col = "#99999940", lwd = 2,lty=2)
  abline(h = hl, col = "#99999940", lwd = 2,lty=2)
  text(x, y, labs[cond], font = 2, col = coly[cond], cex = 0.9)

  mtext("Satisficed",3,line=0.6)

  temp <- data.frame(B50, LTY,Satisfice=c("Fail","Pass")[as.numeric(cond)+1])
  row.names(temp) <- MSEobj@MPs[1:MSEobj@nMPs]
  ord=order(temp[,1]*temp[,2],decreasing=T)
  temp[ord,]
  on.exit(par(op))

}


#' Deparment of Fisheries and Oceans biomass quantile plot
#'
#' A plot of biomass relative to BMSY quantiles over projected years 
#'
#' @param MSEobj An MSE object of class MSE produced by DLMtool function runMSE
#' @param maxcol Integer how many columns for panel plots?
#' @param qcol A color, the quantile coloration
#' @param lcol A color, the mean B/BMSY line
#' @param curyr The current calendar year
#' @param quants A vector 2 long for the quantiles e.g. 0.1 and 0.9 for the 10th and 90th quantiles
#' @param addline Should two individual example simulations be added to the plot?
#' @author T. Carruthers
#' @export DFO_quant
DFO_quant<-function(MSEobj,maxcol=6,qcol=rgb(0.4,0.8,0.95), lcol= "dodgerblue4",curyr=2018,quants=c(0.1,0.9),addline=T){
  
  if(is.na(maxcol))maxcol=ceiling(length(MSEobj@MPs)/0.5) # defaults to portrait 1:2
  MPs<-MSEobj@MPs
  nMPs<-length(MPs)
  yrs<-curyr+(1:MSEobj@proyears)
  
  plots<-split(1:nMPs, ceiling(seq_along(1:nMPs)/maxcol))
  
  nr<-length(plots)*2
  nc<-maxcol
  
  mat<-array(0,c(nc,nr*1.5))
  ind<-floor(0.5+(1:nr)*1.5)
  mat[,ind]<-1:(nr*nc)
  mat<-t(mat)
  ht<-rep(0.2,nr*1.5)
  ht[ind]<-1
  layout(mat,heights=ht)
  op<-par(mai=c(0.3,0.25,0.01,0.01),omi=c(0.5,0.6,0.05,0.05))
  
  B_BMSY<-MSEobj@B_BMSY
  Yd<-MSEobj@C/ MSEobj@OM$RefY
  
  Blims <- c(0,quantile(B_BMSY,0.95))
  Ylims<- c(0,quantile(Yd,0.95))
  
  plotquant<-function(x,p=c(0.1,0.9),yrs,qcol,lcol,addline=T){
    ny<-length(yrs)
    qs<-apply(x,2,quantile,p=p)
    polygon(c(yrs,yrs[ny:1]),c(qs[1,],qs[2,ny:1]),border=NA,col=qcol)
    
    if(addline)for(i in 1:2)lines(yrs,x[i,],col=lcol,lty=i)
    lines(yrs,apply(x,2,quantile,p=0.5),lwd=2,col="white")
  }
  
  for(pp in 1:length(plots)){
    
    toplot<-unlist(plots[pp])
    nt<-length(toplot)
    
    for(i in toplot){
      plot(range(yrs),Blims,col="white")
      add_zones_bar_horiz(textpos=curyr+0.3*MSEobj@proyears)
      axis(2)
      plotquant(B_BMSY[,i,],p=quants,yrs,qcol,lcol)
      mtext(MSEobj@MPs[i],3,line=0.2,font=2)
      if(i==toplot[1])mtext("B/BMSY",2,line=2.5)
    }
    if(nt<maxcol)for(i in 1:(maxcol-nt))plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="y label", xlab="x lablel",axes=F)
    
    for(i in toplot){
      plot(range(yrs),Ylims,col="white")
      plotquant(Yd[,i,],p=quants,yrs,qcol,lcol)
      if(i==toplot[1])mtext("Rel. Yd.",2,line=2.3)
    }
    if(nt<maxcol)for(i in 1:(maxcol-nt))plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="y label", xlab="x lablel",axes=F)
    
  }
  
  mtext("Projection Year",1,line=0.7,outer=T)
  on.exit(par(op))
}



#' COSEWIC forward projection plot
#'
#' Projection of biomass under three scenarios: no catch, FMSY fishing and status quo fishing
#' This plot if for an MSE object created from runMSE with the argument MPs=c("NFref","FMSYref","curE")
#'
#' @param MSEobj An object of class MSE created from runMSE() with the argument MPs=c("NFref","FMSYref","curE")
#' @param syear Starting year of the projection for graphing purposes
#' @return A plot
#' @author T. Carruthers
#' @export COSEWIC_plot
COSEWIC_plot<-function(MSEobj,syear=2015){

  if(sum(is.na(MSEobj@OM$Blow))>0)stop("At least one Blow calculation failed, try runMSE again with argument CalcBlow=TRUE")
          # ZeroC  FMSYref   curE
  cols<-c("green","orange","blue")

  cols[3]<-makeTransparent(cols[3])

  maxage<-dim(MSEobj@SSB_hist)[2]
  nyears<-dim(MSEobj@SSB_hist)[3]
  proyears<-dim(MSEobj@B)[3]
  refyears<-1:proyears
  irefyears<-proyears:1
  yrs<-refyears+syear-1
  hyrs<-((syear-nyears+1):syear)-1
  nsim<-dim(MSEobj@SSB_hist)[1]
  nMPs<-length(MSEobj@MPs)

  op<- par(mfcol=c(2,2),mai=c(0.45,0.8,0.1,0.01),omi=c(0.25,0.02,0.6,0.01))

  labs<-paste0("(",letters[1:4],")")
  ylabline<-2.8
  labline<-0.2

  # temporary fix prior to adding SSB projections to the MSE object
  SSBrat<-mean(apply(MSEobj@SSB_hist[,,nyears,],1,sum) / MSEobj@B[,2,1])
  SSBproj<-MSEobj@B*SSBrat
  SSBproj_rel<-SSBproj/MSEobj@OM$Blow

  # --- Historical plots -----------------

  SSB_hist<-apply(MSEobj@SSB_hist,c(1,3),sum)
  SSBMSYref<- SSB_hist/(MSEobj@OM$BMSY_B0*SSB_hist[,1])

  Cq<-apply(SSBMSYref,2,quantile,p=c(0.025,0.5,0.975))
  ylim<-c(0,6)

  plot(range(hyrs),range(ylim),col='white',xlab="",ylab="")
  abline(h=0.5,lty=2,lwd=1.5)
  abline(h=1,lwd=1.5)

  polygon(c(hyrs,hyrs[nyears:1]),
          c(Cq[1,1:nyears],Cq[3,nyears:1]),
          col=cols[3],border=F)
  lines(hyrs,Cq[2,1:nyears])

  #mtext(MSEobjnam[i],3,line=0.7)

  mtext(labs[1],3,line=labline,adj=0.02,cex=0.8)
  mtext("Hist. SSB / SSBMSY",2,line=ylabline)

  # --- Worm plots -----------------------

  # Absolute

  Cq<-apply(SSBproj,2:3,quantile,p=c(0.025,0.5,0.975))
  ylims<-c(0,max(Cq))

  plot(range(yrs),range(ylims),col='white',xlab="",ylab="")
  abline(h=0,lty=2,lwd=1.5)

  polygon(c(yrs,yrs[proyears:1]),
          c(Cq[1,3,refyears],Cq[3,3,irefyears]),
          col=cols[3],border=FALSE)
  lines(yrs,Cq[2,3,])

  for(MP in 1:2){
    lines(yrs,c(Cq[1,MP,refyears]),col=cols[MP],lwd=1.5)
    lines(yrs,c(Cq[3,MP,refyears]),col=cols[MP],lwd=1.5)
  }

  mtext("Proj. SSB",2,line=ylabline)
  mtext(labs[2],3,line=labline,adj=0.02,cex=0.8)



  # As a fraction of Blow -------------------------

  Cq<-apply(log(SSBproj_rel),2:3,quantile,p=c(0.025,0.5,0.975))
  ylims<-c(-2,max(Cq))

  plot(range(yrs),range(ylims),col='white',xlab="",ylab="")
  abline(h=0,lty=2,lwd=1.5)

  polygon(c(yrs,yrs[proyears:1]),
          c(Cq[1,3,refyears],Cq[3,3,irefyears]),
          col=cols[3],border=FALSE)
  lines(yrs,Cq[2,3,])

  for(MP in 1:2){
    lines(yrs,c(Cq[1,MP,refyears]),col=cols[MP],lwd=1.5)
    lines(yrs,c(Cq[3,MP,refyears]),col=cols[MP],lwd=1.5)
  }

  mtext("Proj. SSB / SSBlow (log)",2,line=ylabline)
  mtext(labs[3],3,line=labline,adj=0.02,cex=0.8)


  # --- probability plot -----------------

  Plim<-apply(SSBproj<MSEobj@OM$Blow,2:3,sum)/nsim*100
  matplot(x=yrs,t(Plim),ylim=c(0,10),type='l',lty=1,col=cols,xlab="",ylab="")

  mtext("Prob SSB < Blow",2,line=ylabline)
  mtext(labs[4],3,line=labline,adj=0.02,cex=0.8)

  legend('topleft',legend=c("No Catch","FMSY","Status quo"),bty='n',text.col=cols)

  mtext("Year",1,line=0.15,outer=T)
  
  mtext("Biomass projection plot relative to Blow for "
             , 3, line=2, outer = T, font = 2, cex = 0.8)
  
  mtext(MSEobj@Name,line=0.6, 3, outer = T, font = 2, cex = 0.8)
  on.exit(par(op))

}

#' Subset an OM cpars slot 
#'
#' Subset the custom parameters of an operating model 
#'
#' @param OM An object of class OM
#' @param sims A logical vector OM@nsim long of simulations to either retain (TRUE) or remove (FALSE)
#' @return An object of class OM
#' @author T. Carruthers
#' @export SubCpars
SubCpars<-function(OM,sims){
  
  for(i in 1:length(OM@cpars)){
    
    OM@nsim<-sum(sims)
    
    if(class(OM@cpars[[i]])=="matrix"){
      OM@cpars[[i]]<-OM@cpars[[i]][sims,]
    }else if(class(OM@cpars[[i]])=="array"){
      OM@cpars[[i]]<-OM@cpars[[i]][sims,,]
    }else{    
      OM@cpars[[i]]<-OM@cpars[[i]][sims]
    }
    
  }
  
  OM
  
}


# =================== Reports =========================================================

#' Create a standard DFO MSE report 
#'
#' Provides performacne plots typical in the assessment of Canadian fish stocks.  
#'
#' @param MSEobj An object of class MSE
#' @param output_file The directory and filename you wish to use for the report e.g. "C:/temp/myMSEreport.html"
#' @param author The person who made this report
#' @param title The title of the report
#' @author T. Carruthers
#' @export DFO_report
DFO_report<-function(MSEobj,output_file=NA,author="Author not specified",title=NA){
  
  if(is.na(output_file))output_file=paste0(getwd(),"/DFO MSE report.html")
  if(is.na(title))title=paste0("DFO Generic MSE Report for",MSEobj@Name)
  params<-new('list')
  params$title<-title
  params$subtitle<-"A prototype MSE performance evaluation"
  params$author<-author
  params$MSEobj<-MSEobj
  rmarkdown::render(input=system.file("DFO_generic.Rmd", package="DLMtool"), output_file=output_file,params=params)
  
}


#' Create a standard DFO performance table 
#'
#' P_Cr_S is the probability of being in the critical zone in the first 10 projected years
#' P_Ct_S is the probability of being in the cautious zone in the first 10 projected years
#' P_H_S is the probability of being in the healthy zone in the first 10 projected years
#' POF_S is the probability of overfishing in the first 10 projected years
#' STY is the mean yield relative to FMSY management over the first 10 projected years
#' P_Cr_L is the probability of being in the critical zone in the last 10 projected years
#' P_Ct_L is the probability of being in the cautious zone in the last 10 projected years
#' P_H_L is the probability of being in the healthy zone in the last 10 projected years
#' POF_L is the probability of overfishing in the last 10 projected years
#' LTY is the mean yield relative to FMSY management over the last 10 projected years
#' AAVY is the average annual variability in yield over the whole projection phrased as a CV percentage
#' P_Reb is the probability the stock has rebuilt to over BMSY in 2 mean generation times
#'
#' @param MSEobj An object of class MSE
#' @param rnd The number of significant figures for rounding. 
#' @author T. Carruthers
#' @export DFO_tab
DFO_tab<-function(MSEobj,rnd=0){
  
  shortterm<-1:min(10,MSEobj@proyears)
  P_Cr_S<-round(apply(MSEobj@B_BMSY[,,shortterm]<0.4,2,mean)*100,rnd)
  P_Ct_S<-round(apply(MSEobj@B_BMSY[,,shortterm]>0.4 & MSEobj@B_BMSY[,,shortterm]<0.8,2,mean)*100,rnd)
  P_H_S<-round(apply(MSEobj@B_BMSY[,,shortterm]>0.8,2,mean)*100,rnd)
  
  longterm<-max(1,MSEobj@proyears-10):MSEobj@proyears
  P_Cr_L<-round(apply(MSEobj@B_BMSY[,,longterm]<0.4,2,mean)*100,rnd)
  P_Ct_L<-round(apply(MSEobj@B_BMSY[,,longterm]>0.4 & MSEobj@B_BMSY[,,longterm]<0.8 ,2,mean)*100,rnd)
  P_H_L<-round(apply(MSEobj@B_BMSY[,,longterm]>0.8,2,mean)*100,rnd)
  
  refY<-MSEobj@OM$RefY
  STY<-round(apply(MSEobj@C[,,shortterm]/refY,2,mean)*100,rnd)
  LTY<-round(apply(MSEobj@C[,,longterm]/refY,2,mean)*100,rnd)
  
  POF_S<-round(apply(MSEobj@F_FMSY[,,shortterm]>1,2,mean)*100,rnd)
  POF_L<-round(apply(MSEobj@F_FMSY[,,longterm]>1,2,mean)*100,rnd)
  
  MGT2<-ceiling(MSEobj@OM$MGT*2)
  MGT2[MGT2<3]=3
  Bind<-cbind(as.matrix(expand.grid(1:MSEobj@nsim,1:MSEobj@nMPs)),rep(MGT2,MSEobj@nMPs))
  Bmat<-array(MSEobj@B_BMSY[Bind],c(MSEobj@nsim,MSEobj@nMPs))
  P_Reb<-round(apply(Bmat>1,2,mean)*100,rnd)
  
  y1 <- 1:(MSEobj@proyears - 1)
  y2 <- 2:MSEobj@proyears
  AAVY <- round(apply((((MSEobj@C[, , y1] - MSEobj@C[, , y2])/MSEobj@C[, , y2])^2)^0.5, 2, mean, na.rm = T)*100,rnd)
  
  MP<-MSEobj@MPs
  
  tab<-data.frame(MP,P_Cr_S, P_Ct_S, P_H_S, POF_S, STY, P_Cr_L, P_Ct_L, P_H_L, POF_L, LTY, AAVY, P_Reb)
  
  tab<-tab[order(tab$LTY,decreasing=T),]
  tab
  
}

#' A formatted version of the standard DFO performance plot, color coded by thresholds
#'
#' Crit_S is the probability of being in the critical zone in the first 10 projected years
#' Caut_S is the probability of being in the cautious zone in the first 10 projected years
#' Health_S is the probability of being in the healthy zone in the first 10 projected years
#' OvFish_S is the probability of overfishing in the first 10 projected years
#' Yield_S is the mean yield relative to FMSY management over the first 10 projected years
#' Crit is the probability of being in the critical zone in the last 10 projected years
#' Caut is the probability of being in the cautious zone in the last 10 projected years
#' Health is the probability of being in the healthy zone in the last 10 projected years
#' OvFish is the probability of overfishing in the last 10 projected years
#' Yield is the mean yield relative to FMSY management over the last 10 projected years
#' AAVY is the average annual variability in yield over the whole projection phrased as a CV percentage
#' Reb is the probability the stock has rebuilt to over BMSY in 2 mean generation times
#'
#' @param Ptab1 A DFO performance table made by DFO_tab()
#' @param thresh A vector of thresholds for each column Health, Yield and Reb are 'greater than threshold' conditions 
#' @author T. Carruthers
#' @importFrom dplyr %>% mutate
#' @importFrom kableExtra cell_spec kable_styling column_spec add_header_above
#' 
#' @export 
DFO_tab_formatted<-function(Ptab1,thresh=c(10,     40,     50,    50,    50,  10,     40,     50,    50,    50,  20,   50)){
  #                                       P_Cr_S, P_Ct_S, P_H_S, POF_S, STY, P_Cr_L, P_Ct_L, P_H_L, POF_L, LTY, AAVY, P_Reb
  
  # save(Ptab1,file="Ptab1")
  MPs<-as.character(Ptab1$MP)
  Data <- DLMtool::SimulatedData
  runMPs <- applyMP(Data, MPs, reps = 2, nsims=1, silent=TRUE)
  recs <- runMPs[[1]]
  type <- matrix(0, nrow=length(MPs),ncol=4) # TAC TAE SL MPA
  for (mm in seq_along(recs)) {
    type[mm,1] <- as.integer(length(recs[[mm]]$TAC) > 0)
    type[mm,2] <- as.integer(length(recs[[mm]]$Effort)>0)
    type[mm,3] <- as.integer(length(recs[[mm]]$LR5)>0)
    type[mm,4] <- as.integer(!is.na(recs[[mm]]$Spatial[1,1]))
  }
  totneeded<-apply(type,1,sum)
  
  MP_Type<-rep("TAC",length(MPs))
  MP_Type[type[,2]==1]<-"TAE"
  MP_Type[type[,3]==1]<-"SzLim"
  MP_Type[type[,4]==1]<-"MPA"
  MP_Type[totneeded>1]<-"Mixed"
  
  Ptab2<-Ptab1 #[,1:ncol(Ptab1)]
  Ptab2<-cbind(Ptab2[,1],MP_Type,Ptab2[,2:ncol(Ptab2)])
  MP <- Crit_S <- Caut_S <- Health_S <- OvFish_S <- Yield_S <- Crit <- Caut <-  Health <- OvFish <- Reb <- NULL # hack for CRAN checks
  names(Ptab2)<-c("MP","MP_Type","Crit_S","Caut_S","Health_S","OvFish_S","Yield_S","Crit","Caut","Health","OvFish","Yield","AAVY","Reb")
  
  #       P_Cr_S,               P_Ct_S,                 P_H_S,                 POF_S,                  STY,  AAVY, P_Reb
  PIsmet<-Ptab2[,3]<thresh[1] & Ptab2[,4] < thresh[2] & Ptab2[,5] >thresh[3] & Ptab2[,6] < thresh[4] & Ptab2[,7] > thresh[5] &
          Ptab2[,8]<thresh[6] & Ptab2[,9] < thresh[7] & Ptab2[,10]>thresh[8] & Ptab2[,11] < thresh[9]& Ptab2[,12]> thresh[10] &
          Ptab2[,13]< thresh[11] & Ptab2[,14] > thresh[12]
  MPcols<-rep('green',length(MPs))
  MPcols[!PIsmet]<-'red'
 
  # Rankings
  ord<-order(Ptab2$Yield,decreasing = T)
  Ptab2<-Ptab2[ord,]
  MPcols<-MPcols[ord]
  
  Ptab2 %>%
    mutate(
      #MP = row.names(.),
      MP =  cell_spec(MP, "html", color = MPcols, bold = T),
      MP_Type =  cell_spec(MP_Type, "html", bold = T),
      Crit_S = ifelse(Crit_S < thresh[1],
                      cell_spec(Crit_S, "html", color = "green",italic = T),
                      cell_spec(Crit_S, "html", color = "red", italic = T)),
      Caut_S = ifelse(Caut_S < thresh[2],
                      cell_spec(Caut_S, "html", color = "green", italic = T),
                      cell_spec(Caut_S, "html", color = "red", italic = T)),
      Health_S = ifelse(Health_S >= thresh[3],
                     cell_spec(Health_S, "html", color = "green", italic = T),
                     cell_spec(Health_S, "html", color = "red", italic = T)),
      OvFish_S = ifelse(OvFish_S < thresh[4],
                       cell_spec(OvFish_S, "html", color = "green", italic = T),
                       cell_spec(OvFish_S, "html", color = "red", italic = T)),
      Yield_S = ifelse(Yield_S >= thresh[5],
                      cell_spec(Yield_S, "html", color = "green", italic = T),
                      cell_spec(Yield_S, "html", color = "red", italic = T)),
      Crit = ifelse(Crit < thresh[6],
                      cell_spec(Crit, "html", color = "green", italic = T),
                      cell_spec(Crit, "html", color = "red", italic = T)),
      Caut = ifelse(Caut < thresh[7],
                      cell_spec(Caut, "html", color = "green", italic = T),
                      cell_spec(Caut, "html", color = "red", italic = T)),
      Health = ifelse(Health > thresh[8],
                   cell_spec(Health, "html", color = "green", italic = T),
                   cell_spec(Health, "html", color = "red", italic = T)),
      OvFish = ifelse(OvFish < thresh[9],
                       cell_spec(OvFish, "html", color = "green", italic = T),
                       cell_spec(OvFish, "html", color = "red", italic = T)),
      Yield = ifelse(Yield >= thresh[10],
                     cell_spec(Yield, "html", color = "green", italic = T),
                     cell_spec(Yield, "html", color = "red", italic = T)),
      AAVY = ifelse(AAVY < thresh[11],
                       cell_spec(AAVY, "html", color = "green", italic = T),
                       cell_spec(AAVY, "html", color = "red", italic = T)),
      Reb = ifelse(Reb >= thresh[12],
                     cell_spec(Reb, "html", color = "green", italic = T),
                     cell_spec(Reb, "html", color = "red", italic = T))
    )%>%
    #select(everything())%>%
    knitr::kable("html", escape = F,align = "c") %>%
    kable_styling("striped", full_width = F)%>%
    column_spec(5, width = "3cm")  %>%
    add_header_above(c(" ", " ","First 10 years of MSE projection" = 5, "Last 10 years of MSE projection" = 5,
                       "",""))
  
}

