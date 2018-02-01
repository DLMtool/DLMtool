
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
DFO_proj <- function(MSEobj,maxplot=3) {
  
  nsim<-MSEobj@nsim
  nMPs<-MSEobj@nMPs
  proyears<-MSEobj@proyears
  plotorg<-split(1:nMPs, ceiling(seq_along(1:nMPs)/maxplot))
  op<-par(mfrow=c(maxplot,2),mai=c(0.7,0.8,0.05,0.1),omi=c(0.01,0.01,0.5,0.01))
  
  for(i in 1:nMPs){
  
    DFO_Kobe(Br=MSEobj@B_BMSY[,i,proyears],Fr=MSEobj@F_FMSY[,i,proyears])
    legend('top',MSEobj@MPs[i],text.font=2,bty='n')
    DFO_Kobe_TS(Brel=MSEobj@B_BMSY[,i,],Frel=MSEobj@F_FMSY[,i,],labs=c("Current","Projection"))
    legend('top',MSEobj@MPs[i],text.font=2,bty='n')
    
    if(i==max(plotorg[[ceiling(i/maxplot)]]))  mtext(c("Status at end of projection (Projected)",
                                                       "Projected time series"),3,at=c(0.27,0.77),line=0.8,outer=T,font=2)
    
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
#' @author T. Carruthers
#' @export DFO_plot
DFO_plot<-function(MSEobj){
 
  par(mai=c(1,1,0.02,0.02))
  yend <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
  POF<-apply(MSEobj@F_FMSY[,,yend],2,mean,na.rm=T)
  
  POFed<-apply(MSEobj@B_BMSY[,,yend] ,2,mean, na.rm = T)
  
  col<-makeTransparent(c("red","dark green","blue","orange","black"),99)
  plot(POFed,POF,col="white",xlab="",ylab="",main="",axes=F)

  add_zones(textpos=quantile(POF,0.95))
  xs<-pretty(seq(min(POFed),max(POFed),length.out=8))
  ys<-pretty(seq(min(POF),max(POF),length.out=8))
  axis(1,xs,xs)
  axis(2,ys,ys)
  text(POFed,POF,MSEobj@MPs,col=col,font=2,cex=0.9)
  
  mtext("B/BMSY",1,line=2.5)
  mtext("F/FMSY",2,line=2.5)
  
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


DFO_Kobe<-function(Br,Fr){
  
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
  
  mtext("B/BMSY",1,line=2.5)
  mtext("F/FMSY",2,line=2.5)
  
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


add_zones<-function(textpos){

  cols<-c("grey86","grey96","white",
          "grey84","grey94","grey97")
  
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
