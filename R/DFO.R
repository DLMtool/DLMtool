
# ==========================================================================
# === Tools for COSEWIC designation ========================================
# ==========================================================================

# Generic tradeoffplot
tradeoffplot <- function(x, y, xlab, ylab, labs, cex, vl, hl) {
  adjj <- c(0.7, 1.3)
  XLim <- c(min(c(-10, min(x, na.rm = T) * adjj)), max(c(max(x, na.rm = T) *
                                                           adjj, 110)))
  YLim <- c(min(c(-10, min(y, na.rm = T) * adjj)), max(c(max(y, na.rm = T) *
                                                           adjj, 110)))
  coly <- rep(c("#0000ff95", "#ff000095", "#20ff1095"), 50)[1:length(labs)]
  #coly[labs %in% c("AvC", "curE", "FMSYref")] <- "black"
  col[grepl%in%c("AvC","curE","FMSYref")]<-black
  # plot(NA,xlim=range(x,na.rm=T)*adjj,ylim=range(y,na.rm=T)*adjj,xlab=xlab,ylab=ylab)
  plot(NA, xlim = XLim, ylim = YLim, xlab = xlab, ylab = ylab)
  abline(v = vl, col = "#99999940", lwd = 2)
  abline(h = hl, col = "#99999940", lwd = 2)
  text(x, y, labs, font = 2, col = coly, cex = 0.9)
}


#' Deparment of Fisheries and Oceans default plot 1
#'
#' A preliminary plot for returning trade-offs plots and performance table for
#' total yield, variability in yield, probability of overfishing and likelihood
#' of biomass dropping below 50 per cent BMSY
#'
#'
#' @usage DFO_plot(MSEobj,nam=NA,type=NA,panel=T)
#' @param MSEobj An object of class MSE
#' @param nam Title of plot
#' @param type Plots full range of data if NA. Plots a subset that meet thresholds if not NA.
#' @return A table of performance metrics.
#' @author T. Carruthers
#' @export DFO_plot2
DFO_plot2 <- function(MSEobj, nam = NA, type = NA, panel = T,Bcut=50, Ycut=50) {

  Yd <- rep(NA, MSEobj@nMPs)
  B50 <- rep(NA, MSEobj@nMPs)
  LTY <- rep(NA, MSEobj@nMPs)
  yend <- max(MSEobj@proyears - 4, 1):MSEobj@proyears

  refMP<-match("FMSYref",MSEobj@MPs)
  RefYd<-apply(MSEobj@C[, refMP, yend],1,mean)

 # RefYd <- MSEobj@OM$RefY

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

