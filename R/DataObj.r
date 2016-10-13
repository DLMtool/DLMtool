
# A collection of functions which operate on DLM data object (class DLM_data) 
# which is used as input object for the management procedures. 

# These functions are used either internally in the MSE (runMSE.r) or used 
# manually to apply a MP (or MPs) to a particular data object.

# January 2016
# Tom Carruthers UBC (t.carruthers@fisheries.ubc.ca)
# Adrian Hordyk (a.hordyk@murdoch.edu.au)

DLMdiag <- function(DLM_data,command="available",reps=5,timelimit=1){
  funcs1<-c(avail("DLM_output"),avail("DLM_input"))
  good<-rep(TRUE,length(funcs1))
  report<-rep("Worked fine",length(funcs1))
  test<-new('list')
  timey<-new('list')
  options(show.error.messages = FALSE)
  for(y in 1:length(funcs1)){
    if(class(match.fun(funcs1[y]))=="DLM_output"){
      time1<-Sys.time()
      suppressWarnings({
        setTimeLimit(timelimit*1.5)
        test[[y]]<-try(do.call(funcs1[y],list(x=1,DLM_data=DLM_data,reps=5)),silent=T)
		if (class(test[[y]]) == "list") test[[y]] <-  test[[y]][[1]]
        setTimeLimit(Inf)
      })
    }else{
      time1<-Sys.time()
      suppressWarnings({
        setTimeLimit(timelimit*1.5)
        test[[y]]<-try(do.call(funcs1[y],list(x=1,DLM_data=DLM_data)),silent=T)
        setTimeLimit(Inf)
      })
    }
    time2<-Sys.time()
    timey[[y]]<-time2-time1
    if(class(test[[y]])=="try-error"){
      report[[y]]<-"Insufficient data"
      good[[y]]<-FALSE
    }else if(sum(is.na(test[[y]]))==length(test[[y]])){
      report[[y]]<-"Produced all NA scores"
      good[[y]]<-FALSE
    }
    if(timey[[y]]>timelimit){
      report[[y]]<-"Exceeded the user-specified time limit"
      good[[y]]<-FALSE
    }
  } # end of funcs
  options(show.error.messages = TRUE)
  if(command=="available")return(funcs1[good])
  if(command=="not available")return(cbind(funcs1[!good],report[!good]))
  if(command=="needed")return(needed(DLM_data,funcs=funcs1[!good]))
}

needed<-function(DLM_data,funcs=NA){
  if(is.na(funcs[1]))funcs<-avail("DLM_output")
  slots<-slotNames('DLM_data')
  rr <- try(slot(DLM_data, "Misc"), silent=TRUE)
  if (class(rr) == "try-error") DLM_data@Misc <- list()
  
  slotnams<-paste("DLM_data@",slotNames(DLM_data),sep="")
  repp<-rep("",length(funcs))
  DLM_data@Misc <- list()
  for(i in 1:length(funcs)){
    temp<-format(match.fun(funcs[i]))
    temp<-paste(temp[1:(length(temp))],collapse=" ")
    rec<-""
    for(j in 1:length(slotnams)){
      if(grepl(slotnams[j],temp)&NAor0(slot(DLM_data,slots[j])))rec<-c(rec,slots[j])
    }
    repp[i]<-paste(funcs[i],": ",paste(rec[2:length(rec)],collapse=", "),sep="")
  }
  repp
}

parallelMPs<-function(x,DLM_data,reps,MPs,ss) sapply(ss,MPs[x],DLM_data,reps=reps)

OneRep<-function(DLM_data){
  DLM_data@CV_Cat=DLM_data@CV_Dt=DLM_data@CV_AvC=DLM_data@CV_Ind=DLM_data@CV_Mort=DLM_data@CV_FMSY_M=DLM_data@CV_BMSY_B0=DLM_data@CV_Cref=DLM_data@CV_Bref=DLM_data@CV_Iref=DLM_data@CV_Rec=DLM_data@CV_Dep=DLM_data@CV_Abun=DLM_data@CV_L50=DLM_data@CV_vbK=DLM_data@CV_vbLinf=DLM_data@CV_vbt0=DLM_data@CV_LFC=DLM_data@CV_LFS=DLM_data@CV_wla=DLM_data@CV_wlb=DLM_data@CV_steep=DLM_data@sigmaL=tiny
  DLM_data
}

plotOFL<-function(DLM_data,xlims=NA,perc=0.5){

  cols<-rep(c('black','red','green','blue','orange','brown','purple','dark grey','violet','dark red','pink','dark blue','grey'),4)
  ltys<-rep(1:4,each=13)
  
  funcs<-DLM_data@MPs
  nMPs<-length(funcs)
   
  if(is.na(xlims[1])|length(xlims)!=2){
    xlims<-quantile(DLM_data@TAC,c(0.005,0.90),na.rm=T)
    if(xlims[1]<0)xlims[1]<-0
  }
  if(!NAor0(DLM_data@Ref)){
    if(xlims[1]>DLM_data@Ref)xlims[1]<-max(0,0.98*DLM_data@Ref)
    if(xlims[2]<DLM_data@Ref)xlims[2]<-1.02*DLM_data@Ref
    if(xlims[2]>DLM_data@Ref*2)xlims[2]<-2*DLM_data@Ref
  }
  ylims<-c(0,1)

  plot(NA,NA,xlim=xlims,ylim=ylims,main="",xlab="",ylab="",col="white",lwd=3,type="l")
  abline(h=0)
  if(!NAor0(DLM_data@Ref)){
    abline(v=DLM_data@Ref,col="light grey",lwd=3)
    if(!NAor0(DLM_data@Ref_type[1]))legend('bottomright',DLM_data@Ref_type,text.col="grey",bty='n')
  }
    
  for(m in 1:nMPs){
      
    if(sum(!is.na(DLM_data@TAC[m,,1]))>10){  # only plot if there are sufficient non-NA TAC samples
      x<-density(DLM_data@TAC[m,,1],from=0,na.rm=T)$x
      y<-density(DLM_data@TAC[m,,1],from=0,na.rm=T)$y
      y<-y/max(y)
      lines(x,y,col=cols[m])
    }else{
       print(paste("Method ",funcs[m]," produced too many NA TAC values for plotting densities",sep=""))
    }
    if(!is.na(perc[1]))abline(v=quantile(DLM_data@TAC[m,,1],p=perc,na.rm=T),col=cols[m],lty=2)
  }
    
  cind<-1:nMPs
  legend('topright',funcs,text.col=cols[cind],col=cols[cind],lty=1,bty='n',cex=0.75)
  
  mtext(paste("OFL (",DLM_data@Units,")",sep=""),1,outer=F,line=2.6)
  mtext(paste("Standardized relative frequency",sep=""),2,outer=F,line=2.6)
  #mtext(paste("OFL calculation for ",DLM_data@Name,sep=""),3,outer=F,line=1)
  
}

# Primary functions
Can <- function(DLM_data,timelimit=1)DLMdiag(DLM_data,"available",timelimit=timelimit)
Cant <- function(DLM_data,timelimit=1)DLMdiag(DLM_data,"not available",timelimit=timelimit)
Needed <- function(DLM_data,timelimit=1)DLMdiag(DLM_data,"needed",timelimit=timelimit)

# A function that determines the inputs for a given data-limited method of 
# class DLM_output and then analyses the sensitivity of TAC estimates to 
# marginal differences in each input. 
Sense <-function(DLM_data,MP,nsense=6,reps=100,perc=c(0.05,0.5,0.95),ploty=T){

  DLM_data2<-DLM_data
  nm <-deparse(substitute(DLM_data2))
  refTAC<-quantile(getTAC(DLM_data2,MP,reps)[[1]],perc,na.rm=T)
 
  DLM_data<-DLM_data2
  reqs<-Required(MP)#read.csv(paste(getwd(),"/Data/Data requirements.csv",sep=""),header=T)
  ind<-(1:nrow(reqs))[reqs[,match(MP,names(reqs))]=="Y"]
  #for(i in 1:length(reqs))
  
 
  slotsCV<-slotNames('DLM_data')[grep("CV_",slotNames('DLM_data'))]  
  slots<-rep("",length(slotsCV))
  for(i in 1:length(slotsCV))slots[i]<-substr(slotsCV[i],4,nchar(slotsCV[i]))
   
  ind<-slots%in%unlist(strsplit(reqs[2],", "))
  slots<-slots[ind]
  slotsCV<-slotsCV[ind]
  sname<-slots
  nslots<-length(slots)

  nrep<-nslots*nsense
  DLM_data<-replic8(DLM_data,nrep)
  pss<-seq(0,1,length.out=nsense+2)[2:(nsense+1)]
  vals<-array(NA,dim=c(nslots,nsense))

  for(i in 1:nslots){
    ind<-(((i-1)*nsense+1):(i*nsense))
    mn<-attr(DLM_data,slots[i])[1]
    cv<-attr(DLM_data,slotsCV[i])[1]*2 # twice the CV of the variable specified in the DLM object
    if(class(attr(DLM_data,slots[i]))=='numeric'){
      if(mn>0){
        attr(DLM_data,slots[i])[ind]<-qlnorm(pss,mconv(mn,cv*mn),sdconv(mn,cv*mn))
        vals[i,]<-qlnorm(pss,mconv(mn,cv*mn),sdconv(mn,cv*mn))
      }else{
        attr(DLM_data,slots[i])[ind]<--qlnorm(pss,mconv(-mn,cv*-mn),sdconv(-mn,cv*-mn))
        vals[i,]<--qlnorm(pss,mconv(-mn,cv*-mn),sdconv(-mn,cv*-mn))
      }
    }else{
      cv<-attr(DLM_data,slotsCV[i])[1]
      attr(DLM_data,slots[i])[ind,]<-attr(DLM_data,slots[i])[ind,]*qlnorm(pss,mconv(1,cv),sdconv(1,cv))
      vals[i,]<-qlnorm(pss,mconv(1,cv),sdconv(1,cv))
    }
  }

  TACa<-getTAC(DLM_data,MPs=MP,reps=reps)[[1]]
  TACa<-apply(TACa,3,quantile,p=perc,na.rm=T)
  LB<-((1:nslots)-1)*4+1
  UB<-(1:nslots)*4
  sense<-matrix(data=NA,nrow=4*nslots,ncol=nsense+1)
  
  for(i in 1:nslots){
    ind<-((i-1)*nsense+1):(i*nsense)
    dat<-TACa[,ind]

    sense[LB[i],2:(nsense+1)]<-vals[i,]
    sense[(LB[i]+1):UB[i],2:(nsense+1)]<-dat
    sense[LB[i],1]<-slots[i]
    sense[(LB[i]+1):UB[i],1]<-perc
  }
  
  DLM_data2@Sense<-sense

  if(ploty){
   ylimy<-range(TACa)
   #dev.new2(width=10,height=0.5+3*ceiling(nslots/2))
   par(mfrow=c(ceiling(nslots/2),2),mai=c(0.4,0.4,0.01,0.01),omi=c(0.4,0.4,0.4,0.01))
   for(i in 1:nslots){
     ind<-(((i-1)*nsense+1):(i*nsense))
     dat<-TACa[,ind]
     xlimy<-range(vals[i,])
     plot(xlimy,rep(refTAC[2],2),ylim=ylimy,xlim=xlimy,type='l',col="#99999960",main="",xlab="",ylab="")
     abline(h=refTAC[c(1,3)],col="#99999960",lty=2)
     abline(v=slot(DLM_data2,slots[i]),col="#99999960",lty=2)
     lines(vals[i,],dat[2,],col="red",lwd=1.5)
     lines(vals[i,],dat[1,],col="red",lty=2,lwd=1.5)
     lines(vals[i,],dat[3,],col="red",lty=2,lwd=1.5)
     legend('top',legend=sname[i],text.col='blue',bty='n')
   }

   mtext(paste("Output control (",DLM_data@Units,")",sep=""),2,outer=T,line=0.5)
   mtext("Parameter / variable input level",1,outer=T,line=0.5)
   mtext(paste("Sensitivity analysis for ",DLM_data@Name,": ",MP,sep=""),3,outer=T,line=0.5)
  }
  #assign(nm,DLM2,envir=.GlobalEnv)
  DLM_data2
}

# Replicates position 1 data to multiple positions for sensitivity testing etc
replic8<-function(DLM_data,nrep){

  slotnam<-slotNames(DLM_data)
  slotnam<-slotnam[slotnam!="Ref"&slotnam!="OM"&slotnam!="MaxAge"&slotnam!="CAL_bins"&slotnam!="Year"]
  
  for(sl in 1:length(slotnam)){
    slt<-attr(DLM_data,slotnam[sl])
    if(class(slt)=='matrix'){
      attr(DLM_data,slotnam[sl])<-matrix(rep(slt,each=nrep),nrow=nrep,ncol=ncol(slt))
    }else if(class(slt)=='numeric'){
      attr(DLM_data,slotnam[sl])<-rep(slt,nrep)
    }else if(class(slt)=='array'){
      attr(DLM_data,slotnam[sl])<-array(rep(slt,each=nrep),dim=c(nrep,dim(slt)[2:3]))
    }
  }
  DLM_data
}

# A function that returns the stochastic TAC recommendations from a vector of 
# data-limited MPs (DLM_output) given a data-limited data object DLM_data
TAC<-function(DLM_data,MPs=NA,reps=100,maxlines=6,perc=NA,xlims=NA,timelimit=1){

  nm <-deparse(substitute(DLM_data))
  PosMPs<-Can(DLM_data,timelimit=timelimit)
  PosMPs<-PosMPs[PosMPs%in%avail("DLM_output")]
  DLM_data@PosMPs<-PosMPs
  if(!is.na(MPs[1]))DLM_data@MPs<-MPs[MPs%in%PosMPs]
  if(is.na(MPs[1]))DLM_data@MPs<-PosMPs
  funcs<-DLM_data@MPs

  if(length(funcs)==0){
    stop("None of the methods 'MPs' are possible given the data available")
  }else{
    temp <- getTAC(DLM_data,MPs=funcs,reps) 
	TACa<- temp[[1]]
	DLM_data <- temp[[2]]
    DLM_data@TAC<-TACa
    return(DLM_data)
    #assign(nm,DLM,envir=.GlobalEnv)
  }

}

getTAC<-function(DLM_data,MPs=NA,reps=100){

  nsims<-length(DLM_data@Mort)
  nMPs<-length(MPs)
  TACa<-array(NA,dim=c(nMPs,reps,nsims))

  if(!sfIsRunning()|(nMPs<8&nsims<8)){
    for(ff in 1:nMPs){
      temp <- sapply(1:nsims,MPs[ff],DLM_data=DLM_data,reps=reps)
	  if (mode(temp) == "numeric") TACa[ff,,] <- temp
	  if (mode(temp) == "list") {
  	    TACa[ff,,] <- unlist(temp[1,])
	    for (x in 1:nsims) DLM_data@Misc[[x]] <- temp[2,x][[1]]
	  }	
    }
  }else{
    sfExport(list=c("DLM_data"))
    if(nsims<8){
      sfExport(list=c("MPs","reps"))
      for(ss in 1:nsims){
	    temp <<- t(sfSapply(1:length(MPs),parallelMPs,DLM_data=DLM_data,reps=reps,MPs=MPs,ss=ss))
		if (mode(temp) == "numeric") TACa[,,ss] <- temp
		if (mode(temp) == "list") {
		  Lens <- unlist(lapply(temp, length))
		  for (X in 1:length(Lens)) {
		    Classes <- unlist(lapply(temp[,X][[1]], class))
			if (length(unique(Classes)) == 1) { # no Misc object
			  TACa[X,,ss] <- unlist(temp[,X])
			} else { # a Misc object is include 
			  ind <- which(Classes == "list")
			  TACa[X,,ss] <- unlist(temp[,X][[1]][1:(ind-1),])
			  DLM_data@Misc[[ss]] <- temp[,X][[1]][ind,]
			}
		  }
	    } 
      }
    }else{
      for(ff in 1:nMPs){
	    temp <- sfSapply(1:nsims,MPs[ff],DLM_data=DLM_data,reps=reps)
		if (mode(temp) == "numeric") TACa[ff,,] <- temp
	    if (mode(temp) == "list") {
  	      TACa[ff,,] <- unlist(temp[1,])
	      for (x in 1:nsims) DLM_data@Misc[[x]] <- temp[2,x][[1]]
	    }
      }
    }
  }
  for(ff in 1:nMPs){
    if(sum(is.na(TACa[ff,,]))>sum(!is.na(TACa[ff,,]))){  # only plot if there are sufficient non-NA TAC samples
      print(paste("Method ",MPs[ff]," produced greater than 50% NA values",sep=""))
    }
  }
  out <- list(TACa, DLM_data)
  return(out)
}

Sam<-function(DLM_data,MPs=NA,reps=100,maxlines=10,perc=0.5){
  nm <-deparse(substitute(DLM))
  DLM_data@PosMPs<-MPs
  funcs<-DLM_data@PosMPs
  nMPs<-length(funcs)
  DLM_data@MPs<-funcs
  temp <- getTAC(DLM_data,MPs=funcs,reps)
  TACa<- temp[[1]]
  DLM_data <- temp[[2]]
  nsim<-length(DLM_data@Mort)
  ref<-array(rep(DLM_data@Ref,nMPs),c(nsim,nMPs))
  TACm<-apply(TACa,c(3,1),quantile,p=perc,na.rm=T)
  TACbias<-(TACm-ref)/ref *100
  POF<-round(apply(TACbias>0,2,sum)/length(DLM_data@Mort)*100,1)
  DLM_data@TAC<-TACa
  DLM_data@TACbias<-TACbias
  DLM_data
}


# Input Control Functions 
# Wrapper function for input control methods 
runInMP <- function(DLM_data,MPs=NA,reps=100) {

  nsims<-length(DLM_data@Mort)
  nMPs<-length(MPs)
  # len <- 4 + DLM_data@MaxAge
  len <- 8
  InC <- array(NA, dim=c(len, nsims, nMPs))
  if(!sfIsRunning()|(nMPs<8&nsims<8)){
    for(ff in 1:nMPs){
      temp <- sapply(1:nsims,MPs[ff],DLM_data=DLM_data,reps=reps)
	  if (mode(temp) == "numeric") {
	    Nrow <- nrow(temp)
		if (Nrow < len) {
		  dif <- ( len - Nrow) 
		  temp <- rbind(temp, matrix(NA, nrow=dif, ncol=ncol(temp)))
		}	
	    InC[,,ff] <- temp
	  }	
	  if (mode(temp) == "list") {
	    temp2 <- unlist(temp[1,])
		temp2 <- matrix(temp2, ncol=nsims)
		Nrow <- nrow(temp2)
		if (Nrow < len) {
		  dif <- ( len - Nrow) 
		  temp2 <- rbind(temp2, matrix(NA, nrow=dif, ncol=ncol(temp2)))
		}
  	    InC[,,ff] <- temp2
	    for (x in 1:nsims) DLM_data@Misc[[x]] <- temp[2,x][[1]]
	  }	
    }
  }else{
    sfExport(list=c("DLM_data"))
    if(nsims<8){
      sfExport(list=c("MPs","reps"))
      for(ss in 1:nsims){
	    temp <- t(sfSapply(1:length(MPs),parallelMPs,DLM_data=DLM_data,reps=reps,MPs=MPs,ss=ss))
		if (mode(temp) == "numeric") {
	      Nrow <- nrow(temp)
		  if (Nrow < len) {
		    dif <- ( len - Nrow) 
		    temp <- rbind(temp, matrix(NA, nrow=dif, ncol=ncol(temp)))
		  }	
	      InC[,,ff] <- temp
	    }	
		if (mode(temp) == "list") {
		  Lens <- unlist(lapply(temp, length))
		  ind <- which(Lens > 1) # these have Misc objects
		  for (X in 1:length(Lens)) {
		    if (any(X == ind)) {
			  temp2 <- unlist(temp[[1]][1,X])
			  temp2 <- matrix(temp2, ncol=nsims)
			  Nrow <- nrow(temp2)
		      if (Nrow < len) {
		       dif <- ( len - Nrow) 
		       temp2 <- rbind(temp2, matrix(NA, nrow=dif, ncol=ncol(temp2)))
		      }
			  InC[,ss,X] <- temp2
			  DLM_data@Misc[[ss]] <- temp[[1]][2,X][[1]]
			} else {
			  InC[,ss,X] <- unlist(temp[,X])
			}
		  }
	    } 
      }
    }else{
      for(ff in 1:nMPs){
	    temp <- sfSapply(1:nsims,MPs[ff],DLM_data=DLM_data,reps=reps)
		if (mode(temp) == "numeric") {
		  Nrow <- nrow(temp)
		  if (Nrow < len) {
		    dif <- (len - Nrow) 
		    temp <- rbind(temp, matrix(NA, nrow=dif, ncol=ncol(temp)))
		  }	
	      InC[,,ff] <- temp
		} 
	    if (mode(temp) == "list") {
		  temp2 <- unlist(temp[1,])
		  temp2 <- matrix(temp2, ncol=nsims)
		  Nrow <- nrow(temp2)
		  if (Nrow < len) {
		    dif <- ( len - Nrow) 
		    temp2 <- rbind(temp2, matrix(NA, nrow=dif, ncol=ncol(temp2)))
		  }
  	      InC[,,ff] <- temp2
	      for (x in 1:nsims) DLM_data@Misc[[x]] <- temp[2,x][[1]]
	    }
      }
    }
  }

  out <- list(InC, DLM_data)
  return(out)
}


# Plot TAC 
boxplot.DLM_data <- function(x, outline=FALSE, ...) {
  x <- DLM_data
  if (class(DLM_data) != "DLM_data") stop("Object must be of class 'DLM_data'")
  if (all(is.na(dim(DLM_data@TAC[,,1])))) stop("No TAC data found")
  tacs <- t(DLM_data@TAC[,,1])
  MPs <- DLM_data@MPs 
  ind <- grep("ref", MPs)
  if (length(ind) > 0) {
    tacs <- tacs[,-ind]
    MPs <- MPs[-ind]
  }
  ord <- order(apply(tacs, 2, median, na.rm=TRUE))
  MPs <- MPs[ord]
  tacs <- tacs[,ord]
  
  cols <- rainbow(30)
  ymax <- quantile(apply(tacs, 2, quantile, 0.99, na.rm=TRUE), 0.99)
  ymin <- quantile(apply(tacs, 2, quantile, 0.01, na.rm=TRUE), 0.01)
  ylim <- c(ymin, ymax)
  par(mfrow=c(1,1), oma=c(2,4,1,0), mar=c(3,3,0,0))
  boxplot(tacs, names=MPs, las=1, col=cols, outline=outline, frame=FALSE, ylim=ylim, 
          horizontal=TRUE, ...)
  
  mtext(paste("TAC (",DLM_data@Units,")",sep=""), side=1, outer=T,line=0.5, cex=1.25)
  mtext(side=2, "Management Procedures", outer=TRUE, line=3, cex=1.25)
  mtext(paste("TAC calculation for ",DLM_data@Name,sep=""),3,outer=T,line=-0.5, cex=1.25)
  
  data.frame(MP=MPs, Median=round(apply(tacs, 2, median, na.rm=TRUE),2),
             SD=round(apply(tacs, 2, sd, na.rm=TRUE),2), Units=DLM_data@Units)
  
}

# run input control method 
Input <- function(DLM_data, MPs=NA, reps=100, timelimit=10, CheckMPs=TRUE) {
  print("Checking which MPs can be run")
  flush.console()
  if(CheckMPs) PosMPs <- Can(DLM_data, timelimit=timelimit)
  if(!CheckMPs) PosMPs <- MPs
  PosMPs <- PosMPs[PosMPs%in%avail("DLM_input")]
  if(!is.na(MPs[1]))DLM_data@MPs<-MPs[MPs%in%PosMPs]
  if(is.na(MPs[1]))DLM_data@MPs<-PosMPs
  funcs<-DLM_data@MPs

  if(length(funcs)==0){
    stop("None of the methods 'MPs' are possible given the data available")
  }else{
    Out <- matrix(NA, nrow=length(funcs), ncol=6)
	colnames(Out) <- c("Effort", "Area 1",  "Area 2", "SL50", "SL95", "UpperLimit")
	rownames(Out) <- funcs
    for (mm in 1:length(funcs)) {
	  print(paste("Running", mm, "of", length(funcs), "-", funcs[mm]))
	  flush.console()
      runIn <- runInMP(DLM_data,MPs=funcs[mm], reps=reps)[[1]][,,1]
	  Out[mm,] <- runIn[2:7]
	  Out[,4:6] <- round(Out[,4:6], 2)
	}
  }
  Out 

}

