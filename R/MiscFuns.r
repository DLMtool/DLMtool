# DLMtool Supporting Functions 
# January 2016
# Tom Carruthers UBC (t.carruthers@fisheries.ubc.ca)
# Adrian Hordyk (a.hordyk@murdoch.edu.au)

# Various short functions that are used internally in the DLMtool.
# Functions are only used internally, and not directly accessible or available 
# to users of the package and therefore do not appear in the help manual.

# Generic class finder
avail<-function(classy){
  chk <- "package:DLMtooldev" %in% search() 
  if (chk) { # development version
    return(unique(c(ls('package:DLMtooldev')[unlist(lapply(ls('package:DLMtooldev'),
      getclass,classy=classy))], ls(envir=.GlobalEnv)[unlist(
	  lapply(ls(envir=.GlobalEnv),getclass,classy=classy))]))) 
  } else {	  
    return(unique(c(ls('package:DLMtool')[unlist(lapply(ls('package:DLMtool'),
    getclass,classy=classy))], ls(envir=.GlobalEnv)[unlist(
	lapply(ls(envir=.GlobalEnv),getclass,classy=classy))])))
  }	
}


# A function that finds all methods in the environment and searches the function
# text for slots in the DLM data object
Required <- function(funcs=NA){
  if(is.na(funcs[1]))funcs<-c(avail("DLM_output"),avail("DLM_input"))
  slots<-slotNames('DLM_data')
  slotnams<-paste("DLM_data@",slotNames('DLM_data'),sep="")
  repp<-rep("",length(funcs))

  for(i in 1:length(funcs)){
    temp<-format(match.fun(funcs[i]))
    temp<-paste(temp[1:(length(temp))],collapse=" ")
    rec<-""
    for(j in 1:length(slotnams))if(grepl(slotnams[j],temp))rec<-c(rec,slots[j])
    if(length(rec)>1)repp[i]<-paste(rec[2:length(rec)],collapse=", ")
  }
  cbind(funcs,repp,deparse.level=0)
}

# A way of locating where the package was installed so you can find example 
# data files and code etc.
DLMDataDir<-function(stock=NA){
  chk <- "package:DLMtooldev" %in% search() 
  if (chk) { # dev version 
    if(is.na(stock)){
      return(paste(searchpaths()[match("package:DLMtooldev",search())],"/",sep=""))
    }else{
      return(paste(searchpaths()[match("package:DLMtooldev",search())],"/",stock,".csv",sep=""))
    }
  } else {
    if(is.na(stock)){
      return(paste(searchpaths()[match("package:DLMtool",search())],"/",sep=""))
    }else{
      return(paste(searchpaths()[match("package:DLMtool",search())],"/",stock,".csv",sep=""))
    }
  }  
}

# What MPs may be run (best case scenario) for various data-availability 
#  scenarios?
Fease<-function(feaseobj,outy="table"){
  
  if(class(feaseobj)!="DLM_fease")stop("Incorrect format: you need an object of class DLM_fease")
  
  sloty<-c("Cat","Ind","AvC","Dt","Rec","CAA","CAL","Mort","L50","L95","vbK",
           "vbLinf","vbt0","wla","wlb","steep","LFC","LFS","Cref","Bref","Iref","Dep","Abun")
  
  type<-c("Catch","Index","Catch","Index","Recruitment_index","Catch_at_age","Catch_at_length",
          "Natural_mortality_rate","Maturity_at_length","Maturity_at_length","Growth","Growth","Growth",
          "Length_weight_conversion","Length_weight_conversion","Stock_recruitment_relationship",
          "Fleet_selectivity","Fleet_selectivity","Target_catch","Target_biomass","Target_index",
          "Index","Abundance") 
  
  ncases<-length(feaseobj@Case)
  slots<-slotNames(feaseobj)
  ns<-length(slots)
  ftab<-array(TRUE,c(ns-2,ncases))
  for(j in 3:ns)ftab[j-2,]<-as.logical(as.numeric(slot(feaseobj,slots[j])))
  
  req<-Required()
  nMPs<-nrow(req)
  gridy<-array("",c(nMPs,ncases))
  for(i in 1:ncases){
    types<-slotNames(feaseobj)[3:17][ftab[,i]]
    slots<-sloty[type%in%types]
    for(m in 1:nMPs){
      brec<-unlist(strsplit(req[m,2],", "))
      brec<-brec[grep("CV_",brec,invert=T)] #remove CV dependencies (we think we can guess these...)
      brec<-brec[brec!="Year"&brec!="MaxAge"&brec!="FMSY_M"&brec!="BMSY_B0"&brec!="t"&brec!="OM"&brec!="MPrec"&brec!="CAL_bins"]
      nr<-length(brec) 
      if(nr==0){
        gridy[m,i]<-"Yes"
      }else{ 
        cc<-0
        for(r in 1:nr){ #loop over requirements
          if(brec[r]%in%slots)cc<-cc+1
        }
        if(cc==nr)gridy[m,i]<-"Yes"
      }
    }
  }
  gridy<-as.data.frame(gridy)
  row.names(gridy)=req[,1]
  names(gridy)=feaseobj@Case
  if(outy=="table")return(gridy)
  if(outy!="table"&class(outy)!="numeric")return(req[,1][gridy[,1]=="Yes"])
  if(class(outy)=="numeric"){
    if(outy<(ncases+1)){
      return(req[,1][gridy[,as.integer(outy)]=="Yes"])
    }else{
      return(req[,1][gridy[,1]=="Yes"])
    }  
  }
  
}


# Internal functions used in runMSE 
getmov<-function(x,Prob_staying,Frac_area_1){
  test<-optim(par=c(0,0,0),movfit,method="L-BFGS-B",lower=rep(-6,3),upper=rep(6,3),prb=Prob_staying[x],frac=Frac_area_1[x])
  mov<-array(c(test$par[1],test$par[2],0,test$par[3]),dim=c(2,2))
  mov<-exp(mov)
  mov/array(apply(mov,1,sum),dim=c(2,2))
}

movfit<-function(par,prb,frac){
  mov<-array(c(par[1],par[2],0,par[3]),dim=c(2,2))
  mov<-exp(mov)
  mov<-mov/array(apply(mov,1,sum),dim=c(2,2))
  dis<-c(frac,1-frac)
  for(i in 1:100)dis<-apply(array(dis,c(2,2))*mov,2,sum)
  (log(mov[1,1])-log(prb))^2+(log(frac)-log(dis[1]))^2
}

getq<-function(x,dep,Find,Perr,Marray,hs,Mat_age,Wt_age,R0,V,nyears,maxage,mov,Spat_targ,SRrel,aR,bR){
  opt<-optimize(qopt,log(c(0.0075,15)),depc=dep[x],Fc=Find[x,],Perrc=Perr[x,],
                     Mc=Marray[x,],hc=hs[x],Mac=Mat_age[x,],Wac=Wt_age[x,,],
                     R0c=R0,Vc=V[x,,],nyears=nyears,maxage=maxage,movc=mov[x,,],
                     Spat_targc=Spat_targ[x],SRrelc=SRrel[x],aRc=aR[x,],bRc=bR[x,])
  return(exp(opt$minimum))
}

qopt<-function(lnq,depc,Fc,Perrc,Mc,hc,Mac,Wac,R0c,Vc,nyears,maxage,movc,Spat_targc,SRrelc,aRc,bRc,opt=T){
  qc<-exp(lnq)
  nareas<-nrow(movc)
  #areasize<-c(asizec,1-asizec)
  idist<-rep(1/nareas,nareas)
  for(i in 1:300)idist<-apply(array(idist,c(2,2))*movc,2,sum)

  N<-array(exp(-Mc[1]*((1:maxage)-1))*R0c,dim=c(maxage,nareas))*array(rep(idist,each=maxage),dim=c(maxage,nareas))
  SSN<-Mac*N   # Calculate initial spawning stock numbers
  Biomass<-N*Wac[,1]
  SSB<-SSN*Wac[,1]                               # Calculate spawning stock biomass

  B0<-sum(Biomass)
  R0a<-idist*R0c
  SSB0<-apply(SSB,2,sum)
  SSBpR<-SSB0/R0a                              # Calculate spawning stock biomass per recruit

  for(y in 1:nyears){
    # set up some indices for indexed calculation
    targ<-(apply(Vc[,y]*Biomass,2,sum)^Spat_targc)/mean(apply(Vc[,y]*Biomass,2,sum)^Spat_targc)
    FMc<-array(qc*Fc[y]*Vc[,y],dim=c(maxage,nareas))*array(rep(targ,each=maxage),dim=c(maxage,nareas))                                           # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
    Zc<-FMc+Mc[y]
    N[2:maxage,]<-N[1:(maxage-1),]*exp(-Zc[1:(maxage-1),])         # Total mortality
    if(SRrelc==1){
      N[1,]<-Perrc[y]*(0.8*R0a*hc*apply(SSB,2,sum))/(0.2*SSBpR*R0a*(1-hc)+(hc-0.2)*apply(SSB,2,sum))  # Recruitment assuming regional R0 and stock wide steepness
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
    SSB<-SSN*Wac[,y]
    Biomass<-N*Wac[,y]
    SBiomass<-SSN*Wac[,y]
    #print(sum(Biomass))
  } # end of year
  return((log(depc)-log(sum(SBiomass)/sum(SSB0)))^2)
}


## Operating Model Functions ---------------------------------------------------
# These functions are used to manually specify, choose, or estimate various 
# parameters of the Operating Model.
# The functions typically take OM object (or Stock or Fleet) and return the same
# object with the relevant parameters populated.

# A highly dubious means of getting very uncertain estimates of current stock 
# biomass and (equilibrium) fishing mortality rate from growth, natural 
# mortality rate, recruitment and fishing selectivity.
ML2D<-function(OM,ML,nsim=100,ploty=T,Dlim=c(0.05,0.6)){
  
  maxage<-OM@maxage
  M<-runif(nsim,OM@M[1],OM@M[2]) # Natural mortality rate
  h<-runif(nsim,OM@h[1],OM@h[2]) # Steepness
  Linf<-runif(nsim,OM@Linf[1],OM@Linf[2]) # Maximum length
  K<-runif(nsim,OM@K[1],OM@K[2]) # Maximum growth rate
  t0<-runif(nsim,OM@t0[1],OM@t0[2]) # Theorectical length at age zero

  LFS<-runif(nsim,OM@LFS[1],OM@LFS[2])*mean(OM@L50)
  AFS<-L2A(t0,Linf,K,LFS,maxage)

  L5<-runif(nsim,OM@L5[1],OM@L5[2])*mean(OM@L50)
  age05<-L2A(t0,Linf,K,L5,maxage)

  Vmaxage<-runif(nsim,OM@Vmaxlen[1],OM@Vmaxlen[2])#runif(BT_fleet@Vmaxage[1],BT_fleet@Vmaxage[2]) # selectivity of oldest age class

  LM<-runif(nsim,OM@L50[1],OM@L50[2])
  AM<-L2A(t0,Linf,K,LM,maxage)

  # age at maturity
  a<-OM@a # length-weight parameter a
  b<-OM@b # length-weight parameter b

  mod<-AFS          # the age at modal (or youngest max) selectivity
  deriv<-getDNvulnS(mod,age05,Vmaxage,maxage,nsim)           # The vulnerability schedule
  vuln<-deriv[[1]]

  Agearray<-array(rep(1:maxage,each=nsim),c(nsim,maxage))
  mat<-1/(1+exp((AM-(Agearray))/(AM*0.1)))  # Maturity at age array

  nyears<-100
  #bootfun<-function(dat,ind)mean(dat[ind])
  #MLo<-boot(MLt,bootfun,nsim)
  #ML<-MLo$t
  out<-CSRA(M,h,Linf,K,t0,AM,a,b,vuln,mat,ML=rep(ML,nsim),NA,NA,maxage,nyears)
  cond<-out[,1]>Dlim[1]&out[,1]<Dlim[2]&out[,2]<2.5 # Stock levels are unlikely to be above 80% unfished, F is unlikely to be above 2.5
  
  if(ploty){
    par(mfrow=c(1,2))
    plot(density(out[cond,1],from=0,adj=0.4),main="Depletion")
    plot(density(out[cond,2],from=0,adj=0.4),main="Fishing mortality rate")
    OM@D<-quantile(out[cond,1],c(0.05,0.95))
  }

  OM
}

# Composition stock reduction analysis 
CSRA<-function(M,h,Linf,K,t0,AM,a,b,vuln,mat,ML,CAL,CAA,maxage,nyears){ 
  nsim<-length(M)
  Dep<-rep(NA,nsim)
  Fm<-rep(NA,nsim)
  for(i in 1:nsim){
    fit<-optimize(CSRAfunc,log(c(0.0001,5)),Mc=M[i],hc=h[i],maxage,nyears,Linfc=Linf[i],Kc=K[i],t0c=t0[i],AMc=AM[i],
                  ac=a,bc=b,vulnc=vuln[i,],matc=mat[i,],MLc=ML[i],CAL=NA,CAA=NA,opt=T)
    
    
    out<-CSRAfunc(fit$minimum,Mc=M[i],hc=h[i],maxage,nyears,Linfc=Linf[i],Kc=K[i],t0c=t0[i],AMc=AM[i],
                  ac=a,bc=b,vulnc=vuln[i,],matc=mat[i,],MLc=ML[i],CAL=NA,CAA=NA,opt=3)
    
    Dep[i]<-out[1]
    Fm[i]<-out[2]
    
    
  }
  cbind(Dep,Fm)
}

# The function that CSRA operates on
CSRAfunc<-function(lnF,Mc,hc,maxage,nyears,AFSc,AFCc,Linfc,Kc,t0c,AMc,ac,bc,vulnc,matc,MLc,CAL,CAA,opt=T,meth="ML"){
  
  Fm<-exp(lnF)
  Fc<-vulnc*Fm
  Lac<-Linfc*(1-exp(-Kc*((1:maxage)-t0c)))
  Wac<-ac*Lac^bc
  N<-exp(-Mc*((1:maxage)-1))
  SSN<-matc*N                                 # Calculate initial spawning stock numbers
  Biomass<-N*Wac
  SSB<-SSN*Wac                               # Calculate spawning stock biomass
  
  B0<-sum(Biomass)
  SSB0<-sum(SSB)
  SSN0<-SSN
  SSBpR<-sum(SSB)                             # Calculate spawning stock biomass per recruit
  SSNpR<-SSN
  Zc<-Fc+Mc
  CN<-array(NA,dim=c(nyears,maxage))
  HR<-rep(0,maxage)
  pen<-0
  for(y in 1:nyears){
    VB<-Biomass*vulnc*exp(-Mc)
    CN[y,]<-N*(1-exp(-Zc))*(Fc/Zc)
    N[2:maxage]<-N[1:(maxage-1)]*exp(-Zc[1:(maxage-1)])         # Total mortality
    N[1]<-(0.8*hc*sum(SSB))/(0.2*SSBpR*(1-hc)+(hc-0.2)*sum(SSB))  # Recruitment assuming regional R0 and stock wide steepness
    Biomass<-N*Wac
    SSN<-N*matc
    SSB<-SSN*Wac
  } # end of year
  
  pred<-sum((CN[nyears,]*Lac))/sum(CN[nyears,])
  fobj<-(pred-MLc)^2 # Currently a least squares estimator. Probably not worth splitting hairs WRT likelihood functions!
  if(opt==1){return(fobj)
  }else{c(sum(SSB)/sum(SSB0),Fm)
  }
}

# Stochastic inverse growth curve used to back-calculate age at first capture from length at first capture
getAFC<-function(t0c,Linfc,Kc,LFC,maxage){ 
  nsim<-length(t0c)
  agev<-c(0.0001,1:maxage)
  agearray<-matrix(rep(agev,each=nsim),nrow=nsim)
  Larray<-Linfc*(1-exp(-Kc*(agearray-t0c)))
  matplot(agev,t(Larray),type='l')
  abline(h=LFC,col="#ff000030",lwd=2)
  AFC<-(log(1-(LFC/Linfc))/-Kc)+t0c
  abline(v=AFC,col="#0000ff30",lwd=2)
  AFC
}  

L2A<-function(t0c,Linfc,Kc,Len,maxage){ 
  nsim<-length(t0c)
  agev<-c(0.0001,1:maxage)
  agearray<-matrix(rep(agev,each=nsim),nrow=nsim)
  Larray<-Linfc*(1-exp(-Kc*(agearray-t0c)))
  matplot(agev,t(Larray),type='l')
  abline(h=Len,col="#ff000030",lwd=2)
  age<-(log(1-(Len/Linfc))/-Kc)+t0c
  abline(v=age,col="#0000ff30",lwd=2)
  age
}  


# Sketch trends in historical fishing mortality --------------------------------
# Takes Fleet object, runs Sketch function for user to specify historical effort
# Then returns Fleet object with Effort objects populated with output of Sketch
ChooseEffort <- function(FleetObj, Years=NULL) {
  nyears <- FleetObj@nyears
  runSketch <- SketchFun(nyears, Years)
  FleetObj@EffYears <- runSketch[,1]
  FleetObj@EffLower <- runSketch[,2]
  FleetObj@EffUpper <- runSketch[,3]
  return(FleetObj)
}

# Helper functions for above
identifyPch <- function(x, y = NULL, n = length(x), pch = 19, ...)
{
    xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
    sel <- rep(FALSE, length(x)); res <- integer(0)
    while(sum(sel) < n) {
        ans <- identify(x[!sel], y[!sel], n = 1, plot = FALSE, ...)
        if(!length(ans)) break
        ans <- which(!sel)[ans]
        points(x[ans], y[ans], pch = pch)
        sel[ans] <- TRUE
		res <- c(res, ans)
    }
    out <- cbind(x[res], y[res])
	out <- out[order(out[,1]),]
	return(out)
}
# function to allow users to pick points off a plotted grid
SketchFun <- function(nyears=NULL, Years=NULL) {
  
  if (length(Years) == 0 & length(nyears) == 0) stop()
  if (length(Years) > 0) { 
    nyears <- length(Years)
	years <- Years
  }	
  if (length(Years) == 0) years <- 1:nyears
  par(mfrow=c(1,1), mar=c(5, 4, 5, 2))
  ys <- seq(from=0, to=1, by=0.05)
  years1 <- seq(from=years[1], by=2, to=max(years))
  years2 <- seq(from=years[2], by=2, to=max(years))
  grd1 <- expand.grid(years1, ys)
  grd2 <- expand.grid(years2, ys)
  grd3 <- expand.grid(years, ys)
  Xs <- grd3[,1]
  Ys <- grd3[,2]

  plot(grd1[,1], grd1[,2], col="black", pch=19, cex=0.2, xlab="Years", ylab="Variable", cex.lab=1.5)
  points(grd2[,1], grd2[,2], col="darkgrey", pch=19, cex=0.2)
  
  mtext(side=3, "Right Click to Finish. Escape to Quit.", xpd=NA, cex=1.25)
  line1 <- "Use mouse to select points on the grid"
  line2 <- "First and last year must be selected."
  line3 <- "Select two points in a single year to represent range of uncertainty"
  par(xpd=TRUE) 
  text(years[1],par("usr")[4]+0.15*(par("usr")[4]-par("usr")[3]), line1,cex=1, pos=4) 
  text(years[1],par("usr")[4]+0.1125*(par("usr")[4]-par("usr")[3]), line2,cex=1, pos=4) 
  text(years[1],par("usr")[4]+0.075*(par("usr")[4]-par("usr")[3]),line3,cex=1, pos=4)  
  par(xpd=FALSE)
  message(line1,"\n", line2, "\n", line3, "\n")
  flush.console()
  
  par()
  out <- NULL
  out <- identifyPch(x=Xs, y=Ys, tolerance=0.1)
  while(is.null(dim(out))) {
    message("Must choose more than one point")
	flush.console()
    out <- identifyPch(x=Xs, y=Ys, tolerance=0.1)
  }
  while(min(out[,1]) != years[1]) {
    message("Choose point(s) for first year (usually 0)")
	flush.console()
	dat <- rbind(out, identifyPch(x=Xs, y=Ys, tolerance=0.1))
	out <- dat[order(dat[,1]),]
  }
  while(max(out[,1]) != years[length(years)]) {
    message("Choose point(s) for last year (nyear)")
	flush.console()
	dat <- rbind(out, identifyPch(x=Xs, y=Ys, tolerance=0.1))
	out <- dat[order(dat[,1]),]
  }
  ord <- order(out[,1], out[,2])
  out <- out[ord,]
  yrs <- unique(out[,1])
  mat <- matrix(NA, nrow=length(yrs), ncol=3)
  mat[,1] <- yrs
  ind <- which(!duplicated(out[,1]))
  mat[,2:3] <- out[ind,2]
  for (X in seq_along(yrs)) {
    chk <- out[,1] %in% yrs[X]
	ind <- range(which(chk))
    if(sum(chk) > 1) {
	  mat[X,2:3] <- out[ind,2] 
	}
  }

  lines(mat[,1], mat[,2])
  lines(mat[,1], mat[,3])

  colnames(mat) <- c("Years", "Lower", "Upper")
  return(mat)
}



# Sketch Historical Selectivity Patterns ---------------------------------------
ChooseSelect <- function(Fleet, Stock=NULL, FstYr=NULL, SelYears=NULL) {
  
  LastYr <- as.numeric(format(Sys.Date(), format="%Y"))
  if (is.null(FstYr)) {
    message("*****************")
    message("Enter first historical year")
    message("Note: *nyears* will be specified from this value")
    FstYr <- as.numeric(readline("First Historical Year: "))
	message("\n")
  }
  if (is.null(SelYears)) {
    message("Enter each selectivity break point year, seperated by a comma")
    message("Note: break points are the years where selectivity pattern changed")
    message("Note: break points must be within year range")
    inString <- readline("Enter each selectivity break point year, seperated by a comma: ")
    options(warn=-1)
    SelYears <- as.numeric(unlist(strsplit(inString, ",")))
    if (is.na(SelYears))  SelYears <- as.numeric(unlist(strsplit(inString, " ")))
    options(warn=0)
  }	
  SelYears <- sort(SelYears)
  if (SelYears[1] != FstYr) SelYears <- c(FstYr, SelYears)
  message("Break Points Years are: ")
  print(SelYears)
  if (length(SelYears) < 2) stop("Must be more than one year")
  if (max(SelYears) > LastYr) stop("Must specify historical year")
  if (min(SelYears) < FstYr) stop("Year before first year")
  flush.console()
  Selnyears <- length(SelYears)
 
  Years <- FstYr:LastYr #SelYears[1]:SelYears[length(SelYears)]
  Fleet@nyears <- length(Years)
  ind <- round((Range(SelYears, Max=LastYr, Min=FstYr)) * Fleet@nyears,0) +1
  ind[length(ind)] <- max(ind) - 1  
  Fleet@AbsSelYears <- SelYears 
  Fleet@SelYears <- ind
  
  Fleet@L5 <- matrix(0, nrow=Selnyears, ncol=2)
  Fleet@LFS <- matrix(0, nrow=Selnyears, ncol=2)
  Fleet@Vmaxlen <- matrix(0, nrow=Selnyears, ncol=2)
  
  # if(is.null(Stock)) Stock <- NA
  set.par <- par(no.readonly=TRUE)
  message("Select selectivity points on plot")
  flush.console()
  for (N in 1:Selnyears) {
    BlankPlot(Stock=Stock, Yr=SelYears[N], N=N)
    L5Out <- ChooseL5()
    Fleet@L5[N,] <- sort(L5Out[,1])
    LFSout <- ChooseLFS(L5Out)
    Fleet@LFS[N,] <- sort(LFSout[,1])
    Vmaxout <- ChooseVmaxlen()
    Fleet@Vmaxlen[N,] <- sort(Vmaxout[,2])
    polygon(x=c(0, max(Fleet@L5[N,]), max(Fleet@LFS[N,]), 3, 
     rev(c(0, min(Fleet@L5[N,]), min(Fleet@LFS[N,]), 3))),
	 y= c(0, 0.05, 1, min(Fleet@Vmaxlen[N,]),
	 rev(c(0, 0.05, 1, max(Fleet@Vmaxlen[N,])))), col="grey")
    par(ask=TRUE)
  }	
  par(set.par)
  ChckSelect(Fleet, Stock)
  Fleet
}

# Rough Plot of Historical Selectivity Patterns --------------------------------
ChckSelect <- function(Fleet, Stock=NULL) {
 if (length(Fleet@SelYears) < 1) stop("No break points in selectivity pattern")
 n <- length(Fleet@SelYears)
 if (n < 4) par(mfrow=c(n, 1), mar=c(4,4,1,1), oma=c(2,3,1,1), bty="l")
 if(n >= 4) {
   N <- ceiling(n/2)
   par(mfrow=c(N, 2), , mar=c(4,4,1,1), oma=c(2,3,1,1), bty="l")
 }
 for (X in 1:n) {
   plot(c(0,3), c(0,1), type="n", xlab="", ylab="")
   if(length(Fleet@AbsSelYears) > 0) title(Fleet@AbsSelYears[X])
   if(length(Fleet@AbsSelYears) == 0) title(Fleet@SelYears[X])
   polygon(x=c(0, max(Fleet@L5[X,]), max(Fleet@LFS[X,]), 3, 
     rev(c(0, min(Fleet@L5[X,]), min(Fleet@LFS[X,]), 3))),
	 y= c(0, 0.05, 1, min(Fleet@Vmaxlen[X,]),
	 rev(c(0, 0.05, 1, max(Fleet@Vmaxlen[X,])))), col="grey")
  lines(c(1,1), c(0,1), lty=3)
  text(1.1, 0.2, "L50", cex=1.25)
}
  mtext(side=2, outer=TRUE, "Selectivity", cex=1.25, xpd=NA)
  mtext(side=1, outer=TRUE, "Relative Length", cex=1.25, xpd=NA)
}
  
# Supporting
# BlankPlot
BlankPlot <- function(Stock=NULL, Yr=NULL, N=NULL) {
  Max <- 3 
  AxCex <- 1.3
  By <- 0.05 
  par(mfrow=c(1,1), mai=c(2, 1, .5, .3), oma=c(1,1,1,1))
  plot(c(0,3), c(0,1), type="n", xlab="", ylab="", axes=FALSE)
  mtext(side=2, line=3, "Selectivity", cex=AxCex)
  axis(side=2)
  Xax <- seq(from=0, to=Max-By, by=2*By)
  axis(side=1, at=Xax)
  axis(side=1, at=c(2, Max), label=c("", "Lmax"), xpd=NA)
  mtext(side=1, line=3.5, "Relative Length", cex=AxCex)
  axis(side=1, at=1, line=1.5, label="L50")
  if (!is.null(Stock) & class(Stock) == "Stock") {
    L50 <- mean(Stock@L50) # mean length at maturity
    MatAx <- L50 * Xax
    axis(side=1, line=5.5, at=Xax, labels=MatAx)
    axis(side=1, line=5.5, at=c(2, Max), label=c("", "Lmax"), xpd=NA)
    mtext(side=1, line=8.5, "Approx. Length", cex=AxCex)
    axis(side=1, at=1, line=6.5, label="Mean L50")
  } 
  if (N == 1) {
    title(paste("Choose selectivity points for Year", Yr, "(First Year)"))
  } else {
    title(paste("Choose selectivity points for Year", Yr))
  }	
}
# Choose L5 
ChooseL5 <- function() {
  By <- 0.05
  Xs <-seq(from=0, to=1, by=By)
  Ys <- rep(0.05, length(Xs))
  points(Xs, Ys, col="gray", cex=0.5)
  text(0.5, 0.2, "Choose two points for L5")
  L5out <- identifyPch(x=Xs, y=Ys, tolerance=0.1, n=2)
  L5out
}
# Choose LFS 
ChooseLFS <- function(L5out) {
  Max <- 3
  By <- 0.05
  Xs <-seq(from=max(L5out[,1]), to=Max, by=By)
  Ys <- rep(1, length(Xs))
  points(Xs, Ys, col="gray", cex=0.5)
  text(1.5, 0.5, "Choose two points for LFS")
  LFSout <- identifyPch(x=Xs, y=Ys, tolerance=0.1, n=2)
  LFSout
}
# Choose Vmaxlen
ChooseVmaxlen <- function(L5out) {
  Max <- 3 
  Ys <- seq(from=0, to=1, by=0.05) 
  Xs <- rep(Max, length(Ys))
  points(Xs, Ys, col="gray", cex=0.5)
  text(2, 0.8, "Choose two points for selectivity\n at maximum length")
  Vmaxout <- identifyPch(x=Xs, y=Ys, tolerance=0.1, n=2)
  Vmaxout
}
	









