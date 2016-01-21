# # MSE code

# # Define classes for operating models
# setClass("OM",representation(Name="character",nyears="numeric",maxage="numeric",R0="numeric",M="numeric",
                              # Msd="numeric",Mgrad="numeric",h="numeric",SRrel="numeric",Linf="numeric",K="numeric",t0="numeric",
                              # Ksd="numeric",Kgrad="numeric",Linfsd="numeric",Linfgrad="numeric",recgrad="numeric",
                              # a="numeric",b="numeric",D="numeric",
                              # Size_area_1="numeric",Frac_area_1="numeric",Prob_staying="numeric",
                              # Source="character", L50="numeric", L50_95="numeric", L5="numeric", LFS="numeric",
							  # Vmaxlen="numeric", beta="numeric", Spat_targ="numeric", Fsd="numeric",
							  # EffYears="numeric", EffLower="numeric", EffUpper="numeric", 
							  # # Fgrad="numeric", 
							  # qinc="numeric",qcv="numeric",AC="numeric", Cobs="numeric",Cbiascv="numeric",
							  # CAA_nsamp="numeric",CAA_ESS="numeric",
                              # CAL_nsamp="numeric",CAL_ESS="numeric",CALcv="numeric",
                              # Iobs="numeric",Perr="numeric",
                              # Mcv="numeric",Kcv="numeric",t0cv="numeric",Linfcv="numeric",
                              # LFCcv="numeric",
                              # LFScv="numeric",
                              # B0cv="numeric",FMSYcv="numeric",FMSY_Mcv="numeric",BMSY_B0cv="numeric",
                              # LenMcv="numeric",rcv="numeric",
                              # Dbiascv="numeric",Dcv="numeric",Btbias="numeric",Btcv="numeric",
                              # Fcurbiascv="numeric",Fcurcv="numeric",hcv="numeric",
                              # Icv="numeric",maxagecv="numeric",
                              # Reccv="numeric",Irefcv="numeric",Crefcv="numeric",Brefcv="numeric"))


# setMethod("initialize", "OM", function(.Object,Stock,Fleet,Observation){

  # if(class(Stock)!='Stock')print(paste('Could not build operating model:',deparse(substitute(Stock)),'not of class Stock'))
  # if(class(Fleet)!='Fleet')print(paste('Could not build operating model:',deparse(substitute(Fleet)),'not of class Fleet'))
  # if(class(Observation)!='Observation')print(paste('Could not build operating model:',deparse(substitute(Observation)),'not of class Observation'))
  # if(class(Stock)!='Stock'|class(Fleet)!='Fleet'|class(Observation)!='Observation')stop()

  # .Object@Name<-paste("Stock:",Stock@Name,"  Fleet:",Fleet@Name,"  Observation model:",Observation@Name,sep="")
  # # Now copy the values for stock, fleet and observation slots to same slots in the Sim object
  # Sslots<-slotNames(Stock)
  # for(i in 2:length(Sslots))slot(.Object,Sslots[i])<-slot(Stock,Sslots[i])
  # Fslots<-slotNames(Fleet)
  # for(i in 2:length(Fslots))slot(.Object,Fslots[i])<-slot(Fleet,Fslots[i])
  # Oslots<-slotNames(Observation)
  # for(i in 2:length(Oslots))slot(.Object,Oslots[i])<-slot(Observation,Oslots[i])
  # .Object

# })

# setClass("Stock",representation(Name="character",maxage="numeric",R0="numeric",M="numeric",
                # Msd="numeric",Mgrad="numeric",h="numeric",SRrel="numeric",Linf="numeric",K="numeric",t0="numeric",
                # Ksd="numeric",Kgrad="numeric",Linfsd="numeric",Linfgrad="numeric",recgrad="numeric",
                # a="numeric",b="numeric",D="numeric",Perr="numeric", Size_area_1="numeric",Frac_area_1="numeric",
				        # Prob_staying="numeric",AC="numeric",
				        # L50="numeric", L50_95="numeric",Source="character"))


# setMethod("initialize", "Stock", function(.Object,file=NA){
  # if (!is.na(file)) {
    # if (file.exists(file)) {
      # dat <- read.csv(file,header=F,colClasses="character") # read 1st sheet
      # dname<-dat[,1]
      # dat<-dat[,2:ncol(dat)]
      
      # .Object@Name<-dat[match("Name", dname),1]
      # .Object@maxage<-as.numeric(dat[match("maxage",dname),1])
      # .Object@R0<-as.numeric(dat[match("R0",dname),1])
      # .Object@M<-as.numeric(dat[match("M",dname),1:2])
      # .Object@Msd<-as.numeric(dat[match("Msd",dname),1:2])
      # .Object@Mgrad<-as.numeric(dat[match("Mgrad",dname),1:2])
      # .Object@h<-as.numeric(dat[match("h",dname),1:2])
      # .Object@SRrel<-as.numeric(dat[match("SRrel",dname),1])
      # .Object@Linf<-as.numeric(dat[match("Linf",dname),1:2])
      # .Object@K<-as.numeric(dat[match("K",dname),1:2])
      # .Object@t0<-as.numeric(dat[match("t0",dname),1:2])
      # .Object@Ksd<-as.numeric(dat[match("Ksd",dname),1:2])
      # .Object@Kgrad<-as.numeric(dat[match("Kgrad",dname),1:2])
      # .Object@Linfsd<-as.numeric(dat[match("Linfsd",dname),1:2])
      # .Object@Linfgrad<-as.numeric(dat[match("Linfgrad",dname),1:2])
      # .Object@recgrad<-as.numeric(dat[match("recgrad",dname),1:2])
      # .Object@a<-as.numeric(dat[match("a",dname),1])
      # .Object@b<-as.numeric(dat[match("b",dname),1])
      # .Object@D<-as.numeric(dat[match("D",dname),1:2])
      # .Object@Perr<-as.numeric(dat[match("Perr",dname),1:2])
      # .Object@AC<-as.numeric(dat[match("AC",dname),1:2])
      # .Object@Size_area_1<-as.numeric(dat[match("Size_area_1",dname),1:2])
      # .Object@Frac_area_1<-as.numeric(dat[match("Frac_area_1",dname),1:2])
      # .Object@Prob_staying<-as.numeric(dat[match("Prob_staying",dname),1:2])
      # .Object@L50<-as.numeric(dat[match("L50", dname),1:2])
      # .Object@L50_95<-as.numeric(dat[match("L50_95", dname),1:2])
      # .Object@Source<-dat[match("Source", dname),1]
	# } else {
	  # message("File doesn't exist")
	# }
   # }	
  # .Object

# })

# setClass("Fleet",representation(Name="character",nyears="numeric", Spat_targ="numeric",
       	 # Fsd="numeric", qinc="numeric",qcv="numeric", 
		 # EffYears="numeric", EffLower="numeric", EffUpper="numeric",
		 # L5="numeric", LFS="numeric",Vmaxlen="numeric"))

# # Added EffYears etc 
# setMethod("initialize", "Fleet", function(.Object,file=NA){
  # if (!is.na(file)) {
    # if (file.exists(file)) {
      # dat <- read.csv(file,header=F,colClasses="character") # read 1st sheet
      # dname<-dat[,1]
      # dat<-dat[,2:ncol(dat)]
      
      # .Object@Name<-dat[match("Name", dname),1]
      # .Object@nyears <-as.numeric(dat[match("nyears",dname),1])
      # .Object@Spat_targ<-as.numeric(dat[match("Spat_targ",dname),1:2])
      # .Object@Fsd<-as.numeric(dat[match("Fsd",dname),1:2])
      # # .Object@Fgrad<-as.numeric(dat[match("Fgrad",dname),1:2])
      # nEffYears <- ncol(dat[match("EffYears",dname),])
	  # chk <- as.numeric(dat[match("EffYears",dname),1:nEffYears])
	  # ind <- which(!is.na(chk))
	  # nEffYears <- length(ind)
	  
      # .Object@EffYears <-as.numeric(dat[match("EffYears",dname),1:nEffYears])
      # .Object@EffLower <-as.numeric(dat[match("EffLower",dname),1:nEffYears])
      # .Object@EffUpper <-as.numeric(dat[match("EffUpper",dname),1:nEffYears])
      
      # .Object@qinc<-as.numeric(dat[match("qinc",dname),1:2])
      # .Object@qcv<-as.numeric(dat[match("qcv",dname),1:2])
      # .Object@L5<-as.numeric(dat[match("L5",dname),1:2])
      # .Object@LFS<-as.numeric(dat[match("LFS",dname),1:2])
      # .Object@Vmaxlen<-as.numeric(dat[match("Vmaxlen",dname),1:2])
    # } else {
	  # message("File doesn't exist")
	# }
   # }	
 	
  # .Object
# })


# setClass("lmmodel",representation(Name="character",models="list"))

# setMethod("initialize", "lmmodel", function(.Object,Name,models){
  
  # .Object@Name<-Name
  # .Object@models<-models
  # .Object
  
# })



# setClass("Observation",representation(Name="character",LenMcv="numeric",
                # Cobs="numeric",Cbiascv="numeric",CAA_nsamp="numeric",CAA_ESS="numeric",
                # CAL_nsamp="numeric",CAL_ESS="numeric",CALcv="numeric",
                # Iobs="numeric",Mcv="numeric",Kcv="numeric",t0cv="numeric",Linfcv="numeric",
                # LFCcv="numeric",LFScv="numeric",B0cv="numeric",
                # FMSYcv="numeric",FMSY_Mcv="numeric",BMSY_B0cv="numeric",
                # rcv="numeric", Dbiascv="numeric",Dcv="numeric",
                # Btbias="numeric",Btcv="numeric",Fcurbiascv="numeric",Fcurcv="numeric",
                # hcv="numeric",Icv="numeric",maxagecv="numeric",Reccv="numeric",
                # Irefcv="numeric",Crefcv="numeric",Brefcv="numeric",beta="numeric"))


# setMethod("initialize", "Observation", function(.Object,file=NA){
  # if (!is.na(file)) {
    # if (file.exists(file)) {
      # dat <- read.csv(file,header=F,colClasses="character") # read 1st sheet
      # dname<-dat[,1]
      # dat<-dat[,2:ncol(dat)]
      
      # .Object@Name<-dat[match("Name", dname),1]
      # .Object@LenMcv<-as.numeric(dat[match("LenMcv",dname),1])
      # .Object@Cobs<-as.numeric(dat[match("Cobs",dname),1:2])
      # .Object@Cbiascv<-as.numeric(dat[match("Cbiascv",dname),1])
      # .Object@CAA_nsamp<-as.numeric(dat[match("CAA_nsamp",dname),1:2])
      # .Object@CAA_ESS<-as.numeric(dat[match("CAA_ESS",dname),1:2])
      # .Object@CAL_nsamp<-as.numeric(dat[match("CAA_nsamp",dname),1:2])
      # .Object@CAL_ESS<-as.numeric(dat[match("CAA_ESS",dname),1:2])
      # .Object@CALcv<-as.numeric(dat[match("CALcv",dname),1:2])
      # .Object@Iobs<-as.numeric(dat[match("Iobs",dname),1:2])
      # .Object@Mcv<-as.numeric(dat[match("Mcv",dname),1])
      # .Object@Kcv<-as.numeric(dat[match("Kcv",dname),1])
      # .Object@t0cv<-as.numeric(dat[match("t0cv",dname),1])
      # .Object@Linfcv<-as.numeric(dat[match("Linfcv",dname),1])
      # .Object@LFCcv<-as.numeric(dat[match("LFCcv",dname),1])
      # .Object@LFScv<-as.numeric(dat[match("LFScv",dname),1])
      # .Object@B0cv<-as.numeric(dat[match("B0cv",dname),1])
      # .Object@FMSYcv<-as.numeric(dat[match("FMSYcv",dname),1])
      # .Object@FMSY_Mcv<-as.numeric(dat[match("FMSY_Mcv",dname),1])
      # .Object@BMSY_B0cv<-as.numeric(dat[match("BMSY_B0cv",dname),1])
      # .Object@rcv<-as.numeric(dat[match("rcv",dname),1])
      # .Object@Dbiascv<-as.numeric(dat[match("Dbiascv",dname),1])
      # .Object@Dcv<-as.numeric(dat[match("Dcv",dname),1:2])
      # .Object@Btbias<-as.numeric(dat[match("Btbias",dname),1:2])
      # .Object@Btcv<-as.numeric(dat[match("Btcv",dname),1:2])
      # .Object@Fcurbiascv<-as.numeric(dat[match("Fcurbiascv",dname),1])
      # .Object@Fcurcv<-as.numeric(dat[match("Fcurcv",dname),1:2])
      # .Object@hcv<-as.numeric(dat[match("hcv",dname),1])
      # .Object@Icv<-as.numeric(dat[match("Icv",dname),1])
      # .Object@maxagecv<-as.numeric(dat[match("maxagecv",dname),1])
      # .Object@Reccv<-as.numeric(dat[match("Reccv",dname),1:2])
      # .Object@Irefcv<-as.numeric(dat[match("Irefcv",dname),1])
      # .Object@Crefcv<-as.numeric(dat[match("Crefcv",dname),1])
      # .Object@Brefcv<-as.numeric(dat[match("Brefcv",dname),1])
      # .Object@beta<-as.numeric(dat[match("beta",dname),1:2])
	# } else {
	  # message("File doesn't exist")
	# }
  # }	
  # .Object

# })


# setClass("MSE",representation(Name="character",nyears="numeric",proyears="numeric",nMPs="numeric",MPs="character",
                              # nsim="numeric",OM="data.frame",Obs="data.frame",B_BMSY="array",
                              # F_FMSY="array",B="array",FM="array",C="array",TAC="array",SSB_hist="array",CB_hist="array",FM_hist="array"))

# setMethod("initialize", "MSE", function(.Object,Name,nyears,proyears,nMPs,MPs,
                                                # nsim,OMtable,Obs,B_BMSYa,F_FMSYa,Ba,FMa,Ca,TACa,SSB_hist,CB_hist,FM_hist){
  # .Object@Name<-Name
  # .Object@nyears <-nyears
  # .Object@proyears<-proyears
  # .Object@nMPs<-nMPs
  # .Object@MPs<-MPs
  # .Object@nsim<-nsim
  # .Object@OM<-OMtable
  # .Object@Obs<-Obs
  # .Object@B_BMSY<-B_BMSYa
  # .Object@F_FMSY<-F_FMSYa
  # .Object@B<-Ba
  # .Object@FM<-FMa
  # .Object@C<-Ca
  # .Object@TAC<-TACa
  # .Object@SSB_hist<-SSB_hist
  # .Object@CB_hist<-CB_hist
  # .Object@FM_hist<-FM_hist
  # .Object
# })

# sampy<-function(x) sample(x,1,prob=!is.na(x))



# runMSE <- function(OM="1", MPs=NA, nsim=48, proyears=28, interval=4, pstar=0.5,
                   # maxF=0.8, timelimit=1, reps=1, custompars=0){
  
  # print("Loading operating model")
  # flush.console()
  # if(class(OM)!="OM")stop("You must specify an operating model")
  # #if(!sfIsRunning())stop("You must initialize snowfall functions sfInit() see ??DLMtool")
  
  # nyears <- OM@nyears  # number of  historical years
  # maxage <- OM@maxage  # maximum age (no plus group)
  
  # dep <- runif(nsim,OM@D[1],OM@D[2])  # sample from the range of user-specified depletion (Bcurrent/B0)
  # Esd <- runif(nsim,OM@Fsd[1],OM@Fsd[2]) # interannual variability in fishing effort (log normal sd)
  
  # Deriv <- getEffhist(Esd, nyears, EffYears=OM@EffYears, EffLower=OM@EffLower, EffUpper=OM@EffUpper) # Historical fishing effort
  # Find <- Deriv[[1]] # Calculate fishing effort rate
  # dFfinal <- Deriv[[2]] # Final gradient in fishing effort yr-1 
  
  # dep[order(dFfinal)]<-dep[order(dep,decreasing=T)] # robustifies 
   
  # # Deriv <- getFhist(nsim,Esd,nyears,dFmin=OM@Fgrad[1],dFmax=OM@Fgrad[2],bb=0.5)     # Calculate fishing effort rate
  # # Find <- Deriv[[1]]     # Calculate fishing effort rate
  # # dFfinal <- Deriv[[2]]  # Final gradient in fishing effort yr-1 
 
  # # Sample operating model parameters ===========================================================
  # procsd<-runif(nsim,OM@Perr[1],OM@Perr[2])         # Process error standard deviation
  # AC<-runif(nsim,OM@AC[1],OM@AC[2])    # auto correlation parameter for recruitment deviations recdev(t)<-AC*recdev(t-1)+(1-AC)*recdev_proposed(t)
  
  # M<-runif(nsim,OM@M[1],OM@M[2]) # natural mortality rate 	
  # Msd<-runif(nsim,OM@Msd[1],OM@Msd[2]) # sample inter annual variability in M from specified range
  # Mgrad<-runif(nsim,OM@Mgrad[1],OM@Mgrad[2]) # sample gradient in M (M y-1)
  # hs<-runif(nsim,OM@h[1],OM@h[2]) # sample of recruitment compensation (steepness - fraction of unfished recruitment at 20% of unfished biomass)
  # Linf<-runif(nsim,OM@Linf[1],OM@Linf[2]) # sample of asymptotic length
  # Linfsd<-runif(nsim,OM@Linfsd[1],OM@Linfsd[2]) # sample of interannual variability in Linf
  # Linfgrad<-runif(nsim,OM@Linfgrad[1],OM@Linfgrad[2]) # sample of gradient in Linf (Linf y-1)
  # recgrad<-runif(nsim,OM@recgrad[1],OM@recgrad[2]) # gradient in recent recruitment
  # K<-runif(nsim,OM@K[1],OM@K[2])  # now predicted by a log-linear model
  # Ksd<-runif(nsim,OM@Ksd[1],OM@Ksd[2])#runif(nsim,OM@Ksd[1],OM@Ksd[2])# sd is already added in the linear model prediction
  # Kgrad<-runif(nsim,OM@Kgrad[1],OM@Kgrad[2]) # gradient in Von-B K parameter (K y-1)
  # t0<-runif(nsim,OM@t0[1],OM@t0[2]) # a sample of theoretical age at length zero
  # lenM <- array(runif(nsim*50,OM@L50[1],OM@L50[2]),c(nsim,50)) # length at 50% maturity
  # lenM[lenM/Linf>0.8]<-NA
  # lenM<-apply(lenM,1,function(x)x[!is.na(x)][1])
  # len95 <- array(lenM + runif(nsim*50,OM@L50_95[1],OM@L50_95[2]),c(nsim,50)) # length at 95% maturity
  # len95[len95/Linf>0.9]<-NA
  # len95<-apply(len95,1,function(x)x[!is.na(x)][1])
  
  # if (max(OM@L5) > 1) {
    # message("L5 set too high (maximum value of 1). \nDefaulting to L5 = 1, with no variability")
    # OM@L5[OM@L5 > 1] <- 1 
  # }
  
  # L5 <- runif(nsim, OM@L5[1], OM@L5[2]) * lenM        # length at 0.05% selectivity ascending
  # LFS <- runif(nsim, OM@LFS[1], OM@LFS[2]) * lenM     # first length at 100% selection
  # # LFS <- array(runif(nsim*50, OM@LFS[1], OM@LFS[2]) * lenM,dim=c(nsim,50))     # first length at 100% selection
 
  # ind <- which(LFS/Linf>1, arr.ind=T)
  # if (length(ind) > 0) {
    # message("LFS too high (LFS > Linf) in some cases. \nDefaulting to LFS = Linf for the affected simulations")
    # LFS[ind] <- Linf[ind]
  # }
   
  # # LFS[LFS/Linf>1]<-NA
  # # LFS<-apply(LFS,1,function(x)x[!is.na(x)][1])
  
  # Vmaxlen <- runif(nsim,OM@Vmaxlen[1],OM@Vmaxlen[2])   # selectivity at maximum length (change variable name later if needed)
  # Spat_targ<-runif(nsim,OM@Spat_targ[1],OM@Spat_targ[2]) # spatial targetting Ba^targetting param 
  # Frac_area_1<-runif(nsim,OM@Frac_area_1[1],OM@Frac_area_1[2]) # sampled fraction of unfished biomass in area 1 (its a two area model by default)
  # Prob_staying<-runif(nsim,OM@Prob_staying[1],OM@Prob_staying[2]) # sampled probability of individuals staying in area 1 among years
  # Size_area_1<-runif(nsim,OM@Size_area_1[1],OM@Size_area_1[2]) # currently redundant parameter for the habitat area size of area 1
  
  # # Sample observation error model paramters ===============================================================
  # Csd<-runif(nsim,OM@Cobs[1],OM@Cobs[2])                                    # Sampled catch observation error (lognormal sd)
  # Cbias<-rlnorm(nsim,mconv(1,OM@Cbiascv),sdconv(1,OM@Cbiascv))              # Sampled catch bias (log normal sd)
  # CAA_nsamp<-ceiling(runif(nsim,OM@CAA_nsamp[1],OM@CAA_nsamp[2]))                     # Number of catch-at-age observations
  # CAA_ESS<-ceiling(runif(nsim,OM@CAA_ESS[1],OM@CAA_ESS[2]))                                          # Effective sample size
  # CAL_nsamp<-runif(nsim,OM@CAL_nsamp[1],OM@CAL_nsamp[2])                              # Observation error standard deviation for single catch at age by area
  # CAL_ESS<-ceiling(runif(nsim,OM@CAL_ESS[1],OM@CAL_ESS[2]))                                          # Effective sample size
  # CALcv<-runif(nsim,OM@CALcv[1],OM@CALcv[2])                                    # Observation error standard deviation for single catch at age by area
  # betas<-exp(runif(nsim,log(OM@beta[1]),log(OM@beta[2]))) # the sampled hyperstability / hyperdepletion parameter beta>1 (hyperdepletion) beta<1 (hyperstability)
  # Isd<-runif(nsim,OM@Iobs[1],OM@Iobs[2])                    # Abundance index observation error (log normal sd)
  # Derr<-runif(nsim,OM@Dcv[1],OM@Dcv[2])
  # Dbias<-rlnorm(nsim,mconv(1,OM@Dbiascv),sdconv(1,OM@Dbiascv)) # sample of depletion bias
  # Mbias<-rlnorm(nsim,mconv(1,OM@Mcv),sdconv(1,OM@Mcv))         # sample of M bias
  # FMSY_Mbias<-rlnorm(nsim,mconv(1,OM@FMSY_Mcv),sdconv(1,OM@FMSY_Mcv)) # sample of FMSY/M bias
  # ntest<-20                               # number of trials
  
  # lenMbias<-rlnorm(nsim,mconv(1,OM@LenMcv),sdconv(1,OM@LenMcv))      # sample of length at maturity bias - assume same error as age based maturity
  # LFCbias<-rlnorm(nsim,mconv(1,OM@LFCcv),sdconv(1,OM@LFCcv))        # sample of length at first capture bias
  # LFSbias<-rlnorm(nsim,mconv(1,OM@LFScv),sdconv(1,OM@LFScv))        # sample of length at full selection bias
  # Aerr<-runif(nsim,OM@Btcv[1],OM@Btcv[2])
  # Abias<-exp(runif(nsim,log(OM@Btbias[1]),log(OM@Btbias[2])))#rlnorm(nsim,mconv(1,OM@Btbiascv),sdconv(1,OM@Btbiascv))    # smaple of current abundance bias
  # Kbias<-rlnorm(nsim,mconv(1,OM@Kcv),sdconv(1,OM@Kcv))              # sample of von B. K parameter bias
  # t0bias<-rlnorm(nsim,mconv(1,OM@t0cv),sdconv(1,OM@t0cv))           # sample of von B. t0 parameter bias
  # Linfbias<-rlnorm(nsim,mconv(1,OM@Linfcv),sdconv(1,OM@Linfcv))     # sample of von B. maximum length bias
  # Irefbias<-rlnorm(nsim,mconv(1,OM@Irefcv),sdconv(1,OM@Irefcv))     # sample of bias in reference (target) abundance index
  # Crefbias<-rlnorm(nsim,mconv(1,OM@Crefcv),sdconv(1,OM@Crefcv))     # sample of bias in reference (target) catch index
  # Brefbias<-rlnorm(nsim,mconv(1,OM@Brefcv),sdconv(1,OM@Brefcv))     # sample of bias in reference (target) biomass index
  # Recsd<-runif(nsim,OM@Reccv[1],OM@Reccv[2])                        # Recruitment deviation 
  
  # # Sample fishing efficiency parameters =======================================================
  # qinc<-runif(nsim,OM@qinc[1],OM@qinc[2])
  # qcv<-runif(nsim,OM@qcv[1],OM@qcv[2])                 # interannual variability in catchability
  
 # # dat<-as.data.frame(cbind(procsd,AC,M,Msd,Mgrad,hs,Linf,Linfsd,Linfgrad,recgrad,K,Ksd,Kgrad,t0,lenM,len95,L5,LFS,
  # #           Vmaxlen,Spat_targ,Frac_area_1,Prob_staying,Size_area_1,Csd,Cbias,CAA_nsamp,CAA_ESS,CALcv,betas,
  # #           Isd,Derr,Dbias,Mbias,FMSY_Mbias,lenMbias,LFCbias,LFSbias,Aerr,Abias,Kbias,
  # #           t0bias,Linfbias,Irefbias,Crefbias,Brefbias,Recsd,qinc,qcv))
  
  # #save(dat,file="F:/DLM/Operating models/Other/custompars")
 
  # # Sample custom parameters ===================================================================
  # if(sum(custompars)!=0){
    # if(nrow(custompars)<nsim){
      # ind<-sample(nrow(custompars),nsim,replace=T)
    # }else{
      # ind<-sample(nrow(custompars),nsim,replace=F)
    # }
    # for(i in 1:ncol(custompars))assign(names(custompars)[i],custompars[ind,i])
  # }

  # procmu<--0.5*(procsd)^2 # adjusted log normal mean
  # Perr<-array(rnorm((nyears+proyears)*nsim,rep(procmu,nyears+proyears),rep(procsd,nyears+proyears)),c(nsim,nyears+proyears))
  # for(y in 2:(nyears+proyears))Perr[,y]<-AC*Perr[,y-1]+Perr[,y]*(1-AC*AC)^0.5#2#AC*Perr[,y-1]+(1-AC)*Perr[,y] # apply a pseudo AR1 autocorrelation to rec devs (log space)
  # Perr <-exp(Perr) # normal space (mean 1 on average)
  # R0 <- OM@R0 # Initial recruitment
  
  # Marray<-gettempvar(M,Msd,Mgrad,nyears+proyears,nsim) # M by sim and year according to gradient and inter annual variability
  # SRrel<-rep(OM@SRrel,nsim) # type of Stock-recruit relationship. 1=Beverton Holt, 2=Ricker
  # Linfarray<-gettempvar(Linf,Linfsd,Linfgrad,nyears+proyears,nsim) # Linf array
  # Karray<-gettempvar(K,Ksd,Kgrad,nyears+proyears,nsim) # the K array
  
  # Agearray<-array(rep(1:maxage,each=nsim),dim=c(nsim,maxage))   # Age array
  # Len_age<-array(NA,dim=c(nsim,maxage,nyears+proyears)) # Length at age array
  # ind<-as.matrix(expand.grid(1:nsim,1:maxage,1:(nyears+proyears))) # an index for calculating Length at age
  # Len_age[ind]<-Linfarray[ind[,c(1,3)]]*(1-exp(-Karray[ind[,c(1,3)]]*(Agearray[ind[,1:2]]-t0[ind[,1]])))
  # Wt_age<-array(NA,dim=c(nsim,maxage,nyears+proyears)) # Weight at age array
  # Wt_age[ind] <- OM@a*Len_age[ind]^OM@b                  # Calculation of weight array
  
  # ageM <- -((log(1-lenM/Linf))/K) + t0 # calculate ageM from L50 and growth parameters (non-time-varying)
  # ageM[ageM < 1] <- 1 # age at maturity must be at least 1 
  # age95 <- -((log(1-len95/Linf))/K) + t0 	
  # age95[age95 < 1] <- 1.5 # must be greater than 0 and ageM
  
  # ageMsd <- sapply(1:nsim,getroot,ageM,age95)
  # ageMarray <- array(ageM,dim=c(nsim,maxage)) # Age at maturity array
  # Mat_age <- 1/(1+exp((ageMarray-(Agearray))/(ageMarray*ageMsd)))  # Maturity at age array
  # mod <- -((log(1-LFS/Linf))/K) + t0 # the age at modal (or youngest max) selectivity
  # age05 <- -((log(1-L5/Linf))/K) + t0# the highest age at %5 selectivity
  
  # V <- array(NA, dim=c(nsim, maxage, nyears+proyears))
   # for (Yr in 1:(nyears+proyears)) { # selectivity pattern for all years (possibly updated in projection)
   # V[,,Yr] <- t(sapply(1:nsim, SelectFun, L5, LFS, Vmaxlen, Linfs=Linfarray[,Yr], Lens=Len_age[,,Yr]))
  # }
  
  # Asize<-cbind(Size_area_1,1-Size_area_1)
  
   
  # print("Optimizing for user-specified movement")  # Print a progress update
  # flush.console()                                  # refresh the console
  
  # if(sfIsRunning()){ # if the cluster is initiated 
    # sfExport(list=c("Frac_area_1","Prob_staying")) # export some of the new arrays and ...
    # mov<-array(t(sfSapply(1:nsim,getmov,Frac_area_1=Frac_area_1,Prob_staying=Prob_staying)),dim=c(nsim,2,2)) # numerically determine movement probability parameters to match Prob_staying and Frac_area_1
  # }else{ # no cluster initiated
    # mov<-array(t(sapply(1:nsim,getmov,Frac_area_1=Frac_area_1,Prob_staying=Prob_staying)),dim=c(nsim,2,2)) # numerically determine movement probability parameters to match Prob_staying and Frac_area_1
  # }
  
  # nareas<-2  # default is a two area model
  # N<-array(NA,dim=c(nsim,maxage,nyears,nareas))        # stock numbers array
  # Biomass<-array(NA,dim=c(nsim,maxage,nyears,nareas))  # stock biomass array
  # VBiomass<-array(NA,dim=c(nsim,maxage,nyears,nareas)) # vulnerable biomass array
  
  # SSN<-array(NA,dim=c(nsim,maxage,nyears,nareas)) # spawning stock numbers array
  # SSB<-array(NA,dim=c(nsim,maxage,nyears,nareas)) # spawning stock biomass array
  # FM<-array(NA,dim=c(nsim,maxage,nyears,nareas))  # fishing mortality rate array
  # Z<-array(NA,dim=c(nsim,maxage,nyears,nareas))   # total mortlaity rate array
  
  # Agearray<-array(rep(1:maxage,each=nsim),dim=c(nsim,maxage))   # Age array
  # surv<-exp(-Marray[,1])^(Agearray-1)                           # Survival array
  # Nfrac<-surv*Mat_age                                           # predicted Numbers of mature ages
  # initdist<-as.matrix(cbind(Frac_area_1,1-Frac_area_1))         # Get the initial spatial distribution of each simulated population
  
  # R0a<-R0*initdist                                              # Unfished recruitment by area
  
  # SAYR<-as.matrix(expand.grid(1:nareas,1,1:maxage,1:nsim)[4:1]) # Set up some array indexes sim (S) age (A) year (Y) region/area (R)
  # SAY<-SAYR[,1:3]
  # SA<-SAYR[,1:2]
  # SR<-SAYR[,c(1,4)]
  # S<-SAYR[,1]
  # SY<-SAYR[,c(1,3)]
  
  # SSN[SAYR]<-Nfrac[SA]*R0*initdist[SR]                           # Calculate initial spawning stock numbers
  # N[SAYR]<-R0*surv[SA]*initdist[SR]                              # Calculate initial stock numbers
  # Biomass[SAYR]<-N[SAYR]*Wt_age[SAY]                             # Calculate initial stock biomass
  # SSB[SAYR]<-SSN[SAYR]*Wt_age[SAY]                               # Calculate spawning stock biomass
  # VBiomass[SAYR]<-Biomass[SAYR]*V[SAY]                            # Calculate vunerable biomass
  # SSN0<-apply(SSN[,,1,],c(1,3),sum)                              # Calculate unfished spawning stock numbers
  # SSB0<-apply(SSB[,,1,],1,sum)                                   # Calculate unfished spawning stock numbers
  # SSBpR<-SSB0/R0                                                 # Spawning stock biomass per recruit
  # SSB0a<-apply(SSB[,,1,],c(1,3),sum)                             # Calculate unfished spawning stock numbers
  # bR<-log(5*hs)/(0.8*SSB0a)                                      # Ricker SR params
  # aR<-exp(bR*SSB0a)/SSBpR                                        # Ricker SR params
  
  # print("Optimizing for user-specified depletion")               # Print a progress update
  # flush.console()                                                # update console
  
  # if(sfIsRunning()){
    # sfExport(list=c("dep","Find","Perr","Marray","hs","Mat_age","Wt_age","R0","V","nyears","maxage","SRrel","aR","bR"))
    # qs<-sfSapply(1:nsim,getq,dep,Find,Perr,Marray,hs,Mat_age,Wt_age,R0,V,nyears,maxage,mov,Spat_targ,SRrel,aR,bR) # find the q that gives current stock depletion
  # }else{
    # qs <- sapply(1:nsim,getq,dep,Find,Perr,Marray,hs,Mat_age,Wt_age,R0,V,nyears,maxage,mov,Spat_targ,SRrel,aR,bR) # find the q that gives current stock depletion
  # }
  
  # # Check that depletion target is reached
  # HighQ <- which(qs> 13)
  # if (length(HighQ) > 0) { # If q has hit bound, re-sample depletion and try again. Tries 10 times 
    # # and then alerts user
    # Err <- TRUE
    # Nsec <- 10
    # Nprob <- length(HighQ)
    # message(Nprob," simulations have final biomass that is not close to specified depletion. \n qs have hit bounds.\n ")
    # message("Re-sampling depletion and trying again")
    # #This is likely going to cause problems later on, such as infinite Fs")
    # #message("Check life history gradients")
    # # if (length(HighQ) == 1) HighQ <- rep(HighQ,2)
    # #PlotLHs(HighQ) # plot life history matrices	
    # # Attempt again with different depletions
    # count <- 0
    # while (Err & count < 10) {
      # count <- count + 1
      # message("Attempt ", count)
      # flush.console()
      # Nprob <- length(HighQ)
      # dep[HighQ] <- runif(Nprob,OM@D[1],OM@D[2])
      # if(sfIsRunning()){
        # sfExport(list=c("dep","Find","Perr","Marray","hs","Mat_age","Wt_age","R0","V","nyears","maxage","SRrel","aR","bR"))
        # qs[HighQ]<-sfSapply(HighQ,getq,dep,Find,Perr,Marray,hs,Mat_age,Wt_age,R0,V,nyears,maxage,mov,Spat_targ,SRrel,aR,bR) # find the q that gives current stock depletion
      # }else{
        # qs[HighQ] <- sapply(HighQ,getq,dep,Find,Perr,Marray,hs,Mat_age,Wt_age,R0,V,nyears,maxage,mov,Spat_targ,SRrel,aR,bR) # find the q that gives current stock depletion
      # }
      # HighQ <- which(qs> 13)
      # if (length(HighQ) == 0) Err <- FALSE
    # }
    # if (!Err) {
      # cat ("Success \n")
      # cat("q range = ", range(qs),"\n")
      # flush.console()
    # }
    # if (Err) { # still a problem
      # print(qs[HighQ])
      # cat ("qs still very high \n")
      # cat ("Press ESC to quit now, or will continue in", Nsec, "seconds \n")
      # flush.console()
      # for (xx in 1:Nsec) {
        # Sys.sleep(1)
        # message(Nsec+1-xx)
        # flush.console()
      # }
    # }  
  # }
  
  # print("Calculating historical stock and fishing dynamics")     # Print a progress update
  # flush.console()                                                # update console
  
  # fishdist<-(apply(VBiomass[,,1,],c(1,3),sum)^Spat_targ)/apply(apply(VBiomass[,,1,],c(1,3),sum)^Spat_targ,1,mean)  # spatial preference according to spatial biomass
  # FM[SAYR]<-qs[S]*Find[SY]*V[SAY]*fishdist[SR]                    # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
  # Z[SAYR]<-FM[SAYR]+Marray[SY]                                   # Total mortality rate                 
  
  # for(y in 1:(nyears-1)){
    # # set up some indices for indexed calculation
    # SAYR<-as.matrix(expand.grid(1:nareas,y,1:maxage,1:nsim)[4:1]) # Set up some array indexes sim (S) age (A) year (Y) region/area (R)
    # SAY1R<-as.matrix(expand.grid(1:nareas,y+1,1:maxage,1:nsim)[4:1])
    # SAY<-SAYR[,1:3]
    # SA<-SAYR[,1:2]
    # SR<-SAYR[,c(1,4)]
    # S<-SAYR[,1]
    # SY<-SAYR[,c(1,3)]
    # SY1<-SAY1R[,c(1,3)]
    # indMov<-as.matrix(expand.grid(1:nareas,1:nareas,y+1,1:maxage,1:nsim)[5:1]) # Movement master index
    # indMov2<-indMov[,c(1,2,3,4)]                                               # Movement from index
    # indMov3<-indMov[,c(1,4,5)]                                                 # Movement to index
    
    # if(SRrel[1]==1){
      # N[,1,y+1,]<-Perr[,y]*(0.8*R0a*hs*apply(SSB[,,y,],c(1,3),sum))/(0.2*SSBpR*R0a*(1-hs)+(hs-0.2)*apply(SSB[,,y,],c(1,3),sum))  # Recruitment assuming regional R0 and stock wide steepness
    # }else{ # most transparent form of the Ricker uses alpha and beta params
      # N[,1,y+1,]<-Perr[,y+nyears]*aR*apply(SSB[,,y,],c(1,3),sum)*exp(-bR*apply(SSB[,,y,],c(1,3),sum))
    # }
    
    # fishdist<-(apply(VBiomass[,,y,],c(1,3),sum)^Spat_targ)/apply(apply(VBiomass[,,y,],c(1,3),sum)^Spat_targ,1,mean)   # spatial preference according to spatial biomass
    # FM[SAY1R]<-qs[S]*Find[SY1]*V[SAY]*fishdist[SR]                           # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
    # Z[SAY1R]<-FM[SAY1R]+Marray[SY]                                          # Total mortality rate
    # N[,2:maxage,y+1,]<-N[,1:(maxage-1),y,]*exp(-Z[,1:(maxage-1),y,])        # Total mortality
    # temp<-array(N[indMov2]*mov[indMov3],dim=c(nareas,nareas,maxage,nsim))   # Move individuals
    # N[,,y+1,]<-apply(temp,c(4,3,1),sum)
    # Biomass[SAY1R]<-N[SAY1R]*Wt_age[SAY]                                    # Calculate biomass
    # VBiomass[SAY1R]<-Biomass[SAY1R]*V[SAY]                                   # Calculate vulnerable biomass
    # SSN[SAY1R]<-N[SAY1R]*Mat_age[SA]                                        # Calculate spawning stock numbers
    # SSB[SAY1R]<-SSN[SAY1R]*Wt_age[SAY]                                      # Calculate spawning stock biomass
    
  # } # end of year
  
  # CN<-apply(N*(1-exp(-Z))*(FM/Z),c(1,3,2),sum)  # Catch in numbers
  # CN[is.na(CN)]<-0
  # CB<-Biomass*(1-exp(-Z))*(FM/Z)                                            # Catch in biomass
  
  # Cbiasa<-array(Cbias,c(nsim,nyears+proyears))                              # Bias array
  # Cerr<-array(rlnorm((nyears+proyears)*nsim,mconv(1,rep(Csd,(nyears+proyears))),sdconv(1,rep(Csd,nyears+proyears))),c(nsim,nyears+proyears)) # composite of bias and observation error
  # Cobs<-Cbiasa[,1:nyears]*Cerr[,1:nyears]*apply(CB,c(1,3),sum)              # Simulated observed catch (biomass)
  
  # CAA<-array(NA,dim=c(nsim,nyears,maxage))                                  # Catch  at age array
  # cond<-apply(CN,1:2,sum,na.rm=T)<1                                         # this is a fix for low sample sizes. If CN is zero across the board a single fish is caught in age class of model selectivity (dumb I know)
  # fixind<-as.matrix(cbind(expand.grid(1:nsim,1:nyears),rep(ceiling(mod),nyears))) # more fix
  # CN[fixind[cond,]]<-1                                                      # puts a catch in the most vulnerable age class
  # for(i in 1:nsim)for(j in 1:nyears)CAA[i,j,]<-ceiling(-0.5+rmultinom(1,CAA_nsamp[i],CN[i,j,])*CAA_nsamp[i]/CAA_ESS[i]) # a multinomial observation model for catch-at-age data
  
  # LatASD <- Len_age * 0.1 # This is currently fixed to cv of 10%
  # MaxBin <- ceiling(max(Linfarray) + 2 * max(LatASD))
  # binWidth <- ceiling(0.05 * MaxBin)
  # CAL_bins <- seq(from=0, to=MaxBin+binWidth, by=binWidth) 
  # CAL_binsmid <- seq(from=0.5*binWidth, by=binWidth, length=length(CAL_bins)-1)
  # nCALbins <- length(CAL_binsmid)
  
  # CAL <- array(NA,dim=c(nsim,nyears,nCALbins))                                # the catch at length array
  # LFC <- rep(NA,nsim) # length at first capture
  
  # for(i in 1:nsim){
    # for(j in 1:nyears){
      # tempCN<-rmultinom(1, size=CAL_ESS[i], prob=CN[i,j,])
      # #ages <- rep(1:maxage,tempCN)+runif(sum(tempCN),-0.5,0.5)          # sample expected age
      # lens <- unlist(sapply(1:maxage, function (X) rnorm(tempCN[X],  Len_age[i,X,j], LatASD[i,X,j])))
      # lens[lens > (max(Linfarray) + 2 * max(LatASD))|lens>max(CAL_bins)] <- max(Linfarray) + 2 * max(LatASD) # truncate at 2 sd 
      # CAL[i,j,] <- hist(lens,CAL_bins,plot=F)$counts                       # assign to bins
      # LFC[i] <- min(c(lens,LFC[i]),na.rm=T)                                # get the smallest CAL observation
      # #CAL[i,j,] <- ceiling(rmultinom(1, size=ESS[i], prob=tempCAL)*CAL_nsamp[i]*CAL_ESS[i]-0.5) # could replace with Dirichlet distribution
    # }
  # }
  
  # Ierr<-array(rlnorm((nyears+proyears)*nsim,mconv(1,rep(Isd,nyears+proyears)),sdconv(1,rep(Isd,nyears+proyears))),c(nsim,nyears+proyears))
  # II<-(apply(Biomass,c(1,3),sum)*Ierr[,1:nyears])^betas     # apply hyperstability / hyperdepletion
  # II<-II/apply(II,1,mean)                                   # normalize
  
  # print("Calculating MSY reference points")                 # Print a progress update
  # flush.console()                                           # update the console
  
  # if(sfIsRunning()){
    # sfExport(list=c("Marray","hs","Mat_age","Wt_age","R0","V","nyears","maxage")) # export some newly made arrays to the cluster

    # MSYrefs<-sfSapply(1:nsim,getFMSY,Marray,hs,Mat_age,Wt_age,R0,V=V[,,nyears],maxage,nyears,proyears=200,Spat_targ,mov,SRrel,aR,bR) # optimize for MSY reference points
	
  # }else{

    # MSYrefs<-sapply(1:nsim,getFMSY,Marray,hs,Mat_age,Wt_age,R0,V=V[,,nyears],maxage,nyears,proyears=200,Spat_targ,mov,SRrel,aR,bR) # optimize for MSY reference points
  # }
    
  # MSY<-MSYrefs[1,]  # record the MSY results (Vulnerable)
  # FMSY<-MSYrefs[2,] # instantaneous FMSY  (Vulnerable)
  # BMSY<-(MSY/(1-exp(-FMSY))) # Biomass at MSY (Vulnerable)
  # BMSY_B0<-MSYrefs[3,] # SSBMSY relative to unfished (SSB)
  # BMSY_B0bias<-array(rlnorm(nsim*ntest,mconv(1,OM@BMSY_B0cv),sdconv(1,OM@BMSY_B0cv)),dim=c(nsim,ntest)) # trial samples of BMSY relative to unfished
  
  # print("Calculating reference yield - best fixed F strategy") # Print a progress update
  # flush.console()                                              # update the console
  
  # if(sfIsRunning()){ # Numerically optimize for F that provides highest long term yield
    # RefY<-sfSapply(1:nsim,getFref,Marray=Marray,Wt_age=Wt_age,Mat_age=Mat_age,Perr=Perr,N_s=N[,,nyears,],SSN_s=SSN[,,nyears,],
                   # Biomass_s=Biomass[,,nyears,],VBiomass_s=VBiomass[,,nyears,],SSB_s=SSB[,,nyears,],
                   # Vn=V[,,nyears],hs=hs,R0a=R0a,nyears=nyears,proyears=proyears,nareas=nareas,maxage=maxage,mov=mov,SSBpR=SSBpR,
                   # aR=aR,bR=bR,SRrel=SRrel)
  # }else{
    # RefY<-sapply(1:nsim,getFref,Marray=Marray,Wt_age=Wt_age,Mat_age=Mat_age,Perr=Perr,N_s=N[,,nyears,],SSN_s=SSN[,,nyears,],
                 # Biomass_s=Biomass[,,nyears,],VBiomass_s=VBiomass[,,nyears,],SSB_s=SSB[,,nyears,],
                 # Vn=V[,,nyears],hs=hs,R0a=R0a,nyears=nyears,proyears=proyears,nareas=nareas,maxage=maxage,mov=mov,SSBpR=SSBpR,
                 # aR=aR,bR=bR,SRrel=SRrel) 
  # }

 
  # Depletion<-(apply(Biomass[,,nyears,],1,sum)/apply(Biomass[,,1,],1,sum))#^betas   # apply hyperstability / hyperdepletion
  # #cbind(dep,Depletion)
  # FMSY_M<-FMSY/M                      # ratio of true FMSY to natural mortality rate M
  # # LFS<-Linf*(1-exp(-K*(mod-t0)))      # Length at full selection
  # A<-apply(VBiomass[,,nyears,],1,sum) # Abundance
  # OFLreal<-A*FMSY                     # the true simulated Over Fishing Limit
  # Recerr<-array(rlnorm((nyears+proyears)*nsim,mconv(1,rep(Recsd,(nyears+proyears))),sdconv(1,rep(Recsd,nyears+proyears))),c(nsim,nyears+proyears))
  
  # test<-array(BMSY_B0*BMSY_B0bias,dim=c(nsim,ntest)) # the simulated observed BMSY_B0 
  # indy<-array(rep(1:ntest,each=nsim),c(nsim,ntest))  # index
  # indy[test>0.9]<-NA                                 # interval censor
  # BMSY_B0bias<-BMSY_B0bias[cbind(1:nsim,apply(indy,1,min,na.rm=T))] # sample such that BMSY_B0<90%
  
  # I3<-apply(Biomass,c(1,3),sum)^betas     # apply hyperstability / hyperdepletion
  # I3<-I3/apply(I3,1,mean)                 # normalize index to mean 1
  # Iref<-apply(I3[,1:5],1,mean)*BMSY_B0    # return the real target abundance index corresponding to BMSY
  
  # hsim<-rep(NA,nsim)                      # simulate values in steepness 
  # cond<-hs>0.6
  # hsim[cond]<-0.2+rbeta(sum(hs>0.6),alphaconv((hs[cond]-0.2)/0.8,(1-(hs[cond]-0.2)/0.8)*OM@hcv),betaconv((hs[cond]-0.2)/0.8,(1-(hs[cond]-0.2)/0.8)*OM@hcv))*0.8
  # hsim[!cond]<-0.2+rbeta(sum(hs<0.6),alphaconv((hs[!cond]-0.2)/0.8,(hs[!cond]-0.2)/0.8*OM@hcv),betaconv((hs[!cond]-0.2)/0.8,(hs[!cond]-0.2)/0.8*OM@hcv))*0.8
  # hbias<-hsim/hs                          # back calculate the simulated bias
  
  # DLM_data<-new('DLM_data',stock="MSE")             # create a blank DLM data object
  # if(reps==1)DLM_data<-OneRep(DLM_data)             # make stochastic variables certain for only one rep
  # DLM_data<-replic8(DLM_data,nsim)                  # make nsim sized slots in the DLM data object
  # DLM_data@Name<-OM@Name
  # DLM_data@Year<-1:nyears
  # DLM_data@Cat<-Cobs
  # DLM_data@Ind<-II
  # DLM_data@Rec<-apply(N[,1,,],c(1,2),sum)*Recerr[,1:nyears] 
  # DLM_data@t<-rep(nyears,nsim)
  # DLM_data@AvC<-apply(Cobs,1,mean)
  # DLM_data@Dt<-Dbias*Depletion*rlnorm(nsim,mconv(1,Derr),sdconv(1,Derr))
  # DLM_data@Mort<-M*Mbias
  # DLM_data@FMSY_M<-FMSY_M*FMSY_Mbias
  # DLM_data@BMSY_B0<-BMSY_B0*BMSY_B0bias
  # DLM_data@Cref<-MSY*Crefbias
  # DLM_data@Bref<-BMSY*Brefbias
  # DLM_data@Iref<-Iref*Irefbias
  # DLM_data@LFC<-LFC*LFCbias
  # DLM_data@LFS<-LFS*LFSbias
  # DLM_data@CAA<-CAA
  # DLM_data@Dep<-Dbias*Depletion*rlnorm(nsim,mconv(1,Derr),sdconv(1,Derr))
  # DLM_data@Abun<-A*Abias*rlnorm(nsim,mconv(1,Aerr),sdconv(1,Aerr))
  # DLM_data@vbK<-K*Kbias
  # DLM_data@vbt0<-t0*t0bias
  # DLM_data@vbLinf<-Linf*Linfbias
  # DLM_data@L50 <- lenM * lenMbias
  # DLM_data@L95 <- len95 * lenMbias
  # DLM_data@L95[DLM_data@L95>0.9*DLM_data@vbLinf]<-0.9*DLM_data@vbLinf[DLM_data@L95>0.9*DLM_data@vbLinf] # Set a hard limit on ratio of L95 to Linf
  # DLM_data@L50[DLM_data@L50>0.9*DLM_data@L95]<-0.9*DLM_data@L95[DLM_data@L50>0.9*DLM_data@L95] # Set a hard limit on ratio of L95 to Linf
  # DLM_data@steep<-hs*hbias
  # DLM_data@CAL_bins<-CAL_bins
  # DLM_data@CAL<-CAL
  # MLbin<-(CAL_bins[1:(length(CAL_bins)-1)]+CAL_bins[2:length(CAL_bins)])/2
  # temp<-CAL*rep(MLbin,each=nsim*nyears)
  # DLM_data@ML<-apply(temp,1:2,sum)/apply(CAL,1:2,sum) 
  # DLM_data@Lc<-array(MLbin[apply(CAL,1:2,which.max)],dim=c(nsim,nyears))
  # nuCAL<-CAL
  # for(i in 1:nsim)for(j in 1:nyears)nuCAL[i,j,1:match(max(1,DLM_data@Lc[i,j]),MLbin)]<0
  # temp<-nuCAL*rep(MLbin,each=nsim*nyears)
  # DLM_data@Lbar<-apply(temp,1:2,sum)/apply(nuCAL,1:2,sum)
  # DLM_data@MaxAge<-maxage
  # DLM_data@Units<-"unitless"
  # DLM_data@Ref<-OFLreal
  # DLM_data@Ref_type<-'Simulated OFL'
  # DLM_data@wla<-rep(OM@a,nsim)
  # DLM_data@wlb<-rep(OM@b,nsim)
  # DLM_data@OM<-as.data.frame(cbind(RefY,M,Depletion,A,BMSY_B0,FMSY_M,Mgrad,Msd,procsd,Esd,dFfinal,MSY,qinc,qcv,
                                   # FMSY,Linf,K,t0,hs,Linfgrad,Kgrad,Linfsd,recgrad,Ksd,ageM,
                                   # L5,LFS,Vmaxlen,LFC,OFLreal,Spat_targ,Frac_area_1,Prob_staying,AC)) # put all the operating model parameters in one table
  
  # DLM_data@Obs<-as.data.frame(cbind(Cbias,Csd,CAA_nsamp,CAA_ESS,CAL_nsamp,CAL_ESS,Isd,Dbias,Derr,Mbias,FMSY_Mbias,BMSY_B0bias,
                                    # lenMbias,LFCbias,LFSbias,Abias,Aerr,Kbias,t0bias,Linfbias,hbias,Irefbias,Crefbias,Brefbias,betas))  # put all the observation error model parameters in one table
  
  # DLM_data@LHYear<-OM@nyears # Last historical year is nyears (for fixed MPs)
  # DLM_data@MPrec<-Cobs[,nyears]
  # #assign("DLM_data",DLM_data,envir=.GlobalEnv) # for debugging fun
  
  # # Run projections ===========================================================================
  # qmu<--0.5*qcv^2                                      # Mean
  # qvar<-array(exp(rnorm(proyears*nsim,rep(qmu,proyears),rep(qcv,proyears))),c(nsim,proyears)) # Variations in interannual variation
  # FinF<-Find[,nyears]
  
  # print("Determining available methods")  # print an progress report
  # flush.console()                         # update the console
  
  # PosMPs <- Can(DLM_data, timelimit=timelimit)                 # list all the methods that could be applied to a DLM data object
  # # print(PosMPs)
  # if(is.na(MPs[1])) {
    # MPs<-PosMPs      # if the user does not supply an argument MPs run the MSE or all available methods
    # print("No MPs specified: running all available")
    # #print(MPs)
  # }	
  # if(!is.na(MPs[1]))MPs<-MPs[MPs%in%PosMPs] # otherwise run the MSE for all methods that are deemed possible
  # if(length(MPs)==0)stop('MSE stopped: no viable methods \n\n') # if none of the user specied methods are possible stop the run
  
  # nMP <- length(MPs)                    # the total number of methods used
  
  # MSElist<-list(DLM_data)[rep(1,nMP)]        # create a data object for each method (they have identical historical data and branch in projected years)
  
  # B_BMSYa<-array(NA,dim=c(nsim,nMP,proyears))  # store the projected B_BMSY
  # F_FMSYa<-array(NA,dim=c(nsim,nMP,proyears))  # store the projected F_FMSY
  # Ba<-array(NA,dim=c(nsim,nMP,proyears))       # store the projected Biomass
  # FMa<-array(NA,dim=c(nsim,nMP,proyears))      # store the projected fishing mortality rate
  # Ca<-array(NA,dim=c(nsim,nMP,proyears))       # store the projected catch
  # TACa<-array(NA,dim=c(nsim,nMP,proyears))     # store the projected TAC recommendation
  
  # for(mm in 1:nMP){    # MSE Loop over methods
    
    # print(paste(mm,"/",nMP," Running MSE for ",MPs[mm],sep=""))  # print a progress report
    # flush.console()                                                  # update the console
    
    # # projection arrays
    # N_P<-array(NA,dim=c(nsim,maxage,proyears,nareas))
    # Biomass_P<-array(NA,dim=c(nsim,maxage,proyears,nareas))
    # VBiomass_P<-array(NA,dim=c(nsim,maxage,proyears,nareas))
    # SSN_P<-array(NA,dim=c(nsim,maxage,proyears,nareas))
    # SSB_P<-array(NA,dim=c(nsim,maxage,proyears,nareas))
    # FM_P<-array(NA,dim=c(nsim,maxage,proyears,nareas))
    # FM_nospace<-array(NA,dim=c(nsim,maxage,proyears,nareas)) # stores prospective F before reallocation to new areas
    # FML<-array(NA,dim=c(nsim,nareas)) # last apical F
    # Z_P<-array(NA,dim=c(nsim,maxage,proyears,nareas))
    # CB_P<-array(NA,dim=c(nsim,maxage,proyears,nareas))
    
    # # indexes
    # SAYRL<-as.matrix(expand.grid(1:nsim,1:maxage,nyears,1:nareas))   # Final historical year
    # SAYRt<-as.matrix(expand.grid(1:nsim,1:maxage,1+nyears,1:nareas)) # Trajectory year
    # SAYR<-as.matrix(expand.grid(1:nsim,1:maxage,1,1:nareas))
    # SYt<-SAYRt[,c(1,3)]
    # SAYt<-SAYRt[,1:3]
    # SR<-SAYR[,c(1,4)]
    # SA1<-SAYR[,1:2]
    # S1<-SAYR[,1]
    # SY1<-SAYR[,c(1,3)]
    # SAY1<-SAYR[,1:3]
    
    # SYA<-as.matrix(expand.grid(1:nsim,1,1:maxage))         # Projection year
    # SY<-SYA[,1:2]
    # SA<-SYA[,c(1,3)]
    # SAY<-SYA[,c(1,3,2)]
    # S<-SYA[,1]
    
    # if(SRrel[1]==1){
      # N_P[,1,1,]<-Perr[,nyears]*(0.8*R0a*hs*apply(SSB[,,nyears,],c(1,3),sum))/(0.2*SSBpR*R0a*(1-hs)+(hs-0.2)*apply(SSB[,,nyears,],c(1,3),sum))  # Recruitment assuming regional R0 and stock wide steepness
    # }else{ # most transparent form of the Ricker uses alpha and beta params
      # N_P[,1,1,]<-Perr[,nyears]*aR*apply(SSB[,,nyears,],c(1,3),sum)*exp(-bR*apply(SSB[,,nyears,],c(1,3),sum))
    # }
    # indMov<-as.matrix(expand.grid(1:nareas,1:nareas,1,1:maxage,1:nsim)[5:1])
    # indMov2<-indMov[,c(1,2,3,4)]
    # indMov3<-indMov[,c(1,4,5)]
    
    # N_P[,2:maxage,1,]<-N[,1:(maxage-1),nyears,]*exp(-Z[,1:(maxage-1),nyears,])        # Total mortality
    # temp<-array(N_P[indMov2]*mov[indMov3],dim=c(nareas,nareas,maxage,nsim))   # Move individuals
    # N_P[,,1,]<-apply(temp,c(4,3,1),sum)
    # Biomass_P[SAYR]<-N_P[SAYR]*Wt_age[SAY1]                                    # Calculate biomass
    # VBiomass_P[SAYR]<-Biomass_P[SAYR]*V[SAYt]                                  # Calculate vulnerable biomass
    # SSN_P[SAYR]<-N_P[SAYR]*Mat_age[SA1]                                        # Calculate spawning stock numbers
    # SSB_P[SAYR]<-SSN_P[SAYR]*Wt_age[SAY1]
    # FML<-apply(FM[,,nyears,],c(1,3),max)
    
      
    # if(class(match.fun(MPs[mm]))=="DLM_output"){
      
      # TACused <- apply(Sam(MSElist[[mm]],MPs=MPs[mm],perc=pstar,reps=reps)@TAC,3,quantile,p=pstar,na.rm=T)
      # TACa[,mm,1]<-TACused
      # fishdist<-(apply(VBiomass_P[,,1,],c(1,3),sum)^Spat_targ)/apply(apply(VBiomass_P[,,1,],c(1,3),sum)^Spat_targ,1,mean)   # spatial preference according to spatial biomass
      # CB_P[SAYR]<-Biomass_P[SAYR]*(1-exp(-V[SAYt]*fishdist[SR]))      # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space
      # temp<-CB_P[,,1,]/apply(CB_P[,,1,],1,sum)   # how catches are going to be distributed
      # CB_P[,,1,]<-TACused*temp           # debug - to test distribution code make TAC = TAC2, should be identical
      # temp<-CB_P[SAYR]/(Biomass_P[SAYR]*exp(-Marray[SYt]/2)) # Pope's approximation
      # temp[temp>(1-exp(-maxF))]<-1-exp(-maxF)
      # FM_P[SAYR]<--log(1-temp)
      
    # }else{ # input control
      
      # inc<-sapply(1:nsim,MPs[mm],DLM_data=MSElist[[mm]])
      # Ai<-inc[1,]
      # Ei<-inc[2,]
      # Si<-t(inc[3:4,])
      # Vi<-t(inc[5:(4+maxage),])
      # y<-1
      
      # if(sum(Si!=1)==0){ # if there is no spatial closure
        # if(sum(!is.na(Vi[1,1]))==0){ # if no vulnerability schedule is specified
          # newVB<-apply(VBiomass_P[,,y,],c(1,3),sum) # vulnerability isn't changed
          # fishdist<-(newVB^Spat_targ)/apply(newVB^Spat_targ,1,mean)   # spatial preference according to spatial biomass
          # FM_P[SAYR]<-FinF[S1]*Ei[S1]*V[SAYt]*fishdist[SR]*qvar[SY1]*qs[S1]*(1+qinc[S1]/100)^y   # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
        # }else{
          # newVB<-apply(VBiomass_P[,,y,]*Vi[SA1],c(1,3),sum) # vulnerability modified
          # fishdist<-(newVB^Spat_targ)/apply(newVB^Spat_targ,1,mean)   # spatial preference according to spatial biomass
          # FM_P[SAYR]<-FinF[S1]*Ei[S1]*Vi[SA1]*fishdist[SR]*qvar[SY1]*qs[S1]*(1+qinc[S1]/100)^y   # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
        # }
      # }else{  # A spatial closure
        # if(sum(!is.na(Vi[1,1]))==0){ # if no vulnerability schedule is specified
          # newVB<-apply(VBiomass_P[,,y,],c(1,3),sum) # vulnerability isn't changed
          # fishdist<-(newVB^Spat_targ)/apply(newVB^Spat_targ,1,mean)   # spatial preference according to spatial biomass
          # Emult<-1+((2/apply(fishdist*Si,1,sum))-1)*Ai  # allocate effort to new area according to fraction allocation Ai
          # FM_P[SAYR]<-FinF[S1]*Ei[S1]*V[SAYt]*Si[SR]*fishdist[SR]*Emult[S1]*qvar[SY1]*qs[S1]^(1+qinc[S1]/100)^y 
        # }else{
          # newVB<-apply(VBiomass_P[,,y,]*Vi[SA1],c(1,3),sum) # vulnerability modified
          # fishdist<-(newVB^Spat_targ)/apply(newVB^Spat_targ,1,mean)   # spatial preference according to spatial biomass
          # Emult<-1+((2/apply(fishdist*Si,1,sum))-1)*Ai  # allocate effort to new area according to fraction allocation Ai
          # FM_P[SAYR]<-FinF[S1]*Ei[S1]*Vi[SA1]*Si[SR]*fishdist[SR]*Emult[S1]*qvar[SY1]*qs[S1]^(1+qinc[S1]/100)^y 
        # } # vulnerability specified
      # }  # spatial closure specified  
    # }   # input control  
    
    # CB_P[SAYR]<-Biomass_P[SAYR]*(1-exp(-FM_P[SAYR])) 
    # Z_P[SAYR]<-FM_P[SAYR]+Marray[SYt] 
   
    # upyrs<-1+(0:(floor(proyears/interval)-1))*interval  # the years in which there are updates (every three years)
    # cat(".")
    # flush.console()
    
    # for(y in 2:proyears){
      
      # cat(".")
      # flush.console()
      # if(class(match.fun(MPs[mm]))=="DLM_output") TACa[,mm,y]<-TACused
      # SAYRt<-as.matrix(expand.grid(1:nsim,1:maxage,y+nyears,1:nareas)) # Trajectory year
      # SAYt<-SAYRt[,1:3]
      # SYt<-SAYRt[,c(1,3)]
      # SAY1R<-as.matrix(expand.grid(1:nsim,1:maxage,y-1,1:nareas))
      # SAYR<-as.matrix(expand.grid(1:nsim,1:maxage,y,1:nareas))
      # SY<-SAYR[,c(1,3)]
      # SA<-SAYR[,1:2]
      # S1<-SAYR[,1]
      
      # SAY<-SAYR[,1:3]
      # S<-SAYR[,1]
      # SR<-SAYR[,c(1,4)]
      # SA2YR<-as.matrix(expand.grid(1:nsim,2:maxage,y,1:nareas))
      # SA1YR<-as.matrix(expand.grid(1:nsim,1:(maxage-1),y-1,1:nareas))
      # indMov<-as.matrix(expand.grid(1:nareas,1:nareas,y,1:maxage,1:nsim)[5:1])
      # indMov2<-indMov[,c(1,2,3,4)]
      # indMov3<-indMov[,c(1,4,5)]
      
      # N_P[SA2YR]<-N_P[SA1YR]*exp(-Z_P[SA1YR])         # Total mortality
      # if(SRrel[1]==1){
        # N_P[,1,y,]<-Perr[,y+nyears]*(0.8*R0a*hs*apply(SSB_P[,,y-1,],c(1,3),sum))/(0.2*SSBpR*R0a*(1-hs)+(hs-0.2)*apply(SSB_P[,,y-1,],c(1,3),sum))  # Recruitment assuming regional R0 and stock wide steepness
      # }else{ # most transparent form of the Ricker uses alpha and beta params
        # N_P[,1,y,]<-Perr[,y+nyears]*aR*apply(SSB_P[,,y-1,],c(1,3),sum)*exp(-bR*apply(SSB_P[,,y-1,],c(1,3),sum))
      # }
      
      # temp<-array(N_P[indMov2]*mov[indMov3],dim=c(nareas,nareas,maxage,nsim))  # Move individuals
      # N_P[,,y,]<-apply(temp,c(4,3,1),sum)
      
      # Biomass_P[SAYR]<-N_P[SAYR]*Wt_age[SAYt]                                    # Calculate biomass
      # VBiomass_P[SAYR]<-Biomass_P[SAYR]*V[SAYt]                       # Calculate vulnerable biomass
      # SSN_P[SAYR]<-N_P[SAYR]*Mat_age[SA]                                       # Calculate spawning stock numbers
      # SSB_P[SAYR]<-SSN_P[SAYR]*Wt_age[SAYt]                                      # Calculate spawning stock biomass
      
      # if(y%in%upyrs){  # rewrite the DLM object and run the TAC function
        # yind<-upyrs[match(y,upyrs)-1]:(upyrs[match(y,upyrs)]-1)
        # CNtemp<-array(N_P[,,yind,]*exp(Z_P[,,yind,])*(1-exp(-Z_P[,,yind,]))*(FM_P[,,yind,]/Z_P[,,yind,]),c(nsim,maxage,interval,nareas))
        # CBtemp<-array(Biomass_P[,,yind,]*exp(Z_P[,,yind,])*(1-exp(-Z_P[,,yind,]))*(FM_P[,,yind,]/Z_P[,,yind,]),c(nsim,maxage,interval,nareas))
        # CNtemp[is.na(CNtemp)]<-1e-20
        # CBtemp[is.na(CNtemp)]<-1e-20
        # CNtemp<-apply(CNtemp,c(1,3,2),sum)
        # CNtemp[CNtemp==0]<-CNtemp[CNtemp==0]+tiny
        # Cobs<-Cbiasa[,nyears+yind]*Cerr[,nyears+yind]*apply(CBtemp,c(1,3),sum)
        # Cobs[is.na(Cobs)]<-tiny
        # Recobs<-Recerr[,nyears+yind]*apply(array(N_P[,1,yind,],c(nsim,interval,nareas)),c(1,2),sum)
        
        # CAA<-array(NA,dim=c(nsim,interval,maxage))                                  # Catch  at age array
        # for(i in 1:nsim)for(j in 1:interval)CAA[i,j,]<-ceiling(-0.5+rmultinom(1,CAA_nsamp[i],CNtemp[i,j,])*CAA_nsamp[i]/CAA_ESS[i]) # a multinomial observation model for catch-at-age data
        
        # CAL <- array(NA,dim=c(nsim,interval,nCALbins))                                # the catch at length array
        # CNtemp[is.na(CNtemp)]<-0
        # for(i in 1:nsim){
          # for(j in 1:interval){
            # yy<-yind[j]
            # tempCN<-rmultinom(1, size=CAL_ESS[i], prob=CNtemp[i,j,])
            # #ages <- rep(1:maxage,tempCN)+runif(sum(tempCN),-0.5,0.5)          # sample expected age
            # lens <- unlist(sapply(1:maxage, function (X) rnorm(tempCN[X],  Len_age[i,X,yy+nyears], LatASD[i,X,yy+nyears])))
            # lens[lens > (max(Linfarray) + 2 * max(LatASD))|lens>max(CAL_bins)] <- max(Linfarray) + 2 * max(LatASD) # truncate at 2 sd 
            # CAL[i,j,] <- hist(lens,CAL_bins,plot=F)$counts                       # assign to bins
          # }
        # }
        
        # I2<-cbind(apply(Biomass,c(1,3),sum),apply(Biomass_P,c(1,3),sum)[,1:(y-1)])*Ierr[,1:(nyears+(y-1))]^betas
        # I2[is.na(I2)]<-tiny
        # I2<-I2/apply(I2,1,mean)
        
        # Depletion<-apply(Biomass_P[,,y,],1,sum)/apply(Biomass[,,1,],1,sum)
        # A<-apply(VBiomass_P[,,y,],1,sum)
        # A[is.na(A)]<-tiny
        # OFLreal<-A*FMSY
        
        # # assign all the new data
        # MSElist[[mm]]@OM$A<-A
        # MSElist[[mm]]@Year<-1:(nyears+y-1)
        # MSElist[[mm]]@Cat<-cbind(MSElist[[mm]]@Cat,Cobs)
        # MSElist[[mm]]@Ind<-I2
        # MSElist[[mm]]@Rec<-cbind(MSElist[[mm]]@Rec,Recobs)
        # MSElist[[mm]]@t<-rep(nyears+y,nsim)
        # MSElist[[mm]]@AvC<-apply(MSElist[[mm]]@Cat,1,mean)
        # MSElist[[mm]]@Dt<-Dbias*Depletion*rlnorm(nsim,mconv(1,Derr),sdconv(1,Derr))
        # oldCAA<-MSElist[[mm]]@CAA
        # MSElist[[mm]]@CAA<-array(0,dim=c(nsim,nyears+y-1,maxage))
        # MSElist[[mm]]@CAA[,1:(nyears+y-interval-1),]<-oldCAA
        # MSElist[[mm]]@CAA[,nyears+yind,]<-CAA
        # MSElist[[mm]]@Dep<-Dbias*Depletion*rlnorm(nsim,mconv(1,Derr),sdconv(1,Derr))
        # MSElist[[mm]]@Abun<-A*Abias*rlnorm(nsim,mconv(1,Aerr),sdconv(1,Aerr))
        # MSElist[[mm]]@CAL_bins<-CAL_bins
        # oldCAL<-MSElist[[mm]]@CAL
        # MSElist[[mm]]@CAL<-array(0,dim=c(nsim,nyears+y-1,nCALbins))
        # MSElist[[mm]]@CAL[,1:(nyears+y-interval-1),]<-oldCAL
        # MSElist[[mm]]@CAL[,nyears+yind,]<-CAL[,1:interval,]
        
        # temp<-CAL*rep(MLbin,each=nsim*interval)
        # MSElist[[mm]]@ML<-cbind(MSElist[[mm]]@ML,apply(temp,1:2,sum)/apply(CAL,1:2,sum)) 
        # MSElist[[mm]]@Lc<-cbind(MSElist[[mm]]@Lc,array(MLbin[apply(CAL,1:2,which.max)],dim=c(nsim,interval)))
        # nuCAL<-CAL
        # for(i in 1:nsim)for(j in 1:interval)nuCAL[i,j,1:match(max(1,MSElist[[mm]]@Lc[i,j]),MLbin)]<0
        # temp<-nuCAL*rep(MLbin,each=nsim*interval)
        # MSElist[[mm]]@Lbar<-cbind(MSElist[[mm]]@Lbar,apply(temp,1:2,sum)/apply(nuCAL,1:2,sum))
        
        # MSElist[[mm]]@Ref<-OFLreal
        # MSElist[[mm]]@Ref_type<-'Simulated OFL'
        
        # #assign("DLM_data",MSElist[[mm]],envir=.GlobalEnv) # for debugging fun
        
        # if(class(match.fun(MPs[mm]))=="DLM_output"){
          # TACused<-apply(Sam(MSElist[[mm]],MPs=MPs[mm],perc=pstar,reps=reps)@TAC,3,quantile,p=pstar,na.rm=TRUE)
          # NAs <- which(is.na(TACused))
          # if (length(NAs) >0 ) { # robustifying TAC setting!
            # TACused[NAs] <- TACa[NAs,mm,y-1] #
            # if (!exists("store")) store <- list()
            # store <- append(store, c(MPs[mm], NAs))
          # }
          # TACa[,mm,y]<-TACused
          # MSElist[[mm]]@MPrec<-TACused
          # fishdist<-(apply(VBiomass_P[,,y,],c(1,3),sum)^Spat_targ)/apply(apply(VBiomass_P[,,y,],c(1,3),sum)^Spat_targ,1,mean)   # spatial preference according to spatial biomass
          # CB_P[SAYR]<-Biomass_P[SAYR]*(1-exp(-V[SAYt]*fishdist[SR]))      # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space
          # temp<-CB_P[,,y,]/apply(CB_P[,,y,],1,sum)   # how catches are going to be distributed
          # CB_P[,,y,]<-TACused*temp           # debug - to test distribution code make TAC = TAC2, should be identical
          # temp<-CB_P[SAYR]/(Biomass_P[SAYR]*exp(-Marray[SYt]/2)) # Pope's approximation
          # temp[temp>(1-exp(-maxF))]<-1-exp(-maxF)
          # FM_P[SAYR]<--log(1-temp)
          
        # }else{
          
          # inc<-sapply(1:nsim,MPs[mm],DLM_data=MSElist[[mm]])
          # Ai<-inc[1,]
          # Ei<-inc[2,]
          # Si<-t(inc[3:4,])
          # Vi<-t(inc[5:(4+maxage),])
          # MSElist[[mm]]@MPrec<-Ei
          
          # if(sum(Si!=1)==0){ # if there is no spatial closure
            # if(sum(!is.na(Vi[1,1]))==0){ # if no vulnerability schedule is specified
              # newVB<-apply(VBiomass_P[,,y,],c(1,3),sum) # vulnerability isn't changed
              # fishdist<-(newVB^Spat_targ)/apply(newVB^Spat_targ,1,mean)   # spatial preference according to spatial biomass
              # FM_P[SAYR]<-FinF[S1]*Ei[S1]*V[SAYt]*fishdist[SR]*qvar[SY]*qs[S1]*(1+qinc[S1]/100)^y   # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
            # }else{
              # newVB<-apply(VBiomass_P[,,y,]*Vi[SA],c(1,3),sum) # vulnerability modified
              # fishdist<-(newVB^Spat_targ)/apply(newVB^Spat_targ,1,mean)   # spatial preference according to spatial biomass
              # FM_P[SAYR]<-FinF[S1]*Ei[S1]*Vi[SY]*fishdist[SR]*qvar[SY]*qs[S1]*(1+qinc[S1]/100)^y   # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
            # }
          # }else{  # A spatial closure
            # if(sum(!is.na(Vi[1,1]))==0){ # if no vulnerability schedule is specified
              # newVB<-apply(VBiomass_P[,,y,],c(1,3),sum) # vulnerability isn't changed
              # fishdist<-(newVB^Spat_targ)/apply(newVB^Spat_targ,1,mean)   # spatial preference according to spatial biomass
              # Emult<-1+((2/apply(fishdist*Si,1,sum))-1)*Ai  # allocate effort to new area according to fraction allocation Ai
              # FM_P[SAYR]<-FinF[S1]*Ei[S1]*V[SAYt]*Si[SR]*fishdist[SR]*Emult[S1]*qvar[SY]*qs[S1]*(1+qinc[S1]/100)^y 
            # }else{
              # newVB<-apply(VBiomass_P[,,y,]*Vi[SA],c(1,3),sum) # vulnerability modified
              # fishdist<-(newVB^Spat_targ)/apply(newVB^Spat_targ,1,mean)   # spatial preference according to spatial biomass
              # Emult<-1+((2/apply(fishdist*Si,1,sum))-1)*Ai  # allocate effort to new area according to fraction allocation Ai
              # FM_P[SAYR]<-FinF[S1]*Ei[S1]*V[SAYt]*Si[SR]*fishdist[SR]*Emult[S1]*qvar[SY]*qs[S1]*(1+qinc[S1]/100)^y 
            # } #vuln not changed
          # }   # spatial closure
          
        # }   # input or output control 
        # CB_P[SAYR]<-Biomass_P[SAYR]*(1-exp(-FM_P[SAYR])) 
        # Z_P[SAYR]<-FM_P[SAYR]+Marray[SYt] 
        
      # }else{ # not an update yr
        
        # if(class(match.fun(MPs[mm]))=="DLM_output"){
          # CB_P[SAYR]<-Biomass_P[SAYR]*(1-exp(-fishdist[SR]*V[SAYt]))      # ignore magnitude of effort or q increase (just get distribution across age and fishdist across space
          # temp<-CB_P[,,y,]/apply(CB_P[,,y,],1,sum)   # how catches are going to be distributed
          # CB_P[,,y,]<-TACused*temp           # debug - to test distribution code make TAC = TAC2, should be identical
          # temp<-CB_P[SAYR]/(Biomass_P[SAYR]*exp(-Marray[SYt]/2)) # Pope's approximation
          # temp[temp>(1-exp(-maxF))]<-1-exp(-maxF)
          # FM_P[SAYR]<--log(1-temp)
        # }else{ #input contol
          # FM_P[SAYR]<-FM_P[SAY1R]*qvar[SY]*(1+qinc[S1]/100)^y  # add fishing efficiency changes and variability
        # }
        # CB_P[SAYR]<-Biomass_P[SAYR]*(1-exp(-FM_P[SAYR])) 
        # Z_P[SAYR]<-FM_P[SAYR]+Marray[SYt]
      # } # not an update year
    # } # end of year
    
    
    
    # B_BMSYa[,mm,]<-apply(Biomass_P,c(1,3),sum)/BMSY
    # F_FMSYa[,mm,]<-(-log(1-apply(CB_P,c(1,3),sum)/(apply(CB_P,c(1,3),sum)+apply(VBiomass_P,c(1,3),sum))))/FMSY
    # Ba[,mm,]<-apply(Biomass_P,c(1,3),sum)
    # FMa[,mm,]<--log(1-apply(CB_P,c(1,3),sum)/(apply(CB_P,c(1,3),sum)+apply(VBiomass_P,c(1,3),sum)))
    # Ca[,mm,]<-apply(CB_P,c(1,3),sum)
    # cat("\n")
  # }    # end of mm methods
  
  # new('MSE',Name=OM@Name,nyears,proyears,nMP,MPs,nsim,OMtable=DLM_data@OM,DLM_data@Obs,B_BMSYa,F_FMSYa,Ba,FMa,Ca,TACa,SSB_hist=SSB,CB_hist=CB,FM_hist=FM)
  
# }

# Tplot<-function(MSEobj,nam=NA){
  # FMSYr<-quantile(MSEobj@F_FMSY,c(0.001,0.90),na.rm=T)
  # BMSYr<-quantile(MSEobj@B_BMSY,c(0.001,0.975),na.rm=T)

  # colsse<-rainbow(100,start=0,end=0.36)[1:100]
  # colB<-rep(colsse[100],ceiling(BMSYr[2]*100))
  # colB[1:100]<-colsse
  # colB<-makeTransparent(colB,60)
  # colsse<-rainbow(200,start=0,end=0.36)[200:1]
  # colF<-rep(colsse[200],ceiling(FMSYr[2]*100))
  # colF[1:200]<-colsse
  # colF<-makeTransparent(colF,60)

  # Yd<-rep(NA,MSEobj@nMPs)
  # P10<-rep(NA,MSEobj@nMPs)
  # P50<-rep(NA,MSEobj@nMPs)
  # P100<-rep(NA,MSEobj@nMPs)
  # POF<-rep(NA,MSEobj@nMPs)
  # yind<-max(MSEobj@proyears-4,1):MSEobj@proyears
  # RefYd<-MSEobj@OM$RefY

  # for(mm in 1:MSEobj@nMPs){
    # Yd[mm]<-round(mean(apply(MSEobj@C[,mm,yind],1,mean,na.rm=T)/RefYd,na.rm=T)*100,1)
  # #cbind(MSEobj@C[,mm,yind],unlist(MSEobj@OM$MSY))
    # POF[mm]<-round(sum(MSEobj@F_FMSY[,mm,]>1,na.rm=T)/prod(dim(MSEobj@F_FMSY[,mm,]),na.rm=T)*100,1)
    # P10[mm]<-round(sum(MSEobj@B_BMSY[,mm,]<0.1,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
    # P50[mm]<-round(sum(MSEobj@B_BMSY[,mm,]<0.5,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
    # P100[mm]<-round(sum(MSEobj@B_BMSY[,mm,]<1,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
  # }
  
  # #dev.new2(width=7,height=7)
  # par(mfrow=c(2,2),mai=c(0.9,1,0.1,0.1),omi=c(0.1,0.1,0.4,0))

  # tradeoffplot(POF,Yd,"Prob. of overfishing (%)", "Relative yield",MSEobj@MPs[1:MSEobj@nMPs],vl=50,hl=100)
  # tradeoffplot(P100,Yd,"Prob. biomass < BMSY (%)", "Relative yield",MSEobj@MPs[1:MSEobj@nMPs],vl=50,hl=100)
  # tradeoffplot(P50,Yd,"Prob. biomass < 0.5BMSY (%)", "Relative yield",MSEobj@MPs[1:MSEobj@nMPs],vl=50,hl=100)
  # tradeoffplot(P10,Yd,"Prob. biomass < 0.1BMSY (%)", "Relative yield",MSEobj@MPs[1:MSEobj@nMPs],vl=50,hl=100)

  # if(is.na(nam))mtext(deparse(substitute(MSEobj)),3,outer=T,line=0.3,font=2)
  # if(!is.na(nam) & !is.character(nam))mtext(MSEobj@Name,3,outer=T,line=0.3,font=2)
  # if(!is.na(nam) & is.character(nam))mtext(nam,3,outer=T,line=0.3,font=2)
# }

# Tplot2<-function(MSEobj,nam=NA){
  # LTY<-rep(NA,MSEobj@nMPs)
  # STY<-rep(NA,MSEobj@nMPs)
  # VY<-rep(NA,MSEobj@nMPs)
  # B10<-rep(NA,MSEobj@nMPs)
  # yend<-max(MSEobj@proyears-4,1):MSEobj@proyears
  # ystart<-1:5
  # RefYd<-MSEobj@OM$RefY
  # y1<-1:(MSEobj@proyears-1)
  # y2<-2:MSEobj@proyears
  # for(mm in 1:MSEobj@nMPs){
    # LTY[mm]<-round(sum(MSEobj@C[,mm,yend]/RefYd>0.5,na.rm=T)/(MSEobj@nsim*length(yend)),3)*100
    # STY[mm]<-round(sum(MSEobj@C[,mm,ystart]/RefYd>0.5,na.rm=T)/(MSEobj@nsim*length(ystart)),3)*100
    # AAVY<-apply(((MSEobj@C[,mm,y1]-MSEobj@C[,mm,y2])^2)^0.5,1,mean,na.rm=T)/apply(MSEobj@C[,mm,y2],1,mean,na.rm=T)
    # VY[mm]<-round(sum(AAVY<0.1,na.rm=T)/MSEobj@nsim,3)*100
    # B10[mm]<-round(sum(MSEobj@B_BMSY[,mm,]>0.1,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,])),3)*100
  # }
  # par(mfrow=c(1,2),mai=c(1.5,1.5,0.1,0.1),omi=c(0.1,0.1,0.4,0))
  # tradeoffplot(STY,LTY,"P(Short term yield over half FMSY)","P(Long term yield over half FMSY)",MSEobj@MPs[1:MSEobj@nMPs],vl=1,hl=1)
  # tradeoffplot(B10,VY,"P(Biomass > 0.1 BMSY)", "P(CV in yield less than 0.1)",MSEobj@MPs[1:MSEobj@nMPs],vl=1,hl=1)
  # if(is.na(nam))mtext(deparse(substitute(MSEobj)),3,outer=T,line=0.3,font=2)
  # if(!is.na(nam))mtext(MSEobj@Name,3,outer=T,line=0.3,font=2)
# }


# NOAA_plot<-function(MSEobj,nam=NA,type=NA,panel=T){
  
  # Yd<-rep(NA,MSEobj@nMPs)
  # B50<-rep(NA,MSEobj@nMPs)
  # PNOF<-rep(NA,MSEobj@nMPs)
  # LTY<-rep(NA,MSEobj@nMPs)
  # STY<-rep(NA,MSEobj@nMPs)
  # VY<-rep(NA,MSEobj@nMPs)
  
  # y1<-1:(MSEobj@proyears-1)
  # y2<-2:MSEobj@proyears
  
  # yend<-max(MSEobj@proyears-4,1):MSEobj@proyears
  
  # RefYd<-MSEobj@OM$RefY
  
  # for(mm in 1:MSEobj@nMPs){
    
    # PNOF[mm]<-round(sum(MSEobj@F_FMSY[,mm,]<1,na.rm=T)/prod(dim(MSEobj@F_FMSY[,mm,]),na.rm=T)*100,1)
    # B50[mm]<-round(sum(MSEobj@B_BMSY[,mm,]>0.5,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
    # LTY[mm]<-round(sum(MSEobj@C[,mm,yend]/RefYd>0.5,na.rm=T)/(MSEobj@nsim*length(yend)),3)*100
    # AAVY<-apply((((MSEobj@C[,mm,y1]-MSEobj@C[,mm,y2])/MSEobj@C[,mm,y2])^2)^0.5,1,mean,na.rm=T) 
    # VY[mm]<-round(sum(AAVY<0.15,na.rm=T)/MSEobj@nsim,3)*100

  # }
  
  # #dev.new2(width=7,height=7)
  # if(panel)par(mfrow=c(1,2),mai=c(1.5,1.5,0.1,0.1),omi=c(0.1,0.1,0.4,0))
  
  # if(is.na(type)){
    # tradeoffplot(PNOF,LTY,"Prob. of not overfishing (%)", "Long-term yield ",MSEobj@MPs[1:MSEobj@nMPs],vl=50,hl=100)
    # tradeoffplot(B50,VY,"Prob. biomass above half BMSY (%)", "Prob. AAVY less than 15%",MSEobj@MPs[1:MSEobj@nMPs],vl=80,hl=50)
  # }else{
    # tradeoffplot3(PNOF,LTY,"Prob. of not overfishing (%)", "Long-term yield",MSEobj@MPs[1:MSEobj@nMPs],vl=50,hl=100,xlim=c(45,105),ylim=c(0,105))
    # tradeoffplot3(B50,VY,"Prob. biomass above half BMSY (%)", "Prob. AAVY less than 15%",MSEobj@MPs[1:MSEobj@nMPs],vl=80,hl=50,xlim=c(75,105),ylim=c(45,105))
  # }
  
  # #if(is.na(nam))mtext(deparse(substitute(MSEobj)),3,outer=T,line=0.3,font=2)
  # #if(!is.na(nam) & !is.character(nam))mtext(MSEobj@Name,3,outer=T,line=0.3,font=2)
  # #if(!is.na(nam) & is.character(nam))mtext(nam,3,outer=T,line=0.3,font=2)
  
  
  # temp<-data.frame(PNOF,B50,LTY,VY)
  # row.names(temp)<-MSEobj@MPs[1:MSEobj@nMPs]
  # temp
  
# }



# Pplot<-function(MSEobj,nam=NA){
  
  # FMSYr<-quantile(MSEobj@F_FMSY,c(0.001,0.90),na.rm=T)
  # BMSYr<-quantile(MSEobj@B_BMSY,c(0.001,0.975),na.rm=T)
  
  # colsse<-rainbow(100,start=0,end=0.36)[1:100]
  # colB<-rep(colsse[100],ceiling(BMSYr[2]*100))
  # colB[1:100]<-colsse
  # colB<-makeTransparent(colB,60)
  # colsse<-rainbow(200,start=0,end=0.36)[200:1]
  # colF<-rep(colsse[200],ceiling(FMSYr[2]*100))
  # colF[1:200]<-colsse
  # colF<-makeTransparent(colF,60)
  
  # Yd<-rep(NA,MSEobj@nMPs)
  # P10<-rep(NA,MSEobj@nMPs)
  # P50<-rep(NA,MSEobj@nMPs)
  # P100<-rep(NA,MSEobj@nMPs)
  # POF<-rep(NA,MSEobj@nMPs)
  # yind<-max(MSEobj@proyears-4,1):MSEobj@proyears
  # RefYd<-MSEobj@OM$RefY
  
  # for(mm in 1:MSEobj@nMPs){
    # Yd[mm]<-round(mean(apply(MSEobj@C[,mm,yind],1,mean,na.rm=T)/RefYd,na.rm=T)*100,1)
    # #cbind(MSEobj@C[,mm,yind],unlist(MSEobj@OM$MSY))
    # POF[mm]<-round(sum(MSEobj@F_FMSY[,mm,]>1,na.rm=T)/prod(dim(MSEobj@F_FMSY[,mm,]),na.rm=T)*100,1)
    # P10[mm]<-round(sum(MSEobj@B_BMSY[,mm,]<0.1,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
    # P50[mm]<-round(sum(MSEobj@B_BMSY[,mm,]<0.5,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
    # P100[mm]<-round(sum(MSEobj@B_BMSY[,mm,]<1,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
  # }
    
  # nr<-ceiling(MSEobj@nMPs/8)
  # nc<-ceiling(MSEobj@nMPs/nr)
  # nr<-nr*2
  # MSEcols<-c('red','green','blue','orange','brown','purple','dark grey','violet','dark red','pink','dark blue','grey')
  # temp<-array(0,c(nr*2+(nr/2-1),nc*2))
  # i<-0
  # for(c in 1:nc){
    # for(r in 1:nr){
      # i<-i+1
      # temp[(ceiling(r/2)-1)+(1:2)+(r-1)*2,(1:2)+(c-1)*2]<-((c-1)*nr)+r
    # }
  # }
  # par(mfcol=c(nr,nc),mai=c(0.2,0.35,0.3,0.01),omi=c(0.5,0.4,0.4,0.05))
  # layout(temp)
  # #dev.new2(width=nc*3,height=nr*3)
  # #
  # lwdy<-2.5

  # for(mm in 1:MSEobj@nMPs){
    # plot(MSEobj@F_FMSY[1,mm,],ylim=FMSYr,col=colF[ceiling(mean(MSEobj@F_FMSY[1,mm,],na.rm=T)*100)],type='l',lwd=lwdy)
    # for(i in 1:MSEobj@nsim)lines(MSEobj@F_FMSY[i,mm,],col=colF[ceiling(mean(MSEobj@F_FMSY[i,mm,],na.rm=T)*100)],lwd=lwdy)
    # abline(h=100,col="grey",lwd=3)
    # mtext(MSEobj@MPs[mm],3,outer=F,line=0.6)
    # legend('topright',c(paste(POF[mm],"% POF",sep=""),
                 # paste(Yd[mm],"% FMSY yield",sep="")),bty='n',cex=0.8)
    # if(mm%in%(1:(nr/2)))mtext("F/FMSY",2,line=2.5,outer=F)
    # abline(h=1,col=makeTransparent("grey",30),lwd=2.5)
  
    # plot(MSEobj@B_BMSY[1,mm,],ylim=BMSYr,col=colB[ceiling(MSEobj@B_BMSY[1,mm,MSEobj@proyears]*100)],type='l',lwd=lwdy)
    # for(i in 1:MSEobj@nsim)lines(MSEobj@B_BMSY[i,mm,],col=colB[ceiling(MSEobj@B_BMSY[i,mm,MSEobj@proyears]*100)],lwd=lwdy)
    # abline(h=100,col="grey",lwd=3)
    # legend('topright',c(paste(P100[mm],"% < BMSY",sep=""),
                 # paste(P50[mm],"% < 0.5BMSY",sep=""),
                 # paste(P10[mm],"% < 0.1BMSY",sep="")),bty='n',cex=0.8)
    # if(mm%in%(1:(nr/2)))mtext("B/BMSY",2,line=2.5,outer=F)
    # abline(h=1,col=makeTransparent("grey",30),lwd=2.5)
  
  # }
  # mtext("Projection year",1,outer=T,line=1.2)
  # if(is.na(nam))mtext(deparse(substitute(MSEobj)),3,outer=T,line=0.3,font=2)
  # if(!is.na(nam))mtext(MSEobj@Name,3,outer=T,line=0.3,font=2)
# }



# Kplot<-function(MSEobj,maxsim=60,nam=NA){
  # nr<-floor((MSEobj@nMPs)^0.5)
  # nc<-ceiling((MSEobj@nMPs)/nr)
    
  # FMSYr<-quantile(MSEobj@F_FMSY,c(0.001,0.90),na.rm=T)
  # BMSYr<-quantile(MSEobj@B_BMSY,c(0.001,0.975),na.rm=T)
    
  # #dev.new2(width=nc*3,height=nr*3.6)
  # par(mfrow=c(nr,nc),mai=c(0.45,0.45,0.45,0.01),omi=c(0.45,0.3,0.35,0.01))
  
  # colsse<-rainbow(MSEobj@proyears,start=0.63,end=0.95)[1:MSEobj@proyears]
  # colsse<-makeTransparent(colsse,95)
  
  # for(mm in 1:MSEobj@nMPs){
    # plot(c(MSEobj@B_BMSY[1,mm,1],MSEobj@B_BMSY[1,mm,2]),
         # c(MSEobj@F_FMSY[1,mm,1],MSEobj@F_FMSY[1,mm,2]),xlim=BMSYr,ylim=FMSYr,
         # col=colsse[1],type='l')
    
    # OO<-round(sum(MSEobj@B_BMSY[,mm,MSEobj@proyears]<1&MSEobj@F_FMSY[,mm,MSEobj@proyears]>1,na.rm=T)/MSEobj@nsim*100,1)
    # OU<-round(sum(MSEobj@B_BMSY[,mm,MSEobj@proyears]>1&MSEobj@F_FMSY[,mm,MSEobj@proyears]>1,na.rm=T)/MSEobj@nsim*100,1)
    # UO<-round(sum(MSEobj@B_BMSY[,mm,MSEobj@proyears]<1&MSEobj@F_FMSY[,mm,MSEobj@proyears]<1,na.rm=T)/MSEobj@nsim*100,1)
    # UU<-round(sum(MSEobj@B_BMSY[,mm,MSEobj@proyears]>1&MSEobj@F_FMSY[,mm,MSEobj@proyears]<1,na.rm=T)/MSEobj@nsim*100,1)
    
    # #alp<-80
    # #polygon(c(1,-1000,-1000,1),c(1,1,1000,1000),col=makeTransparent("orange",alp),border=makeTransparent("orange",alp))
    # #polygon(c(1,1000,1000,1),c(1,1,1000,1000),col=makeTransparent("yellow",alp),border=makeTransparent("yellow",alp))
    # #polygon(c(1,-1000,-1000,1),c(1,1,-1000,-1000),col=makeTransparent("yellow",alp),border=makeTransparent("yellow",alp))
    # #polygon(c(1,1000,1000,1),c(1,1,-1000,-1000),col=makeTransparent("green",alp),border=makeTransparent("yellow",alp))
    
    
    # abline(h=1,col="grey",lwd=3)
    # abline(v=1,col="grey",lwd=3)
    # #abline(v=c(0.1,0.5),col="grey",lwd=2)
    # rng<-1:min(maxsim,MSEobj@nsim)
    # for(i in rng){
      # for(y in 1:(MSEobj@proyears-1)){
        # lines(c(MSEobj@B_BMSY[i,mm,y],MSEobj@B_BMSY[i,mm,y+1]),
              # c(MSEobj@F_FMSY[i,mm,y],MSEobj@F_FMSY[i,mm,y+1]),
              # col=colsse[y],lwd=1.6)
      # }
    # }
    
    # points(MSEobj@B_BMSY[rng,mm,1],MSEobj@F_FMSY[rng,mm,1],pch=19,cex=0.8,col=colsse[1])
    # points(MSEobj@B_BMSY[rng,mm,MSEobj@proyears],MSEobj@F_FMSY[rng,mm,MSEobj@proyears],pch=19,cex=0.8,col=colsse[MSEobj@proyears])
    
    # if(mm==1)legend('right',c("Start","End"),bty='n',text.col=c(colsse[1],colsse[MSEobj@proyears]),pch=19,col=c(colsse[1],colsse[MSEobj@proyears]))
    # legend('topleft',paste(OO,"%",sep=""),bty='n',text.font=2)
    # legend('topright',paste(OU,"%",sep=""),bty='n',text.font=2)
    # legend('bottomleft',paste(UO,"%",sep=""),bty='n',text.font=2)
    # legend('bottomright',paste(UU,"%",sep=""),bty='n',text.font=2)
    
    # mtext(MSEobj@MPs[mm],3,line=0.45)
  # }
  # mtext("B/BMSY",1,outer=T,line=1.4)
  # mtext("F/FMSY",2,outer=T,line=0.2)
  # if(is.na(nam))mtext(deparse(substitute(MSEobj)),3,outer=T,line=0.25,font=2)
  # if(!is.na(nam))mtext(MSEobj@Name,3,outer=T,line=0.25,font=2)
# }

# # Plotting code for MSE object
# setMethod("plot",
  # signature(x = "MSE"),
  # function(x){
  # MSEobj<-x
  # Pplot(MSEobj)
  # Kplot(MSEobj)
  # Tplot(MSEobj)
# })

# # Value of information analysis
# VOI<-function(MSEobj,ncomp=6,nbins=8,maxrow=8,Ut=NA,Utnam="Utility"){# Value of information  ====================================================================
  # objnam<-deparse(substitute(MSEobj))
  # nsim<-MSEobj@nsim
 
  # if(is.na(Ut[1])){
    # Ut<-array(NA,c(nsim,MSEobj@nMPs))
    # yind<-max(MSEobj@proyears-4,1):MSEobj@proyears
    # RefYd<-MSEobj@OM$RefY
  
    # for(mm in 1:MSEobj@nMPs){
      # Ut[,mm]<-apply(MSEobj@C[,mm,yind],1,mean,na.rm=T)/RefYd*100
    # #POF[,mm]<-apply(MSEobj@F_FMSY[,mm,]>1,1,sum)/MSEobj@proyears
    # #P10[,mm]<-apply(MSEobj@B_BMSY[,mm,]<0.1,1,sum)/MSEobj@proyears
    # }
    # Utnam<-"Long-term yield relative to MSY (%)"
  # }

  # MPs<-MSEobj@MPs
  # nMPs<-MSEobj@nMPs
 
  # onlycor<-c("RefY","A","MSY","Linf","t0","OFLreal","Spat_targ")
  # vargood<-(apply(MSEobj@OM,2,sd)/(apply(MSEobj@OM,2,mean)^2)^0.5)>0.005
  # # MSEobj@OM<- MSEobj@OM[,(!names(MSEobj@OM)%in%onlycor)&vargood]
  # MSEobj@OM[,which((!names(MSEobj@OM)%in%onlycor)&vargood)]  
  # OMp<-apply(MSEobj@OM,2,quantile,p=seq(0,1,length.out=nbins+1))
  # Obsp<-apply(MSEobj@Obs,2,quantile,p=seq(0,1,length.out=nbins+1))
  # OMv<-array(NA,c(nMPs,ncol(MSEobj@OM),nbins))
  # Obsv<-array(NA,c(nMPs,ncol(MSEobj@Obs),nbins))
  
  # for(mm in 1:nMPs){
    # for(j in 1:nbins){
      # for(i in 1:ncol(MSEobj@OM)){
        # cond<-MSEobj@OM[,i]>OMp[j,i]&MSEobj@OM[,i]<OMp[j+1,i]
        # OMv[mm,i,j]<-mean(Ut[cond,mm],na.rm=T)
      # }
      # for(i in 1:ncol(MSEobj@Obs)){
        # cond<-MSEobj@Obs[,i]>Obsp[j,i]&MSEobj@Obs[,i]<Obsp[j+1,i]
        # Obsv[mm,i,j]<-mean(Ut[cond,mm],na.rm=T)
      # }
    # }
  # }
  
  # # -- Operating model variables
  # OMs<-apply(OMv,1:2,sd,na.rm=T)
  # OMstr<-array("",c(nMPs*2,ncomp+1))
   
  # for(mm in 1:nMPs){
    # ind<-order(OMs[mm,],decreasing=T)[1:ncomp]
    # OMstr[1+(mm-1)*2,1]<-MPs[mm]
    # OMstr[1+(mm-1)*2,2:(1+ncomp)]<-names(MSEobj@OM[ind])
    # OMstr[2+(mm-1)*2,2:(1+ncomp)]<-round(OMs[mm,ind],2)
  # }
  # OMstr<-data.frame(OMstr)
  # names(OMstr)<-c("MP",1:ncomp)
  
  
  # # -- Observation model variables
  # slots<-c( "Cat",  "Cat","AvC",  "AvC","CAA",      "CAA",    "CAL",      "CAL",    "Ind","Dep",  "Dep", "Dt",   "Dt", "Mort", "FMSY_M",    "BMSY_B0",     "L50",      "L95",    "LFC",    "LFS",    "Abun",  "Abun","vbK",  "vbt0",  "vbLinf",  "Steep","Iref",    "Cref",    "Bref")
  # Obsnam<-c("Cbias","Csd","Cbias","Csd","CAA_nsamp","CAA_ESS","CAL_nsamp","CAL_ESS","Isd","Dbias","Derr","Dbias","Derr","Mbias","FMSY_Mbias","BMSY_B0bias", "lenMbias","lenMbias","LFCbias","LFSbias","Abias","Aerr","Kbias","t0bias","Linfbias","hbias","Irefbias","Crefbias","Brefbias")
  # Obss<-apply(Obsv,1:2,sd,na.rm=T)
  # Obsstr<-array("",c(nMPs*2,ncomp+1))
  # for(mm in 1:nMPs){
    # relobs<-Obsnam[slots%in%unlist(strsplit(Required(MPs[mm])[,2],split=", "))]
    # ind<-(1:ncol(MSEobj@Obs))[match(relobs,names(MSEobj@Obs))]
    # pos<-names(MSEobj@Obs)[ind]# possible observation processes
    # maxy<-min(max(1,length(pos)),ncomp,na.rm=T)
    # ind2<-order(Obss[mm,ind],decreasing=T)[1:maxy]
    # Obsstr[1+(mm-1)*2,1]<-MPs[mm]
    # Obsstr[1+(mm-1)*2,2:(1+maxy)]<-pos[ind2]
    # Obsstr[2+(mm-1)*2,2:(1+maxy)]<-round(Obss[mm,ind][ind2],2)
  # }
  # Obsstr<-data.frame(Obsstr)
  # names(Obsstr)<-c("MP",1:ncomp)

  # ncols<-40
  # #colsse<-makeTransparent(rainbow(ncols,start=0,end=0.36),95)[ncols:1]
  # colsse<-makeTransparent(rainbow(ncols,start=0,end=0.36),90)[ncols:1]
  # minsd<-0
  # maxsd<-max(OMs,na.rm=T)
  # coly<-ceiling(OMs/maxsd*ncols)
 
  # # Operating model variables
  # mbyp <- split(1:nMPs, ceiling(1:nMPs/maxrow))
  # ylimy=c(0,max(OMv,na.rm=T)*1.2)              
  
  
  # for(pp in 1:length(mbyp)){
    
    # par(mfrow=c(length(mbyp[[pp]]),ncomp),mai=c(0.15,0.1,0.15,0.05),omi=c(0.1,0.9,0.3,0.05))
   
    # for(mm in mbyp[[pp]]){
      # for(cc in 1:ncomp){
        # rind<-(mm-1)*2+1
        # y<-Ut[,mm]
        # cind<-match(OMstr[rind,1+cc],names(MSEobj@OM))
        # x<-MSEobj@OM[,cind]
        # plot(x,y,col="white",axes=F,ylim=ylimy)
        # axis(1,pretty(OMp[,cind]),pretty(OMp[,cind]),cex.axis=0.8,padj=-1.5)
        # abline(v=OMp[,cind],col="#99999960")
        # points(x,y,col=colsse[coly[mm,cind]],pch=19,cex=0.8)
        # x2<-(OMp[1:nbins,cind]+OMp[2:(nbins+1),cind])/2
        # y2<-OMv[mm,cind,]
        # lines(x2,y2)
        # legend('bottomright',legend=round(OMs[mm,cind],2),bty='n',cex=0.8)
        # legend('topleft',legend=OMstr[rind,1+cc],bty='n',cex=0.85)
        # if(cc==1){ 
          # mtext(MPs[mm],2,font=2,outer=F,cex=0.8,line=2)
          # ytick<-pretty(seq(ylimy[1],ylimy[2]*1.3,length.out=10))
          # axis(2,ytick,ytick,cex.axis=0.8)
        # } # only first column
      # } # parameters (columns)
    # } # MPs (rows)
    
    # mtext(Utnam,2,outer=T,cex=0.9,line=3.5)
    # mtext(paste("Operating model parameters: ",objnam,"@OM",sep=""),3,outer=T,font=2,cex=0.9)
    
  # } # Plots
  
  # # Observation model values
  
  # ylimy=c(0,max(Obsv,na.rm=T)*1.2)              
  # minsd<-0
  # maxsd<-max(Obss)
  # coly<-ceiling(Obss/maxsd*ncols)
  
  # if(sum(is.na(Obsstr)|Obsstr=="")<(ncomp+1)*nMPs*2-nMPs){ # only if there is data to plot
  
  # for(pp in 1:length(mbyp)){
    
    # par(mfrow=c(length(mbyp[[pp]]),ncomp),mai=c(0.15,0.1,0.15,0.05),omi=c(0.1,0.9,0.3,0.05))
    
    # for(mm in mbyp[[pp]]){
      # rind<-(mm-1)*2+1
      # npres<-sum(Obsstr[rind+1,]!="")
      # for(cc in 1:ncomp){
        # if(!is.na(npres)&cc<(npres+1)){
          # y<-Ut[,mm]
          # cind<-match(Obsstr[rind,1+cc],names(MSEobj@Obs))
          # x<-MSEobj@Obs[,cind]
          # plot(x,y,col="white",axes=F,ylim=ylimy)
          # axis(1,pretty(Obsp[,cind]),pretty(Obsp[,cind]),cex.axis=0.8,padj=-2)
          # abline(v=Obsp[,cind],col="#99999960")
          # points(x,y,col=colsse[coly[mm,cind]],pch=19,cex=0.8)
          # x2<-(Obsp[1:nbins,cind]+Obsp[2:(nbins+1),cind])/2
          # y2<-Obsv[mm,cind,]
          # lines(x2,y2)
          # legend('bottomright',legend=round(Obss[mm,cind],2),bty='n',cex=0.8)
          # legend('topleft',legend=Obsstr[rind,1+cc],bty='n',cex=0.75)
          # if(cc==1){ 
            # mtext(MPs[mm],2,font=2,outer=F,cex=0.6,line=2)
            # ytick<-pretty(seq(ylimy[1],ylimy[2]*1.3,length.out=10))
            # axis(2,ytick,ytick,cex.axis=0.8)
          # } # only first column
        # }else{
          # plot(0,type='n',axes=FALSE,ann=FALSE)
          # if(cc==1){ 
            # mtext(MPs[mm],2,font=2,outer=F,cex=0.6,line=2)
          # } # only first column
        # }  
      # } # parameters (columns)
    # } # MPs (rows)
    
    # mtext(Utnam,2,outer=T,cex=0.9,line=3.5)
    # mtext(paste("Observation model parameters: ",objnam,"@Obs",sep=""),3,outer=T,font=2,cex=0.9)
   
  # } # Plots
  # } # if there is data to plot
  
  # list(OMstr,Obsstr)
  
# } # VOI

# tradeoffplot<-function(x,y,xlab,ylab,labs,cex,vl,hl){
   # adjj<-c(0.7,1.3)
   # XLim <- c(min(c(-10, min(x,na.rm=T)*adjj)), max(c(max(x,na.rm=T)*adjj, 110)))
   # YLim <- c(min(c(-10, min(y,na.rm=T)*adjj)), max(c(max(y,na.rm=T)*adjj, 110)))
   # coly<-rep(c('#0000ff95','#ff000095','#20ff1095'),50)[1:length(labs)]
   # coly[labs%in%c("AvC","curE","FMSYref")]<-'black'
   # # plot(NA,xlim=range(x,na.rm=T)*adjj,ylim=range(y,na.rm=T)*adjj,xlab=xlab,ylab=ylab)
   # plot(NA,xlim=XLim,ylim=YLim,xlab=xlab,ylab=ylab)
   # abline(v=vl,col="#99999940",lwd=2)
   # abline(h=hl,col="#99999940",lwd=2)
   # text(x,y,labs,font=2,col=coly,cex=0.9)
# }

# tradeoffplot2<-function(x,y,xlab,ylab,cex=1,vl,hl,coly,leg){
  # adjj<-c(0.7,1.3)
  # plot(NA,xlim=range(x,na.rm=T)*adjj,ylim=range(y,na.rm=T)*adjj,xlab=xlab,ylab=ylab)
  # abline(v=vl,col="grey",lwd=2)
  # abline(h=hl,col="grey",lwd=2)
  # for(m in 1:nrow(x))points(x[m,],y[m,],col=makeTransparent(coly[m],50),pch=19,cex=cex)
  # if(!is.na(leg[1]))legend('topright',legend=leg,text.col=coly,bty='n')
# }

# tradeoffplot3<-function(x,y,xlab,ylab,labs,cex,vl,hl,xlim,ylim){
  # coly<-rep(c('#0000ff95','#ff000095','#20ff1095'),10)
  # plot(NA,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
  # abline(v=vl,col="#99999940",lwd=2)
  # abline(h=hl,col="#99999940",lwd=2)
  # text(x,y,labs,font=2,col=coly,cex=1)
# }


# makeTransparent<-function(someColor, alpha=100){
  # newColor<-col2rgb(someColor)
  # apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    # blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
# }

# getmov<-function(x,Prob_staying,Frac_area_1){
  # test<-optim(par=c(0,0,0),movfit,method="L-BFGS-B",lower=rep(-6,3),upper=rep(6,3),prb=Prob_staying[x],frac=Frac_area_1[x])
  # mov<-array(c(test$par[1],test$par[2],0,test$par[3]),dim=c(2,2))
  # mov<-exp(mov)
  # mov/array(apply(mov,1,sum),dim=c(2,2))
# }

# movfit<-function(par,prb,frac){
  # mov<-array(c(par[1],par[2],0,par[3]),dim=c(2,2))
  # mov<-exp(mov)
  # mov<-mov/array(apply(mov,1,sum),dim=c(2,2))
  # dis<-c(frac,1-frac)
  # for(i in 1:100)dis<-apply(array(dis,c(2,2))*mov,2,sum)
  # (log(mov[1,1])-log(prb))^2+(log(frac)-log(dis[1]))^2
# }

# getq<-function(x,dep,Find,Perr,Marray,hs,Mat_age,Wt_age,R0,V,nyears,maxage,mov,Spat_targ,SRrel,aR,bR){
  # opt<-optimize(qopt,log(c(0.0075,15)),depc=dep[x],Fc=Find[x,],Perrc=Perr[x,],
                     # Mc=Marray[x,],hc=hs[x],Mac=Mat_age[x,],Wac=Wt_age[x,,],
                     # R0c=R0,Vc=V[x,,],nyears=nyears,maxage=maxage,movc=mov[x,,],
                     # Spat_targc=Spat_targ[x],SRrelc=SRrel[x],aRc=aR[x,],bRc=bR[x,])
  # return(exp(opt$minimum))
# }


# qopt<-function(lnq,depc,Fc,Perrc,Mc,hc,Mac,Wac,R0c,Vc,nyears,maxage,movc,Spat_targc,SRrelc,aRc,bRc,opt=T){
  # qc<-exp(lnq)
  # nareas<-nrow(movc)
  # #areasize<-c(asizec,1-asizec)
  # idist<-rep(1/nareas,nareas)
  # for(i in 1:300)idist<-apply(array(idist,c(2,2))*movc,2,sum)

  # N<-array(exp(-Mc[1]*((1:maxage)-1))*R0c,dim=c(maxage,nareas))*array(rep(idist,each=maxage),dim=c(maxage,nareas))
  # SSN<-Mac*N   # Calculate initial spawning stock numbers
  # Biomass<-N*Wac[,1]
  # SSB<-SSN*Wac[,1]                               # Calculate spawning stock biomass

  # B0<-sum(Biomass)
  # R0a<-idist*R0c
  # SSB0<-apply(SSB,2,sum)
  # SSBpR<-SSB0/R0a                              # Calculate spawning stock biomass per recruit

  # for(y in 1:nyears){
    # # set up some indices for indexed calculation
    # targ<-(apply(Vc[,y]*Biomass,2,sum)^Spat_targc)/mean(apply(Vc[,y]*Biomass,2,sum)^Spat_targc)
    # FMc<-array(qc*Fc[y]*Vc[,y],dim=c(maxage,nareas))*array(rep(targ,each=maxage),dim=c(maxage,nareas))                                           # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
    # Zc<-FMc+Mc[y]
    # N[2:maxage,]<-N[1:(maxage-1),]*exp(-Zc[1:(maxage-1),])         # Total mortality
    # if(SRrelc==1){
      # N[1,]<-Perrc[y]*(0.8*R0a*hc*apply(SSB,2,sum))/(0.2*SSBpR*R0a*(1-hc)+(hc-0.2)*apply(SSB,2,sum))  # Recruitment assuming regional R0 and stock wide steepness
    # }else{
      # N[1,]<- aRc*apply(SSB,2,sum)*exp(-bRc*apply(SSB,2,sum)) 
    # }
      
    # #print(N[1])
    # indMov<-as.matrix(expand.grid(1:nareas,1:nareas,1:maxage)[3:1])
    # indMov2<-indMov[,1:2]
    # indMov3<-indMov[,2:3]
    # temp<-array(N[indMov2]*movc[indMov3],dim=c(nareas,nareas,maxage))
    # N<-apply(temp,c(3,1),sum)
    # SSN<-N*Mac
    # SSB<-SSN*Wac[,y]
    # Biomass<-N*Wac[,y]
    # SBiomass<-SSN*Wac[,y]
    # #print(sum(Biomass))
  # } # end of year
  # return((log(depc)-log(sum(SBiomass)/sum(SSB0)))^2)
# }

# getinitdist<-function(tol,mov,indMain){
  # init<-array(1/nareas,dim=c(nsim,nareas))
  # ind4<-as.matrix(cbind(rep(1:nsim,each=nareas*nareas),indMain[,2]))
  # i<-0
  # delta<-1
  # #for(i in 1:100){
  # while(delta > tol){
    # i<-i+1
    # trial<-init
    # temp<-array(init[ind4]*mov[indMain],dim=c(nareas,nareas,nsim))
    # init<-apply(temp,c(3,1),sum)
    # delta<-max((trial-init)^2)
  # }
  # print(paste("Converged in ",i," iterations"))
  # init
# }

# getFMSY<-function(x,Marray,hs,Mat_age,Wt_age,R0,V,maxage,nyears,proyears,Spat_targ,mov,SRrel,aR,bR){
  # opt<-optimize(FMSYopt,log(c(0.001,5)),
                     # Mc=Marray[x,nyears],hc=hs[x],Mac=Mat_age[x,],Wac=Wt_age[x,,nyears],
                     # R0c=R0,Vc=V[x,],maxage=maxage,nyears=nyears,proyears=proyears,Spat_targc=Spat_targ[x],movc=mov[x,,],SRrelc=SRrel[x],aRc=aR[x,],bRc=bR[x,],Opt=T)
				
  # return(FMSYopt(opt$minimum,
                     # Mc=Marray[x,nyears],hc=hs[x],Mac=Mat_age[x,],Wac=Wt_age[x,,nyears],
                     # R0c=R0,Vc=V[x,],maxage=maxage,nyears=nyears,proyears=proyears,Spat_targc=Spat_targ[x],movc=mov[x,,],SRrelc=SRrel[x],aRc=aR[x,],bRc=bR[x,],Opt=F)
         # )
# }

# FMSYopt<-function(lnF,Mc,hc,Mac,Wac,R0c,Vc,maxage,nyears,proyears,Spat_targc,movc,SRrelc,aRc,bRc,Opt=T){

  # FMSYc<-exp(lnF)
  # nareas<-nrow(movc)
  # #areasize<-c(asizec,1-asizec)
  # idist<-rep(1/nareas,nareas)
  # for(i in 1:100)idist<-apply(array(idist,c(2,2))*movc,2,sum)

  # N<-array(exp(-Mc*((1:maxage)-1))*R0c,dim=c(maxage,nareas))*array(rep(idist,each=maxage),dim=c(maxage,nareas))
  # SSN<-Mac*N   # Calculate initial spawning stock numbers
  # Biomass<-N*Wac
  # VBiomass<-Biomass*Vc
  # SSB<-SSN*Wac                              # Calculate spawning stock biomass
  
  # B0 <- sum(Biomass)
  # VB0<-sum(VBiomass)
  # R0a<-idist*R0c
  # SSB0<-apply(SSB,2,sum)
  # SSBpR<-SSB0/R0a

  # N<-N/2                              # Calculate spawning stock biomass per recruit
  # SSN<-Mac*N   # Calculate initial spawning stock numbers
  # Biomass<-N*Wac
  # SSB<-SSN*Wac                              # Calculate spawning stock biomass

  # for(y in 1:nyears){
    # # set up some indices for indexed calculation
    # dis<-apply(Vc*Biomass,2,sum)/sum(Vc*Biomass)
    # targ<-(dis^Spat_targc)/mean(dis^Spat_targc)
    # FMc<-array(FMSYc*Vc,dim=c(maxage,nareas))*array(rep(targ,each=maxage),dim=c(maxage,nareas))                                           # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
    # Zc<-FMc+Mc
    # CN<-N*(1-exp(-Zc))*(FMc/Zc)
    # CB<-CN*Wac

    # N[2:maxage,]<-N[1:(maxage-1),]*exp(-Zc[1:(maxage-1),])         # Total mortality
    # if(SRrelc==1){
      # N[1,]<-(0.8*R0a*hc*apply(SSB,2,sum))/(0.2*SSBpR*R0a*(1-hc)+(hc-0.2)*apply(SSB,2,sum))  # Recruitment assuming regional R0 and stock wide steepness
    # }else{
      # N[1,]<- aRc*apply(SSB,2,sum)*exp(-bRc*apply(SSB,2,sum)) 
    # }
    # #print(N[1])
    # N[1,]<-apply(array(N[1,],c(2,2))*movc,2,sum)
    # SSN<-N*Mac
    # SSB<-SSN*Wac
    # Biomass<-N*Wac
    # VBiomass<-Biomass*Vc
    # #print(sum(Biomass))
  # } # end of year
 
  # CBc<-sum(CB)
  # if(Opt){
    # return(-CBc)
  # }else{
    # return(c(CBc,-log(1-(CBc/(sum(VBiomass)+CBc))),sum(SSB)/sum(SSB0)))
  # }
# }

# Sam<-function(DLM_data,MPs=NA,reps=100,maxlines=10,perc=0.5){
  # nm <-deparse(substitute(DLM))
  # DLM_data@PosMPs<-MPs
  # funcs<-DLM_data@PosMPs
  # nMPs<-length(funcs)
  # DLM_data@MPs<-funcs
  # TACa<-getTAC(DLM_data,MPs=funcs,reps)
  # nsim<-length(DLM_data@Mort)
  # ref<-array(rep(DLM_data@Ref,nMPs),c(nsim,nMPs))
  # TACm<-apply(TACa,c(3,1),quantile,p=perc,na.rm=T)
  # TACbias<-(TACm-ref)/ref *100
  # POF<-round(apply(TACbias>0,2,sum)/length(DLM_data@Mort)*100,1)
  # DLM_data@TAC<-TACa
  # DLM_data@TACbias<-TACbias
  # DLM_data
# }


# getFhist<-function(nsim,Esd,nyears,dFmin,dFmax,bb){

  # ne<-nsim*3                                                         # Number of simulated effort datasets
  # dEfinal<-runif(ne,dFmin,dFmax)#(exp(rnorm(ne,mean=demu,sd=desd))-1)*6               # Sample the final gradient in effort
  # a<-(dEfinal-bb)/nyears                                         # Derive slope to get there from intercept
  # a<-array(a,dim=c(ne,nyears))                                  # Slope array
  # bb<-array(bb,dim=c(ne,nyears))                                  # Intercept array
  # x<-array(rep(1:nyears,each=ne),dim=c(ne,nyears))              # Year array
  # dE<-a*x+bb                                                     # Change in effort
  # # E<-array(NA,dim=c(ne,nyears))                                 # Define total effort array
  # # E[,1]<-dE[,1]
  # # for(y in 2:nyears){
    # # E[,y]<-apply(dE[,1:y],1,sum)
  # # }
  # # E<-E/array(apply(E,1,mean),dim=c(ne,nyears))                  # Standardise Effort to average 1
  # E2 <- t(apply(dE, 1,cumsum))  								# Define total effort array
  # E2 <- E2/array(apply(E2,1,mean),dim=c(ne,nyears))             # Standardise Effort to average 1
  # E <- E2 # 
  # cond <- apply(E,1,min)>0
  # pos<-(1:ne)[cond]
  # pos<-pos[1:nsim]
  # #environment("dEfinal")<-asNamespace('DLMtool')#assign("dFfinal",dEfinal[pos],envir=.GlobalEnv)
  
  # E<-E[pos,]                                 # Sample only those without negative effort
  # Emu<--0.5*Esd^2
  # Eerr<-array(exp(rnorm(nyears*nsim,rep(Emu,nyears),rep(Esd,nyears))),c(nsim,nyears))
  # outy<-new('list')
  # outy[[1]]<-E*Eerr
  # outy[[2]]<-dEfinal[pos]
  # outy
# }

# densnorm<-function(sd1){   # difference in density from 0.05 given a standard deviation sd1 (sd_asc) and age at maximum vulnerability modo
  # (0.05-(dnorm(0,mod[i],sd1)/dnorm(mod[i],mod[i],sd1)))^2
# }

# densnormasc<-function(sd1,age_05,mody){
  # (0.05-(dnorm(age_05,mody,sd1)/dnorm(mody,mody,sd1)))^2
# }

# getsdasc<-function(sm,age05,mod){
  # optimize(densnormasc,interval=c(0.5,100),age_05=age05[sm],mody=mod[sm])$minimum
# }

# densnormdesc<-function(sd2,V_maxage,maxy,mody){
  # (V_maxage-(dnorm(maxy,mody,sd2)/dnorm(mody,mody,sd2)))^2
# }

# getsddesc<-function(sm,Vmaxage,maxage,mod){
  # optimize(densnormdesc,interval=c(0.5,10000),V_maxage=Vmaxage[sm],maxy=maxage,mody=mod[sm])$minimum
# }

# getDNvulnS<-function(mod,age05,Vmaxage,maxage,nsim){
  # sd_asc<-sapply(1:nsim,getsdasc,age05=age05,mod=mod)
  # sd_desc<-sapply(1:nsim,getsddesc,Vmaxage=Vmaxage,maxage=maxage,mod=mod)
  # V<-array(NA,dim=c(nsim,maxage))
  # for(i in 1:nsim){
    # V[i,1:ceiling(mod[i])]<-dnorm(1:ceiling(mod[i]),mod[i],sd_asc[i])
    # V[i,(1+ceiling(mod[i])):maxage]<-dnorm((1+ceiling(mod[i])):maxage,mod[i],sd_desc[i])
    # V[i,(1+ceiling(mod[i])):maxage]<-V[i,(1+ceiling(mod[i])):maxage]/V[i,1+ceiling(mod[i])]#/V[i,floor(mod[i])+1]
    # V[i,1:ceiling(mod[i])]<-V[i,1:ceiling(mod[i])]/dnorm(mod[i],mod[i],sd_asc[i])#,mod[i],sd_asc[i])#V[i,floor(mod[i])]

  # }
  # outy<-new('list')
  # outy[[1]]<-V
  # outy[[2]]<-mod-1.18*sd_asc
  # outy
# }

# gettempvar<-function(targ,targsd,targgrad,nyears,nsim){   # creates a time series per simulation that has gradient grad and random normal walk wiht sigma
  # mutemp<--0.5*targsd^2
  # temp<-array(1,dim=c(nsim,nyears))
  # for(i in 2:nyears){
    # temp[,i]<-temp[,i]*exp(rnorm(nsim,mutemp,targsd))
  # }
  # yarray<-array(rep((1:nyears)-1,each=nsim),dim=c(nsim,nyears))
  # temp<-temp*(1+targgrad/100)^yarray
  # targ*temp/apply(temp,1,mean)
# }


# getFref<-function(x,Marray,Wt_age,Mat_age,Perr,N_s,SSN_s, Biomass_s,VBiomass_s,SSB_s,
                  # Vn,hs,R0a,nyears,proyears,nareas,maxage,mov,SSBpR,aR,bR,SRrel){
  
  # opt<-optimize(doprojPI,log(c(0.001,10)),
                # Mvec=Marray[x,(nyears+1):(nyears+proyears)],Wac=Wt_age[x,,(nyears+1):(nyears+proyears)],Mac=Mat_age[x,],
                # Pc=Perr[x,(nyears+1):(nyears+proyears)],N_c=N_s[x,,],SSN_c=SSN_s[x,,],Biomass_c=Biomass_s[x,,],
                # VBiomass_c=VBiomass_s[x,,],SSB_c=SSB_s[x,,],Vc=Vn[x,],hc=hs[x],R0ac=R0a[x,],proyears,nareas,maxage,movc=mov[x,,],SSBpRc=SSBpR[x],aRc=aR[x,],bRc=bR[x,],SRrelc=SRrel[x])
  # # print(exp(opt$minimum))
  # return(-opt$objective)
  
# }

# doprojPI<-function(lnF,Mvec,Wac,Mac,Pc,N_c,SSN_c,Biomass_c,VBiomass_c,SSB_c,Vc,hc,R0ac,proyears,nareas,maxage,movc,SSBpRc,aRc,bRc,SRrelc){
  
  # FF<-exp(lnF)
  
  # N_P<-array(NA,dim=c(maxage,proyears,nareas))
  # Biomass_P<-array(NA,dim=c(maxage,proyears,nareas))
  # VBiomass_P<-array(NA,dim=c(maxage,proyears,nareas))
  # SSN_P<-array(NA,dim=c(maxage,proyears,nareas))
  # SSB_P<-array(NA,dim=c(maxage,proyears,nareas))
  # FM_P<-array(NA,dim=c(maxage,proyears,nareas))
  # Z_P<-array(NA,dim=c(maxage,proyears,nareas))
  # CB_P<-rep(NA,proyears)
  
  # AYR<-as.matrix(expand.grid(1:maxage,1,1:nareas))
  # YA<-as.matrix(expand.grid(1,1:maxage))         # Projection year
  # Y<-YA[,1]
  # A<-YA[,2]
  # AY<-YA[,c(2,1)]
  
  # N_P[AYR]<-N_c#[AYRL]
  # SSN_P[AYR]<-SSN_c#SSN[AYRL]
  # Biomass_P[AYR]<-Biomass_c#[AYRL]
  # VBiomass_P[AYR]<-VBiomass_c#[AYRL]
  # SSB_P[AYR]<-SSB_c#[AYRL]
  
  # FM_P[AYR]<-FF*Vc[A]
  # Z_P[AYR]<-FM_P[A]+Mvec[Y]
  
  # for(y in 2:proyears){
    
    # AY1R<-as.matrix(expand.grid(1:maxage,y-1,1:nareas))
    # AYR<-as.matrix(expand.grid(1:maxage,y,1:nareas))
    # Y<-AYR[,2]
    # A<-AYR[,1]
    # AY<-AYR[,1:2]
    # R<-AYR[,3]
    # A2YR<-as.matrix(expand.grid(2:maxage,y,1:nareas))
    # A1YR<-as.matrix(expand.grid(1:(maxage-1),y-1,1:nareas))
    # A1Y<-as.matrix(expand.grid(1:(maxage-1),y-1))
    
    # indMov<-as.matrix(expand.grid(1:nareas,1:nareas,y,1:maxage)[4:1])
    # indMov2<-indMov[,c(1,2,3)]
    # indMov3<-indMov[,c(3,4)]
    
    # N_P[A2YR]<-N_P[A1YR]*exp(-Z_P[A1Y])         # Total mortality
    
    # if(SRrelc==1){
      # N_P[1,y,]<-Pc[y]*(0.8*R0ac*hc*apply(SSB_P[,y-1,],2,sum))/(0.2*SSBpRc*R0ac*(1-hc)+(hc-0.2)*apply(SSB_P[,y-1,],2,sum))  # Recruitment assuming regional R0 and stock wide steepness
    # }else{
      # N_P[1,y,]<-Pc[y]*aRc*apply(SSB_P[,y-1,],2,sum)*exp(-bRc*apply(SSB_P[,y-1,],2,sum))  
    # }
    
    # temp<-array(N_P[indMov2]*movc[indMov3],dim=c(nareas,nareas,maxage))  # Move individuals
    # N_P[,y,]<-apply(temp,c(3,1),sum)
    
    # Biomass_P[AYR]<-N_P[AYR]*Wac[AY]                                    # Calculate biomass
    # VBiomass_P[AYR]<-Biomass_P[AYR]*Vc[A]                       # Calculate vulnerable biomass
    # SSN_P[AYR] <-N_P[AYR]*Mac[A]                                       # Calculate spawning stock numbers
    # SSB_P[AYR]<-SSN_P[AYR]*Wac[AY] # Calculate spawning stock biomass
    # FM_P[AYR]<-FF*Vc[A]
    # Z_P[AYR]<-FM_P[AYR]+Mvec[Y]
    # CNtemp<-N_P[,y,]*exp(Z_P[,y,])*(1-exp(-Z_P[,y,]))*(FM_P[,y,]/Z_P[,y,])
    # CB_P[y]<-sum(Biomass_P[,y,]*exp(Z_P[,y,])*(1-exp(-Z_P[,y,]))*(FM_P[,y,]/Z_P[,y,]))
	
	# # CB_P[y] <- sum(CNtemp*Wac[AY]) 
  # } # end of year
  # return(-mean(CB_P[(proyears-min(4,(proyears-1))):proyears],na.rm=T))
  
# }


  
  
# # Function to generate length composition for single simulation and single year
# GenLenFun <- function(NatAGTG, LenatAgeGTG, LenBin, LenMid) {
  # Nbins <- length(LenMid)
  # SizeComp <- rep(0, Nbins)
  # for (L in 1:length(LenMid)) {
    # temp <- NatAGTG
    # ind <- LenatAgeGTG <= LenBin[L+1] & LenatAgeGTG > LenBin[L]
    # temp[!ind] <- 0
    # SizeComp[L] <- SizeComp[L] + sum(temp)
  # }
  # return(SizeComp) 
# }

# comp<-function(MSEobj,MPs=NA){
  
  # if(is.na(MPs))MPs<-MSEobj@MPs
  # notm<-MPs[!(MPs%in%MSEobj@MPs)]
  # canm<-MPs[MPs%in%MSEobj@MPs]
  # if(length(notm)>0)print(paste("Methods",paste(notm,collapse=", "),"were not carried out in MSE",deparse(substitute(MSEobj)),sep=" "))
  
  # if(length(canm)==0)stop(paste('None of the methods you specified were carried out in the MSE', deparse(substitute(MSEobj)),sep=""))
  
  # if(length(canm)>4){
    # print(paste('A maximum of four methods can be compared at once. Plotting first four:',paste(canm[1:4],collapse=", "),sep=" "))
    # canm<-canm[1:4]
  # } 
  
  # mind<-match(canm,MSEobj@MPs)
  # nm<-length(mind)
  # nsim<-MSEobj@nsim
  # proyears<-MSEobj@proyears
  
  # Yd<-array(NA,c(nm,nsim))
  # P10<-array(NA,c(nm,nsim))
  # P50<-array(NA,c(nm,nsim))
  # P100<-array(NA,c(nm,nsim))
  # POF<-array(NA,c(nm,nsim))
  # yind<-max(MSEobj@proyears-4,1):MSEobj@proyears
  # RefYd<-MSEobj@OM$RefY
  
  # for(m in 1:nm){
    # mm<-mind[m]
    # Yd[m,]<-round(apply(MSEobj@C[,mm,yind],1,mean,na.rm=T)/RefYd*100,1)
    # POF[m,]<-round(apply(MSEobj@F_FMSY[,mm,]>1,1,sum,na.rm=T)/proyears*100,1)
    # P10[m,]<-round(apply(MSEobj@B_BMSY[,mm,]<0.1,1,sum,na.rm=T)/proyears*100,1)
    # P50[m,]<-round(apply(MSEobj@B_BMSY[,mm,]<0.5,1,sum,na.rm=T)/proyears*100,1)
    # P100[m,]<-round(apply(MSEobj@B_BMSY[,mm,]<1,1,sum,na.rm=T)/proyears*100,1)
  # }
  
  # MSEcols<-c('red','green','blue','orange')
  
  # #dev.new2(width=7,height=7)
  # par(mfrow=c(2,2),mai=c(0.85,0.7,0.1,0.1),omi=rep(0.01,4))
  
  # tradeoffplot2(POF,Yd,"Prob. of overfishing (%)", "Relative yield",vl=50,hl=100,coly=MSEcols,leg=NA)
  # tradeoffplot2(P100,Yd,"Prob. biomass < BMSY (%)", "Relative yield",vl=50,hl=100,coly=MSEcols,leg=canm)
  # tradeoffplot2(P50,Yd,"Prob. biomass < 0.5BMSY (%)", "Relative yield",vl=50,hl=100,coly=MSEcols,leg=NA)
  # tradeoffplot2(P10,Yd,"Prob. biomass < 0.1BMSY (%)", "Relative yield",vl=50,hl=100,coly=MSEcols,leg=NA)
  
# }

# setMethod("summary",
          # signature(object = "MSE"),
          # function(object){            

    # MSEobj<-object      
    # nm<-MSEobj@nMPs
    # nsim<-MSEobj@nsim
    # proyears<-MSEobj@proyears
    
    # Yd<-P10<-P50<-P100<-POF<-LTY<-STY<-VY<-array(NA,c(nm,nsim))
    
    # yind<-max(MSEobj@proyears-4,1):MSEobj@proyears
    # RefYd<-MSEobj@OM$RefY
    # yend<-max(MSEobj@proyears-9,1):MSEobj@proyears
    # ystart<-1:10
    # y1<-1:(MSEobj@proyears-1)
    # y2<-2:MSEobj@proyears
    
    # for(m in 1:nm){
      # Yd[m,]<-round(apply(MSEobj@C[,m,yind],1,mean,na.rm=T)/RefYd*100,1)
      # POF[m,]<-round(apply(MSEobj@F_FMSY[,m,]>1,1,sum,na.rm=T)/proyears*100,1)
      # P10[m,]<-round(apply(MSEobj@B_BMSY[,m,]<0.1,1,sum,na.rm=T)/proyears*100,1)
      # P50[m,]<-round(apply(MSEobj@B_BMSY[,m,]<0.5,1,sum,na.rm=T)/proyears*100,1)
      # P100[m,]<-round(apply(MSEobj@B_BMSY[,m,]<1,1,sum,na.rm=T)/proyears*100,1)
      # LTY[m]<-round(sum(MSEobj@C[,m,yend]/RefYd>0.5)/(MSEobj@nsim*length(yend))*100,1)
      # STY[m]<-round(sum(MSEobj@C[,m,ystart]/RefYd>0.5)/(MSEobj@nsim*length(ystart))*100,1)
      # AAVY<-apply(((MSEobj@C[,m,y1]-MSEobj@C[,m,y2])^2)^0.5,1,mean)/apply(MSEobj@C[,m,y2],1,mean)
      # VY[m]<-round(sum(AAVY<0.1)/MSEobj@nsim*100,1)
    # }
    # nr<-2
    # out<-cbind(MSEobj@MPs,round(apply(Yd,1,mean,na.rm=T),nr),round(apply(Yd,1,sd,na.rm=T),nr),
                             # round(apply(POF,1,mean,na.rm=T),nr),round(apply(POF,1,sd,na.rm=T),nr),
                             # round(apply(P10,1,mean,na.rm=T),nr),round(apply(P10,1,sd,na.rm=T),nr),
                             # round(apply(P50,1,mean,na.rm=T),nr),round(apply(P50,1,sd,na.rm=T),nr),
                             # round(apply(P100,1,mean,na.rm=T),nr),round(apply(P100,1,sd,na.rm=T),nr),
                             # round(apply(LTY,1,mean,na.rm=T),nr),
                             # round(apply(STY,1,mean,na.rm=T),nr),
                             # round(apply(VY,1,mean,na.rm=T),nr))
    # out<-as.data.frame(out)
    # names(out)<-c("MP","Yield","stdev","POF","stdev ","P10","stdev",
                  # "P50","stdev","P100","stdev","LTY","STY","VY")
    # out[,1]<-as.character(out[,1])
    # for(i in 2:ncol(out))out[,i]<-as.numeric(as.character(out[,i]))
    # out
  # })



# CheckConverg <- function(MSEobj, thresh=2, Plot=TRUE) {
  # nm<-MSEobj@nMPs
  # nsim<-MSEobj@nsim
  # proyears<-MSEobj@proyears
  
  # Yd <- CumlYd <- array(NA,c(nm,nsim))
  # P10 <- CumlP10 <- array(NA,c(nm,nsim))
  # P50 <- CumlP50 <- array(NA,c(nm,nsim))
  # P100 <- CumlP100 <- array(NA,c(nm,nsim))
  # POF <- CumlPOF <- array(NA,c(nm,nsim))
  # yind <-max(MSEobj@proyears-4,1):MSEobj@proyears
  # RefYd<-MSEobj@OM$RefY
  
  # for(m in 1:nm){
    # Yd[m,] <-round(apply(MSEobj@C[,m,yind],1,mean,na.rm=T)/RefYd*100,1)
    # POF[m,] <-round(apply(MSEobj@F_FMSY[,m,]>1,1,sum,na.rm=T)/proyears*100,1)
    # P10[m,] <-round(apply(MSEobj@B_BMSY[,m,]<0.1,1,sum,na.rm=T)/proyears*100,1)
    # P50[m,] <-round(apply(MSEobj@B_BMSY[,m,]<0.5,1,sum,na.rm=T)/proyears*100,1)
    # P100[m,] <-round(apply(MSEobj@B_BMSY[,m,]<1,1,sum,na.rm=T)/proyears*100,1)
    # CumlYd[m,] <- cumsum(Yd[m,]) / seq_along(Yd[m,])#/ mean(Yd[m,], na.rm=TRUE) 
    # CumlPOF[m,] <- cumsum(POF[m,]) / seq_along(POF[m,])# / mean(POF[m,], na.rm=TRUE)
    # CumlP10[m,] <- cumsum(P10[m,]) / seq_along(P10[m,])# / mean(P10[m,], na.rm=TRUE)
    # CumlP50[m,] <- cumsum(P50[m,]) / seq_along(P50[m,])# / mean(P50[m,], na.rm=TRUE)
    # CumlP100[m,] <- cumsum(P100[m,]) / seq_along(P100[m,])# / mean(P100[m,], na.rm=TRUE)
  # }
  
  # # CumlYd[is.nan(CumlYd)] <- 1
  # # CumlPOF[is.nan(CumlPOF)] <- 1
  # # CumlP10[is.nan(CumlP10)] <- 1
  # # CumlP50[is.nan(CumlP50)] <- 1
  # # CumlP100[is.nan(CumlP100)] <- 1
  # if (Plot) {
    # par(mfrow=c(2,3), cex.axis=1.5, cex.lab=1.7, oma=c(1,1,0,0), mar=c(5,5,1,1), bty="l")
    # matplot(t(CumlYd), type="l", xlab="Iteration", ylab="Rel. Yield")
    # matplot(t(CumlPOF), type="l", xlab="Iteration", ylab="Prob. F/FMSY > 1")
    # matplot(t(CumlP10), type="l", xlab="Iteration", ylab="Prob. B/BMSY < 0.1")
    # matplot(t(CumlP50), type="l", xlab="Iteration", ylab="Prob. B/BMSY < 0.5")
    # matplot(t(CumlP100), type="l", xlab="Iteration", ylab="Prob. B/BMSY < 1")
  # }
  
  # Chk <- function(X) {
    # # checks if difference in
    # # last 10 iterations is greater than thresh
    # L <- length(X)
    # Y <- 1: min(nsim, 10)
    
    # # return(mean(abs((X[L-1:10] - X[L]))/X[L], na.rm=TRUE) > thresh)
    # return(mean(abs((X[L-Y] - X[L])), na.rm=TRUE) > thresh)
  # }
  
  # NonCon <- sort(unique(c(which(apply(CumlYd, 1, Chk)), which(apply(CumlPOF, 1, Chk)),
                          # which(apply(CumlP10, 1, Chk)), which(apply(CumlP50, 1, Chk)), 	
                          # which(apply(CumlP100, 1, Chk)))))
  
  # if (length(NonCon) == 1) NonCon <- rep(NonCon, 2)
  # if (length(NonCon) > 0) {
    # if (Plot) {
      # plot(c(0,1), c(0,1), type="n", axes=FALSE, xlab="", ylab="")
      # text(0.5, 0.5, "Some MPs have not converged", cex=1)
      # #ev.new()
      # par(mfrow=c(2,3), cex.axis=1.5, cex.lab=1.7, oma=c(1,1,0,0), mar=c(5,5,1,1), bty="l")
      # matplot(t(CumlYd[NonCon,]), type="b", xlab="Iteration", ylab="Rel. Yield", lwd=2)
      # matplot(t(CumlPOF[NonCon,]), type="b", xlab="Iteration", ylab="Prob. F/FMSY > 1", lwd=2)
      # matplot(t(CumlP10[NonCon,]), type="b", xlab="Iteration", ylab="Prob. B/BMSY < 0.1", lwd=2)
      # matplot(t(CumlP50[NonCon,]), type="b", xlab="Iteration", ylab="Prob. B/BMSY < 0.5", lwd=2)
      # matplot(t(CumlP100[NonCon,]), type="b", xlab="Iteration", ylab="Prob. B/BMSY < 1", lwd=2)
      # legend(nsim*1.25, 50, legend=MSEobj@MPs[NonCon], col=1:length(NonCon), bty="n", 
             # xpd=NA, lty=1, lwd=2, cex=1.25)
    # }  
    
    # message("Some MPs may not have converged in ", nsim, " iterations (threshold = ", 
            # thresh, "%)")
    # message("MPs are: ", paste(MSEobj@MPs[NonCon], " "))
    # message("MPs #: ", paste(NonCon, " "))
    # return(data.frame(Num=NonCon, MP=MSEobj@MPs[NonCon]))
  # }
  # if (length(NonCon) == 0) {
    # if (Plot) {
      # plot(c(0,1), c(0,1), type="n", axes=FALSE, xlab="", ylab="")
      # text(0.5, 0.5, "All MPs converged", cex=1)
    # }	
    # message("All MPs appear to have converged in ", nsim, " iterations (threshold = ", 
            # thresh, "%)")
  # }

# }

# TwoSidedFun <- function(L1, s1, s2, Lens) {
  # Sl <- rep(0, length(Lens))
  # Sl[Lens < L1] <- exp( -((Lens[Lens < L1] - L1)^2)/(2*s1^2))
  # Sl[Lens >=L1 ] <- exp( -((Lens[Lens >= L1] - L1)^2)/(2*s2^2))
  # return(Sl)
# }

# getroot<-function(X,ageM,age95)uniroot(getSlopeFun, interval=c(0.0001, 5), age50=ageM[X], age95=age95[X])$root
# getSlope1 <- function(tst, L1, L0.05 ) (0.05 - TwoSidedFun(L1=L1, s1=tst, s2=1000, L0.05 ))^2
# getSlope2 <- function(tst, L1, s1, Linf, MaxSel) (MaxSel - TwoSidedFun(L1=L1, s1=s1, s2=tst, Linf))^2
# getSlopeFun <- function(SD, age50, age95) 0.95 - (1/(1+exp((age50-age95)/(age50*SD))))

# SelectFun <- function(i, SL0.05, SL1, MaxSel, Linfs, Lens) {
  # s1 <- optimise(getSlope1, interval=c(0, 100), L1=SL1[i], L0.05=SL0.05[i])$minimum
  # s2 <- optimise(getSlope2, interval=c(0, 1000), L1=SL1[i], s1=s1, Linf=Linfs[i], MaxSel=MaxSel[i])$minimum 
  # TwoSidedFun(L1=SL1[i], s1=s1, s2=s2, Lens=Lens[i,])
# }

# SelectFunGTG <- function(i, SL0.05, SL1, MaxSel, Linfs, LenGTG) {
  # s1 <- optimise(getSlope1, interval=c(0, 100), L1=SL1[i], L0.05=SL0.05[i])$minimum
  # s2 <- optimise(getSlope2, interval=c(0, 1000), L1=SL1[i], s1=s1, Linf=Linfs[i], MaxSel=MaxSel[i])$minimum 
  # NGTG <- dim(LenGTG)[1]
  # t(sapply(1:NGTG, function (X) TwoSidedFun(L1=SL1[i], s1=s1, s2=s2, Lens=LenGTG[X,i,])))
# }

# FitSelect <- function(Pars, V, Linf, Lens) {
  # SL0.05 <- (Pars[1])
  # SL1 <- (Pars[2])
  # MaxSel <- (Pars[3])
  # Lens <- t(as.matrix(Lens))
  # SS <- sum((V-SelectFun(1, SL0.05, SL1, MaxSel, Linf, Lens))^2)
  # return(SS)
# }


# range01 <- function(x) (x-min(x))/(max(x)-min(x)) # function to standardize to minimum and maximum

# getEffhist <- function(Esd, nyears, EffYears, EffLower, EffUpper) { # Simulate historical fishing effort trajectory
  # if(length(EffLower) == length(EffUpper) & length(EffUpper) == length(EffYears)) {
    # nsim <- length(Esd) # get nsim 
    # refYear <- floor(range01(EffYears+0.5) * nyears) + 1 # standardize years 
    # refYear[length(refYear)] <- nyears # first year is year 1 
	# Effs <- mapply(runif, n=nsim, min=EffLower, max=EffUpper) # sample Effort
    # if (nsim > 1) {
	  # Effs <- Effs / apply(Effs, 1, max) # standardize effort
	  # effort <- t(sapply(1:nsim, function (x) approx(x=refYear, y = Effs[x,],  method = "linear", n=nyears)$y)) # linear interpolation
	# }  
	# if (nsim == 1) {
	  # Effs <- Effs/max(Effs)
	  # effort <- approx(x=refYear, y = Effs,  method = "linear", n=nyears)$y
 	# }
    
    # Emu <- -0.5*Esd^2
    # Eerr <-array(exp(rnorm(nyears*nsim,rep(Emu,nyears),rep(Esd,nyears))),c(nsim,nyears)) # calc error
    # out <- NULL
    # eff <-effort * Eerr # add error 
	# out[[1]] <- eff
	# out[[2]] <- (effort[,nyears]-effort[,nyears-4])/5
    # return(out)
  # } else {
    # message("Input vectors of effort years and bounds not of same length")
    # return(NULL)
  # }
# }

# # helper function for above
# identifyPch <- function(x, y = NULL, n = length(x), pch = 19, ...)
# {
    # xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
    # sel <- rep(FALSE, length(x)); res <- integer(0)
    # while(sum(sel) < n) {
        # ans <- identify(x[!sel], y[!sel], n = 1, plot = FALSE, ...)
        # if(!length(ans)) break
        # ans <- which(!sel)[ans]
        # points(x[ans], y[ans], pch = pch)
        # sel[ans] <- TRUE
		# res <- c(res, ans)
    # }
    # out <- cbind(x[res], y[res])
	# out <- out[order(out[,1]),]
	# return(out)
# }







# ############################
# # These require help files #
# ############################

# # Generic Sketch function to allow users to pick points off a plotted grid
# SketchFun <- function(nyears=NULL, Years=NULL) {
  
  # if (length(Years) == 0 & length(nyears) == 0) stop()
  # if (length(Years) > 0) { 
    # nyears <- length(Years)
	# years <- Years
  # }	
  # if (length(Years) == 0) years <- 1:nyears
  # par(mfrow=c(1,1), mar=c(5, 4, 5, 2))
  # ys <- seq(from=0, to=1, by=0.05)
  # years1 <- seq(from=years[1], by=2, to=max(years))
  # years2 <- seq(from=years[2], by=2, to=max(years))
  # grd1 <- expand.grid(years1, ys)
  # grd2 <- expand.grid(years2, ys)
  # grd3 <- expand.grid(years, ys)
  # Xs <- grd3[,1]
  # Ys <- grd3[,2]

  # plot(grd1[,1], grd1[,2], col="black", pch=19, cex=0.2, xlab="Years", ylab="Variable", cex.lab=1.5)
  # points(grd2[,1], grd2[,2], col="darkgrey", pch=19, cex=0.2)
  
  # mtext(side=3, "Right Click to Finish. Escape to Quit.", xpd=NA, cex=1.25)
  # line1 <- "Use mouse to select points on the grid"
  # line2 <- "First and last year must be selected."
  # line3 <- "Select two points in a single year to represent range of uncertainty"
  # par(xpd=TRUE) 
  # text(years[1],par("usr")[4]+0.15*(par("usr")[4]-par("usr")[3]), line1,cex=1, pos=4) 
  # text(years[1],par("usr")[4]+0.1125*(par("usr")[4]-par("usr")[3]), line2,cex=1, pos=4) 
  # text(years[1],par("usr")[4]+0.075*(par("usr")[4]-par("usr")[3]),line3,cex=1, pos=4)  
  # par(xpd=FALSE)
  # message(line1,"\n", line2, "\n", line3, "\n")
  # flush.console()
  
  # par()
  # out <- NULL
  # out <- identifyPch(x=Xs, y=Ys, tolerance=0.1)
  # while(is.null(dim(out))) {
    # message("Must choose more than one point")
	# flush.console()
    # out <- identifyPch(x=Xs, y=Ys, tolerance=0.1)
  # }
  # while(min(out[,1]) != years[1]) {
    # message("Choose point(s) for first year (usually 0)")
	# flush.console()
	# dat <- rbind(out, identifyPch(x=Xs, y=Ys, tolerance=0.1))
	# out <- dat[order(dat[,1]),]
  # }
  # while(max(out[,1]) != years[length(years)]) {
    # message("Choose point(s) for last year (nyear)")
	# flush.console()
	# dat <- rbind(out, identifyPch(x=Xs, y=Ys, tolerance=0.1))
	# out <- dat[order(dat[,1]),]
  # }
  # ord <- order(out[,1], out[,2])
  # out <- out[ord,]
  # yrs <- unique(out[,1])
  # mat <- matrix(NA, nrow=length(yrs), ncol=3)
  # mat[,1] <- yrs
  # ind <- which(!duplicated(out[,1]))
  # mat[,2:3] <- out[ind,2]
  # for (X in seq_along(yrs)) {
    # chk <- out[,1] %in% yrs[X]
	# ind <- range(which(chk))
    # if(sum(chk) > 1) {
	  # mat[X,2:3] <- out[ind,2] 
	# }
  # }

  # lines(mat[,1], mat[,2])
  # lines(mat[,1], mat[,3])

  # colnames(mat) <- c("Years", "Lower", "Upper")
  # return(mat)
# }

# # Function for help users specify historical effort
# # Takes Fleet object, runs Sketch function for user to specify historical effort
# # Then returns Fleet object with Effort objects populated with output of Sketch
# ChooseEffort <- function(FleetObj, Years=NULL) {
  # nyears <- FleetObj@nyears
  # runSketch <- SketchFun(nyears, Years)
  # FleetObj@EffYears <- runSketch[,1]
  # FleetObj@EffLower <- runSketch[,2]
  # FleetObj@EffUpper <- runSketch[,3]
  # return(FleetObj)
# }

# # Subset MSE 
# # Function to subset MSE object by MPs or simulations
# # returns a subsetted MSE object 
# # MPs can be either numeric, or MP name 
# # Simulations can be numeric, or logical 
# Sub <- function(MSEobj, MPs=NULL, sims=NULL) {
  # Class <- class(MPs)
  # if(Class == "NULL") subMPs <- MSEobj@MPs
  # if(Class == "integer" | Class == "numeric") subMPs <- MSEobj@MPs[as.integer(MPs)]
  # if(Class == "character") subMPs <- MPs
  # SubMPs <- which(MSEobj@MPs %in% subMPs)
  # not <- (subMPs %in% MSEobj@MPs) # Check for MPs misspelled
  # ind <- which(not == FALSE)
  # newMPs <- MSEobj@MPs[SubMPs]
  # if (length(SubMPs) < 1) stop("MPs not found")
  # if (length(ind > 0)) {
    # message("Warning: MPs not found - ", paste0(subMPs[ind], " "))
	# message("Subsetting by MPs: ", paste0(newMPs, " "))
  # }
  
  # ClassSims <- class(sims)
  # if (ClassSims == "NULL") SubIts <- 1:MSEobj@nsim
  # if (ClassSims == "integer" | ClassSims == "numeric") {
    # sims <- 1:min(MSEobj@nsim, max(sims))
	# SubIts <- as.integer(sims)
  # }	
  # if (ClassSims == "logical") SubIts <- which(sims)

  # SubF <- MSEobj@F_FMSY[SubIts,SubMPs,]
  # SubB <- MSEobj@B_BMSY[SubIts,SubMPs,]
  # OutOM <- MSEobj@OM[SubIts,]
  
  # SubResults <- new('MSE',Name=MSEobj@Name, nyears=MSEobj@nyears, proyears=MSEobj@proyears, nMPs=length(SubMPs),
	# MPs=newMPs, nsim=length(SubIts), OMtable=OutOM, Obs=MSEobj@Obs[SubIts,], B_BMSYa=SubB, F_FMSYa=SubF, Ba=MSEobj@B[SubIts,SubMPs,], 
	# FMa=MSEobj@FM[SubIts,SubMPs,], Ca=MSEobj@C[SubIts,SubMPs,], TACa=MSEobj@TAC[SubIts,SubMPs,], SSB_hist=MSEobj@SSB_hist[SubIts,,,],
	# CB_hist=MSEobj@CB_hist[SubIts,,,], FM_hist=MSEobj@FM_hist[SubIts,,,])
  
 # return(SubResults)
# }

# DOM <- function(MSEobj, MPtg=NA) {
  # if (any(is.na(MPtg))) MPtg <- MSEobj@MPs 
  # proyears<-MSEobj@proyears
  # nMP <- MSEobj@nMPs
  # ind <- which(MSEobj@MPs %in%  MPtg)
  # MPr <- which(!(MSEobj@MPs %in%  MPtg))
  # yind<-max(MSEobj@proyears-4,1):MSEobj@proyears
  # y1<-1:(MSEobj@proyears-1)
  # y2<-2:MSEobj@proyears
  # Mat <- matrix(0, nrow=length(MPtg), ncol=nMP)
  # rownames(Mat) <- MPtg
  # colnames(Mat) <- MSEobj@MPs
  # POF <- P100 <- YieldMat <- IAVmat <- Mat
  # for (X in 1:length(MPtg)) {
    # # Overfishing (F > FMSY)
    # ind1 <- as.matrix(expand.grid(1:nsim, ind[X], 1:proyears))
    # ind2 <- as.matrix(expand.grid(1:nsim,  1:nMP, 1:proyears))
    # t1 <- apply(array(MSEobj@F_FMSY[ind1] > 1, dim=c(nsim, 1, proyears)), c(1,2), sum, na.rm=TRUE) 
    # t2 <- apply(array(MSEobj@F_FMSY[ind2] > 1, dim=c(nsim, nMP, proyears)), c(1,2), sum, na.rm=TRUE)
    # POF[X,] <- round(apply(matrix(rep(t1, nMP), nrow=nsim) < t2, 2, sum) / nsim * 100, 0)
	# # B < BMSY
    # t1 <- apply(array(MSEobj@B_BMSY[ind1] < 1, dim=c(nsim, 1, proyears)), c(1,2), sum, na.rm=TRUE) 
    # t2 <- apply(array(MSEobj@B_BMSY[ind2] < 1, dim=c(nsim, nMP, proyears)), c(1,2), sum, na.rm=TRUE) 
    # P100[X,] <- round(apply(matrix(rep(t1, nMP), nrow=nsim) < t2, 2, sum, na.rm=TRUE) / nsim * 100, 0)
    # # Relative yield in last 5 years
    # ind1 <- as.matrix(expand.grid(1:nsim, ind[X], yind))
    # ind2 <- as.matrix(expand.grid(1:nsim,  1:nMP, yind))
    # t1 <- apply(array(MSEobj@C[ind1], dim=c(nsim, 1, length(yind))), c(1,2), sum, na.rm=TRUE)
    # t2 <- apply(array(MSEobj@C[ind2], dim=c(nsim, nMP, length(yind))), c(1,2), sum, na.rm=TRUE)
    # YieldMat[X,] <- round(apply(matrix(rep(t1, nMP), nrow=nsim) > t2, 2, sum, na.rm=TRUE) / nsim * 100, 0)
    # # interannual variation in catch 
    # ind1 <- as.matrix(expand.grid(1:nsim, ind[X], y1))
    # ind2 <- as.matrix(expand.grid(1:nsim, ind[X], y2)) 
    # AAVY1 <- apply(array(((MSEobj@C[ind1]-MSEobj@C[ind2])^2)^0.5, dim=c(nsim, 1, length(y1))),1,mean,na.rm=T)/apply(array(MSEobj@C[ind2], dim=c(nsim, 1, length(y1))),1,mean,na.rm=T) 
    # ind1 <- as.matrix(expand.grid(1:nsim, 1:nMP, y1))
    # ind2 <- as.matrix(expand.grid(1:nsim, 1:nMP, y2)) 
    # AAVY2 <- apply(array(((MSEobj@C[ind1]-MSEobj@C[ind2])^2)^0.5, dim=c(nsim, nMP, length(y1))),c(1,2),mean,na.rm=T)/apply(array(MSEobj@C[ind2], dim=c(nsim, nMP, length(y1))),c(1,2),mean,na.rm=T)
    # IAVmat[X,] <- round(apply(matrix(rep(AAVY1, nMP), nrow=nsim) < AAVY2, 2, sum, na.rm=TRUE) / nsim * 100, 0)
  # }
  # out <- list() 
  # out$POF <- POF
  # out$P100 <- P100
  # out$Yd <- YieldMat
  # out$AAVY <- IAVmat
  # return(out)
# }




