

makeData <- function(Biomass, CBret, Cret, N, SSB, VBiomass, StockPars, 
                     FleetPars, ObsPars, ImpPars, RefPoints,
                     ErrList, OM, SampCpars, initD, silent=FALSE) {
  
  if(!silent) message("Simulating observed data")
  
  Name <- OM@Name
  nyears <- OM@nyears
  proyears <- OM@proyears
  nsim <- OM@nsim 
  nareas <- StockPars$nareas
  reps <- OM@reps
  
  Data <- new("Data", stock = "MSE")  # create a blank DLM data object
  if (reps == 1) Data <- OneRep(Data)  # make stochastic variables certain for only one rep
  Data <- replic8(Data, nsim)  # make nsim sized slots in the DLM data object
  
  Data@Name <- Name
  Data@Year <- 1:nyears
  
  # --- Observed catch ----
  # Simulated observed retained catch (biomass)
  Cobs <- ErrList$Cbiasa[, 1:nyears] * ErrList$Cerr[, 1:nyears] * apply(CBret, c(1, 3), sum)  
  Data@Cat <- Cobs 
  
  # --- Index of total abundance ----
  # Index of abundance from total biomass - beginning of year before fishing
  # apply hyperstability / hyperdepletion
  II <- (apply(Biomass, c(1, 3), sum)^ObsPars$betas) * ErrList$Ierr[, 1:nyears]  
  II <- II/apply(II, 1, mean)  # normalize
  Data@Ind <- II # index of total abundance
  
  # --- Real Indices ----
  if (!is.null(SampCpars$Data)) {
    # real data has been provided 
    types <- SampCpars$Data@Type
    chk <- types %in% c("Biomass", "VBiomass", "SpBiomass")
    if (any(!chk)) 
      if (!is.na(types)) 
        stop("Invalid index type in cpars$Data@Type. \nValid types are: 'Biomass', 'VBiomass', 'SpBiomass'", call.=FALSE)
    if (all(is.na(types))) {
      if (!silent) message("Data object provided in cpars but Data@Type is empty. Ignoring Data")
    } else {
      indices <- SampCpars$Data@RInd[1,,]
      if (all(is.na(indices))) {
        if (!silent) message("Data object provided in cpars but Data@RInd is empty. Ignoring Data")
      } else {
        if (!silent) message("Conditioning simulated RInd on cpars$Data@RInd")
        # add Real indices to Data object 
        Data@Type <- SampCpars$Data@Type
        RInd <- SampCpars$Data@RInd[1,,]
        temp <- array(RInd, dim=c(nrow(RInd), ncol(RInd), nsim))
        Data@RInd <- aperm(temp, c(3,1,2))
        
        # Calculate implied observation error
        sim.indices <- list(apply(Biomass, c(1, 3), sum),
                            apply(VBiomass, c(1, 3), sum),
                            apply(SSB, c(1, 3), sum))
        
        ind.type <- c('Biomass', 'VBiomass', 'SpBiomass')
        out.list <- list()
        for (type in Data@Type) {
          out.list[[type]] <- sapply(1:nsim, indfitwrap, type=type, 
                                  sim.indices=sim.indices, Data=Data)
        }   
      }
    }
  }
  SpBeta <- out.list$SpBiomass[1,] %>% unlist()
  SpAC <- out.list$SpBiomass[2,] %>% unlist()
  SpSD <- out.list$SpBiomass[3,] %>% unlist()
  

  sim <- 2
  Res <- rlnorm(nyears, 0, SpSD[sim])
  for(y in 2:nyears)Res[y]<-SpAC[sim]*Res[y-1]+Res[y]*(1-SpAC[sim]*SpAC[sim])^0.5
  index <- sim.indices[[3]]
  stindex <- index[sim,]/mean(index[sim,])
  index2 <- (stindex^SpBeta[sim]) * Res
  index2 <- (index2/mean(index2)) * mean(index[sim,])
  
  plot(index2, ylim=c(0, max(c(index2,index[sim,] ))), type="l", col="green")
  lines(Data@RInd[sim, 3, ], col="blue")
  
  lines(index[sim,], col="red")
  
  
  ## UP TO HERE #####
  
  ErrList$Ierr <- array(rlnorm((nyears + proyears) * nsim, 
                               mconv(1, rep(Isd, nyears + proyears)), 
                               sdconv(1, rep(Isd, nyears + proyears))), 
                        c(nsim, nyears + proyears))
  

  indfitwrap <- function(x, type, sim.indices, Data, plot=FALSE) {
    sim.index <- sim.indices[[match(type, ind.type)]][x,]
    obs.ind <- Data@RInd[x,match(type, Data@Type),]
    Year <- Data@Year
    indfit(sim.index,obs.ind, Year, plot, lcex=0.8)
  }
  
  indfit <- function(sim.index,obs.ind, Year, plot=FALSE, lcex=0.8){
    
    lcs<-function(x){
      x1<-x/mean(x) # rescale to mean 1
      x2<-log(x1)     # log it
      x3<-x2-mean(x2) # mean 0
      x3
    }
    getbeta<-function(beta,x,y)sum((y-x^beta)^2)
    
    sim.index <- lcs(sim.index) # log space conversion of standardized simulated index
    obs.ind <- lcs(obs.ind) # log space conversion of standardized observed ind
  
    if(plot){
      par(mfrow=c(1,2),mai=c(0.7,0.5,0.05,0.01),omi=c(0.01,0.2,0.01,0.01))
      plot(exp(sim.index),exp(obs.ind),xlab="",ylab="",pch=19,col=rgb(0,0,0,0.5))
      mtext("Model estimate",1,line=2.2)
      mtext("Index",2,outer=T,line=0)
    }
    
    opt<-optimize(getbeta,x=exp(sim.index),y=exp(obs.ind),interval=c(0.1,10))
    res<-exp(obs.ind)-(exp(sim.index)^opt$minimum)
    ac<-acf(res,plot=F)$acf[2,1,1] # lag-1 autocorrelation
    
    res2<-obs.ind-sim.index                  # linear, without hyperdepletion / hyperstability
    ac2<-acf(res2,plot=F)$acf[2,1,1] # linear AC
    
    if(plot){
      SSBseq<-seq(min(exp(sim.index)),max(exp(sim.index)),length.out=1000)
      lines(SSBseq,SSBseq^opt$minimum,col='#0000ff90',pch=19)
      legend('bottomright',legend=round(c(sum((obs.ind-sim.index)^2),opt$objective),3),text.col=c("black","blue"),bty='n',title="SSQ",cex=lcex)
      legend('topleft',legend=round(opt$minimum,3),text.col="blue",bty='n',title='Hyper-stability, beta',cex=lcex)
      legend('left',legend=round(cor(sim.index,obs.ind),3),bty='n',title='Correlation',cex=lcex)
      
      plot(Year,sim.index,ylab="",xlab="",ylim=range(c(obs.ind,sim.index)),type="l")
      mtext("Year",1,line=2.2)
      points(Year,obs.ind,col='#ff000090',pch=19)
      legend('topleft',legend=round(ac,3),text.col="red",bty='n',title="Lag 1 autocorrelation",cex=lcex)
      legend('bottomleft',legend=round(sd(res),3),text.col="red",bty='n',title="Residual StDev",cex=lcex)
      legend('topright',legend=c("Model estimate","Index"),text.col=c("black","red"),bty='n',cex=lcex)
    }

    data.frame(beta=opt$minimum,AC=ac,sd=sd(exp(obs.ind)/(exp(sim.index)^opt$minimum)),
               cor=cor(sim.index,obs.ind),AC2=ac2,sd2=sd(obs.ind-sim.index))
    
    # list(stats=data.frame(beta=opt$minimum,AC=ac,sd=sd(exp(obs.ind)/(exp(sim.index)^opt$minimum)),
    #                       cor=cor(sim.index,obs.ind),AC2=ac2,sd2=sd(obs.ind-sim.index)),
    #      mult=exp(obs.ind)/(exp(sim.index)^opt$minimum))
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # --- Index of recruitment ----
  Data@Rec <- apply(N[, 1, , ], c(1, 2), sum) * ErrList$Recerr[, 1:nyears] 
  Data@t <- rep(nyears, nsim) # number of years of data
  
  # --- Average catch ----
  Data@AvC <- apply(Cobs, 1, mean) # average catch over all years
  
  # --- Depletion ----
  # observed depletion
  Depletion <- apply(SSB[,,nyears,],1,sum)/RefPoints$SSB0 # current depletion
  Data@Dt <- ObsPars$Dbias * Depletion * 
    rlnorm(nsim, mconv(1, ObsPars$Derr), sdconv(1, ObsPars$Derr))
  
  Data@Dep <- ObsPars$Dbias * Depletion * 
    rlnorm(nsim, mconv(1, ObsPars$Derr), sdconv(1, ObsPars$Derr))  
  
  # --- Life-history parameters ----
  Data@vbLinf <- StockPars$Linfarray[,nyears] * ObsPars$Linfbias # observed vB Linf
  Data@vbK <- StockPars$Karray[,nyears] * ObsPars$Kbias # observed vB K
  Data@vbt0 <- StockPars$t0array[,nyears] * ObsPars$t0bias # observed vB t0
  Data@Mort <- StockPars$Marray[,nyears] * ObsPars$Mbias # natural mortality
  Data@L50 <- StockPars$L50array[,nyears] * ObsPars$lenMbias # observed length at 50% maturity
  Data@L95 <- StockPars$L95array[,nyears] * ObsPars$lenMbias # observed length at 95% maturity
  Data@L95[Data@L95 > 0.9 * Data@vbLinf] <- 0.9 * Data@vbLinf[Data@L95 > 0.9 * Data@vbLinf]  # Set a hard limit on ratio of L95 to Linf
  Data@L50[Data@L50 > 0.9 * Data@L95] <- 0.9 * Data@L95[Data@L50 > 0.9 * Data@L95]  # Set a hard limit on ratio of L95 to Linf
  Data@LenCV <- StockPars$LenCV # variablity in length-at-age - no error at this time
  Data@sigmaR <-  StockPars$procsd # observed sigmaR - assumed no obs error
  Data@MaxAge <- StockPars$maxage # maximum age - no error - used for setting up matrices only
  
  # Observed steepness values 
  hs <- StockPars$hs
  if (!is.null(OM@cpars[['hsim']])) {
    hsim <- SampCpars$hsim
    hbias <- hsim/hs  # back calculate the simulated bias
    if (OM@hbiascv == 0) hbias <- rep(1, nsim) 
    ObsPars$hbias <- hbias 
  } else {
    hsim <- rep(NA, nsim)  
    cond <- hs > 0.6
    hsim[cond] <- 0.2 + rbeta(sum(hs > 0.6), alphaconv((hs[cond] - 0.2)/0.8, (1 - (hs[cond] - 0.2)/0.8) * OM@hbiascv), 
                              betaconv((hs[cond] - 0.2)/0.8,  (1 - (hs[cond] - 0.2)/0.8) * OM@hbiascv)) * 0.8
    hsim[!cond] <- 0.2 + rbeta(sum(hs <= 0.6), alphaconv((hs[!cond] - 0.2)/0.8,  (hs[!cond] - 0.2)/0.8 * OM@hbiascv), 
                               betaconv((hs[!cond] - 0.2)/0.8, (hs[!cond] - 0.2)/0.8 * OM@hbiascv)) * 0.8
    hbias <- hsim/hs  # back calculate the simulated bias
    if (OM@hbiascv == 0) hbias <- rep(1, nsim) 
    ObsPars$hbias <- hbias
  }
  Data@steep <- hs * ObsPars$hbias # observed steepness
  
  # --- Reference points ----
  # Simulate observation error in BMSY/B0 
  ntest <- 20  # number of trials  
  BMSY_B0bias <- array(rlnorm(nsim * ntest, 
                              mconv(1, OM@BMSY_B0biascv), sdconv(1, OM@BMSY_B0biascv)), 
                       dim = c(nsim, ntest))  # trial samples of BMSY relative to unfished  
  test <- array(RefPoints$SSBMSY_SSB0 * BMSY_B0bias, dim = c(nsim, ntest))  # the simulated observed BMSY_B0 
  indy <- array(rep(1:ntest, each = nsim), c(nsim, ntest))  # index
  indy[test > max(0.9, max(RefPoints$SSBMSY_SSB0))] <- NA  # interval censor
  BMSY_B0bias <- BMSY_B0bias[cbind(1:nsim, apply(indy, 1, min, na.rm = T))]  # sample such that BMSY_B0<90%
  ObsPars$BMSY_B0bias <- BMSY_B0bias
  
  Data@FMSY_M <- RefPoints$FMSY_M * ObsPars$FMSY_Mbias # observed FMSY/M
  Data@BMSY_B0 <- RefPoints$SSBMSY_SSB0 * ObsPars$BMSY_B0bias # observed BMSY/B0
  Data@Cref <- RefPoints$MSY * ObsPars$Crefbias # Catch reference - MSY with error
  Data@Bref <- RefPoints$VBMSY * ObsPars$Brefbias # Vuln biomass ref - VBMSY with error
  
  # Generate values for reference SBMSY/SB0
  # should be calculated from unfished - won't be correct if initD is set
  I3 <- apply(Biomass, c(1, 3), sum)^ObsPars$betas  # apply hyperstability / hyperdepletion
  I3 <- I3/apply(I3, 1, mean)  # normalize index to mean 1
  if (!is.null(initD)) {
    b1 <- apply(Biomass, c(1, 3), sum)
    b2 <- matrix(RefPoints$BMSY, nrow=nsim, ncol=nyears)
    ind <- apply(abs(b1/ b2 - 1), 1, which.min) # find years closest to BMSY
    Iref <- diag(I3[1:nsim,ind])  # return the real target abundance index closest to BMSY
  } else {
    Iref <- apply(I3[, 1:5], 1, mean) * RefPoints$SSBMSY_SSB0  # return the real target abundance index corresponding to BMSY
  }
  Data@Iref <- Iref * ObsPars$Irefbias # index reference with error
  
  # --- Abundance ----
  # Calculate vulnerable and spawning biomass abundance --
  M_array <- array(0.5*StockPars$M_ageArray[,,nyears], dim=c(nsim, StockPars$maxage, nareas))
  A <- apply(VBiomass[, , nyears, ] * exp(-M_array), 1, sum) # Abundance (mid-year before fishing)
  Asp <- apply(SSB[, , nyears, ] * exp(-M_array), 1, sum)  # Spawning abundance (mid-year before fishing)
  OFLreal <- A * (1-exp(-RefPoints$FMSY))  # the true simulated Over Fishing Limit
  
  Data@Abun <- A * ObsPars$Abias * 
    rlnorm(nsim, mconv(1, ObsPars$Aerr), sdconv(1, ObsPars$Aerr)) # observed vulnerable abundance
  Data@SpAbun <- Asp * ObsPars$Abias * 
    rlnorm(nsim, mconv(1, ObsPars$Aerr), sdconv(1, ObsPars$Aerr)) # spawing abundance
  
  # --- Catch-at-age ----
  Data@CAA <- simCAA(nsim, nyears, StockPars$maxage, Cret, ObsPars$CAA_ESS, ObsPars$CAA_nsamp) 
  
  # --- Catch-at-length ----
  vn <- apply(N, c(1,2,3), sum) * FleetPars$retA[,,1:nyears] # numbers at age in population that would be retained
  vn <- aperm(vn, c(1,3, 2))
  
  CALdat <- simCAL(nsim, nyears, StockPars$maxage, ObsPars$CAL_ESS, 
                   ObsPars$CAL_nsamp, StockPars$nCALbins, StockPars$CAL_binsmid, 
                   vn, FleetPars$retL, StockPars$Linfarray, 
                   StockPars$Karray, StockPars$t0array, StockPars$LenCV)
  
  Data@CAL_bins <- StockPars$CAL_bins
  Data@CAL <- CALdat$CAL # observed catch-at-length
  Data@ML <- CALdat$ML # mean length
  Data@Lc <- CALdat$Lc # modal length 
  Data@Lbar <- CALdat$Lbar # mean length above Lc 
  
  Data@LFC <- CALdat$LFC * ObsPars$LFCbias # length at first capture
  Data@LFS <- FleetPars$LFS[nyears,] * ObsPars$LFSbias # length at full selection
  
  # --- Previous Management Recommendations ----
  Data@MPrec <- apply(CBret, c(1, 3), sum) # catch in last year
  Data@MPeff <- rep(1, nsim) # effort in last year = 1 
  
  # --- Store OM Parameters ----
  # put all the operating model parameters in one table
  ind <- which(lapply(StockPars, length) == nsim)
  stock <- as.data.frame(StockPars[ind])
  stock$Fdisc <- NULL
  ind <- which(lapply(FleetPars, length) == nsim)
  fleet <- as.data.frame(FleetPars[ind])
  
  ind <- which(lapply(ImpPars, length) == nsim)
  imp <- as.data.frame(ImpPars[ind])
  refs <- RefPoints %>% select('MSY', 'FMSY', 'SSBMSY_SSB0', 'BMSY_B0',
                               'UMSY', 'FMSY_M', 'RefY', 'Blow', 'MGT', 'SSB0')
  
  OMtable <- data.frame(stock, fleet, imp, refs, ageM=StockPars$ageM[,nyears], 
                     L5=FleetPars$L5[nyears, ], LFS=FleetPars$LFS[nyears, ], 
                     Vmaxlen=FleetPars$Vmaxlen[nyears, ],
                     LR5=FleetPars$LR5[nyears,], LFR=FleetPars$LFR[nyears,], 
                     Rmaxlen=FleetPars$Rmaxlen[nyears,], 
                     DR=FleetPars$DR[nyears,], OFLreal)
                 
  OMtable <- OMtable[,order(names(OMtable))]
  Data@OM <- OMtable
  
  # --- Store Obs Parameters ----
  ObsTable <- as.data.frame(ObsPars)
  ObsTable <- ObsTable[,order(names(ObsTable))]
  Data@Obs <- ObsTable # put all the observation error model parameters in one table
  
  # --- Misc ----
  Data@Units <- "unitless"
  Data@Ref_type <- "Simulated OFL"
  Data@wla <- rep(StockPars$a, nsim)
  Data@wlb <- rep(StockPars$b, nsim)
  Data@nareas <- nareas
  Data@Ref <- OFLreal 
  Data@LHYear <- nyears  # Last historical year is nyears (for fixed MPs)
  Data@Misc <- vector("list", nsim)
  
  Data
}


updateData <- function(Data, OM, MPCalcs, Effort, Biomass, Biomass_P, CB_Pret, 
                       N_P, SSB_P, VBiomass_P, RefPoints, ErrList, FMSY_P, retA_P, 
                       retL_P, StockPars, FleetPars, ObsPars, 
                       upyrs, interval, y=2, 
                       mm=1, Misc) {
  
  yind <- upyrs[match(y, upyrs) - 1]:(upyrs[match(y, upyrs)] - 1) # index
  
  nyears <- OM@nyears
  proyears <- OM@proyears
  nsim <- OM@nsim 
  nareas <- StockPars$nareas
  reps <- OM@reps
  
  Data@Year <- 1:(nyears + y - 1)
  Data@t <- rep(nyears + y, nsim)
  
  # --- Simulate catches ---- 
  CBtemp <- CB_Pret[, , yind, , drop=FALSE] # retained catch-at-age
  CNtemp <- retA_P[,,yind+nyears, drop=FALSE] * 
    apply(N_P[,,yind,, drop=FALSE], c(1,2,3), sum) # retained age structure
  CBtemp[is.na(CBtemp)] <- tiny
  CBtemp[!is.finite(CBtemp)] <- tiny
  CNtemp[is.na(CNtemp)] <- tiny
  CNtemp[!is.finite(CNtemp)] <- tiny
  CNtemp <- aperm(CNtemp, c(1,3,2))
  
  # --- Observed catch ----
  # Simulated observed retained catch (biomass)
  Cobs <- ErrList$Cbiasa[, nyears + yind] * ErrList$Cerr[, nyears + yind] * 
    apply(CBtemp, c(1, 3), sum, na.rm = TRUE)
  Data@Cat <- cbind(Data@Cat, Cobs) 
  
  # --- Index of total abundance ----
  I2 <- (cbind(apply(Biomass, c(1, 3), sum), 
               apply(Biomass_P, c(1, 3), sum)[, 1:(y - 1)])^ObsPars$betas) * 
    ErrList$Ierr[, 1:(nyears + (y - 1))]
  
  I2[is.na(I2)] <- tiny
  I2 <- I2/apply(I2, 1, mean)
  Data@Ind <- I2
  
  # --- Index of recruitment ----
  Recobs <- ErrList$Recerr[, nyears + yind] * apply(array(N_P[, 1, yind, ], 
                                                          c(nsim, interval[mm], nareas)),
                                                    c(1, 2), sum)
  Data@Rec <- cbind(Data@Rec, Recobs)
  
  # --- Average catch ----
  Data@AvC <- apply(Data@Cat, 1, mean)
  
  # --- Depletion ----
  Depletion <- apply(SSB_P[, , y, ], 1, sum)/RefPoints$SSB0 
  Depletion[Depletion < tiny] <- tiny
  Data@Dt <- ObsPars$Dbias * Depletion * rlnorm(nsim, mconv(1, ObsPars$Derr), sdconv(1, ObsPars$Derr))
  Data@Dep <- ObsPars$Dbias * Depletion * rlnorm(nsim, mconv(1, ObsPars$Derr), sdconv(1, ObsPars$Derr))
  
  # --- Update life-history parameter estimates for current year ----
  Data@vbLinf <- StockPars$Linfarray[,nyears+y] * ObsPars$Linfbias # observed vB Linf
  Data@vbK <- StockPars$Karray[,nyears+y] * ObsPars$Kbias # observed vB K
  Data@vbt0 <- StockPars$t0array[,nyears+y] * ObsPars$t0bias # observed vB t0
  Data@Mort <- StockPars$Marray[,nyears+y] * ObsPars$Mbias # natural mortality
  Data@L50 <- StockPars$L50array[,nyears+y] * ObsPars$lenMbias # observed length at 50% maturity
  Data@L95 <- StockPars$L95array[,nyears+y] * ObsPars$lenMbias # observed length at 95% maturity
  Data@L95[Data@L95 > 0.9 * Data@vbLinf] <- 0.9 * Data@vbLinf[Data@L95 > 0.9 * Data@vbLinf]  # Set a hard limit on ratio of L95 to Linf
  Data@L50[Data@L50 > 0.9 * Data@L95] <- 0.9 * Data@L95[Data@L50 > 0.9 * Data@L95]  # Set a hard limit on ratio of L95 to Linf
  
  
  # --- Abundance ----
  # Calculate vulnerable and spawning biomass abundance --
  M_array <- array(0.5*StockPars$M_ageArray[,,nyears+y], dim=c(nsim, StockPars$maxage, nareas))
  A <- apply(VBiomass_P[, , y, ] * exp(-M_array), 1, sum) # Abundance (mid-year before fishing)
  Asp <- apply(SSB_P[, , y, ] * exp(-M_array), 1, sum)  # Spawning abundance (mid-year before fishing)
  Data@Abun <- A * ObsPars$Abias * rlnorm(nsim, mconv(1, ObsPars$Aerr), sdconv(1, ObsPars$Aerr))
  Data@SpAbun <- Asp * ObsPars$Abias * rlnorm(nsim, mconv(1, ObsPars$Aerr), sdconv(1, ObsPars$Aerr))
  Data@Ref <- A * (1 - exp(-FMSY_P[,mm,y])) 
  
  # --- Catch-at-age ----
  # previous CAA
  oldCAA <- Data@CAA
  Data@CAA <- array(0, dim = c(nsim, nyears + y - 1, StockPars$maxage))
  Data@CAA[, 1:(nyears + y - interval[mm] - 1), ] <- oldCAA[, 1:(nyears + y - interval[mm] - 1), ] 
  # update CAA
  CAA <- simCAA(nsim, yrs=length(yind), StockPars$maxage, Cret=CNtemp, ObsPars$CAA_ESS, ObsPars$CAA_nsamp)
  Data@CAA[, nyears + yind, ] <- CAA
  
  

  # --- Catch-at-length ----
  oldCAL <- Data@CAL
  Data@CAL <- array(0, dim = c(nsim, nyears + y - 1, StockPars$nCALbins))
  Data@CAL[, 1:(nyears + y - interval[mm] - 1), ] <- oldCAL[, 1:(nyears + y - interval[mm] - 1), ]
  
  CAL <- array(NA, dim = c(nsim, interval[mm], StockPars$nCALbins))  
  vn <- (apply(N_P[,,,], c(1,2,3), sum) * retA_P[,,(nyears+1):(nyears+proyears)]) # numbers at age that would be retained
  vn <- aperm(vn, c(1,3,2))
  
  CALdat <- simCAL(nsim, nyears=length(yind), StockPars$maxage, ObsPars$CAL_ESS, 
                   ObsPars$CAL_nsamp, StockPars$nCALbins, StockPars$CAL_binsmid, 
                   vn=vn[,yind,, drop=FALSE], retL=retL_P[,,nyears+yind, drop=FALSE],
                   Linfarray=StockPars$Linfarray[,nyears + yind, drop=FALSE],  
                   Karray=StockPars$Karray[,nyears + yind, drop=FALSE], 
                   t0array=StockPars$t0array[,nyears + yind,drop=FALSE],
                   LenCV=StockPars$LenCV)

  Data@CAL[, nyears + yind, ] <- CALdat$CAL # observed catch-at-length
  Data@ML <- cbind(Data@ML, CALdat$ML) # mean length
  Data@Lc <- cbind(Data@Lc, CALdat$Lc) # modal length 
  Data@Lbar <- cbind(Data@Lbar, CALdat$Lbar) # mean length above Lc 
  
  Data@LFC <- CALdat$LFC * ObsPars$LFCbias # length at first capture
  Data@LFS <- FleetPars$LFS[nyears+y,] * ObsPars$LFSbias # length at full selection

  # --- Previous Management Recommendations ----
  Data@MPrec <- MPCalcs$TACrec # last MP  TAC recommendation
  Data@MPeff <- Effort[, mm, y-1] # last recommended effort
  
  Data@Misc <- Misc
  
  Data
}


