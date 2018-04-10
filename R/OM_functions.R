#' Forces correlation among operating model parameters for M, K, Linf and L50
#'
#' @description Uses typical correlations among estimated parameters to generate realistic samples for natural mortality rate (M), growth rate (K), maximum length (Linf) and length at 50% maturity (L50), these are placed in the cpars slot
#' @param OM An operating model object with M, growth, stock-recruitment and maturity parameters specified.
#' @param nsim The number of simulated values to create (note that OM@nsim will be used preferentially). 
#' @param plot Should the sampled parameters and distributions be plotted?
#' @return An object of class OM with a populated (or appended) cpars slot
#' @author T. Carruthers (Canadian DFO grant)
#' @importFrom mvtnorm rmvnorm
#' @export 
#' @examples
#' testOM<-ForceCor(testOM)
ForceCor<-function(OM,nsim=48,plot=T){
  
  if("nsim"%in%slotNames(OM))nsim<-OM@nsim
  if("seed"%in%slotNames(OM))set.seed(OM@seed)
  
  OM@nsim<-nsim
  
  colline=makeTransparent('blue',60)
  lwdline=4
  histcol='black'
  # Estimation cross correlation is typically -0.9 for Linf and K
  sigma<-matrix(c(1,-0.9,-0.9,1),nrow=2)
  
  # Other Parameters correlation  from meta analysis (technically incorrect in the sapce of other v
  #                  M     K    F_linf    Linf
  sigma=matrix(c(1.000,  0.721, -0.568, -0.721,
                 0.721,  1.000, -0.107, -0.910,
                 -0.568, -0.107, 1.000,  0.407,
                 -0.721, -0.910, 0.407,  1.000),nrow=4)
  
  
  means<-c(mean(OM@M),mean(OM@K),mean(OM@L50),mean(OM@Linf))
  
  sim<-as.data.frame(mvtnorm::rmvnorm(nsim,rep(0,4),sigma))
  sim<-sim-rep(apply(sim,2,mean),each=nsim)#mean 0
  sim<-sim/rep(apply(sim,2,sd),each=nsim)# cv=1
  cvs<-c(OM@M[2]-OM@M[1],OM@K[2]-OM@K[1],OM@L50[2]-OM@L50[1],OM@Linf[2]-OM@Linf[1])/(1.96*2)/means
  sim<-exp(sim*rep(cvs,each=nsim))
  sim<-sim/rep(apply(sim,2,mean),each=nsim)
  sim<-sim*rep(means,each=nsim)
  
  if(plot){ 
    par(mfrow=c(4,4),mai=c(0.3,0.3,0.4,0.05),omi=c(0.02,0.02,0.3,0.02))
    
    labs<-c("M","K","L50","Linf")
    bounds<-matrix(c(OM@M,OM@K,OM@L50,OM@Linf),nrow=2)
    
    for(i in 1:4){
      
      for(j in 1:4){
        
        if(i == j){
          
          if(i==1){
            
            hist(sim[,1],main="Natural mortality rate (M)",col=histcol,border='white',xlab="",axes=F)
            axis(1)
            abline(v=OM@M,col=colline,lwd=lwdline)
            
          }else if(i==2){
            
            hist(sim[,2],main="Growth rate (K)",col=histcol,border='white',xlab="",axes=F)
            axis(1)
            abline(v=OM@K,col=colline,lwd=lwdline)
            
            
          }else if(i==4){
            hist(sim[,4],main="Maximum length (Linf)",col=histcol,border='white',xlab="",axes=F)
            axis(1)
            abline(v=OM@Linf,col=colline,lwd=lwdline)

          }else{
            
            hist(sim[,3],main="Length at 50% maturity (L50)",col=histcol,border='white',xlab="",axes=F)
            axis(1)
            abline(v=OM@L50,col=colline,lwd=lwdline)
            
          }
          
        }else{ # not positive diagonal
          
          plot(sim[,j],sim[,i],axes=F,col="white")
          polygon(bounds[c(1,1,2,2),j],bounds[c(1,2,2,1),i],col=colline,border="white")
          points(sim[,j],sim[,i],pch=19)
          
          axis(1)
          axis(2)
          
        }
        
      }
    }
    
    mtext("Sampled parameters and cross-correlations",3,font=2,outer=T)
    
  } # if plot
  
  names(sim)<-c("M","K","L50","Linf")
  
  OM@cpars$M<-sim$M
  OM@cpars$K<-sim$K
  OM@cpars$L50<-sim$L50
  OM@cpars$Linf<-sim$Linf
  
  OM
  
}

#' Impute correlated life-history parameters
#' 
#' Impute correlated life-history parameters by imputing K and L50 values from M and Linf
#'
#' @param OM An object of class 'OM'
#' @param samp.m Numeric. A multiple of `OM@sim`. The number of initial samples to generate 
#' before filtering to ensure predicted values are within ranges specified in OM (e.g., `OM@M`)
#' @param plot Logical. Should the plot be produced?
#' @param ign.bounds Logical. Should the OM bounds be ignored if it is not possible to generate
#' enough predicted values within the specified bounds
#'
#' @return An OM with `OM@cpars` populated with `OM@nsim` samples of M, K, Linf and L50
#' @author A. Hordyk
#' @export
#' @examples 
#' testOM<-ImputeLH(testOM)
#' 
#' @importFrom Amelia amelia
#' @importFrom dplyr bind_rows %>%
#'
#'
ImputeLH <- function(OM, samp.m=100, plot=TRUE, ign.bounds=TRUE) {
  set.seed(OM@seed)
  DBdata <- DLMtool::LHdata %>% select(M, K, relL, MK)
  if (class(OM) != "OM") stop('Object must be class "OM"', call. = FALSE)
  
  inM <- inK <- inLinf <- inL50 <- NULL
  
  for (nm in c("M", "K", "Linf", "L50")) {
    if (length(OM@cpars[[nm]])>0) {
      tempval <- OM@cpars[[nm]]
    } else {
      tempval <- slot(OM, nm)
    }
    assign(paste0('in', nm), tempval)
  }
  
  parlist <- list(M=inM, K=inK, Linf=inLinf, L50=inL50)
  means <- lapply(parlist, mean)
  
  Ms <- Ks <- Linfs <- L50s <- NULL
  # Generate samples of known parameters - keep cpars
  nsamp <- OM@nsim * samp.m
  
  ncpars <- max(unlist(lapply(parlist, length)))
  if (ncpars < nsamp) {
    sampRows <- sample(1:ncpars, nsamp, replace=TRUE)  
  } else {
    sampRows <- sample(1:ncpars, nsamp, replace=FALSE) 
  }
  
  incpars <- data.frame(M=TRUE, K=TRUE, Linf=TRUE, L50=TRUE)
  for (nm in c("M", "K", "Linf", "L50")) {
    if (length(parlist[[nm]]) == 2) {
      varsd <- (max(parlist[[nm]]) - mean(parlist[[nm]]))/2
      varmean <- mean(parlist[[nm]])
      var <- rnorm(nsamp*100, varmean, varsd)
      var <- var[var>=parlist[[nm]][1] & var<=parlist[[nm]][2]][1:nsamp]
      assign(paste0(nm, "s"), var)
      incpars[nm] <- FALSE
    } else {
      var <- as.numeric(unlist((OM@cpars[nm])))[sampRows]
      assign(paste0(nm, "s"), var)
      incpars[nm] <- TRUE
    }
  }

  if (rowSums(incpars) == 4) stop("M, K, Linf, and L50 all in OM@cpars. Nothing to impute!", call.=FALSE)
  indata_all <- data.frame(M=Ms, K=Ks, L50=L50s, Linf=Linfs)
  if (any(indata_all$L50>indata_all$Linf)) stop("L50 is greater than Linf")
  indata <- data.frame(M=Ms, K=Ks)
  
  ind <- apply(indata, 2, mean) == 0
  indata[,ind] <- NA # drop any unknown parameters (e.g OM@Linf = c(0,0)= unknown)
 
  if (!all(is.na(indata$M)) & !incpars$K) {
    indata$K <- NA # if M is present, impute K 
    message("Imputing K from M")
  }
  if (!all(is.na(indata$K)) & !incpars$M) {
    indata$M <- NA # if K is present, impute M 
    message("Imputing M from K")
  }
  if (incpars$K & incpars$M & !(incpars$L50 | incpars$Linf)) {
    message('M and K in cpars. Imputing L50/Linf')
    indata$MK <- Ms/Ks
  }
  if ((incpars$L50 & incpars$Linf & !(incpars$K | incpars$M))) {
    message('L50 and Linf in cpars. Imputing M/K')
    indata$relL <- L50s/Linfs
  }
  
  # if (!all(is.na(indata$Linf)) & !incpars$L50) indata$L50 <- NA # if Linf is present, impute L50 
  # if (!all(is.na(indata$L50)) & !incpars$Linf) indata$Linf <- NA # if L50 is present, impute Linf
  
  ind <- colSums(apply(indata, 2, is.na)) == 0
  indata <- indata[,ind, drop=FALSE]
  if(length(indata$MK) > 0) message("M/K provided. Imputing L50/Linf")
  if (ncol(indata) < 1) stop("At least one parameter but be provided - M or K", call.=FALSE)
  
  # Transform to approximate multivariate normal
  tpow <- 1/3
  tdata <- as.data.frame(DBdata^tpow)
  # result <- MVN::mvn(data = tdata, mvnTest = "mardia", univariatePlot = "histogram")
  
  logs <- NULL
  bounds <- matrix(c(1, 0.025, 0.8, # bounds for M
                     2, 0.025, 1.5, # bounds for K 
                     3, 0.2, 0.8, #  bounds for L50/Linf
                     4, 0.2, 4), # bounds for M/K
                   nrow=ncol(DBdata), ncol=3, byrow=TRUE)
  if (!is.null(logs)) bounds[logs, 2:3] <- log(bounds[logs, 2:3] )
  bounds[,2:3] <- bounds[,2:3]^tpow # transform bounds
  
  tindata <- as.data.frame(indata^tpow)
  
  # impute missing values
  sink("temp")
  mod <- Amelia::amelia(alldata, bounds=bounds, logs=logs, max.resample=5000, verbose=FALSE) 
  sink()
  unlink("temp")
  
  tMeanEsts <- Reduce("+", mod$imputations) / length(mod$imputations)
  tMeanEsts <- tMeanEsts[(nrow(DBdata)+1):nrow(alldata),]
  MeanEsts <- as.data.frame(tMeanEsts^(1/tpow))
  
  
  
  if (is.null(indata$K)) MeanEsts$K <- MeanEsts$M / MeanEsts$MK
  if (is.null(indata$M)) MeanEsts$M <- MeanEsts$MK * MeanEsts$K 
  MeanEsts$Linf <- indata_all$Linf # MeanEsts$relL * MeanEsts$Linf
  MeanEsts$L50 <- indata_all$Linf * MeanEsts$relL 
  
  MeanEsts$relL <- NULL
  MeanEsts$MK <- NULL
  
  
  Out <- MeanEsts

  if (ign.bounds) {
    ok <- array(TRUE, dim=dim(MeanEsts[,1:4]))
    cols <- 1:4
    for (xx in cols) {
      rng <- range(parlist[[xx]])
      ok[,xx] <- MeanEsts[,xx] >= rng[1] & MeanEsts[,xx] <= rng[2]
    }
    chk <- apply(ok, 2, prod)
    oksamps <- chk == 1  
    if (sum(chk) !=4) {
      ind <- which(!oksamps)
      message('Some predicted samples are outside specified bounds in OM: ', paste0(names(parlist)[ind], ", "), 
              "\nIgnoring bounds for these parameters.")
    }
  }
  
  
  
  # filter to satisfy all bounds
  if (!ign.bounds){
    maxcount <- 50
    chk <- sum(apply(ok, 1, prod)==1)
    count <- 0
    while (chk < OM@nsim & count < maxcount) {
      count <- count + 1
      oksamps <- as.data.frame(t(apply(ok, 2, sum)))  
      ind <- which.min(oksamps)
      MeanEsts[,ind] <- indata_all[,ind]
      cols2 <- cols[!cols %in% ind]
      ok <- array(TRUE, dim=dim(MeanEsts[,1:4]))
      for (xx in cols2) {
        rng <- range(parlist[[xx]])
        if (mean(rng)!=0) {
          ok[,xx] <- MeanEsts[,xx] >= rng[1] & MeanEsts[,xx] <= rng[2]  
        }
        
      }
      warning('Could not generate sufficient predicted samples within specified bounds for ', names(parlist)[ind], 
              ".\nTry increase 'samp.m' and re-run. \nIgnoring predicted values for this parameter and using specified bounds")
      chk <- sum(apply(ok, 1, prod)==1)
    }
    Out <- MeanEsts[apply(ok, 1, prod)==1,]
  }
  
  
  Out <- Out[1:OM@nsim,]
  
  if(plot){ 
    par(mfrow=c(4,4),mai=c(0.3,0.3,0.4,0.05),omi=c(0.02,0.02,0.3,0.02))
    
    colline=makeTransparent('blue',60)
    lwdline=4
    histcol='black'
    labs<-c("M","K","Linf","L50")
    bounds<-matrix(unlist(lapply(parlist, range)),nrow=2)
    
    for(i in 1:4){
      rng <- range(parlist[[i]])
      rng2 <- range(Out[,i])
      rng2[1] <- min(c(min(rng), min(rng2)))
      rng2[2] <- max(c(max(rng), max(rng2)))
      for(j in 1:4){
        
        if(i == j){
          
          if(i==1){
            
            hist(Out[,1],main="Natural mortality rate (M)",col=histcol,border='white',xlab="",axes=F, xlim=rng2)
            axis(1)
            abline(v=rng,col=colline,lwd=lwdline)
            
          }else if(i==2){
            
            hist(Out[,2],main="Growth rate (K)",col=histcol,border='white',xlab="",axes=F, xlim=rng2)
            axis(1)
            abline(v=rng,col=colline,lwd=lwdline)
            
            
          }else if(i==3){
            
            hist(Out[,3],main="Asymptotic length (Linf)",col=histcol,border='white',xlab="",axes=F, xlim=rng2)
            axis(1)
            abline(v=rng,col=colline,lwd=lwdline)
            
            
          }else{
            
            hist(Out[,4],main="Length at 50% maturity (L50)",col=histcol,border='white',xlab="",axes=F, xlim=rng2)
            axis(1)
            abline(v=rng,col=colline,lwd=lwdline)
            
          }
          
        }else{ # not positive diagonal
          
          plot(Out[,j],Out[,i],axes=F,col="white")
          polygon(bounds[c(1,1,2,2),j],bounds[c(1,2,2,1),i],col=colline,border="white")
          points(Out[,j],Out[,i],pch=19)
          
          axis(1)
          axis(2)
          
        }
        
      }
    }
  }
  
  OM@cpars$M <- Out$M 
  OM@cpars$K <- Out$K
  OM@cpars$Linf <- Out$Linf
  OM@cpars$L50 <- Out$L50
  
  OM@M <- c(0,0)
  OM@K <- c(0,0)
  OM@Linf <- c(0,0)
  OM@L50 <- c(0,0)
  OM
  
}




#' Replace an existing Stock, Fleet, Obs, or Imp object 
#' 
#' A function that replaces a Stock, Fleet, Obs, or Imp object from an 
#' OM with one from another OM. Mainly used for internal functions.
#' 
#' @param OM An operating model object (class OM) which will be updated with a sub-model from another OM
#' @param from The OM object from which the sub-model is being taken
#' @param Sub A character string specifying what object type to replace
#' "Stock", "Fleet", "Obs" or "Imp" (default is all four which is probably not what you want to do)
#' @param Quiet Should the function not return a text message
#' @return An object of class OM
#' @author A. Hordyk
#' @examples 
#' \dontrun{
#' OM <- Replace(OM, fromOM, "Stock")
#' }
#' 
#' @export 
Replace <- function(OM, from, Sub=c("Stock", "Fleet", "Obs", "Imp"),Quiet=F) {
  if (class(OM) =="character") OM <- get(OM)
  if (class(from) !="OM") fromOM <- get(from)
  if (class(OM) !="OM") stop("OM must be of class OM ", call.=FALSE)
  if (class(from) !="OM") stop("''from' must be of class OM ", call.=FALSE)
  Sub <- match.arg(Sub, several.ok=TRUE)
  
  Stock <- SubOM(OM, "Stock")
  Fleet <- SubOM(OM, "Fleet")
  Obs <- SubOM(OM, "Obs")
  Imp <- SubOM(OM, "Imp")
  
  if(!Quiet)message("Replacing sub-models:", paste0(" ", Sub))
  for (x in 1:length(Sub)) {
    assign(Sub[x], SubOM(from, Sub[x]))
  }
  
  outOM <- new("OM", Stock, Fleet, Obs, Imp) 
  
  OMsl <- slotNames('OM')
  allSl <- c(slotNames('Stock'), slotNames('Fleet'), slotNames('Obs'), slotNames('Imp'))
  repsl <- OMsl[!OMsl %in% allSl]
  for (sl in repsl) slot(outOM, sl) <- slot(OM, sl)
  slot(outOM, 'Name') <- slot(OM, 'Name')
  
  outOM 
} 


#' Subset a Stock, Fleet, Obs, or Imp object from an OM object
#' 
#' A function that strips out a Stock, Fleet, Obs, or Imp object from a 
#' complete OM object. Mainly used for internal functions.
#' 
#' @param OM An operating model object (class OM)
#' @param Sub A character string specifying what object type to strip out
#' "Stock", "Fleet", "Obs", or "Imp"
#' @return An object of class Stock, Fleet, Obs, or Imp
#' @author A. Hordyk
#' @examples 
#' Stock <- SubOM(DLMtool::testOM, "Stock")
#' class(Stock)
#' @export 
SubOM <- function(OM, Sub=c("Stock", "Fleet", "Obs", "Imp")) {
  if (class(OM) !="OM") stop("OM must be of class OM ", call.=FALSE)
  Sub <- match.arg(Sub)
  temp <- new(Sub)
  
  slots <- slotNames(temp)
  for (X in seq_along(slots)) 
    slot(temp, slots[X]) <- slot(OM, slots[X]) 
  
  colon <- gregexpr(":", temp@Name)
  space <- gregexpr("  ", temp@Name)
  ind <- switch(Sub, Stock=1, Fleet=2, Obs=3, Imp=4)
  
  if (ind < 4) temp@Name <- substr(temp@Name, colon[[1]][ind]+1, space[[1]][ind]-1)
  if (ind == 4) temp@Name <- substr(temp@Name, colon[[1]][ind]+1, nchar(temp@Name))
  
  temp 
}
