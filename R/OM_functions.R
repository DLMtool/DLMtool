#' Forces correlation among operating model parameters for M, K, Linf and L50
#'
#' @description Uses typical correlations among estimated parameters to generate realistic samples for natural mortality rate (M), growth rate (K), maximum length (Linf) and length at 50% maturity (L50), these are placed in the cpars slot
#' @param OM An operating model object with M, growth, stock-recruitment and maturity parameters specified.
#' @param nsim The number of simulated values to create (note that OM@nsim will be used preferentially). 
#' @param plot Should the sampled parameters and distributions be plotted?
#' @return An object of class OM with a populated (or appended) cpars slot
#' @author T. Carruthers (Canadian DFO grant)
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

#' Predict missing life-history parameters
#' 
#' Predict missing life-history based on taxonomic information and 
#' hierachical model fitted to FishBase life-history parameters
#' 
#' @describeIn LH2OM Predict missing life-history and populate `OM@cpars`
#'
#' @param OM An object of class 'OM'
#' @param dist Character. Should parameters be sampled from a uniform (`unif`) or
#' normal (`norm`) distribution?
#' @param filterMK Logical. Should the predicted M and K parameters be filtered within the range specied in `inpars`or `OM`? 
#' e.g. `OM@M` and `OM@K`. Empty slots or slots with all values of 0 are considered unknown.
#' @param plot Logical. Should the plot be produced?
#' @param Class Optional higher order taxonomic information
#' @param Order Optional higher order taxonomic information
#' @param Family Optional higher order taxonomic information
#' @param msg Logical. Should messages be printed?
#' 
#' @return LH2OM: An OM with `OM@cpars` populated with `OM@nsim` samples of M, K, Linf and L50
#' @author A. Hordyk
#' @references Thorson, J. T., S. B. Munch, J. M. Cope, and J. Gao. 2017. 
#' Predicting life history parameters for all fishes worldwide. Ecological Applications. 27(8): 2262--2276  
#' @source \url{https://github.com/James-Thorson/FishLife/} 
#' @export
#' @examples
#' myOM<-LH2OM(DLMtool::testOM)
#' 
#' # drawing known parameters from normal distribution
#' myOM <- LH2OM(DLMtool::testOM, dist='norm')
#' 

#' 
LH2OM <- function(OM, dist=c("unif", "norm"), filterMK=FALSE, plot=TRUE,
                  Class = "predictive", Order = "predictive", 
                  Family = "predictive", msg=TRUE, db=DLMtool::LHdatabase) {
  if (class(OM) != 'OM') stop("OM must be class 'OM'")
  dist <- match.arg(dist)
  set.seed(OM@seed)
  if (length(OM@nsim)<1) OM@nsim <- 48
  if (length(OM@cpars)>0) {
    cnames <- names(OM@cpars)
    if (any(c("Linf", "L50", "M", "K") %in% cnames)) {
      message("Life-history parameters already in OM@cpars.\nReturning original OM")
      return(OM)
    }
  }
  
  sls <- c("Linf", "L50", "K", "M")
  for (sl in sls) {
    slval <- slot(OM, sl)
    if(any(is.na(slval)) | all(slval==0) | length(slval)<1) {
      assign(sl, NA)
    } else {
      assign(sl, slval)
    }
  }
  
  if (length(OM@M)>2) {
    message("Age-dependant M has been set in OM. (length(OM@M) > 2)\nReturning original OM")
    return(OM)
  }
  
  Genus <- unlist(strsplit(OM@Species, " "))[1]
  Species <- unlist(strsplit(OM@Species, " "))[2]
  if (is.na(Genus) || nchar(Genus)<1) Genus <- "predictive"
  if (is.na(Species) || nchar(Species)<1) Species <- "predictive"
  
  Out <- predictLH(inpars=list(Linf=Linf, L50=L50, K=K, M=M), 
                   Genus, Species, nsamp=OM@nsim, db=db, dist=dist, 
                   filterMK=filterMK, plot=plot, Class=Class, Order=Order, Family=Family, msg=msg)
  if (is.null(Out)) {
    message('Could not complete prediction. Returning original OM')
    return(OM)
  }
  OM@Linf <- c(0,0)
  OM@L50 <- c(0,0)
  OM@M <- c(0,0)
  OM@K <- c(0,0)
  
  OM@cpars$Linf <- Out$Linf
  OM@cpars$M <- Out$M
  OM@cpars$K <- Out$K
  OM@cpars$L50 <- Out$L50
  OM
}

#' @describeIn LH2OM Predict missing life-history based on taxonomic information and 
#' hierachical model fitted to FishBase life-history parameters
#' @param inpars A named list with lower and upper bounds of provided parameters: *Linf*, *L50*,
#' *K* and *M* (must be length 2). Unknown or missing parameters should not be included. For example, an empty list assumes
#' that all four life history parameters are unknown and need to be estimated. See Details below for more information.
#' @param Genus Character string specifying the Genus name. Optional. Default is 'predictive'  
#' @param Species Character string specifying the Species name. Optional. Default is 'predictive'. If full species name 
#' (Genus + Species) is not found if FishLife database (based on FishBase) higher order taxonomony will be used (e.g., Family)
#' for the predictions. 
#' @param nsamp The number of samples to return
#' @param db Database from FishLife model with fitted model results
#' 
#' @details The model predicts missing life-history parameters based on provided parameters and taxonomic information.
#' If both *M* and *K* are provided in `inpars` or  `OM`, *K* values are predicted and predictions filtered 
#' so that resulting *K* values are within bounds specified in `inpars$K` or `OM@K` (see `filterMK`).
#' 
#' If both *Linf* and *L50* are provided in `inpars` or `OM`, *L50* values are predicted and values in `inpars$L50` 
#' or `OM@L50` are ignored.
#'
#' @return predictLH: A data.frame with `nsamp` rows with `Linf`, `L50`, `K`, and `M` values.
#' @export
#'
#' @examples
#' # predict life-history parameters and return a data frame
#' 
#' # predict all life-history parameters 
#' Predicts <- predictLH(list(), "Katsuwonus", "pelamis")
#' head(Predicts)
#' 
#' # predict L50 from Linf, and M and K
#' Predicts <- predictLH(list(Linf=c(90, 95)), "Katsuwonus", "pelamis")
#' 
#' # predict L50 and K 
#' Predicts <- predictLH(list(Linf=c(90, 95), M=c(0.8, 0.9)), "Katsuwonus", "pelamis")
#' 
#' # predict L50 and K sampling Linf and M from normal distribution 
#' Predicts <- predictLH(list(Linf=c(90, 95), M=c(0.8, 0.9)), "Katsuwonus", "pelamis", dist='norm')
#' 
predictLH <- function(inpars=list(), Genus="predictive", Species="predictive", nsamp=100, 
                      db=DLMtool::LHdatabase, dist=c("unif", "norm"), 
                      filterMK=TRUE, plot=TRUE, Class = "predictive", Order = "predictive", 
                      Family = "predictive", msg=TRUE) {
  
  dist <- match.arg(dist)
  
  inpars_1 <- inpars
  # checks 
  names <- names(inpars)
  lens <- lapply(inpars, length) == 2
  valnames <- c("Linf", "L50",  "K",  "M") 
  if (!prod(names %in% valnames)) stop("invalid names in inpars. Valid names are: ", paste0(valnames, " "))
  
  notnames <- !valnames %in% names
  
  chkNA <- as.logical(unlist(lapply(lapply(inpars_1, is.na), prod)))
  notnames[chkNA] <- TRUE
  
  inpars_1[valnames[!valnames %in% names]] <- NA
  
  # names(inpars_1) <- valnames
  if (any(notnames)) {
    for(x in seq_along(valnames)) {
      if (notnames[x]) {
        nm <- valnames[x]
        inpars[[nm]] <- NA
        if (msg) message("Predicting ", nm) 
      }
    } 
  }
  
  if (prod(valnames[1:2] %in% names)) {
    if (lens['Linf']) {
      if (prod(valnames[1:2] %in% names)) {
        if (lens['Linf']) {
          if (msg) message("Predicting L50 from Linf")  
          lens['L50'] <- FALSE
          inpars$L50 <- NA
        }
        if (lens['L50']) {
          if (msg) message("Predicting Linf from L50")  
          lens['Linf'] <- FALSE
          inpars$Linf <- NA
        }
      }
      lens['L50'] <- FALSE
      inpars$L50 <- NA
    }
    if (lens['L50']) {
      if (msg) message("Predicting Linf from L50")  
      lens['Linf'] <- FALSE
      inpars$Linf <- NA
    }
  } 
  
  if (prod(valnames[3:4] %in% names)) {
    if (lens['M']) {
      if (msg) message("Predicting K from M")  
      lens['K'] <- FALSE
      inpars$K <- NA
    }
    if (lens['K']) {
      if (msg) message("Predicting M from K")  
      lens['M'] <- FALSE
      inpars$M <- NA
    }
  }
  multi <- 100
  filterM <- filterK <- FALSE
  if (prod(c("K", "M") %in% names) & filterMK & !(is.na(inpars_1$K) || is.na(inpars_1$M))) {
    if (all(is.na(inpars$M))) {
      filterM <- TRUE
      if (msg) message("Filtering predicted M within bounds: ", paste0(inpars_1$M, " "))
    }
    if (all(is.na(inpars$K))) {
      filterK <- TRUE
      if (msg) message("Filtering predicted K within bounds: ", paste0(inpars_1$K, " "))
    } 
    multi <- 500 
  }
  
  # get predictions from FishLife 
  taxa <- gettaxa(Class, Order, Family, Genus, Species, msg=msg)
  if (is.null(taxa)) return(NULL)
  if (class(db) != "list") stop("db must be database list from FishLife", call.=FALSE)
  Which <- grep(taxa, db$ParentChild_gz[,'ChildName'])
  mu <- db$ParHat$beta_gj[Which,]
  covar <- db$Cov_gjj[Which,,]
  names(mu) <- gsub("Loo", "Linf", names(mu))
  names(mu) <- gsub("Lm", "L50", names(mu))
  
  sampvals <- exp(as.data.frame(MASS::mvrnorm(nsamp*multi, mu, covar))) %>% dplyr::select(Linf, K, M, L50)
  sampvals$relLm <- sampvals$L50/sampvals$Linf 
  sampvals$MK <- sampvals$M/sampvals$K

  outpars <- inpars
  # sample input parameters
  
  if (dist=="unif") {
    for (x in names(outpars)) {
      if(!all(is.na(outpars[[x]]))) outpars[[x]] <- myrunif(nsamp*multi, min(outpars[[x]]), max(outpars[[x]]))
    }
  } else {
    for (x in names(outpars)) {
      if(!all(is.na(outpars[[x]]))) {
        varsd <- (max(outpars[[x]]) - mean(outpars[[x]]))/2
        varmean <- mean(outpars[[x]])
        outpars[[x]] <- rnorm(nsamp*multi, varmean, varsd)
      }
    }
  }
  
  # generate predicted parameters 
  missing <- lapply(lapply(outpars, is.na), prod) ==1
  missnm <- names(outpars)[missing]
  for (x in missnm) {
    if (x == "L50") {
      outpars$L50 <- outpars$Linf * sampvals$relLm
    }
    if (x == "Linf") {
      outpars$Linf <- outpars$L50 / sampvals$relLm
    }
    if (x == "M") {
      outpars$M <- outpars$K * sampvals$MK
    }
    if (x == "K") {
      outpars$K <- outpars$M / sampvals$MK
    }
  }
  

  # Still missing
  missing <- lapply(lapply(outpars, is.na), prod) ==1
  missnm <- names(outpars)[missing]
  for (x in missnm) {
    outpars[[x]] <- sampvals[[x]]
  }
  
  Out <- as.data.frame(do.call("cbind", outpars))
  
  # drop any samples where L50 > 0.95 Linf 
  ind <- Out$L50 > 0.95*Out$Linf
  Out <- Out[!ind,]
  
  if(filterK) {
    ind <- Out$K > min(inpars_1$K) & Out$K < max(inpars_1$K)
    if (sum(ind)<2) {
      warning('No samples of K within bounds: ', paste(as.character(inpars_1$K), collapse=" "), 
              "\nIgnoring bounds on K")
    } else {
      Out <- Out[ind,]
    }
    
  }
  if (filterM) {
    ind <- Out$M > min(inpars_1$M) & Out$M < max(inpars_1$M)
    if (sum(ind)<2) {
      warning('No samples of M within bounds: ', paste(as.character(inpars_1$M), collapse=" "), 
              "\nIgnoring bounds on M")
    } else {
      Out <- Out[ind,]
    }
  }
  
  if(nrow(Out) < nsamp) {
    warning("Could not generate ", nsamp, ' samples within specified bounds. Sampling with replacement')
    rows <- sample(1:nrow(Out), nsamp, replace=TRUE)
    Out <- Out[rows,]
  } else {
    Out <- Out[1:nsamp,]  
  }
  
  Out <- Out %>% dplyr::select(valnames)
  
  if(plot){ 
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(mfrow=c(4,4),mai=c(0.3,0.3,0.4,0.05),omi=c(0.02,0.02,0.3,0.02))
    
    colline=DLMtool::makeTransparent('blue',60)
    lwdline=4
    histcol='black'
    if (length(inpars_1)>1) {
      dobounds <- TRUE
      bounds<-matrix(NA,nrow=2, ncol=4)  
      count <- 0
      for (nm in valnames) {
        count <- count + 1
        tempval <- inpars_1[[nm]]
        if (length(tempval)>0) bounds[,count] <- range(tempval)
      }
    } else {
      dobounds <- FALSE
    }
    
    
    for(i in 1:4){
      if (length(inpars_1[[valnames[i]]])>1) {
        rng <- range(inpars_1[[valnames[i]]])  
      } else {
        rng <- c(NA, NA)
      }
      
      rng2 <- range(Out[,i])
      # options(warn=-1)
      if (!all(is.na(rng))) {
        rng2[1] <- min(c(min(rng, na.rm=TRUE), min(rng2, na.rm=TRUE)), na.rm=TRUE)
        rng2[2] <- max(c(max(rng, na.rm=TRUE), max(rng2, na.rm=TRUE)), na.rm=TRUE)
      }

      # options(warn=0)
      
      for(j in 1:4){
        if(i == j){
          if(i==1){
            hist(Out[,1],main="Asymptotic length (Linf)",col=histcol,border='white',xlab="",axes=F, xlim=rng2, ylab="")
            axis(1)
            abline(v=rng,col=colline,lwd=lwdline)
          }else if(i==2){
            hist(Out[,2],main="Length at 50% maturity (L50)",col=histcol,border='white',xlab="",axes=F, xlim=rng2, ylab="")
            axis(1)
            abline(v=rng,col=colline,lwd=lwdline)
          }else if(i==3){
            hist(Out[,3],main="Growth rate (K)",col=histcol,border='white',xlab="",axes=F, xlim=rng2, ylab="")
            axis(1)
            abline(v=rng,col=colline,lwd=lwdline)
          }else{
            hist(Out[,4],main="Natural mortality rate (M)",col=histcol,border='white',xlab="",axes=F, xlim=rng2, ylab="")
            axis(1)
            abline(v=rng,col=colline,lwd=lwdline)
          }
        }else{ # not positive diagonal
          
          plot(Out[,j],Out[,i],axes=F,col="white", xlab="", ylab="")
          if (dobounds) polygon(bounds[c(1,1,2,2),j],bounds[c(1,2,2,1),i],col=colline,border="white")
          points(Out[,j],Out[,i],pch=19)
          
          axis(1)
          axis(2)
        }
        
      }
    }
  }
  Out
}


# ImputeLH <- function(OM, samp.m=100, plot=TRUE, ign.bounds=TRUE) {
#   set.seed(OM@seed)
#   relL <- MK <- M <- K <- NULL
#   DBdata <- DLMtool::LHdata
#   DBdata <- DBdata%>% dplyr::select(M, K, relL, MK)
#   
#   if (class(OM) != "OM") stop('Object must be class "OM"', call. = FALSE)
#   
#   inM <- inK <- inLinf <- inL50 <- NULL
#   
#   for (nm in c("M", "K", "Linf", "L50")) {
#     if (length(OM@cpars[[nm]])>0) {
#       tempval <- OM@cpars[[nm]]
#     } else {
#       tempval <- slot(OM, nm)
#     }
#     assign(paste0('in', nm), tempval)
#   }
#   
#   parlist <- list(M=inM, K=inK, Linf=inLinf, L50=inL50)
#   means <- lapply(parlist, mean)
#   
#   Ms <- Ks <- Linfs <- L50s <- NULL
#   # Generate samples of known parameters - keep cpars
#   nsamp <- OM@nsim * samp.m
#   
#   ncpars <- max(unlist(lapply(parlist, length)))
#   if (ncpars < nsamp) {
#     sampRows <- sample(1:ncpars, nsamp, replace=TRUE)  
#   } else {
#     sampRows <- sample(1:ncpars, nsamp, replace=FALSE) 
#   }
#   
#   incpars <- data.frame(M=TRUE, K=TRUE, Linf=TRUE, L50=TRUE)
#   for (nm in c("M", "K", "Linf", "L50")) {
#     if (length(parlist[[nm]]) == 2) {
#       varsd <- (max(parlist[[nm]]) - mean(parlist[[nm]]))/2
#       varmean <- mean(parlist[[nm]])
#       var <- rnorm(nsamp*100, varmean, varsd)
#       var <- var[var>=parlist[[nm]][1] & var<=parlist[[nm]][2]][1:nsamp]
#       assign(paste0(nm, "s"), var)
#       incpars[nm] <- FALSE
#     } else {
#       var <- as.numeric(unlist((OM@cpars[nm])))[sampRows]
#       assign(paste0(nm, "s"), var)
#       incpars[nm] <- TRUE
#     }
#   }
# 
#   if (rowSums(incpars) == 4) stop("M, K, Linf, and L50 all in OM@cpars. Nothing to impute!", call.=FALSE)
#   indata_all <- data.frame(M=Ms, K=Ks, L50=L50s, Linf=Linfs)
#   if (any(indata_all$L50>indata_all$Linf)) stop("L50 is greater than Linf")
#   indata <- data.frame(M=Ms, K=Ks)
#   
#   ind <- apply(indata, 2, mean) == 0
#   indata[,ind] <- NA # drop any unknown parameters (e.g OM@Linf = c(0,0)= unknown)
#  
#   if (!all(is.na(indata$M)) & !incpars$K) {
#     indata$K <- NA # if M is present, impute K 
#     message("Imputing K from M")
#   }
#   if (!all(is.na(indata$K)) & !incpars$M) {
#     indata$M <- NA # if K is present, impute M 
#     message("Imputing M from K")
#   }
#   if (incpars$K & incpars$M & !(incpars$L50 | incpars$Linf)) {
#     message('M and K in cpars. Imputing L50/Linf')
#     indata$MK <- Ms/Ks
#   }
#   if ((incpars$L50 & incpars$Linf & !(incpars$K | incpars$M))) {
#     message('L50 and Linf in cpars. Imputing M/K')
#     indata$relL <- L50s/Linfs
#   }
#   
#   # if (!all(is.na(indata$Linf)) & !incpars$L50) indata$L50 <- NA # if Linf is present, impute L50 
#   # if (!all(is.na(indata$L50)) & !incpars$Linf) indata$Linf <- NA # if L50 is present, impute Linf
#   
#   ind <- colSums(apply(indata, 2, is.na)) == 0
#   indata <- indata[,ind, drop=FALSE]
#   if(length(indata$MK) > 0) message("M/K provided. Imputing L50/Linf")
#   if (ncol(indata) < 1) stop("At least one parameter but be provided - M or K", call.=FALSE)
#   
#   # Transform to approximate multivariate normal
#   tpow <- 1/3
#   tdata <- as.data.frame(DBdata^tpow)
#   # result <- MVN::mvn(data = tdata, mvnTest = "mardia", univariatePlot = "histogram")
#   
#   logs <- NULL
#   bounds <- matrix(c(1, 0.025, 0.8, # bounds for M
#                      2, 0.025, 1.5, # bounds for K 
#                      3, 0.2, 0.8, #  bounds for L50/Linf
#                      4, 0.2, 4), # bounds for M/K
#                    nrow=ncol(DBdata), ncol=3, byrow=TRUE)
#   if (!is.null(logs)) bounds[logs, 2:3] <- log(bounds[logs, 2:3] )
#   bounds[,2:3] <- bounds[,2:3]^tpow # transform bounds
#   
#   tindata <- as.data.frame(indata^tpow)
#   
#   alldata <- bind_rows(tdata, tindata)
#   # impute missing values
#   sink("temp")
#   mod <- Amelia::amelia(alldata, bounds=bounds, logs=logs, max.resample=5000, verbose=FALSE) 
#   sink()
#   unlink("temp")
#   
#   tMeanEsts <- Reduce("+", mod$imputations) / length(mod$imputations)
#   tMeanEsts <- tMeanEsts[(nrow(DBdata)+1):nrow(alldata),]
#   MeanEsts <- as.data.frame(tMeanEsts^(1/tpow))
#   
#   
#   
#   if (is.null(indata$K)) MeanEsts$K <- MeanEsts$M / MeanEsts$MK
#   if (is.null(indata$M)) MeanEsts$M <- MeanEsts$MK * MeanEsts$K 
#   MeanEsts$Linf <- indata_all$Linf # MeanEsts$relL * MeanEsts$Linf
#   MeanEsts$L50 <- indata_all$Linf * MeanEsts$relL 
#   
#   MeanEsts$relL <- NULL
#   MeanEsts$MK <- NULL
#   
#   
#   Out <- MeanEsts
# 
#   if (ign.bounds) {
#     ok <- array(TRUE, dim=dim(MeanEsts[,1:4]))
#     cols <- 1:4
#     for (xx in cols) {
#       rng <- range(parlist[[xx]])
#       ok[,xx] <- MeanEsts[,xx] >= rng[1] & MeanEsts[,xx] <= rng[2]
#     }
#     chk <- apply(ok, 2, prod)
#     oksamps <- chk == 1  
#     if (sum(chk) !=4) {
#       ind <- which(!oksamps)
#       message('Some predicted samples are outside specified bounds in OM: ', paste0(names(parlist)[ind], ", "), 
#               "\nIgnoring bounds for these parameters.")
#     }
#   }
#   
#   
#   
#   # filter to satisfy all bounds
#   if (!ign.bounds){
#     maxcount <- 50
#     chk <- sum(apply(ok, 1, prod)==1)
#     count <- 0
#     while (chk < OM@nsim & count < maxcount) {
#       count <- count + 1
#       oksamps <- as.data.frame(t(apply(ok, 2, sum)))  
#       ind <- which.min(oksamps)
#       MeanEsts[,ind] <- indata_all[,ind]
#       cols2 <- cols[!cols %in% ind]
#       ok <- array(TRUE, dim=dim(MeanEsts[,1:4]))
#       for (xx in cols2) {
#         rng <- range(parlist[[xx]])
#         if (mean(rng)!=0) {
#           ok[,xx] <- MeanEsts[,xx] >= rng[1] & MeanEsts[,xx] <= rng[2]  
#         }
#         
#       }
#       warning('Could not generate sufficient predicted samples within specified bounds for ', names(parlist)[ind], 
#               ".\nTry increase 'samp.m' and re-run. \nIgnoring predicted values for this parameter and using specified bounds")
#       chk <- sum(apply(ok, 1, prod)==1)
#     }
#     Out <- MeanEsts[apply(ok, 1, prod)==1,]
#   }
#   
#   
#   Out <- Out[1:OM@nsim,]
#   
#   if(plot){ 
#     par(mfrow=c(4,4),mai=c(0.3,0.3,0.4,0.05),omi=c(0.02,0.02,0.3,0.02))
#     
#     colline=makeTransparent('blue',60)
#     lwdline=4
#     histcol='black'
#     labs<-c("M","K","Linf","L50")
#     bounds<-matrix(unlist(lapply(parlist, range)),nrow=2)
#     
#     for(i in 1:4){
#       rng <- range(parlist[[i]])
#       rng2 <- range(Out[,i])
#       rng2[1] <- min(c(min(rng), min(rng2)))
#       rng2[2] <- max(c(max(rng), max(rng2)))
#       for(j in 1:4){
#         
#         if(i == j){
#           
#           if(i==1){
#             
#             hist(Out[,1],main="Natural mortality rate (M)",col=histcol,border='white',xlab="",axes=F, xlim=rng2)
#             axis(1)
#             abline(v=rng,col=colline,lwd=lwdline)
#             
#           }else if(i==2){
#             
#             hist(Out[,2],main="Growth rate (K)",col=histcol,border='white',xlab="",axes=F, xlim=rng2)
#             axis(1)
#             abline(v=rng,col=colline,lwd=lwdline)
#             
#             
#           }else if(i==3){
#             
#             hist(Out[,3],main="Asymptotic length (Linf)",col=histcol,border='white',xlab="",axes=F, xlim=rng2)
#             axis(1)
#             abline(v=rng,col=colline,lwd=lwdline)
#             
#             
#           }else{
#             
#             hist(Out[,4],main="Length at 50% maturity (L50)",col=histcol,border='white',xlab="",axes=F, xlim=rng2)
#             axis(1)
#             abline(v=rng,col=colline,lwd=lwdline)
#             
#           }
#           
#         }else{ # not positive diagonal
#           
#           plot(Out[,j],Out[,i],axes=F,col="white")
#           polygon(bounds[c(1,1,2,2),j],bounds[c(1,2,2,1),i],col=colline,border="white")
#           points(Out[,j],Out[,i],pch=19)
#           
#           axis(1)
#           axis(2)
#           
#         }
#         
#       }
#     }
#   }
#   
#   OM@cpars$M <- Out$M 
#   OM@cpars$K <- Out$K
#   OM@cpars$Linf <- Out$Linf
#   OM@cpars$L50 <- Out$L50
#   
#   OM@M <- c(0,0)
#   OM@K <- c(0,0)
#   OM@Linf <- c(0,0)
#   OM@L50 <- c(0,0)
#   OM
#   
# }




#' Replace an existing Stock, Fleet, Obs, or Imp object 
#' 
#' A function that replaces a Stock, Fleet, Obs, or Imp object from an 
#' OM with one from another object.
#' 
#' @param OM An operating model object (class OM) which will be updated with a sub-model from another OM
#' @param from An object of class `OM`, `Stock`, `Fleet`, `Obs`, or `Imp` to be replace the values in `OM`
#' @param Sub A character string specifying what object type to replace (only used if `from` is class `OM`)
#' "Stock", "Fleet", "Obs" or "Imp" (default is all four which is probably not what you want to do)
#' @param Name Character. Name for the new OM object (`OM@Name`)
#' @param silent Should messages be printed?
#' @return An object of class OM
#' @author A. Hordyk
#' @examples 
#' # Replace Stock 
#' OM <- DLMtool::testOM
#' OM2 <- Replace(OM, Blue_shark)
#' 
#' # Replace Fleet 
#' OM <- DLMtool::testOM
#' OM2 <- Replace(OM, Generic_DecE)
#' 
#' # Replace Fleet from another OM 
#' OM1 <- new("OM", Albacore, Generic_DecE, Perfect_Info, Overages)
#' OM2 <- new("OM", Blue_shark, Generic_IncE, Generic_Obs, Perfect_Imp)
#' OM1a <- Replace(OM1, OM2, "Fleet")
#' 
#' @export 
Replace <- function(OM, from,Sub=c("Stock", "Fleet", "Obs", "Imp"),  Name=NULL, silent=FALSE) {
  if (class(OM) =="character") OM <- get(OM)
  if (class(OM) !="OM") stop("OM must be of class OM ", call.=FALSE)
  if (class(from) =="character") from <- get(from)
  if (!class(from) %in% c("OM", "Stock", "Fleet", "Obs", "Imp")) 
    stop("from must be class `OM`, `Stock`, `Fleet`, `Obs`, or `Imp`", call.=FALSE)
  
  Stock <- SubOM(OM, "Stock")
  Fleet <- SubOM(OM, "Fleet")
  Obs <- SubOM(OM, "Obs")
  Imp <- SubOM(OM, "Imp")
  
  if (class(from) == "OM") {
    Sub <- match.arg(Sub, several.ok=TRUE)
    if (length(Sub)==4) warning("Replacing all OM components. Probably not what you want to do ...")
    
    if(!silent) message("Replacing sub-models:", paste0(" ", Sub))
    for (x in 1:length(Sub)) {
      assign(Sub[x], SubOM(from, Sub[x]))
    }
  
  } else {
    if(!silent) message("Replacing sub-model: ", class(from))
    assign(class(from), from)
  }
  outOM <- new("OM", Stock, Fleet, Obs, Imp) 
  

  OMsl <- slotNames('OM')
  allSl <- c(slotNames('Stock'), slotNames('Fleet'), slotNames('Obs'), slotNames('Imp'))
  repsl <- OMsl[!OMsl %in% allSl]
  for (sl in repsl) slot(outOM, sl) <- slot(OM, sl)
  
  if (is.null(Name)) {
    slot(outOM, 'Name') <- paste0('REPLACED -- ', slot(OM, 'Name'))
  } else {
    slot(outOM, 'Name') <- Name  
  }
  
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

# modified from FishLife::Search_species
# https://github.com/James-Thorson/FishLife
gettaxa <- function(Class = "predictive", Order = "predictive", 
                    Family = "predictive", Genus="predictive", Species="predictive",
                    ParentChild_gz=DLMtool::LHdatabase$ParentChild_gz,
                    msg=TRUE) {

  Match = 1:nrow(rfishbase::fishbase)
  if (Class != "predictive") 
    Match = Match[which(tolower(rfishbase::fishbase$Class[Match]) == tolower(Class))]
  if (Order != "predictive") 
    Match = Match[which(tolower(rfishbase::fishbase$Order[Match]) == tolower(Order))]
  if (Family != "predictive") 
    Match = Match[which(tolower(rfishbase::fishbase$Family[Match]) == tolower(Family))]
  if (Genus != "predictive") 
    Match = Match[which(tolower(rfishbase::fishbase$Genus[Match]) == tolower(Genus))]
  if (Species != "predictive") 
    Match = Match[which(tolower(rfishbase::fishbase$Species[Match]) == tolower(Species))]
  
  full_taxonomy <- c(Class, Order, Family, Genus, Species)
  spIn <- trimws(paste(gsub("predictive", "", full_taxonomy), collapse=" "))
  if (length(Match) == 0) {
    if(msg) message(spIn, ' not found in FishBase database')
    Class <- Order <- Family <- Genus <- Species <- "predictive"
  }
  full_taxonomy <- c(Class, Order, Family, Genus, Species)
  
  if (!all(Species == "predictive")) {
    if (length(unique(rfishbase::fishbase[Match, "Species"])) != 1) 
      stop("inputs are not unique")
    if (length(unique(rfishbase::fishbase[Match, "Species"])) == 1) 
      full_taxonomy[5] = unique(rfishbase::fishbase[Match, "Species"])[1]
  }
  if (!all(c(Species, Genus) == "predictive")) {
    if (length(unique(rfishbase::fishbase[Match, "Genus"])) != 1) 
      stop("inputs are not unique")
    if (length(unique(rfishbase::fishbase[Match, "Genus"])) == 1) 
      full_taxonomy[4] = unique(rfishbase::fishbase[Match, "Genus"])[1]
  }
  if (!all(c(Species, Genus, Family) == "predictive")) {
    if (length(unique(rfishbase::fishbase[Match, "Family"])) != 1) 
      stop("inputs are not unique")
    if (length(unique(rfishbase::fishbase[Match, "Family"])) == 1) 
      full_taxonomy[3] = unique(rfishbase::fishbase[Match, "Family"])[1]
  }
  if (!all(c(Species, Genus, Family, Order) == "predictive")) {
    if (length(unique(rfishbase::fishbase[Match, "Order"])) != 1) 
      stop("inputs are not unique")
    if (length(unique(rfishbase::fishbase[Match, "Order"])) == 1) 
      full_taxonomy[2] = unique(rfishbase::fishbase[Match, "Order"])[1]
  }
  if (!all(c(Species, Genus, Family, Order, Class) == "predictive")) {
    if (length(unique(rfishbase::fishbase[Match, "Class"])) != 1) 
      stop("inputs are not unique")
    if (length(unique(rfishbase::fishbase[Match, "Class"])) == 1) 
      full_taxonomy[1] = unique(rfishbase::fishbase[Match, 
                                                    "Class"])[1]
  }
  match_taxonomy = full_taxonomy
  Count = 1
  Group = NA
  while (is.na(Group)) {
    Group = match(paste(tolower(match_taxonomy), collapse = "_"), 
                  tolower(ParentChild_gz[, "ChildName"]))
    if (is.na(Group)) {
      match_taxonomy[length(match_taxonomy) - Count + 
                       1] = "predictive"
      Count = Count + 1
    }
  }
  
  fullname <- gsub("_", " ", ParentChild_gz[Group,  "ChildName"])
  ind <- !grepl("predictive", strsplit(fullname, " ")[[1]])
  if (all(!ind)) {
    if (msg) message("Predicting from all species in FishBase")
  } else if (any(!ind)) {
    if (msg) message("Closest match: ", fullname)
  } else {
    if (msg) message("Species match: ", fullname)
  }
  

  match_taxonomy = unique(as.character(Add_predictive(ParentChild_gz[Group, "ChildName"])))
  match_taxonomy
}

# from FishLife::Search_species
# https://github.com/James-Thorson/FishLife
Add_predictive = function(char_vec) {
  return_vec = char_vec
  for (i in 1:length(return_vec)) {
    vec = strsplit(as.character(return_vec[i]), "_")[[1]]
    return_vec[i] = paste(c(vec, rep("predictive", 5 - length(vec))), collapse = "_")
  }
  return(return_vec)
}


