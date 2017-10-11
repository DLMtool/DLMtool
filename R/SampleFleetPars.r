
#' Sample Fleet Parameters
#'
#' @param Fleet An object of class 'Fleet' or class 'OM' 
#' @param Stock An object of class 'Stock' or a list of sampled Stock parameters. 
#' Ignored if 'Fleet' is class 'OM'
#' @param nsim Number of simulations. Ignored if 'Fleet' is class 'OM'
#' @param nyears Number of historical years. Ignored if 'Fleet' is class 'OM'
#' @param proyears Number of projection years. Ignored if 'Fleet' is class 'OM'
#' @param cpars Optional named list of custom parameters. Ignored if 'Fleet' is class 'OM'
#'
#' @return A named list of sampled Fleet parameters
#' @export
#'
SampleFleetPars <- function(Fleet, Stock=NULL, nsim=NULL, nyears=NULL, proyears=NULL, cpars=NULL) {
  if (class(Fleet) != "Fleet" & class(Fleet) != "OM") 
    stop("First argument must be class 'Fleet' or 'OM'")
  
  if (class(Fleet) != "OM" & class(Stock) != "Stock" & class(Stock) != "list") 
    stop("Must provide 'Stock' object", call.=FALSE)
  
  # Get custom pars if they exist
  if (class(Fleet) == "OM" && length(Fleet@cpars) > 0 && is.null(cpars)) 
    cpars <- SampleCpars(Fleet@cpars, Fleet@nsim)  # custom parameters exist in Stock/OM object
  if (length(cpars) > 0) { # custom pars exist - assign to function environment 
    Names <- names(cpars)
    for (X in 1:length(Names)) assign(names(cpars)[X], cpars[[X]])
  }
  
  if (class(Fleet) == "OM") {
    nsim <- Fleet@nsim
    nyears <- Fleet@nyears 
    proyears <- Fleet@proyears
    StockPars <- SampleStockPars(Fleet, nsim, nyears, proyears, cpars)
    for (X in 1:length(StockPars)) assign(names(StockPars)[X], StockPars[[X]])
  }
  
  Fleet <- updateMSE(Fleet) # update to add missing slots with default values
  if (class(Stock) == "Stock") {
    Stock <- updateMSE(Stock) # update to add missing slots with default values
    # Sample Stock Pars - need some to calculate selectivity at age and length  
    StockPars <- SampleStockPars(Stock, nsim, nyears, proyears, cpars)
    for (X in 1:length(StockPars)) assign(names(StockPars)[X], StockPars[[X]])
  } 
  if (class(Stock) == "list") for (X in 1:length(Stock)) assign(names(Stock)[X], Stock[[X]])

  Fleetout <- list()
  
  # == Sample Historical Fishing Effort =====
  if (!exists("Esd", inherits = FALSE)) Esd <- runif(nsim, Fleet@Esd[1], Fleet@Esd[2])  # interannual variability in fishing effort (log normal sd)
  if (!exists("EffLower", inherits = FALSE)) EffLower <- Fleet@EffLower
  if (!exists("EffUpper", inherits = FALSE)) EffUpper <- Fleet@EffUpper 
  if (!exists("EffYears", inherits = FALSE)) EffYears <- Fleet@EffYears
  
  # if (max(EffYears) != nyears && max(EffYears) != 1) stop("Maximum EffYears (", max(EffYears), ") not equal to nyears (", nyears, ")")
  
  if (!exists("Find", inherits = FALSE)) {
    if (any(is.na(EffLower)) || any(is.na(EffUpper)) || any(is.na(EffYears))) {
      message("NAs in EffLower, EffUpper, or EffYears")
      Find <- matrix(NA, nsim, nyears)
      Deriv <- list(Find, rep(NA, nsim))
    } else {
      Deriv <- getEffhist(Esd, nyears, EffYears = EffYears, EffLower = EffLower, EffUpper = EffUpper)  # Historical fishing effort  
      Find <- Deriv[[1]]  # Calculate fishing effort rate    
    }

  }   

  if (!exists("dFfinal", inherits = FALSE)) {
    if (exists("Deriv", inherits = FALSE)) {
      dFfinal <- Deriv[[2]]  # Final gradient in fishing effort yr-1 
    } else {
      dFfinal <- rep(NA, nsim)
    }
  }
  if (any(dim(Find) != c(nsim, nyears))) stop("Find must be matrix with dimensions: nsim (", nsim, "), nyears (", nyears, ") but is: ", paste(dim(Find), ""))
  
  Fleetout$Esd <- Esd
  Fleetout$Find <- Find
  Fleetout$dFfinal <- dFfinal
  
  # === Spatial Targetting ====
  if (!exists("Spat_targ", inherits = FALSE))  Spat_targ <- runif(nsim, Fleet@Spat_targ[1], Fleet@Spat_targ[2])  # spatial targetting Ba^targetting param 
  
  Fleetout$Spat_targ <- Spat_targ
  
  # === Sample fishing efficiency parameters ====
  if (!exists("qinc", inherits = FALSE)) qinc <- runif(nsim, Fleet@qinc[1], Fleet@qinc[2])
  if (!exists("qcv", inherits = FALSE)) qcv <- runif(nsim, Fleet@qcv[1], Fleet@qcv[2])  # interannual variability in catchability
  
  # === Simulate future variability in fishing efficiency ====
  qmu <- -0.5 * qcv^2  # Mean
  if (!exists("qvar", inherits = FALSE)) qvar <- array(exp(rnorm(proyears * nsim, rep(qmu, proyears), rep(qcv, proyears))), c(nsim, proyears))  # Variations in interannual variation
  FinF <- Find[, nyears]  # Effort in final historical year
  
  Fleetout$qinc <- qinc
  Fleetout$qcv <- qcv
  Fleetout$qvar <- qvar
  Fleetout$FinF <- FinF
  
  # ==== Sample selectivity parameters ====
  Selnyears <- length(Fleet@SelYears)
  # are selectivity parameters relative to size at maturity?
  chk <- class(Fleet@isRel)
  if (length(Fleet@isRel) < 1) 
    Fleet@isRel <- "true"
  if (chk == "character") {
    chkRel <- tolower(Fleet@isRel)
    if (chkRel == "true" | Fleet@isRel == "1") multi <- L50
    if (chkRel == "false" | Fleet@isRel == "0")multi <- 1
  }
  if (chk == "numeric") {
    if (Fleet@isRel == 1) multi <- L50
    if (Fleet@isRel == 0) multi <- 1
  }
  if (!exists("L5", inherits = FALSE)) L5 <- runif(nsim, Fleet@L5[1], Fleet@L5[2]) * multi  # length at 0.05% selectivity ascending
  if (!exists("LFS", inherits = FALSE)) LFS <- runif(nsim, Fleet@LFS[1], Fleet@LFS[2]) * multi  # first length at 100% selection
  if (!exists("Vmaxlen", inherits = FALSE)) Vmaxlen <- runif(nsim, Fleet@Vmaxlen[1], Fleet@Vmaxlen[2])  # selectivity at maximum length
  
  Vmaxlen[Vmaxlen<=0] <- tiny
  L5s <- LFSs <- Vmaxlens <- NULL  # initialize 
  
  if (Selnyears > 1) {   # change of selectivity in historical years 
    # length at 0.05% selectivity ascending
    L5s <- mapply(runif, n = nsim, min = Fleet@L5Lower, max = Fleet@L5Upper) * multi
    # first length at 100% selection
    LFSs <- mapply(runif, n = nsim, min = Fleet@LFSLower, max = Fleet@LFSUpper) *  multi
    # selectivity at maximum length
    Vmaxlens <- mapply(runif, n = nsim, min = Fleet@VmaxLower, max = Fleet@VmaxUpper)
  } else {
    L5s <- LFSs <- Vmaxlens <- NA
  }
  Fleetout$L5s <- L5s
  Fleetout$LFSs <- LFSs
  Fleetout$Vmaxlens <- Vmaxlens
  
  
  # == Calculate Selectivity at Length ====
  nCALbins <- length(CAL_binsmid)
  SLarray <- array(NA, dim=c(nsim, nCALbins, nyears+proyears)) # Selectivity-at-length 

  if (exists("V", inherits=FALSE)) { # V has been passed in with custompars
    if(dim(V)[3] != proyears + nyears) V<-abind::abind(V,array(V[,,nyears],c(nsim,maxage,proyears)),along=3) # extend future Vulnerabiliy according to final historical vulnerability
    # assign L5, LFS and Vmaxlen - dodgy loop 
    # could calculate length at 5% selectivity from vB
    L5 <- matrix(NA, nrow = nyears + proyears, ncol = nsim)
    LFS <- matrix(NA, nrow = nyears + proyears, ncol = nsim)
    Vmaxlen <- matrix(NA, nrow = nyears + proyears, ncol = nsim)
    
    for (yr in 1:(nyears+proyears)) {
      for (s in 1:nsim) {
        ind <- min(which(V[s,,yr] >=0.05))
        L5[yr, s] <- Len_age[s, ind, yr]
        ind2 <- min(which(V[s,,yr] >=0.50))
        if (ind2 == ind) ind2 <- ind + 1
        LFS[yr, s] <- Len_age[s, ind2, yr]
        Vmaxlen[yr, s] <- V[s, maxage, yr]
        # SLarray[s,, yr] <- SelectFun(s, SL0.05=L5[yr, ], SL1=LFS[yr, ], MaxSel=Vmaxlen[yr, ], 
                                     # maxlens=Len_age[, maxage, nyears], Lens=CAL_binsmid)
      }
      SLarray[,, yr] <- t(sapply(1:nsim, getsel, lens=CAL_binsmidMat, L5[yr,], LFS[yr,], Vmaxlen[yr,], Linf) )
    }
    
  }
  
  # == Calculate Selectivity at Age and Length ====
  CAL_binsmidMat <- matrix(CAL_binsmid, nrow=nsim, ncol=length(CAL_binsmid), byrow=TRUE)
  
  
  if (!exists("V", inherits=FALSE)) { # don't run if V has been passed in with custompars 
    if (Selnyears <= 1) {    
      L5 <- matrix(L5, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
      LFS <- matrix(LFS, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
      Vmaxlen <- matrix(Vmaxlen, nrow = nyears + proyears, ncol = nsim, byrow = TRUE) 
      
      # ind <- which(LFS/matrix(Linf, nrow = proyears + nyears, ncol = nsim, byrow = TRUE) > 1, arr.ind = T)
      # if (length(ind) > 0) {
        # message("LFS too high (LFS > Linf) in some cases. \nDefaulting to LFS = 0.9 Linf for the affected simulations")
        # LFS[ind] <- Linf[ind[, 2]] * 0.9
      # } 
      
      
      # Calculate selectivity-at-age  curve 
      V <- array(NA, dim = c(nsim, maxage, nyears + proyears)) 
      # s1 <- sapply(1:nsim, function(i) optimize(getSlope1, interval = c(0, 1e+05), 
      #                                           LFS = LFS[1, i], L0.05 = L5[1,i])$minimum)	
      # if (all(Vmaxlen >= 0.99)) s2 <- rep(1E5, nsim)
      # Vmaxlen[Vmaxlen ==0] <- 0.001 # fix for when Vmaxlen == 0 
      # if (!all(Vmaxlen >= 0.99)) 
      #   s2 <- sapply(1:nsim, function(i) optimize(getSlope2, interval = c(0, 1e+05), 
      #                                             LFS = LFS[1,i], s1=s1[i], maxlen=maxlen[i], 
      #                                             MaxSel=Vmaxlen[1, i])$minimum)
    
      srs <- (Linf - LFS[1,]) / ((-log(Vmaxlen[1,drop=FALSE],2))^0.5) # selectivity parameters are constant for all years
      sls <- (LFS[1,] - L5[1, ]) /((-log(0.05,2))^0.5)
      
      # Calculate selectivity at length class 
      
      if (nsim>1) SelLength <- t(sapply(1:nsim, getsel, lens=CAL_binsmidMat, lfs=LFS[1, ], sls=sls, srs=srs))
      if (nsim == 1) SelLength <- getsel(1, lens=CAL_binsmidMat, lfs=LFS[1, ], sls=sls, srs=srs)
    
      for (yr in 1:(nyears+proyears)) {
        
        # Calculate selectivity at age class 
        # V[ , , yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[1,i], s1[i], s2[i], lens=Len_age[i,,yr])))
        if(nsim>1) V[ , , yr] <- t(sapply(1:nsim, getsel, lens=Len_age[,,yr], lfs=LFS[1,], sls=sls, srs=srs))
        
        if(nsim == 1) V[ , , yr] <- getsel(x=1, lens=t(matrix(Len_age[,,yr])), lfs=LFS[1,], sls=sls, srs=srs)
        # Calculate selectivity at length class 
        # SLarray[,, yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[1,i], s1[i], s2[i], lens=CAL_binsmid)))
        SLarray[,, yr] <- SelLength
      }	 
      
    }
    
    if (Selnyears > 1) {
      # More than one break point in historical selection pattern
      L5 <- matrix(0, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
      LFS <- matrix(0, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
      Vmaxlen <- matrix(0, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
      SelYears <- Fleet@SelYears
      
      ind <- which(LFSs/ matrix(Linf, nrow=nsim, ncol=Selnyears) > 1, arr.ind = T)
      if (length(ind) > 0) {
        message("LFS too high (LFS > Linf) in some cases. \nDefaulting to LFS = 0.9 Linf for the affected simulations")
        LFSs[ind] <- Linf[ind[, 1]] * 0.9
      }     
      

      # Calculate selectivity-at-age  curve 
      V <- array(NA, dim = c(nsim, maxage, nyears + proyears))     
      
      for (X in 1:(Selnyears - 1)) {	
        bkyears <- SelYears[X]:SelYears[X + 1]
        L5[bkyears, ] <- matrix(rep((L5s[, X]), length(bkyears)), ncol = nsim, byrow = TRUE)
        LFS[bkyears, ] <- matrix(rep((LFSs[, X]), length(bkyears)), ncol = nsim, byrow = TRUE)
        Vmaxlen[bkyears, ] <- matrix(rep((Vmaxlens[, X]), length(bkyears)), ncol = nsim, byrow = TRUE)
        
        srs <- (Linf - LFS[bkyears[1],]) / ((-log(Vmaxlen[bkyears[1],],2))^0.5) #
        sls <- (LFS[bkyears[1],] - L5[bkyears[1], ]) /((-log(0.05,2))^0.5)
      
        # Calculate selectivity at length class 
        SelLength <- t(sapply(1:nsim, getsel, lens=CAL_binsmidMat, lfs=LFS[bkyears[1],], sls=sls, srs=srs))
        
        
        # s1 <- sapply(1:nsim, function(i) optimize(getSlope1, interval = c(0, 1e+05), 
        #                                           LFS = LFSs[i, X], L0.05 = L5s[i, X])$minimum)
        # s2 <- sapply(1:nsim, function(i) optimize(getSlope2, interval = c(0, 1e+05), 
        #                                           LFS = LFSs[i, X], s1=s1[i], maxlen=maxlen[i], 
        #                                           MaxSel=Vmaxlens[i, X])$minimum)	
        for (yr in bkyears) {
          # V[ , , yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[yr, i], s1[i], s2[i], lens=Len_age[i,,yr])))
          # SLarray[,, yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[1,i], s1[i], s2[i], lens=CAL_binsmid)))   
          V[ , , yr] <-  t(sapply(1:nsim, getsel, lens=Len_age[,,yr], lfs=LFS[yr,], sls=sls, srs=srs))
          SLarray[,, yr] <- SelLength 
        }
      }
      
      restYears <- max(SelYears):(nyears + proyears)
      L5[restYears, ] <- matrix(rep((L5s[, Selnyears]), length(restYears)), ncol = nsim, byrow = TRUE)
      LFS[restYears, ] <- matrix(rep((LFSs[, Selnyears]), length(restYears)), ncol = nsim, byrow = TRUE)
      Vmaxlen[restYears, ] <- matrix(rep((Vmaxlens[, Selnyears]), length(restYears)), ncol = nsim, byrow = TRUE)
      
      # s1 <- sapply(1:nsim, function(i) optimize(getSlope1, interval = c(0, 1e+05), 
      #                                           LFS = LFSs[i, Selnyears], L0.05 = L5s[i, Selnyears])$minimum)
      # s2 <- sapply(1:nsim, function(i) optimize(getSlope2, interval = c(0, 1e+05), 
      #                                           LFS = LFSs[i, Selnyears], s1=s1[i], maxlen=maxlen[i], 
      #                                           MaxSel=Vmaxlens[i, Selnyears])$minimum)	
    
      srs <- (Linf - LFS[restYears[1],]) / ((-log(Vmaxlen[restYears[1],],2))^0.5) #
    
      sls <- (LFS[restYears[1],] - L5[restYears[1], ]) /((-log(0.05,2))^0.5)
      
      # Calculate selectivity at length class 
      SelLength <- t(sapply(1:nsim, getsel, lens=CAL_binsmidMat, lfs=LFS[restYears[1],], sls=sls, srs=srs))
      
      for (yr in restYears) { 
        # V[ , , restYears] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[yr, i], s1[i], s2[i], lens=Len_age[i,,yr])))		
        # SLarray[,, yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[1,i], s1[i], s2[i], lens=CAL_binsmid))) 
        V[ , , yr] <- t(sapply(1:nsim, getsel, lens=Len_age[,,yr], lfs=LFS[yr,], sls=sls, srs=srs))
        SLarray[,, yr] <- SelLength
        
      }	 
    }
  } # end of 'if V exists'
  
  if (any((dim(V) != c(nsim, maxage, proyears+nyears)))) 
    stop("V must have dimensions: nsim (", nsim,") maxage (", maxage, 
         ") proyears+nyears (", proyears+nyears, ") \nbut has ", 
         dim(V)[1], " ", dim(V)[2], " ", dim(V)[3], call.=FALSE)
  
  
  # == Sample Retention Parameters ====
  if(!exists("LR5", inherits = FALSE)) LR5 <- runif(nsim, min(Fleet@LR5), max(Fleet@LR5)) * multi
  if(!exists("LFR", inherits = FALSE)) LFR <- runif(nsim, min(Fleet@LFR), max(Fleet@LFR)) * multi
  if(!exists("Rmaxlen", inherits = FALSE)) Rmaxlen <- runif(nsim, min(Fleet@Rmaxlen), max(Fleet@Rmaxlen))
  if(!exists("DR", inherits = FALSE)) DR <- runif(nsim, min(Fleet@DR), max(Fleet@DR))
    
  if (any(LR5 > LFR)) stop('LR5 is greater than LFR', call.=FALSE)
  # == Calculate Retention Curve ====
  Rmaxlen[Rmaxlen<=0] <- tiny 
  LR5 <- matrix(LR5, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
  LFR <- matrix(LFR, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
  Rmaxlen <- matrix(Rmaxlen, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)
  DR <- matrix(DR, nrow = nyears + proyears, ncol = nsim, byrow = TRUE)

  
  # s1 <- sapply(1:nsim, function(i) optimize(getSlope1, interval = c(0, 1e+05), 
  #                                           LFS = LFR[1, i], L0.05 = LR5[1,i])$minimum)
  # if (all(Rmaxlen >= 0.99)) s2 <- rep(1E5, nsim)
  # if (!all(Rmaxlen >= 0.99)) 
  #   s2 <- sapply(1:nsim, function(i) optimize(getSlope2, interval = c(0, 1e+05), 
  #                                             LFS = LFR[1,i], s1=s1[i], maxlen=maxlen[i], 
  #                                             MaxSel=Rmaxlen[1, i])$minimum)
  
  srs <- (Linf - LFR[1,]) / ((-log(Rmaxlen[1,],2))^0.5) # selectivity parameters are constant for all years
  sls <- (LFR[1,] - LR5[1,]) /((-log(0.05,2))^0.5)
  
  RetLength <- t(sapply(1:nsim, getsel, lens=CAL_binsmidMat, lfs=LFR[1,], sls=sls, srs=srs))
  
  if (!exists("retA", inherits=FALSE)) {
    retA <- array(NA, dim = c(nsim, maxage, nyears + proyears)) # retention at age
    for (yr in 1:(nyears+proyears)) {
      # Calculate retention at age class 
      # retA[ , , yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFR[yr,i], s1[i], s2[i], lens=Len_age[i,,yr])))
      if (nsim>1) retA[ , , yr] <- t(sapply(1:nsim, getsel, lens=Len_age[,,yr], lfs=LFR[1,], sls=sls, srs=srs))
      if (nsim == 1) retA[ , , yr] <- getsel(1, lens=t(matrix(Len_age[,,yr])), lfs=LFR[1,], sls=sls, srs=srs)
      
    } 
  } else {
    # check dimensions 
    if (any((dim(retA) != c(nsim, maxage, proyears+nyears)))) 
      stop("retA must have dimensions: nsim (", nsim,") maxage (", maxage, 
           ") proyears+nyears (", proyears+nyears, ") \nbut has ", 
           dim(retA)[1], " ", dim(retA)[2], " ", dim(retA)[3], call.=FALSE) 
  }
  
  # 
  if (!exists("retL", inherits=FALSE)) {
    retL <- array(NA, dim = c(nsim, nCALbins, nyears + proyears)) # retention at length
    for (yr in 1:(nyears+proyears)) {
      # Calculate retention at length class 
      # retL[,, yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFR[yr,i], s1[i], s2[i], lens=CAL_binsmid)))   
      retL[,, yr] <- RetLength
    }
  } else {
    # check dimensions 
    if (any((dim(retL) != c(nsim, nCALbins, proyears+nyears)))) 
      stop("retL must have dimensions: nsim (", nsim,") nCALbins (", nCALbins, 
           ") proyears+nyears (", proyears+nyears, ") \nbut has ", 
           dim(retL)[1], " ", dim(retL)[2], " ", dim(retL)[3], call.=FALSE) 
  }
 

  V2 <- V
  SLarray2 <- SLarray
  
  # Apply general discard rate 
  dr <- aperm(abind::abind(rep(list(DR), maxage), along=3), c(2,3,1))
  retA <- (1-dr) * retA

  dr <- aperm(abind::abind(rep(list(DR), nCALbins), along=3), c(2,3,1))
  retL <- (1-dr) * retL
  
  # update realized vulnerablity curve with retention and dead discarded fish 
  Fdisc_array1 <- array(Fdisc, dim=c(nsim, maxage, nyears+proyears))
  V <- V * (retA + (1-retA)*Fdisc_array1) # Realised selection at age
  
  Fdisc_array2 <- array(Fdisc, dim=c(nsim, nCALbins, nyears+proyears))
  SLarray <- SLarray2 * (retL + (1-retL)*Fdisc_array2) # Realised selection at length
  
  # Realised Retention curves
  retA <- retA * V2
  retL <- retL * SLarray2
  
  Fleetout$Fdisc <- Fdisc
  Fleetout$Fdisc_array1 <- Fdisc_array1
  Fleetout$Fdisc_array2 <- Fdisc_array2
  Fleetout$LR5 <- LR5  
  Fleetout$LFR <- LFR 
  Fleetout$Rmaxlen <- Rmaxlen
  Fleetout$DR <- DR
  
  Fleetout$retA <- retA  # retention-at-age array - nsim, maxage, nyears+proyears
  Fleetout$retL <- retL  # retention-at-length array - nsim, nCALbins, nyears+proyears
  
  Fleetout$L5 <- L5  
  Fleetout$LFS <- LFS 
  Fleetout$Vmaxlen <- Vmaxlen 
  Fleetout$V <- V  # realized vulnerability-at-age
  Fleetout$SLarray <- SLarray # realized vulnerability-at-length
  Fleetout$V2 <- V2 # original vulnerablity-at-age curve 
  Fleetout$SLarray2 <- SLarray2 # original vulnerablity-at-length curve 
  
  Fleetout 
}
