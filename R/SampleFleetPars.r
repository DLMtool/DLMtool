
#' Sample Fleet Parameters
#'
#' @param Fleet An object of class 'Fleet' or class 'OM' 
#' @param Stock An object of class 'Stock' or a list of sampled Stock parameters. Ignored if 'Fleet' is class 'OM'
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
  
  if (class(Fleet) != "OM" & class(Stock) != "Stock" & class(Stock) != "list") stop("Must provide 'Stock' object", call.=FALSE)
  
  # Get custom pars if they exist
  if (class(Fleet) == "OM" && length(Fleet@cpars) > 0 && is.null(cpars)) cpars <- SampleCpars(Fleet@cpars, Fleet@nsim)  # custom parameters exist in Stock/OM object
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
  
  if (!exists("Find", inherits = FALSE)) {
    Deriv <- getEffhist(Esd, nyears, EffYears = EffYears, EffLower = EffLower, EffUpper = EffUpper)  # Historical fishing effort  
    Find <- Deriv[[1]]  # Calculate fishing effort rate
    dFfinal <- Deriv[[2]]  # Final gradient in fishing effort yr-1 
  } else {
    dFfinal <- NA
  }
  if (any(dim(Find) != c(nsim, nyears))) stop("Find must be matrix with dimensions: nsim, nyears")
  
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
  
  L5s <- LFSs <- Vmaxlens <- NULL  # initialize 
  
  if (Selnyears > 1) {   # change of selectivity in historical years 
    # length at 0.05% selectivity ascending
    L5s <- mapply(runif, n = nsim, min = Fleet@L5Lower, max = Fleet@L5Upper) * multi
    # first length at 100% selection
    LFSs <- mapply(runif, n = nsim, min = Fleet@LFSLower, max = Fleet@LFSUpper) *  multi
    # selectivity at maximum length
    Vmaxlens <- mapply(runif, n = nsim, min = Fleet@VmaxLower, max = Fleet@VmaxUpper)
  }
  
  
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
        SLarray[s,, yr] <- SelectFun(s, SL0.05=L5[yr, ], SL1=LFS[yr, ], MaxSel=Vmaxlen[yr, ], 
                                     maxlens=Len_age[, maxage, nyears], Lens=CAL_binsmid)
      }
    }
    
  }
  
  # == Calculate Selectivity at Age and Length ====
  
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
      s1 <- sapply(1:nsim, function(i) optimize(getSlope1, interval = c(0, 1e+05), 
                                                LFS = LFS[1, i], L0.05 = L5[1,i])$minimum)	
      if (all(Vmaxlen >= 0.99)) s2 <- rep(1E5, nsim)
      if (!all(Vmaxlen >= 0.99)) 
        s2 <- sapply(1:nsim, function(i) optimize(getSlope2, interval = c(0, 1e+05), 
                                                  LFS = LFS[1,i], s1=s1[i], maxlen=maxlen[i], 
                                                  MaxSel=Vmaxlen[1, i])$minimum)
      for (yr in 1:(nyears+proyears)) {
        # Calculate selectivity at age class 
        V[ , , yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[1,i], s1[i], s2[i], lens=Len_age[i,,yr])))
        # Calculate selectivity at length class 
        SLarray[,, yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[1,i], s1[i], s2[i], lens=CAL_binsmid)))   
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
        
        s1 <- sapply(1:nsim, function(i) optimize(getSlope1, interval = c(0, 1e+05), 
                                                  LFS = LFSs[i, X], L0.05 = L5s[i, X])$minimum)
        s2 <- sapply(1:nsim, function(i) optimize(getSlope2, interval = c(0, 1e+05), 
                                                  LFS = LFSs[i, X], s1=s1[i], maxlen=maxlen[i], 
                                                  MaxSel=Vmaxlens[i, X])$minimum)	
        for (yr in bkyears) {
          V[ , , yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[yr, i], s1[i], s2[i], lens=Len_age[i,,yr])))
          SLarray[,, yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[1,i], s1[i], s2[i], lens=CAL_binsmid)))   		 
        }
      }
      
      restYears <- max(SelYears):(nyears + proyears)
      L5[restYears, ] <- matrix(rep((L5s[, Selnyears]), length(restYears)), ncol = nsim, byrow = TRUE)
      LFS[restYears, ] <- matrix(rep((LFSs[, Selnyears]), length(restYears)), ncol = nsim, byrow = TRUE)
      Vmaxlen[restYears, ] <- matrix(rep((Vmaxlens[, Selnyears]), length(restYears)), ncol = nsim, byrow = TRUE)
      
      s1 <- sapply(1:nsim, function(i) optimize(getSlope1, interval = c(0, 1e+05), 
                                                LFS = LFSs[i, Selnyears], L0.05 = L5s[i, Selnyears])$minimum)
      s2 <- sapply(1:nsim, function(i) optimize(getSlope2, interval = c(0, 1e+05), 
                                                LFS = LFSs[i, Selnyears], s1=s1[i], maxlen=maxlen[i], 
                                                MaxSel=Vmaxlens[i, Selnyears])$minimum)	
      for (yr in restYears) { 
        V[ , , restYears] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[yr, i], s1[i], s2[i], lens=Len_age[i,,yr])))		
        SLarray[,, yr] <- t(sapply(1:nsim, function(i) TwoSidedFun(LFS[1,i], s1[i], s2[i], lens=CAL_binsmid))) 
      }	 
    }
  } # end of 'if V exists'
  
  if (any((dim(V) != c(nsim, maxage, proyears+nyears)))) 
    stop("V must have dimensions: nsim (", nsim,") maxage (", maxage, 
         ") proyears+nyears (", proyears+nyears, ") \nbut has ", 
         dim(V)[1], " ", dim(V)[2], " ", dim(V)[3], call.=FALSE)
  
  
  
  
  Fleetout$L5 <- L5  
  Fleetout$LFS <- LFS 
  Fleetout$Vmaxlen <- Vmaxlen 
  Fleetout$V <- V 
  Fleetout$SLarray <- SLarray 
  
  Fleetout 
}