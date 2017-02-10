# DLMtool Internal Functions 
# Various short functions that are used internally in the DLMtool.
# Functions are only used internally, and not directly accessible or
# available to users of the package and therefore do not appear in the
# help manual.


# Misc. Defintions
# -------------------------------------------------------------
utils::globalVariables(c("R0", "Mdb", "mod", "i", "nareas", "nsim", "dFfinal"))
tiny <- 1e-15  # define tiny variable
proportionMat <- TL <- Wa <- SurvWeiMat <- r <- lx <- logNormDensity <- sumlogNormDen <- NULL
proportionMat = vector()


#' Calculate FMSY and related metrics 
#' 
#' @param x internal parameter
#' @param Marray internal parameter
#' @param hs internal parameter
#' @param Mat_age internal parameter
#' @param Wt_age internal parameter
#' @param R0 internal parameter
#' @param V internal parameter
#' @param maxage internal parameter
#' @param nyears internal parameter
#' @param proyears internal parameter
#' @param Spat_targ internal parameter
#' @param mov internal parameter
#' @param SRrel internal parameter
#' @param aR internal parameter
#' @param bR internal parameter
#' 
#' @keywords internal
#' @export getFMSY
getFMSY <- function(x, Marray, hs, Mat_age, Wt_age, R0, V, maxage, nyears, 
    proyears, Spat_targ, mov, SRrel, aR, bR) {
    opt <- optimize(FMSYopt, log(c(0.001, 5)), Mc = Marray[x, nyears], 
        hc = hs[x], Mac = Mat_age[x, ], Wac = Wt_age[x, , nyears], R0c = R0[x], 
        Vc = V[x, ], maxage = maxage, nyears = nyears, proyears = proyears, 
        Spat_targc = Spat_targ[x], movc = mov[x, , ], SRrelc = SRrel[x], 
        aRc = aR[x, ], bRc = bR[x, ], Opt = T)
    return(FMSYopt(opt$minimum, Mc = Marray[x, nyears], hc = hs[x], Mac = Mat_age[x, 
        ], Wac = Wt_age[x, , nyears], R0c = R0[x], Vc = V[x,  ], maxage = maxage, 
        nyears = nyears, proyears = proyears, Spat_targc = Spat_targ[x], 
        movc = mov[x, , ], SRrelc = SRrel[x], aRc = aR[x, ], bRc = bR[x, 
            ], Opt = F))
}

#' Calculate FMSY and related metrics using Rcpp code 
#' 
#' @param x internal parameter
#' @param Marray internal parameter
#' @param hs internal parameter
#' @param Mat_age internal parameter
#' @param Wt_age internal parameter
#' @param R0 internal parameter
#' @param V internal parameter
#' @param maxage internal parameter
#' @param nyears internal parameter
#' @param proyears internal parameter
#' @param Spat_targ internal parameter
#' @param mov internal parameter
#' @param SRrel internal parameter
#' @param aR internal parameter
#' @param bR internal parameter
#' 
#' @keywords internal
#' @export getFMSY2	
getFMSY2 <- function(x, Marray, hs, Mat_age, Wt_age, R0, V, maxage, nyears, 
    proyears, Spat_targ, mov, SRrel, aR, bR, Control=1) {
    opt <- optimize(projOpt_cpp, log(c(0.001, 8)),
		Mc = Marray[x, nyears], hc = hs[x], Mac = Mat_age[x, ], Wac = Wt_age[x, , nyears], R0c = R0[x], 
        Vc = V[x, ,nyears], nyears = nyears, maxage = maxage, movc = mov[x, , ], Spat_targc = Spat_targ[x],
        SRrelc = SRrel[x], aRc = aR[x, ], bRc = bR[x, ], proyears = proyears, Control=Control)

	MSY <- -opt$objective 
	MSYs <- projOpt_cpp(lnIn = opt$minimum, 
		Mc = Marray[x, nyears], hc = hs[x], Mac = Mat_age[x, ], Wac = Wt_age[x, , nyears], R0c = R0[x], 
        Vc = V[x, ,nyears], nyears = nyears, maxage = maxage, movc = mov[x, , ], Spat_targc = Spat_targ[x],
        SRrelc = SRrel[x], aRc = aR[x, ], bRc = bR[x, ], proyears = proyears, Control=2)
    SSB_MSY <- MSYs[1]				
    B_MSY <- MSYs[2] 
    V_BMSY <- MSYs[3]
		
    F_MSYv <- -log(1 - (MSY/(V_BMSY+MSY))) 
	F_MSYb <- -log(1 - (MSY/(SSB_MSY+MSY))) 
    return(c(MSY = MSY, FMSY = F_MSYv, SSB = SSB_MSY, B = B_MSY, VB=V_BMSY, F_MSYb=F_MSYb))				
}

#' Internal function FMSY and related metrics 
#' 
#' @param lnF internal parameter
#' @param Mc internal parameter
#' @param hc internal parameter
#' @param Mac internal parameter
#' @param Wac internal parameter
#' @param R0c internal parameter
#' @param Vc internal parameter
#' @param maxage internal parameter
#' @param nyears internal parameter
#' @param proyears internal parameter
#' @param Spat_targc internal parameter
#' @param movc internal parameter
#' @param SRrelc internal parameter
#' @param aRc internal parameter
#' @param bRc internal parameter
#' @param Opt internal parameter
#' 
#' @keywords internal
#' @export FMSYopt
FMSYopt <- function(lnF, Mc, hc, Mac, Wac, R0c, Vc, maxage, nyears, proyears, 
    Spat_targc, movc, SRrelc, aRc, bRc, Opt = T) {
    
    FMSYc <- exp(lnF)
    nareas <- nrow(movc)
    # areasize<-c(asizec,1-asizec)
    idist <- rep(1/nareas, nareas)
    for (i in 1:300) idist <- apply(array(idist, c(2, 2)) * movc, 2, sum)
    
    N <- array(exp(-Mc * ((1:maxage) - 1)) * R0c, dim = c(maxage, nareas)) * 
        array(rep(idist, each = maxage), dim = c(maxage, nareas))
    SSN <- Mac * N  # Calculate initial spawning stock numbers
    Biomass <- N * Wac
    VBiomass <- Biomass * Vc
    SSB <- SSN * Wac  # Calculate spawning stock biomass
    
    B0 <- sum(Biomass)
    VB0 <- sum(VBiomass)
    R0a <- idist * R0c
    SSB0 <- apply(SSB, 2, sum)
    SSBpR <- SSB0/R0a
    
    N <- N/2  # Calculate spawning stock biomass per recruit
    SSN <- Mac * N  # Calculate initial spawning stock numbers
    Biomass <- N * Wac
    SSB <- SSN * Wac  # Calculate spawning stock biomass
    
    for (y in 1:nyears) {
        # set up some indices for indexed calculation
        dis <- apply(Vc * Biomass, 2, sum)/sum(Vc * Biomass)
        targ <- (dis^Spat_targc)/mean(dis^Spat_targc)
        FMc <- array(FMSYc * Vc, dim = c(maxage, nareas)) * array(rep(targ, 
            each = maxage), dim = c(maxage, nareas))  # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
        Zc <- FMc + Mc
        CN <- N * (1 - exp(-Zc)) * (FMc/Zc)
        CB <- CN * Wac
        
        N[2:maxage, ] <- N[1:(maxage - 1), ] * exp(-Zc[1:(maxage - 1), 
            ])  # Total mortality
        if (SRrelc == 1) {
            N[1, ] <- (0.8 * R0a * hc * apply(SSB, 2, sum))/(0.2 * SSBpR * 
                R0a * (1 - hc) + (hc - 0.2) * apply(SSB, 2, sum))  # Recruitment assuming regional R0 and stock wide steepness
        } else {
            N[1, ] <- aRc * apply(SSB, 2, sum) * exp(-bRc * apply(SSB, 
                2, sum))
        }
        # print(N[1])
        N[1, ] <- apply(array(N[1, ], c(2, 2)) * movc, 2, sum)
        SSN <- N * Mac
        SSB <- SSN * Wac
        Biomass <- N * Wac
        VBiomass <- Biomass * Vc
        # print(sum(Biomass))
    }  # end of year
    
    CBc <- sum(CB)
    if (Opt) {
        return(-CBc)
    } else {
        return(c(MSY = CBc, FMSY = -log(1 - (CBc/(sum(VBiomass) + CBc))), 
            SSB = sum(SSB), SSB_SSB0 = sum(SSB)/sum(SSB0), B = sum(N * 
                Wac), B_B0 = sum(N * Wac)/B0))
    }
}



getFhist <- function(nsim, Esd, nyears, dFmin, dFmax, bb) {
    
    ne <- nsim * 3  # Number of simulated effort datasets
    dEfinal <- runif(ne, dFmin, dFmax)  #(exp(rnorm(ne,mean=demu,sd=desd))-1)*6               # Sample the final gradient in effort
    a <- (dEfinal - bb)/nyears  # Derive slope to get there from intercept
    a <- array(a, dim = c(ne, nyears))  # Slope array
    bb <- array(bb, dim = c(ne, nyears))  # Intercept array
    x <- array(rep(1:nyears, each = ne), dim = c(ne, nyears))  # Year array
    dE <- a * x + bb  # Change in effort
    # E<-array(NA,dim=c(ne,nyears)) # Define total effort array
    # E[,1]<-dE[,1] for(y in 2:nyears){ E[,y]<-apply(dE[,1:y],1,sum) }
    # E<-E/array(apply(E,1,mean),dim=c(ne,nyears)) # Standardise Effort to
    # average 1
    E2 <- t(apply(dE, 1, cumsum))  # Define total effort array
    E2 <- E2/array(apply(E2, 1, mean), dim = c(ne, nyears))  # Standardise Effort to average 1
    E <- E2  # 
    cond <- apply(E, 1, min) > 0
    pos <- (1:ne)[cond]
    pos <- pos[1:nsim]
    # environment('dEfinal')<-asNamespace('DLMtool')#assign('dFfinal',dEfinal[pos],envir=.GlobalEnv)
    
    E <- E[pos, ]  # Sample only those without negative effort
    Emu <- -0.5 * Esd^2
    Eerr <- array(exp(rnorm(nyears * nsim, rep(Emu, nyears), rep(Esd, nyears))), 
        c(nsim, nyears))
    outy <- new("list")
    outy[[1]] <- E * Eerr
    outy[[2]] <- dEfinal[pos]
    outy
}


#' Internal Get Reference F 
#' 
#' @param x internal parameter
#' @param Marray internal parameter
#' @param Wt_age internal parameter
#' @param Mat_age internal parameter
#' @param Perr internal parameter
#' @param N_s internal parameter
#' @param SSN_s internal parameter
#' @param N_s internal parameter
#' @param Biomass_s internal parameter
#' @param VBiomass_s internal parameter
#' @param SSB_s internal parameter
#' @param Vn internal parameter
#' @param hs internal parameter
#' @param R0a internal parameter
#' @param nyears internal parameter
#' @param proyears internal parameter
#' @param nareas internal parameter
#' @param maxage internal parameter
#' @param mov internal parameter
#' @param SSBpR internal parameter
#' @param aR internal parameter
#' @param bR internal parameter
#' @param SRrel internal parameter
#' @param Spat_targ internal parameter
#' 
#' @keywords internal
#' @export getFref
getFref <- function(x, Marray, Wt_age, Mat_age, Perr, N_s, SSN_s, Biomass_s, 
    VBiomass_s, SSB_s, Vn, hs, R0a, nyears, proyears, nareas, maxage, mov, 
    SSBpR, aR, bR, SRrel, Spat_targ) {
    
    opt <- optimize(doprojPI, log(c(0.001, 10)), Mvec = Marray[x, (nyears + 1):(nyears + proyears)], 
	  Wac = Wt_age[x, , (nyears + 1):(nyears + proyears)], Mac = Mat_age[x, ], 
	    Pc = Perr[x, (nyears + 1):(nyears + proyears)], N_c = N_s[x, , ], 
		SSN_c = SSN_s[x, , ], Biomass_c = Biomass_s[x, , ], 
		VBiomass_c = VBiomass_s[x, , ], SSB_c = SSB_s[x, , ], Vc = Vn[x, , ], 
		hc = hs[x], R0ac = R0a[x, ], proyears, nareas, maxage, movc = mov[x, , ], 
		SSBpRc = SSBpR[x], aRc = aR[x, ], bRc = bR[x, ], SRrelc = SRrel[x], 
        spat_targ = Spat_targ[x])
    # print(exp(opt$minimum))
    return(-opt$objective)
    
}


#' Internal Get Reference F using Rcpp 
#' 
#' @param x internal parameter
#' @param Marray internal parameter
#' @param Wt_age internal parameter
#' @param Mat_age internal parameter
#' @param Perr internal parameter
#' @param N_s internal parameter
#' @param SSN_s internal parameter
#' @param N_s internal parameter
#' @param Biomass_s internal parameter
#' @param VBiomass_s internal parameter
#' @param SSB_s internal parameter
#' @param Vn internal parameter
#' @param hs internal parameter
#' @param R0a internal parameter
#' @param nyears internal parameter
#' @param proyears internal parameter
#' @param nareas internal parameter
#' @param maxage internal parameter
#' @param mov internal parameter
#' @param SSBpR internal parameter
#' @param aR internal parameter
#' @param bR internal parameter
#' @param SRrel internal parameter
#' @param Spat_targ internal parameter
#' 
#' @keywords internal
#' @export getFref2
getFref2 <- function(x, Marray, Wt_age, Mat_age, Perr, N_s, SSN_s, Biomass_s, 
    VBiomass_s, SSB_s, Vn, hs, R0a, nyears, proyears, nareas, maxage, mov, 
    SSBpR, aR, bR, SRrel, Spat_targ) {
    
    opt <- optimize(doprojPI_cpp, log(c(0.001, 10)), Mvec = Marray[x, (nyears + 1):(nyears + proyears)], 
	  Wac = Wt_age[x, , (nyears + 1):(nyears + proyears)], Mac = Mat_age[x, ], 
	    Pc = Perr[x, (nyears + 1):(nyears + proyears)], N_c = N_s[x, , 1,], 
		SSN_c = SSN_s[x, , 1, ], Biomass_c = Biomass_s[x, , 1, ], 
		VBiomass_c = VBiomass_s[x, , 1, ], SSB_c = SSB_s[x, , 1, ], Vc = Vn[x, , ], 
		hc = hs[x], R0ac = R0a[x, ], proyears, nareas, maxage, movc = mov[x, , ], 
		SSBpRc = SSBpR[x], aRc = aR[x, ], bRc = bR[x, ], SRrelc = SRrel[x], 
        Spat_targc = Spat_targ[x])
    # print(exp(opt$minimum))
    return(-opt$objective)
    
}

			
#' Internal Projection Function
#' 
#' @param lnF internal parameter
#' @param Mvec internal parameter
#' @param Wac internal parameter
#' @param Pc internal parameter
#' @param N_c internal parameter
#' @param SSN_c internal parameter
#' @param Biomass_c internal parameter
#' @param VBiomass_c internal parameter
#' @param SSB_c internal parameter
#' @param Vc internal parameter
#' @param hc internal parameter
#' @param R0ac internal parameter
#' @param proyears internal parameter
#' @param nareas internal parameter
#' @param maxage internal parameter
#' @param movc internal parameter
#' @param SSBpRc internal parameter
#' @param aRc internal parameter
#' @param bRc internal parameter
#' @param SRrelc internal parameter
#' @param spat_targ internal parameter
#' 
#' @keywords internal
#' @export doprojPI
doprojPI <- function(lnF, Mvec, Wac, Mac, Pc, N_c, SSN_c, Biomass_c, VBiomass_c, 
    SSB_c, Vc, hc, R0ac, proyears, nareas, maxage, movc, SSBpRc, aRc, bRc, 
    SRrelc, spat_targ) {
    
    FF <- exp(lnF)
    
    N_P <- array(NA, dim = c(maxage, proyears, nareas))
    Biomass_P <- array(NA, dim = c(maxage, proyears, nareas))
    VBiomass_P <- array(NA, dim = c(maxage, proyears, nareas))
    SSN_P <- array(NA, dim = c(maxage, proyears, nareas))
    SSB_P <- array(NA, dim = c(maxage, proyears, nareas))
    FM_P <- array(NA, dim = c(maxage, proyears, nareas))
    Z_P <- array(NA, dim = c(maxage, proyears, nareas))
    CB_P <- rep(NA, proyears)
    
    # AYR<-as.matrix(expand.grid(1:maxage,1,1:nareas))
    AYR <- as.matrix(expand.grid(1:maxage, 1, 1:nareas))
    YA <- as.matrix(expand.grid(1, 1:maxage))  # Projection year
    Y <- YA[, 1]
    A <- YA[, 2]
    AY <- YA[, c(2, 1)]
    R <- matrix(AYR[, 3])
    
    N_P[AYR] <- N_c  #[AYRL]
    SSN_P[AYR] <- SSN_c  #SSN[AYRL]
    Biomass_P[AYR] <- Biomass_c  #[AYRL]
    VBiomass_P[AYR] <- VBiomass_c  #[AYRL]
    SSB_P[AYR] <- SSB_c  #[AYRL]
    
    fishdist <- (apply(VBiomass_P[, 1, ], 2, sum)^spat_targ)/mean(apply(VBiomass_P[, 
        1, ], 2, sum)^spat_targ)  # spatial preference according to spatial biomass
    fishdist <- matrix(fishdist, ncol = nareas)
    FM_P[AYR] <- FF * Vc[AY] * fishdist[R]
    # FM_P[AYR]<-FF*Vc[A]
    Z_P[AYR] <- FM_P[AYR] + Mvec[Y]
    
    for (y in 2:proyears) {
        AY1R <- as.matrix(expand.grid(1:maxage, y - 1, 1:nareas))
        AYR <- as.matrix(expand.grid(1:maxage, y, 1:nareas))
        Y <- AYR[, 2]
        A <- AYR[, 1]
        AY <- AYR[, 1:2]
        R <- AYR[, 3]
        A2YR <- as.matrix(expand.grid(2:maxage, y, 1:nareas))
        A1YR <- as.matrix(expand.grid(1:(maxage - 1), y - 1, 1:nareas))
        A1Y <- as.matrix(expand.grid(1:(maxage - 1), y - 1))
        
        indMov <- as.matrix(expand.grid(1:nareas, 1:nareas, y, 1:maxage)[4:1])
        indMov2 <- indMov[, c(1, 2, 3)]
        indMov3 <- indMov[, c(3, 4)]
        
        N_P[A2YR] <- N_P[A1YR] * exp(-Z_P[A1YR])  # Total mortality
        if (SRrelc == 1) {
            N_P[1, y, ] <- Pc[y] * (0.8 * R0ac * hc * apply(SSB_P[, y - 1, ], 2, sum))/
			  (0.2 * SSBpRc * R0ac * (1 - hc) + (hc - 0.2) * apply(SSB_P[, y - 1, ], 2, sum))  # Recruitment assuming regional R0 and stock wide steepness
        } else {
            N_P[1, y, ] <- Pc[y] * aRc * apply(SSB_P[, y - 1, ], 2, sum) * 
                exp(-bRc * apply(SSB_P[, y - 1, ], 2, sum))
        }
        
        temp <- array(N_P[indMov2] * movc[indMov3], dim = c(nareas, nareas,  maxage))  # Move individuals		
        N_P[, y, ] <- apply(temp, c(3, 1), sum)
        
        Biomass_P[AYR] <- N_P[AYR] * Wac[AY]  # Calculate biomass
        VBiomass_P[AYR] <- Biomass_P[AYR] * Vc[AY]  # Calculate vulnerable biomass
        SSN_P[AYR] <- N_P[AYR] * Mac[A]  # Calculate spawning stock numbers
        SSB_P[AYR] <- SSN_P[AYR] * Wac[AY]  # Calculate spawning stock biomass
        
        fishdist <- (apply(VBiomass_P[, y, ], 2, sum)^spat_targ)/mean(apply(VBiomass_P[, 
            y, ], 2, sum)^spat_targ)  # spatial preference according to spatial biomass
        FM_P[AYR] <- FF * Vc[AY] * fishdist[R]
        # FM_P[AYR]<-FF*Vc[A]
        
        Z_P[AYR] <- FM_P[AYR] + Mvec[Y]
        CNtemp <- N_P[, y, ] * exp(Z_P[, y, ]) * (1 - exp(-Z_P[, y, ])) * 
            (FM_P[, y, ]/Z_P[, y, ])
        CB_P[y] <- sum(Biomass_P[, y, ] * exp(Z_P[, y, ]) * (1 - exp(-Z_P[, 
            y, ])) * (FM_P[, y, ]/Z_P[, y, ]))
        
        CNtemp <- (FM_P[, y, ]/Z_P[, y, ] * N_P[, y, ] * (1 - exp(-Z_P[, 
            y, ])))
        CB_P[y] <- sum(FM_P[, y, ]/Z_P[, y, ] * Biomass_P[, y, ] * (1 - 
            exp(-Z_P[, y, ])))
        
        # temp <- sum(CNtemp*Wac[AY]) print(c(CB_P[y], temp))
        
        # CB_P[y] <- sum(CNtemp*Wac[AY])
    }  # end of year
    
    # plot(CB_P, type='l', ylim=c(0, max(CB_P, na.rm=TRUE))) Biomass_P[,1,]
    # apply(N_P[,1:6,1:2], c(1,2), sum)
    
    return(-mean(CB_P[(proyears - min(4, (proyears - 1))):proyears], na.rm = T))  
}


# Calculate initial spatial distribution
getinitdist <- function(tol, mov, indMain) {
    init <- array(1/nareas, dim = c(nsim, nareas))
    ind4 <- as.matrix(cbind(rep(1:nsim, each = nareas * nareas), indMain[, 
        2]))
    i <- 0
    delta <- 1
    # for(i in 1:100){
    while (delta > tol) {
        i <- i + 1
        trial <- init
        temp <- array(init[ind4] * mov[indMain], dim = c(nareas, nareas, 
            nsim))
        init <- apply(temp, c(3, 1), sum)
        delta <- max((trial - init)^2)
    }
    print(paste("Converged in ", i, " iterations"))
    init
}

# Demographic model
demofn <- function(log.r, M, amat, sigma, K, Linf, to, hR, maxage, a, b) demographic2(log.r, 
    M, amat, sigma, K, Linf, to, hR, maxage = maxage, a, b)$epsilon

demographic2 = function(log.r, M, amat, sigma, K, Linf, to, hR, maxage, 
    a, b) {
    # switch on and off to use either S or m in MC simulations
    r = exp(log.r)
    lx = exp(-M)^((1:maxage) - 1)  #survivorship
    logNormDensity = (dnorm(x = log((1:maxage)), mean = log(amat), sd = sigma))/(1:maxage)  #Maturity ogive calculation
    logNormDensity[1] = 0
    sumlogNormDen = sum(logNormDensity)
    NormalisedMaturity = logNormDensity/sumlogNormDen
    proportionMat[1] = NormalisedMaturity[1]
    for (i in 2:maxage) proportionMat[i] = proportionMat[i - 1] + NormalisedMaturity[i]
    TL = Linf * (1 - exp(-K * ((1:maxage) - to)))  #length at age
    Wa = a * TL^b  #wegith at age
    SurvWeiMat = lx * Wa * proportionMat  #survivorship X weight X maturity
    SBPR = sum(SurvWeiMat)  #Spawner biomass per recruit
    RPS = 1/(SBPR * (1 - hR)/(4 * hR))  # Beverton Holt
    # RPS=(5*hR)^(5/4)/SBPR # Ricker Recruitment per spawner biomass
    RPF = Wa * proportionMat * RPS  #Recruits per female
    Lotka = lx * RPF * exp(-(1:maxage) * r)
    sumLotka = sum(Lotka)
    epsilon = (1 - sumLotka)^2  #objective function
    return(list(epsilon = epsilon, r = r))
}

# Various functions for determing selectivity curve from input
# parameters difference in density from 0.05 given a standard deviation
# sd1 (sd_asc) and age at maximum vulnerability modo
densnorm <- function(sd1) {
    (0.05 - (dnorm(0, mod[i], sd1)/dnorm(mod[i], mod[i], sd1)))^2
}

densnormasc <- function(sd1, age_05, mody) {
    (0.05 - (dnorm(age_05, mody, sd1)/dnorm(mody, mody, sd1)))^2
}

getsdasc <- function(sm, age05, mod) {
    optimize(densnormasc, interval = c(0.5, 100), age_05 = age05[sm], mody = mod[sm])$minimum
}

densnormdesc <- function(sd2, V_maxage, maxy, mody) {
    (V_maxage - (dnorm(maxy, mody, sd2)/dnorm(mody, mody, sd2)))^2
}

getsddesc <- function(sm, Vmaxage, maxage, mod) {
    optimize(densnormdesc, interval = c(0.5, 10000), V_maxage = Vmaxage[sm], 
        maxy = maxage, mody = mod[sm])$minimum
}

getDNvulnS <- function(mod, age05, Vmaxage, maxage, nsim) {
    sd_asc <- sapply(1:nsim, getsdasc, age05 = age05, mod = mod)
    sd_desc <- sapply(1:nsim, getsddesc, Vmaxage = Vmaxage, maxage = maxage, 
        mod = mod)
    V <- array(NA, dim = c(nsim, maxage))
    for (i in 1:nsim) {
        V[i, 1:ceiling(mod[i])] <- dnorm(1:ceiling(mod[i]), mod[i], sd_asc[i])
        V[i, (1 + ceiling(mod[i])):maxage] <- dnorm((1 + ceiling(mod[i])):maxage, 
            mod[i], sd_desc[i])
        V[i, (1 + ceiling(mod[i])):maxage] <- V[i, (1 + ceiling(mod[i])):maxage]/V[i, 
            1 + ceiling(mod[i])]  #/V[i,floor(mod[i])+1]
        V[i, 1:ceiling(mod[i])] <- V[i, 1:ceiling(mod[i])]/dnorm(mod[i], 
            mod[i], sd_asc[i])  #,mod[i],sd_asc[i])#V[i,floor(mod[i])]
        
    }
    outy <- new("list")
    outy[[1]] <- V
    outy[[2]] <- mod - 1.18 * sd_asc
    outy
}



#' Internal Two sided selectivity curve
#'
#' @param LFS first length at fulll selection 
#' @param s1 ascending slope
#' @param s2 descending slope 
#' @param lens vector of lengths 
#' @keywords internal
#' @export TwoSidedFun
TwoSidedFun <- function(LFS, s1, s2, lens) {
    Sl <- rep(0, length(lens))
    Sl[lens < LFS] <- exp(-((lens[lens < LFS] - LFS)^2)/(2 * s1^2))
    Sl[lens >= LFS] <- exp(-((lens[lens >= LFS] - LFS)^2)/(2 * s2^2))
    return(Sl)
}

#' Internal function to calculate ascending slope of selectivity curve 
#'
#' @param s1 slope 1 
#' @param LFS first length at fulll selection 
#' @param L0.05 length at 5 percent selection 
#' @keywords internal
#' @export getSlope1
getSlope1 <- function(s1, LFS, L0.05) 
  (0.05 - TwoSidedFun(LFS, s1 = s1, s2 = 1E5, lens=L0.05))^2


#' Internal function to calculate slope
#'
#' @param s2 desceding slope 
#' @param LFS length one
#' @param s1 ascending slope 
#' @param maxlen length of oldest age class
#' @param MaxSel selectivity of maxlen 
#' @keywords internal
#' @export getSlope2
getSlope2 <- function(s2, LFS, s1, maxlen, MaxSel) 
  (MaxSel - TwoSidedFun(LFS, s1, s2, maxlen))^2

 
  
#' Selectivity at length function
#'
#' @param i index 
#' @param SL0.05 length at 5 percent selection
#' @param SL1 length at full selection
#' @param MaxSel Maximum selectivity
#' @param maxlens maximum length 
#' @param Lens vector of lengths 
#' @keywords internal
#' @export SelectFun
SelectFun <- function(i, SL0.05, SL1, MaxSel, maxlens, Lens) {
  s1 <- optimise(getSlope1, interval = c(0, 1e+06), LFS = SL1[i], L0.05 = SL0.05[i])$minimum
  s2 <- optimise(getSlope2, interval = c(0, 1e+06), LFS = SL1[i], s1 = s1, maxlen = maxlens[i], 
    MaxSel = MaxSel[i])$minimum
  if (is.vector(Lens)) TwoSidedFun(LFS = SL1[i], s1 = s1, s2 = s2, lens = Lens)  #nsim = 1
  else TwoSidedFun(LFS = SL1[i], s1 = s1, s2 = s2, lens = Lens[i, ])  #nsim > 1
}


#' Calculate slope from ageM and age95 
#'
#' @param X index 
#' @param ageM age at maturity
#' @param age95 age at 95 percent maturity
#' @keywords internal
#' @export getroot
getroot <- function(X, ageM, age95) 
  uniroot(getSlopeFun, interval = c(1e-04, 5), age50 = ageM[X], age95 = age95[X])$root
  
#' Internal function to calculate slope
#'
#' @param SD standard deviation 
#' @param age50 age at maturity
#' @param age95 age at 95 percent maturity
#' @keywords internal
#' @export getSlopeFun
getSlopeFun <- function(SD, age50, age95) 
  0.95 - (1/(1 + exp((age50 - age95)/(age50 * SD))))
  
  
# Selectivity at length function for GTG model
SelectFunGTG <- function(i, SL0.05, SL1, MaxSel, Linfs, LenGTG) {
    s1 <- optimise(getSlope1, interval = c(0, 100), L1 = SL1[i], L0.05 = SL0.05[i])$minimum
    s2 <- optimise(getSlope2, interval = c(0, 1000), L1 = SL1[i], s1 = s1, 
        Linf = Linfs[i], MaxSel = MaxSel[i])$minimum
    NGTG <- dim(LenGTG)[1]
    t(sapply(1:NGTG, function(X) TwoSidedFun(SL1[i], s1, s2, LenGTG[X, i, ])))
}

# Obj value for opt routine
FitSelect <- function(Pars, V, Linf, Lens) {
    SL0.05 <- (Pars[1])
    SL1 <- (Pars[2])
    MaxSel <- (Pars[3])
    Lens <- t(as.matrix(Lens))
    SS <- sum((V - SelectFun(1, SL0.05, SL1, MaxSel, Linf, Lens))^2)
    return(SS)
}


# Growth-Type-Group Functions (not currently used - Jan 2016)
GenLenFun <- function(NatAGTG, LenatAgeGTG, LenBin, LenMid) {
    Nbins <- length(LenMid)
    SizeComp <- rep(0, Nbins)
    for (L in 1:length(LenMid)) {
        temp <- NatAGTG
        ind <- LenatAgeGTG <= LenBin[L + 1] & LenatAgeGTG > LenBin[L]
        temp[!ind] <- 0
        SizeComp[L] <- SizeComp[L] + sum(temp)
    }
    return(SizeComp)
}






