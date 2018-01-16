
## General Supporting Functions ####

iVB <- function(t0, K, Linf, L) {
  max(1, ((-log(1 - L/Linf))/K + t0))  # Inverse Von-B
}

getr <- function(x, Data, Mvec, Kvec, Linfvec, t0vec, hvec, maxage,
                 r_reps = 100) {
  r <- rep(NA, r_reps)
  for (i in 1:r_reps) {
    log.r = log(0.3)

    opt = optimize(demofn, lower = log(1e-04), upper = log(1.4), M = Mvec[i],
                   amat = iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],
                              Data@L50[x]), sigma = 0.2, K = Kvec[i], Linf = Linfvec[i],
                   to = t0vec[i], hR = hvec[i], maxage = maxage, a = Data@wla[x],
                   b = Data@wlb[x])
    # demographic2(opt$minimum,M[x],ageM[x],0.2,K[x],Linf,t0,steepness[x],maxage,a,b)$r
    r[i] <- exp(opt$minimum)
  }
  r
}

#' TAC Filter
#'
#' Filters vector of TAC recommendations by replacing negatives with NA and
#' and values beyond five standard deviations from the mean as NA
#'
#' @param TAC A numeric vector of TAC recommendations
#' @author T. Carruthers
#' @export
TACfilter <- function(TAC) {
  TAC[TAC < 0] <- NA  # Have to robustify due to R optmization problems.. work in progress.
  TAC[TAC > (mean(TAC, na.rm = T) + 5 * stats::sd(TAC, na.rm = T))] <- NA  # remove very large TAC samples
  return(as.numeric(TAC))
}

## Catch curve function ####
CC <- function(x, Data, reps = 100) {
  ny <- dim(Data@CAA)[2]
  CAA <- apply(Data@CAA[x, max(ny - 2, 1):ny, ], 2, sum)  # takes last two years as the sample (or last year if there is only one)
  maxageobs <- length(CAA)
  AFS <- which.max(CAA)
  AFS[AFS > (maxageobs - 3)] <- maxageobs - 3  # provides at least three datapoints

  nS <- ceiling(sum(CAA)/2)
  y <- log(CAA[AFS:maxageobs]/sum(CAA[AFS:maxageobs], na.rm = T))
  xc <- 1:length(y)
  y[y == "-Inf"] <- NA
  mod <- lm(y ~ xc)
  chk <- sum(is.na(coef(mod)))  # check if model failed
  if (chk) {
    return(NA)
  } else {
    coefs <- summary(mod, weights = CAA[AFS:maxageobs])$coefficients[2, 1:2]
    coefs[is.nan(coefs)] <- tiny
    return(-rnorm(reps, coefs[1], coefs[2]))
  }
}


## DBSRA supporting functions ####

fn <- function(n, BMSY_K) {
  # optimizer to find parameter n according to sampled BMSY/B0 (theta)
  thetapred <- n^(-1/(n - 1))
  (BMSY_K - thetapred)^2
}

getn <- function(BMSY_K) {
  # wrapper for n finder
  optimize(fn, c(0.01, 6), BMSY_K = BMSY_K)$minimum  #get the optimum
}

gety <- function(n) (n^(n/(n - 1)))/(n - 1)  # More DBSRA code: get the y parameter for n

prodPTF <- function(depletion, n, MSY) {
  # Pella-Tomlinson production function required for DB-SRA
  y <- (n^(n/(n - 1)))/(n - 1)
  MSY * y * depletion - MSY * y * depletion^n
}

DBSRAopt <- function(lnK, C_hist, nys, Mdb, FMSY_M, BMSY_K, Bt_K, adelay) {
  # the optimization for B0 given DBSRA assumptions
  Kc <- exp(lnK)
  n <- getn(BMSY_K)
  g <- gety(n)
  FMSY <- FMSY_M * Mdb
  UMSY <- (FMSY/(FMSY + Mdb)) * (1 - exp(-(FMSY + Mdb)))
  MSY <- Kc * BMSY_K * UMSY
  # Bjoin rules from Dick & MacCall 2011
  Bjoin_K <- 0.5
  if (BMSY_K < 0.3)
    Bjoin_K <- 0.5 * BMSY_K
  if (BMSY_K > 0.3 & BMSY_K < 0.5)
    Bjoin_K <- 0.75 * BMSY_K - 0.075
  Bjoin <- Bjoin_K * Kc
  PBjoin <- prodPTF(Bjoin_K, n, MSY)
  cp <- (1 - n) * g * MSY * (Bjoin^(n - 2)) * Kc^-n
  Bc <- rep(NA, nys)
  Bc[1] <- Kc
  obj <- 0
  for (yr in 2:nys) {
    yref <- max(1, yr - adelay)
    if (Bc[yref] > Bjoin | BMSY_K > 0.5) {
      Bc[yr] <- Bc[yr - 1] + g * MSY * (Bc[yref]/Kc) - g * MSY *
        (Bc[yref]/Kc)^n - C_hist[yr - 1]
    } else {
      Bc[yr] <- Bc[yr - 1] + Bc[yref] * ((PBjoin/Bjoin) + cp * (Bc[yref] -
                                                                  Bjoin)) - C_hist[yr - 1]
    }
    if (Bc[yr] < 0)
      obj <- obj + log(-Bc[yr])
    Bc[yr] <- max(1e-06, Bc[yr])
  }
  obj + ((Bc[nys]/Kc) - Bt_K)^2
}  # end of DBSRA optimization function


## Delay-Difference supporting functions ####
DD_R <- function(params, opty, So_DD, Alpha_DD, Rho_DD, ny_DD, k_DD, wa_DD, E_hist,
                 C_hist, UMSYprior) {
  UMSY_DD = exp(params[1])
  MSY_DD = exp(params[2])
  q_DD = exp(params[3])
  SS_DD = So_DD * (1 - UMSY_DD)  # Initialise for UMSY, MSY and q leading.
  Spr_DD = (SS_DD * Alpha_DD/(1 - SS_DD) + wa_DD)/(1 - Rho_DD * SS_DD)
  DsprDu_DD = (Alpha_DD + Spr_DD * (1 + Rho_DD - 2 * Rho_DD * SS_DD))/((1 - Rho_DD * SS_DD) * (1 - SS_DD))
  DsprDu_DD = DsprDu_DD + Alpha_DD * SS_DD/((1 - Rho_DD * SS_DD) * (1 - SS_DD)^2) - Spr_DD/(1 - SS_DD)
  DsprDu_DD = -So_DD * DsprDu_DD
  Arec_DD = 1/(((1 - UMSY_DD)^2) * (Spr_DD + UMSY_DD * DsprDu_DD))
  Brec_DD = UMSY_DD * (Arec_DD * Spr_DD - 1/(1 - UMSY_DD))/MSY_DD
  Spr0_DD = (So_DD * Alpha_DD/(1 - So_DD) + wa_DD)/(1 - Rho_DD * So_DD)
  Ro_DD = (Arec_DD * Spr0_DD - 1)/(Brec_DD * Spr0_DD)
  Bo_DD = Ro_DD * Spr0_DD
  No_DD = Ro_DD/(1 - So_DD)

  B_DD <- rep(NA, ny_DD + 1)
  N_DD <- rep(NA, ny_DD + 1)
  R_DD <- rep(NA, ny_DD + k_DD)
  Cpred_DD <- rep(NA, ny_DD)

  B_DD[1] = Bo_DD
  N_DD[1] = No_DD
  R_DD[1:k_DD] = Ro_DD

  for (tt in 1:ny_DD) {

    Surv_DD = So_DD * exp(-q_DD * E_hist[tt])
    Cpred_DD[tt] = B_DD[tt] * (1 - exp(-q_DD * E_hist[tt]))
    Sp_DD = B_DD[tt] - Cpred_DD[tt]
    R_DD[tt + k_DD] = Arec_DD * Sp_DD/(1 + Brec_DD * Sp_DD)
    B_DD[tt + 1] = Surv_DD * (Alpha_DD * N_DD[tt] + Rho_DD * B_DD[tt]) +
      wa_DD * R_DD[tt + 1]
    N_DD[tt + 1] = Surv_DD * N_DD[tt] + R_DD[tt + 1]

  }
  tiny <- 1e-15
  Cpred_DD[Cpred_DD < tiny] <- tiny

  if (opty == 1) {
    # The following conditions must be met for positive values
    # of Arec_DD and Brec_DD, respectively:
    # umsy * DsprDu + Spr_DD > 0 and Arec_DD * Spr_DD * (1 - UMSY_DD) - 1 > 0
    # Thus, create a likelihood penalty of 100 if either condition is not met
    umsy_penalty <- ifelse(Spr_DD + UMSY_DD * DsprDu_DD > 0, 0, 100)
    alpha_penalty <- ifelse(Arec_DD * Spr_DD * (1 - UMSY_DD) - 1 > 0, 0, 100)
    
    test <- dnorm(log(Cpred_DD), log(C_hist), 0.25, log = T)
    test2 <- dlnorm(UMSY_DD, log(UMSYprior[1]), UMSYprior[2], log = T)
    test[is.na(test)] <- -1000
    test[test == (-Inf)] <- -1000
    if (is.na(test2) | test2 == -Inf | test2 == Inf)
      test2 <- 1000
    return(-sum(test, test2) + umsy_penalty + alpha_penalty)  # return objective function
  } else if (opty == 2) {
    # return MLE TAC estimate
    UMSY_DD * B_DD[ny_DD]
  } else if (opty == 3) {
    B_DD[tt + 1]/Bo_DD
  } else {
    cbind(C_hist, Cpred_DD)  # return observations vs predictions
  }
}

## Mean Length supporting functions ####

bheq <- function(K, Linf, Lc, Lbar) {
  K * (Linf - Lbar)/(Lbar - Lc)
}

bhnoneq <- function(year, mlen, ss, K, Linf, Lc, nbreaks, styrs, stZ) {
  mlen[mlen <= 0 | is.na(mlen)] <- -99
  ss[ss <= 0 | is.na(ss) | mlen == -99] <- 0
  stpar <- c(stZ, styrs)

  # results <-
  # optim(stpar,bhnoneq_LL,method='BFGS',year=year,Lbar=mlen,ss=ss,
  # nbreaks=nbreaks,K=K,Linf=Linf,Lc=Lc,control=list(maxit=1e6))
  results <- try(optim(stpar, bhnoneq_LL, method = "Nelder-Mead", year = year,
                       Lbar = mlen, ss = ss, nbreaks = nbreaks, K = K, Linf = Linf, Lc = Lc,
                       control = list(maxit = 1e+06), hessian = FALSE), silent=TRUE)
  if (class(results) == "try-error") {
    return(FALSE)
  } else return(results$par[1:(nbreaks + 1)])
}



MLne <- function(x, Data, Linfc, Kc, ML_reps = 100, MLtype = "dep") {
  year <- 1:dim(Data@CAL)[2]
  nlbin <- ncol(Data@CAL[x, , ])
  nlyr <- nrow(Data@CAL[x, , ])
  mlbin <- (Data@CAL_bins[1:nlbin] + Data@CAL_bins[2:(nlbin + 1)])/2
  nbreaks <- 1
  Z <- matrix(NA, nrow = ML_reps, ncol = nbreaks + 1)
  Z2 <- rep(NA, ML_reps)
  # temp<-apply(Data@CAL[x,,],2,sum) Lc<-mlbin[which.max(temp)] #
  # modal length

  # dd <- dim(Data@CAL[x,,]) curLen <- Data@CAL[x,dd[1],] Lc <-
  # mlbin[which.max(curLen)] Lc <- Data@LFS[,x] Lc <- Lc[length(Lc)]
  Lc <- Data@LFS[x]

  for (i in 1:ML_reps) {
    mlen <- rep(NA, length(year))
    ss <- ceiling(apply(Data@CAL[x, , ], 1, sum)/2)
    if (MLtype == "dep") {
      for (y in 1:length(year)) {
        if (sum(Data@CAL[x, y, ] > 0) > 0.25 * length(Data@CAL[x, y, ])) {
          temp2 <- sample(mlbin, ceiling(sum(Data@CAL[x, y, ])/2), replace = T, prob = Data@CAL[x, y, ])
          mlen[y] <- mean(temp2[temp2 >= Lc], na.rm = TRUE)
        }
      }

      fitmod <- bhnoneq(year = year, mlen = mlen, ss = ss, K = Kc[i], Linf = Linfc[i],
                        Lc = Lc, nbreaks = nbreaks, styrs = ceiling(length(year) * ((1:nbreaks)/(nbreaks + 1))),
                        stZ = rep(Data@Mort[x], nbreaks + 1))
      if (all(fitmod == FALSE)) {
        Z[i, ] <- NA
      } else Z[i, ] <- fitmod
    } else {

      # ind<-(which.min(((Data@CAL_bins-Data@LFS[x])^2)^0.5)-1):(length(Data@CAL_bins)-1)
      for (y in 1:length(year)) {
        if (sum(Data@CAL[x, y, ] > 0) > 0.25 * length(Data@CAL[x, y, ])) {
          temp2 <- sample(mlbin, ceiling(sum(Data@CAL[x, y, ])/2), replace = T, prob = Data@CAL[x, y, ])
          mlen[y] <- mean(temp2[temp2 >= Lc], na.rm = TRUE)
        }
      }
      mlen <- mean(mlen[(length(mlen) - 2):length(mlen)], na.rm = TRUE)
      Z2[i] <- bheq(K = Kc[i], Linf = Linfc[i], Lc = Lc, Lbar = mlen)

    }
  }
  # Z <- Z[,ncol(Z)] # last estimate of Z? Z needs to be vector reps long
  if (MLtype == "F") {
    Z2[Z2<0] <- NA
    return(Z2)
  }
  if (MLtype == "dep") return(Z)
}


## SP supporting functions ####
SPSRAopt <- function(lnK, dep, r, Ct, PE) {
  nyears <- length(Ct)
  B <- rep(NA, nyears)
  B[1] <- exp(lnK)
  OBJ <- 0
  for (y in 2:nyears) {
    if ((B[y - 1] - Ct[y - 1]) < 0)
      OBJ <- OBJ + (B[y - 1] - Ct[y - 1])^2
    B[y] <- max(0.01, B[y - 1] - Ct[y - 1])
    B[y] <- B[y] + r * B[y] * (1 - B[y]/B[1]) * PE[y]
  }
  return(OBJ + ((B[nyears]/B[1]) - dep)^2)
}


## SRA supporting functions ####
SRAfunc <- function(lnR0c, Mc, hc, maxage, LFSc, LFCc, Linfc, Kc, t0c,
                    AMc, ac, bc, Catch, CAA, opt = 1) {

  ny <- length(Catch)
  AFC <- log(1 - min(0.99, LFCc/Linfc))/-Kc + t0c
  AFS <- log(1 - min(0.99, LFSc/Linfc))/-Kc + t0c
  if (AFC >= 0.7 * maxage)
    AFC <- 0.7 * maxage
  if (AFS >= 0.9 * maxage)
    AFS <- 0.9 * maxage
  KES <- max(2, ceiling(mean(c(AFC, AFS))))
  vul <- rep(1, maxage)
  vul[1:(KES - 1)] <- 0
  Mac <- rep(1, maxage)
  Mac[1:max(1, floor(AMc))] <- 0
  Lac <- Linfc * (1 - exp(-Kc * ((1:maxage) - t0c)))
  Wac <- ac * Lac^bc
  R0c <- exp(lnR0c)
  N <- exp(-Mc * ((1:maxage) - 1)) * R0c
  SSN <- Mac * N  # Calculate initial spawning stock numbers
  Biomass <- N * Wac
  SSB <- SSN * Wac  # Calculate spawning stock biomass

  B0 <- sum(Biomass)
  SSB0 <- sum(SSB)
  SSN0 <- SSN
  SSBpR <- sum(SSB)/R0c  # Calculate spawning stock biomass per recruit
  SSNpR <- SSN/R0c

  CN <- array(NA, dim = c(ny, maxage))
  HR <- rep(0, maxage)
  pen <- 0
  for (y in 1:ny) {
    # set up some indices for indexed calculation
    VB <- Biomass[KES:maxage] * exp(-Mc)
    CB <- Catch[y] * VB/sum(VB)
    testHR <- CB[1]/VB[1]
    if (testHR > 0.8)
      pen <- pen + (testHR - 0.8)^2
    HR[KES:maxage] <- min(testHR, 0.8)
    FMc <- -log(1 - HR)  # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
    Zc <- FMc + Mc

    CN[y, ] <- N * (1 - exp(-Zc)) * (FMc/Zc)
    N[2:maxage] <- N[1:(maxage - 1)] * exp(-Zc[1:(maxage - 1)])  # Total mortality
    N[1] <- (0.8 * R0c * hc * sum(SSB))/(0.2 * SSBpR * R0c * (1 - hc) +
                                           (hc - 0.2) * sum(SSB))  # Recruitment assuming regional R0 and stock wide steepness
    # print(N[1])
    Biomass <- N * Wac
    SSN <- N * Mac
    SSB <- SSN * Wac

  }  # end of year

  CN[CN < 0] <- 0  # stop any negative catches
  syear <- ny - dim(CAA)[1] + 1
  pred <- CN[syear:ny, ]
  pred <- pred/array(apply(pred, 1, sum), dim = c(dim(CAA)[1], maxage))

  fobj <- pen - sum(log(pred + tiny) * CAA, na.rm = T)
  if (opt == 1) {
    return(fobj)
  } else if (opt == 2) {
    return(sum(Biomass))
  } else if (opt == 3) {
    sum(SSB)/sum(SSB0)
  }
  # CBc<-sum(CB)
}

SRAFMSY <- function(lnFMc, Mc, hc, maxage, LFSc, LFCc, Linfc, Kc, t0c,
                    AMc, ac, bc, opt = T) {

  FMc <- exp(lnFMc)
  ny <- 100
  AFC <- log(1 - min(0.99, LFCc/Linfc))/-Kc + t0c
  AFS <- log(1 - min(0.99, LFSc/Linfc))/-Kc + t0c
  if (AFC >= 0.7 * maxage)
    AFC <- 0.7 * maxage
  if (AFS >= 0.9 * maxage)
    AFS <- 0.9 * maxage

  KES <- max(2, ceiling(mean(c(AFC, AFS))))
  vul <- rep(1, maxage)
  vul[1:(KES - 1)] <- 0
  Mac <- rep(1, maxage)
  Mac[1:max(1, floor(AMc))] <- 0
  Lac <- Linfc * (1 - exp(-Kc * ((1:maxage) - t0c)))
  Wac <- ac * Lac^bc
  R0c <- 1
  N <- exp(-Mc * ((1:maxage) - 1)) * R0c
  SSN <- Mac * N  # Calculate initial spawning stock numbers
  Biomass <- N * Wac
  SSB <- SSN * Wac  # Calculate spawning stock biomass

  B0 <- sum(Biomass)
  SSB0 <- sum(SSB)
  SSN0 <- SSN
  SSBpR <- sum(SSB)/R0c  # Calculate spawning stock biomass per recruit
  SSNpR <- SSN/R0c

  N <- N/2
  SSN <- Mac * N  # Calculate initial spawning stock numbers
  Biomass <- N * Wac
  SSB <- SSN * Wac

  for (y in 1:ny) {
    # set up some indices for indexed calculation Fishing mortality rate
    # determined by effort, catchability, vulnerability and spatial
    # preference according to biomass
    Zc <- FMc * vul + Mc
    CN <- N * (1 - exp(-Zc)) * (FMc/Zc)
    CB <- CN * Wac
    Biomass <- N * Wac
    N[2:maxage] <- N[1:(maxage - 1)] * exp(-Zc[1:(maxage - 1)])  # Total mortality
    N[1] <- (0.8 * R0c * hc * sum(SSB))/(0.2 * SSBpR * R0c * (1 - hc) +
                                           (hc - 0.2) * sum(SSB))  # Recruitment assuming regional R0 and stock wide steepness
    # print(N[1])
    SSN <- N * Mac
    SSB <- SSN * Wac

  }  # end of year

  if (opt) {
    return(-sum(CB))
  } else {
    return(FMc)
  }
}


## Sampling steepness parameter (h)  ####
#
# Code from Quang Huynh that fixes the bug where h is sometimes sampled > 1 or < 0.2

#' Sample steepness given mean and cv
#'
#' @param n number of samples
#' @param mu mean h
#' @param cv cv of h
#'
#' @author Q. Huynh
#'
sample_steepness2 <- function(n, mu, cv) {

  if(n == 1) return(mu)
  else {
    sigma <- mu * cv
    
    mu.beta.dist <- (mu - 0.2)/0.8
    sigma.beta.dist <- sigma/0.8
    
    beta.par <- derive_beta_par(mu.beta.dist, sigma.beta.dist)
    
    h.transformed <- rbeta(n, beta.par[1], beta.par[2])
    
    h <- 0.8 * h.transformed + 0.2
    h[h > 0.99] <- 0.99
    h[h < 0.2] <- 0.2
    
    return(h)
  }
  
}


#' This function reduces the CV by 5 per cent until steepness values can be sampled without error
#'
#'
#' @param mu mean h
#' @param sigma sd of h
#'
#' @author Q. Huynh
#'
derive_beta_par <- function(mu, sigma) {

  a <- alphaconv(mu, sigma)
  b <- betaconv(mu, sigma)

  if(a <= 0 || b <= 0) {
    sigma <- 0.95 * sigma
    Recall(mu, sigma)
  }
  else return(c(a, b))

}


## VPA supporting functions ####
VPAopt = function(theta, Cat, yt, S, maxage, wa, pmat, opt = T, usewat = F) {

  Uterm <- exp(theta[1])/(1 + exp(theta[1]))
  minagecom = 1
  minagesel = ceiling(maxage * 0.66)
  avg_yrsel <- max(2, (min(5, floor(dim(Cat)[1]/2))))

  sig = exp(theta[2])
  tiny = 1e-10
  n = dim(Cat)[1]
  A = dim(Cat)[2]
  Nat = matrix(NA, n, A)  # Numbers-at-age matrix
  Ut = rep(NA, length = n)
  ai = 1:(A - 2)
  va = c(plogis(1:(A - 4), 2, 0.2), rep(1, 4))  # Initial values at the terminal selectivy
  Ut[n] = Uterm

  for (j in 1:15) {
    # Numerical convergence to terminal F print(Ut)
    Nat[n, ] = Cat[n, ]/(Uterm * va)  # Initialize the terminal year

    for (i in (n - 1):1) {
      Nat[i, ai] = Nat[i + 1, ai + 1]/S + Cat[i, ai]
      Nat[i, A - 1] = (Nat[i + 1, A]/S + Cat[i, A - 1] + Cat[i, A]) *
        (Cat[i, A - 1]/(Cat[i, A - 1] + Cat[i, A] + tiny))
      Nat[i, A] = (Nat[i + 1, A]/S + Cat[i, A - 1] + Cat[i, A]) *
        (Cat[i, A]/(Cat[i, A - 1] + Cat[i, A] + tiny))
    }
    # modify this parameters if need it ##

    # minagesel = 8
    minagecom = 1

    Ut = rowSums(Cat[, minagesel:(A - 1)])/rowSums(Nat[, minagesel:(A -
                                                                      1)])  # Exploitation rate for fully recruited fish
    # Ut[n] = 0.4
    vat = Cat/Nat/Ut  # relative vulnerablility at age
    va = colMeans(vat[(n - avg_yrsel):(n - minagecom), ])  # update terminal vul
    va[minagesel:A] = 1

  }

  vat[is.na(vat)] <- 1
  Ut[n] = Uterm
  if (usewat == T)
    vbt = rowSums((Nat * vat) * wa) else vbt = (Nat * vat) %*% wa
  if (usewat == T)
    bt = as.vector(rowSums(Nat * wa)) else bt = as.vector(Nat %*% wa)
  fec = pmat * wa
  zt = log(yt/vbt)
  epsilon = zt - mean(zt)
  if (usewat == T)
    ssb = as.vector(rowSums(Nat * fec)) else ssb = as.vector(Nat %*% fec)
  predcpue = exp(mean(zt)) * vbt  ### check again if bt or vbt
  cpue_q = yt/exp(mean(zt))
  qhat = exp(mean(epsilon))

  lnl = sum(dnorm(epsilon, mean = 0, sd = sig, log = T))

  if (opt) {
    return(-lnl)
  } else {
    # ss = sum(epsilon^2) lnl = 0.5*n*log(ss)
    return(list(Uterm = Uterm, va = va, rt = Nat[, 1], ssb = ssb, yt = yt,
                vbt = vbt, cpue_q = cpue_q, Nat = Nat, vat = vat, Ut = Ut,
                bt = bt, predcpue = predcpue, epsilon = epsilon/sig, lnl = lnl,
                qhat = qhat, minagesel = minagesel, minagecom = minagecom,
                avg_yrsel = avg_yrsel))
  }


}

VPAFMSY <- function(lnFMc, Mc, hc, maxage, vul, Linfc, Kc, t0c, AMc, ac,
                    bc, opt = T, ny = 50) {

  FMc <- exp(lnFMc)

  Mac <- rep(1, maxage)
  Mac[1:max(1, floor(AMc))] <- 0
  Lac <- Linfc * (1 - exp(-Kc * ((1:maxage) - t0c)))
  Wac <- ac * Lac^bc
  R0c <- 1
  N <- exp(-Mc * ((1:maxage) - 1)) * R0c
  SSN <- Mac * N  # Calculate initial spawning stock numbers
  Biomass <- N * Wac
  SSB <- SSN * Wac  # Calculate spawning stock biomass

  B0 <- sum(Biomass)
  SSB0 <- sum(SSB)
  SSN0 <- SSN
  SSBpR <- sum(SSB)/R0c  # Calculate spawning stock biomass per recruit
  SSNpR <- SSN/R0c

  N <- N/2
  SSN <- Mac * N  # Calculate initial spawning stock numbers
  Biomass <- N * Wac
  SSB <- SSN * Wac

  for (y in 1:ny) {
    # set up some indices for indexed calculation Fishing mortality rate
    # determined by effort, catchability, vulnerability and spatial
    # preference according to biomass
    Zc <- FMc * vul + Mc
    CN <- N * (1 - exp(-Zc)) * (FMc/Zc)
    CB <- CN * Wac
    Biomass <- N * Wac
    N[2:maxage] <- N[1:(maxage - 1)] * exp(-Zc[1:(maxage - 1)])  # Total mortality
    N[1] <- (0.8 * R0c * hc * sum(SSB))/(0.2 * SSBpR * R0c * (1 - hc) +
                                           (hc - 0.2) * sum(SSB))  # Recruitment assuming regional R0 and stock wide steepness
    # print(N[1])
    SSN <- N * Mac
    SSB <- SSN * Wac

  }  # end of year

  if (opt) {
    return(-sum(CB))
  } else {
    return(FMc)
  }
}


## YPR supporting functions ####
# Yield per recruit estimate of FMSY Meaghan Bryan 2013
YPRopt = function(Linfc, Kc, t0c, Mdb, a, b, LFS, maxage, reps = 100) {

  nf <- 200
  frates <- seq(0, 3, length.out = nf)
  Winf = a * Linfc^b
  rat <- LFS/Linfc
  rat[rat > 0.8] <- 0.8  # need to robustify this for occasionally very high samples of LFS
  tc = log(1 - rat)/-Kc + t0c
  tc = round(tc, 0)
  tc[tc < 1] <- 1
  tc[tc > maxage] <- maxage

  vul <- array(0, dim = c(reps, maxage))
  mat <- array(0, dim = c(reps, maxage))
  lx <- array(NA, dim = c(reps, maxage))
  lxo <- array(NA, dim = c(reps, maxage))

  ypr <- array(NA, dim = c(reps, nf))
  sbpr <- array(NA, dim = c(reps, nf))
  sbpr.ratio <- array(NA, dim = c(reps, nf))
  sbpr.dif <- array(NA, dim = c(reps, nf))

  f.max <- array(NA, dim = c(reps, maxage))

  # average weight at age - follow von Bertalanffy growth
  age <- array(rep(1:maxage, each = reps), dim = c(reps, maxage))
  la <- Linfc * (1 - exp(-Kc * ((age - t0c))))
  wa <- a * la^b

  # vulnerability schedule - assumes knife-edge vulnerability, where all
  # individuals age tc to maxage are fully vulnerbale all individulas
  # less than age tc are not vulnerable
  for (i in 1:reps) {
    if (tc[i] > 0)
      vul[i, tc[i]:maxage] <- 1
    if (tc[i] > 1)
      mat[i, max(1, tc[i] - 1):maxage] <- 1
  }

  lx[, 1] <- 1
  lxo[, 1] <- 1
  for (k in 1:nf) {
    for (i in 2:maxage) {
      lx[, i] = lx[, i - 1] * exp(-(Mdb + vul[, i - 1] * frates[k]))
      lxo[, i] = lx[, i] * exp(-Mdb)
    }
    phi_vb = apply(lx * wa * vul, 1, sum)
    sbpro = apply(lxo * wa * mat, 1, sum)

    ypr[, k] = (1 - exp(-frates[k])) * phi_vb
    sbpr[, k] = apply(lx * wa * mat, 1, sum)
    sbpr.ratio[, k] = sbpr[, k]/sbpro
    sbpr.dif[, k] = abs(sbpr.ratio[, k] - 0.3)  #hard code comparison ratio
  }

  # frates[apply(ypr,1,which.max)] Fmaxypr

  # More code that derived F0.1 in 'per recruit analysis.R' (Meaghan
  # Bryan)
  slope.origin = (ypr[, 2] - ypr[, 1])/(frates[2] - frates[1])
  slope.10 = round(0.1 * slope.origin, 2)

  slope = array(NA, dim = dim(ypr))  #vector(length=length(ypr))
  slope[, 1] = slope.origin
  for (i in 3:ncol(ypr)) {
    slope[, i - 1] = round((ypr[, i] - ypr[, i - 1])/(frates[i] - frates[i - 1]), 2)
  }
  dif = abs(slope - slope.10)
  dif[is.na(dif)] <- 1e+11
  frates[apply(dif, 1, which.min)]  #frates[which.min(dif)]
}

