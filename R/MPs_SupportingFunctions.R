##### MP plotting #####

plotAvC <- function(x, Data, meanC, histCatch, yr.ind, lwd=3, cex.lab=1.25) {
  op <- par(no.readonly = TRUE)
  on.exit(op)
  par(mfrow=c(1,1))
  plot(c(Data@Year[yr.ind], Data@Year[max(yr.ind)]+1), c(histCatch,NA), type="l", 
       xlab="Year", ylab=paste0("Catch (", Data@Units, ")"), lwd=lwd, bty="l", las=1, cex.lab=cex.lab)
  abline(v=Data@LHYear, lty=2, col="darkgray") #
  text(Data@LHYear, max(histCatch, na.rm=TRUE)*0.9, "Last Historical Year", pos=2, xpd=NA)
  lines(c(min(Data@Year), Data@LHYear), rep(mean(Data@Cat[x,yr.ind]),2), lty=2) #
  text(quantile(Data@Year, 0.1), meanC*1.1, pos=4, "Average Historical Catch")
  boxplot(Rec@TAC, add=TRUE, at=max(Data@Year)+1, col="grey", width=1, outline=TRUE, axes=FALSE)
  text(max(Data@Year)+1, quantile(Rec@TAC, 0.05), "TAC", col="black", pos=2)
}

plotDCAC <- function(x, Data, dcac, yrs, lwd=3, cex.lab=1.25) {
  op <- par(no.readonly = TRUE)
  on.exit(op)
  par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(5,4,1,4))
  yr.lst <- max(yrs)
  ylim <- c(0, max(c(Data@Cat[x,1:yr.lst], dcac)))
  plot(c(Data@Year[yrs], Data@Year[max(yrs)]+1:3), c(Data@Cat[x,1:yr.lst],NA, NA, NA), type="l", 
       xlab="Year", ylab=paste0("Catch (", Data@Units, ")"), lwd=lwd, bty="l", las=1, cex.lab=cex.lab,
       ylim=ylim)
  abline(v=Data@LHYear, lty=2, col="darkgray") #
  
  text(Data@LHYear, max(Data@Cat[x,1:yr.lst], na.rm=TRUE)*0.9, "Last Historical Year", pos=2, xpd=NA)
  lines(c(min(Data@Year), Data@LHYear), rep(mean(Data@Cat[x,1:yr.lst]),2), lty=2) #
  text(quantile(Data@Year, 0.1), mean(Data@Cat[x,1:yr.lst])*1.1, pos=4, "Average Historical Catch")
}

plotDCACadd <- function(TAC, Data, Bt_K) {
  boxplot(TAC, add=TRUE, at=max(Data@Year)+1, col="darkgrey", width=1, outline=TRUE, axes=FALSE)
  text(max(Data@Year)+1, quantile(TAC, 0.05), "TAC", col="black", pos=2)
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(new = T)
  plot(c(1, max(Data@Year)+3), c(0,1), type="n", axes=FALSE, xlab="", ylab="")
  quants <- quantile(Bt_K, c(0.025, 0.5, 0.975))
  points(max(Data@Year)+3, quants[2], pch=16, col="blue", cex=1.5)
  lines(c(max(Data@Year)+3, max(Data@Year)+3), c(quants[1], quants[3]), col="blue")
  axis(side=4, las=1, col="blue", labels=FALSE)
  at = axTicks(4)
  mtext(side = 4, text = at, at = at, col = "blue", line = 1, las=1)
  mtext(side=4, "Depletion (median + 95 percentiles)", line=3, cex=1.25, col="blue")
}

plotBK <- function(DF) {
  DF2 <- DF %>% filter(vars %in% c("Lc/Linf", "K", "Fmax"))
  p1 <- ggplot(DF2, aes(x=vars, y=vals)) + geom_boxplot() + 
    theme_classic() + expand_limits(y=0) + labs(x="", y='Values')
  DF3 <- DF %>% filter(!vars %in% c("Lc/Linf", "K", "Fmax"))
  p2 <- ggplot(DF3, aes(x=vars, y=vals)) + geom_boxplot() + 
    theme_classic() + expand_limits(y=0) + labs(x="", y='Values')
  
  gridExtra::grid.arrange(p1, p2, nrow=2)
}


plotCompSRA <- function(runCompSRA, TAC) {
  op <- par(no.readonly = TRUE)
  on.exit(op)
  
  CAA <-runCompSRA$CAA
  CAA <- CAA/apply(CAA, 1, sum)
  nsamps <- nrow(CAA)
  ages <- 1:ncol(CAA)
  nreps <- length(runCompSRA$pred)
  
  nplots <- nsamps + 2 
  
  ncol <- ceiling(sqrt(nplots))
  nrow <- ceiling(nplots/ncol)
  par(mfrow=c(nrow, ncol), oma=c(2,2,3,2))
  
  for (x in 1:nsamps) {
    ylim <- c(0, max(CAA[x,], max(unlist(runCompSRA$pred))))
    plot(ages, CAA[x,], type="l", lwd=3, bty="n", xlab="Age", ylab="Frequency", ylim=ylim)
    for (r in 1:nreps) matplot(ages, runCompSRA$pred[[r]][x,], add=TRUE, type="l")
  }
  mtext("Catch-at-age (+ fitted)", side=3, outer=TRUE)
  
  ylim <- c(0, max(c(runCompSRA$Bt_K, runCompSRA$FMSY)))
  boxplot(runCompSRA$Bt_K, runCompSRA$FMSY, ylim=ylim, las=1, names=c("Depletion", "FMSY"))
  
  ylim <- c(0, max(c(runCompSRA$Ac, TAC)))
  boxplot(runCompSRA$Ac, TAC, ylim=ylim, las=1, names=c("Abundance", "TAC"))
}
# default plotting options
leg.pos <- col1 <- col2 <- col3 <- col4 <- pt.cex <- tex.cex <- cex.lab <- lwd <- leg.post <- NULL
MP.plot <- new.env()
MP.plot$col1 <- "#2B2E6E"
MP.plot$col2 <- "#9F812D"
MP.plot$col3 <- "black"
MP.plot$col4 <- "darkgray"
MP.plot$pt.cex <- 2
MP.plot$tex.cex <- 1.1
MP.plot$cex.lab <- 1.25
MP.plot$lwd <-3
MP.plot$leg.pos <- "topleft"







Itarget_p <- function(...) {
 
  inlist <- list(...)

  # default plotting options
  for (nm in names(MP.plot)) {
    if (!nm %in% names(inlist)) {
      assign(nm, MP.plot[[nm]])
    } else {
      assign(nm, inlist[[nm]])
    }
  }


  ylim <- c(0, max(c(1, I0, Itarget, Data@Ind[x, ])))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow=c(1,2))
  plot(Data@Year, Data@Ind[x, ], ylim=ylim, type='l', bty="l", 
       xlab="Year", ylab="Index of Total Abundance", las=1, cex.lab=cex.lab) 
  abline(v=Data@LHYear, lty=3, col=col4)
  text(Data@LHYear, min(ylim)*1.15, "Last Historical Year", col=col4, pos=2)
  
  lines(Data@Year[ind3], rep(Iave, length(ind3)), col=col2, lwd=lwd)
  lines(Data@Year[ind], rep(I0, length(Data@Year[ind])), col=col2, lwd=lwd-1)
  
  lines(Data@Year[ind], rep(Irecent, length(ind)), col=col1, lwd=lwd)
  points(max(Data@Year[ind]), Itarget, col=col3, cex=pt.cex, pch=16)
  

  legend(leg.pos, lty=c(1,1,1,NA), pch=c(NA, NA, NA, 16), col=c(col1, col2, col2, col3), lwd=c(lwd, lwd, lwd-1, NA),
         legend=c(as.expression(bquote("I"[recent]~(mean~last~ .(yrsmth) ~ "years"))),
                  as.expression(bquote("I"[average]~(mean~last~ .(yrsmth*2) ~ "historical years"))),
                  as.expression(bquote("I"[0]~(0.8~"I"[average]))),
                  as.expression(bquote("I"[target]~(.(Imulti)~"I"[average])))),
         bty="n", pt.cex=pt.cex, cex=tex.cex, xpd=TRUE)
  
  
  
  ylim <- range(c(C_dat,TAC, Data@Cat[x, ]))

  plot(c(Data@Year, max(Data@Year)+1), c(Data@Cat[x, ], NA), ylim=ylim, type="l", bty="l", 
       xlab="Year", ylab="Catch", las=1, cex.lab=cex.lab)
  lines(Data@Year[ind2], C_dat, col=col1, lwd=lwd-1)
  
  abline(v=Data@LHYear, lty=3, col=col4)
  text(Data@LHYear, min(ylim)*1.1, "Last Historical Year", col=col4, pos=2)
  lines(Data@Year[ind2], rep(mean(C_dat), length(ind2)), col=col1, lwd=lwd)
  if ((1-xx) ==1) {
    legend(leg.pos, lty=c(1, NA, NA), pch=c(NA, 16), lwd=lwd, col=c(col1, "black"), 
           legend =c(as.expression(bquote("C"[average]~(last~ .(yrsmth) ~ "historical years"))),
                     "TAC"),
           bty="n", pt.cex=pt.cex, cex=tex.cex, xpd=TRUE) 
  } else {
    points(mean(Data@Year[ind2]), mean(TACstar), cex=pt.cex, pch=16, col=col1)
    legend(leg.pos, lty=c(1, NA, NA), pch=c(NA, 16, 16), lwd=lwd, col=c(col1, col1, "black"), 
           legend =c(as.expression(bquote("C"[average]~(last~ .(yrsmth) ~ "historical years"))),
                     as.expression(bquote("C"[star]~(.(1-xx)~"C"[average]))),
                     "TAC"),
           bty="n", pt.cex=pt.cex, cex=tex.cex, xpd=TRUE) 
  }
  
  points(max(Data@Year)+1, mean(TAC), cex=pt.cex, pch=16)
  if (!all(TAC/mean(TAC) == 1)) {
    boxplot(TAC, at=max(Data@Year)+1, add=TRUE, xpd=NA, axes=FALSE, outline = FALSE)  
  }
}



## General Supporting Functions ####





#' Inverse von Bertalanffy 
#'
#' Calculate the age given length from the vB equation
#' @param t0 Hypothetical age when length is 0
#' @param K Growth coefficient
#' @param Linf Asymptotic length
#' @param L Length
#' @export
#' @keywords internal 
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




## Delay-Difference supporting functions ####
DD_R <- function(params, opty, So_DD, Alpha_DD, Rho_DD, ny_DD, k_DD, wa_DD, E_hist,
                 C_hist, UMSYprior) {
  UMSY_DD = 1/(1 + exp(-params[1])) # Logit transform to constrain u between 0-1
  MSY_DD = exp(params[2])
  q_DD = exp(params[3])
  SS_DD = So_DD * (1 - UMSY_DD)  # Initialise for UMSY, MSY and q leading.
  Spr_DD = (SS_DD * Alpha_DD/(1 - SS_DD) + wa_DD)/(1 - Rho_DD * SS_DD)
  DsprDu_DD = ((Alpha_DD + Spr_DD * (1 + Rho_DD - 2 * Rho_DD * SS_DD))/((1 - Rho_DD * SS_DD) * (1 - SS_DD)) + 
    Alpha_DD * SS_DD/((1 - Rho_DD * SS_DD) * (1 - SS_DD)^2) - Spr_DD/(1 - SS_DD)) * -So_DD
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
    umsy_penalty <- ifelse(Spr_DD + UMSY_DD * DsprDu_DD > 0, 0, UMSY_DD * 100)
    alpha_penalty <- ifelse(Arec_DD * Spr_DD * (1 - UMSY_DD) - 1 > 0, 0, UMSY_DD * 100)
    
    sigma <- sqrt(sum((log(C_hist) - log(Cpred_DD))^2)/ny_DD) # Analytical solution
    
    test <- dnorm(log(C_hist), log(Cpred_DD), sigma, log = T)
    test2 <- dbeta(UMSY_DD, UMSYprior[1], UMSYprior[2], log = T)
    test[is.na(test)] <- -1000
    test[test == (-Inf)] <- -1000
    if (is.na(test2) | test2 == -Inf | test2 == Inf)
      test2 <- 1000
    return(-sum(test, test2) + umsy_penalty + alpha_penalty)  # return objective function
  } else if (opty == 2) {
    # return MLE TAC estimate
    UMSY_DD * B_DD[ny_DD + 1]
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
  if (AFC >= 0.7 * maxage) AFC <- 0.7 * maxage
  if (AFS >= 0.9 * maxage) AFS <- 0.9 * maxage
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
    VB <- Biomass[KES:maxage] * exp(-Mc)
    CB <- Catch[y] * VB/sum(VB)
    testHR <- CB[1]/VB[1]
    if (testHR > 0.8)  pen <- pen + (testHR - 0.8)^2
    HR[KES:maxage] <- min(testHR, 0.8)
    FMc <- -log(1 - HR)  # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
    Zc <- FMc + Mc

    CN[y, ] <- N * (1 - exp(-Zc)) * (FMc/Zc)
    N[2:maxage] <- N[1:(maxage - 1)] * exp(-Zc[1:(maxage - 1)])  # Total mortality
    N[1] <- (0.8 * R0c * hc * sum(SSB))/(0.2 * SSBpR * R0c * (1 - hc) +
                                           (hc - 0.2) * sum(SSB))  # Recruitment assuming regional R0 and stock wide steepness
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
  } 
  if (opt == 2) {
    return(list(B=sum(Biomass), D=sum(SSB)/sum(SSB0), pred=pred))
  } 
 
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
# VPAopt = function(theta, Cat, yt, S, maxage, wa, pmat, opt = T, usewat = F) {
# 
#   Uterm <- exp(theta[1])/(1 + exp(theta[1]))
#   minagecom = 1
#   minagesel = ceiling(maxage * 0.66)
#   avg_yrsel <- max(2, (min(5, floor(dim(Cat)[1]/2))))
# 
#   sig = exp(theta[2])
#   tiny = 1e-10
#   n = dim(Cat)[1]
#   A = dim(Cat)[2]
#   Nat = matrix(NA, n, A)  # Numbers-at-age matrix
#   Ut = rep(NA, length = n)
#   ai = 1:(A - 2)
#   va = c(plogis(1:(A - 4), 2, 0.2), rep(1, 4))  # Initial values at the terminal selectivy
#   Ut[n] = Uterm
# 
#   for (j in 1:15) {
#     # Numerical convergence to terminal F print(Ut)
#     Nat[n, ] = Cat[n, ]/(Uterm * va)  # Initialize the terminal year
# 
#     for (i in (n - 1):1) {
#       Nat[i, ai] = Nat[i + 1, ai + 1]/S + Cat[i, ai]
#       Nat[i, A - 1] = (Nat[i + 1, A]/S + Cat[i, A - 1] + Cat[i, A]) *
#         (Cat[i, A - 1]/(Cat[i, A - 1] + Cat[i, A] + tiny))
#       Nat[i, A] = (Nat[i + 1, A]/S + Cat[i, A - 1] + Cat[i, A]) *
#         (Cat[i, A]/(Cat[i, A - 1] + Cat[i, A] + tiny))
#     }
#     # modify this parameters if need it ##
# 
#     # minagesel = 8
#     minagecom = 1
# 
#     Ut = rowSums(Cat[, minagesel:(A - 1)])/rowSums(Nat[, minagesel:(A -
#                                                                       1)])  # Exploitation rate for fully recruited fish
#     # Ut[n] = 0.4
#     vat = Cat/Nat/Ut  # relative vulnerablility at age
#     va = colMeans(vat[(n - avg_yrsel):(n - minagecom), ])  # update terminal vul
#     va[minagesel:A] = 1
# 
#   }
# 
#   vat[is.na(vat)] <- 1
#   Ut[n] = Uterm
#   if (usewat == T)
#     vbt = rowSums((Nat * vat) * wa) else vbt = (Nat * vat) %*% wa
#   if (usewat == T)
#     bt = as.vector(rowSums(Nat * wa)) else bt = as.vector(Nat %*% wa)
#   fec = pmat * wa
#   zt = log(yt/vbt)
#   epsilon = zt - mean(zt)
#   if (usewat == T)
#     ssb = as.vector(rowSums(Nat * fec)) else ssb = as.vector(Nat %*% fec)
#   predcpue = exp(mean(zt)) * vbt  ### check again if bt or vbt
#   cpue_q = yt/exp(mean(zt))
#   qhat = exp(mean(epsilon))
# 
#   lnl = sum(dnorm(epsilon, mean = 0, sd = sig, log = T))
# 
#   if (opt) {
#     return(-lnl)
#   } else {
#     # ss = sum(epsilon^2) lnl = 0.5*n*log(ss)
#     return(list(Uterm = Uterm, va = va, rt = Nat[, 1], ssb = ssb, yt = yt,
#                 vbt = vbt, cpue_q = cpue_q, Nat = Nat, vat = vat, Ut = Ut,
#                 bt = bt, predcpue = predcpue, epsilon = epsilon/sig, lnl = lnl,
#                 qhat = qhat, minagesel = minagesel, minagecom = minagecom,
#                 avg_yrsel = avg_yrsel))
#   }
# 
# 
# }
# 
# VPAFMSY <- function(lnFMc, Mc, hc, maxage, vul, Linfc, Kc, t0c, AMc, ac,
#                     bc, opt = T, ny = 50) {
# 
#   FMc <- exp(lnFMc)
# 
#   Mac <- rep(1, maxage)
#   Mac[1:max(1, floor(AMc))] <- 0
#   Lac <- Linfc * (1 - exp(-Kc * ((1:maxage) - t0c)))
#   Wac <- ac * Lac^bc
#   R0c <- 1
#   N <- exp(-Mc * ((1:maxage) - 1)) * R0c
#   SSN <- Mac * N  # Calculate initial spawning stock numbers
#   Biomass <- N * Wac
#   SSB <- SSN * Wac  # Calculate spawning stock biomass
# 
#   B0 <- sum(Biomass)
#   SSB0 <- sum(SSB)
#   SSN0 <- SSN
#   SSBpR <- sum(SSB)/R0c  # Calculate spawning stock biomass per recruit
#   SSNpR <- SSN/R0c
# 
#   N <- N/2
#   SSN <- Mac * N  # Calculate initial spawning stock numbers
#   Biomass <- N * Wac
#   SSB <- SSN * Wac
# 
#   for (y in 1:ny) {
#     # set up some indices for indexed calculation Fishing mortality rate
#     # determined by effort, catchability, vulnerability and spatial
#     # preference according to biomass
#     Zc <- FMc * vul + Mc
#     CN <- N * (1 - exp(-Zc)) * (FMc/Zc)
#     CB <- CN * Wac
#     Biomass <- N * Wac
#     N[2:maxage] <- N[1:(maxage - 1)] * exp(-Zc[1:(maxage - 1)])  # Total mortality
#     N[1] <- (0.8 * R0c * hc * sum(SSB))/(0.2 * SSBpR * R0c * (1 - hc) +
#                                            (hc - 0.2) * sum(SSB))  # Recruitment assuming regional R0 and stock wide steepness
#     # print(N[1])
#     SSN <- N * Mac
#     SSB <- SSN * Wac
# 
#   }  # end of year
# 
#   if (opt) {
#     return(-sum(CB))
#   } else {
#     return(FMc)
#   }
# }
# 

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

