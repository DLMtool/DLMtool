##### MP plotting #####

AvC_plot <- function(x, Data, Rec, meanC, histCatch, yr.ind, lwd=3, cex.lab=1.25) {
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
  text(max(Data@Year)+1, quantile(Rec@TAC, 0.05, na.rm=TRUE), "TAC", col="black", pos=2)
}



DCAC_plot <- function(x, Data, dcac, TAC, Bt_K, yrs, lwd=3, cex.lab=1.25) {
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

  boxplot(TAC, add=TRUE, at=max(Data@Year)+1, col="darkgrey", width=1, outline=TRUE, axes=FALSE)
  text(max(Data@Year)+1, quantile(TAC, 0.95, na.rm=TRUE), "TAC", col="black", pos=3)

  par(new = T)
  plot(c(1, max(Data@Year)+3), c(0,1), type="n", axes=FALSE, xlab="", ylab="")
  quants <- quantile(Bt_K, c(0.025, 0.5, 0.975), na.rm=TRUE)
  points(max(Data@Year)+3, quants[2], pch=16, col="blue", cex=1.5)
  lines(c(max(Data@Year)+3, max(Data@Year)+3), c(quants[1], quants[3]), col="blue")
  axis(side=4, las=1, col="blue", labels=FALSE)
  at = axTicks(4)
  mtext(side = 4, text = at, at = at, col = "blue", line = 1, las=1)
  mtext(side=4, "Depletion (median + 95 percentiles)", line=3, cex=1.25, col="blue")
}

BK_plot <- function(DF) {

  DF2 <- DF %>% dplyr::filter(vars %in% c("Lc/Linf", "K", "Fmax"))
  p1 <- ggplot(DF2, aes(x=vars, y=vals)) + geom_boxplot() + 
    theme_classic() + expand_limits(y=0) + labs(x="", y='Values')
  DF3 <- DF %>% dplyr::filter(!vars %in% c("Lc/Linf", "K", "Fmax"))
  p2 <- ggplot(DF3, aes(x=vars, y=vals)) + geom_boxplot() + 
    theme_classic() + expand_limits(y=0) + labs(x="", y='Values')
  
  gridExtra::grid.arrange(p1, p2, nrow=2)
}


CompSRA_plot <- function(runCompSRA, TAC) {
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


DBSRA_plot <- function(runDBSRA, Data, TAC) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow=c(2,2))
  Btrend <- t(runDBSRA$Btrend)
  B0s <- Btrend[1,]
  relB <- Btrend/matrix(B0s, nrow=nrow(Btrend), ncol=ncol(Btrend), byrow=TRUE)
  Years <- Data@Year
  matplot(Years, Btrend, type='l', ylim=c(0, max(Btrend)), bty="n", xlab="Year", ylab=paste("Biomass (", Data@Units, ")"), las=1)
  matplot(Years, relB, type='l', ylim=c(0, max(relB)), bty="n", xlab="Year", ylab='B/B0', las=1)
  if (!is.null(runDBSRA$hcr)) {
    abline(h=runDBSRA$hcr[1], lty=2, col="gray")
    abline(h=runDBSRA$hcr[2], lty=3, col="gray")
  }
  plot(c(Years, max(Years+1)), c(runDBSRA$C_hist, NA), type='l', lwd=2,  bty="n", 
       xlab="Year", ylab=paste("Catch (", Data@Units, ")"), las=1, ylim=c(0, max(c(runDBSRA$C_hist, TAC), na.rm=TRUE)))
  if (all(round(TAC / mean(TAC, na.rm=TRUE),1) ==1 )) {
    points(max(Years)+1, mean(TAC, na.rm=TRUE), pch=16, cex=2, col="blue")
    text(max(Years)+1, mean(TAC, na.rm=TRUE), "TAC", pos=1, col="blue")
  } else {
    boxplot(TAC, add=TRUE, at=max(Years)+1, col="grey", width=1, outline=TRUE, 
            axes=FALSE)
    text(max(Years)+1, quantile(TAC, 0.95, na.rm=TRUE), "TAC", pos=3, col="blue")
  }
  
  df <- data.frame(B_B0=runDBSRA$Bt_Kstore, BMSY_B0=runDBSRA$BMSY_K_Mstore,
                   FMSY_M=runDBSRA$FMSY_Mstore)
  boxplot(df, las=1)
}




DD_plot <- function(x, runDD, Data, TAC) {
  C_hist <- runDD$C_hist 
  I_hist <- runDD$I_hist
  E_hist <- runDD$E_hist
  B_DD <- runDD$B_DD
  dep <- runDD$dep
  Cpred_DD <- runDD$Cpredict
  Year <- runDD$Year
  
  B0est <- matrix(B_DD[1,], nrow=nrow(B_DD), ncol=ncol(B_DD), byrow=TRUE)
  relB <- B_DD/B0est
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow=c(2,2))
  
  Years <- c(Year, max(Year)+1)
  
  matplot(Years, B_DD, type='l', ylim=c(0, max(B_DD)), bty="n", xlab="Year", ylab=paste0("Estimated Biomass (", Data@Units, ")"), las=0)
  matplot(Years, relB, type='l', ylim=c(0, max(relB)), bty="n", xlab="Year", ylab='B/B0', las=1)
  if (!is.null(runDD$hcr)) {
    abline(h=runDD$hcr[1], lty=2, col="gray")
    abline(h=runDD$hcr[2], lty=3, col="gray")
  }
  
  plot(Year, E_hist, bty="n", type="l", lwd=2, ylab="Standardized Index", xlab='Year')
  lines(Year, I_hist, lwd=2, lty=2)
  legend("topleft", bty="n", lwd=2, lty=1:2, legend=c("Effort", "Index"))
  
  
  plot(Years, c(C_hist, NA), type='l', bty="n", 
       xlab="Year", ylab=paste("Catch (", Data@Units, ")"), las=0, ylim=c(0, max(c(C_hist, TAC, Cpred_DD))), lwd=3)
  matplot(Year, Cpred_DD, type="l", lwd=1, add=TRUE)
  
  if (all(round(TAC / mean(TAC),1) ==1 )) {
    points(max(Years), mean(TAC), pch=16, cex=2, col="blue")
    text(max(Years), mean(TAC), "TAC", pos=1, col="blue")
  } else {
    boxplot(TAC, add=TRUE, at=max(Years), col="grey", width=1, outline=TRUE, 
            axes=FALSE)
    text(max(Years), quantile(TAC, 0.95, na.rm=TRUE), "TAC", pos=3, col="blue")
  }
}


DynF_plot <- function(C_dat, C_hist, TAC, yrsmth, B_dat, B_hist, Data, SP_hist, 
                      ind, ind1, G_new, Frat,newF, years) {
  op <- par(no.readonly = TRUE)
  on.exit(op)
  
  par(mfrow=c(2,2))
  
  
  plot(years, exp(B_dat), pch=16, bty="n", 
       xlab=paste0("Year (last ", yrsmth, " years)"), ylab=paste0("Abundance (", Data@Units, ")"))
  lines(years, B_hist)
  legend("topright", bty="n", lty=1, legend="Smoothed Abundance")
  
  plot(B_hist[ind1], SP_hist, type="p", pch=16, bty="n", 
       xlab=paste0("Biomass (last ", yrsmth-1, " years)"), ylab='Surplus Production')
  
  ylim <- c(0, max(c(TAC,(exp(C_dat)))))
  years <- c(years, max(years)+1)
  plot(years, c(exp(C_dat), NA), pch=16, bty="n", ylim=ylim,
       xlab=paste0("Year (last ", yrsmth, " years)"), ylab=paste0("Catch (", Data@Units, ")"))
  lines(years, c(C_hist, NA))
  legend("topright", bty="n", lty=1, legend="Smoothed Catch")
  if (all(round(TAC / mean(TAC),1) ==1 )) {
    points(max(years), mean(TAC), pch=16, cex=2, col="blue")
    text(max(years), mean(TAC), "TAC", pos=1, col="blue")
  } else {
    boxplot(TAC, add=TRUE, at=max(years), col="blue", width=1, outline=TRUE, 
            axes=FALSE)
    text(max(years), quantile(TAC, 0.95, na.rm=TRUE), "TAC", pos=3, col="blue")
  }
  
  boxplot(data.frame(SP_Gradient=G_new, Fmsy=Frat, updatedF=newF), bty="l")
}

Fadapt_plot <- DynF_plot





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
