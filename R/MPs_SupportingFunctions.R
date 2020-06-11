##### MP plotting #####

AvC_plot <- function(x, Data, Rec, meanC, histCatch, yr.ind, lwd=3, cex.lab=1.25) {
  op <- par(no.readonly = TRUE)
  on.exit(op)
  par(mfrow=c(1,1))
  plot(c(Data@Year[yr.ind], Data@Year[max(yr.ind)]+1), c(histCatch,NA), type="l", 
       xlab="Year", ylab=paste0("Catch (", Data@Units, ")"), lwd=lwd, bty="l", las=1, cex.lab=cex.lab)
  abline(v=Data@LHYear[1], lty=2, col="darkgray") #
  text(Data@LHYear[1], max(histCatch, na.rm=TRUE)*0.9, "Last Historical Year", pos=2, xpd=NA)
  lines(c(min(Data@Year), Data@LHYear[1]), rep(mean(Data@Cat[x,yr.ind]),2), lty=2) #
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
  abline(v=Data@LHYear[1], lty=2, col="darkgray") #
  
  text(Data@LHYear[1], max(Data@Cat[x,1:yr.lst], na.rm=TRUE)*0.9, "Last Historical Year", pos=2, xpd=NA)
  lines(c(min(Data@Year), Data@LHYear[1]), rep(mean(Data@Cat[x,1:yr.lst]),2), lty=2) #
  text(quantile(Data@Year, 0.1), mean(Data@Cat[x,1:yr.lst])*1.1, pos=4, "Average Historical Catch")

  boxplot(TAC, add=TRUE, at=max(Data@Year)+1, col="darkgrey", width=1, outline=TRUE, axes=FALSE)
  text(max(Data@Year)+1, quantile(TAC, 0.95, na.rm=TRUE), "TAC", col="black", pos=3)

  par(new = T)
  plot(c(1, max(Data@Year)+3), c(0,1), type="n", axes=FALSE, xlab="", ylab="")
  quants <- quantile(Bt_K, c(0.025, 0.5, 0.975), na.rm=TRUE)
  points(max(Data@Year)+3, quants[2], pch=16, col="blue", cex=1.5)
  lines(c(max(Data@Year)+3, max(Data@Year)+3), c(quants[1], quants[3]), col="blue")
  axis(side=4, las=1, col="blue", labels=FALSE)
  at = graphics::axTicks(4)
  mtext(side = 4, text = at, at = at, col = "blue", line = 1, las=1)
  mtext(side=4, "Depletion (median + 95 percentiles)", line=3, cex=1.25, col="blue")
}

BK_plot <- function(DF) {
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package \"gridExtra\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  vars <- vals <- NULL # R check hack
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




DD_plot <- function(x, runDD, Data, TAC=NULL, Eff=NULL) {
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
  
  if (!is.null(TAC)) {
    plot(Years, c(C_hist, NA), type='l', bty="n", 
         xlab="Year", ylab=paste("Catch (", Data@Units, ")"), las=0, ylim=c(0, max(c(C_hist, TAC, Cpred_DD), na.rm=TRUE)), lwd=3)
    matplot(Year, Cpred_DD, type="l", lwd=1, add=TRUE)
    
    
    if (all(round(TAC / mean(TAC, na.rm=TRUE),1) ==1 )) {
      points(max(Years), mean(TAC, na.rm=TRUE), pch=16, cex=2, col="blue")
      text(max(Years), mean(TAC, na.rm=TRUE), "TAC", pos=1, col="blue")
    } else {
      boxplot(TAC, add=TRUE, at=max(Years), col="grey", width=1, outline=TRUE, 
              axes=FALSE)
      text(max(Years), quantile(TAC, 0.95, na.rm=TRUE), "TAC", pos=3, col="blue")
    }
  } else {
    Years <- c(Data@LHYear[1],Data@LHYear+1)
    Eff <- c(Data@MPeff[x], Eff)
    plot(Years, Eff, type="b", ylab="Effort", axes=FALSE)
    axis(side=1, at=Years)
    axis(side=2)
    abline(v=Years[1], col="lightgray", lty=2)
    text(Years[1], Eff[1], pos=4, "Previous Year")
    text(Years[2], Eff[2], pos=2, "Next Year")
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
  
  ylim <- c(0, max(c(TAC,(exp(C_dat))), na.rm=TRUE))
  years <- c(years, max(years)+1)
  plot(years, c(exp(C_dat), NA), pch=16, bty="n", ylim=ylim,
       xlab=paste0("Year (last ", yrsmth, " years)"), ylab=paste0("Catch (", Data@Units, ")"))
  lines(years, c(C_hist, NA))
  legend("topright", bty="n", lty=1, legend="Smoothed Catch")
  if (all(round(TAC / mean(TAC, na.rm=TRUE),1) ==1 )) {
    points(max(years), mean(TAC, na.rm=TRUE), pch=16, cex=2, col="blue")
    text(max(years), mean(TAC, na.rm=TRUE), "TAC", pos=1, col="blue")
  } else {
    boxplot(TAC, add=TRUE, at=max(years), col="blue", width=1, outline=TRUE, 
            axes=FALSE)
    text(max(years), quantile(TAC, 0.95, na.rm=TRUE), "TAC", pos=3, col="blue")
  }
  
  boxplot(data.frame(SP_Gradient=G_new, Fmsy=Frat, updatedF=newF), bty="l")
}

EtargetLopt_plot <- function(x, rec, Data, ind,Lopt) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow=c(1,2))
  
  ylim <- range(c(Data@ML[ind], Lopt))
  plot(Data@Year[ind], Data@ML[ind], xlab="Year", ylab="Mean Length", bty="l",
       pch=16, ylim=ylim)
  abline(h=mean(Data@ML[ind], na.rm=TRUE), lty=2)
  text(mean(Data@Year[ind]), mean(Data@ML[ind]), pos=3, "Mean")
  
  abline(h=Lopt, lty=3)
  text(mean(Data@Year[ind]), Lopt, pos=3, "Lopt", xpd=NA)
  
  
  Years <- c(Data@LHYear[1],Data@LHYear+1)
  Eff <- c(Data@MPeff[x], rec@Effort)
  plot(Years, Eff, type="b", ylab="Effort", axes=FALSE)
  axis(side=1, at=Years)
  axis(side=2)
  abline(v=Years[1], col="lightgray", lty=2)
  text(Years[1], Eff[1], pos=4, "Previous Year")
  text(Years[2], Eff[2], pos=2, "Next Year")
}


Fadapt_plot <- DynF_plot

Fdem_plot <- function(runFdem, Data) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow=c(1,2))
  boxplot(cbind(runFdem$Ac, runFdem$TAC), names=c("Abundance", "TAC"), ylab=Data@Units)
  if (all(round(mean(runFdem$FMSY)/runFdem$FMSY,1)==1)) {
    fmsy <- mean(runFdem$FMSY)
    boxplot(fmsy, ylab=expression("F"[MSY]))
  } else {
    boxplot(runFdem$FMSY, ylab=expression("F"[MSY]))  
  }
  
}


Fratio_plot <- function(x, Data, TAC, runFrat) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  if ("Bt_K" %in% names(runFrat)) {
    par(mfrow=c(1,3))
    incBt_K <- TRUE
  } else {
    par(mfrow=c(1,2))  
    incBt_K <- FALSE
  }
  
  boxplot(cbind(runFrat$Abun, TAC), names=c("Abundance", "TAC"), ylab=Data@Units)
  if (all(round(mean(runFrat$Frat)/runFrat$FMSY,1)==1)) {
    fmsy <- mean(runFrat$Frat)
    boxplot(fmsy, ylab=expression("F"[MSY]))
  } else {
    boxplot(runFrat$Frat, ylab=expression("F"[MSY]))  
  }
  if (incBt_K) {
    if (all(round(mean(runFrat$Bt_K)/runFrat$Bt_K,1)==1)) {
      Bt_K <- mean(runFrat$Bt_K)
      boxplot(Bt_K,ylab="Depletion")  
    } else {
      boxplot(runFrat$Bt_K, ylab="Depletion")  
    }
  }
}


GB_CC_plot <- function(x, Catrec, TAC, Data) {
  ylim <- range(c(TAC, Catrec))
  tt <- boxplot(TAC, ylab=paste0("TAC (", Data@Units, ")"), ylim=ylim)
  points(1, Data@Cref[x], pch=16, col="orange", cex=2)
  text(1, Data@Cref[x], pch=16, col="orange", "Cref", pos=2)
  
  points(1, Catrec, pch=16, col="blue", cex=2)
  text(1, Catrec, pch=16, col="blue", "Last Catch", pos=2)
  
}

GB_slope_plot <- function(Data, ind, I_hist, MuC, TAC, Islp) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow=c(1,3))
  yrs <- Data@Year[ind]
  plot(yrs, I_hist, xlab="Year", ylab="Index of Abundance", type="l", lwd=2, bty="l", las=1)
  boxplot(Islp, ylab="log Index slope")
  boxplot(cbind(MuC, TAC), names=c("Last Catch", "TAC"), ylab=Data@Units)
}



GB_target_plot <- function(Itarg, Irec, I0, Data, Catrec, TAC) {
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow=c(1,2))
  
  ylim <- range(c(Itarg, Irec, I0))
  boxplot(Itarg, ylim=ylim, ylab="Index target")
  points(1, Irec, pch=16, col="blue", cex=2)
  text(1, Irec, "Irec", col="blue", cex=1.25, pos=2)
  
  points(1, I0, pch=16, col="orange", cex=2)
  text(1, I0, "I0", col="orange", cex=1.25, pos=2)
  
  
  ylim <- range(c(Catrec, TAC))
  boxplot(cbind(Catrec, TAC), ylim=ylim, ylab=Data@Units, names=c("Last Catch", "TAC"))
  
}



Gcontrol_plot <- function(years, ind1, yrsmth, SP_new, SP_hist, B_dat, B_hist, C_dat, C_hist, TAC, Data) {
  op <- par(no.readonly = TRUE)
  on.exit(op)
  
  par(mfrow=c(2,2))
  
  plot(years, exp(B_dat), pch=16, bty="n", 
       xlab=paste0("Year (last ", yrsmth, " years)"), ylab=paste0("Abundance (", Data@Units, ")"))
  lines(years, B_hist)
  legend("topright", bty="n", lty=1, legend="Smoothed Abundance")
  
  plot(B_hist[ind1], SP_hist, type="p", pch=16, bty="n", 
       xlab=paste0("Biomass (last ", yrsmth-1, " years)"), ylab='Surplus Production')
  
  ylim <- c(0, max(c(TAC,(exp(C_dat))), na.rm=TRUE))
  years <- c(years, max(years)+1)
  plot(years, c(exp(C_dat), NA), pch=16, bty="n", ylim=ylim,
       xlab=paste0("Year (last ", yrsmth, " years)"), ylab=paste0("Catch (", Data@Units, ")"))
  lines(years, c(C_hist, NA))
  legend("topright", bty="n", lty=1, legend="Smoothed Catch")
  if (all(round(TAC / mean(TAC, na.rm=TRUE),1) ==1 )) {
    points(max(years), mean(TAC, na.rm=TRUE), pch=16, cex=2, col="blue")
    text(max(years), mean(TAC, na.rm=TRUE), "TAC", pos=1, col="blue")
  } else {
    boxplot(TAC, add=TRUE, at=max(years), col="blue", width=1, outline=TRUE, 
            axes=FALSE)
    text(max(years), quantile(TAC, 0.95, na.rm=TRUE), "TAC", pos=3, col="blue")
  }
  
  boxplot(SP_new, bty="l", ylab="Predicted Surplus Production")
}



ICI_plot <- function(Years, Index, ci.low, ci.high, TAC, Cat, Data) {
  op <- par(no.readonly = TRUE)
  on.exit(op)
  par(mfrow=c(1,2))
  plot(Years, Index, type="l", bty="l", lwd=2, xlab="Years", ylab="Index")
  
  lines(Years, rep(mean(ci.low, na.rm=TRUE), length(Years)), lty=2)
  text(quantile(Years,0.05), mean(ci.low, na.rm=TRUE), pos=1, "CI_low")
  lines(Years, rep(mean(ci.high, na.rm=TRUE), length(Years)), lty=2)
  text(quantile(Years, 0.05), mean(ci.high, na.rm=TRUE), pos=3, "CI_high")
  
  ylim <- range(c(TAC, Cat))
  boxplot(TAC, col="grey", width=1, outline=TRUE, ylab=paste0("TAC (", Data@Units, ")"), ylim=ylim)
  points(Cat, pch=16, cex=1.5, col="blue", xpd=NA)
  text(1, Cat, "Last Catch", pos=2, col="blue")
  
}



Iratio_plot <- function(Data, I.num, ind.num, I.den, ind.den, alpha, TAC, Cat) {
  op <- par(no.readonly = TRUE)
  on.exit(op)
  par(mfrow=c(1,3))
  plot(Data@Year, Data@Ind, xlab="Year", ylab="Index", bty="l", lwd=2, type="l")
  
  lines(Data@Year[ind.num], rep(mean(I.num), length(ind.num)), lty=2, col='blue')
  lines(Data@Year[ind.den], rep(mean(I.den), length(ind.den)), lty=3, col='blue')
  
  boxplot(alpha, ylab="alpha", col="grey")
  
  ylim <- range(c(TAC, Cat))
  boxplot(TAC, col="grey", width=1, outline=TRUE, ylab=paste0("TAC (", Data@Units, ")"), ylim=ylim)
  points(Cat, pch=16, cex=1.5, col="blue", xpd=NA)
  text(1, Cat, "Last Catch", pos=2, col="blue")
  
}



Islope_plot <- function(runIslope, Data) {
  op <- par(no.readonly = TRUE)
  on.exit(op)
  par(mfrow=c(1,2))
  plot(runIslope$Years, log(runIslope$I_hist), xlab="Year", ylab="log Index", bty="l", type="l", lwd=2)
  
  ylim <- range(c(runIslope$C_dat, runIslope$TACstar, runIslope$TAC))
  Years1 <- c(runIslope$Years, max(runIslope$Years)+1)
  plot(Years1, c(runIslope$C_dat, NA), xlab="Year", ylab=paste0("Catch (", Data@Units, ")"), 
       bty="l", type="l", lwd=2, ylim=ylim)
  points(max(runIslope$Years), mean(runIslope$TACstar, na.rm=TRUE), col="blue", cex=2, pch=16)
  text(max(runIslope$Years), mean(runIslope$TACstar, na.rm=TRUE), col="blue",'TAC*', pos=2, cex=1.5)
  
  boxplot(at=max(Years1), runIslope$TAC, add=TRUE, axes=FALSE, col="gray")
  text(max(Years1), quantile(runIslope$TAC, 0.95), pos=3, "TAC", col="black", cex=1.5, xpd=NA)
  
}

ITe_plot <- function(x, Data, rec, yrsmth) {
  ind <- max(1, (length(Data@Year) - yrsmth + 1)):length(Data@Year)
  deltaI <- mean(Data@Ind[x, ind])/Data@Iref[x]
  
  op <- par(no.readonly = TRUE)
  on.exit(op)
  par(mfrow=c(1,2), oma=c(1,1,1,1), mar=c(5,4,1,4))
  
  ylim <- range(c(Data@Ind[x, ind],Data@Iref[x]))
  plot(Data@Year[ind], Data@Ind[x, ind], ylim=ylim, type="b", bty="l",
       xlab="Year", ylab="Index", lwd=2)
  abline(h= mean(Data@Ind[x, ind]), lty=2)
  text(median(Data@Year[ind]), mean(Data@Ind[x, ind]), "Mean Index", pos=3)
  abline(h=Data@Iref[x], lty=3)
  text(median(Data@Year[ind]), Data@Iref[x], "Target Index", pos=3)     
  
  plot(c(max(Data@Year), max(Data@Year)+1), c(Data@MPeff[x], rec@Effort), 
       type="b", xlab="Year", ylab="Effort", bty="l", lwd=2)
  
}

Lratio_BHI_plot <- function(mlbin, CAL, LSQ, Lref, Data, x, TAC, Cc, yrsmth) {
  op <- par(no.readonly = TRUE)
  on.exit(op)
  par(mfrow=c(1,2))
  plot(mlbin, CAL, type="l", xlab="Length", ylab=paste0("Count (last ", yrsmth, " years)"),
       bty="l", lwd=2)
  abline(v=mean(LSQ), lty=3)
  text(mean(LSQ), quantile(CAL, 0.75), pos=4, "Mean length")
  
  abline(v=mean(Lref), lty=3, col="blue")
  text(mean(Lref), quantile(CAL, 0.95), pos=4, "Reference length", col="blue")
  
  
  ylim <- range(c(Data@Cat[x,], TAC, Cc ))
  
  plot(c(Data@Year, max(Data@Year)+1), c(Data@Cat[x,],NA), type="l", xlab="Year", 
       ylab=paste0("Catch (", Data@Units, ")"),
       lwd=2,  bty="l", ylim=ylim)
  
  boxplot(TAC, col="blue", axes=FALSE, at=max(Data@Year)+1, add=TRUE)
  text(max(Data@Year)+1, quantile(TAC, 0.95, na.rm=TRUE), "TAC", col="blue")
  points(max(Data@Year), mean(Cc, na.rm=TRUE), cex=2, pch=16, col="green")
  text(max(Data@Year), mean(Cc, na.rm=TRUE), "last Catch", col="green", pos=1, xpd=NA)
  
}



MCD_plot <- function(Data, AvC, Bt_K, TAC) {
  op <- par(no.readonly = TRUE)
  on.exit(op)
  par(mfrow=c(1,3))
  boxplot(Bt_K, ylab=paste0("Depletion"), col="gray")
  ylim <- range(c(AvC, TAC), na.rm=TRUE)
  boxplot(AvC, ylab=paste0("Mean catch (", Data@Units, ")"), ylim=ylim, col="gray")
  boxplot(TAC, ylab=paste0("TAC (", Data@Units, ")"), ylim=ylim, col="gray")
  
}




Rcontrol_plot <- function(rsamp, ind, G_new, B_hist, SP_hist, SP_mu, B_dat, Data, yind, C_hist, TAC, TACa) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow=c(2,3))
  
  
  
  plot(Data@Year[ind], B_hist, type="l", bty="l", xlab="Year",
       ylab="Smoothed Biomass", lwd=2)
  points(Data@Year[ind], exp(B_dat))
  
  ylim <- range(c(SP_hist, SP_mu))
  plot(c(Data@Year[yind], max(Data@Year[yind])+1), c(SP_hist, NA), type="l", 
       bty="l", xlab="Year", ylab="Surplus Production", lwd=2, ylim=ylim)
  points(max(Data@Year[yind])+1, mean(SP_mu), pch=16, cex=1.5, col='blue')
  text(max(Data@Year[yind])+1, mean(SP_mu), pos=2, "predicted SP", col="blue", xpd=NA)
  
  
  ylim <- range(c(C_hist, TAC, TACa))
  plot(c(Data@Year[ind], max(Data@Year[ind])+1), c(C_hist, NA), type="l", 
       bty="l", xlab="Year", ylab=paste0("Catch (", Data@Units, ")"), lwd=2, ylim=ylim)
  points(max(Data@Year[ind])+1, median(TACa, na.rm=TRUE), cex=1.5, pch=16, col="red")
  text(max(Data@Year[ind])+1, median(TACa, na.rm=TRUE), "TAC_init", col="red", pos=2)
  boxplot(at=max(Data@Year[ind])+1, TAC, col='blue', axes=FALSE, add=TRUE)
  text(max(Data@Year[ind])+1, quantile(TAC, 0.95, na.rm=TRUE), pos=2, "TAC_adj", col="blue", xpd=NA)
  
  boxplot(rsamp, ylab="intrinsic rate of increase", col="gray")
  boxplot(G_new, ylab="G", col="gray")
  
}



SPSRA_plot <- function(runSPSRA, Data, x) {
  op <- par(no.readonly = TRUE)
  on.exit(op)
  par(mfrow=c(3,2), oma=c(1,1,1,1), mar=c(5,4,1,4))
  
  if(all(round(runSPSRA$Ksamp/mean(runSPSRA$Ksamp, na.rm=TRUE),2) == 1)) {
    boxplot(mean(runSPSRA$Ksamp, na.rm=TRUE), ylab="Intrinsic rate of increase")
  } else {
    boxplot(runSPSRA$Ksamp, ylab=paste0("Unfished biomass (", Data@Units, ")")) 
  }
  if(all(round(runSPSRA$dep/mean(runSPSRA$dep, na.rm=TRUE),2) == 1)) {
    boxplot(mean(runSPSRA$dep, na.rm=TRUE), ylab="Depletion")
  } else {
    boxplot(runSPSRA$dep, ylab="Depletion")
  }
  if(all(round(runSPSRA$rsamp/mean(runSPSRA$rsamp, na.rm=TRUE),2) == 1)) {
    boxplot(mean(runSPSRA$rsamp, na.rm=TRUE), ylab="Intrinsic rate of increase")
  } else {
    boxplot(runSPSRA$rsamp, ylab="Intrinsic rate of increase")  
  }
  if(all(round(runSPSRA$MSY/mean(runSPSRA$MSY, na.rm=TRUE),2) == 1)) {
    boxplot(mean(runSPSRA$MSY, na.rm=TRUE), ylab="MSY")
  } else {
    boxplot(runSPSRA$MSY, ylab="MSY")
  }
 
  
  TAC <- runSPSRA$TAC
  
  if(all(round(TAC/mean(TAC, na.rm=TRUE),2) == 1)) {
    boxplot(mean(TAC, na.rm=TRUE), ylab="TAC")
  } else {
    boxplot(TAC, ylab="TAC")
  }
  
  ylim <- range(c(Data@Cat[x,], TAC), na.rm=TRUE)
  plot(c(Data@Year, max(Data@Year)+1), c(Data@Cat[x,], NA), xlab="Year", ylab=paste0('Catch (', Data@Units, ')'),
       bty="l", type="l", lwd=2, ylim=ylim)
  boxplot(TAC, axes=FALSE, add=TRUE, at=max(Data@Year)+1, col="blue", width=2)
  
  
}


YPR_plot <- function(runYPR, Data,reps) {
  frates <- runYPR$frates
  ypr <- runYPR$ypr
  dif <- runYPR$dif
  F0.1 <- runYPR$F0.1
  Ac <- runYPR$Ac
  TAC <- runYPR$TAC
  
  op <- par(mfrow=c(1,3))
  on.exit(par(op, no.readonly = TRUE))
  matplot(frates, t(ypr), type="l", col=1:reps, xlab="F", ylab="Yield-per-recruit")
  savey <- NULL
  for (r in 1:reps) {
    savey[r] <- ypr[r,which.min(dif[r,])]
    points(F0.1[r], savey[r], col=r, pch=16)
  }
  
  text(F0.1[which.max(savey)], ypr[ which.max(savey),which.min(dif[ which.max(savey),])], "F0.1", pos=2, col=which.max(savey))
  text(F0.1[which.min(savey)], ypr[ which.min(savey),which.min(dif[ which.min(savey),])], "F0.1", pos=2, col=which.min(savey))
  
  boxplot(Ac, ylab=paste0("Abundance (", Data@Units, ")"))
  boxplot(TAC, ylab=paste0("TAC (", Data@Units, ")"))
}


size_lim_plot <- function(x, Data, Rec) {
  Val <- Var <- NULL # cran check hacks
  Linf <- Data@vbLinf[x]
  Lens <- 1:Linf
  
  LR5 <- Rec@LR5
  LFR <- Rec@LFR 
  Rmaxlen <- Rec@Rmaxlen
  if (length(Rmaxlen)<1) Rmaxlen <-1
  HS <- Rec@HS 
  L5 <- Rec@L5
  LFS <- Rec@LFS 
  Vmaxlen <- Rec@Vmaxlen
  if (length(Vmaxlen)<1) Vmaxlen <-1
  
  Ret <- Sel <- rep(NA, length(Lens))
  if (length(LR5)>0 && length(LFR)>0) {
    srs <- (Linf - LFR) / ((-log(Rmaxlen,2))^0.5)
    srs[!is.finite(srs)] <- Inf
    sls <- (LFR - LR5) /((-log(0.05,2))^0.5)
    Ret <- getsel(1,lens=Lens, lfs=LFR, sls, srs)
  }
  if (length(HS)>0) Ret[Lens>=HS] <- 0 
  
  if (length(L5)>0 && length(LFS)>0) {
    srs <- (Linf - LFS) / ((-log(Vmaxlen,2))^0.5)
    srs[!is.finite(srs)] <- Inf
    sls <- (LFS - L5) /((-log(0.05,2))^0.5)
    Sel <- getsel(1,lens=Lens, lfs=LFS, sls, srs)
  }
  
  df <- data.frame(Lens=rep(Lens,2),Val=c(Ret, Sel), 
                   Var=rep(c("Retention", "Selectivity"), each=length(Lens)))
  
  p1 <- ggplot2::ggplot(data=subset(df, !is.na(Val)), ggplot2::aes(x=Lens, y=Val, color=Var)) + 
    ggplot2::geom_line(size=1.2) +
    ggplot2::theme_classic() +
    ggplot2::labs(x="Length", y="Proportion", color="")
  print(p1)
  
}

curE_plot <- function(x, rec, Data) {
  
  Years <- c(Data@LHYear[1],Data@LHYear+1)
  Eff <- c(Data@MPeff[x], rec@Effort)
  plot(Years, Eff, type="b", ylab="Effort", axes=FALSE)
  axis(side=1, at=Years)
  axis(side=2)
  abline(v=Years[1], col="lightgray", lty=2)
  text(Years[1], Eff[1], pos=4, "Previous Year")
  text(Years[2], Eff[2], pos=2, "Next Year")
  
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
#' Age-based Catch Curve
#'
#' @param x Iteration number
#' @param Data An object of class `Data`
#' @param reps Number of reps 
#' @param plot Logical. Show the plot?
#'
#' @return A vector of length `reps` of samples of the negative slope of the catch-curve (Z)
#' @export
#'
#' @examples
#' CC(1, DLMtool::SimulatedData, plot=TRUE)
CC <- function(x, Data, reps = 100, plot=FALSE) {
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
    
    if (plot) {
      op <- par(mfrow=c(1,1), no.readonly = TRUE)
      on.exit(par(op))
      plot(xc, y, xlab="Age", ylab="log N", bty="n", pch=16, cex=1.2)
      lines(xc[1:length(predict(mod))], predict(mod))
      text(median(xc[1:length(predict(mod))]), median(predict(mod)), paste0("Z =", round(-coefs[1],2)), pos=4)
    }
    
    return(-rnorm(reps, coefs[1], coefs[2]))
  }
}




## Delay-Difference supporting functions ####


## Mean Length supporting functions ####

bheq <- function(K, Linf, Lc, Lbar) {
  K * (Linf - Lbar)/(Lbar - Lc)
}

#' @importFrom Rcpp evalCpp
bhnoneq <- function(year, mlen, ss, K, Linf, Lc, nbreaks, styrs, stZ) {
  mlen[mlen <= 0 | is.na(mlen)] <- -99
  ss[ss <= 0 | is.na(ss)] <- 0 
  ss[mlen == -99] <- 0
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
  Lc <- Data@Lc[x] # Data@LFS[x]

  for (i in 1:ML_reps) {
    # mlen <- rep(NA, length(year))
    ss <- ceiling(apply(Data@CAL[x, , ], 1, sum)/2)
    if (MLtype == "dep") {
      mlen <- Data@Lbar[x,]
      # for (y in 1:length(year)) {
      #   if (sum(Data@CAL[x, y, ] > 0) > 0.25 * length(Data@CAL[x, y, ])) {
      #     temp2 <- sample(mlbin, ceiling(sum(Data@CAL[x, y, ])/2), replace = T, prob = Data@CAL[x, y, ])
      #     mlen[y] <- mean(temp2[temp2 >= Lc], na.rm = TRUE)
      #   }
      # }
      # 
      fitmod <- bhnoneq(year = year, mlen = mlen, ss = ss, K = Kc[i], Linf = Linfc[i],
                        Lc = Lc, nbreaks = nbreaks, styrs = ceiling(length(year) * ((1:nbreaks)/(nbreaks + 1))),
                        stZ = rep(Data@Mort[x], nbreaks + 1))
      if (all(fitmod == FALSE)) {
        Z[i, ] <- NA
      } else Z[i, ] <- fitmod
    } else {

      # ind<-(which.min(((Data@CAL_bins-Data@LFS[x])^2)^0.5)-1):(length(Data@CAL_bins)-1)
      # for (y in 1:length(year)) {
      #   if (sum(Data@CAL[x, y, ] > 0) > 0.25 * length(Data@CAL[x, y, ])) {
      #     temp2 <- sample(mlbin, ceiling(sum(Data@CAL[x, y, ])/2), replace = T, prob = Data@CAL[x, y, ])
      #     mlen[y] <- mean(temp2[temp2 >= Lc], na.rm = TRUE)
      #   }
      # }
      # mlen <- mean(mlen[(length(mlen) - 2):length(mlen)], na.rm = TRUE)
      mlen <- Data@Lbar[x,]
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
#' @export 
#' @keywords internal 
#' @describeIn sample_steepness2 sample steepness values
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
#' @describeIn sample_steepness2 derive beta parameter
#' @export
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
