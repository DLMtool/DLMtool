SampleImpPars <- function(Imp, nsim, nyears, proyears) {
  if (class(Imp) != "Imp" & class(Imp) != "Imp") 
    stop("First argument must be class 'Imp' or 'Imp'")
  
  # === Sample implementation error parameters ====
  TACSD <- runif(nsim, Imp@TACSD[1], Imp@TACSD[2])  # Sampled TAC error (lognormal sd)
  TACFrac <- runif(nsim, Imp@TACFrac[1], Imp@TACFrac[2])  # Sampled TAC fraction (log normal sd)
  
  ESD <- runif(nsim, Imp@ESD[1], Imp@ESD[2])  # Sampled Effort error (lognormal sd)
  EFrac <- runif(nsim, Imp@EFrac[1], Imp@EFrac[2])  # Sampled Effort fraction (log normal sd)
  
  SizeLimSD<-runif(nsim,Imp@SizeLimSD[1],Imp@SizeLimSD[2])
  SizeLimFrac<-runif(nsim,Imp@SizeLimFrac[1],Imp@SizeLimFrac[2])
  
  DiscMort<-runif(nsim,Imp@DiscMort[1],Imp@DiscMort[2])
  
  
}