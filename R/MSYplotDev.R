# 
# # OM <- testOM
# # silent=TRUE
# # sim=1; yr=1
# OM <- tinyErr(OM)
# Hist <- runMSE(OM, Hist=TRUE)
# 
# Hist@Ref$MSY[sim]
# Hist@Ref$SSB0[sim]
# 
# 
# 
# 
# MSYFuns <- function(OM, sim=1, yr=1, silent=TRUE) {
#   if (class(OM) != "OM") stop('Object must be class `OM`', call.=FALSE)
#   
#   nsim <- OM@nsim
#   nyears <- OM@nyears 
#   if (yr > nyears) {
#     message('yr is greater than `OM@nyears`. Defaulting to `yr=OM@nyears`')
#     yr <- OM@nyears
#   }
#  
#   proyears <- 2
#   SampCpars <- list()
#   SampCpars <- if(length(OM@cpars)>0) SampleCpars(OM@cpars, nsim, msg=!silent)
#   
#   # Stock Parameters & assign to function environment
#   StockPars <- SampleStockPars(OM, nsim, nyears, proyears, SampCpars, msg=!silent)
#   for (X in 1:length(StockPars)) assign(names(StockPars)[X], StockPars[[X]])
#   
#   # Fleet Parameters & assign to function environment
#   FleetPars <- SampleFleetPars(SubOM(OM, "Fleet"), Stock=StockPars, nsim, 
#                                nyears, proyears, cpars=SampCpars)
# 
#   M_at_Age <- StockPars$M_ageArray[sim, , yr]
#   Wt_at_Age <- StockPars$Wt_age[sim,, yr]
#   Mat_at_Age <- StockPars$Mat_age[sim,, yr]
#   V_at_Age <- FleetPars$V[sim,, yr]
#   
#   maxage <- OM@maxage
#   R0x <- StockPars$R0[sim]
#   SRrelx <- StockPars$SRrel[sim]
#   hx <- StockPars$hs[sim]
#   
#   Uvec <- seq(0, 1, by=.025)
#   DF <- sapply(log(Uvec), 
#          MSYCalcs, M_at_Age, Wt_at_Age, Mat_at_Age, V_at_Age, 
#                        maxage, R0x, SRrelx, hx, opt=2, plusgroup=0)
#   DF <- t(DF)
#   DF <- data.frame(DF)
#   
#   DF2 <- data.frame(Ages=1:maxage, N.Mortality=M_at_Age, Weight=Wt_at_Age,
#                     Selectivity=V_at_Age, Maturity=Mat_at_Age)
#   
#   DF2 <- tidyr::gather(DF2, "key", "value", 2:5)
#   
#   p1 <- ggplot2::ggplot(DF2 %>% dplyr::filter(key %in% c("Maturity", "Selectivity")), 
#                   ggplot2::aes(x=Ages, y=value, color=key)) + 
#     ggplot2::geom_line() +
#     ggplot2::geom_point() + 
#     ggplot2::theme_classic() +
#     ggplot2::labs(x="Age (Year)", y="Probability", col="Key")
#   
# 
#   ggplot2::ggplot(DF, ggplot2::aes(x=F, y=Yield)) + ggplot2::geom_line()
#   ggplot2::ggplot(DF, ggplot2::aes(x=F, y=SB)) + ggplot2::geom_line()
#   
#   
#   names(DF)
#   plot(Uvec, DF$Yield)
#   plot(Uvec, DF$SB_SB0)
#   plot(Uvec, DF$B_B0)
#   plot(Uvec, DF$VB_VB0)
#   
#   plot(DF$F, DF$Yield)
#   plot(Uvec, DF$RelRec)
#   plot(DF$SB_SB0, DF$RelRec/R0x)
#   abline(h=hx, lty=2)
#   abline(v=0.2, lty=2)
#   
#   plot(DF$F, Uvec)
#   
# }
