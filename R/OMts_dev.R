library(DLMtool); library(dplyr)
OM <- testOM

DLMextra()
library(DLMextra)
OMs <- avail("OM")
dflist <- list(); n <- 0
for (om in OMs) {
  n <- n+1
  OM <- get(om)
  dflist[[n]] <- data.frame(OM=om, maxage=OM@maxage)
}
DF <- do.call("rbind", dflist) %>% as.data.frame()

DF %>% arrange(maxage)


OMts <- function(OM, ts=4, silent=FALSE) {
  if(class(OM) != "OM") stop("Object must be class `OM`", call.=FALSE)
  
  OMout <- OM 
  # Natural mortality
  OMout@M <- OM@M/ts
  OMout@M2 <- OM@M2/ts
  
  # Maximum age
  OMout@maxage <- OM@maxage*ts
  
  # Growth Parameters
  OMout@K <- OM@K/ts
  OMout@t0 <- OM@t0 * ts
  
  # Effort Trends
  OMout@EffYears <- OM@EffYears * ts
  
  # Interval & Years
  OMout@interval <- OM@interval
  OMout@nyears <- OM@nyears * ts
  OMout@proyears <- OM@proyears * ts
  
  # Cpars
  if (length(OM@cpars$M_at_Length) >0) {
    OMout@cpars$M_at_Length <- OM@cpars$M_at_Length/ts
  }

  if (length(OM@cpars$ageM) >0) {
    OMout@cpars$ageM <- OM@cpars$ageM*ts
  }
  if (length(OM@cpars$age95) >0) {
    OMout@cpars$age95 <- OM@cpars$age95*ts
  }
  
  # At-age - can't deal with for now
  nms <- c('M_ageArray', 'Mat_age', 'Len_age', "Wt_age", "V", "retA", "Find", 'Perr',
           "Marray", "Karray")
  for (nm in nms) {
    if (!is.null(OM@cpars[[nm]]))
      stop("Cannot convert cpars$", nm, " to different time-step. Specify externally to converted OM", call.=FALSE)
  }
   
  # Message 
  if (!silent) message("Converting from annual time-step to ", ts, ' time-steps per year')
  OMout
}

# Test Capelin_GSL_DFO - short-lived
# test with curE

OM <- Lesser_Amberjack_GOM_NOAA
OM@cpars$Len_age <- OM@cpars$Wt_age <- OM@cpars$Find <- OM@cpars$Perr <- OM@cpars$V <- 
  OM@cpars$Marray <- NULL
OM2 <- OMts(OM)

nsim <- 400
OM@nsim <- nsim
OM2@nsim <- nsim

MSE1 <- runMSE(OM, MP=c("FMSYref", "curE"), ntrials=300, parallel = TRUE, fracD=0.0001)
MSE2 <- runMSE(OM2, MP=c("FMSYref", "curE"), ntrials=300, parallel = TRUE, fracD=0.0001)


TradePlot(MSE1)
TradePlot(MSE2)


# No error test











