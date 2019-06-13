library(DLMtool); library(dplyr)
# OM <- testOM
# 
# DLMextra()
# library(DLMextra)
# OMs <- avail("OM")
# 

# saveRDS(Jonah_Crab_LFA34_DFO, "saveOM.rdata")

def.args <- DLMtool:::dev.mode(); for (nm in names(def.args)) assign(nm, def.args[[nm]])



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
  OMout@nyears <- OM@nyears 
  OMout@proyears <- OM@proyears 

  # Cpars
  if (length(OM@cpars) > 0) {
    stop("Cannot convert OM with cpars")
  }
  
  # if (length(OM@cpars$M_at_Length) >0) {
  #   OMout@cpars$M_at_Length <- OM@cpars$M_at_Length/ts
  # }
  # 
  # if (length(OM@cpars$ageM) >0) {
  #   OMout@cpars$ageM <- OM@cpars$ageM*ts
  # }
  # if (length(OM@cpars$age95) >0) {
  #   OMout@cpars$age95 <- OM@cpars$age95*ts
  # }
  # 
  # # At-age - can't deal with for now
  # nms <- c('M_ageArray', 'Mat_age', 'Len_age', "Wt_age", "V", "retA", "Find", 'Perr',
  #          "Marray", "Karray")
  # for (nm in nms) {
  #   if (!is.null(OM@cpars[[nm]]))
  #     stop("Cannot convert cpars$", nm, " to different time-step. Specify externally to converted OM", call.=FALSE)
  # }

  # Message
  if (!silent) message("Converting from annual time-step to ", ts, ' time-steps per year')
  OMout@cpars$nts <- ts # add a new slot here
  OMout
}

#

OM_annual <- readRDS("saveOM.rdata")
OM_annual@interval <- 1
OM_annual@K <- c(0.11, 0.15)
OM_annual@cpars <- list()

OM_annual@nsim <- 3
# OM_annual <- tinyErr(OM_annual)
OM_annual@AC <- c(0,0)

OM <- OMts(OM_annual, 4)

MSE1 <- runMSE(OM, Hist=TRUE)

MSE2 <- runMSE(OM_annual, Hist=TRUE)


nsim <- 400
OM@nsim <- nsim
OM2@nsim <- nsim

MSE1 <- runMSE(OM, MP=c("FMSYref", "curE"), ntrials=300, parallel = TRUE, fracD=0.0001)
MSE2 <- runMSE(OM_annual, MP=c("FMSYref", "curE"), ntrials=300, parallel = TRUE, fracD=0.0001)


TradePlot(MSE1)
TradePlot(MSE2)



OM@cpars$Len_age <- OM@cpars$Wt_age <- OM@cpars$Find <- OM@cpars$Perr <- OM@cpars$V <-
  OM@cpars$Marray <- NULL

# No error test

