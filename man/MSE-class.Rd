\name{MSE-class}
\Rdversion{1.1}
\docType{class}
\alias{MSE-class}
%%\alias{initialize,MSE-method}
%%\alias{plot,MSE-method}

\title{Class \code{"MSE"}}
\description{
A Management Strategy Evaluation object that contains information about simulation conditions and performance of data-limited methods
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MSE", Name, nyears, proyears, nMPs, MPs, nsim, OMtable, Obs, B_BMSYa, F_FMSYa, Ba, FMa, Ca, OFLa)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{Name}:}{Name of the MSE object}
    \item{\code{nyears}:}{The number of years for the historical simulation}
    \item{\code{proyears}:}{The number of years for the projections - closed loop simulations}
    \item{\code{nMPs}:}{Number of management procedures simulation tested}
    \item{\code{MPs}:}{The names of the MPs that were tested}
    \item{\code{nsim}:}{Number of simulations}
    \item{\code{OM}:}{A table of nsim rows with a column for each sampled parameter of the operating model\cr
     \itemize{
     \item RefY: reference yield, the highest long-term yield (mean over last five years of projection) obtained from a fixed F strategy. This is a useful reference point for framing performance of MPs because it standardizes for starting point and future productivity. 
     \item M: instantaneous natural mortality rate
     \item Depletion: stock depletion (biomass / unfished biomass) in the final historical year (prior to projection)
     \item A: abundance (biomass) updated in each management update of projection
     \item BMSY_B0: most productive stock size relative to unfished
     \item FMSY_M: fishing mortality rate divided by natural mortality rate
     \item Mgrad: mean average percentage gradient in natural mortality rate (percentage per time step)
     \item Msd: interannual variability in natural mortality rate (lognormal CV)
     \item procsd: process error - CV in log-normal recruitment deviations
     \item Esd: interannual variability in historical effort (fishing mortality rate)
     \item dFfinal: gradient in fishing mortality rate over final five years of the historical simulation
     \item MSY: Maximum Sustainable Yield
     \item qinc: mean percentage increase in fishing efficiency (catchability) in projected years (input controls only)
     \item qcv: interannual variability in future fishing efficiency (catchability) in projected years (input controls only)
     \item CALcv: variability in lengths at age around the growth curve (normal CV)
     \item FMSY: Fishing mortality rate at Maximum Sustainable Yield
     \item Linf: maximum length (von Bertalanffy Linf parameter)
     \item K: maximum growth rate (von Bertalanffy K parameter)
     \item t0: theoretical length at age zero (von Bertalanffy t0 parameter)
     \item hs: steepness of the stock recruitment relationship (the fraction of unfished recruitment at a fifth of unfished stock levels)
     \item Linfgrad: mean gradient in maximum length (per cent per time step)
     \item Kgrad: mean gradient in maximum growth rate (per cent per time step)
     \item Linfsd: interannual variability in maximum length (log normal CV)
     \item recgrad: gradient in recruitment strength (age 1 population numbers) over last 10 years of historical simulations
     \item Ksd: interannual variability in maximum growth rate (log normal CV)
     \item ageM: age at 50 per cent maturity
     \item LFS: length at full selection (the shortest length class where fishery selectivity is 100 per cent)
     \item age05: the age at 5 percent selectivity (ascending limb of selectivity curve)
     \item Vmaxage: the selectivity of the oldest age class (controls dome shape of selectivity curve)
     \item LFC: length at first capture, the smallest length that can be caught by the gear
     \item OFLreal: the true simulated Over Fishing Limit (FMSY x biomass) updated in each management update of the projection
     \item Spat_targ: spatial targetting parameter, fishing mortality rate across areas is proportional to vulnerable biomass raised to the power of this number. 
     \item Frac_area_1: the fraction of unfished biomass inhabiting area 1 (can be seen as fraction of habitat in area 1 or relative size of area 1)
     \item Prob_staying: the probability that individuals in area 1 remain there between time-steps
     \item AC: autocorrelation in recruitment
    }}
    \item{\code{Obs}:}{A table of nsim rows with a column for each sampled parameter of the observation model\cr
     \itemize{
     \item Cbias: bias in observed catches
     \item Csd: observation error in observed catches (lognormal CV)
     \item CAA_nsamp: the number of catch-at-age observations per time step
     \item CAA_ESS: the effective sample size of multinomial catch-at-age observation model (number of independent draws)
     \item CAL_nsamp: the number of catch-at-length observations per time step
     \item CAL_ESS: the effective sample size of multinomial catch-at-length observation model (number of independent draws)
     \item Isd: observation error in relative abundance index (lognormal CV)
     \item Dbias: bias in observed stock depletion (also applies to depletion Dt for DCAC)
     \item Mbias: bias in observed natural mortality rate
     \item FMSY_Mbias: bias in ratio of FMSY to natural mortality rate
     \item BMSY_B0bias: bias in ratio of most productive stock size relative to unfished
     \item AMbias: bias in age at 50 per cent maturity
     \item LFCbias: bias in length at first capture
     \item LFSbias: bias in length at full selection
     \item Abias: bias in observed current absolute stock biomass
     \item Kbias: bias in maximum growth rate (von Bertalanffy K parameter)
     \item t0bias: bias in theoretical length at age zero (von Bertalanffy t0 parameter)
     \item Linfbias: bias in maximum length (von Bertalanffy Linf parameter)
     \item hbias: bias in observed steepness of the stock recruitment relationship
     \item Irefbias: bias in abundance index corresponding to BMSY stock levels
     \item Crefbias: bias in MSY prediction (target or reference catch)
     \item Brefbias: bias in BMSY stock levels (target or reference biomass levels)}}
    \item{\code{B_BMSY}:}{Stored biomass relative to BMSY over the projection (an array with dimensions nsim, nMPs, proyears) }
    \item{\code{F_FMSY}:}{Stored fishing mortality rate relative to FMSY over the projection (an array with dimensions nsim, nMPs, proyears) }
    \item{\code{B}:}{Stored stock biomass over the projection (an array with dimensions nsim, nMPs, proyears)}
    \item{\code{FM}:}{Stored fishing mortality rate over the projection (an array with dimensions nsim, nMPs, proyears)}
    \item{\code{C}:}{Stored catches (taken) over the projection (an array with dimensions nsim, nMPs, proyears) }
    \item{\code{TAC}:}{Stored Total Allowable Catch (prescribed) over the projection (an array with dimensions nsim, nMPs, proyears)(note that this is NA for input controls) }
    \item{\code{SSB_hist}:}{Stored historical spawning stock biomass (historical simulations - an array with dimensions nsim, nages, nyears, nareas)}
    \item{\code{CB_hist}:}{Stored historical catches in weight (historical simulations - an array with dimensions nsim, nages, nyears, nareas)}
    \item{\code{FM_hist}:}{Stored historical fishing mortality rate (historical simulations - an array with dimensions nsim, nages, nyears, nareas)}

  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "MSE")}: ... }
    \item{plot}{\code{signature(x = "MSE")}: ... }
	 }
}
\author{
T. Carruthers
}
\keyword{classes}
