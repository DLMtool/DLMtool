\name{LBSPRSim}
\alias{LBSPRSim}

\title{Per-Recruit Length-Based Spawning Potential Ratio Simulation Model}
\description{
Function that takes stock (biology) and fleet (selectivity and fishing mortality) parameters and returns expected equilibrum size strucure of the catch and spawning potential ratio. 
}
\usage{
LBSPRSim(StockPars, FleetPars, SizeBins = NULL, P = 0.001, Nage = 101)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{StockPars}{
The life-history parameters of the stock.
}
  \item{FleetPars}{
The relative fishing mortality (F/M) and the selectivity-at-length parameters.
}
  \item{SizeBins}{
Optional - specifies the maximum length bin (ToSize) and the length increment of the binned length data (Linc). If empty, Linc = 5 and ToSize = Linf + MaxSD * SDLinf
}
  \item{P}{
Percentage of initial cohort still alive at maximum age.
}
  \item{Nage}{
Maximum age of generic age-strucured model. Not sensitive, unless set too low (e.g., < 50). 
}
}
\references{
Hordyk, A.R., Ono, K., Sainsbury, K.J., Loneragan, N., and Prince, J.D. 2015. Some explorations of the life history ratios to describe length composition, spawning-per-recruit, and the spawning potential ratio. ICES J. Mar. Sci. 72: 204-216.

Hordyk, A.R., Ono, K., Valencia, S.R., Loneragan, N.R., and Prince, J.D. 2015. A novel length-based empirical estimation method of spawning potential ratio (SPR), and tests of its performance, for small-scale, data-poor fisheries. ICES J. Mar. Sci. 72: 217-231.
}
\author{
A. Hordyk
}
