\name{DBSRA_40}
\alias{DBSRA_40}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Depletion-Based Stock Reduction Analysis assuming 40 per cent stock depletion
}
\description{
DBSRA assuming that current stock depletion is exactly 40 per cent of unfished stock levels. 
}
\usage{
DBSRA_40(x, DLM_data, reps = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{
A position in a data-limited methods data object
}
  \item{DLM_data}{
A data-limited methods data object
}
  \item{reps}{
The number of stochastic samples of the TAC recommendation
}
}
\references{
Dick, E.J., MacCall, A.D., 2010. Estimates of sustainable yield for 50 data-poor stocks in the Pacific Coast groundfish fishery management plan. Technical memorandum. Southwest fisheries Science Centre, Santa Cruz, CA. National Marine Fisheries Service, National Oceanic and Atmospheric Administration of the U.S. Department of Commerce. NOAA-TM-NMFS-SWFSC-460.
}
\author{
T. Carruthers.
}
\note{
A 40 percent assumption for current depletion is more or less the most optimistic state for a stock (ie very close to BMSY/B0 for many stocks). 
}
