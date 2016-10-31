\name{DCAC_40}
\alias{DCAC_40}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Depletion Corrected Average Catch assuming 40 per cent stock depletion}
\description{
DCAC assuming that current stock biomass is exactly 40 per cent of unfished levels. 
}
\usage{
DCAC_40(x, DLM_data, reps = 100)
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
MacCall, A.D., 2009. Depletion-corrected average catch: a simple formula for estimating sustainable yields in data-poor situations. ICES J. Mar. Sci. 66, 2267-2271.
}
\author{
T. Carruthers
}
\note{
The 40 percent depletion assumption doesn't really affect DCAC that much as it already makes TAC recommendations that are quite MSY-like. 
}
