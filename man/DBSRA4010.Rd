\name{DBSRA4010}
\alias{DBSRA4010}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Depletion-Based Stock Reduction Analysis paired with 40-10 harvest control rule
}
\description{
User prescribed BMSY/B0, M, FMSY/M are used to find B0 and therefore the OFL by back-constructing the stock to match a user specified level of stock depletion (OFL = M * FMSY/M * depletion* B0). In this method DBSRA is paried with the 40-10 rule that throttles back the OFL to zero at 10 percent of unfished biomass. 
}
\usage{
DBSRA4010(x, DLM_data, reps = 100)
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
Dick, E.J., MacCall, A.D., 2011. Depletion-Based Stock Reduction Analysis: A catch-based method for determining sustainable yields for data-poor fish stocks. Fish. Res. 110, 331-341.
}
\author{
T. Carruthers
}
