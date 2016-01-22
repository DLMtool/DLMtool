\name{GB_slope}
\alias{GB_slope}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Geromont and Butterworth index slope Harvest Control Rule
}
\description{
An MP similar to SBT1 that modifies a time-series of catch recommendations and aims for a stable catch rates.
}
\usage{
GB_slope(x, DLM_data, reps = 100, yrsmth = 5, lambda = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A position in data-limited methods data object
}
  \item{DLM_data}{
A data-limited methods data object
}
  \item{reps}{
The number of TAC samples
}
  \item{yrsmth}{
Number of years for evaluating slope in relative abundance index
}
  \item{lambda}{
A gain parameter
}
}
\details{
Note that this is my interpretation of their approach and is now stochastic. Currently it is generalized and is not 'tuned' to more detailed assessment data which might explain why in some cases it leads to stock declines. 
}
\references{
Geromont, H.F. and Butterworth, D.S. 2014. Complex assessment or simple management procedures for efficient fisheries management: a comparative study. ICES J. Mar. Sci. doi:10.1093/icesjms/fsu017
}
\author{
T. Carruthers
}
