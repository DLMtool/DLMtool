\name{Fdem_CC}
\alias{Fdem_CC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Demographic FMSY method using catch-curve analysis to estimate recent Z}
\description{
FMSY is calculated as r/2 from a demographic r prior method, current abudnance is estimated from naive catch curve analysis.
}
\usage{
Fdem_CC(x, DLM_data, reps = 100, Fmin=0.005)
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
  \item{Fmin}{
The minimum fishing mortality rate derived from the catch-curve analysis
}
}
\references{
McAllister, M.K., Pikitch, E.K., and Babcock, E.A. 2001. Using demographic methods to construct Bayesian priors for the intrinsic rate of increase in the Schaefer model and implications for stock rebuilding. Can. J. Fish. Aquat. Sci. 58: 1871-1890.
}
\author{
T. Carruthers
}
