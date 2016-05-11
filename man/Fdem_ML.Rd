\name{Fdem_ML}
\alias{Fdem_ML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Demographic FMSY method that uses mean length data to estimate recent Z
}
\description{
Demographic F (r/2) method using the mean length estimator to calculate current abundance. 
}
\usage{
Fdem_ML(x, DLM_data, reps = 100)
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
}
\references{
McAllister, M.K., Pikitch, E.K., and Babcock, E.A. 2001. Using demographic methods to construct Bayesian priors for the intrinsic rate of increase in the Schaefer model and implications for stock rebuilding. Can. J. Fish. Aquat. Sci. 58: 1871-1890.
}
\author{
T. Carruthers
}
\note{
The mean length extension was programmed by Gary Nelson as part of his excellent R package 'fishmethods'
}
