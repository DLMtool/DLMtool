\name{DCAC_ML}
\alias{DCAC_ML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Depletion-Based Stock Reduction Analysis using mean-length estimator of current depletion
}
\description{
DCAC that uses the mean length estimator to calculate current stock depletion. 
}
\usage{
DCAC_ML(x, DLM_data, reps = 100)
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
The mean length extension was programmed by Gary Nelson as part of his excellent R package 'fishmethods'
}
