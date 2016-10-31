\name{DBSRA_ML}
\alias{DBSRA_ML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Depletion-Based Stock Reduction Analysis using mean length estimator of stock depletion}
\description{
DBSRA using the mean length estimator to calculate current stock depletion. 
}
\usage{
DBSRA_ML(x, DLM_data, reps = 100)
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
The number of stochastic samples of the quota recommendation
}
}
\references{
Dick, E.J., MacCall, A.D., 2011. Depletion-Based Stock Reduction Analysis: A catch-based method for determining sustainable yields for data-poor fish stocks. Fish. Res. 110, 331-341.
}
\author{
T. Carruthers
}
\note{
The mean length extension was programmed by Gary Nelson as part of his excellent R package 'fishmethods'
}
