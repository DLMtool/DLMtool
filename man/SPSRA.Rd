\name{SPSRA}
\alias{SPSRA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Surplus Production Stock Reduction Analysis
}
\description{
A surplus production equivalent of DB-SRA that uses a demographically derived prior for intrinsic rate of increase (McAllister method, below)
}
\usage{
SPSRA(x, DLM_data, reps = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
A position in a data-limited methods data object
}
\item{DLM_data}{
A data-limited methods data object (class DLM)
}
  \item{reps}{
The number of samples of the TAC taken for the calculation of the quota
}

}
\references{
McAllister, M.K., Pikitch, E.K., and Babcock, E.A. 2001. Using demographic methods to construct Bayesian priors for the intrinsic rate of increase in the Schaefer model and implications for stock rebuilding. Can. J. Fish. Aquat. Sci. 58: 1871-1890.
}
\author{
T. Carruthers
}
