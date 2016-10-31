\name{MCD}
\alias{MCD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mean Catch Depletion}
\description{
A simple average catch-depletion MP that was included to demonstrate just how informative an estimate of current stock depletion can be. TAC=2*D*AvC
}
\usage{
MCD(x, DLM_data, reps = 100)
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
\value{
A numeric vector of TAC recommendations
}
\author{
T. Carruthers
}

