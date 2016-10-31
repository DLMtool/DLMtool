\name{MCD4010}
\alias{MCD4010}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mean Catch Depletion}
\description{
A simple average catch-depletion MP linked to a 40-10 harvest controle rule that was included to demonstrate just how informative an estimate of current stock depletion can be. TAC=d(1-d)AvC
}
\usage{
MCD4010(x, DLM_data, reps = 100)
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

