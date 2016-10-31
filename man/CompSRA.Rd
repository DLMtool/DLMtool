\name{CompSRA}
\alias{CompSRA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Age-composition-based estimate of current stock depletion given constant Z linked to an FMSY estimate to provide OFL}
\description{
Estimates an OFL based on a Stock Reduction analysis fitted to current age-composition data. Knife-edge vulnerability at age at maturity allows for an FMSY estimate. OFL=FMSY*F/C
}
\usage{
CompSRA(x, DLM_data, reps = 100)
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
The number of stochastic samples of the TAC.
}
}
\author{
T. Carruthers
}
\note{
Given a fixed historical F, What level of depletion gives you this length composition? 
}
