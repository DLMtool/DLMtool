\name{BK_ML}
\alias{BK_ML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Beddington and Kirkwood life-history analysis with mean-length estimator of current abundance}
\description{
Uses an approximation to FMSY based on length at first capture and an estimate of current abundance based on a mean-length estimator. 
}
\usage{
BK_ML(x, DLM_data, reps = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Position in a data-limited methods data object
}
  \item{DLM_data}{
A data-limited methods data object (class DLM_data)
}
  \item{reps}{
The number of samples of the TAC recommendation
}
}
\references{
Beddington, J.R., Kirkwood, G.P., 2005. The estimation of potential yield and stock status using life history parameters. Philos. Trans. R. Soc. Lond. B Biol. Sci. 360, 163-170.
}
\author{
T. Carruthers
}
\note{
The mean length extension was programmed by Gary Nelson as part of his excellent R package 'fishmethods'
}
