\name{BK_CC}
\alias{BK_CC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Beddington and Kirkwood life-history method combined with catch curve analysis
}
\description{
Calculates an OFL using a catch curve estimate of current F and an approximation of FMSY based on length at first capture. 
}
\usage{
BK_CC(x, DLM_data, reps = 100, Fmin=0.005)
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
  \item{Fmin}{
The minimum fishing mortality rate that is derived from the catch-curve (interval censor)
}
}

\references{
Beddington, J.R., Kirkwood, G.P., 2005. The estimation of potential yield and stock status using life history parameters. Philos. Trans. R. Soc. Lond. B Biol. Sci. 360, 163-170.
}
\author{
T. Carruthers
}
