\name{BK}
\alias{BK}
\title{Beddington and Kirkwood life-history MP (simple version)}
\description{
Sets an OFL according to current abundance and an approximation of FMSY based on length at first capture. 
}
\usage{
BK(x, DLM_data, reps = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A position in a data-limited methods data object.
}
  \item{DLM_data}{
A data-limited methods data object.
}
  \item{reps}{
The number of stochastic samples of the TAC recommendation
}
}
\references{
Beddington, J.R., Kirkwood, G.P., 2005. The estimation of potential yield and stock status using life history parameters. Philos. Trans. R. Soc. Lond. B Biol. Sci. 360, 163-170.
}
\author{
T. Carruthers.
}
\note{
This is the simple version of the BK MP. The paper has a more complex approach that might work better. 
}
