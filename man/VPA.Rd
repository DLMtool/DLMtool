\name{VPA}
\alias{VPA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust Virtual Population Analysis 
}
\description{
A simple 2 parameter (q, terminal F) VPA. Note that this is an early version that needs more testing.
}
\usage{
VPA(x, DLM_data, reps = 100)
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
\value{
A numeric vector of TAC recommendations
}
\references{
Method based on VPA of Carl Walters and Roberto Licandeo.
}
\author{
C. Walters (Model), R. Licandeo (R code), T. Carruthers (DLMtool implementation)
}
\note{
This is an early version of the VPA and seems not to be working correctly at present. More testing and development is needed. 
}
