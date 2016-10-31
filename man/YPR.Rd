\name{YPR}
\alias{YPR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Yield Per Recruit analysis to get FMSY proxy F01}
\description{
A simple yield per recruit approximation to FMSY (F01) which is the position of the ascending YPR curve for which dYPR/dF = 0.1(dYPR/d0)
}
\usage{
YPR(x, DLM_data, reps = 100)
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
The number of samples of the TAC
}
}
\value{
A numeric vector of TAC samples
}
\references{
Beverton and Holt. 1954. 
}
\author{
Meaghan Bryan and Tom Carruthers
}
\note{
Based on the code of Meaghan Bryan
}
