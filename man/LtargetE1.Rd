\name{LtargetE1}
\alias{LtargetE1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A management procedure that incrementally adjusts the Effort to reach a target mean length in catches.
}
\description{
A effort based version of the least biologically precautionary of four target length MPs proposed by Geromont and Butterworth 2014.
}
\usage{
LtargetE1(x, DLM_data, reps = 100, yrsmth = 5, xx = 0, xL = 1.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A position in data-limited methods data object
}
  \item{DLM_data}{
A data-limited methods data object
}
  \item{reps}{
The number of samples
}
  \item{yrsmth}{
Years over which to smooth recent estimates of surplus production  
}
\item{xx}{
Parameter controlling the fraction of mean catch to start using in first year
}
\item{xL}{
Parameter controlling the magnitude of the target mean length of catches relative to average length in catches.
}
}
\value{
A numeric vector of input controls
}
\author{
T. Carruthers
}