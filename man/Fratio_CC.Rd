\name{Fratio_CC}
\alias{Fratio_CC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A data-limited method that uses FMSY/M ratio and a naive catch-curve estimate of recent Z}
\description{
Calculates the OFL based on a fixed ratio of FMSY to M and a catch curve estimate of current stock size.
}
\usage{
Fratio_CC(x, DLM_data, reps = 100, Fmin = 0.005)
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
The number of TAC samples
}
\item{Fmin}{
Minimum current fishing mortality rate for the catch-curve analysis
}
}
\references{
Gulland, J.A., 1971. The fish resources of the ocean. Fishing News Books, West Byfleet, UK. 

Martell, S., Froese, R., 2012. A simple method for estimating MSY from catch and resilience. Fish Fish. doi: 10.1111/j.1467-2979.2012.00485.x.
}
\author{
T. Carruthers
}
