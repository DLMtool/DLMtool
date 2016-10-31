\name{Fratio4010}
\alias{Fratio4010}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{An FMSY/M ratio method paired with the 40-10 rule}
\description{
Calculates the OFL based on a fixed ratio of FMSY to M multiplied by a current estimate of abundance. In this method DBSRA is paired with the 40-10 rule that throttles back the OFL to zero at 10 percent of unfished biomass.
}
\usage{
Fratio4010(x, DLM_data, reps = 100)
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
}

\references{
Gulland, J.A., 1971. The fish resources of the ocean. Fishing News Books, West Byfleet, UK. 

Martell, S., Froese, R., 2012. A simple method for estimating MSY from catch and resilience. Fish Fish. doi: 10.1111/j.1467-2979.2012.00485.x.
}
\author{
T. Carruthers
}
