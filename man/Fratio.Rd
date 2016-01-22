\name{Fratio}
\alias{Fratio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An FMSY/M ratio method
}
\description{
Calculates the OFL based on a fixed ratio of FMSY to M multiplied by a current estimate of abundance. 
}
\usage{
Fratio(x, DLM_data, reps = 100)
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
The number of samples of the TAC recommendation
}
}
\details{
A simple method that tends to outperform many other approaches alarmingly often even when current biomass is relatively poorly known. The low stock crash potential is largely due to the quite large difference between Fmax and FMSY for most stocks. 
}
\value{
A numeric vector of TAC recommendations
}
\references{
Gulland, J.A., 1971. The fish resources of the ocean. Fishing News Books, West Byfleet, UK. 

Martell, S., Froese, R., 2012. A simple method for estimating MSY from catch and resilience. Fish Fish. doi: 10.1111/j.1467-2979.2012.00485.x.
}
\author{
T. Carruthers
}
