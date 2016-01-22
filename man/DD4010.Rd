\name{DD4010}
\alias{DD4010}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Delay - Difference Stock Assessment with UMSY and MSY leading coupled with a 40-10 harvest control rule
}
\description{
A simple delay-difference assessment that estimates the OFL using a time-series of catches and a relative abundance index. In this version of the DD MP a 40-10 rule is imposed over the OFL recommendation.  
}
\usage{
DD4010(x, DLM_data, reps = 100)
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
Method based on equations of Carl Walters
}
\author{
T. Carruthers
}
