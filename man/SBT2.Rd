\name{SBT2}
\alias{SBT2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{SBT complex MP}
\description{
An MP that makes incremental adjustments to TAC recommendations based on index levels relative to target levels (BMSY/B0) and catch levels relative to target levels (MSY)
}
\usage{
SBT2(x, DLM_data, reps = 100, epsB=0.25,epsR=0.75,tauR=5,tauB=7,gamma=1)
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
  \item{epsB}{
Control parameter
}
  \item{epsR}{
Control parameter
}
  \item{tauR}{
Control parameter
}
  \item{tauB}{
Control parameter
}
  \item{gamma}{
Control parameter
}
}
\details{
This isn't exactly the same as the proposed methods and is stochastic in this implementation. The method doesn't tend to work too well under many circumstances possibly due to the lack of 'tuning' that occurs in the real SBT assessment environment. You could try asking Rich Hillary at CSIRO about this approach. 
}

\references{
http://www.ccsbt.org/site/recent_assessment.php
}
\author{
T. Carruthers
}
