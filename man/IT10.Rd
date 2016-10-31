\name{IT10}
\alias{IT10}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Index Target 10}
\description{
An index target MP where the TAC is modified according to current index levels (mean index over last 5 years) relative to a target level. Maximum annual changes are 10 per cent. 
}
\usage{
IT10(x, DLM_data, reps = 100,yrsmth=5,mc=0.1)
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
The number of stochastic samples of the quota recommendation
}
\item{yrsmth}{
The number of historical years over which to average the index
}
\item{mc}{
The maximum fractional change in the TAC among years.  
}
}
\value{
A numeric vector of TAC recommendations
}
\author{
T. Carruthers
}

