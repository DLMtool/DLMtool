\name{Gcontrol}
\alias{Gcontrol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
G-control MP
}
\description{
A harvest control rule proposed by Carl Walters that uses trajectory in inferred surplus production to make upward/downward adjustments to TAC recommendations
}
\usage{
Gcontrol(x, DLM_data, reps = 100, yrsmth = 10, gg = 2, glim = c(0.5, 2))
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
The number of quota samples
}
\item{yrsmth}{
The number of years over which to smooth catch and biomass data
}
  \item{gg}{
A gain parameter
}
  \item{glim}{
A constraint limiting the maximum level of change in quota recommendations
}
}
\references{
Made-up for this package. Carruthers et al. 2015. Performance of Simple Management Procedures. 
}
\author{
C. Walters and T. Carruthers
}
