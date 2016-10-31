\name{Rcontrol}
\alias{Rcontrol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Harvest Control Rule using prior for intrinsic rate of increase}
\description{
An MP proposed by Carl Walters that modifies TACs according to trends in apparent surplus production that includes information from a demographically derived prior for intrinsic rate of increase
}
\usage{
Rcontrol(x, DLM_data, reps = 100, yrsmth = 10, gg = 2, glim = c(0.5, 2))
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
The number of years for smoothing catch and biomass data
}
  \item{gg}{
A gain parameters
}
  \item{glim}{
Limits for the change in TAC among years
}
}
\references{
Made-up for this package.
}
\author{
C. Walters and T. Carruthers
}
