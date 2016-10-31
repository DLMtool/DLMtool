\name{Rcontrol2}
\alias{Rcontrol2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MP using prior for intrinsic rate of increase with a quadratic approximation to surplus production}
\description{
An MP proposed by Carl Walters that modifies quotas according to trends in apparent surplus production that includes information from a demographically derived prior for intrinsic rate of increase. This is different from Rcontrol because it includes a quadratic approximation of recent trend in surplus production given biomass
}
\usage{
Rcontrol2(x, DLM_data, reps = 100, yrsmth = 10, gg = 2, glim = c(0.5, 2))
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
