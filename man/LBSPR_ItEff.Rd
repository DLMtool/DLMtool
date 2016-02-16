\name{LBSPR_ItEff}
\alias{LBSPR_ItEff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Length-based SPR model with HCR that iteratively adjusts Effort
}
\description{
Iteratively adjusts Effort based on distance between estimated and target SPR (40\%), and slope of recent SPR estimates. 
}
\usage{
LBSPR_ItEff(x, DLM_data, yrsmth = 1, reps = reps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Simulation number
}
  \item{DLM_data}{
DLM_data object
}
  \item{yrsmth}{
Number of years to smooth length data - not currently used
}

  \item{reps}{
Not currently used
}
}

\author{
A. Hordyk
}
