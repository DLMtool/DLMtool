\name{LBSPR_ItTAC}
\alias{LBSPR_ItTAC}

\title{
Length-based SPR model with HCR that iteratively adjusts TAC
}
\description{
Iteratively adjusts TAC based on distance between estimated and target SPR (40\%), and slope of recent SPR estimates.
}
\usage{
LBSPR_ItTAC(x, DLM_data, yrsmth = 1, reps = reps)
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
