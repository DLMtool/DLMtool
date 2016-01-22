\name{GB_target}
\alias{GB_target}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Geromont and Butterworth target CPUE and catch MP
}
\description{
An MP similar to SBT2 that modifies a time-series of catch recommendations and aims for target catch rate and catch level based on BMSY/B0 and MSY, respectively.
}
\usage{
GB_target(x, DLM_data, reps = 100, w = 0.5)
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
}  \item{w}{
A gain parameter
}
}
\details{
Note that this is my interpretation of their MP and is now stochastic. Currently it is generalized and is not 'tuned' to more detailed assessment data which might explain why in some cases it leads to stock declines. 
}
\references{
Geromont, H.F. and Butterworth, D.S. 2014. Complex assessment or simple management procedures for efficient fisheries management: a comparative study. ICES J. Mar. Sci. doi:10.1093/icesjms/fsu017
}
\author{
T. Carruthers
}
