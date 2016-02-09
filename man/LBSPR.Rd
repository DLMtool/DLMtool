\name{LBSPR}
\alias{LBSPR}
\title{
Apply LBSPR model to time-series of catch-at-length
}
\description{
Apply LBSPR model to time-series of catch-at-length data, and return smoothed estimates of SPR, and relative fishing mortality (F/M). 
}
\usage{
LBSPR(x, DLM_data, yrsmth = 1, perc = pstar, reps = reps)
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
Number of years to smooth length data - currently not used
}
  \item{perc}{
Currently not used
}
  \item{reps}{
Currently not used
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Hordyk, A.R., Ono, K., Sainsbury, K.J., Loneragan, N., and Prince, J.D. 2015. Some explorations of the life history ratios to describe length composition, spawning-per-recruit, and the spawning potential ratio. ICES J. Mar. Sci. 72: 204-216.

Hordyk, A.R., Ono, K., Valencia, S.R., Loneragan, N.R., and Prince, J.D. 2015. A novel length-based empirical estimation method of spawning potential ratio (SPR), and tests of its performance, for small-scale, data-poor fisheries. ICES J. Mar. Sci. 72: 217-231.
}
\author{
A. Hordyk 
}
