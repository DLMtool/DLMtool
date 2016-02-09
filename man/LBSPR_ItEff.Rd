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
LBSPR_ItEff(x, DLM_data, yrsmth = 1, perc = pstar, reps = reps)
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
  \item{perc}{
Not currently used
}
  \item{reps}{
Not currently used
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
%% ~put references to the literature/web site here ~
}
\author{
A. Hordyk
}
