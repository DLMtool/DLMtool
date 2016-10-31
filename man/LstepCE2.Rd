\name{LstepCE2}
\alias{LstepCE2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A management procedure that incrementally adjusts the Effort according to the mean length of recent catches.}
\description{
A effort-based version of one of the four adaptive length-based MPs proposed by Geromont and Butterworth 2014.
}
\usage{
LstepCE2(x, DLM_data, reps = 100, yrsmth = 5, xx = 0, stepsz = 0.1, 
  llim = c(0.96, 0.98, 1.05))
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
The number of samples
}
  \item{yrsmth}{
Years over which to smooth recent estimates of surplus production  
}
\item{xx}{
Parameter controlling the fraction of mean catch to start using in first year
}
\item{stepsz}{
Parameter controlling the size of the effort update increment.
}
\item{llim}{
A vector of length reference points that determine the conditions for increasing, maintaining or reducing the effort.
}
}
\value{
A numeric vector of input controls
}
\author{
T. Carruthers
}
