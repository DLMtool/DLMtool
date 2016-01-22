\name{LstepCC1}
\alias{LstepCC1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A management procedure that incrementally adjusts the TAC according to the mean length of recent catches.
}
\description{
The least biologically precautionary of four adaptive length-based MPs proposed by Geromont and Butterworth 2014. Tested by Carruthers et al. 2015
}
\usage{
LstepCC1(x, DLM_data, reps = 100, yrsmth = 5, xx=0, stepsz=0.05, llim=c(0.96,0.98,1.05))
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
Years over which to smooth recent estimates of surplus production  
}
\item{xx}{
Parameter controlling the fraction of mean catch to start using in first year
}
\item{stepsz}{
Parameter controlling the size of the TAC update increment.
}
\item{llim}{
A vector of length reference points that determine the conditions for increasing, maintaining or reducing the TAC.
}
}
\details{
Tested by Carruthers et al. 2015. 
}
\value{
A numeric vector of TAC recommendations
}
\references{
Carruthers et al. 2015. Performance evaluation of simple management procedures. Fish and Fisheries. In press.

Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for data-poor fisheries; forecasting with few data. ICES J. Mar. Sci. doi:10.1093/icesjms/fst232
}
\author{
T. Carruthers
}
