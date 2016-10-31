\name{Ltarget1}
\alias{Ltarget1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A management procedure that incrementally adjusts the TAC to reach a target mean length in catches.}
\description{
The least biologically precautionary of four target length MPs proposed by Geromont and Butterworth 2014. Tested by Carruthers et al. 2015
}
\usage{
Ltarget1(x, DLM_data, reps = 100, yrsmth = 5, xx=0, xL=1.05)
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
\item{xL}{
Parameter controlling the magnitude of the target mean length of catches relative to average length in catches.
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
