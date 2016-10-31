\name{ItargetE4}
\alias{ItargetE4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A management procedure that incrementally adjusts the Effort to reach a target CPUE / relative abundance index}
\description{
An effort-based version of the most biologically precautionary of two index/CPUE target MPs proposed by Geromont and Butterworth 2014.
}
\usage{
ItargetE4(x, DLM_data, reps = 100, yrsmth = 5, xx = 0, Imulti = 2.5)
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
\item{Imulti}{
Parameter controlling how much larger target CPUE / index is compared with recent levels.
}
}
\details{
Tested by Carruthers et al. 2015. 
}
\value{
A numeric vector of input controls
}
\references{
Carruthers et al. 2015. Performance evaluation of simple management procedures. Fish and Fisheries. In press.

Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for data-poor fisheries; forecasting with few data. ICES J. Mar. Sci. doi:10.1093/icesjms/fst232
}
\author{
T. Carruthers
}
