\name{L95target}
\alias{L95target}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A management procedure that adjusts the TAC up/down from reference (target) level (that is a fraction of mean recent premanagement catches) to reach a target mean length of fish caught.
}
\description{
This MP is based on Ltarget1 proposed by Geromont and Butterworth 2014, but here the target and limit mean lengths are based on the length at maturity distribution rather than an arbitrary multiplicative of the mean length. 
}
\usage{
L95target(x, DLM_data, reps = 100, yrsmth = 5, buffer=0)
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
Years over which to calculate the mean historical catch
}
\item{buffer}{
Parameter controlling the fraction of mean catch to set the reference (or target) TAC level - acts as a precautionary buffer 
}
}
\value{
A numeric vector of TAC recommendations
}
\references{
Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for data-poor fisheries; forecasting with few data. ICES J. Mar. Sci. doi:10.1093/icesjms/fst232
}
\author{
HF Geromont
}
