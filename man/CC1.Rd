\name{CC1}
\alias{CC1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Constant catch management procedure of Geromont and Butterworth (2014)}
\description{
The TAC is the average catch over last yrsmth years.
}
\usage{
CC1(x, DLM_data, reps = 100, yrsmth = 5, xx=0)
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
Years over which to calculate mean catches  
}
\item{xx}{
Parameter controlling the TAC. Mean catches are multiplied by (1-xx)
}
}
\details{
This is one of four constant catch rules of Geromont and Butterworth 2014.
}
\value{
A numeric vector of TAC recommendations
}
\references{
Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for data-poor fisheries; forecasting with few data. ICES J. Mar. Sci. doi:10.1093/icesjms/fst232
}
\author{
T. Carruthers
}
