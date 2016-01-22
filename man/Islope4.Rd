\name{Islope4}
\alias{Islope4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A management procedure that incrementally adjusts the TAC to maintain a constant CPUE or relative abundance index.
}
\description{
The most biologically precautionary of two constant index / CPUE methods proposed by Geromont and Butterworth 2014. Tested by Carruthers et al. 2015
}
\usage{
Islope4(x, DLM_data, reps = 100, yrsmth = 5, lambda=0.2,xx=0.4)
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
\item{lambda}{
A gain parameter controlling the speed in update in TAC.
}
\item{xx}{
Parameter controlling the fraction of mean catch to start using in first year
}

}
\details{
Tested by Carruthers et al. 2015. 
}
\value{
A numeric vector of quota recommendations
}
\references{
Carruthers et al. 2015. Performance evaluation of simple management procedures. Fish and Fisheries. In press.

Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for data-poor fisheries; forecasting with few data. ICES J. Mar. Sci. doi:10.1093/icesjms/fst232
}
\author{
T. Carruthers
}
