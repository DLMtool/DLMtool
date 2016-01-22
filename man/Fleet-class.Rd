\name{Fleet-class}
\Rdversion{1.1}
\docType{class}
\alias{Fleet-class}
%%\alias{initialize,Fleet-method}

\title{Class \code{"Fleet"}}
\description{
The component of the operating model that controls fishing dynamics
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Fleet", OM)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{Name}:}{Name of the Fleet object}
    \item{\code{nyears}:}{The number of years for the historical simulation}
    \item{\code{Spat_targ}:}{Distribution of fishing in relation to spatial biomass: F is proportional to B^Spat_targ (uniform distribution) }
    \item{\code{LFS}:}{Shortest length that is fully vulnerable to fishing (uniform distribution)}
    \item{\code{L5}:}{Shortest length corresponding ot 5 percent vulnerability (uniform distribution)}
    \item{\code{Vmaxlen}:}{The vulnerability of the longest (oldest) fish (uniform distribution)}
    \item{\code{Fsd}:}{Inter-annual variability in fishing mortality rate}
    \item{\code{EffYears}:}{Vector of verticies, years at which to simulate varying relative effort}
    \item{\code{EffLower}:}{Lower bound on relative effort corresponding to EffYears (uniform distribution)}
    \item{\code{EffYears}:}{Uppper bound on relative effort corresponding to EffYears (uniform distribution)}
    \item{\code{qinc}:}{Average percentage change in fishing efficiency (uniform distribution)(applicable only to forward projection and input controls)}
    \item{\code{qcv}:}{Inter-annual variability in fishing efficiency (uniform distribution)(applicable only to forward projection and input controls)}
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "Fleet")}: ... }
	 }
}
\author{
T. Carruthers
}
\examples{
showClass("Fleet")
}
\keyword{classes}
