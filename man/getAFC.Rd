\name{getAFC}
\alias{getAFC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate age at first capture from length at first capture and growth}
\description{
As title.  
}
\usage{
getAFC(t0c,Linfc,Kc,LFC,maxage)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t0c}{
    A vector of theoretical age at length zero (von Bertalanffy growth)
  }
  \item{Linfc}{
    A vector of maximum length (von Bertalanffy growth)
  }
  \item{Kc}{
    A vector of maximum growth rate (von Bertalanffy growth)
  }
  \item{LFC}{
    A vector of length at first capture
  }
  \item{maxage}{
    Maximum age
  }
}
\author{
T. Carruthers
}

