\name{L2A}
\alias{L2A}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Length to age conversion}
\description{
Simple deterministic length to age conversion given inverse von Bertalanffy growth.
}
\usage{
L2A(t0c,Linfc,Kc,Len,maxage)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t0c}{
Theoretical age at length zero
}
  \item{Linfc}{
Maximum length
}
  \item{Kc}{
Maximum growth rate
}
  \item{Len}{
Length
} 
  \item{maxage}{
Maximum age
}
}
\value{
An age (vector of ages, matrix of ages) corresponding with Len
}
\author{
T. Carruthers
}

