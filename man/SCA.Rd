\name{SCA}
\alias{SCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust Statistical Catch at Age using Stock Synthesis 3 (currently not working)
}
\description{
A simple SCA with one fleet (that has single CPUE index and catch-age composition data for every historical year) which estimates annual recruitment deviations and logistic selectivity. 
}
\usage{
SCA(x, DLM_data, reps = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A position in a data-limited methods data object
}
  \item{DLM_data}{
A data-limited methods data object
}
  \item{reps}{
The number of stochastic samples of the TAC recommendation
}
}
\value{
A numeric vector of TAC recommendations
}
\references{
Stock Sythesis 3. Methot, R. 2012. NOAA Fisheries Toolbox: http://nft.nefsc.noaa.gov/Stock_Synthesis_3.htm
}
\author{
T. Carruthers (DLMtool implementation)
}
\note{
This is an early version of SCA and currently doesn't work. The code is put up here to aid in collaboration to sort this out!
}
