\name{DDe}
\alias{DDe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Effort control version of DD - Delay - Difference Stock Assessment with UMSY and MSY leading}
\description{
A simple delay-difference assessment that estimates and recommends FMSY using a time-series of catches and a relative abundance index. 
}
\usage{
DDe(x, DLM_data, reps = 100)
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
\references{
Method based on equations of Carl Walters (bug him with questions and expect colourful responses) 
}
\author{
T. Carruthers
}
\note{
This DD model is observation error only and has does not estimate process error (recruitment deviations). Similar to many other assessment models it depends on a whole host of dubious assumptions such as temporally stationary productivity and proportionality between the abundance index and real abundance. Unsurprisingly the extent to which these assumptions are violated tends to be the biggest driver of performance for this method. 
}
