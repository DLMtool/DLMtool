\name{MPStats}
\alias{MPStats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Statistics for MP Performance
}
\description{
Function calculates probabilities and other statistics for a range of performance metrics
}
\usage{
MPStats(MSEobj, 
  PMRefs = list(B_BMSY = 0.5, B_B0 = 0.2, F_FMSY = 1, AAVY = 30, AAVE = 30), 
  lastYrs = 10, UseMean = TRUE, msg = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MSEobj}{
An object of class MSE
}
  \item{PMRefs}{
A list of reference points for the performance metrics (must be named)
}
  \item{lastYrs}{
The last number of years in the projection to calculate the statistics
}
  \item{UseMean}{
Logical. Calculate mean (TRUE) or median (FALSE)?
}
  \item{msg}{
Logical. Print out messages?
}
}
\author{
A. Hordyk
}