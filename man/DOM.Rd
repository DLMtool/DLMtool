\name{DOM}
\alias{DOM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
How dominant is an MP?
}
\description{
The DOM function examines how consistently an MP outperforms another. For example DCAC might provide higher yield than AvC on average but outperforms AvC in less than half of simulations. 
}
\usage{
DOM(MSEobj, MPtg=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MSEobj}{
An object of class 'MSE'
}
  \item{MPtg}{
A character vector of management procedures for cross examination
}
}
\value{
A matrix of performance comparisons length(MPtg) rows by MSE@nMPs columns
}
\author{
A. Hordyk
}

