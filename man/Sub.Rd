\name{Sub}
\alias{Sub}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subset MSE object by management procedure (MP) or simulation. 
}
\description{
 Subset the MSE object by particular MPs (either MP number or name), or 
 particular simulations, or a subset of the projection years (e.g., 1: < projection years). 
}
\usage{
Sub(MSEobj, MPs=NULL, sims=NULL, years=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MSEobj}{
A MSE object.
}
\item{MPs}{
A vector MPs names or MP numbers to subset the MSE object. 
 Defaults to all MPs.
}
  \item{sims}{
A vector of simulation numbers to subset the MSE object. Can also be a 
logical vector. Defaults to all simulations.
}
\item{years}{
A numeric vector of projection years. Should start at 1 and increase by one to some value equal or less than the total number of projection years.
}
}
\author{
A. Hordyk
}
