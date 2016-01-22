\name{Sub}
\alias{Sub}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subset MSE object by management procedure (MP) or simulation. 
}
\description{
 Subset the MSE object by particular MPs (either MP number or name), or 
 particular simulations. 
}
\usage{
Sub(MSEobj, MPs=NULL, sims=NULL)
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
A vector of simulation numbers to subset the MSE object. can also be a 
logical vector. Defaults to all simulations.
}
}
\author{
A. Hordyk
}
