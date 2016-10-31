\name{VOI}
\alias{VOI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate Value Of Information}
\description{
A function that relates operating model parameters and parameters of the observation model to yield (by default). A user can also specific their own utility values (Ut) which is arranged in a matrix of nsim rows and nMP columns. 
}
\usage{
VOI(MSEobj, ncomp = 6, nbins = 8, maxrow = 8, Ut = NA, Utnam = "Utility")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MSEobj}{An object of class MSE}
  \item{ncomp}{Maximum number of variables to examine per MP}
  \item{nbins}{Number of percentile bins for sampled parameters of the operating model or observation model, which is used for calculating variability in utility across the sampled range of each parameter}
  \item{maxrow}{maximum number of MPs per plot}
  \item{Ut}{A matrix of user-specified utility values of nsim rows and nMPs columns}
  \item{Utnam}{The name of the utility measure for plotting}
}
\author{
T. Carruthers
}
