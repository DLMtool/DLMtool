\name{VOI2}
\alias{VOI2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate Value Of Information 2}
\description{
A function that relates operating model parameters and parameters of the observation model to relative yield (yield over last 5 years of projection relative to a 'best F' scenario that maximizes yield).
}
\usage{
VOI2(MSEobj, ncomp = 6, nbins = 4, Ut = NA, Utnam = "yield", lay = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MSEobj}{An object of class MSE}
  \item{ncomp}{Maximum number of observation variables to examine per MP}
  \item{nbins}{Number of bins for sampled observation variables used for calculating variability in utility across the sampled range of each parameter}
  \item{Ut}{A matrix of user-specified utility values of nsim rows and nMPs columns}
  \item{Utnam}{The name of the utility measure for plotting}
  \item{lay}{Controls whether labels are in lay terms or not}
}
\author{
T. Carruthers
}
\note{
VOI2 assumes that relative cost for each type of improvement in data is linearly related to the number of samples (e.g. nCAAobs) or square function of improved precision and bias e.g.: relative cost= 1/(newCV/oldCV)^2  
}
