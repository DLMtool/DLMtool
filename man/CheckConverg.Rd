\name{CheckConverg}
\alias{CheckConverg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MSE convergence diagnostic}
\description{
Have I undertaken enough simulations (nsim)? Has my MSE converged on stable (reliable) peformance metrics?
}
\usage{
CheckConverg(MSEobj,thresh=2, Plot=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MSEobj}{
An object of class 'MSE'
}
\item{thresh}{
The convergence threshold (percentage). If mean perforamnce metrics are within thresh percent of the second to last interation, the MSE can be considered to have converged.
}
\item{Plot}{
Should figures be plotted?
}
}
\author{
A. Hordyk
}

