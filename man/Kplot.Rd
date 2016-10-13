\name{Kplot}
\alias{Kplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
KOBE plot: a projection by projection plot of F/FMSY and B/BMSY
}
\description{
A standard KOBE plot by each method that also shows the percentage of methods that ended up in each quadrant.
}
\usage{
Kplot(MSEobj,maxsim=60, MPs=NA, sims=NULL, maxMP=9, nam=NA, cex.leg=1.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MSEobj}{
An object of class MSE
}
\item{maxsim}{
Maximum number of simulations (lines) to plot on each panel. 
}
\item{MPs}{
Optional subset MSE object by MP
}
\item{sims}{
Optional subset MSE object by simulation
}
\item{maxMP}{
Maximum number of MPs to include in plot
}
\item{nam}{
The name of the plot
}
\item{cex.leg}{
Size of legend
}
}
\author{
T. Carruthers with some additions from A. Hordyk
}
\note{
Apologies for the nauseating shading. 
}
