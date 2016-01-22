\name{ML2D}
\alias{ML2D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Depletion and F estimation from mean length of catches
}
\description{
A highly dubious means of getting very uncertain estimates of current stock biomass and (equilibrium) fishing mortality rate from growth, natural mortality rate, recruitment and fishing selectivity. 
}
\usage{
ML2D(OM,ML,nsim=100,ploty=T,Dlim=c(0.05,0.6))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{OM}{
An object of class 'OM'
}
  \item{ML}{
A estimate of current mean length of catches
}
  \item{nsim}{
Number of simulations
}
  \item{ploty}{
Produce a plot of depletion and F
} 
  \item{Dlim}{
Limits on the depletion that is returned as a fraction of unfished biomass. 
}
}
\value{
A table of nsim rows and 2 columns (depletion, fishing mortality rate)
}
\author{
T. Carruthers
}

