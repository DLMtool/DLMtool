\name{movfit}
\alias{movfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Optimization function that returns the squared difference between user specified and calculated movement parameters. }
\description{
The user specifies the probability of staying in the same area and spatial heterogeneity (both in the unfished state). This function returns the squared difference between these values and those produced by the three logit movement model.  
}
\usage{
movfit(par,prb,frac)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{
Three parameters in the logit space that control the four probabilities of moving between 2 areas
}
  \item{prb}{
User specified probability that individuals in area 1 remain in that area (unfished conditions)
}
  \item{frac}{
User specified fraction of individuals found in area 1 (unfished conditions)
}
}
\details{
This is paired with getmov to find the correct movement model. 
}
\author{
T. Carruthers
}
