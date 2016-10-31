\name{getq}
\alias{getq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Optimization function that find the catchability (q where F=qE) value required to get to user-specified stock depletion (current biomass / unfished biomass)}
\description{
The user specifies the level of stock depleiton. This funciton takes the derived effort trajectories and finds the catchabiltiy to get the stock there.   
}
\usage{
getq(x,dep,Find,Perr,Marray,hs,
     Mat_age,Wt_age,R0,V,nyears,
     maxage,mov,Spat_targ,SRrel,
     aR,bR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
}
  \item{dep}{
}
  \item{Find}{
}
  \item{Perr}{
}
  \item{Marray}{
}
  \item{hs}{
}
  \item{Mat_age}{
}
  \item{Wt_age}{
}
  \item{R0}{
}
  \item{V}{
}
  \item{nyears}{
}
  \item{maxage}{
}
  \item{mov}{
}
  \item{Spat_targ}{
}
  \item{SRrel}{
}
  \item{aR}{
}
  \item{bR}{
}
}
\details{
Paired with qopt. 
}
\author{
T. Carruthers
}
