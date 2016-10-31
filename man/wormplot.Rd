\name{wormplot}
\alias{wormplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Biomass wormplot}
\description{
A worm plot for plotting the likelihood of meeting biomass targets in future years.
}
\usage{
wormplot(MSEobj, Bref=0.5, LB=0.25, UB=0.75)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MSEobj}{ Object of class MSE, output of the runMSE function}
\item{Bref}{The reference fraction of BMSY (to evaluate the probability of exceeding this level)}
\item{LB}{The lower bound probability that seperates red (bad) and yellow (O.K.) colored segments}
\item{UB}{The upper bound probability that seperates yellow (O.K.) and green (good) colored segments}
}
\details{
Returns a matrix of nMPs rows and proyears columns which is the fraction of simulations for which biomass was above Bref. 
}

\author{
T. Carruthers
}

