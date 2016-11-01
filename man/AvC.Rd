\name{AvC}
\alias{AvC}
\title{An MP where the TAC is average historical catch}
\description{
A simple average catch MP that is included to demonstrate a 'status quo' management option.
}
\usage{
AvC(x, DLM_data, reps = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A position in a data-limited methods data object
}
  \item{DLM_data}{
A data-limited methods data object
}
  \item{reps}{
The number of stochastic samples of the quota recommendation
}
}
\value{
A numeric vector of TAC recommendations
}
\author{
T. Carruthers
}

