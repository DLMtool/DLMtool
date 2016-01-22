\name{ITM}
\alias{ITM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Index Target based on natural mortality rate
}
\description{
An index target MP where the TAC is modified according to current index levels (mean index over last yrsmth years) relative to a target level. Maximum fractional annual changes are mc. 
mc=(5+M*25)/100
yrsmth=4*(1/M)^(0.25)
}
\usage{
ITM(x, DLM_data, reps = 100)
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

