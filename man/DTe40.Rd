\name{DTe40}
\alias{DTe40}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Effort searching MP aiming for 40 per cent stock depletion}
\description{
A very simple MP that modifies effort to reach 40 percent stock depletion
}
\usage{
DTe40(x, DLM_data, reps = 100, alpha=0.4, LB=0.9, UB=1.1)
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
The number of stochastic samples of the TAC recommendation
}
  \item{alpha}{
The target level of depletion
}
\item{LB}{
The lowest permitted factor of previous fishing effort
}
\item{UB}{
The highest permitted factor of previous fishing effort
}
}

\author{
T. Carruthers
}

