\name{NFref}
\alias{NFref}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
No Fishing Reference MP
}
\description{
A reference MP that sets annual catch to zero (or very close to it). Used for looking at variability in stock with no fishing.
}
\usage{
NFref(x, DLM_data, reps = 100)
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
A TAC of 0.01
}
\author{
A. Hordyk
}
