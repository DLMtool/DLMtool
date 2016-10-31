\name{matlenlim}
\alias{matlenlim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A data-limited method in which fishing vulnerability is set according to the maturity curve}
\description{
An example of the implementation of input controls in the DLM toolkit, where selectivity-at-length is set equivalent to maturity-at-length
}
\usage{
matlenlim(x, DLM_data, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A position in a data-limited methods object
}
  \item{DLM_data}{
A data-limited methods object
}
\item{...}{Optional additional arguments that are ignored. Note arguments \code{reps} or \code{...} are required for all input controls }
}

\value{
A vector of input control recommendations, with values for length at first capture and full selection
}
\references{
Made-up for this package
}
\author{
T. Carruthers
}
