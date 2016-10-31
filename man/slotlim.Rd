\name{slotlim}
\alias{slotlim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{An data-limited method which sets a slot limit}
\description{
An example of the implementation of input controls in the DLM toolkit, where selectivity-at-length is set using a slot limit; that is, a minimum and maximum legal length.  The maximum limit is set here, quite arbitrarily, as the 75th percentile between the new minimum legal length and the estimated asymptotic length.  
}
\usage{
slotlim(x, DLM_data, ...)
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
A vector of input control recommendations, with values for length at first capture, full selection, and maximum size limit in the 5th, 6th, and 7th elements of the vector
}
\references{
Made-up for this package
}
\author{
A. Hordyk
}
