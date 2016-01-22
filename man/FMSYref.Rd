\name{FMSYref}
\alias{FMSYref}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A reference FMSY method (uses perfect information about FMSY)
}
\description{
FMSY is taken from the operating model stored at DLM@OM$FMSY
}
\usage{
FMSYref(x, DLM_data, reps = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A position in data-limited methods data object
}
  \item{DLM_data}{
A data-limited methods data object
}
  \item{reps}{
The number of TAC samples
}
}
\details{
Note that you can out-performm this MP even though it has perfect information of FMSY and current abundance. The requirment for fixed F is actually quite strict and is by no means the upper limit in terms of yield. Don't panic if your method beats this one for yield, especially for short-lived species of high temporal variability in productivity!
}
\author{
T. Carruthers
}
