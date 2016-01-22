\name{FMSYref75}
\alias{FMSYref75}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A reference FMSY method that fishes at three quarters of FMSY (uses perfect information about FMSY)
}
\description{
FMSY is taken from the operating model stored at DLM@OM$FMSY
}
\usage{
FMSYref75(x, DLM_data, reps = 100)
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
Note that you can out-performm this method easily. The requirment for fixed F is actually quite strict and is by no means the upper limit in terms of yield. Don't panic if your method beats this one for yield!

Interesting that the reduction in yield is no way near commensurate with the reduction in F as predicted by a yield curve and expressed in the pretty good yield theory. 
}
\author{
T. Carruthers
}
