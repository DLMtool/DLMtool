\name{EDCAC}
\alias{EDCAC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extra Depletion Corrected Average Catch
}
\description{
Essentially DCAC multiplied by 2*depletion and divided by BMSY/B0 (Bpeak)
}
\usage{
EDCAC(x, DLM_data, reps = 100)
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
}
\references{
MacCall, A.D., 2009. Depletion-corrected average catch: a simple formula for estimating sustainable yields in data-poor situations. ICES J. Mar. Sci. 66, 2267-2271.
Harford W. and Carruthers, T. 2015. A simple modification to MacCAll's depletion-corrected average catch. In draft, intended for Fish. Bull. 
}
\author{
W. Harford and T. Carruthers
}

