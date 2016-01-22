\name{DCAC}
\alias{DCAC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Depletion Corrected Average Catch
}
\description{
A method of calculating an MSY proxy (FMSY * BMSY and therefore the OFL at most productive stock size) based on average catches accounting for the windfall catch that got the stock down to BMSY levels.
}
\usage{
DCAC(x, DLM_data, reps = 100)
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
}
\author{
T. Carruthers
}
\note{
It's probably worth noting that DCAC TAC recomemndations do not tend to zero as depletion tends to zero. It adjusts for depletion only in calculating historical average catch. It follows that at stock levels much below BMSY, DCAC tends to chronically overfish. 
}
