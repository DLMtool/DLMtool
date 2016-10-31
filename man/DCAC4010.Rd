\name{DCAC4010}
\alias{DCAC4010}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Depletion Corrected Average Catch paired with the 40-10 rule}
\description{
A method of calculating an MSY proxy (FMSY * BMSY and therefore the OFL at most productive stock size) based on average catches accounting for the windfall catch that got the stock down to BMSY levels. In this method DCAC is paired with the 40-10 rule that throttles back the OFL to zero at 10 percent of unfished stock size (the OFL is not subject to downward adjustment above 40 percent unfished)
}
\usage{
DCAC4010(x, DLM_data, reps = 100)
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
DCAC can overfish below BMSY levels. The 40-10 harvest control rule largely resolves this problem providing an MP with surprisingly good performance even at low stock levels. 
}
