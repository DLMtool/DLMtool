\name{DBSRA}
\alias{DBSRA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Depletion-Based Stock Reduction Analysis
}
\description{
User prescribed BMSY/B0, M, FMSY/M are used to find B0 and therefore the OFL by back-constructing the stock to match a user specified level of stock depletion (OFL = M * FMSY/M * depletion* B0).
}
\usage{
DBSRA(x, DLM_data, reps = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A position in a data-limited methods object.
}
  \item{DLM_data}{
A data-limited methods object.
}
  \item{reps}{
The number of samples of the TAC (OFL) recommendation.
}
}
\details{
You specify a range of stock depletion and, given historical catches DB-SRA calculates what unfished biomass must have been to get you here given samples for M, FMSY relative to M and also BMSY relative to Bunfished. 
}
\value{
A vector of TAC (OFL) values.
}
\references{
Dick, E.J., MacCall, A.D., 2011. Depletion-Based Stock Reduction Analysis: A catch-based method for determining sustainable yields for data-poor fish stocks. Fish. Res. 110, 331-341.
}
\author{
T. Carruthers
}
\note{
This is set up to return the OFL (FMSY * current biomass).

You may  have noticed that you -the user- specify three of the factors that make the quota recommendation. So this can be quite a subjective method.

Also the DB-SRA method of this package isn't exactly the same as the original method of Dick and MacCall (2011) because it has to work for simulated depletions above BMSY/B0 and even on occasion over B0. Also it doesn't have the modification for flatfish life histories that has previously been applied by Dick and MacCall.
}

