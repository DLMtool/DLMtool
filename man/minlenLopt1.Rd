\name{minlenLopt1}
\alias{minlenLopt1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This input control sets the minimum length of fish caught to a fraction of the length that maximises the biomass, Lopt. 
}
\description{
This aim of this simple MP is restrict the catch of small fish to rebuild the stock biomass towards the optimal length, Lopt, expressed in terms of the growth parameters Lopt=b/(M/k+b) (Hordyk et al. (2014)
}
\usage{
minlenLopt1(x, DLM_data, reps = 100, buffer = 0.1)
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
\item{buffer}{
Parameter controlling the fraction of Lopt to set the minimum length of fish caught: minlen=Lopt*(0.7+buffer).
}
}
\value{
The length at first caprture, LFC, and length at full selectivity
}
\references{
Hordyk, A., Ono, K., Sainsbury, K., Loneragan, N., and J. Prince. 2014. Some explorations of the life history ratios to describe length composition, spawning-per-recruit, and the spawning potential ratio ICES Journal of Marine Science, doi:10.1093/icesjms/fst235.	
}
\author{
HF Geromont
}
