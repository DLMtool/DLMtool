\name{EtargetLopt}
\alias{EtargetLopt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Effort MP: adjust effort up/down if mean length above/below Ltarget 
}
\description{
Effort MP: adjust effort up/down if mean length above/below Ltarget 
}
\usage{
EtargetLopt(x, DLM_data, reps = 100, yrsmth=3)
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
\item{yrsmth}{
Number of years to calculate average length
}
}
\value{
An adjustment for fishing effort
}
\author{
HF Geromont
}
