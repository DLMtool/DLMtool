\name{runInMP}
\alias{runInMP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Runs input control MPs on a DLM_data object.
}
\description{
Function runs a MP (or MPs) of class 'DLM_input' and returns a list: input control recommendation(s) in element 1 and DLM_data object in element 2.
}
\usage{
runInMP(DLM_data, MPs = NA, reps = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DLM_data}{
 A object of class DLM_data
}
  \item{MPs}{
A vector of MPs of class 'DLM_input'
}
  \item{reps}{
Number of stochastic repititions - often not used in input control MPs.
}
}

\author{
A. Hordyk
}


