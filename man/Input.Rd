\name{Input}
\alias{Input}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to run a set of input control methods
}
\description{
Runs a set of input control methods are returns the output in a single table
}
\usage{
Input(DLM_data, MPs = NA, reps = 100, timelimit = 10, CheckMPs = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DLM_data}{
 A DLM_data object
}
  \item{MPs}{
 A list of input MPs, if NA all available input MPs are run
}
  \item{reps}{
 Number of repetitions (for those methods that use them) 
}
  \item{timelimit}{
 Maximum timelimit to run MP (in seconds)
}
  \item{CheckMPs}{
 Logical, the Can function is run if this is TRUE
}
}

\author{
A. Hordyk
}

