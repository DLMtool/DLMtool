\name{Can}
\alias{Can}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{What data-limited methods can be applied to this DLM_data object?}
\description{
An diagnostic tool that looks up the slot requirements of each method and compares this to the data available to limit the analysis to methods that have the correct data, do not produce errors and run within a time limit. Time limit is the maximum time taken to carry out five reps (stochastic samples) of a given method and is in units of seconds. 
}
\usage{
Can(DLM_data, timelimit = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DLM_data}{
A data-limited methods data object (class DLM_data)
}
  \item{timelimit}{
The maximum time (seconds) taken for a method to undertake 10 reps (this filters out methods that are too slow)
}
}

