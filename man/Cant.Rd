\name{Cant}
\alias{Cant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
What methods can't be applied to this DLM data object
}
\description{
The methods that don't have sufficient data, lead to errors or don't run in time along with a list of their data requirments. 
}
\usage{
Cant(DLM_data, timelimit = 1)
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
