\name{Fease}
\alias{Fease}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MP feasibility diagnostic}
\description{
What MPs may be run (best case scenario) for various data-availability scenarios?
}
\usage{
Fease(feaseobj,outy="table")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{feaseobj}{
An object of class 'DLM_fease'
}
\item{outy}{
Determines whether you would like a full table or some column of the table for a specific case of the feasibility object. When set equal to table, the full table is produced. When set equal to an integer number the names of MPs that are feasible for that case are returned. 
}
}
\author{
T. Carruthers
}

