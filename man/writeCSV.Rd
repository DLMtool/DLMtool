\name{writeCSV}
\alias{writeCSV}

\title{Internal function to write CSVs for objects}
\description{
Used internally in the DLMtool package to write CSV files from an existing DLMtool object
}
\usage{
writeCSV(inobj, tmpfile = NULL, objtype = c("Stock", "Fleet", "Observation", "DLM_data", 
  "OM", "DLM_fease"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inobj}{
A object of class Stock, Fleet, Observation, DLM_data, OM, or DLM_fease
}
  \item{tmpfile}{
The full file path and name for the saved CSV file
}
  \item{objtype}{
The class corresonding to the \code{inobj}
}
}
\author{
A. Hordyk
}
