\name{makePerf}
\alias{makePerf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert a OM object to one without observation or process error
}
\description{
Takes an existing OM object and converts it to one without any observation error, and 
very little process error.  Used for debugging and testing that MPs perform as expected 
under perfect conditions.
}
\usage{
makePerf(OMin, except = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{OMin}{
An object of class \code{OM}
}
  \item{except}{
An optional vector of slot names in the OM that will not be changed (not tested perfectly so watch out!)
}
}

\value{
A new \code{OM} object
}

\author{
A. Hordyk
}

