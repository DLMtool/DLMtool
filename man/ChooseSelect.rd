\name{ChooseSelect}
\alias{ChooseSelect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Manually choose the historical selectivity pattern.
}
\description{
 Input the first historical year, and all years where selectivity pattern changed (separated by comma).  
 
 Interactive plot which allows users to specify a range for the length at 5\% and full selection (LFS), as well as selectivity at maximum length for each year.  
 Produces a simple plot which shows the range in selectivity pattern for each break-point year. Selectivity-at-length is fixed in between break-point years.
 Note that this function replaces 'nyears' in the Fleet object with the value defined here (FstYr:current year).
}
\usage{
ChooseSelect(Fleet, Stock=NULL, FstYr=NULL, SelYears=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Fleet}{
  A fleet object.
  }
  \item{Stock}{
  Optional Stock object. If provided, average length-at-maturity is included on plot for reference.
  }
  \item{FstYr}{
  Optional value for first historical year. If empty, user must specify the year in console.   
  }
  \item{SelYears}{
  Optional vector of values for each year where selectivity pattern changed. If empty, user must specify the years in console (comma separated).   
  }
}
\author{
A. Hordyk
}


