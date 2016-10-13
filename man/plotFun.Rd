\name{plotFun}
\alias{plotFun}

\title{
Print out plotting functions
}
\description{
This function prints out the available plotting functions for objects of class MSE or DLM_data
}
\usage{
plotFun(class = c("MSE", "DLM_data"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{class}{
Character string. Prints out the plotting functions for objects of this class.
}
}

\author{
A. Hordyk
}
\note{
Basically the function looks for any functions in the DLMtool that have the word `plot` in them. 
There is a chance that some plotting functions are missed. Let us know if you find any and we will add them.
}