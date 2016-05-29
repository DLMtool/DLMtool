\name{VOIplot}
\alias{VOIplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Yet another Value of Information Plot
}
\description{
A function that relates parameters of the observation model and the operating model parameters  
to yield. 
}
\usage{
VOIplot(MSEobj, MPs = NA, nvars = 5, nMP = 4, Par = c("Obs", "OM"), 
  doPlot = TRUE, incStat = FALSE, availMP = NULL, acceptMP = NULL, 
  incNames = TRUE, labcex = 0.8, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MSEobj}{
An object of class MSE
}
  \item{MPs}{
The MPs to plot. If NA it will plot the first nMP from MSEobj
}
  \item{nvars}{
The number of observation or operating model parameters to plot (number of columns)
}
  \item{nMP}{
The maximum number of MPs to plot (number of rows)
}
  \item{Par}{
Plot Operating Model (OM) or Observation (Obs) parameters?
}
  \item{doPlot}{
Output the plot?
}
  \item{incStat}{
Include a print out of statistic describing the curviness of the line?
}
  \item{availMP}{
Optional character string of MPs that are available. These names are colored black
}
  \item{acceptMP}{
Optional character string of MPs that are acceptable. These names are colored green if they are also in availMP
}
  \item{incNames}{
Include the names?
}
  \item{labcex}{
Character size of the label
}
  \item{\dots}{
Optional additional arguments
}
}

\value{
A list of all the information included in the plot
}
\author{
A. Hordyk
}
