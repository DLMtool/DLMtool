\name{Cplot}
\alias{Cplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the median biomass and yield relative to last historical year}
\description{
Compare median biomass and yield in first year and last 5 years of projection
}
\usage{
Cplot(MSEobj, MPs = NA, lastYrs = 5, XMin = NULL, YMin = NULL, ShowLabs = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MSEobj}{
An object of class MSE
}
  \item{MPs}{
Optional subset by MP
}
  \item{lastYrs}{
Last number of years of projection to calculate median
}
  \item{XMin}{
Optional minimum for the x-axis
}
  \item{YMin}{
Optional minimum for the y-axis
}
  \item{ShowLabs}{
Logical. Show the MP labels? Otherwise only plot points
}
}

\value{
Invisibly returns a data frame containing information shown in the plot
}
\author{
A. Hordyk
}