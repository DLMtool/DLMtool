\name{boxplot.DLM_data}
\alias{boxplot.DLM_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Boxplot of TAC recommendations}
\description{
Boxplot of TAC recommendations from a DLM_data object
}
\usage{
\method{boxplot}{DLM_data}(x, outline = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class MSE
}
  \item{outline}{
Logical. Include outliers in plot?
}
  \item{...}{
Optional additional arguments passed to \code{boxplot}
}
}
\value{
Returns a data frame containing the information shown in the plot
}
\author{
A. Hordyk
}