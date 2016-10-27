\name{Splot}
\alias{Splot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scatter plot of B/BMSY or B/B0 and F/FMSY
}
\description{
Scatter plot of B/BMSY or B/B0 and F/FMSY 
}
\usage{
Splot(MSEobj = NULL, MPs = NA, All = TRUE, Var = c("B_BMSY", "B_B0"), 
  lastYrs = 10, Fref = 1, BMSYref = 1, B0ref = 0.4, cex.MP = 1, Fbg = FALSE, 
  Bbg = FALSE, Props = FALSE, TP = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MSEobj}{
An object of class MSE
}
  \item{MPs}{
Optional subset by MP 
}
  \item{All}{
Logical. Plot all points or just the mean?
}
  \item{Var}{
What to plot on the x-axis: \code{B_BMSY} or \code{B_B0}
}
  \item{lastYrs}{
Last number of years in projection to calculate statistics
}
  \item{Fref}{
Location of F statistic reference line
}
  \item{BMSYref}{
Location of B_MSY statistic reference line
}
  \item{B0ref}{
Location of B_0 statistic reference line
}
  \item{cex.MP}{
size of MP label
}
  \item{Fbg}{
Logical. Include background colours for F-statistic?
}
  \item{Bbg}{
Logical. Include background colours for B-statistic?
}
  \item{Props}{
Logical. Display the proportion of points in each quadrant?
}
  \item{TP}{
Logical. Use transparent colours?
}
}
\author{
A. Hordyk
}