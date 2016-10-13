\name{Pplot}
\alias{Pplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A projection by projection plot of F/FMSY and B/BMSY 
}
\description{
A shorter version of the plot method for MSEs that just shows the projected trends in stock status and over exploitation
}
\usage{
Pplot(MSEobj, YVar=c("B_BMSY", "F_FMSY"), MPs=NA, sims=NULL, 
  traj=c("all", "quant"), quants=c(0.1, 0.9), incquant=TRUE, quantcol="lightgray",
  RefYield=c("lto", "curr"), LastYr=TRUE,  maxMP=6, alpha=60, cex.axis=1.35, 
  cex.lab=1.4,  YLab=NULL, incMP=TRUE, MPcex=1.4, incLeg=TRUE, cex.leg=1.5, 
  legPos="topleft", yline=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MSEobj}{
An object of class MSE
}
  \item{YVar}{
   What to plot on the y-axis? Options are: \code{c("B_B0", "B_BMSY", "F_FMSY", "Yield")}
}
  \item{MPs}{
  Optional subset by MP
}
  \item{sims}{
  Optional subset by simulation
}
  \item{traj}{
   Plot all projections (\code{all}) or only quantiles (\code{quant})
}
  \item{quants}{
  Numeric vector of length 2 specifying the quantiles (e.g., 10th and 90th. Median is always included)
}
  \item{incquant}{
  Logical. Include the quantiles or only plot median?
}
  \item{quantcol}{
  Colour of the quantile polygon
}
  \item{RefYield}{
   Should yield be relative to long-term optimum (\code{lto}) or last historical year (\code{curr})
}
  \item{LastYr}{
   Logical. Include the last historical year in the yield projections?
}
  \item{maxMP}{
  Maximum number of MPs to plot
}
  \item{alpha}{
  Alpha for transparency of lines
}
  \item{cex.axis}{
  Size of axis text
}
  \item{cex.lab}{
  Size of axis label
}
  \item{YLab}{
  Optional label for y-axis
}
  \item{incMP}{
  Logical. Include name of MP?
}
  \item{MPcex}{
  Size of MP label
}
  \item{incLeg}{
  Logical. Include a legend?
}
  \item{cex.leg}{
  Size of legend text
}
  \item{legPos}{
  Legend position
}
  \item{yline}{
   Optional horizontal line
}
  \item{...}{
   Additional arguments to be passed to plotting functions
}
}
\author{
T. Carruthers & A.Hordyk
}
