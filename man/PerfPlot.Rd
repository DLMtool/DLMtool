\name{PerfPlot}
\alias{PerfPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Joint Satisficing Trade-off Plot  
}
\description{
Calculates and plots the probability of each method meeting the individual performance metrics, 
as well as joint satisficing; that is, the probabilitity that a method meets the minimum requirements 
for the  performance metrics.
}
\usage{
PerfPlot(MSEobj, PMs = c("B_BMSY", "F_FMSY", "AAVY"), PLim = 50, lastYrs = 10, 
  BmsyRef = 0.5, B0Ref = 0.2, FRef = 1.25, maxVar = 15, 
  ShowCols = TRUE, ShowLabs = TRUE, AvailMPs = NULL, 
  bplot = TRUE, maxmp = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MSEobj}{
Object of class MSE, output of the runMSE function
}
  \item{PMs}{
Character vector of performance metrics.  Choices are \code{"B_BMSY", "B_B0", "F_FMSY", "AAVY"}, and defaults are \code{"B_BMSY", "F_FMSY", "AAVY"}.
}
  \item{PLim}{
Minimum probability threshold for meeting the performance criteria as a percentage. Default is 50\%.
}
  \item{lastYrs}{
The final number of years of the projection period to summarize results over. 
}
  \item{BmsyRef}{
Reference level for BMSY, in proportion, i.e., 0.5 = 0.5BMSY
}
  \item{B0Ref}{
Reference level for B0, in proportion, i.e., 0.2 = 0.2B0
}
  \item{FRef}{
Reference level for F relative to FMSY, in proportion, i.e., 1.25 = 1.25FMSY
}
  \item{maxVar}{
Reference for average annual variability in yield in percent
}
  \item{ShowCols}{
Logical to specify if MP labels are shown
}
  \item{ShowLabs}{
Logical to specify if background colors are shown
}
  \item{AvailMPs}{
vector of MPs that *could* be applied to the fishery, i.e., sufficient data exists. These a plotted with different symbol
}
  \item{bplot}{
Logical to specify if the barplots should be created
}
  \item{maxmp}{
Maximum number of methods to include in each barplot
}
}
\details{
Returns a list containing the probability of each method meeting the individual performance metrics, and the joint satisficing of all the performance metrics with average yield in final \code{lastYrs} of the projection period. 
}
\author{
A. Hordyk
}
