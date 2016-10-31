\name{TradePlot}
\alias{TradePlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generic Trade-off Plot}
\description{
Creates a trade-off plot (up to four panels) of built-in performance metrics.  
}
\usage{
TradePlot(MSEobj, XAxis=c("Overfishing", "Biomass:BMSY"), 
	YAxis=c("Long-term Yield", "AnnualVar"), XThresh=c(30, 80), YThresh=c(0,50),
	maxVar=15, BmsyRef=0.5, B0Ref=0.2, AvailMPs=NULL, ShowLabs=FALSE, 
	ShowCols=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MSEobj}{ Object of class MSE, output of the runMSE function}
\item{XAxis}{Character string describing the performance metrics for the x-axis (or x-axes if vector; max 4). Must be chosen for list of existing PMs and same length as YAxis. See \code{PMs}}
\item{YAxis}{Character string describing the performance metrics for the y-axis (or y-axes if vector; max 4). Must be chosen for list of existing PMs and same length as XAxis. See \code{PMs}}
\item{XThresh}{Minimum threshold values in percent (i.e., 50 = 50\%) for the x-axes (must be same length as XAxis)}
\item{YThresh}{Minimum threshold values in percent (i.e., 50 = 50\%) for the y-axes (must be same length as YAxis)}
\item{maxVar}{Reference for average annual variability in yield in percent}
\item{BmsyRef}{Reference level of BMSY, in proportion, i.e., 0.5 = 0.5BMSY}
\item{B0Ref}{Reference level of B0, in proportion, i.e., 0.2 = 0.2B0}
\item{AvailMPs}{vector of MPs that *could* be applied to the fishery, i.e., sufficient data exists. These a plotted with different symbol}
\item{ShowLabs}{Logical to specify if MP labels are shown}
\item{ShowCols}{Logical to specify if background colors are shown}
}
\details{
Returns a list containing the names of performance metrics that meet the minimum performance metrics for each trade-off, and ranks the MPs by increasing distance from the top-right corner. 
}

\author{
A. Hordyk
}

