\name{NOAA_plot}
\alias{NOAA_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
National Oceanographic and Atmospheric Administration default plot 1
}
\description{
A preliminary plot for returning trade-offs plots and performance table for total yield, variability in yield, probability of overfishing and likelihood of biomass dropping below 50 per cent BMSY
}
\usage{
NOAA_plot(MSEobj,nam=NA,type=NA,panel=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MSEobj}{
An object of class MSE
}
\item{nam}{
Title of plot
}
\item{type}{
Plots full range of data if NA. Plots a subset that meet thresholds if not NA. 
}
\item{panel}{
Should a two panel plot be made or should plots be made in sequence.
}
}
\value{
A table of performance metrics.
}
\author{
T. Carruthers
}

