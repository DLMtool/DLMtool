\name{Sense}
\alias{Sense}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sensitivity analysis
}
\description{
A function that determines the inputs for a given data-limited method of class DLM_output and then analyses the sensitivity of TAC estimates to marginal differences in each input. The range used for sensitivity is based on the user-specified CV for that input (e.g. CV_Mort, Mort)
}
\usage{
Sense(DLM_data, MP, nsense = 6, reps = 100, perc = c(0.05, 0.5, 0.95), ploty = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DLM_data}{
A data-limited methods data object (class DLM_data)
}
  \item{MP}{
A character string representing an MP applied in calculating the TAC recommendations in the DLM object
}
  \item{nsense}{
The number of points over which to calculate the TAC (resolution)
}
  \item{reps}{
The number of samples of the quota taken for the calculation of the TAC
}
  \item{perc}{
The percentile of the sample TAC
}
  \item{ploty}{
A logical switch, (T/F, should a plot be drawn?)
}
}
\author{
T. Carruthers
}
