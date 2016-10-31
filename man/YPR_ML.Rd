\name{YPR_ML}
\alias{YPR_ML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Yield Per Recruit analysis to get FMSY proxy F01 paired with a mean-length estimate of current stock size}
\description{
A simple yield per recruit approximation to FMSY (F01) which is the position of the ascending YPR curve for which dYPR/dF = 0.1(dYPR/d0)
A mean-length estimate of recent Z is used to infer current abundance
}
\usage{
YPR_ML(x, DLM_data, reps = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A position in a data-limited methods data object
}
  \item{DLM_data}{
A data-limited methods data object
}
  \item{reps}{
The number of samples of the TAC
}
}
\author{
Meaghan Bryan and T. Carruthers
}
\note{
The mean length extension was programmed by Gary Nelson as part of his excellent R package 'fishmethods'
}
