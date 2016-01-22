\name{YPR_CC}
\alias{YPR_CC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Yield Per Recruit analysis to get FMSY proxy F01 paired to a naive catch curve estimate of recent Z
}
\description{
A simple yield per recruit approximation to FMSY (F01) which is the position of the ascending YPR curve for which dYPR/dF = 0.1(dYPR/d0)
A naive catch-curve analysis is used to determine recent Z which given M (Mort) gives F and thus abundance = Ct/(1-exp(-F))
}
\usage{
YPR_CC(x, DLM_data, reps = 100, Fmin=0.005)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
A position in a data-limited methods data object
}
\item{DLM_data}{
A data-limited methods data object (class DLM)
}
  \item{reps}{
The number of samples of the TAC
}
  \item{Fmin}{
The minimum fishing mortality rate inferred from the catch-curve analysis
}
}
\author{
Meaghan Bryan and T. Carruthers
}
