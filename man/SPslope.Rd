\name{SPslope}
\alias{SPslope}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Slope in surplus production MP}
\description{
A management procedure that makes incremental adjustments to TAC recommendations based on the apparent trend in recent surplus production. Based on the theory of Mark Maunder (IATTC)
}
\usage{
SPslope(x, DLM_data, reps = 100, yrsmth = 4, alp = c(0.9, 1.1), bet = c(1.5, 0.9))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A position in data-limited methods data object
}
  \item{DLM_data}{
A data-limited methods data object
}
  \item{reps}{
The number of quota samples
}
  \item{yrsmth}{
Years over which to smooth recent estimates of surplus production  
}
\item{alp}{
Condition for modifying the DLM_data (bounds on change in abundance)
}
  \item{bet}{
Limits for how much the DLM_data can change among years
}
}
\details{
Note that this isn't exactly what Mark has previously suggested and is stochastic in this implementation. 
}
\value{
A numeric vector of DLM_data recommendations
}
\references{
http://www.iattc.org/Meetings/Meetings2014/MAYSAC/PDFs/SAC-05-10b-Management-Strategy-Evaluation.pdf
}
\author{
T. Carruthers
}
