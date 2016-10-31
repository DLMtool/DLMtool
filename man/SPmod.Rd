\name{SPmod}
\alias{SPmod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Surplus production based catch-limit modifier}
\description{
An MP that makes incremental adjustments to TAC recommendations based on the apparent trend in surplus production. Based on the theory of Mark Maunder (IATTC)
}
\usage{
SPmod(x, DLM_data, reps = 100, alp = c(0.8, 1.2), bet = c(0.8, 1.2))
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
\item{alp}{
Condition for modifying the TAC (bounds on change in abundance)
}
  \item{bet}{
Limits for how much the TAC can change among years
}
}
\details{
Note that this isn't exactly what Mark has previously suggested and is stochastic in this implementation. 
}
\value{
A numeric vector of TAC recommendations
}
\references{
http://www.iattc.org/Meetings/Meetings2014/MAYSAC/PDFs/SAC-05-10b-Management-Strategy-Evaluation.pdf
}
\author{
T. Carruthers
}
