\name{OptFun}
\alias{OptFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Objective function for LBSPR methods
}
\description{
Used internally to define objective function for LBSPR methods.
}
\usage{
OptFun(tryFleetPars, LenDat, StockPars, SizeBins = NULL, mod = c("GTG", "LBSPR"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tryFleetPars}{
Vector of relative fishing mortality (F/M) and selectivity-at-length parameters.
}
  \item{LenDat}{
Binned length data
}
  \item{StockPars}{
Life history parameters of stock.
}
  \item{SizeBins}{
Information on the length bins.
}
  \item{mod}{
Optional for alternative models - only "LBSPR" currently used.
}
}

\author{
A. Hordyk 
}
