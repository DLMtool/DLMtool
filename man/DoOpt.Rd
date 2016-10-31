\name{DoOpt}
\alias{DoOpt}

\title{Optimization routine for LBSPR methods}
\description{
Internal optimization routine. 
}
\usage{
DoOpt(StockPars, LenDat, SizeBins = NULL, mod = c("GTG", "LBSPR"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{StockPars}{
Life history parameters of stock.
}
  \item{LenDat}{
Binned length data
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
