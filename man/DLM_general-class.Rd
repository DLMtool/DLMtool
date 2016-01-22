\name{DLM_general-class}
\Rdversion{1.1}
\docType{class}
\alias{DLM_general-class}

\title{Class \code{"DLM_general"}}
\description{
An object for storing general toolkit data. The data are stored in the right format in the slot data. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("DLM_general", stock)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{Name}:}{The name of the data}
    \item{\code{data}:}{The data correctly formated}
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "DLM_general")}: ... }
    \item{plot}{\code{signature(x = "DLM_general")}: ... }
    \item{summary}{\code{signature(object = "DLM_general")}: ... }
	 }
}

\author{
T. Carruthers
}
\examples{
newdata<-new('DLM_general')
}
\keyword{classes}
