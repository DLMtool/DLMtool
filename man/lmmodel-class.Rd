\name{lmmodel-class}
\Rdversion{1.1}
\docType{class}
\alias{lmmodel-class}

\title{Class \code{"lmmodel"}}
\description{
An object for storing fitted linear model objects in this case the relationship between M, age-at-maturity and the von B. K parameter.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("lmmodel", stock)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{Name}:}{The Name of the list of linear models}
    \item{\code{models}:}{A list of fitted linear models}
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "lmmodel")}: ... }
	 }
}
\author{
T. Carruthers
}
\examples{
newdata<-new('lmmodel',"Name",new('list'))
}
\keyword{classes}
