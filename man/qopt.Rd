\name{qopt}
\alias{qopt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Internal optimization function that find the catchability (q where F=qE) value required to get to user-specified stock depletion (current biomass / unfished biomass)}
\description{
The user specifies the level of stock depleiton. This funciton takes the derived effort trajectories and finds the catchabiltiy to get the stock there.   
}
\usage{
qopt(lnq,depc,Fc,Perrc,Mc,hc,Mac,
     Wac,R0c,Vc,nyears,maxage,movc,
     Spat_targc,SRrelc,aRc,bRc,opt=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lnq}{
}
  \item{depc}{
}
  \item{Fc}{
}
  \item{Perrc}{
}
  \item{Mc}{
}
  \item{hc}{
}
  \item{Mac}{
}
  \item{Wac}{
}
  \item{R0c}{
}
  \item{Vc}{
}
  \item{nyears}{
}
  \item{maxage}{
}
  \item{movc}{
}
  \item{Spat_targc}{
}
  \item{SRrelc}{
}
  \item{aRc}{
}
  \item{bRc}{
}
  \item{opt}{
}
}
\details{
Paired with qopt. 
}
\author{
T. Carruthers
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
