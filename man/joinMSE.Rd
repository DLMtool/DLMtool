\name{joinMSE}
\alias{joinMSE}
\title{
Join multiple MSE objects together
}
\description{
Joins two or more MSE objects together. MSE objects must have identical number of historical years, projection years, and management procedures.
}
\usage{
joinMSE(MSEobj1 = NULL, MSEobj2 = NULL, MSEobjs = NULL)
}
\arguments{
  \item{MSEobj1}{
The first object of class \code{MSE} 
}
  \item{MSEobj2}{
A second object of class \code{MSE} that has identical operating model and MPs as \code{MSEobj1}
}
  \item{MSEobjs}{
A character vector with the names of the MSE objects. Must all have identical operating model and MPs. \code{MSEobj1} and \code{MSEobj2} are ignored if \code{MSEobjs} is provided.
}
}
\value{
An object of class \code{MSE} 
}
\author{
A. Hordyk
}

