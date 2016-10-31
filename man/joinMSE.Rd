\name{joinMSE}
\alias{joinMSE}
\title{Join multiple MSE objects together}
\description{
Joins two or more MSE objects together. MSE objects must have identical number of historical years, and projection years. 
}
\usage{
joinMSE(MSEobjs = NULL)
}
\arguments{
  \item{MSEobjs}{
A list of MSE objects. Must all have identical operating model and MPs. MPs which don't appear in all MSE objects will be dropped.
}
}
\value{
An object of class \code{MSE} 
}
\author{
A. Hordyk
}

