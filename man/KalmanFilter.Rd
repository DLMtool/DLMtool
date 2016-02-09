\name{KalmanFilter}
\alias{KalmanFilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kalman filter and Rauch-Tung-Striebel smoother
}
\description{
Kalman filter to predict new points and smoother for time-series. 
}
\usage{
KalmanFilter(RawEsts, R = 1, Q = 0.1, Int = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RawEsts}{
Vector of numeric values to be filtered and smoothed.
}
  \item{R}{
Variance of sampling noise
}
  \item{Q}{
Variance of random walk increments
}
  \item{Int}{
Covariance of initial uncertainty
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
A. Hordyk
}

