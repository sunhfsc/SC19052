% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rwC}
\alias{rwC}
\title{A random walk Metropolis sampler using Rcpp}
\usage{
rwC(sigma, x0, N)
}
\arguments{
\item{sigma}{The standard deviation in the normal distribution (double)}

\item{x0}{The initial point (double)}

\item{N}{The length of the chain (int)}
}
\value{
The random numbers and the number of rejections (list)
}
\description{
A random walk Metropolis sampler for generating the standard Laplace distribution
}
\examples{
\dontrun{
N<-1000
sigma<-0.05
x0<-25
rwC(sigma,x0,N)
}
}
