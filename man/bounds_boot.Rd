% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multipleversions_code_110316.R
\name{bounds_boot}
\alias{bounds_boot}
\title{Bootstrapping upper and lower bounds}
\usage{
bounds_boot(Y0, Y1, D0, D1, n_boot)
}
\arguments{
\item{Y0}{observed outcomes in the Z = 0 arm}

\item{Y1}{observed outcomes in the Z = 1 arm}

\item{D0}{observed treatment types in the Z = 0 arm}

\item{D1}{observed treatment types in the Z = 1 arm}

\item{n_boot}{number of bootstrap replications to perform}
}
\value{
an n_boot x 4 matrix of upper and lower bounds estimates
from the bootstrap samples
}
\description{
Performs non-parameteric bootstrap of upper and lower bounds
}

