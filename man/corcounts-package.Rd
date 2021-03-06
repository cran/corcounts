\name{corcounts-package}
\alias{corcounts-package}
\alias{corcounts}
\docType{package}
\title{
Generate correlated count random variables
}
\description{
Sample high-dimensional correlated count random variables
with approximate prespecified Pearson correlation and exact margins.
}
\details{
\tabular{ll}{
Package: \tab corcounts\cr
Type: \tab Package\cr
Version: \tab 1.4\cr
Date: \tab 2009-11-13\cr
License: \tab GPL (>= 3)\cr
LazyLoad: \tab yes\cr
}
Specify the marginal distributions and parameters and the desired correlation matrix and run 'rcounts()'.

In order to allow for regression, i.e. for individual parameters for each cluster, run 'rcounts.reg()'.
}
\author{
Maintainer: Vinzenz Erhardt <erhardt@ma.tum.de>
}
\references{
Erhardt, V., Czado, C. (2009)
A method for approximately sampling high-dimensional count variables with prespecified Pearson correlation.
Submitted for publication.

Czado, C., Erhardt, V., Min, A., Wagner, S. (2007) 
Zero-inflated generalized Poisson models with regression effects on the mean, 
dispersion and zero-inflation level applied to patent outsourcing rates.
Statistical Modelling 7 (2), 125-153.
}
\keyword{ package }
