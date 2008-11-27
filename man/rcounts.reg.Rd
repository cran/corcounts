\name{rcounts.reg}
\alias{rcounts.reg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate correlated count random variables with individual parameters for each cluster }
\description{
'rcounts.reg' is used to sample high-dimensional correlated count random variables
with approximate prespecified Pearson correlation and exact margins.
}
\usage{
rcounts.reg(N, margins, mu, phi=matrix(NA,N,length(margins)), omega=matrix(NA,N,length(margins)), psi=matrix(NA,N,length(margins)), corstr, corpar, conv=0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{ number of observations to be generated per margin (should be at least 500).}
  \item{margins}{ Vector of margin tokens. Its length T is the dimension. See details. }
  \item{mu}{ Matrix of dimension N x T of means for the Poisson, GP, ZIP, ZIGP and NB margins. }
  \item{phi}{ Matrix of dimension N x T of dispersion parameters for the GP, and ZIGP margins. For Poisson, ZIP and NB margins, an 'NA' can be provided. }
  \item{omega}{ Matrix of dimension N x T of zero-inflation parameters for the ZIP and ZIGP margins. For Poisson, GP and NB margins, an 'NA' can be provided. }
  \item{psi}{ Matrix of dimension N x T of size parameters for the NB margins. For Poisson, GP, ZIP and ZIGP margins, an 'NA' can be provided. }
  \item{corstr}{ Correlation structure. Can be 'ex' for exchangeable, 'AR1' for AR(1) and 'unstr' for unstructured. }
  \item{corpar}{ Correlation parameter. Scalar correlation for 'ex' and 'AR1' and matrix of dimension TxT for 'unstr'. }
  \item{conv}{ Convergence criterion }
}
\details{
The entries in 'margins' can be specified as 'Poi' for Poisson, 'GP' for generalized Poisson, 'ZIP' for zero-inflated Poisson, 'ZIGP' for zero-inflated generalized Poisson and 'NB' for negative-binomial.

NOTE: there is a tradeoff between too small N (decreasing accuracy of the resulting correlation) and too high N (dramatically increasing computation time).
}
\value{
The function will return a matrix of counts of dimension N x T.
}
\author{ Vinzenz Erhardt }
\examples{
N <- 500

# bivariate example
margins <- c("ZIGP","GP")
mu <-    matrix(runif(N*2,10,20),N,2)
phi <-   matrix(runif(N*2,1,3),N,2)
omega <- matrix(c(runif(N,0,.3),rep(NA,N)),N,2)
corstr <- "ex"
corpar <- .5
Y <- rcounts.reg(N=N, margins=margins, mu=mu, phi=phi, omega=omega, corstr=corstr, corpar=corpar)
cor(Y)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ distribution }
