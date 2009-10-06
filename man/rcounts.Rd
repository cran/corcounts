\name{rcounts}
\alias{rcounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate correlated count random variables }
\description{
'rcounts' is used to sample high-dimensional correlated count random variables
with approximate prespecified Pearson correlation and exact margins.
}
\usage{
rcounts(N, margins, mu, phi=rep(NA,length(margins)), 
        omega=rep(NA,length(margins)), psi=rep(NA,length(margins)), 
        corstr, corpar, conv=0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{ number of observations to be generated per margin (should be at least 500).}
  \item{margins}{ Vector of margin tokens. Its length T is the dimension. See details. }
  \item{mu}{ Vector of length T of means for the Poisson, GP, ZIP, ZIGP and NB margins. }
  \item{phi}{ Vector of length T of dispersion parameters for the GP, and ZIGP margins. For Poisson, ZIP and NB margins, an 'NA' can be provided. }
  \item{omega}{ Vector of length T of zero-inflation parameters for the ZIP and ZIGP margins. For Poisson, GP and NB margins, an 'NA' can be provided. }
  \item{psi}{ Vector of length T of size parameters for the NB margins. For Poisson, GP, ZIP and ZIGP margins, an 'NA' can be provided. }
  \item{corstr}{ Correlation structure. Can be 'ex' for exchangeable, 'AR1' for AR(1) and 'unstr' for unstructured. }
  \item{corpar}{ Correlation parameter. Scalar correlation for 'ex' and 'AR1' and matrix of dimension TxT for 'unstr'. }
  \item{conv}{ Convergence criterion }
}
\details{
The entries in 'margins' can be specified as 'Poi' for Poisson, 'GP' for generalized Poisson, 'ZIP' for zero-inflated Poisson, 'ZIGP' for zero-inflated generalized Poisson and 'NB' for negative-binomial.
}
\value{
The function will return a matrix of counts of dimension N x T.
}
\author{ Vinzenz Erhardt }
\examples{
N <- 5000

# high precision in dimension 2
margins <- c("ZIGP","GP")
mu <-    c(10,  15)
phi <-   c(1.5, 3.5)
omega <- c(.25, NA)
corstr <- "ex"
corpar <- .5
Y <- rcounts(N=N, margins=margins, mu=mu, phi=phi, omega=omega, corstr=corstr, 
             corpar=corpar, conv=0.0001)
cor(Y)


# five-dimensional examples
margins <- c("ZIGP","GP","Poi","NB","ZIP")
mu <-    c(10,  25,  12, 20, 28)
phi <-   c(1.5, 2,   NA, NA, NA)
omega <- c(.25, NA,  NA, NA, .2)
psi <-   c(NA,  NA,  NA, 7,  NA)

# Exchangeable structure with correlation of 0.5
corstr <- "ex"
corpar <- .5
Y <- rcounts(N=N, margins=margins, mu=mu, phi=phi, omega=omega, psi=psi, 
             corstr=corstr, corpar=corpar)
cor(Y)

# AR(1) structure with correlation of corr(Y(t1), Y(t2)) = .8 ^ |t1-t2|
corstr <- "AR1"
corpar <- .8
Y <- rcounts(N=N, margins=margins, mu=mu, phi=phi, omega=omega, psi=psi, 
             corstr=corstr, corpar=corpar)
cor(Y)

# Unstructured correlation. Create random symmetric positive definite 
# matrix using function 'unstructured'
corstr <- "unstr"
corpar <- unstructured(5)
corpar
Y <- rcounts(N=N, margins=margins, mu=mu, phi=phi, omega=omega, psi=psi, 
             corstr=corstr, corpar=corpar)
cor(Y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ distribution }
