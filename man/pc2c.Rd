\name{pc2c}
\alias{pc2c}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate a corresponding correlation matrix to partial correlations }
\description{
'c2pc' is used to calculate the corresponding correlation matrix of dimension T times T out of partial correlations.
}
\usage{
pc2c(Theta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Theta}{ A T times T matrix with partial correlations. See details. }
}
\details{
The partial correlations in Theta have to be specified as

Theta =

12........13........14..........15...........16

..........23|1......24|1........25|1.........26|1
          
....................34|12.......35|12........36|12     ...
                    
................................45|123.......46|123
                                
.............................................56|1234
                                             

                        ...

and may be NA elsewhere. Theta has to be of dimension T times T.

This routine only calculates partial correlations conditional on 1, 12, 123, 1234, etc.. Partial correlations conditional on other margins can be obtained by a permutation of margins.
}
\value{
A symmetric positive definite correlation matrix of dimension T times T.
}
\author{ Vinzenz Erhardt }
\examples{
# create random uniform(0,1) partial correlations in dimension 8
dimension <- 8
Theta <- matrix(NA,dimension,dimension)
for (i in 2:dimension) {
  for (j in 1:(i-1)) {
    Theta[j,i] <- runif(1,-1,1)
  }
}
Theta

# calculate corresponding correlation matrix
C <- pc2c(Theta)
C

# transform back to partial correlations
c2pc(C)

# equivalence with original Theta
Theta - c2pc(C)
}
\keyword{ multivariate }