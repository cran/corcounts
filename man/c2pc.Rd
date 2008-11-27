\name{c2pc}
\alias{c2pc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate partial correlations from a correlation matrix }
\description{
'c2pc' is used to calculate partial correlations from a correlation matrix.
}
\usage{
c2pc(Cin)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Cin}{ A symmetric positive definite correlation matrix. }
}
\details{
If you obtain values not in [-1,1], your correlation matrix is not positive definite.

This routine only calculates partial correlations conditional on 1, 12, 123, 1234, etc.. Partial correlations conditional on other margins can be obtained by a permutation of margins.
}
\value{
The partial correlations calculated will be

Theta =

12........13........14..........15...........16

..........23|1......24|1........25|1.........26|1
          
....................34|12.......35|12........36|12     ...
                    
................................45|123.......46|123
                                
.............................................56|1234
                                             

                        ...
}
\author{ Vinzenz Erhardt }
\examples{
# create 8 dimensional symmetric positive correlation matrix with random entries
Cin <- unstructured(8)
Cin

Theta <- c2pc(Cin)
Theta

# transform Theta back to obtain the correlation matrix
pc2c(Theta)

# identical with Cin
pc2c(Theta) - Cin
}
\seealso{ Package 'corpcor' calculates partial correlations conditional on ALL other margins. }
\keyword{ multivariate }