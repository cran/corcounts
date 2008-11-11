\name{unstructured}
\alias{unstructured}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate random correlation matrix }
\description{ 'unstructured' generates a random correlation matrix of dimension T with random entries. To ensure positive definiteness, a matrix of partial correlations with random entries uniform on [-0.9, 0.9] will be generated an the corresponding correlation matrix be calculated from it using a bijective recursive relation between them. }
\usage{
unstructured(dimension)
}
\arguments{
  \item{dimension}{ Dimension T of the correlation matrix. }
}
\value{
A correlation matrix of dimension T x T.}
\examples{
unstructured(10)
}
\keyword{ multivariate }