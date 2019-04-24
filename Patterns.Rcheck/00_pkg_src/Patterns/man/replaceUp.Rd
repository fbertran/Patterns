\name{replaceUp}
\alias{replaceUp}
\title{
Replace matrix values triangular upper part and by band for the lower part.
}
\description{
F matrices utility function.
}
\usage{
replaceUp(a, b, k)
}
\arguments{
  \item{a}{The matrix to be replaced}
  \item{b}{The matrix with the replacement values}
  \item{k}{The extend of the replacement: 0 (upper part only), 1 (upper part and first extra diagonal), in general an entry is replaced if (row(a) - col(a)) <= k}
}
\value{
A matrix (same size as a) 
}
\examples{
a=matrix(1:9,3,3)
b=matrix(1,3,3)
replaceUp(a,b,0)
replaceUp(a,b,1)
replaceUp(a,b,2)
}
\keyword{manip}
