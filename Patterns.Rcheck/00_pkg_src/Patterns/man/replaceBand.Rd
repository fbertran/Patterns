\name{replaceBand}
\alias{replaceBand}
\title{
Replace matrix values by band.
}
\description{
F matrices utility function.
}
\usage{
replaceBand(a, b, k)
}
\arguments{
  \item{a}{The matrix to be replaced}
  \item{b}{The matrix with the replacement values}
  \item{k}{The extend of the replacement: 0 (diagonal only), 1 (diagonal and first extra diagonal), in general an entry is replaced if abs(row(a) - col(a)) <= k}
}
\value{
A matrix (same size as a) 
}
\examples{
a=matrix(1:9,3,3)
b=matrix(0,3,3)
replaceBand(a,b,0)
replaceBand(a,b,1)
replaceBand(a,b,2)
}
\keyword{manip}
