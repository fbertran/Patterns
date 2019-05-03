\name{replaceDown}
\alias{replaceDown}
\title{
Replace matrix values triangular lower part and by band for the upper part.
}
\description{
F matrices utility function.
}
\usage{
replaceDown(a, b, k)
}
\arguments{
  \item{a}{The matrix to be replaced}
  \item{b}{The matrix with the replacement values}
  \item{k}{The extend of the replacement: 0 (lower part and diagonal only), 1 (lower part and first extra diagonal), in general an entry is replaced if -(row(a) - col(a)) <= k}
}
\value{
A matrix (same size as a) 
}
\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}

\examples{
a=matrix(1:9,3,3)
b=matrix(1,3,3)
replaceDown(a,b,0)
replaceDown(a,b,1)
replaceDown(a,b,2)
}
\keyword{manip}
