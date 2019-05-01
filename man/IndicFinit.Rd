\name{IndicFinit}
\alias{IndicFinit}
\title{
Create initial F matrices using specific intergroup actions for network inference.
}
\description{
This is an helper function to create initial values F matrices for networks.
}
\usage{
IndicFinit(sqF, ngrp, Indic, low.trig = TRUE)
}
\arguments{
  \item{sqF}{Size of an F cell}
  \item{ngrp}{Number of groups}
  \item{Indic}{Matrix to specify where there is an interaction from one group to another}
  \item{low.trig}{Fill the lower trigonal matrices with ones}
}
\value{
An array of size (sqF, sqF, ngrp).
}
\examples{
IndicFinit(3, 2, matrix(1,2,2)-diag(2))
}
\keyword{models}
