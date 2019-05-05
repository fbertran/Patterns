\name{IndicFshape}
\alias{IndicFshape}
\title{
Create F matrices using specific intergroup actions for network inference.
}
\description{
This is an helper function to create values F matrices using specific intergroup actions for network inference.
}
\usage{
IndicFshape(sqF, ngrp, Indic)
}
\arguments{
  \item{sqF}{Size of an F cell}
  \item{ngrp}{Number of groups}
  \item{Indic}{Matrix to specify where there is an interaction from one group to another}
}
\value{
An array of size (sqF, sqF, ngrp).
}
\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}

\examples{
IndicFshape(3, 2, matrix(1,2,2)-diag(2))
}
\keyword{models}
