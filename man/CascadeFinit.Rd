\name{CascadeFinit}
\alias{CascadeFinit}
\title{
Create initial F matrices for cascade networks inference.
}
\description{
This is an helper function to create initial values F matrices for cascade networks.
}
\usage{
CascadeFinit(sqF, ngrp, low.trig = TRUE)
}
\arguments{
  \item{sqF}{Size of an F cell}
  \item{ngrp}{Number of groups}
  \item{low.trig}{Fill the lower trigonal matrices with ones}
}
\value{
An array of sizes c(sqF, sqF, ngrp).
}
\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}

\examples{
CascadeFinit(3,2)
CascadeFinit(4,3)
plotF(CascadeFinit(4,3),choice = "F")
CascadeFinit(3,2,low.trig=FALSE)
CascadeFinit(4,3,low.trig=FALSE)
plotF(CascadeFinit(4,3,low.trig=FALSE),choice = "F")
}
\keyword{models}