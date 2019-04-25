\name{CascadeFshape}
\alias{CascadeFshape}
\title{
Create F matrices shaped for cascade networks inference.
}
\description{
This is an helper function to create F matrices with special shape used for cascade networks.
}
\usage{
CascadeFshape(sqF, ngrp)
}
\arguments{
  \item{sqF}{Size of an F cell}
  \item{ngrp}{Number of groups}
}
\value{
An array of sizes c(sqF, sqF, ngrp).
}

\examples{
CascadeFshape(3,2)
CascadeFshape(4,3)
}

\keyword{models}
