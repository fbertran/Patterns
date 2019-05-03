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
\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}

\examples{
CascadeFshape(3,2)
plotF(CascadeFshape(3,2),choice = "Fshape")
CascadeFshape(4,3)
plotF(CascadeFshape(4,3),choice = "Fshape")
}

\keyword{models}
