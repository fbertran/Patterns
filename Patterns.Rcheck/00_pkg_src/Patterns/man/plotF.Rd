\name{plotF}
\alias{plotF}
\title{Plot functions for the F matrices.}
\description{
The graphical output will differ according to the option used.
}
\usage{
plotF(x, choice = "Fshape", nround = 2, pixmap.color = terrain.colors(20))
}

\arguments{
  \item{x}{The F matrix.}
  \item{choice}{A string: either "F", "Fpixmap", "Fshape", or "Fshapepixmap"}
  \item{nround}{An integer. For numerical F matrices only. The number of decimal numbers to display.}
  \item{pixmap.color}{For pixmap plots.}
}

\value{Nothing.}

\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}


\examples{
#For numerical/inferred F matrices
plotF(CascadeFinit(4,4),choice="F", nround=1)
plotF(CascadeFinit(4,4),choice="Fpixmap")

#For theoritical F matrices
plotF(CascadeFshape(4,4),choice="Fshape")
plotF(CascadeFshape(4,4),choice="Fshapepixmap")
}

\keyword{dplot}

