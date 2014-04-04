\name{geneNeighborhood}
\alias{geneNeighborhood}
\alias{geneNeighborhood-methods}
\alias{geneNeighborhood,network-method}
\title{Find the neighborhood of a set of nodes.}

\description{
Find the neighborhood of a set of nodes.
}
\usage{
	geneNeighborhood(net,targets,...)
}

\arguments{
 \item{net}{a network object}
 \item{targets}{a vector containing the set of nodes}
 \item{...}{Optional arguments. See plot options. }
}

\value{
  The neighborhood of the targeted genes. 
}

\references{
Vallat, L., Kemper, C. A., Jung, N., Maumy-Bertrand, M., Bertrand, F., Meyer, N., ... & Bahram, S. (2013). Reverse-engineering the genetic circuitry of a cancer cell with predicted intervention in chronic lymphocytic leukemia. \emph{Proceedings of the National Academy of Sciences}, 110(2), 459-464.
}

\author{
Nicolas Jung, Bertrand Frederic, Myriam Maumy-Bertrand.
}



\examples{
	#See vignette
}

\keyword{methods}
