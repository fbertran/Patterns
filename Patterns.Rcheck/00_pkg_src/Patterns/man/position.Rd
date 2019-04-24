\name{position}
\alias{position}
\title{
Retrieve network position for consistent plotting.
}
\description{
Utility function to plot networks.
}
\usage{
position(net, ...)
}
\arguments{
  \item{net}{Network}
  \item{\dots}{Additionnal parameters
}
}
\value{
Matrix with as many rows as the number of edges of network and three columns (name, xcoord, ycoord).
}
\examples{
data(network)
position(network)
}
\keyword{dplot}
