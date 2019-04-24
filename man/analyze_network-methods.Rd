\name{analyze_network}
\alias{analyze_network}
\alias{analyze_network-methods}
\alias{analyze_network,network-method}

\title{Analysing the network}
\description{
Calculates some indicators for each node in the network.
}
\usage{
	analyze_network(Omega,nv,...)
}

\arguments{
 \item{Omega}{a network object}
 \item{nv}{the level of cutoff at which the analysis should be done}	
 \item{...}{Optional arguments}
}

\value{
  A matrix containing, for each node, its betweenness,its degree, its output, its closeness.
}

\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}


\examples{
	data(network)
#	analyze_network(network,nv=0)
}

\keyword{methods}




