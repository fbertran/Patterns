\name{evolution}
\alias{evolution}
\alias{evolution-methods}
\alias{evolution,network-method}

\title{See the evolution of the network with change of cutoff}
\description{
See the evolution of the network with change of cutoff
}
\usage{
	evolution(net,list_nv,\dots)
}

\arguments{
 \item{net}{a network object}
 \item{list_nv}{a vector of cutoff at which the network should be shown}
 \item{\dots}{Optionnal arguments:
 	\describe{ \item{gr}{a vector giving the group of each gene}
  \item{color.vertex}{a vector giving the color of each node}
  \item{fix}{logical, should the position of the node in the network be calculated once at the beginning ? Defaut to TRUE.}
  \item{taille}{vector giving the size of the plot. Default to c(2000,1000)}
}
}}

\value{ A HTML page with the evolution of the network.}


\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}



\examples{
	data(network)
	sequence<-seq(0,0.2,length.out=20)
	evolution(network,sequence)
}

\keyword{methods}

