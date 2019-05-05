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

\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}



\examples{
data(Selection)
data(infos)
#Find probesets for EGR1
pbst_EGR1 = infos[infos$hgnc_symbol=="EGR1", "affy_hg_u133_plus_2"]

gene_IDs = infos[match(Selection@name, infos$affy_hg_u133_plus_), "hgnc_symbol"]

data(network)
#A nv value can chosen using the cutoff function
nv=.11 
EGR1<-which(is.element(Selection@name,pbst_EGR1))
P<-position(network,nv=nv)

geneNeighborhood(network,targets=EGR1,nv=nv,ini=P,
label_v=gene_IDs)
}

\keyword{methods}
