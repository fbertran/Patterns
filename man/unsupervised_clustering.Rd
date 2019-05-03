\name{unsupervised_clustering}
\alias{unsupervised_clustering}
\alias{unsupervised_clustering,micro_array,numeric,numeric-method}
\title{
Cluster a micro_array object: performs the clustering.
}
\description{
Based on soft clustering performed by the Mfuzz package.
}
\usage{
unsupervised_clustering(M1, clust, mestim, ...)
}
\arguments{
  \item{M1}{Object of micro_array class.}
  \item{clust}{Number of clusters.}
  \item{mestim}{Fuzzification parameter.}
  \item{\dots}{Additional parameters.}
}
\value{
An object of class micro_array with the group slot updated by groups deduced from the soft clustering result.
}
\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}

\examples{
if(require(CascadeData)){
data(micro_S, package="CascadeData")
M<-as.micro_array(micro_S[1:100,],1:4,6)
mc<-unsupervised_clustering_auto_m_c(M)
MwithGrp=unsupervised_clustering(M, 8, mc$m, screen=NULL, heatmap=FALSE, new.window = FALSE)
# Other options
unsupervised_clustering(M, 8, mc$m, screen=c(3,3), heatmap=TRUE, new.window = FALSE)
# Plot the clusters
plot(MwithGrp)
}
}
\keyword{cluster}
