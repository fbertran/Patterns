\name{unsupervised_clustering_auto_m_c}
\alias{unsupervised_clustering_auto_m_c}
\alias{unsupervised_clustering_auto_m_c,micro_array-method}
\title{
Cluster a micro_array object: determine optimal fuzzification parameter and number of clusters.
}
\description{
Based on soft clustering performed by the Mfuzz package.
}
\usage{
unsupervised_clustering_auto_m_c(M1, ...)
}
\arguments{
  \item{M1}{Object of micro_array class.}
  \item{\dots}{Additional parameters.}
}
\value{
\item{m}{Estimate of the optimal fuzzification parameter.}
\item{c}{Estimate of the optimal number of clusters.}
\item{csearch}{More result from the cselection function of the Mfuzz package}
}
\examples{
if(require(CascadeData)){
data(micro_S, package="CascadeData")
M<-as.micro_array(micro_S[1:100,],1:4,6)
mc<-unsupervised_clustering_auto_m_c(M)
}
}
\keyword{cluster}
