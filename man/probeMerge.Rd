\name{probeMerge}
\alias{probeMerge}
\alias{probeMerge,micro_array-method}
\title{
Function to merge probesets
}
\description{
Used to collapse probesets using the collapseRows function of the WGCNA package
}
\usage{
probeMerge(x, ...)
}
\arguments{
  \item{x}{Microarray}
  \item{\dots}{Additionnal parameters to the collapseRows function of the WGCNA package}
}
\value{
Formal class 'micro_array' [package "Patterns"] with 7 slots
}
\examples{
if(require(CascadeData)){
data(micro_S)
D<-as.micro_array(micro_S[1:2000,],1:4,6)
D@gene_ID<-jetset::scores.hgu133plus2[D@name,"EntrezID"]
PM <- probeMerge(D)
}
}
\keyword{manip}
