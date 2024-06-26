\name{clustInference}
\alias{clustInference}
\alias{clustInference-methods}
\alias{clustInference,micro_array,numeric-method}
\title{
A function to explore a dataset and cluster its rows.
}
\description{
Based on soft clustering performed by the Mfuzz package.
}
\usage{
clustInference(microarray, vote.index, ...)
}
\arguments{
  \item{microarray}{A microarray to cluster}
  \item{vote.index}{Option for cluster attribution}
  \item{\dots}{
    Optional arguments: 	  	
    \describe{ 
      \item{new.window}{Boolean. New X11 window for plots. Defaults to FALSE.}
    }
  }
}
\value{A list of two elements:
  \item{res.matrix}{A data.frame of nrows(microarray) observations of 3 variables (name, cluster, maj.vote.index).}
  \item{prop.matrix}{Additionnal info.}
}

\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}

\examples{
library(Patterns)
if(require(CascadeData)){
data(micro_S, package="CascadeData")
D<-Patterns::as.micro_array(micro_S[1:20,],1:4,6)
b<-Patterns::clustInference(D,0.5)
b
}
}
\keyword{cluster}
