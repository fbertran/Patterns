\name{clustExploration}
\alias{clustExploration}
\alias{clustExploration-methods}
\alias{clustExploration,micro_array-method}
\title{
A function to explore a dataset and cluster its rows.
}
\description{
Based on soft clustering performed by the Mfuzz package.
}
\usage{
clustExploration(microarray,...)
}
\arguments{
  \item{microarray}{A microarray to cluster}
  \item{\dots}{
    Optional arguments: 	  	
    \describe{ 
      \item{new.window}{Boolean. New X11 window for plots. Defaults to FALSE.}
    }
  }
}
\value{
A data.frame of nrows(microarray) observations of 3 variables (name, cluster, maj.vote.index).
}

\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}

\examples{
library(Patterns)
if(require(CascadeData)){
data(micro_S, package="CascadeData")
D<-Patterns::as.micro_array(micro_S[1:100,],1:4,6)
a<-clustExploration(D)
a
}
}
\keyword{cluster}
