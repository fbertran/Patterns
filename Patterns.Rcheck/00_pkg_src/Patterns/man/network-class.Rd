\name{network-class}
\Rdversion{1.1}
\docType{class}
\alias{network-class}


\title{Class \code{"network"}}
\description{
2254
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("network", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{network}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{name}:}{Object of class \code{"vector"} ~~ }
    \item{\code{F}:}{Object of class \code{"array"} ~~ }
    \item{\code{convF}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{convO}:}{Object of class \code{"vector"} ~~ }
    \item{\code{time_pt}:}{Object of class \code{"vector"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{analyze_network}{\code{signature(Omega = "network")}: ... }
    \item{cutoff}{\code{signature(Omega = "network")}: ... }
    \item{evolution}{\code{signature(net = "network")}: ... }
    \item{geneNeighborhood}{\code{signature(net = "network")}: ... }
    \item{plot}{\code{signature(x = "network", y = "ANY")}: ... }
    \item{plot}{\code{signature(x = "network", y = "micro_array")}: ... }
    \item{position}{\code{signature(net = "network")}: ... }
    \item{print}{\code{signature(x = "network")}: ... }
	 }
}

\examples{
showClass("network")
}
\keyword{classes}
