\name{micro_array-class}
\Rdversion{1.1}
\docType{class}
\alias{micro_array-class}


\title{Class \code{"micro_array"}}
\description{
The Class
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("micro_array", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{microarray}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{name}:}{Object of class \code{"vector"} ~~ }
    \item{\code{gene_ID}:}{Object of class \code{"vector"} ~~ }
    \item{\code{group}:}{Object of class \code{"vector"} ~~ }
    \item{\code{start_time}:}{Object of class \code{"vector"} ~~ }
    \item{\code{time}:}{Object of class \code{"vector"} ~~ }
    \item{\code{subject}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{dim}{\code{signature(x = "micro_array")}: ... }
    \item{genePeakSelection}{\code{signature(M1 = "micro_array", M2 = "micro_array", peak = "numeric")}: ... }
    \item{geneSelection}{\code{signature(M1 = "micro_array", M2 = "micro_array", tot.number = "numeric")}: ... }
    \item{head}{\code{signature(x = "micro_array")}: ... }
    \item{inference}{\code{signature(M = "micro_array")}: ... }
    \item{plot}{\code{signature(x = "micro_array", y = "ANY")}: ... }
    \item{plot}{\code{signature(x = "network", y = "micro_array")}: ... }
    \item{predict}{\code{signature(object = "micro_array")}: ... }
    \item{print}{\code{signature(x = "micro_array")}: ... }
    \item{summary}{\code{signature(object = "micro_array")}: ... }
    \item{unionMicro}{\code{signature(M1 = "micro_array", M2 = "micro_array")}: ... }
	 }
}
\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}

\examples{
showClass("micro_array")
}
\keyword{classes}
