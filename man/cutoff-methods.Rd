\name{cutoff}
\alias{cutoff}
\alias{cutoff-methods}
\alias{cutoff,network-method}
\title{Choose the best cutoff}

\description{
Allows estimating the best cutoff. For a sequence of cutoff, the p value corresponding 
}
\usage{
	cutoff(Omega,...)
}

\arguments{
 \item{Omega}{a network object}
 \item{...}{Optional arguments: 
 	\describe{
    \item{sequence}{a vector corresponding to the sequence of cutoffs that will be tested.}
 	\item{x_min}{an integer ; only values over x_min are further retained for performing the test.}	
       }

 }
}

\value{
  A list containing two objects :
  \item{p.value}{the p values corresponding to the sequence of cutoff}
  \item{p.value.inter}{the smoothed p value vector, using the loess function}
}

\references{
Vallat, L., Kemper, C. A., Jung, N., Maumy-Bertrand, M., Bertrand, F., Meyer, N., ... & Bahram, S. (2013). Reverse-engineering the genetic circuitry of a cancer cell with predicted intervention in chronic lymphocytic leukemia. \emph{Proceedings of the National Academy of Sciences}, 110(2), 459-464.
}

\author{
Nicolas Jung, Bertrand Frederic, Myriam Maumy-Bertrand.
}


\examples{
\donttest{
		data(network)
		cutoff(network)
		#See vignette for more details
}
}

\keyword{methods}

