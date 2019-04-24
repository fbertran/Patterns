\name{compare}
\alias{compare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Some basic criteria of comparison between actual and inferred network.
}
\description{
Allows comparison between actual and inferred network.}
\usage{
compare(Net, Net_inf, nv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Net}{
A network object containing the actual network.
}
  \item{Net_inf}{
A network object containing the inferred network.
}
  \item{nv}{
A number that indicates at which level of cutoff the comparison should be done.
}
}

\value{
A vector containing : sensibility, predictive positive value, and the F-score
}

\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}

\examples{
	#See vignette
}

