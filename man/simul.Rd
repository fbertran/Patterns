\name{simul}
\alias{simul}
\alias{M}
\alias{Net}
\alias{Net_inf_PL}
\docType{data}
\title{
Simulated data for examples.
}
\description{
Three objects: M, microarray, Net, network, and Net_inf_PL, the reverse-engineered network with the `Patterns` package using the fitfun="LASSO" default function and a cascade network setting.
}
\usage{data("simul")}
\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}

\examples{
data(simul)
head(M)
str(M)
str(Net)
str(Net_inf_PL)
}
\keyword{datasets}