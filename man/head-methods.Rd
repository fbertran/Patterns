\name{head}
\alias{methods}
\alias{head-methods}
\alias{head,ANY-method}
\alias{head,micro_array-method}
\title{Overview of a micro_array object}
\description{
Overview of a micro_array object.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "ANY")}}{
Gives an overview.
}

\item{\code{signature(x = "micro_array")}}{
Gives an overview.
}
}}

\examples{
  if(require(CascadeData)){
	data(micro_US)
	micro_US<-as.micro_array(micro_US,time=c(60,90,210,390),subject=6)
	head(micro_US)
	}
}

\keyword{methods}

