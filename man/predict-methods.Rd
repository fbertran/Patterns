\name{predict}
\alias{predict}
\alias{predict-methods}
\alias{predict,ANY-method}
\alias{predict,micro_array-method}
\title{  Methods for Function \code{predict}  }
\description{
Prediction for knock out intervention in the gene regulatory network.
}

\usage{
	predict(object,Omega,nv=0,target=NULL)
}

\arguments{
 \item{object}{a micro_array object}
 \item{Omega}{a network object.}
 	\item{nv}{ [=0] numeric ; the level of the cutoff}	
 	 	\item{targets}{ [NULL] vector ; which genes are knocked out ?}
       }

\keyword{methods}

