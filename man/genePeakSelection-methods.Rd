\name{genePeakSelection}
\alias{genePeakSelection}
\alias{geneSelection}
\alias{genePeakSelection-methods}
\alias{geneSelection-methods}
\alias{geneSelection,micro_array,numeric-method}
\alias{genePeakSelection,micro_array,numeric-method}
\title{Methods for selecting genes }
\description{
Find the neighborhood of a set of nodes. Please refer to the vignette for more details.
}
\usage{
	geneSelection(M1,tot.number,...)
	genePeakSelection(M1,pic,... )
}

\arguments{
 \item{M1}{a micro_array object. The stimulated measurements.}     
 \item{tot.number}{an integer. The number of selected genes. [only for geneSelection]}
\item{pic}{interger. At which time points measurements should the genes be selected [optionnal for geneSelection].}

 \item{...}{Optional arguments: \describe{
  \item{M2}{a micro_array object. The unstimulated measurements.}
 		 		\item{data_log}{logical (default to TRUE) ; should data be logged ?}
 		\item{wanted.patterns}{a matrix with wanted patterns [only for geneSelection].}
 		\item{forbidden.patterns}{a matrix with forbidden patterns [only for geneSelection].}
 		\item{durPic}{vector of size 2 (default to c(1,1)) ; the first elements gives the length of the pic at the left, the second at the right. [only for genePeakSelection]}
 		\item{abs_val}{ logical (default to TRUE) ; should genes be selected on the basis of their absolute value expression ? [only for genePeakSelection]}
 		\item{alpha_diff}{float ; the risk level }
 		
 }}
}

\value{
  A micro_array object.  
}

\examples{
	#See vignette.
}

\keyword{methods}
