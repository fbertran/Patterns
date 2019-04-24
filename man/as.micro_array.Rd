\name{as.micro_array}
\alias{as.micro_array}

\title{
Coerce a matrix into a micro_array object.
}
\description{
Coerce a matrix into a micro_array object.
}
\usage{
as.micro_array(M, time, subject, name_probe = NULL, gene_ID = NULL)
}

\arguments{
  \item{M}{
A matrix. Contains the microarray measurements. Should of size N * K, with N the number of genes and K=T*P with T the number of time points, and P the number of subjects. This matrix should be created using cbind(M1,M2,...) with M1 a N*T matrix with the measurements for patient 1, M2 a N*T matrix with the measurements for patient 2.
}
  \item{time}{
A vector. The time points measurements.
}
  \item{subject}{
The number of subjects. 
}
  \item{name_probe}{
Vector with the row names of the micro_array. 
}
  \item{gene_ID}{
Vector with the gene IDs of the row names of the micro array. 
}
}

\value{A micro_array object.}

\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}




\examples{
if(require(CascadeData)){
	data(micro_US, package="CascadeData")
	micro_US<-as.micro_array(micro_US[1:100,],time=c(60,90,210,390),subject=6)
	plot(micro_US)
	}
}
