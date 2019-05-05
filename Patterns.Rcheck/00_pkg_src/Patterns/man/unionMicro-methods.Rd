\name{unionMicro-methods}
\docType{methods}
\alias{unionMicro}
\alias{unionMicro-methods}
\alias{unionMicro,micro_array,micro_array-method}
\alias{unionMicro,list,ANY-method}
\title{Makes the union between two micro_array objects.}
\description{
Makes the union between two micro_array objects.
}
\section{Methods}{
\describe{

\item{\code{signature(M1 = "micro_array", M2 = "micro_array")}}{
Returns a micro_array object which is the union of M1 and M2. 
}

\item{\code{signature(M1 = "list", M2 = "ANY")}}{
Returns a micro_array object which is the union of the elements of M1. 
}
}}
\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}

\examples{
if(require(CascadeData)){
data(micro_S, package="CascadeData")
#Create another microarray object with 100 genes
Mbis<-M<-as.micro_array(micro_S[1:100,],1:4,6)
#Rename the 100 genes
Mbis@name<-paste(M@name,"bis")
rownames(Mbis@microarray) <- Mbis@name
#Union (merge without duplicated names) of the two microarrays. 
str(unionMicro(M,Mbis))
}
}
\keyword{methods}

