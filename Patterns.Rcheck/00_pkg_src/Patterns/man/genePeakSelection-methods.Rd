\name{genePeakSelection}
\alias{genePeakSelection}
\alias{geneSelection}
\alias{genePeakSelection-methods}
\alias{geneSelection-methods}
\alias{geneSelection,micro_array,numeric-method}
\alias{geneSelection,micro_array,micro_array,numeric-method}
\alias{geneSelection,list,list,numeric-method}
\alias{genePeakSelection,micro_array,numeric-method}
\title{Methods for selecting genes }
\description{
Selection of differentially expressed genes.
}
\usage{
	geneSelection(x,y,tot.number,...)
	genePeakSelection(x,peak,...)
}   
\arguments{
\item{x}{either a micro_array object or a list of micro_array objects. In the first case, the micro_array object represents the stimulated measurements. In the second case, the control unstimulated data (if present) should be the first element of the list.}     
\item{y}{either a micro_array object or a list of strings. In the first case, the micro_array object represents the stimulated measurements. In the second case, the list is the way to specify the contrast: \describe{
\item{First element:}{ condition, condition&time or pattern. The condition specification is used when the overall is to compare two conditions. 
The condition&time specification is used when comparing two conditions at two precise time points.
The pattern specification allows to decide which time point should be differentially expressed.}
\item{Second element:}{a vector of length 2. The two conditions which should be compared. If a condition is used as control, it should be the first element of the vector. However, if this control is not measured throught time, 
the option cont=TRUE should be used.}
\item{Third element:}{depends on the first element. 
It is no needed if condition has been specified. 
If condition&time has been specified, then this is a vector containing the time point at which the comparison should be done.
If pattern has been specified, then this is a vector of 0 and 1 of length T, where T is the number of time points. The time points with desired differential expression are provided with 1. 
}}}      
 \item{tot.number}{an integer. The number of selected genes. If tot.number <0 all differentially genes are selected. If tot.number > 1, tot.number is the maximum of diffenrtially genes that will be selected. 
 If 0<tot.number<1, tot.number represents the proportion of diffenrentially genes that are selected.}
\item{peak}{interger. At which time points measurements should the genes be selected [optionnal for geneSelection].}
 \item{...}{Optional arguments: \describe{
  \item{M2}{a micro_array object. The unstimulated measurements.}
 		 		\item{data_log}{logical (default to TRUE) ; should data be logged ?}
 		\item{wanted.patterns}{a matrix with wanted patterns [only for geneSelection].}
 		\item{forbidden.patterns}{a matrix with forbidden patterns [only for geneSelection].}
 		\item{durPeak}{vector of size 2 (default to c(1,1)) ; the first elements gives the length of the peak at the left, the second at the right. [only for genePeakSelection]}
 		\item{abs_val}{ logical (default to TRUE) ; should genes be selected on the basis of their absolute value expression ? [only for genePeakSelection]}
 		\item{alpha_diff}{float ; the risk level }
 }}
}

\value{
  A micro_array object.  
}

\author{
Nicolas Jung, Frédéric Bertrand , Myriam Maumy-Bertrand.
}


\references{
Jung, N., Bertrand, F., Bahram, S., Vallat, L., and Maumy-Bertrand, M. (2014). Cascade: a R-package to study, predict and simulate the diffusion of a signal through a temporal gene network. \emph{Bioinformatics}, btt705. 
	
Vallat, L., Kemper, C. A., Jung, N., Maumy-Bertrand, M., Bertrand, F., Meyer, N., ... & Bahram, S. (2013). Reverse-engineering the genetic circuitry of a cancer cell with predicted intervention in chronic lymphocytic leukemia. \emph{Proceedings of the National Academy of Sciences}, 110(2), 459-464.
}

\examples{
\donttest{
  if(require(CascadeData)){
	data(micro_US)
	micro_US<-as.micro_array(micro_US,time=c(60,90,210,390),subject=6)
	data(micro_S)
	micro_S<-as.micro_array(micro_S,time=c(60,90,210,390),subject=6)

  #Basically, to find the 50 more significant expressed genes you will use:
  Selection_1<-geneSelection(x=micro_S,y=micro_US,
  tot.number=50,data_log=TRUE)
  summary(Selection_1)
  
  #If we want to select genes that are differentially 
  #at time t60 or t90 :
  Selection_2<-geneSelection(x=micro_S,y=micro_US,tot.number=30,
  wanted.patterns=
  rbind(c(0,1,0,0),c(1,0,0,0),c(1,1,0,0)))
  summary(Selection_2)

  #To select genes that have a differential maximum of expression at a specific time point.
  
  Selection_3<-genePeakSelection(x=micro_S,y=micro_US,peak=1,
  abs_val=FALSE,alpha_diff=0.01)
  summary(Selection_3)
  }

  if(require(CascadeData)){
data(micro_US)
micro_US<-as.micro_array(micro_US,time=c(60,90,210,390),subject=6)
data(micro_S)
micro_S<-as.micro_array(micro_S,time=c(60,90,210,390),subject=6)
#Genes with differential expression at t1
Selection1<-geneSelection(x=micro_S,y=micro_US,20,wanted.patterns= rbind(c(1,0,0,0)))
#Genes with differential expression at t2
Selection2<-geneSelection(x=micro_S,y=micro_US,20,wanted.patterns= rbind(c(0,1,0,0)))
#Genes with differential expression at t3
Selection3<-geneSelection(x=micro_S,y=micro_US,20,wanted.patterns= rbind(c(0,0,1,0)))
#Genes with differential expression at t4
Selection4<-geneSelection(x=micro_S,y=micro_US,20,wanted.patterns= rbind(c(0,0,0,1)))
#Genes with global differential expression 
Selection5<-geneSelection(x=micro_S,y=micro_US,20)

#We then merge these selections:
Selection<-unionMicro(list(Selection1,Selection2,Selection3,Selection4,Selection5))
print(Selection)

#Prints the correlation graphics Figure 4:
summary(Selection,3)

##Uncomment this code to retrieve geneids.
#library(org.Hs.eg.db)
#
#ff<-function(x){substr(x, 1, nchar(x)-3)}
#ff<-Vectorize(ff)
#
##Here is the function to transform the probeset names to gene ID.
#
#library("hgu133plus2.db")
#
#probe_to_id<-function(n){  
#x <- hgu133plus2SYMBOL
#mp<-mappedkeys(x)
#xx <- unlist(as.list(x[mp]))
#genes_all = xx[(n)]
#genes_all[is.na(genes_all)]<-"unknown"
#return(genes_all)
#}
#Selection@name<-probe_to_id(Selection@name)
  }
	}
}
\keyword{methods}
