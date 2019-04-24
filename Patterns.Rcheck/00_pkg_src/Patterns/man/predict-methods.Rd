\name{predict}
\alias{predict}
\alias{predict-methods}
\alias{predict,ANY-method}
\alias{predict,micro_array-method}
\title{  Methods for Function \code{predict}  }
\description{
 Prediction of the gene expressions after a knock-out experience 
}
\usage{
	predict(object,...)
}

\arguments{
 \item{object}{a micro_array object}
 \item{...}{Other arguments: 
 	\describe{
    \item{Omega}{a netowork object.}
 	\item{nv}{ [=0] numeric ; the level of the cutoff}	
 	 	\item{targets}{ [NULL] vector ; which genes are knocked out ?}
       }
  }
}

\examples{
\donttest{
data(Selection)
data(network)
#A nv value can chosen using the cutoff function
nv=.11
EGR1<-which(match(Selection@gene_ID,"EGR1")==1)
P<-position(network,nv=nv)

#We predict gene expression modulations within the network if EGR1 is experimentaly knocked-out. 
prediction_ko5<-predict(Selection,network@network,nv=nv,targets=EGR1)

#Then we plot the results. Here for example we see changes at time point t2:
plot(prediction_ko5,time=2,ini=P,label_v=Selection@name)
}
}
\keyword{methods}

