\name{compare-methods}
\docType{methods}
\alias{compare-methods}
\alias{compare}
\alias{compare,network,network,numeric-method}
\title{
Some basic criteria of comparison between actual and inferred network.
}
\description{
Allows comparison between actual and inferred network.}
\section{Methods}{
\describe{

\item{\code{signature(Net = "network", Net_inf = "network", nv = "numeric")}}{
\describe{
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
}                        

}
}


\value{
A vector containing : sensitivity, predictive positive value, the usual F-score (2*ppv*sens/(sppvpe+sens)), the 1/2 ponderated Fscore ((1+0.5^2)*ppv*sens/(ppv/4+sens)) and the 2 ponderated Fscore ((1+2^2)*ppv*sens/(ppv*4+sens)).
}

\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}

\examples{
data(simul)

#Comparing true and inferred networks
Crit_values=NULL

#Here are the cutoff level tested
test.seq<-seq(0,max(abs(Net_inf_PL@network*0.9)),length.out=200)
for(u in test.seq){
	Crit_values<-rbind(Crit_values,Patterns::compare(Net,Net_inf_PL,u))
}
matplot(test.seq,Crit_values,type="l",ylab="Criterion value",xlab="Cutoff level",lwd=2)
legend(x="topleft", legend=colnames(Crit_values), lty=1:5,col=1:5,ncol=2,cex=.9)
}

