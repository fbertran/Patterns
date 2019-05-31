\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,micro_array,ANY-method}
\alias{plot,micropredict,ANY-method}
\alias{plot,network,ANY-method}
\title{Plot}
\description{
Considering the class of the argument which is passed to plot, the graphical output differs. 
}
\section{Methods}{
\describe{

\item{\code{signature(x = "micro_array", y = "ANY",...)}}{
\describe{
 \item{x}{a micro\_array object}
 \item{list_nv}{a vector of cutoff at which the network should be shown}
 }
}                        
\item{\code{signature(x = "network", y = "ANY",...)}}{
\describe{
 \item{x}{a network object}
 \item{\dots}{Optionnal arguments:
 	\describe{ \item{gr}{a vector giving the group of each gene}
 	\item{choice}{what graphic should be plotted: either "F" (for a representation of the matrices F) or "network".}
  \item{nv}{the level of cutoff. Defaut to 0.}
  \item{ini}{using the ``position'' function, you can fix the position of the nodes}
  \item{color.vertex}{a vector defining the color of the vertex}
 \item{ani}{vector giving the size of the plot. Default to c(2000,1000). The animation can only be created in the working directory. See the help page of the animation method.}
 \item{video}{if ani is TRUE and video is TRUE, the animation result is a GIF video}
 \item{label_v}{vector defining the vertex labels}
\item{legend.position}{position of the legend}
 \item{frame.color}{color of the frames}
\item{label.hub}{logical ; if TRUE only the hubs are labeled}                                     
 \item{edge.arrow.size}{size of the arrows ; default to 0.7}
\item{edge.thickness}{edge thickness ; default to 1.}   
}
}}}

\item{\code{signature(x = "micropredict", y = "ANY",...)}}{
\describe{
 \item{x}{a micropredict object}
 \item{\dots}{Optionnal arguments: see plot for network}
}}
}
}
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
\keyword{methods}

