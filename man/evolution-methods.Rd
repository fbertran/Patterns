\name{evolution}
\alias{evolution}
\alias{evolution-methods}
\alias{evolution,network-method}

\title{See the evolution of the network with change of cutoff}
\description{
See the evolution of the network with change of cutoff
}
\usage{
	evolution(net,list_nv,\dots)
}

\arguments{
 \item{net}{a network object}
 \item{list_nv}{a vector of cutoff at which the network should be shown}
 \item{\dots}{Optionnal arguments:
 	\describe{ \item{gr}{a vector giving the group of each genee. Defaults to NULL}
  \item{color.vertex}{a vector giving the color of each nodee. Defaults to NULL}
  \item{color.edge}{a vector giving the color of each edge. Defaults to NULL}
  \item{fix}{logical, should the position of the node in the network be calculated once at the beginning ? Defaults to TRUE.}
  \item{size}{vector giving the size of the plot. Defaults to c(2000,1000)}
  \item{label_v}{vector giving the labels of each vertex. Defaults to 1:dim(net@network)[1]}
  \item{legend}{string giving the position of the legend. Defaults to "topleft"}
  \item{frame.color}{string giving the color of the frame of the plot. Defaults to "black"}
  \item{label.hub}{label hubs. Defaults to FALSE}
  \item{outdir}{Directory to save the animation. No default value since it must be specified by the user.}
  \item{type.ani}{Type of animation. Defaults to "html"}
}
}}

\details{
Several types of outputs are available using the type.ani option.
\itemize{ 
\item{html}
\item{latex (requires latex)}
\item{swf (requires swftools)}
\item{video (requires ffmpeg)}
\item{gif}
\item{manual_gif}
}
}

\value{ A HTML page with the evolution of the network.}


\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}



\examples{
\donttest{
data(network)
sequence<-seq(0,0.2,length.out=20)

destdir = tempdir()

#We set the working directory to a special directory only in the 
# case of online documentation creation. This should not be run 
# by a package user.
try(if(in_pkgdown()){destdir = 
  "~/Github/Patterns/docs/reference/evolution/"})
	
#Example of use of the evolution method with an html output.
evolution(network,sequence,type.ani = "html",outdir=destdir)

#Example of use of the evolution method with an animated gif output.
try(if(in_pkgdown()){evolution(network,sequence,type.ani = "gif",outdir=destdir)})
  
#The .gif creation is only working in the user working directory and not in temp directory.
\dontrun{
evolution(network,sequence,type.ani = "gif",outdir=getwd())
}
}
}

\keyword{methods}

