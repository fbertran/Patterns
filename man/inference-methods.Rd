\name{inference}
\alias{inference}
\alias{inference-methods}
\alias{inference,micro_array-method}
\title{ Reverse-engineer the network}
\description{
Reverse-engineer the network.
}
\usage{
	inference(M,\dots)
}

\arguments{
 \item{M}{a micro_array object.}
  \item{\dots}{Optional arguments: 	 
  \describe{
\item{M}{Data}
\item{tour.max=30}{tour.max + 1 = maximal number of steps.}
\item{g=function(x){1/x}}{After each step, the new solution is choosen as (the old solution + g(x) * the new solution)/(1+g(x)) where x is the number of steps.}
\item{conv=0.001}{Convergence criterion.}
\item{cv.subjects=TRUE}{Subjectwise cross validation: should the cross validation be done by removing the subject one by one?}
\item{nb.folds=NULL}{Relevant only if no subjectwise cross validation (i.e. cv.subjects=FALSE). The number of folds in cross validation.}
\item{eps=10^-5}{Threshold for rounding coefficients to 0 (i.e. machine zero).}
\item{type.inf="iterative"}{"iterative" or "noniterative" : should the algorithm be computed iteratively or only for one step? For highly homogeneous clusters, the "noniterative" option is suffisant.}
\item{Fshape=NULL}{Shape of the F matrix.}
\item{Finit=NULL}{Init values of the F matrix.}
\item{Omega=NULL}{Init values for the Omega matrix.}
\item{fitfun="LASSO"}{Function to infer the Omega matrix at each step.}
\item{use.Gram=TRUE}{Optional parameter for the lasso in the `lars` package.}
\item{error.stabsel=0.05}{Optional parameter for the stability selection algorithm in the `c060` package.}
\item{pi_thr.stabsel=0.6}{Optional parameter for the stability selection algorithm in the `c060` package.}
\item{priors=NULL}{A priori weights for the links between genes.}
\item{mc.cores=getOption("mc.cores", 2L)}{Number of cores.}
\item{intercept.stabpath=TRUE}{Use intercept in stability selection models?}
\item{steps.seq=.95}{Optional parameter for the SelectBoost algorithm in the `SelectBoost` package.}
\item{limselect=.95}{Optional parameter for the SelectBoost algorithm in the `SelectBoost` package.}
\item{use.parallel=TRUE}{Use parallel computing?}
}
}
}

\details{
The fitting built-in fitting functions (`fitfun`) provided with the `Patterns` package are :
\describe{
\item{LASSO}{from the `lars` package (default value)}
\item{LASSO2}{from the `glmnet` package}
\item{SPLS}{from the `spls` package}
\item{ELASTICNET}{from the `elasticnet` package}
\item{stability.c060}{from the `c060` package implementation of stability selection}
\item{stability.c060.weighted}{a new weighted version of the `c060` package implementation of stability selection}
\item{robust}{lasso from the `lars` package with light random Gaussian noise added to the explanatory variables}
\item{selectboost.weighted}{a new weighted version of the `selectboost` package implementation of the selectboost algorithm to look for the more stable links against resampling that takes into account the correlated structure of the predictors. If no weights are provided, equal weigths are for all the variables (=non weighted case).}
}
}


\value{
A network object.  
}


\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}


\examples{
	#data(micro_US)
	#inference(micro_US)
	#See vignette for more details
}

\keyword{methods}

