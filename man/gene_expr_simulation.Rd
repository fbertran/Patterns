\name{gene_expr_simulation}
\alias{gene_expr_simulation}
\alias{gene_expr_simulation-methods}
\alias{gene_expr_simulation,network-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulates microarray data based on a given network. 
}
\description{
Simulates microarray data based on a given network. 
}
\usage{
gene_expr_simulation(network,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{network}{A network object.}
\item{...}{\describe{
  \item{time_label}{a vector containing the time labels.}
  \item{subject}{the number of subjects}
  \item{peak_level}{the mean level of peaks. }
	\item{act_time_group}{ [NULL] vector ; at which time the groups (defined by sort(unique(group))) are activated ?}
}}
}

\value{
A micro_array object.
}

\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}



\examples{
data(Net)
set.seed(1)

#We simulate gene expressions according to the network Net
Msim<-Patterns::gene_expr_simulation(
	network=Net,
	time_label=rep(1:4,each=25),
	subject=5,
	peak_level=200)
head(Msim)
}
