\name{network_random}
\alias{network_random}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates a network.
}
\description{
Generates a network.
}
\usage{
network_random(nb, time_label, exp, init, regul, min_expr, max_expr, casc.level)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nb}{
Integer. The number of genes.
}
  \item{time_label}{
Vector. The time points measurements.
}
  \item{exp}{
The exponential parameter, as in the barabasi.game function in igraph package.
}
  \item{init}{
The attractiveness of the vertices with no adjacent edges. See barabasi.game function.
}
  \item{regul}{
A vector mapping each gene with its number of regulators.
}
  \item{min_expr}{
Minimum of strength of a non-zero link
}
  \item{max_expr}{
Maximum of strength of a non-zero link
}
  \item{casc.level}{
...}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A network object. 
}
\references{
Vallat, L., Kemper, C. A., Jung, N., Maumy-Bertrand, M., Bertrand, F., Meyer, N., ... & Bahram, S. (2013). Reverse-engineering the genetic circuitry of a cancer cell with predicted intervention in chronic lymphocytic leukemia. \emph{Proceedings of the National Academy of Sciences}, 110(2), 459-464.
}

\author{
Nicolas Jung, Bertrand Frederic, Myriam Maumy-Bertrand.
}


\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
	#See vignette
}