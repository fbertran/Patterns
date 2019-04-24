\name{CLL}
\alias{CLL}
\docType{data}
\title{
Expression data from healthy and malignant (chronic lymphocytic leukemia, CLL) human B-lymphocytes after B-cell receptor stimulation (GSE 39411 dataset)
}
\description{
B-cells were negatively selected from healthy donors and previously untreated CLL patients. BCR stimulated and unstimulated control B-cells were treated at four time points after stimulation for total RNA extraction and hybridization on Affymetrix microarrays.
}
\usage{data("CLL")}
\format{
  The format is:
 chr "CLL"
}
\details{
Three different cell populations (6 healthy B-lymphocytes, 6 leukemic CLL B-lymphocyte of indolent form and 5 leukemic CLL B-lymphocyte of aggressive form) were stimulated in vitro with an anti-IgM antibody, activating the B-cell receptor (BCR). We analyzed the gene expression at 4 time points (60, 90, 210 and 390 minutes). Each gene expression measurement is performed both in stimulated cells and in control unstimulated cells.
For one aggressive CLL case, we silenced expression of DUSP1 by transfecting DUSP1-specific RNAi and, as a control, transfected cells with a non-targeting RNAi. We then stimulated the BCR of these cells and analyzed the gene expression at the same time points in stimulated cells and in control unstimulated cells.
}
\source{
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39411
}
\references{
Vallat, L., Kemper, C. A., Jung, N., Maumy-Bertrand, M., Bertrand, F., Meyer, N., … Bahram, S. (2013). Reverse-engineering the genetic circuitry of a cancer cell with predicted intervention in chronic lymphocytic leukemia. Proceedings of the National Academy of Sciences of the United States of America, 110(2), 459–464. 
}
\examples{
data(CLL)
str(CLL)
}
\keyword{datasets}
