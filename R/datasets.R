#' Expression data from healthy and malignant (chronic lymphocytic leukemia,
#' CLL) human B-lymphocytes after B-cell receptor stimulation (GSE 39411
#' dataset)
#' 
#' B-cells were negatively selected from healthy donors and previously
#' untreated CLL patients. BCR stimulated and unstimulated control B-cells were
#' treated at four time points after stimulation for total RNA extraction and
#' hybridization on Affymetrix microarrays.
#' 
#' The dataset provided with package is the first five lines of the full
#' dataset. The full dataset can be downloaded from the github repository of
#' the package
#' (https://raw.githubusercontent.com/fbertran/Patterns/master/add_data/CLL.RData).
#' 
#' Three different cell populations (6 healthy B-lymphocytes, 6 leukemic CLL
#' B-lymphocyte of indolent form and 5 leukemic CLL B-lymphocyte of aggressive
#' form) were stimulated in vitro with an anti-IgM antibody, activating the
#' B-cell receptor (BCR). We analyzed the gene expression at 4 time points (60,
#' 90, 210 and 390 minutes). Each gene expression measurement is performed both
#' in stimulated cells and in control unstimulated cells. For one aggressive
#' CLL case, we silenced expression of DUSP1 by transfecting DUSP1-specific
#' RNAi and, as a control, transfected cells with a non-targeting RNAi. We then
#' stimulated the BCR of these cells and analyzed the gene expression at the
#' same time points in stimulated cells and in control unstimulated cells.
#' 
#' @name CLL
#' @docType data
#' @format The format is: chr "CLL"
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @references Vallat, L., Kemper, C. A., Jung, N., Maumy-Bertrand, M.,
#' Bertrand, F., Meyer, N., … Bahram, S. (2013). Reverse-engineering the
#' genetic circuitry of a cancer cell with predicted intervention in chronic
#' lymphocytic leukemia. Proceedings of the National Academy of Sciences of the
#' United States of America, 110(2), 459–464.
#' @source https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39411
#' @keywords datasets
#' @examples
#' 
#' data(CLL)
#' str(CLL)
#' 
#' \donttest{
#' CLLfile <- "https://raw.githubusercontent.com/fbertran/Patterns/master/add_data/CLL.RData"
#' load(CLLfile)
#' 
#' str(CLL)
#' }
#' 
NULL

#' Details on some probesets of the affy_hg_u133_plus_2 platform.
#' 
#' Dataset with information on the affy_hg_u133_plus_2 platform such as
#' probeset name (affy_hg_u133_plus_2), ensembl_gene_id, entrezgene,
#' hgnc_symbol, chromosome_name, start_position, end_position and band.
#' 
#' Data.frame with 8859 rows and 8 variables.
#' 
#' @name infos
#' @docType data
#' @format The format is: chr "infos"
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords datasets
#' @examples
#' 
#' data(infos)
#' 
NULL

#' Simulated microarray.
#' 
#' Simulated M, microarray.
#' 
#' 
#' @name M
#' @docType data
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords datasets
#' @examples
#' 
#' data(M)
#' head(M)
#' str(M)
#' 
NULL

#' Reverse-engineered network of the M and Net simulated data.
#' 
#' The reverse-engineered network with the `Patterns` package using the
#' fitfun="LASSO" default function and a cascade network setting.
#' 
#' 
#' @name Net_inf_PL
#' @docType data
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords datasets
#' @examples
#' 
#' data(Net_inf_PL)
#' str(Net_inf_PL)
#' 
NULL

#' Simulated network for examples.
#' 
#' Simulated network.
#' 
#' 
#' @name Net
#' @docType data
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords datasets
#' @examples
#' 
#' data(Net)
#' str(Net)
#' 
NULL

#' A example of an inferred network (4 groups case).
#' 
#' This dataset is a network example with 102 nodes, 4 times and 4 groups.
#' 
#' A network class object [package "Patterns"] with 6 slots.
#' 
#' @name network
#' @docType data
#' @format The format is: chr "network"
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords datasets
#' @examples
#' 
#' data(network)
#' str(network)
#' plot(network)
#' 
NULL

#' A example of an inferred cascade network (2 groups case).
#' 
#' This dataset is a cascade network example with 53 nodes, 4 times and 2
#' groups.
#' 
#' A network class object [package "Patterns"] with 6 slots.
#' 
#' @name network2gp
#' @docType data
#' @format The format is: chr "network2gp"
#' @keywords datasets
#' @examples
#' 
#' data(network2gp)
#' str(network2gp)
#' plot(network2gp)
#' 
NULL

#' A example of an inferred cascade network (4 groups case).
#' 
#' This dataset is a cascade network example with 102 nodes, 4 times and 4
#' groups.
#' 
#' A network class object [package "Patterns"] with 6 slots.
#' 
#' @name networkCascade
#' @docType data
#' @format The format is: chr "networkCascade"
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords datasets
#' @examples
#' 
#' data(networkCascade)
#' str(networkCascade)
#' plot(networkCascade)
#' 
NULL

#' Selection of genes.
#' 
#' 20 (at most) genes with differential expression at t1, 20 (at most) genes
#' with differential expression at t2, 20 (at most) genes with differential
#' expression at t3, 20 (at most) genes with differential expression at t4 et
#' 20 (at most) genes with global differential expression were selected.
#' 
#' 
#' @name Selection
#' @docType data
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords datasets
#' @examples
#' 
#' data(Selection)
#' head(Selection)
#' summary(Selection,3)
#' 
NULL





