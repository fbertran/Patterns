#' @keywords internal
#' @aliases Patterns-package Patterns NULL
#'
#' @references F. Bertrand, I. Aouadi, N. Jung, R. Carapito, L. Vallat, S. Bahram, M. Maumy-Bertrand (2020). SelectBoost: a general algorithm to enhance the performance of variable selection methods, \emph{Bioinformatics}, \doi{10.1093/bioinformatics/btaa855}.
#' 
#' C. Schleiss, [...], M. Maumy-Bertrand, S. Bahram, F. Bertrand, and L. Vallat. (2021). Temporal multiomic modelling reveals a B-cell receptor proliferative program in chronic lymphocytic leukemia. \emph{Leukemia}.
#'
"_PACKAGE"


#' @importFrom grDevices col2rgb
#' @importFrom grDevices colorRamp
#' @importFrom grDevices dev.cur
# #' @importFrom grDevices dev.new
#' @importFrom grDevices grey
#' @importFrom grDevices rainbow
#' @importFrom grDevices rgb
#' @importFrom grDevices terrain.colors
#' @importFrom graphics abline
#' @importFrom graphics hist
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics matplot
#' @importFrom graphics par
#' @importFrom graphics rect
#' @importFrom graphics text
#' @importFrom graphics grid
#' @importFrom stats aggregate
#' @importFrom stats coef
#' @importFrom stats cor
#' @importFrom stats hclust
#' @importFrom stats lm
#' @importFrom stats loess
#' @importFrom stats model.matrix
#' @importFrom stats na.omit
#' @importFrom stats quantile
#' @importFrom stats rbinom
#' @importFrom stats reshape
#' @importFrom stats rmultinom
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats var
#' @importFrom stats wilcox.test
#' @importFrom methods new
#' @importFrom utils sessionInfo
#' @importFrom VGAM rlaplace
#' @importFrom VGAM rpareto
#' @importFrom VGAM zeta
#' @importFrom abind abind
#' @importFrom tnet as.tnet
#' @importFrom tnet betweenness_w
#' @importFrom tnet closeness_w
#' @importFrom tnet degree_w
#' @importFrom limma contrasts.fit
#' @importFrom limma eBayes
#' @importFrom limma lmFit
#' @importFrom limma makeContrasts
#' @importFrom limma topTable
#' @importFrom igraph get.edgelist 
#' @importFrom igraph graph.adjacency
#' @importFrom igraph layout.fruchterman.reingold
#' @importFrom igraph neighborhood
#' @importFrom igraph vcount
#' @importFrom igraph V
#' @importFrom nnls nnls
#' @importFrom plotrix color2D.matplot
#' 
#' @importFrom cluster agnes
#' @importFrom Mfuzz standardise
#' @importFrom Mfuzz mestimate
#' @importFrom Mfuzz cselection
#' @importFrom Mfuzz Dmin
#' @importFrom Mfuzz overlap.plot
#' @importFrom Mfuzz mfuzz.plot
#' @importFrom Mfuzz mfuzz
#' @importFrom e1071 cmeans
#' @importFrom SelectBoost fastboost
#' @importFrom SelectBoost group_func_2
#' @importFrom gplots heatmap.2
#' @importFrom jetset jscores
#' @exportMethod analyze_network
#' @exportMethod clustExploration
#' @exportMethod clustInference
#' @exportMethod compare
#' @exportMethod cutoff
#' @exportMethod dim
#' @exportMethod evolution
#' @exportMethod gene_expr_simulation
#' @exportMethod geneNeighborhood
#' @exportMethod genePeakSelection
#' @exportMethod geneSelection
#' @exportMethod head
#' @exportMethod inference
#' @exportMethod probeMerge
#' @exportMethod plot
#' @exportMethod position
#' @exportMethod predict
#' @exportMethod show
#' @exportMethod summary
#' @exportMethod unionOmics
#  #' @exportMethod unionNGseq
#  
#' @exportClass omics_array
# #' @exportClass nextgen_seq
#' @exportClass omics_network
#' @exportClass omics_predict
# #' @exportClass ngseqpredict
#' 
#' @export network_random
#' @export unsupervised_clustering_auto_m_c
#' @export unsupervised_clustering
#' @export as.omics_array
#' @export plotF
# #' @export "as.nextgen_seq
#' @export replaceBand
#' @export replaceUp
#' @export replaceDown
#' @export IndicFinit
#' @export IndicFshape
#' @export CascadeFshape
#' @export CascadeFinit
#' 
NULL

setGeneric("analyze_network",package="Patterns",def = function(Omega,nv,...){standardGeneric("analyze_network")})
setGeneric("clustExploration",package="Patterns",def = function(omicsarray,...){standardGeneric("clustExploration")})
setGeneric("clustInference",package="Patterns",def = function(omicsarray,vote.index,...){standardGeneric("clustInference")})
setGeneric("compare",package="Patterns",def = function(Net,Net_inf,nv){standardGeneric("compare")})
setGeneric("cutoff",package="Patterns",def = function(Omega,... ){standardGeneric("cutoff")})
setGeneric("evolution",package="Patterns",def = function(net,list_nv,... ){standardGeneric("evolution")})
setGeneric("gene_expr_simulation",package="Patterns",def = function(omics_network,...){standardGeneric("gene_expr_simulation")})
#setGeneric("gene_counts_simulation",package="Patterns",def = function(network,...){standardGeneric("gene_counts_simulation")})
setGeneric("geneNeighborhood",package="Patterns",def = function(net,targets,... ){standardGeneric("geneNeighborhood")})
setGeneric("genePeakSelection",package="Patterns",def = function(x,peak,... ){standardGeneric("genePeakSelection")})
setGeneric("geneSelection",package="Patterns",def = function(x,y,tot.number,... ){standardGeneric("geneSelection")})
setGeneric("inference",package="Patterns",def = function(M,... ){standardGeneric("inference")})
setGeneric("position",package="Patterns",def = function(net,... ){standardGeneric("position")})
#setGeneric("inferenceCascade",package="Patterns",def = function(M,... ){standardGeneric("inferenceCascade")})
#setGeneric("predict",package="Patterns",def = function(object,...){standardGeneric("predict")})
setGeneric("probeMerge",package="Patterns",def = function(x,...){standardGeneric("probeMerge")})
setGeneric("unionOmics",package="Patterns",def = function(M1,M2 ){standardGeneric("unionOmics")})
#setGeneric("unionNGseq",package="Patterns",def = function(C1,C2 ){standardGeneric("unionNGseq")})
setGeneric("unsupervised_clustering_auto_m_c",package="Patterns",def=function(M1,... ){standardGeneric("unsupervised_clustering_auto_m_c")})
setGeneric("unsupervised_clustering",package="Patterns",def=function(M1,clust,mestim,... ){standardGeneric("unsupervised_clustering")})

#' Class \code{"omics_array"}
#' 
#' The \code{"omics_array"} class
#' 
#' 
#' @name omics_array-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("omics_array", ...)}.
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords classes
#' @examples
#' 
#' showClass("omics_array")
#' 
setClass(Class = "omics_array",
         slots=c(omicsarray="matrix",                        
                        name="vector",
                        gene_ID="vector",
                        group=c("vector",NULL),
                        start_time=c("vector",NULL),
                        time=c("vector",NULL),
                        subject="numeric"),
         prototype = prototype(group = 0,start_time=0),
         validity=function(object){
           if(dim(object@omicsarray)[2] != length(object@time)*object@subject){
             stop("[Error: ]Number of colomns must be equal to the number of time points * the number of subject")
           }
           if(dim(object@omicsarray)[1] != length(object@name)&&length(object@name)!=0){
             stop("[Error: ] Length of the vector of names must equal to the number of genes")
           }
           if(dim(object@omicsarray)[1] != length(object@group)&&length(object@group)!=1){
             print(object@group)
             stop("[Error: ] Length of the vector of group must equal to the number of genes or null")
           }
           if(dim(object@omicsarray)[1] != length(object@start_time)&&length(object@start_time)!=1){
             stop("[Error: ] Length of the vector of starting time must equal to the number of genes or null")
           }
           if(object@subject<1){
             stop("[Error: ] There must be at least one subject")
           }	
         }
)

#' Class \code{"omics_network"}
#' 
#' The \code{"omics_network"} class
#' 
#' 
#' @name omics_network-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("omics_network", ...)}.
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords classes
#' @examples
#' 
#' showClass("omics_network")
#' 
setClass(
  Class = "omics_network",
  slots=c(omics_network = "matrix",
    name = "vector",
    F = "array",
    convF = "matrix",
    convO = "vector",
    time_pt = "vector"
  )
)

#' Class \code{"omics_predict"}
#' 
#' The \code{"omics_predict"} class
#' 
#' 
#' @name omics_predict-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("omics_predict", ...)}.
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords classes
#' @examples
#' 
#' showClass("omics_predict")
#' 
setClass(Class = "omics_predict",
         slots=c(omicsarray_unchanged="omics_array"
                        ,omicsarray_changed="omics_array"
                        ,omicsarray_predict="omics_array"
                        ,nv="numeric"
                        ,omics_network="omics_network"
                        ,targets="numeric")
)




