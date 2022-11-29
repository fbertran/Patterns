#' \code{Show} methods
#' 
#' Methods for generic function \code{show}
#' 
#' 
#' @name show-methods
#' @aliases show-methods show,ANY-method show,omics_array-method
#' show,omics_network-method
#' @docType methods
#' @section Methods: \describe{
#' \item{list("signature(object = \"ANY\")")}{ }
#' \item{list("signature(object = \"omics_array\")")}{ Print an object of class omics_array }
#' \item{list("signature(object = \"omics_network\")")}{ Print an object of class omics_network } 
#' }
#' 
#' @param object an object of class omics-array or omics_network
#' 
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords methods
#' 
setMethod("show","omics_array",function(object)
{
  cat(paste("This is a omics_array S4 class. It contains : \n - (@omicsarray) a matrix of dimension ",
            dim(object@omicsarray)[1],"*",dim(object@omicsarray)[2],
            "\n          .... [omics expressions (or abundancies)] \n - (@name) a vector of length ",
            length(object@name),
            " .... [omics names (i.e. gene or peptides names)] \n","\n - (@gene_ID) a vector of length ",
            length(object@gene_ID),
            " .... [omics ID] \n","- (@group) a vector of length ",
            length(object@group),
            " .... [groups for omics] \n","- (@start_time) a vector of length ",
            length(object@start_time),
            "\n          .... [first differential expression for omics] \n",
            "- (@time)a vector of length ",
            length(object@time),
            " .... [time points]\n",
            "- (@subject) an integer  .... [number of subject]")) 	
  invisible(object)
})

#' Overview of a omics_array object
#' 
#' Overview of a omics_array object.
#' 
#' 
#' @aliases methods head-methods head,ANY-method head,omics_array-method
#' @section Methods: \describe{
#' 
#' \item{list("signature(x = \"ANY\")")}{ Gives an overview. }
#' 
#' \item{list("signature(x = \"omics_array\")")}{ Gives an overview. } }
#' @param x an object of class `omics_array`.
#' @param ... additional parameters
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords methods
#' @examples
#' 
#'  if(require(CascadeData)){
#' 	data(micro_US)
#' 	micro_US<-as.omics_array(micro_US,time=c(60,90,210,390),subject=6)
#' 	head(micro_US)
#' 	}
#' 
setMethod("head",signature="omics_array",function(x,...)
{
  K<-dim(x@omicsarray)[2]
  K<-min(K,3)
  return(list(omicsarray=head(x@omicsarray[,1:K]),
     name=head(x@name),
     gene_ID=head(x@gene_ID),
     group=head(x@group),
     start_time=head(x@start_time),
     time=x@time,
     subject=x@subject))
}
)


#' \code{Summary} methods
#' 
#' Methods for function \code{summary}
#' 
#' 
#' @name summary-methods
#' @aliases summary-methods summary,ANY-method summary,omics_array-method
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(object = \"ANY\")")}{ %% ~~describe this method here~~
#' }
#' 
#' \item{list("signature(object = \"omics_array\")")}{ %% ~~describe this
#' method here~~ } }
#' 
#' @param object an object of class omics-array
#' @param nb.graph (optionnal) choose the graph to plot. Displays all graphs by default.
#' @param ... additional parameters.
#' 
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords methods
#' 
setMethod("summary","omics_array",function(object,nb.graph=NULL,...)
{
  require(cluster)
  print(summary(object@omicsarray))
  G<-object@omicsarray
  colnames(G)<-paste(rep(paste("T",object@time),object@subject), as.character(rep(paste("subject",1:object@subject),each=length(object@time))))
  z <- cor(G)
  require(lattice)
  #rownames(z)<-NULL
  if(is.null(nb.graph) || nb.graph==1 ){
    print(lattice::levelplot(z,aspect="iso",xlab="",ylab="", 
                    scales=list(x=list(rot=90)), 
                    ylab.right = "Level of correlation",
                    par.settings = list(layout.widths = list(axis.key.padding = 0,
                                                             ylab.right = 2))))
  }
  #if(is.null(nb.graph)){dev.new()}
  G<-object@omicsarray
  colnames(G)<-paste(rep(paste("T",object@time),object@subject), as.character(rep(paste("subject",1:object@subject),each=length(object@time))))
  w<-cluster::agnes(t(G))[1]$order
  G<-G[,w]
  z <- cor(G)
  
  if(is.null(nb.graph) || nb.graph==2 ){
    print(lattice::levelplot(z,aspect="iso",xlab="",ylab="", 
                    scales=list(x=list(rot=90)), 
                    ylab.right = "Level of correlation",
                    par.settings = list(layout.widths = list(axis.key.padding = 0,
                                                             ylab.right = 2))))
  }
  if(dim(object@omicsarray)[1]<1000){
    #if(is.null(nb.graph)){dev.new()}
    
    R<-object@omicsarray
    w<-cluster::agnes(R)[1]$order
    R<-R[w,]
    z <- cor(t(R))
    
    if(is.null(nb.graph) || nb.graph==3 ){
      print(lattice::levelplot(z,aspect="iso",xlab="",ylab="", 
                      scales=list(x=list(rot=90)), 
                      ylab.right = "Level of correlation",
                      par.settings = list(layout.widths = list(axis.key.padding = 0,
                                                               ylab.right = 2))))
    }
  }}
)

#' Dimension of the data
#' 
#' Dimension of the data
#' 
#' 
#' @name dim
#' @aliases dim dim-methods dim,omics_array-method
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(x = \"omics_array\")")}{ Gives the dimension of the
#' matrix of measurements. } }
#' @param x an object of class `omics_array`.
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords methods
setMethod("dim","omics_array",function(x)
{
  return(dim(x@omicsarray))
})




#' Plot
#' 
#' Considering the class of the argument which is passed to plot, the graphical
#' output differs.
#' 
#' 
#' @name plot-methods
#' @aliases plot-methods plot,omics_array,ANY-method
#' plot,omics_predict,ANY-method plot,omics_network,ANY-method
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(x = \"omics_array\", y = \"ANY\",...)")}{ \describe{
#' \item{x}{a omics_array object} \item{list_nv}{a vector of cutoff at which
#' the network should be shown} } } \item{list("signature(x = \"omics_network\", y =
#' \"ANY\",...)")}{ \describe{ \item{x}{a omics_network object}
#' \item{list()}{Optionnal arguments: \describe{ \item{gr}{a vector giving the
#' group of each gene} \item{choice}{what graphic should be plotted: either "F"
#' (for a representation of the matrices F) or "network".} \item{nv}{the level
#' of cutoff. Defaut to 0.} \item{ini}{using the ``position'' function, you can
#' fix the position of the nodes} \item{color.vertex}{a vector defining the
#' color of the vertex} \item{ani}{vector giving the size of the plot. Default
#' to c(2000,1000). The animation can only be created in the working directory.
#' See the help page of the animation method.} \item{video}{if ani is TRUE and
#' video is TRUE, the animation result is a GIF video} \item{label_v}{vector
#' defining the vertex labels} \item{legend.position}{position of the legend}
#' \item{frame.color}{color of the frames} \item{label.hub}{logical ; if TRUE
#' only the hubs are labeled} \item{edge.arrow.size}{size of the arrows ;
#' default to 0.7} \item{edge.thickness}{edge thickness ; default to 1.} } }}}
#' 
#' \item{list("signature(x = \"omics_predict\", y = \"ANY\",...)")}{ \describe{
#' \item{x}{a omics_predict object} \item{list()}{Optional arguments: see plot
#' for omics_network} }} }
#' 
#' @param x a omics_array object, a omics_network object or a omics_predict object
#' @param y optional and not used if x is an appropriate structure
#' @param gr a vector giving the group of each gene 
#' @param choice what graphic should be plotted: either "F"
#' (for a representation of the matrices F) or "network".
#' @param nv the level of cutoff. Defaut to `0`.
#' @param ini using the ``position'' function, you can
#' fix the position of the nodes.
#' @param color.vertex a vector defining the color of the vertex.
#' @param ani animated plot?
#' @param size vector giving the size of the plot. Default to `c(2000,1000)`.
#' @param video if ani is TRUE and video is TRUE, the result of the animation is saved as an animated GIF. 
#' @param label_v vector defining the vertex labels.
#' @param legend.position position of the legend.
#' @param frame.color color of the frames.
#' @param label.hub logical ; if TRUE only the hubs are labeled.
#' @param edge.arrow.size size of the arrows ; default to 0.7.
#' @param edge.thickness edge thickness ; default to 1.
#' @param weight.node nodes weighting. Defaults to `NULL`.
#' @param horiz landscape? Defaults to `TRUE`.
#' @param time sets the time for plot of the prediction. Defaults to `NULL`
#' @param color.edge color of the edges
#' @param nround number of digits to display
#' @param ani.img.name name of image file  for animations
#' @param ani.imgdir name of the image directory for animations
#' @param ani.htmlfile name of the html file for animations
#' @param outdir name of the outdir for animations
#' @param ani.group.legend legend for animations
#' @param layout layout of the graphs
#' @param alpha transparency of the graphs
#' @param pixmap.color color for pixmap graphs
#' @param ... additional parameters
#' 
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords methods
#' @examples
#' 
#' if(require(CascadeData)){
#' data(micro_US, package="CascadeData")
#' micro_US<-as.omics_array(micro_US[1:100,],time=c(60,90,210,390),subject=6)
#' plot(micro_US)
#' }
#' 
setMethod("plot","omics_array",function(x,...)
{
  
  requireNamespace("lattice")
  #requireNamespace("grid")	
  xs<-t(x@omicsarray)
  
  rownames(xs)<-1:dim(xs)[1]
  ys<-x@time
  YS<-rep(ys,x@subject)
  suj<- rep(paste("Subject",1:x@subject,sep=" "),each=length(ys))
  ID<-paste("ID",1:dim(xs)[1],sep="")
  cclus<-unique(x@group)
  cclus<-cclus[order(cclus)]
  U<-data.frame(xs,suj,ys)
  V<- reshape(U,idvar="ID",varying=list(1:dim(xs)[2]), v.names = "conc", direction = "long")
  
  if(length(unique(x@group))>1){
    gr<-rep(paste("Cluster",x@group,sep=" "),each=x@subject*length(x@time))
    #if(!attr(dev.cur(),"names")=="pdf"){dev.new()}
    print(lattice::xyplot(V$conc~V$ys|V$suj,as.table=TRUE,xlab="Time",ylab="Gene Expression",group=rep(1:(x@subject*dim(xs)[2]),each=length(x@time)),type="l",scales=list(x=list(relation="free",at=x@time),y=list(relation="free")),col=rep(x@group,each=x@subject),key=list(
      space="right",
      lines=list(type="l",col=cclus),
      text=list(text=paste("Cluster",as.character(cclus)))
    )))
    #if(!attr(dev.cur(),"names")=="pdf"){dev.new()}
    print(lattice::xyplot(V$conc~V$ys|gr,as.table=TRUE,xlab="Time",ylab="Gene Expression",group=rep(1:(x@subject*dim(xs)[2]),each=length(x@time)),type="l",scales=list(x=list(relation="free",at=x@time),y=list(relation="free")),col=rep(1:x@subject,dim(xs)[2]),key=list(
      space="right",
      lines=list(type="l",col=1:x@subject),
      text=list(text=paste("Subject",as.character(1:x@subject)))
    )))
    for(i in 1:x@subject){
      ss<-V$suj
      sss<-ss==paste("Subject",i,sep=" ")
      #if(!attr(dev.cur(),"names")=="pdf"){dev.new()}
      print(lattice::xyplot(V$conc[sss]~V$ys[sss]|gr[sss],as.table=TRUE,xlab="Time",ylab="Gene Expression",type="l",main=paste("Subject",i,sep=" "),group=rep(1:(x@subject*dim(xs)[2]),each=length(x@time))[sss],scales=list(x=list(relation="free",at=x@time),y=list(relation="free"))))
    }
  }
  else{
    lattice::xyplot(V$conc~V$ys|V$suj,xlab="Time",as.table=TRUE,group=rep(1:(x@subject*dim(xs)[2]),each=length(x@time)),ylab="Gene Expression",scales=list(x=list(relation="free",at=x@time),y=list(relation="free")),type="l",col="black")
  }
}
)


#' Function to merge probesets
#' 
#' Used to collapse probesets using the collapseRows function of the WGCNA
#' package
#' 
#' 
#' @aliases probeMerge probeMerge,omics_array-method
#' @param x omicsarray
#' @param \dots Additionnal parameters to the collapseRows function of the
#' WGCNA package
#' @return Formal class 'omics_array' [package "Patterns"] with 7 slots
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords manip
#' @examples
#' 
#' if(require(CascadeData)){
#' data(micro_S)
#' D<-as.omics_array(micro_S[1:2000,],1:4,6)
#' D@gene_ID<-jetset::scores.hgu133plus2[D@name,"EntrezID"]
#' PM <- probeMerge(D)
#' }
#' 
setMethod("probeMerge","omics_array",function(x,...)
{
  #requireNamespace('WGCNA',quietly = TRUE)
  probeID<-x@name
  geneID<-x@gene_ID
  M<-x@omicsarray
  
  NAsPositions<-which(is.na(geneID))
  geneID[NAsPositions]<-paste("Unknown",1:length(NAsPositions))
  
  selec<-WGCNA::collapseRows(datET=M,rowGroup=geneID,rowID=probeID,...)
  
  M1<-new("omics_array"
          ,omicsarray=selec$datETcollapsed
          ,name=selec$group2row[,2]
          ,gene_ID=selec$group2row[,1]
          ,time=x@time
          ,subject=x@subject
          ,group=x@group[selec$selectedRow]
          ,start_time=x@start_time[selec$selectedRow])
  return(M1)
})



#' Cluster a omics_array object: determine optimal fuzzification parameter and
#' number of clusters.
#' 
#' Based on soft clustering performed by the Mfuzz package.
#' 
#' 
#' @aliases unsupervised_clustering_auto_m_c
#' unsupervised_clustering_auto_m_c,omics_array-method
#' @param M1 Object of omics_array class.
#' @param clust [NULL] Number of clusters.
#' @param mestim [NULL] Fuzzification parameter.
#' @param M2 [NULL] Object of omics_array class,
#' @param data_log [TRUE] Should data be logged?
#' @param screen [NULL] Specify `screen` parameter.
#' @param crange [NULL] Specify `crange` parameter.
#' @param repeats [NULL] Specify `repeats` parameter.
#' @param cselect [TRUE] Estimate `cselect` parameter.
#' @param dminimum [FALSE] Estimate `dminimum` parameter.
#' 
#' @return \item{m}{Estimate of the optimal fuzzification parameter.}
#' \item{c}{Estimate of the optimal number of clusters.} \item{csearch}{More
#' result from the cselection function of the Mfuzz package}
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords cluster
#' @examples
#' 
#' if(require(CascadeData)){
#' data(micro_S, package="CascadeData")
#' M<-as.omics_array(micro_S[1:100,],1:4,6)
#' mc<-unsupervised_clustering_auto_m_c(M)
#' }
#' 
setMethod(f="unsupervised_clustering_auto_m_c", 	
          signature=c("omics_array"),
          definition=function(M1,
                              clust=NULL,
                              mestim=NULL,
                              M2=NULL,
                              data_log=TRUE,
                              screen=NULL,
                              crange=NULL,
                              repeats=NULL,
                              cselect=TRUE,
                              dminimum=FALSE){
            require("Mfuzz", quietly = TRUE, warn.conflicts = FALSE)
            if(is.null(crange)){crange=seq(4,32,4)}
            if(is.null(repeats)){repeats=5}
            if(is.null(M2)){
              if(data_log==TRUE){
                M<-log(M1@omicsarray)
              }
              else{
                M<-(M1@omicsarray)
              }
            } else {
              if(data_log==TRUE){
                M1_mic<-log(M1@omicsarray)
                M2_mic<-log(M2@omicsarray)
              }
              else{
                M1_mic<-M1@omicsarray
                M2_mic<-M2@omicsarray
              }
              M=M1_mic-M2_mic
            } 
            colnames(M)<-NULL
            Em<-Biobase::ExpressionSet(M)
            Em.s <- Mfuzz::standardise(Em)
            if(is.null(mestim)){mestim <- Mfuzz::mestimate(Em.s)
            }
            if(cselect){
              if(is.null(clust)){
                tmp <- Mfuzz::cselection(Em.s,m=mestim,crange=crange,repeats=repeats,visu=TRUE)
                clust <- ifelse(which.max(rowSums(t(tmp)-crange)<0)==1,crange[dim(tmp)[2]],crange[which.max(rowSums(t(tmp)-crange)<0)-1])
              } else {tmp=NA}
            }
            if(dminimum){
              if(is.null(clust)){
                tmp <- Mfuzz::Dmin(Em.s,m=mestim,crange=crange,repeats=repeats,visu=TRUE)
                clust <- NA
              } else {tmp=NA}
            }
            return(list(m=mestim,c=clust,csearch=tmp))
          }
)



#' Cluster a omics_array object: performs the clustering.
#' 
#' Based on soft clustering performed by the Mfuzz package.
#' 
#' 
#' @aliases unsupervised_clustering
#' unsupervised_clustering,omics_array,numeric,numeric-method
#' @param M1 Object of omics_array class.
#' @param clust Number of clusters.
#' @param mestim Fuzzification parameter.
#' @param M2 [NULL] Object of omics_array class,
#' @param data_log [TRUE] Should data be logged?
#' @param screen [NULL] Specify `mfrow` parameter.
#' @param heatmap [TRUE] Plot heatmaps?
#' @param new.window [TRUE] Use new window?
#' 
#' @return An object of class omics_array with the group slot updated by groups
#' deduced from the soft clustering result.
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords cluster
#' @examples
#' 
#' if(require(CascadeData)){
#' data(micro_S, package="CascadeData")
#' M<-as.omics_array(micro_S[51:100,],1:4,6)
#' mc<-unsupervised_clustering_auto_m_c(M)
#' MwithGrp=unsupervised_clustering(M, 4, mc$m, screen=NULL, heatmap=FALSE, new.window = FALSE)
#' # Other options
#' unsupervised_clustering(M, 4, mc$m, screen=c(2,2), heatmap=TRUE, new.window = FALSE)
#' # Plot the clusters
#' plot(MwithGrp)
#' }
#' 
setMethod(f="unsupervised_clustering", 	
          signature=c("omics_array","numeric","numeric"),
          definition=function(M1,
                              clust,
                              mestim,
                              M2=NULL,
                              data_log=TRUE,
                              screen=NULL,
                              heatmap=TRUE,
                              new.window=TRUE){
            require("Mfuzz", quietly = TRUE, warn.conflicts = FALSE)
            if(is.null(M2)){
              if(data_log==TRUE){
                M<-log(M1@omicsarray)
              }
              else{
                M<-(M1@omicsarray)
              }
            } else {
              if(data_log==TRUE){
                M1_mic<-log(M1@omicsarray)
                M2_mic<-log(M2@omicsarray)
              }
              else{
                M1_mic<-(M1@omicsarray)
                M2_mic<-(M2@omicsarray)
              }
              M=M1_mic-M2_mic
            } 
            colnames(M)<-NULL
            Em<-Biobase::ExpressionSet(M)
            Em.s <- Mfuzz::standardise(Em)
            cl <- Mfuzz::mfuzz(Em.s,c=clust,m=mestim)
            Mfuzz::overlap.plot(cl,over =Mfuzz::overlap(cl), thres = 0.05)
            if(!attr(dev.cur(),"names")=="pdf"){
              if(is.null(screen)){Mfuzz::mfuzz.plot(Em.s,cl,new.window=new.window)
              } else {Mfuzz::mfuzz.plot(Em.s,cl,mfrow=screen,new.window=new.window)}} else {
                if(is.null(screen)){Mfuzz::mfuzz.plot(Em.s,cl,new.window=FALSE)
                } else {Mfuzz::mfuzz.plot(Em.s,cl,mfrow=screen,new.window=FALSE)}}  
            if(heatmap){
              require(gplots,quietly = TRUE)
              #if(!attr(dev.cur(),"names")=="pdf"){dev.new()}
              gplots::heatmap.2(Biobase::exprs(Em.s), dendrogram = 'row', Colv = FALSE, col = gplots::greenred(75), 
                        key = FALSE, keysize = 1.0, symkey = FALSE, density.info = 'none',
                        trace = 'none', colsep = rep(1:10), sepcolor = 'white', sepwidth = 0.05,
                        hclustfun = function (c){hclust(c, method = 'average')}, 
                        labRow = NA, cexCol = 1)
            }
            MM<-M1
            MM@group<-cl$cluster
            return(MM)
          }
)

#' Methods for selecting genes
#' 
#' Selection of differentially expressed genes.
#' 
#' 
#' @name geneSelection
#' @aliases genePeakSelection geneSelection genePeakSelection-methods
#' geneSelection-methods geneSelection,omics_array,numeric-method
#' geneSelection,omics_array,omics_array,numeric-method
#' geneSelection,list,list,numeric-method
#' genePeakSelection,omics_array,numeric-method
#' @param x either a omics_array object or a list of omics_array objects. In
#' the first case, the omics_array object represents the stimulated
#' measurements. In the second case, the control unstimulated data (if present)
#' should be the first element of the list.
#' @param y either a omics_array object or a list of strings. In the first
#' case, the omics_array object represents the stimulated measurements. In the
#' second case, the list is the way to specify the contrast: \describe{
#' \item{First element:}{ condition, condition&time or pattern. The condition
#' specification is used when the overall is to compare two conditions.  The
#' condition&time specification is used when comparing two conditions at two
#' precise time points. The pattern specification allows to decide which time
#' point should be differentially expressed.} \item{Second element:}{a vector
#' of length 2. The two conditions which should be compared. If a condition is
#' used as control, it should be the first element of the vector. However, if
#' this control is not measured throught time, the option cont=TRUE should be
#' used.} \item{Third element:}{depends on the first element.  It is no needed
#' if condition has been specified.  If condition&time has been specified, then
#' this is a vector containing the time point at which the comparison should be
#' done. If pattern has been specified, then this is a vector of 0 and 1 of
#' length T, where T is the number of time points. The time points with desired
#' differential expression are provided with 1.  }}
#' @param tot.number an integer. The number of selected genes. If tot.number <0
#' all differentially genes are selected. If tot.number > 1, tot.number is the
#' maximum of diffenrtially genes that will be selected.  If 0<tot.number<1,
#' tot.number represents the proportion of diffenrentially genes that are
#' selected.
#' @param peak interger. At which time points measurements should the genes be
#' selected [optionnal for geneSelection].
#' @param data_log logical (default to TRUE); should data be logged ?
#' @param wanted.patterns a matrix with wanted patterns [only for geneSelection].
#' @param forbidden.patterns a matrix with forbidden patterns [only for geneSelection].
#' @param durPeak vector of size 2 (default to c(1,1)) ; the first elements gives the length of the peak at
#' the left, the second at the right. [only for genePeakSelection]
#' @param abs_val logical (default to TRUE) ; should genes be selected on the
#' basis of their absolute value expression ? [only for genePeakSelection]
#' @param alpha_diff float; the risk level
#' @param alpha float; the risk level. Default to `alpha=0.05`
#' @param Design the design matrix of the experiment. Defaults to `NULL`.
#' @param lfc log fold change value used in limma's `topTable`. Defaults to 0.
#' @param cont use contrasts. Defaults to `FALSE`. 
#' @param f.asso function used to assess the association between the genes. 
#' The default value `NULL` implies the use of the usual `mean` function. 
#' @param return.diff [FALSE] if TRUE then the function returns the stimulated expression of the differentially expressed genes 
#' @return A omics_array object.
#' @author Frédéric Bertrand , Myriam Maumy-Bertrand.
#' @keywords methods
#' @examples
#' 
#' \donttest{
#'   if(require(CascadeData)){
#' 	data(micro_US)
#' 	micro_US<-as.omics_array(micro_US,time=c(60,90,210,390),subject=6)
#' 	data(micro_S)
#' 	micro_S<-as.omics_array(micro_S,time=c(60,90,210,390),subject=6)
#' 
#'   #Basically, to find the 50 more significant expressed genes you will use:
#'   Selection_1<-geneSelection(x=micro_S,y=micro_US,
#'   tot.number=50,data_log=TRUE)
#'   summary(Selection_1)
#'   
#'   #If we want to select genes that are differentially 
#'   #at time t60 or t90 :
#'   Selection_2<-geneSelection(x=micro_S,y=micro_US,tot.number=30,
#'   wanted.patterns=
#'   rbind(c(0,1,0,0),c(1,0,0,0),c(1,1,0,0)))
#'   summary(Selection_2)
#' 
#'   #To select genes that have a differential maximum of expression at a specific time point.
#'   
#'   Selection_3<-genePeakSelection(x=micro_S,y=micro_US,peak=1,
#'   abs_val=FALSE,alpha_diff=0.01)
#'   summary(Selection_3)
#'   }
#' 
#'   if(require(CascadeData)){
#' data(micro_US)
#' micro_US<-as.omics_array(micro_US,time=c(60,90,210,390),subject=6)
#' data(micro_S)
#' micro_S<-as.omics_array(micro_S,time=c(60,90,210,390),subject=6)
#' #Genes with differential expression at t1
#' Selection1<-geneSelection(x=micro_S,y=micro_US,20,wanted.patterns= rbind(c(1,0,0,0)))
#' #Genes with differential expression at t2
#' Selection2<-geneSelection(x=micro_S,y=micro_US,20,wanted.patterns= rbind(c(0,1,0,0)))
#' #Genes with differential expression at t3
#' Selection3<-geneSelection(x=micro_S,y=micro_US,20,wanted.patterns= rbind(c(0,0,1,0)))
#' #Genes with differential expression at t4
#' Selection4<-geneSelection(x=micro_S,y=micro_US,20,wanted.patterns= rbind(c(0,0,0,1)))
#' #Genes with global differential expression 
#' Selection5<-geneSelection(x=micro_S,y=micro_US,20)
#' 
#' #We then merge these selections:
#' Selection<-unionOmics(list(Selection1,Selection2,Selection3,Selection4,Selection5))
#' print(Selection)
#' 
#' #Prints the correlation graphics Figure 4:
#' summary(Selection,3)
#' 
#' ##Uncomment this code to retrieve geneids.
#' #library(org.Hs.eg.db)
#' #
#' #ff<-function(x){substr(x, 1, nchar(x)-3)}
#' #ff<-Vectorize(ff)
#' #
#' ##Here is the function to transform the probeset names to gene ID.
#' #
#' #library("hgu133plus2.db")
#' #
#' #probe_to_id<-function(n){  
#' #x <- hgu133plus2SYMBOL
#' #mp<-mappedkeys(x)
#' #xx <- unlist(as.list(x[mp]))
#' #genes_all = xx[(n)]
#' #genes_all[is.na(genes_all)]<-"unknown"
#' #return(genes_all)
#' #}
#' #Selection@name<-probe_to_id(Selection@name)
#'   }
#' 	}
#' 
setMethod(
  f = "geneSelection",
  signature = c("omics_array", "omics_array", "numeric"),
  definition = function(x,
                        y,
                        tot.number,
                        data_log = TRUE,
                        wanted.patterns = NULL,
                        forbidden.patterns = NULL,
                        peak = NULL,
                        alpha = 0.05,
                        Design = NULL,
                        lfc = 0) {
    warnings(
      "The use of the geneSelection with signature c(omics_array,omics_array) is depreciated",
      immediate. = TRUE
    )
    if (!is.null(sessionInfo()$otherPkgs$limma$Version)) {
      BBB <- strsplit(sessionInfo()$otherPkgs$limma$Version, "[.]")
    } else {
      if (!is.null(sessionInfo()$loadedOnly$limma$Version)) {
        BBB <- strsplit(sessionInfo()$loadedOnly$limma$Version, "[.]")
      } else {
        stop("Please install the limma BioConductor package")
      }
    }
    
    if (!(BBB[[1]][1] > 3 ||
          (BBB[[1]][1] == 3 && BBB[[1]][2] > 18) ||
          (BBB[[1]][1] == 3 &&
           BBB[[1]][2] == 18 && BBB[[1]][3] >= 13)))
    {
      stop("Upgrade your version of Limma (>= 3.18.13)")
    }
    
    
    indic <- 0
    M1 <- x
    M2 <- y
    if (is.null(M2)) {
      indic <- 1
      M2 <- M1
    }
    
    require(limma)
    
    if (data_log == TRUE) {
      M1_mic <- log(M1@omicsarray)
      M2_mic <- log(M2@omicsarray)
    } else{
      M1_mic <- (M1@omicsarray)
      M2_mic <- (M2@omicsarray)
    }
    
    if (is.null(rownames(M1_mic))) {
      rownames(M1_mic) <- paste("probe ", 1:dim(M1_mic)[1])
    }
    if (is.null(rownames(M2_mic))) {
      rownames(M2_mic) <- paste("probe ", 1:dim(M2_mic)[1])
    }
    
    if (indic == 1) {
      M2_mic <- M2_mic * 0
    }
    
    
    
    colnames(M1_mic) <-
      paste(rep("US", length(M1@time) * M1@subject),
            rep(M1@time, M1@subject),
            sep = "")
    colnames(M2_mic) <-
      paste(rep("S", length(M2@time) * M2@subject), rep(M2@time, M2@subject), sep =
              "")
    #rownames(M1_mic)<- paste("probe",rownames(M1_mic) )
    #rownames(M2_mic)<-  paste("probe",rownames(M2_mic) )
    M <- cbind(M1_mic, M2_mic)
    
    T <- length(M1@time)
    
    #Construction de la matrice de design
    
    if (is.null(Design)) {
      design <- t(matrix(rep(0, (
        M1@subject + M2@subject
      ) * T * (T * 2)), 2 * T))
      
      for (i in 1:(T)) {
        design[which(colnames(M) %in% paste("US", M1@time[i], sep = "")), i] <-
          1
        design[which(colnames(M) %in% paste("S", M2@time[i], sep =
                                              "")), i + T] <- 1
      }
      
      
      vnom <-
        paste(c(
          paste(rep("US_time", length(M1@time)), M1@time[1:(length(M1@time))], sep =
                  ""),
          paste(rep("S_time", length(M1@time)), M1@time[1:(length(M1@time))], sep =
                  "")
        ), sep = "")
      colnames(design) <- vnom
    } else{
      design <- Design$X
    }
    #block<-rep(1:(M1@subject*2),each=T)
    #dupcor <- duplicateCorrelation(M,design,block=block)
    #model<-lmFit(M,design,block=block)
    model <- lmFit(M, design)
    if (is.null(Design)) {
      diff <- paste(vnom[1:T], vnom[1:T + length(M1@time)], sep = "-")
      contr <-
        makeContrasts(contrasts = diff, levels = vnom)
    } else{
      contr <- Design$contr
    }
    model2 <- contrasts.fit(model, contr)
    model.test <- eBayes(model2)
    
    
    
    p.val.ind <- model.test$p.value
    rownames(p.val.ind) <- rownames(M1_mic)
    
    p.val.all <-
      topTable(model.test,
               p.value = alpha,
               number = dim(M)[1],
               lfc = lfc)
    kkkk <- 0
    if (is.null(p.val.all$ID)) {
      kkkk <- 1
      ID <- row.names(p.val.all)
      p.val.all <- cbind(ID, p.val.all)
      f_nom <- function(nom) {
        which(x@name %in% nom)
      }
      f_nom <- Vectorize(f_nom)
      row.names(p.val.all) <-
        unlist(f_nom(as.character(p.val.all$ID)))
    }
    if (tot.number > 0 && tot.number <= 1) {
      tot.number <- round(tot.number * dim(p.val.all)[1])
    }
    
    if (kkkk <- 0) {
      ind.class <- rownames(p.val.all)
      p.val.ind.class <- p.val.ind[(ind.class), ]
    } else{
      ind.class <- rownames(p.val.all)
      p.val.ind.class <- p.val.ind[as.numeric(ind.class), ]
      
    }
    f.p.val <- function(x) {
      if (x < alpha) {
        return(1)
      }
      else{
        return(0)
      }
    }
    
    p.val.ind.class.cat <-
      apply(p.val.ind.class, c(1, 2), f.p.val)
    choix <-
      cbind(p.val.ind.class.cat, rep(0, dim(p.val.ind.class.cat)[1]))
    ch <- dim(choix)[2]
    f_test_patt <- function(x, pat) {
      if (sum(abs(x - pat)) == 0) {
        return(1)
      }
      else{
        return(0)
      }
    }
    
    if (!is.null(forbidden.patterns)) {
      for (i in 1:dim(forbidden.patterns)[1]) {
        f_test2 <- function(x) {
          f_test_patt(x, forbidden.patterns[i, ])
        }
        S <- apply(choix[, 1:(ch - 1)], 1, f_test2)
        choix <- choix[S == 0, ]
      }
    }
    
    
    if (!is.null(wanted.patterns)) {
      sel <- rep(0, dim(choix)[1])
      choix[, ch] <- choix[, ch] * 0
      for (i in 1:dim(wanted.patterns)[1]) {
        f_test2 <- function(x) {
          f_test_patt(x, wanted.patterns[i, 1:T])
        }
        sel <- sel + apply(choix[, 1:(ch - 1)], 1, f_test2)
        
      }
      for (j in 1:length(sel)) {
        if ((sum(choix[1:(j - 1), ch]) < tot.number ||
             tot.number < 0) && sel[j] >= 1) {
          choix[j, ch] <- 1
        }
      }
    }
    
    
    f_test_peak <- function(x, peak) {
      x[peak]
    }
    
    if (!is.null(peak)) {
      f_test2 <- function(x) {
        f_test_peak(x, peak)
      }
      
      # if(tot.number>0){
      # S<-apply(choix[,1:(ch-1)],1,f_test2)
      # for(j in 1:length(S)){
      # if(sum(S[1:j])<=tot.number && S[j]==1){choix[j,ch]<-1}
      # }
      # }
      # else{
      
      choix[, ch] <- apply(choix[, 1:(ch - 1)], 1, f_test2)
      
      #}
      
    }
    
    
    if (tot.number > 0 && is.null(wanted.patterns)) {
      i <- 1
      PN <- dim(choix)[1]
      PN <- 1:PN
      
      if (is.null(peak)) {
        while (sum(choix[, ch]) < tot.number) {
          choix[i, ch] <- 1
          i <- i + 1
        }
        choix[PN >= i, ch] <- 0
        choix <- choix[choix[, ch] == 1, ]
      }
      else{
        choix <- choix[choix[, ch] == 1, ]
        PN <- dim(choix)[1]
        choix <- choix[1:min(PN, tot.number), ]
      }
      
    }
    
    if ((!is.null(wanted.patterns))) {
      choix <- choix[choix[, ch] == 1, ]
    }
    if (!is.null(peak) &&
        tot.number <= 0) {
      choix <- choix[choix[, ch] == 1, ]
    }
    
    temps.prem <- function(x) {
      u <- 0
      i <- 1
      while (u == 0) {
        if (x[i] == 1) {
          return(i)
        }
        else
          i = i + 1
      }
    }
    
    
    R <- apply(choix, 1, temps.prem)
    n <- length(R)
    MM1 <- M1_mic[rownames(choix), ] - M2_mic[rownames(choix), ]
    
    
    r3 <- function(choix) {
      sortie <- NULL
      for (i in 1:length(choix)) {
        sortie <- c(sortie, which(rownames(M1_mic) %in% rownames(choix)[i]))
      }
      sortie
    }
    
    
    
    M <-
      new(
        "omics_array",
        omicsarray = MM1,
        name = M1@name[r3(choix)],
        gene_ID = M1@gene_ID[r3(choix)],
        time = M1@time,
        subject = M1@subject,
        group = as.vector(R),
        start_time = as.vector(R)
      )
    
    
    
    return(M)
    
  }
  
)

#' @rdname geneSelection
setMethod(
  f = "geneSelection",
  signature = c("list", "list", "numeric"),
  definition = function(x,
                        y,
                        tot.number,
                        data_log = TRUE,
                        alpha = 0.05,
                        cont = FALSE,
                        lfc = 0,
                        f.asso = NULL,
                        return.diff=FALSE) {
    #here we check if the version of Limma is ok
    
    if (!is.null(sessionInfo()$otherPkgs$limma$Version)) {
      BBB <- strsplit(sessionInfo()$otherPkgs$limma$Version, "[.]")
    } else {
      if (!is.null(sessionInfo()$loadedOnly$limma$Version)) {
        BBB <- strsplit(sessionInfo()$loadedOnly$limma$Version, "[.]")
      } else {
        stop("Please install the limma BioConductor package")
      }
    }
    
    if (!(BBB[[1]][1] > 3 ||
          (BBB[[1]][1] == 3 && BBB[[1]][2] > 18) ||
          (BBB[[1]][1] == 3 &&
           BBB[[1]][2] == 18 && BBB[[1]][3] >= 13)))
    {
      stop("Upgrade your version of Limma (>= 3.18.13)")
    }
    
    
    M <- x
    contrast <- y
    
    M_mic <- M
    require(limma)
    n <- length(M)
    Time <- length(M[[2]]@time)
    
    if (cont == TRUE && is.null(f.asso)) {
      f.asso <- "mean"
      
    }
    
    Subj <- (M[[2]]@subject)
    
    
    if (data_log == TRUE) {
      for (i in 1:n) {
        M_mic[[i]] <- log(M[[i]]@omicsarray) / log(2)
      }
    } else{
      for (i in 1:n) {
        M_mic[[i]] <- (M[[i]]@omicsarray)
      }
    }
    
    
    
    
    #colN<-function(N,i){
    #			colnames(N)<-paste(rep("Time ",Time*Subj),rep(M[[1]]@time,Subj), sep="")
    #			return((N))
    #			}
    #	M_mic<-rapply(M_mic,colN,how="list")
    #
    #
    #
    Ma <- abind::abind(M_mic, along = 2)
    T <- Time
    
    
    #Construction de la matrice de design
    
    
    condition <-
      as.factor(paste("condition", rep(1:n, unlist(lapply(
        M, ncol
      ))), sep = ""))
    
    if (cont == FALSE) {
      timeT <-
        paste("time", rep(1:T, sum(unlist(
          lapply(M, function(x) {
            x@subject
          })
        ))), sep = "")
    } else{
      timeT <-
        c(rep("time0", ncol(M[[1]])), paste("time", rep(1:T, sum(
          unlist(lapply(M[2:length(M)], function(x) {
            x@subject
          }))
        )), sep = ""))
      
    }
    
    Fac <- as.factor(paste(condition, timeT, sep = "."))
    
    if (contrast[[1]] != "condition") {
      formule <-  ~ -1 + Fac
      design <- model.matrix(formule)
      colnames(design) <- levels(Fac)
    } else{
      formule <-  ~ -1 + condition
      design <- model.matrix(formule)
      colnames(design) <- levels(condition)
    }
    
    model <- lmFit(Ma, design)
    
    
    if (contrast[[1]] == "condition") {
      contrastM <-
        makeContrasts(
          contrasts = paste("condition", contrast[[2]][1], "-condition",
                            contrast[[2]][2],
                            sep = ""),
          levels = design
        )
    }
    #if(contrast[[1]]=="time"){
    #
    #			if(contrast[[2]][1]!=1){
    #				contrastM[contrast[[2]][1]+n]<-contrastM[contrast[[2]][1]+n]+1
    #			}else{
    #				contrastM[1:n]<-contrastM[1:n]+1
    #				contrastM[(n+1):(n+T-1)]<-contrastM[(n+1):(n+T-1)]-1
    #				}
    #
    #			if(contrast[[2]][2]!=1){
    #				contrastM[contrast[[2]][2]+n]<-contrastM[contrast[[2]][2]+n]-1
    #			}else{
    #				contrastM[1:n]<-contrastM[1:n]-1
    #				contrastM[(n+1):(n+T-1)]<-contrastM[(n+1):(n+T-1)]+1
    #				}
    #
    #
    #		}
    #
    
    if (contrast[[1]] == "patterns" &&
        (cont == FALSE || contrast[[2]][1] != 1)) {
      coma = NULL
      for (j in 1:Time) {
        coma <-
          c(
            coma,
            paste(
              "condition",
              contrast[[2]][1],
              ".time",
              j,
              "-condition",
              contrast[[2]][2],
              ".time",
              j,
              sep = ""
            )
          )
      }
      contrastM <- makeContrasts(contrasts = coma, levels = design)
      
    }
    
    if (contrast[[1]] == "patterns" &&
        (cont == TRUE  && contrast[[2]][1] == 1)) {
      coma = NULL
      for (j in 1:Time) {
        coma <-
          c(
            coma,
            paste(
              "condition",
              contrast[[2]][1],
              ".time",
              0,
              "-condition",
              contrast[[2]][2],
              ".time",
              j,
              sep = ""
            )
          )
      }
      contrastM <- makeContrasts(contrasts = coma, levels = design)
      
    }
    
    if (contrast[[1]] == "condition&time") {
      coma <-
        paste(
          "condition",
          contrast[[2]][1],
          ".time",
          contrast[[3]][1],
          "-condition",
          contrast[[2]][2],
          ".time",
          contrast[[3]][2],
          sep = ""
        )
      contrastM <- makeContrasts(contrasts = coma, levels = design)
      
    }
    
    
    
    model2 <- contrasts.fit(model, -contrastM)
    model.test <- eBayes(model2)
    
    p.val.ind <- model.test$p.value
    
    
    
    
    nb.tot <- length(M[[1]]@name)
    
    
    
    if (contrast[[1]] == "patterns") {
      tableT <- matrix(0, nb.tot, Time + 1)
      row.names(tableT) <- 1:nb.tot
      tableT[, Time + 1] <- model.test$F.p.value
      for (j in 1:Time) {
        p.val.all <-
          topTable(
            model.test,
            p.value = alpha,
            number = nb.tot,
            lfc = lfc,
            coef = j
          )
        
        if (is.null(p.val.all$ID)) {
          ID <- row.names(p.val.all)
          p.val.all <- cbind(ID, p.val.all)
          f_nom <- function(nom) {
            which(x[[1]]@name %in% nom)
          }
          f_nom <- Vectorize(f_nom)
          row.names(p.val.all) <- f_nom(p.val.all$ID)
        }
        tableT[as.numeric(row.names(p.val.all)), j] <- 1
        
      }
      
      f.test <- function(x1) {
        rep <- FALSE
        for (i in 1:nrow(contrast[[3]])) {
          if (sum(x1 == contrast[[3]][i, ]) == Time) {
            rep <- TRUE
          }
        }
        return(rep)
      }
      
      B <- apply(tableT[, 1:Time], 1, f.test)
      tableT <- tableT[B, ]
      tableT <- tableT[which(tableT[, Time + 1] < alpha), ]
      tableT <- tableT[order(tableT[, Time + 1]), ]
      nb.ret <- nrow(tableT)
      if (tot.number >= 1 && nb.ret > tot.number) {
        nb.ret <- tot.number
      }
      
      if (tot.number < 1 && tot.number > 0) {
        nb.ret <- round(nb.ret * tot.number)
        
      }
      K1 <-
        M_mic[[contrast[[2]][1]]][as.numeric(row.names(tableT[1:nb.ret, ])), ]
      K2 <-
        M_mic[[contrast[[2]][2]]][as.numeric(row.names(tableT[1:nb.ret, ])), ]
      
      
      
    } else{
      p.val.all <- topTable(model.test,
                            p.value = alpha,
                            number = nb.tot,
                            lfc = lfc)
      if (dim(p.val.all)[1] == 0) {
        print("The selection is empty")
      } else{
        print("The selection is not empty")
      }
      
      if (is.null(p.val.all$ID)) {
        ID <- row.names(p.val.all)
        p.val.all <- cbind(ID, p.val.all)
        f_nom <- function(nom) {
          which(x[[1]]@name %in% nom)
        }
        f_nom <- Vectorize(f_nom)
        row.names(p.val.all) <- f_nom(p.val.all$ID)
      }
      nb.ret <- dim(p.val.all)[1]
      
      if (tot.number >= 1 && nb.ret > tot.number) {
        nb.ret <- tot.number
      }
      
      if (tot.number < 1 && tot.number > 0) {
        nb.ret <- round(nb.ret * tot.number)
        
      }
      
      p.val.all <- p.val.all[1:nb.ret, ]
      
      
      K1 <-
        M_mic[[contrast[[2]][1]]][as.character(p.val.all$ID), ]
      K2 <-
        M_mic[[contrast[[2]][2]]][as.character(p.val.all$ID), ]
    }
    
    
    if (!is.null(f.asso) && cont == TRUE) {
      K1 <- apply(K1, 1, f.asso)
    }
    if (!is.null(f.asso) && cont == FALSE) {
      for (i in 1:Time) {
        K1[i, ] <- apply(K1[seq(i, ncol(K1), by = Time), ], 1, f.asso)
      }
    }
    if (sum(dim(K1) == dim(K2)) == 2) {
      if(return.diff){
        print(
          "This function returns the stimulated expression of the differentially expressed genes"
        )
        MM1 <- K2 - K1 #-K1
      } else {
        print(
          "This function returns the stimulated expression"
        )
        MM1 <- K2
      }
    } else{
      warning("The number of patient is not equal. This function returns the stimulated expression")
      MM1 <- K2
    }
    
    if (contrast[[1]] == "patterns") {
      f_gr <- function(x) {
        min(which(x == 1))
      }
      
      group = apply(tableT[1:nb.ret, 1:Time], 1, f_gr)
      
      
    } else{
      group = rep(0, nb.ret)
    }
    
    head(x[[1]])
    if (length(x[[1]]@gene_ID > 2) & length(x[[1]]@name > 2)) {
      azert <- cbind(x[[1]]@name, x[[1]]@gene_ID)
      row.names(azert) <- x[[1]]@name
      head(azert)
      gene_ID <- azert[row.names(MM1), 2]
    } else{
      gene_ID <- 0
    }
    M <-
      new(
        "omics_array",
        omicsarray = MM1,
        name = row.names(MM1),
        gene_ID = gene_ID,
        time = M[[contrast[[2]][2]]]@time,
        subject = Subj,
        group = group,
        start_time = group
      )
    
    
    return(M)
    
  }
  
)


#' @rdname geneSelection
setMethod(f="genePeakSelection", 
          signature=c("omics_array","numeric"),
          definition=function(x,peak,y=NULL,data_log=TRUE,durPeak=c(1,1),abs_val=TRUE,alpha_diff=0.05){
            
            M1<-x
            M2<-y
            Select<-geneSelection(M1,tot.number=-1,M2,data_log=data_log,peak=peak)
            
            if(abs_val==FALSE){
              M<-Select@omicsarray
              sel<-rep(0,length(M[,1]))
              comp<-M[,(1:M1@subject-1)*length(M1@time)+peak]
              test1<-function(y){
                if(wilcox.test(y[1:M1@subject],y[(M1@subject+1):(2*M1@subject)], alternative="less")$p.value<alpha_diff){
                  return(1)
                }
                else{
                  return(0)
                }
              }
              uu<-0
              if((peak-durPeak[1])>0){
                for(i in 1:(peak-durPeak[1])){
                  comp1<-M[,(1:M1@subject-1)*length(M1@time)+i]
                  comp2<-cbind(comp1,comp)
                  sel<-sel+apply(comp2,1,test1)
                  uu<-uu+1				
                }
              }
              
              if((peak+durPeak[2])<=length(M1@time)){
                for(i in (peak+durPeak[2]):length(M1@time)){
                  comp1<-M[,(1:M1@subject-1)*length(M1@time)+i]
                  comp2<-cbind(comp1,comp)
                  sel<-sel+apply(comp2,1,test1)	
                  uu<-uu+1			
                }
              }
              
              
              
              #print(uu)
              N1<-Select@name[sel==uu]
              
              M<--Select@omicsarray
              sel<-rep(0,length(M[,1]))
              comp<-M[,(1:M1@subject-1)*length(M1@time)+peak]
              test1<-function(y){
                if(wilcox.test(y[1:M1@subject],y[(M1@subject+1):(2*M1@subject)], alternative="less")$p.value<alpha_diff){
                  return(1)
                }
                else{
                  return(0)
                }
              }
              uu<-0
              if((peak-durPeak[1])>0){
                for(i in 1:(peak-durPeak[1])){
                  comp1<-M[,(1:M1@subject-1)*length(M1@time)+i]
                  comp2<-cbind(comp1,comp)
                  sel<-sel+apply(comp2,1,test1)
                  uu<-uu+1				
                }
              }
              
              if((peak+durPeak[2])<=length(M1@time)){
                for(i in (peak+durPeak[2]):length(M1@time)){
                  comp1<-M[,(1:M1@subject-1)*length(M1@time)+i]
                  comp2<-cbind(comp1,comp)
                  sel<-sel+apply(comp2,1,test1)	
                  uu<-uu+1			
                }
              }
              
              N11<-Select@name[sel==uu]
              #Select2<-geneSelection(M1,M2,tot.number=-1,data_log=data_log,wanted.patterns=t(as.matrix(c(0,1,0,0,1000))))
              #N2<-Select2@name
              
              N<-c(N1,N11)
              Mi<-new("omics_array",omicsarray=Select@omicsarray[which(Select@name %in% N ),],name=N,gene_ID=Select@gene_ID[which(Select@name %in% N )],time=M1@time,subject=M1@subject,start_time=Select@start_time[which(Select@name %in% N)],group=rep(peak,length(N)))			
              return(Mi)
            } 		
            else{
              M<-abs(Select@omicsarray)
              sel<-rep(0,length(M[,1]))
              comp<-M[,(1:M1@subject-1)*length(M1@time)+peak]
              test1<-function(y){
                if(wilcox.test(y[1:M1@subject],y[(M1@subject+1):(2*M1@subject)], alternative="less")$p.value<0.05){
                  return(1)
                }
                else{
                  return(0)
                }
              }
              uu<-0
              if((peak-durPeak[1])>0){
                for(i in 1:(peak-durPeak[1])){
                  comp1<-M[,(1:M1@subject-1)*length(M1@time)+i]
                  comp2<-cbind(comp1,comp)
                  sel<-sel+apply(comp2,1,test1)
                  uu<-uu+1				
                }
              }
              
              if((peak+durPeak[2])<=length(M1@time)){
                for(i in (peak+durPeak[2]):length(M1@time)){
                  comp1<-M[,(1:M1@subject-1)*length(M1@time)+i]
                  comp2<-cbind(comp1,comp)
                  sel<-sel+apply(comp2,1,test1)	
                  uu<-uu+1			
                }
              }
              
              N<-Select@name[sel==uu]
              #Select2<-geneSelection(M1,M2,tot.number=-1,data_log=data_log,wanted.patterns=t(as.matrix(c(0,1,0,0,1000))))
              #N2<-Select2@name
              
              
              Mi<-new("omics_array",omicsarray=Select@omicsarray[which(Select@name %in% N ),],name=N,gene_ID=Select@gene_ID[which(Select@name %in% N )],time=M1@time,subject=M1@subject,group=rep(peak,length(N)),start_time=Select@start_time[which(Select@name %in% N)])	
              
              
            }
            
            
            
          }
)


#' Makes the union between two omics_array objects.
#' 
#' Makes the union between two omics_array objects.
#' 
#' 
#' @name unionOmics-methods
#' @aliases unionOmics unionOmics-methods
#' unionOmics,omics_array,omics_array-method unionOmics,list,ANY-method
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(M1 = \"omics_array\", M2 = \"omics_array\")")}{
#' Returns a omics_array object which is the union of M1 and M2.  }
#' 
#' \item{list("signature(M1 = \"list\", M2 = \"ANY\")")}{ Returns a omics_array
#' object which is the union of the elements of M1.  } }
#' 
#' @param M1 a omics-array or a list of omics-arrays
#' @param M2 a omics-array or nothing if M1 is a list of omics-arrays
#' 
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords methods
#' @examples
#' 
#' if(require(CascadeData)){
#' data(micro_S, package="CascadeData")
#' #Create another omicsarray object with 100 genes
#' Mbis<-M<-as.omics_array(micro_S[1:100,],1:4,6)
#' #Rename the 100 genes
#' Mbis@name<-paste(M@name,"bis")
#' rownames(Mbis@omicsarray) <- Mbis@name
#' #Union (merge without duplicated names) of the two omicsarrays. 
#' str(unionOmics(M,Mbis))
#' }
#' 
setMethod(f="unionOmics", 
          signature=c("omics_array","omics_array"),
          definition=function(M1,M2){
            
            
            
            nom1<-rownames(M1@omicsarray)
            nom2<-rownames(M2@omicsarray)
            corres<-cbind(c(M1@name,M2@name),c(nom1,nom2))
            corres<-unique(unique(corres,MARGIN=2))
            corres2<-cbind(c(M1@gene_ID,M2@gene_ID),c(nom1,nom2))
            corres2<-unique(unique(corres2,MARGIN=2))
            NOM<-unique(c(nom1,nom2))
            n<-length(NOM)
            m1<-M1@omicsarray[which(nom1 %in% NOM),]
            NOM2<-NOM[-which(nom1 %in% NOM)]
            m2<-M2@omicsarray[which(nom2 %in% NOM2),]		
            M<-rbind(m1,m2)
            
            gr1<-M1@group[which(nom1 %in% NOM)]
            gr2<-M2@group[which(nom2 %in% NOM2)]
            gr<-c(gr1,gr2)
            str1<-M1@start_time[which(nom1 %in% NOM)]
            str2<-M2@start_time[which(nom2 %in% NOM2)]
            str<-c(str1,str2)
            rep<-new("omics_array",omicsarray=M,name=corres[,1],gene_ID=corres2[,1],time=M1@time,subject=M1@subject,group=gr,start_time=str)
            return(rep)
          }
          
) 

setMethod(f="unionOmics", 
          signature=c("list","ANY"),
          definition=function(M1,M2){
            
            rep<-unionOmics(M1[[1]],M1[[1]])
            if(length(M1)>1){	
              for(i in 2:length(M1)){
                rep<-unionOmics(rep,M1[[i]])
                
              }
            }
            return(rep)
          }
          
)



#' A function to explore a dataset and cluster its rows.
#' 
#' Based on soft clustering performed by the Mfuzz package.
#' 
#' 
#' @aliases clustExploration clustExploration-methods
#' clustExploration,omics_array-method
#' @param omicsarray A omicsarray to cluster
#' @param new.window Boolean. New X11 window for plots. Defaults to FALSE.
#' @return A data.frame of nrows(omicsarray) observations of 3 variables (name,
#' cluster, maj.vote.index).
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords cluster
#' @examples
#' 
#' library(Patterns)
#' if(require(CascadeData)){
#' data(micro_S, package="CascadeData")
#' D<-Patterns::as.omics_array(micro_S[1:100,],1:4,6)
#' a<-clustExploration(D)
#' a
#' }
#' 
setMethod(f="clustExploration", 
          signature=c("omics_array"),
          definition=function(omicsarray,new.window=FALSE){
  
  require("Mfuzz", quietly = TRUE, warn.conflicts = FALSE)
  T<-length(omicsarray@time)
  P<-omicsarray@subject
  M1<-omicsarray@omicsarray
  
  M<-M1[,1:T]
  for(i in 2:P ){
    M<-rbind(M,M1[,1:T+(i-1)*T])
  }
  
  M<-t(scale(t(M),scale=FALSE))
  
  M<-M/sqrt(diag((M)%*%t(M)))
  
  lambda_max<-max(eigen(M%*%t(M)/dim(M)[1])$values)
  
  
  #for the choice of m, I recommend the lecture of :
  # http://ccc.inaoep.mx/~ariel/2012/Analysis%20of%20parameter%20selections%20for%20fuzzy%20c-means.pdf
  
  
  
  if( lambda_max <0.5){
    m_max<-min(1/(1-2*lambda_max),2.5)
  }else{
    m_max<-2.5
  }
  rownames(M)<-NULL
  M_exp<- Biobase::ExpressionSet(M)
  
  indic<-"continue"
  k<-1
  
  while(indic=="continue"){
    k<-k+1
    cl <- Mfuzz::mfuzz(M_exp, c = k,m=m_max)
    M_clust<-matrix(cl$cluster,nrow(M1),P)
    mv<-apply(M_clust,1,majority_vote)  
    mind<-apply(M_clust,1,majority_indice)  
    cont<-mean(mind/P)
    if(cont<0.4 & k>2){indic="stop"}else{
      aM_clust<-M_clust
      amv<-mv
      amind<-mind
      acont<-cont
    }
    
  }
  k<-k-1
  cl <- Mfuzz::mfuzz(M_exp, c = k,m=m_max)
  if(k<5){
    Mfuzz::mfuzz.plot(M_exp,cl=cl,mfrow=c(1,k), new.window = new.window)
  }else{
    Mfuzz::mfuzz.plot(M_exp,cl=cl,mfrow=c(ceiling(sqrt(k)),ceiling(sqrt(k))), new.window = new.window)
  }
  
  
  return(data.frame(name=omicsarray@name,cluster=amv,maj.vote.index=amind))
  
  
}  
)



#' A function to explore a dataset and cluster its rows.
#' 
#' Based on soft clustering performed by the Mfuzz package.
#' 
#' 
#' @aliases clustInference clustInference-methods
#' clustInference,omics_array,numeric-method
#' @param omicsarray A omicsarray to cluster
#' @param vote.index Option for cluster attribution
#' @param new.window Boolean. New X11 window for plots. Defaults to FALSE.
#' @return A list of two elements: \item{res.matrix}{A data.frame of
#' nrows(omicsarray) observations of 3 variables (name, cluster,
#' maj.vote.index).} \item{prop.matrix}{Additionnal info.}
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords cluster
#' @examples
#' 
#' library(Patterns)
#' if(require(CascadeData)){
#' data(micro_S, package="CascadeData")
#' D<-Patterns::as.omics_array(micro_S[1:20,],1:4,6)
#' b<-Patterns::clustInference(D,0.5)
#' b
#' }
#' 
setMethod(f="clustInference", 
          signature=c("omics_array","numeric"),
          definition=function(omicsarray,vote.index,new.window=FALSE){  
  require("Mfuzz", quietly = TRUE, warn.conflicts = FALSE)
            
  T<-length(omicsarray@time)
  P<-omicsarray@subject
  M1<-omicsarray@omicsarray
  
  M<-M1[,1:T]
  for(i in 2:P ){
    M<-rbind(M,M1[,1:T+(i-1)*T])
  }
  
  M<-t(scale(t(M),scale=FALSE))
  
  M<-M/sqrt(diag((M)%*%t(M)))
  
  lambda_max<-max(eigen(M%*%t(M)/dim(M)[1])$values)
  
  
  
  if( lambda_max <0.5){
    m_max<-min(1/(1-2*lambda_max),2.5)
  }else{
    m_max<-2.5
  }
  rownames(M)<-NULL
  M_exp<- Biobase::ExpressionSet(M)
  
  indic<-"continue"
  k<-1
  within_error_hard<-NULL
  while(indic=="continue"){
    k<-k+1
    cl<-Mfuzz::mfuzz(M_exp, c = k,m=m_max)
    we<-cl$withinerror
    for(f in 1:50){
      CL <- Mfuzz::mfuzz(M_exp, c = k,m=m_max)
      wee<-CL$withinerror
      if(wee<we){
        cl<-CL
        we<-wee
      }
      
    }
    
    M_clust<-matrix(cl$cluster,nrow(M1),P)
    mv<-apply(M_clust,1,majority_vote)  
    mind<-apply(M_clust,1,majority_indice)  
    cont<-mean(mind/P)
    within_error_hard<-c(within_error_hard,mean((M-cl$centers[cl$cluster,])^2))
    
    if((cont<vote.index|| min(table(mv))<(ceiling(((T)^2/P))+1)) & k>2){indic="stop"}else{
      aM_clust<-M_clust
      amv<-mv
      amind<-mind
      acont<-cont
      acl<-cl
    }
    
  }
  
  prop.clust<-rep(0,nrow(M))
  Amv<-rep(amv,P)
  
  for(kk in 1:nrow(M)){
    
    prop.clust[kk]<-acl$member[kk,Amv[kk]]
    
    
  }
  
  
  within_error_hard<-within_error_hard[-length(within_error_hard)]
  prop.matrix<-matrix(prop.clust,nrow(M1),P)
  k<-k-1
  cl <- Mfuzz::mfuzz(M_exp, c = k,m=m_max)
  plot(2:(length(within_error_hard)+1),within_error_hard)
  
  
  
  
  if(k<5){
    Mfuzz::mfuzz.plot(M_exp,cl=cl,mfrow=c(1,k), new.window = new.window)
  }else{
    Mfuzz::mfuzz.plot(M_exp,cl=cl,mfrow=c(ceiling(sqrt(k)),ceiling(sqrt(k))), new.window = new.window)
  }
  
  
  return(res.matrix=list(data.frame(name=omicsarray@name,cluster=amv,maj.vote.index=amind),prop.matrix=prop.matrix))
  
  
}  
)
