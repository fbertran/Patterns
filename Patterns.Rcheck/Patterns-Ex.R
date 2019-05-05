pkgname <- "Patterns"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "Patterns-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('Patterns')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("CLL")
### * CLL

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CLL
### Title: Expression data from healthy and malignant (chronic lymphocytic
###   leukemia, CLL) human B-lymphocytes after B-cell receptor stimulation
###   (GSE 39411 dataset)
### Aliases: CLL
### Keywords: datasets

### ** Examples

data(CLL)
str(CLL)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CLL", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CascadeFinit")
### * CascadeFinit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CascadeFinit
### Title: Create initial F matrices for cascade networks inference.
### Aliases: CascadeFinit
### Keywords: models

### ** Examples

CascadeFinit(3,2)
CascadeFinit(4,3)
plotF(CascadeFinit(4,3),choice = "F")
CascadeFinit(3,2,low.trig=FALSE)
CascadeFinit(4,3,low.trig=FALSE)
plotF(CascadeFinit(4,3,low.trig=FALSE),choice = "F")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CascadeFinit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CascadeFshape")
### * CascadeFshape

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CascadeFshape
### Title: Create F matrices shaped for cascade networks inference.
### Aliases: CascadeFshape
### Keywords: models

### ** Examples

CascadeFshape(3,2)
plotF(CascadeFshape(3,2),choice = "Fshape")
CascadeFshape(4,3)
plotF(CascadeFshape(4,3),choice = "Fshape")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CascadeFshape", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("IndicFinit")
### * IndicFinit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: IndicFinit
### Title: Create initial F matrices using specific intergroup actions for
###   network inference.
### Aliases: IndicFinit
### Keywords: models

### ** Examples

IndicFinit(3, 2, matrix(1,2,2)-diag(2))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("IndicFinit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("IndicFshape")
### * IndicFshape

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: IndicFshape
### Title: Create F matrices using specific intergroup actions for network
###   inference.
### Aliases: IndicFshape
### Keywords: models

### ** Examples

IndicFshape(3, 2, matrix(1,2,2)-diag(2))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("IndicFshape", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Selection")
### * Selection

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Selection
### Title: Selection of genes.
### Aliases: Selection
### Keywords: datasets

### ** Examples

data(Selection)
head(Selection)
summary(Selection,3)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Selection", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("analyze_network-methods")
### * analyze_network-methods

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: analyze_network
### Title: Analysing the network
### Aliases: analyze_network analyze_network-methods
###   analyze_network,network-method
### Keywords: methods

### ** Examples

data(network)
analyze_network(network,nv=0)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("analyze_network-methods", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("as.micro_array")
### * as.micro_array

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: as.micro_array
### Title: Coerce a matrix into a micro_array object.
### Aliases: as.micro_array

### ** Examples

if(require(CascadeData)){
	data(micro_US, package="CascadeData")
	micro_US<-as.micro_array(micro_US[1:100,],time=c(60,90,210,390),subject=6)
	plot(micro_US)
	}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("as.micro_array", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("clustExploration-methods")
### * clustExploration-methods

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: clustExploration
### Title: A function to explore a dataset and cluster its rows.
### Aliases: clustExploration clustExploration-methods
###   clustExploration,micro_array-method
### Keywords: cluster

### ** Examples

library(Patterns)
if(require(CascadeData)){
data(micro_S, package="CascadeData")
D<-Patterns::as.micro_array(micro_S[1:100,],1:4,6)
a<-clustExploration(D)
a
}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("clustExploration-methods", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("clustInference-methods")
### * clustInference-methods

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: clustInference
### Title: A function to explore a dataset and cluster its rows.
### Aliases: clustInference clustInference-methods
###   clustInference,micro_array,numeric-method
### Keywords: cluster

### ** Examples

library(Patterns)
if(require(CascadeData)){
data(micro_S, package="CascadeData")
D<-Patterns::as.micro_array(micro_S[1:100,],1:4,6)
b<-clustInference(D,0.5)
b
}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("clustInference-methods", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compare-methods")
### * compare-methods

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compare-methods
### Title: Some basic criteria of comparison between actual and inferred
###   network.
### Aliases: compare-methods compare compare,network,network,numeric-method

### ** Examples

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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compare-methods", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cutoff-methods")
### * cutoff-methods

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cutoff
### Title: Choose the best cutoff
### Aliases: cutoff cutoff-methods cutoff,network-method
### Keywords: methods

### ** Examples




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cutoff-methods", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("evolution-methods")
### * evolution-methods

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: evolution
### Title: See the evolution of the network with change of cutoff
### Aliases: evolution evolution-methods evolution,network-method
### Keywords: methods

### ** Examples




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("evolution-methods", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("geneNeighborhood-methods")
### * geneNeighborhood-methods

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: geneNeighborhood
### Title: Find the neighborhood of a set of nodes.
### Aliases: geneNeighborhood geneNeighborhood-methods
###   geneNeighborhood,network-method
### Keywords: methods

### ** Examples

data(Selection)
data(infos)
#Find probesets for EGR1
pbst_EGR1 = infos[infos$hgnc_symbol=="EGR1", "affy_hg_u133_plus_2"]

gene_IDs = infos[match(Selection@name, infos$affy_hg_u133_plus_), "hgnc_symbol"]

data(network)
#A nv value can chosen using the cutoff function
nv=.11 
EGR1<-which(is.element(Selection@name,pbst_EGR1))
P<-position(network,nv=nv)

geneNeighborhood(network,targets=EGR1,nv=nv,ini=P,
label_v=gene_IDs)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("geneNeighborhood-methods", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("genePeakSelection-methods")
### * genePeakSelection-methods

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: genePeakSelection
### Title: Methods for selecting genes
### Aliases: genePeakSelection geneSelection genePeakSelection-methods
###   geneSelection-methods geneSelection,micro_array,numeric-method
###   geneSelection,micro_array,micro_array,numeric-method
###   geneSelection,list,list,numeric-method
###   genePeakSelection,micro_array,numeric-method
### Keywords: methods

### ** Examples




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("genePeakSelection-methods", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gene_expr_simulation")
### * gene_expr_simulation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gene_expr_simulation
### Title: Simulates microarray data based on a given network.
### Aliases: gene_expr_simulation gene_expr_simulation-methods
###   gene_expr_simulation,network-method

### ** Examples

data(simul)
set.seed(1)

#We simulate gene expressions according to the network Net
Msim<-Patterns::gene_expr_simulation(
	network=Net,
	time_label=rep(1:4,each=25),
	subject=5,
	peak_level=200)
head(Msim)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gene_expr_simulation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("head-methods")
### * head-methods

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: head
### Title: Overview of a micro_array object
### Aliases: methods head-methods head,ANY-method head,micro_array-method
### Keywords: methods

### ** Examples

  if(require(CascadeData)){
	data(micro_US)
	micro_US<-as.micro_array(micro_US,time=c(60,90,210,390),subject=6)
	head(micro_US)
	}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("head-methods", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("inference-methods")
### * inference-methods

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: inference
### Title: Reverse-engineer the network
### Aliases: inference inference-methods inference,micro_array-method
### Keywords: methods

### ** Examples




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("inference-methods", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("infos")
### * infos

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: infos
### Title: Details on some probesets of the affy_hg_u133_plus_2 platform.
### Aliases: infos
### Keywords: datasets

### ** Examples

data(infos)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("infos", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("micro_array-class")
### * micro_array-class

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: micro_array-class
### Title: Class '"micro_array"'
### Aliases: micro_array-class
### Keywords: classes

### ** Examples

showClass("micro_array")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("micro_array-class", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("micropred-class")
### * micropred-class

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: micropredict-class
### Title: Class '"micropred"'
### Aliases: micropredict-class
### Keywords: classes

### ** Examples

showClass("network")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("micropred-class", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("network-class")
### * network-class

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: network-class
### Title: Class '"network"'
### Aliases: network-class
### Keywords: classes

### ** Examples

showClass("network")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("network-class", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("network")
### * network

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: network
### Title: A example of an inferred network (4 groups case).
### Aliases: network
### Keywords: datasets

### ** Examples

data(network)
str(network)
plot(network)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("network", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("network2gp")
### * network2gp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: network2gp
### Title: A example of an inferred cascade network (2 groups case).
### Aliases: network2gp
### Keywords: datasets

### ** Examples

data(network2gp)
str(network2gp)
plot(network2gp)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("network2gp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("networkCascade")
### * networkCascade

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: networkCascade
### Title: A example of an inferred cascade network (4 groups case).
### Aliases: networkCascade
### Keywords: datasets

### ** Examples

data(networkCascade)
str(networkCascade)
plot(networkCascade)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("networkCascade", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("network_random")
### * network_random

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: network_random
### Title: Generates a network.
### Aliases: network_random

### ** Examples

set.seed(1)
Net<-network_random(
	nb=100,
	time_label=rep(1:4,each=25),
	exp=1,
	init=1,
	regul=round(rexp(100,1))+1,
	min_expr=0.1,
	max_expr=2,
	casc.level=0.4
	)
plot(Net)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("network_random", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot-methods")
### * plot-methods

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot-methods
### Title: Plot
### Aliases: plot-methods plot,micro_array,ANY-method
###   plot,micropredict,ANY-method plot,network,ANY-method
### Keywords: methods

### ** Examples

if(require(CascadeData)){
data(micro_US, package="CascadeData")
micro_US<-as.micro_array(micro_US[1:100,],time=c(60,90,210,390),subject=6)
plot(micro_US)
}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot-methods", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotF")
### * plotF

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotF
### Title: Plot functions for the F matrices.
### Aliases: plotF
### Keywords: dplot

### ** Examples

#For numerical/inferred F matrices
plotF(CascadeFinit(4,4),choice="F", nround=1)
plotF(CascadeFinit(4,4),choice="Fpixmap")

#For theoritical F matrices
plotF(CascadeFshape(4,4),choice="Fshape")
plotF(CascadeFshape(4,4),choice="Fshapepixmap")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotF", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("position")
### * position

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: position
### Title: Retrieve network position for consistent plotting.
### Aliases: position
### Keywords: dplot

### ** Examples

data(network)
position(network)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("position", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predict-methods")
### * predict-methods

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predict
### Title: Methods for Function 'predict'
### Aliases: predict predict-methods predict,ANY-method
###   predict,micro_array-method
### Keywords: methods

### ** Examples




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predict-methods", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("probeMerge")
### * probeMerge

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: probeMerge
### Title: Function to merge probesets
### Aliases: probeMerge probeMerge,micro_array-method
### Keywords: manip

### ** Examples

if(require(CascadeData)){
data(micro_S)
D<-as.micro_array(micro_S[1:2000,],1:4,6)
D@gene_ID<-jetset::scores.hgu133plus2[D@name,"EntrezID"]
PM <- probeMerge(D)
}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("probeMerge", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("replaceBand")
### * replaceBand

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: replaceBand
### Title: Replace matrix values by band.
### Aliases: replaceBand
### Keywords: manip

### ** Examples

a=matrix(1:9,3,3)
b=matrix(0,3,3)
replaceBand(a,b,0)
replaceBand(a,b,1)
replaceBand(a,b,2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("replaceBand", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("replaceDown")
### * replaceDown

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: replaceDown
### Title: Replace matrix values triangular lower part and by band for the
###   upper part.
### Aliases: replaceDown
### Keywords: manip

### ** Examples

a=matrix(1:9,3,3)
b=matrix(1,3,3)
replaceDown(a,b,0)
replaceDown(a,b,1)
replaceDown(a,b,2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("replaceDown", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("replaceUp")
### * replaceUp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: replaceUp
### Title: Replace matrix values triangular upper part and by band for the
###   lower part.
### Aliases: replaceUp
### Keywords: manip

### ** Examples

a=matrix(1:9,3,3)
b=matrix(1,3,3)
replaceUp(a,b,0)
replaceUp(a,b,1)
replaceUp(a,b,2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("replaceUp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simul")
### * simul

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simul
### Title: Simulated data for examples.
### Aliases: simul M Net Net_inf_PL
### Keywords: datasets

### ** Examples

data(simul)
head(M)
str(M)
str(Net)
str(Net_inf_PL)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simul", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("unionMicro-methods")
### * unionMicro-methods

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: unionMicro-methods
### Title: Makes the union between two micro_array objects.
### Aliases: unionMicro unionMicro-methods
###   unionMicro,micro_array,micro_array-method unionMicro,list,ANY-method
### Keywords: methods

### ** Examples

if(require(CascadeData)){
data(micro_S, package="CascadeData")
#Create another microarray object with 100 genes
Mbis<-M<-as.micro_array(micro_S[1:100,],1:4,6)
#Rename the 100 genes
Mbis@name<-paste(M@name,"bis")
rownames(Mbis@microarray) <- Mbis@name
#Union (merge without duplicated names) of the two microarrays. 
str(unionMicro(M,Mbis))
}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("unionMicro-methods", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("unsupervised_clustering")
### * unsupervised_clustering

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: unsupervised_clustering
### Title: Cluster a micro_array object: performs the clustering.
### Aliases: unsupervised_clustering
###   unsupervised_clustering,micro_array,numeric,numeric-method
### Keywords: cluster

### ** Examples

if(require(CascadeData)){
data(micro_S, package="CascadeData")
M<-as.micro_array(micro_S[1:100,],1:4,6)
mc<-unsupervised_clustering_auto_m_c(M)
MwithGrp=unsupervised_clustering(M, 8, mc$m, screen=NULL, heatmap=FALSE, new.window = FALSE)
# Other options
unsupervised_clustering(M, 8, mc$m, screen=c(3,3), heatmap=TRUE, new.window = FALSE)
# Plot the clusters
plot(MwithGrp)
}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("unsupervised_clustering", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("unsupervised_clustering_auto_m_c")
### * unsupervised_clustering_auto_m_c

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: unsupervised_clustering_auto_m_c
### Title: Cluster a micro_array object: determine optimal fuzzification
###   parameter and number of clusters.
### Aliases: unsupervised_clustering_auto_m_c
###   unsupervised_clustering_auto_m_c,micro_array-method
### Keywords: cluster

### ** Examples

if(require(CascadeData)){
data(micro_S, package="CascadeData")
M<-as.micro_array(micro_S[1:100,],1:4,6)
mc<-unsupervised_clustering_auto_m_c(M)
}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("unsupervised_clustering_auto_m_c", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
