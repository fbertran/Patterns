### R code from vignette source 'Patterns.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: Patterns.Rnw:27-29
###################################################
options(width=60)
options(continue="   ")


###################################################
### code chunk number 2: Patterns.Rnw:77-78 (eval = FALSE)
###################################################
## install.packages("name_of_the_package")


###################################################
### code chunk number 3: Patterns.Rnw:81-83 (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("name_of_the_package")


###################################################
### code chunk number 4: Patterns.Rnw:90-91
###################################################
library(Patterns)


###################################################
### code chunk number 5: Patterns.Rnw:99-101
###################################################
data(micro_S)
data(micro_US)


###################################################
### code chunk number 6: Patterns.Rnw:110-111
###################################################
colnames(micro_S)


###################################################
### code chunk number 7: Patterns.Rnw:115-117
###################################################
micro_S<-as.micro_array(micro_S,time=c(60,90,210,390),subject=6)
micro_US<-as.micro_array(micro_US,time=c(60,90,210,390),subject=6)


###################################################
### code chunk number 8: Patterns.Rnw:123-124
###################################################
print(micro_S)


###################################################
### code chunk number 9: Patterns.Rnw:129-130
###################################################
head(micro_S)


###################################################
### code chunk number 10: Patterns.Rnw:147-148
###################################################
Selection<-geneSelection(M1=micro_S,M2=micro_US,tot.number=50,data_log=TRUE)


###################################################
### code chunk number 11: Patterns.Rnw:161-162
###################################################
summary(Selection)


###################################################
### code chunk number 12: Patterns.Rnw:168-169
###################################################
summary(Selection,nb.graph=2)


###################################################
### code chunk number 13: Patterns.Rnw:172-173
###################################################
summary(Selection,nb.graph=3)


###################################################
### code chunk number 14: Patterns.Rnw:187-192
###################################################
#If we want to select genes that are differentially 
#at time t60 or t90 :
Selection<-geneSelection(M1=micro_S,M2=micro_US,tot.number=30,
wanted.patterns=
rbind(c(0,1,0,0),c(1,0,0,0),c(1,1,0,0)))


###################################################
### code chunk number 15: Patterns.Rnw:199-201
###################################################
Selection<-genePicSelection(M1=micro_S,M2=micro_US,1,
abs_val=FALSE,alpha_diff=0.01)


###################################################
### code chunk number 16: Patterns.Rnw:206-227
###################################################
#Select early genes (t1 or t2)
Selection1<-geneSelection(M1=micro_S,M2=micro_US,20,
wanted.patterns=
rbind(c(0,1,0,0),c(1,0,0,0),c(1,1,0,0)))

#Section genes with first significant differential 
#expression at t1:

Selection2<-geneSelection(M1=micro_S,M2=micro_US,20,
pic=1)

#Section genes with first significant differential 
#expression at t2:

Selection3<-geneSelection(M1=micro_S,M2=micro_US,20,
pic=2)

#Select later genes (t3 or t4)
Selection4<-geneSelection(M1=micro_S,M2=micro_US,50,
wanted.patterns=
rbind(c(0,0,1,0),c(0,0,0,1),c(1,1,0,0)))


###################################################
### code chunk number 17: Patterns.Rnw:232-236
###################################################
Selection<-unionMicro(Selection1,Selection2)
Selection<-unionMicro(Selection,Selection3)
Selection<-unionMicro(Selection,Selection4)
print(Selection)


###################################################
### code chunk number 18: Patterns.Rnw:238-240
###################################################
#Prints the correlation graphics Figure 4:
summary(Selection,3)


###################################################
### code chunk number 19: Patterns.Rnw:245-246
###################################################
summary(Selection,3)


###################################################
### code chunk number 20: Patterns.Rnw:253-255
###################################################
Selection2gp<-unionMicro(Selection1,Selection2)
Selection2gp<-unionMicro(Selection2gp,Selection3)


###################################################
### code chunk number 21: Patterns.Rnw:346-349 (eval = FALSE)
###################################################
## network<-inference(Selection)
## networkCascade<-inference(Selection,Finit=CascadeFinit(4,4),Fshape=CascadeFshape(4,4))
## network2gp<-inference(Selection2gp)


###################################################
### code chunk number 22: Patterns.Rnw:351-352
###################################################
load(system.file("extdata", "networks.Rdata", package = "Patterns"))


###################################################
### code chunk number 23: Patterns.Rnw:357-362 (eval = FALSE)
###################################################
## plot(network,choice="F")
## plot(networkCascade,choice="F")
## plot(network2gp,choice="F")
## plot(network,choice="network",gr=Selection@group)
## plot(network2gp,choice="network",gr=Selection2gp@group)


###################################################
### code chunk number 24: Patterns.Rnw:368-369
###################################################
plot(network,choice="F")


###################################################
### code chunk number 25: Patterns.Rnw:372-373
###################################################
plot(networkCascade,choice="F")


###################################################
### code chunk number 26: Patterns.Rnw:376-377
###################################################
plot(network2gp,choice="F")


###################################################
### code chunk number 27: Patterns.Rnw:380-381
###################################################
plot(network,choice="network",gr=Selection@group)


###################################################
### code chunk number 28: Patterns.Rnw:384-385
###################################################
plot(network2gp,choice="network",gr=Selection2gp@group)


###################################################
### code chunk number 29: Patterns.Rnw:397-399 (eval = FALSE)
###################################################
## evolution(network,seq(0,0.4,by=0.01),gr=Selection@group,fix=TRUE)
## evolution(network,seq(0,0.4,by=0.01),gr=Selection@group,fix=FALSE)


###################################################
### code chunk number 30: Patterns.Rnw:413-415
###################################################
evol_cutoff<-cutoff(network)
nv<-0.07


###################################################
### code chunk number 31: Patterns.Rnw:423-424
###################################################
plot(network,choice="network",gr=Selection@group,nv=0.07)


###################################################
### code chunk number 32: Patterns.Rnw:427-428
###################################################
plot(evol_cutoff$sequence,evol_cutoff$p.value.inter,type="l",xlab="cutoff",ylab="p.value")


###################################################
### code chunk number 33: Patterns.Rnw:451-453
###################################################
analyze<-analyze_network(network,nv)
head(analyze)


###################################################
### code chunk number 34: Patterns.Rnw:460-461 (eval = FALSE)
###################################################
## plot(network,nv=nv,gr=Selection@group,ani=TRUE)


###################################################
### code chunk number 35: Patterns.Rnw:469-474 (eval = FALSE)
###################################################
## P<-position(network,nv=nv)
## #plotting the network with the group coloring:
## plot(network,nv=nv,gr=Selection@group,ini=P)
## #plotting the network without the group coloring:
## plot(network,nv=nv,ini=P)


###################################################
### code chunk number 36: Patterns.Rnw:476-477
###################################################
P<-position(network,nv=nv)


###################################################
### code chunk number 37: Patterns.Rnw:488-492 (eval = FALSE)
###################################################
## geneNeighborhood(network,targets=16,nv=nv,ini=P,
## 	label.hub=TRUE,label_v=Selection@name)
## #label.hub: only hubs vertex should have a name
## #label_v: name of the vertex


###################################################
### code chunk number 38: Patterns.Rnw:498-500
###################################################
prediction_ko16<-predict(Selection,network,nv=nv,targets=16)
prediction_ko16_2gp<-predict(Selection2gp,network2gp,nv=nv,targets=16)


###################################################
### code chunk number 39: Patterns.Rnw:504-507 (eval = FALSE)
###################################################
## #We plot the results ; here for example we see changes at time point t2
## plot(prediction_ko16,time=2,ini=P,label.hub=TRUE,label_v=Selection@name)
## plot(prediction_ko16_2gp,time=2,ini=P,label.hub=TRUE,label_v=Selection2gp@name)


###################################################
### code chunk number 40: Patterns.Rnw:512-513
###################################################
geneNeighborhood(network,targets=16,nv=nv,ini=P,label.hub=TRUE,label_v=Selection@name)


###################################################
### code chunk number 41: Patterns.Rnw:516-517
###################################################
plot(prediction_ko16,time=2,ini=P,label.hub=TRUE,label_v=Selection@name)


###################################################
### code chunk number 42: Patterns.Rnw:520-521
###################################################
plot(prediction_ko16_2gp,time=2,ini=P,label.hub=TRUE,label_v=Selection2gp@name)


###################################################
### code chunk number 43: keep.source
###################################################
#We set the seed to make the results reproducible 
set.seed(1)

#We create a random scale free network
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

#We change F matrices
ngrp<-4
T<-4 
F<-array(0,c(T,T,ngrp*ngrp))
for(i in 1:(ngrp*ngrp)){diag(F[,,i])<-1}
F[,,2]<-F[,,2]*0.2
F[2,1,2]<-1
F[3,2,2]<-1
F[,,4]<-F[,,2]*0.3
F[3,1,4]<-1
F[,,5]<-F[,,2]
Net@F<-F

#We simulate gene expression according to the network Net
M<-gene_expr_simulation(
	network=Net,
	time_label=rep(1:4,each=25),
	subject=5,
	level_pic=200)


###################################################
### code chunk number 44: keep.source
###################################################
load(system.file("extdata", "Net_inf2.RData", package = "Patterns"))
Net_inf<-Net_inf2


###################################################
### code chunk number 45: keep.source (eval = FALSE)
###################################################
## #We infer the new network
## Net_inf<-inference(M)


###################################################
### code chunk number 46: keep.source
###################################################
#Comparing true and inferred networks
F_score<-rep(0,200)
#Here are the cutoff level tested
test.seq<-seq(0,max(abs(Net_inf@network*0.9)),length.out=200)

u<-0
for(i in test.seq){
	u<-u+1
	F_score[u]<-compare(Net,Net_inf,i)[3]
}


###################################################
### code chunk number 47: keep.source
###################################################
#Choosing the cutoff
cut.seq<-cutoff(Net_inf)
plot(cut.seq$sequence,cut.seq$p.value.inter)


###################################################
### code chunk number 48: Patterns.Rnw:605-607
###################################################
plot(cut.seq$sequence,cut.seq$p.value.inter,type="l",xlab="cutoff",ylab="p.value")
abline(v=0.24,col="red")


###################################################
### code chunk number 49: Patterns.Rnw:610-612
###################################################
plot(test.seq,F_score,type="b",xlab="cutoff",ylab="Fscore")
abline(v=0.24,col="red")


