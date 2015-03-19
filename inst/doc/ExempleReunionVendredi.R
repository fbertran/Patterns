### R code from vignette source 'Patterns.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Patterns.Rnw:28-30
###################################################
options(width=60)
options(continue="   ")


###################################################
### code chunk number 2: Patterns.Rnw:78-79 (eval = FALSE)
###################################################
## install.packages("name_of_the_package")


###################################################
### code chunk number 3: Patterns.Rnw:82-84 (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("name_of_the_package")


###################################################
### code chunk number 4: Patterns.Rnw:91-92
###################################################
library(Patterns)


###################################################
### code chunk number 5: Patterns.Rnw:100-102
###################################################
data(micro_S)
data(micro_US)


###################################################
### code chunk number 6: Patterns.Rnw:111-112
###################################################
colnames(micro_S)


###################################################
### code chunk number 7: Patterns.Rnw:116-118
###################################################
micro_S<-as.micro_array(micro_S,time=c(60,90,210,390),subject=6)
micro_US<-as.micro_array(micro_US,time=c(60,90,210,390),subject=6)


###################################################
### code chunk number 8: Patterns.Rnw:124-125
###################################################
print(micro_S)


###################################################
### code chunk number 9: Patterns.Rnw:130-131
###################################################
head(micro_S)


###################################################
### code chunk number 10: Patterns.Rnw:148-149
###################################################
Selection<-geneSelection(x=micro_S,y=micro_US,tot.number=50,data_log=TRUE)


###################################################
### code chunk number 11: Patterns.Rnw:162-163
###################################################
summary(Selection)


###################################################
### code chunk number 12: Patterns.Rnw:169-170
###################################################
summary(Selection,nb.graph=2)


###################################################
### code chunk number 13: Patterns.Rnw:173-174
###################################################
summary(Selection,nb.graph=3)


###################################################
### code chunk number 14: Patterns.Rnw:188-193
###################################################
#If we want to select genes that are differentially 
#at time t60 or t90 :
Selection<-geneSelection(x=micro_S,y=micro_US,tot.number=30,
wanted.patterns=
rbind(c(0,1,0,0),c(1,0,0,0),c(1,1,0,0)))


###################################################
### code chunk number 16: Patterns.Rnw:207-228
###################################################
#Select early genes (t1 or t2)
Selection1<-geneSelection(x=micro_S,y=micro_US,20,
wanted.patterns=
rbind(c(0,1,0,0),c(1,0,0,0),c(1,1,0,0)))
Selection1@gene_ID<-Selection1@name
#Section genes with first significant differential 
#expression at t1:

Selection2<-geneSelection(x=micro_S,y=micro_US,20,
pic=1)
Selection2@gene_ID<-Selection2@name
#Section genes with first significant differential 
#expression at t2:

Selection3<-geneSelection(x=micro_S,y=micro_US,20,
pic=2)
Selection3@gene_ID<-Selection3@name
#Select later genes (t3 or t4)
Selection4<-geneSelection(x=micro_S,y=micro_US,50,
wanted.patterns=
rbind(c(0,0,1,0),c(0,0,0,1),c(1,1,0,0)))
Selection4@gene_ID<-Selection4@name

###################################################
### code chunk number 17: Patterns.Rnw:233-237
###################################################
Selection<-unionMicro(Selection1,Selection2)
Selection<-unionMicro(Selection,Selection3)
Selection<-unionMicro(Selection,Selection4)
print(Selection)


###################################################
### code chunk number 18: Patterns.Rnw:239-241
###################################################
#Prints the correlation graphics Figure 4:
summary(Selection,3)


###################################################
### code chunk number 19: Patterns.Rnw:246-247
###################################################
summary(Selection,3)


###################################################
### code chunk number 20: Patterns.Rnw:254-256
###################################################
Selection2gp<-unionMicro(Selection1,Selection2)
Selection2gp<-unionMicro(Selection2gp,Selection3)


###################################################
### code chunk number 21: Patterns.Rnw:347-350 (eval = FALSE)
###################################################
network<-inference(Selection,type.inf="noniterative")
networkCascade<-inference(Selection,Finit=CascadeFinit(4,4),Fshape=CascadeFshape(4,4),type.inf="noniterative")
network2gp<-inference(Selection2gp,type.inf="noniterative")


###################################################
### code chunk number 22: Patterns.Rnw:352-353
###################################################
#load(system.file("extdata", "networks.Rdata", package = "Patterns"))


###################################################
### code chunk number 23: Patterns.Rnw:358-363 (eval = FALSE)
###################################################
## plot(network,choice="F")
## plot(networkCascade,choice="F")
## plot(network2gp,choice="F")
## plot(network,choice="network",gr=Selection@group)
## plot(network2gp,choice="network",gr=Selection2gp@group)


###################################################
### code chunk number 24: Patterns.Rnw:369-370
###################################################

plot(network,choice="F")


###################################################
### code chunk number 25: Patterns.Rnw:373-374
###################################################
plot(networkCascade,choice="F")


###################################################
### code chunk number 26: Patterns.Rnw:377-378
###################################################
plot(network2gp,choice="F")


###################################################
### code chunk number 27: Patterns.Rnw:381-382
###################################################
plot(network,choice="network",gr=Selection@group)


###################################################
### code chunk number 28: Patterns.Rnw:385-386
###################################################
plot(network2gp,choice="network",gr=Selection2gp@group)


###################################################
### code chunk number 29: Patterns.Rnw:398-400 (eval = FALSE)
###################################################
evolution(network,seq(0,0.4,by=0.01),gr=Selection@group,fix=TRUE)
evolution(network,seq(0,0.4,by=0.01),gr=Selection@group,fix=FALSE)


###################################################
### code chunk number 30: Patterns.Rnw:414-416
###################################################
#evol_cutoff<-cutoff(network)
nv<-0.07


###################################################
### code chunk number 31: Patterns.Rnw:424-425
###################################################
plot(network,choice="network",gr=Selection@group,nv=0.07)


###################################################
### code chunk number 32: Patterns.Rnw:428-429
###################################################
#plot(evol_cutoff$sequence,evol_cutoff$p.value.inter,type="l",xlab="cutoff",ylab="p.value")


###################################################
### code chunk number 33: Patterns.Rnw:452-454
###################################################
#analyze<-analyze_network(network,nv)
#head(analyze)


###################################################
### code chunk number 34: Patterns.Rnw:461-462 (eval = FALSE)
###################################################
## plot(network,nv=nv,gr=Selection@group,ani=TRUE)


###################################################
### code chunk number 35: Patterns.Rnw:470-475 (eval = FALSE)
###################################################
## P<-position(network,nv=nv)
## #plotting the network with the group coloring:
## plot(network,nv=nv,gr=Selection@group,ini=P)
## #plotting the network without the group coloring:
## plot(network,nv=nv,ini=P)


###################################################
### code chunk number 36: Patterns.Rnw:477-478
###################################################
P<-position(network,nv=nv)


###################################################
### code chunk number 37: Patterns.Rnw:489-493 (eval = FALSE)
###################################################
## geneNeighborhood(network,targets=16,nv=nv,ini=P,
## 	label.hub=TRUE,label_v=Selection@name)
## #label.hub: only hubs vertex should have a name
## #label_v: name of the vertex


###################################################
### code chunk number 38: Patterns.Rnw:499-501
###################################################
# prediction_ko16<-predict(Selection,network,nv=nv,targets=16)
# prediction_ko16_2gp<-predict(Selection2gp,network2gp,nv=nv,targets=16)


###################################################
### code chunk number 39: Patterns.Rnw:505-508 (eval = FALSE)
###################################################
## #We plot the results ; here for example we see changes at time point t2
## plot(prediction_ko16,time=2,ini=P,label.hub=TRUE,label_v=Selection@name)
## plot(prediction_ko16_2gp,time=2,ini=P,label.hub=TRUE,label_v=Selection2gp@name)


###################################################
### code chunk number 40: Patterns.Rnw:513-514
###################################################
geneNeighborhood(network,targets=16,nv=nv,ini=P,label.hub=TRUE,label_v=Selection@name)


###################################################
### code chunk number 41: Patterns.Rnw:517-518
###################################################
#plot(prediction_ko16,time=2,ini=P,label.hub=TRUE,label_v=Selection@name)


###################################################
### code chunk number 42: Patterns.Rnw:521-522
###################################################
#plot(prediction_ko16_2gp,time=2,ini=P,label.hub=TRUE,label_v=Selection2gp@name)


###################################################
### code chunk number 43: keep.source
###################################################
#We set the seed to make the results reproducible 
# set.seed(1)
# 
# #We create a random scale free network
# Net<-network_random(
# 	nb=100,
# 	time_label=rep(1:4,each=25),
# 	exp=1,
# 	init=1,
# 	regul=round(rexp(100,1))+1,
# 	min_expr=0.1,
# 	max_expr=2,
# 	casc.level=0.4
# 	)
# 
# #We change F matrices
# ngrp<-4
# T<-4 
# F<-array(0,c(T,T,ngrp*ngrp))
# for(i in 1:(ngrp*ngrp)){diag(F[,,i])<-1}
# F[,,2]<-F[,,2]*0.2
# F[2,1,2]<-1
# F[3,2,2]<-1
# F[,,4]<-F[,,2]*0.3
# F[3,1,4]<-1
# F[,,5]<-F[,,2]
# Net@F<-F
# 
# #We simulate gene expression according to the network Net
# M<-gene_expr_simulation(
# 	network=Net,
# 	time_label=rep(1:4,each=25),
# 	subject=5,
# 	level_pic=200)
# 
# 
# ###################################################
# ### code chunk number 44: keep.source
# ###################################################
# load(system.file("extdata", "Net_inf2.RData", package = "Patterns"))
# Net_inf<-Net_inf2
# 
# 
# ###################################################
# ### code chunk number 45: keep.source (eval = FALSE)
# ###################################################
# ## #We infer the new network
# ## Net_inf<-inference(M)
# 
# 
# ###################################################
# ### code chunk number 46: keep.source
# ###################################################
# #Comparing true and inferred networks
# F_score<-rep(0,200)
# #Here are the cutoff level tested
# test.seq<-seq(0,max(abs(Net_inf@network*0.9)),length.out=200)
# 
# u<-0
# for(i in test.seq){
# 	u<-u+1
# 	F_score[u]<-compare(Net,Net_inf,i)[3]
# }
# 
# 
# ###################################################
# ### code chunk number 47: keep.source
# ###################################################
# #Choosing the cutoff
# cut.seq<-cutoff(Net_inf)
# plot(cut.seq$sequence,cut.seq$p.value.inter)
# 
# 
# ###################################################
# ### code chunk number 48: Patterns.Rnw:606-608
# ###################################################
# plot(cut.seq$sequence,cut.seq$p.value.inter,type="l",xlab="cutoff",ylab="p.value")
# abline(v=0.24,col="red")
# 
# 
# ###################################################
# ### code chunk number 49: Patterns.Rnw:611-613
# ###################################################
# plot(test.seq,F_score,type="b",xlab="cutoff",ylab="Fscore")
# abline(v=0.24,col="red")
# 
# 
