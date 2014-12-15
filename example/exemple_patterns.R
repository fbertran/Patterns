

#Exemple d'utilisation de patterns

###############################
data(CLL)

require(biomaRt)

affyids=c("202763_at","209310_s_at","207500_at")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
infos<-getBM(attributes=c("affy_hg_u133_plus_2","ensembl_gene_id","entrezgene","hgnc_symbol","chromosome_name","start_position","end_position","band"),
        filters = "affy_hg_u133_plus_2", values = CLL[,2], mart = ensembl,uniqueRows=TRUE, checkFilters = TRUE)
                 


hea_US<-CLL[,which((1:48)%%8<5&(1:48)%%8>0)+2]
hea_S<-CLL[,which(!((1:48)%%8<5&(1:48)%%8>0))+2]

agg_US<-CLL[,which((1:40)%%8<5&(1:40)%%8>0)+98]
agg_S<-CLL[,which(!((1:40)%%8<5&(1:40)%%8>0))+98]

m_hea_US<-as.micro_array(hea_US,c(60,90,210,390),6,name=CLL[,1],gene_ID=CLL[,2])
m_hea_S<- as.micro_array(hea_S,c(60,90,210,390),6,name=CLL[,1],gene_ID=CLL[,2])
  
m_agg_US<-as.micro_array((agg_US),c(60,90,210,390),5,name=CLL[,1],gene_ID=CLL[,2])
m_agg_S<- as.micro_array((agg_S),c(60,90,210,390),5,name=CLL[,1],gene_ID=CLL[,2])

selection1<-geneSelection(list(m_agg_US,m_agg_S),list("condition&time",c(1,2),c(1,1)),-1,alpha=0.1)
selection2<-geneSelection(list(m_agg_US,m_agg_S),list("condition&time",c(1,2),c(1,1)+1),-1,alpha=0.1)
selection3<-geneSelection(list(m_agg_US,m_agg_S),list("condition&time",c(1,2),c(1,1)+2),-1,alpha=0.05)
selection4<-geneSelection(list(m_agg_US,m_agg_S),list("condition&time",c(1,2),c(1,1)+3),-1,alpha=0.05)

selection<-unionMicro(list(selection1,selection2,selection3,selection4))

length(selection@gene_ID)
save(list=c("selection","infos"),file="travail.RData")

################################################################



require(biomaRt)
data(travail)


#################################################################

#Etude des facteurs de transcription


library(XML)

doc <- readHTMLTable(
  doc=htmlParse("http://www.bioguo.org/AnimalTFDB/BrowseAllTF.php?spe=Homo_sapiens",encoding = "UTF-8")
)

TF<-unique(as.character(doc[[5]][,4]))
TF<-TF[order(TF)]
TF<-TF[-1]


#TF contient la liste des gene ID des FT

tfs<-which(selection@gene_ID %in% TF)

matplot(t(selection@microarray[tfs,]),type="l",lty=1)
kk<-kmeans((selection@microarray[tfs,]),10)
matplot(t(kk$centers),type="l",lty=1)

TFi<-function(x) length(which(TF %in% x))



n<-80
kre<-kmeans(selection@microarray,n)
kre
lll<-split(selection@name,kre$cluster)


require("clusterProfiler")
require("AnnotationFuncs")
require(org.Hs.eg.db)
require(DCGL)
pp<-list()

for(k in 1:n){
  print(k)
pp[[k]]<-translate(lll[[k]],from=org.Hs.egSYMBOL2EG,simplify=TRUE)
# GOs[[k]]<-enrichGO(pp, organism = "human", ont = "MF", pvalueCutoff = 0.05,
#          pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 5,
#          readable = FALSE)

}

names(pp)<-paste("X",1:n,sep="")
test<-compareCluster(pp,fun="enrichGO", organism="human", pvalueCutoff=0.15)

plot(test)






TFu<-(unlist(lapply(pp,TFi)))
TFy<-unlist(lapply(pp,length))

plot(TFu/TFy)
plot(TFu)
sum(TFu)

entrez<-translate(selection@gene_ID,from=org.Hs.egSYMBOL2EG,simplify=TRUE)

geneName<-translate(entrez[which(TF %in% entrez)],from=org.Hs.egSYMBOL,simplify=TRUE)

which(selection@gene_ID %in% "EGR1")

