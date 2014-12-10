

#Exemple d'utilisation de patterns


data(CLL)


hea_US<-CLL[,which((1:48)%%8<5&(1:48)%%8>0)+2]
hea_S<-CLL[,which(!((1:48)%%8<5&(1:48)%%8>0))+2]


agg_US<-CLL[,which((1:40)%%8<5&(1:40)%%8>0)+98]
agg_S<-CLL[,which(!((1:40)%%8<5&(1:40)%%8>0))+98]


m_hea_US<-as.micro_array(hea_US,c(60,90,210,390),6,name=CLL[,2],gene_ID=CLL[,1])
m_hea_S<- as.micro_array(hea_S,c(60,90,210,390),6,name=CLL[,2],gene_ID=CLL[,1])
  
m_agg_US<-as.micro_array(agg_US,c(60,90,210,390),5,name=CLL[,2],gene_ID=CLL[,1])
m_agg_S<- as.micro_array(agg_S,c(60,90,210,390),5,name=CLL[,2],gene_ID=CLL[,1])




geneSelection()


