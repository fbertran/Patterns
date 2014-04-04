setClass(Class = "micro_array",
representation(microarray="matrix",name="vector",group=c("vector",NULL),start_time=c("vector",NULL),time=c("vector",NULL),subject="numeric"),
prototype = prototype(group = 0,start_time=0),
validity=function(object){

	if(dim(object@microarray)[2] != length(object@time)*object@subject){

				stop("[Error: ]Number of colomns must be equal to the number of time points * the number of subject")
		}
	if(dim(object@microarray)[1] != length(object@name)&&length(object@name)!=0){

				stop("[Error: ] Length of the vector of names must equal to the number of genes")
		}

	if(dim(object@microarray)[1] != length(object@group)&&length(object@group)!=1){

				print(object@group)
				stop("[Error: ] Length of the vector of group must equal to the number of genes or null")
		}
	
	if(dim(object@microarray)[1] != length(object@start_time)&&length(object@start_time)!=1){


				stop("[Error: ] Length of the vector of starting time must equal to the number of genes or null")
		}
		
	if(length(object@name)!=length(unique(object@name))){

				stop("[Error: ] Names must be unique")
		}
	if(object@subject<1){

				stop("[Error: ] There must be at least one subject")
		}	

	}

)

setMethod("head","micro_array",function(x,...)
{
	cat("The matrix :")
	cat("\n")
	cat("\n")
	K<-dim(x@microarray)[2]
	K<-min(K,3)
	
	print(head(x@microarray[,1:K]))
	cat("...")
	cat("\n")
	cat("\n")
	cat("Vector of names :")
	cat("\n")
	print(head(x@name))
	cat("...")
	cat("\n")
	cat("Vector of group :")
	cat("\n")
	print(head(x@group))
	cat("...")
	cat("\n")
	cat("Vector of starting time :")
	cat("\n")
	print(head(x@start_time))
	cat("...")
	cat("\n")
	cat("Vector of time :")
	cat("\n")
	print(x@time)
	cat("\n")
	cat("Number of subject :")
	cat("\n")
	print(x@subject)
}
)


setMethod("print","micro_array",function(x,...)
{
cat(paste("This is a micro_array S4 class. It contains : \n - (@microarray) a matrix of dimension ",dim(x@microarray)[1],"*",dim(x@microarray)[2],"\n          .... [gene expressions] \n - (@name) a vector of length ",length(x@name)," .... [gene names] \n","- (@group) a vector of length ",length(x@group)," .... [groups for genes] \n","- (@start_time) a vector of length ",length(x@start_time),"\n          .... [first differential expression for genes] \n","- (@time)a vector of length ",length(x@time)," .... [time points]\n","- (@subject) an integer  .... [number of subject]")) 	
})


setMethod("summary","micro_array",function(object,nb.graph=NULL,...)
{
	require(cluster)
	print(summary(object@microarray))
	G<-object@microarray
	colnames(G)<-paste(rep(paste("T",object@time),object@subject), as.character(rep(paste("subject",1:object@subject),each=length(object@time))))
	z <- cor(G)
	require(lattice)
	rownames(z)<-NULL
	if(is.null(nb.graph) || nb.graph==1 ){
	print(levelplot(z,xaxt="n",aspect="iso",xlab="",ylab=""))}
	if(is.null(nb.graph)){if(!attr(dev.cur(),"names")=="pdf"){dev.new()}}
	G<-object@microarray
	colnames(G)<-paste(rep(paste("T",object@time),object@subject), as.character(rep(paste("subject",1:object@subject),each=length(object@time))))
	w<-agnes(t(G))[1]$order
	G<-G[,w]
	z <- cor(G)
	rownames(z)<-NULL
	if(is.null(nb.graph) || nb.graph==2 ){
	print(levelplot(z,xaxt="n",aspect="iso",xlab="",ylab=""))}
	
	if(dim(object@microarray)[1]<1000){
		if(is.null(nb.graph)){if(!attr(dev.cur(),"names")=="pdf"){dev.new()}}
		
	R<-object@microarray
	w<-agnes(R)[1]$order
	R<-R[w,]
	z <- cor(t(R))
	rownames(z)<-NULL
	if(is.null(nb.graph) || nb.graph==3 ){
	print(levelplot(z,xaxt="n",aspect="iso",xlab="",ylab=""))}
	}
	})

setMethod("dim","micro_array",function(x)
{
	return(dim(x@microarray))
})
setMethod("plot","micro_array",function(x,y,...)
{
		
require(lattice)
	
xs<-t(x@microarray)

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
if(!attr(dev.cur(),"names")=="pdf"){dev.new()}
print(xyplot(V$conc~V$ys|V$suj,xlab="Time",ylab="Gene Expression",group=rep(1:(x@subject*dim(xs)[2]),each=length(x@time)),type="l",scales="free",col=rep(x@group,each=x@subject),key=list(
space="right",
lines=list(type="l",col=cclus),
text=list(text=paste("Cluster",as.character(cclus)))
)))
if(!attr(dev.cur(),"names")=="pdf"){dev.new()}
print(xyplot(V$conc~V$ys|gr,xlab="Time",ylab="Gene Expression",group=rep(1:(x@subject*dim(xs)[2]),each=length(x@time)),type="l",scales="free",col=rep(1:x@subject,dim(xs)[2]),key=list(
space="right",
lines=list(type="l",col=1:x@subject),
text=list(text=paste("Subject",as.character(1:x@subject)))
)))
for(i in 1:x@subject){
ss<-V$suj
sss<-ss==paste("Subject",i,sep=" ")
if(!attr(dev.cur(),"names")=="pdf"){dev.new()}
print(xyplot(V$conc[sss]~V$ys[sss]|gr[sss],xlab="Time",ylab="Gene Expression",type="l",main=paste("Subject",i,sep=" "),group=rep(1:(x@subject*dim(xs)[2]),each=length(x@time))[sss],scales="free"))
}
}
else{
	xyplot(V$conc~V$ys|V$suj,xlab="Time",group=rep(1:(x@subject*dim(xs)[2]),each=length(x@time)),ylab="Gene Expression",scales="free",type="l",col="black")
	}
}
)
 

setMethod(f="unsupervised_clustering_auto_m_c", 	
	signature=c("micro_array"),
	definition=function(M1,clust=NULL,mestim=NULL,M2=NULL,data_log=TRUE,screen=NULL,crange=NULL,repeats=NULL,cselect=TRUE,dminimum=FALSE){
require("Mfuzz")
if(is.null(crange)){crange=seq(4,32,4)}
if(is.null(repeats)){repeats=5}
if(is.null(M2)){
  if(data_log==TRUE){
			M<-log(M1@microarray)
	}
	else{
			M<-(M1@microarray)
	}
} else {
  if(data_log==TRUE){
			M1_mic<-log(M1@microarray)
			M2_mic<-log(M2@microarray)
	}
	else{
			M1_mic<-M1@microarray
			M2_mic<-M2@microarray
	}
  M=M1_mic-M2_mic
} 
  colnames(M)<-NULL
  Em<-ExpressionSet(M)
  Em.s <- standardise(Em)
if(is.null(mestim)){mestim <- mestimate(Em.s)
  }
if(cselect){
  if(is.null(clust)){
  tmp <- cselection(Em.s,m=mestim,crange=crange,repeats=repeats,visu=TRUE)
  clust <- ifelse(which.max(rowSums(t(tmp)-crange)<0)==1,crange[dim(tmp)[2]],crange[which.max(rowSums(t(tmp)-crange)<0)-1])
  } else {tmp=NA}
}
if(dminimum){
  if(is.null(clust)){
  tmp <- Dmin(Em.s,m=mestim,crange=crange,repeats=repeats,visu=TRUE)
  clust <- NA
  } else {tmp=NA}
}
  return(list(m=mestim,c=clust,csearch=tmp))
}
)



setMethod(f="unsupervised_clustering", 	
	signature=c("micro_array","numeric","numeric"),
	definition=function(M1,clust,mestim,M2=NULL,data_log=TRUE,screen=NULL,heatmap=TRUE){
require("Mfuzz")
if(is.null(M2)){
  if(data_log==TRUE){
			M<-log(M1@microarray)
	}
	else{
			M<-(M1@microarray)
	}
} else {
  if(data_log==TRUE){
			M1_mic<-log(M1@microarray)
			M2_mic<-log(M2@microarray)
	}
	else{
			M1_mic<-(M1@microarray)
			M2_mic<-(M2@microarray)
	}
  M=M1_mic-M2_mic
} 
  colnames(M)<-NULL
  Em<-ExpressionSet(M)
  Em.s <- standardise(Em)
  cl <- mfuzz(Em.s,c=clust,m=mestim)
  overlap.plot(cl,over =overlap(cl), thres = 0.05)
  if(!attr(dev.cur(),"names")=="pdf"){
  if(is.null(screen)){mfuzz.plot(Em.s,cl)
  } else {mfuzz.plot(Em.s,cl,mfrow=screen)}} else {
  if(is.null(screen)){mfuzz.plot(Em.s,cl,new.window=FALSE)
  } else {mfuzz.plot(Em.s,cl,mfrow=screen,new.window=)}}  
  if(heatmap){
  library(gplots) 
  if(!attr(dev.cur(),"names")=="pdf"){dev.new()}
  heatmap.2(exprs(Em.s), dendrogram = 'row', Colv = FALSE, col = greenred (75), 
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





 
 
setMethod(f="geneSelection", 
	signature=c("micro_array","numeric"),
	definition=function(M1,tot.number,M2=NULL,data_log=TRUE,wanted.patterns=NULL,forbidden.patterns=NULL,pic=NULL,alpha=0.05){
    indic<-0
  if(is.null(M2)){
  indic<-1
  M2<-M1
  }

		require(limma)		

		if(data_log==TRUE){
			M1_mic<-log(M1@microarray)
			M2_mic<-log(M2@microarray)
		}
		else{
			M1_mic<-(M1@microarray)
			M2_mic<-(M2@microarray)
		}
      if(indic==1){ M2_mic<-M2_mic*0}

		colnames(M1_mic)<-paste(rep("US",length(M1@time)*M1@subject),rep(M1@time,M1@subject), sep="")
		colnames(M2_mic)<-paste(rep("S",length(M2@time)*M2@subject),rep(M2@time,M2@subject), sep="")
		rownames(M1_mic)<- M1@name
		rownames(M2_mic)<- M1@name
    M<-cbind(M1_mic,M2_mic)

		T<-length(M1@time)
		
		#Construction de la matrice de design

		design<-t(matrix(rep(0,(M1@subject+M2@subject)*T*(T*2)),2*T))
		
		for(i in 1:(T)){
			design[which(colnames(M)%in%paste("US",M1@time[i],sep="")),i]<-1
			design[which(colnames(M)%in%paste("S",M2@time[i],sep="")),i+T]<-1	
		}


		vnom<-paste(c(paste(rep("US_time",length(M1@time)),M1@time[1:(length(M1@time))],sep=""),paste(rep("S_time",length(M1@time)),M1@time[1:(length(M1@time))],sep="")),sep="")
		colnames(design)<-vnom
		
		#block<-rep(1:(M1@subject*2),each=T)
		#dupcor <- duplicateCorrelation(M,design,block=block)
		#model<-lmFit(M,design,block=block)
		model<-lmFit(M,design)
		diff<-paste(vnom[1:T],vnom[1:T+length(M1@time)],sep="-")
		contr<-makeContrasts(contrasts=diff,levels=vnom)
		model2<-contrasts.fit(model,contr)
		model.test<-eBayes(model2)
		

		
p.val.ind<-model.test$p.value
rownames(p.val.ind)<-M1@name
p.val.all<-topTable(model.test,p.value=0.05,number=dim(M)[1])

ind.class<-rownames(p.val.all)
p.val.ind.class<-p.val.ind[ind.class,]
		
		f.p.val<-function(x){
			if(x<0.05){return(1)}
			else{return(0)}
		}

		p.val.ind.class.cat<-apply(p.val.ind.class,c(1,2),f.p.val)
		choix<-cbind(p.val.ind.class.cat,rep(0,dim(p.val.ind.class.cat)[1]))
		ch<-dim(choix)[2]
		f_test_patt<-function(x,pat){
				if(sum(abs(x-pat))==0){
					return(1)
				}
				else{
					return(0)		
				}}

		if(!is.null(forbidden.patterns)){			
			
			for(i in 1:dim(forbidden.patterns)[1]){
				f_test2<-function(x){f_test_patt(x,forbidden.patterns[i,])}
				S<-apply(choix[,1:(ch-1)],1,f_test2)
				choix<-choix[S==0,]
			}
		}


		if(!is.null(wanted.patterns)){
			sel<-rep(0,dim(choix)[1])
			choix[,ch]<-choix[,ch]*0
			for(i in 1:dim(wanted.patterns)[1]){
				f_test2<-function(x){f_test_patt(x,wanted.patterns[i,1:T])}
				sel<-sel+apply(choix[,1:(ch-1)],1,f_test2)

			}				
			for(j in 1:length(sel)){
					if((sum(choix[1:(j-1),ch])<tot.number || tot.number<0) && sel[j]>=1){choix[j,ch]<-1}
				}			
		}


		f_test_pic<-function(x,pic){
				x[pic]
				}

		if(!is.null(pic)){
				f_test2<-function(x){f_test_pic(x,pic)}
				
				# if(tot.number>0){
				# S<-apply(choix[,1:(ch-1)],1,f_test2)
				# for(j in 1:length(S)){
					# if(sum(S[1:j])<=tot.number && S[j]==1){choix[j,ch]<-1}
				# }
				# }
				# else{
					
					choix[,ch]<-apply(choix[,1:(ch-1)],1,f_test2)
					
					#}
				
			}
		
		
		if(tot.number>0 && is.null(wanted.patterns) ){
			i<-1
			PN<-dim(choix)[1]
			PN<-1:PN
			
			if(is.null(pic)){
			while(sum(choix[,ch])<tot.number){
				
				choix[i,ch]<-1
				i<-i+1
			}
			choix[PN>=i,ch]<-0
			choix<-choix[choix[,ch]==1,]
		}
		else{
			choix<-choix[choix[,ch]==1,]
			PN<-dim(choix)[1]
			choix<-choix[1:min(PN,tot.number),]
			}
		
		}

		if((!is.null(wanted.patterns))  ){choix<-choix[choix[,ch]==1,]}
		if(!is.null(pic)&&tot.number<=0){choix<-choix[choix[,ch]==1,]}
	
		temps.prem<-function(x){
			u<-0
			i<-1
			while(u==0){
				if(x[i]==1){return(i)}
				else
				i=i+1
			}
		}
		
		
		R<-apply(choix,1,temps.prem)
		n<-length(R)
		MM1<-M1_mic[rownames(choix),]-M2_mic[rownames(choix),]
		
		M<-new("micro_array",microarray=MM1,name=rownames(choix),time=M1@time,subject=M1@subject,group=as.vector(R),start_time=as.vector(R))
		
		
			
			return(M)

	}

)

setMethod(f="genePicSelection", 
	signature=c("micro_array","numeric"),
	definition=function(M1,pic,M2=NULL,data_log=TRUE,durPic=c(1,1),abs_val=TRUE,alpha_diff=0.05){
			
			
			Select<-geneSelection(M1,tot.number=-1,M2,data_log=data_log,pic=pic)
			
			if(abs_val==FALSE){
			M<-Select@microarray
			sel<-rep(0,length(M[,1]))
			comp<-M[,(1:M1@subject-1)*length(M1@time)+pic]
			test1<-function(y){
				if(wilcox.test(y[1:M1@subject],y[(M1@subject+1):(2*M1@subject)], alternative="less")$p.value<alpha_diff){
					return(1)
				}
				else{
					return(0)
					}
				}
				uu<-0
			if((pic-durPic[1])>0){
			for(i in 1:(pic-durPic[1])){
				comp1<-M[,(1:M1@subject-1)*length(M1@time)+i]
				comp2<-cbind(comp1,comp)
				sel<-sel+apply(comp2,1,test1)
				uu<-uu+1				
			}
			}
		
			if((pic+durPic[2])<=length(M1@time)){
			for(i in (pic+durPic[2]):length(M1@time)){
				comp1<-M[,(1:M1@subject-1)*length(M1@time)+i]
				comp2<-cbind(comp1,comp)
				sel<-sel+apply(comp2,1,test1)	
				uu<-uu+1			
			}
			}
			
			
			
			#print(uu)
			N1<-Select@name[sel==uu]
			
			M<--Select@microarray
			sel<-rep(0,length(M[,1]))
			comp<-M[,(1:M1@subject-1)*length(M1@time)+pic]
			test1<-function(y){
				if(wilcox.test(y[1:M1@subject],y[(M1@subject+1):(2*M1@subject)], alternative="less")$p.value<alpha_diff){
					return(1)
				}
				else{
					return(0)
					}
				}
				uu<-0
			if((pic-durPic[1])>0){
			for(i in 1:(pic-durPic[1])){
				comp1<-M[,(1:M1@subject-1)*length(M1@time)+i]
				comp2<-cbind(comp1,comp)
				sel<-sel+apply(comp2,1,test1)
				uu<-uu+1				
			}
			}
		
			if((pic+durPic[2])<=length(M1@time)){
			for(i in (pic+durPic[2]):length(M1@time)){
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
			Mi<-new("micro_array",microarray=Select@microarray[which(Select@name %in% N ),],name=N,time=M1@time,subject=M1@subject,start_time=Select@start_time[which(Select@name %in% N)],group=rep(pic,length(N)))			
			return(Mi)
			} 		
		else{
			M<-abs(Select@microarray)
			sel<-rep(0,length(M[,1]))
			comp<-M[,(1:M1@subject-1)*length(M1@time)+pic]
			test1<-function(y){
				if(wilcox.test(y[1:M1@subject],y[(M1@subject+1):(2*M1@subject)], alternative="less")$p.value<0.05){
					return(1)
				}
				else{
					return(0)
					}
				}
				uu<-0
			if((pic-durPic[1])>0){
			for(i in 1:(pic-durPic[1])){
				comp1<-M[,(1:M1@subject-1)*length(M1@time)+i]
				comp2<-cbind(comp1,comp)
				sel<-sel+apply(comp2,1,test1)
				uu<-uu+1				
			}
			}
		
			if((pic+durPic[2])<=length(M1@time)){
			for(i in (pic+durPic[2]):length(M1@time)){
				comp1<-M[,(1:M1@subject-1)*length(M1@time)+i]
				comp2<-cbind(comp1,comp)
				sel<-sel+apply(comp2,1,test1)	
				uu<-uu+1			
			}
			}
			
			N<-Select@name[sel==uu]
			#Select2<-geneSelection(M1,M2,tot.number=-1,data_log=data_log,wanted.patterns=t(as.matrix(c(0,1,0,0,1000))))
			#N2<-Select2@name
			
		
			Mi<-new("micro_array",microarray=Select@microarray[which(Select@name %in% N ),],name=N,time=M1@time,subject=M1@subject,group=rep(pic,length(N)),start_time=Select@start_time[which(Select@name %in% N)])	
			
			
		}
		
		
		
		}
		)
		

setMethod(f="unionMicro", 
	signature=c("micro_array","micro_array"),
	definition=function(M1,M2){
		
		nom1<-M1@name
		nom2<-M2@name
		NOM<-unique(c(nom1,nom2))
		n<-length(NOM)
		m1<-M1@microarray[which(nom1 %in% NOM),]
		NOM2<-NOM[-which(nom1 %in% NOM)]
		m2<-M2@microarray[which(nom2 %in% NOM2),]		
		M<-rbind(m1,m2)
		
		gr1<-M1@group[which(nom1 %in% NOM)]
		gr2<-M2@group[which(nom2 %in% NOM2)]
		gr<-c(gr1,gr2)
		str1<-M1@start_time[which(nom1 %in% NOM)]
		str2<-M2@start_time[which(nom2 %in% NOM2)]
		str<-c(str1,str2)
		rep<-new("micro_array",microarray=M,name=NOM,time=M1@time,subject=M1@subject,group=gr,start_time=str)
		return(rep)
				}
		
		) 
 