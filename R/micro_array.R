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
           
           #	if(length(object@name)!=length(unique(object@name))){ #Disappeared in Cascade 1.03
           #
           #				stop("[Error: ] Names must be unique")
           #		}
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
  #rownames(z)<-NULL
  if(is.null(nb.graph) || nb.graph==1 ){
    print(levelplot(z,aspect="iso",xlab="",ylab="", 
                    scales=list(x=list(rot=90)), 
                    ylab.right = "Level of correlation",
                    par.settings = list(layout.widths = list(axis.key.padding = 0,
                                                             ylab.right = 2))))
  }
  if(is.null(nb.graph)){dev.new()}
  G<-object@microarray
  colnames(G)<-paste(rep(paste("T",object@time),object@subject), as.character(rep(paste("subject",1:object@subject),each=length(object@time))))
  w<-agnes(t(G))[1]$order
  G<-G[,w]
  z <- cor(G)
  
  if(is.null(nb.graph) || nb.graph==2 ){
    print(levelplot(z,aspect="iso",xlab="",ylab="", 
                    scales=list(x=list(rot=90)), 
                    ylab.right = "Level of correlation",
                    par.settings = list(layout.widths = list(axis.key.padding = 0,
                                                             ylab.right = 2))))
  }
  if(dim(object@microarray)[1]<1000){
    if(is.null(nb.graph)){dev.new()}
    
    R<-object@microarray
    w<-agnes(R)[1]$order
    R<-R[w,]
    z <- cor(t(R))
    
    if(is.null(nb.graph) || nb.graph==3 ){
      print(levelplot(z,aspect="iso",xlab="",ylab="", 
                      scales=list(x=list(rot=90)), 
                      ylab.right = "Level of correlation",
                      par.settings = list(layout.widths = list(axis.key.padding = 0,
                                                               ylab.right = 2))))
    }
  }}
)

setMethod("dim","micro_array",function(x)
{
  return(dim(x@microarray))
})
setMethod("plot","micro_array",function(x,y,...)
{
  
  require(lattice)
  require(grid)	
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
    print(xyplot(V$conc~V$ys|V$suj,as.table=TRUE,xlab="Time",ylab="Gene Expression",group=rep(1:(x@subject*dim(xs)[2]),each=length(x@time)),type="l",scales=list(x=list(relation="free",at=x@time),y=list(relation="free")),col=rep(x@group,each=x@subject),key=list(
      space="right",
      lines=list(type="l",col=cclus),
      text=list(text=paste("Cluster",as.character(cclus)))
    )))
    if(!attr(dev.cur(),"names")=="pdf"){dev.new()}
    print(xyplot(V$conc~V$ys|gr,as.table=TRUE,xlab="Time",ylab="Gene Expression",group=rep(1:(x@subject*dim(xs)[2]),each=length(x@time)),type="l",scales=list(x=list(relation="free",at=x@time),y=list(relation="free")),col=rep(1:x@subject,dim(xs)[2]),key=list(
      space="right",
      lines=list(type="l",col=1:x@subject),
      text=list(text=paste("Subject",as.character(1:x@subject)))
    )))
    for(i in 1:x@subject){
      ss<-V$suj
      sss<-ss==paste("Subject",i,sep=" ")
      if(!attr(dev.cur(),"names")=="pdf"){dev.new()}
      print(xyplot(V$conc[sss]~V$ys[sss]|gr[sss],as.table=TRUE,xlab="Time",ylab="Gene Expression",type="l",main=paste("Subject",i,sep=" "),group=rep(1:(x@subject*dim(xs)[2]),each=length(x@time))[sss],scales=list(x=list(relation="free",at=x@time),y=list(relation="free"))))
    }
  }
  else{
    xyplot(V$conc~V$ys|V$suj,xlab="Time",as.table=TRUE,group=rep(1:(x@subject*dim(xs)[2]),each=length(x@time)),ylab="Gene Expression",scales=list(x=list(relation="free",at=x@time),y=list(relation="free")),type="l",col="black")
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
          signature=c("micro_array","micro_array","numeric"),
          definition=function(x,y,tot.number,data_log=TRUE,wanted.patterns=NULL,forbidden.patterns=NULL,pic=NULL,alpha=0.05,Design=NULL,lfc=0){
            
            BBB<-strsplit(sessionInfo()[5]$otherPkgs$limma$Version,"[.]")
            
            if( !(BBB[[1]][1]>3 || (BBB[[1]][1]==3 && BBB[[1]][2]>18) || 
                    (BBB[[1]][1]==3 && BBB[[1]][2]==18 && BBB[[1]][3]>=13 ) ))
            {stop("Upgrade your version of Limma (>= 3.18.13)")}
            
            
            indic<-0
            M1<-x
            M2<-y
            if(is.null(M2)){
              indic<-1
              M2<-M1
            }
            
            require(limma)		
            
            if(data_log==TRUE){
              M1_mic<-log(M1@microarray)
              M2_mic<-log(M2@microarray)
            } else{
              M1_mic<-(M1@microarray)
              M2_mic<-(M2@microarray)
            }
            
            if(is.null(rownames(M1_mic))){rownames(M1_mic)<-paste("probe ",1:dim(M1_mic)[1])}
            if(is.null(rownames(M2_mic))){rownames(M2_mic)<-paste("probe ",1:dim(M2_mic)[1])}
            
            if(indic==1){ M2_mic<-M2_mic*0}
            
            
            
            colnames(M1_mic)<-paste(rep("US",length(M1@time)*M1@subject),rep(M1@time,M1@subject), sep="")
            colnames(M2_mic)<-paste(rep("S",length(M2@time)*M2@subject),rep(M2@time,M2@subject), sep="")
            #rownames(M1_mic)<- paste("probe",rownames(M1_mic) ) 
            #rownames(M2_mic)<-  paste("probe",rownames(M2_mic) ) 
            M<-cbind(M1_mic,M2_mic)
            
            T<-length(M1@time)
            
            #Construction de la matrice de design
            
            if(is.null(Design)){
              design<-t(matrix(rep(0,(M1@subject+M2@subject)*T*(T*2)),2*T))
              
              for(i in 1:(T)){
                design[which(colnames(M)%in%paste("US",M1@time[i],sep="")),i]<-1
                design[which(colnames(M)%in%paste("S",M2@time[i],sep="")),i+T]<-1	
              }
              
              
              vnom<-paste(c(paste(rep("US_time",length(M1@time)),M1@time[1:(length(M1@time))],sep=""),paste(rep("S_time",length(M1@time)),M1@time[1:(length(M1@time))],sep="")),sep="")
              colnames(design)<-vnom
            }else{
              design<-Desing$X		
            }
            #block<-rep(1:(M1@subject*2),each=T)
            #dupcor <- duplicateCorrelation(M,design,block=block)
            #model<-lmFit(M,design,block=block)
            model<-lmFit(M,design)
            if(is.null(Design)){
              diff<-paste(vnom[1:T],vnom[1:T+length(M1@time)],sep="-")
              contr<-makeContrasts(contrasts=diff,levels=vnom)  		}else{
                contr<-Design$contr
              }
            model2<-contrasts.fit(model,contr)
            model.test<-eBayes(model2)
            
            
            
            p.val.ind<-model.test$p.value
            rownames(p.val.ind)<-rownames(M1_mic) 
            
            p.val.all<-topTable(model.test,p.value=alpha,number=dim(M)[1],lfc=lfc)
            kkkk<-0
            if(is.null(p.val.all$ID)){
              kkkk<-1
              ID<-row.names(p.val.all)
              p.val.all<-cbind(ID,p.val.all)
              f_nom<-function(nom){
                which(x@name %in% nom)			
              }
              f_nom<-Vectorize(f_nom)
              row.names(p.val.all)<-unlist(f_nom(as.character(p.val.all$ID)))
            }
            if(tot.number>0 && tot.number<=1){
              tot.number<-round(tot.number*dim(p.val.all)[1])
            }
            
            if(kkkk<-0){
              ind.class<-rownames(p.val.all)
              p.val.ind.class<-p.val.ind[(ind.class),]
            }else{
              ind.class<-rownames(p.val.all)
              p.val.ind.class<-p.val.ind[as.numeric(ind.class),]
              
            }	
            f.p.val<-function(x){
              if(x<alpha){return(1)}
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
            
            
            r3<-function(choix){
              sortie<-NULL
              for(i in 1:length(choix)){
                sortie<-c(sortie,which(rownames(M1_mic)%in%rownames(choix)[i]))
              }
              sortie
            }
            
            
            
            M<-new("micro_array",microarray=MM1,name=M1@name[r3(choix)],time=M1@time,subject=M1@subject,group=as.vector(R),start_time=as.vector(R))
            
            
            
            return(M)
            
          }
          
)


setMethod(f="geneSelection", 
          signature=c("list","list","numeric"),
          definition=function(x,y,tot.number,data_log=TRUE,alpha=0.05,cont=FALSE,lfc=0,f.asso=NULL){
            
            
            BBB<-strsplit(sessionInfo()[5]$otherPkgs$limma$Version,"[.]")
            
            if( !(BBB[[1]][1]>3 || (BBB[[1]][1]==3 && BBB[[1]][2]>18) || 
                    (BBB[[1]][1]==3 && BBB[[1]][2]==18 && BBB[[1]][3]>=13 ) ))
            {stop("Upgrade your version of Limma (>= 3.18.13)")}
            
            
            M<-x
            contrast<-y
            
            M_mic<-M
            require(limma)		
            n<-length(M)
            Time<-length(M[[2]]@time)
            
            if(cont==TRUE && is.null(f.asso) ){
              
              f.asso<-"mean"
              
            }
            
            Subj<-(M[[2]]@subject)
            
            
            if(data_log==TRUE){
              
              for(i in 1:n){
                M_mic[[i]]<-log(M[[i]]@microarray)
              }
            } else{
              for(i in 1:n){
                M_mic[[i]]<-(M[[i]]@microarray)
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
            Ma<-abind(M_mic,along=2)
            T<-Time
            
            
            #Construction de la matrice de design
            
            
            condition<-as.factor(paste("condition",rep(1:n,unlist(lapply(M,ncol))),sep=""))
            
            if(cont==FALSE){
              timeT<-paste("time",rep(1:T,sum(unlist(lapply(M,function(x){x@subject})))
              ),sep="") }else{
                timeT<-c(rep("time0",ncol(M[[1]]) ),paste("time",rep(1:T,sum(unlist(lapply( M[2:length(M)],function(x){x@subject})))
                ),sep=""))
                
              }
            
            Fac<-as.factor( paste(condition,timeT,sep="."))
            
            if(contrast[[1]]!="condition"){
              formule<-~-1+Fac
              design<-model.matrix(formule)
              colnames(design)<-levels(Fac)
            }else{
              formule<-~-1+condition
              design<-model.matrix(formule)
              colnames(design)<-levels(condition)
            }
            
            model<-lmFit(Ma,design)
            
            
            if(contrast[[1]]=="condition"){
              contrastM<-makeContrasts(contrasts=paste("condition",contrast[[2]][1],"-condition",
                                                       contrast[[2]][2],
                                                       sep=""),levels=design)
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
            
            if(contrast[[1]]=="patterns" && (cont==FALSE || contrast[[2]][1] != 1)){
              coma=NULL	
              for(j in 1:Time){	
                coma<-c(coma,paste("condition", contrast[[2]][1],".time", j,"-condition", contrast[[2]][2],".time", j,sep=""))
              }
              contrastM<-makeContrasts(contrasts=coma,levels=design)
              
            }
            
            if(contrast[[1]]=="patterns" && (cont==TRUE  && contrast[[2]][1] == 1)){
              coma=NULL	
              for(j in 1:Time){	
                coma<-c(coma,paste("condition", contrast[[2]][1],".time", 0,"-condition", contrast[[2]][2],".time", j,sep=""))
              }
              contrastM<-makeContrasts(contrasts=coma,levels=design)
              
            }
            
            if(contrast[[1]]=="condition&time"){
              
              
              coma<-paste("condition", contrast[[2]][1],".time", contrast[[3]][1],"-condition", contrast[[2]][2],".time", contrast[[3]][2],sep="")
              contrastM<-makeContrasts(contrasts=coma,levels=design)
              
            } 
            
            
            
            model2<-contrasts.fit(model,-contrastM)
            model.test<-eBayes(model2)
            
            p.val.ind<-model.test$p.value
            
            
            
            
            nb.tot<-length(M[[1]]@name)
            
            
            
            if(contrast[[1]]=="patterns"){
              tableT<-matrix(0,nb.tot,Time+1)	
              row.names(tableT)<-1:nb.tot
              tableT[,Time+1]<-model.test$F.p.value
              for(j in 1:Time){
                
                p.val.all<-topTable(model.test,p.value=alpha,number=nb.tot,lfc=lfc,coef=j)
                
                if(is.null(p.val.all$ID)){
                  ID<-row.names(p.val.all)
                  p.val.all<-cbind(ID,p.val.all)
                  f_nom<-function(nom){
                    which(x[[1]]@name %in% nom)			
                  }
                  f_nom<-Vectorize(f_nom)
                  row.names(p.val.all)<-f_nom(p.val.all$ID)
                }
                tableT[as.numeric(row.names(p.val.all)),j]<-1
                
              }
              
              f.test<-function(x1){
                rep<-FALSE
                for(i in 1:nrow(contrast[[3]])){
                  if(sum(x1==contrast[[3]][i,])==Time){
                    rep<-TRUE
                  }
                }
                return(rep)
              }
              
              B<-apply(tableT[,1:Time],1,f.test)
              tableT<-tableT[B,]
              tableT<-tableT[which(tableT[,Time+1]<alpha),]
              tableT<-tableT[order(tableT[,Time+1]),]
              nb.ret<-nrow(tableT)
              if(tot.number>=1 && nb.ret>tot.number){
                nb.ret<-tot.number
              }
              
              if(tot.number<1 && tot.number>0){
                nb.ret<-round(nb.ret*tot.number)
                
              }
              K1<-M_mic[[contrast[[2]][1]]][as.numeric(row.names(tableT[1:nb.ret,])),]
              K2<-M_mic[[contrast[[2]][2]]][as.numeric(row.names(tableT[1:nb.ret,])),]
              
              
              
            }else{
              
              p.val.all<-topTable(model.test,p.value=alpha,number=nb.tot,lfc=lfc)
              if(is.null(p.val.all$ID)){
                ID<-row.names(p.val.all)
                p.val.all<-cbind(ID,p.val.all)
                f_nom<-function(nom){
                  which(x[[1]]@name %in% nom)			
                }
                f_nom<-Vectorize(f_nom)
                row.names(p.val.all)<-f_nom(p.val.all$ID)
              }
              nb.ret<-dim(p.val.all)[1]
              
              if(tot.number>=1 && nb.ret>tot.number){
                nb.ret<-tot.number
              }
              
              if(tot.number<1 && tot.number>0){
                nb.ret<-round(nb.ret*tot.number)
                
              }
              
              p.val.all<-p.val.all[1:nb.ret,]
              
              
              K1<-M_mic[[contrast[[2]][1]]][p.val.all$ID,]
              K2<-M_mic[[contrast[[2]][2]]][p.val.all$ID,]
            }
            
            
            if(!is.null(f.asso) && cont==TRUE){
              K1<-apply(K1,1,f.asso)
            }
            if(!is.null(f.asso) && cont==FALSE){
              for(i in 1:Time){
                K1[i,]<-apply(K1[seq(i,ncol(K1),by=Time),],1,f.asso)
              }
            }
            if(sum(dim(K1)==dim(K2))==2){
            MM1<-K2 -K1
            }else{
              warning("The number of patient is not equal. This function returns the stimulated expression (instead of log fold change)")
              MM1<-K2
            }
            
            if(contrast[[1]]=="patterns"){
              
              f_gr<-function(x){ min(which(x==1))}
              
              group=apply(tableT[1:nb.ret,1:Time],1,f_gr)
              
              
            } else{
              group=rep(0,nb.ret)
            }
            
            M<-new("micro_array",microarray=MM1,name=row.names(MM1),time=M[[contrast[[2]][2]]]@time,subject=Subj,group=group,start_time=group)
            
            
            return(M)
            
          }
          
)



setMethod(f="genePeakSelection", 
          signature=c("micro_array","numeric"),
          definition=function(x,pic,y=NULL,data_log=TRUE,durPic=c(1,1),abs_val=TRUE,alpha_diff=0.05){
            
            M1<-x
            M2<-y
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
            
            
            
            nom1<-rownames(M1@microarray)
            nom2<-rownames(M2@microarray)
            corres<-cbind(c(M1@name,M2@name),c(nom1,nom2))
            corres<-unique(unique(corres,MARGIN=2))
            
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
            rep<-new("micro_array",microarray=M,name=corres[,1],time=M1@time,subject=M1@subject,group=gr,start_time=str)
            return(rep)
          }
          
) 


setMethod(f="unionMicro", 
          signature=c("list","ANY"),
          definition=function(M1,M2){
            
            rep<-unionMicro(M1[[1]],M1[[1]])
            if(length(M1)>1){	
              for(i in 2:length(M1)){
                rep<-unionMicro(rep,M1[[i]])
                
              }
            }
            return(rep)
          }
          
)  



setMethod(f="clustExploration", 
          signature=c("micro_array"),
          definition=function(microarray){
  
  require(Mfuzz)
  T<-length(microarray@time)
  P<-microarray@subject
  M1<-microarray@microarray
  
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
  M_exp<- ExpressionSet(M)
  
  indic<-"continue"
  k<-1
  
  while(indic=="continue"){
    k<-k+1
    cl <- mfuzz(M_exp, c = k,m=m_max)
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
  cl <- mfuzz(M_exp, c = k,m=m_max)
  if(k<5){
    mfuzz.plot(M_exp,cl=cl,mfrow=c(1,k))
  }else{
    mfuzz.plot(M_exp,cl=cl,mfrow=c(ceiling(sqrt(k)),ceiling(sqrt(k))))
  }
  
  
  return(data.frame(name=microarray@name,cluster=amv,maj.vote.index=amind))
  
  
}  
)




setMethod(f="clustInference", 
          signature=c("micro_array","numeric"),
          definition=function(microarray,vote.index){  
  require(Mfuzz)
  T<-length(microarray@time)
  P<-microarray@subject
  M1<-microarray@microarray
  
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
  M_exp<- ExpressionSet(M)
  
  indic<-"continue"
  k<-1
  within_error_hard<-NULL
  while(indic=="continue"){
    k<-k+1
    cl<-mfuzz(M_exp, c = k,m=m_max)
    we<-cl$withinerror
    for(f in 1:50){
      CL <- mfuzz(M_exp, c = k,m=m_max)
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
  cl <- mfuzz(M_exp, c = k,m=m_max)
  plot(2:(length(within_error_hard)+1),within_error_hard)
  
  
  
  
  if(k<5){
    mfuzz.plot(M_exp,cl=cl,mfrow=c(1,k))
  }else{
    mfuzz.plot(M_exp,cl=cl,mfrow=c(ceiling(sqrt(k)),ceiling(sqrt(k))))
  }
  
  
  return(list(data.frame(name=microarray@name,cluster=amv,maj.vote.index=amind),prop.matrix=prop.matrix))
  
  
}  
)