setGeneric("unsupervised_clustering_auto_m_c",package="Mfuzz",def=function(M1,... ){standardGeneric("unsupervised_clustering_auto_m_c")})
setGeneric("unsupervised_clustering",package="Mfuzz",def=function(M1,clust,mestim,... ){standardGeneric("unsupervised_clustering")})
setGeneric("geneSelection",package="Patterns",def = function(x,y,tot.number,... ){standardGeneric("geneSelection")})
setGeneric("genePeakSelection",package="Patterns",def = function(x,pic,... ){standardGeneric("genePeakSelection")})
setGeneric("unionMicro",package="Limma",def = function(M1,M2 ){standardGeneric("unionMicro")})
setGeneric("unionNGseq",package="Limma",def = function(C1,C2 ){standardGeneric("unionNGseq")})
setGeneric("position",package="igraph",def = function(net,... ){standardGeneric("position")})
setGeneric("geneNeighborhood",package="igraph",def = function(net,targets,... ){standardGeneric("geneNeighborhood")})
setGeneric("evolution",package="igraph",def = function(net,list_nv,... ){standardGeneric("evolution")})
setGeneric("inferenceCascade",package="Patterns",def = function(M,... ){standardGeneric("inferenceCascade")})
setGeneric("inference",package="Patterns",def = function(M,... ){standardGeneric("inference")})
setGeneric("cutoff",package="Patterns",def = function(Omega,... ){standardGeneric("cutoff")})
setGeneric("analyze_network",package="Patterns",def = function(Omega,nv,...){standardGeneric("analyze_network")})
setGeneric("predict",package="Patterns",def = function(object,...){standardGeneric("predict")})
setGeneric("gene_expr_simulation",def = function(network,...){standardGeneric("gene_expr_simulation")})
setGeneric("gene_counts_simulation",def = function(network,...){standardGeneric("gene_counts_simulation")})
setGeneric("clustExploration",def = function(microarray){standardGeneric("clustExploration")})
setGeneric("clustInference",def = function(microarray,vote.index){standardGeneric("clustInference")})
setGeneric("probeMerge",def = function(x,...){standardGeneric("probeMerge")})



as.micro_array<-function(M,time,subject,name_probe=NULL,gene_ID=NULL){
  if(is.null(row.names(M))){row.names(M)<-paste("gene",1:dim(M)[1])}
  if(!is.null(name_probe)){row.names(M)<-name_probe}
  g<-0
  if(!is.null(gene_ID)) g<-gene_ID
  return(new("micro_array",microarray=as.matrix(M),
             name=row.names(M),
             gene_ID=g,
             time=time,
             subject=subject,
             group=0,
             start_time=0))    
}

as.nextgen_seq<-function(C,time,subject){
  if(is.null(row.names(C))){row.names(C)<-paste("gene",1:dim(M)[1])}
  return(new("nextgen_seq",nextgenseq=as.matrix(C),name=row.names(C),time=time,subject=subject,group=0,start_time=0))	
}


lasso_reg<-function(M,Y,K,eps,priors){
  require(lars)
  cvlars<-try(cv.lars(t(M),(Y),intercept=FALSE,K=K,plot.it=FALSE,eps=10^-5,use.Gram=FALSE))
  n<-try(cvlars$index[which.min(cvlars$cv)+max(1,which.min(cvlars$cv[which.min(cvlars$cv):length(cvlars$cv)]<=min(cvlars$cv)+cvlars$cv.error[which.min(cvlars$cv)])-1)-1])
  model<-try(lars(t(M),(Y),intercept=FALSE,eps=10^-5,use.Gram=FALSE))
  repu<-try(coef(model,s=n,mode="fraction"))
  if(!is.vector(repu)){repu<-rep(0,dim(M)[1])}
  return(repu)
}


lasso_reg2<-function(M,Y,eps,foldid=foldid,priors,nfolds){
  require(glmnet)
  foldid=rep(1:nfolds,each=ncol(M)/nfolds)
  if(is.null(priors)) priors<-rep(1,nrow(M))
  cvlars<-try(cv.glmnet(t(M),Y,nfolds=nfolds,foldid=foldid,penalty.factor=priors,intercept=FALSE))
  n<-which(cvlars$cvm==min(cvlars$cvm))
  lambda<-cvlars$lambda[n]
  print(lambda)
  model<-glmnet(t(M),Y,intercept=FALSE,thresh=10^(-5),penalty.factor=priors,lambda=lambda)
  return( as.vector(coef(model)[-1]))
}

spls_reg<-function(M,Y,K,eps){
  require(spls)
  cvspls<-try(cv.spls(t(M),(Y),fold=K,K=1:10,eta = seq(0.1,0.9,0.1),plot.it=FALSE))
  model<-try(spls(t(M),(Y),eta=cvspls$eta.opt,K=cvspls$K.opt,eps=10^-5))
  repu<-try(coef(model))
  if(!is.vector(repu)){repu<-rep(0,dim(M)[1])}
  return(repu)
}


enet_reg<-function(M,Y,K,eps){
  require(elasticnet)
  M<-t(M)
  colnames(M)<-1:dim(M)[2]	  
  M2<-(M[,which(apply(M,2,sum)!=0)])
  cvenetP<-try(cv.enet((M2),(Y),lambda=0.05,s=seq(0,1,length=100),mode="fraction",intercept=FALSE,K=K,plot.it=FALSE,eps=10^-5))
  cvenet<-cvenetP
  n<-try(cvenetP$s[which.min(cvenet$cv)+max(1,which.min(cvenet$cv[which.min(cvenet$cv):length(cvenet$cv)]<=min(cvenet$cv)+cvenet$cv.error[which.min(cvenet$cv)])-1)-1])	  
  modelP<-try(enet((M2),(Y),intercept=FALSE,eps=10^-5))
  repu2 <- predict(modelP, s=n, type="coef", mode="fraction")	  
  repu<-rep(0,dim(M)[2])
  #    print(sum((apply(M,1,sum)!=0)))
  repu[as.numeric(names(repu2$coefficients))]<-repu2$coefficients
  if(!is.vector(repu)){repu<-rep(0,dim(M)[2])}  
  return(repu)
}


F_f<-function(F,x){
  
  for(i in 1:dim(F)[1]){
    F[col(F)==row(F)-(i-1)]<-x[i]
  }
  
  
  #for(i in 1:dim(F)[1]){
  #if(sum(abs(F[i,]))!=0){
  #F[i,]<-F[i,]/sum(abs(F[i,]))
  #}
  #}
  #for(i in  1:dim(F)[1]){
  #if(sum(F[i,])==0){F[i,i]<-1}
  #}
  return(F)
}

sumabso<-function(x){
  d<-sum(abs(x))
  if(d==0){
    d<-1
  }
  return(d)
}

expo<-function(x,hub){
  sum(x>=hub)/sum(x>0)
}

choice_cutoff<-function(O,nb,eps,hub,plot.g=FALSE,prop.hub=NULL){
  O<-abs(O)
  #plus petit omega sup a eps
  minO<-min(O[O>eps])
  #mediane des omegas
  qq<-quantile(O[O>eps],0.50)
  #max des omega
  maxO<-max(O)
  #cutoffs
  sequence<-seq(minO,maxO,length.out=nb)
  
  Mcut<-array(0,c(dim(O)[1],nb))
  #pour chaque valeur du cutoff on calcule 
  for(i in 1:nb){
    Mcut[,i]<-apply(O>sequence[i],1,sum)
  }
  
  expoh<-function(x){expo(x,hub)}
  hh<-apply(Mcut[,1:(dim(Mcut)[2]-1)],2,expoh)
  if(plot.g==TRUE){
    matplot(t(Mcut),type="l")
    
    if(!attr(dev.cur(),"names")=="pdf"){dev.new()}
    plot(sequence[which(hh>0)],hh[hh>0],type="l",xlab="cut off",ylab="Proportion of hubs")
  }
  
  auto<-min(which(hh==min(hh[hh>0])))
  if(!is.null(prop.hub)){
    if(prop.hub>min(hh[hh>0])){
      auto<-min(which(hh==min(hh[hh>prop.hub])))
      #abline(h=prop.hub)
    }
    
  }
  
  #print(Mcut[,auto])
  return(sequence[auto])
  
}

choice_cutoff_final<-function(O,nb,eps,hub,plot.g=FALSE,prop.hub){
  
  if(length(hub)>1){
    choix<-NULL
    for(i in hub){
      choix<-c(choix,choice_cutoff(O,nb,eps,i,prop.hub=prop.hub)	)
    }
    plot(hub,choix,type="l")
    lines(hub,predict(loess(choix ~ hub, span=0.75)),col="red")
    return(choix)
  }
  else{
    choix<-NULL
    for(i in prop.hub){
      choix<-c(choix,choice_cutoff(O,nb,eps,hub,prop.hub=i)	)
    }
    plot(prop.hub,choix,type="l")
    lines(prop.hub,predict(loess(choix ~ prop.hub, span=0.75)),col="red")
    return(choix)
  }
}

#simulations

network_random<-function(nb,time_label,exp,init,regul,min_expr,max_expr,casc.level){
  
  net<-matrix(0,nb,nb)
  net2<-net
  while(!identical(regul[which(time_label != 1)] , apply(net,2,sum)[which(time_label != 1)])) {
    for(i in 1:nb){
      if(time_label[i] !=1){
        
        if(rbinom(1,1,1-casc.level)==1){	
          reg<-which(time_label<time_label[i])}
        else{
          reg<-which(time_label==(time_label[i]-1))
        }
        if(length(reg)!=0 && sum(sum(net[,i])<regul[i])){
          pb<-apply(net[reg,],1,sum)^exp
          pb<-(pb+init)/(sum(pb+init))
          r<-rmultinom(1, 1, pb)
          net[reg[which(r==1)],i]<-1
          net2[reg[which(r==1)],i]<-runif(1,min_expr,max_expr)*(-1)^rbinom(1,1,0.5)
        }
      }
    }
  }
  length(unique(time_label))->T
  F<-CascadeFinit(T,T)
  N<-new("network",network=net2,name=paste("gene",1:nb),F=F,convO=0,convF=matrix(0,1,1),time_pt=1:length(unique(time_label)))
  return(N)
}

compare<-function(Net,Net_inf,nv){
  N1<-Net@network
  N2<-Net_inf@network
  N1[abs(N1)>0]<-1
  N1[abs(N1)<=0]<-0
  N2[abs(N2)>nv]<-1
  N2[abs(N2)<=nv]<-0
  Nb<-sum(N1)
  sens<-sum((N1-2*N2)==-1)/sum(N1==1)
  spe<-sum((N1-2*N2)==-1)/sum(N2==1)
  Fscore<-2*sens*spe/(sens+spe)
  return(c(sens,spe,Fscore))
  
}

# Band text matrices

replaceBand <- function(a, b, k) {
  swap <- abs(row(a) - col(a)) <= k
  a[swap] <- b[swap]
  a
}

replaceUp <- function(a, b, k) {
  swap <- (row(a) - col(a)) <= k
  a[swap] <- b[swap]
  a
}

replaceDown <- function(a, b, k) {
  swap <- -(row(a) - col(a)) <= k
  a[swap] <- b[swap]
  a
}

# Cascade Finit and Fshape generator

CascadeFshape=function(sqF,ngrp){
  nF=ngrp*ngrp
  Fshape<-array("0",c(sqF,sqF,nF))
  for(ii in 1:nF){   
    if((ii-1)%/%ngrp+1>=ifelse(ii%%ngrp==0,ngrp,ii%%ngrp)){
      Fshape[,,ii]<-"0"
    } else {
      lchars <- paste("a",c(4,1,2,3),sep="")
      tempFshape<-matrix("0",sqF,sqF)
      for(bb in 0:(sqF-1)){
        tempFshape<-replaceDown(tempFshape,matrix(lchars[bb+1],sqF,sqF),-bb)
      }
      tempFshape <- replaceBand(tempFshape,matrix("0",sqF,sqF),0)
      Fshape[,,ii]<-tempFshape
    }
  }
  return(Fshape)
}


CascadeFinit=function(sqF,ngrp){
  nF=ngrp*ngrp
  Finit<-array(0,c(sqF,sqF,nF))
  for(ii in 1:nF){   
    if((ii-1)%/%ngrp+1>=ifelse(ii%%ngrp==0,ngrp,ii%%ngrp)){
      Finit[,,ii]<-0
    } else {
      #       Finit[,,ii]<-lower.tri(matrix(1,4,4))*matrix(1,4,4)
      Finit[,,ii]<-cbind(rbind(rep(0,sqF-1),diag(1,sqF-1)),rep(0,sqF))
    }
  }
  return(Finit)
}

IndicFshape <- function(sqF,ngrp,Indic){
  nF=ngrp*ngrp
  Fshape<-array("0",c(sqF,sqF,nF))
  for(ii in 1:ngrp){
    for(jj in 1:ngrp){
      if(!Indic[ii,jj]){
        Fshape[,,(ii-1)*ngrp+jj]<-"0"
      } else {
        lchars <- paste("a",c(4,1,2,3),sep="")
        tempFshape<-matrix("0",sqF,sqF)
        for(bb in 0:(sqF-1)){
          tempFshape<-replaceDown(tempFshape,matrix(lchars[bb+1],sqF,sqF),-bb)
        }
        tempFshape <- replaceBand(tempFshape,matrix("0",sqF,sqF),0)
        Fshape[,,(ii-1)*ngrp+jj]<-tempFshape
      }
    }
  }
  return(Fshape)
}

IndicFinit <- function(sqF,ngrp,Indic,low.trig=TRUE){
  nF=ngrp*ngrp
  Finit<-array(0,c(sqF,sqF,nF))
  for(ii in 1:ngrp){
    for(jj in 1:ngrp){
      if(!Indic[ii,jj]){
        Finit[,,(ii-1)*ngrp+jj]<-0
      } else {
        if(low.trig){
          Finit[,,(ii-1)*ngrp+jj]<-lower.tri(matrix(1,4,4))*matrix(1,4,4)
        } else {
          Finit[,,(ii-1)*ngrp+jj]<-cbind(rbind(rep(0,sqF-1),diag(1,sqF-1)),rep(0,sqF))
        }
      }
    }
  }
  return(Finit)
}

majority_vote<-function(x){
  
  which.max(tabulate(x))
}  

majority_indice<-function(x){
  
  max(tabulate(x))
}    



boost<-function(X,Y,corr=0.8,B=100,normalize=TRUE,eps=10^(-4)){
  func<-"selec_meth"
  selec_meth<-function(X,Y){
    
    N<-dim(X)[1]      
    modd<-lars(X,Y)
    tau2<-var(Y)
    index<-which.min(N * log(2 * pi * tau2) + apply((predict(modd,newx=X)$fit-Y)^2,2,sum)/tau2 + log(N) * (apply(predict(modd,type="coef")$coef!=0,1,sum)+1))
    predict(modd,s=index,type="coef")$coefficients
    
  }
    group<-function(X,c0){
      
      Xcor<-abs(cor(X))
      Xcor[Xcor<c0]<-0
      Xcor[Xcor!=0]<-1
      
      dete<-function(x){which(x == 1)}
      res<-apply(Xcor,2,dete)
      return(res)
      
      
    }
  
  
  group2<-function(x)group(x,1)
  
  
  #require(snowfall)
  #require(lars)
  
  
  #if(multi==FALSE){
  #if(cpu>1){
  #sfInit( parallel=TRUE, cpus=cpu )
  #}else{
  #sfInit( parallel=FALSE, cpus=cpu )	
  #	}
  #sfLibrary(movMF)
  #sfLibrary(lars)
  #sfLibrary(msgps)
  
  require(movMF)
  require(lars)
  require(msgps)
  
  #sfLibrary(snowfall)}
  
  if(normalize==TRUE){
    
    X<-X-matrix(rep(apply(X,2,mean),dim(X)[1]),dim(X)[1],dim(X)[2],byrow=TRUE)	
    X<-t(t(X)/sqrt(diag(t(X)%*%X)))
    
  }
  
  Correlation_matrice<-t(X)%*%X
  Correlation_indice<-Correlation_matrice*0	
  Correlation_indice[which((Correlation_matrice)>=0)]<-1		
  Correlation_indice[which(-(Correlation_matrice)>=0)]<--1	
  diag(Correlation_indice)<-1
  
  groups<-group2(X)
  
  #nb_correlated_variables<-apply(abs(Correlation_indice),2,sum)
  #		#cat("Number of correlated variables : \n")
  
  
  Boot<-matrix(rep(0,dim(X)[2]*B),B,dim(X)[2])
  
  func_passage1<-function(x){
    
    Xpass<-matrix(0,dim(X)[1],dim(X)[1]-1)
    Xpass[col(Xpass)>row(Xpass)]<-1
    diag(Xpass)<-1
    Xpass[col(Xpass)==(row(Xpass)-1)]<--(1:(dim(X)[1]-1))
    Xpass<-t(t(Xpass)/sqrt(diag(t(Xpass)%*%Xpass)))
    return(solve(t(Xpass)%*%Xpass)%*%t(Xpass)%*%x)
  }	
  
  func_passage2<-function(x){
    
    Xpass<-matrix(0,dim(X)[1],dim(X)[1]-1)
    Xpass[col(Xpass)>row(Xpass)]<-1
    diag(Xpass)<-1
    Xpass[col(Xpass)==(row(Xpass)-1)]<--(1:(dim(X)[1]-1))
    Xpass<-t(t(Xpass)/sqrt(diag(t(Xpass)%*%Xpass)))
    return(apply( matrix(rep(x,dim(X)[1]),dim(X)[1],dim(X)[1]-1,byrow=TRUE)*Xpass,1,sum))
  }
  
  Xret<-array(0,c(dim(X)[1],dim(X)[2],B))
  Xb<-X
  
  
  
  simul1<-function(j){
    
    if(length(groups[[j]])>=2){
      
      indice<-groups[[j]]
      corr_set<-t(t(X[,indice])*Correlation_indice[j,indice])
      corr_set2<-apply(corr_set,2,func_passage1)
      BB<-coef(movMF(t(corr_set2),1))$theta
      newv<-rmovMF(1,BB)
      newv<-func_passage2(newv)
    }else{
      
      newv<-X[,j]
      
    }
    
    
    return(newv)			
  }
  
  simul1<-Vectorize(simul1)
  
  simul2<-function(){simul1(1:(dim(X)[2]))}
  
  Xret<-replicate(B,simul2())
  
  
  select<-function(k){
    rr<-eval(parse(text=paste(func,"(Xret[,,",k,"],Y)",sep="")))
    
    return(rr)	
  }
  
  
  Boot<-lapply(1:B,select)
  
  Boot<-matrix(unlist(Boot),B,dim(X)[2],byrow=TRUE)
  
  Fs<-apply(abs(Boot)>eps,2,sum)
  
  
  return(Fs/B)
}
