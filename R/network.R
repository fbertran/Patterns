setClass(Class = "network",
representation(network="matrix",name="vector",F="array",convF="matrix",convO="vector",time_pt="vector")
)



setMethod("print","network",function(x,...){
	cat(paste("This is a S4 class with : \n - (@network) a matrix of dimension ",dim(x@network)[1],"*",dim(x@network)[2]," .... [the network] \n - (@name) a vector of length ",length(x@name)," .... [gene names] \n","- (@F) a array of dimension ",dim(x@F)[1],"*",dim(x@F)[2],"*",dim(x@F)[3]," .... [F matrices] \n","- (@convF) a matrix of dimension ",dim(x@convF)[1],"*",dim(x@convF)[2]," .... [convergence (L1 norm) of array F] \n","- (@convO)a vector of length ",length(x@convO)," .... [convergence (L1 norm) of matrix Omega]\n","- (@time_pt) an vector of length",length(x@time_pt),"  .... [time points]")) 
})



setMethod("analyze_network","network",function(Omega,nv){
	require(tnet)
	O<-Omega@network
	Omega<-Omega@network
    O[abs(O)<=nv]<-0
   	G<-graph.adjacency(O,weighted=TRUE)
   	if(.Platform$OS.type=="unix"){
    get.edgelist(G)+1->Q
    }else{
    get.edgelist(G)->Q	
    	}
    weight<-rep(0,dim(Q)[1])
	for(i in 1:dim(Q)[1]){weight[i]<-Omega[Q[i,1],Q[i,2]]}
	R<-cbind(Q,abs(weight))
   	G<-as.tnet(R)
   	C<-cbind(betweenness_w(G)[,1:2], degree_w(G)[,2:3],closeness_w(G,gconly=FALSE)[,2])
	colnames(C)<-c("node","betweenness","degree","output","closeness")
	return(C)
}
)



setMethod("evolution","network",function(net,list_nv,gr=NULL,color.vertex=NULL,fix=TRUE,gif=TRUE,taille=c(2000,1000)){

Omega<-net
	
	if(gif==TRUE){
		require(animation)
		ani.options(ani.height=taille[2],ani.width=taille[1],outdir = getwd())
		saveHTML({
			
			POS<-position(net,nv=list_nv[1])
		
	for(i in list_nv){
		
		if(fix==TRUE){
		plot(Omega,nv=i,gr=gr,ini=POS,color.vertex=color.vertex)
		}
		else{
		 	plot(Omega,nv=i,gr=gr,color.vertex=color.vertex)
		
		}
    text(-1,1,round(i,3))
		
	}
			
		
		})
		
		
	}
	
	else{
		
		
		par(ask=TRUE)
	POS<-position(net,nv=list_nv[1])
	if(is.null(gr)){gr<-rep(1,dim(Omega@network)[1])} 
	for(i in list_nv){
		plot(Omega,nv=i,gr=gr,ini=POS,color.vertex=color.vertex)
		
	}

	par(ask=FALSE)	
		
		}
	
	}
	
	
)

setMethod("position","network",function(net,nv=0){
	require(igraph)
    O<-net@network
  	Omega<-net@network
    O[abs(O)<=nv]<-0
    O[abs(O)>nv]<-1
      
        nom<-1:dim(O)[1]
    enle<-which(apply(O,1,sum)+apply(O,2,sum)==0)
 	if(length(enle)!=0){
 	nom<-nom[-enle]
 	O<-O[-enle,-enle]
 	 
 	 }
       
    G<-graph.adjacency(O,weighted=TRUE)
    	if(.Platform$OS.type=="unix"){
    get.edgelist(G)+1->Q
    }else{
    get.edgelist(G)->Q	
    	}
    Q[,1]<-nom[Q[,1]]
    Q[,2]<-nom[Q[,2]]
      
   # dev.new(width=150,heigth=300,xlim=c(min(L[,1]),max(L[,1])))
    #par(mar=c(0,0,0,0), oma=c(0,0,0,0),mai=c(0,0,0,0))
	L<-layout.fruchterman.reingold(G,xmin=0,xmax=150,ymin=0,ymax=300)
	return(cbind(nom,L))
}

)




setMethod("plot","network",function(x,y,choice="network",nv=0,gr=NULL,ini=NULL,color.vertex=NULL,weight.node=NULL,ani=FALSE,taille=c(2000,1000),label_v=1:dim(x@network)[1],horiz=TRUE,legend.position="topleft",frame.color="black",label.hub=FALSE,nround=2,...)
{
	
if(choice=="F"){
  library(plotrix)
	F<-x@F
  nF<-dim(F)[3]
  ngrp=sqrt(dim(F)[3])
  ymax<-max(F)
  couleur<-rainbow(ngrp)
#  par(mfrow=c(ngrp,ngrp))
FF=NULL
for(i in 1:ngrp){
FFa=NULL
for(j in 1:ngrp){
FFa=cbind(FFa,round(F[,,(i-1)*ngrp+j],nround))
#par(mai=c(0.1,0.1,0.1,0.1))
#color2D.matplot(x=round(F[,,(i-1)*ngrp+j],2),cs1=c(0,1,1),cs2=c(1,1,0),cs3=0,show.values=TRUE,axes=FALSE,main="",xlab="",ylab="",show.legend=TRUE)
}
FF=rbind(FF,FFa)
}
par(mar=c(0,0,0,0), oma=c(0,0,0,0),mai=c(0,0,0,0))
color2D.matplot(x=FF,cs1=c(0,.5,1),cs2=c(.5,1,0),cs3=c(1,2,0),show.values=nround,axes=FALSE,main="",xlab="",ylab="",show.legend=FALSE)
abline(h=4*(1:(ngrp-1)),lwd=3)
abline(v=4*(1:(ngrp-1)),lwd=3)
}
if(choice=="network"){
	
	couleur<-0	
	O<-x@network
  	Omega<-x@network
    O[abs(O)<=nv]<-0
    O[abs(O)>nv]<-1
	F<-x@F
  nF<-dim(F)[3]
  ngrp=sqrt(dim(F)[3])
    
	
	
	require(igraph)
	if(is.null(gr)){gr<-rep(1,dim(O)[1])}
	
	if(is.null(color.vertex)){
    	  	color.vertex<-rainbow(length(unique(gr)))
    	  	couleur<-1
    	}
    	  
	
	
	nom<-1:dim(O)[1]
    enle<-which(apply(O,1,sum)+apply(O,2,sum)==0)
 	if(length(enle)!=0){
 	nom<-nom[-enle]
 	
 	O<-O[-enle,-enle]
 	 
 	 if(!is.null(weight.node)){weight.node<-weight.node[-enle]}
 	}
 	
 	if(!is.null(ini)){
 	ini<-ini[which(ini[,1] %in% nom),]
    }
    nom2<-nom
 	nom2<-label_v[nom]
    if(!is.null(weight.node)){
    	size<-weight.node
    }    else{
    	size<-apply(O,1,sum)
    	}
   
    size<-size/max(size)*10
    size[size<2]<-2 
   	
   	if(label.hub==TRUE){nom2[which(size<3)]<-NA}
   	
   	  
   	G<-graph.adjacency(O,weighted=TRUE)
   	
    	if(.Platform$OS.type=="unix"){
    get.edgelist(G)+1->Q
    }else{
    get.edgelist(G)->Q	
    	}
    Q[,1]<-nom[Q[,1]]
    Q[,2]<-nom[Q[,2]]
    size.ed<-rep(0,dim(Q)[1])
  
    for(i in 1:dim(Q)[1]){
    size.ed[i]<-abs(Omega[Q[i,1],Q[i,2]])
    }
    fpl<-function(x){3*exp(8*x)/(14+exp(8*x))}
    
    size.ed<-fpl(size.ed)
    size.ed<-size.ed

	

	
		if(ani==TRUE){
			if(is.null(gr)){error("Need of groups")}
			T<-length(x@time_pt)
			require(animation)
			ani.options(ani.height=taille[2],ani.width=taille[1],outdir = getwd())
			
			if(is.null(ini)){
					
  					ini<-position(x,nv)[]
				}
			
			
				#ini<-ini[which(ini[,1] %in% nom),]
				ini<-ini[,2:3]
			gr2<-gr[nom]
			
			saveHTML({
				
					
										
 for(i in c(0.999,0.99,0.98,0.97,0.96,0.95)){											plot(G,layout=ini,vertex.size=size,edge.width=size.ed,vertex.color=grey(i),vertex.label.cex=1,edge.color=grey(i),edge.arrow.size=0.2,asp=0,vertex.label=nom2,vertex.frame.color=frame.color)
				}
for(i in 1:(ngrp)){
	
Ver.col<-rep(grey(0.95),length(nom))	
Edge.col<-rep(grey(0.95),dim(Q)[1])
mms<-20
PP<-1:mms
for(p in PP){	
	
	
coul<-col2rgb(color.vertex[i])
indi<-which(gr2==i)
Ver.col[indi]<-rgb(coul[1],coul[2],coul[3],alpha=(p-1)*255/(mms), max = 255)


for(k in 1:ngrp){
coul<-col2rgb(color.vertex[k])	
indi3<-which(gr[Q[,1]]==k & gr[Q[,2]]==i )

Edge.col[indi3]<-rgb(coul[1],coul[2],coul[3],alpha=255-(p-1)*255/(mms), max = 255)
}

if(i!=ngrp){
for(k in 1:ngrp){
	for(k2 in (1:ngrp)[-k]){
coul<-col2rgb(color.vertex[k])	
indi3<-which(gr[Q[,1]]==k & gr[Q[,2]]==k2 )




al<-round(mms/(k2-k))
ald<-min(1+al*(i-k),mms)
if(k2==(i+1)){
	alf<-mms
}
else{
alf<-min(al*(i-k+1),mms)}
alseq<-seq(25,255,length.out=mms)

alphaseq2<-seq(alseq[ald],alseq[alf],length.out=mms)

Edge.col[indi3]<-rgb(coul[1],coul[2],coul[3],alpha=alphaseq2[p], max = 255)
}}}

plot(G,layout=ini,vertex.size=size,edge.width=size.ed,vertex.color=Ver.col,vertex.label.cex=1,edge.color=Edge.col,edge.arrow.size=0.2,asp=0,vertex.label=nom2,vertex.frame.color=frame.color)
if(length(unique(gr[nom])) >1){
legend(legend.position,horiz=horiz,pch=1,col=color.vertex[unique(gr[nom])[order(unique(gr[nom]))]],legend=paste("Cluster",unique(gr[nom])[order(unique(gr[nom]))]))} 
	}
	

	
}
				
				
				})
		}	
		else{
	
	
    	par(mar=c(0,0,0,0), oma=c(0,0,0,0),mai=c(0,0,0,0))
    
       

    

  if(is.null(ini)){
  	L<-position(x,nv)[,2:3]
  
 if(couleur==0){plot(G,layout=L,vertex.size=size,edge.width=size.ed,vertex.color=color.vertex[nom],vertex.label.cex=1,edge.color=color.vertex[gr[Q[,2]]],edge.arrow.size=0.2,asp=0,vertex.label=nom2,vertex.frame.color=frame.color)}  
 else{  plot(G,layout=L,vertex.size=size,edge.width=size.ed,vertex.color=color.vertex[gr[nom]],vertex.label.cex=1,edge.color=color.vertex[gr[Q[,2]]],edge.arrow.size=0.2,asp=0,vertex.label=nom2,vertex.frame.color=frame.color)}
 
 if(length(unique(gr[nom])) >1){
legend(legend.position,horiz=horiz,pch=1,col=color.vertex[unique(gr[nom])[order(unique(gr[nom]))]],legend=paste("Cluster",unique(gr[nom])[order(unique(gr[nom]))])) }
	}
	else{
		
		L<-ini[,2:3]
if(couleur==0){plot(G,layout=L,vertex.size=size,edge.width=size.ed,vertex.color=color.vertex[nom],vertex.label.cex=1,edge.color=color.vertex[gr[Q[,2]]],edge.arrow.size=0.2,asp=0,vertex.label=nom2,vertex.frame.color=frame.color)}  
 else{  plot(G,layout=L,vertex.size=size,edge.width=size.ed,vertex.color=color.vertex[gr[nom]],vertex.label.cex=1,edge.color=color.vertex[gr[Q[,2]]],edge.arrow.size=0.2,asp=0,vertex.label=nom2,vertex.frame.color=frame.color)}
 
if(length(unique(gr[nom])) >1){
		legend(legend.position,horiz=horiz,pch=1,col=color.vertex[unique(gr[nom])[order(unique(gr[nom]))]],legend=paste("Cluster",unique(gr[nom])[order(unique(gr[nom]))])) }
	
		L<-ini
		}
	
	}

}	
})

setMethod("geneNeighborhood","network",function(net,targets,nv=0,order=length(net@time_pt)-1,label_v=NULL,ini=NULL,frame.color="white",label.hub=FALSE){
	require(igraph)
	O<-net@network
  	Omega<-net@network
    O[abs(O)<=nv]<-0
    O[abs(O)>nv]<-1
	G<-graph.adjacency(O,weighted=TRUE)
	couleur<-colorRamp(c("red","blue"))
	color<-rep(grey(0.92),length(V(G)))
	if(is.null(label_v)){label_v<-1:dim(O)[1]}
	K<-vector("list",order)
	for(j in order:1){
	
		if(.Platform$OS.type=="unix"){
    N<-neighborhood(G, order=j, nodes=targets-1, mode=c( "out"))
    }else{
    N<-neighborhood(G, order=j, nodes=targets, mode=c( "out"))
    	}
	
	
	K[[j]]<-N
	gg<-length(targets)
		
	for(i in 1:gg){
		if(.Platform$OS.type=="unix"){
    color[N[[i]]+1]<-rgb(couleur((j-1)*1/(order-1)),alpha=255,max=255)
    }else{
    color[N[[i]]+1]<-rgb(couleur((j-1)*1/(order-1)),alpha=255,max=255)
    	}
		
	}
	}
	color[targets]<-"green"
	plot(net,nv=nv,ini=ini,color.vertex=color,frame.color=frame.color,label.hub=label.hub,label_v=label_v)
	legend("topright",legend=paste("order",1:order),pch=1,col=rgb(couleur((1:order-1)*1/(order-1)),max=255))
	return(K)
	}
	)
	

setMethod("cutoff","network",function(Omega,sequence=NULL,x_min=NULL){
	
plfit<-function(x=rpareto(1000,10,2.5),method="limit",value=c(),finite=FALSE,nowarn=FALSE,nosmall=FALSE){
      #init method value to NULL	
      vec <- c() ; sampl <- c() ; limit <- c(); fixed <- c()
   ###########################################################################################
   #
   #  test and trap for bad input
   #
      switch(method, 
        range = vec <- value,
        sample = sampl <- value,
        limit = limit <- value,
        fixed = fixed <- value,
        argok <- 0)
      
      if(exists("argok")){stop("(plfit) Unrecognized method")}
   
      if( !is.null(vec) && (!is.vector(vec) || min(vec)<=1 || length(vec)<=1) ){
        print(paste("(plfit) Error: ''range'' argument must contain a vector > 1; using default."))
        vec <- c()
      }
      if( !is.null(sampl) && ( !(sampl==floor(sampl)) ||  length(sampl)>1 || sampl<2 ) ){
        print(paste("(plfit) Error: ''sample'' argument must be a positive integer > 2; using default."))
        sample <- c()
      }
      if( !is.null(limit) && (length(limit)>1 || limit<1) ){
        print(paste("(plfit) Error: ''limit'' argument must be a positive >=1; using default."))
        limit <- c()
      }
      if( !is.null(fixed) && (length(fixed)>1 || fixed<=0) ){
        print(paste("(plfit) Error: ''fixed'' argument must be a positive >0; using default."))
        fixed <- c() 
      }   
   
   #  select method (discrete or continuous) for fitting and test if x is a vector
      fdattype<-"unknow"
      if( is.vector(x,"numeric") ){ fdattype<-"real" }
      if( all(x==floor(x)) && is.vector(x) ){ fdattype<-"integer" }
      if( all(x==floor(x)) && min(x) > 1000 && length(x) > 100 ){ fdattype <- "real" }
      if( fdattype=="unknow" ){ stop("(plfit) Error: x must contain only reals or only integers.") }
   
   #
   #  end test and trap for bad input
   #
   ###########################################################################################
   
   ###########################################################################################
   #
   #  estimate xmin and alpha in the continuous case
   #
      if( fdattype=="real" ){
   
       xmins <- sort(unique(x))
       xmins <- xmins[-length(xmins)]
   
       if( !is.null(limit) ){
         xmins <- xmins[xmins>=limit]
       } 
       if( !is.null(fixed) ){
         xmins <- fixed
       }
       if( !is.null(sampl) ){ 
         xmins <- xmins[unique(round(seq(1,length(xmins),length.out=sampl)))]
       }
    
       dat <- rep(0,length(xmins))
       z   <- sort(x)
       
       for( xm in 1:length(xmins) ){
         xmin <- xmins[xm]
         z    <- z[z>=xmin]
         n    <- length(z)
         # estimate alpha using direct MLE
         a    <- n/sum(log(z/xmin))
         # truncate search if nosmall is selected
         if( nosmall ){
          if((a-1)/sqrt(n) > 0.1){
              dat <- dat[1:(xm-1)]
              print(paste("(plfit) Warning : xmin search truncated beyond",xmins[xm-1]))
              break
          }
         }
         # compute KS statistic
         cx   <- c(0:(n-1))/n
         cf   <- 1-(xmin/z)^a
         dat[xm] <- max(abs(cf-cx))
       }
    
       D     <- min(dat)
       xmin  <- xmins[min(which(dat<=D))]
       z     <- x[x>=xmin]
       n     <- length(z)
       alpha <- 1 + n/sum(log(z/xmin))
   
       if( finite ){
         alpha <- alpha*(n-1)/n+1/n # finite-size correction
       }
      
   
      }
   #
   #  end continuous case
   #
   ###########################################################################################
   
   ###########################################################################################
   #
   #  estimate xmin and alpha in the discrete case
   #
      if( fdattype=="integer" ){
   
       if( is.null(vec) ){ vec<-seq(1.5,3.5,.01) } # covers range of most practical scaling parameters
       zvec <- zeta(vec)
   
       xmins <- sort(unique(x))
       xmins <- xmins[-length(xmins)]
   
       if( !is.null(limit) ){
         limit <- round(limit)
         xmins <- xmins[xmins>=limit]
       } 
   
       if( !is.null(fixed) ){
         xmins <- fixed
       }
   
       if( !is.null(sampl) ){ 
         xmins <- xmins[unique(round(seq(1,length(xmins),length.out=sampl)))]
       }
   
       if( is.null(xmins) || length(xmins) < 2){
         stop("(plfit) error: x must contain at least two unique values.")
       }
   
       if(length(which(xmins==0) > 0)){
         stop("(plfit) error: x must not contain the value 0.")
       }
   
       xmax <- max(x)
        dat <- matrix(0,nrow=length(xmins),ncol=2)
          z <- x
       for( xm in 1:length(xmins) ){
         xmin <- xmins[xm]
         z    <- z[z>=xmin]
         n    <- length(z)
         # estimate alpha via direct maximization of likelihood function
         #  vectorized version of numerical calculation
         # matlab: zdiff = sum( repmat((1:xmin-1)',1,length(vec)).^-repmat(vec,xmin-1,1) ,1);
         if(xmin==1){
           zdiff <- rep(0,length(vec))
         }else{
          zdiff <- apply(rep(t(1:(xmin-1)),length(vec))^-t(kronecker(t(array(1,xmin-1)),vec)),2,sum)
         }
         # matlab: L = -vec.*sum(log(z)) - n.*log(zvec - zdiff);
         L <- -vec*sum(log(z)) - n*log(zvec - zdiff);
         I <- which.max(L)
         # compute KS statistic
         fit <- cumsum((((xmin:xmax)^-vec[I])) / (zvec[I] - sum((1:(xmin-1))^-vec[I])))
         cdi <- cumsum(hist(z,c(min(z)-1,(xmin+.5):xmax,max(z)+1),plot=FALSE)$counts/n)
         dat[xm,] <- c(max(abs( fit - cdi )),vec[I])
       }
       D     <- min(dat[,1])
       I     <- which.min(dat[,1])
       xmin  <- xmins[I]
       n     <- sum(x>=xmin)
       alpha <- dat[I,2]
   
       if( finite ){
         alpha <- alpha*(n-1)/n+1/n # finite-size correction
       }
     
      }
   #
   #  end discrete case
   #
   ###########################################################################################
   
   #  return xmin, alpha and D in a list
      return(list(xmin=xmin,alpha=alpha,D=D))
   }


sequence_test<-sequence

if(is.null(x_min)){
	x_min<-round(dim(Omega@network)[1]*0.02)
	
}
if(is.null(sequence_test)){
	sequence_test<-seq(0,min(max(abs(Omega@network*0.9)),0.4),length.out=200)
	}
	
#install_github('poweRlaw', 'csgillespie')
	
cc<-0
u<-0
f<-rep(0,length(sequence_test))
for(cc in sequence_test){
	u<-u+1
O<-Omega@network

 O[abs(O)>cc]<-1
  O[abs(O)<=cc]<-0
  liste<-apply(O,1,sum)+1
  if(length(which(liste>round(x_min)+1))<4){break}

xmin = plfit(liste[which(liste!=1)],value=1)$D
f[u]<-xmin
}

p<-which(sequence_test>0.25)[1]

f<-f[1:p]
K<-1:length(f)
f2<-f
f<-predict(loess(f~K))
par(mar=c(5,4,0.2,0.2))
plot(sequence_test[1:p],f,xlab="cutoff sequence",ylab="scale-freeness test p-value",type="l",lwd=2,ylim=c(min(f)-0.02,max(f)))
abline(h=0.1,col="red")
rect(0.10,0.10,0.15,1,col=rgb(0,1,0,0.5))
rect(0.05,0.10,0.10,1,col=rgb(0.5,1,0,0.4))
rect(-0.5,0.10,0.05,1,col=rgb(1,0,0,0.4))
rect(0.15,0.10,0.5,1,col=rgb(0.5,1,0,0.4))
abline(h=0.1,col="red",lwd=3)
legend("bottomright",lty=c(1,0,0,0),col=c("red",rgb(0,1,0,0.5),rgb(0.5,1,0,0.4),
rgb(1,0,0,0.4)),pch=c(NA,15,15,15),legend=c("p-value=0.1:  above this line, scale-freeness may be assumed","best area of choice (determined by simulation)","less recommended area of choice (determined by simulation)","area of choice to be avoided (determined by simulation)"),
lwd=c(3,5,5,5),cex=1)

return(list(p.value=f2,p.value.inter=f,sequence=sequence_test[1:p])	)
	
}

)	



	
#setMethod("cutoff","network",function(Omega,sequence=NULL,x_min=NULL){
#	
#plfit<-function(x=rpareto(1000,10,2.5),method="limit",value=c(),finite=FALSE,nowarn=FALSE,nosmall=FALSE){
#      #init method value to NULL	
#      vec <- c() ; sampl <- c() ; limit <- c(); fixed <- c()
#   ###########################################################################################
#   #
#   #  test and trap for bad input
#   #
#      switch(method, 
#        range = vec <- value,
#        sample = sampl <- value,
#        limit = limit <- value,
#        fixed = fixed <- value,
#        argok <- 0)
#      
#      if(exists("argok")){stop("(plfit) Unrecognized method")}
#   
#      if( !is.null(vec) && (!is.vector(vec) || min(vec)<=1 || length(vec)<=1) ){
#        print(paste("(plfit) Error: ''range'' argument must contain a vector > 1; using default."))
#        vec <- c()
#      }
#      if( !is.null(sampl) && ( !(sampl==floor(sampl)) ||  length(sampl)>1 || sampl<2 ) ){
#        print(paste("(plfit) Error: ''sample'' argument must be a positive integer > 2; using default."))
#        sample <- c()
#      }
#      if( !is.null(limit) && (length(limit)>1 || limit<1) ){
#        print(paste("(plfit) Error: ''limit'' argument must be a positive >=1; using default."))
#        limit <- c()
#      }
#      if( !is.null(fixed) && (length(fixed)>1 || fixed<=0) ){
#        print(paste("(plfit) Error: ''fixed'' argument must be a positive >0; using default."))
#        fixed <- c() 
#      }   
#   
#   #  select method (discrete or continuous) for fitting and test if x is a vector
#      fdattype<-"unknow"
#      if( is.vector(x,"numeric") ){ fdattype<-"real" }
#      if( all(x==floor(x)) && is.vector(x) ){ fdattype<-"integer" }
#      if( all(x==floor(x)) && min(x) > 1000 && length(x) > 100 ){ fdattype <- "real" }
#      if( fdattype=="unknow" ){ stop("(plfit) Error: x must contain only reals or only integers.") }
#   
#   #
#   #  end test and trap for bad input
#   #
#   ###########################################################################################
#   
#   ###########################################################################################
#   #
#   #  estimate xmin and alpha in the continuous case
#   #
#      if( fdattype=="real" ){
#   
#       xmins <- sort(unique(x))
#       xmins <- xmins[-length(xmins)]
#   
#       if( !is.null(limit) ){
#         xmins <- xmins[xmins>=limit]
#       } 
#       if( !is.null(fixed) ){
#         xmins <- fixed
#       }
#       if( !is.null(sampl) ){ 
#         xmins <- xmins[unique(round(seq(1,length(xmins),length.out=sampl)))]
#       }
#    
#       dat <- rep(0,length(xmins))
#       z   <- sort(x)
#       
#       for( xm in 1:length(xmins) ){
#         xmin <- xmins[xm]
#         z    <- z[z>=xmin]
#         n    <- length(z)
#         # estimate alpha using direct MLE
#         a    <- n/sum(log(z/xmin))
#         # truncate search if nosmall is selected
#         if( nosmall ){
#          if((a-1)/sqrt(n) > 0.1){
#              dat <- dat[1:(xm-1)]
#              print(paste("(plfit) Warning : xmin search truncated beyond",xmins[xm-1]))
#              break
#          }
#         }
#         # compute KS statistic
#         cx   <- c(0:(n-1))/n
#         cf   <- 1-(xmin/z)^a
#         dat[xm] <- max(abs(cf-cx))
#       }
#    
#       D     <- min(dat)
#       xmin  <- xmins[min(which(dat<=D))]
#       z     <- x[x>=xmin]
#       n     <- length(z)
#       alpha <- 1 + n/sum(log(z/xmin))
#   
#       if( finite ){
#         alpha <- alpha*(n-1)/n+1/n # finite-size correction
#       }
#      
#   
#      }
#   #
#   #  end continuous case
#   #
#   ###########################################################################################
#   
#   ###########################################################################################
#   #
#   #  estimate xmin and alpha in the discrete case
#   #
#      if( fdattype=="integer" ){
#   
#       if( is.null(vec) ){ vec<-seq(1.5,3.5,.01) } # covers range of most practical scaling parameters
#       zvec <- zeta(vec)
#   
#       xmins <- sort(unique(x))
#       xmins <- xmins[-length(xmins)]
#   
#       if( !is.null(limit) ){
#         limit <- round(limit)
#         xmins <- xmins[xmins>=limit]
#       } 
#   
#       if( !is.null(fixed) ){
#         xmins <- fixed
#       }
#   
#       if( !is.null(sampl) ){ 
#         xmins <- xmins[unique(round(seq(1,length(xmins),length.out=sampl)))]
#       }
#   
#       if( is.null(xmins) || length(xmins) < 2){
#         stop("(plfit) error: x must contain at least two unique values.")
#       }
#   
#       if(length(which(xmins==0) > 0)){
#         stop("(plfit) error: x must not contain the value 0.")
#       }
#   
#       xmax <- max(x)
#        dat <- matrix(0,nrow=length(xmins),ncol=2)
#          z <- x
#       for( xm in 1:length(xmins) ){
#         xmin <- xmins[xm]
#         z    <- z[z>=xmin]
#         n    <- length(z)
#         # estimate alpha via direct maximization of likelihood function
#         #  vectorized version of numerical calculation
#         # matlab: zdiff = sum( repmat((1:xmin-1)',1,length(vec)).^-repmat(vec,xmin-1,1) ,1);
#         if(xmin==1){
#           zdiff <- rep(0,length(vec))
#         }else{
#          zdiff <- apply(rep(t(1:(xmin-1)),length(vec))^-t(kronecker(t(array(1,xmin-1)),vec)),2,sum)
#         }
#         # matlab: L = -vec.*sum(log(z)) - n.*log(zvec - zdiff);
#         L <- -vec*sum(log(z)) - n*log(zvec - zdiff);
#         I <- which.max(L)
#         # compute KS statistic
#         fit <- cumsum((((xmin:xmax)^-vec[I])) / (zvec[I] - sum((1:(xmin-1))^-vec[I])))
#         cdi <- cumsum(hist(z,c(min(z)-1,(xmin+.5):xmax,max(z)+1),plot=FALSE)$counts/n)
#         dat[xm,] <- c(max(abs( fit - cdi )),vec[I])
#       }
#       D     <- min(dat[,1])
#       I     <- which.min(dat[,1])
#       xmin  <- xmins[I]
#       n     <- sum(x>=xmin)
#       alpha <- dat[I,2]
#   
#       if( finite ){
#         alpha <- alpha*(n-1)/n+1/n # finite-size correction
#       }
#     
#      }
#   #
#   #  end discrete case
#   #
#   ###########################################################################################
#   
#   #  return xmin, alpha and D in a list
#      return(list(xmin=xmin,alpha=alpha,D=D))
#   }
#
#
#sequence_test<-sequence
#
#if(is.null(x_min)){
#	x_min<-round(dim(Omega@network)[1]*0.02)
#	
#}
#if(is.null(sequence_test)){
#	sequence_test<-seq(0,min(max(abs(Omega@network*0.9)),0.4),length.out=2000)
#	}
#	
##install_github('poweRlaw', 'csgillespie')
#	
#cc<-0
#u<-0
#f<-rep(0,length(sequence_test))
#for(cc in sequence_test){
#	u<-u+1
#O<-Omega@network
#
# O[abs(O)>cc]<-1
#  O[abs(O)<=cc]<-0
#  liste<-apply(O,1,sum)+1
#  if(length(which(liste>round(x_min)+1))<4){break}
#
#xmin = plfit(liste[which(liste!=1)],value=1)$D
#f[u]<-xmin
#}
#
#f<-f[1:u]
#K<-1:length(f)
#f2<-f
#f<-predict(loess(f~K))
#plot(sequence_test[1:u],f,xlab="cutoff sequence",ylab="p value")
#return(list(p.value=f2,p.value.inter=f,sequence=sequence_test[1:u])	)
#	
#}
#
#)	
#


