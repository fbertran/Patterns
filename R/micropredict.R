setClass(Class = "micropredict",
representation(microarray_unchanged="micro_array",microarray_changed="micro_array",microarray_predict="micro_array",nv="numeric",network="network",targets="numeric")	)

setMethod("plot",c("micropredict"),function(x,time=NULL,label_v=NULL,frame.color="white",ini=NULL,label.hub=FALSE){
	net<-x@network
	micro<-predict(x@microarray_unchanged,x@network,nv=x@nv)@microarray_predict
	nv<-x@nv
	micro_pred<-x@microarray_predict
	targets<-x@targets
	O<-net@network
  	Omega<-net@network
    O[abs(O)<=nv]<-0
    O[abs(O)>nv]<-1
    if(is.null(label_v)){label_v<-1:dim(O)[1]}
	
	G<-graph.adjacency(O,weighted=TRUE)
	
	couleur1<-colorRamp(c(grey(0.95),"blue"))
	couleur2<-colorRamp(c(grey(0.95),"red"))
	color<-rep(grey(0.0999),length(V(G)))
	
	if(is.null(ini)){
	P<-position(net,nv)
	}else{
		P<-ini
		}
		
	if(is.null(time)){
		K<-2:(length(micro@time))
		}
		else{K<-time}
		for(time in K){
	if(!attr(dev.cur(),"names")=="pdf"){dev.new()}
		if(time==length(micro@time)){time<-0}
		sup_pred<-1:dim(micro@microarray)[2]
	sup_pred<-sup_pred[sup_pred%%length(micro@time)==time]
	M<-apply(micro@microarray[,sup_pred]-micro_pred@microarray[,sup_pred],1,mean)
	
	maxi<-round(max(abs(M[-targets])),2)
	if(max(abs(M))!=0){M<-M/max(abs(M[-targets]))}
	long<-1:length(M)
	long<-long[-targets]
	for(j in long){
		
		if(M[j]>=0){
			color[j]<-rgb(couleur2(M[j]),alpha=255,max=255)
					}
					else{
						
						color[j]<-rgb(couleur1(-M[j]),alpha=255,max=255)

					}
		
	}
	if(!is.null(targets)){color[targets]<-"green"}
	plot(net,nv=nv,color.vertex=color,ini=P,label_v=label_v,frame.color=frame.color,label.hub=label.hub)
couleur3<-colorRamp(c("blue","grey","red"))
nb.col<-300
coll<-rgb(couleur3(1:nb.col/nb.col),max=255)

rect(seq(-1,-0.5,length.out=nb.col)[1:(nb.col-1)],-1,seq(-1,-0.5,length.out=nb.col)[2:nb.col],-0.95,border="transparent",col=coll)
text(-1,-1.05,-maxi)
text(-0.51,-1.05,maxi)
text(-0.751,-1.05,"0")	
if(time==0){time<-length(micro@time)}
 text(-0.8,1.1,paste("Time point prediction =",time, sep=" "))
		}
		
}
)



