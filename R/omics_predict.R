#' @rdname plot-methods
setMethod("plot"
          ,c("omics_predict")
          ,function(x
                    ,time=NULL
                    ,label_v=NULL
                    ,frame.color="white"
                    ,ini=NULL
                    ,label.hub=FALSE
                    ,edge.arrow.size=0.7
                    ,edge.thickness=1){
            net<-x@omics_network
            omics<-x@omicsarray_changed
            nv<-x@nv
            omics_pred<-x@omicsarray_predict
            targets<-x@targets
            O<-net@omics_network
            Omega<-net@omics_network
            O[abs(O)<nv]<-0
            O[abs(O)>=nv]<-1
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
              K<-2:(length(omics@time))
            }
            else{K<-time}
            for(time in K){
              #if(length(K)>1){
                #if(!attr(dev.cur(),"names")=="pdf"){dev.new()}
              #}
              if(time==length(omics@time)){time<-0}
              sup_pred<-1:dim(omics@omicsarray)[2]
              sup_pred<-sup_pred[sup_pred%%length(omics@time)==time]
              M<-apply(omics@omicsarray[,sup_pred]-omics_pred@omicsarray[,sup_pred],1,mean)
              
              maxi<-round(max(abs(M[-targets])),2)
              if(max(abs(M))!=0){M<-M/max(abs(M[-targets]))}
              long<-1:length(M)
              long<-long[-targets]
              for(j in long){
                
                if(M[j]>0){
                  color[j]<-rgb(couleur2(max(M[j],0.2)),alpha=255,maxColorValue=255)
                }
                else{
                  color[j]<-rgb(couleur1(max(-M[j],0.2)),alpha=255,maxColorValue=255)
                  
                }
                if(M[j]==0){
                  color[j]<- rgb(couleur2(0),alpha=255,maxColorValue=255)
                }
                
              }
              if(!is.null(targets)){color[targets]<-"green"}
              plot(net,nv=nv,color.vertex=color,ini=P,label_v=label_v,frame.color=frame.color,label.hub=label.hub,
                   edge.arrow.size=edge.arrow.size,edge.thickness=edge.thickness)
              couleur3<-colorRamp(c("blue","grey","red"))
              nb.col<-300
              coll<-rgb(couleur3(1:nb.col/nb.col),maxColorValue=255)
              
              rect(seq(-1,-0.5,length.out=nb.col)[1:(nb.col-1)],-1,seq(-1,-0.5,length.out=nb.col)[2:nb.col],-0.95,border="transparent",col=coll)
              text(-1,-1.05,-maxi)
              text(-0.51,-1.05,maxi)
              text(-0.751,-1.05,"0")	
              if(time==0){time<-length(omics@time)}
              text(-0.8,1.1,paste("Time point prediction =",time, sep=" "))
            }
            
          }
)
