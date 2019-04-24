setMethod("predict"
          ,c("micro_array")
          ,function(object
                    ,Omega
                    ,act_time_group=NULL
                    ,nv=0
                    ,targets=NULL
                    ,adapt=TRUE
          ){
#            require(magic)
            
            micro<-object
            if(!is.null(targets)){
              micro@microarray[targets,]<-0
            }
            if(!is.null(act_time_group)){
              stop("Cluster activation times must be provided in the act_time_group numeric vector.")
            }
            #groups
            groupe<-micro@group
            #measurements
            M<-micro@microarray
            #number of timepoints
            T<-length(micro@time)
            #gene groups
            gr<-micro@group 

            #group vector
            ngrp<-length(unique(gr))

            if(all(micro@gene_ID==0)){
              gene<-1:length(groupe)
            } else {
              gene<-micro@gene_ID
            }
            gene2<-gene
            #First timepoints for all the subjects
            supp<-seq(1,T*object@subject,T)
            #All the timepoints
            supp2<-1:(T*object@subject)
            #All timepoints except the first one
            supp2<-supp2[-supp]
            #Removing silenced genes
            if(!is.null(targets)){
              gene<-gene[-targets]
            }
            #Links
            O<-Omega@network
            #Cutoff
            O[abs(O)<nv]<-0
            #F matrix
            F<-Omega@F
            #Encore le microarray
            microP<-micro
            #predictors
            sup_pred<-rep(1:T,micro@subject)+rep(seq(0,T*(micro@subject-1),T),each=T)
            
            
            #still micro_array
            micro2<-micro
            #F matrix index
            u<-0

            #loop on the T timepoints to iteratively fill the expression data
            #we need to select the gene groups that were activated before time peak
            #act_time_group
#
#            for(peak in 2:(T)){
#  
            for(grpjj in 1:ngrp){
              #On garde les groupes predicteurs
#              IND<-which(groupe[gene2]%in%1:(peak-1))
              
              IND<-which(groupe[gene]%in%(1:ngrp)[-grpjj])
              grIND<-groupe[IND]
              #On silence les targets
              if(!is.null(targets)){micro2@microarray[targets,]<-micro@microarray[targets,]}
              #Genes predicteurs aux temps 1..T
              pred<-micro2@microarray[IND,sup_pred]
              #print(pred[IND==72,])
              #Pour chaque groupe sauf le grpjj
              for(k in (1:ngrp)[-grpjj]) {
                #On prend les indices du kieme groupe
                ind<-which(grIND %in% k)
                #Fonction produit de la matrice F avec un vecteur
                f<-function(x){(F[,,grpjj+(k-1)*ngrp]%*%(x))}#/(sqrt(sum((F[,,grpjj+(k-1)*ngrp]%*%(x))^2)))}
                #Pour chaque sujet
                for(i in 1:micro@subject){
                  pred[ind,1:T+(i-1)*T]<-t(apply(pred[ind,1:T+(i-1)*T, drop = FALSE],1,f))
                }
              }
              pred[is.na(pred)]<-0
              #Indices des genes du groupe d'interet
              IND2<-which(groupe[gene]==(grpjj))
              #Pour chaque gene du groupe reponse
              for(j in gene[IND2]){
                #predj<-(pred)*O[IND,j]
                #Pour chaque gene on cherche les genes du groupe predicteur qui ont un lien non nul
                predj<-pred[O[IND,j]!=0,]
                #Si il y en a
                if(length(predj)!=0){
                  #On isole les valeurs finales dans Y (reponse)
                  Y<-micro@microarray[j,sup_pred]
                  if(adapt==TRUE){
                    if(!is.null(dim(predj))){
                      mm<-lm(Y~t(predj)-1)
                    }
                    else{
                      mm<-lm(Y~(predj)-1)
                    }
                    #On remplace par la valeur inferee
                    micro2@microarray[j,sup_pred]<-predict(mm)
                    #On update les coefficients de Omega
                    O[IND,j][O[IND,j]!=0]<-coef(mm)[]
                  }
                  #Si il n'y en a pas
                }
                else{
                  #On prend la somme ponderee des sorties 
                  #des genes exprimes aux autres temps
                  predj<-apply((pred)*O[IND,j],2,sum)	
                  micro2@microarray[j,sup_pred]<-predj[sup_pred]
                }
                
              }
            }
            #Pour pouvoir faire un plot
            micro33<-micro2
            if(!is.null(targets)){ 
              pppp<- unique(unlist(geneNeighborhood(Omega,targets,nv,graph=FALSE)))
              genes3<-gene2[-pppp]                           
            }            else{
              genes3<-gene2
            }
            if(!is.null(targets)){
              micro33@microarray[ genes3,]<-micro@microarray[genes3,]
            }		
            
            if(is.null(targets)){
              targets<- -1
            }
            subjects<-object@subject
            times<-object@time
            ntimes<-length(times)
            patients<-paste(rep("P",subjects*ntimes),rep(1:subjects,each=ntimes),sep="")
            temps<-paste(rep("T",subjects*ntimes),rep(times,subjects),sep="")
            indicateurs<-paste(patients,temps,sep="")
            expr<-rep("log(S/US)",subjects*ntimes)
            nomscol<-paste(expr,":",indicateurs)
            #nomscol<-paste(patients,timestring1)
            
            colnames(object@microarray)<-nomscol
            colnames(micro@microarray)<-nomscol
            colnames(micro33@microarray)<-nomscol
            return(new("micropredict"
                       ,microarray_unchanged=object
                       ,microarray_changed=micro
                       ,microarray_predict=micro33
                       ,nv=nv
                       ,network=Omega#Avant =network
                       ,targets=targets
            )
            )	
          }
)



setMethod(f="inference"
          ,signature=c("micro_array")
          ,definition=function(M
                               ,tour.max=30
                               ,g=function(x){1/x}
                               ,conv=0.001
                               ,cv.subjects=TRUE
                               ,nb.folds=NULL
                               ,eps=10^-5
                               ,type.inf="iterative"
                               ,Fshape=NULL
                               ,Finit=NULL
                               ,Omega=NULL
                               ,fitfun="LASSO"
                               ,use.Gram=TRUE
                               ,error.stabsel=0.05
                               ,pi_thr.stabsel=0.6
                               ,priors=NULL
                               ,mc.cores=getOption("mc.cores", 2L)
                               ,intercept.stabpath=TRUE
                               ,steps.seq=.95
                               ,limselect=.95
                               ,use.parallel=TRUE
                               ){

            #Fshape=NULL
            #Finit=NULL
            #
            #M<-Selection
            #tour.max=30
            #g=function(x){1/x}
            #conv=0.001
            #cv.subjects=TRUE
            #nb.folds=NULL
            #eps=10^-5
            #type.inf="iterative"
            #Omega=NULL
            #TempEnforce=FALSE
            
            #Package requis
            require(nnls)
            #Quelques indicateurs
            mat<-M@microarray
            
            if(is.null(priors)) priors<-matrix(1,nrow(mat),nrow(mat))
            if(!is.matrix(priors)) stop("priors should be a matrix")
            if(!prod(dim(priors) == rep(nrow(mat),2))==1) stop("priors should have the same dimension than omega")
            #cat(str(priors))
            
            #La matrice contenant les donnees
            gr<-M@group 
            N<-dim(mat)[1] 
            # Le vecteur des groupes
            ngrp<-length(unique(gr))
            #Nombre de genes
            T<-length(unique(M@time)) 
            #Nombre de mesures par patient (souvent nombre de temps de mesures)
            sqF<-length(unique(M@time)) 
            #Nombre de patients
            P<-M@subject
            #On compte toutes les matrices m?mes celles sur la diagonales. F est form?e de T*T blocs
            nF<-ngrp*ngrp 	
            charslist=vector("list",nF)
            
            #La condition suivante determine le nombre de folds 
            # pour la cross validation
            if(is.null(nb.folds)){
              K<-sqF-1
            }            else{
              K<-nb.folds
            }
            
            if(is.null(Finit)){
              #Initialisation des matrices Fabinit par defaut
              # F[,,1] correspond à la matrice F11
              # F[,,2] correspond à la matrice F12
              # F[,,T+1] correspond à la matrice F21
              # F[,,T+2] correspond à la matrice F22
              # F[,,2T] correspond à la matrice F2T
              # ... (toutes les matrices y compris les nulles sont indicées F11 par exemple)
              #Lexicographic order for matrices 11 -> 1T -> 21 -> 2T -> 32 -> .... TT
              #row number is (ii-1)%/%T+1
              #column number is ii%%T if ii%%T!=0 and T if ii%%T==0
              Finit<-array(0,c(sqF,sqF,nF))	
              for(ii in 1:nF){    
                if((ii%%(ngrp+1))==1){
                  Finit[,,ii]<-0
                } else {
                  Finit[,,ii]<-cbind(rbind(rep(0,sqF-1),diag(1,sqF-1)),rep(0,sqF))+rbind(cbind(rep(0,sqF-1),diag(1,sqF-1)),rep(0,sqF))
                }
              }
            } else {
              if(dim(Finit)[3]!=nF){stop("Wrong number of Finit matrices: ",dim(Finit)[3]," instead of ",nF,".",sep="")}
              if((dim(Finit)[1]!=sqF)|(dim(Finit)[2]!=sqF)){stop("Finit matrices must be squared of order ",sqF,".",sep="")}
            }
            F <- Finit
            if(is.null(Fshape)){
              #Initialisation des matrices Fabshape par defaut
              #J'ai prefere construire une matrice F tri dimensionnelle laquelle contient toutes les matrices Fab
              # F[,,1] correspond a la matrice F11
              # F[,,2] correspond a la matrice F12
              # F[,,T+1] correspond a la matrice F21
              # F[,,T+2] correspond a la matrice F22
              # F[,,2T] correspond a la matrice F2T
              # ... (toutes les matrices y compris les nulles sont indicées F11 par exemple)
              #Lexicographic order for matrices 11 -> 1T -> 21 -> 2T -> 32 -> .... TT
              #row number is (ii-1)%/%T+1
              #column number is ii%%T if ii%%T!=0 and T if ii%%T==0
              Fshape<-array("0",c(sqF,sqF,nF)) 
              for(ii in 1:nF){  
                if((ii%%(ngrp+1))==1){
                  Fshape[,,ii]<-"0"
                } else {
                  lchars <- paste("a",1:(2*sqF-1),sep="")
                  tempFshape<-matrix("0",sqF,sqF)
                  for(bb in (-sqF+1):(sqF-1)){
                    tempFshape<-replaceUp(tempFshape,matrix(lchars[bb+sqF],sqF,sqF),-bb)
                  }
                  tempFshape <- replaceBand(tempFshape,matrix("0",sqF,sqF),0)
                  Fshape[,,ii]<-tempFshape
                }
              }
            } else {
              if(dim(Fshape)[3]!=nF){stop("Wrong number of Fshape matrices: ",dim(Fshape)[3]," instead of ",nF,".",sep="")}
              if((dim(Fshape)[1]!=sqF)|(dim(Fshape)[2]!=sqF)){stop("Fshape matrices must be squared of order ",sqF,".",sep="")}
            }
            
            #Initialisation de la matrice Omega
            
            if(is.null(Omega)){
              Omega<-array(0,c(N,N))
            } else 
            {
              if(!((dim(Omega)[1]==N)|(Omega[2]==N))){stop(paste("The Omega matrix must be squared of order N=",N,".",sep=""))}
              if(!(all(Omega==0|Omega==1))){stop("The Omega coordinates must all be 0s or 1s")}
              #if(any(Omega<0|Omega>1)){stop("The Omega coordinates must all lie between 0 and 1")}
            }
            
            
            #Initialisation des deux indicateurs de convergence
            convF<-rep(mean(F^2),nF)
            convO<-mean(mat^2)
            
            #Support : correspond aux colonnes de mat qui servent pour la prediction
            #Attention le modele se sert que des temps 1 a T-1 pour chaque patient pour les predicteurs
            
            sup_pred<-rep(1:sqF,P)+rep(seq(0,sqF*(P-1),sqF),each=sqF)
            
            #Initialisation du nombre de tours	
            tour<-1
            
            #Si on veut faire la version non iterative
            # il faut deux passages : dans le premier on infere
            # F et dans le second on infere omega
            if(type.inf=="noniterative"){
              tour.max<-2
            }
            
            #L'algorithme commence ici
            while((tour<= tour.max && convO[length(convO)]>conv) || tour<=2){ 
              #Condition d'arret : soit nombre de tour max atteint, 
              #soit la convergence de la matrice Omega est suffisante
              
              cat(paste("We are at step : ",tour))
              cat("\n") 
              #Pour montrer que l'algorithme est en train de calculer
              
              OmegaS<-Omega 
              #OmegaS comme sauvegarde ; necessaire pour calculer 
              # la convergence
              
              u<-0 
              #u a un role essentiel, puisque qu'il determine 
              # laquelle des matrices Fab est indicee. 
              # Se referer plus haut pour en connaitre l'ordre
              cat("Computing Group (out of ",ngrp,") : ",sep="")
              for(grpjj in 1:ngrp){ 
                cat("\n",grpjj);#if(grpjj < ngrp){cat(" ... ")}
                #Comprendre ici que grpjj correspond au groupe du gene REPONSE.
                #Ici nous cherchons les genes possiblement predicteurs.
                #Nous recuperons le groupe des individus possiblement predicteurs.			
                
                IND<-which(gr %in% (1:ngrp)[-(grpjj)])
                #Ici nous cherchons les genes possiblement 
                # PREDICTEURS.
                #En effet, le gene reponse etant de groupe grpjj, 
                # un predicteur ne peut etre du groupe grpjj
                
                grIND<-gr[IND] 			
                #Nous recuperons le groupe des individus possiblement
                # predicteurs.
                #Agir ici avec Omega binaire
                pred<-mat[IND,sup_pred]		
                #Nous creons la matrice des predicteurs 
                #Ou ici avec Omega ponderation
                
                for(k in 1:ngrp) {
                  #cat(paste("(",k,")",sep=""));
                  #Cette boucle sert a transformer les
                  #predicteurs en fonction des matrices Fab
                  #Le groupe de la variable reponse
                  # etant ngrp, les groupes des predicteurs 
                  # sont tous les groupes
                  
                  ind<-which(grIND %in% k) 
                  #On regarde successivement, et dans 
                  # l'ordre, tous les groupes possibles			
                  #						u<-u+1 
                  #u est initialise a 0. 
                  #La premiere fois, u vaut donc 1, grpjj 1, et k =1. Donc F[,,u]=F11.
                  #La deuxieme fois qu'on arrive ici, grpjj est 1, et k vaut 2 F[,,u]=F21
                  #La troisieme fois grpjj=1 k=3 donc F[,,u]=F31
                  # ... c'est bien l'ordre dans lequel nous avons range les F			
                  #       grpjj+(k-1)*ngrp  grpjj=1 k=1 u=1 F11
                  #       grpjj+(k-1)*ngrp  grpjj=1 k=2 u=5 F21
                  
                  f<-function(xf){(F[,, grpjj+(k-1)*ngrp]%*%xf)} 
                  #On construit une fonction generique
                  
                  for(i in 1:P){			
                    pred[ind,1:(sqF)+(i-1)*(sqF)]<-t(apply(pred[ind,1:(sqF)+(i-1)*(sqF)],1,f)) 
                    #La transformation est faite ici.
                  }
                }
                #En sortant de la boucle avec k, pred 
                #contient les genes predicteurs correctement 
                #transformes. 				
                
                pred[is.na(pred)]<-0 
                #Ceci est une securite, au cas ou une matrice F
                #deviendrait nulle
                
                #On construit ci dessous la matrice des vecteurs
                # reponse
                Y<-mat[which(gr %in% grpjj),sup_pred]
                Omega[IND, which(gr %in% grpjj)]<-Omega[IND, which(gr %in% grpjj)]*0
                
                #Nous allons passer au Lasso
                if(fitfun=="LASSO2"){
                  
                 priors2<-priors[IND,which(gr %in% grpjj)]
                 Y2<-cbind(1:nrow(Y),Y)
                  if(norm(pred,type="F")>eps){     
                    fun_lasso2<-function(x){cat(".");
                      lasso_reg2(pred,x[-1],nfolds=P,foldid=rep(1:P,each=ncol(pred)/P),priors=priors2[,x[1]])}
                    options(show.error.messages = FALSE)
                    Omega[IND, which(gr %in% grpjj)]<-apply(Y2,1,fun_lasso2)
                    options(show.error.messages = TRUE)
                  }
                }
                
                
                if(fitfun=="LASSO"){
                  if(norm(pred,type="F")>eps){     
                    if(cv.subjects==TRUE){
                      cv.folds1=function(n,folds){
                      split(1:dim(pred)[2]
                            ,rep(1:P,each=dim(pred)[2]/P))}
                    #                  cv.fun.name="cv.folds1"
                    } else {
                      cv.folds1=lars::cv.folds
                    #                    cv.fun.name="lars::cv.folds"
                    }
                  #                cat(cv.fun.name)
                  fun_lasso<-function(x){lasso_reg(pred,x,K=K,eps,cv.fun=cv.folds1
                                                 #,cv.fun.name=cv.fun.name
                  )} 
                  # retenir<-lars::cv.folds 
                  #Ceci permet de changer une fonction interne de lars
                  # qui s'occupe de la validation croisee. 
                  Omega[IND, which(gr %in% grpjj)]<-apply(Y,1,fun_lasso)
                  }
                }
                if(fitfun=="SPLS"){
                  if(norm(pred,type="F")>eps){     
                    if(cv.subjects==TRUE){
                      cv.folds1=function(n,folds){
                        split(1:dim(pred)[2],rep(1:P,each=dim(pred)[2]/P))
                        }
                    } else {
                      cv.folds1=function(n, folds){return(split(sample(1:n), rep(1:folds, length = n)))}
                      }
                      fun_spls<-function(x){spls_reg(pred,x,K=K,eps)} 
                      Omega[IND, which(gr %in% grpjj)]<-apply(Y,1,fun_spls)
                  }
                }
                if(fitfun=="ELASTICNET"){
                  if(norm(pred,type="F")>eps){     
                    if(cv.subjects==TRUE){
                      cv.folds1=function(n,folds){
                        split(1:dim(pred)[2],rep(1:P,each=dim(pred)[2]/P))
                        }
                    }else{
                      cv.folds1=lars::cv.folds
                      }
                      fun_enet<-function(x){lasso_reg(pred,x,K=K,eps,cv.fun=cv.folds1)} 
                      Omega[IND, which(gr %in% grpjj)]<-apply(Y,1,fun_enet)
                  }
                }
                if(fitfun=="stability.c060"){
                  require(c060);#cat(".")
                  cat("mc.cores=",mc.cores,sep="")
                  
                  fun_stab<-function(g,mc.cores=mc.cores,intercept.stabpath=intercept.stabpath){
                    
                    if(sum(pred)==0){
                      return(rep(0,nrow(pred)))
                      cat(".")               
                    }else{
                  LL=rep(0,nrow(pred));error.inf=TRUE
                  #cat(intercept.stabpath)
                  try({essai<-c060::stabpath(g,t(pred),mc.cores=mc.cores,intercept=intercept.stabpath);
                  varii<-c060::stabsel(essai,error=error.stabsel,pi_thr=pi_thr.stabsel)$stable;
                  lambda<-c060::stabsel(essai,error=error.stabsel,pi_thr=pi_thr.stabsel)$lambda;
#                  L<-lars(t(pred),g,use.Gram=use.Gram)
#                  LL<-predict(L,s=lambda,mode="lambda",type="coef")$coefficients
                  L<-glmnet::glmnet(t(pred),g,intercept=intercept.stabpath);
                  LL<-as.matrix(predict(L,s=lambda,type="coef"))[-1,1]})
                  try({LL[-varii]<-0;
                  error.inf=FALSE})
                  if(error.inf){cat("!")} else {cat(".")}
                  if(!is.vector(LL)){LL<-rep(0,nrow(pred))}
                  return(LL)
                  }
                  }
                
                  options(show.error.messages = FALSE)
                  Omega[IND, which(gr %in% grpjj)]<-apply(Y,1,fun_stab,mc.cores=mc.cores,intercept.stabpath=intercept.stabpath)
                  options(show.error.messages = TRUE)
                }
                
                
                
                
                
                if(fitfun=="stability.c060.weighted"){
                  require(c060);#cat(".")
                  cat("mc.cores=",mc.cores,sep="")
                  priors2<-priors[IND,which(gr %in% grpjj)]
                  Y2<-cbind(1:nrow(Y),Y)
                  
                  fun_stab_weighted<-function(g,mc.cores=mc.cores,intercept.stabpath=intercept.stabpath,penalty.factor=penalty.factor){

                      if(sum(pred)==0){
                      return(rep(0,nrow(pred)))
                        cat(".")               
                    }else{
                      stabpath <- function (y, x, size = 0.632, steps = 100, penalty.factor = rep(1,ncol(x)), weakness=1, mc.cores = getOption("mc.cores", 
                                                                                                                                               2L), ...) 
                      {
                        fit <- glmnet::glmnet(x, y, ...)
                        if (class(fit)[1] == "multnet" | class(fit)[1] == "lognet") 
                        y <- as.factor(y)
                        p <- ncol(x)
                        subsets <- sapply(1:steps, function(v) {
                          sample(1:nrow(x), nrow(x) * size)
                        })
                        if (.Platform$OS.type != "windows") {
                          res <- parallel::mclapply(1:steps, mc.cores = mc.cores, glmnet.subset.weighted, 
                                                    subsets, x, y, lambda = fit$lambda, penalty.factor, weakness, p, 
                                                    ...)
                        }
                        else {
                          cl <- parallel::makePSOCKcluster(mc.cores)
                          parallel::clusterExport(cl, c("glmnet", "drop0"))
                          res <- parallel::parLapply(cl, 1:steps, glmnet.subset.weighted, subsets, 
                                           x, y, lambda = fit$lambda, penalty.factor, weakness, p, ...)
                          parallel::stopCluster(cl)
                        }
                        #str(res)
                        res <- res[unlist(lapply(lapply(res, dim), function(x) x[2] == 
                                                   dim(res[[1]])[2]))]
                        x <- as.matrix(res[[1]])
                        qmat <- matrix(ncol = ncol(res[[1]]), nrow = length(res))
                        qmat[1, ] <- colSums(as.matrix(res[[1]]))
                        for (i in 2:length(res)) {
                          qmat[i, ] <- colSums(as.matrix(res[[i]]))
                          x <- x + as.matrix(res[[i]])
                        }
                        x <- x/length(res)
                        qs <- colMeans(qmat)
                        out <- list(fit = fit, x = x, qs = qs)
                        class(out) <- "stabpath"
                        return(out)
                      }
                      
                      glmnet.subset.weighted <- function (index, subsets, x, y, lambda, penalty.factor, weakness, p, ...) 
                      {
                        if (length(dim(y)) == 2 | class(y) == "Surv") {
                          glmnet::glmnet(x[subsets[, index], ], y[subsets[, index], ], 
                                         lambda = lambda, penalty.factor = penalty.factor*1/runif(p, weakness, 1), ...)$beta != 0
                        }
                        else {
                          if (is.factor(y) & length(levels(y)) > 2) {
                            temp <- glmnet::glmnet(x[subsets[, index], ], y[subsets[, 
                                                                                    index]], lambda = lambda, penalty.factor = penalty.factor*1/runif(p, weakness, 1), ...)[[2]]
                            temp <- lapply(temp, as.matrix)
                            Reduce("+", lapply(temp, function(x) x != 0))
                          }
                          else {
                            glmnet::glmnet(x[subsets[, index], ], y[subsets[, index]], 
                                           lambda = lambda, penalty.factor = penalty.factor*1/runif(p, weakness, 1), ...)$beta != 0
                          }
                        }
                      }
                      
#                      res <- stabpath(y,x,penalty.factor=c(rep(1,999),rep(0,1)),mc.cores=2)
                      #assignInNamespace("stabpath",,ns="c060")    
                      #assignInNamespace("glmnet.subset",, ns="c060") 
                      #assignInNamespace("glmnet.subset.weighted",, ns="c060")    
                      
                      rm(varii)
                      LL=rep(0,nrow(pred));error.comp=TRUE;error.inf=TRUE
                      #cat(intercept.stabpath)
                      try({
                      essai<-stabpath(g[-1],t(pred),mc.cores=mc.cores,intercept=intercept.stabpath,penalty.factor=priors2[,g[1]]);
                      varii<-c060::stabsel(essai,error=error.stabsel,pi_thr=pi_thr.stabsel)$stable;
                      lambda<-c060::stabsel(essai,error=error.stabsel,pi_thr=pi_thr.stabsel)$lambda;
                      #                  L<-lars(t(pred),g,use.Gram=use.Gram)
                      #                  LL<-predict(L,s=lambda,mode="lambda",type="coef")$coefficients
                      L<-glmnet::glmnet(t(pred),g[-1],intercept=intercept.stabpath,penalty.factor=priors2[,g[1]]);
                      LL<-as.matrix(predict(L,s=lambda,type="coef"))[-1,1];
                      error.comp=FALSE
                      })
                      if(error.comp){cat("!",geterrmessage(),"\n")} else {cat("")}
                      try({LL[-varii]<-0;error.inf=FALSE})
                      if(!error.comp&error.inf){cat("!",geterrmessage(),"\n")} else {cat(".")}
                      if(!is.vector(LL)){LL<-rep(0,nrow(pred))}
                      return(LL)
                    }
                  }
                  
                  options(show.error.messages = FALSE)
                  Omega[IND, which(gr %in% grpjj)]<-apply(Y2,1,fun_stab_weighted,mc.cores=mc.cores,intercept.stabpath=intercept.stabpath)
                  options(show.error.messages = TRUE)
                }
                
                
                
                if(fitfun=="robust"){
                  #require(movMF)
                  require(lars)
                  #require(msgps)
                  fun_robust<-function(g){
                    if(sum(pred)==0){
                      return(rep(0,nrow(pred)))
                    }else{
                      essai<-boost(t(pred)+rnorm(prod(dim(pred)),0,0.001),g)  
                      varii<-which(essai==1)
                      lambda<-0
                      L<-lars::lars(t(pred),g)
                      
                      LL<-predict(L,s=lambda,mode="lambda",type="coef")$coefficients
                      LL[-varii]<-0
                      return(LL)
                    }
                  }
                  
                  Omega[IND, which(gr %in% grpjj)]<-apply(Y,1,fun_robust)
                  
                }  
                
                if(fitfun=="selectboost.weighted"){
                  #require(movMF)
                  require(lars)
                  #require(msgps)
                  requireNamespace("SelectBoost")
                  priors2<-priors[IND,which(gr %in% grpjj)]
                  Y2<-cbind(1:nrow(Y),Y)

                  fun_robust_weighted<-function(g,mc.cores=mc.cores, steps.seq = steps.seq, limselect = limselect, use.parallel = use.parallel){
                    if(norm(pred,type="F")<=eps){     
                      return(rep(0,nrow(pred)))
                      cat(".")               
                    }else{
                      if(exists("varii")){rm(varii)}
                      LL=rep(0,nrow(pred));error.comp=TRUE;error.inf=TRUE
                      #cat(intercept.stabpath)
                      try({
                      essai<-suppressMessages(SelectBoost::fastboost(t(pred),g[-1],SelectBoost::group_func_2,SelectBoost::lasso_cv_glmnet_min_weighted,corrfunc="crossprod",normalize=TRUE, B=100, use.parallel=use.parallel, ncores=mc.cores,c0lim=FALSE, steps.seq = steps.seq, priors=priors2[,g[1]]))
                      varii<-which(essai>=limselect)

                      resultat<-glmnet::cv.glmnet(t(pred),g[-1],nfolds=10,penalty.factor=priors2[,g[1]])
                      LL<-predict(resultat,s="lambda.min",type="coef")[-1,1]
                      error.comp=FALSE
                      })
                      if(error.comp){cat("!",geterrmessage(),"\n")} else {cat("")}
                      try({LL[-varii]<-0;error.inf=FALSE})
                      if(!error.comp&error.inf){cat("!",geterrmessage(),"\n")} else {cat(".")}
                      if(!is.vector(LL)){LL<-rep(0,nrow(pred))}
                      return(LL)
                    }
                  }
                  Omega[IND, which(gr %in% grpjj)]<-apply(Y2,1,fun_robust_weighted,mc.cores=mc.cores,steps.seq=.95,limselect=.95,use.parallel=use.parallel)
                  }  
                
                
                                
              }
              cat("\n")
              #fin de la boucle for avec peak ; 
              #la matrice omega est inferee
              
              co<-apply(Omega,2,sumabso)
              Omega<-t(t(Omega)/co)
              
              if(tour!=1 && type.inf=="iterative"){
                Omega<-(g(tour)*Omega+OmegaS)/(1+g(tour)) 
                #On prend seulement une partie de l'innovation
              }
              
              
              convO<-c(convO,mean(abs(Omega-OmegaS)))
              
              if( type.inf=="iterative"){
                cat(paste("The convergence of the network is (L1 norm) :", round(convO[length(convO)],5)))	
                cat("\n")	
              }
              uuu<-0
              sauvF<-F
              
              if(tour==1 && type.inf=="noniterative"){
                Omega<-Omega*0+1 
                #Tous les pr?dicteurs sont mis egaux a 1
                # dans le cadre de l'inf?rence non iterative
              }
              
              #Maintenant, il s'agit d'estimer la matrice F ;
              # on reprend la meme maniere de faire que pour Omega
              #Nous annotons les differences
              
              for(grpjj in 1:ngrp){
                
                
                # En effet, le gene reponse etant de groupe grpjj, 
                #un predicteur ne peut etre que de groupe different
                #Important dans patterns : etant donne qu'il apparaisse dans 
                # les memes equations Fij, i=1...T, i<>j doivent etre estimees
                # en meme temps
                
                IND<-which(gr %in% (1:ngrp)[-(grpjj)])
                grIND<-gr[IND]
                sup_pred<-rep(1:sqF,P)+rep(seq(0,sqF*(P-1),sqF),each=sqF)
                pred<-(mat[IND,sup_pred])
                IND2<-which(gr %in% grpjj)
                
                charslist=vector("list",nF)
                Xf<-NULL
                for(i in (1:ngrp)[-grpjj]){ 
                  #Cette  boucle permet de creer la matrice des 
                  # predicteurs selon une forme pratique
                  X<-NULL
                  suma<-function(x){sum(abs(x))}
                  f<-Vectorize(function(x){
                    #Multiplie la matrice des pr??dicteurs du groupe i
                    # par les omega correspondants aux genes du groupe
                    # dont les indices sont dans x et on somme les 
                    # valeurs absolues
                    apply(pred[which(grIND==i),]*Omega[IND[which(grIND==i)],x],2,suma)})
                  Xa<-(f(IND2))
                  Xb<-NULL
                  FF=Fshape[,,(i-1)*ngrp+grpjj]
                  #        cat((i-1)*sqF+grpjj)
                  #        cat("\n")
                  #        cat(FF)
                  #        cat("\n")
                  chars=sort(setdiff(unique(as.vector(FF)),"0"))
                  if(length(chars)>0){
                    charslist[[(i-1)*ngrp+grpjj]]<-chars
                    for(p in 1:P){
                      Q<-NULL
                      q<-as.vector(Xa[1:sqF+(p-1)*sqF,])
                      for(cc in 1:length(chars)){
                        inds=FF==chars[cc]
                        ff<-function(xff){return(inds%*%xff)}
                        Q<-cbind(Q,unlist(tapply(q,factor(rep(1:length(IND2),rep(sqF,length(IND2)))),ff)))
                      }
                      X<-rbind(X,Q)
                    }
                    #        cat(Q)
                    #        cat("\n")
                    #        cat(X)
                    #        cat("\n")
                    Xf<-cbind(Xf,X)
                  }
                }
                
                if(!is.null(Xf)){
                  Y<-c(t(mat[IND2,sup_pred]))
                  pond<-rep(0,P)
                  coeffi<-array(0,c(P,dim(Xf)[2]))
                  #  rm(model)
                  for(pat in 1:P){
                    support<-1:length(Y)
                    enl<-(1:(length(Y)/(P))+(pat-1)*(length(Y)/(P)))
                    support<-support[-enl]                                                             
                    model<-nnls(abs(Xf[support,]),abs(Y[support]))
                    pond[pat]<-1/(mean((Xf[enl,]%*%coef(model)-Y[enl])^2) )
                    coeffi[pat,]<-coef(model)
                  }
                  
                  model<-apply(coeffi*( pond/sum(pond)),2,sum)
                }
                
                ncoeff=0
                for(jj in 1:ngrp){
                  if(!is.null(charslist[[grpjj+(jj-1)*ngrp]])){
                    charsjj=charslist[[grpjj+(jj-1)*ngrp]]
                    TempFF <- matrix(0,sqF,sqF)
                    FF=Fshape[,,grpjj+(jj-1)*ngrp]
                    for(aa in 1:length(charsjj)){
                      ncoeff=ncoeff+1
                      inds=FF==charsjj[aa]
                      TempFF=TempFF+inds*model[ncoeff]         
                    }
                    F[,,grpjj+(jj-1)*ngrp]<-TempFF
                    #        cat(TempFF)
                    #        cat("\n")
                    #        rm(TempFF)
                  }
                }
              } 
              #fin de la boucle peak
              
              if(tour==1 && type.inf=="noniterative"){
                Omega<-Omega*0 
                #Tous les pr?dicteurs sont remis egaux a 0
                # dans le cadre de l'inf?rence non iterative
              }
              
              
              if( type.inf=="iterative"){
                F<-(g(tour)*F+sauvF)/(1+g(tour))
              }
              cc<-rep(0,nF)
              for(i in 1:nF){cc[i]<-mean(abs((F[,,i]/sum(F[,,i])-sauvF[,,i]/sum(sauvF[,,i]))))}
              convF<-cbind(convF,cc)
              tour<-tour+1
            }
            
            if(type.inf=="iterative"){
              plot(convO[-1],type="l")
              if(!attr(dev.cur(),"names")=="pdf"){dev.new()}
              matplot(t(convF),type="l")
            }
            else{
              F<-sauvF   
            }
            result<-new("network"
                        ,network=Omega
                        ,name=M@name
                        ,F=F
                        ,convF=convF
                        ,convO=convO
                        ,time_pt=M@time
            )
            return(result)
          }
)


setMethod("gene_expr_simulation"
          ,"network"
          ,function(network
                    ,time_label=1:4
                    ,subject=5
                    ,peak_level=100
                    ,act_time_group=1:4
          ){
            require(VGAM)
            
            N<-network@network
            M<-matrix(0,dim(network@network)[1],length(unique(time_label))*subject)
            T<-length(unique(time_label))
            gene1<-which(time_label==1)
            supp<-seq(1,dim(M)[2],by=length(unique(time_label)))
            M[gene1,supp]<-VGAM::rlaplace(length(supp)*length(gene1),peak_level,peak_level*0.9)*(-1)^rbinom(length(supp)*length(gene1),1,0.5)
            supp<-(1:dim(M)[2])[-supp]
            M[gene1,supp]<-VGAM::rlaplace(length(supp)*length(gene1),0,peak_level*0.3) 
            
            
            for(i in 2:T){
              
              genei<-which(time_label==i)
              supp<-seq(1,dim(M)[2],by=length(unique(time_label)))
              M[genei,supp]<-VGAM::rlaplace(length(supp)*length(genei),0,peak_level*0.3)
              for(j in genei){
                for( t in 2:T){
                  M[j,supp+t-1]<-apply(N[,j]*M[,supp+(t-2)],2,sum) + 	rnorm(length(supp+t),0,50)	
                  
                }			
                
              }
              
              
            }
            
            MM<-as.micro_array(M,1:length(unique(time_label)),subject)
            MM@group<-time_label
            
            
            
            G<-Patterns::predict(MM,Omega=network,act_time_group=act_time_group)@microarray_predict #Before Cascade 1.03, we have G<-predict(MM,network,adapt=FALSE)@microarray_predict
            supp<-seq(1,dim(M)[2],by=length(unique(time_label)))
            G@microarray[,supp]<-M[, supp]	
            return(G)
          }
)
