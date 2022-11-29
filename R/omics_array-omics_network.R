#' Methods for Function \code{predict}
#' 
#' Prediction of the gene expressions after a knock-out experience for cascade
#' networks.
#' 
#' The plot of prediction of knock down experiments (i.e. targets<>NULL) is
#' still in beta testing for the moment.
#' 
#' 
#' @aliases predict predict-methods predict,ANY-method
#' predict,omics_array-method
#' @param object a omics_array object.
#' @param Omega a omics_network object.
#' @param act_time_group [NULL] vector; at which time the groups (defined by sort(unique(group))) are activated ?
#' @param nv [=0] numeric ; the level of the cutoff
#' @param targets [NULL] vector ; which genes are knocked out ?
#' @param adapt [TRUE] boolean; do not raise an error if used with vectors 
#' 
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords methods
#' @examples
#' 
#' \donttest{
#' data(Selection)
#' data(infos)
#' pbst_NR4A1 = infos[infos$hgnc_symbol=="NR4A1", "affy_hg_u133_plus_2"]
#' pbst_EGR1 = infos[infos$hgnc_symbol=="EGR1", "affy_hg_u133_plus_2"]
#' gene_IDs = infos[match(Selection@name, infos$affy_hg_u133_plus_), "hgnc_symbol"]
#' 
#' data(networkCascade)
#' #A nv value can chosen using the cutoff function
#' nv = .02
#' NR4A1<-which(is.element(Selection@name,pbst_NR4A1))
#' EGR1<-which(is.element(Selection@name,pbst_EGR1))
#' P<-position(networkCascade,nv=nv)
#' 
#' #We predict gene expression modulations within the network if NR4A1 is experimentaly knocked-out. 
#' prediction_ko5_NR4A1<-predict(Selection,networkCascade,nv=nv,targets=NR4A1,act_time_group=1:4)
#' 
#' #Then we plot the results. Here for example we see changes at time points t2, t3 ans t4:
#' plot(prediction_ko5_NR4A1,time=2:4,ini=P,label_v=gene_IDs)
#' 
#' #We predict gene expression modulations within the network if EGR1 is experimentaly knocked-out. 
#' prediction_ko5_EGR1<-predict(Selection,networkCascade,nv=nv,targets=EGR1,act_time_group=1:4)
#' 
#' #Then we plot the results. Here for example we see changes at time point t2, t3 ans t4:
#' plot(prediction_ko5_EGR1,time=2:4,ini=P,label_v=gene_IDs)
#' }
#' 
setMethod("predict"
          , c("omics_array")
          , function(object
                     ,
                     Omega
                     ,
                     act_time_group = NULL
                     ,
                     nv = 0
                     ,
                     targets = NULL
                     ,
                     adapt = TRUE) {
            #            require(magic)
            
            omics <- object
            if (!is.null(targets)) {
              omics@omicsarray[targets, ] <- 0
            }
            if (is.null(act_time_group)) {
              stop("Cluster activation times must be provided in the act_time_group numeric vector.")
            }
            #groups
            groupe <- omics@group
            #measurements
            #M <- omics@omicsarray
            #number of timepoints
            T <- length(omics@time)
            #gene groups
            gr <- omics@group
            
            #group vector
            ngrp <- length(unique(gr))
            vgrp <- sort(unique(gr))
            
            if (all(omics@gene_ID == 0)) {
              gene <- 1:length(groupe)
            } else {
              gene <- omics@gene_ID
            }
            gene2 <- gene
            #Naming groupe vector for easy retrieving of group membership
            names(groupe) <- gene2
            #First timepoints for all the subjects
            supp <- seq(1, T * object@subject, T)
            #All the timepoints
            supp2 <- 1:(T * object@subject)
            #All timepoints except the first one
            supp2 <- supp2[-supp]
            #Removing silenced genes
            if (!is.null(targets)) {
              gene <- gene[-targets]
            }
            #Links
            O <- Omega@omics_network
            #Cutoff
            O[abs(O) < nv] <- 0
            colnames(O) <- gene2
            rownames(O) <- gene2
            #F matrix
            F <- Omega@F
            omicsP <- omics
            #predictors
            sup_pred <-
              rep(1:T, omics@subject) + rep(seq(0, T * (omics@subject - 1), T), each =
                                              T)
            
            
            omics2 <- omics
            #F matrix index
            u <- 0
            
            #loop on the T timepoints to iteratively fill the expression data
            #we need to select the gene groups that were activated before time peak
            #act_time_group
            #
            
            for (peak in 2:(T)) {
              for (grpjj in vgrp[act_time_group == peak]) {
                IND <- which(groupe[gene2] %in% vgrp[act_time_group < peak])
                grIND <- groupe[IND]
                if (!is.null(targets)) {
                  omics2@omicsarray[targets, ] <- omics@omicsarray[targets, ]
                }
                pred <- omics2@omicsarray[IND, sup_pred]
                for (k in (1:ngrp)[-grpjj]) {
                  ind <- which(grIND %in% k)
                  f <-
                    function(x) {
                      (F[, , grpjj + (k - 1) * ngrp] %*% (x))
                    }
                  for (i in 1:omics@subject) {
                    pred[ind, 1:T + (i - 1) * T] <-
                      t(apply(pred[ind, 1:T + (i - 1) * T, drop = FALSE], 1, f))
                  }
                }
                pred[is.na(pred)] <- 0
                IND2 <- which(groupe[gene] == (grpjj))
                for (j in gene[IND2]) {
                  predj <- pred[O[IND, j] != 0, ]
                  if (length(predj) != 0) {
                    Y <- omics@omicsarray[j, sup_pred]
                    if (adapt == TRUE) {
                      if (!is.null(dim(predj))) {
                        mm <- lm(Y ~ t(predj) - 1)
                      }
                      else{
                        mm <- lm(Y ~ (predj) - 1)
                      }
                      omics2@omicsarray[j, sup_pred] <- predict(mm)
                      O[IND, j][O[IND, j] != 0] <- coef(mm)[]
                    }
                  }
                  else{
                    predj <- apply((pred) * O[IND, j], 2, sum)
                    omics2@omicsarray[j, sup_pred] <- predj[sup_pred]
                  }
                  
                }
              }
            }
            omics33 <- omics2
            if (!is.null(targets)) {
              pppp <-
                unique(unlist(geneNeighborhood(Omega, targets, nv, graph = FALSE)))
              genes3 <- gene2[-pppp]
            }            else{
              genes3 <- gene2
            }
            if (!is.null(targets)) {
              omics33@omicsarray[genes3, ] <- omics@omicsarray[genes3, ]
            }
            
            if (is.null(targets)) {
              targets <- -1
            }
            subjects <- object@subject
            times <- object@time
            ntimes <- length(times)
            patients <-
              paste(rep("P", subjects * ntimes),
                    rep(1:subjects, each = ntimes),
                    sep = "")
            temps <-
              paste(rep("T", subjects * ntimes), rep(times, subjects), sep = "")
            indicateurs <- paste(patients, temps, sep = "")
            expr <- rep("log(S/US)", subjects * ntimes)
            nomscol <- paste(expr, ":", indicateurs)
            
            colnames(object@omicsarray) <- nomscol
            colnames(omics@omicsarray) <- nomscol
            colnames(omics33@omicsarray) <- nomscol
            return(
              new(
                "omics_predict"
                ,
                omicsarray_unchanged = object
                ,
                omicsarray_changed = omics
                ,
                omicsarray_predict = omics33
                ,
                nv = nv
                ,
                omics_network = Omega
                ,
                targets = targets
              )
            )
          })



#' Reverse-engineer the network
#' 
#' Reverse-engineer the network.
#' 
#' The fitting built-in fitting functions (`fitfun`) provided with the
#' `Patterns` package are : \describe{ \item{LASSO}{from the `lars` package
#' (default value)} \item{LASSO2}{from the `glmnet` package} \item{SPLS}{from
#' the `spls` package} \item{ELASTICNET}{from the `elasticnet` package}
#' \item{stability.c060}{from the `c060` package implementation of stability
#' selection} \item{stability.c060.weighted}{a new weighted version of the
#' `c060` package implementation of stability selection} \item{robust}{lasso
#' from the `lars` package with light random Gaussian noise added to the
#' explanatory variables} \item{selectboost.weighted}{a new weighted version of
#' the `selectboost` package implementation of the selectboost algorithm to
#' look for the more stable links against resampling that takes into account
#' the correlated structure of the predictors. If no weights are provided,
#' equal weigths are for all the variables (=non weighted case).} }
#' 
#' The weights are viewed as a penalty factors in the penalized regression
#' model: it is a number that multiplies the lambda value in the minimization
#' problem to allow differential shrinkage, [Friedman et al.
#' 2010](https://web.stanford.edu/~hastie/Papers/glmnet.pdf), equation 1 page
#' 3. If equal to 0, it implies no shrinkage, and that variable is always
#' included in the model. Default is 1 for all variables. Infinity means that
#' the variable is excluded from the model. Note that the weights are rescaled
#' to sum to the number of variables.
#' 
#' @name inference
#' @aliases inference inference-methods inference,omics_array-method
#' @param M a omics_array object.
#' @param tour.max [30] tour.max + 1 = maximal number of steps.
#' @param g After each step, the new solution is choosen as (the
#' old solution + g(x) * the new solution)/(1+g(x)) where x is the number of
#' steps. Defaults to `g=function(x) 1/x`
#' @param conv [0.001] Convergence criterion.
#' @param cv.subjects [TRUE] Subjectwise cross validation: should the cross validation be done by removing the subject one by one?
#' @param nb.folds [NULL] Relevant only if no subjectwise cross validation (i.e. cv.subjects=FALSE). The number of folds in cross validation.
#' @param eps [10^-5] Threshold for rounding coefficients to 0 (i.e. machine zero).
#' @param type.inf ["iterative"] "iterative" or "noniterative" : should the algorithm be computed iteratively or only for one step? For highly homogeneous clusters, the "noniterative" option is suffisant.
#' @param Fshape [NULL] Shape of the F matrix.
#' @param Finit [NULL] Init values of the F matrix.
#' @param Omega [NULL] Init values for the Omega matrix.
#' @param fitfun ["LASSO"] Function to infer the Omega matrix at each step.
#' @param use.Gram [TRUE] Optional parameter for the lasso in the `lars` package.
#' @param error.stabsel [0.05] Optional parameter for the stability selection algorithm in the `c060` package.
#' @param pi_thr.stabsel [0.6] Optional parameter for the stability selection algorithm in the `c060` package.
#' @param priors [NULL] A priori weights for
#' the links between the actors. 0 means that an actor is always included in
#' the predictive model, 1 is a neutral weighting and +infinity that the actor
#' is never used in the model. For a given predictive model, the weighting
#' vector is normalized so that its sum is equal to the number of predictors in
#' the model.
#' @param mc.cores [getOption("mc.cores", 2L)] Number of cores.
#' @param intercept.stabpath [TRUE] Use intercept in stability selection models?
#' @param steps.seq [.95] Optional parameter for the SelectBoost algorithm in the `SelectBoost` package. 
#' @param limselect [.95] Optional parameter for the SelectBoost algorithm in the `SelectBoost` package.
#' @param use.parallel [TRUE] Use parallel computing?
#' @param verbose [TRUE] Info on the completion of the fitting process
#' @param show.error.messages [FALSE] Should the error messages of the Omega estimating function be returned?
#' @return A omics_network object.
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @keywords methods
#' @examples
#' 
#' \donttest{
#' #With simulated data, default shaped F matrix and default LASSO from the lars package
#' #as fitting function
#' data(M)
#' infM <- inference(M)
#' str(infM)
#' plot(infM, choice="F", nround=0)
#' plot(infM, choice="F", nround=1)
#' 
#' #With simulated data, cascade network shaped F matrix (1 group per time measurement case) 
#' #and default LASSO from the lars package as fitting function
#' infMcasc <- inference(M, Finit=CascadeFinit(4,4), Fshape=CascadeFshape(4,4))
#' str(infMcasc)
#' plot(infMcasc, choice="F", nround=0)
#' plot(infMcasc, choice="F", nround=1)
#' 
#' #With selection of genes from GSE39411
#' data(Selection)
#' infSel <- inference(Selection, Finit=CascadeFinit(4,4), Fshape=CascadeFshape(4,4))
#' str(infSel)
#' str(infSel)
#' plot(infSel, choice="F", nround=0)
#' plot(infSel, choice="F", nround=1)
#' }
#' 
setMethod(f="inference"
          ,signature=c("omics_array")
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
                               ,verbose=TRUE
                               ,show.error.messages = FALSE
                               ){

            require(nnls, quietly = TRUE, warn.conflicts = FALSE);on.exit(unloadNamespace("package:nnls"))
            mat<-M@omicsarray
            
            if(is.null(priors)) priors<-matrix(1,nrow(mat),nrow(mat))
            if(!is.matrix(priors)) stop("priors should be a matrix")
            if(!prod(dim(priors) == rep(nrow(mat),2))==1) stop("priors should have the same dimension than omega")

            gr<-M@group 
            N<-dim(mat)[1] 
            ngrp<-length(unique(gr))
            T<-length(unique(M@time)) 
            sqF<-length(unique(M@time)) 
            P<-M@subject
            nF<-ngrp*ngrp 	
            charslist=vector("list",nF)
            
            if(is.null(nb.folds)){
              K<-sqF-1
            }            else{
              K<-nb.folds
            }
            
            if(is.null(Finit)){
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
            
            if(is.null(Omega)){
              Omega<-array(0,c(N,N))
            } else 
            {
              if(!((dim(Omega)[1]==N)|(Omega[2]==N))){stop(paste("The Omega matrix must be squared of order N=",N,".",sep=""))}
              if(!(all(Omega==0|Omega==1))){stop("The Omega coordinates must all be 0s or 1s")}
              #if(any(Omega<0|Omega>1)){stop("The Omega coordinates must all lie between 0 and 1")}
            }
            
            convF<-rep(mean(F^2),nF)
            convO<-mean(mat^2)
            
            sup_pred<-rep(1:sqF,P)+rep(seq(0,sqF*(P-1),sqF),each=sqF)
            
            tour<-1
            
            if(type.inf=="noniterative"){
              tour.max<-2
            }
            
            while((tour<= tour.max && convO[length(convO)]>conv) || tour<=2){ 

              if(verbose){
              cat(paste("We are at step : ",tour))
              cat("\n") 
              }
              OmegaS<-Omega 

              u<-0 
              if(verbose){
                cat("Computing Group (out of ",ngrp,") : ",sep="")
              }
              for(grpjj in 1:ngrp){ 
                if(verbose){
                  cat("\n",grpjj)
                }

                IND<-which(gr %in% (1:ngrp)[-(grpjj)])

                grIND<-gr[IND] 			
                pred<-mat[IND,sup_pred]		

                for(k in 1:ngrp) {

                  ind<-which(grIND %in% k) 

                  f<-function(xf){(F[,, grpjj+(k-1)*ngrp]%*%xf)} 
                  #generic function
                  
                  for(i in 1:P){			
                    pred[ind,1:(sqF)+(i-1)*(sqF)]<-t(apply(pred[ind,1:(sqF)+(i-1)*(sqF)],1,f)) 
                    #transform is done here.
                  }
                }
                #predictors gens are now in pred and correctly transfromed
                
                pred[is.na(pred)]<-0 
                #Protection in case of a F matrix becoming null.
                
                #Responses' matrix
                Y<-mat[which(gr %in% grpjj),sup_pred]
                Omega[IND, which(gr %in% grpjj)]<-Omega[IND, which(gr %in% grpjj)]*0
                
                if(fitfun=="LASSO2"){
                  
                 priors2<-priors[IND,which(gr %in% grpjj)]
                 Y2<-cbind(1:nrow(Y),Y)
                  if(norm(pred,type="F")>eps){
                    save_show.error.messages = options()$show.error.messages
                    options(show.error.messages = show.error.messages)
                    on.exit(options(show.error.messages = save_show.error.messages))
                    if(cv.subjects==TRUE){
                      fun_lasso2<-function(x){if(verbose){cat(".")};
                        lasso_reg2(pred,x[-1],foldid=rep(1:P,each=ncol(pred)/P),priors=priors2[,x[1]])}
                    } else {
                      fun_lasso2<-function(x){if(verbose){cat(".")};
                        lasso_reg2(pred,x[-1],foldid=sample(rep(1:K,length=ncol(pred))),priors=priors2[,x[1]])}
                    }
                    Omega[IND, which(gr %in% grpjj)]<-apply(Y2,1,fun_lasso2)
                    options(show.error.messages = save_show.error.messages)
                  }
                }
                
                
                if(fitfun=="LASSO"){
                  if(norm(pred,type="F")>eps){     
                    save_show.error.messages = options()$show.error.messages
                    options(show.error.messages = show.error.messages)
                    on.exit(options(show.error.messages = save_show.error.messages))
                    if(cv.subjects==TRUE){
                      cv.folds1=function(n,folds){
                      split(1:dim(pred)[2]
                            ,rep(1:P,each=dim(pred)[2]/P))}
                    } else {
                      cv.folds1=lars::cv.folds
                    }
                  fun_lasso<-function(x){if(verbose){cat(".")};lasso_reg(pred,x,K=K,eps,cv.fun=cv.folds1
                                                 #,cv.fun.name=cv.fun.name
                  )} 
                  Omega[IND, which(gr %in% grpjj)]<-apply(Y,1,fun_lasso)
                  options(show.error.messages = save_show.error.messages)
                  }
                }
                if(fitfun=="SPLS"){
                  if(norm(pred,type="F")>eps){     
                    save_show.error.messages = options()$show.error.messages
                    options(show.error.messages = show.error.messages)
                    on.exit(options(show.error.messages = save_show.error.messages))
                    if(cv.subjects==TRUE){
                      cv.folds1=function(n,folds){
                        split(1:dim(pred)[2],rep(1:P,each=dim(pred)[2]/P))
                        }
                    } else {
                      cv.folds1=function(n, folds){return(split(sample(1:n), rep(1:folds, length = n)))}
                      }
                      fun_spls<-function(x){if(verbose){cat(".")};spls_reg(pred,x,K=K,eps,cv.fun=cv.folds1)} 
                      Omega[IND, which(gr %in% grpjj)]<-apply(Y,1,fun_spls)
                      options(show.error.messages = save_show.error.messages)
                  }
                }
                if(fitfun=="ELASTICNET"){
                  if(norm(pred,type="F")>eps){     
                    save_show.error.messages = options()$show.error.messages
                    options(show.error.messages = show.error.messages)
                    on.exit(options(show.error.messages = save_show.error.messages))
                    if(cv.subjects==TRUE){
                      cv.folds1=function(n,folds){
                        split(1:dim(pred)[2],rep(1:P,each=dim(pred)[2]/P))
                        }
                    }else{
                      cv.folds1=lars::cv.folds
                      }
                      fun_enet<-function(x){if(verbose){cat(".")};enet_reg(pred,x,K=K,eps,cv.fun=cv.folds1)} 
                      Omega[IND, which(gr %in% grpjj)]<-apply(Y,1,fun_enet)
                      options(show.error.messages = save_show.error.messages)
                  }
                }
                if(fitfun=="stability.c060"){
                  require(c060);
                  if(verbose){cat("mc.cores=",mc.cores,sep="")}
                  
                  fun_stab<-function(g,mc.cores=mc.cores,intercept.stabpath=intercept.stabpath){
                    
                    if(sum(pred)==0){
                      return(rep(0,nrow(pred)))
                      if(verbose){cat(".")}               
                    }else{
                  LL=rep(0,nrow(pred));error.inf=TRUE
                  try({essai<-c060::stabpath(g,t(pred),mc.cores=mc.cores,intercept=intercept.stabpath);
                  respath<-c060::stabsel(essai,error=error.stabsel,pi_thr=pi_thr.stabsel);
                  varii<-respath$stable;
                  lambda<-respath$lambda;
                  L<-glmnet::glmnet(t(pred),g,intercept=intercept.stabpath);
                  LL<-as.matrix(predict(L,s=lambda,type="coef"))[-1,1]})
                  try({LL[-varii]<-0;
                  error.inf=FALSE})
                  if(verbose){if(error.inf&options()$show.error.messages){cat("!")} else {cat(".")}}
                  if(!is.vector(LL)){LL<-rep(0,nrow(pred))}
                  return(LL)
                  }
                  }
                
                  save_show.error.messages = options()$show.error.messages
                  options(show.error.messages = show.error.messages)
                  on.exit(options(show.error.messages = save_show.error.messages))
                  Omega[IND, which(gr %in% grpjj)]<-apply(Y,1,fun_stab,mc.cores=mc.cores,intercept.stabpath=intercept.stabpath)
                  options(show.error.messages = save_show.error.messages)
                }
                
                
                
                
                
                if(fitfun=="stability.c060.weighted"){
                  require(c060);
                  if(verbose){cat("mc.cores=",mc.cores," ",sep="")}
                  priors2<-priors[IND,which(gr %in% grpjj)]
                  Y2<-cbind(1:nrow(Y),Y)
                  
                  fun_stab_weighted<-function(g,mc.cores=mc.cores,intercept.stabpath=intercept.stabpath,penalty.factor=penalty.factor){

                      if(sum(pred)==0){
                      return(rep(0,nrow(pred)))
                        if(verbose){cat(".")}               
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
                        if (length(dim(y)) == 2 | inherits(y, "Surv")) {
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
                      
                      suppressWarnings(rm(varii))
                      LL=rep(0,nrow(pred));error.comp=TRUE;error.inf=TRUE
                      try({
                      essai<-stabpath(g[-1],t(pred),mc.cores=mc.cores,intercept=intercept.stabpath,penalty.factor=priors2[,g[1]]);
                      respath<-c060::stabsel(essai,error=error.stabsel,pi_thr=pi_thr.stabsel);
                      varii<-respath$stable;
                      lambda<-respath$lambda;
                      L<-glmnet::glmnet(t(pred),g[-1],intercept=intercept.stabpath,penalty.factor=priors2[,g[1]]);
                      LL<-as.matrix(predict(L,s=lambda,type="coef"))[-1,1];
                      error.comp=FALSE
                      })
                      if(verbose){if(error.comp&options()$show.error.messages){cat("!",geterrmessage(),"\n")} else {cat("")}}
                      try({LL[-varii]<-0;error.inf=FALSE})
                      if(verbose){if(!error.comp&error.inf&options()$show.error.messages){cat("!",geterrmessage(),"\n")} else {cat(".")}}
                      if(!is.vector(LL)){LL<-rep(0,nrow(pred))}
                      return(LL)
                    }
                  }
                  
                  save_show.error.messages = options()$show.error.messages
                  options(show.error.messages = show.error.messages)
                  on.exit(options(show.error.messages = save_show.error.messages))
                  Omega[IND, which(gr %in% grpjj)]<-apply(Y2,1,fun_stab_weighted,mc.cores=mc.cores,intercept.stabpath=intercept.stabpath)
                  options(show.error.messages = save_show.error.messages)
                }
                
                
                
                if(fitfun=="robust"){
                  require(lars)
                  fun_robust<-function(g){
                    if(sum(pred)==0){
                      return(rep(0,nrow(pred)))
                      if(verbose){cat(".")}               
                    }else{
                      essai<-robustboost(t(pred)+rnorm(prod(dim(pred)),0,0.001),g)  
                      varii<-which(essai==1)
                      lambda<-0
                      L<-lars::lars(t(pred),g)
                      
                      LL<-predict(L,s=lambda,mode="lambda",type="coef")$coefficients
                      LL[-varii]<-0
                      if(verbose){cat(".")}               
                      return(LL)
                    }
                  }
                  
                  save_show.error.messages = options()$show.error.messages
                  options(show.error.messages = show.error.messages)
                  on.exit(options(show.error.messages = save_show.error.messages))
                  Omega[IND, which(gr %in% grpjj)]<-apply(Y,1,fun_robust)
                  options(show.error.messages = save_show.error.messages)
                }  
                
                if(fitfun=="selectboost.weighted"){
                  requireNamespace("SelectBoost");on.exit(unloadNamespace("package:SelectBoost"))
                  priors2<-priors[IND,which(gr %in% grpjj)]
                  Y2<-cbind(1:nrow(Y),Y)
                  
                  if(cv.subjects==TRUE){
                    folds_id_glmnet=rep(1:P,each=ncol(pred)/P)
                  } else {
                    folds_id_glmnet=sample(rep(1:K,length=ncol(pred)))
                  }
                  
                  lasso_cv_glmnet_min_weighted_Patterns <- 
                    function(X,Y,priors){
                      requireNamespace("glmnet")
                      if(is.null(priors)) priors<-rep(1,ncol(X))
                      resultat<-glmnet::cv.glmnet(X,Y,foldid=folds_id_glmnet,penalty.factor=priors)
                      coefvec<-try(as.vector(coef(resultat,s="lambda.min")[-1]))
                      # if(!is.vector(coefvec)){repu<-rep(0,ncol(X))}
                      if(!is.vector(coefvec)){coefvec<-rep(0,ncol(X))}
                      return(coefvec)
                    }
                  
                  
                  fun_selectboost_weighted<-function(g,mc.cores=mc.cores, steps.seq = steps.seq, limselect = limselect, use.parallel = use.parallel){
                    if(norm(pred,type="F")<=eps){     
                      return(rep(0,nrow(pred)))
                      if(verbose){cat(".")}               
                    }else{
                      if(exists("varii")){rm(varii)}
                      LL=rep(0,nrow(pred));error.comp=TRUE;error.inf=TRUE
                      #if(verbose){cat(intercept.stabpath)}
                      try({
                      essai<-suppressWarnings(SelectBoost::fastboost(t(pred),g[-1],SelectBoost::group_func_2,lasso_cv_glmnet_min_weighted_Patterns,corrfunc="crossprod",normalize=TRUE, B=100, use.parallel=use.parallel, ncores=mc.cores,c0lim=FALSE, steps.seq = steps.seq, priors=priors2[,g[1]]))
                      varii<-which(essai>=limselect)
                      
                      resultat<-suppressWarnings(glmnet::cv.glmnet(t(pred),g[-1],foldid=folds_id_glmnet,penalty.factor=priors2[,g[1]]))
                      LL<-predict(resultat,s="lambda.min",type="coef")[-1,1]
                      error.comp=FALSE
                      })
                      if(verbose){if(error.comp&options()$show.error.messages){cat("!",geterrmessage(),"\n")}}
                      try({LL[-varii]<-0;error.inf=FALSE})
                      if(verbose){if(!error.comp&error.inf&options()$show.error.messages){cat("!",geterrmessage(),"\n")} else {cat(".")}}
                      if(!is.vector(LL)){LL<-rep(0,nrow(pred))}
                      return(LL)
                    }
                  }
                  save_show.error.messages = options()$show.error.messages
                  options(show.error.messages = show.error.messages)
                  on.exit(options(show.error.messages = save_show.error.messages))
                  Omega[IND, which(gr %in% grpjj)]<-apply(Y2,1,fun_selectboost_weighted,mc.cores=mc.cores,steps.seq=.95,limselect=.95,use.parallel=use.parallel)
                  options(show.error.messages = save_show.error.messages)
                }  
                
                
                                
              }
              if(verbose){cat("\n")}

              co<-apply(Omega,2,sumabso)
              Omega<-t(t(Omega)/co)
              
              if(tour!=1 && type.inf=="iterative"){
                Omega<-(g(tour)*Omega+OmegaS)/(1+g(tour)) 
              }
              
              
              convO<-c(convO,mean(abs(Omega-OmegaS)))
              
              if(verbose){
                if( type.inf=="iterative"){
                  cat(paste("The convergence of the network is (L1 norm) :", round(convO[length(convO)],5)))	
                  cat("\n")	
                }
              }
              uuu<-0
              sauvF<-F
              
              if(tour==1 && type.inf=="noniterative"){
                Omega<-Omega*0+1 
              }
              
              for(grpjj in 1:ngrp){

                IND<-which(gr %in% (1:ngrp)[-(grpjj)])
                grIND<-gr[IND]
                sup_pred<-rep(1:sqF,P)+rep(seq(0,sqF*(P-1),sqF),each=sqF)
                pred<-(mat[IND,sup_pred])
                IND2<-which(gr %in% grpjj)
                
                charslist=vector("list",nF)
                Xf<-NULL
                for(i in (1:ngrp)[-grpjj]){ 
                  X<-NULL
                  suma<-function(x){sum(abs(x))}
                  f<-Vectorize(function(x){
                    apply(pred[which(grIND==i),]*Omega[IND[which(grIND==i)],x],2,suma)})
                  Xa<-(f(IND2))
                  Xb<-NULL
                  FF=Fshape[,,(i-1)*ngrp+grpjj]
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
                    Xf<-cbind(Xf,X)
                  }
                }
                
                if(!is.null(Xf)){
                  Y<-c(t(mat[IND2,sup_pred]))
                  pond<-rep(0,P)
                  coeffi<-array(0,c(P,dim(Xf)[2]))
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
                  }
                }
              } 

              if(tour==1 && type.inf=="noniterative"){
                Omega<-Omega*0 
              }
              
              
              if(type.inf=="iterative"){
                F<-(g(tour)*F+sauvF)/(1+g(tour))
              }
              cc<-rep(0,nF)
              for(i in 1:nF){cc[i]<-mean(abs((F[,,i]/sum(F[,,i])-sauvF[,,i]/sum(sauvF[,,i]))))}
              convF<-cbind(convF,cc)
              tour<-tour+1
            }
            
            if(type.inf=="iterative"){
              plot(convO[-1],type="l")
              #if(!attr(dev.cur(),"names")=="pdf"){dev.new()}
              matplot(t(convF),type="l")
            }
            else{
              F<-sauvF   
            }
            result<-new("omics_network"
                        ,omics_network=Omega
                        ,name=M@name
                        ,F=F
                        ,convF=convF
                        ,convO=convO
                        ,time_pt=M@time
            )
            return(result)
          }
)


#' Simulates omicsarray data based on a given network.
#' 
#' Simulates omicsarray data based on a given network.
#' 
#' 
#' @aliases gene_expr_simulation gene_expr_simulation-methods
#' gene_expr_simulation,omics_network-method
#' @param omics_network A omics_network object.
#' @param time_label a vector containing the time labels.
#' @param subject the number of subjects
#' @param peak_level the mean level of peaks.
#' @param act_time_group  [NULL] vector ; at which time the groups (defined by sort(unique(group))) are activated ?
#' @return A omics_array object.
#' @author Bertrand Frederic, Myriam Maumy-Bertrand.
#' @examples
#' 
#' data(Net)
#' set.seed(1)
#' 
#' #We simulate gene expressions according to the network Net
#' Msim<-Patterns::gene_expr_simulation(
#' 	omics_network=Net,
#' 	time_label=rep(1:4,each=25),
#' 	subject=5,
#' 	peak_level=200)
#' head(Msim)
#' 
setMethod("gene_expr_simulation"
          ,"omics_network"
          ,function(omics_network
                    ,time_label=1:4
                    ,subject=5
                    ,peak_level=100
                    ,act_time_group=1:4
          ){
            require(VGAM)
            
            N<-omics_network@omics_network
            M<-matrix(0,dim(omics_network@omics_network)[1],length(unique(time_label))*subject)
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
            
            MM<-as.omics_array(M,1:length(unique(time_label)),subject)
            MM@group<-time_label
            
            
            
            G<-Patterns::predict(MM,Omega=omics_network,act_time_group=act_time_group)@omicsarray_predict
            supp<-seq(1,dim(M)[2],by=length(unique(time_label)))
            G@omicsarray[,supp]<-M[, supp]	
            return(G)
          }
)
