setGeneric("analyze_network",package="Patterns",def = function(Omega,nv,...){standardGeneric("analyze_network")})
setGeneric("clustExploration",package="Patterns",def = function(microarray,...){standardGeneric("clustExploration")})
setGeneric("clustInference",package="Patterns",def = function(microarray,vote.index,...){standardGeneric("clustInference")})
setGeneric("compare",package="Patterns",def = function(Net,Net_inf,nv){standardGeneric("compare")})
setGeneric("cutoff",package="Patterns",def = function(Omega,... ){standardGeneric("cutoff")})
setGeneric("evolution",package="Patterns",def = function(net,list_nv,... ){standardGeneric("evolution")})
setGeneric("gene_expr_simulation",package="Patterns",def = function(network,...){standardGeneric("gene_expr_simulation")})
#setGeneric("gene_counts_simulation",package="Patterns",def = function(network,...){standardGeneric("gene_counts_simulation")})
setGeneric("geneNeighborhood",package="Patterns",def = function(net,targets,... ){standardGeneric("geneNeighborhood")})
setGeneric("genePeakSelection",package="Patterns",def = function(x,peak,... ){standardGeneric("genePeakSelection")})
setGeneric("geneSelection",package="Patterns",def = function(x,y,tot.number,... ){standardGeneric("geneSelection")})
setGeneric("inference",package="Patterns",def = function(M,... ){standardGeneric("inference")})
setGeneric("position",package="Patterns",def = function(net,... ){standardGeneric("position")})
#setGeneric("inferenceCascade",package="Patterns",def = function(M,... ){standardGeneric("inferenceCascade")})
#setGeneric("predict",package="Patterns",def = function(object,...){standardGeneric("predict")})
setGeneric("probeMerge",package="Patterns",def = function(x,...){standardGeneric("probeMerge")})
setGeneric("unionMicro",package="Patterns",def = function(M1,M2 ){standardGeneric("unionMicro")})
#setGeneric("unionNGseq",package="Patterns",def = function(C1,C2 ){standardGeneric("unionNGseq")})
setGeneric("unsupervised_clustering_auto_m_c",package="Patterns",def=function(M1,... ){standardGeneric("unsupervised_clustering_auto_m_c")})
setGeneric("unsupervised_clustering",package="Patterns",def=function(M1,clust,mestim,... ){standardGeneric("unsupervised_clustering")})



as.micro_array <-
  function(M,
           time,
           subject,
           name_probe = NULL,
           gene_ID = NULL,
           group=0,
           start_time=0) {
    if (is.null(row.names(M))) {
      row.names(M) <- paste("gene", 1:dim(M)[1])
    }
    if (!is.null(name_probe)) {
      row.names(M) <- name_probe
    }
    g <- 0
    if (!is.null(gene_ID))
      g <- gene_ID
    return(
      new(
        "micro_array",
        microarray = as.matrix(M),
        name = row.names(M),
        gene_ID = g,
        time = time,
        subject = subject,
        group = group,
        start_time = start_time
      )
    )
  }

as.nextgen_seq <- function(C, time, subject) {
  if (is.null(row.names(C))) {
    row.names(C) <- paste("gene", 1:dim(C)[1])
  }
  return(
    new(
      "nextgen_seq",
      nextgenseq = as.matrix(C),
      name = row.names(C),
      time = time,
      subject = subject,
      group = 0,
      start_time = 0
    )
  )
}


cv.lars1 <- function (x, y, K = 10, index, trace = FALSE, plot.it = TRUE, 
                      se = TRUE, type = c("lasso", "lar", "forward.stagewise", 
                                          "stepwise"), mode = c("fraction", "step"), cv.fun
                      #, cv.fun.name
                      , ...) 
{
  #  requireNamespace("lars")
  #  if(trace){cat(cv.fun.name)}
  type = match.arg(type)
  if (missing(mode)) {
    mode = switch(type, lasso = "fraction", lar = "step", 
                  forward.stagewise = "fraction", stepwise = "step")
  }
  else mode = match.arg(mode)
  all.folds <- cv.fun(length(y), K)
  #  if(trace){cat(all.folds[[1]],"\n")}
  if (missing(index)) {
    index = seq(from = 0, to = 1, length = 100)
    if (mode == "step") {
      fit = lars::lars(x, y, type = type, ...)
      nsteps = nrow(fit$beta)
      maxfold = max(sapply(all.folds, length))
      nsteps = min(nsteps, length(y) - maxfold)
      index = seq(nsteps)
    }
  }
  residmat <- matrix(0, length(index), K)
  for (i in seq(K)) {
    omit <- all.folds[[i]]
    fit <- lars::lars(x[-omit, , drop = FALSE], y[-omit], trace = trace, 
                      type = type, ...)
    fit <- lars::predict.lars(fit, x[omit, , drop = FALSE], mode = mode, 
                              s = index)$fit
    if (length(omit) == 1) 
      fit <- matrix(fit, nrow = 1)
    residmat[, i] <- apply((y[omit] - fit)^2, 2, mean)
    if(trace){cat("\n CV Fold", i, "\n\n")}
  }
  cv <- apply(residmat, 1, mean)
  cv.error <- sqrt(apply(residmat, 1, var)/K)
  object <- list(index = index, cv = cv, cv.error = cv.error, 
                 mode = mode)
  if (plot.it) 
    lars::plotCVLars(object, se = se)
  invisible(object)
}


lasso_reg<-function(M,Y,K,eps,priors,cv.fun=lars::cv.folds
                    #,cv.fun.name="lars::cv.folds"
){
  #  require(lars)
  #  cat("lasso_reg",cv.fun.name,"\n")
  cvlars1<-try(cv.lars1(t(M),(Y),intercept=FALSE,K=K,plot.it=FALSE,eps=10^-5,cv.fun=cv.fun
                      #, cv.fun.name=cv.fun.name
  ))
  n<-try(cvlars1$index[which.min(cvlars1$cv)+max(1,which.min(cvlars1$cv[which.min(cvlars1$cv):length(cvlars1$cv)]<=min(cvlars1$cv)+cvlars1$cv.error[which.min(cvlars1$cv)])-1)-1])
  model<-try(lars::lars(t(M),(Y),intercept=FALSE,eps=10^-5,use.Gram=FALSE))
  repu<-try(lars::coef.lars(model,s=n,mode="fraction"))
  if(!is.vector(repu)){repu<-rep(0,dim(M)[1])}
  return(repu)
}


lasso_reg_old<-function(M,Y,K,eps,priors){
  cvlars<-try(lars::cv.lars(t(M),(Y),intercept=FALSE,K=K,plot.it=FALSE,eps=10^-5,use.Gram=FALSE))
  n<-try(cvlars$index[which.min(cvlars$cv)+max(1,which.min(cvlars$cv[which.min(cvlars$cv):length(cvlars$cv)]<=min(cvlars$cv)+cvlars$cv.error[which.min(cvlars$cv)])-1)-1])
  model<-try(lars::lars(t(M),(Y),intercept=FALSE,eps=10^-5,use.Gram=FALSE))
  repu<-try(lars::coef.lars(model,s=n,mode="fraction"))
  if(!is.vector(repu)){repu<-rep(0,dim(M)[1])}
  return(repu)
}


lasso_reg2<-function(M,Y,eps,foldid=foldid,priors){
  #require(glmnet)
  if(is.null(priors)) priors<-rep(1,nrow(M))
  cvglmnet<-try(glmnet::cv.glmnet(t(M),Y,foldid=foldid,penalty.factor=priors,intercept=FALSE))
#  n<-try(which(cvglmnet$cvm==min(cvglmnet$cvm)))
#  lambda<-try(cvglmnet$lambda[n])
#  try(print(lambda))
#  model<-try(glmnet(t(M),Y,intercept=FALSE,thresh=10^(-5),penalty.factor=priors,lambda=lambda))
#  coef(cv.fit,s="lambda.min")
#  repu<-try(as.vector(coef(model)[-1]))
  repu<-try(as.vector(glmnet::coef.cv.glmnet(cvglmnet,s="lambda.min")[-1]))
  if(!is.vector(repu)){repu<-rep(0,dim(M)[1])}
  return(repu)
}


cv.spls1 <- function (x, y, fold = 10, K, eta, kappa = 0.5, select = "pls2", 
                      fit = "simpls", scale.x = TRUE, scale.y = FALSE, plot.it = TRUE, cv.fun
                      #, cv.fun.name
                      , verbose=FALSE, ...) 
{
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  ip <- c(1:p)
  y <- as.matrix(y)
  q <- ncol(y)
  correctp1 <- function (x, y, eta, K, kappa, select, fit) 
  {
    if (min(eta) < 0 | max(eta) >= 1) {
      if (max(eta) == 1) {
        stop("eta should be strictly less than 1!")
      }
      if (length(eta) == 1) {
        stop("eta should be between 0 and 1!")
      }
      else {
        stop("eta should be between 0 and 1! \n  Choose appropriate range of eta!")
      }
    }
    if (max(K) > ncol(x)) {
      stop("K cannot exceed the number of predictors! Pick up smaller K!")
    }
    if (max(K) >= nrow(x)) {
      stop("K cannot exceed the sample size! Pick up smaller K!")
    }
    if (min(K) <= 0 | !all(K%%1 == 0)) {
      if (length(K) == 1) {
        stop("K should be a positive integer!")
      }
      else {
        stop("K should be a positive integer! \n  Choose appropriate range of K!")
      }
    }
    if (kappa > 0.5 | kappa < 0) {
      warning("kappa should be between 0 and 0.5! kappa=0.5 is used.")
      kappa <- 0.5
    }
    if (select != "pls2" & select != "simpls") {
      warning("Invalid PLS algorithm for variable selection: pls2 algorithm is used.")
      select <- "pls2"
    }
    fits <- c("simpls", "kernelpls", "widekernelpls", "oscorespls")
    if (!any(fit == fits)) {
      warning("Invalid PLS algorithm for model fitting: simpls algorithm is used.")
      fit <- "simpls"
    }
    list(K = K, eta = eta, kappa = kappa, select = select, fit = fit)
  }
  type <- correctp1(x, y, eta, K, kappa, select, fit)
  eta <- type$eta
  K <- type$K
  kappa <- type$kappa
  select <- type$select
  fit <- type$fit
  foldi <- cv.fun(n=n, folds=fold)
  mspemat <- matrix(0, length(eta), length(K))
  for (i in 1:length(eta)) {
    if(verbose){cat(paste("eta =", eta[i], "\n"))}
    mspemati <- matrix(0, fold, length(K))
    for (j in 1:fold) {
      omit <- foldi[[j]]
      object <- spls::spls(x[-omit, , drop = FALSE], y[-omit, 
                                                 , drop = FALSE], eta = eta[i], kappa = kappa, 
                     K = max(K), select = select, fit = fit, scale.x = scale.x, 
                     scale.y = scale.y, trace = FALSE)
      newx <- x[omit, , drop = FALSE]
      newx <- scale(newx, object$meanx, object$normx)
      betamat <- object$betamat
      for (k in K) {
        pred <- newx %*% betamat[[k]] + matrix(1, nrow(newx), 
                                               1) %*% object$mu
        mspemati[j, (k - min(K) + 1)] <- mean(apply((y[omit, 
                                                       ] - pred)^2, 2, mean))
      }
    }
    mspemat[i, ] <- apply(mspemati, 2, mean)
  }
  minpmse <- min(mspemat)
  rownames(mspemat) <- eta
  colnames(mspemat) <- K
  mspecol <- apply(mspemat, 2, min)
  msperow <- apply(mspemat, 1, min)
  K.opt <- min(K[mspecol == minpmse])
  eta.opt <- max(eta[msperow == minpmse])
  if(verbose){cat(paste("\nOptimal parameters: eta = ", eta.opt, ", ", sep = ""))
    cat(paste("K = ", K.opt, "\n", sep = ""))
  }
  if (plot.it) {
    spls::heatmap.spls(t(mspemat), xlab = "K", ylab = "eta", main = "CV MSPE Plot", 
                 coln = 16, as = "n")
  }
  rownames(mspemat) <- paste("eta=", eta)
  colnames(mspemat) <- paste("K =", K)
  cv <- list(mspemat = mspemat, eta.opt = eta.opt, K.opt = K.opt)
  invisible(cv)
}



spls_reg<-function(M,Y,K,eps,cv.fun=cv.fun
                   #, cv.fun.name=cv.fun.name
                   ){
                   repu=M[,1]
                   repu[]<-0
                   VarM=apply(M,1,var)
                   NonZeroVarM=!(VarM<eps)
                   NonZeroM<-M[NonZeroVarM,]
                   cvspls<-try(cv.spls1(t(NonZeroM),Y,fold=K,K=1:10,eta = seq(0.1,0.9,0.1),plot.it=FALSE,cv.fun=cv.fun
                                        #, cv.fun.name=cv.fun.name
                   ))
                   model<-try(spls::spls(t(NonZeroM),Y,eta=cvspls$eta.opt,K=cvspls$K.opt,eps=10^-5))
                   tempcoefs<-try(spls::coef.spls(model)[,1])
                   repu[names(tempcoefs)] <- tempcoefs
                   if(!is.vector(repu)){repu<-rep(0,dim(M)[1])}
                   return(repu)
                   }



spls_reg_old<-function(M,Y,K,eps){
                   cvspls<-try(spls::cv.spls(t(M),(Y),fold=K,K=1:10,eta = seq(0.1,0.9,0.1),plot.it=FALSE))
                   model<-try(spls::spls(t(M),(Y),eta=cvspls$eta.opt,K=cvspls$K.opt,eps=10^-5))
                   repu<-try(spls::coef.spls(model))
                   if(!is.vector(repu)){repu<-rep(0,dim(M)[1])}
                   return(repu)
                   }


cv.enet1 <- function (x, y, K = 10, lambda, s, mode, trace = FALSE, plot.it = TRUE, 
                      se = TRUE, cv.fun
                      #, cv.fun.name
                      , ...)  
{
  all.folds <- cv.fun(length(y), K)
  residmat <- matrix(0, length(s), K)
  for (i in seq(K)) {
    omit <- all.folds[[i]]
    fit <- elasticnet::enet(x[-omit, ], y[-omit], lambda = lambda, ...)
    fit <- elasticnet::predict.enet(fit, x[omit, , drop = FALSE], s = s, mode = mode)$fit
    if (length(omit) == 1) {
      fit <- matrix(fit, nrow = 1)
    }
    residmat[, i] <- apply((y[omit] - fit)^2, 2, mean)
    if(trace){cat("\n CV Fold", i, "\n\n")}
  }
  cv <- apply(residmat, 1, mean)
  cv.error <- sqrt(apply(residmat, 1, var)/K)
  object <- list(s = s, cv = cv, cv.error = cv.error)
  if (plot.it) {
    plot(s, cv, type = "b", xlab = mode, ylim = range(cv, 
                                                      cv + cv.error, cv - cv.error))
    if (se) 
      lars::error.bars(s, cv + cv.error, cv - cv.error, width = 1/length(s))
  }
  invisible(object)
}


enet_reg<-function(M,Y,K,eps,cv.fun=cv.fun
                   #, cv.fun.name=cv.fun.name
                   ){
  # require(elasticnet)
  M<-t(M)
  colnames(M)<-1:dim(M)[2]	  
  M2<-(M[,which(apply(M,2,sum)!=0)])
  cvenetP<-try(cv.enet1(M2,Y,lambda=0.05,s=seq(0,1,length=100),mode="fraction",intercept=FALSE,K=K,plot.it=FALSE,
                        eps=10^-5, cv.fun=cv.fun))
  cvenet<-cvenetP
  n<-try(cvenetP$s[which.min(cvenet$cv)+max(1,which.min(cvenet$cv[which.min(cvenet$cv):length(cvenet$cv)]<=min(cvenet$cv)+cvenet$cv.error[which.min(cvenet$cv)])-1)-1])	  
  modelP<-try(elasticnet::enet((M2),(Y),intercept=FALSE,eps=10^-5))
  repu2 <- elasticnet::predict.enet(modelP, s=n, type="coef", mode="fraction")	  
  repu<-rep(0,dim(M)[2])
  #    print(sum((apply(M,1,sum)!=0)))
  repu[as.numeric(names(repu2$coefficients))]<-repu2$coefficients
  if(!is.vector(repu)){repu<-rep(0,dim(M)[2])}  
  return(repu)
}


enet_reg_old<-function(M,Y,K,eps){
  # require(elasticnet)
  M<-t(M)
  colnames(M)<-1:dim(M)[2]	  
  M2<-(M[,which(apply(M,2,sum)!=0)])
  cvenetP<-try(elasticnet::cv.enet((M2),(Y),lambda=0.05,s=seq(0,1,length=100),mode="fraction",intercept=FALSE,K=K,plot.it=FALSE,eps=10^-5))
  cvenet<-cvenetP
  n<-try(cvenetP$s[which.min(cvenet$cv)+max(1,which.min(cvenet$cv[which.min(cvenet$cv):length(cvenet$cv)]<=min(cvenet$cv)+cvenet$cv.error[which.min(cvenet$cv)])-1)-1])	  
  modelP<-try(elasticnet::enet((M2),(Y),intercept=FALSE,eps=10^-5))
  repu2 <- elasticnet::predict.enet(modelP, s=n, type="coef", mode="fraction")	  
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
    
    #if(!attr(dev.cur(),"names")=="pdf"){dev.new()}
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


CascadeFinit=function(sqF,ngrp,low.trig=TRUE){
  nF=ngrp*ngrp
  Finit<-array(0,c(sqF,sqF,nF))
  for(ii in 1:nF){   
    if((ii-1)%/%ngrp+1>=ifelse(ii%%ngrp==0,ngrp,ii%%ngrp)){
      Finit[,,ii]<-0
    } else {
      if(low.trig){
        Finit[,,ii]<-lower.tri(matrix(1,sqF,sqF))*matrix(1,sqF,sqF)
      } else {
        Finit[,,ii]<-cbind(rbind(rep(0,sqF-1),diag(1,sqF-1)),rep(0,sqF))
      }
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
          Finit[,,(ii-1)*ngrp+jj]<-lower.tri(matrix(1,sqF,sqF))*matrix(1,sqF,sqF)
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


plotF <- function(x
                  ,
                  choice = "Fshape"
                  ,
                  nround = 2
                  ,
                  pixmap.color = terrain.colors(20))
{
  if (choice == "F") {
    requireNamespace("plotrix", quietly = TRUE)
    F <- x
    sizeF <- dim(F)[1]
    nF <- dim(F)[3]
    ngrp = sqrt(dim(F)[3])
    ymax <- max(F)
    coloring <- rainbow(ngrp)

    FF = NULL
    for (i in 1:ngrp) {
      FFa = NULL
      for (j in 1:ngrp) {
        FFa = cbind(FFa, round(F[, , (i - 1) * ngrp + j], nround))
      }
      FF = rbind(FF, FFa)
    }
    par(
      mar = c(0, 0, 0, 0),
      oma = c(0, 0, 0, 0),
      mai = c(0, 0, 0, 0)
    )
    plotrix::color2D.matplot(
      x = FF,
      cs1 = c(0, .5, 1),
      cs2 = c(.5, 1, 0),
      cs3 = c(1, 2, 0),
      show.values = nround,
      axes = FALSE,
      main = "",
      xlab = "",
      ylab = "",
      show.legend = FALSE
    )
    abline(h = sizeF * (1:(ngrp - 1)), lwd = 3, col="grey50")
    abline(v = sizeF * (1:(ngrp - 1)), lwd = 3, col="grey50")
  }
  if (choice == "Fpixmap") {
    F <- x
    nF <- dim(F)[3]
    ngrp = sqrt(dim(F)[3])
    ymax <- max(F)
    coloring <- rainbow(ngrp)

    FF = NULL
    for (i in 1:ngrp) {
      FFa = NULL
      for (j in 1:ngrp) {
        FFa = cbind(FFa, round(F[, , (i - 1) * ngrp + j], nround))
      }
      FF = rbind(FF, FFa)
    }
    par(
      mar = c(0, 0, 0, 0),
      oma = c(0, 0, 0, 0),
      mai = c(0, 0, 0, 0)
    )
    x.temp <- pixmap::pixmapGrey(data = FF, cellres = c(2, 2))
    plot(x.temp)
    rm(x.temp)
  }
  if (choice == "Fshape") {
    requireNamespace("plotrix", quietly = TRUE)
    Fshape <- x
    sizeF <- dim(Fshape)[1]
    nF <- dim(Fshape)[3]
    ngrp = sqrt(dim(Fshape)[3])
    ymax <- max(Fshape)
    nround = 0
    
    Fshape[Fshape == "0"] <- NA
    
    allcoefFtheo = na.omit(unique(as.vector(Fshape)))
    ncoefFtheo = length(allcoefFtheo)
    coloring <- rainbow(ncoefFtheo)
    
    Ftemp = array(NA, c(nrow(Fshape), ncol(Fshape), nF))
    for (coefFtheo in 1:ncoefFtheo) {
      Ftemp[Fshape == allcoefFtheo[coefFtheo]] <- coefFtheo
    }
    
    FF = NULL
    FFc = NULL
    for (i in 1:ngrp) {
      FFa = NULL
      FFcol = NULL
      for (j in 1:ngrp) {
        FFa = cbind(FFa, Ftemp[, , (i - 1) * ngrp + j])
        FFcol = cbind(FFcol, apply(Ftemp[, , (i - 1) * ngrp + j],c(1,2),function(x) {coloring[x]}))
      }
      FF = rbind(FF, FFa)
      FFc = rbind(FFc, FFcol)
    }
    
    FFa[is.na(FFa)]<-0
    FFcol[is.na(FFcol)]<-"#000000"
    
    par(
      mar = c(0, 0, 0, 0),
      oma = c(0, 0, 0, 0),
      mai = c(0, 0, 0, 0)
    )
    plotrix::color2D.matplot(
      x = FF,
      cs1 = c(0, .5, 1),
      cs2 = c(.5, 1, 0),
      cs3 = c(1, 2, 0),
      show.values = nround,
      axes = FALSE,
      main = "",
      xlab = "",
      ylab = "",
      show.legend = FALSE,
      cellcolors = FFc
    )
    abline(h = sizeF * (1:(ngrp - 1)), lwd = 3, col="grey50")
    abline(v = sizeF * (1:(ngrp - 1)), lwd = 3, col="grey50")
  }
  if (choice == "Fshapepixmap") {
    Fshape <- x
    nF <- dim(Fshape)[3]
    ngrp = sqrt(dim(Fshape)[3])
    ymax <- max(Fshape)
    
    allcoefFtheo = unique(as.vector(Fshape))
    ncoefFtheo = length(allcoefFtheo)
    coloring <- rainbow(ncoefFtheo)
    
    Ftemp = array(NA, c(nrow(Fshape), ncol(Fshape), nF))
    for (coefFtheo in 1:ncoefFtheo) {
      Ftemp[Fshape == allcoefFtheo[coefFtheo]] <- coefFtheo
    }
    
    FF = NULL
    for (i in 1:ngrp) {
      FFa = NULL
      for (j in 1:ngrp) {
        FFa = cbind(FFa, round(Ftemp[, , (i - 1) * ngrp + j], nround))
      }
      FF = rbind(FF, FFa)
    }
    par(
      mar = c(0, 0, 0, 0),
      oma = c(0, 0, 0, 0),
      mai = c(0, 0, 0, 0)
    )
    x.temp <- pixmap::pixmapGrey(data = FF, cellres = c(2, 2))
    plot(x.temp)
    rm(x.temp)
  }
}
















robustboost <- function(X,
                  Y,
                  corr = 0.8,
                  B = 100,
                  normalize = TRUE,
                  eps = 10 ^ (-4)) {
  func <- "selec_meth"
  selec_meth <- function(X, Y) {
    N <- dim(X)[1]
    Nc <- dim(X)[2]
    modd <- lars::lars(X, Y)
    tau2 <- var(Y)
    index <-
      which.min(N * log(2 * pi * tau2) + apply((
        lars::predict.lars(modd, newx = X)$fit - Y
      ) ^ 2, 2, sum) / tau2 + log(N) * (apply(
        lars::predict.lars(modd, type = "coef")$coef != 0, 1, sum
      ) + 1))
    lars::predict.lars(modd, s = index, type = "coef")$coefficients
    
  }
  group <- function(X, c0) {
    Xcor <- abs(cor(X))
    Xcor[Xcor < c0] <- 0
    Xcor[Xcor != 0] <- 1
    
    dete <- function(x) {
      which(x == 1)
    }
    res <- apply(Xcor, 2, dete)
    return(res)
    
    
  }
  
  group2 <- function(x)
    group(x, 1)
  
  N <- dim(X)[1]
  Nc <- dim(X)[2]
  
  if (normalize == TRUE) {
    X <-
      X - matrix(rep(apply(X, 2, mean), N), N, Nc, byrow =
                   TRUE)
    X <- t(t(X) / sqrt(diag(t(X) %*% X)))
    
  }
  
  Correlation_matrice <- t(X) %*% X
  Correlation_indice <- Correlation_matrice * 0
  Correlation_indice[which((Correlation_matrice) >= 0)] <- 1
  Correlation_indice[which(-(Correlation_matrice) >= 0)] <- -1
  diag(Correlation_indice) <- 1
  
  groups <- group2(X)
  
  Boot <- matrix(rep(0, Nc * B), B, Nc)

  Xpass <- matrix(0, N, N - 1)
  Xpass[col(Xpass) > row(Xpass)] <- 1
  diag(Xpass) <- 1
  Xpass[col(Xpass) == (row(Xpass) - 1)] <- -(1:(N - 1))
  Xpass <- t(t(Xpass) / sqrt(diag(t(Xpass) %*% Xpass)))
  olsXpass <- solve(t(Xpass) %*% Xpass) %*% t(Xpass)
  nXpass <- t(t(Xpass) / sqrt(diag(t(Xpass) %*% Xpass)))
  
  func_passage1 <- function(x, olsXpass = olsXpass) {
    return(olsXpass %*% x)
  }
  
  func_passage2 <- function(x, N = N, nXpass = nXpass) {
    return(apply(matrix(
      rep(x, N), N, N - 1, byrow = TRUE
    ) * nXpass, 1, sum))
  }
  
  Xret <- array(0, c(N, Nc, B))
  Xb <- X
  
  
  
  simul1 <- function(j) {
    if (length(groups[[j]]) >= 2) {
      indice <- groups[[j]]
      corr_set <- t(t(X[, indice]) * Correlation_indice[j, indice])
      corr_set2 <- apply(corr_set, 2, func_passage1)
      BB <- coef(movMF::movMF(t(corr_set2), 1))$theta
      newv <- movMF::rmovMF(1, BB)
      newv <- func_passage2(newv)
    } else{
      newv <- X[, j]
      
    }
    
    
    return(newv)
  }
  
  simul1 <- Vectorize(simul1)
  
  simul2 <- function() {
    simul1(1:(Nc))
  }
  
  Xret <- replicate(B, simul2())
  
  
  select <- function(k) {
    rr <- eval(parse(text = paste(func, "(Xret[,,", k, "],Y)", sep = "")))
    
    return(rr)
  }
  
  
  Boot <- lapply(1:B, select)
  
  Boot <- matrix(unlist(Boot), B, Nc, byrow = TRUE)
  
  Fs <- apply(abs(Boot) > eps, 2, sum)
  
  
  return(Fs / B)
}
