

###generate from logistic models:
gendata.logistic2 <- function(n0,n,pi.true,logistic.fun){
  n1 <- n-n0
  z <- sample(2,n1, replace = TRUE, prob = c(1-pi.true,pi.true))#z=1:Y=0;z=2:Y=1
  y <- z-1
  Y.ind0 <- which(y==0);Y.ind1 <- which(y==1)
  size.g <- n0+length(Y.ind1);size.h <- length(Y.ind0)#sample sizes for densities g & h
  if(size.h==0){stop("sample size for the density h is 0!")}
  
  #generate size.g data with density g
  num.Y1 <- 0; Xall.G <- NULL
  while(num.Y1<size.g){
    newX <- c(runif(2),rbinom(1,1,0.6),sample.int(5,1))
    newY <- rbinom(n=1,size=1,prob=1/(1+exp(-logistic.fun(newX))))
    if(newY==1){
      Xall.G <- rbind(Xall.G,matrix(data=newX,nrow=1))
      num.Y1 <- num.Y1+1 
    }
  }
  dist.info <- list(mean.g=apply(Xall.G,2,mean),cov.g=cov(Xall.G))
  
  #generate size.h data with density h
  num.Y0 <- 0; Xall.H <- NULL
  while (num.Y0<size.h) {
    newX <- c(runif(2),rbinom(1,1,0.6),sample.int(5,1))
    newY <- rbinom(n=1,size=1,prob=1/(1+exp(-logistic.fun(newX))))
    if(newY==0){
      Xall.H <- rbind(Xall.H,matrix(data=newX,nrow=1))
      num.Y0 <- num.Y0+1 
    }
  }
  dist.info$mean.h <- apply(Xall.H,2,mean)
  dist.info$cov.h <- cov(Xall.H)
  
  pdata <- Xall.G[1:n0,,drop=F]
  udata <- rbind(Xall.H,Xall.G[(n0+1):(nrow(Xall.G)),,drop=F])
  y.udata <- c(rep(0,size.h),rep(1,length(Y.ind1)))
  
  nameX <- NULL
  for(j in 1:ncol(Xall.G)){nameX <- c(nameX,paste0("X",j))}
  colnames(pdata) <- colnames(udata) <- nameX
  pdata <- as.data.frame(pdata);udata <- as.data.frame(udata)
  return(list(pdata=pdata,udata=udata,y.udata=y.udata,dist.emp.info=dist.info))
}



###generate formula using b-splines in function gam 
gen.formula.bs2 <- function(dimC,dimL){
  out <- paste0("resp~s(X",1,",bs='bs',k=k,m=c(ord,ord.pen))")
  if(dimC>1){
    for(j in 2:dimC){
      out <- paste0(out,"+",paste0("s(X",j,",bs='bs',k=k,m=c(ord,ord.pen))"))
    }
  }
  if(dimL>0){
    for(j in 1:dimL){
      out <- paste0(out,"+",paste0("X",j+dimC))
    }
  }
  return(out)
}


###initial values 

#method3: linear fit
init.lin <- function(pdata,udata,pi.true,seed0){
  set.seed(seed0)
  pert.pi <- runif(1,min=-0.1,max=0.1) 
  pi.guess <- pi.true+pert.pi
  emout <- EM.addDR.linear(pdata,udata,maxit=1000,thres=1e-4,
                           ini.Y=pi.guess,ini.coef=NULL)#we do not specify initial coef here
  return(list(pi=emout$pi.est,linpred=emout$linpred,
              linpred.labeled=emout$linpred.labeled,conv=emout$conv))
}


#method4: initial 1-step gam fit
init.gam <- function(pdata,udata,pi.true,seed0,
                     nu=0.01,k=10,ord=3,ord.pen=2,
                     dimC=2,dimL=2){#nu and k can be reset as those in the main EM procedure
  set.seed(seed0)
  pert.pi <- runif(1,min=-0.1,max=0.1)#truncnorm::rtruncnorm(1,a=-0.2,b=0.2)
  pi.est <- pi.true+pert.pi
  n0 <- nrow(pdata);n1 <- nrow(udata);n <- n1+n0
  resp <- c(rep(1,n),rep(0,n1))#construct responses
  X <- rbind(pdata,udata,udata)#expand X
  w <- c(rep(1,n0),rep(pi.est,n1),rep(1-pi.est,n1)) #construct weights
  fm1 <- as.formula(gen.formula.bs2(dimC,dimL))
  gamout0 <- gam(formula = fm1,family = binomial(link = "logit"), 
                 data = X, weights = w, sp=rep(nu,ncol(X)),control = list(scalePenalty=F))
  conv <- gamout0$converged
  linpred <- gamout0$linear.predictors[(n0+1):n]
  linpred.labeled <- gamout0$linear.predictors[1:n0]
  return(list(pi=pi.est,linpred=linpred,linpred.labeled=linpred.labeled,conv=conv))
}


#method5: assuming no mixture, density ratio is of linear form
init.nomix.lin <- function(pdata,udata,pi.true,seed0){
  set.seed(seed0)
  pert.pi <- runif(1,min=-0.1,max=0.1)#truncnorm::rtruncnorm(1,a=-0.2,b=0.2)
  pi.est <- pi.true+pert.pi
  
  n0 <- nrow(pdata);n1 <- nrow(udata);n <- n1+n0
  resp <- c(rep(1,n0),rep(0,n1))#Z
  X <- as.matrix(rbind(pdata,udata))
  glmout0 <- glm(resp~X,family = binomial(link = "logit"))
  linpred <- glmout0$linear.predictors
  return(list(pi=pi.est,linpred=linpred[(n0+1):n],linpred.labeled=linpred[1:n0],conv=T))
}


#method6: assuming no mixture, density ratio is of additive form
init.nomix.add <- function(pdata,udata,pi.true,seed0,
                           nu=0.01,k=5,ord=3,ord.pen=2,
                           dimC=2,dimL=2){#nu and k can be reset as those in the main EM procedure
  set.seed(seed0)
  pert.pi <- runif(1,min=-0.1,max=0.1)#truncnorm::rtruncnorm(1,a=-0.2,b=0.2)
  pi.est <- pi.true+pert.pi
  
  n0 <- nrow(pdata);n1 <- nrow(udata);n <- n1+n0
  resp <- c(rep(1,n0),rep(0,n1))#Z 
  X <- rbind(pdata,udata)
  fm1 <- as.formula(gen.formula.bs2(dimC,dimL))
  gamout0 <- gam(formula = fm1,family = binomial(link = "logit"), 
                 data = X, sp=rep(nu,ncol(X)),control = list(scalePenalty=F))
  conv <- gamout0$converged
  linpred <- gamout0$linear.predictors[(n0+1):n]
  linpred.labeled <- gamout0$linear.predictors[1:n0]
  return(list(pi=pi.est,linpred=linpred,linpred.labeled=linpred.labeled,conv=conv))
}




obsLL.pen <- function(n0,n1,pi.est,linpred,linpred.labeled,gamout,dimL){
  out <- obsLL.pen0(n0,n1,pi.est,linpred,linpred.labeled)
  
  p <- length(gamout$smooth)#dimension of continuous X 
  p.para <- 1+dimL
  coef.all <- gamout$coefficients[-c(1:p.para)]
  coef.mat <- matrix(coef.all,ncol=p,byrow = F)
  quadpen <- numeric(p)
  for(j in 1:p){
    quadpen[j] <- (gamout$smooth[[j]]$sp)*(coef.mat[,j])%*%(gamout$smooth[[j]]$S[[1]])%*%(coef.mat[,j])
  }
  out <- out-sum(quadpen)
  
  return(out)
}

#a function allowing all 4 initials above (only one can be chosen)
EM.addDR.initcomp2 <- function(seed0,pdata,udata,maxit=1000,thres=1e-4,
                              ini.method=3,pi.true=0.4,
                              k=10,ord=3,ord.pen=2,nu=0.01,eval.df=NULL,
                              init=NULL,dimC=2,dimL=2){#newly added: init can be specified! 
  n0 <- nrow(pdata);n1 <- nrow(udata);n <- n1+n0
  pi.est <- error.pi <- numeric(length = 0)#record the estimated pi
  obsLL <- obsLL.incre <- numeric(length = 0)#record the observed log-lik
  resp <- c(rep(1,n),rep(0,n1))#construct responses
  X <- rbind(pdata,udata,udata)#expand X
  conv <- vector(mode = "logical", length =0)#convergence status
  
  ###initialization: 
  if(is.null(init)){
    if(ini.method==3){
      init <- init.lin(pdata,udata,pi.true=pi.true,seed0=seed0)
    }else if(ini.method==4){
      init <- init.gam(pdata,udata,pi.true=pi.true,seed0=seed0,
                       nu=nu,k=k,ord=ord,ord.pen=ord.pen,dimC=dimC,dimL=dimL)
    }else if(ini.method==5){
      init <- init.nomix.lin(pdata,udata,pi.true=pi.true,seed0=seed0)
    }else if(ini.method==6){
      init <- init.nomix.add(pdata,udata,pi.true=pi.true,seed0=seed0,
                             nu=nu,k=k,ord=ord,ord.pen=ord.pen,dimC=dimC,dimL=dimL)
    }
  }
  linPred <- init$linpred
  pi.est[1] <- init$pi
  conv[1] <- init$conv
  obsLL[1] <- -Inf
  
  ###EM algorithm
  for(r in 1:maxit){
    #E-step: 
    Y.work <- pi.est[r]/(exp(-linPred)*(n0/n1+pi.est[r])+pi.est[r])
    
    #M-step: 
    pi.est[r+1] <- mean(Y.work) #update pi
    if(pi.est[r+1]<0.045){stop(paste("too small pi occurs at iteration: ",r))}
    
    w <- c(rep(1,n0),Y.work,1-Y.work) #update weights
    fm1 <- as.formula(gen.formula.bs2(dimC=dimC,dimL=dimL))
    gamout  <- gam(formula = fm1,family = binomial(link = "logit"), 
                   data = X, weights = w, sp=rep(nu,ncol(X)),control = list(scalePenalty=F))
    conv[r+1] <- gamout$converged
    linPred <- gamout$linear.predictors[(n0+1):n]
    obsLL[r+1] <- obsLL.pen(n0,n1,pi.est[r+1],linPred,
                            gamout$linear.predictors[1:n0],gamout = gamout,
                            dimL=dimL)
    
    #check & break
    if(identical(obsLL[r+1],NaN)){
      stop(paste("NaN observed likelihood occurs at iteration: ",r))
    }
    obsLL.incre[r+1] <- obsLL[r+1]-obsLL[r]
    error.pi[r+1] <- abs((pi.est[r+1]-pi.est[r]))
    if(obsLL.incre[r+1]<thres && error.pi[r+1]<5*thres){break}
  }
  
  info.df <- data.frame(pi=pi.est,re.err.pi=error.pi,
                        obsLL=obsLL,obsLL.incre=obsLL.incre,conv=conv)
  ###classification for udata
  class.label <- pi.est[r+1]/(exp(-linPred)*(n0/n1+pi.est[r+1])+pi.est[r+1])
  
  ###functions estimation
  #set evaluation points
  if(is.null(eval.df)){
    eval.df <- matrix(nrow=101,ncol=dimC)
    colnames(eval.df) <- names(pdata)[1:dimC]#Note: continuous rvs should be put at first
    for(k in 1:dimC){
      xk <- c(pdata[,k],udata[,k])
      eval.df[,k] <- seq(from=min(xk),to=max(xk),length.out=101)
    }
    eval.df <- as.data.frame(eval.df)
  }
  
  #information cri
  obsLL0 <- obsLL.pen0(n0,n1,pi.est[r+1],linPred,gamout$linear.predictors[1:n0])#no penalty obs loglik
  obsLL0 <- -2*obsLL0
  edf <- sum(gamout$edf)
  
  #pass init value to bootstrap
  init.fin <- list(pi=pi.est[r+1],conv=ifelse(r==maxit,F,T),linpred=linPred,
                   linpred.labeled=gamout$linear.predictors[1:n0])

  
  return(list(conv=ifelse(r==maxit,F,T),
              piEst=pi.est[r+1],df=info.df,gamObj=gamout,
              label.pred=class.label,obsLL.fin=obsLL[r+1],
              aic= obsLL0+2*edf,bic=obsLL0+log(n)*edf,
              eval.pts=eval.df,
              eval.terms=predict(gamout,type="terms",newdata = eval.df),
              init=init.fin))
}



#choose best initials
EM.addDR.init4.2 <- function(seed0,pdata,udata,maxit=1000,thres=1e-4,
                             pi.true=0.4,k=10,ord=3,ord.pen=2,nu=0.01,each=5,
                             eval.df=NULL,dimC=2,dimL=2){
  #using 4 methods to generate initials 
  set.seed(seed0)
  seed.seq <- sample.int(10000000,each)
  
  n0 <- nrow(pdata);n1 <- nrow(udata);n <- n1+n0
  pi.est <- error.pi <- numeric(length = 0)#record the estimated pi
  obsLL <- obsLL.incre <- numeric(length = 0)#record the observed log-lik
  resp <- c(rep(1,n),rep(0,n1))#construct responses
  X <- rbind(pdata,udata,udata)#expand X
  
  
  test.l.4method <- list()#record the best for each initial method
  test.l.4method.val <- rep(-Inf,4)#record the best obsLL for each initial method
  for(methodNo in 3:4){
    test.l <- list();test.l.val <- rep(-Inf,each)
    for(j in 1:each){
      test <- try(EM.addDR.initcomp2(seed.seq[j],pdata,udata,maxit=maxit,thres=thres,
                                    ini.method=methodNo,pi.true=pi.true,
                                    k=k,ord=ord,ord.pen=ord.pen,nu=nu,eval.df=eval.df,
                                    init=NULL,dimC=dimC,dimL=dimL),
                  silent = T)
      if(class(test)!="try-error"){
        test.l[[j]] <- test
        test.l.val[j] <- test$obsLL.fin
      }
    }
    if(length(test.l)>0){
      id <- which.max(test.l.val)
      test.l.4method[[methodNo]] <- test.l[[id]]
      test.l.4method.val[methodNo] <- test.l.val[id]
    }
    
  }
  
  best.method <- which.max(test.l.4method.val)
  if(best.method!=1){
    test <- test.l.4method[[best.method]] 
    test$best.init.method <- best.method
    return(test)
  }else{
    return(list(conv=F))
  }
  
}




AIC.addDR2 <- function(nu.seq,seed,pdata,udata,maxit=1000,thres=1e-4,
                       pi.true=0.4,k=10,ord=3,ord.pen=2,eval.df=NULL,
                       initinput=NULL,dimC=2,dimL=2,each=5){
  aicout <- rep(NA,length(nu.seq))
  bicout <- rep(NA,length(nu.seq));edf <- rep(NA,length(nu.seq))
  pi.all <- rep(NA,length(nu.seq))
  emout.l <- list()
  
  if(is.null(initinput)){
    ###no initial input: try all initials for each tuning parameter
    for(j in 1:length(nu.seq)){
      emout <- EM.addDR.init4.2(seed,pdata,udata,maxit=maxit,thres=thres,
                                pi.true=pi.true,k=k,ord=ord,ord.pen=ord.pen,nu=nu.seq[j],
                                each=each,eval.df = eval.df,dimC=dimC,dimL=dimL)
      if(isTRUE(emout$conv)){
        aicout[j] <- emout$aic
        bicout[j] <- emout$bic; edf[j] <- sum(emout$gamObj$edf)
        pi.all[j] <- emout$piEst
        emout.l[[j]] <- emout
      }
    }
    
  }else{
    ###initial input: try one given initial for each tuning parameter
    for(j in 1:length(nu.seq)){
      emout <- try(EM.addDR.initcomp2(seed0=NULL,pdata,udata,maxit=maxit,thres=thres,
                                     ini.method=NULL,pi.true=NULL,
                                     k=k,ord=ord,ord.pen=ord.pen,nu=nu.seq[j],
                                     eval.df=eval.df,init=initinput,dimC=dimC,dimL=dimL),
                   silent = T)
      if(class(emout)!="try-error"){
        if(isTRUE(emout$conv)){
          aicout[j] <- emout$aic
          bicout[j] <- emout$bic; edf[j] <- sum(emout$gamObj$edf)
          pi.all[j] <- emout$piEst
          emout.l[[j]] <- emout
        }
      }
    }
    
  }
  id <- which.min(aicout)
  dfall <- data.frame(nu=nu.seq,aic=aicout,bic=bicout,piest=pi.all,edf=edf)
  if(length(id)>0){
    emout <- emout.l[[id]]
    out <- list(tuneid=id,tune=nu.seq[id],gamObj=emout$gamObj,
                eval.pts=emout$eval.pts,eval.terms=emout$eval.terms,
                conv=emout$conv, initMethd=emout$best.init.method,pi=emout$piEst)
    out$label.pred <- ifelse(emout$label.pred>0.5,1,0)
    out$init <- emout$init#pass linpred as future init
    out$df <- emout$df
    out$icall <- dfall
  }else{
    out <- list(conv=F)
  }
  return(out)
}




getest <- function(test,dimC,dimL,pdata){
  #test is obtained by 'AIC.addDR2'
  if(dimL>0){
    fakedata <- pdata[1,][,-c(1:dimC),drop=F]
    fakedata <- do.call(rbind,replicate(nrow(test$eval.pts), fakedata, simplify = FALSE))
    eval.pts <- cbind(test$eval.pts,fakedata)
    eval.terms <- predict(test$gamObj,type="terms",newdata = eval.pts)
    eval.terms <- eval.terms[,-c(1:dimL)]
    lincoef <- test$gamObj$coefficient[-1][1:dimL]
  }else{
    eval.terms <- predict(gamout0,type="terms",newdata = eval.df)
    lincoef <- NA
  }
  return(list(pi=test$pi,lincoef=lincoef,eval.terms=eval.terms,eval.pts=test$eval.pts))
}




###np bootstrap
npBoot.aic.core <- function(run,init,
                            nu.seq,seed,pdata,udata,maxit=1000,thres=1e-4,
                            pi.true=0.4,k=10,ord=3,ord.pen=2,eval.df=NULL,
                            dimC=2,dimL=2,each=1){
  #re-sample:
  set.seed(run)
  pdata.selec.id <- sample(c(1:nrow(pdata)),nrow(pdata),T)
  pdata.selec <- pdata[pdata.selec.id,]
  udata.selec.id <- sample(c(1:nrow(udata)),nrow(udata),T)
  udata.selec <- udata[udata.selec.id,]
  #modify initial
  getinit <- init
  getinit$linpred <- getinit$linpred[udata.selec.id]
  getinit$linpred.labeled <- getinit$linpred.labeled[pdata.selec.id]
  #aic tuning selection
  res <- AIC.addDR2(nu.seq=nu.seq,seed=seed,pdata.selec,udata.selec,maxit=maxit,thres=thres,
                    pi.true=pi.true,k=k,ord=ord,ord.pen=ord.pen,eval.df=eval.df,
                    initinput=getinit,dimC=dimC,dimL=dimL,each=each)
  if(isTRUE(res$conv)){
    return(list(pi=res$pi,tuneid=res$tuneid,icall=res$icall))
  }else{return(list(pi=NA,tuneid=NA))}
}




naive.est <- function(pdata,udata,nu=0.01,k=10,ord=3,ord.pen=2,
                      dimC=2,dimL=2){
  names.seq <- NULL
  for(j in 1:ncol(pdata)){names.seq <- c(names.seq,paste0("X",j))}
  names(pdata) <- names.seq;names(udata) <- names.seq
  
  n0 <- nrow(pdata);n1 <- nrow(udata);n <- n1+n0
  resp <- c(rep(1,n0),rep(0,n1))#Z 
  X <- rbind(pdata,udata)
  fm1 <- as.formula(gen.formula.bs2(dimC,dimL))
  gamout0 <- gam(formula = fm1,family = binomial(link = "logit"), 
                 data = X, sp=rep(nu,ncol(X)),control = list(scalePenalty=F))
  conv <- gamout0$converged
  
  eval.df <- matrix(nrow=101,ncol=dimC)
  colnames(eval.df) <- names(pdata)[1:dimC]#Note: continuous rvs should be put at first
  for(k in 1:dimC){
    xk <- c(pdata[,k],udata[,k])
    qq <- quantile(xk,probs =c(0.05,0.95))
    eval.df[,k] <- seq(from=qq[1],to=qq[2],length.out=101)
  }
  eval.df <- as.data.frame(eval.df)
  
  if(dimL>0){
    fakedata <- pdata[1,][,-c(1:dimC),drop=F]
    fakedata <- do.call(rbind,replicate(nrow(eval.df), fakedata, simplify = FALSE))
    newdata <- cbind(eval.df,fakedata)
    eval.terms <- predict(gamout0,type="terms",newdata = newdata)
    eval.terms <- eval.terms[,-c(1:dimL)]
  }else{
    eval.terms <- predict(gamout0,type="terms",newdata = eval.df)
  }
  
  return(list(conv=conv,eval.pts=eval.df,eval.terms=eval.terms,gamObj=gamout0))
}
