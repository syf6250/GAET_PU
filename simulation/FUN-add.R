############some basic functions for the proposed method

######generate formula using b-splines in function gam 
gen.formula.bs <- function(dim.X){
  out <- paste0("resp~s(X",1,",bs='bs',k=k,m=c(ord,ord.pen))")
  if(dim.X>1){
    for(j in 2:dim.X){
      out <- paste0(out,"+",paste0("s(X",j,",bs='bs',k=k,m=c(ord,ord.pen))"))
    }
  }
  return((out))
}
######generate knots according to data quantiles.
gen.q.knots <- function(x,k,ord){
  #k: number of B-spline basis
  #ord: order of B-spline
  #generate inner knots:
  kts <- quantile(x,probs = seq(0,1,length.out=k-ord+1),names=F)
  #generate outer knots:
  dif <- diff(range(x))/(k-ord)
  for(jj in 1:ord){
    kts <- c(kts[1]-dif,kts,kts[length(kts)]+dif)
  }
  return(kts)
}

######observed log-likelihood (penalized)
obsLL.pen <- function(n0,n1,pi.est,linpred,linpred.labeled,gamout){
  out <- obsLL.pen0(n0,n1,pi.est,linpred,linpred.labeled)
  
  p <- length(gamout$smooth)#dimension of X(since each coordinate is smoothed individually)
  coef.all <- gamout$coefficients[-1]
  coef.mat <- matrix(coef.all,ncol=p,byrow = F)
  quadpen <- numeric(p)
  for(j in 1:p){
    quadpen[j] <- (gamout$smooth[[j]]$sp)*(coef.mat[,j])%*%(gamout$smooth[[j]]$S[[1]])%*%(coef.mat[,j])
  }
  out <- out-sum(quadpen)
  
  return(out)
}



############initializations


###method1: true value
init0 <- function(pdata,udata,pi.true,PrY1,logistic.fun){
  n0 <- nrow(pdata);n1 <- nrow(udata)
  shift <- log((n0/n1+pi.true)/(1-pi.true))
  linpred <- apply(udata, 1, logistic.fun)+log((1-PrY1)/PrY1)#true log density ratio
  linpred <- linpred+shift
  linpred.labeled <- apply(pdata, 1, logistic.fun)+log((1-PrY1)/PrY1)+shift
  return(list(pi=pi.true,linpred=linpred,linpred.labeled=linpred.labeled,
              conv=T))
}


###method2: true value + small perturbation 
init0.ptb <- function(pdata,udata,pi.true,seed0,logistic.fun){
  p <- ncol(pdata)
  set.seed(seed0)
  pert <- truncnorm::rtruncnorm(p,a=-0.5,b=0.5); pert.pi <- runif(1,min=-0.2,max=0.2)#truncnorm::rtruncnorm(1,a=-0.2,b=0.2)
  #pert <- runif(dim.X,min=-0.5,max=0.5); pert.pi <- runif(1,min=-0.2,max=0.2)
  logistic.fun2 <- function(x){
    logistic.fun(x)+as.numeric(crossprod(pert,x))
  }
  #monte carlo method to see Pr(Y=1) (population level)
  sum <- 0
  for(j in 1:500000){
    sum <- sum + 1/(1+exp(-logistic.fun2(runif(p))))
  }
  return(init0(pdata,udata,pi.true+pert.pi,sum/500000,logistic.fun2))
}


###method3: linear fit
init.lin <- function(pdata,udata,pi.true,seed0){
  set.seed(seed0)
  pert.pi <- runif(1,min=-0.2,max=0.2)#truncnorm::rtruncnorm(1,a=-0.2,b=0.2)
  pi.guess <- pi.true+pert.pi
  emout <- EM.addDR.linear(pdata,udata,maxit=1000,thres=1e-4,
                           ini.Y=pi.guess,ini.coef=NULL)
  return(list(pi=emout$pi.est,linpred=emout$linpred,
              linpred.labeled=emout$linpred.labeled,conv=emout$conv))
}


###method4: initial 1-step gam fit
init.gam <- function(pdata,udata,pi.true,seed0,
                     nu=0.1,k=5,ord=3,ord.pen=2){
  set.seed(seed0)
  pert.pi <- runif(1,min=-0.2,max=0.2)#truncnorm::rtruncnorm(1,a=-0.2,b=0.2)
  pi.est <- pi.true+pert.pi
  n0 <- nrow(pdata);n1 <- nrow(udata);n <- n1+n0
  resp <- c(rep(1,n),rep(0,n1))#construct responses
  X <- rbind(pdata,udata,udata)#expand X
  w <- c(rep(1,n0),rep(pi.est,n1),rep(1-pi.est,n1)) #construct weights
  fm1 <- as.formula(gen.formula.bs(ncol(pdata)))
  gamout0 <- gam(formula = fm1,family = binomial(link = "logit"), 
                 data = X, weights = w, sp=rep(nu,ncol(X)),control = list(scalePenalty=F))
  conv <- gamout0$converged
  linpred <- gamout0$linear.predictors[(n0+1):n]
  linpred.labeled <- gamout0$linear.predictors[1:n0]
  return(list(pi=pi.est,linpred=linpred,linpred.labeled=linpred.labeled,conv=conv))
}



###method5: assuming no mixture, density ratio is of linear form
init.nomix.lin <- function(pdata,udata,pi.true,seed0){
  set.seed(seed0)
  pert.pi <- runif(1,min=-0.2,max=0.2)#truncnorm::rtruncnorm(1,a=-0.2,b=0.2)
  pi.est <- pi.true+pert.pi
  
  n0 <- nrow(pdata);n1 <- nrow(udata);n <- n1+n0
  resp <- c(rep(1,n0),rep(0,n1))#Z
  X <- as.matrix(rbind(pdata,udata))
  glmout0 <- glm(resp~X,family = binomial(link = "logit"))
  linpred <- glmout0$linear.predictors
  return(list(pi=pi.est,linpred=linpred[(n0+1):n],linpred.labeled=linpred[1:n0],conv=T))
}



###method6: assuming no mixture, density ratio is of additive form
init.nomix.add <- function(pdata,udata,pi.true,seed0,
                           nu=0.1,k=5,ord=3,ord.pen=2){
  set.seed(seed0)
  pert.pi <- runif(1,min=-0.2,max=0.2)#truncnorm::rtruncnorm(1,a=-0.2,b=0.2)
  pi.est <- pi.true+pert.pi
  
  n0 <- nrow(pdata);n1 <- nrow(udata);n <- n1+n0
  resp <- c(rep(1,n0),rep(0,n1))#Z 
  X <- rbind(pdata,udata)
  fm1 <- as.formula(gen.formula.bs(ncol(pdata)))
  gamout0 <- gam(formula = fm1,family = binomial(link = "logit"), 
                 data = X, sp=rep(nu,ncol(X)),control = list(scalePenalty=F))
  conv <- gamout0$converged
  linpred <- gamout0$linear.predictors[(n0+1):n]
  linpred.labeled <- gamout0$linear.predictors[1:n0]
  return(list(pi=pi.est,linpred=linpred,linpred.labeled=linpred.labeled,conv=conv))
}








###a function allowing all initials above (only one can be chosen)
EM.addDR.initcomp <- function(seed0,pdata,udata,maxit=1000,thres=1e-4,
                              ini.method=1,pi.true=0.4,PrY1=0.424,logistic.fun,
                              k=10,ord=3,ord.pen=2,nu=0.1,eval.df=NULL,
                              init=NULL){
  n0 <- nrow(pdata);n1 <- nrow(udata);n <- n1+n0
  pi.est <- error.pi <- numeric(length = 0)#record the estimated pi
  obsLL <- obsLL.incre <- numeric(length = 0)#record the observed log-lik
  resp <- c(rep(1,n),rep(0,n1))#construct responses
  X <- rbind(pdata,udata,udata)#expand X
  conv <- vector(mode = "logical", length =0)#convergence status
  
  ###initialization: 
  if(is.null(init)){
    if(ini.method==1){
      init <- init0(pdata,udata,pi.true=pi.true,PrY1=PrY1,logistic.fun=logistic.fun)
    }else if(ini.method==2){
      init <- init0.ptb(pdata,udata,pi.true=pi.true,seed0=seed0,logistic.fun=logistic.fun)
    }else if(ini.method==3){
      init <- init.lin(pdata,udata,pi.true=pi.true,seed0=seed0)
    }else if(ini.method==4){
      init <- init.gam(pdata,udata,pi.true=pi.true,seed0=seed0,
                       nu=nu,k=k,ord=ord,ord.pen=ord.pen)
    }else if(ini.method==5){
      init <- init.nomix.lin(pdata,udata,pi.true=pi.true,seed0=seed0)
    }else if(ini.method==6){
      init <- init.nomix.add(pdata,udata,pi.true=pi.true,seed0=seed0,
                             nu=nu,k=k,ord=ord,ord.pen=ord.pen)
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
    w <- c(rep(1,n0),Y.work,1-Y.work) #update weights
    fm1 <- as.formula(gen.formula.bs(ncol(pdata)))
    gamout  <- gam(formula = fm1,family = binomial(link = "logit"), 
                   data = X, weights = w, sp=rep(nu,ncol(X)),control = list(scalePenalty=F))
    conv[r+1] <- gamout$converged
    linPred <- gamout$linear.predictors[(n0+1):n]
    obsLL[r+1] <- obsLL.pen(n0,n1,pi.est[r+1],linPred,
                            gamout$linear.predictors[1:n0],gamout = gamout)
    
    #check & break
    if(identical(obsLL[r+1],NaN)){
      stop(paste("NaN observed likelihood occurs at iteration: ",r))
    }
    obsLL.incre[r+1] <- obsLL[r+1]-obsLL[r]
    error.pi[r+1] <- abs((pi.est[r+1]-pi.est[r])/pi.est[r])
    if(obsLL.incre[r+1]<thres && error.pi[r+1]<10*thres){break}
  }
  
  info.df <- data.frame(pi=pi.est,re.err.pi=error.pi,
                        obsLL=obsLL,obsLL.incre=obsLL.incre,conv=conv)
  ###classification for udata
  class.label <- pi.est[r+1]/(exp(-linPred)*(n0/n1+pi.est[r+1])+pi.est[r+1])
  
  ###functions estimation
  #set evaluation points
  if(is.null(eval.df)){
    eval.df <- matrix(nrow=101,ncol=ncol(pdata))
    colnames(eval.df) <- names(pdata)
    for(k in 1:ncol(pdata)){
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
  init.fin <- list(linpred=linPred)
  #init.fin$linpred.labeled <- gamout$linear.predictors[1:n0]
  
  return(list(conv=ifelse(r==maxit,F,T),
              piEst=pi.est[r+1],df=info.df,gamObj=gamout,
              label.pred=class.label,obsLL.fin=obsLL[r+1],
              aic= obsLL0+2*edf,bic=obsLL0+log(n)*edf,
              eval.pts=eval.df,edf=edf,
              eval.terms=predict(gamout,type="terms",newdata = eval.df),
              init=init.fin))
}





############functions that supports multiple initializations
###a function allowing four initial methods together
EM.addDR.init4 <- function(seed0,pdata,udata,maxit=1000,thres=1e-4,
                           pi.true=0.4,k=10,ord=3,ord.pen=2,nu=0.01,each=5,eval.df=NULL){
  #using 4 methods to generate initials 
  set.seed(seed0)
  seed.seq <- sample.int(10000000,each)
  
  n0 <- nrow(pdata);n1 <- nrow(udata);n <- n1+n0
  pi.est <- error.pi <- numeric(length = 0)#record the estimated pi
  obsLL <- obsLL.incre <- numeric(length = 0)#record the observed log-lik
  resp <- c(rep(1,n),rep(0,n1))#construct responses
  X <- rbind(pdata,udata,udata)#expand X
  
  
  test.l.4method <- list()#record the best for each initial method
  test.l.4method.val <- rep(-Inf,6)#record the best obsLL for each initial method
  for(methodNo in 3:6){
    test.l <- list();test.l.val <- rep(-Inf,each)
    for(j in 1:each){
      test <- try(EM.addDR.initcomp(seed.seq[j],pdata,udata,maxit=maxit,thres=thres,
                                    ini.method=methodNo,pi.true=pi.true,
                                    PrY1=NULL,logistic.fun=NULL,
                                    k=k,ord=ord,ord.pen=ord.pen,nu=nu,eval.df=eval.df),
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



######AIC 
AIC.addDR <- function(nu.seq,seed,pdata,udata,maxit=1000,thres=1e-4,
                      pi.true=0.4,k=10,ord=3,ord.pen=2,eval.df=NULL,
                      Fvalout=T,initinput=NULL){
  aicout <- rep(NA,length(nu.seq))
  emout.l <- list()
  edf <- rep(NA,length(nu.seq))
  pi.all <- rep(NA,length(nu.seq))
  
  if(is.null(initinput)){
    ###no initial input: try all initials for each tuning parameter
    for(j in 1:length(nu.seq)){
      emout <- EM.addDR.init4(seed,pdata,udata,maxit=maxit,thres=thres,
                              pi.true=pi.true,k=k,ord=ord,ord.pen=ord.pen,nu=nu.seq[j],
                              each=5,eval.df = eval.df)
      if(isTRUE(emout$conv)){
        aicout[j] <- emout$aic
        edf[j] <- sum(emout$gamObj$edf)
        pi.all[j] <- emout$piEst
        emout.l[[j]] <- emout
      }
    }
    
  }else{
    ###initial input: try one given initial for each tuning parameter
    for(j in 1:length(nu.seq)){
      emout <- try(EM.addDR.initcomp(seed0=NULL,pdata,udata,maxit=maxit,thres=thres,
                                     ini.method=NULL,pi.true=NULL,PrY1=NULL,logistic.fun=NULL,
                                     k=k,ord=ord,ord.pen=ord.pen,nu=nu.seq[j],
                                     eval.df=eval.df,init=initinput),silent = T)
      if(class(emout)!="try-error"){
        if(isTRUE(emout$conv)){
          aicout[j] <- emout$aic
          edf[j] <- sum(emout$gamObj$edf)
          pi.all[j] <- emout$piEst
          emout.l[[j]] <- emout
        }
      }
    }
  }
  
  id <- which.min(aicout)
  if(length(id)>0){
    emout <- emout.l[[id]]
    out <- list(tuneid=id,tune=nu.seq[id],
                conv=emout$conv, initMethd=emout$best.init.method,pi=emout$piEst)
    out$label.pred <- ifelse(emout$label.pred>0.5,1,0)
    out$df <- emout$df
    out$icall <- data.frame(nu=nu.seq,aic=aicout,piest=pi.all,edf=edf)
    out$init <- emout$init#pass linpred as future init
    out$gamObj <- emout$gamObj
    if(isTRUE(Fvalout)){
      out$eval.terms <- emout$eval.terms
    }
  }else{
    out <- list(conv=F)
  }
  return(out)
}


