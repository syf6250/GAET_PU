############functions for parametric linear density ratio model

############parametric linear model
######obtain initial values for parametric linear fit:
EM.addDR.linear.init <- function(pdata,udata,ini.Y=0.1,ini.coef=NULL){
  ###NOTE: initial \pi should be given (ini.Y); otherwise it is set as 0.1
  n0 <- nrow(pdata);n1 <- nrow(udata);n <- n1+n0
  data <- rbind(pdata,udata)
  
  if(!is.null(ini.coef)){#if coefficients are provided already 
    if(length(ini.coef)==(ncol(pdata)+1)){
      out.coef <- as.numeric(ini.coef)
      linpred <- as.matrix(data)%*%out.coef[-1]+out.coef[1]
      linpred <- as.numeric(linpred)+log((n0/n1+ini.Y)/(1-ini.Y))#linear predictor of logistic model
    }else{
      stop("The dimension of ini.coef should match dimension of covariates +1")
    }
  }else{#when only an initial guess of \pi is given; (0.1 by default)
    resp <- c(rep(1,n),rep(0,n1))#construct responses
    X <- as.matrix(rbind(pdata,udata,udata))#expand X
    w <- c(rep(1,n0),rep(ini.Y,n1),rep(1-ini.Y,n1)) #construct weights
    glmout0 <- glm(resp~X,family = binomial(link = "logit"),weights = w)#initial estimation 
    out.coef <- as.numeric(coef(glmout0))
    out.coef[1] <- out.coef[1]-log((n0/n1+ini.Y)/(1-ini.Y))
    linpred <- glmout0$linear.predictors
  }
  return(list(pi=ini.Y,coef=out.coef,
              linpred=linpred[(n0+1):n],linpred.labeled=linpred[1:n0]))
}

######EM: parametric linear method
EM.addDR.linear <- function(pdata,udata,maxit=500,thres=1e-4,
                            ini.Y=0.1,ini.coef=NULL){
  n0 <- nrow(pdata);n1 <- nrow(udata);n <- n1+n0
  pi.est <- error.pi <- numeric(length = 0)#record the estimated pi
  obsLL <- obsLL.incre <- numeric(length = 0)#record the observed log-lik
  resp <- c(rep(1,n),rep(0,n1))#construct responses
  X <- as.matrix(rbind(pdata,udata,udata))#expand X
  
  
  ###initialization: 
  init <- EM.addDR.linear.init(pdata,udata,ini.Y=ini.Y,ini.coef=ini.coef)
  linPred <- init$linpred
  pi.est[1] <- init$pi
  #compute obs log-lik
  obsLL[1] <- obsLL.pen0(n0,n1,pi.est[1],init$linpred,init$linpred.labeled)
  
  
  ###EM algorithm
  for(r in 1:maxit){
    #E-step: 
    Y.work <- pi.est[r]/(exp(-linPred)*(n0/n1+pi.est[r])+pi.est[r])
    
    #M-step: 
    pi.est[r+1] <- mean(Y.work) #update pi
    w <- c(rep(1,n0),Y.work,1-Y.work) #update weights
    glmout <- glm(resp~X,family = binomial(link = "logit"),weights = w)
    linPred <- glmout$linear.predictors[(n0+1):n]
    obsLL[r+1] <- obsLL.pen0(n0,n1,pi.est[r+1],linPred,
                             glmout$linear.predictors[1:n0])
    
    #check & break
    if(identical(obsLL[r+1],NaN)){
      #stop(paste("NaN observed likelihood occurs at iteration: ",r))
      return(list(conv=F,pi.est=pi.est[r+1],NaNoccur=T,obsLL.fin=obsLL[r],
                  label.pred=pi.est[r+1]/(exp(-linPred)*(n0/n1+pi.est[r+1])+pi.est[r+1])))
    }
    obsLL.incre[r+1] <- obsLL[r+1]-obsLL[r]
    error.pi[r+1] <- abs((pi.est[r+1]-pi.est[r]))#/pi.est[r])
    if(obsLL.incre[r+1]<thres && error.pi[r+1]<5*thres){break}
  }
  ###classification for udata
  class.label <- pi.est[r+1]/(exp(-linPred)*(n0/n1+pi.est[r+1])+pi.est[r+1])
  
  info.df <- data.frame(pi=pi.est,re.err.pi=error.pi,
                        obsLL=obsLL,obsLL.incre=obsLL.incre)
  coef.est <- glmout$coefficients 
  coef.est[1] <- coef.est[1]-log((n0/n1+pi.est[r+1])/(1-pi.est[r+1]))
  linpred <- glmout$linear.predictors
  
  return(list(conv=ifelse(r==maxit,F,T),
              df=info.df,glmObj=glmout,NaNoccur=F,
              coef.est=coef.est,pi.est=pi.est[r+1],
              linpred=linpred[(n0+1):n],linpred.labeled=linpred[1:n0],
              label.pred=class.label,obsLL.fin=obsLL[r+1]))
}
######observed log-likelihood 
obsLL.pen0 <- function(n0,n1,pi.est,linpred,linpred.labeled){
  ratio <- n0/(n0+n1*pi.est)
  part1 <- log(ratio/(1+exp(-linpred.labeled)))
  exptem <- exp(-linpred)
  part2 <- log((exptem+1-ratio)/(1+exptem))
  return(sum(part1)+sum(part2))
}





EM.addDR.linear.minit <- function(seed0,pdata,udata,maxit=500,thres=1e-4,
                                  pi.true=0.4,ini.coef=NULL,each=5){
  set.seed(seed0)
  pert.pi <- runif(5,min=-0.1,max=0.1)
  pi.seq <- pi.true+pert.pi
  
  test.l <- list();test.l.val <- rep(-Inf,each)
  test.l.val2 <- rep(-Inf,each);test.l.2 <- list()#for NaN cases
  for(j in 1:each){
    paraout <- try(EM.addDR.linear(pdata,udata,maxit=maxit,thres=thres,
                                   ini.Y=pi.seq[j],ini.coef=ini.coef),silent = T)
    # if(class(paraout)!="try-error"){
    #   test.l[[j]] <- paraout
    #   test.l.val[j] <- paraout$obsLL.fin
    # }
    if(isFALSE(paraout$NaNoccur)){
      test.l[[j]] <- paraout
      test.l.val[j] <- paraout$obsLL.fin
    }else{
      test.l.2[[j]] <- paraout
      test.l.val2[j] <- paraout$obsLL.fin
    }
  }
  
  if(length(test.l)>0){
    id <- which.max(test.l.val)
    return(test.l[[id]])
  }else{
    #return(list(conv=F))
    id <- which.max(test.l.val2)
    return(test.l.2[[id]])
  }
}

