############functions for generating PU data for simulation

###generate data directly from g and h 
gendata.normal <- function(n0,n,pi.true,
                           mean1=rep(1,5),mean0=rep(0,5),
                           sigma1=diag(1,5),sigma0=diag(1,5)){
  n1 <- n-n0
  z <- sample(2,n1, replace = TRUE, prob = c(1-pi.true,pi.true))#z=1:Y=0;z=2:Y=1
  y <- z-1
  Y.ind0 <- which(y==0);Y.ind1 <- which(y==1)
  pdata <- mvtnorm::rmvnorm(n0,mean=mean1,sigma=sigma1)
  udata <- matrix(nrow=n1,ncol=length(mean0))
  udata[Y.ind0,] <- mvtnorm::rmvnorm(length(Y.ind0),mean=mean0,sigma=sigma0)
  udata[Y.ind1,] <- mvtnorm::rmvnorm(length(Y.ind1),mean=mean1,sigma=sigma1)
  nameX <- NULL
  for(j in 1:length(mean1)){nameX <- c(nameX,paste0("X",j))}
  colnames(pdata) <- colnames(udata) <- nameX
  pdata <- as.data.frame(pdata);udata <- as.data.frame(udata)
  return(list(pdata=pdata,udata=udata,y.udata=y))
}


gendata.logistic <- function(n0,n,pi.true,logistic.fun,dim.X){
  n1 <- n-n0
  z <- sample(2,n1, replace = TRUE, prob = c(1-pi.true,pi.true))#z=1:Y=0;z=2:Y=1
  y <- z-1
  Y.ind0 <- which(y==0);Y.ind1 <- which(y==1)
  size.g <- n0+length(Y.ind1);size.h <- length(Y.ind0)#sample sizes for densities g & h
  if(size.h==0){stop("sample size for the density h is 0!")}
  
  #generate size.g data with density g
  num.Y1 <- 0; Xall.G <- NULL
  while(num.Y1<size.g){
    newX <- runif(dim.X)
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
    newX <- runif(dim.X)
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
  for(j in 1:dim.X){nameX <- c(nameX,paste0("X",j))}
  colnames(pdata) <- colnames(udata) <- nameX
  pdata <- as.data.frame(pdata);udata <- as.data.frame(udata)
  return(list(pdata=pdata,udata=udata,y.udata=y.udata,dist.emp.info=dist.info))
}



