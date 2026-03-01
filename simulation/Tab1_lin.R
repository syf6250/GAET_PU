###simulation for linear method in Table 1

library(doParallel)
source("FUN-gendata.R")
source("FUN-para.R")

pi.true <- 0.4
RUN <- 2000
prop <- 5
n1seq <- c(250,500,750,1250)
n0seq <- prop*n1seq
nseq <- n0seq+n1seq

para_Linear <- function(run,n0,n,pi.true=0.4,logistic.fun,dim.X=5){
  seed <- 77*run+7
  set.seed(seed)
  alldata <- gendata.logistic(n0,n,pi.true=pi.true,logistic.fun,dim.X=dim.X)
  pdata <- alldata$pdata;udata <- alldata$udata
  
  test <- EM.addDR.linear.minit(seed,pdata,udata,maxit=500,thres=1e-4,
                                pi.true=pi.true,ini.coef=NULL,each=5)
  out <- list()
  if(!is.null(test$pi.est)){
    out$conv <- test$conv
    out$pi <- test$pi.est
    label.pred <- ifelse(test$label.pred>0.5,1,0)
    label.true <- alldata$y.udata
    out$classperf <- data.frame(FP=length(which(label.pred-label.true==1))/length(which(label.true==0)),
                                FN=length(which(label.pred-label.true==-1))/length(which(label.true==1)),
                                miserr=length(which(label.pred-label.true!=0))/length(label.true))
  }else{out$conv <- F}
  return(out)
}

sum.paracomp <- function(res.l,pi.true=0.4){
  conv.status <- sapply(res.l, function(x){x$conv})
  conv.id <- which(conv.status==T)
  out <- list(conv.prop=length(conv.id)/length(res.l))
  piest <- sapply(res.l[conv.id], function(x){x$pi})
  out$pi <- data.frame(bias=mean(piest)-pi.true,sd=sd(piest),
                       mse=(mean(piest)-pi.true)^2+(sd(piest))^2)
  class.df <- do.call(rbind,lapply(res.l[conv.id], function(x){x$classperf}))
  out$classerr <- apply(class.df,2,mean)
  return(out)
}


###Setting 1
logistic.fun <- function(x){
  -16+6*sum(x)
  #-14+sum((x-0.3)^2)*24
  #-13+6*x[1]+24*(x[2]-0.3)^2+1/(x[3]+0.1)-5*cos(5*x[4])+2*exp(4*x[5]-2)
}

ptm <- proc.time()
for(jj in 1:length(nseq)){
  n <- nseq[jj]
  n0 <- n0seq[jj]
  cluster <- makeCluster(min(97,detectCores()-2))
  registerDoParallel(cluster)
  res.l <- list()
  res.l <- foreach(run = 1:RUN)%dopar%{
    para_Linear(run,n0=n0,n=n,pi.true=pi.true,
                logistic.fun=logistic.fun,dim.X=5)
  }
  stopCluster(cl = cluster)
  print(sum.paracomp(res.l,pi.true=pi.true))
}
proc.time() - ptm



###Setting 2
logistic.fun <- function(x){
  #-16+6*sum(x)
  -14+sum((x-0.3)^2)*24
  #-13+6*x[1]+24*(x[2]-0.3)^2+1/(x[3]+0.1)-5*cos(5*x[4])+2*exp(4*x[5]-2)
}

ptm <- proc.time()
for(jj in 1:length(nseq)){
  n <- nseq[jj]
  n0 <- n0seq[jj]
  cluster <- makeCluster(min(97,detectCores()-2))
  registerDoParallel(cluster)
  res.l <- list()
  res.l <- foreach(run = 1:RUN)%dopar%{
    para_Linear(run,n0=n0,n=n,pi.true=pi.true,
                logistic.fun=logistic.fun,dim.X=5)
  }
  stopCluster(cl = cluster)
  print(sum.paracomp(res.l,pi.true=pi.true))
}
proc.time() - ptm



###Setting 3
logistic.fun <- function(x){
  #-16+6*sum(x)
  #-14+sum((x-0.3)^2)*24
  -13+6*x[1]+24*(x[2]-0.3)^2+1/(x[3]+0.1)-5*cos(5*x[4])+2*exp(4*x[5]-2)
}

ptm <- proc.time()
for(jj in 1:length(nseq)){
  n <- nseq[jj]
  n0 <- n0seq[jj]
  cluster <- makeCluster(min(97,detectCores()-2))
  registerDoParallel(cluster)
  res.l <- list()
  res.l <- foreach(run = 1:RUN)%dopar%{
    para_Linear(run,n0=n0,n=n,pi.true=pi.true,
                logistic.fun=logistic.fun,dim.X=5)
  }
  stopCluster(cl = cluster)
  print(sum.paracomp(res.l,pi.true=pi.true))
}
proc.time() - ptm

