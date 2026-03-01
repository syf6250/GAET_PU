library(mgcv)
library(doParallel)
source("FUN-gendata.R")
source("FUN-para.R")
source("FUN-add.R")

###settings
prop <- 5
n1seq <- c(250,500,750,1250)
n0seq <- prop*n1seq
nseq <- n0seq+n1seq
pi.true <- 0.4
dim.X <- 5
logistic.fun <- function(x){
  -14+sum((x-0.3)^2)*24
}



#main function
aic.addEM <- function(run,n0,n,nu.seq,
                      pi.true=0.4,logistic.fun=logistic.fun,dim.X=5){
  seed <- 77*run+7
  set.seed(seed)
  alldata <- gendata.logistic(n0,n,pi.true=pi.true,logistic.fun,dim.X=dim.X)
  pdata <- alldata$pdata;udata <- alldata$udata
  ###function evaluation
  eval.pts <- seq(from=0,to=1,length.out=101)
  eval.df <- matrix(rep(eval.pts,ncol(pdata)),nrow=101,ncol=ncol(pdata))
  colnames(eval.df) <- names(pdata);eval.df <- as.data.frame(eval.df)
  
  
  aicout <- rep(NA,length(nu.seq))#numeric(length(nu.seq))
  edf <- rep(NA,length(nu.seq))
  emout.l <- list()
  for(j in 1:length(nu.seq)){
    emout <- EM.addDR.init4(seed,pdata,udata,maxit=1000,thres=1e-4,
                            pi.true=pi.true,k=10,ord=3,ord.pen=2,nu=nu.seq[j],
                            each=5,eval.df = eval.df)
    if(isTRUE(emout$conv)){
      aicout[j] <- emout$aic
      edf[j] <- emout$edf
      emout.l[[j]] <- emout
    }
  }
  id <- which.min(aicout)
  if(length(id)>0){
    emout <- emout.l[[id]]
    out <- list(tuneid=id,tune=nu.seq[id],
                conv=emout$conv, initMethd=emout$best.init.method,pi=emout$piEst)
    label.pred <- ifelse(emout$label.pred>0.5,1,0)
    label.true <- alldata$y.udata
    out$classperf <- data.frame(FP=length(which(label.pred-label.true==1))/length(which(label.true==0)),
                                FN=length(which(label.pred-label.true==-1))/length(which(label.true==1)),
                                miserr=length(which(label.pred-label.true!=0))/length(label.true))
    out$eval.terms <- emout$eval.terms
    out$aic.df <- data.frame(aic=aicout,edf=edf)
  }else{
    out <- list(conv=F)
  }
  return(out)
}

sum.aic.addEM <- function(res.l,pi.true=0.4,RUN,nu.seq){
  allout <- matrix(nrow=1,ncol=7)
  colnames(allout) <- c("conv","bias","sd","mse","FP","FN","Err")
  
  conv.status <- do.call(rbind,lapply(res.l, function(x){x$conv}))
  allout[,1] <- apply(conv.status,2,mean)
  
  edfall <- do.call(rbind,lapply(res.l, function(x){x$aic.df$edf}))
  edfmean <- apply(edfall,2,mean,na.rm=T)
  
  piest <- do.call(rbind,lapply(res.l, function(x){x$pi}))
  allout[,2] <- apply(piest,2,mean,na.rm=T)-pi.true
  allout[,3] <- apply(piest,2,sd,na.rm=T)
  allout[,4] <- allout[,2]^2+allout[,3]^2
  
  FP <- do.call(rbind,lapply(res.l, function(x){x$classperf[,1]}))
  FN <- do.call(rbind,lapply(res.l, function(x){x$classperf[,2]}))
  Err <- do.call(rbind,lapply(res.l, function(x){x$classperf[,3]}))
  allout[,5] <- apply(FP,2,mean,na.rm=T)
  allout[,6] <- apply(FN,2,mean,na.rm=T)
  allout[,7] <- apply(Err,2,mean,na.rm=T)
  
  init.table <- do.call(rbind,lapply(res.l, function(x){x$initMethd}))
  
  tune.table <- table(factor(sapply(res.l,function(x){x$tuneid}),levels = 1:length(nu.seq)))
  tune.df <- data.frame(tune=nu.seq,frq=as.numeric(tune.table/RUN))
  
  eval.l <- lapply(res.l, function(x){x$eval.terms})
  
  return(list(summary=allout,pi.all=piest,Err.all=Err,init.table=init.table,
              tune.all=sapply(res.l,function(x){x$tuneid}),
              tune.df=tune.df,eval.l=eval.l,edfmean=edfmean))
}





###simulation
RUN <- 2000
nu.seq.ori <- 10^seq(from=-7.5,to=-5.5,by=0.5)
perf.all <- NULL

ptm <- proc.time()
for(jj in 1:length(nseq)){
  n <- nseq[jj]
  n0 <- n0seq[jj]
  nu.seq <- nu.seq.ori*n0
  
  cluster <- makeCluster(min(97,detectCores()-2))
  registerDoParallel(cluster)
  res.l <- list()
  res.l <- foreach(run = 1:RUN,.packages = c("mgcv"))%dopar%{
    aic.addEM(run,n0=n0,n=n,nu.seq=nu.seq,
              pi.true=pi.true,logistic.fun=logistic.fun,dim.X=dim.X)
  }
  stopCluster(cl = cluster)
  
  sumout <- sum.aic.addEM(res.l,pi.true=pi.true,RUN,nu.seq)
  save(sumout,file=paste0("n",jj,"_aicsum_DRquad.Rdata"))
  perf.all <- rbind(perf.all,sumout$summary)
}
proc.time() - ptm
print(perf.all)




