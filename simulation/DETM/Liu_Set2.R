library(glmnet)
library(doParallel)
source("FUN-gendata.R")
source("PUEM/R/DETMfull.R")
source("PUEM/R/DETMnull.R")
source("puem_detm_wrapper_simple_conv.R")

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

unify.fun <- function(run,n0,n,pi.true=0.4,logistic.fun=logistic.fun,dim.X=5){
  seed <- 77*run+7
  set.seed(seed)
  alldata <- gendata.logistic(n0,n,pi.true=pi.true,logistic.fun,dim.X=dim.X)
  pdata <- alldata$pdata;udata <- alldata$udata
  
  fit <- puem_detm_analyze(posdata = pdata,negdata = udata,pi_init = pi.true,
                           pi_input = pi.true,pihalf = TRUE,standardize = TRUE)
  
  if(isTRUE(fit$convergence)){
    label.pred <- fit$predictions_unlabeled$predicted_label
    label.true <- alldata$y.udata
    df <- data.frame(FP=length(which(label.pred-label.true==1))/length(which(label.true==0)),
                     FN=length(which(label.pred-label.true==-1))/length(which(label.true==1)),
                     miserr=length(which(label.pred-label.true!=0))/length(label.true))
    return(list(conv=T,pi=fit$pi_estimation,cover=fit$pi_input_check$covered_by_CI,classperf=df,
                CI=fit$pi_confidence_interval$ci))
  }else{
    return(list(conv=F,pi=NA,cover=NA,classperf=data.frame(FP=NA,FN=NA,miserr=NA),CI=rep(NA,2)))
  }
}
sum.unify.fun <- function(res.l,pi.true=0.4){
  allout <- matrix(nrow=1,ncol=8)
  colnames(allout) <- c("conv","bias","sd","mse","FP","FN","Err","cover.prop")
  
  conv.status <- do.call(rbind,lapply(res.l, function(x){x$conv}))
  allout[,1] <- apply(conv.status,2,mean)
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
  allout[,8] <- apply(do.call(rbind,lapply(res.l, function(x){x$cover})),2,mean,na.rm=T)
  
  piall <- cbind(do.call(rbind,lapply(res.l, function(x){x$CI})),piest)
  colnames(piall) <- c("lower","upper","point")
  return(list(summary=allout,piall=piall))
}


RUN <- 2000
perf.all <- NULL

ptm <- proc.time()
for(jj in 1:length(nseq)){
  n <- nseq[jj]
  n0 <- n0seq[jj]
  ncores <- 120
  cluster <- makeCluster(ncores)
  registerDoParallel(cluster)
  res.l <- list()
  res.l <- foreach(run = 1:RUN,.packages = c("glmnet"))%dopar%{
    unify.fun(run,n0=n0,n=n,
              pi.true=pi.true,logistic.fun=logistic.fun,dim.X=dim.X)
  }
  stopCluster(cl = cluster)
  sumout <- sum.unify.fun(res.l,pi.true=pi.true)
  save(sumout, file = paste0("n", jj, "_Liu_DRquad_RUN", RUN, ".Rdata"))
  perf.all <- rbind(perf.all,sumout$summary)
}
elapsed <- proc.time() - ptm

print(elapsed)
print(perf.all)





