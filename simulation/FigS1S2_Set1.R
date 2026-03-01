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
  -16+6*sum(x)
}
nu.seq.ori <- 10^seq(from=-6,to=0,by=1)


#main function
fix.addEM <- function(run,n0,n,nu.seq,
                      pi.true=0.4,logistic.fun=logistic.fun,dim.X=5){
  seed <- 77*run+7
  set.seed(seed)
  alldata <- gendata.logistic(n0,n,pi.true=pi.true,logistic.fun,dim.X=dim.X)
  pdata <- alldata$pdata;udata <- alldata$udata
  ###function evaluation
  eval.pts <- seq(from=0,to=1,length.out=101)
  eval.df <- matrix(rep(eval.pts,ncol(pdata)),nrow=101,ncol=ncol(pdata))
  colnames(eval.df) <- names(pdata);eval.df <- as.data.frame(eval.df)
  
  
  aicout <- edf <- rep(NA,length(nu.seq))
  piest <- rep(NA,length(nu.seq))
  FP <- FN <- miserr <- rep(NA,length(nu.seq))
  conv <- rep(T,length(nu.seq))
  #emout.l <- list()
  for(j in 1:length(nu.seq)){
    emout <- EM.addDR.init4(seed,pdata,udata,maxit=1000,thres=1e-4,
                            pi.true=pi.true,k=10,ord=3,ord.pen=2,nu=nu.seq[j],
                            each=5,eval.df = eval.df)
    if(isTRUE(emout$conv)){
      aicout[j] <- emout$aic
      edf[j] <- emout$edf
      #emout.l[[j]] <- emout
      piest[j] <- emout$piEst
      
      label.pred <- ifelse(emout$label.pred>0.5,1,0)
      label.true <- alldata$y.udata
      FP[j] <- length(which(label.pred-label.true==1))/length(which(label.true==0))
      FN[j] <- length(which(label.pred-label.true==-1))/length(which(label.true==1))
      miserr[j] <- length(which(label.pred-label.true!=0))/length(label.true)
    }else{
      conv[j] <- F
    }
  }
  
  return(list(conv=conv,piest=piest,aic=aicout,edf=edf,FP=FP,FN=FN,Err=miserr))
}

sum.fix.addEM <- function(res.l,pi.true=0.4,nu.seq){
  allout <- matrix(nrow=length(nu.seq),ncol=9)
  colnames(allout) <- c("conv","bias","sd","mse","FP","FN","Err","edf","aic")
  
  conv.status <- do.call(rbind,lapply(res.l, function(x){x$conv}))
  allout[,1] <- apply(conv.status,2,mean)
  edfall <- do.call(rbind,lapply(res.l, function(x){x$edf}))
  allout[,8] <- apply(edfall,2,mean,na.rm=T)
  aicall <- do.call(rbind,lapply(res.l, function(x){x$aic}))
  allout[,9] <- apply(aicall,2,mean,na.rm=T)
  
  piall <- do.call(rbind,lapply(res.l, function(x){x$piest}))
  allout[,2] <- apply(piall,2,mean,na.rm=T)-pi.true
  allout[,3] <- apply(piall,2,sd,na.rm=T)
  allout[,4] <- allout[,2]^2+allout[,3]^2
  
  FP <- do.call(rbind,lapply(res.l, function(x){x$FP}))
  FN <- do.call(rbind,lapply(res.l, function(x){x$FN}))
  Err <- do.call(rbind,lapply(res.l, function(x){x$Err}))
  allout[,5] <- apply(FP,2,mean,na.rm=T)
  allout[,6] <- apply(FN,2,mean,na.rm=T)
  allout[,7] <- apply(Err,2,mean,na.rm=T)
  
  return(allout)
}





###simulation
RUN <- 2000
allout <- list()
mse.t <- err.t <- NULL

ptm <- proc.time()
for(jj in 1:length(nseq)){
  n <- nseq[jj]
  n0 <- n0seq[jj]
  nu.seq <- nu.seq.ori*n0
  
  cluster <- makeCluster(min(120,detectCores()-2))
  registerDoParallel(cluster)
  res.l <- list()
  res.l <- foreach(run = 1:RUN,.packages = c("mgcv"))%dopar%{
    fix.addEM(run,n0=n0,n=n,nu.seq=nu.seq,
              pi.true=pi.true,logistic.fun=logistic.fun,dim.X=dim.X)
  }
  stopCluster(cl = cluster)
  
  sumout <- sum.fix.addEM(res.l,pi.true=pi.true,nu.seq)
  allout[[jj]] <- sumout
  #openxlsx2::write_xlsx(sumout,paste0("fixperf_prop5_DRlin2_n",n,".xlsx"))
  mse.t <- rbind(mse.t,sumout[,4]);err.t <- rbind(err.t,sumout[,7])
  #print(sumout[,4])
}
proc.time() - ptm
print(mse.t);print(err.t)

#save(allout,file="sensi_DRlin_Run2000.Rdata")
###plot:
mse.ratio <- t(apply(mse.t,1,function(x){x/min(x)}))
err.ratio <- t(apply(err.t,1,function(x){x/min(x)}))
save(mse.ratio,file="sensi_DRlin_Run2000_mseR.Rdata")
save(err.ratio,file="sensi_DRlin_Run2000_errR.Rdata")




