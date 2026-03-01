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


#main function
aic.addEM.bt <- function(run,n0,n,nu.seq,B=500,
                         pi.true=0.4,logistic.fun=logistic.fun,dim.X=5,btune=T){
  seed <- 77*run+7
  set.seed(seed)
  alldata <- gendata.logistic(n0,n,pi.true=pi.true,logistic.fun,dim.X=dim.X)
  pdata <- alldata$pdata;udata <- alldata$udata
  ###function evaluation
  eval.pts <- seq(from=0,to=1,length.out=101)
  eval.df <- matrix(rep(eval.pts,ncol(pdata)),nrow=101,ncol=ncol(pdata))
  colnames(eval.df) <- names(pdata);eval.df <- as.data.frame(eval.df)
  
  ###aic
  aicres <- AIC.addDR(nu.seq,seed,pdata,udata,maxit=1000,thres=1e-4,
                      pi.true=pi.true,k=10,ord=3,ord.pen=2,eval.df=eval.df,
                      Fvalout=F,initinput=NULL)
  if(isFALSE(aicres$conv)){return(list(conv=F))}
  #collect results:
  out <- list(tuneid=aicres$tuneid,tune=nu.seq[aicres$tuneid],
              conv=aicres$conv, initMethd=aicres$initMethd,pi=aicres$pi)
  label.pred <- aicres$label.pred
  label.true <- alldata$y.udata
  out$classperf <- data.frame(FP=length(which(label.pred-label.true==1))/length(which(label.true==0)),
                              FN=length(which(label.pred-label.true==-1))/length(which(label.true==1)),
                              miserr=length(which(label.pred-label.true!=0))/length(label.true))
  
  
  ###bootstrap
  piall.seq <- rep(NA,B)
  if(isFALSE(btune)){
    tune.work <- out$tune
  }else{tune.work <- nu.seq}
  for(jj in 1:B){
    set.seed(jj)
    pdata.selec.id <- sample(c(1:nrow(pdata)),nrow(pdata),T)
    pdata.selec <- pdata[pdata.selec.id,]
    udata.selec.id <- sample(c(1:nrow(udata)),nrow(udata),T)
    udata.selec <- udata[udata.selec.id,]
    
    getinit <- list(linpred=aicres$init$linpred[udata.selec.id],
                    conv=aicres$conv,pi=aicres$pi)
    bootout <- AIC.addDR(tune.work,seed,pdata.selec,udata.selec,maxit=1000,thres=1e-4,
                         pi.true=pi.true,k=10,ord=3,ord.pen=2,eval.df=eval.df,
                         Fvalout=F,initinput=getinit)
    if(isTRUE(bootout$conv)){
      piall.seq[jj] <- bootout$pi
    }
  }
  
  ###summarize
  info.bt.pi <- matrix(c(quantile(piall.seq,probs = c(0.025,0.975,0.05,0.95),na.rm=T),
                         mean(piall.seq,na.rm = T),sd(piall.seq,na.rm = T)),nrow=1)
  colnames(info.bt.pi) <- c("L1","U1","L2","U2","mean","se")
  out$piBoot <- as.data.frame(info.bt.pi)
  
  piall.seq <- log(piall.seq/(1-piall.seq))
  info.bt.pi <- matrix(c(quantile(piall.seq,probs = c(0.025,0.975,0.05,0.95),na.rm=T),
                         mean(piall.seq,na.rm = T),sd(piall.seq,na.rm = T)),nrow=1)
  colnames(info.bt.pi) <- c("L1","U1","L2","U2","mean","se")
  out$piBoot.logit <- as.data.frame(info.bt.pi)
  return(out)
}

#summarize the results
sum.bt.addEM <- function(res.l,pi.true=0.4,RUN,nu.seq){
  allout <- matrix(nrow=1,ncol=7)
  colnames(allout) <- c("conv","bias","sd","mse","FP","FN","Err")
  
  conv.status <- sapply(res.l, function(x){x$conv})
  allout[,1] <- mean(conv.status)
  
  piest <- sapply(res.l, function(x){x$pi})
  allout[,2] <- mean(piest,na.rm=T)-pi.true
  allout[,3] <- sd(piest,na.rm=T)
  allout[,4] <- allout[,2]^2+allout[,3]^2
  
  FP <- sapply(res.l, function(x){x$classperf[,1]});allout[,5] <- mean(FP,na.rm=T)
  FN <- sapply(res.l, function(x){x$classperf[,2]});allout[,6] <- mean(FN,na.rm=T)
  Err <- sapply(res.l, function(x){x$classperf[,3]});allout[,7] <- mean(Err,na.rm=T)
  
  init.methods <- sapply(res.l, function(x){x$initMethd})
  
  tune.table <- table(factor(sapply(res.l,function(x){x$tuneid}),levels = 1:length(nu.seq)))
  tune.df <- data.frame(tune=nu.seq,frq=as.numeric(tune.table/RUN))
  
  bootpi <- do.call(rbind,lapply(res.l, function(x){x$piBoot}))
  bootpi2 <- do.call(rbind,lapply(res.l, function(x){x$piBoot.logit}))
  
  return(list(summary=allout,pi.all=piest,pi.all.logit=log(piest/(1-piest)),
              Err.all=Err,init.methods=init.methods,
              tune.df=tune.df,bootpi=bootpi,bootpi.logit=bootpi2))
}
#compute coverage rate
cover.prop <- function(pi.all,bootpi,pi.true=0.4,debias=F){
  #piall is sumout$pi.all or sumout$pi.all.logit
  #bootpi is sumout$bootpi or sumout$bootpi.logit
  cover.mat <- matrix(nrow=length(pi.all),ncol=6)
  colnames(cover.mat) <- c("L95M1","L90M1","L95M2","L90M2","L95M3","L90M3")
  
  if(isFALSE(debias)){
    tominus <- 0
  }else{
    tominus <- bootpi$mean-pi.all
  }
  
  ll <- bootpi$L1-tominus; uu <- bootpi$U1-tominus
  cover.mat[,1] <- ifelse(ll<=pi.true & uu>=pi.true,T,F)
  ll <- bootpi$L2-tominus; uu <- bootpi$U2-tominus
  cover.mat[,2] <- ifelse(ll<=pi.true & uu>=pi.true,T,F)
  
  ll <- 2*(pi.all)-bootpi$U1-tominus
  uu <- 2*(pi.all)-bootpi$L1-tominus
  cover.mat[,3] <- ifelse(ll<=pi.true & uu>=pi.true,T,F)
  ll <- 2*(pi.all)-bootpi$U2-tominus
  uu <- 2*(pi.all)-bootpi$L2-tominus
  cover.mat[,4] <- ifelse(ll<=pi.true & uu>=pi.true,T,F)
  
  ll <- pi.all-tominus-(bootpi$se)*1.96
  uu <- pi.all-tominus+(bootpi$se)*1.96
  cover.mat[,5] <- ifelse(ll<=pi.true & uu>=pi.true,T,F)
  ll <- pi.all-tominus-(bootpi$se)*1.64
  uu <- pi.all-tominus+(bootpi$se)*1.64
  cover.mat[,6] <- ifelse(ll<=pi.true & uu>=pi.true,T,F)
  
  return(apply(cover.mat, 2, mean))
}



#simulation
RUN <- 2000
nu.seq.ori <- 10^seq(from=-5,to=-1,by=1)
tosave <- list()
perf.all <- NULL
cover <- cover.logit <- NULL

ptm <- proc.time()
for(kk in 1:length(nseq)){
  n <- nseq[kk]
  n0 <- n0seq[kk]
  nu.seq <- nu.seq.ori*n0
  
  cluster <- makeCluster(min(95,detectCores()-2))
  registerDoParallel(cluster)
  res.l <- list()
  res.l <- foreach(run = 1:RUN,.packages = c("mgcv"))%dopar%{
    aic.addEM.bt(run,n0=n0,n=n,nu.seq=nu.seq,B=500,
                 pi.true=pi.true,logistic.fun=logistic.fun,dim.X=dim.X,btune = F)
  }
  stopCluster(cl = cluster)
  
  sumout <- sum.bt.addEM(res.l,pi.true=pi.true,RUN,nu.seq)
  tosave[[kk]] <- sumout
  perf.all <- rbind(perf.all,sumout$summary)
  cover <- rbind(cover, cover.prop(sumout$pi.all,sumout$bootpi,pi.true=pi.true))
  cover.logit <- rbind(cover.logit, cover.prop(sumout$pi.all.logit,sumout$bootpi.logit,pi.true=log(pi.true/(1-pi.true))))
}
proc.time() - ptm
print(cover.logit[,3])













