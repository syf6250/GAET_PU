library(mgcv)
library(gss)
library(fda)
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
zeta <- 0#Set different values for zeta (when zeta==0, the result is the empirical size)
logistic.fun <- function(x){
  -14+sum((x[2:5]-0.3)^2)*24+(x[1]-0.3)^2*24*zeta
}
base.f <- function(x){(x-0.3)^2*24*zeta}

quad <- gauss.quad(size = 20, interval = c(0,1))
nodes <- quad$pt
weights <- quad$wt
bsval <- fourier(nodes,11,1)


###main function
EM.btfun <- function(run,n0,n,nu.seq,B,nodes,weights,bsval,
                     pi.true,logistic.fun,base.f,dim.X,btune=F){
  seed <- 77*run+7
  set.seed(seed)
  alldata <- gendata.logistic(n0,n,pi.true=pi.true,logistic.fun=logistic.fun,dim.X=dim.X)
  pdata <- alldata$pdata;udata <- alldata$udata
  ###function evaluation
  eval.df <- matrix(rep(nodes,ncol(pdata)),nrow=length(nodes),ncol=ncol(pdata))
  colnames(eval.df) <- names(pdata);eval.df <- as.data.frame(eval.df)
  
  ###aic
  aicres <- AIC.addDR(nu.seq,seed,pdata,udata,maxit=1000,thres=1e-4,
                      pi.true=pi.true,k=10,ord=3,ord.pen=2,eval.df=eval.df,
                      Fvalout=T,initinput=NULL)
  if(isFALSE(aicres$conv)){return(list(conv=F,tuneid=NA,initMethd=NA,
                                       bsumsq=NA,bBoot=rep(NA,4),testout05=NA,testout01=NA))}
  
  out <- list(tuneid=aicres$tuneid,tune=nu.seq[aicres$tuneid],
              conv=aicres$conv,initMethd=aicres$initMethd,pi=aicres$pi)
  #shift to obtain u1 est!
  X <- rbind(pdata,udata,udata)
  shift <- apply(base.f(X), 2, mean)
  eval.terms2 <- aicres$eval.terms+matrix(rep(shift,nrow(aicres$eval.terms)),
                                          nrow=nrow(aicres$eval.terms),byrow=T)#modified function estimates
  bsumsq <- crossprod(as.numeric(weights%*%(bsval*as.numeric(eval.terms2[,1]))))
  out$bsumsq <- as.numeric(bsumsq)
  best.tem <- eval.terms2[,1]
  
  ###bootstrap
  if(isFALSE(btune)){
    tune.work <- out$tune
  }else{tune.work <- nu.seq}
  
  b0sumsq.seq <- rep(NA,B)
  for(jj in 1:B){
    set.seed(jj)
    pdata.selec.id <- sample(c(1:nrow(pdata)),nrow(pdata),T)
    pdata.selec <- pdata[pdata.selec.id,]
    udata.selec.id <- sample(c(1:nrow(udata)),nrow(udata),T)
    udata.selec <- udata[udata.selec.id,]
    
    getinit <- list(linpred=aicres$init$linpred[udata.selec.id],conv=aicres$conv,pi=aicres$pi)
    bootout <- AIC.addDR(tune.work,seed,pdata.selec,udata.selec,maxit=1000,thres=1e-4,
                         pi.true=pi.true,k=10,ord=3,ord.pen=2,eval.df=eval.df,
                         Fvalout=T,initinput=getinit)
    if(isTRUE(bootout$conv)){
      X.selec <- rbind(pdata.selec,udata.selec,udata.selec)
      shift <- apply(base.f(X.selec), 2, mean)
      eval.terms2 <- bootout$eval.terms+matrix(rep(shift,nrow(bootout$eval.terms)),
                                               nrow=nrow(bootout$eval.terms),byrow=T)
      diff.eval <- eval.terms2[,1]-best.tem
      b0sumsq.seq[jj] <- crossprod(as.numeric(weights%*%(bsval*as.numeric(diff.eval))))
    }
  }
  out$b0sumsq.seq <- b0sumsq.seq
  Bootres <- c(quantile(b0sumsq.seq,probs = c(0.95,0.99),na.rm=T),
               mean(b0sumsq.seq,na.rm=T),sd(b0sumsq.seq,na.rm=T))
  names(Bootres)[3] <- 'mean';names(Bootres)[4] <- 'sd'
  out$bBoot <- Bootres
  out$testout05 <- ifelse(out$bsumsq>as.numeric(Bootres[1]),T,F)
  out$testout01 <- ifelse(out$bsumsq>as.numeric(Bootres[2]),T,F)
  return(out)
}

sum.EM.btfun <- function(res.l,pi.true=0.4,RUN,nu.seq){
  allout <- matrix(nrow=1,ncol=6)
  colnames(allout) <- c("conv","bias","sd","mse","power05","power01")
  
  conv.status <- sapply(res.l, function(x){x$conv})
  allout[,1] <- mean(conv.status)
  piest <- sapply(res.l, function(x){x$pi})
  allout[,2] <- mean(piest,na.rm=T)-pi.true
  allout[,3] <- sd(piest,na.rm=T)
  allout[,4] <- allout[,2]^2+allout[,3]^2
  allout[,5] <- sum(sapply(res.l, function(x){x$testout05}),na.rm = T)
  allout[,6] <- sum(sapply(res.l, function(x){x$testout01}),na.rm = T)
  
  init.methods <- sapply(res.l, function(x){x$initMethd})
  
  tune.table <- table(factor(sapply(res.l,function(x){x$tuneid}),levels = 1:length(nu.seq)))
  tune.df <- data.frame(tune=nu.seq,frq=as.numeric(tune.table/RUN))
  return(list(summary=allout,init.methods=init.methods,tune.df=tune.df,
              bsumsq=sapply(res.l, function(x){x$bsumsq}),bBoot=lapply(res.l, function(x){x$bBoot})))
}


###simulations
RUN <- 2000
nu.seq.ori <- 10^seq(from=-7.5,to=-5.5,by=0.5)
allres <- list()
perf.all <- NULL;tune.df <- NULL

ptm <- proc.time()
for(kk in 1:length(nseq)){
  n <- nseq[kk]
  n0 <- n0seq[kk]
  nu.seq <- nu.seq.ori*n0
  
  cluster <- makeCluster(min(97,detectCores()-2))
  registerDoParallel(cluster)
  res.l <- list()
  res.l <- foreach(run = 1:RUN,.packages = c("mgcv"))%dopar%{
    EM.btfun(run,n0=n0,n=n,nu.seq=nu.seq,B=500,nodes=nodes,weights=weights,bsval=bsval,
             pi.true=pi.true,logistic.fun=logistic.fun,base.f=base.f,dim.X=dim.X,btune=F)
  }
  stopCluster(cl = cluster)
  sumout <- sum.EM.btfun(res.l,pi.true=pi.true,RUN=RUN,nu.seq=nu.seq)
  allres[[kk]] <- sumout
  perf.all <- rbind(perf.all,sumout$summary)
  tune.df <- rbind(tune.df,sumout$tune.df$frq)
  #print(perf.all)
  #print(tune.df)
  print(sumout$summary)
}
#save(allres,file="sumBTfun-quad_scale_zeta0.Rdata")
proc.time() - ptm

perf.all[,5]/RUN

