library(mgcv)
library(fastDummies)
library(doParallel)
library(ggplot2)
source("FUN-para.R")
source("FUN-real.R")

raw <- read.csv("wilt_training.csv")
raw <- rbind(raw,read.csv("wilt_testing.csv"))
raw <- raw[,c(2:6,1)]
for(j in 1:(ncol(raw)-1)){
  raw[,j] <- (raw[,j]-min(raw[,j]))/(max(raw[,j])-min(raw[,j]))
}

dataY0 <- subset(raw,class=='w')#Y=0
dataY1 <- subset(raw,class=='n')#Y=1


prob.seq <- c(0.7,0.8,0.9) #=P(Z=1|Y=1) 
#pi=P(Y=1|Z=0)=P(Z=0|Y=1)P(Y=1)/{P(Z=0|Y=1)P(Y=1)+P(Y=0)}=(1-prob)P(Y=1)/{(1-prob)P(Y=1)+P(Y=0)}
pi.true.seq <- (1-prob.seq)*nrow(dataY1)/((1-prob.seq)*nrow(dataY1)+nrow(dataY0))


simu.fun.bank <- function(run,prob,pi.true,dataY1,dataY0,nu.seq,dimC,dimL,B=500){
  ###data
  set.seed(7*run)
  z.y1 <- rbinom(n=nrow(dataY1),size = 1,prob = prob)
  pdata <- dataY1[which(z.y1==1),]
  udata <- rbind(dataY0,dataY1[which(z.y1==0),])
  label.true <- ifelse(udata$class=='w',0,1)
  pdata$class <- NULL; udata$class <- NULL
  names.seq <- NULL
  for(j in 1:ncol(pdata)){names.seq <- c(names.seq,paste0("X",j))}
  names(pdata) <- names.seq;names(udata) <- names.seq
  #mean(label.true) #not exactly pi.true
  pdata <- as.data.frame(pdata);udata <- as.data.frame(udata)
  PUratio <- nrow(pdata)/nrow(udata)
  
  out <- list(prob=prob,pi.true=pi.true,PUratio=PUratio)
  ######method-linear
  testL <- EM.addDR.linear.minit(77, pdata, udata,maxit=1000,thres=1e-4,
                                 pi.true=pi.true,ini.coef=NULL,each=5)
  out$pi.L <- testL$pi.est
  label.pred <- ifelse(testL$label.pred>0.5,1,0)
  out$class.L <- data.frame(FP=length(which(label.pred-label.true==1))/length(which(label.true==0)),
                            FN=length(which(label.pred-label.true==-1))/length(which(label.true==1)),
                            miserr=length(which(label.pred-label.true!=0))/length(label.true))
  out$convL <- testL$conv
  
  
  ######method-additive 
  test <- AIC.addDR2(nu.seq,7777,pdata,udata,maxit=1000,thres=1e-4,
                     pi.true=pi.true,k=10,ord=3,ord.pen=2,eval.df=NULL,
                     initinput=NULL,dimC=dimC,dimL=dimL,each=5)
  if(isTRUE(test$conv)){
    out$icall <- test$icall;out$tuneid <- test$tuneid;out$pi <- test$pi;out$initmethod <- test$initMethd
    initA <- test$init
    label.pred <- ifelse(test$label.pred>0.5,1,0)
    out$class <- data.frame(FP=length(which(label.pred-label.true==1))/length(which(label.true==0)),
                            FN=length(which(label.pred-label.true==-1))/length(which(label.true==1)),
                            miserr=length(which(label.pred-label.true!=0))/length(label.true))
  }else{
    initA <- NULL
    out$tuneid <- NA;out$pi <- NA;out$initmethod <- NA
    out$class <- data.frame(FP=NA,FN=NA,miserr=NA)
  }
  
  
  ######bootstrap interval (to be done; must use initA!)
  piall <- rep(NA,B)
  if(!is.null(initA)){
    for(run in 1:B){
      #re-sample:
      set.seed(run)
      pdata.selec.id <- sample(c(1:nrow(pdata)),nrow(pdata),T)
      pdata.selec <- pdata[pdata.selec.id,]
      udata.selec.id <- sample(c(1:nrow(udata)),nrow(udata),T)
      udata.selec <- udata[udata.selec.id,]
      #modify initial
      getinit <- initA#additive fit as initials
      getinit$linpred <- getinit$linpred[udata.selec.id]
      getinit$linpred.labeled <- getinit$linpred.labeled[pdata.selec.id]
      res <- AIC.addDR2(nu.seq=nu.seq[test$tuneid],#no tuning selection
                        seed=NULL,pdata.selec,udata.selec,
                        maxit=1000,thres=1e-4,
                        pi.true=NULL,k=10,ord=3,ord.pen=2,eval.df=test$eval.pts,
                        initinput=getinit,dimC=dimC,dimL=dimL,each=NULL)
      if(isTRUE(res$conv)){piall[run] <- res$pi}
    }
    out$conv.Boot <- 1-sum(is.na(piall))/B
    
    quant <- quantile(piall,probs = c(0.025,0.975),na.rm = T)
    out$pi.boot <- data.frame(q025=quant[1],q975=quant[2],mean=mean(piall,na.rm=T),sd=sd(piall,na.rm=T))
    interval.low <- 2*(test$pi)-quant[2];interval.up <- 2*(test$pi)-quant[1]
    out$cover <- ifelse(interval.low<=pi.true && interval.up>=pi.true,T,F)
    
    piall <- log(piall/(1-piall))
    quant <- quantile(piall,probs = c(0.025,0.975),na.rm = T)
    out$pi.boot.logit <- data.frame(q025=quant[1],q975=quant[2],mean=mean(piall,na.rm=T),sd=sd(piall,na.rm=T))
    pi2 <- log(test$pi/(1-test$pi));pi.true2 <- log(pi.true/(1-pi.true))
    interval.low <- 2*pi2-quant[2];interval.up <- 2*pi2-quant[1]
    out$cover.logit <- ifelse(interval.low<=pi.true2 && interval.up>=pi.true2,T,F)
  }else{
    out$conv.Boot <- NA
    out$pi.boot <- data.frame(q025=NA,q975=NA,mean=NA,sd=NA)
    out$cover <- NA
  }
  return(out)
}





nu.seq <- 10^seq(from=-10.5,to=-6.5,by=1)
RUN <- 500

ptm <- proc.time()
allres <- list()
for(kk in 1:length(prob.seq)){
  
  cluster <- makeCluster(min(95,detectCores()-2))
  registerDoParallel(cluster)
  out.l <- list()
  out.l <- foreach(run = 1:RUN,.packages = c("mgcv"))%dopar%{
    simu.fun.bank(run,prob=prob.seq[kk],pi.true=pi.true.seq[kk],
                  dataY1=dataY1,dataY0=dataY0,nu.seq=nu.seq,dimC=5,dimL=0,B=500) 
  }
  stopCluster(cl = cluster)
  allres[[kk]] <- out.l
}
proc.time() - ptm# 

save(allres,file="wilt_P789.Rdata")



###summary
perf.l <- list()
for(kk in 1:length(prob.seq)){
  piestL <- sapply(allres[[kk]], function(x){x$pi.L})
  piest <- sapply(allres[[kk]], function(x){x$pi})
  classL <- do.call(rbind,lapply(allres[[kk]], function(x){x$class.L}))
  names(classL) <- c("FP.L","FN.L","miserr.L")
  classA <- do.call(rbind,lapply(allres[[kk]], function(x){x$class}))
  convL <- sapply(allres[[kk]], function(x){x$convL})
  cover.res <- sapply(allres[[kk]], function(x){x$cover})
  cover.res2 <- sapply(allres[[kk]], function(x){x$cover.logit})
  convBoot <- sapply(allres[[kk]], function(x){x$conv.Boot})
  
  perf.l[[kk]] <- list(tune=table(sapply(allres[[kk]], function(x){x$tuneid})),#table of tuning selection
                       est=data.frame(noconvL=(RUN-sum(convL))/RUN,noconv=sum(is.na(piest))/RUN,
                                      pi.true=pi.true.seq[kk],meanL=mean(piestL,na.rm=T),sdL=sd(piestL,na.rm=T),
                                      mean=mean(piest,na.rm=T),sd=sd(piest,na.rm=T),
                                      cover=sum(cover.res,na.rm = T)/RUN,
                                      convBoot=mean(convBoot,na.rm=T)),
                       classification=c(apply(classL,2,mean,na.rm=T),apply(classA,2,mean,na.rm=T)) )
}
aa <- do.call(rbind,lapply(perf.l, function(x){x$est}))
#lapply(perf.l, function(x){x$tune})
(aa$meanL-aa$pi.true)*100;(aa$mean-aa$pi.true)*100#bias
aa$sdL*100;aa$sd*100#se
data.frame(mseL=(aa$meanL-aa$pi.true)^2+(aa$sdL)^2,mse=(aa$mean-aa$pi.true)^2+(aa$sd)^2)*100#mse
do.call(rbind,lapply(perf.l, function(x){x$classification}))#classfication
sapply(perf.l, function(x){x$est$cover})#empirical coverage probabilities


