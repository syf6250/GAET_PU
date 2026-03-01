library(mgcv)
library(fastDummies)
library(doParallel)
library(gss)
library(fda)
source("FUN-para.R")
source("FUN-real.R")

#set various kk for different pi (Table S6): kk=1,2,3
kk <- 1

raw <- read.table("spambase.data",sep=',')
raw <- raw[,c(6,26,46,56,5,7,8,16,17,52,53,57,25,27,37,45,58)]
#names(raw) <- c('over','hpl','edu','capmax','V58')
for(j in 1:(ncol(raw)-1)){
  raw[,j] <- log(0.1+raw[,j])
  raw[,j] <- (raw[,j]-min(raw[,j]))/(max(raw[,j])-min(raw[,j]))
}
dataY0 <- subset(raw,V58==1)#Y=0
dataY1 <- subset(raw,V58==0)#Y=1


prob.seq <- c(0.7,0.8,0.9) #=P(Z=1|Y=1) 
#pi=P(Y=1|Z=0)=P(Z=0|Y=1)P(Y=1)/{P(Z=0|Y=1)P(Y=1)+P(Y=0)}=(1-prob)P(Y=1)/{(1-prob)P(Y=1)+P(Y=0)}
pi.true.seq <- (1-prob.seq)*nrow(dataY1)/((1-prob.seq)*nrow(dataY1)+nrow(dataY0))
prob <- prob.seq[kk];pi.true <- pi.true.seq[kk]


###design integral range and evaluation points(outside the loop)
nodes.mat <- weights.mat <-  matrix(nrow=20,ncol=4)
bsval.l <- list()
for(j in 1:4){
  qt <- c(0,1)#quantile(raw[,j],prob=c(0.05,0.95))
  quad <- gauss.quad(size = 20, interval = c(qt[1],qt[2]))
  nodes.mat[,j] <- quad$pt
  weights.mat[,j] <- quad$wt
  bsval.l[[j]] <-  fourier(quad$pt,11,diff(range(qt)))
}

###summarize tunings
load("spam_P789.Rdata")
allres <- allres[[kk]]
nu.seq <- 10^seq(from=-4.5,to=-2.5,by=0.5)
tuneid.seq <- sapply(allres, function(x){x$tuneid})

#density ratio model under true labels (store 'true' functions)
naive.out <- naive.est(dataY1[,c(1:16)],dataY0[,c(1:16)],nu=0.0001,k=6,ord=3,ord.pen=2,
                       dimC=4,dimL=12)


simu.spam.funbt <- function(run,prob,pi.true,dataY1,dataY0,nu.seq,dimC=4,dimL=12,B=500,
                            nodes.mat,weights.mat,bsval.l,tuneid.seq,naive.out){
  ###generalize data
  set.seed(7*run)
  z.y1 <- rbinom(n=nrow(dataY1),size = 1,prob = prob)
  pdata <- dataY1[which(z.y1==1),]
  udata <- rbind(dataY0,dataY1[which(z.y1==0),])
  label.true <- ifelse(udata$V58==1,0,1)
  pdata$V58 <- NULL; udata$V58 <- NULL
  names.seq <- NULL
  for(j in 1:ncol(pdata)){names.seq <- c(names.seq,paste0("X",j))}
  names(pdata) <- names.seq;names(udata) <- names.seq
  pdata <- as.data.frame(pdata);udata <- as.data.frame(udata)
  PUratio <- nrow(pdata)/nrow(udata)
  out <- list(prob=prob,pi.true=pi.true,PUratio=PUratio)
  
  ###design integral range and evaluation points(inside the loop) 
  colnames(nodes.mat) <- colnames(weights.mat) <- names(pdata)[1:dimC]
  fakedata <- pdata[1,][,-c(1:dimC),drop=F]
  fakedata <- do.call(rbind,replicate(nrow(nodes.mat), fakedata, simplify = FALSE))
  eval.pts <- cbind(nodes.mat,fakedata)
  eval.df <- as.data.frame(eval.pts)
  
  
  ###method
  test <- AIC.addDR2(nu.seq[tuneid.seq[run]],7777,pdata,udata,maxit=1000,thres=1e-4,
                     pi.true=pi.true,k=10,ord=3,ord.pen=2,eval.df=eval.df,
                     initinput=NULL,dimC=dimC,dimL=dimL,each=5)
  if(isTRUE(test$conv)){
    out$pi <- test$pi
    X <- rbind(pdata,udata,udata)
    shift <- apply(predict(naive.out$gamObj,type="terms",newdata = X), 2, mean)[-c(1:dimL)]
    eval.terms2 <- test$eval.terms[,-c(1:dimL)]+matrix(rep(shift,nrow(test$eval.terms)),
                                                       nrow=nrow(test$eval.terms),byrow=T)
    out$eval.terms <- eval.terms2#save this output for future plot
    
    bsumsq <- numeric(dimC)
    for(k in 1:dimC){
      bsumsq[k] <- crossprod(as.numeric(weights.mat[,k]%*%(bsval.l[[k]]*as.numeric(eval.terms2[,k]))))
    }
    out$bsumsq <- bsumsq
    
    ###bootstrap
    b0sumsq.all <- matrix(nrow=B,ncol=dimC)
    for(runB in 1:B){
      #re-sample:
      set.seed(runB)
      pdata.selec.id <- sample(c(1:nrow(pdata)),nrow(pdata),T)
      pdata.selec <- pdata[pdata.selec.id,]
      udata.selec.id <- sample(c(1:nrow(udata)),nrow(udata),T)
      udata.selec <- udata[udata.selec.id,]
      #modify initial
      getinit <- test$init#additive fit as initials
      getinit$linpred <- getinit$linpred[udata.selec.id]
      getinit$linpred.labeled <- getinit$linpred.labeled[pdata.selec.id]
      #refit
      res <- AIC.addDR2(nu.seq=test$tune,#no tuning selection
                        seed=NULL,pdata.selec,udata.selec,maxit=1000,thres=1e-4,
                        pi.true=NULL,k=10,ord=3,ord.pen=2,eval.df=eval.df,
                        initinput=getinit,dimC=dimC,dimL=dimL,each=NULL)
      if(isTRUE(res$conv)){
        X.selec <- rbind(pdata.selec,udata.selec,udata.selec)
        shift <- apply(predict(naive.out$gamObj,type="terms",newdata = X.selec),2,mean)[-c(1:dimL)]
        eval.terms2 <- res$eval.terms[,-c(1:dimL)]+matrix(rep(shift,nrow(res$eval.terms)),
                                                          nrow=nrow(res$eval.terms),byrow=T)
        diff.eval <- eval.terms2-out$eval.terms
        for(k in 1:dimC){
          b0sumsq.all[runB,k] <- crossprod(as.numeric(weights.mat[,k]%*%(bsval.l[[k]]*as.numeric(diff.eval[,k]))))
        }
      }
    }
    out$b0sumsq.all <- b0sumsq.all
    quant95 <- apply(b0sumsq.all,2,quantile,na.rm=T,probs=0.95)
    out$sig <- bsumsq>=quant95
    out$quant95 <- quant95
  }
  return(out)
}


RUN <- 500
ptm <- proc.time()
cluster <- makeCluster(min(95,detectCores()-2))
registerDoParallel(cluster)
output <- list()
output <- foreach(run = 1:RUN,.packages = c("mgcv"))%dopar%{
  simu.spam.funbt(run,prob=prob,pi.true=pi.true,dataY1=dataY1,dataY0=dataY0,
                  nu.seq=nu.seq,dimC=4,dimL=12,B=500,
                  nodes.mat=nodes.mat,weights.mat=weights.mat,bsval.l=bsval.l,
                  tuneid.seq=tuneid.seq,naive.out=naive.out)
}
stopCluster(cl = cluster)
proc.time() - ptm#

apply(sapply(output,function(x){x$sig}),1,mean)








