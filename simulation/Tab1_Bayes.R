###simulation for Bayes method in Table 1

library(doParallel)
source("FUN-gendata.R")
RUN <- 2000

paracomp0_Linear_Bayes <- function(run,n0,n,pi.true=0.4,logistic.fun,dim.X=5,PrY1=0.512){
  set.seed(77*run+7)
  alldata <- gendata.logistic(n0,n,pi.true=pi.true,logistic.fun,dim.X=dim.X)
  pdata <- alldata$pdata;udata <- alldata$udata
  
  logi <- apply(udata, 1, logistic.fun)
  label.out <- pi.true*(1-PrY1)/(pi.true*(1-PrY1)+(1-pi.true)*PrY1*exp(-logi))
  label.pred <- ifelse(label.out>0.5,1,0)
  label.true <- alldata$y.udata
  classperf <- data.frame(FP=length(which(label.pred-label.true==1))/length(which(label.true==0)),
                          FN=length(which(label.pred-label.true==-1))/length(which(label.true==1)),
                          miserr=length(which(label.pred-label.true!=0))/length(label.true))
  return(classperf)
}

#Setting 1
logistic.fun <- function(x){
  -16+6*sum(x)
  #-14+sum((x-0.3)^2)*24
  #-13+6*x[1]+24*(x[2]-0.3)^2+1/(x[3]+0.1)-5*cos(5*x[4])+2*exp(4*x[5]-2)
}
# #monte carlo method to see Pr(Y=1) (population level)
# sum <- 0
# for(j in 1:5000000){
#   sum <- sum + 1/(1+exp(-logistic.fun(runif(5))))
# }
# sum/5000000#about  0.409 (linear)/0.509(quad)/0.47(general)


ptm <- proc.time()
cluster <- makeCluster(min(97,detectCores()-2))
registerDoParallel(cluster)
res.l <- list()
res.l <- foreach(run = 1:RUN)%dopar%{
                   paracomp0_Linear_Bayes(run,n0=6250,n=7500,pi.true=0.4,
                                          logistic.fun=logistic.fun,dim.X=5,PrY1=0.409)
                 }
stopCluster(cl = cluster)
proc.time() - ptm
print(apply(do.call(rbind,res.l),2,mean))




#Setting 2
logistic.fun <- function(x){
  #-16+6*sum(x)
  -14+sum((x-0.3)^2)*24
  #-13+6*x[1]+24*(x[2]-0.3)^2+1/(x[3]+0.1)-5*cos(5*x[4])+2*exp(4*x[5]-2)
}
ptm <- proc.time()
cluster <- makeCluster(min(97,detectCores()-2))
registerDoParallel(cluster)
res.l <- list()
res.l <- foreach(run = 1:RUN)%dopar%{
  paracomp0_Linear_Bayes(run,n0=6250,n=7500,pi.true=0.4,
                         logistic.fun=logistic.fun,dim.X=5,PrY1=0.509)
}
stopCluster(cl = cluster)
proc.time() - ptm
print(apply(do.call(rbind,res.l),2,mean))



#Setting 3
logistic.fun <- function(x){
  #-16+6*sum(x)
  #-14+sum((x-0.3)^2)*24
  -13+6*x[1]+24*(x[2]-0.3)^2+1/(x[3]+0.1)-5*cos(5*x[4])+2*exp(4*x[5]-2)
}
ptm <- proc.time()
cluster <- makeCluster(min(97,detectCores()-2))
registerDoParallel(cluster)
res.l <- list()
res.l <- foreach(run = 1:RUN)%dopar%{
  paracomp0_Linear_Bayes(run,n0=6250,n=7500,pi.true=0.4,
                         logistic.fun=logistic.fun,dim.X=5,PrY1=0.47)
}
stopCluster(cl = cluster)
proc.time() - ptm
print(apply(do.call(rbind,res.l),2,mean))



