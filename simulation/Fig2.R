source("FUN-gendata.R")
library(reshape2)
library(ggplot2)
library(latex2exp)
#need to use outputs of 'Tab1_Set1.R', 'Tab1_Set2.R', 'Tab1_Set3.R'
#suppose now they are in the folder "results"
path1 <- "results/"

prop <- 5
n1seq <- c(250,500,750)
n0seq <- prop*n1seq
nseq <- n0seq+n1seq

dim.X <- 5
RUN <- 2000
eval.pts <- seq(from=0,to=1,length.out=101)

#function id to be drawn for all settings
fid <- c(3,3,3)


###setting 2: 
base.f <- function(x){(x-0.3)^2*24}#corresponds to fid[2]!
logistic.fun <- function(x){
  -14+sum((x-0.3)^2)*24
}

#produce points to be plotted
mean.l <- list()
for(k in 1:3){
  load(paste0(path1,"n",k,"_aicsum_DRquad.Rdata"))
  allout <- sumout
  n0 <- n0seq[k]
  n <- nseq[k]
  
  shift.mat <- NULL
  for(run in 1:RUN){
    seed <- 77*run+7
    set.seed(seed)
    alldata <- gendata.logistic(n0=n0,n=n,pi.true=0.4,logistic.fun,dim.X=dim.X)
    pdata <- alldata$pdata;udata <- alldata$udata
    X <- rbind(pdata,udata,udata)
    shift.mat <- rbind(shift.mat,apply(base.f(X), 2, mean))#SHIFTs of true function
  }
  
  feval.l <- list()
  for(run in 1:RUN){
    feval0 <- allout$eval.l[[run]]
    feval.l[[run]] <- feval0+matrix(rep(shift.mat[run,],101),nrow=101,byrow=T)
  }
  
  feval.l2 <- sapply(feval.l, function(x){x[,fid[2]]})
  mean.l[[k]] <- apply(feval.l2, 1, mean)
}


df <- data.frame(t=eval.pts,mean1=mean.l[[1]],mean2=mean.l[[2]],mean3=mean.l[[3]])
df$true <- base.f(df$t)
colnames(df)[2:4] <- c('n=1500','n=3000','n=4500')
df <- df[,c(1,5,2,3,4)]
df2 <- reshape2::melt(df,id.vars = 1,variable.name = 'case')

pquad <- ggplot(data=df2,aes(x=t,y=value,group=case))+
  geom_line(aes(linetype=case, color=case),linewidth=1)+
  ggtitle("Setting 2")+ylab(TeX('$u_3(x_3)$'))+xlab(TeX('$x_3$'))+
  scale_color_manual(values=c('grey','#00BA38','#F8766D','blue'))+
  scale_linetype_manual(values=c("solid","dotdash", "dotted","dashed"))+
  theme(plot.title = element_text(hjust = 0.5),
        legend.title=element_blank(),
        legend.position = c(0.2,0.7),
        legend.key.width = unit(1.2, "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 17))




###setting 3: 
base.f <- function(x){1/(x+0.1)}#corresponds to fid[3]!
logistic.fun <- function(x){
  -13+6*x[1]+24*(x[2]-0.3)^2+1/(x[3]+0.1)-5*cos(5*x[4])+2*exp(4*x[5]-2)
}

#produce points to be plotted
mean.l <- list()
for(k in 1:3){
  load(paste0(path1,"n",k,"_aicsum_DRgen.Rdata"))
  allout <- sumout
  n0 <- n0seq[k]
  n <- nseq[k]
  
  shift.mat <- NULL
  for(run in 1:RUN){
    seed <- 77*run+7
    set.seed(seed)
    alldata <- gendata.logistic(n0=n0,n=n,pi.true=0.4,logistic.fun,dim.X=dim.X)
    pdata <- alldata$pdata;udata <- alldata$udata
    X <- rbind(pdata,udata,udata)
    shift.mat <- rbind(shift.mat,apply(base.f(X), 2, mean))#This is only true for fid[3] component!
  }
  
  feval.l <- list()
  for(run in 1:RUN){
    feval0 <- allout$eval.l[[run]]
    feval.l[[run]] <- feval0+matrix(rep(shift.mat[run,],101),nrow=101,byrow=T)
  }
  
  feval.l2 <- sapply(feval.l, function(x){x[,fid[3]]})
  mean.l[[k]] <- apply(feval.l2, 1, mean)
}

df <- data.frame(t=eval.pts,mean1=mean.l[[1]],mean2=mean.l[[2]],mean3=mean.l[[3]])
df$true <- base.f(df$t)
colnames(df)[2:4] <- c('n=1500','n=3000','n=4500')
df <- df[,c(1,5,2,3,4)]
df2 <- reshape2::melt(df,id.vars = 1,variable.name = 'case')

pgen <- ggplot(data=df2,aes(x=t,y=value,group=case))+
  geom_line(aes(linetype=case, color=case),linewidth=1)+
  ggtitle("Setting 3")+ylab(TeX('$u_3(x_3)$'))+xlab(TeX('$x_3$'))+
  scale_color_manual(values=c('grey','#00BA38','#F8766D','blue'))+
  scale_linetype_manual(values=c("solid","dotdash", "dotted","dashed"))+
  theme(plot.title = element_text(hjust = 0.5),
        legend.title=element_blank(),
        legend.position = c(0.7,0.7),
        legend.key.width = unit(1.2, "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 17))

p1 <- gridExtra::grid.arrange(pquad,pgen,nrow=1)
ggsave(filename='funplot.pdf',plot=p1,device="pdf")
