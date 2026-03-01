library(mgcv)
source("FUN-gendata.R")
library(reshape2)
library(ggplot2)
library(latex2exp)
#need to use outputs of 'Tab1_Set1.R', 'Tab1_Set2.R', 'Tab1_Set3.R'
#suppose now they are in the folder "results"
path1 <- "results/"

prop <- 5
n1seq <- c(250,500,750,1250)
n0seq <- prop*n1seq
nseq <- n0seq+n1seq

dim.X <- 5
rr <- 0.15
pirange <- c(0.4-rr, 0.4+rr)

for(nid in 1:4){
  ###setting 2: 
  load(paste0(path1,"n",nid,"_aicsum_DRquad.Rdata"))
  df <- data.frame(pi=sumout$pi.all)
  p2 <- ggplot(data=df,aes(x=pi)) + geom_histogram(color="black", fill="white")+
    geom_vline(aes(xintercept=mean(pi)),color="blue", linetype="dashed", size=1)+
    geom_vline(aes(xintercept=0.4),color="#FF6666", linetype="solid", size=1)+
    scale_x_continuous(limits = pirange) + 
    ggtitle("Setting 2")+xlab(TeX('$pi$'))+
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size = 16))
  
  
  ###setting 3: 
  load(paste0(path1,"n",nid,"_aicsum_DRgen.Rdata"))
  df <- data.frame(pi=sumout$pi.all)
  p3 <- ggplot(data=df,aes(x=pi)) + geom_histogram(color="black", fill="white")+
    geom_vline(aes(xintercept=mean(pi)),color="blue", linetype="dashed", size=1)+
    geom_vline(aes(xintercept=0.4),color="#FF6666", linetype="solid", size=1)+
    scale_x_continuous(limits = pirange) + 
    ggtitle("Setting 3")+xlab(TeX('$pi$'))+
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size = 16))
  
  
  ###setting 1: 
  load(paste0(path1,"n",nid,"_aicsum_DRlin.Rdata"))
  df <- data.frame(pi=sumout$pi.all)
  p1 <- ggplot(data=df,aes(x=pi)) + geom_histogram(color="black", fill="white")+
    geom_vline(aes(xintercept=mean(pi)),color="blue", linetype="dashed", size=1)+
    geom_vline(aes(xintercept=0.4),color="#FF6666", linetype="solid", size=1)+
    scale_x_continuous(limits = pirange) + 
    ggtitle("Setting 1")+xlab(TeX('$pi$'))+
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size = 16))
  
  ppp <- gridExtra::grid.arrange(p1,p2,p3,nrow=1)
  ggsave(filename=paste0('hist2000-n',nseq[nid],'.pdf'),plot=ppp,device="pdf")
}


