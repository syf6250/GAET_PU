library(mgcv)
library(doParallel)
library(gss)
library(fda)
source("FUN-para.R")
source("FUN-real.R")
library(disdat)
po <- disPo("SWI")
bg <- disBg("SWI")
library(doParallel)
source("FUN-para.R")
source("FUN-testreal.R")
po <- subset(po,spid=='swi28')
#nrow(unique(po[, c("x","y")]))==nrow(po)

#apply(po,2,function(x){length(unique(x))})
pdata <- po[,c(10:11,13,15:17,19,7:9,12,14,18)]
udata <- bg[,c(10:11,13,15:17,19,7:9,12,14,18)]

maxall <- minall <- rep(NA,ncol(pdata))
for(i in 1:ncol(pdata)){
  all <- c(pdata[,i],udata[,i])
  maxall[i] <- max(all);minall[i] <- min(all)
  pdata[,i] <- (pdata[,i]-minall[i])/(maxall[i]-minall[i])
  udata[,i] <- (udata[,i]-minall[i])/(maxall[i]-minall[i])
}
#apply(pdata,2,range);apply(udata,2,range)
names.seq <- NULL
for(i in 1:ncol(pdata)){names.seq <- c(names.seq,paste0("X",i))}
names(pdata) <- names.seq;names(udata) <- names.seq
n1 <- nrow(udata);n0 <- nrow(pdata);n <- n0+n1



######method-linear
testL <- EM.addDR.linear.minit(77, pdata, udata,maxit=1000,thres=1e-4,
                               pi.true=0.4,ini.coef=NULL,each=5)
#testL$conv
testL$pi.est



######method-proposed
nu.seq <- 10^seq(from=-4.5,to=-3.5,by=0.25)
eval.pts <- seq(from=0.05,to=0.95,by=0.01)
eval.df <- matrix(nrow=length(eval.pts),ncol=7)#==dimC
colnames(eval.df) <- colnames(pdata)[1:7]#==dimC
for(j in 1:7){#==dimC
  eval.df[,j] <- eval.pts
}
eval.df <- data.frame(eval.df)


###to determine nu.seq
ptm <- proc.time()
res.l <- list()
cluster <- makeCluster(min(length(nu.seq),detectCores()-2))
registerDoParallel(cluster)
res.l <- foreach(run = 1:length(nu.seq),.packages = c("mgcv"))%dopar%{
  AIC.addDR2(nu.seq[run],7777,pdata,udata,maxit=1000,thres=1e-4,
             pi.true=0.4,k=10,ord=3,ord.pen=2,eval.df=eval.df,
             initinput=NULL,dimC=7,dimL=6,each=5)
}
stopCluster(cl = cluster)
proc.time() - ptm#s
#save(res.l,file="treePOswi28_dimC7dimL6_list_-4.5_-3.5_noinit.Rdata")
#sapply(res.l, function(x){x$conv})
aicall <- do.call(c,lapply(res.l, function(x){x$icall$aic}))#AIC tuning selection
#aicall
#do.call(c,lapply(res.l, function(x){x$icall$edf}))
tuneid <- which.min(aicall)
test <- res.l[[tuneid]]
test$pi#pi estimation


#############bootstrap CI
B <- 2000
###design integral range and evaluation points (outside the loop)
nodes.mat <- weights.mat <-  matrix(nrow=20,ncol=7)#ncol==dimC
bsval.l <- list()
for(j in 1:7){#==dimC
  qt <- quantile(pdata[,j],probs=c(0.05,0.95))
  quad <- gauss.quad(size = 20, interval = qt)
  nodes.mat[,j] <- quad$pt
  weights.mat[,j] <- quad$wt
  bsval.l[[j]] <-  fourier(quad$pt,11,diff(range(qt)))
}
colnames(nodes.mat) <- colnames(weights.mat) <- names(pdata)[1:7]#==dimC
eval.df.test <- as.data.frame(nodes.mat)


ptm <- proc.time()
cluster <- makeCluster(min(120,detectCores()-2))
registerDoParallel(cluster)
res.l <- list()
res.l <- foreach(run = 1:B,.packages = c("mgcv"))%dopar%{
  npBoot.aic.core(run,init=test$init,full.obj=test$gamObj,
                  nu.seq=nu.seq[tuneid],seed=777,pdata,udata,maxit=1000,thres=1e-4,
                  pi.true=NULL,k=10,ord=3,ord.pen=2,eval.df=eval.df,
                  dimC=7,dimL=6,each=NULL,eval.df.test=eval.df.test)
}
stopCluster(cl = cluster)
proc.time() - ptm#
save(res.l,file="boot_treePOswi28_dimC7dimL6_list_-4.5_-3.5_notune.Rdata")


###collect results
sumCI <- getCIest(test,res.l,dimC=7,dimL=6,pdata,
                  eval.df.test=eval.df.test,weights.mat=weights.mat,bsval.l=bsval.l,B=B)
#sumCI$conv.prop
sumCI$CI.pi['95',c(3:4)]#confidence interval
sumCI$sigtest#T:refuse; F: accept
sumCI$pval#p-value



#############plot
library(ggplot2)
X <- rbind(pdata,udata,udata)

j=1
plotdomain <- c(10:85)
x.pts <- (test$eval.pts[,j]*(maxall[j]-minall[j])+minall[j])[plotdomain]
u.val <- test$eval.terms[plotdomain,j]
intcpt <- -testL$coef.est[j+1]*mean(X[,j])
u.val.Lin <- (test$eval.pts[plotdomain,j])*testL$coef.est[j+1]+intcpt
u.ci.upper <- sumCI$CI.fun95[[j]][plotdomain,4]
u.ci.lower <- sumCI$CI.fun95[[j]][plotdomain,3]
df <- data.frame(x = x.pts,u = u.val,u_lin = u.val.Lin,
                 u_ci_lower = u.ci.lower,u_ci_upper = u.ci.upper)

uplot.ddeg <- ggplot(df, aes(x = x)) +
  #geom_ribbon(aes(ymin = u_ci_lower, ymax = u_ci_upper),alpha = 0.2) +
  # curve 1: u.val
  geom_line(aes(y = u, color = "GAET", linetype = "GAET"),linewidth = 1) +
  # curve 2: u.val.Lin
  geom_line(aes(y = u_lin, color = "linear", linetype = "linear"),linewidth = 1) +
  scale_color_manual(values = c("GAET" = "#1f77b4", "linear" = "#d62728")) +
  scale_linetype_manual(values = c("GAET" = "solid", "linear" = "dashed")) +
  labs(x = "ddeg", y = expression(hat(u)(ddeg)),color = NULL, linetype = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = c(0.7,0.2),legend.key.width  = unit(1.5, "cm"),
        legend.title = element_text(size = 14),legend.text  = element_text(size = 10),
        axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
        axis.text.x  = element_text(size = 10),axis.text.y  = element_text(size = 10))


j=2
plotdomain <- c(10:88)
x.pts <- (test$eval.pts[,j]*(maxall[j]-minall[j])+minall[j])[plotdomain]
u.val <- test$eval.terms[plotdomain,j]
intcpt <- -testL$coef.est[j+1]*mean(X[,j])
u.val.Lin <- (test$eval.pts[plotdomain,j])*testL$coef.est[j+1]+intcpt
u.ci.upper <- sumCI$CI.fun95[[j]][plotdomain,4]
u.ci.lower <- sumCI$CI.fun95[[j]][plotdomain,3]
df <- data.frame(x = x.pts,u = u.val,u_lin = u.val.Lin,
                 u_ci_lower = u.ci.lower,u_ci_upper = u.ci.upper)
uplot.nutri <- ggplot(df, aes(x = x)) +
  #geom_ribbon(aes(ymin = u_ci_lower, ymax = u_ci_upper),alpha = 0.2) +
  # curve 1: u.val
  geom_line(aes(y = u, color = "GAET", linetype = "GAET"),linewidth = 1) +
  # curve 2: u.val.Lin
  geom_line(aes(y = u_lin, color = "linear", linetype = "linear"),linewidth = 1) +
  scale_color_manual(values = c("GAET" = "#1f77b4", "linear" = "#d62728")) +
  scale_linetype_manual(values = c("GAET" = "solid", "linear" = "dashed")) +
  labs(x = "nutri", y = expression(hat(u)(nutri)),color = NULL, linetype = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position="none",
        axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
        axis.text.x  = element_text(size = 10),axis.text.y  = element_text(size = 10))


j=3
plotdomain <- c(11:80)
x.pts <- (test$eval.pts[,j]*(maxall[j]-minall[j])+minall[j])[plotdomain]
u.val <- test$eval.terms[plotdomain,j]
intcpt <- -testL$coef.est[j+1]*mean(X[,j])
u.val.Lin <- (test$eval.pts[plotdomain,j])*testL$coef.est[j+1]+intcpt
u.ci.upper <- sumCI$CI.fun95[[j]][plotdomain,4]
u.ci.lower <- sumCI$CI.fun95[[j]][plotdomain,3]
df <- data.frame(x = x.pts,u = u.val,u_lin = u.val.Lin,
                 u_ci_lower = u.ci.lower,u_ci_upper = u.ci.upper)

uplot.precyy <- ggplot(df, aes(x = x)) +
  #geom_ribbon(aes(ymin = u_ci_lower, ymax = u_ci_upper),alpha = 0.2) +
  # curve 1: u.val
  geom_line(aes(y = u, color = "GAET", linetype = "GAET"),linewidth = 1) +
  # curve 2: u.val.Lin
  geom_line(aes(y = u_lin, color = "linear", linetype = "linear"),linewidth = 1) +
  scale_color_manual(values = c("GAET" = "#1f77b4", "linear" = "#d62728")) +
  scale_linetype_manual(values = c("GAET" = "solid", "linear" = "dashed")) +
  labs(x = "precyy", y = expression(hat(u)(precyy)),color = NULL, linetype = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position="none",
        axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
        axis.text.x  = element_text(size = 10),axis.text.y  = element_text(size = 10))


j=4
plotdomain <- c(1:75)
x.pts <- (test$eval.pts[,j]*(maxall[j]-minall[j])+minall[j])[plotdomain]
u.val <- test$eval.terms[plotdomain,j]
intcpt <- -testL$coef.est[j+1]*mean(X[,j])
u.val.Lin <- (test$eval.pts[plotdomain,j])*testL$coef.est[j+1]+intcpt
u.ci.upper <- sumCI$CI.fun95[[j]][plotdomain,4]
u.ci.lower <- sumCI$CI.fun95[[j]][plotdomain,3]
df <- data.frame(x = x.pts,u = u.val,u_lin = u.val.Lin,
                 u_ci_lower = u.ci.lower,u_ci_upper = u.ci.upper)

uplot.slope <- ggplot(df, aes(x = x)) +
  #geom_ribbon(aes(ymin = u_ci_lower, ymax = u_ci_upper),alpha = 0.2) +
  # curve 1: u.val
  geom_line(aes(y = u, color = "GAET", linetype = "GAET"),linewidth = 1) +
  # curve 2: u.val.Lin
  geom_line(aes(y = u_lin, color = "linear", linetype = "linear"),linewidth = 1) +
  scale_color_manual(values = c("GAET" = "#1f77b4", "linear" = "#d62728")) +
  scale_linetype_manual(values = c("GAET" = "solid", "linear" = "dashed")) +
  labs(x = "slope", y = expression(hat(u)(slope)),color = NULL, linetype = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position="none",
        axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
        axis.text.x  = element_text(size = 10),axis.text.y  = element_text(size = 10))


j=5
plotdomain <- c(20:90)
x.pts <- (test$eval.pts[,j]*(maxall[j]-minall[j])+minall[j])[plotdomain]
u.val <- test$eval.terms[plotdomain,j]
intcpt <- -testL$coef.est[j+1]*mean(X[,j])
u.val.Lin <- (test$eval.pts[plotdomain,j])*testL$coef.est[j+1]+intcpt
u.ci.upper <- sumCI$CI.fun95[[j]][plotdomain,4]
u.ci.lower <- sumCI$CI.fun95[[j]][plotdomain,3]
df <- data.frame(x = x.pts,u = u.val,u_lin = u.val.Lin,
                 u_ci_lower = u.ci.lower,u_ci_upper = u.ci.upper)

uplot.sradyy <- ggplot(df, aes(x = x)) +
  #geom_ribbon(aes(ymin = u_ci_lower, ymax = u_ci_upper),alpha = 0.2) +
  # curve 1: u.val
  geom_line(aes(y = u, color = "GAET", linetype = "GAET"),linewidth = 1) +
  # curve 2: u.val.Lin
  geom_line(aes(y = u_lin, color = "linear", linetype = "linear"),linewidth = 1) +
  scale_color_manual(values = c("GAET" = "#1f77b4", "linear" = "#d62728")) +
  scale_linetype_manual(values = c("GAET" = "solid", "linear" = "dashed")) +
  labs(x = "sradyy", y = expression(hat(u)(sradyy)),color = NULL, linetype = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position="none",
        axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
        axis.text.x  = element_text(size = 10),axis.text.y  = element_text(size = 10))

j=6
plotdomain <- c(5:90)
x.pts <- (test$eval.pts[,j]*(maxall[j]-minall[j])+minall[j])[plotdomain]
u.val <- test$eval.terms[plotdomain,j]
intcpt <- -testL$coef.est[j+1]*mean(X[,j])
u.val.Lin <- (test$eval.pts[plotdomain,j])*testL$coef.est[j+1]+intcpt
u.ci.upper <- sumCI$CI.fun95[[j]][plotdomain,4]
u.ci.lower <- sumCI$CI.fun95[[j]][plotdomain,3]
df <- data.frame(x = x.pts,u = u.val,u_lin = u.val.Lin,
                 u_ci_lower = u.ci.lower,u_ci_upper = u.ci.upper)

uplot.swb <- ggplot(df, aes(x = x)) +
  #geom_ribbon(aes(ymin = u_ci_lower, ymax = u_ci_upper),alpha = 0.2) +
  # curve 1: u.val
  geom_line(aes(y = u, color = "GAET", linetype = "GAET"),linewidth = 1) +
  # curve 2: u.val.Lin
  geom_line(aes(y = u_lin, color = "linear", linetype = "linear"),linewidth = 1) +
  scale_color_manual(values = c("GAET" = "#1f77b4", "linear" = "#d62728")) +
  scale_linetype_manual(values = c("GAET" = "solid", "linear" = "dashed")) +
  labs(x = "swb", y = expression(hat(u)(swb)),color = NULL, linetype = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position="none",
        axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
        axis.text.x  = element_text(size = 10),axis.text.y  = element_text(size = 10))


j=7
plotdomain <- c(15:90)
x.pts <- (test$eval.pts[,j]*(maxall[j]-minall[j])+minall[j])[plotdomain]
u.val <- test$eval.terms[plotdomain,j]
intcpt <- -testL$coef.est[j+1]*mean(X[,j])
u.val.Lin <- (test$eval.pts[plotdomain,j])*testL$coef.est[j+1]+intcpt
u.ci.upper <- sumCI$CI.fun95[[j]][plotdomain,4]
u.ci.lower <- sumCI$CI.fun95[[j]][plotdomain,3]
df <- data.frame(x = x.pts,u = u.val,u_lin = u.val.Lin,
                 u_ci_lower = u.ci.lower,u_ci_upper = u.ci.upper)

uplot.topo <- ggplot(df, aes(x = x)) +
  #geom_ribbon(aes(ymin = u_ci_lower, ymax = u_ci_upper),alpha = 0.2) +
  # curve 1: u.val
  geom_line(aes(y = u, color = "GAET", linetype = "GAET"),linewidth = 1) +
  # curve 2: u.val.Lin
  geom_line(aes(y = u_lin, color = "linear", linetype = "linear"),linewidth = 1) +
  scale_color_manual(values = c("GAET" = "#1f77b4", "linear" = "#d62728")) +
  scale_linetype_manual(values = c("GAET" = "solid", "linear" = "dashed")) +
  labs(x = "topo", y = expression(hat(u)(topo)),color = NULL, linetype = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position="none",
        axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
        axis.text.x  = element_text(size = 10),axis.text.y  = element_text(size = 10))

p7 <- gridExtra::grid.arrange(uplot.ddeg,uplot.nutri,uplot.precyy,
                              uplot.slope,uplot.sradyy,uplot.swb,uplot.topo,ncol=2)

ggsave("treePO_uplot.pdf",p7,device = "pdf")