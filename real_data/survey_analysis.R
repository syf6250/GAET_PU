library(mgcv)
library(fastDummies)
library(ggplot2)
library(doParallel)
library(gss)
library(fda)
source("FUN-para.R")
source("FUN-real.R")

#############data
raw <- haven::read_dta("ZA5688_v6-0-0.dta")#downloaded from https://search.gesis.org/research_data/ZA5688
morale_df <- data.frame(a1=raw$qe20_1,a2=raw$qe20_2,a3=raw$qe20_3,a4=raw$qe20_4,
                        a5=raw$qe20_5,a6=raw$qe20_6,a7=raw$qe20_7)
Tax_morale <- apply(morale_df, 1, mean,na.rm=T)
rawdata <- data.frame(age=raw$d11,tax_morale=Tax_morale,gender=raw$d10,#urban=raw$d25,
                      occupation=raw$d15a_r1,financial_problem=raw$d60,detection_risk=raw$qe3,
                      y=as.numeric(raw$qe14))
###delete missing cases
pid <- which(rawdata$y==1)
rawdata$y[-pid] <- 0
id2del <- which(is.na(rawdata$detection_risk)==T)
rawdata <- rawdata[-id2del,]
id2del <- which(is.na(rawdata$financial_problem)==T)
rawdata <- rawdata[-id2del,]
id2del <- which(is.na(rawdata$tax_morale)==T)
rawdata <- rawdata[-id2del,]
#sum(complete.cases(rawdata))==nrow(rawdata)

###generate dummy
rawdata$age <- as.numeric(rawdata$age)
rawdata$gender <- factor(rawdata$gender,levels=c("1", "2"),labels=c("M","F"))
rawdata$occupation <- factor(rawdata$occupation,levels=c("1","2","3"),labels=c("Self_employed","Employed","Unemployed"))
rawdata$financial_problem <- factor(rawdata$financial_problem,levels=c("1","2","3"),labels=c("Most","Occasional","None"))
rawdata$detection_risk <- factor(rawdata$detection_risk,levels=c("1","2","3","4"),labels=c("High","high","small","Small"))
rawdata2 <- fastDummies::dummy_cols(rawdata,
                                    select_columns = c("gender","occupation","financial_problem","detection_risk"),
                                    remove_selected_columns=T,remove_most_frequent_dummy=T)

###standardization
agerange <- c(min(rawdata2$age),max(rawdata2$age))
rawdata2$age <- (rawdata2$age-min(rawdata2$age))/(max(rawdata2$age)-min(rawdata2$age))
rawdata2$tax_morale <- (rawdata2$tax_morale-min(rawdata2$tax_morale))/(max(rawdata2$tax_morale)-min(rawdata2$tax_morale))
###split
raw.pdata <- subset(rawdata2,y==1)[,-c(3)]#1145
raw.udata <- subset(rawdata2,y==0)[,-c(3)]#23086
###set up
n1 <- nrow(raw.udata);n0 <- nrow(raw.pdata);n <- n0+n1
names.seq <- NULL
for(j in 1:ncol(raw.pdata)){names.seq <- c(names.seq,paste0("X",j))}
pdata <- raw.pdata;names(pdata) <- names.seq
udata <- raw.udata;names(udata) <- names.seq 






#############estimation
######method-linear
testL <- EM.addDR.linear.minit(77, pdata, udata,maxit=1000,thres=1e-4,
                               pi.true=0.2,ini.coef=NULL,each=5)
#testL$conv
testL$pi.est#


######method-proposed
nu.seq <- 10^seq(from=-5,to=-3,by=0.5)
eval.pts <- seq(from=0.05,to=0.95,by=0.01)
eval.df <- data.frame(X1=eval.pts)

ptm <- proc.time()
res.l <- list()
cluster <- makeCluster(min(length(nu.seq),detectCores()-2))
registerDoParallel(cluster)
res.l <- foreach(run = 1:length(nu.seq),.packages = c("mgcv"))%dopar%{
  AIC.addDR2(nu.seq[run],7777,pdata,udata,maxit=1000,thres=1e-4,
             pi.true=0.2,k=10,ord=3,ord.pen=2,eval.df=eval.df,
             initinput=NULL,dimC=1,dimL=9,each=5)
}
stopCluster(cl = cluster)
proc.time() - ptm#
#save(res.l,file="survey_-5_-3_noinit.Rdata")
#sapply(res.l, function(x){x$conv})
aicall <- do.call(c,lapply(res.l, function(x){x$icall$aic}))#AIC tuning selection
#aicall
tuneid <- which.min(aicall)
test <- res.l[[tuneid]]
test$pi#estimated pi
test$pi*(n1/n)+n0/n#prevalence
X <- rbind(pdata,udata,udata)



#############bootstrap CI
B <- 2000
###design integral range and evaluation points (outside the loop)
nodes.mat <- weights.mat <-  matrix(nrow=20,ncol=1)#ncol==dimC
bsval.l <- list()
for(j in 1:1){#==dimC
  qt <- c(0,0.9)
  quad <- gauss.quad(size = 20, interval = qt)
  nodes.mat[,j] <- quad$pt
  weights.mat[,j] <- quad$wt
  bsval.l[[j]] <-  fourier(quad$pt,11,diff(range(qt)))
}
colnames(nodes.mat) <- colnames(weights.mat) <- names(pdata)[1:1]#==dimC
eval.df.test <- as.data.frame(nodes.mat)



ptm <- proc.time()
cluster <- makeCluster(min(120,detectCores()-2))
registerDoParallel(cluster)
res.l <- list()
res.l <- foreach(run = 1:B,.packages = c("mgcv"))%dopar%{
  npBoot.aic.core(run,init=test$init,full.obj=test$gamObj,
                  nu.seq=nu.seq[tuneid],seed=777,pdata,udata,maxit=1000,thres=1e-4,
                  pi.true=NULL,k=10,ord=3,ord.pen=2,eval.df=eval.df,
                  dimC=1,dimL=9,each=NULL,eval.df.test=eval.df.test)
}
stopCluster(cl = cluster)
proc.time() - ptm# 
#save(res.l,file="boot_survey_-5_-3_notune_shortrange.Rdata")

###collect results 
sumCI <- getCIest(test,res.l,dimC=1,dimL=9,pdata,
                  eval.df.test=eval.df.test,weights.mat=weights.mat,bsval.l=bsval.l,B=B)
#sumCI$conv.prop
sumCI$CI.pi['95',c(3,4)]#confidence interval
sumCI$sigtest#T:refuse; F: accept
sumCI$pval#p-value



#############function plot
plotdomain <- c(1:80)
age.pts <- (test$eval.pts*diff(agerange)+agerange[1])[plotdomain,1]
u.val <- test$eval.terms[plotdomain,1]
u.ci.upper <- sumCI$CI.fun95[[1]][plotdomain,4]
u.ci.lower <- sumCI$CI.fun95[[1]][plotdomain,3]
intcpt <- -testL$coef.est[2]*mean(X$X1)
u.val.Lin <- (test$eval.pts[plotdomain,1])*testL$coef.est[2]+intcpt
df <- data.frame(
  age = age.pts,
  u = u.val,
  u_ci_lower = u.ci.lower,
  u_ci_upper = u.ci.upper,
  u_lin = u.val.Lin
)
uplot <- ggplot(df, aes(x = age)) +
  #geom_ribbon(aes(ymin = u_ci_lower, ymax = u_ci_upper),alpha = 0.2) +
  geom_line(aes(y = u, color = "GAET", linetype = "GAET"),
            linewidth = 1) +
  geom_line(aes(y = u_lin, color = "linear", linetype = "linear"),
            linewidth = 1) +
  scale_color_manual(values = c("GAET" = "#1f77b4", "linear" = "#d62728")) +
  scale_linetype_manual(values = c("GAET" = "solid", "linear" = "dashed")) +
  labs(x = "age", y = expression(hat(u)(age)),
       color = NULL, linetype = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = c(0.2,0.3),legend.key.width  = unit(1.5, "cm"),
        legend.title = element_text(size = 16),legend.text  = element_text(size = 12),
        axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16),
        axis.text.x  = element_text(size = 12),axis.text.y  = element_text(size = 12))
ggsave("survey_uplot.pdf",uplot,device = "pdf")

