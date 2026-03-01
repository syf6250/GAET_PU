###normal test for estimations of pi
library(ggplot2)
pval.mat <- matrix(nrow=3,ncol=4)

#need to use outputs of 'Tab1_Set1.R', 'Tab1_Set2.R', 'Tab1_Set3.R'
#suppose now they are in the folder "results"
path1 <- "results/"

#setting 1

pi.list <- list()
for(j in 1:4){
  load(paste0(path1,"n",j,"_aicsum_DRlin.Rdata"))
  pi.list[[j]] <- sumout$pi.all
}
sapply(pi.list, function(x){shapiro.test(x)$p.value})
pval.mat[1,] <- sapply(pi.list, function(x){shapiro.test(x)$p.value})
#qqnorm(pi.list[[1]])


qq1 <-  ggplot(data.frame(x = pi.list[[1]]), aes(sample = x)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "Setting 1 (n = 1500)",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  )+theme(plot.title = element_text(hjust = 0.5))#+coord_cartesian(xlim = c(-2.5, 2.5))
qq2 <-  ggplot(data.frame(x = pi.list[[2]]), aes(sample = x)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "Setting 1 (n = 3000)",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  )+theme(plot.title = element_text(hjust = 0.5))#+coord_cartesian(xlim = c(-2.5, 2.5))
qq3 <-  ggplot(data.frame(x = pi.list[[3]]), aes(sample = x)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "Setting 1 (n = 4500)",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  )+theme(plot.title = element_text(hjust = 0.5))#+coord_cartesian(xlim = c(-2.5, 2.5))
qq4 <-  ggplot(data.frame(x = pi.list[[4]]), aes(sample = x)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "Setting 1 (n = 7500)",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  )+theme(plot.title = element_text(hjust = 0.5))#+coord_cartesian(xlim = c(-2.5, 2.5))

qq.all <- gridExtra::grid.arrange(qq1,qq2,qq3,qq4,nrow=1)
ggsave(filename=paste0('Set',1,'qqplot','.pdf'),plot=qq.all,device="pdf")






#setting 2

pi.list <- list()
for(j in 1:4){
  load(paste0(path1,"n",j,"_aicsum_DRquad.Rdata"))
  pi.list[[j]] <- sumout$pi.all
}
sapply(pi.list, function(x){shapiro.test(x)$p.value})
pval.mat[2,] <- sapply(pi.list, function(x){shapiro.test(x)$p.value})
#qqnorm(pi.list[[1]])


qq1 <-  ggplot(data.frame(x = pi.list[[1]]), aes(sample = x)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "Setting 2 (n = 1500)",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  )+theme(plot.title = element_text(hjust = 0.5))#+coord_cartesian(xlim = c(-2.5, 2.5))
qq2 <-  ggplot(data.frame(x = pi.list[[2]]), aes(sample = x)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "Setting 2 (n = 3000)",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  )+theme(plot.title = element_text(hjust = 0.5))#+coord_cartesian(xlim = c(-2.5, 2.5))
qq3 <-  ggplot(data.frame(x = pi.list[[3]]), aes(sample = x)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "Setting 2 (n = 4500)",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  )+theme(plot.title = element_text(hjust = 0.5))#+coord_cartesian(xlim = c(-2.5, 2.5))
qq4 <-  ggplot(data.frame(x = pi.list[[4]]), aes(sample = x)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "Setting 2 (n = 7500)",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  )+theme(plot.title = element_text(hjust = 0.5))#+coord_cartesian(xlim = c(-2.5, 2.5))

qq.all <- gridExtra::grid.arrange(qq1,qq2,qq3,qq4,nrow=1)
ggsave(filename=paste0('Set',2,'qqplot','.pdf'),plot=qq.all,device="pdf")








#setting 3

pi.list <- list()
for(j in 1:4){
  load(paste0(path1,"n",j,"_aicsum_DRgen.Rdata"))
  pi.list[[j]] <- sumout$pi.all
}
sapply(pi.list, function(x){shapiro.test(x)$p.value})
pval.mat[3,] <- sapply(pi.list, function(x){shapiro.test(x)$p.value})
#qqnorm(pi.list[[1]])


qq1 <-  ggplot(data.frame(x = pi.list[[1]]), aes(sample = x)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "Setting 3 (n = 1500)",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  )+theme(plot.title = element_text(hjust = 0.5))#+coord_cartesian(xlim = c(-2.5, 2.5))
qq2 <-  ggplot(data.frame(x = pi.list[[2]]), aes(sample = x)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "Setting 3 (n = 3000)",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  )+theme(plot.title = element_text(hjust = 0.5))#+coord_cartesian(xlim = c(-2.5, 2.5))
qq3 <-  ggplot(data.frame(x = pi.list[[3]]), aes(sample = x)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "Setting 3 (n = 4500)",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  )+theme(plot.title = element_text(hjust = 0.5))#+coord_cartesian(xlim = c(-2.5, 2.5))
qq4 <-  ggplot(data.frame(x = pi.list[[4]]), aes(sample = x)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "Setting 3 (n = 7500)",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  )+theme(plot.title = element_text(hjust = 0.5))#+coord_cartesian(xlim = c(-2.5, 2.5))

qq.all <- gridExtra::grid.arrange(qq1,qq2,qq3,qq4,nrow=1)
ggsave(filename=paste0('Set',3,'qqplot','.pdf'),plot=qq.all,device="pdf")



print(pval.mat)
