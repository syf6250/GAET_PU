library(reshape2)
library(ggplot2)
#need to use outputs of 'FigS1S2_Set1.R', 'FigS1S2_Set2.R', 'FigS1S2_Set3.R'
#suppose now they are in the folder "results"
path1 <- "results/"

prop <- 5
n1seq <- c(250,500,750,1250)
n0seq <- prop*n1seq
nseq <- n0seq+n1seq


###Setting 2
load(paste0(path1,"sensi_DRquad_Run2000_mseR.Rdata"))
k_vals <- -10:-6

df <- melt(mse.ratio[,1:5], varnames = c("row", "col"), value.name = "value")
df$n <- nseq[df$row]
df$k <- k_vals[df$col]
df$n <- as.factor(df$n)
p.mse.2 <- ggplot(df, aes(x = k, y = value, group = n, color = n,
                          linetype = n, shape = n)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(x = expression(italic("l")), y = "mse") +ggtitle("Setting 2")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5, size = 19, face = "bold"),
                   legend.position = c(0.3, 0.7),
                   legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
                   legend.key = element_rect(fill = "white"),legend.key.width  = unit(1.5, "cm"),
                   legend.title = element_text(size = 18),legend.text  = element_text(size = 15),
                   axis.title.x = element_text(size = 18),axis.title.y = element_text(size = 18),
                   axis.text.x  = element_text(size = 15),axis.text.y  = element_text(size = 15))


###Setting 3
load(paste0(path1,"sensi_DRgen_Run2000_mseR.Rdata"))
k_vals <- -10:-6

df <- melt(mse.ratio[,1:5], varnames = c("row", "col"), value.name = "value")
df$n <- nseq[df$row]
df$l <- k_vals[df$col]
df$n <- as.factor(df$n)
p.mse.3 <- ggplot(df, aes(x = l, y = value, group = n, color = n,
                          linetype = n, shape = n)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(x = expression(italic("l")), y = "mse") +ggtitle("Setting 3")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5, size = 19, face = "bold"),
                   legend.position = c(0.3, 0.7),
                   legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
                   legend.key = element_rect(fill = "white"),legend.key.width  = unit(1.5, "cm"),
                   legend.title = element_text(size = 18),legend.text  = element_text(size = 15),
                   axis.title.x = element_text(size = 18),axis.title.y = element_text(size = 18),
                   axis.text.x  = element_text(size = 15),axis.text.y  = element_text(size = 15))


###Setting 1
load(paste0(path1,"sensi_DRlin_Run2000_mseR.Rdata"))
k_vals <- -4:0

df <- melt(mse.ratio[,3:7], varnames = c("row", "col"), value.name = "value")
df$n <- nseq[df$row]
df$k <- k_vals[df$col]
df$n <- as.factor(df$n)
p.mse.1 <- ggplot(df, aes(x = k, y = value, group = n, color = n,
                          linetype = n, shape = n)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(x = expression(italic("l")), y = "mse") +ggtitle("Setting 1")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5, size = 19, face = "bold"),
                   legend.position = c(0.7, 0.7),
                   legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
                   legend.key = element_rect(fill = "white"),legend.key.width  = unit(1.5, "cm"),
                   legend.title = element_text(size = 18),legend.text  = element_text(size = 15),
                   axis.title.x = element_text(size = 18),axis.title.y = element_text(size = 18),
                   axis.text.x  = element_text(size = 15),axis.text.y  = element_text(size = 15))

p.mse.all <- gridExtra::grid.arrange(p.mse.1,p.mse.2,p.mse.3,nrow=1)
ggsave("sensi_mse.pdf",p.mse.all,device = "pdf")




###Err
###Setting 2
load(paste0(path1,"sensi_DRquad_Run2000_errR.Rdata"))
k_vals <- -9:-5

df <- melt(err.ratio[,2:6], varnames = c("row", "col"), value.name = "value")
df$n <- nseq[df$row]
df$k <- k_vals[df$col]
df$n <- as.factor(df$n)
p.err.2 <- ggplot(df, aes(x = k, y = value, group = n, color = n,
                          linetype = n, shape = n)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(x = expression(italic("l")), y = "Err") +ggtitle("Setting 2")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5, size = 19, face = "bold"),
                   legend.position = c(0.3, 0.7),
                   legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
                   legend.key = element_rect(fill = "white"),legend.key.width  = unit(1.5, "cm"),
                   legend.title = element_text(size = 18),legend.text  = element_text(size = 15),
                   axis.title.x = element_text(size = 18),axis.title.y = element_text(size = 18),
                   axis.text.x  = element_text(size = 15),axis.text.y  = element_text(size = 15))


###Setting 3
load(paste0(path1,"sensi_DRgen_Run2000_errR.Rdata"))
k_vals <- -9:-5

df <- melt(err.ratio[,2:6], varnames = c("row", "col"), value.name = "value")
df$n <- nseq[df$row]
df$k <- k_vals[df$col]
df$n <- as.factor(df$n)
p.err.3 <- ggplot(df, aes(x = k, y = value, group = n, color = n,
                          linetype = n, shape = n)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(x = expression(italic("l")), y = "Err") +ggtitle("Setting 3")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5, size = 19, face = "bold"),
                   legend.position = c(0.3, 0.7),
                   legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
                   legend.key = element_rect(fill = "white"),legend.key.width  = unit(1.5, "cm"),
                   legend.title = element_text(size = 18),legend.text  = element_text(size = 15),
                   axis.title.x = element_text(size = 18),axis.title.y = element_text(size = 18),
                   axis.text.x  = element_text(size = 15),axis.text.y  = element_text(size = 15))

###Setting 1
load(paste0(path1,"sensi_DRlin_Run2000_errR.Rdata"))
k_vals <- -4:0

df <- melt(err.ratio[,3:7], varnames = c("row", "col"), value.name = "value")
df$n <- nseq[df$row]
df$k <- k_vals[df$col]
df$n <- as.factor(df$n)
p.err.1 <- ggplot(df, aes(x = k, y = value, group = n, color = n,
                          linetype = n, shape = n)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(x = expression(italic("l")), y = "Err") +ggtitle("Setting 1")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5, size = 19, face = "bold"),
                   legend.position = c(0.7, 0.7),
                   legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
                   legend.key = element_rect(fill = "white"),legend.key.width  = unit(1.5, "cm"),
                   legend.title = element_text(size = 18),legend.text  = element_text(size = 15),
                   axis.title.x = element_text(size = 18),axis.title.y = element_text(size = 18),
                   axis.text.x  = element_text(size = 15),axis.text.y  = element_text(size = 15))

p.err.all <- gridExtra::grid.arrange(p.err.1,p.err.2,p.err.3,nrow=1)
ggsave("sensi_err.pdf",p.err.all,device = "pdf")
