#need to use outputs of 'Tab1_Set1.R', 'Tab1_Set2.R', 'Tab1_Set3.R'
#suppose now they are in the folder "results"
path1 <- "results/"
out <- matrix(nrow=12,ncol=4)


for(j in 1:4){
  load(paste0(path1,"n",j,"_aicsum_DRlin.Rdata"))
  out[j,] <- table(factor(sumout$init.table,levels = 3:6))
}

for(j in 1:4){
  load(paste0(path1,"n",j,"_aicsum_DRquad.Rdata"))
  out[j+4,] <- table(factor(sumout$init.table,levels = 3:6))
}

for(j in 1:4){
  load(paste0(path1,"n",j,"_aicsum_DRgen.Rdata"))
  out[j+8,] <- table(factor(sumout$init.table,levels = 3:6))
}

print(out/2000)


