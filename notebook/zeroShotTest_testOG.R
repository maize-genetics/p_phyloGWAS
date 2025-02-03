rm(list = ls())
# Use reticulate to import numpy
Sys.setenv(RETICULATE_PYTHON = "/usr/bin/python")
library(reticulate)
library(Biostrings)
library(parallel)

softmax <- function(logit) {
  p <- exp(logit) / sum(exp(logit))
  return(p)
}

py_config()
np <- import("numpy")

# load numpy
npz <- np$load("/workdir/hackathon_202406/output/finalTestOG_ESM_logit.npz")
test_list <- lapply(npz$files, function(x) npz[[x]])
names(test_list) <- npz$files

# softmax
prob = mclapply(test_list,function(x) {
  tmp = t(apply(x,1,softmax))
  colnames(tmp) = limma::strsplit2("LAGVSERTIDPKQNFYMHWC","")
  return(tmp)
},mc.cores = 5)

# load fasta
testfasta = readAAStringSet("/workdir/hackathon_202406/output/finalTestOGs.aa.fa")

# zero shot score calc.
zeroshot = c()
for(j in 1:length(prob)) { # loop through sequences
  if(j%%100==0) print(j)
  tmp = strsplit(as.character(testfasta[[j]]),"")[[1]]
  refProb = c()
  secProb = c()
  for(i in 1:ifelse(length(tmp)<=999,length(tmp),999)){ # loop through sites (max of 999)
    if(!tmp[i]%in%limma::strsplit2("LAGVSERTIDPKQNFYMHWC","")) next
    refProb = c(refProb,prob[[j]][i,tmp[i]]) # P(obs)
    secProb = c(secProb,prob[[j]][i,which.max(prob[[j]][i,colnames(prob[[j]])!=tmp[i]])]) # P.max(others)
  }
  zeroshot = c(zeroshot,mean(log(refProb/secProb))) # take the mean log ratio of the ps
}
names(zeroshot) = names(prob)
zeroshot_scaled = unlist(tapply(zeroshot,limma::strsplit2(names(zeroshot),":")[,1], function(x) (x-mean(x,na.rm =T))/sd(x,na.rm = T)))
names(zeroshot_scaled) = names(prob)

outTab = data.frame(OG = limma::strsplit2(names(prob),split = ":")[,1],
                    assemblyID = limma::strsplit2(names(prob),split = ":")[,2],
                    zeroShotScore = zeroshot, scaledZeroShotScore = zeroshot_scaled)

write.table(outTab, "/workdir/sh2246/p_phyloGWAS/output/finalTestOGs_ESM2ZeroShot.txt",sep = "\t",quote = F,row.names = T,col.names = T)

omega = read.table("/workdir/sh2246/p_phyloGWAS/output/allOG_dNdS_test.txt")

rownames(omega) = omega[,2]

common_col = intersect(names(zeroshot),omega[,2])

cor(zeroshot[common_col],omega[common_col,15],use = "complete")
cor(zeroshot[common_col],omega[common_col,16],use = "complete")
cor(zeroshot[common_col],omega[common_col,17],use = "complete")
cor(zeroshot[common_col],omega[common_col,18],use = "complete")

cor(zeroshot[common_col],omega[common_col,15],use = "complete",method = "spearman")
cor(zeroshot[common_col],omega[common_col,16],use = "complete",method = "spearman")
cor(zeroshot[common_col],omega[common_col,17],use = "complete",method = "spearman")
cor(zeroshot[common_col],omega[common_col,18],use = "complete",method = "spearman")

# cor(zeroshot_scaled[common_col],omega[common_col,15],use = "complete",method = "spearman")
# cor(zeroshot_scaled[common_col],omega[common_col,16],use = "complete",method = "spearman")
# cor(zeroshot_scaled[common_col],omega[common_col,17],use = "complete",method = "spearman")
# cor(zeroshot_scaled[common_col],omega[common_col,18],use = "complete",method = "spearman")

library(LSD)
par(mfrow = c(1,3))
heatscatter(zeroshot[common_col],omega[common_col,15], xlab = "avg. zero-shot score",ylab = "dS",main = "",ylim = c(0,9))
text(8,8, "rho = -0.01")
heatscatter(zeroshot[common_col],omega[common_col,16], xlab = "avg. zero-shot score",ylab = "dN",main = "")
text(8,0.8, "rho = -0.38")
heatscatter(zeroshot[common_col],omega[common_col,17], xlab = "avg. zero-shot score",ylab = "dN/dS",main = "",ylim = c(0,6))
text(8,5.5, "rho = -0.41")


pheatmap::pheatmap(cor(cbind(zeroshot[common_col],omega[common_col,3:18]),use = "complete",method = "spearman"),cluster_rows = F,cluster_cols = F)


length(unique(limma::strsplit2(names(which(zeroshot[common_col]<2&omega[common_col,16]<0.2)),":")[,1]))


       