rm(list = ls())
# Use reticulate to import numpy
Sys.setenv(RETICULATE_PYTHON = "/usr/bin/python")
library(reticulate)
library(Biostrings)

py_config()
np <- import("numpy")

args = commandArgs(trailingOnly=TRUE)

inputLogit = args[1]
inputFASTA = args[2]
ouutput = args[3]

# load numpy
npz <- np$load(inputLogit)
prob <- lapply(npz$files, function(x) npz[[x]])
names(prob) <- npz$files

# load fasta
testfasta = readAAStringSet(inputFASTA)

# zero shot score calc.
zeroshot = c()
for(j in 1:length(prob)) { # loop through sequences
  if(j%%100==0) print(j)
  tmp = strsplit(as.character(testfasta[[j]]),"")[[1]]
  refProb = c()
  secProb = c()
  for(i in 1:length(tmp)){ # loop through sites (max of 999) fixed by Jingjing
    if(!tmp[i]%in%limma::strsplit2("LAGVSERTIDPKQNFYMHWC","")) next
    refProb = c(refProb,prob[[j]][i,tmp[i]]) # P(obs)
    secProb = c(secProb,prob[[j]][i,which.max(prob[[j]][i,colnames(prob[[j]])!=tmp[i]])]) # P.max(others)
  }
  zeroshot = c(zeroshot,mean(log(refProb/secProb))) # take the mean log ratio of the ps
}
names(zeroshot) = names(prob)
zeroshot_scaled = unlist(tapply(zeroshot,limma::strsplit2(names(zeroshot),":")[,1], function(x) (x-mean(x,na.rm =T))/sd(x,na.rm = T)))
names(zeroshot_scaled) = names(prob)

outTab = data.frame(seqID = names(prob),zeroShotScore = zeroshot, scaledZeroShotScore = zeroshot_scaled)

write.table(outTab, ouutput,sep = "\t",quote = F,row.names = T,col.names = T)
