# Rscript to calculate distance matrix based on MSA
# author: Sheng-Kai Hsu
# date created: 2023.09.11
# date last edited: 2023.09.11


rm(list=ls())
library(ape)
library(limma)
library(parallel)

args <- commandArgs(TRUE)

if (args[1] == "--help") {
  message("################################ HELP ################################\n
          ## Written by Sheng-Kai Hsu\n## Arguments :\n
          # --input (required)\n\n
          # --reference (required) \n\n
          #--output path to output, default is current directory (required)\n
          \n######################################################################")
  q("no")
}

outDir = getwd()
i = -1
while (i < I(length(args)-1)) {
  i = i + 2
  if (args[i] %in% c("--input")) {
    inDir = args[i + 1]
  }else if (args[i] %in% c("--output")) {
    outDir = args[i + 1]
  }
  else{
    message("\n### ERROR ###\n\n Argument not recognized:")
    print(args[i])
    message("\n### ERROR ###\n")
    q("no")
  }
}

# dat=read.dna(inDir,format = "fasta")
# dat = dat[sort(rownames(dat)),]
sampleID = gsub(".txt","",strsplit2(inDir,"/")[,ncol(strsplit2(inDir,"/"))])
dat = read.table(inDir,header = T)
refIdx = grepl("B73",colnames(dat))
print(which(refIdx))
if (any(refIdx)){
  d  = dat[,which(refIdx)[1]]
  d = cbind(d, limma::strsplit2(rownames(d),":"))
  d = d[d[,4]=="0",]
  d = d[,-4]
  write.table(d,paste0(outDir,sampleID,"_toB73.txt"),quote = F,sep = "\t")
  
}
