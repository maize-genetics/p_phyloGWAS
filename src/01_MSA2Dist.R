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

dat=read.dna(inDir,format = "fasta")
sampleID = gsub(".fa","",strsplit2(inDir,"/")[,ncol(strsplit2(inDir,"/"))])
# dat = dat[sort(rownames(dat)),]
d = dist.dna(dat,pairwise.deletion = T,as.matrix = T,model = "K81",gamma = T)
write.table(d,paste0(outDir,sampleID,"_Dist.txt"),quote = F,sep = "\t")
