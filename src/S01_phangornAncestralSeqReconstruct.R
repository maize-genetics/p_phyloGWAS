# R function to read in orthogroup alignment and reconstruct ancestral sequence
# author: Sheng-Kai Hsu
# date created: 2022.11.23
# date last edited: 2023.02.01
# gap state added

rm(list=ls())
library(phangorn)


args <- commandArgs(TRUE)

if (args[1] == "--help") {
  message("################################ HELP ################################\n
          ## Written by Sheng-Kai Hsu\n## Arguments :\n
          # --input (required)\n\n
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

# load data
dat=read.phyDat(inDir,format = "fasta",type = "USER",levels=c("a", "c", "g", "t", "-") )

# dist estimate & tree construct
datD=dist.ml(dat,model = "F81")
datTree=njs(datD)

#fit ML model
datPML=pml(datTree,dat)

#ancestral sequence reconstruction
datASR=ancestral.pml(datPML,type = "ml")

rootASP=t(datASR[[getRoot(datTree)]]) # root ancestral sequence
rootAS=c("A","C","G","T","-")[apply(rootASP,2,which.max)]

#add Ancestral sequence to MSA
datOut=dat
datOut$MRCA=as.numeric(as.factor(rootAS))

#output
write.phyDat(datOut,outDir,format = "fasta")
