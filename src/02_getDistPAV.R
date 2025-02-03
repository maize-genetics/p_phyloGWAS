# Rscript to calculate distance from outgroups
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
  }else if (args[i] %in% c("--reference")) {
    ID = args[i + 1]
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

d_flist=list.files(inDir,full.names = T)
dList = mclapply(d_flist,function(x) read.table(x,header = T,sep = "\t",row.names = 1),mc.cores = 30)
names(dList)=gsub("_Dist.txt","",strsplit2(d_flist,split = "/")[,ncol(strsplit2(d_flist,split = "/"))])

allTaxa = rownames(dList[[1]])
# 
# dList.fill = mclapply(dList,function(x) {
#   tmp = x
#   colnames(tmp) = rownames(tmp)
#   missTaxa = setdiff(allTaxa,rownames(x))
#   if (length(missTaxa)!=0){
#     tmp = cbind(tmp,matrix(rep(1,nrow(tmp)*length(missTaxa)),nrow = nrow(tmp),ncol = length(missTaxa)))
#     tmp = rbind(tmp,matrix(rep(1,ncol(tmp)*length(missTaxa)),ncol = ncol(tmp),nrow = length(missTaxa)))
#     rownames(tmp)[-(1:nrow(x))] = missTaxa
#     colnames(tmp)[-(1:ncol(x))] = missTaxa 
#   } else tmp = tmp 
#   tmp = tmp[sort(rownames(tmp)),sort(rownames(tmp))]
#   return(tmp)
# },mc.cores = 30)

dist2Ref = sapply(dList,function(x) x[,as.numeric(ID)])
dist2Ref_filtered = dist2Ref[,apply(dist2Ref,2,function(x) sum(is.na(x)))>2]

dist2Ref_filtered[is.na(dist2Ref_filtered)] = 1

brPAV = t(apply(dist2Ref_filtered,2,function(x){
  avg = mean(x[x>0&x<1],na.rm = T)
  sdev = sd(x[x>0&x<1], na.rm = T)
  return(as.numeric(x < avg+1*sdev)) # branch length over x sd away of the avg were considered loss of function
}))

dim(brPAV)
length(allTaxa)
colnames(brPAV) = allTaxa
write.table(brPAV,paste0(outDir,"distPAV.txt"),
            quote = F,col.names = T,row.names = T,sep="\t")
