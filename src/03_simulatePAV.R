# script to simulate trait-associated PAV at different frequency
# author: Sheng-Kai Hsu
# date created: 2023.09.13
# date last edited: 2023.09.13


rm(list=ls())
library(ape)
library(limma)

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
  if (args[i] %in% c("--treeDir")) {
  inDir = args[i + 1]
  }else if (args[i] %in% c("--pMatrix")) {
    pDir = args[i + 1]
  }else if (args[i] %in% c("--traitID")) {
    traitID = args[i + 1]
  }else if (args[i] %in% c("--freq")) {
    freq = args[i + 1]
  }else if (args[i] %in% c("--replicate")) {
    r = args[i + 1]
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
# define functions
branchSampling=function(tree,N,targetID){
  tre=makeNodeLabel(tree,prefix = "")
  allsubTreeLabel=lapply(tre$node.label,function(x) extract.clade(tre,node = x)$tip.label)
  targetNode=tre$node.label[sapply(allsubTreeLabel,function(x) all(x%in%targetID))]
  targetNodeBL=tre$edge.length[tre$edge[,2]%in%(as.numeric(targetNode)+length(tre$tip.label))]
  targetTip=which(tre$tip.label%in%targetID)
  targetTipBL=tre$edge.length[tre$edge[,2]%in%targetTip]
  target=c(targetNode,tre$tip.label[targetTip])
  targetBL=c(targetNodeBL,targetTipBL)
  D=c(rpois(5e5,length(target)*.2),rpois(1e5,length(target)))
  d=sample(1+D,N)
  idx=c()
  for (j in 1:length(d)){
    mutEdge=unique(sample(target,d[j],prob = targetBL/sum(targetBL),replace = T))
    # identify the keep substree below the selected node (extract.clade() or keep.tip())
    mutTaxa=c()
    for (i in 1:length(mutEdge)){
      if (mutEdge[i]%in%tre$tip.label[targetTip])  mutTaxa=c(mutTaxa,mutEdge[i])
      else mutTaxa=c(mutTaxa,extract.clade(tre,node = mutEdge[i])$tip.label)
    }
    idx=rbind(idx,as.numeric(tre$tip.label[targetTip]%in%mutTaxa))
  }
  return(idx)
}
mutPAVDelCoalescent=function(p,n,targetID,freq,simGeno){
  out=matrix(1,nrow = n,ncol = nrow(p))
  colnames(out) = rownames(p)
  cap=ceiling(freq*length(targetID))
  # print(cap)
  samplingPool=simGeno[rowSums(simGeno!=0)==cap,]
  idx=list()
  for (j in 1:n){
    idx[[j]]=colnames(simGeno)[samplingPool[sample(1:nrow(samplingPool),1),]!=0]
  }
  for (i in 1:n){
    out[i,idx[[i]]]=0 
  }
  return(out)
}

# load data
tree = read.tree(inDir)
pMatrix = read.csv(pDir,header = T,row.names = 1)  
# tree = read.tree("p_phyloGWAS/data/speciesTrees_finaltaxa_median4dDz_tabasco4000.newick")  
# pMatrix = read.csv("p_phyloGWAS/data/ePC_all_data_310Taxa_09.12.23.csv",header = T,row.names = 1)  
# lambda = 2
# N = 1e5  
# n = 10
freq = as.numeric(freq)
r =as.numeric(r)
# filter
pMatrix.filtered = pMatrix[rownames(pMatrix)%in%tree$tip.label,]
targetID = rownames(pMatrix.filtered[pMatrix.filtered[,traitID]<quantile(pMatrix.filtered[,traitID],0.1),])

mutGeno=branchSampling(tree,1e4,targetID)
colnames(mutGeno)=tree$tip.label[tree$tip.label%in%targetID]

for(i in 1:r) {
  out = mutPAVDelCoalescent(pMatrix.filtered,10,targetID,freq,simGeno = mutGeno)
  rownames(out) = paste0("simPAV",0:9)
  
  #output
  write.table(out,paste0(outDir,"_",i,".txt"),quote = F,sep = "\t",row.names = T,col.names = T)
}

