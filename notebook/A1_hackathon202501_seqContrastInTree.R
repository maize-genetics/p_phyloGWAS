args = commandArgs(trailingOnly=TRUE)

# load libraries
library(ape)
library(phytools)
library(phangorn)
library(testthat)
library(limma)
library(data.table)

#read Labeled tree
tree <- read.tree(args[1])
tree <- read.tree("Dropbox/postDoc/projects/p_phyloGWAS/output/RAxML_Labeled_bestTree.OG0000923_Relax.tree")

#tree naming
treeRenamed = tree
treeRenamed$tip.label = gsub('\\{Foreground\\}|\\{Background\\}','',tree$tip.label)
treeRenamed = makeNodeLabel(treeRenamed)

# read MSA
aa=c("-","*","a","c","d","e","f","g","h","i","k","l","m","n","p","q","r","s","t","v","w","y")
dat=read.phyDat(args[2],format = "fasta",type = "USER",levels= aa)
dat=read.phyDat("Dropbox/postDoc/projects/p_phyloGWAS/output/OG0000923.aa.aligned.fa",format = "fasta",type = "USER",levels= aa)
dat = dat[grep(":0$",names(dat)),]
names(dat) = gsub("\\.|-","_",strsplit2(names(dat),":")[,2])
dat_filled = dat[treeRenamed$tip.label,]
dat_filled = as.data.frame(dat_filled)
dat_filled[,which(!treeRenamed$tip.label%in%names(dat))] = rep("-",1446)
dat_filled = as.phyDat(dat_filled,type = "USER",levels= aa)

#ancestral sequence reconstruction
datPML=pml(treeRenamed,dat_filled)
datASR=ancestral.pml(datPML,type = "ml",return = "phyDat")

# identified foreground node and their corresponding ancestors
foregroundNode = grep("Foreground",c(tree$tip.label,tree$node.label))
foregroundPairsIdx = cbind(Ancestors(tree,foregroundNode,type = "parent"),foregroundNode)
foregroundPairsSeq = cbind(apply(as.character(datASR[foregroundPairsIdx[,1],]),1,paste,collapse = ""),
                           apply(as.character(datASR[foregroundPairsIdx[,2],]),1,paste,collapse = ""))
foregroundPairsSeq = gsub("-","",foregroundPairsSeq)
colnames(foregroundPairsSeq) = c("seq1","seq2")
foregroundPairsSeq[foregroundPairsSeq==""] = NA
foregroundPairsSeq = na.omit(foregroundPairsSeq)

backgroundNode = sample(grep("Background",c(tree$tip.label,tree$node.label)))
if (length(backgroundNode)>length(foregroundNode)) backgroundNode = sample(backgroundNode,length(foregroundNode))
backgroundPairsIdx = cbind(Ancestors(tree,backgroundNode,type = "parent"),backgroundNode)
backgroundPairsSeq = cbind(apply(as.character(datASR[backgroundPairsIdx[,1],]),1,paste,collapse = ""),
                           apply(as.character(datASR[backgroundPairsIdx[,2],]),1,paste,collapse = ""))
backgroundPairsSeq = gsub("-","",backgroundPairsSeq)
colnames(backgroundPairsSeq) = c("seq1","seq2")
backgroundPairsSeq[backgroundPairsSeq==""] = NA
backgroundPairsSeq = na.omit(backgroundPairsSeq)

fwrite(toupper(foregroundPairsSeq),"Dropbox/postDoc/projects/p_phyloGWAS/output/fgPairs.test.tsv",sep = "\t")
fwrite(toupper(backgroundPairsSeq),"Dropbox/postDoc/projects/p_phyloGWAS/output/bgPairs.test.tsv",sep = "\t")

