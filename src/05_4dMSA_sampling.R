#Phylogenetic tree construction for BNI taxa
#author: Sheng-Kai Hsu
#created: 2023.02.28
#last edited: 2024.03.18
#note 2024/03/18: finalize analysis

rm(list=ls())

library(ape)
library(limma)
library(parallel)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
# arg 1: metadata/translateTable; arg 2: path to 4d MSA directory; arg 3: out path

# load metadata
metadata=read.csv(args[1],header = T)
metadata$Sample=gsub(".fasta|.fa","",metadata$genome_path)
metadata$spTaxa=paste0(metadata$Plant_ID,": ",metadata$species.name)

#load 4d msa per transcript
fasta.4d.list=list.files(args[2],".4d.fa$",full.names = T)
fasta.4d.aln=mclapply(fasta.4d.list,function(x) read.FASTA(x),mc.cores = 10)
fasta.4d.alnm=mclapply(fasta.4d.aln,function(x) as.matrix(x),mc.cores = 10)

# select 5k random and filter for sufficient coverage
numTaxa=sapply(fasta.4d.alnm,nrow)
set.seed(123)
idx = sample(which(numTaxa>27),5000)
fasta.4d.alnm3 = fasta.4d.alnm[idx]
fasta.4d.alnm3 = mclapply(fasta.4d.alnm3,function(x) {rownames(x)[grep("Pavag",rownames(x))] = "REF"; return(x)},mc.cores = 30)

# concat MSA for species tree construction
allTaxa = rownames(fasta.4d.alnm3[[which(sapply(fasta.4d.alnm3,nrow)==39)[1]]])
fasta.4d.alnm3 = mclapply(fasta.4d.alnm3,function(x) {
  missTaxa = setdiff(allTaxa,rownames(x))
  tmp = matrix(rep("-",ncol(x)*length(missTaxa)),nrow = length(missTaxa),ncol = ncol(x),byrow = T)
  rownames(tmp) = missTaxa
  tmp = as.matrix(as.DNAbin(tmp))
  out = rbind(x,tmp)
  if (ncol(out) > 20) out = out[allTaxa,sample(1:ncol(out),20)]
  return(out)
},mc.cores = 30)

fasta.4d.aln.merged=do.call(cbind, fasta.4d.alnm3)

# output
write.FASTA(fasta.4d.aln.merged,args[3])
