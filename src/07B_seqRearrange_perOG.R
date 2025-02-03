# R function to rearrange input fasta per taxa to per gene
# author: Sheng-Kai Hsu
# created: 2024.03.25
# last edited: 2024.04.19
# note: rewrite from the jupyter notebook

rm(list=ls())
library(Biostrings)
library(parallel)
library(limma)
library(rtracklayer)

args = commandArgs(trailingOnly=TRUE)
# arg 1: path to fasta; arg 2: OGID; arg 3: query AA; arg 4: out path

faPerTaxa_flist = list.files(args[1],pattern = "*.fa",full.names = T)
faPerTaxa_list = mclapply(faPerTaxa_flist,function(x) readDNAStringSet(x),mc.cores = 2)

names(faPerTaxa_list)=gsub(".fa","",strsplit2(faPerTaxa_flist,"/")[,ncol(strsplit2(faPerTaxa_flist,"/"))])

# OGID=read.table("/workdir/sh2246/p_evolBNI/output/Pv-672_allTranscriptID.txt",header = F)[,1]

query = readAAStringSet(args[3])
OGID = args[2]

# N = length(faPerTaxa_list)+1

faPerOG = DNAStringSet(lapply(faPerTaxa_list,function(x) {
  tmp=x[[OGID]]
  if (is.null(tmp)) tmp = DNAString("-")
  return(tmp)
}))

dir.create(args[4])
writeXStringSet(faPerOG,paste0(args[4],OGID,".fasta"))

