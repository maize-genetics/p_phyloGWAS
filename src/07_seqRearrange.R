# R function to rearrange input fasta per taxa to per gene
# author: Sheng-Kai Hsu
# created: 2024.03.25
# last edited: 2024.03.25
# note: rewrite from the jupyter notebook

rm(list=ls())
library(Biostrings)
library(parallel)
library(limma)
library(rtracklayer)

args = commandArgs(trailingOnly=TRUE)
# arg 1: path to fasta; arg 2: threads; arg 3: query AA; arg 4: out path

faPerTaxa_flist = list.files(args[1],pattern = "*.fa",full.names = T)
faPerTaxa_list = mclapply(faPerTaxa_flist,function(x) readDNAStringSet(x),mc.cores = args[2])

names(faPerTaxa_list)=gsub(".fa","",strsplit2(faPerTaxa_flist,"/")[,ncol(strsplit2(faPerTaxa_flist,"/"))])

# OGID=read.table("/workdir/sh2246/p_evolBNI/output/Pv-672_allTranscriptID.txt",header = F)[,1]

query = readAAStringSet(args[3])
OGID = names(query)

# N = length(faPerTaxa_list)+1

faPerOG_list = list()
for (i in OGID){
  faPerOG_list[[i]] = DNAStringSet(lapply(faPerTaxa_list,function(x) {
    tmp=x[[i]]
    if (is.null(tmp)) tmp = DNAString("-")
    return(tmp)
  }))
  # faPerOG_list[[i]][[N]]=pvCDS[[i]]
  # names(faPerOG_list[[i]])[N] = i
}

dir.create(args[4])
for (i in OGID){
  writeXStringSet(faPerOG_list[[i]],paste0(args[4],i,".fasta"))
}
