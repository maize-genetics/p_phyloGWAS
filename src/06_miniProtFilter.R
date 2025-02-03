# R function to filter miniProt alignment result to retain one most "functional" copy per query
# author: Sheng-Kai Hsu
# Date created: 2023.05.02
# Date last edited: 2023.05.02
# input includes  (1) miniProt gff
#                 (2) query aa fasta
#                 (3) unfiltered CDS fasta
#                 (4) reelProt score 
#                 (5) method: p or r; filter for primary or best reelProt score
#                 (6) coverage cutoff
# output includes (1) a info tab. (MPID query start end range coverage  rank  reelProtScore)
#                 (2) filtered fasta
#                 (3) filtered mRNA gff (optional)

#### load package ####
rm(list=ls())
library(rtracklayer)
library(limma)
library(Biostrings)

####read arguments#### 
args <- commandArgs(TRUE)

if (args[1] == "--help") {
  message("################################ HELP ################################\n
          ## Written by Sheng-Kai Hsu\n
          ## Arguments :\n
          ## --gffDir <Dir to miniProt gff3 file> \n
          ## --queryFa <query a.a fasta> \n 
          ## --inFa <unfiltered CDS fasta> \n 
          ## --reelProt <reelProt score for each CDS> \n
          ## --method <p or r; filter for primary or best reelProt score> \n
          ## --cov <query coverage cutoff> \n
          ## --outTab <output directory for info Tab.; optional> \n
          ## --outFa <output directory for filtered fasta> \n
          ## --outGff <output directory for filtered mRNA gffl> \n
          ######################################################################")
  q("no")
}

i = -1
while (i < I(length(args)-1)) {
  i = i + 2
  if (args[i] %in% c("--gffDir")) {
    gffDir = args[i + 1]
  }else if (args[i] %in% c("--queryFa")) {
    queryFaDir = args[i + 1]
  }else if (args[i] %in% c("--inFa")) {
    inFaDir = args[i + 1]
  }else if (args[i] %in% c("--reelProt")) {
    reelProtScore = args[i + 1]
  }else if (args[i] %in% c("--method")) {
    method = args[i + 1]
  }else if (args[i] %in% c("--cov")) {
    coverage = args[i + 1]
  }else if(args[i]%in%c("--outTab")) {
    outTabDir = args[i+1]
  }else if(args[i]%in%c("--outFa")) {
    outFaDir = args[i+1]
  }else if(args[i]%in%c("--outGff")) {
    outGffDir = args[i+1]
  }
  else{
    message("\n### ERROR ###\n\n Argument not recognized:")
    print(args[i])
    message("\n### ERROR ###\n")
    q("no")
  }
}

# gffDir = "output/01_Ab-Traiperm_572-DRAFT-PanAnd-1.0.gff"
# queryFaDir = "output/Pv-672_v3.0_peptide.fa"
# inFaDir = "output/01_Ab-Traiperm_572-DRAFT-PanAnd-1.0.fa"
# reelProtScore = NULL
# method = "p"
# coverage = 0.9
# outTabDir = "output/miniProtFilterTest.txt"
# outFaDir = "output/miniProtFilterTest.fa"
# outGffDir = "output/miniProtFilterTest.gff"
#### read inputs ####
gff=readGFF(gffDir)
queryFa = readAAStringSet(queryFaDir)
inFa = readDNAStringSet(inFaDir)

#### filtering ####
#query length
queryFa_length=width(queryFa)
names(queryFa_length)=names(queryFa)

# gff parsing 
gff.mrna = gff[gff[,3]=="mRNA",]

queryStat = strsplit2(gff.mrna$Target," ")
queryStat = as.data.frame(queryStat)
colnames(queryStat) = c("query","start","end")
queryStat$range = as.numeric(queryStat[,3])-as.numeric(queryStat[,2])+1

# query coverage
qcov=c()
for (i in unique(queryStat[,1])){
  qcov=c(qcov,queryStat[queryStat[,1]==i,4]/queryFa_length[i])
}
queryStat$coverage = qcov

queryStat$rank = gff.mrna$Rank
queryStat$MPID = gff.mrna$ID 
# reelgene score to be added

idx1 = queryStat$coverage > coverage
if (method == "p") idx2 = gff.mrna$Rank==1
# if (method == "r") idx2 = gff.mrna$Rank==1 # update needed with reelGene score

outFa = inFa[queryStat$MPID[idx1&idx2]]
names(outFa) = queryStat[idx1&idx2,1]

gff.mrna.filtered = gff.mrna[idx1&idx2,]

filteredID = gff.mrna.filtered$ID

gff.all.filtered = gff[gff$ID %in% filteredID | as.character(gff$Parent) %in%filteredID,]
#### output ####
write.table(queryStat,outTabDir,col.names = T,row.names = F,quote = F,sep = "\t")
writeXStringSet(outFa,outFaDir)
export.gff3(gff.all.filtered,outGffDir)
