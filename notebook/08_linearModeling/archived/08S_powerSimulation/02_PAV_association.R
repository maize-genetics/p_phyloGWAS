# script for PAV association study
# author: Sheng-Kai Hsu
# date created: 2022.08.04
# date last edited: 2023.01.03
# note: single PAV to 10 PAV 

rm(list=ls())
home.dir="~/Dropbox/postDoc/projects/p_phyloGWAS/"
setwd(home.dir)
library(ape)
library(limma)
library(tidyverse)
require(reshape2)
require(plyr)
require(rMVP)      # R version of the MVP package
require(milorGWAS) # mixed logistic regression, see Slack pan_and for reference
source("src/S01_phyloK.R")

####data input and process####
gff_genePos=read.delim("data/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.mRNA.gff3",header = F)
rownames(gff_genePos)=substr(gff_genePos$V9,4,23)
perenniality=read.csv("data/Perennitality_Numeric_Phenotype_by_Taxa.csv",header = T,row.names = 1)
PAV=read.table("data/100ASM_minAlign_pt8_maxNM_pt2_pav.txt",header = T)
inputDist=read.table("data/allRepTranscript4dSites_distMatrix_v2.txt",header = T)
inputTree=read.tree("data/allRepTranscript4dSites_speciesTree.newick")

#get taxa intersect
infTaxa=rownames(inputDist)

#phenotype
perennialityUsed=perenniality[infTaxa,]
perennialityUsed[is.na(perennialityUsed)]=1
perennialityUsed=data.frame(sample=infTaxa,pheno=perennialityUsed)

anno_col=data.frame(life_history=perennialityUsed$pheno)
rownames(anno_col)=perennialityUsed$sample

#genotype
colnames(PAV)=strsplit2(colnames(PAV),".final.contig")[,1]
colnames(PAV)=strsplit2(colnames(PAV),"_megahit")[,1]
colnames(PAV)[grep("PanAnd",colnames(PAV))]=sort(infTaxa[grep("PanAnd",infTaxa)])
PAV$'B73_REF_FULL'=1
PAVUsed=PAV[rownames(PAV)%in%rownames(gff_genePos[grep("chr",gff_genePos$V1),]),infTaxa]

map_panand2 = data.frame(gene=rownames(gff_genePos),chr = gff_genePos[,1],POS=gff_genePos[,4])
map_panand2=map_panand2[grep("chr",map_panand2$chr),]
map_panand2$gene = as.character(map_panand2$gene)
map_panand2 = map_panand2 %>% filter(gene %in% rownames(PAVUsed))
dim(map_panand2)

PAVUsed=PAVUsed[map_panand2$gene,]
PAVfiltered=PAVUsed[apply(PAVUsed,1,sum)>=9,]

map_panand3=map_panand2[map_panand2$gene%in%rownames(PAVfiltered),]

# genoTabLC=data.frame("rs#"=rownames(PAVfiltered),alleles="A/T",
#                      chrom="chr1",pos=1:dim(PAVfiltered)[1],
#                      strand=NA,"assembly#"=NA,center=NA,protLSID=NA,assayLSID=NA,panel=NA,QCcode=NA)
# PAVrecoded=apply(PAVfiltered,c(1,2),function(x) ifelse(x,"A","T"))
# genoTab=as.data.frame(t(cbind(genoTabLC,PAVrecoded)))

# write.table(t(PAVfiltered),"output/PAVTable_filtered.txt",quote = F,row.names = T,col.names = T,sep = "\t")

####kinship calculation####
K1=(1-inputDist)*2
K2=phyloK(inputTree) 
#write.table(K2, "output/kinshipMatrix_SBL.txt",row.names = T,col.names = T,quote = F)
K3=read.delim("output/kinshipMatrix_IBS.txt",header = F,row.names = 1,skip = 3)
colnames(K3)=rownames(K3)
####PAV_simulation####
mutPAVDel=function(PAVtab,perennialTab,transcriptID,freq=c(1,0)){
  out=PAVtab
  idx=replicate(length(transcriptID),c(sample(colnames(PAVtab)[perennialTab==0],size = ceiling(freq[1]*sum(perennialTab==0)),replace = F),
                                       sample(colnames(PAVtab)[perennialTab==1],size = ceiling(freq[2]*sum(perennialTab==1)),replace = F)))
  for (i in 1:length(transcriptID)){
    out[transcriptID,idx[,i]]=0 
  }
  return(out)
} #to do: write this function outside and do unit testing
# mutPAVDel(PAVtab = PAVfiltered,perennialTab = perennialityUsed$pheno,transcriptID = sample(rownames(PAVfiltered)[which(rowSums(PAVfiltered)==90)],10,replace = F),freq=c(0.5,0.01))

GWAS=function(P,G,K,map,nPC=0,method="GLM"){
  # preparing the inputs (PAV, K matrix and map) according to MVP needs
  
  genotype_PAV     = as.big.matrix(G) # dimension: PAV x samples
  genotype_Kmatrix = as.big.matrix(K)           # dimension: samples x samples
  
  # perenniality  # data.frame, dimension: samples x (sample id + traits)
  # map_panand2   # data.frame, dimension: PAVs x 3 (gene, chr and position)
  
  #head(map_panand2)
  #map_panand2$chr = as.numeric(gsub(map_panand2$chr,pattern = 'chr',replacement = '')) # ignore it. does not matter
  
  GWAS_perenniality <- rMVP::MVP(
    # Data Inputs 
    phe  = P, 
    geno = genotype_PAV,
    K    = genotype_Kmatrix,
    map  = map,
    
    # Method 
    method=method, 
    
    # Number of Principal components for each mode ---
    nPC.GLM = nPC, nPC.MLM = nPC,
    
    # Threshold 
    permutation.threshold = F, threshold = 0.05,
    
    # Outputs and computing parameters
    file.output=F, ncpus=8)
}

set.seed(1000)
resList0=list()
resList1=list()
resList2=list()
resList3=list()
mut=list()
for (freq in as.character(c(100,90,80,70,60,50,40,30,20,10,1))){
  mutTargets=replicate(100,sample(rownames(PAVfiltered)[which(rowSums(PAVfiltered)==90)],10,replace = F))
  mut[[freq]]=mutTargets
#  rownames(PAVUsed)[mutTargets]
  
  PAVMut=apply(mutTargets,2,function(x) {
    mutPAVDel(PAVtab = PAVfiltered,perennialTab = perennialityUsed$pheno,transcriptID = x,freq=c(as.numeric(freq)/100,0.01))
  })
  resList0[[freq]]=lapply(PAVMut,function(x) GWAS(perennialityUsed,x,K1,map_panand3))
  resList1[[freq]]=lapply(PAVMut,function(x) GWAS(perennialityUsed,x,K1,map_panand3,nPC = 2,method = c("GLM","MLM")))
  resList2[[freq]]=lapply(PAVMut,function(x) GWAS(perennialityUsed,x,K2,map_panand3,nPC = 2,method = c("GLM","MLM")))
  resList3[[freq]]=lapply(PAVMut,function(x) GWAS(perennialityUsed,x,K3,map_panand3,nPC = 2,method = c("GLM","MLM")))
  
}
names(mut)=c(100,90,80,70,60,50,40,30,20,10,1)

saveRDS(resList0,"output/resList0_v2.rds")
saveRDS(resList1,"output/resList1_v2.rds")
saveRDS(resList2,"output/resList2_v2.rds")
saveRDS(resList3,"output/resList3_v2.rds")

blRes0=GWAS(perennialityUsed,PAVfiltered,K1,map_panand3)
blRes1=GWAS(perennialityUsed,PAVfiltered,K1,map_panand3,nPC = 2,method = c("GLM","MLM"))
blRes2=GWAS(perennialityUsed,PAVfiltered,K2,map_panand3,nPC = 2,method = c("GLM","MLM"))
blRes3=GWAS(perennialityUsed,PAVfiltered,K3,map_panand3,nPC = 2,method = c("GLM","MLM"))

dpTab0=c()
dpTab1=c()
dpTab2=c()
dpTab3=c()
sTab0=c()
sTab1=c()
sTab2=c()
sTab3=c()
for (j in as.character(c(100,90,80,70,60,50,40,30,20,10,1))){
  dP0=c()
  dP1=c()
  dP2=c()
  dP3=c()
  ap0=c()
  ap1=c()
  ap2=c()
  ap3=c()
  for(i in 1:ncol(mut[[j]])){
    dP0=rbind(dP0,sum(p.adjust(resList0[[j]][[i]][[2]][,3],"bonferroni")[rownames(PAVfiltered)%in%mut[[j]][,i]]<0.05,na.rm=T))
    dP1=rbind(dP1,colSums(sapply(resList1[[j]][[i]][-c(1,4)],function(x) p.adjust(x[,3],"bonferroni")[rownames(PAVfiltered)%in%mut[[j]][,i]]<0.05),na.rm=T))
    dP2=rbind(dP2,colSums(sapply(resList2[[j]][[i]][-c(1,4)],function(x) p.adjust(x[,3],"bonferroni")[rownames(PAVfiltered)%in%mut[[j]][,i]]<0.05),na.rm=T))
    dP3=rbind(dP3,colSums(sapply(resList3[[j]][[i]][-c(1,4)],function(x) p.adjust(x[,3],"bonferroni")[rownames(PAVfiltered)%in%mut[[j]][,i]]<0.05),na.rm=T))
    ap0=rbind(ap0,sum(p.adjust(resList0[[j]][[i]][[2]][,3],"bonferroni")<0.05,na.rm=T) )
    ap1=rbind(ap1,sapply(c("glm.results","mlm.results"),function(x) sum(p.adjust(resList1[[j]][[i]][[x]][,3],"bonferroni")<0.05,na.rm=T) ))
    ap2=rbind(ap2,sapply(c("glm.results","mlm.results"),function(x) sum(p.adjust(resList2[[j]][[i]][[x]][,3],"bonferroni")<0.05,na.rm=T) ))
    ap3=rbind(ap3,sapply(c("glm.results","mlm.results"),function(x) sum(p.adjust(resList3[[j]][[i]][[x]][,3],"bonferroni")<0.05,na.rm=T) ))
  }
  dpTab0=rbind(dpTab0,mean(dP0/10,na.rm=T))
  dpTab1=rbind(dpTab1,apply(dP1/10,2,mean,na.rm = T))
  dpTab2=rbind(dpTab2,apply(dP2/10,2,mean,na.rm = T))
  dpTab3=rbind(dpTab3,apply(dP3/10,2,mean,na.rm = T))
  sTab0=rbind(sTab0,sum(dP0/ap0,na.rm=T)/nrow(ap0))
  sTab1=rbind(sTab1,colSums(dP1/ap1,na.rm=T)/nrow(ap1))
  sTab2=rbind(sTab2,colSums(dP2/ap2,na.rm=T)/nrow(ap2))
  sTab3=rbind(sTab3,colSums(dP3/ap3,na.rm=T)/nrow(ap3))
}

rownames(dpTab0)=c(100,90,80,70,60,50,40,30,20,10,1)
rownames(dpTab1)=c(100,90,80,70,60,50,40,30,20,10,1)
rownames(dpTab2)=c(100,90,80,70,60,50,40,30,20,10,1)
rownames(dpTab3)=c(100,90,80,70,60,50,40,30,20,10,1)
rownames(sTab0)=c(100,90,80,70,60,50,40,30,20,10,1)
rownames(sTab1)=c(100,90,80,70,60,50,40,30,20,10,1)
rownames(sTab2)=c(100,90,80,70,60,50,40,30,20,10,1)
rownames(sTab3)=c(100,90,80,70,60,50,40,30,20,10,1)



png("output/powerSimulation_10PAV.png",width = 8.7*1.5, height = 8.7*1.5,units = "cm",pointsize = 8,res = 600)
par(mfrow=c(3,1))
barplot(cbind(dpTab0,dpTab1),beside = T,xlab = "Frequency of deletion in annual taxa (%)",ylab = "Power",legend.text = T,
        names.arg = c("NM","GLM","MLM"),main = "K (Distance)")
barplot(cbind(dpTab0,dpTab2),beside = T,xlab = "Frequency of deletion in annual taxa (%)",ylab = "Power",
        names.arg = c("NM","GLM","MLM"),main = "K (Shared Branch Length)")
barplot(cbind(dpTab0,dpTab3),beside = T,xlab = "Frequency of deletion in annual taxa (%)",ylab = "Power",
        names.arg = c("NM","GLM","MLM"),main = "K (IBS)")
dev.off()

png("output/powerSimulation_10PAV_sensitivity.png",width = 8.7*1.5, height = 8.7*1.5,units = "cm",pointsize = 8,res = 600)
par(mfrow=c(3,1))
barplot(cbind(sTab0,sTab1),beside = T,xlab = "Frequency of deletion in annual taxa (%)",ylab = "Sensitivity",legend.text = T,
        names.arg = c("NM","GLM","MLM"),main = "K (Distance)",ylim=c(0,1))
barplot(cbind(sTab0,sTab2),beside = T,xlab = "Frequency of deletion in annual taxa (%)",ylab = "Sensitivity",
        names.arg = c("NM","GLM","MLM"),main = "K (Shared Branch Length)",ylim=c(0,1))
barplot(cbind(sTab0,sTab3),beside = T,xlab = "Frequency of deletion in annual taxa (%)",ylab = "Sensitivity",
        names.arg = c("NM","GLM","MLM"),main = "K (IBS)",ylim=c(0,1))
dev.off()

png("output/powerSimulation_10PAV_qqplot_freq100.png",width = 8.7*1.5, height = 8.7*1.5,units = "cm",pointsize = 8,res = 600)
par(mfrow=c(3,3))
qqman::qq(resList0$`100`[[1]]$glm.results[,3],col=ifelse(rownames(PAVfiltered[!is.na(resList0$`100`[[1]]$glm.results[,3]),])[order(-log(na.omit(resList0$`100`[[1]]$glm.results[,3])),decreasing = T)]%in%mut[["100"]][,1],"red","black"))
qqman::qq(resList1$`100`[[1]]$glm.results[,3],col=ifelse(rownames(PAVfiltered[!is.na(resList1$`100`[[1]]$glm.results[,3]),])[order(-log(na.omit(resList1$`100`[[1]]$glm.results[,3])),decreasing = T)]%in%mut[["100"]][,1],"red","black"))
qqman::qq(resList1$`100`[[1]]$mlm.results[,3],col=ifelse(rownames(PAVfiltered[!is.na(resList1$`100`[[1]]$mlm.results[,3]),])[order(-log(na.omit(resList1$`100`[[1]]$mlm.results[,3])),decreasing = T)]%in%mut[["100"]][,1],"red","black"))

qqman::qq(resList0$`100`[[1]]$glm.results[,3],col=ifelse(rownames(PAVfiltered[!is.na(resList0$`100`[[1]]$glm.results[,3]),])[order(-log(na.omit(resList0$`100`[[1]]$glm.results[,3])),decreasing = T)]%in%mut[["100"]][,1],"red","black"))
qqman::qq(resList2$`100`[[1]]$glm.results[,3],col=ifelse(rownames(PAVfiltered[!is.na(resList2$`100`[[1]]$glm.results[,3]),])[order(-log(na.omit(resList2$`100`[[1]]$glm.results[,3])),decreasing = T)]%in%mut[["100"]][,1],"red","black"))
qqman::qq(resList2$`100`[[1]]$mlm.results[,3],col=ifelse(rownames(PAVfiltered[!is.na(resList2$`100`[[1]]$mlm.results[,3]),])[order(-log(na.omit(resList2$`100`[[1]]$mlm.results[,3])),decreasing = T)]%in%mut[["100"]][,1],"red","black"))

qqman::qq(resList0$`100`[[1]]$glm.results[,3],col=ifelse(rownames(PAVfiltered[!is.na(resList0$`100`[[1]]$glm.results[,3]),])[order(-log(na.omit(resList0$`100`[[1]]$glm.results[,3])),decreasing = T)]%in%mut[["100"]][,1],"red","black"))
qqman::qq(resList3$`100`[[1]]$glm.results[,3],col=ifelse(rownames(PAVfiltered[!is.na(resList3$`100`[[1]]$glm.results[,3]),])[order(-log(na.omit(resList3$`100`[[1]]$glm.results[,3])),decreasing = T)]%in%mut[["100"]][,1],"red","black"))
qqman::qq(resList3$`100`[[1]]$mlm.results[,3],col=ifelse(rownames(PAVfiltered[!is.na(resList3$`100`[[1]]$mlm.results[,3]),])[order(-log(na.omit(resList3$`100`[[1]]$mlm.results[,3])),decreasing = T)]%in%mut[["100"]][,1],"red","black"))
dev.off()


png("output/powerSimulation_10PAV_qqplot_freq10.png",width = 8.7*1.5, height = 8.7*1.5,units = "cm",pointsize = 8,res = 600)
par(mfrow=c(3,3))
qqman::qq(resList0$`10`[[1]]$glm.results[,3],col=ifelse(rownames(PAVfiltered[!is.na(resList0$`10`[[1]]$glm.results[,3]),])[order(-log(na.omit(resList0$`10`[[1]]$glm.results[,3])),decreasing = T)]%in%mut[["10"]][,1],"red","black"))
qqman::qq(resList1$`10`[[1]]$glm.results[,3],col=ifelse(rownames(PAVfiltered[!is.na(resList1$`10`[[1]]$glm.results[,3]),])[order(-log(na.omit(resList1$`10`[[1]]$glm.results[,3])),decreasing = T)]%in%mut[["10"]][,1],"red","black"))
qqman::qq(resList1$`10`[[1]]$mlm.results[,3],col=ifelse(rownames(PAVfiltered[!is.na(resList1$`10`[[1]]$mlm.results[,3]),])[order(-log(na.omit(resList1$`10`[[1]]$mlm.results[,3])),decreasing = T)]%in%mut[["10"]][,1],"red","black"))

qqman::qq(resList0$`10`[[1]]$glm.results[,3],col=ifelse(rownames(PAVfiltered[!is.na(resList0$`10`[[1]]$glm.results[,3]),])[order(-log(na.omit(resList0$`10`[[1]]$glm.results[,3])),decreasing = T)]%in%mut[["10"]][,1],"red","black"))
qqman::qq(resList2$`10`[[1]]$glm.results[,3],col=ifelse(rownames(PAVfiltered[!is.na(resList2$`10`[[1]]$glm.results[,3]),])[order(-log(na.omit(resList2$`10`[[1]]$glm.results[,3])),decreasing = T)]%in%mut[["10"]][,1],"red","black"))
qqman::qq(resList2$`10`[[1]]$mlm.results[,3],col=ifelse(rownames(PAVfiltered[!is.na(resList2$`10`[[1]]$mlm.results[,3]),])[order(-log(na.omit(resList2$`10`[[1]]$mlm.results[,3])),decreasing = T)]%in%mut[["10"]][,1],"red","black"))

qqman::qq(resList0$`10`[[1]]$glm.results[,3],col=ifelse(rownames(PAVfiltered[!is.na(resList0$`10`[[1]]$glm.results[,3]),])[order(-log(na.omit(resList0$`10`[[1]]$glm.results[,3])),decreasing = T)]%in%mut[["10"]][,1],"red","black"))
qqman::qq(resList3$`10`[[1]]$glm.results[,3],col=ifelse(rownames(PAVfiltered[!is.na(resList3$`10`[[1]]$glm.results[,3]),])[order(-log(na.omit(resList3$`10`[[1]]$glm.results[,3])),decreasing = T)]%in%mut[["10"]][,1],"red","black"))
qqman::qq(resList3$`10`[[1]]$mlm.results[,3],col=ifelse(rownames(PAVfiltered[!is.na(resList3$`10`[[1]]$mlm.results[,3]),])[order(-log(na.omit(resList3$`10`[[1]]$mlm.results[,3])),decreasing = T)]%in%mut[["10"]][,1],"red","black"))
dev.off()

####association using rMVP####
path = paste0("output/test_",2)
dir.create(path)
setwd(path)

# require(AGHmatrix) I did this test to check if the issue was in the Kmatrix. It seems it is not.
#K_matrix92 = Gmatrix(SNPmatrix = data.matrix(PAV_matrix_noscaf),ploidy = 2)

# preparing the inputs (PAV, K matrix and map) according to MVP needs

genotype_PAV     = as.big.matrix(PAVfiltered) # dimension: PAV x samples
genotype_Kmatrix = as.big.matrix(K2)           # dimension: samples x samples

# perenniality  # data.frame, dimension: samples x (sample id + traits)
# map_panand2   # data.frame, dimension: PAVs x 3 (gene, chr and position)

#head(map_panand2)
#map_panand2$chr = as.numeric(gsub(map_panand2$chr,pattern = 'chr',replacement = '')) # ignore it. does not matter

GWAS_perenniality <- rMVP::MVP(
  # Data Inputs 
  phe  = perennialityUsed, 
  geno = genotype_PAV,
  K    = genotype_Kmatrix,
  map  = map_panand3,
  
  # Method 
  method=c("GLM", "MLM", "FarmCPU"), 
  
  # Bin parameters 
  # method.bin = 'FaST-LMM', 
  #  bin.selection = seq(1, 100, 1),
  
  # Number of Principal components for each mode ---
  nPC.FarmCPU = 2,nPC.GLM = 2,nPC.MLM = 2,
  
  # Threshold 
  permutation.threshold = T, threshold = 0.05, maxLoop=3,
  
  # Outputs and computing parameters
  file.output=T, ncpus=1)

GWAS_perenniality$farmcpu.results
pheatmap(PAVfiltered[which(GWAS_perenniality$farmcpu.results[,3]<1e-7),],cluster_rows = F,cluster_cols = T,annotation_col = anno_col)


####rTASSEL####
pheno=readPhenotypeFromDataFrame(phenotypeDF = perennialityUsed,taxaID = "taxa")
geno=readGenotypeTableFromPath("output/PAVTable_filtered.hmp.txt")



a=1:5
pm=matrix(c(0,2,1,1,1,
            2,0,1,1,1,
            1,1,0,1,1,
            1,1,1,0,1,
            1,1,1,1,0),5,5)
pm=pm/sum(pm)
sample(a,2,prob = pm)
