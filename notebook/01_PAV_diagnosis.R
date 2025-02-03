# script for PAV diagnosis
# author: Sheng-Kai Hsu
# date created: 2022.08.03
# date last edited: 2022.08.03

rm(list=ls())

library(limma)
library(scales)
gff=read.delim("~/Dropbox/postDoc/projects/p_phyloGWAS/data/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.mRNA.gff3",header=F)
res=read.csv("~/Dropbox/postDoc/projects/p_phyloGWAS/output/perenniality.GLM.csv",header = T)
PAV=read.delim("~/Dropbox/postDoc/projects/p_phyloGWAS/data/100ASM_minAlign_pt8_maxNM_pt2_pav.txt",header=T,row.names = 1)
rownames(gff)=substr(gff$V9,4,23)

pres=res
pres$chr=as.numeric(strsplit2(res$chr,"chr")[,2])
pres$POS=gff[pres$gene,4]
colnames(pres)[6]="p"
pres$p=-log10(res$perenniality.GLM)
qqman::manhattan(na.omit(pres[pres$chr==9,]),chr = "chr",bp = "POS",p = "p",snp = "gene",logp = F)

pres2=pres
pres2$chr=1
pres2$POS=1:dim(pres2)[1]
qqman::manhattan(na.omit(pres2),chr = "chr",bp = "POS",p = "p",snp = "gene",logp = F)

head(PAV[pres$gene[pres$chr==10],])

png("Dropbox/postDoc/projects/p_phyloGWAS/output/PAV_freqComparison.png",width = 8.7,height = 8.7*1.5,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(3,1))
hist(apply(PAV[pres$gene[pres$chr==1],],1,sum),freq = F,col=alpha("green",0.5),main="chr1",xlab="# of taxa",ylim=c(0,0.035))
hist(apply(PAV[pres$gene[pres$chr==9],],1,sum),freq = F,col=alpha("orange",0.5),main="chr9",xlab="# of taxa",ylim=c(0,0.035))
hist(apply(PAV[pres$gene[pres$chr==10],],1,sum),freq = F,col=alpha("blue",0.5),main="chr10",xlab="# of taxa",ylim=c(0,0.035))
dev.off()
pheatmap::pheatmap(PAV[pres$gene[pres$chr==1][1:1000],],show_rownames = F,show_colnames = F,cluster_rows = F,cluster_cols = F,breaks = c(0,0.5,1),color = c("blue","red"),
                   filename = "Dropbox/postDoc/projects/p_phyloGWAS/output/heatmap_PAV_chr1.png",width = 8.7,height = 8.7,unit="cm",main = "chr1")
pheatmap::pheatmap(PAV[pres$gene[pres$chr==9][1:1000],],show_rownames = F,show_colnames = F,cluster_rows = F,cluster_cols = F,breaks = c(0,0.5,1),color = c("blue","red"),
                   filename = "Dropbox/postDoc/projects/p_phyloGWAS/output/heatmap_PAV_chr9.png",width = 8.7,height = 8.7,unit="cm",main = "chr9")

pheatmap::pheatmap(PAV[pres$gene[pres$chr==9][1001:2000],],show_rownames = F,show_colnames = F,cluster_rows = F,cluster_cols = F,breaks = c(0,0.5,1),color = c("blue","red"),
                   main = "last N")
