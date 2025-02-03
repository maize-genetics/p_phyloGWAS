rm(list=ls())
library(limma)

scaf_flist = list.files("/workdir/sh2246/p_phyloGWAS/output/haploidizationRun",pattern = "*goodScafName.txt",recursive = T,full.names = T)
scaf_flist = scaf_flist[grep("PanAnd|NAM|JGI",scaf_flist)]
d_flist = list.files("/workdir/sh2246/p_phyloGWAS/output/haploidizationRun",pattern = "*dist.tab.txt",recursive = T,full.names = T)
d_flist = d_flist[grep("PanAnd|NAM|JGI",d_flist)]

scafList = lapply(scaf_flist,function(x){read.table(x,header = F)})
dList = lapply(d_flist,function(x){read.table(x,header = F)[,1]})


genomeID = strsplit2(scaf_flist,"/")[,7]

cutoff = data.frame(genome = genomeID,cutoff = 0.025)
cutoff[substr(cutoff[,1],1,2)%in%c("Ac","Ir","Rt","Rr","Sm","Te"),2]=0.00 # haploid assemblies or difficult...
cutoff[substr(cutoff[,1],1,2)%in%c("Hp"),2]=0.0045
cutoff[substr(cutoff[,1],1,2)%in%c("Vc"),2]=0.005
cutoff[substr(cutoff[,1],1,2)%in%c("Ab"),2]=0.0075
cutoff[substr(cutoff[,1],1,2)%in%c("Ud","Ss","Sn","Et","Cc"),2]=0.01
cutoff[substr(cutoff[,1],1,2)%in%c("Hc"),2]=0.013
cutoff[substr(cutoff[,1],1,2)%in%c("Bl","Ag"),2]=0.015
cutoff[substr(cutoff[,1],1,2)%in%c("Zn"),2]=0.06
write.table(cutoff,"/workdir/sh2246/p_phyloGWAS/output/haploidizationRun/distCutoff.txt",quote = F,sep = "\t",row.names = F,col.names = F)

for (i in c(1:42)){
  names(dList[[i]]) = rev(scafList[[i]][-1,1])
}

for (i in c(1:42)){
  barplot(dList[[i]],las = 2, ylab = "",ylim = c(0,.2),main = genomeID[i])
  abline(h = cutoff[i,2],col ="red",lwd = 2)
}

for (i in c(1:42)){
  hist(dList[[i]],las = 2, ylab = "",xlim = c(0,.2),main = genomeID[i],breaks = seq(0,1,0.001))
  abline(v = cutoff[i,2],col ="red",lwd = 2)
}

dir.create("/workdir/sh2246/p_phyloGWAS/output/haploidizationRun/distBarPlot/")
for (i in c(1:42)){
  png(paste0("/workdir/sh2246/p_phyloGWAS/output/haploidizationRun/distBarPlot/",i,".png"),
      width = 8.7*1.5,height = 8.7,units = "cm",res = 300,pointsize = 8)
  barplot(dList[[i]],las = 2, ylab = "",ylim = c(0,.2),main = genomeID[i])
  abline(h = cutoff[i,2],col ="red",lwd = 2)
  dev.off()
}

####dot plots####
library(ggplot2)
dir.create("/workdir/sh2246/p_phyloGWAS/output/anchorwave/dotplot/")
dir.create("/workdir/sh2246/p_phyloGWAS/output/haploidizationRun/keepScaf/")
for (i in genomeID[c(1:16,18:42)]){
  print(i)
  tmpIdx = which(genomeID%in%i)
  purged_scaf = names(which(dList[[tmpIdx]]<cutoff[tmpIdx,2]))
  # a=read.table(paste0('/workdir/sh2246/p_phyloGWAS/output/anchorwave/',i,'.anchors'), head=T)
  # a$refChr=factor(a$refChr, levels=1:10)
  # a$queryChr=factor(a$queryChr, levels=rev(unique(a$queryChr[order(a$refChr,a$referenceEnd)])))
  # a_subset = a[a$queryChr%in%scafList[[tmpIdx]][,1],]
  # a_subset$keep = "1"
  # a_subset$keep[a_subset$queryChr%in%purged_scaf] = "2"
  # # whole genome
  # dot.p <- ggplot(a_subset, aes(x=referenceStart, y=queryStart)) + 
  #   geom_point(size=0.5, aes(color=keep)) + 
  #   facet_grid(queryChr~refChr, scales='free', space='free', switch = "y") + 
  #   scale_color_manual(values=c('#0072B7','#F59A23')) + theme_light() +
  #   theme(axis.title.x.bottom = element_text(face = "bold"), axis.title.y.left = element_text(face="bold")) +
  #   theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  #   theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid = element_blank()) +
  #   theme(axis.title=element_text(size=44)) +
  #   theme(strip.text.y.left = element_text(angle = 0)) +
  #   theme(legend.position = "none") +
  #   theme(strip.text.x.top = element_text(size=24)) +
  #   xlab("") +
  #   ylab("")
  # 
  # ggsave(paste0("/workdir/sh2246/p_phyloGWAS/output/anchorwave/dotplot/",i,".png"), 
  #        device = 'png', dpi=200, plot = dot.p, width = 27.5, height = 27.5, units = "in")
  write.table(setdiff(scafList[[tmpIdx]][,1],purged_scaf),
              paste0("/workdir/sh2246/p_phyloGWAS/output/haploidizationRun/keepScaf/",i,".txt"),
              quote = F,sep = "\t",row.names = F,col.names = F)
  #dot.p
}


for (i in genomeID[c(1,6,20,38)]){
  for (j in 1:10){
    tmpIdx = which(genomeID%in%i)
    pruged_scaf = names(which(dList[[tmpIdx]]<.015))
    a=read.table(paste0('/workdir/sh2246/p_phyloGWAS/output/anchorwave/',i,'.anchors'), head=T)
    a$refChr=factor(a$refChr, levels=1:10)
    a$queryChr=factor(a$queryChr, levels=rev(unique(a$queryChr[order(a$refChr,a$referenceEnd)])))
    a_subset = a[a$queryChr%in%scafList[[tmpIdx]][,1]&a$refChr%in%j,]
    a_subset$keep = "1"
    a_subset$keep[a_subset$queryChr%in%pruged_scaf] = "2"
    # whole genome
    dot.p <- ggplot(a_subset, aes(x=referenceStart, y=queryStart)) + 
      geom_point(size=0.1, aes(color=keep)) + 
      facet_grid(queryChr~refChr, scales='free', space='free', switch = "y") + 
      scale_color_manual(values=c('#0072B7','#F59A23')) + theme_light() +
      theme(axis.title.x.bottom = element_text(face = "bold"), axis.title.y.left = element_text(face="bold")) +
      theme(axis.text = element_blank(), axis.ticks = element_blank()) +
      theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid = element_blank()) +
      theme(axis.title=element_text(size=44)) +
      theme(strip.text.y.left = element_text(angle = 0,size = 3)) +
      theme(legend.position = "none") +
      theme(strip.text.x.top = element_text(size=24)) +
      xlab("") +
      ylab("")
    
    print(dot.p)
    }
  
  # ggsave(paste0("/workdir/sh2246/p_phyloGWAS/output/anchorwave/dotplot/",i,".png"), 
  #        device = 'png', dpi=200, plot = dot.p, width = 27.5, height = 17, units = "in")
}  


#### Hc and Hp: special case with undistinguishable distance; how about chr-specific cutoff####
# Hc - keep two copies per chromosome (tetraploid)
i = genomeID[10]
tmpIdx = which(genomeID%in%i)
# purged_scaf = names(which(dList[[tmpIdx]]<cutoff[tmpIdx,2]))
a=read.table(paste0('/workdir/sh2246/p_phyloGWAS/output/anchorwave/',i,'.anchors'), head=T)
a$refChr=factor(a$refChr, levels=1:10)
a$queryChr=factor(a$queryChr, levels=rev(unique(a$queryChr[order(a$refChr,a$referenceEnd)])))
a_subset = a[a$queryChr%in%scafList[[tmpIdx]][,1],]

cutoff2 = c(0.012,.017,.012,.013,.0215,.015,.02,.024,.018,.02)
purged_scaf = c()
for (j in 1:10){
  qChr = as.character(unique(a_subset[a_subset$refChr==j,]$queryChr))
  # print(j)
  dList[[tmpIdx]][qChr]
  purged_scaf = c(purged_scaf,names(which(dList[[tmpIdx]][qChr]<cutoff2[j])))
}
purged_scaf = unique(purged_scaf)
a_subset$keep = "1"
a_subset$keep[a_subset$queryChr%in%purged_scaf] = "2"
# whole genome

dot.p <- ggplot(a_subset, aes(x=referenceStart, y=queryStart)) + 
  geom_point(size=0.5, aes(color=keep)) + 
  facet_grid(queryChr~refChr, scales='free', space='free', switch = "y") + 
  scale_color_manual(values=c('#0072B7','#F59A23')) + theme_light() +
  theme(axis.title.x.bottom = element_text(face = "bold"), axis.title.y.left = element_text(face="bold")) +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid = element_blank()) +
  theme(axis.title=element_text(size=44)) +
  theme(strip.text.y.left = element_text(angle = 0)) +
  theme(legend.position = "none") +
  theme(strip.text.x.top = element_text(size=24)) +
  xlab("") +
  ylab("")
ggsave(paste0("/workdir/sh2246/p_phyloGWAS/output/anchorwave/dotplot/",i,"_v2.png"), 
       device = 'png', dpi=200, plot = dot.p, width = 27.5, height = 27.5, units = "in")
write.table(setdiff(scafList[[tmpIdx]][,1],purged_scaf),
            paste0("/workdir/sh2246/p_phyloGWAS/output/haploidizationRun/keepScaf/",i,"_v2.txt"),
            quote = F,sep = "\t",row.names = F,col.names = F)


#Hp - highly possible a autohexaploid... - keep one copy per chromosome
i = genomeID[11]
tmpIdx = which(genomeID%in%i)
# purged_scaf = names(which(dList[[tmpIdx]]<cutoff[tmpIdx,2]))
a=read.table(paste0('/workdir/sh2246/p_phyloGWAS/output/anchorwave/',i,'.anchors'), head=T)
a$refChr=factor(a$refChr, levels=1:10)
a$queryChr=factor(a$queryChr, levels=rev(unique(a$queryChr[order(a$refChr,a$referenceEnd)])))
a_subset = a[a$queryChr%in%scafList[[tmpIdx]][,1],]

cutoff2 = rep(0.02, 10)
# cutoff2[9] = 0.012

purged_scaf = c()
for (j in 1:10){
  qChr = as.character(unique(a_subset[a_subset$refChr==j,]$queryChr))
  # print(j)
  #barplot(dList[[tmpIdx]][qChr],main = j,ylim = c(0,.02))
  purged_scaf = c(purged_scaf,names(which(dList[[tmpIdx]][qChr]<cutoff2[j])))
}
purged_scaf = unique(purged_scaf)
purged_scaf = purged_scaf[!purged_scaf%in%"scaf_55"] # manually add chr9 right arm back due to potential mis-assembly of scaf_12
a_subset$keep = "1"
a_subset$keep[a_subset$queryChr%in%purged_scaf] = "2"
# whole genome
# a_subset = a_subset[a_subset$refChr==9,]
dot.p <- ggplot(a_subset, aes(x=referenceStart, y=queryStart)) + 
  geom_point(size=0.5, aes(color=keep)) + 
  facet_grid(queryChr~refChr, scales='free', space='free', switch = "y") + 
  scale_color_manual(values=c('#0072B7','#F59A23')) + theme_light() +
  theme(axis.title.x.bottom = element_text(face = "bold"), axis.title.y.left = element_text(face="bold")) +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid = element_blank()) +
  theme(axis.title=element_text(size=44)) +
  theme(strip.text.y.left = element_text(angle = 0)) +
  theme(legend.position = "none") +
  theme(strip.text.x.top = element_text(size=24)) +
  xlab("") +
  ylab("")
ggsave(paste0("/workdir/sh2246/p_phyloGWAS/output/anchorwave/dotplot/",i,"_v2.png"), 
       device = 'png', dpi=200, plot = dot.p, width = 27.5, height = 27.5, units = "in")
write.table(setdiff(scafList[[tmpIdx]][,1],purged_scaf),
            paste0("/workdir/sh2246/p_phyloGWAS/output/haploidizationRun/keepScaf/",i,"_v2.txt"),
            quote = F,sep = "\t",row.names = F,col.names = F)
