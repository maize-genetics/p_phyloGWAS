# Compute envPCs (PCA over per-species environmental-feature quantiles) from the raw
# per-occurrence environmental data pulled by src/09_pulling_envData.r, and write out the
# derived pipeline outputs consumed downstream (envData_707Poaceae*, HyPhy adaptation lists,
# loading diagnostics).
# author: Sheng-Kai Hsu
# (ported from notebook/06_envirotyping/06B_visualizationEnvAdapt.ipynb, cells 0-43,
#  as part of the 06A/06B restructuring; invoked from 06A_spCoordEnvData.sh)

library(ape)
library(magrittr)
library(reshape2)
library(plyr)
library(pheatmap)

PHYLOGWAS_ROOT <- Sys.getenv("PHYLOGWAS_ROOT", unset = "/workdir/sh2246/p_phyloGWAS")

source(file.path(PHYLOGWAS_ROOT, "src/process_synthetic_fun.R"))

#load metadata (assemblyID to species name)
metadata = read.delim(file.path(PHYLOGWAS_ROOT, "data/Poaceae_metadata_filtered_2025.08.28.tsv"),header = T)
metadata$spTaxa = paste(metadata[,1],metadata[,3],sep = ":")
colnames(metadata)[3] = "latest_name"
rownames(metadata) = metadata$assemblyID

# load tree
spTre = read.tree(file.path(PHYLOGWAS_ROOT, "output/PoaceaeTree_angiosperm353_astral_filtered_20250407.nwk"))
spTre.rooted = root(spTre,"ASM1935983v1",resolve.root = T)
spTre.rooted$edge.length[is.na(spTre.rooted$edge.length)] = 0.1

# load raw env. data
envTrait_range = read.delim(file.path(PHYLOGWAS_ROOT, "output/metadataFormalOut/formal_envData_20240820.txt"),na.strings = "NA")

envTrait_range_merged = merge(metadata[,c(1,3)],envTrait_range,by = "latest_name")

# per-species dominant Koppen climate class (KG3), needed by 06B's KG3/tree-overlay sections -
# derived here since it's the only downstream consumer of the raw per-occurrence env data
KG3_class = tapply(envTrait_range_merged$KG3,envTrait_range_merged$assemblyID,function(x) as.numeric(names(table(x))[which.max(table(x))]))
write.table(data.frame(assemblyID = names(KG3_class), KG3_class = KG3_class),
            file.path(PHYLOGWAS_ROOT, "output/KG3_perSpecies_20250804.txt"),
            quote = F, row.names = F, sep = "\t")

exampleDat = envTrait_range_merged$bio01_Annual_Mean_Temperature[envTrait_range_merged$latest_name%in%"Andropogon gerardii"]
png(file.path(PHYLOGWAS_ROOT, "output/figure/suppFig/FigS1_envPCPipeline_a2.png"),
    width = 4,height = 4,units = "cm",res = 600, pointsize = 6)
par(mar = c(3,3,2,2))
plot(density(exampleDat,na.rm = T),main = expression(italic("Andropogon gerardi")),
     xlab = "Bio01 Annual Mean Temperature",mgp = c(2,.7,0),col = "grey30", ylim = c(0,.12))
abline(v = quantile(exampleDat,c(.1,.5,.9),na.rm = T),col = c("blue","black","red"))
text(quantile(exampleDat,c(.1,.5,.9)),0.11,labels =  paste0(c(10,50,90),"%"),
     col = c("blue","black","red"),pos = 4,offset = .1,cex = .8)
dev.off()

# for a few features in GSDE: 156 ==NA
problemFeature = c()
for(i in grep("GSDE",colnames(envTrait_range_merged))){
    if(max(envTrait_range_merged[,i],na.rm = T)==156) problemFeature = c(problemFeature,colnames(envTrait_range_merged)[i])
}

envTrait_range_merged[,problemFeature][envTrait_range_merged[,problemFeature]==156] = NA

computing_quantitave_features <-
envTrait_range_merged %>%  # combined coordinates from BIEN and GBIF
  reshape2:::melt(
    id.vars=c("latest_name","assemblyID",
              "decimalLatitude","decimalLongitude","Elevation_m")) %>%
  plyr::ddply(.(assemblyID , variable), plyr::summarise,
              quan10 = quantile(value,0.10, na.rm=T),            # bottom 10%
              quan50 = quantile(value,0.50, na.rm=T),            # median
              quan90 = quantile(value,0.90, na.rm=T)) %>%        # top 90%
  reshape2:::melt(variable.name='quantile',value.name = 'q_value') %>%
  reshape2::acast(assemblyID ~ variable+quantile,value.var = 'q_value')

ePCs <-
  computing_quantitave_features %>%
  process_synthetic(n.synthetic = 40)

saveRDS(ePCs,file.path(PHYLOGWAS_ROOT, 'output/ePC_20250804.rds'))

png(file.path(PHYLOGWAS_ROOT, "output/figure/suppFig/FigS1_envPCPipeline_a3.png"),
    width = 4,height = 4,unit = "cm",pointsize = 6,res = 600)
par(mar = c(3,3,1,1),mgp = c(2,.7,0))
plot(ePCs$variance.explained.by.eigen[1:10,2],xlab = 'envPCs',ylab = 'Variance explained (%)')
dev.off()

commonID = intersect(spTre.rooted$tip.label,rownames(ePCs$environmental.features))
commonID = intersect(commonID,metadata$assemblyID)
eTraits_filtered = as.data.frame(ePCs$environmental.features[commonID,])
spTre.rooted.filtered = keep.tip(spTre.rooted,commonID)

write.tree(spTre.rooted.filtered,file.path(PHYLOGWAS_ROOT, "output/PoaceaeTree_angiosperm353_astral_filtered_withEnvData_20250804.nwk"))

write.table(eTraits_filtered,file.path(PHYLOGWAS_ROOT, "output/envData_707Poaceae_20250804.txt"),
            quote = F, col.names = T, row.names = T,sep = "\t")

eTraits_filtered_out = cbind(species_name = metadata[commonID,]$latest_name,eTraits_filtered)

write.table(eTraits_filtered_out,file.path(PHYLOGWAS_ROOT, "output/envData_707Poaceae_withSpName_20250804.txt"),
            quote = F, col.names = T, row.names = T,sep = "\t")

### write out assembly ID with peripheral envPC values for HyPhy analysis
envPC1_coldTolSp = rownames(eTraits_filtered)[eTraits_filtered$envPC_1<quantile(eTraits_filtered$envPC_1,0.3)]
envPC1_warmTolSp = rownames(eTraits_filtered)[eTraits_filtered$envPC_1>quantile(eTraits_filtered$envPC_1,0.7)]
write.table(envPC1_coldTolSp,file.path(PHYLOGWAS_ROOT, "output/coldAdaptedAssemblies.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
write.table(envPC1_warmTolSp,file.path(PHYLOGWAS_ROOT, "output/warmAdaptedAssemblies.txt"),quote = F,sep = "\t",row.names = F,col.names = F)

envPC2_droughtTolSp = rownames(eTraits_filtered)[eTraits_filtered$envPC_2<quantile(eTraits_filtered$envPC_2,0.3)]
envPC2_wetTolSp = rownames(eTraits_filtered)[eTraits_filtered$envPC_2>quantile(eTraits_filtered$envPC_2,0.7)]
write.table(envPC2_droughtTolSp,file.path(PHYLOGWAS_ROOT, "output/droughtAdaptedAssemblies.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
write.table(envPC2_wetTolSp,file.path(PHYLOGWAS_ROOT, "output/wetAdaptedAssemblies.txt"),quote = F,sep = "\t",row.names = F,col.names = F)

envPC3_sandTolSp = rownames(eTraits_filtered)[eTraits_filtered$envPC_3<quantile(eTraits_filtered$envPC_3,0.3)]
envPC3_clayTolSp = rownames(eTraits_filtered)[eTraits_filtered$envPC_3>quantile(eTraits_filtered$envPC_3,0.7)]
write.table(envPC3_sandTolSp,file.path(PHYLOGWAS_ROOT, "output/sandAdaptedAssemblies.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
write.table(envPC3_clayTolSp,file.path(PHYLOGWAS_ROOT, "output/clayAdaptedAssemblies.txt"),quote = F,sep = "\t",row.names = F,col.names = F)

# loading analysis
env_metadata = read.table(file.path(PHYLOGWAS_ROOT, "data/env_metadata.txt"),sep ="\t",header = T)
ePC_loading_labels = data.frame(categories = rep(env_metadata[,3],each = 3))
rownames(ePC_loading_labels) = paste(rep(env_metadata[,1],each = 3),c("quan10","quan50","quan90"),sep = "_")
envp = pheatmap(na.omit(ePCs$variable.correlation[order(ePC_loading_labels$categories,decreasing = T),1:10]),cluster_rows = F,cluster_cols = F,show_rownames = F,annotation_row = ePC_loading_labels)

png(file.path(PHYLOGWAS_ROOT, "output/envPC_loading_20250804.png"),width = 16,height = 16, units = "cm",res = 600,pointsize = 8)
envp
dev.off()

envp2 = pheatmap(na.omit(ePCs$variable.correlation[order(ePC_loading_labels$categories,decreasing = T),1:10]),cluster_rows = F,cluster_cols = F,show_rownames = T,annotation_row = ePC_loading_labels,fontsize_row = 4)
png(file.path(PHYLOGWAS_ROOT, "output/envPC_loading_withName_20250804.png"),width = 20,height = 30, units = "cm",res = 600,pointsize = 8)
envp2
dev.off()

write.table(ePCs$variable.correlation,file.path(PHYLOGWAS_ROOT, "output/envPC_load.txt"),sep = "\t")

png(file.path(PHYLOGWAS_ROOT, "output/figure/suppFig/FigS1_envPCPipeline_c.png"),width = 15.85,height = 8.7*.5,
    unit = "cm",res = 600, pointsize = 8)
par(mar = c(3,10,1,.5),mfrow = c(1,5))
for (i in 1:5){
    tmp = ePCs$variable.correlation[,i]
    names(tmp) = rownames(ePCs$variable.correlation)
    idx = tapply(tmp,ePC_loading_labels$categories,function(x) x[rank(abs(x),na.last = F)>=(length(x)-2)])[-1]
    names(idx) = NULL
    idx = unlist(idx)
    par(mgp = c(3,.1,0))
    barplot(rev(abs(idx)),main = "",horiz = T,las =2,xlim = c(0,max(abs(idx))*1.01),cex.names = .5,cex.axis = .5,
            col = rep(scales::hue_pal()(6)[c(2,4,3,5,1)],each = 3),
            density = ifelse(rev(idx)>0,NA,50),names.arg = rev(names(idx)),xaxt = 'n')
    par(mgp = c(3,.75,0))
    axis(1,cex.axis = .5,las = 2)
}
dev.off()
