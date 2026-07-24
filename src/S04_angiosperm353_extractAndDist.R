# Extract angiosperm353 per-gene alignments from the gap-stripped per-OG CDS MSAs
# (05A_treeConstruction's own gap-strip step) and compute a K81 genetic-distance matrix.
# Feeds 05A_treeConstruction's RAxML gene-tree step (output/geneTree_angiosperm353/*.fa).
# author: Sheng-Kai Hsu
# (ported from notebook/05_phylotreeConstruction/05B_neutralPhylogenyVisualization.ipynb,
#  cells 0-15, as part of the 05A/05B restructuring)

library(ape)
library(limma)
library(parallel)

PHYLOGWAS_ROOT <- Sys.getenv("PHYLOGWAS_ROOT", unset = "/workdir/sh2246/p_phyloGWAS")

# load msa per transcript of the angiosperm 353 loci
angio353 = read.table(file.path(PHYLOGWAS_ROOT, "output/angiosperm353ToOGName.txt"),header = F)[,1]
fasta.list = list.files(file.path(PHYLOGWAS_ROOT, "output/CDSMSAPerOG_gs/"), full.names = TRUE)
fasta.list = fasta.list[substr(basename(fasta.list),1,9)%in%angio353]
fasta.aln=mclapply(fasta.list,function(x) read.FASTA(x),mc.cores = 5)
fasta.alnm=mclapply(fasta.aln,function(x) as.matrix(x),mc.cores = 2)

fasta.alnm2 = mclapply(fasta.alnm,function(x) x[grep(":0$",rownames(x)),],mc.cores = 30)

fasta.alnm2 = mclapply(fasta.alnm2,function(x) {
    tmp = x
    rownames(tmp) = strsplit2(rownames(x),":")[,2]
    return(tmp)
},mc.cores = 30)

dir.create(file.path(PHYLOGWAS_ROOT, "output/geneTree_angiosperm353/"))
for (i in seq_along(fasta.alnm2)){
    write.FASTA(fasta.alnm2[[i]],file.path(PHYLOGWAS_ROOT, paste0("output/geneTree_angiosperm353/gene_",i,".fa")))
}

# genetic distance estimate based on 353 loci
allTaxa = Reduce("union",lapply(fasta.alnm2,rownames))

fasta.alnm3 = mclapply(fasta.alnm2,function(x) {
    missTaxa = setdiff(allTaxa,rownames(x))
    tmp = matrix(rep("-",ncol(x)*length(missTaxa)),nrow = length(missTaxa),ncol = ncol(x),byrow = T)
    rownames(tmp) = missTaxa
    tmp = as.matrix(as.DNAbin(tmp))
    out = rbind(x,tmp)
    if (ncol(out) > 20) out = out[allTaxa,sample(1:ncol(out),20)]
    return(out)
},mc.cores = 30)

fasta.aln.merged=do.call(cbind, fasta.alnm3)

d.angio353 =dist.dna(fasta.aln.merged,pairwise.deletion = T,model = "K81",as.matrix = T)

d.angio353[d.angio353==Inf] = NA
d.angio353[d.angio353>=0.3] = NA

write.table(d.angio353,file.path(PHYLOGWAS_ROOT, "output/angiosperm353_geneticDistance.txt"),col.names = T,row.names = F,quote = F,sep = "\t")
