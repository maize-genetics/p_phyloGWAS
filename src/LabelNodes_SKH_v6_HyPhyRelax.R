args = commandArgs(trailingOnly=TRUE)

pkgTest <- function(x)
  {
    if (!require(x,character.only = TRUE))
    {
      install.packages(x,dep=TRUE)
        if(!require(x,character.only = TRUE)) stop("Package not found")
    }
  }

pkgTest("ape")

library(ape)
library(phytools)
library(treeio)
library(phangorn)

tree <- read.tree(args[1])
tree <- root(tree, args[2])
Targets <- read.table(args[3])[,1]

# tree <- read.tree("Dropbox/postDoc/projects/p_phyloGWAS/output/RAxML_bestTree.OG0000928.tree")
# tree <- root(tree, "ASM1935983v1")
# Targets <- read.table("Dropbox/postDoc/projects/p_phyloGWAS/output/OG0000928.target")[,1]

# drop tips not in MSA
#msa = read.FASTA(args[4])
# msa = read.FASTA("Dropbox/postDoc/projects/p_phyloGWAS/output/subset_clean.fa")
#tipToKeep = tree$tip.label[gsub("-|\\.","_",tree$tip.label)%in%names(msa)]
#tree = keep.tip(tree, tipToKeep)

# identify tip and node to label
Targets = which(tree$tip.label%in%Targets)
# t0=Sys.time()
toLabel = c()
while (length(Targets)>0){
  tmp.sis = getSisters(tree,Targets[1])
  tmp.sis.child = c()
  if (length(tmp.sis)>=2) {
    tmp.sis.child = unlist(lapply(tmp.sis,function(x) {
      if(x > length(tree$tip.label)) unlist(Descendants(tree,x,type = "tips"))
      }))
  } else if(tmp.sis > length(tree$tip.label)) tmp.sis.child = unlist(Descendants(tree,tmp.sis,type = "tips"))
  
  if(any(!tmp.sis.child%in%Targets)&any(!tmp.sis%in%Targets)) {
    toLabel = c(toLabel,Targets[1])
    Targets = setdiff(Targets,Targets[1])
  } else {
    Targets = c(Targets,getMRCA(tree,c(Targets[1],tmp.sis)))
    Targets = setdiff(Targets,c(Targets[1],tmp.sis))
  }
}
# t1=Sys.time()
# t1-t0

#label the tree
treeLabeled = tree 
treeLabeled$node.label=rep("{Background}",treeLabeled$Nnode)
treeLabeled$node.label[toLabel[toLabel>length(tree$tip.label)]-length(tree$tip.label)]="{Foreground}"
# treeLabeled$node.label[is.na(treeLabeled$node.label)]="{Background}"
treeLabeled$tip.label <- paste0(tree$tip.label,"{Background}")
treeLabeled$tip.label[toLabel[toLabel<=length(tree$tip.label)]] <- paste0(tree$tip.label[toLabel[toLabel<=length(tree$tip.label)]],"{Foreground}")
treeLabeled$tip.label = gsub("-|\\.","_",treeLabeled$tip.label)

# write.tree(treeLabeled,"p_panAndOGASR/paml_run/test.tre")
write.tree(treeLabeled,paste0(args[5]))
# png(paste0("RAxML_Labeled_bestTree.ungapped_",args[1],".png"), width=800, height=800)
# plot(treeLabeled, main = paste0("RAxML_bestTree.ungapped_",args[1]))
# add.scale.bar(cex = 0.8, font = 2, col = "red")
# dev.off()
