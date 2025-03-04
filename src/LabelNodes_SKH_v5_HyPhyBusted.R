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

# tree <- read.tree("Dropbox/postDoc/projects/p_panAndOGASR/output/Pavag02G345800.aln.gs.fw.l0.fa_dir/RAxML_bestTree.ungapped_Pavag02G345800.aln.gs.fw.l0.fa")
# tree <- root(tree, "Pavag02G345800.1.v3.1")
# Targets <- read.table("Dropbox/postDoc/projects/p_panAndOGASR/output/Pavag02G345800.aln.gs.fw.l0.fa_dir/Pavag02G345800.target")[,1]

Targets = which(tree$tip.label%in%Targets)
# t0=Sys.time()
toLabel = c()
while (length(Targets)>0){
  tmp.sis = getSisters(tree,Targets[1])
  if (length(tmp.sis)>=2) {
    tmp.sis = unlist(lapply(tmp.sis,function(x) {
      if(x > length(tree$tip.label)) unlist(Descendants(tree,x,type = "tips"))
      }))
  } else if(tmp.sis > length(tree$tip.label)) tmp.sis = unlist(Descendants(tree,tmp.sis,type = "tips"))
  
  if(any(!tmp.sis%in%Targets)) {
    toLabel = c(toLabel,Targets[1])
    Targets = setdiff(Targets,Targets[1])
  } else {
    Targets = c(Targets,getMRCA(tree,c(Targets[1],tmp.sis)))
    Targets = setdiff(Targets,c(Targets[1],tmp.sis))
  }
}
# t1=Sys.time()
# t1-t0

treeLabeled = tree 
treeLabeled$node.label=rep("",treeLabeled$Nnode)
treeLabeled$node.label[toLabel[toLabel>length(tree$tip.label)]-length(tree$tip.label)]="{Foreground}"
treeLabeled$tip.label[toLabel[toLabel<=length(tree$tip.label)]] <- paste0(tree$tip.label[toLabel[toLabel<=length(tree$tip.label)]],"{Foreground}")
treeLabeled$tip.label = gsub("-|\\.","_",treeLabeled$tip.label)

# write.tree(treeLabeled,"p_panAndOGASR/paml_run/test.tre")
write.tree(treeLabeled,paste0(args[4]))
# png(paste0("RAxML_Labeled_bestTree.ungapped_",args[1],".png"), width=800, height=800)
# plot(treeLabeled, main = paste0("RAxML_bestTree.ungapped_",args[1]))
# add.scale.bar(cex = 0.8, font = 2, col = "red")
# dev.off()
