# R function to derive phylogenetic K matrix from a tree based on the shared branch length
# author: Sheng-Kai Hsu
# date created: 2022.07.29
# date last edited: 2022.08.05
library(ape)

phyloK=function(tree){
  if (class(tree)!="phylo") stop(print("tree not in correct class"))
  mrca_idx=mrca(tree)
  k=apply(mrca_idx,c(1,2),function(x) {
    anc_b_idx=c()
    ancNode=x
    while(ancNode!=length(tree$tip.label)+1){
      anc_b_idx=c(anc_b_idx,which(tree$edge[,2]==ancNode))
      ancNode=tree$edge[,1][tree$edge[,2]==ancNode]
    }
    out=sum(tree$edge.length[anc_b_idx])
    max_brlen=max(node.depth.edgelength(tree))
    out=2*out/max_brlen
    return(out)
  })
}
