args = commandArgs(trailingOnly=TRUE)

# pkgTest <- function(x)
#   {
#     if (!require(x,character.only = TRUE))
#     {
#       install.packages(x,dep=TRUE)
#         if(!require(x,character.only = TRUE)) stop("Package not found")
#     }
#   }
# 
# pkgTest("ape")

# load libraries
library(ape)
library(phytools)
library(phangorn)
library(testthat)

# define function that identify tip and node to label
getNodeToLabel = function(Targets,tree){
  toLabel = c()
  while (length(Targets)>0){
    focalTarget = max(Targets)
    tmp.sis = getSisters(tree,focalTarget)
    tmp.sis.child = c()
    if (length(tmp.sis)==0){
      toLabel = c(focalTarget)
      Targets = c()
    } else if (length(tmp.sis)>=2) {
      tmp.sis.child = unlist(lapply(tmp.sis,function(x) {
        if(x > length(tree$tip.label)) unlist(Descendants(tree,x,type = "tips"))
      }))
      if(any(!tmp.sis.child%in%Targets)) {
        toLabel = c(toLabel,focalTarget)
        Targets = setdiff(Targets,focalTarget)
      } else {
        Targets = c(Targets,getMRCA(tree,c(focalTarget,tmp.sis)))
        Targets = setdiff(Targets,c(focalTarget,tmp.sis))
      }
    } else if(tmp.sis > length(tree$tip.label) & !tmp.sis%in%Targets) {
      tmp.sis.child = unlist(Descendants(tree,tmp.sis,type = "tips"))
      if(any(!tmp.sis.child%in%Targets)) {
        toLabel = c(toLabel,focalTarget)
        Targets = setdiff(Targets,focalTarget)
      } else {
        Targets = c(Targets,getMRCA(tree,c(focalTarget,tmp.sis)))
        Targets = setdiff(Targets,c(focalTarget,tmp.sis,tmp.sis.child))
      }
    } else{
      if(!tmp.sis%in%Targets) {
        toLabel = c(toLabel,focalTarget)
        Targets = setdiff(Targets,focalTarget)
      } else {
        Targets = c(Targets,getMRCA(tree,c(focalTarget,tmp.sis)))
        Targets = setdiff(Targets,c(focalTarget,tmp.sis))
      }
    }
  }
  return(sort(toLabel))
}

#### unit test ####
if ("unitTest"%in%args){
  set.seed(123)
  testTree = rtree(10)
  testTree$edge.length[3:4] = 0
  testTree = makeNodeLabel(testTree)
  
  print("test: no tip in target")
  testTarget = c()
  toLabelTest = getNodeToLabel(testTarget,testTree)
  toLabelTruth = c()
  test_that("two sister nodes (both in target)",{
    expect_equal(toLabelTest,toLabelTruth)
  })
  
  print("test: all tips in target")
  testTarget = c(1:10)
  toLabelTest = getNodeToLabel(testTarget,testTree)
  toLabelTruth = c(11)
  test_that("two sister nodes (both in target)",{
    expect_equal(toLabelTest,toLabelTruth)
  })
  
  print("test: all branches of interest at tips")
  testTarget = c(1,4,7)
  toLabelTest = getNodeToLabel(testTarget,testTree)
  toLabelTruth = c(1,4,7)
  test_that("test: all branches of interest at tips",{
    expect_equal(toLabelTest,toLabelTruth)
  })
  
  print("test: all branches of interest none at tips")
  testTarget = c(1,2,3,4,5)
  toLabelTest = getNodeToLabel(testTarget,testTree)
  toLabelTruth = c(12,16)
  test_that("test: all branches of interest none at tips",{
    expect_equal(toLabelTest,toLabelTruth)
  })
  
  print("test: full clade all in target (& two sister nodes in the search)")
  testTarget = c(1,2,3,4)
  toLabelTest = getNodeToLabel(testTarget,testTree)
  toLabelTruth = c(4,12)
  test_that("full clade all in target",{
    expect_equal(toLabelTest,toLabelTruth)
  })
  
  print("test: a clade but one tip (& two sister nodes in the search)")
  testTarget = c(2,3,4)
  toLabelTest = getNodeToLabel(testTarget,testTree)
  toLabelTruth = c(2,3,4)
  test_that("a clade but one tip",{
    expect_equal(toLabelTest,toLabelTruth)
  })
  print("unit test done")
  q()
}

#### deployment ####
# load data
tree <- read.tree(args[1])
rootName = tree$tip.label[grep(args[2],tree$tip.label)[1]]
tree <- root(tree, rootName)
TargetsName <- read.table(args[3])[,1]

# tree <- read.tree("Dropbox/postDoc/projects/p_phyloGWAS/output/RAxML_bestTree.OG0000928.tree")
# tree <- root(tree, "ASM1935983v1")
# TargetsName <- read.table("Dropbox/postDoc/projects/p_phyloGWAS/output/OG0000928.target")[,1]

# drop tips not in MSA
# msa = read.FASTA(args[4])
# # msa = read.FASTA("Dropbox/postDoc/projects/p_phyloGWAS/output/subset_clean.fa")
# tipToKeep = tree$tip.label[gsub("-|\\.","_",tree$tip.label)%in%names(msa)]
# tree = keep.tip(tree, tipToKeep)
Targets = which(tree$tip.label%in%TargetsName)
toLabel = getNodeToLabel(Targets,tree)

# validation testing
print("validation test:")
print("Capture all targets?")
test_that("all descendent in targets",{
  expect_equal(sort(unlist(Descendants(tree,toLabel,type = "tips"))),Targets)
})
print("No nested nodes/tips?")
test_that("upper most node",{
  # expect_true(all(sapply(toLabelTest,function(x) all(unlist(Descendants(testTree,x,type = "children")) <= x))))
  expect_false(any(sapply(toLabel,function(x) any(setdiff(unlist(Descendants(tree,x,type = "all")),Targets) %in% toLabel))))
})

#label the tree
treeLabeled = tree 
treeLabeled$node.label=rep("{Background}",treeLabeled$Nnode)
treeLabeled$node.label[toLabel[toLabel>length(tree$tip.label)]-length(tree$tip.label)]="{Foreground}"
treeLabeled$tip.label <- paste0(tree$tip.label,"{Background}")
treeLabeled$tip.label[toLabel[toLabel<=length(tree$tip.label)]] <- paste0(tree$tip.label[toLabel[toLabel<=length(tree$tip.label)]],"{Foreground}")
treeLabeled$tip.label = gsub("-|\\.","_",treeLabeled$tip.label)

write.tree(treeLabeled,paste0(args[4]))
