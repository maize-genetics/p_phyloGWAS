# R function to wrap get4d.msa in rphast for parallelization
# author: Sheng-Kai Hsu
# date created: 2022.08.25
# date last edited: 2022.08.25

library(rphast)

args <- commandArgs(TRUE)

if (args[1] == "--help") {
  message("################################ HELP ################################\n## Written by Sheng-Kai Hsu\n## Arguments :\n# --input (required)\n\n#--output path to output, default is current directory (required)\n#--feature feature gff file\n\n######################################################################")
  q("no")
}

outDir = getwd()
i = -1
while (i < I(length(args)-1)) {
  i = i + 2
  if (args[i] %in% c("--input")) {
    inDir = args[i + 1]
  }else if (args[i] %in% c("--output")) {
    outDir = args[i + 1]
  }else if (args[i] %in% c("--feature")) {
    featDir = args[i + 1]
  }
  else{
    message("\n### ERROR ###\n\n Argument not recognized:")
    print(args[i])
    message("\n### ERROR ###\n")
    q("no")
  }
}

testMSA=read.msa(inDir)
feat=read.feat(featDir)
out=get4d.msa(testMSA,feat)

write.msa(out,file = outDir,format = "FASTA") 
