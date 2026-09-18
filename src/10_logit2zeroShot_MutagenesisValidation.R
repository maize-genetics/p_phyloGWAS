rm(list = ls())
# Use reticulate to import numpy. Any python with numpy works; override with
#   export RETICULATE_PYTHON=/path/to/python
# before calling Rscript if /usr/bin/python lacks numpy on your machine.
if (Sys.getenv("RETICULATE_PYTHON") == "") Sys.setenv(RETICULATE_PYTHON = "/usr/bin/python")
library(reticulate)
library(Biostrings)

py_config()
np <- import("numpy")

args = commandArgs(trailingOnly=TRUE)

inputLogit = args[1]
inputFASTA = args[2]
output = args[3]
# regex used to derive the gene/group key from a seqID, for per-group z-scaling.
# default matches At mutant naming ("<gene>_mut_N" / "<gene>_neutralMut_N") and WT
# FASTA headers (space-delimited); override for datasets with a different scheme.
groupRegex = ifelse(length(args) >= 4, args[4], "_| ")

AA20 = limma::strsplit2("LAGVSERTIDPKQNFYMHWC","")[1,]

# load numpy
npz <- np$load(inputLogit)
prob <- lapply(npz$files, function(x) npz[[x]])
names(prob) <- npz$files

# load fasta
testfasta = readAAStringSet(inputFASTA)
# seqID -> sequence, keyed the same way `4_ESM_logits.py` keys the npz (first
# whitespace-delimited token of the FASTA header, i.e. Biopython's record.id)
fastaID = sub("\\s.*$", "", names(testfasta))
seqByID = setNames(as.character(testfasta), fastaID)

# zero shot score calc.
# score(seq) = mean over standard-AA sites of log( P(obs) / max_{a != obs} P(a) )
zeroshot = rep(NA_real_, length(prob))
names(zeroshot) = names(prob)
for(j in seq_along(prob)) { # loop through sequences
  if(j%%100==0) print(j)
  seqID = names(prob)[j]
  seqStr = seqByID[[seqID]] # match logits to sequence BY NAME, not position
  if (is.null(seqStr)) {
    warning("No FASTA sequence found for logit entry '", seqID, "' - skipping")
    next
  }
  tmp = strsplit(seqStr,"")[[1]]
  P = prob[[j]]
  colnames(P) = AA20

  siteIdx = which(tmp %in% AA20)
  siteIdx = siteIdx[siteIdx <= nrow(P)] # guard against length mismatch
  obsCol = match(tmp[siteIdx], AA20)

  refProb = P[cbind(siteIdx, obsCol)] # P(obs)
  # P.max(others): mask out the observed column per site, then take the row max
  Pmasked = P[siteIdx, , drop = FALSE]
  Pmasked[cbind(seq_along(siteIdx), obsCol)] = -Inf
  secProb = apply(Pmasked, 1, max)

  zeroshot[j] = mean(log(refProb/secProb)) # mean log ratio of the ps
}

groupKey = limma::strsplit2(names(zeroshot), groupRegex)[,1]
zeroshot_scaled = ave(zeroshot, groupKey, FUN = function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))
names(zeroshot_scaled) = names(zeroshot)

outTab = data.frame(seqID = names(prob), zeroShotScore = zeroshot, scaledZeroShotScore = zeroshot_scaled)

write.table(outTab, output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
