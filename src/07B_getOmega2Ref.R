# # R script for dn/ds ratio calculation for each exon/gene
# author: Sheng-Kai Hsu
# date created: 2022.11.20
# date last edited: 2024.08.15
# fix the issue of negative omega
# update a 07B version for omega to reference only
# update to take care of Ns as gap
# bug fix for edge cases (ref length not divisble by 3 & seq with < 2 sensible codon)
rm(list=ls())
library(Biostrings)
library(MSA2dist)
library(limma)
library(doParallel)

# patches to fix bug on negative dn/ds
codonmat2pnps <- function(codonmat){
  stopifnot("Error: input needs to be a codonmat of 2 columns"=
              dim(codonmat)[2] == 2)
  seq1_name <- colnames(codonmat)[1]
  seq2_name <- colnames(codonmat)[2]
  count_codons <- dim(codonmat)[1]
  # modified by Sheng-Kai Hsu to distinguish codon insertion and base insertion
  emptyCodon_idx = Reduce(union,apply(codonmat, 2,function(x) grep("---", x)))
  count_empty <- length(emptyCodon_idx)
  if(count_empty > 0){
    codonmat <- codonmat[-emptyCodon_idx, , drop = FALSE]
  }
  # Sheng-Kai's edit end #
  insertions_idx <- Reduce(union,apply(codonmat, 2,function(x) grep("-", x)))
  count_insertions <- length(insertions_idx)
  if(count_insertions > 0){
    codonmat <- codonmat[-insertions_idx, , drop = FALSE]
  }
  codonnumber <- apply(codonmat, 2, MSA2dist::codon2numberTCAG)
  Ns_idx <- unique(unlist(apply(codonnumber, 2, function(x) which(is.na(x)))))
  count_Ns <- length(Ns_idx)
  if(count_Ns > 0){
    codonmat <- codonmat[-Ns_idx, , drop = FALSE]
    codonnumber <- codonnumber[-Ns_idx, , drop = FALSE]
  }
  SA_Nei <- sum(MSA2dist::GENETIC_CODE_TCAG[codonmat[, 1], 4])
  SB_Nei <- sum(MSA2dist::GENETIC_CODE_TCAG[codonmat[, 2], 4])
  identical_codons_idx <- which(codonnumber[, 1] == codonnumber[, 2])
  identical_codons <- length(identical_codons_idx)
  if(identical_codons > 0){
    codonmat <- codonmat[-identical_codons_idx, , drop = FALSE]
    codonnumber <- codonnumber[-identical_codons_idx, , drop = FALSE]
  }
  ## At this point, codonA and codonB are "real" codons (no N's or -'s)
  ## but are not identical
  syn_codons <- 0
  nonsyn_codons <- 0
  if(nrow(codonmat) > 0){
    syn_nonsyn_codons <- apply(codonmat, 1,
                               function(x) MSA2dist::compareCodons(x[1], x[2]))
    syn_codons <- syn_codons + sum(syn_nonsyn_codons[1, ])
    nonsyn_codons <- nonsyn_codons + sum(syn_nonsyn_codons[2, ])
  }
  count_ambiguous_codons <- count_insertions + count_Ns
  count_compared_codons <- count_codons - count_ambiguous_codons
  potential_syn <- ((SA_Nei / 3) + (SB_Nei / 3)) / 2
  potential_nonsyn <- (3 * count_compared_codons) - potential_syn
  ps <- syn_codons / potential_syn
  pn <- nonsyn_codons / potential_nonsyn
  ds <- ifelse(ps < 0.75 , ((-3 / 4) * log(1 - (4 * (ps / 3)))), NA)
  dn <- ifelse(pn < 0.75, ((-3 / 4) * log(1 - (4 * (pn / 3)))), NA)
  dnds <- NA
  if(dn != 0 & ds != 0 & !is.na(dn) & !is.na(ds)){ dnds <- dn/ds}
  pnps <- NA
  if(pn != 0 & ps != 0 & !is.na(pn) & !is.na(ps)){ pnps <- pn/ps}
  codonmat_out <- setNames(c(seq1_name, seq2_name, count_codons,
                         count_compared_codons, count_ambiguous_codons, count_insertions,
                         count_Ns, syn_codons, nonsyn_codons, potential_syn, potential_nonsyn,
                         ps, pn, pnps, ds, dn, dnds),
                       c("seq1", "seq2", "Codons", "Compared", "Ambigiuous", "Indels", "Ns",
                         "Sd", "Sn", "S", "N", "ps", "pn", "pn/ps", "ds", "dn", "dn/ds"))
  attr(codonmat_out, "class") <- "pnps"
  return(codonmat_out)
}

dnastring2kaks_modified <- function(cds,
                           model = "Li",
                           threads = 1,
                           isMSA = TRUE,
                           sgc = "1",
                           reference,
                           verbose = FALSE,
                           ...) {
  stopifnot("Error: input needs to be a DNAStringSet" = methods::is(cds, "DNAStringSet"))
  stopifnot("Error: either choose model 'Li' or 'NG86' or KaKs_Calculator2 model" =
              model %in% c("Li", "NG86", "NG", "LWL", "LPB", "MLWL", "MLPB", "GY", "YN", "MYN", "MS", "MA", "GNG", "GLWL", "GLPB", "GMLWL", "GMLPB", "GYN", "GMYN"))
  
  # Identify the reference sequence
  ref_index <- if (is.numeric(reference)) {
    reference
  } else {
    tmp = grep(reference, names(cds))
    tmp[which.max(width(cds[tmp]))]
  }
  
  if (is.na(ref_index) || ref_index > length(cds)) {
    stop("Error: Invalid reference sequence specified")
  }
  
  cds.names <- names(cds)
  cds.names <- gsub(" ", "_", gsub(" $", "", gsub("\\s+", " ", cds.names)))
  names(cds) <- cds.names
  
  if (model == "Li") {
    if (isMSA) {
      OUT <- seqinr::kaks(dnastring2aln(cds))
      OUT.ka <- as.matrix(OUT$ka)
      OUT.ks <- as.matrix(OUT$ks)
      OUT.vka <- as.matrix(OUT$vka)
      OUT.vks <- as.matrix(OUT$vks)
      
      ka_df <- data.frame(seq2 = colnames(OUT.ka), ka = OUT.ka[ref_index, ])
      ks_df <- data.frame(seq2 = colnames(OUT.ks), ks = OUT.ks[ref_index, ])
      vka_df <- data.frame(seq2 = colnames(OUT.vka), vka = OUT.vka[ref_index, ])
      vks_df <- data.frame(seq2 = colnames(OUT.vks), vks = OUT.vks[ref_index, ])
      
      OUT <- Reduce(function(x, y) merge(x, y, by = "seq2"), list(ka_df, ks_df, vka_df, vks_df))
      OUT <- OUT %>% dplyr::mutate(seq1 = cds.names[ref_index])
      OUT <- as.data.frame(OUT)
      attr(OUT, "model") <- "Li"
      attr(OUT, "align") <- "FALSE"
      attr(OUT, "MSA2dist.class") <- "dnastring2kaks"
      return(OUT)
    } else {
      if (.Platform$OS.type == "windows") {
        cl <- parallel::makeCluster(threads)
      } else {
        cl <- parallel::makeForkCluster(threads)
      }
      doParallel::registerDoParallel(cl)
      
      OUT <- foreach(j = seq_len(length(cds)), .combine = rbind, .packages = c('foreach')) %dopar% {
        if (j != ref_index) {
          res <- unlist(seqinr::kaks(MSA2dist::dnastring2aln(MSA2dist::cds2codonaln(cds[ref_index], cds[j], ...))))
          data.frame(Comp1 = ref_index, Comp2 = j, t(as.data.frame(res)))
        }
      }
      parallel::stopCluster(cl)
      OUT <- as.data.frame(OUT)
      attr(OUT, "model") <- "Li"
      attr(OUT, "align") <- "TRUE"
      attr(OUT, "MSA2dist.class") <- "dnastring2kaks"
      return(OUT)
    }
  } else if (model == "NG86") {
    if (isMSA) {
      if (.Platform$OS.type == "windows") {
        cl <- parallel::makeCluster(threads)
      } else {
        cl <- parallel::makeForkCluster(threads)
      }
      doParallel::registerDoParallel(cl)
      
      codonmat <- MSA2dist::dnastring2codonmat(cds)
      OUT <- foreach(j = seq_len(ncol(codonmat)), .combine = rbind, .packages = c('foreach')) %dopar% {
        # if (j != ref_index) {
          res <- codonmat2pnps(codonmat[, c(ref_index, j)])
          data.frame(Comp1 = ref_index, Comp2 = j, t(as.data.frame(res)))
        # }
      }
      parallel::stopCluster(cl)
      OUT <- as.data.frame(OUT)
      attr(OUT, "model") <- "NG86"
      attr(OUT, "align") <- "FALSE"
      attr(OUT, "MSA2dist.class") <- "dnastring2kaks"
      return(OUT)
    } else {
      if (.Platform$OS.type == "windows") {
        cl <- parallel::makeCluster(threads)
      } else {
        cl <- parallel::makeForkCluster(threads)
      }
      doParallel::registerDoParallel(cl)
      
      OUT <- foreach(j = seq_len(length(cds)), .combine = rbind, .packages = c('foreach')) %dopar% {
        # if (j != ref_index) {
          res <- MSA2dist::codonmat2pnps(MSA2dist::dnastring2codonmat(MSA2dist::cds2codonaln(cds[ref_index], cds[j], ...)))
          data.frame(Comp1 = ref_index, Comp2 = j, seq1 = cds.names[ref_index], seq2 = cds.names[j], t(res))
        # }
      }
      parallel::stopCluster(cl)
      OUT <- as.data.frame(OUT)
      attr(OUT, "model") <- "NG86"
      attr(OUT, "align") <- "TRUE"
      attr(OUT, "MSA2dist.class") <- "dnastring2kaks"
      return(OUT)
    }
  } else {
    if (isMSA) {
      OUT <- rcpp_KaKs(cdsstr = as.character(cds), sgc = sgc, method = model, verbose = verbose)
      OUT <- as.data.frame(t(tibble::column_to_rownames(tidyr::as_tibble(setNames(stringr::str_split(OUT$results_vec, "\t"), OUT$results_names)) %>%
                                                          tibble::add_column(rownames = OUT$rownames), "rownames")))
      attr(OUT, "model") <- model
      attr(OUT, "align") <- "FALSE"
      attr(OUT, "MSA2dist.class") <- "dnastring2kaks"
      return(OUT)
    } else {
      if (.Platform$OS.type == "windows") {
        cl <- parallel::makeCluster(threads)
      } else {
        cl <- parallel::makeForkCluster(threads)
      }
      doParallel::registerDoParallel(cl)
      
      OUT <- foreach(j = seq_len(length(cds)), .combine = rbind, .packages = c('foreach')) %dopar% {
        # if (j != ref_index) {
          tmp_out <- rcpp_KaKs(cdsstr = as.character(MSA2dist::cds2codonaln(cds[ref_index], cds[j], ...)), sgc = sgc, method = model)
          tmp_out <- as.data.frame(t(tibble::column_to_rownames(tidyr::as_tibble(setNames(stringr::str_split(tmp_out$results_vec, "\t"), tmp_out$results_names)) %>%
                                                                  tibble::add_column(rownames = tmp_out$rownames), "rownames")))
          tmp_out["Comp1"] <- ref_index
          tmp_out["Comp2"] <- j
          tmp_out
        # }
      }
      parallel::stopCluster(cl)
      OUT <- as.data.frame(OUT)
      attr(OUT, "model") <- model
      attr(OUT, "align") <- "TRUE"
      attr(OUT, "MSA2dist.class") <- "dnastring2kaks"
      return(OUT)
    }
  }
}


args <- commandArgs(TRUE)

if (args[1] == "--help") {
  message("################################ HELP ################################\n
          ## Written by Sheng-Kai Hsu\n## Arguments :\n
          # --input (required)\n\n
          #--output path to output, default is current directory (required)\n
          # --ref name of the reference sequence (required) \n
          \n######################################################################")
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
  }else if (args[i] %in% c("--ref")) {
    ref = args[i + 1]
  }
  else{
    message("\n### ERROR ###\n\n Argument not recognized:")
    print(args[i])
    message("\n### ERROR ###\n")
    q("no")
  }
}

#loading
# inDir = "~/Dropbox/postDoc/projects/p_phyloGWAS/output/OG0005032_mafft.fa"
# ref = "ASM1935983v1"

dat=readDNAStringSet(inDir,format = "fasta")

# gap-stripping
refIdx <- if (is.numeric(ref)) {
  ref
} else {
  tmp = grep(ref, names(dat))
  tmp[which.max(width(dat[tmp]))]
}
print(paste("Reference:",names(dat)[refIdx]))
dat = dnastring2dnabin(dat)
gapInRef = which(as.character(dat[refIdx,]) == "-")
if (length(gapInRef)>0) dat = dat[,-gapInRef] else dat = dat
dat = dnabin2dnastring(dat)

# if(!grepl("^ATG",dat[length(dat)])) dat = reverseComplement(dat)

# find non-empty taxa
# nonemptyIdx=which(sapply(dat,function(x) grepl("[ATCG]",x)))
secondFilterIdx=!apply(dnastring2codonmat(dat),2,function(x) sum(!grepl("-",x))<=1)

# calculate query coverage
nonGapLength = sapply(dat,function(x) sum(strsplit2(x,"")!='-'&strsplit2(x,"")!='N'))
queryCov = nonGapLength/nonGapLength[refIdx]

if(width(dat[refIdx])%%3==0){
  # dnds calculation
  dnds=dnastring2kaks_modified(dat[nonGapLength>5&secondFilterIdx,],model="NG86", threads=1,reference = ref)
  
  # output
  out = cbind(dnds[,-c(1:2)],queryCov[nonGapLength>5&secondFilterIdx])
  write.table(out,outDir,quote = F,sep = "\t",col.names = F,row.names = F)
}


