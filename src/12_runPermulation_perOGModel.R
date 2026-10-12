library(data.table)
library(tidyverse)
library(parallel)
library(ape)
library(geiger)
library(taxize)
library(dplyr)

library(asreml)
library(asremlPlus)

fit_model_v4 <- function(OG_id, data, Kmat,traitVec) {
  #   library(asreml)
  # gc()
  # Cfixed/gammaPar must be set before fitting, not left to infoCriteria()'s lazy
  # refit: when Cfixed is NULL, infoCriteria(IClikelihood="full") triggers
  # asreml::update.asreml() to recompute it, and that refit's call re-evaluation
  # breaks in this pipeline's forked-worker context (surfaces as either
  # "'vpc.char' is not an exported object from 'namespace:asreml'" or, depending on
  # the model formula, an "object not found" scoping error) -- confirmed via direct
  # reproduction across asreml 4.2.0.302, 4.2.0.355, and two separate conda envs, so
  # this is not a version-pin issue. Setting Cfixed up front avoids the refit
  # entirely. See revision-analyses.md / environment-r.yml for the full story.
  asreml.options(verbose = FALSE, Cfixed = TRUE, gammaPar = FALSE)
  # print(OG_id)
  filtered_data <- subset(data, OG == OG_id)
  filtered_data <- filtered_data[,c("assemblyID","PAV","scaled.dn.ds","scaledESM", "scaledPlantCAD", "PMS",
                                    "compdup", "queryCov","pacbio")]
  filtered_data = merge(filtered_data,traitVec,by = "assemblyID")
  tmp = filtered_data

  if (length(unique(filtered_data$PAV))==2){
    fullFM = as.formula(paste0(colnames(traitVec)[2], "~ PAV + PAV:scaled.dn.ds + PAV:scaledESM + PAV:scaledPlantCAD + PAV:PMS + compdup + PAV:queryCov + pacbio"))
    reducedFM = as.formula(paste0(colnames(traitVec)[2], "~ compdup + queryCov + pacbio"))
    
    model_full <- asreml(fixed = fullFM,
                         random = ~ vm(assemblyID, Kmat) , ai.sing = F, data = filtered_data)
    # model_reduced1 <- asreml(fixed = fullFM,
    #                          ai.sing = F, data = filtered_data)
    filtered_data <- tmp
    model_reduced <- asreml(fixed = reducedFM,
                             random = ~ vm(assemblyID, Kmat) , ai.sing = F, data = filtered_data)  
    
    totalVar = var(filtered_data[,ncol(filtered_data)])
    modelSummary <- summary(model_full)
    IC_full = infoCriteria(model_full,IClikelihood = "full")
    IC_reduced = infoCriteria(model_reduced,IClikelihood = "full")
    logLik_full <- IC_full$loglik
    logLik_reduced <- IC_reduced$loglik
    BIC_full <- IC_full$BIC
    BIC_reduced <- IC_reduced$BIC
    
    LRT_statistic <- -2 * (logLik_reduced - logLik_full)
    p_value <- pchisq(LRT_statistic, df = model_reduced$nedf-model_full$nedf, lower.tail = FALSE)
    # gene_explained <-  (model_reduced2$sigma2-model_full$sigma2)/totalVar
    phylo_explained <- modelSummary$varcomp$component[1]/sum(modelSummary$varcomp$component)
    
    modelWald <- wald.asreml(model_full)

    partialCoeff <- as.data.frame(model_full$coefficients$fixed)[paste0("PAV_1",c("",":scaled.dn.ds",":scaledESM",":scaledPlantCAD",":PMS_1")),1]
    partialPval <- as.data.frame(modelWald)[paste0("PAV",c("",":scaled.dn.ds",":scaledESM",":scaledPlantCAD",":PMS")),4]
    out = c(OG_id,BIC_full,BIC_reduced,LRT_statistic,p_value,phylo_explained,totalVar,partialCoeff,partialPval)
    names(out) = c("OG","BIC_full","BIC_reduced","LR","p","phyloVAE","totalVar",
                   paste("partialCoeff",c("PAV","dNdS","ESM2","plantCad","PMS"),sep = "_"),
                   paste("partialP",c("PAV","dNdS","ESM2","plantCad","PMS"),sep = "_"))
    
  } else {
    fullFM = as.formula(paste0(colnames(traitVec)[2], "~ scaled.dn.ds + scaledESM + scaledPlantCAD + PMS + compdup + queryCov + pacbio"))
    reducedFM = as.formula(paste0(colnames(traitVec)[2], "~ compdup + queryCov + pacbio"))
    
    model_full <- asreml(fixed = fullFM,
                         random = ~ vm(assemblyID, Kmat) , ai.sing = F, data = filtered_data)
    # model_reduced1 <- asreml(fixed = fullFM,
    #                          ai.sing = F, data = filtered_data)
    filtered_data <- tmp
    model_reduced <- asreml(fixed = reducedFM,
                            random = ~ vm(assemblyID, Kmat) , ai.sing = F, data = filtered_data)  
    
    totalVar = var(filtered_data[,ncol(filtered_data)])
    modelSummary <- summary(model_full)
    IC_full = infoCriteria(model_full,IClikelihood = "full")
    IC_reduced = infoCriteria(model_reduced,IClikelihood = "full")
    logLik_full <- IC_full$loglik
    logLik_reduced <- IC_reduced$loglik
    BIC_full <- IC_full$BIC
    BIC_reduced <- IC_reduced$BIC
    
    LRT_statistic <- -2 * (logLik_reduced - logLik_full)
    p_value <- pchisq(LRT_statistic, df = model_reduced$nedf-model_full$nedf, lower.tail = FALSE)
    # gene_explained <-  (model_reduced2$sigma2-model_full$sigma2)/totalVar
    phylo_explained <- modelSummary$varcomp$component[1]/sum(modelSummary$varcomp$component)
    
    modelWald <- wald.asreml(model_full)
    partialCoeff <- as.data.frame(model_full$coefficients$fixed)[c("PAV","scaled.dn.ds","scaledESM","scaledPlantCAD","PMS_1"),1]
    partialPval <- as.data.frame(modelWald)[c("PAV","scaled.dn.ds","scaledESM","scaledPlantCAD","PMS"),4]
    out = c(OG_id,BIC_full,BIC_reduced,LRT_statistic,p_value,phylo_explained,totalVar,
            partialCoeff,partialPval)
  }
  names(out) = c("OG","BIC_full","BIC_reduced","LR","p","phyloVAE","totalVar",
                 paste("partialCoeff",c("PAV","dNdS","ESM2","plantCad","PMS"),sep = "_"),
                 paste("partialP",c("PAV","dNdS","ESM2","plantCad","PMS"),sep = "_"))
  
  
  #     return(model_full)
  if(model_full$converge) return(out)
  else return(NULL)
}
fit_model_v4_GLM <- function(OG_id, data, Kmat,traitVec) {
  #   library(asreml)
  # gc()
  # Cfixed/gammaPar must be set before fitting, not left to infoCriteria()'s lazy
  # refit: when Cfixed is NULL, infoCriteria(IClikelihood="full") triggers
  # asreml::update.asreml() to recompute it, and that refit's call re-evaluation
  # breaks in this pipeline's forked-worker context (surfaces as either
  # "'vpc.char' is not an exported object from 'namespace:asreml'" or, depending on
  # the model formula, an "object not found" scoping error) -- confirmed via direct
  # reproduction across asreml 4.2.0.302, 4.2.0.355, and two separate conda envs, so
  # this is not a version-pin issue. Setting Cfixed up front avoids the refit
  # entirely. See revision-analyses.md / environment-r.yml for the full story.
  asreml.options(verbose = FALSE, Cfixed = TRUE, gammaPar = FALSE)
  # print(OG_id)
  filtered_data <- subset(data, OG == OG_id)
  filtered_data <- filtered_data[,c("assemblyID","PAV","scaled.dn.ds","scaledESM", "scaledPlantCAD", "PMS",
                                    "compdup", "queryCov","pacbio")]
  filtered_data = merge(filtered_data,traitVec,by = "assemblyID")
  tmp = filtered_data
  
  if (length(unique(filtered_data$PAV))==2){
    fullFM = as.formula(paste0(colnames(traitVec)[2], "~ PAV + PAV:scaled.dn.ds + PAV:scaledESM + PAV:scaledPlantCAD + PAV:PMS + compdup + PAV:queryCov + pacbio"))
    reducedFM = as.formula(paste0(colnames(traitVec)[2], "~ compdup + queryCov + pacbio"))
    
    model_full <- asreml(fixed = fullFM, ai.sing = F, data = filtered_data)
    # model_reduced1 <- asreml(fixed = fullFM,
    #                          ai.sing = F, data = filtered_data)
    filtered_data <- tmp
    model_reduced <- asreml(fixed = reducedFM, ai.sing = F, data = filtered_data)  
    
    totalVar = var(filtered_data[,ncol(filtered_data)])
    modelSummary <- summary(model_full)
    IC_full = infoCriteria(model_full,IClikelihood = "full")
    IC_reduced = infoCriteria(model_reduced,IClikelihood = "full")
    logLik_full <- IC_full$loglik
    logLik_reduced <- IC_reduced$loglik
    BIC_full <- IC_full$BIC
    BIC_reduced <- IC_reduced$BIC
    
    LRT_statistic <- -2 * (logLik_reduced - logLik_full)
    p_value <- pchisq(LRT_statistic, df = model_reduced$nedf-model_full$nedf, lower.tail = FALSE)
    # gene_explained <-  (model_reduced2$sigma2-model_full$sigma2)/totalVar
    phylo_explained <- modelSummary$varcomp$component[1]/sum(modelSummary$varcomp$component)
    
    modelWald <- wald.asreml(model_full)
    
    partialCoeff <- as.data.frame(model_full$coefficients$fixed)[paste0("PAV_1",c("",":scaled.dn.ds",":scaledESM",":scaledPlantCAD",":PMS_1")),1]
    partialPval <- as.data.frame(modelWald)[paste0("PAV",c("",":scaled.dn.ds",":scaledESM",":scaledPlantCAD",":PMS")),4]
    out = c(OG_id,BIC_full,BIC_reduced,LRT_statistic,p_value,phylo_explained,totalVar,partialCoeff,partialPval)
    names(out) = c("OG","BIC_full","BIC_reduced","LR","p","phyloVAE","totalVar",
                   paste("partialCoeff",c("PAV","dNdS","ESM2","plantCad","PMS"),sep = "_"),
                   paste("partialP",c("PAV","dNdS","ESM2","plantCad","PMS"),sep = "_"))
    
  } else {
    fullFM = as.formula(paste0(colnames(traitVec)[2], "~ scaled.dn.ds + scaledESM + scaledPlantCAD + PMS + compdup + queryCov + pacbio"))
    reducedFM = as.formula(paste0(colnames(traitVec)[2], "~ compdup + queryCov + pacbio"))
    
    model_full <- asreml(fixed = fullFM, ai.sing = F, data = filtered_data)
    # model_reduced1 <- asreml(fixed = fullFM,
    #                          ai.sing = F, data = filtered_data)
    filtered_data <- tmp
    model_reduced <- asreml(fixed = reducedFM, ai.sing = F, data = filtered_data)  
    
    totalVar = var(filtered_data[,ncol(filtered_data)])
    modelSummary <- summary(model_full)
    IC_full = infoCriteria(model_full,IClikelihood = "full")
    IC_reduced = infoCriteria(model_reduced,IClikelihood = "full")
    logLik_full <- IC_full$loglik
    logLik_reduced <- IC_reduced$loglik
    BIC_full <- IC_full$BIC
    BIC_reduced <- IC_reduced$BIC
    
    LRT_statistic <- -2 * (logLik_reduced - logLik_full)
    p_value <- pchisq(LRT_statistic, df = model_reduced$nedf-model_full$nedf, lower.tail = FALSE)
    # gene_explained <-  (model_reduced2$sigma2-model_full$sigma2)/totalVar
    phylo_explained <- modelSummary$varcomp$component[1]/sum(modelSummary$varcomp$component)
    
    modelWald <- wald.asreml(model_full)
    partialCoeff <- as.data.frame(model_full$coefficients$fixed)[c("PAV","scaled.dn.ds","scaledESM","scaledPlantCAD","PMS_1"),1]
    partialPval <- as.data.frame(modelWald)[c("PAV","scaled.dn.ds","scaledESM","scaledPlantCAD","PMS"),4]
    out = c(OG_id,BIC_full,BIC_reduced,LRT_statistic,p_value,phylo_explained,totalVar,
            partialCoeff,partialPval)
  }
  names(out) = c("OG","BIC_full","BIC_reduced","LR","p","phyloVAE","totalVar",
                 paste("partialCoeff",c("PAV","dNdS","ESM2","plantCad","PMS"),sep = "_"),
                 paste("partialP",c("PAV","dNdS","ESM2","plantCad","PMS"),sep = "_"))
  
  
  #     return(model_full)
  if(model_full$converge) return(out)
  else return(NULL)
}

PHYLOGWAS_ROOT <- Sys.getenv("PHYLOGWAS_ROOT", unset = "/workdir/sh2246/p_phyloGWAS")
MC_CORES <- as.integer(Sys.getenv("MC_CORES", unset = "35"))

args = commandArgs(trailingOnly=TRUE)
# arg 1: path to trait table; arg 2: path to master table; arg 3: path to phyloK; arg 4: path to tree;
# arg 5: name of response variable (e.g. envPC_1);
# arg 6: output directory; arg 7: permulation or not (TRUE or FALSE)
# arg 8: GLM or not (TRUE or FALSE)

# dir.trait = file.path(PHYLOGWAS_ROOT, "output/envData_707Poaceae_20250804.txt")
# dir.masterTab = file.path(PHYLOGWAS_ROOT, "output/masterDataTable_PAVFill_20250901_test.txt")
# dir.phyloK = file.path(PHYLOGWAS_ROOT, 'output/phyloK_728Poaceae_astral_20250407.txt')
# dir.tree = file.path(PHYLOGWAS_ROOT, "output/PoaceaeTree_angiosperm353_astral_filtered_20250407.nwk")
# responseVar = "envPC_1"
# dir.out = file.path(PHYLOGWAS_ROOT, "output/finalModels/")
# permulation = TRUE

dir.trait = args[1]
dir.masterTab = args[2]
dir.phyloK = args[3]
dir.tree = args[4]
responseVar = args[5]
dir.out = args[6]
permulation = as.logical(args[7])
GLM = as.logical(args[8])

if(!dir.exists(dir.out)) dir.create(dir.out)
dir.out = paste0(dir.out,"/",responseVar)
if(!dir.exists(dir.out)) dir.create(dir.out)
outDirRes = paste0(dir.out,"/","ASREML_res.txt")
outDirResGLM = paste0(dir.out,"/","ASREML_GLMres.txt")
outDirPerm = paste0(dir.out,"/","permulation")

print("start loading data")
envData = read.table(dir.trait,header = T)
ePCIdx = grep("PC",colnames(envData))
envDf = cbind(rownames(envData),envData[ePCIdx])
colnames(envDf)[1] = "assemblyID"
print("trait data loaded")

testData_merged = data.table::fread(dir.masterTab,header = T,data.table = F,nThread = 35)
testData_merged$PMS = as.factor(testData_merged$PMS)
testData_merged$PAV = as.factor(testData_merged$PAV)
testData_merged$assemblyID = as.factor(testData_merged$assemblyID)
testData_merged$pacbio = as.factor(testData_merged$pacbio)
testData_merged$groups2 = as.factor(testData_merged$groups2)
print("master table loaded")

phyloKMat <- read.table(dir.phyloK)
colnames(phyloKMat) <- rownames(phyloKMat)
print("phyloK loaded")

spTre = read.tree(dir.tree)
spTre = keep.tip(spTre,rownames(envDf))
spTre.rooted = root(spTre,"ASM1935983v1")
spTre.rooted$edge.length[is.na(spTre.rooted$edge.length)] = 0.1 # to fill tip branch length
spTre.rooted.resolved = ape::multi2di(spTre.rooted)
print("tree loaded")

# trait vector
trait_test = envDf[,c("assemblyID",responseVar)]
trait_testVec = trait_test[,2]
names(trait_testVec) = rownames(envDf)
print(head(trait_test))
print("modeling start")
OG_list = unique(testData_merged$OG)
# OG_list = unique(testRes_envPC1_v2$OG[testRes_envPC1_v2$p<0.001])
# set.seed(123)
# OG_list = sample(testRes_envPC1_v2$OG,1000)
# mappingFile = read.table("/workdir/sh2246/p_phyloGWAS_archived/output/OGToPv_mapping_v2.txt")
# mappingFile2 = read.table(file.path(PHYLOGWAS_ROOT, "data/OGToZm_mapping_v2.txt"))
# colnames(mappingFile) = c("PvID","OG")
# colnames(mappingFile2) = c("ZmID","OG")
# mappingFileMerged = merge(mappingFile,mappingFile2,by = "OG",all = T)

asreml.options(maxit = 40,verbose = F,trace=F, Cfixed = TRUE, gammaPar = FALSE)
test = try(fit_model_v4("OG0000327",data = testData_merged,Kmat = phyloKMat,traitVec = trait_test),silent = T)
print(test)
test = try(fit_model_v4_GLM("OG0000327",data = testData_merged,Kmat = phyloKMat,traitVec = trait_test),silent = T)
print(test)

# Fixed set of column names fit_model_v4()/fit_model_v4_GLM() always produce on a
# successful fit (see their names(out) <- ... assignments above). Used to validate and
# combine per-OG results below.
#
# This used to be inferred at runtime from a single hardcoded reference OG
# ("OG0000327")'s fit outcome (length(test), from the smoke test above) -- fragile,
# and confirmed to silently corrupt the whole output table when that one OG doesn't
# behave the way the full-707-taxon run happened to: on an Andropogoneae-only subset
# (483 taxa), OG0000327 fit differently and the length-based filter kept the FAILED
# fits instead of the successful ones, producing a garbled table of error-message text
# under auto-generated V1..V19613 column names (see revision-analyses.md).
#
# Every OG now gets a row in the output regardless of fit outcome: a failed fit (an R
# error caught by try(), or fit_model_v4's own internal non-convergence check
# returning NULL) gets NA for every result column and converged=FALSE, rather than
# being silently dropped (previous intended behavior) or corrupting the whole table
# (the actual bug). Failed fits are kept, not filtered out, per author's request --
# convergence failure is itself an informative outcome, not something to hide.
RESULT_COLS = c("OG","BIC_full","BIC_reduced","LR","p","phyloVAE","totalVar",
               paste("partialCoeff",c("PAV","dNdS","ESM2","plantCad","PMS"),sep = "_"),
               paste("partialP",c("PAV","dNdS","ESM2","plantCad","PMS"),sep = "_"))

combineResults = function(rawList, ogList) {
  converged = sapply(rawList, function(x) !is.null(x) && !inherits(x, "try-error") && length(x) == length(RESULT_COLS))
  out = matrix(NA_character_, nrow = length(ogList), ncol = length(RESULT_COLS), dimnames = list(NULL, RESULT_COLS))
  if (any(converged)) out[converged, ] = t(simplify2array(rawList[converged]))
  out[, "OG"] = ogList
  out = as.data.frame(out, stringsAsFactors = FALSE)
  for (col in setdiff(RESULT_COLS, "OG")) out[[col]] = as.numeric(out[[col]])
  out$converged = converged
  cat(sprintf("  %d / %d OGs converged (%.1f%%)\n", sum(converged), length(ogList), 100 * mean(converged)))
  out
}

if(!file.exists(outDirRes)){
  print("perOG modeling start")
  asreml.options(maxit = 40,verbose = F, Cfixed = TRUE, gammaPar = FALSE)
  testRes = mclapply(OG_list,function(x) try(fit_model_v4(x,testData_merged,phyloKMat,trait_test),
                                             silent = T),mc.cores = MC_CORES)
  testRes = combineResults(testRes, OG_list)
  # testRes = merge(testRes,mappingFileMerged,by = "OG",all.x = T)
  write.table(testRes,outDirRes,quote = F,sep = "\t",row.names = F)
  print("perOG modeling done and saved")
}

if(!file.exists(outDirResGLM)&GLM){
  print("perOG GLM modeling start")
  asreml.options(maxit = 40,verbose = F, Cfixed = TRUE, gammaPar = FALSE)
  testResGLM = mclapply(OG_list,function(x) try(fit_model_v4_GLM(x,testData_merged,phyloKMat,trait_test),
                                             silent = T),mc.cores = MC_CORES)
  testResGLM = combineResults(testResGLM, OG_list)
  # testResGLM = merge(testResGLM,mappingFileMerged,by = "OG",all.x = T)
  write.table(testResGLM,outDirResGLM,quote = F,sep = "\t",row.names = F)
  print("perOG GLM modeling done and saved")
}

if(permulation){
  testRes = read.table(outDirRes,header =T)
  topOG_list = unique(testRes$OG[!is.na(testRes$p) & testRes$p<0.001])
} else {
  # `} else stop` (bare reference to the stop function, not a call -- it never actually halts)
  # used to be here; harmless on its own, but the guard below used to read `permulation &
  # length(topOG_list)>1` with `&` (not `&&`), which evaluates both sides eagerly regardless
  # of permulation's value, so with permulation=FALSE and topOG_list never assigned this died
  # with "object 'topOG_list' not found" -- confirmed via a real run (revision-analyses.md's
  # no-permulation reruns). topOG_list is now defined either way, and the guard below uses
  # `&&` so it short-circuits and skips permulation cleanly when permulation is FALSE.
  topOG_list = character(0)
}

# emp_p is only trustworthy when enough of the 1000 permutation reps actually converged for
# that OG -- see the emp_p computation below. Overridable since 100 is a judgment call, not a
# hard statistical requirement.
MIN_PERM_CONVERGED <- as.integer(Sys.getenv("MIN_PERM_CONVERGED", unset = "100"))

# print(topOG_list)
if(permulation && length(topOG_list)>1){
  if(!dir.exists(outDirPerm)) dir.create(outDirPerm)
  
  print("permulation start")
  envPC_perm = list()
  set.seed(123)
  for (i in 1:1000){
    # rate estimates and simulation
    rate_test = ratematrix(spTre.rooted.resolved, trait_testVec)
    trait_sim = sim.char(spTre.rooted.resolved,rate_test,nsim = 1)
    #     trait_simulated = cbind(trait_simulated,trait_sim[,,1])
    envPC_perm[[i]] = data.frame("assemblyID"=rownames(envDf),"trait" = sort(trait_testVec)[rank(trait_sim[,,1])])
  }
  # print(head(envPC1_perm[[1]]))
  
  testRes_permList = list()
  for(i in 1:1000){
    testRes_permList[[i]] = mclapply(topOG_list,function(x) try(fit_model_v4(x,testData_merged,phyloKMat,envPC_perm[[i]]),silent = T),mc.cores = MC_CORES)
    testRes_permList[[i]] = combineResults(testRes_permList[[i]], topOG_list)
    # testRes_envPC1_permList[[i]] = merge(testRes_envPC1_permList[[i]],mappingFileMerged,by = "OG",all.x = T)
    write.table(testRes_permList[[i]],paste0(outDirPerm,"/ASREML_res_perm_",i,".txt"),
                quote = F,sep = "\t",row.names = F)
  }
  testRes_permMat = sapply(testRes_permList,function(x){
    tmp = x
    rownames(tmp) = tmp$OG
    tmp = tmp[topOG_list,]
    return(tmp[,5])
  })
  rownames(testRes_permMat) = topOG_list
  isTopOG = !is.na(testRes$p) & testRes$p < 0.001  # NA-safe version of testRes$p < 0.001 -- see topOG_list above
  testRes_permMat = cbind(testRes$p[isTopOG],
                                     testRes_permMat[topOG_list,])
  # emp_p = proportion of permutation reps more extreme than the real fit, among permutation
  # reps that actually converged for that OG -- NOT divided by the fixed 1000 permutation
  # count. Dividing by 1000 regardless of how many permutations converged (the previous
  # behavior) silently deflates emp_p -- i.e. makes OGs look more significant than warranted
  # -- whenever permulation convergence is low. That rarely mattered on the original full-
  # 707-taxon runs (~100% permulation convergence there), but on reduced-taxon-count subsets
  # (see revision-analyses.md) permulation convergence for topOG_list OGs can drop to single
  # digits out of 1000, since these OGs were specifically selected for having a strong signal
  # under the REAL trait -- permuting destroys that signal, and the remaining fit is harder to
  # converge with fewer taxa. n_perm_converged is recorded so low-reliability emp_p values are
  # visible, and emp_p is set to NA outright when fewer than MIN_PERM_CONVERGED permutations
  # converged: an emp_p computed from a handful of samples isn't trustworthy at the emp_p <=
  # 0.001 candidate-calling threshold this pipeline uses elsewhere.
  n_perm_converged = apply(testRes_permMat[,-1,drop = FALSE],1,function(x) sum(!is.na(x)))
  emp_p = apply(testRes_permMat,1,function(x) {
    n_conv = sum(!is.na(x[-1]))
    if (n_conv == 0) return(NA_real_)
    sum(x[-1] < x[1], na.rm = TRUE) / n_conv
  })
  emp_p[n_perm_converged < MIN_PERM_CONVERGED] = NA
  testRes$emp_p = NA
  testRes$emp_p[isTopOG] = emp_p
  testRes$n_perm_converged = NA_integer_
  testRes$n_perm_converged[isTopOG] = n_perm_converged
  write.table(testRes,paste0(dir.out,"/","ASREML_res_empPadded.txt"),quote = F,sep = "\t",row.names = F)
}

