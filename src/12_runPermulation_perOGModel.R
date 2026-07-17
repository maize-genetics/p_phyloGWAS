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
  asreml.options(verbose = FALSE)
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
  asreml.options(verbose = FALSE)
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

args = commandArgs(trailingOnly=TRUE)
# arg 1: path to trait table; arg 2: path to master table; arg 3: path to phyloK; arg 4: path to tree;
# arg 5: name of response variable (e.g. envPC_1);
# arg 6: output directory; arg 7: permulation or not (TRUE or FALSE)
# arg 8: GLM or not (TRUE or FALSE)

# dir.trait = "/workdir/sh2246/p_phyloGWAS/output/envData_707Poaceae_20250804.txt"
# dir.masterTab = "/workdir/sh2246/p_phyloGWAS/output/masterDataTable_PAVFill_20250901_test.txt"
# dir.phyloK = '/workdir/sh2246/p_phyloGWAS/output/phyloK_728Poaceae_astral_20250407.txt'
# dir.tree = "/workdir/sh2246/p_phyloGWAS/output/PoaceaeTree_angiosperm353_astral_filtered_20250407.nwk"
# responseVar = "envPC_1"
# dir.out = "/workdir/sh2246/p_phyloGWAS/output/finalModels/"
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
# mappingFile2 = read.table("/workdir/sh2246/p_phyloGWAS/data/OGToZm_mapping_v2.txt")
# colnames(mappingFile) = c("PvID","OG")
# colnames(mappingFile2) = c("ZmID","OG")
# mappingFileMerged = merge(mappingFile,mappingFile2,by = "OG",all = T)

asreml.options(maxit = 40,verbose = F,trace=F)
test = try(fit_model_v4("OG0000327",data = testData_merged,Kmat = phyloKMat,traitVec = trait_test),silent = T)
test = try(fit_model_v4("OG0000327",data = testData_merged,Kmat = phyloKMat,traitVec = trait_test),silent = T)
print(test)
test = try(fit_model_v4_GLM("OG0000327",data = testData_merged,Kmat = phyloKMat,traitVec = trait_test),silent = T)
print(test)

if(!file.exists(outDirRes)){
  print("perOG modeling start")
  asreml.options(maxit = 40,verbose = F)
  testRes = mclapply(OG_list,function(x) try(fit_model_v4(x,testData_merged,phyloKMat,trait_test),
                                             silent = T),mc.cores = 35)
  testRes = t(simplify2array(testRes[sapply(testRes,length)==length(test)]))
  testRes = as.data.frame(testRes)
  # testRes = merge(testRes,mappingFileMerged,by = "OG",all.x = T)
  write.table(testRes,outDirRes,quote = F,sep = "\t",row.names = F)
  print("perOG modeling done and saved")
}

if(!file.exists(outDirResGLM)&GLM){
  print("perOG GLM modeling start")
  asreml.options(maxit = 40,verbose = F)
  testRes = mclapply(OG_list,function(x) try(fit_model_v4_GLM(x,testData_merged,phyloKMat,trait_test),
                                             silent = T),mc.cores = 35)
  testRes = t(simplify2array(testRes[sapply(testRes,length)==length(test)]))
  testRes = as.data.frame(testRes)
  # testRes = merge(testRes,mappingFileMerged,by = "OG",all.x = T)
  write.table(testRes,outDirResGLM,quote = F,sep = "\t",row.names = F)
  print("perOG GLM modeling done and saved")
}

if(permulation){
  testRes = read.table(outDirRes,header =T)
  topOG_list = unique(testRes$OG[testRes$p<0.001])
} else stop

# print(topOG_list)
if(permulation&length(topOG_list)>1){
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
    testRes_permList[[i]] = mclapply(topOG_list,function(x) try(fit_model_v4(x,testData_merged,phyloKMat,envPC_perm[[i]]),silent = T),mc.cores = 35)
    testRes_permList[[i]] = t(simplify2array(testRes_permList[[i]][sapply(testRes_permList[[i]],length)==length(test)]))
    testRes_permList[[i]] = as.data.frame(testRes_permList[[i]])
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
  testRes_permMat = cbind(testRes$p[testRes$p < 0.001],
                                     testRes_permMat[topOG_list,])
  emp_p = apply(testRes_permMat,1,function(x) sum(x[-1]<x[1],na.rm = T))/ncol(testRes_permMat[,-1])
  testRes$emp_p = NA
  testRes$emp_p[testRes$p < 0.001] = emp_p
  write.table(testRes,paste0(dir.out,"/","ASREML_res_empPadded.txt"),quote = F,sep = "\t",row.names = F)
}

