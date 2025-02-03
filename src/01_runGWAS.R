library(asreml)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(parallel)
install.packages("qqman")
library(qqman)

allOG_fast <- read.table('data/allOGs_withHeaders_dnds.txt', header = TRUE)
metadata <- read.csv('data/Poaceae_metadata_2024.06.12.csv')
phyloK <- read.table('data/phyloK_800Poaceae_astral.txt')
tabasco <- read.csv('data/tabasco_Poaceae800.csv')

allOG_fast <- allOG_fast %>% separate_wider_delim(seq2,":", names = c("OG","assemblyID","seqNumber"))
allOG_fast <- merge(allOG_fast, metadata, by = "assemblyID")
allOG_fast <- merge(allOG_fast, tabasco, by = "assemblyID")
allOG_fast$OG <- gsub("_R_","", allOG_fast$OG)


forFiltering <- allOG_fast %>% group_by(OG) %>% count(assemblyID) %>% filter(n > 12) %>% summarise(unique(OG))
allOG_fast <- allOG_fast %>% filter(!(OG %in% forFiltering$OG))

testData <-  allOG_fast %>% filter(lifeHistory == "annual" | lifeHistory == "perennial") %>% filter(!is.na(ds)) %>% filter(!is.na(dn)) %>% filter(!is.na(dn.ds)) %>% filter(!is.na(queryCov))

testData$assemblyID <- as.factor(testData$assemblyID)
testData$lifeHistory <- as.factor(testData$lifeHistory)

testData$lifeHistoryNumeric <- as.numeric(testData$lifeHistory)

testData <- testData %>% filter(assemblyID %in% rownames(phyloK))

colnames(phyloK) <- rownames(phyloK)

testData$compdup <- testData$complete + testData$duplicated

OG_list <- unique(testData$OG)

fit_model <- function(OG_id, data, Kmat) {
  library(asreml)
  asreml.options(verbose = FALSE)
  filtered_data <- subset(data, OG == OG_id)
  #filtered_data <- filtered_data %>% group_by(assemblyID) %>% slice_max(queryCov, n = 1, with_ties = FALSE) %>% ungroup()
  phyloK_OG <- Kmat[rownames(Kmat) %in% filtered_data$assemblyID,]
  phyloK_OG <- phyloK_OG[,colnames(phyloK_OG) %in% filtered_data$assemblyID]
  model <- asreml(fixed = lifeHistoryNumeric ~ dn.ds, random = ~ vm(assemblyID, phyloK_OG) + compdup, ai.sing = TRUE, data = filtered_data)
  #model <- asreml(fixed = lifeHistoryNumeric ~ dn.ds, random = ~compdup, ai.sing = TRUE, data = filtered_data)
  modelSummary <- summary(model)
  loglik <- modelSummary$loglik
  tabasco_explained <- modelSummary$varcomp$component[1]
  phylo_explained <- modelSummary$varcomp$component[2]
  modelWald <- wald.asreml(model)
  dnds_pval <- modelWald[2,4]
  dnds_explained <- modelWald[2,2]/sum(modelWald[,2])
  return(c(OG_id, dnds_pval, loglik, tabasco_explained, phylo_explained, dnds_explained))
}


num_cores <- 6
cl <- makeCluster(num_cores)

# Export necessary objects and functions to the cluster
clusterExport(cl, list("fit_model", "testData", "OG_list", "phyloK"))

# Run the model fitting in parallel
results <- parLapply(cl, OG_list, function(id) try(fit_model(id, testData, phyloK), silent = TRUE))

# Stop the cluster after execution
stopCluster(cl)

# Combine the results into a matrix
results_matrix <- do.call(rbind, lapply(results, function(x) if(!inherits(x, "try-error")) x else c(NA, NA, NA, NA, NA,NA)))
colnames(results_matrix) <- c("OG_id", "dnds_pval", "loglik", "tabasco_explained", "phylo_explained", "dnds_explained")

# Convert to data frame for easier handling if needed
results_df <- as.data.frame(results_matrix, stringsAsFactors = FALSE)
results_df$dnds_pval <- as.numeric(results_df$dnds_pval)

write.table(results_df, "full_dnds_results.txt", sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

plot.new()
png(filename = "manhattan.png", width = 20, height = 20, units = "cm", res = 600, pointsize = 6)
results_df %>% ggplot(aes(x = OG_id, y = -log10(dnds_pval))) + geom_point()
dev.off()

plot.new()
png(filename = "histogram.png", width = 20, height = 20, units = "cm", res = 600, pointsize = 6)
results_df %>% ggplot(aes(x = -log10(dnds_pval))) + geom_histogram(bins = 50)
dev.off()

plot.new()
png(filename = "qq.png", width = 20, height = 20, units = "cm", res = 600, pointsize = 6)
result_qq <- qq(results_df$dnds_pval, pch = 1)
dev.off()

