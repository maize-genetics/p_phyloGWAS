rm(list = ls())
library(data.table)
library(dplyr)
df = fread("/workdir/sh2246/p_phyloGWAS/output/jonathanOG/matches.tsv",data.table = F)
df_filtered <- df %>%
  group_by(V1) %>%
  slice_min(order_by = V11, with_ties = FALSE) %>%
  ungroup()

df_filtered = df_filtered[df_filtered$V3>70,]
fwrite(df_filtered,"/workdir/sh2246/p_phyloGWAS/output/jonathanOG/OGTranslationTable.txt",sep = "\t")
