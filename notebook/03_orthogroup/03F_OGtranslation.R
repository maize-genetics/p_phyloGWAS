rm(list = ls())
library(data.table)
library(dplyr)
PHYLOGWAS_ROOT <- Sys.getenv("PHYLOGWAS_ROOT", unset = "/workdir/sh2246/p_phyloGWAS")
df = fread(file.path(PHYLOGWAS_ROOT, "output/jonathanOG/matches.tsv"),data.table = F)
df_filtered <- df %>%
  group_by(V1) %>%
  slice_min(order_by = V11, with_ties = FALSE) %>%
  ungroup()

df_filtered = df_filtered[df_filtered$V3>70,]
fwrite(df_filtered,file.path(PHYLOGWAS_ROOT, "output/jonathanOG/OGTranslationTable.txt"),sep = "\t")
