setwd("~/Documents/Etudes/Stage_projet_LOT/CRBE/données/analyses_di")


library(phyloseq)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(dplyr)

ps <- readRDS("donnees_nettoyees.rds")
sample_names(ps) <- gsub("\\.", "-", sub("^X", "", sample_names(ps)))

familles_a_enlever <- c("Moraceae", "Fabaceae") 
ps <- subset_samples(ps, !plant_family %in% familles_a_enlever)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

df_alpha <- estimate_richness(ps, measures = c("Observed", "Shannon")) |>
  cbind(sample_data(ps)) |>
  transform(
    expShannon = exp(Shannon),
    grp = interaction(organ, position, project)
  )

saveRDS(ps, "ps_final.rds")
saveRDS(df_alpha, "df_alpha.rds")
