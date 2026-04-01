setwd("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses")

library(phyloseq)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(dplyr)

ps <- readRDS("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/donnees_nettoyees.rds")
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


sample_data(ps)$Zone <- "Unknown"
sample_data(ps)$Zone[substr(sample_data(ps)$LOT_sampleID, 9, 9) == "1"] <- "Bottom"
sample_data(ps)$Zone[substr(sample_data(ps)$LOT_sampleID, 9, 9) %in% c("2", "3")] <- "top_trunk"
sample_data(ps)$Zone[substr(sample_data(ps)$LOT_sampleID, 9, 9) == "4"] <- "Heart"
sample_data(ps)$Zone[substr(sample_data(ps)$LOT_sampleID, 9, 9) %in% c("5", "6")] <- "periph_canopy"



saveRDS(ps, "~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/ps_final.rds")
saveRDS(df_alpha, "~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/df_alpha.rds")
