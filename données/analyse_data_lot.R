setwd("~/Documents/Etudes/Stage_projet_LOT/CRBE/données")

library(phyloseq)

phyloseq_lot <- readRDS("phyloseq_lot.rds")

unique(sample_data(phyloseq_lot)$project)
unique(sample_data(phyloseq_lot)$organ)
phyloseq_lot2 <- subset_samples(phyloseq_lot, organ != "flower")
unique(sample_data(phyloseq_lot2)$organ)

unique(sample_data(phyloseq_lot2)$position)
phyloseq_lot2_clean <- subset_samples(phyloseq_lot2, position %in% c("endophyte", "epiphyte"))
unique(sample_data(phyloseq_lot2_clean)$position)

unique(sample_data(phyloseq_lot2_clean)$plant_family)
unique(sample_data(phyloseq_lot2_clean)$plant_genus)


plant_df <- data.frame(
  sample = sample_names(phyloseq_lot2_clean),
  plant_genus = sample_data(phyloseq_lot2_clean)$plant_genus,
  plant_family = sample_data(phyloseq_lot2_clean)$plant_family,
  plant_species = sample_data(phyloseq_lot2_clean)$plant_species,
  row.names = NULL
)


tax_df <- as.data.frame(tax_table(phyloseq_lot2_clean))
get_taxa_unique(phyloseq_lot2_clean, "gbr268_Kingdom")
get_taxa_unique(phyloseq_lot2_clean, "gbr268_Phylum")
na.omit(unique(as.data.frame(tax_table(phyloseq_lot2_clean))$gbr268_Phylum))


library(ggplot2)

saveRDS(phyloseq_lot2_clean, file = "donnees_nettoyees.rds")




alpha_meas <- c("Observed", "Chao1", "Shannon", "Simpson")
p <- plot_richness(phyloseq_lot2_clean, x = "project", measures = alpha_meas, color = "organ")
p + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

