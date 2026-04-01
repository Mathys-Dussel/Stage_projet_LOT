setwd("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses")

library(vegan)
library(phyloseq)
library(ggplot2)
library(patchwork)

ps <- readRDS("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/ps_final.rds")

sum(sample_data(ps)$Zone == "Unknown", na.rm = TRUE)

ps_hel <- transform_sample_counts(ps, function(x) sqrt(x / sum(x)))
ps_hel <- prune_samples(sample_sums(ps_hel) > 0, ps_hel)
otu_hel <- as(otu_table(ps_hel), "matrix")
metadata <- as(sample_data(ps_hel), "data.frame")

bray_dist <- vegdist(otu_hel, method = "bray")

# Global 
perm_global <- adonis2(bray_dist ~ project + organ + position, data = metadata, by = "margin", permutations = 99)
nmds_global <- ordinate(ps_hel, method = "NMDS", distance = "bray", k = 3)

p1 <- plot_ordination(ps_hel, nmds_global, color = "project", shape = "organ") +
  stat_ellipse(aes(group = project, fill = project), alpha = 0.1, geom = "polygon") +
  labs(title = "1. Structure Globale (Projet)", 
       subtitle = paste0("R2 Project: ", round(perm_global$R2[1], 3))) +
  theme_bw()

# Compartimentation 
ps_roots <- subset_samples(ps_hel, organ == "root")
nmds_roots <- ordinate(ps_roots, "NMDS", "bray")

p2 <- plot_ordination(ps_roots, nmds_roots, color="position") + 
  geom_line(aes(group = sample_id), alpha = 0.2) +
  labs(title = "2. Filtrage (Racines: Epi vs Endo)") +
  theme_minimal()

# Niche 
ps_epi <- subset_samples(ps_hel, position == "epiphyte")
ps_endo <- subset_samples(ps_hel, position == "endophyte")

p3_epi <- plot_ordination(ps_epi, ordinate(ps_epi, "NMDS", "bray"), color = "organ") +
  stat_ellipse() + labs(title = "3a. Organes (Épiphytes)") + theme_bw()

p3_endo <- plot_ordination(ps_endo, ordinate(ps_endo, "NMDS", "bray"), color = "organ") +
  stat_ellipse() + labs(title = "3b. Organes (Endophytes)") + theme_bw()

# Hôte 
perm_fam <- adonis2(bray_dist ~ plant_family, data = metadata, permutations = 99)

p4 <- plot_ordination(ps_hel, nmds_global, color = "plant_family") +
  stat_ellipse(alpha = 0.05) +
  labs(title = "4. Effet Hôte (Famille)",
       subtitle = paste0("R2 Famille: ", round(perm_fam$R2[1], 3))) +
  theme_bw() + theme(legend.position = "none")

p5 <- plot_ordination(ps_hel, nmds_global, color = "plant_genus") +
  stat_ellipse(alpha = 0.05) +
  labs(title = "4. Effet Hôte (Genre)",
       subtitle = paste0("R2 Genre: ", round(perm_fam$R2[1], 3))) +
  theme_bw() + theme(legend.position = "none")

disper_organ <- betadisper(bray_dist, metadata$organ)
plot(disper_organ)

(p1 | p2) / (p3_epi | p3_endo) / (p4 | p5) + plot_annotation(tag_levels = 'A')




perm_global_fam <- adonis2(bray_dist ~ project + organ + position + plant_family, data = metadata, by = "margin", permutations = 99)



results_list <- list(
  "Global (Project)" = perm_global_fam,
  "Global (Organ)"   = perm_global_fam,
  "Global (Position)" = perm_global_fam,
  "Famille Plante"   = perm_global_fam
)

extract_permanova <- function(perm_obj, row_name) {
  data.frame(
    Facteur = row_name,
    Df = perm_obj[row_name, "Df"],
    R2 = round(perm_obj[row_name, "R2"], 4),
    p_value = perm_obj[row_name, "Pr(>F)"]
  )
}

tab_permanova <- rbind(
  extract_permanova(perm_global_fam, "project"),
  extract_permanova(perm_global_fam, "organ"),
  extract_permanova(perm_global_fam, "position"),
  extract_permanova(perm_global_fam, "plant_family")
)

print(tab_permanova)

write.csv(tab_permanova, "donnees/Tableau_PERMANOVA_Synthese.csv", row.names = FALSE)





############### LOT01 #################

ps_lot01 <- subset_samples(ps_hel, project == "LOT1")
ps_lot01 <- prune_taxa(taxa_sums(ps_lot01) > 0, ps_lot01)

bray_lot01 <- vegdist(otu_table(ps_lot01), method = "bray")
perm_lot01 <- adonis2(bray_lot01 ~ organ + position + Zone + plant_family, data = as(sample_data(ps_lot01), "data.frame"), by = "margin", permutations = 999)

nmds_lot01 <- ordinate(ps_lot01, method = "NMDS", distance = "bray", k = 3)

p_lot1_organ <- plot_ordination(ps_lot01, nmds_lot01, color = "organ") +
  stat_ellipse() +
  labs(title = "LOT1 - Organe") +
  theme_bw()

p_lot1_position <- plot_ordination(ps_lot01, nmds_lot01, color = "position") +
  stat_ellipse() +
  labs(title = "LOT1 - Position") +
  theme_bw()

p_lot1_zone <- plot_ordination(ps_lot01, nmds_lot01, color = "Zone") +
  stat_ellipse() +
  labs(title = "LOT1 - Zone géographique") +
  theme_bw()

p_lot1_family <- plot_ordination(ps_lot01, nmds_lot01, color = "plant_family") +
  stat_ellipse() +
  labs(title = "LOT1 - Famille végétale") +
  theme_bw()

(p_lot1_organ | p_lot1_position) / (p_lot1_zone | p_lot1_family) + plot_annotation(title = "Ordination NMDS - LOT1")






############### LOT02 #################

ps_lot02 <- subset_samples(ps_hel, project == "LOT2")
ps_lot02 <- prune_taxa(taxa_sums(ps_lot02) > 0, ps_lot02)

bray_lot02 <- vegdist(otu_table(ps_lot02), method = "bray")
perm_lot02 <- adonis2(bray_lot02 ~ organ + position + Zone + plant_family, data = as(sample_data(ps_lot02), "data.frame"), by = "margin", permutations = 99)

nmds_lot02 <- ordinate(ps_lot02, method = "NMDS", distance = "bray", k = 3)

p_lot2_organ <- plot_ordination(ps_lot02, nmds_lot02, color = "organ") +
  stat_ellipse() +
  labs(title = "LOT2 - Organe") +
  theme_bw()

p_lot2_position <- plot_ordination(ps_lot02, nmds_lot02, color = "position") +
  stat_ellipse() +
  labs(title = "LOT2 - Position") +
  theme_bw()

p_lot2_zone <- plot_ordination(ps_lot02, nmds_lot02, color = "Zone") +
  stat_ellipse() +
  labs(title = "LOT2 - Zone géographique") +
  theme_bw()

p_lot2_family <- plot_ordination(ps_lot02, nmds_lot02, color = "plant_family") +
  stat_ellipse() +
  labs(title = "LOT2 - Famille végétale") +
  theme_bw()

(p_lot2_organ | p_lot2_position) / (p_lot2_zone | p_lot2_family) + plot_annotation(title = "Ordination NMDS - LOT2")


############### LOT03 #################

ps_lot03 <- subset_samples(ps_hel, project == "LOT3")
ps_lot03 <- prune_taxa(taxa_sums(ps_lot03) > 0, ps_lot03)

bray_lot03 <- vegdist(otu_table(ps_lot03), method = "bray")
perm_lot03 <- adonis2(bray_lot03 ~ organ + position + Zone + plant_family, data = as(sample_data(ps_lot03), "data.frame"), by = "margin", permutations = 99)

nmds_lot03 <- ordinate(ps_lot03, method = "NMDS", distance = "bray", k = 3)

p_lot3_organ <- plot_ordination(ps_lot03, nmds_lot03, color = "organ") +
  stat_ellipse() +
  labs(title = "LOT3 - Organe") +
  theme_bw()

p_lot3_position <- plot_ordination(ps_lot03, nmds_lot03, color = "position") +
  stat_ellipse() +
  labs(title = "LOT3 - Position") +
  theme_bw()

p_lot3_zone <- plot_ordination(ps_lot03, nmds_lot03, color = "Zone") +
  stat_ellipse() +
  labs(title = "LOT3 - Zone géographique") +
  theme_bw()

p_lot3_family <- plot_ordination(ps_lot03, nmds_lot03, color = "plant_family") +
  stat_ellipse() +
  labs(title = "LOT3 - Famille végétale") +
  theme_bw()

(p_lot3_organ | p_lot3_position) / (p_lot3_zone | p_lot3_family) + plot_annotation(title = "Ordination NMDS - LOT3")



############### Racine #################

unique(sample_data(ps_hel)$organ)

ps_racine <- subset_samples(ps_hel, organ == "root")
ps_racine <- prune_taxa(taxa_sums(ps_racine) > 0, ps_racine)

bray_racine <- vegdist(otu_table(ps_racine), method = "bray")
perm_racine <- adonis2(bray_racine ~ project + position + Zone + plant_family, data = as(sample_data(ps_racine), "data.frame"), by = "margin", permutations = 99)
nmds_racine <- ordinate(ps_racine, method = "NMDS", distance = "bray", k = 3)

p_racine_organ <- plot_ordination(ps_racine, nmds_racine, color = "organ") +
  stat_ellipse() +
  labs(title = "Racine - Organe") +
  theme_bw()

p_racine_position <- plot_ordination(ps_racine, nmds_racine, color = "position") +
  stat_ellipse() +
  labs(title = "Racine - Position") +
  theme_bw()

p_racine_zone <- plot_ordination(ps_racine, nmds_racine, color = "Zone") +
  stat_ellipse() +
  labs(title = "Racine - Zone géographique") +
  theme_bw()

p_racine_family <- plot_ordination(ps_racine, nmds_racine, color = "plant_family") +
  stat_ellipse() +
  labs(title = "Racine - Famille végétale") +
  theme_bw()

(p_racine_organ | p_racine_position) / (p_racine_zone | p_racine_family) + plot_annotation(title = "Ordination NMDS - Racine")