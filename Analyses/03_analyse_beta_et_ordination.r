setwd("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses")

library(vegan)
library(phyloseq)
library(ggplot2)
library(patchwork)

ps <- readRDS("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/donnees_nettoyees.rds")

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







results_list <- list(
  "Global (Project)" = perm_global,
  "Global (Organ)"   = perm_global,
  "Global (Position)" = perm_global,
  "Famille Plante"   = perm_fam
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
  extract_permanova(perm_global, "project"),
  extract_permanova(perm_global, "organ"),
  extract_permanova(perm_global, "position"),
  extract_permanova(perm_fam, "plant_family")
)

print(tab_permanova)

write.csv(tab_permanova, "donnees/Tableau_PERMANOVA_Synthese.csv", row.names = FALSE)
