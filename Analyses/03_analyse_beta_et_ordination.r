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
  labs(title = "2. Filtrage (Racines: Epi vs Endo)", subtitle = paste0("R2 Position: ", round(perm_global$R2[3], 3))) +
  theme_minimal()

# Niche 
ps_epi <- subset_samples(ps_hel, position == "epiphyte")
ps_endo <- subset_samples(ps_hel, position == "endophyte")

p3_epi <- plot_ordination(ps_epi, ordinate(ps_epi, "NMDS", "bray"), color = "organ") +
  stat_ellipse() + labs(title = "3a. Organes (Épiphytes)", subtitle = paste0("R2 Organ: ", round(perm_global$R2[2], 3))) + theme_bw()

p3_endo <- plot_ordination(ps_endo, ordinate(ps_endo, "NMDS", "bray"), color = "organ") +
  stat_ellipse() + labs(title = "3b. Organes (Endophytes)", subtitle = paste0("R2 Organ: ", round(perm_global$R2[2], 3))) + theme_bw()

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

# Remove samples with missing values for PERMANOVA variables
meta_lot01 <- as(sample_data(ps_lot01), "data.frame")
valid_samples <- complete.cases(meta_lot01[, c("organ", "position", "Zone", "plant_family", "Morphologie", "Capacité_Nid")])
ps_lot01 <- prune_samples(valid_samples, ps_lot01)

bray_lot01 <- vegdist(otu_table(ps_lot01), method = "bray")
perm_lot01 <- adonis2(bray_lot01 ~ organ + position + Zone + plant_family + Morphologie + Capacité_Nid, data = as(sample_data(ps_lot01), "data.frame"), by = "margin", permutations = 999)

nmds_lot01 <- ordinate(ps_lot01, method = "NMDS", distance = "bray", k = 3)

p_lot1_organ <- plot_ordination(ps_lot01, nmds_lot01, color = "organ") +
  stat_ellipse() +
  labs(title = "LOT1 - Organe", subtitle = paste0("R2 = ", round(perm_lot01["organ", "R2"], 3), ", p = ", perm_lot01["organ", "Pr(>F)"])) +
  theme_bw()

p_lot1_position <- plot_ordination(ps_lot01, nmds_lot01, color = "position") +
  stat_ellipse() +
  labs(title = "LOT1 - Position", subtitle = paste0("R2 = ", round(perm_lot01["position", "R2"], 3), ", p = ", perm_lot01["position", "Pr(>F)"])) +
  theme_bw()

p_lot1_zone <- plot_ordination(ps_lot01, nmds_lot01, color = "Zone") +
  stat_ellipse() +
  labs(title = "LOT1 - Zone", subtitle = paste0("R2 = ", round(perm_lot01["Zone", "R2"], 3), ", p = ", perm_lot01["Zone", "Pr(>F)"])) +
  theme_bw()

p_lot1_family <- plot_ordination(ps_lot01, nmds_lot01, color = "plant_family") +
  stat_ellipse() +
  labs(title = "LOT1 - Famille végétale", subtitle = paste0("R2 = ", round(perm_lot01["plant_family", "R2"], 3), ", p = ", perm_lot01["plant_family", "Pr(>F)"])) +
  theme_bw()

(p_lot1_organ | p_lot1_position) / (p_lot1_zone | p_lot1_family) + plot_annotation(title = "Ordination NMDS - LOT1")




p_org_pos_1=plot_ordination(ps_lot01, nmds_lot01, color = "organ", shape = "position") +
  geom_point(size = 3) +
  scale_shape_manual(values = c("epiphyte" = 1, "endophyte" = 16)) +
  stat_ellipse(aes(group = organ)) +
  labs(title = "LOT1 - Organe", subtitle = paste0("R2 = ", round(perm_lot01["organ", "R2"], 3), ", p = ", perm_lot01["organ", "Pr(>F)"])) +
  theme_bw()


############### LOT02 #################

ps_lot02 <- subset_samples(ps_hel, project == "LOT2")
ps_lot02 <- prune_taxa(taxa_sums(ps_lot02) > 0, ps_lot02)

meta_lot02 <- as(sample_data(ps_lot02), "data.frame")
valid_samples <- complete.cases(meta_lot02[, c("organ", "position", "Zone", "plant_family", "Morphologie", "Capacité_Nid")])
ps_lot02 <- prune_samples(valid_samples, ps_lot02)


bray_lot02 <- vegdist(otu_table(ps_lot02), method = "bray")
# Only keep variables with 2 or more levels
vars_lot02 <- c("organ", "position", "Zone", "plant_family", "Morphologie", "Capacité_Nid")
valid_vars_lot02 <- vars_lot02[sapply(meta_lot02[valid_samples, vars_lot02], function(x) length(unique(x)) > 1)]
formula_lot02 <- as.formula(paste("bray_lot02 ~", paste(valid_vars_lot02, collapse = " + ")))

perm_lot02 <- adonis2(formula_lot02, data = as(sample_data(ps_lot02), "data.frame"), by = "margin", permutations = 99)

nmds_lot02 <- ordinate(ps_lot02, method = "NMDS", distance = "bray", k = 3)

p_lot2_organ <- plot_ordination(ps_lot02, nmds_lot02, color = "organ") +
  stat_ellipse() +
  labs(title = "LOT2 - Organe", subtitle = paste0("R2 = ", round(perm_lot02["organ", "R2"], 3), ", p = ", perm_lot02["organ", "Pr(>F)"])) +
  theme_bw()

p_lot2_position <- plot_ordination(ps_lot02, nmds_lot02, color = "position") +
  stat_ellipse() +
  labs(title = "LOT2 - Position", subtitle = paste0("R2 = ", round(perm_lot02["position", "R2"], 3), ", p = ", perm_lot02["position", "Pr(>F)"])) +
  theme_bw()

p_lot2_zone <- plot_ordination(ps_lot02, nmds_lot02, color = "Zone") +
  stat_ellipse() +
  labs(title = "LOT2 - Zone", subtitle = paste0("R2 = ", round(perm_lot02["Zone", "R2"], 3), ", p = ", perm_lot02["Zone", "Pr(>F)"])) +
  theme_bw()

p_lot2_family <- plot_ordination(ps_lot02, nmds_lot02, color = "plant_family") +
  stat_ellipse() +
  labs(title = "LOT2 - Famille végétale", subtitle = paste0("R2 = ", round(perm_lot02["plant_family", "R2"], 3), ", p = ", perm_lot02["plant_family", "Pr(>F)"])) +
  theme_bw()

(p_lot2_organ | p_lot2_position) / (p_lot2_zone | p_lot2_family) + plot_annotation(title = "Ordination NMDS - LOT2")




p_org_pos_2=plot_ordination(ps_lot02, nmds_lot02, color = "organ", shape = "position") +
  geom_point(size = 3) +
  scale_shape_manual(values = c("epiphyte" = 1, "endophyte" = 16)) +
  stat_ellipse(aes(group = organ)) +
  labs(title = "LOT2 - Organe", subtitle = paste0("R2 = ", round(perm_lot02["organ", "R2"], 3), ", p = ", perm_lot02["organ", "Pr(>F)"])) +
  theme_bw()




############### LOT03 #################

ps_lot03 <- subset_samples(ps_hel, project == "LOT3")
ps_lot03 <- prune_taxa(taxa_sums(ps_lot03) > 0, ps_lot03)

bray_lot03 <- vegdist(otu_table(ps_lot03), method = "bray")
perm_lot03 <- adonis2(bray_lot03 ~ organ + position + Zone + plant_family, data = as(sample_data(ps_lot03), "data.frame"), by = "margin", permutations = 99)

nmds_lot03 <- ordinate(ps_lot03, method = "NMDS", distance = "bray", k = 3)

p_lot3_organ <- plot_ordination(ps_lot03, nmds_lot03, color = "organ") +
  stat_ellipse() +
  labs(title = "LOT3 - Organe", subtitle = paste0("R2 = ", round(perm_lot03["organ", "R2"], 3), ", p = ", perm_lot03["organ", "Pr(>F)"])) +
  theme_bw()

p_lot3_position <- plot_ordination(ps_lot03, nmds_lot03, color = "position") +
  stat_ellipse() +
  labs(title = "LOT3 - Position", subtitle = paste0("R2 = ", round(perm_lot03["position", "R2"], 3), ", p = ", perm_lot03["position", "Pr(>F)"])) +
  theme_bw()

p_lot3_zone <- plot_ordination(ps_lot03, nmds_lot03, color = "Zone") +
  stat_ellipse() +
  labs(title = "LOT3 - Zone", subtitle = paste0("R2 = ", round(perm_lot03["Zone", "R2"], 3), ", p = ", perm_lot03["Zone", "Pr(>F)"])) +
  theme_bw()

p_lot3_family <- plot_ordination(ps_lot03, nmds_lot03, color = "plant_family") +
  stat_ellipse() +
  labs(title = "LOT3 - Famille végétale", subtitle = paste0("R2 = ", round(perm_lot03["plant_family", "R2"], 3), ", p = ", perm_lot03["plant_family", "Pr(>F)"])) +
  theme_bw()

(p_lot3_organ | p_lot3_position) / (p_lot3_zone | p_lot3_family) + plot_annotation(title = "Ordination NMDS - LOT3")



p_org_pos_3=plot_ordination(ps_lot03, nmds_lot03, color = "organ", shape = "position") +
  geom_point(size = 3) +
  scale_shape_manual(values = c("epiphyte" = 1, "endophyte" = 16)) +
  stat_ellipse(aes(group = organ)) +
  labs(title = "LOT3 - Organe", subtitle = paste0("R2 = ", round(perm_lot03["organ", "R2"], 3), ", p = ", perm_lot03["organ", "Pr(>F)"])) +
  theme_bw()



plot_ordination(ps_lot03, nmds_lot03, color = "organ", shape = "position") +
  geom_point(size = 3) +
  scale_shape_manual(values = c("epiphyte" = 1, "endophyte" = 16)) +
  stat_ellipse(aes(group = organ)) +
  labs(title = "LOT3 - Organe", subtitle = paste0("R2 = ", round(perm_lot03["organ", "R2"], 3), ", p = ", perm_lot03["organ", "Pr(>F)"])) +
  theme_bw()


vars_pop <- c("project", "organ", "position")
form_pop <- as.formula(paste("bray_dist ~", paste(vars_pop, collapse = " + ")))
perm_pop <- adonis2(form_pop, data = metadata, by = "margin", permutations = 999)

tab_pop <- data.frame(
  Facteur = vars_pop,
  Df = perm_pop[vars_pop, "Df"],
  R2 = round(perm_pop[vars_pop, "R2"], 4),
  p_value = perm_pop[vars_pop, "Pr(>F)"]
)
rownames(tab_pop) <- NULL

print(perm_pop)
print(tab_pop)

write.csv(tab_pop, "donnees/Tableau_PERMANOVA_project_organ_position.csv", row.names = FALSE)



p1 <- plot_ordination(ps_hel, nmds_global, color = "project", shape = "organ") +
  stat_ellipse(aes(group = project, fill = project), alpha = 0.1, geom = "polygon") +
  labs(title = "Structure Globale (Projet)", 
       subtitle = paste0("R2 = ", tab_pop$R2[1], ", p = ", tab_pop$p_value[1])) +
  theme_bw()

p_org_pos_1=plot_ordination(ps_lot01, nmds_lot01, color = "organ", shape = "position") +
  geom_point(size = 3) +
  scale_shape_manual(values = c("epiphyte" = 1, "endophyte" = 16)) +
  stat_ellipse(aes(group = organ)) +
  labs(title = "LOT1 - Organe x position", subtitle = paste0("R2 = ", tab_pop$R2[2], ", p = ", tab_pop$p_value[2])) +
  theme_bw()

p_org_pos_2=plot_ordination(ps_lot02, nmds_lot02, color = "organ", shape = "position") +
  geom_point(size = 3) +
  scale_shape_manual(values = c("epiphyte" = 1, "endophyte" = 16)) +
  stat_ellipse(aes(group = organ)) +
  labs(title = "LOT2 - Organe x position", subtitle = paste0("R2 = ", tab_pop$R2[3], ", p = ", tab_pop$p_value[3])) +
  theme_bw()

p_org_pos_3=plot_ordination(ps_lot03, nmds_lot03, color = "organ", shape = "position") +
  geom_point(size = 3) +
  scale_shape_manual(values = c("epiphyte" = 1, "endophyte" = 16)) +
  stat_ellipse(aes(group = organ)) +
  labs(title = "LOT3 - Organe x position", subtitle = paste0("R2 = ", tab_pop$R2[4], ", p = ", tab_pop$p_value[4])) +
  theme_bw()

(p1 | p_org_pos_1) / (p_org_pos_2 | p_org_pos_3) 


###########################


ps_lot1_endo <- subset_samples(ps_hel, project == "LOT1" & position == "endophyte")
ps_lot1_endo <- prune_taxa(taxa_sums(ps_lot1_endo) > 0, ps_lot1_endo)

meta_lot1_endo <- as(sample_data(ps_lot1_endo), "data.frame")
valid_samples_endo <- complete.cases(meta_lot1_endo[, "plant_genus", drop = FALSE])
ps_lot1_endo <- prune_samples(valid_samples_endo, ps_lot1_endo)

# Keep the 5 most frequent plant genera
top_5_genera <- names(sort(table(sample_data(ps_lot1_endo)$plant_genus), decreasing = TRUE)[1:5])
ps_lot1_endo <- subset_samples(ps_lot1_endo, plant_genus %in% top_5_genera)
ps_lot1_endo <- prune_taxa(taxa_sums(ps_lot1_endo) > 0, ps_lot1_endo)

bray_lot1_endo <- vegdist(otu_table(ps_lot1_endo), method = "bray")
perm_lot1_endo <- adonis2(bray_lot1_endo ~ plant_genus, data = as(sample_data(ps_lot1_endo), "data.frame"), permutations = 999)

print("Effet plant_genus chez les endophytes du LOT1 (Top 10 genres) :")
print(perm_lot1_endo)

nmds_lot1_endo <- ordinate(ps_lot1_endo, method = "NMDS", distance = "bray", k = 3)
p_lot1_endo_gen <- plot_ordination(ps_lot1_endo, nmds_lot1_endo, color = "plant_genus") +
  stat_ellipse() +
  labs(title = "LOT1 - Endophytes - Genre végétal (Top 5)", 
       subtitle = paste0("R2 = ", round(perm_lot1_endo[1, "R2"], 3), 
                         ", p = ", perm_lot1_endo[1, "Pr(>F)"])) +
  theme_bw()

print(p_lot1_endo_gen)






ps_lot1_epi <- subset_samples(ps_hel, project == "LOT1" & position == "epiphyte")
ps_lot1_epi <- prune_taxa(taxa_sums(ps_lot1_epi) > 0, ps_lot1_epi)

meta_lot1_epi <- as(sample_data(ps_lot1_epi), "data.frame")
valid_samples_epi <- complete.cases(meta_lot1_epi[, "plant_genus", drop = FALSE])
ps_lot1_epi <- prune_samples(valid_samples_epi, ps_lot1_epi)

# Keep the 5 most frequent plant genera
top_5_genera <- names(sort(table(sample_data(ps_lot1_epi)$plant_genus), decreasing = TRUE)[1:5])
ps_lot1_epi <- subset_samples(ps_lot1_epi, plant_genus %in% top_5_genera)
ps_lot1_epi <- prune_taxa(taxa_sums(ps_lot1_epi) > 0, ps_lot1_epi)

bray_lot1_epi <- vegdist(otu_table(ps_lot1_epi), method = "bray")
perm_lot1_epi <- adonis2(bray_lot1_epi ~ plant_genus, data = as(sample_data(ps_lot1_epi), "data.frame"), permutations = 999)

print("Effet plant_genus chez les épiphytes du LOT1 (Top 5 genres) :")
print(perm_lot1_epi)

nmds_lot1_epi <- ordinate(ps_lot1_epi, method = "NMDS", distance = "bray", k = 3)
p_lot1_epi_gen <- plot_ordination(ps_lot1_epi, nmds_lot1_epi, color = "plant_genus") +
  stat_ellipse() +
  labs(title = "LOT1 - Epiphytes - Genre végétal (Top 5)", 
       subtitle = paste0("R2 = ", round(perm_lot1_epi[1, "R2"], 3), 
                         ", p = ", perm_lot1_epi[1, "Pr(>F)"])) +
  theme_bw()

print(p_lot1_epi_gen)











ps_morph <- subset_samples(ps_hel, !is.na(plant_family) & !is.na(Morphologie))
ps_morph <- prune_taxa(taxa_sums(ps_morph) > 0, ps_morph)

bray_morph <- vegdist(otu_table(ps_morph), method = "bray")
perm_morph <- adonis2(bray_morph ~ plant_family + Morphologie, data = as(sample_data(ps_morph), "data.frame"), by = "margin", permutations = 99)

print(perm_morph)

nmds_morph <- ordinate(ps_morph, method = "NMDS", distance = "bray", k = 3)

p_morph_fam <- plot_ordination(ps_morph, nmds_morph, color = "plant_family") +
  stat_ellipse(alpha = 0.5) +
  labs(title = "NMDS - Famille végétal", 
       subtitle = paste0("R2 = ", round(perm_morph["plant_family", "R2"], 3), 
                         ", p = ", perm_morph["plant_family", "Pr(>F)"])) +
  theme_bw()

p_morph_morpho <- plot_ordination(ps_morph, nmds_morph, color = "Morphologie") +
  stat_ellipse(alpha = 0.5) +
  labs(title = "NMDS - Morphologie", 
       subtitle = paste0("R2 = ", round(perm_morph["Morphologie", "R2"], 3), 
                         ", p = ", perm_morph["Morphologie", "Pr(>F)"])) +
  theme_bw()

p_morph_fam | p_morph_morpho




meta_lien <- as(sample_data(ps_hel), "data.frame")
meta_lien <- meta_lien[!is.na(meta_lien$plant_genus) & !is.na(meta_lien$Morphologie), ]

table_genus_morph <- table(meta_lien$plant_genus, meta_lien$Morphologie)
print("Table de contingence : plant_genus vs Morphologie")
print(table_genus_morph)

test_chi2 <- chisq.test(table_genus_morph, simulate.p.value = TRUE)
print(test_chi2)

library(FactoMineR)
library(factoextra)

meta_mca <- as(sample_data(ps_hel), "data.frame")
meta_mca <- meta_mca[complete.cases(meta_mca[, c("plant_genus", "Morphologie", "Capacité_Nid")]), ]
df_mca <- meta_mca[, c("plant_genus", "Morphologie", "Capacité_Nid")]

df_mca[] <- lapply(df_mca, as.factor)

res_mca <- MCA(df_mca, graph = FALSE)

p_mca <- fviz_mca_var(res_mca, 
            repel = TRUE, 
            col.var = "black", 
            shape.var = 15,
            ggtheme = theme_minimal()) +
  labs(title = "Analyse des Correspondances Multiples (ACM)",
     subtitle = "Lien entre Genre végétal, Morphologie et Capacité Nid")

print(p_mca)


df_alpha = readRDS("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/df_alpha.rds")
df_alpha$Morphologie <- as(sample_data(ps), "data.frame")[rownames(df_alpha), "Morphologie"]

head(df_alpha)




p_richesse_morph <- ggplot(subset(df_alpha, !is.na(Morphologie)), aes(x = Morphologie, y = Observed, fill = Morphologie)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  theme_bw() +
  labs(title = "Richesse spécifique en fonction de la Morphologie",
       x = "Morphologie",
       y = "Richesse spécifique (Observée)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p_richesse_morph)

df_alpha$Capacité_Nid <- as(sample_data(ps), "data.frame")[rownames(df_alpha), "Capacité_Nid"]

p_richesse_cap <- ggplot(subset(df_alpha, !is.na(Capacité_Nid)), aes(x = Capacité_Nid, y = Observed, fill = Capacité_Nid)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  theme_bw() +
  labs(title = "Richesse spécifique en fonction de la capacité à stocker de l'eau",
       x = "Capacité Stockage Eau",
       y = "Richesse spécifique (Observée)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p_richesse_cap)

p_richesse_morph | p_richesse_cap






p_expsha_morph <- ggplot(subset(df_alpha, !is.na(Morphologie)), aes(x = Morphologie, y = expShannon, fill = Morphologie)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  theme_bw() +
  labs(title = "Indice de diversité exponentiel de Shannon en fonction de la Morphologie",
       x = "Morphologie",
       y = "Indice de diversité exponentiel de Shannon") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p_expsha_morph)

df_alpha$Capacité_Nid <- as(sample_data(ps), "data.frame")[rownames(df_alpha), "Capacité_Nid"]

p_expsha_cap <- ggplot(subset(df_alpha, !is.na(Capacité_Nid)), aes(x = Capacité_Nid, y = expShannon, fill = Capacité_Nid)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  theme_bw() +
  labs(title = "Indice de diversité exponentiel de Shannon en fonction de la capacité à stocker de l'eau",
       x = "Capacité Stockage Eau",
       y = "Indice de diversité exponentiel de Shannon") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p_expsha_cap)
p_expsha_morph | p_expsha_cap

(p_richesse_morph | p_richesse_cap) / (p_expsha_morph | p_expsha_cap) + plot_annotation(title = "Richesse et diversité en fonction de la Morphologie et de la Capacité à stocker de l'eau")

p_expsha_morph | p_expsha_cap + plot_annotation(title = "Diversité en fonction de la Morphologie et de la Capacité à stocker de l'eau")

