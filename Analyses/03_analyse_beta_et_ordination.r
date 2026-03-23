setwd("~/Documents/Etudes/Stage_projet_LOT/CRBE/données/analyses_di")

ps <- readRDS("ps_final.rds")

library(vegan)

ps_hellinger <- transform_sample_counts(ps, function(x) sqrt(x / sum(x)))

ps_hellinger <- prune_samples(sample_sums(ps_hellinger) > 0, ps_hellinger)

otu_hellinger <- as(otu_table(ps_hellinger), "matrix")

bray_dist <- vegdist(otu_hellinger, method = "bray", na.rm = TRUE)

# jaccard_dist <- vegdist(otu_table > 0, method = "jaccard")


#nmds_bray_K2 <- ordinate(ps_hellinger, method = "NMDS", distance = "bray")

nmds_bray_K3 <- ordinate(ps_hellinger, method = "NMDS", distance = "bray", k = 3)


library(goeveg)

dimcheckMDS(as(otu_table(ps_hellinger), "matrix"), distance = "bray", k = 5, trymax = 10)





bray_dist <- vegdist(as(otu_table(ps_hellinger), "matrix"), method = "bray", na.rm = TRUE)

perm1=adonis2(bray_dist ~ project + organ + position, by = "margin" ,data = as(sample_data(ps_hellinger), "data.frame"), permutations = 99)
perm1


plot_ordination(ps_hellinger, nmds_bray_K3, color = "project", shape = "organ") +
  stat_ellipse(aes(group = project, fill = project), alpha = 0.1, geom = "polygon") +
  annotate("text", x = Inf, y = Inf, 
           label = paste0("PERMANOVA: R²=", round(perm1$R2[1], 2), ", p=", perm1$`Pr(>F)`[1]), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
  theme_bw() +
  labs(title = "NMDS - Bray-Curtis (K=3)")



library(plotly)

nmds_3d_data <- as.data.frame(scores(nmds_bray_K3, display = "sites"))
nmds_3d_data$organ <- sample_data(ps_hellinger)$organ
nmds_3d_data$position <- sample_data(ps_hellinger)$position
nmds_3d_data$project <- sample_data(ps_hellinger)$project

plot_ly(nmds_3d_data, x = ~NMDS1, y = ~NMDS2, z = ~NMDS3, 
  color = ~organ, symbol = ~position, text = ~project,
  type = "scatter3d", mode = "markers", 
  marker = list(size = 5, opacity = 0.8)) %>%
  layout(title = "NMDS 3D - Bray-Curtis (K=3)")




ps_epiphytes <- subset_samples(ps_hellinger, position == "epiphyte")
nmds_epi <- ordinate(ps_epiphytes, method = "NMDS", distance = "bray")

bray_epi <- vegdist(as(otu_table(ps_epiphytes), "matrix"), method = "bray")
perm_epi <- adonis2(bray_epi ~ organ, data = as(sample_data(ps_epiphytes), "data.frame"), permutations = 999)

p_epi=plot_ordination(ps_epiphytes, nmds_epi, color = "organ", shape = "project") +
  stat_ellipse(aes(group = organ, fill = organ), alpha = 0.1, geom = "polygon") +
  annotate("text", x = Inf, y = Inf, 
           label = paste0("PERMANOVA: R²=", round(perm_epi$R2[1], 2), ", p=", perm_epi$`Pr(>F)`[1]), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
  theme_bw() 

perm_epi


ps_endophytes <- subset_samples(ps_hellinger, position == "endophyte")
nmds_endo <- ordinate(ps_endophytes, method = "NMDS", distance = "bray")

bray_endo <- vegdist(as(otu_table(ps_endophytes), "matrix"), method = "bray")
perm_endo <- adonis2(bray_endo ~ organ, data = as(sample_data(ps_endophytes), "data.frame"), permutations = 99)

p_endo=plot_ordination(ps_endophytes, nmds_endo, color = "organ", shape = "project") +
  stat_ellipse(aes(group = organ, fill = organ), alpha = 0.1, geom = "polygon") +
  annotate("text", x = Inf, y = Inf, 
           label = paste0("PERMANOVA: R²=", round(perm_endo$R2[1], 2), ", p=", perm_endo$`Pr(>F)`[1]), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
  theme_bw() 

library(pathwork)
(p_epi + p_endo) + 
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Comparaison des ordinations NMDS (Épiphytes vs Endophytes)")



  nmds_family <- ordinate(ps_hellinger, method = "NMDS", distance = "bray")

  bray_family <- vegdist(as(otu_table(ps_hellinger), "matrix"), method = "bray", na.rm = TRUE)
  df_meta_family <- as(sample_data(ps_hellinger), "data.frame")
  perm_family <- adonis2(bray_family ~ plant_family, data = df_meta_family, permutations = 999)

  plot_ordination(ps_hellinger, nmds_family, color = "plant_family") +
    stat_ellipse(aes(group = plant_family, fill = plant_family), alpha = 0.1, geom = "polygon") +
    annotate("text", x = Inf, y = Inf, 
             label = paste0("PERMANOVA: R²=", round(perm_family$R2[1], 2), ", p=", perm_family$`Pr(>F)`[1]), 
             hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
    theme_bw() +
    labs(title = "NMDS - Taxonomie des plantes (Famille)",
         color = "Famille de plante")




  ps_leaves <- subset_samples(ps_hellinger, tolower(organ) %in% c("leaf", "leaves", "young_leaf", "old_leaf"))

  ps_leaves <- prune_taxa(taxa_sums(ps_leaves) > 0, ps_leaves)

  nmds_leaves <- ordinate(ps_leaves, method = "NMDS", distance = "bray")

  bray_leaves <- vegdist(as(otu_table(ps_leaves), "matrix"), method = "bray", na.rm = TRUE)
  df_meta_leaves <- as(sample_data(ps_leaves), "data.frame")
  perm_leaves <- adonis2(bray_leaves ~ organ + position, data = df_meta_leaves, permutations = 999)

  plot_ordination(ps_leaves, nmds_leaves, color = "organ", shape = "position") +
    stat_ellipse(aes(group = organ, fill = organ), alpha = 0.1, geom = "polygon") +
    annotate("text", x = Inf, y = Inf, 
             label = paste0("Age PERMANOVA: R²=", round(perm_leaves$R2[1], 2), ", p=", perm_leaves$`Pr(>F)`[1]), 
             hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
    theme_bw() +
    labs(title = "NMDS - Succession foliaire (Jeune vs Vieille feuille)",
         color = "Âge",
         shape = "Compartiment (Position)")



metadata_global <- as(sample_data(ps), "data.frame")
modele_lm <- lm(nmds_bray_K3$points[, 1] ~ organ + plant_family+project, data = metadata_global)
summary(modele_lm)
shapiro.test(residuals(modele_lm))

par(mfrow = c(2, 2))
plot(modele_lm)
par(mfrow = c(1, 1))

library(MASS)

adjusted_nmds <- round(nmds_bray_K3$points[, 1] - min(nmds_bray_K3$points[, 1]))
glm_test=glm.nb(adjusted_nmds ~ organ + plant_family+project, data = metadata_global)
summary(glm_test)



metadata_global <- as(sample_data(ps), "data.frame")

dispersion_organ <- betadisper(bray_dist, metadata_global$organ)
permutest(dispersion_organ, permutations = 999)
plot(dispersion_organ, main = "Homogénéité des dispersions multivariées (Organe)")

dispersion_family <- betadisper(bray_dist, metadata_global$plant_family)
permutest(dispersion_family, permutations = 999)
plot(dispersion_family, main = "Homogénéité des dispersions multivariées (Famille de plante)")