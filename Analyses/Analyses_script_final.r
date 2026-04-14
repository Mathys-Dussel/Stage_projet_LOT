setwd("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses")

ps <- readRDS("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/ps_final.rds")
df_alpha <- readRDS("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/df_alpha.rds")

library(phyloseq)
library(iNEXT.3D)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggVennDiagram)
library(dplyr)
library(UpSetR)
library(microbiome)
library(circlize)
library(vegan)
library(tidyr)



############## courbes raréfactions et extrapolations ##############

library(iNEXT.3D)
library(ggplot2)
library(patchwork)

clean_otu <- function(ps_obj, variable) {
  ps_m <- merge_samples(ps_obj, variable)
  tab <- as.matrix(t(otu_table(ps_m)))
  tab_num <- apply(tab, 2, function(x) as.integer(as.character(x)))
  tab_num[is.na(tab_num)] <- 0
  return(tab_num)
}


ps_trim <- filter_taxa(ps, function(x) sum(x) > 500, TRUE)


total_reads_avant <- sum(sample_sums(ps))
total_reads_apres <- sum(sample_sums(ps_trim))
(total_reads_apres / total_reads_avant) * 100



mat_epi <- clean_otu(subset_samples(ps_trim, position == "epiphyte"), "organ")
res_epi <- iNEXT3D(mat_epi, q=c(0,1,2), datatype="abundance", nboot=0, size=seq(10, max(1.5*colSums(mat_epi)), length.out=20))
saveRDS(res_epi, file = "res_epi.rds")

mat_endo <- clean_otu(subset_samples(ps_trim, position == "endophyte"), "organ")
res_endo <- iNEXT3D(mat_endo, q=c(0,1,2), datatype="abundance", nboot=0, size=seq(10, max(1.5*colSums(mat_endo)), length.out=20))
saveRDS(res_endo, file = "res_endo.rds")

mat_proj <- clean_otu(ps_trim, "project")
res_proj <- iNEXT3D(mat_proj, q=c(0,1,2), datatype="abundance", nboot=0, size=seq(10, max(1.5*colSums(mat_proj)), length.out=20))
saveRDS(res_proj, file = "res_proj.rds")






library(dplyr)
library(tidyr)

extraire_tableau_complet <- function(res_obj, categorie) {
  asy_data <- res_obj$TDAsyEst %>%
    select(Assemblage, qTD, TD_asy, s.e.) %>%
    mutate(across(c(TD_asy, s.e.), ~ round(., 3)))
  
  asy_format <- asy_data %>%
    mutate(Value_Display = ifelse(is.na(s.e.), 
                                  as.character(TD_asy), 
                                  paste0(TD_asy, " ± ", s.e.))) %>%
    select(Assemblage, qTD, Value_Display) %>%
    pivot_wider(names_from = qTD, values_from = Value_Display)
  
  equit_data <- asy_data %>%
    select(Assemblage, qTD, TD_asy) %>%
    pivot_wider(names_from = qTD, values_from = TD_asy) %>%
    mutate(Equitability = round(`Shannon diversity` / `Species richness`, 2)) %>%
    select(Assemblage, Equitability)
  
  final_tab <- left_join(asy_format, equit_data, by = "Assemblage") %>%
    mutate(Category = categorie) %>%
    select(Category, everything())
  
  return(final_tab)
}

tab_epi_complet  <- extraire_tableau_complet(res_epi, "Epiphytes")
tab_endo_complet <- extraire_tableau_complet(res_endo, "Endophytes")
tab_proj_complet <- extraire_tableau_complet(res_proj, "Projets")

tableau_final_global <- bind_rows(tab_epi_complet, tab_endo_complet, tab_proj_complet)

print(tableau_final_global)

write.csv(tableau_final_global, "Tableau_Resultats_Synthese_iNEXT.csv", row.names = FALSE)




















subset_q <- function(res, q_val) {
  res_sub <- res
  res_sub$iNextEst$size_based <- res_sub$iNextEst$size_based[res_sub$iNextEst$size_based$Order.q == q_val, ]
  if ("coverage_based" %in% names(res_sub$iNextEst)) {
    res_sub$iNextEst$coverage_based <- res_sub$iNextEst$coverage_based[res_sub$iNextEst$coverage_based$Order.q == q_val, ]
  }
  return(res_sub)
}


p1 <- ggiNEXT3D(res_epi, type = 1, facet.var = "Order.q") + 
          facet_wrap(~Order.q, scales = "free") +
          ggtitle("Epiphyte Diversity (q = 0, 1, 2)") +
          theme_bw() +
          theme(strip.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 14))

p2 <- ggiNEXT3D(res_endo, type = 1, facet.var = "Order.q") + 
          facet_wrap(~Order.q, scales = "free") +
          ggtitle("Endophyte Diversity (q = 0, 1, 2)") +
          theme_bw() +
          theme(strip.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 14))

p3 <- ggiNEXT3D(res_proj, type = 1, facet.var = "Order.q") + 
          facet_wrap(~Order.q, scales = "free") +
          ggtitle("Project Diversity (q = 0, 1, 2)") +
          theme_bw() +
          theme(strip.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 14))






p1 <- p1 + ylab("Epiphyte Diversity")
p2 <- p2 + ylab("Endophyte Diversity")
p3 <- p3 + ylab("Project Diversity")

final_plot <- (p1 / p2 / p3) + 
  plot_layout(guides = "collect") & 
  theme_bw() & 
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) & 
    scale_color_manual(values = c("root" = "brown", "old_leaf" = "darkgreen", "young_leaf" = "lightgreen", "LOT1" = "lightblue", "LOT2" = "#D34949", "LOT3" = "grey")) &
    scale_fill_manual(values = c("root" = "brown", "old_leaf" = "darkgreen", "young_leaf" = "lightgreen", "LOT1" = "lightblue", "LOT2" = "#D34949", "LOT3" = "grey")) & 
  geom_line(linewidth = 0.5) & 
  geom_point(size = 1)

print(final_plot)





############## UpSetR ##############

ps_bin <- transform_sample_counts(ps_rare, function(x) ifelse(x > 0, 1, 0))

meta <- data.frame(sample_data(ps))
meta$category <- paste(meta$project, meta$age, meta$organ, meta$position, sep = " ")

taxa_list <- split(rownames(meta), meta$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

upset(fromList(taxa_list), nsets = length(taxa_list), order.by = "freq", 
      main.bar.color = "forestgreen", sets.bar.color = "steelblue")



meta <- data.frame(sample_data(ps))
meta <- subset(meta, project == "LOT1")
meta$category <- paste(meta$age, meta$organ, meta$position, sep = " ")

taxa_list <- split(rownames(meta), meta$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

upset(fromList(taxa_list), nsets = length(taxa_list), order.by = "freq", 
      main.bar.color = "forestgreen", sets.bar.color = "steelblue")



meta <- data.frame(sample_data(ps))
meta <- subset(meta, project == "LOT2")
meta$category <- paste(meta$age, meta$organ, meta$position, sep = " ")

taxa_list <- split(rownames(meta), meta$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

upset(fromList(taxa_list), nsets = length(taxa_list), order.by = "freq", 
      main.bar.color = "forestgreen", sets.bar.color = "steelblue")




meta <- data.frame(sample_data(ps))
meta <- subset(meta, project == "LOT3")
meta$category <- paste(meta$age, meta$organ, meta$position, sep = " ")

taxa_list <- split(rownames(meta), meta$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

upset(fromList(taxa_list), nsets = length(taxa_list), order.by = "freq", 
      main.bar.color = "forestgreen", sets.bar.color = "steelblue")




library(UpSetR)
library(cowplot)
library(magick)
library(patchwork)

# --- 1. Fonction de capture (pour éviter l'erreur de vecteur) ---
generate_upset_captured <- function(taxa_list, title_text = "") {
  tmp <- tempfile(fileext = ".png")
  # On ajuste la résolution pour que ce soit lisible
  png(tmp, width = 1000, height = 800, res = 120)
  
  print(upset(fromList(taxa_list), 
              nsets = length(taxa_list), 
              order.by = "freq", 
              main.bar.color = "forestgreen", 
              sets.bar.color = "steelblue"))
  
  # Optionnel : ajout d'un titre simple via grid
  grid::grid.text(title_text, x = 0.5, y = 0.95, gp = grid::gpar(fontsize = 15, font = 2))
  
  dev.off()
  
  img <- magick::image_read(tmp)
  return(cowplot::ggdraw() + cowplot::draw_image(img))
}

# --- 2. Préparation des données pour chaque graphique ---

# A. Global (Tous les projets)
meta_all <- data.frame(sample_data(ps))
meta_all$category <- paste(meta_all$project, meta_all$age, meta_all$organ, meta_all$position, sep = " ")
taxa_all <- split(rownames(meta_all), meta_all$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

# B. LOT1
meta_lot1 <- subset(data.frame(sample_data(ps)), project == "LOT1")
meta_lot1$category <- paste(meta_lot1$age, meta_lot1$organ, meta_lot1$position, sep = " ")
taxa_lot1 <- split(rownames(meta_lot1), meta_lot1$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

# C. LOT2
meta_lot2 <- subset(data.frame(sample_data(ps)), project == "LOT2")
meta_lot2$category <- paste(meta_lot2$age, meta_lot2$organ, meta_lot2$position, sep = " ")
taxa_lot2 <- split(rownames(meta_lot2), meta_lot2$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

# D. LOT3
meta_lot3 <- subset(data.frame(sample_data(ps)), project == "LOT3")
meta_lot3$category <- paste(meta_lot3$age, meta_lot3$organ, meta_lot3$position, sep = " ")
taxa_lot3 <- split(rownames(meta_lot3), meta_lot3$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

# --- 3. Génération et capture des 4 plots ---
p1 <- generate_upset_captured(taxa_all, "Global")
p2 <- generate_upset_captured(taxa_lot1, "Project LOT1")
p3 <- generate_upset_captured(taxa_lot2, "Project LOT2")
p4 <- generate_upset_captured(taxa_lot3, "Project LOT3")

# --- 4. Affichage final ---
(p1 + p2) / (p3 + p4)











############## Tableau abondances ##############



library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(pheatmap)
library(tidyr)


ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))


ps_phylum <- tax_glom(ps_rel, taxrank = "gbr268_Phylum", NArm = FALSE)

df_phylum <- psmelt(ps_phylum)



ps_filtered <- filter_taxa(ps, function(x) sum(x > 0) >= (0.05 * length(x)), TRUE)

ps_rel <- transform_sample_counts(ps_filtered, function(x) x / sum(x))


ps_class <- tax_glom(ps_rel, taxrank = "gbr268_Order", NArm = FALSE)

df_class <- psmelt(ps_class)

df_agg <- df_class %>%
    group_by(gbr268_Order, project, organ, position) %>%
    summarise(Abundance = mean(Abundance), .groups = "drop") %>%
    mutate(group = paste(project, organ, position, sep = "_"))

df_heatmap <- df_agg %>%
    dplyr::select(gbr268_Order, group, Abundance) %>%
    pivot_wider(names_from = group, values_from = Abundance, values_fill = 0) %>%
    filter(!is.na(gbr268_Order)) %>%
    as.data.frame()

rownames(df_heatmap) <- df_heatmap$gbr268_Order
df_heatmap$gbr268_Order <- NULL
mat <- as.matrix(df_heatmap)

col_meta <- df_agg %>%
    dplyr::select(group, project, organ, position) %>%
    distinct() %>%
    arrange(match(group, colnames(mat)))

col_ha <- HeatmapAnnotation(
    Project = col_meta$project,
    Organ = col_meta$organ,
    Position = col_meta$position,
    annotation_name_side = "right"
)

row_meta <- df_class %>%
    dplyr::select(gbr268_Order, gbr268_Phylum, gbr268_Class) %>%
    distinct() %>%
    filter(!is.na(gbr268_Order)) %>%
    arrange(match(gbr268_Order, rownames(mat)))

row_ha <- rowAnnotation(
    Phylum = row_meta$gbr268_Phylum,
    Class  = row_meta$gbr268_Class,
    annotation_name_side = "top"
)

color_mapping <- colorRamp2(c(0, max(mat)/2, max(mat)), c("lightblue", "yellow", "red"))
ann_colors = list(
    Project = c("LOT1" = "#7FB3D5", "LOT2" = "#D98880", "LOT3" = "#BDC3C7"),
    Organ = c("root" = "#A93226", "old_leaf" = "#27AD60", "young_leaf" = "#A2D9CE"),
    Position = c("endophyte" = "#FFDAB5", "epiphyte" = "#C0D461")
)

col_ha <- HeatmapAnnotation(
    Project = col_meta$project,
    Organ = col_meta$organ,
    Position = col_meta$position,
    col = ann_colors,
    annotation_name_side = "right"
)

ht <- Heatmap(mat,
    name = "Relative\nabundance",
    col = color_mapping,
    top_annotation = col_ha,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    left_annotation = row_ha,
    row_title = NULL,
    column_title = NULL,
    column_split = col_meta$project,
    heatmap_legend_param = list(direction = "vertical"),
    row_split = row_meta$gbr268_Phylum,
    rect_gp = gpar(col = "black", lwd = 1),
    row_names_side = "right",
    column_names_side = "bottom",
    column_names_rot = -90,
    row_title_rot = 0
)

draw(ht, merge_legend = TRUE)









############## Ordinations ##############

ps <- readRDS("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/ps_final.rds")

library(vegan)

ps_hel <- transform_sample_counts(ps, function(x) sqrt(x / sum(x)))
ps_hel <- prune_samples(sample_sums(ps_hel) > 0, ps_hel)
otu_hel <- as(otu_table(ps_hel), "matrix")
metadata <- as(sample_data(ps_hel), "data.frame")

bray_dist <- vegdist(otu_hel, method = "bray")

perm_global <- adonis2(bray_dist ~ project + organ + position, data = metadata, by = "margin", permutations = 99)
nmds_global <- ordinate(ps_hel, method = "NMDS", distance = "bray", k = 3)

p1 <- plot_ordination(ps_hel, nmds_global, color = "project", shape = "organ") +
  stat_ellipse(aes(group = project, fill = project), alpha = 0.1, geom = "polygon") +
  labs(title = "A. Structure Globale (Projet)", 
       subtitle = paste0("R2 Project: ", round(perm_global$R2[1], 3))) +
  theme_bw()








ps_lot01 <- subset_samples(ps_hel, project == "LOT1")
ps_lot01 <- prune_taxa(taxa_sums(ps_lot01) > 0, ps_lot01)

meta_lot01 <- as(sample_data(ps_lot01), "data.frame")
valid_samples <- complete.cases(meta_lot01[, c("organ", "position")])
ps_lot01 <- prune_samples(valid_samples, ps_lot01)

bray_lot01 <- vegdist(otu_table(ps_lot01), method = "bray")
perm_lot01 <- adonis2(bray_lot01 ~ organ + position, data = as(sample_data(ps_lot01), "data.frame"), by = "margin", permutations = 999)

nmds_lot01 <- ordinate(ps_lot01, method = "NMDS", distance = "bray", k = 3)

p2=plot_ordination(ps_lot01, nmds_lot01, color = "organ", shape = "position") +
  geom_point(size = 3) +
  scale_shape_manual(values = c("epiphyte" = 1, "endophyte" = 16)) +
  stat_ellipse(aes(group = organ)) +
  labs(title = "B. LOT1 - Organe", subtitle = paste0("R2 = ", round(perm_lot01["organ", "R2"], 3), ", p = ", perm_lot01["organ", "Pr(>F)"])) +
  theme_bw()







ps_lot02 <- subset_samples(ps_hel, project == "LOT2")
ps_lot02 <- prune_taxa(taxa_sums(ps_lot02) > 0, ps_lot02)

meta_lot02 <- as(sample_data(ps_lot02), "data.frame")
valid_samples <- complete.cases(meta_lot02[, c("organ", "position")])
ps_lot02 <- prune_samples(valid_samples, ps_lot02)


bray_lot02 <- vegdist(otu_table(ps_lot02), method = "bray")
vars_lot02 <- c("organ", "position")
valid_vars_lot02 <- vars_lot02[sapply(meta_lot02[valid_samples, vars_lot02], function(x) length(unique(x)) > 1)]
formula_lot02 <- as.formula(paste("bray_lot02 ~", paste(valid_vars_lot02, collapse = " + ")))

perm_lot02 <- adonis2(formula_lot02, data = as(sample_data(ps_lot02), "data.frame"), by = "margin", permutations = 99)

nmds_lot02 <- ordinate(ps_lot02, method = "NMDS", distance = "bray", k = 3)


p3=plot_ordination(ps_lot02, nmds_lot02, color = "organ", shape = "position") +
  geom_point(size = 3) +
  scale_shape_manual(values = c("epiphyte" = 1, "endophyte" = 16)) +
  stat_ellipse(aes(group = organ)) +
  labs(title = "C. LOT2 - Organe", subtitle = paste0("R2 = ", round(perm_lot02["organ", "R2"], 3), ", p = ", perm_lot02["organ", "Pr(>F)"])) +
  theme_bw()







ps_lot03 <- subset_samples(ps_hel, project == "LOT3")
ps_lot03 <- prune_taxa(taxa_sums(ps_lot03) > 0, ps_lot03)

bray_lot03 <- vegdist(otu_table(ps_lot03), method = "bray")
perm_lot03 <- adonis2(bray_lot03 ~ organ + position, data = as(sample_data(ps_lot03), "data.frame"), by = "margin", permutations = 99)

nmds_lot03 <- ordinate(ps_lot03, method = "NMDS", distance = "bray", k = 3)

p4=plot_ordination(ps_lot03, nmds_lot03, color = "organ", shape = "position") +
  geom_point(size = 3) +
  scale_shape_manual(values = c("epiphyte" = 1, "endophyte" = 16)) +
  stat_ellipse(aes(group = organ)) +
  labs(title = "D. LOT3 - Organe", subtitle = paste0("R2 = ", round(perm_lot03["organ", "R2"], 3), ", p = ", perm_lot03["organ", "Pr(>F)"])) +
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

  organ_colors <- c("root" = "#794515", "old_leaf" = "forestgreen", "young_leaf" = "#4DCCBD")
  p2 <- p2 + scale_color_manual(values = organ_colors)
  p3 <- p3 + scale_color_manual(values = organ_colors)
  p4 <- p4 + scale_color_manual(values = organ_colors)

(p1 | p2) / (p3 | p4) 



########### Niches de Levins ###########


library(dplyr)
library(phyloseq)
library(ggplot2)
library(ggVennDiagram)
library(microbiome)
library(RColorBrewer)
library(cluster)
library(ape)
library(ggpubr)
library(ggtree)
library(spaa)
library(vegan)
library(ANCOMBC)
library(tidyr)


ps_merged <- merge_samples(ps, "organ")
otu_merged <- as.matrix(otu_table(ps_merged))
if (taxa_are_rows(ps_merged)) { otu_merged <- t(otu_merged) }

levins_raw <- niche.width(otu_merged, method = "levins")
niche_scores <- (as.numeric(levins_raw) - 1) / (3 - 1)
names(niche_scores) <- colnames(levins_raw)

ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
otu_all <- as.matrix(otu_table(ps_rel))
if (taxa_are_rows(ps_rel)) { otu_all <- t(otu_all) }

calc_mean_niche <- function(sample_row, scores) {
  present <- sample_row > 0
  if (sum(present) == 0) return(NA)
  
  sub_scores <- scores[names(sample_row)[present]]
  sub_abund <- sample_row[present]
  
  return(sum(sub_abund * sub_scores, na.rm = TRUE) / sum(sub_abund, na.rm = TRUE))
}

sample_indices <- apply(otu_all, 1, calc_mean_niche, scores = niche_scores)

metadata <- data.frame(sample_data(ps))
final_df <- data.frame(
  Sample = names(sample_indices),
  Levins_Mean = as.numeric(sample_indices)
) %>%
  mutate(organ = metadata$organ[match(Sample, rownames(metadata))])

niche_summary <- final_df %>%
  filter(!is.nan(Levins_Mean) & !is.na(Levins_Mean)) %>%
  group_by(organ) %>%
  summarise(
    Moyenne = mean(Levins_Mean),
    SD = sd(Levins_Mean),
    n = n()
  )

print(niche_summary)



final_df$project <- metadata$project[match(final_df$Sample, rownames(metadata))]
final_df$project <- as.factor(final_df$project)


final_df$organ <- factor(
    final_df$organ,
    levels = c("root", "young_leaf", "old_leaf")
)



p_project_organ <- ggplot(final_df, aes(x = project, y = Levins_Mean, fill = organ)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.6, position = position_dodge(0.8)) +
  geom_point(aes(color = organ), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), alpha = 0.4, size = 1) +
  labs(title = "Niche fongique par projet et organe",
       x = "Projet",
       y = "Indice de Levins Moyen",
       fill = "organ",
       color = "organ") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 11, face = "bold"),
        title = element_text(size = 13, face = "bold"))

p_project_organ_final <- p_project_organ + 
  stat_compare_means(aes(group = organ), method = "kruskal.test", label = "p.signif")

p_project_organ_final <- p_project_organ_final +
    scale_fill_manual(
        values = c(
            "root"       = "#914e27",
            "young_leaf" = "#83d483",
            "old_leaf"   = "#1d5d1d"
        )
    ) +
    scale_color_manual(
        values = c(
            "root"       = "#914e27",
            "young_leaf" = "#83d483",
            "old_leaf"   = "#1d5d1d"
        )
    )

print(p_project_organ_final)



