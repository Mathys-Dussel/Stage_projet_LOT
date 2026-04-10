setwd("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses")

ps <- readRDS("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/ps_final.rds")
df_alpha <- readRDS("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/df_alpha.rds")

library(phyloseq)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggVennDiagram)
library(dplyr)
library(UpSetR)
library(microbiome)
library(circlize)




############## courbes raréfactions et extrapolations ##############







############## UpSetR ##############

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
    Project = c("LOT1" = "lightblue", "LOT2" = "#D34949", "LOT3" = "grey"),
    Organ = c("root" = "brown", "old_leaf" = "darkgreen", "young_leaf" = "lightgreen"),
    Position = c("endophyte" = "orange", "epiphyte" = "pink")
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




############# TreeMap ##############


# install.packages("treemapify")  # Run this once manually if 'treemapify' is not installed
library(treemapify)

df_treemap <- df_class%>%
    filter(organ == "root", position == "epiphyte") %>%
    group_by(project, gbr268_Class) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    mutate(gbr268_Phylum = ifelse(is.na(gbr268_Class), "Unclassified", gbr268_Class))

ggplot(df_treemap, aes(area = Abundance, fill = gbr268_Class, label = gbr268_Class)) +
    treemapify::geom_treemap(color = "white") +
    treemapify::geom_treemap_text(colour = "black", place = "centre", grow = TRUE, reflow = TRUE) +
    theme_bw() +
    theme(legend.position = "right")




df_treemap <- df_class%>%
    filter(organ == "young_leaf", position == "epiphyte") %>%
    group_by(project, gbr268_Class) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    mutate(gbr268_Phylum = ifelse(is.na(gbr268_Class), "Unclassified", gbr268_Class))

ggplot(df_treemap, aes(area = Abundance, fill = gbr268_Class, label = gbr268_Class)) +
    treemapify::geom_treemap(color = "white") +
    treemapify::geom_treemap_text(colour = "black", place = "centre", grow = TRUE, reflow = TRUE) +
    theme_bw() +
    theme(legend.position = "right")




    df_treemap <- df_class%>%
    filter(organ == "old_leaf", position == "epiphyte") %>%
    group_by(project, gbr268_Class) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    mutate(gbr268_Phylum = ifelse(is.na(gbr268_Class), "Unclassified", gbr268_Class))

ggplot(df_treemap, aes(area = Abundance, fill = gbr268_Class, label = gbr268_Class)) +
    treemapify::geom_treemap(color = "white") +
    treemapify::geom_treemap_text(colour = "black", place = "centre", grow = TRUE, reflow = TRUE) +
    theme_bw() +
    theme(legend.position = "right")




    

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
  labs(title = "1. Structure Globale (Projet)", 
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
  labs(title = "LOT1 - Organe", subtitle = paste0("R2 = ", round(perm_lot01["organ", "R2"], 3), ", p = ", perm_lot01["organ", "Pr(>F)"])) +
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
  labs(title = "LOT2 - Organe", subtitle = paste0("R2 = ", round(perm_lot02["organ", "R2"], 3), ", p = ", perm_lot02["organ", "Pr(>F)"])) +
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



(p1 | p2) / (p3 | p4) 



########### Niches de Levins ###########



levins_col <- grep("levins", names(df_alpha), ignore.case = TRUE, value = TRUE)[1]

meta_lev <- as(sample_data(ps), "data.frame")
meta_lev <- cbind(meta_lev, df_alpha)

p_levins_project <- ggplot(meta_lev,
                           aes_string(x = "project", y = levins_col, fill = "organ")) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Projet", y = "Indice de niche de Levins")

p_levins_organ <- ggplot(meta_lev,
                         aes_string(x = "organ", y = levins_col, fill = "project")) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Organe", y = "Indice de niche de Levins")

p_levins_project + p_levins_organ

