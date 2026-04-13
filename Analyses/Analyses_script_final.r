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



ps_trim <- filter_taxa(ps, function(x) sum(x) > 99, TRUE)

ps_merged <- merge_samples(subset_samples(ps_trim, position == "epiphyte"), "organ")
otu_tab <- as.matrix(t(otu_table(ps_merged)))

otu_num <- apply(otu_tab, 2, function(x) as.numeric(as.character(x)))

otu_num[is.na(otu_num)] <- 0

rownames(otu_num) <- taxa_names(ps_merged)
otu_grouped <- as.data.frame(otu_num)


data_list <- as.list(otu_grouped)
data_list <- lapply(data_list, as.integer)
write.table(data_list, "data_list.txt", sep = "\t", row.names = TRUE, col.names = NA)

max_reads <- max(colSums(otu_grouped))
paliers <- seq(10, 3649011, length.out = 20)


res_inext <- iNEXT3D(data = data_list, 
                   diversity = "TD", 
                   q = 0, 
                   datatype = "abundance",
                   nboot = 0,       
                   size = paliers,
                   endpoint = 3649011) 


saveRDS(res_inext, file = "res_inext.rds")

ggiNEXT3D(res_inext, type = 1)












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

ps_trim <- filter_taxa(ps, function(x) sum(x) > 1000, TRUE)


mat_epi <- clean_otu(subset_samples(ps_trim, position == "epiphyte"), "organ")
res_epi <- iNEXT3D(mat_epi, q=c(0,1,2), datatype="abundance", nboot=0, size=seq(10, max(colSums(mat_epi)), length.out=10))

mat_endo <- clean_otu(subset_samples(ps_trim, position == "endophyte"), "organ")
res_endo <- iNEXT3D(mat_endo, q=c(0,1,2), datatype="abundance", nboot=0, size=seq(10, max(colSums(mat_endo)), length.out=10))

mat_proj <- clean_otu(ps_trim, "project")
res_proj <- iNEXT3D(mat_proj, q=c(0,1,2), datatype="abundance", nboot=0, size=seq(10, max(colSums(mat_proj)), length.out=10))


prepare_plot <- function(res, title) {
  p <- ggiNEXT3D(res, type = 1, facet_var = "Order.q") +
    theme_bw() +
    labs(title = title, x = NULL, y = "Diversity") +
    theme(legend.position = "right")
  
  levels(p$data$Order.q) <- c("Richness (q=0)", "Shannon (q=1)", "Simpson (q=2)")
  return(p)
}

p1 <- prepare_plot(res_epi, "Position: Epiphytes")
p2 <- prepare_plot(res_endo, "Position: Endophytes")
p3 <- prepare_plot(res_proj, "Global Projects") + labs(x = "Number of individuals")



final_plot <- p1 / p2 / p3

print(final_plot)

ggsave("Figure_Alpha_Diversity_Complete.pdf", final_plot, width = 14, height = 12)



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



