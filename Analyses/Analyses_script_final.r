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
library(scales)
library(tidyr)

############## courbes raréfactions et extrapolations ##############


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
res_epi <- iNEXT3D(
  mat_epi,
  q = c(0, 1, 2),
  datatype = "abundance",
  nboot = 50,
  size = seq(10, max(1.5 * colSums(mat_epi)), length.out = 40)
)
saveRDS(res_epi, file = "res_epi.rds")

mat_endo <- clean_otu(subset_samples(ps_trim, position == "endophyte"), "organ")
res_endo <- iNEXT3D(
  mat_endo,
  q = c(0, 1, 2),
  datatype = "abundance",
  nboot = 50,
  size = seq(10, max(1.5 * colSums(mat_endo)), length.out = 40)
)
saveRDS(res_endo, file = "res_endo.rds")

mat_proj <- clean_otu(ps_trim, "project")
res_proj <- iNEXT3D(
  mat_proj,
  q = c(0, 1, 2),
  datatype = "abundance",
  nboot = 50,
  size = seq(10, max(1.5 * colSums(mat_proj)), length.out = 40)
)
saveRDS(res_proj, file = "~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees_inext/res_proj.rds")

library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)



res_epi = readRDS( "~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/donnees_inext/res_epi.rds")
res_endo = readRDS("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/donnees_inext/res_endo.rds")
res_proj = readRDS("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/donnees_inext/res_proj.rds")



extraire_tableau_complet <- function(res_obj, categorie) {
asy_data <- res_obj$TDAsyEst %>%
select(Assemblage, qTD, TD_asy, s.e.) %>%
mutate(across(c(TD_asy, s.e.), ~ round(., 3)))

asy_format <- asy_data %>%
mutate(
Value_Display = ifelse(
is.na(s.e.),
as.character(TD_asy),
paste0(TD_asy, " ± ", s.e.)
)
) %>%
select(Assemblage, qTD, Value_Display) %>%
pivot_wider(names_from = qTD, values_from = Value_Display)

equit_data <- asy_data %>%
select(Assemblage, qTD, TD_asy) %>%
pivot_wider(names_from = qTD, values_from = TD_asy) %>%
mutate(
Equitability = round(`Shannon diversity` / `Species richness`, 2)
) %>%
select(Assemblage, Equitability)

final_tab <- left_join(asy_format, equit_data, by = "Assemblage") %>%
mutate(Category = categorie) %>%
select(Category, everything())

return(final_tab)
}


tab_epi_complet <- extraire_tableau_complet(res_epi, "Epiphytes")
tab_endo_complet <- extraire_tableau_complet(res_endo, "Endophytes")
tab_proj_complet <- extraire_tableau_complet(res_proj, "Projets")

tableau_final_global <- bind_rows(tab_epi_complet, tab_endo_complet, tab_proj_complet)

print(tableau_final_global)

write.csv(tableau_final_global,
"~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/donnees_inext/Tableau_Resultats_Synthese_iNEXT.csv",
row.names = FALSE)







# Couleurs pour les deux types de graphiques
custom_colors <- c("root" = "#8B4513", "old_leaf" = "#2E7D32", "young_leaf" = "#81C784")
project_colors <- c("LOT1" = "#1F78B4", "LOT2" = "#E31A1C", "LOT3" = "#33A02C")

plot_final_premium <- function(res_obj, title_label, palette_couleurs) {
  
  df <- as.data.frame(res_obj$TDiNextEst$size_based) %>%
    mutate(
      Order.q = factor(Order.q, levels = c(0, 1, 2), 
                       labels = c("Richness (q=0)", "Shannon (q=1)", "Simpson (q=2)")),
      Method = factor(Method, levels = c("Rarefaction", "Observed", "Extrapolation"))
    )

  ggplot(df, aes(x = m, y = qTD, color = Assemblage, fill = Assemblage)) +
    # OPTIMISATION OMBRES : 
    # 1. Augmentation de l'alpha (0.35)
    # 2. Ajout d'une bordure de ruban très fine (linewidth = 0.1) pour la structure
    geom_ribbon(aes(ymin = qTD.LCL, ymax = qTD.UCL, group = interaction(Assemblage, Order.q)), 
            alpha = 0.5,      # Augmente l'opacité
            color = "black",  # Ajoute un contour pour "marquer" l'ombre
            linewidth = 0.2,  # Très fin pour ne pas alourdir
            linetype = "dotted") +
    
    # Lignes de tendance
    geom_line(aes(linetype = Method), linewidth = 1.1) +
    
    # Points observés plus contrastés
    geom_point(data = filter(df, Method == "Observed"), size = 3.5, shape = 17) +
    
    facet_wrap(~ Order.q, scales = "free", ncol = 3) +
    scale_x_continuous(labels = label_number(scale = 1e-6, suffix = "M")) +
    scale_color_manual(values = palette_couleurs) +
    scale_fill_manual(values = palette_couleurs) +
    scale_linetype_manual(values = c("Rarefaction"="solid", "Extrapolation"="dashed", "Observed"="solid")) +
    
    theme_minimal(base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 14, margin = margin(b = 10)),
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "grey95", color = "grey90"),
      panel.border = element_rect(color = "grey90", fill = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92"), # Grille plus douce
      axis.title = element_text(face = "bold", size = 10),
      legend.key.width = unit(1.5, "cm")
    ) +
    labs(x = "Profondeur de séquençage (Millions de Lectures)", 
         y = "Nombres de Hill (Diversité)", 
         title = title_label)
}

# Génération des graphiques avec leurs palettes respectives
p_epi_final  <- plot_final_premium(res_epi, "A. Communautés Epiphytes", custom_colors)
p_endo_final <- plot_final_premium(res_endo, "B. Communautés Endophytes", custom_colors)
p_proj_final <- plot_final_premium(res_proj, "C. Communautés Projets", project_colors)

# Assemblage final
design_final <- (p_epi_final / p_endo_final / p_proj_final) + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = "bottom", 
        legend.box = "vertical",
        legend.margin = margin(t = 10),
        plot.margin = margin(15, 15, 15, 15))

print(design_final)











############## UpSetR ##############

ps_bin <- transform_sample_counts(ps_rare, function(x) ifelse(x > 0, 1, 0))

meta <- data.frame(sample_data(ps))
meta$category <- paste(
  meta$project,
  meta$age,
  meta$organ,
  meta$position,
  sep = " "
)

taxa_list <- split(rownames(meta), meta$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

upset(
  fromList(taxa_list),
  nsets = length(taxa_list),
  order.by = "freq",
  main.bar.color = "forestgreen",
  sets.bar.color = "steelblue"
)


meta <- data.frame(sample_data(ps))
meta <- subset(meta, project == "LOT1")
meta$category <- paste(meta$age, meta$organ, meta$position, sep = " ")

taxa_list <- split(rownames(meta), meta$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

upset(
  fromList(taxa_list),
  nsets = length(taxa_list),
  order.by = "freq",
  main.bar.color = "forestgreen",
  sets.bar.color = "steelblue"
)


meta <- data.frame(sample_data(ps))
meta <- subset(meta, project == "LOT2")
meta$category <- paste(meta$age, meta$organ, meta$position, sep = " ")

taxa_list <- split(rownames(meta), meta$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

upset(
  fromList(taxa_list),
  nsets = length(taxa_list),
  order.by = "freq",
  main.bar.color = "forestgreen",
  sets.bar.color = "steelblue"
)


meta <- data.frame(sample_data(ps))
meta <- subset(meta, project == "LOT3")
meta$category <- paste(meta$age, meta$organ, meta$position, sep = " ")

taxa_list <- split(rownames(meta), meta$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

upset(
  fromList(taxa_list),
  nsets = length(taxa_list),
  order.by = "freq",
  main.bar.color = "forestgreen",
  sets.bar.color = "steelblue"
)


library(UpSetR)
library(cowplot)
library(magick)
library(patchwork)

generate_upset_captured <- function(taxa_list, title_text = "") {
  tmp <- tempfile(fileext = ".png")
  png(tmp, width = 1000, height = 800, res = 120)

  print(upset(
    fromList(taxa_list),
    nsets = length(taxa_list),
    order.by = "freq",
    main.bar.color = "forestgreen",
    sets.bar.color = "steelblue"
  ))

  grid::grid.text(
    title_text,
    x = 0.5,
    y = 0.95,
    gp = grid::gpar(fontsize = 15, font = 2)
  )

  dev.off()

  img <- magick::image_read(tmp)
  return(cowplot::ggdraw() + cowplot::draw_image(img))
}


meta_all <- data.frame(sample_data(ps))
meta_all$category <- paste(
  meta_all$project,
  meta_all$age,
  meta_all$organ,
  meta_all$position,
  sep = " "
)
taxa_all <- split(rownames(meta_all), meta_all$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

meta_lot1 <- subset(data.frame(sample_data(ps)), project == "LOT1")
meta_lot1$category <- paste(
  meta_lot1$age,
  meta_lot1$organ,
  meta_lot1$position,
  sep = " "
)
taxa_lot1 <- split(rownames(meta_lot1), meta_lot1$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

meta_lot2 <- subset(data.frame(sample_data(ps)), project == "LOT2")
meta_lot2$category <- paste(
  meta_lot2$age,
  meta_lot2$organ,
  meta_lot2$position,
  sep = " "
)
taxa_lot2 <- split(rownames(meta_lot2), meta_lot2$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

meta_lot3 <- subset(data.frame(sample_data(ps)), project == "LOT3")
meta_lot3$category <- paste(
  meta_lot3$age,
  meta_lot3$organ,
  meta_lot3$position,
  sep = " "
)
taxa_lot3 <- split(rownames(meta_lot3), meta_lot3$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

p1 <- generate_upset_captured(taxa_all, "Global")
p2 <- generate_upset_captured(taxa_lot1, "Project LOT1")
p3 <- generate_upset_captured(taxa_lot2, "Project LOT2")
p4 <- generate_upset_captured(taxa_lot3, "Project LOT3")

(p1 + p2) / (p3 + p4)


############## Tableau abondances ##############




library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(matrixStats)

ps_filtered <- filter_taxa(
  ps,
  function(x) sum(x > 0) >= ceiling(0.01 * length(x)),
  TRUE
)

ps_rel <- transform_sample_counts(ps_filtered, function(x) x / sum(x))

ps_order <- tax_glom(ps_rel, taxrank = "gbr268_Order", NArm = FALSE)
df_order <- psmelt(ps_order)

df_agg <- df_order %>%
  group_by(gbr268_Order, gbr268_Phylum, gbr268_Class, project, organ, position) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  mutate(group = paste(project, organ, position, sep = "_"))

df_heatmap <- df_agg %>%
  dplyr::select(gbr268_Order, group, Abundance) %>%
  filter(!is.na(gbr268_Order)) %>%
  group_by(gbr268_Order, group) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = group,
    values_from = Abundance,
    values_fill = list(Abundance = 0)
  ) %>%
  as.data.frame()

rownames(df_heatmap) <- df_heatmap$gbr268_Order
df_heatmap$gbr268_Order <- NULL
mat <- as.matrix(df_heatmap)

pseudo_count <- ifelse(any(mat > 0), min(mat[mat > 0]) / 2, 1e-6)
mat_log <- log10(mat + pseudo_count)

scale_rows <- function(x) {
  row_means <- rowMeans(x, na.rm = TRUE)
  row_sds <- matrixStats::rowSds(x, na.rm = TRUE)
  row_sds[row_sds == 0 | is.na(row_sds)] <- 1
  sweep(sweep(x, 1, row_means, "-"), 1, row_sds, "/")
}

mat_scaled <- scale_rows(mat_log)
mat_scaled[is.na(mat_scaled)] <- 0

col_meta <- df_agg %>%
  dplyr::select(group, project, organ, position) %>%
  distinct() %>%
  arrange(match(group, colnames(mat_scaled)))

row_meta <- df_agg %>%
  dplyr::select(gbr268_Order, gbr268_Phylum, gbr268_Class) %>%
  distinct() %>%
  filter(!is.na(gbr268_Order)) %>%
  arrange(match(gbr268_Order, rownames(mat_scaled)))

ann_colors <- list(
  Project = c("LOT1" = "#7FB3D5", "LOT2" = "#D98880", "LOT3" = "#BDC3C7"),
  Organ = c(
    "root" = "#A93226",
    "old_leaf" = "#27AD60",
    "young_leaf" = "#A2D9CE"
  ),
  Position = c("endophyte" = "#FFDAB5", "epiphyte" = "#C0D461")
)

col_ha <- HeatmapAnnotation(
  Project = col_meta$project,
  Organ = col_meta$organ,
  Position = col_meta$position,
  col = ann_colors,
  annotation_name_side = "right"
)

row_ha <- rowAnnotation(
  Phylum = row_meta$gbr268_Phylum,
  Class = row_meta$gbr268_Class,
  annotation_name_side = "top"
)

max_abs <- max(abs(mat_scaled), na.rm = TRUE)
color_mapping <- colorRamp2(
  c(-max_abs, 0, max_abs),
  c("blue", "yellow", "red")
)
    

ordered_cols <- with(
  col_meta,
  order(
    factor(project, levels = c("LOT1", "LOT2", "LOT3")),
    factor(organ, levels = c("root", "old_leaf", "young_leaf")),
    factor(position, levels = c("epiphyte", "endophyte"))
  )
)

ht <- Heatmap(
  mat_scaled,
  name = "Row Z-score\n(log10 rel. abundance)",
  col = color_mapping,
  top_annotation = col_ha,
  left_annotation = row_ha,
  cluster_rows = TRUE,
  column_order = ordered_cols,
  cluster_columns = FALSE,
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  row_title = NULL,
  column_title = NULL,
  column_split = factor(col_meta$project[ordered_cols], levels = c("LOT1", "LOT2", "LOT3")),
  row_split = row_meta$gbr268_Phylum,
  heatmap_legend_param = list(direction = "vertical"),
  rect_gp = gpar(col = "grey85", lwd = 0.5),
  width = ncol(mat_scaled) * grid::unit(6, "mm"),
  height = nrow(mat_scaled) * grid::unit(4, "mm"),
  row_names_side = "right",
  column_names_side = "bottom",
  column_names_rot = 90,
  row_title_rot = 0
)


draw(ht, merge_legend = TRUE)

output_dir <- "~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

png(file.path(output_dir, "heatmap_orders.png"), width = 7000, height = 10000, res = 800)
draw(ht, merge_legend = TRUE)
dev.off()











############## Ordinations ##############

ps <- readRDS(
  "~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/ps_final.rds"
)

library(vegan)

ps_hel <- transform_sample_counts(ps, function(x) sqrt(x / sum(x)))
ps_hel <- prune_samples(sample_sums(ps_hel) > 0, ps_hel)
otu_hel <- as(otu_table(ps_hel), "matrix")
metadata <- as(sample_data(ps_hel), "data.frame")

bray_dist <- vegdist(otu_hel, method = "bray")

perm_global <- adonis2(
  bray_dist ~ project + organ + position,
  data = metadata,
  by = "margin",
  permutations = 999
)
nmds_global <- ordinate(ps_hel, method = "NMDS", distance = "bray", k = 3)

p1 <- plot_ordination(ps_hel, nmds_global, color = "project", shape = "organ") +
  stat_ellipse(
    aes(group = project, fill = project),
    alpha = 0.1,
    geom = "polygon"
  ) +
  labs(
    title = "A. Structure Globale (Projet)",
    subtitle = paste0("R2 Project: ", round(perm_global$R2[1], 3))
  ) +
  theme_bw()


ps_lot01 <- subset_samples(ps_hel, project == "LOT1")
ps_lot01 <- prune_taxa(taxa_sums(ps_lot01) > 0, ps_lot01)

meta_lot01 <- as(sample_data(ps_lot01), "data.frame")
valid_samples <- complete.cases(meta_lot01[, c("organ", "position")])
ps_lot01 <- prune_samples(valid_samples, ps_lot01)

bray_lot01 <- vegdist(otu_table(ps_lot01), method = "bray")
perm_lot01 <- adonis2(
  bray_lot01 ~ organ + position,
  data = as(sample_data(ps_lot01), "data.frame"),
  by = "margin",
  permutations = 999
)

nmds_lot01 <- ordinate(ps_lot01, method = "NMDS", distance = "bray", k = 3)

p2 = plot_ordination(
  ps_lot01,
  nmds_lot01,
  color = "organ",
  shape = "position"
) +
  geom_point(size = 3) +
  scale_shape_manual(values = c("epiphyte" = 1, "endophyte" = 16)) +
  stat_ellipse(aes(group = organ)) +
  labs(
    title = "B. LOT1 - Organe",
    subtitle = paste0(
      "R2 = ",
      round(perm_lot01["organ", "R2"], 3),
      ", p = ",
      perm_lot01["organ", "Pr(>F)"]
    )
  ) +
  theme_bw()


ps_lot02 <- subset_samples(ps_hel, project == "LOT2")
ps_lot02 <- prune_taxa(taxa_sums(ps_lot02) > 0, ps_lot02)

meta_lot02 <- as(sample_data(ps_lot02), "data.frame")
valid_samples <- complete.cases(meta_lot02[, c("organ", "position")])
ps_lot02 <- prune_samples(valid_samples, ps_lot02)


bray_lot02 <- vegdist(otu_table(ps_lot02), method = "bray")
vars_lot02 <- c("organ", "position")
valid_vars_lot02 <- vars_lot02[sapply(
  meta_lot02[valid_samples, vars_lot02],
  function(x) length(unique(x)) > 1
)]
formula_lot02 <- as.formula(paste(
  "bray_lot02 ~",
  paste(valid_vars_lot02, collapse = " + ")
))

perm_lot02 <- adonis2(
  formula_lot02,
  data = as(sample_data(ps_lot02), "data.frame"),
  by = "margin",
  permutations = 999
)

nmds_lot02 <- ordinate(ps_lot02, method = "NMDS", distance = "bray", k = 3)


p3 = plot_ordination(
  ps_lot02,
  nmds_lot02,
  color = "organ",
  shape = "position"
) +
  geom_point(size = 3) +
  scale_shape_manual(values = c("epiphyte" = 1, "endophyte" = 16)) +
  stat_ellipse(aes(group = organ)) +
  labs(
    title = "C. LOT2 - Organe",
    subtitle = paste0(
      "R2 = ",
      round(perm_lot02["organ", "R2"], 3),
      ", p = ",
      perm_lot02["organ", "Pr(>F)"]
    )
  ) +
  theme_bw()


ps_lot03 <- subset_samples(ps_hel, project == "LOT3")
ps_lot03 <- prune_taxa(taxa_sums(ps_lot03) > 0, ps_lot03)

bray_lot03 <- vegdist(otu_table(ps_lot03), method = "bray")
perm_lot03 <- adonis2(
  bray_lot03 ~ organ + position,
  data = as(sample_data(ps_lot03), "data.frame"),
  by = "margin",
  permutations = 999
)

nmds_lot03 <- ordinate(ps_lot03, method = "NMDS", distance = "bray", k = 3)

p4 = plot_ordination(
  ps_lot03,
  nmds_lot03,
  color = "organ",
  shape = "position"
) +
  geom_point(size = 3) +
  scale_shape_manual(values = c("epiphyte" = 1, "endophyte" = 16)) +
  stat_ellipse(aes(group = organ)) +
  labs(
    title = "D. LOT3 - Organe",
    subtitle = paste0(
      "R2 = ",
      round(perm_lot03["organ", "R2"], 3),
      ", p = ",
      perm_lot03["organ", "Pr(>F)"]
    )
  ) +
  theme_bw()


vars_pop <- c("project", "organ", "position")
form_pop <- as.formula(paste("bray_dist ~", paste(vars_pop, collapse = " + ")))
perm_pop <- adonis2(
  form_pop,
  data = metadata,
  by = "margin",
  permutations = 999
)

tab_pop <- data.frame(
  Facteur = vars_pop,
  Df = perm_pop[vars_pop, "Df"],
  R2 = round(perm_pop[vars_pop, "R2"], 4),
  p_value = perm_pop[vars_pop, "Pr(>F)"]
)
rownames(tab_pop) <- NULL

print(perm_pop)
print(tab_pop)

organ_colors <- c(
  "root" = "#794515",
  "old_leaf" = "forestgreen",
  "young_leaf" = "#4DCCBD"
)
p2 <- p2 + scale_color_manual(values = organ_colors)
p3 <- p3 + scale_color_manual(values = organ_colors)
p4 <- p4 + scale_color_manual(values = organ_colors)

(p1 | p2) / (p3 | p4)


fmt_p <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    ifelse(x < 0.001, "<0.001", formatC(x, format = "f", digits = 3))
  )
}

sig_code <- function(x) {
  dplyr::case_when(
    is.na(x) ~ "",
    x < 0.001 ~ "***",
    x < 0.01 ~ "**",
    x < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

extract_main_effect <- function(adonis_obj, term, panel, figure_title, nmds_obj, ps_obj) {
  adf <- as.data.frame(adonis_obj)
  adf$Term <- rownames(adf)

  if (!term %in% adf$Term) {
    return(
      data.frame(
        Panel = panel,
        Figure = figure_title,
        Effect = term,
        Samples = phyloseq::nsamples(ps_obj),
        Df = NA,
        R2 = NA,
        F = NA,
        p_value = NA_character_,
        Significance = "",
        NMDS_stress = round(nmds_obj$stress, 3),
        stringsAsFactors = FALSE
      )
    )
  }

  row_i <- adf[adf$Term == term, , drop = FALSE]
  p_raw <- row_i[["Pr(>F)"]]

  data.frame(
    Panel = panel,
    Figure = figure_title,
    Effect = term,
    Samples = phyloseq::nsamples(ps_obj),
    Df = row_i$Df,
    R2 = round(row_i$R2, 3),
    F = round(row_i$F, 3),
    p_value = fmt_p(p_raw),
    Significance = sig_code(p_raw),
    NMDS_stress = round(nmds_obj$stress, 3),
    stringsAsFactors = FALSE
  )
}

tableau_ordination_publication <- dplyr::bind_rows(
  extract_main_effect(
    perm_global, "project",
    "A", "Structure globale (Projet)",
    nmds_global, ps_hel
  ),
  extract_main_effect(
    perm_lot01, "organ",
    "B", "LOT1 - Organe",
    nmds_lot01, ps_lot01
  ),
  extract_main_effect(
    perm_lot02, "organ",
    "C", "LOT2 - Organe",
    nmds_lot02, ps_lot02
  ),
  extract_main_effect(
    perm_lot03, "organ",
    "D", "LOT3 - Organe",
    nmds_lot03, ps_lot03
  )
)

print(tableau_ordination_publication)

write.csv(
  tableau_ordination_publication,
  "~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/Tableau_resume_ordinations_publication.csv",
  row.names = FALSE
)

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

# Calcul de la niche de Levins
ps_merged <- merge_samples(ps, "organ")
otu_merged <- as.matrix(otu_table(ps_merged))
if (taxa_are_rows(ps_merged)) {
  otu_merged <- t(otu_merged)
}
levins_raw <- niche.width(otu_merged, method = "levins")
niche_scores <- (as.numeric(levins_raw) - 1) / (3 - 1)
names(niche_scores) <- colnames(levins_raw)

ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
otu_all <- as.matrix(otu_table(ps_rel))
if (taxa_are_rows(ps_rel)) {
  otu_all <- t(otu_all)
}

calc_mean_niche <- function(sample_row, scores) {
  present <- sample_row > 0
  if (sum(present) == 0) return(NA)
  sub_scores <- scores[names(sample_row)[present]]
  sub_abund <- sample_row[present]
  sum(sub_abund * sub_scores, na.rm = TRUE) / sum(sub_abund, na.rm = TRUE)
}

sample_indices <- apply(otu_all, 1, calc_mean_niche, scores = niche_scores)

metadata <- data.frame(sample_data(ps))
final_df <- data.frame(
  Sample = names(sample_indices),
  Levins_Mean = as.numeric(sample_indices)
) %>%
  mutate(
    organ = metadata$organ[match(Sample, rownames(metadata))],
    project = metadata$project[match(Sample, rownames(metadata))]
  )

niche_summary <- final_df %>%
  filter(!is.nan(Levins_Mean) & !is.na(Levins_Mean)) %>%
  group_by(organ) %>%
  summarise(
    Moyenne = mean(Levins_Mean),
    SD = sd(Levins_Mean),
    n = n()
  )

print(niche_summary)

# Tests Kruskal-Wallis
kruskal_organ_by_lot <- final_df %>%
  filter(!is.na(Levins_Mean), !is.na(organ), !is.na(project)) %>%
  group_by(project) %>%
  summarise(
    p_value = kruskal.test(Levins_Mean ~ organ)$p.value,
    .groups = "drop"
  )
print("Kruskal-Wallis par lot (organes):")
print(kruskal_organ_by_lot)

kruskal_lot_by_organ <- final_df %>%
  filter(!is.na(Levins_Mean), !is.na(organ), !is.na(project)) %>%
  group_by(organ) %>%
  summarise(
    p_value = kruskal.test(Levins_Mean ~ project)$p.value,
    .groups = "drop"
  )
print("Kruskal-Wallis par organe (lots):")
print(kruskal_lot_by_organ)

# Post-hoc Dunn test si significatif
if (!requireNamespace("FSA", quietly = TRUE)) install.packages("FSA")
library(FSA)

dunn_organ_by_lot <- final_df %>%
  filter(!is.na(Levins_Mean), !is.na(organ), !is.na(project)) %>%
  group_by(project) %>%
  group_modify(~ {
    if (length(unique(.x$organ)) > 1) {
      dunn <- FSA::dunnTest(Levins_Mean ~ organ, data = .x, method = "bonferroni")
      data.frame(dunn$res, project = unique(.x$project))
    } else {
      data.frame()
    }
  }) %>% ungroup()

dunn_lot_by_organ <- final_df %>%
  filter(!is.na(Levins_Mean), !is.na(organ), !is.na(project)) %>%
  group_by(organ) %>%
  group_modify(~ {
    if (length(unique(.x$project)) > 1) {
      dunn <- FSA::dunnTest(Levins_Mean ~ project, data = .x, method = "bonferroni")
      data.frame(dunn$res, organ = unique(.x$organ))
    } else {
      data.frame()
    }
  }) %>% ungroup()

print("Dunn post-hoc organes par lot :")
print(dunn_organ_by_lot)
print("Dunn post-hoc lots par organe :")
print(dunn_lot_by_organ)

final_df$project <- as.factor(final_df$project)
final_df$organ <- factor(final_df$organ, levels = c("root", "young_leaf", "old_leaf"))

# Préparation des comparaisons pour stat_pvalue_manual
# On prépare une table pour chaque projet avec les comparaisons et p-values
comparisons_organ_by_lot <- dunn_organ_by_lot %>%
  filter(!is.na(P.adj)) %>%
  mutate(
    group1 = gsub(" vs .*", "", Comparison),
    group2 = gsub(".* vs ", "", Comparison),
    y.position = tapply(final_df$Levins_Mean, final_df$project, max, na.rm = TRUE)[project] + 0.05,
    p.adj.signif = symnum(P.adj, corr = FALSE, na = FALSE,
                          cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                          symbols = c("***", "**", "*", "ns"))
  )

# Pour éviter les superpositions, on décale les y.position pour chaque comparaison
comparisons_organ_by_lot <- comparisons_organ_by_lot %>%
  group_by(project) %>%
  mutate(y.position = y.position + row_number() * 0.03) %>%
  ungroup()

# Plot
p_project_organ <- ggplot(
  final_df,
  aes(x = project, y = Levins_Mean, fill = organ)
) +
  geom_boxplot(
    alpha = 0.6,
    outlier.shape = NA,
    width = 0.6,
    position = position_dodge(0.8)
  ) +
  geom_jitter(
    aes(color = organ),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    alpha = 0.5,
    size = 1
  ) +
  labs(
    title = "Niche fongique par projet et organe",
    x = "Projet",
    y = "Indice de Levins Moyen",
    fill = "Organe",
    color = "Organe"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11, face = "bold"),
    title = element_text(size = 13, face = "bold")
  ) +
  scale_fill_manual(
    values = c(
      "root" = "#914e27",
      "young_leaf" = "#83d483",
      "old_leaf" = "#1d5d1d"
    )
  ) +
  scale_color_manual(
    values = c(
      "root" = "#914e27",
      "young_leaf" = "#83d483",
      "old_leaf" = "#1d5d1d"
    )
  )

# Ajout des étoiles de significativité avec stat_pvalue_manual
if (nrow(comparisons_organ_by_lot) > 0) {
  p_project_organ <- p_project_organ +
    stat_pvalue_manual(
      data = comparisons_organ_by_lot,
      label = "p.adj.signif",
      y.position = "y.position",
      xmin = "group1",
      xmax = "group2",
      group = "project",
      tip.length = 0.01,
      step.increase = 0,
      hide.ns = TRUE
    )
}

print(p_project_organ)
