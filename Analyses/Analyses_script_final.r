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


nettoyage_otu <- function(ps_obj, variable) {
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


mat_epi <- nettoyage_otu(subset_samples(ps_trim, position == "epiphyte"), "organ")
res_epi <- iNEXT3D(
  mat_epi,
  q = c(0, 1, 2),
  datatype = "abundance",
  nboot = 50,
  size = seq(10, max(1.5 * colSums(mat_epi)), length.out = 40)
)
saveRDS(res_epi, file = "res_epi.rds")

mat_endo <- nettoyage_otu(subset_samples(ps_trim, position == "endophyte"), "organ")
res_endo <- iNEXT3D(
  mat_endo,
  q = c(0, 1, 2),
  datatype = "abundance",
  nboot = 50,
  size = seq(10, max(1.5 * colSums(mat_endo)), length.out = 40)
)
saveRDS(res_endo, file = "res_endo.rds")

mat_proj <- nettoyage_otu(ps_trim, "project")
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





couleurs_organes <- c("root" = "#8B4513", "old_leaf" = "#2E7D32", "young_leaf" = "#81C784")
couleurs_projets <- c("LOT1" = "#1F78B4", "LOT2" = "#E31A1C", "LOT3" = "#33A02C")

plot_final_premium <- function(res_obj, title_label, palette_couleurs) {
  
  df <- as.data.frame(res_obj$TDiNextEst$size_based) %>%
    mutate(
      Order.q = factor(Order.q, levels = c(0, 1, 2), 
                       labels = c("Richness (q=0)", "Shannon (q=1)", "Simpson (q=2)")),
      Method = factor(Method, levels = c("Rarefaction", "Observed", "Extrapolation"))
    )

  ggplot(df, aes(x = m, y = qTD, color = Assemblage, fill = Assemblage)) +
    geom_ribbon(aes(ymin = qTD.LCL, ymax = qTD.UCL, group = interaction(Assemblage, Order.q)), 
            alpha = 0.5,      
            color = "black",  
            linewidth = 0.2,  
            linetype = "dotted") +
    
    geom_line(aes(linetype = Method), linewidth = 1.1) +
    
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
      panel.grid.major = element_line(color = "grey92"), 
      axis.title = element_text(face = "bold", size = 10),
      legend.key.width = unit(1.5, "cm")
    ) +
    labs(x = "Profondeur de séquençage (Millions de Lectures)", 
         y = "Nombres de Hill (Diversité)", 
         title = title_label)
}

p_epi_final  <- plot_final_premium(res_epi, "A. Communautés Epiphytes", couleurs_organes)
p_endo_final <- plot_final_premium(res_endo, "B. Communautés Endophytes", couleurs_organes)
p_proj_final <- plot_final_premium(res_proj, "C. Communautés Projets", couleurs_projets)

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

assemblage_des_upset <- function(taxa_list, title_text = "") {
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

p1 <- assemblage_des_upset(taxa_all, "Global")
p2 <- assemblage_des_upset(taxa_lot1, "Project LOT1")
p3 <- assemblage_des_upset(taxa_lot2, "Project LOT2")
p4 <- assemblage_des_upset(taxa_lot3, "Project LOT3")

(p1 + p2) / (p3 + p4)





















############## Tableau abondances ##############























library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(matrixStats)

ps_filtrés <- filter_taxa(
  ps,
  function(x) sum(x > 0) >= ceiling(0.01 * length(x)),
  TRUE
)

ps_rel <- transform_sample_counts(ps_filtrés, function(x) x / sum(x))

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




















############# Abondances relatives classes #############















library(scales)



ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
ps_class <- tax_glom(ps_rel, taxrank = "gbr268_Class", NArm = FALSE)
df_class <- psmelt(ps_class)
df_class$organ_position <- interaction(df_class$organ, df_class$position, sep = " - ")

mean_abundances <- df_class %>%
  group_by(gbr268_Class) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE))

low_classes <- mean_abundances$gbr268_Class[mean_abundances$mean_abundance < 0.01]

df_class$Class_grouped <- ifelse(df_class$gbr268_Class %in% low_classes, "Autres", as.character(df_class$gbr268_Class))



df_plot <- df_class %>%
  group_by(plant_family, organ_position, Class_grouped) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = 'drop')

classes_bio <- setdiff(unique(df_plot$Class_grouped), c("NA", "Autres"))

df_plot$Class_grouped <- factor(df_plot$Class_grouped, 
                                levels = c("NA", "Autres", classes_bio))

ggplot(df_plot, aes(x = plant_family, y = Abundance, fill = Class_grouped)) +
  geom_col(position = "fill", color = NA) + 
  facet_wrap(~ organ_position, scales = "free_x", ncol = 3) +
  theme_bw() +
  labs(
    title = "Abondances relatives des classes fongiques",
    x = "Famille de plante", 
    y = "Abondance relative", 
    fill = "Classe fongique"
  ) +
  theme(
    legend.position = "bottom", 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    values = c(
      "NA" = "#707070",      
      "Autres" = "#D3D3D3",  
      setNames(
        hue_pal()(length(classes_bio)),
        classes_bio
      )
    )
  )












  ###################### Analyse des classes les plus abondantes et différences entre groupes ######################















  df_class$leaf_group <- ifelse(df_class$organ == "root", "root", "leaf")

  top_classes <- df_class %>%
    group_by(leaf_group, gbr268_Class) %>%
    summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(gbr268_Class)) %>%
    group_by(leaf_group) %>%
    slice_max(order_by = mean_abundance, n = 2) %>%
    arrange(leaf_group, desc(mean_abundance))

  print(top_classes)

  df_class$leaf_group <- ifelse(df_class$organ == "root", "root", "leaf")

  top_classes_by_group <- df_class %>%
    group_by(leaf_group, position, gbr268_Class) %>%
    summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(gbr268_Class)) %>%
    group_by(leaf_group, position) %>%
    slice_max(order_by = mean_abundance, n = 2) %>%
    arrange(leaf_group, position, desc(mean_abundance))

  print(top_classes_by_group)


  top_classes_leaves <- df_class %>%
    filter(organ %in% c("old_leaf", "young_leaf")) %>%
    group_by(organ, gbr268_Class) %>%
    summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(gbr268_Class)) %>%
    group_by(organ) %>%
    slice_max(order_by = mean_abundance, n = 2) %>%
    arrange(organ, desc(mean_abundance))

  print(top_classes_leaves)


  library(dplyr)

  df_leaves <- df_class %>%
    filter(organ %in% c("young_leaf", "old_leaf")) %>%
    filter(!is.na(gbr268_Class))

  diff_classes <- df_leaves %>%
    group_by(gbr268_Class, organ) %>%
    summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = organ, values_from = mean_abundance, values_fill = 0) %>%
    mutate(diff = young_leaf - old_leaf) %>%
    arrange(desc(abs(diff)))

  print(diff_classes)



  library(dplyr)

  df_leaves <- df_class %>%
    filter(organ %in% c("young_leaf", "old_leaf")) %>%
    filter(!is.na(gbr268_Class))

  diff_classes_position <- df_leaves %>%
    group_by(gbr268_Class, position) %>%
    summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = position, values_from = mean_abundance, values_fill = 0) %>%
    mutate(diff = endophyte - epiphyte) %>%
    arrange(desc(abs(diff)))

  print(diff_classes_position)




  library(dplyr)

  df_roots <- df_class %>%
    filter(organ == "root") %>%
    filter(!is.na(gbr268_Class))

  diff_classes_position_roots <- df_roots %>%
    group_by(gbr268_Class, position) %>%
    summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = position, values_from = mean_abundance, values_fill = 0) %>%
    mutate(diff = endophyte - epiphyte) %>%
    arrange(desc(abs(diff)))

  print(diff_classes_position_roots)

  df_class$leaf_group <- ifelse(df_class$organ == "root", "root", "leaf")

  diff_classes_root_leaf <- df_class %>%
    filter(!is.na(gbr268_Class)) %>%
    group_by(gbr268_Class, leaf_group) %>%
    summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = leaf_group, values_from = mean_abundance, values_fill = 0) %>%
    mutate(diff = root - leaf) %>%
    arrange(desc(abs(diff)))


  print(top_classes)
  print(top_classes_by_group)
  print(top_classes_leaves)
  print(diff_classes)
  print(diff_classes_position)
  print(diff_classes_position_roots)
  print(diff_classes_root_leaf)


  abund_by_family <- df_class %>%
    group_by(plant_family, gbr268_Class) %>%
    summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(gbr268_Class))

  library(pheatmap)
  mat_family_class <- abund_by_family %>%
    tidyr::pivot_wider(names_from = plant_family, values_from = mean_abundance, values_fill = 0) %>%
    as.data.frame()
  rownames(mat_family_class) <- mat_family_class$gbr268_Class
  mat_family_class$gbr268_Class <- NULL
  pheatmap(as.matrix(mat_family_class), 
           cluster_rows = TRUE, cluster_cols = TRUE, 
           main = "Abondance moyenne des classes fongiques par famille de plante")

  kruskal_results <- abund_by_family %>%
    group_by(gbr268_Class) %>%
    summarise(
      p.value = kruskal.test(mean_abundance ~ plant_family)$p.value,
      .groups = "drop"
    ) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
    arrange(p.adj)

  print(kruskal_results)

  signif_classes <- kruskal_results %>% filter(p.adj < 0.05)
  print(signif_classes)




















# %% [Ordinations] 



















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













# %% [Pairwise adonis et tableau de synthèse pour publication] 















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


perm_global <- adonis2(
  bray_dist ~ project + organ + position + plant_family,
  data = metadata,
  by = "margin",
  permutations = 999
)



tab_perm_global <- data.frame(
  Facteur = rownames(perm_global)[1:4],
  Df = perm_global$Df[1:4],
  R2 = round(perm_global$R2[1:4], 4),
  F = round(perm_global$F[1:4], 3),
  p_value = perm_global$`Pr(>F)`[1:4]
)

tab_perm_global$p_value_fmt <- ifelse(
  is.na(tab_perm_global$p_value),
  NA_character_,
  ifelse(tab_perm_global$p_value < 0.001, "<0.001", formatC(tab_perm_global$p_value, format = "f", digits = 3))
)

tab_perm_global$Significance <- dplyr::case_when(
  is.na(tab_perm_global$p_value) ~ "",
  tab_perm_global$p_value < 0.001 ~ "***",
  tab_perm_global$p_value < 0.01 ~ "**",
  tab_perm_global$p_value < 0.05 ~ "*",
  TRUE ~ "ns"
)

print(tab_perm_global)

write.csv(
  tab_perm_global,
  "~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/Tableau_resume_perm_global_publication.csv",
  row.names = FALSE
)












library(pairwiseAdonis)

# %% [Comparaisons filtrées des communautés fongiques (Pairwise Adonis)]

metadata$proj_org_pos <- interaction(metadata$project, metadata$organ, metadata$position, sep = "_", drop = TRUE)
mask <- metadata$proj_org_pos %in% names(which(table(metadata$proj_org_pos) >= 3))
pw_projorgpos <- pairwise.adonis(
  as.dist(as.matrix(bray_dist)[mask, mask]),
  metadata$proj_org_pos[mask],
  p.adjust.m = "fdr"
)
pw_projorgpos <- pw_projorgpos[order(pw_projorgpos$R2, decreasing = TRUE), ]




# %% [Comparaisons par paires des communautés fongiques entre familles de plantes hôtes (Pairwise Adonis)]

pw_plant_family <- pairwise.adonis(bray_dist, metadata$plant_family, p.adjust.m = "fdr")
print(pw_plant_family)



# %% [Comparaisons par paires des communautés fongiques entre familles de plantes hôtes au sein de chaque projet (Pairwise Adonis)]

metadata$family_project <- interaction(metadata$plant_family, metadata$project, drop = TRUE)
fam_proj <- metadata %>%
  count(plant_family, project) %>%
  count(plant_family) %>%
  filter(n >= 2) %>%
  pull(plant_family)
meta_fp <- metadata %>% filter(plant_family %in% fam_proj)

bray_fp <- as.dist(as.matrix(bray_dist)[rownames(meta_fp), rownames(meta_fp)])
meta_fp$family_project <- interaction(meta_fp$plant_family, meta_fp$project, drop = TRUE)
pw_family_project_all <- pairwise.adonis(
  bray_fp,
  meta_fp$family_project,
  p.adjust.m = "fdr"
)
















########### Niches de Levins ###########






















library(dplyr)
library(tidyr)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(spaa)


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
    organ = factor(metadata$organ[match(Sample, rownames(metadata))], 
                   levels = c("root", "young_leaf", "old_leaf")),
    project = as.factor(metadata$project[match(Sample, rownames(metadata))])
  ) %>%
  filter(!is.na(Levins_Mean))





stat_df <- final_df %>%
  group_by(project) %>%
  wilcox_test(Levins_Mean ~ organ) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_y_position(step.increase = 0.1)


  stat_df_clean <- stat_df %>%
    mutate(across(where(is.list), ~ sapply(., function(x) paste(x, collapse = "; "))))
  
  write.csv(
    stat_df_clean,
    "~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/stat_df_levins_niche.csv",
    row.names = FALSE
  )


p <- ggplot(final_df, aes(x = organ, y = Levins_Mean, fill = organ)) +
  geom_boxplot(
    alpha = 0.7,
    outlier.shape = NA,
    width = 0.6
  ) +
  geom_jitter(
    aes(color = organ),
    width = 0.15,
    alpha = 0.4,
    size = 1.2
  ) +
  stat_pvalue_manual(
    stat_df, 
    label = "p.adj.signif", 
    tip.length = 0.02,
    hide.ns = TRUE
  ) +
  facet_wrap(~project, strip.position = "bottom") +
  scale_fill_manual(values = c("root" = "#914e27", "young_leaf" = "#83d483", "old_leaf" = "#1d5d1d")) +
  scale_color_manual(values = c("root" = "#914e27", "young_leaf" = "#83d483", "old_leaf" = "#1d5d1d")) +
  labs(
    title = "Niche fongique par projet et organe",
    x = "Sites",
    y = "Indice de Levins Moyen"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_blank(), # On cache "root/leaf" en bas car on a la légende et les facettes
    axis.ticks.x = element_blank(),
    panel.spacing = unit(1, "lines") # Espace entre les LOTS
  )

print(p)















############################### Autres #########################




















##### Nombre de plantes par famille et par projet
meta <- data.frame(sample_data(ps))
cols_taxo <- c("plant_family")

if (nrow(meta) > 0) {
  meta$plant_id <- sub("^.{6}([0-9]{4}).*$", "\\1", as.character(meta$LOT_sampleID))
} else {
  stop("Le data.frame 'meta' est vide. Vérifiez que l'objet 'ps' contient des données.")
}

df_plantes <- meta %>%
  group_by(across(all_of(cols_taxo))) %>%
  summarise(n_plantes = n_distinct(plant_id), .groups = "drop")

df_plantes_lot <- meta %>%
  group_by(across(all_of(cols_taxo)), project) %>%
  summarise(n_plantes = n_distinct(plant_id), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = project,
    values_from = n_plantes,
    values_fill = 0,
    names_prefix = "n_plantes_"
  )

df_plantes_lot$total <- rowSums(df_plantes_lot[grep("^n_plantes_", names(df_plantes_lot))])
df_plantes_final <- df_plantes_lot


df_nb_plantes_par_famille_projet <- meta %>%
  group_by(plant_family, project) %>%
  summarise(nb_plantes_uniques = n_distinct(plant_id), .groups = "drop")

print(df_nb_plantes_par_famille_projet, n = Inf)

for (lot in c("LOT1", "LOT2", "LOT3")) {
  colname <- paste0("n_plantes_", lot)
  if (!colname %in% names(df_plantes_lot)) {
    df_plantes_lot[[colname]] <- 0
  }
}

cols_final <- c(cols_taxo, paste0("n_plantes_LOT", 1:3), "total")
df_tab_final <- df_plantes_lot[, cols_final, drop = FALSE]

df_tab_final[is.na(df_tab_final)] <- 0

print(df_tab_final, n = Inf)

