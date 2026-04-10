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


######## Tests et visualisations pour les alpha-diversités ########

valeurs_p_kruskal <- lapply(
  df_alpha[c("Observed", "expShannon")], 
  function(metrique) kruskal.test(metrique ~ df_alpha$grp)$p.value
)

plot_sig_heatmap <- function(donnees, metrique, titre) {
  

  resultats_wilcoxon <- pairwise.wilcox.test(
    donnees[[metrique]], 
    donnees$grp, 
    p.adjust.method = "BH"
  )$p.value
  
  donnees_graphique <- resultats_wilcoxon |>
    as.table() |>
    as.data.frame() |>
    na.omit()
  
  ggplot(donnees_graphique, aes(x = Var1, y = Var2, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = cut(Freq, 
                              breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                              labels = c("***", "**", "*", "ns")))) +
    scale_fill_gradient2(low = "firebrick", mid = "red", high = "white", midpoint = 0.05) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = titre, x = "Groupes", y = "Groupes comparés", fill = "Valeur p")
}

plot_sig_heatmap(df_alpha, "Observed", "Significativité - Richesse Observée")

plot_sig_heatmap(df_alpha, "expShannon", "Significativité - ExpShannon")

nb_otus <- nrow(tax_table(ps))
cat("Nombre d'OTUs différents :", nb_otus, "\n")


plot_groupes <- function(df, x, y, facet_by) {
    ggplot(df, aes(x = .data[[x]], y = .data[[y]], fill = .data[[x]])) +
        geom_violin(alpha = 0.7, trim = FALSE) +
        geom_pwc(method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
        facet_wrap(vars(.data[[facet_by]])) + 
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
}


wrap_plots(
    plot_groupes(df_alpha, "project", "Observed", "organ"),
    plot_groupes(df_alpha, "position", "Observed", "organ"),
    plot_groupes(df_alpha, "project", "Observed", "position"),
    ncol = 3
) + plot_annotation(title = "Richesse Spécifique par Catégorie")

wrap_plots(
    plot_groupes(df_alpha, "project", "expShannon", "organ"),
    plot_groupes(df_alpha, "position", "expShannon", "organ"),
    plot_groupes(df_alpha, "project", "expShannon", "position"),
    ncol = 3
) + plot_annotation(title = "ExpShannon par Catégorie")

wrap_plots(
    plot_groupes(df_alpha, "organ", "Observed", "project"),
    plot_groupes(df_alpha, "organ", "Observed", "position"),
    plot_groupes(df_alpha, "position", "Observed", "project"),
    ncol = 3
) + plot_annotation(title = "Richesse Spécifique par Catégorie")

wrap_plots(
    plot_groupes(df_alpha, "organ", "expShannon", "project"),
    plot_groupes(df_alpha, "organ", "expShannon", "position"),
    plot_groupes(df_alpha, "position", "expShannon", "project"),
    ncol = 3
) + plot_annotation(title = "ExpShannon par Catégorie")



# Courbes de raréfaction en fct du nombre de sequences lues







if (!requireNamespace("iNEXT", quietly = TRUE)) {
  install.packages("iNEXT")
}
library(iNEXT)

otu_mat <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)

meta <- data.frame(sample_data(ps))
meta <- meta[!is.na(meta$organ), , drop = FALSE]
otu_mat <- otu_mat[rownames(meta), , drop = FALSE]

org_idx <- split(seq_len(nrow(otu_mat)), meta$organ)

abund_list <- lapply(org_idx, function(i) {
  x <- colSums(otu_mat[i, , drop = FALSE])
  x[x > 0]
})

rich_ine <- iNEXT(abund_list, q = 0, datatype = "abundance")

p_rare <- ggiNEXT(rich_ine, type = 1, facet.var = "site", color.var = "site") +
  labs(
    x = "Nombre de séquences lues",
    y = "Richesse spécifique (q = 0)",
    color = "Organe",
    title = "Courbes de raréfaction / extrapolation de la richesse spécifique par organe"
  ) +
  theme_bw()

print(p_rare)





































library(vegan)
library(ggplot2)
library(dplyr)

get_shannon_acc <- function(ps_obj, target_organ_name) {
  indices <- which(sample_data(ps_obj)$organ == target_organ_name)
  
  otu_full <- as(otu_table(ps_obj), "matrix")
  if (taxa_are_rows(ps_obj)) { otu_full <- t(otu_full) }
  
  otu_sub <- otu_full[indices, , drop = FALSE]
  n_samples <- nrow(otu_sub)
  
  res <- replicate(20, {
    sapply(seq(1, n_samples, length.out = 30), function(n) { 
      n <- round(n)
      if (n == 1) {
        idx <- sample(1:n_samples, 1)
        shan <- diversity(otu_sub[idx, , drop=FALSE], index = "shannon")
      } else {
        idx <- sample(1:n_samples, n)
        combined_otu <- colSums(otu_sub[idx, , drop=FALSE])
        shan <- diversity(combined_otu, index = "shannon")
      }
      return(exp(shan)) 
    })
  })
  
  return(data.frame(
    Samples = seq(1, n_samples, length.out = 30),
    ExpShannon = rowMeans(res),
    SD = apply(res, 1, sd),
    Organe = target_organ_name
  ))
}

organs <- unique(as.character(sample_data(ps)$organ))
df_shannon_acc <- do.call(rbind, lapply(organs, function(x) get_shannon_acc(ps, x)))

ggplot(df_shannon_acc, aes(x = Samples, y = ExpShannon, color = Organe, fill = Organe)) +
  geom_ribbon(aes(ymin = ExpShannon - SD, ymax = ExpShannon + SD), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Diversité Cumulative (Nombre d'espèces effectives)",
       subtitle = "Accumulation de l'Exp(Shannon) : Stabilité de la structure des communautés",
       x = "Nombre d'échantillons (Plantes)",
       y = "Exp(Shannon) global") +
  theme_bw()












# avec la fonction native de phyloseq
alpha_div <- c("Observed", "Chao1", "Shannon", "Simpson")
p <- plot_richness(ps, x = "project", measures = alpha_div, color = "organ")
p + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))





# Diagramme de Venn 

extraire_taxons <- function(physeq, variable) {
  meta <- as.data.frame(sample_data(physeq))
  groupes <- split(rownames(meta), meta[[variable]])
  
  lapply(groupes, function(echantillons) {
    sous_ps <- prune_samples(echantillons, physeq)
    abondances <- taxa_sums(sous_ps)
    names(abondances[abondances > 0])
  })
}

taxons_par_projet <- extraire_taxons(ps, "project")

ggVennDiagram(taxons_par_projet) +
  scale_fill_gradient(low = "white", high = "forestgreen") +
  labs(title = "Taxons partagés entre les projets LOT",
       fill = "Nombre de taxons") 

taxons_par_organe <- extraire_taxons(ps, "organ")

ggVennDiagram(taxons_par_organe) +
  scale_fill_gradient(low = "white", high = "forestgreen") +
  labs(title = "Taxons partagés entre les organes",
       fill = "Nombre de taxons")



# Diagramme d'UpSet
meta <- data.frame(sample_data(ps))
meta$category <- paste(meta$project, meta$age, meta$organ, meta$position, sep = " ")

taxa_list <- split(rownames(meta), meta$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

upset(fromList(taxa_list), nsets = length(taxa_list), order.by = "freq", 
      main.bar.color = "forestgreen", sets.bar.color = "steelblue")

input_upset <- fromList(taxa_list)
input_upset <- fromList(taxa_list)

proj_cols <- c(
  LOT1 = "#1b9e77",
  LOT2 = "#d95f02",
  LOT3 = "#7570b3"
)

# Colonnes correspondant aux catégories des projets
lot1_sets <- grep("^LOT1", colnames(input_upset), value = TRUE)
lot2_sets <- grep("^LOT2", colnames(input_upset), value = TRUE)
lot3_sets <- grep("^LOT3", colnames(input_upset), value = TRUE)

# TRUE si au moins une catégorie du projet est présente dans l'intersection
lot1_query <- function(row, data) any(row[lot1_sets] == 1)
lot2_query <- function(row, data) any(row[lot2_sets] == 1)
lot3_query <- function(row, data) any(row[lot3_sets] == 1)

upset(
  input_upset,
  nsets = ncol(input_upset),
  order.by = "freq",
  queries = list(
    list(
      query      = lot1_query,
      params     = list(),          # params explicitement défini
      color      = proj_cols["LOT1"],
      active     = TRUE,
      query.name = "LOT1"
    ),
    list(
      query      = lot2_query,
      params     = list(),
      color      = proj_cols["LOT2"],
      active     = TRUE,
      query.name = "LOT2"
    ),
    list(
      query      = lot3_query,
      params     = list(),
      color      = proj_cols["LOT3"],
      active     = TRUE,
      query.name = "LOT3"
    )
  ),
  main.bar.color = "grey"
)





meta <- data.frame(sample_data(ps))
meta <- subset(meta, project == "LOT1")
meta$category <- paste(meta$age, meta$organ, meta$position, sep = " ")

taxa_list <- split(rownames(meta), meta$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

p3 =upset(fromList(taxa_list), nsets = length(taxa_list), order.by = "freq", 
      main.bar.color = "forestgreen", sets.bar.color = "steelblue")



meta <- data.frame(sample_data(ps))
meta <- subset(meta, project == "LOT2")
meta$category <- paste(meta$age, meta$organ, meta$position, sep = " ")

taxa_list <- split(rownames(meta), meta$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

p2=upset(fromList(taxa_list), nsets = length(taxa_list), order.by = "freq", 
      main.bar.color = "forestgreen", sets.bar.color = "steelblue")




meta <- data.frame(sample_data(ps))
meta <- subset(meta, project == "LOT3")
meta$category <- paste(meta$age, meta$organ, meta$position, sep = " ")

taxa_list <- split(rownames(meta), meta$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

p1=upset(fromList(taxa_list), nsets = length(taxa_list), order.by = "freq", 
      main.bar.color = "forestgreen", sets.bar.color = "steelblue")





# Diagramme de Chord


categories <- names(taxa_list)
shared_mx <- matrix(0, nrow = length(categories), ncol = length(categories), 
                    dimnames = list(categories, categories))

for (i in seq_along(categories)) {
  for (j in seq_along(categories)) {
    if (i != j) {
      shared_mx[i, j] <- length(intersect(taxa_list[[i]], taxa_list[[j]]))
    }
  }
}

chordDiagram(shared_mx, transparency = 0.3, annotationTrack = c("name", "grid"))




ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

ps_phylum <- tax_glom(ps_rel, taxrank = "gbr268_Phylum", NArm = FALSE)

df_phylum <- psmelt(ps_phylum)

df_chord <- df_phylum %>%
  mutate(organ_position = paste(organ, position, sep = "_")) %>%
  group_by(organ_position, gbr268_Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  filter(Abundance > 0.02) %>%
  na.omit()

circos.clear()
chordDiagram(df_chord, transparency = 0.3, annotationTrack = c("name", "grid"))


ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

ps_class<- tax_glom(ps_rel, taxrank = "gbr268_Class", NArm = FALSE)

df_class <- psmelt(ps_class)

df_chord <- df_class %>%
  mutate(organ_position = paste(organ, position, sep = "_")) %>%
  group_by(organ_position, gbr268_Class) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  filter(Abundance > 0.02) %>%
  na.omit()

circos.clear()
chordDiagram(df_chord, transparency = 0.3, annotationTrack = c("name", "grid"))


install.packages("pheatmap")
library(pheatmap)
library(tidyr)

ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))


ps_phylum <- tax_glom(ps_rel, taxrank = "gbr268_Phylum", NArm = FALSE)

df_phylum <- psmelt(ps_phylum)



install.packages("ComplexHeatmap")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)


ps_filtered <- filter_taxa(ps, function(x) sum(x > 0) >= (0.05 * length(x)), TRUE)

ps_rel <- transform_sample_counts(ps_filtered, function(x) x / sum(x))


ps_class <- tax_glom(ps_rel, taxrank = "gbr268_Class", NArm = FALSE)

df_class <- psmelt(ps_class)

df_agg <- df_class %>%
  group_by(gbr268_Class, project, organ, position) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  mutate(group = paste(project, organ, position, sep = "_"))

df_heatmap <- df_agg %>%
  dplyr::select(gbr268_Class, group, Abundance) %>%
  pivot_wider(names_from = group, values_from = Abundance, values_fill = 0) %>%
  filter(!is.na(gbr268_Class)) %>%
  as.data.frame()

rownames(df_heatmap) <- df_heatmap$gbr268_Class
df_heatmap$gbr268_Class <- NULL
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
  dplyr::select(gbr268_Class, gbr268_Phylum) %>%
  distinct() %>%
  filter(!is.na(gbr268_Class)) %>%
  arrange(match(gbr268_Class, rownames(mat)))

row_ha <- rowAnnotation(
  Phylum = row_meta$gbr268_Phylum,
  annotation_name_side = "top"
)
color_mapping <- colorRamp2(c(0, max(mat)/2, max(mat)), c("white", "yellow", "red"))
ann_colors = list(
  Project = c("LOT1" = "lightblue", "LOT2" = "#D34949", "LOT3" = "grey"), # Ajustez si besoin
  Organ = c("root" = "brown", "old_leaf" = "darkgreen", "young_leaf" = "lightgreen"),
  Position = c("endophyte" = "orange", "epiphyte" = "pink") # Remplacez par vos vraies positions
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
  rect_gp = gpar(col = "cadetblue", lwd = 1),
  row_names_side = "right",
  column_names_side = "bottom",
  column_names_rot = -90,
  row_title_rot = 0
)

draw(ht, merge_legend = TRUE)





















library(phyloseq)
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ape)

ps_fam <- tax_glom(ps, taxrank = "gbr268_Order")

tax_df <- as.data.frame(tax_table(ps_fam))
tax_df[is.na(tax_df)] <- "Unknown"
tax_cols <- c("gbr268_Phylum", "gbr268_Class", "gbr268_Order")
tax_df[tax_cols] <- lapply(tax_df[tax_cols], factor)

tree_fam <- as.phylo(~gbr268_Phylum/gbr268_Class/gbr268_Order, data = tax_df)
tree_fam$tip.label <- taxa_names(ps_fam)
tree_fam$edge.length <- rep(1, nrow(tree_fam$edge))

ps_final <- merge_phyloseq(ps_fam, tree_fam)

ps_avg <- merge_samples(ps_final, "organ")
otu_tab <- as.data.frame(t(otu_table(ps_avg)))
df_heatmap_clean <- df_heatmap %>%
  rename(Organ_Type = organ, LogAbundance = Abundance)

df_total_clean <- df_total %>% 
  rename(TotalAbundance = Total)




tax_data <- as.data.frame(tax_table(ps_final)) %>%
  rownames_to_column("label") %>%
  select(label, gbr268_Phylum, gbr268_Class, gbr268_Order)

library(tidytree)
library(treeio)

obj_tree <- phy_tree(ps_final)
obj_tax  <- as.data.frame(tax_table(ps_final)) %>% 
            rownames_to_column("label") # Indispensable pour la jointure

tree_data_final <- full_join(as_tibble(obj_tree), obj_tax, by = "label") %>% 
                   as.treedata()

p <- ggtree(tree_data_final, layout = "circular", branch.length = "none")

p <- p + aes(color = gbr268_Phylum)

p <- p + 
  geom_fruit(
    data = df_heatmap_clean,
    geom = geom_tile,
    mapping = aes(y = OTU, x = Organ_Type, fill = LogAbundance),
    pwidth = 0.15, 
    offset = 0.05
  ) +
  scale_fill_gradientn(colors = c("white", "#66c2a5", "#fc8d62", "#8da0cb"), name = "Log10 Abund") +
  
  geom_fruit(
    data = df_total_clean,
    geom = geom_bar,
    mapping = aes(y = OTU, x = TotalAbundance),
    stat = "identity",
    orientation = "y",
    fill = "grey70",
    pwidth = 0.1,
    offset = 0.1
  ) +
  theme(legend.position = "right")

print(p)





library(phyloseq)
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ape)
library(tidytree)
library(treeio)
library(ggnewscale)

ps_fam <- tax_glom(ps, taxrank = "gbr268_Order")
tax_df <- as.data.frame(tax_table(ps_fam))
tax_df[is.na(tax_df)] <- "Unknown"
tax_cols <- c("gbr268_Phylum", "gbr268_Class", "gbr268_Order")
tax_df[tax_cols] <- lapply(tax_df[tax_cols], factor)

tree_fam <- as.phylo(~gbr268_Phylum/gbr268_Class/gbr268_Order, data = tax_df)
tree_fam$tip.label <- taxa_names(ps_fam)
tree_fam$edge.length <- rep(1, nrow(tree_fam$edge))

ps_final <- merge_phyloseq(ps_fam, tree_fam)
ps_avg <- merge_samples(ps_final, "organ")
otu_tab <- as.data.frame(t(otu_table(ps_avg)))

df_heatmap_clean <- df_heatmap %>%
  rename(Organ_Type = organ, LogAbundance = Abundance)

df_total_clean <- df_total %>% 
  rename(TotalAbundance = Total)

obj_tree <- phy_tree(ps_final)
obj_tax  <- as.data.frame(tax_table(ps_final)) %>% rownames_to_column("label")
tree_data_final <- full_join(as_tibble(obj_tree), obj_tax, by = "label") %>% as.treedata()

p <- ggtree(tree_data_final, layout = "circular", branch.length = "none", size = 0.2) +
  aes(color = gbr268_Phylum) +
  geom_nodepoint(aes(fill = gbr268_Phylum), shape = 21, size = 1.2, color = "black", stroke = 0.1) +
  geom_tippoint(aes(fill = gbr268_Phylum), shape = 21, size = 1.8, color = "black", stroke = 0.1) +
  scale_color_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  new_scale_fill()
  p <- p + 
    geom_tree(size = 0.4, colour = "grey30")
p <- p + 
  geom_fruit(
    data = df_heatmap_clean,
    geom = geom_tile,
    mapping = aes(y = OTU, x = Organ_Type, fill = LogAbundance),
    pwidth = 0.15, 
    offset = 0.08,
    color = "white",
    size = 0.05
  ) +
  scale_fill_gradientn(colors = c("lightblue", "yellow", "red"), name = "Abundance") +
  new_scale_fill()

p <- p + 
  geom_fruit(
    data = df_total_clean,
    geom = geom_bar,
    mapping = aes(y = OTU, x = TotalAbundance, fill = gbr268_Phylum),
    stat = "identity",
    orientation = "y",
    pwidth = 0.25,
    offset = 0.1,
    axis.params = list(axis = "x", text.size = 1.5, nbreak = 3)
  ) +
  scale_fill_brewer(palette = "Pastel1") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7)
  ) +
  layout_circular() + 
  xlim(-5, NA)

print(p)
