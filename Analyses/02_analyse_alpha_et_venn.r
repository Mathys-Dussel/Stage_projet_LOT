setwd("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses")

ps <- readRDS("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/donnees_nettoyees.rds")
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



library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)

otu_mat <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) { otu_mat <- t(otu_mat) }

out <- rarecurve(otu_mat, step = 500, label = FALSE, tidy = TRUE)

colnames(out) <- c("Site", "Profondeur", "OTUs")

metadata_sub <- data.frame(
  Site = sample_names(ps), 
  Organe = as.character(sample_data(ps)$organ)
)

rare_df <- out %>%
  left_join(metadata_sub, by = "Site")

ggplot(rare_df, aes(x = Profondeur, y = OTUs, group = Site, color = Organe)) +
  geom_line(alpha = 0.15) + 
  scale_color_manual(values = c("root" = "#4DAF4A", "old_leaf" = "#E41A1C", "young_leaf" = "#377EB8")) +
  labs(title = "Courbes de raréfaction (1107 échantillons)",
       x = "Nombre de séquences lues",
       y = "Nombre d'OTUs (Richesse)") +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 2)))






library(vegan)
library(ggplot2)

otu_mat <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) { otu_mat <- t(otu_mat) }

get_specacc <- function(mat, groups, target) {
  sub_mat <- mat[groups == target, ]
  acc <- specaccum(sub_mat, method = "random", permutations = 100)
  
  return(data.frame(
    Samples = acc$sites,
    Richness = acc$richness,
    SD = acc$sd,
    Organe = target
  ))
}

organs <- unique(as.character(sample_data(ps)$organ))
df_acc <- do.call(rbind, lapply(organs, function(x) get_specacc(otu_mat, sample_data(ps)$organ, x)))

ggplot(df_acc, aes(x = Samples, y = Richness, color = Organe, fill = Organe)) +
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Courbe d'accumulation des espèces (SAC)",
       x = "Nombre d'échantillons (Individus)",
       y = "Nombre d'OTUs cumulés") +
  theme_bw()







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
meta$category <- paste(meta$age, meta$organ, meta$position, sep = "_")

taxa_list <- split(rownames(meta), meta$category) |>
  lapply(function(samples) {
    taxa_sums <- taxa_sums(prune_samples(samples, ps))
    names(taxa_sums[taxa_sums > 0])
  })

upset(fromList(taxa_list), nsets = length(taxa_list), order.by = "freq", 
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


