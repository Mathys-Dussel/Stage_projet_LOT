setwd("~/Documents/Etudes/Stage_projet_LOT/CRBE/données/analyses_di")

ps <- readRDS("ps_final.rds")

library(dplyr)
library(phyloseq)
library(ggplot2)

ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

ps_class <- tax_glom(ps_rel, taxrank = "gbr268_Class", NArm = FALSE)

ps_class_filtered <- filter_taxa(ps_class, function(x) mean(x) > 0.02, TRUE)

df_class <- psmelt(ps_class_filtered)

df_class$organ_position <- interaction(df_class$organ, df_class$position, sep = " - ")


ggplot(df_class, aes(x = Sample, y = Abundance, fill = gbr268_Class)) +
    geom_bar(stat = "identity", position = "fill", color = NA) +
    facet_wrap(~ organ_position, scales = "free_x", ncol = 3) +
    theme_bw() +
    labs(title = "Abondances relatives des principales classes fongiques (> 2%) par zone",
         x = "Échantillons ",
         y = "Abondance relative",
         fill = "Classe fongique") +
    theme(legend.position = "bottom",
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())

ggplot(df_class, aes(x = project, y = Abundance, fill = gbr268_Class)) +
    geom_bar(stat = "identity", position = "fill", color = NA) +
    facet_wrap(~ organ_position, scales = "free_x", ncol = 3) +
    theme_bw() +
    labs(title = "Abondances relatives des principales classes fongiques (> 2%) par zone",
         x = "Échantillons ",
         y = "Abondance relative",
         fill = "Classe fongique") +
    theme(legend.position = "bottom",
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())


library(ggVennDiagram)
library(microbiome)
taxons_list <- list(
  Root = core_members(subset_samples(ps, organ == "root"), detection = 0, prevalence = 0.001),
  Young_Leaf = core_members(subset_samples(ps, organ == "young_leaf"), detection = 0, prevalence = 0.001),
  Old_Leaf = core_members(subset_samples(ps, organ == "old_leaf"), detection = 0, prevalence = 0.001)
)

ggVennDiagram(taxons_list, label_alpha = 0) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  labs(title = "Partage des OTUs entre les organes") +
  theme(legend.position = "none")



ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

tax_df <- as.data.frame(tax_table(ps_rel))

new_labels <- tax_df %>%
  mutate(
    Phylum = ifelse(is.na(gbr268_Phylum), "Unclassified", gbr268_Phylum),
    Family = ifelse(is.na(gbr268_Family), "Unclassified", gbr268_Family),
    FullLabel = paste0(Phylum, ": ", Family, " (", rownames(tax_df), ")")
  )

taxa_names(ps_rel) <- new_labels$FullLabel

ps_core <- core(ps_rel, detection = 0.001, prevalence = 0.2)

library(microbiome)
library(RColorBrewer)
library(dplyr)

ps_family <- tax_glom(ps, taxrank = "gbr268_Family", NArm = TRUE)

ps_comp <- microbiome::transform(ps_family, "compositional")

tax_df <- as.data.frame(tax_table(ps_comp))
new_labels <- paste0(tax_df$gbr268_Phylum, ": ", tax_df$gbr268_Family)
taxa_names(ps_comp) <- new_labels

ps_core <- core(ps_comp, detection = 0.001, prevalence = 0.3)

plot_core(ps_core, 
          plot.type = "heatmap", 
          prevalences = seq(0.1, 1, 0.1), 
          detections = c(0.001, 0.005, 0.01, 0.05, 0.1),
          colours = rev(brewer.pal(5, "RdYlBu"))) +
  labs(x = "Seuil de Détection (Abondance Relative %)", 
       y = "Taxonomie (Phylum: Famille)",
       title = "Core Microbiota au niveau Famille") +
  scale_x_discrete(labels = c("0.1%", "0.5%", "1%", "5%", "10%")) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(size = 9)
  )



library(cluster)
library(ape)
library(ggtree)
library(phyloseq)

top_otus <- names(sort(taxa_sums(ps), decreasing = TRUE)[1:100])
ps_top <- prune_taxa(top_otus, ps)

tax_df <- as.data.frame(tax_table(ps_top))

smart_names <- apply(tax_df, 1, function(x) {
  name <- paste(x["gbr268_Genus"], x["gbr268_Species"])
  if (grepl("Unknown", name)) {
    name <- paste0(x["gbr268_Family"], " (unkn. sp.)")
  }
  return(name)
})

taxa_names(ps_top) <- make.unique(smart_names)

tax_df_fact <- as.data.frame(unclass(tax_df), stringsAsFactors = TRUE)
dist_taxo <- daisy(tax_df_fact, metric = "gower")
arbre_taxo <- as.phylo(hclust(dist_taxo, method = "average"))
arbre_taxo$tip.label <- taxa_names(ps_top)
phy_tree(ps_top) <- arbre_taxo

sample_data(ps_top)$dummy <- 1 
df_abundance <- data.frame(Abundance = taxa_sums(ps_top))

ps_top_final <- ps_top 

ggtree(ps_top_final, aes(color = gbr268_Phylum), ladderize = TRUE) +
  geom_tippoint(aes(size = Abundance), alpha = 0.6) + 
  geom_tiplab(size = 2.8, offset = 0.02, fontface = "italic") + 
  scale_size_continuous(range = c(1, 6), name = "Somme des abondances") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.8))) + 
  theme_tree2() +
  theme(legend.position = "right") +
  labs(title = "Top 100 Espèces - Phylogénie Taxonomique",
       subtitle = "La taille des points indique l'abondance totale",
       color = "Phylum")
