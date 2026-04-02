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

setwd("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses")
ps <- readRDS("donnees/ps_final.rds")

# Abondances relatives par classe

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
         x = "Échantillons ", y = "Abondance relative", fill = "Classe fongique") +
    theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggplot(df_class, aes(x = plant_family, y = Abundance, fill = gbr268_Class)) +
    geom_bar(stat = "identity", position = "fill", color = NA) +
    facet_wrap(~ organ_position, scales = "free_x", ncol = 3) +
    theme_bw() +
    labs(title = "Abondances relatives des principales classes fongiques (> 2%) par famille de plante",
       x = "Famille de plante", y = "Abondance relative", fill = "Classe fongique") +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(df_class, aes(x = project, y = Abundance, fill = gbr268_Class)) +
    geom_bar(stat = "identity", position = "fill", color = NA) +
    facet_wrap(~ organ_position, scales = "free_x", ncol = 3) +
    theme_bw() +
    labs(title = "Abondances relatives des principales classes fongiques (> 2%) par zone",
         x = "Projet", y = "Abondance relative", fill = "Classe fongique") +
    theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank())





# Venn Diagramme
taxons_list <- list(
  Root = core_members(subset_samples(ps, organ == "root"), detection = 0, prevalence = 0.001),
  Young_Leaf = core_members(subset_samples(ps, organ == "young_leaf"), detection = 0, prevalence = 0.001),
  Old_Leaf = core_members(subset_samples(ps, organ == "old_leaf"), detection = 0, prevalence = 0.001)
)

ggVennDiagram(taxons_list, label_alpha = 0) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  labs(title = "Partage des OTUs entre les organes") +
  theme(legend.position = "none")

# Heatmap du core microbiota au niveau famille

ps_family <- tax_glom(ps, taxrank = "gbr268_Family", NArm = TRUE)
ps_comp <- microbiome::transform(ps_family, "compositional")
tax_df_core <- as.data.frame(tax_table(ps_comp))
taxa_names(ps_comp) <- paste0(tax_df_core$gbr268_Phylum, ": ", tax_df_core$gbr268_Family)
ps_core <- core(ps_comp, detection = 0.001, prevalence = 0.3)

plot_core(ps_core, plot.type = "heatmap", 
          prevalences = seq(0.1, 1, 0.1), 
          detections = c(0.001, 0.005, 0.01, 0.05, 0.1),
          colours = rev(brewer.pal(5, "RdYlBu"))) +
  labs(x = "Seuil de Détection (%)", y = "Taxonomie (Phylum: Famille)", title = "Core Microbiota au niveau Famille") +
  scale_x_discrete(labels = c("0.1%", "0.5%", "1%", "5%", "10%")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9, face = "bold"))













  library(indicspecies)
  library(tibble)
  library(tidyr)

  otu_mat <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) { otu_mat <- t(otu_mat) }

  metadata_df <- as.data.frame(sample_data(ps))

  indval <- multipatt(otu_mat, metadata_df$plant_family, func = "IndVal.g", duleg = TRUE, control = how(nperm = 99))

  signif_otus <- indval$sign %>%
    rownames_to_column("OTU") %>%
    filter(p.value <= 0.05) %>%
    pivot_longer(cols = starts_with("s."), names_to = "plant_family", values_to = "is_indicator") %>%
    filter(is_indicator == 1) %>%
    mutate(plant_family = gsub("^s\\.", "", plant_family))

  stat_df <- as.data.frame(indval$str) %>%
    rownames_to_column("OTU") %>%
    pivot_longer(
      cols = -OTU, 
      names_to = "plant_family", 
      values_to = "Stat"
    ) %>%
    mutate(plant_family = gsub("^s\\.", "", plant_family))

  tax_info <- as.data.frame(tax_table(ps)) %>% 
    rownames_to_column("OTU") %>% 
    select(OTU, gbr268_Class)

  plot_df <- stat_df %>%
    inner_join(signif_otus, by = c("OTU", "plant_family")) %>%
    left_join(tax_info, by = "OTU") %>%
    mutate(OTU_label = paste0(gbr268_Class, " - ", OTU)) %>%
    select(OTU, OTU_label, plant_family, Stat, gbr268_Class) %>%
    distinct()

  head(stat_df)
  plot_df_agg <- plot_df %>%
    group_by(plant_family, gbr268_Class) %>%
    summarise(Num_OTUs = n(), Mean_Stat = mean(Stat), .groups = "drop")

  ggplot(plot_df_agg, aes(x = plant_family, y = gbr268_Class, size = Num_OTUs, color = Mean_Stat)) +
    geom_point() +
    scale_color_distiller(palette = "YlOrRd", direction = 1) +
    theme_bw() +
    labs(title = "Classes indicatrices (Indicspecies) par Famille de Plante",
         x = "Famille de Plante", 
         y = "Classe Fongique",
         size = "Nombre d'OTUs", 
         color = "IndVal Stat Moyen") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.text.y = element_text(size = 10))




          library(phyloseq)
          library(indicspecies)
          library(ggtree)
          library(tidyverse)
          library(ape)

          ps_filt <- filter_taxa(ps, function(x) sum(x > 0) > (0.05 * nsamples(ps)), TRUE)

          ps_filt <- prune_taxa(names(sort(taxa_sums(ps_filt), TRUE)[1:min(2000, ntaxa(ps_filt))]), ps_filt)

          if(taxa_are_rows(ps_filt)) {
            otu_tab <- as.matrix(t(otu_table(ps_filt)))
          } else {
            otu_tab <- as.matrix(otu_table(ps_filt))
          }

          groups_vec <- as.character(sample_data(ps_filt)$plant_family)

          inv <- multipatt(otu_tab, groups_vec, 
                           func = "IndVal.g", 
                           control = how(nperm = 99), duleg = TRUE) 

          inv_res <- inv$sign %>%
            rownames_to_column("OTU") %>%
            filter(p.value < 0.05) %>%
            mutate(Group_Name = colnames(inv$comb)[index])

          tax_df <- as.data.frame(tax_table(ps_filt)) %>%
            rownames_to_column("OTU") %>%
            filter(OTU %in% inv_res$OTU) %>%
            filter(!is.na(gbr268_Phylum) & !is.na(gbr268_Class) & !is.na(gbr268_Order)) %>%
            filter(gbr268_Phylum != "Unknown" & gbr268_Class != "Unknown" & gbr268_Order != "Unknown") %>%
            mutate(across(c(gbr268_Phylum, gbr268_Class, gbr268_Order), as.factor))

          heatmap_data <- inv_res %>%
            left_join(tax_df, by = "OTU") %>%
            filter(!is.na(gbr268_Order)) %>%
            group_by(gbr268_Order, Group_Name) %>%
            summarise(stat = mean(stat, na.rm = TRUE), .groups = 'drop') %>%
            pivot_wider(names_from = Group_Name, values_from = stat, values_fill = 0) %>%
            column_to_rownames("gbr268_Order")

          tree_structure <- tax_df %>% 
            select(gbr268_Phylum, gbr268_Class, gbr268_Order) %>% 
            distinct()
          tree <- as.phylo(~gbr268_Phylum/gbr268_Class/gbr268_Order, data = tree_structure)

          p <- ggtree(tree, layout = "rectangular") + 
            geom_tiplab(size = 3, align = TRUE, offset = 0.2) 

          p_final <- gheatmap(p, heatmap_data, 
                             offset = 5.0, width = 1.5, color = "black",
                               colnames_angle = 45, colnames_offset_y = -3,
                             font.size = 3) +
            scale_fill_gradient(low = "white", high = "blue", name = "IndVal Moyenne") +
            theme(legend.position = "right", plot.margin = margin(10, 120, 50, 10))

          print(p_final)




            ps_filt <- filter_taxa(ps, function(x) sum(x > 0) > (0.001 * nsamples(ps)), TRUE)

            ps_filt <- prune_taxa(names(sort(taxa_sums(ps_filt), TRUE)[1:min(2000, ntaxa(ps_filt))]), ps_filt)

            if(taxa_are_rows(ps_filt)) {
            otu_tab <- as.matrix(t(otu_table(ps_filt)))
            } else {
            otu_tab <- as.matrix(otu_table(ps_filt))
            }

            groups_vec <- paste(sample_data(ps_filt)$project, 
                      sample_data(ps_filt)$organ, 
                      sample_data(ps_filt)$position, sep = "_")

            inv <- multipatt(otu_tab, groups_vec, 
                     func = "IndVal.g", 
                     control = how(nperm = 99), duleg = TRUE) 

            inv_res <- inv$sign %>%
            rownames_to_column("OTU") %>%
            filter(p.value < 0.05) %>%
            mutate(Group_Name = colnames(inv$comb)[index])

            tax_df <- as.data.frame(tax_table(ps_filt)) %>%
            rownames_to_column("OTU") %>%
            filter(OTU %in% inv_res$OTU) %>%
            filter(!is.na(gbr268_Phylum) & !is.na(gbr268_Class) & !is.na(gbr268_Order)) %>%
            filter(gbr268_Phylum != "Unknown" & gbr268_Class != "Unknown" & gbr268_Order != "Unknown") %>%
            mutate(across(c(gbr268_Phylum, gbr268_Class, gbr268_Order), as.factor))

            heatmap_data <- inv_res %>%
            left_join(tax_df, by = "OTU") %>%
            filter(!is.na(gbr268_Order)) %>%
            group_by(gbr268_Order, Group_Name) %>%
            summarise(stat = mean(stat, na.rm = TRUE), .groups = 'drop') %>%
            pivot_wider(names_from = Group_Name, values_from = stat, values_fill = 0) %>%
            column_to_rownames("gbr268_Order")

            tree_structure <- tax_df %>% 
            select(gbr268_Phylum, gbr268_Class, gbr268_Order) %>% 
            distinct()
            tree <- as.phylo(~gbr268_Phylum/gbr268_Class/gbr268_Order, data = tree_structure)

            p <- ggtree(tree, layout = "rectangular") + 
            geom_tiplab(size = 3, align = TRUE, offset = 1) 

            heatmap_data <- heatmap_data[, sort(colnames(heatmap_data))]

            p_final <- gheatmap(p, heatmap_data, 
                 offset = 5.0, width = 1.5, color = "black",
                   colnames_angle = 90, colnames_offset_y = -5,
                 font.size = 3) +
            scale_fill_gradient(low = "white", high = "blue", name = "IndVal Moyenne") +
            scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
            coord_cartesian(clip = "off") +
            theme(legend.position = "right", plot.margin = margin(10, 120, 10, 10))

            print(p_final)



            p <- ggtree(tree, layout = "dendrogramme") + 
            geom_tiplab(size = 3, align = TRUE, offset = 1) 

print(p)

summary(inv)
summary(inv_res)






            ps_filt <- filter_taxa(ps, function(x) sum(x > 0) > (0.001 * nsamples(ps)), TRUE)

            ps_filt <- prune_taxa(names(sort(taxa_sums(ps_filt), TRUE)[1:min(2000, ntaxa(ps_filt))]), ps_filt)

            if(taxa_are_rows(ps_filt)) {
            otu_tab <- as.matrix(t(otu_table(ps_filt)))
            } else {
            otu_tab <- as.matrix(otu_table(ps_filt))
            }

            groups_vec <- paste(sample_data(ps_filt)$project, 
                      sample_data(ps_filt)$organ, 
                      sample_data(ps_filt)$position, sep = "_")

            inv <- multipatt(otu_tab, groups_vec, 
                     func = "IndVal.g", 
                     control = how(nperm = 99), duleg = TRUE) 

            inv_res <- inv$sign %>%
            rownames_to_column("OTU") %>%
            filter(p.value < 0.05) %>%
            mutate(Group_Name = colnames(inv$comb)[index])

            tax_df <- as.data.frame(tax_table(ps_filt)) %>%
            rownames_to_column("OTU") %>%
            filter(OTU %in% inv_res$OTU) %>%
            filter(!is.na(gbr268_Phylum) & !is.na(gbr268_Class) & !is.na(gbr268_Order) & !is.na(gbr268_Family)) %>%
            filter(gbr268_Phylum != "Unknown" & gbr268_Class != "Unknown" & gbr268_Order != "Unknown"& gbr268_Family != "Unknown") %>%
            mutate(across(c(gbr268_Phylum, gbr268_Class, gbr268_Order, gbr268_Family), as.factor))

            heatmap_data <- inv_res %>%
            left_join(tax_df, by = "OTU") %>%
            filter(!is.na(gbr268_Family)) %>%
            group_by(gbr268_Family, Group_Name) %>%
            summarise(stat = mean(stat, na.rm = TRUE), .groups = 'drop') %>%
            pivot_wider(names_from = Group_Name, values_from = stat, values_fill = 0) %>%
            column_to_rownames("gbr268_Family")

            tree_structure <- tax_df %>% 
            select(gbr268_Phylum, gbr268_Class, gbr268_Order, gbr268_Family) %>% 
            distinct()
            tree <- as.phylo(~gbr268_Phylum/gbr268_Class/gbr268_Order/gbr268_Family, data = tree_structure)

            p <- ggtree(tree, layout = "rectangular") + 
            geom_tiplab(size = 3, align = TRUE, offset = 1) 

            heatmap_data <- heatmap_data[, sort(colnames(heatmap_data))]

            p_final <- gheatmap(p, heatmap_data, 
                 offset = 5.0, width = 1.5, color = "black",
                   colnames_angle = 90, colnames_offset_y = -5,
                 font.size = 3) +
            scale_fill_gradient(low = "white", high = "blue", name = "IndVal Moyenne") +
            scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
            coord_cartesian(clip = "off") +
            theme(legend.position = "right", plot.margin = margin(10, 120, 10, 10))

            print(p_final)

















# Arbre taxonomique des 100 OTUs les plus abondants

top_otus <- names(sort(taxa_sums(ps), decreasing = TRUE)[1:100])
ps_top <- prune_taxa(top_otus, ps)
tax_df_top <- as.data.frame(tax_table(ps_top))

smart_names <- apply(tax_df_top, 1, function(x) {
  name <- paste(x["gbr268_Genus"], x["gbr268_Species"])
  if (grepl("Unknown", name)) name <- paste0(x["gbr268_Family"], " (unkn. sp.)")
  return(name)
})

taxa_names(ps_top) <- make.unique(smart_names)
tax_df_fact <- as.data.frame(unclass(tax_df_top), stringsAsFactors = TRUE)
dist_taxo <- cluster::daisy(tax_df_fact, metric = "gower")
arbre_taxo <- as.phylo(hclust(dist_taxo, method = "average"))
arbre_taxo$tip.label <- taxa_names(ps_top)
phy_tree(ps_top) <- arbre_taxo

abundance_data <- data.frame(label = taxa_names(ps_top), 
                            Abondance = taxa_sums(ps_top))

p_tree <- ggtree(ps_top, aes(color = gbr268_Phylum), ladderize = TRUE) %<+% abundance_data

p_tree +
  geom_tippoint(aes(size = Abondance), alpha = 0.6) + 
  geom_tiplab(size = 2.8, offset = 0.02, fontface = "italic") + 
  scale_size_continuous(range = c(1, 6), name = "Abondance totale") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.8))) + 
  theme_tree2() +
  labs(title = "Top 100 OTU", 
       color = "Phylum")

# ANCOMB2 pour identifier les taxons différentiellement abondants entre les organes

metadata <- as.data.frame(as(sample_data(ps), "data.frame"))
metadata$organ <- factor(metadata$organ, levels = c("old_leaf", "young_leaf", "root"))
sample_data(ps) <- sample_data(metadata)
ps_raw <- prune_samples(sample_sums(ps) > 0, ps)

out_3groups <- ancombc2(data = ps_raw, fix_formula = "organ", p_adj_method = "holm", 
                        group = "organ", struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05)

tax_df_ancom <- as.data.frame(as(tax_table(ps_raw), "matrix"))
tax_df_ancom$taxon <- rownames(tax_df_ancom)
tax_df_ancom$Genre_clean <- ifelse(is.na(tax_df_ancom$gbr268_Genus), 
                                   paste("Uncl.", tax_df_ancom$gbr268_Family), 
                                   as.character(tax_df_ancom$gbr268_Genus))

res_df <- out_3groups$res %>%
  select(taxon, lfc_young = lfc_organyoung_leaf, p_young = p_organyoung_leaf,
         lfc_root = lfc_organroot, p_root = p_organroot)

res_final <- inner_join(res_df, tax_df_ancom, by = "taxon")

top_root <- res_final %>% filter(p_root < 0.05 & lfc_root > 0) %>% arrange(desc(lfc_root)) %>% head(10) %>% mutate(Group = "Racines (vs Vieilles feuilles)", LFC_val = lfc_root)
top_young <- res_final %>% filter(p_young < 0.05 & lfc_young > 0) %>% arrange(desc(lfc_young)) %>% head(10) %>% mutate(Group = "Jeunes feuilles (vs Vieilles feuilles)", LFC_val = lfc_young)
top_old <- res_final %>% filter(p_root < 0.05 & lfc_root < 0 & p_young < 0.05 & lfc_young < 0) %>% arrange(lfc_root) %>% head(10) %>% mutate(Group = "Vieilles feuilles", LFC_val = abs(lfc_root))

top_all <- rbind(top_root, top_young, top_old)

ggplot(top_all, aes(x = reorder(Genre_clean, LFC_val), y = LFC_val, fill = gbr268_Family)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~Group, scales = "free_y", ncol = 1) +
  labs(x = "Genre", y = "Echelle log des abondances", fill = "Famille") +
  theme_bw()



# Niche fongique par organe : Calcul de l'indice de Levins et moyenne pondérée par échantillon

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

p_project <- ggplot(final_df, aes(x = project, y = Levins_Mean, fill = project)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 1, color = "black") +
  labs(title = "Niche fongique par projet",
       x = "Projet",
       y = "Indice de Levins Moyen") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 11, face = "bold"),
        title = element_text(size = 13, face = "bold"))

p_project_final <- p_project + stat_compare_means(method = "kruskal.test", label.y = max(final_df$Levins_Mean, na.rm = TRUE) + 0.05) +
               stat_compare_means(label = "p.signif", method = "wilcox.test", 
                                  ref.group = ".all.", hide.ns = FALSE)

print(p_project_final)


final_df$position <- metadata$position[match(final_df$Sample, rownames(metadata))]
final_df$position <- as.factor(final_df$position)

p_organ_pos <- ggplot(final_df, aes(x = organ, y = Levins_Mean, fill = position)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.6, position = position_dodge(0.8)) +
  geom_point(aes(color = position), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), alpha = 0.4, size = 1) +
  labs(title = "Niche fongique par organe et position",
       x = "Organe",
       y = "Indice de Levins Moyen",
       fill = "Position",
       color = "Position") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        title = element_text(size = 13, face = "bold"))

p_organ_pos_final <- p_organ_pos + 
  stat_compare_means(aes(group = position), method = "wilcox.test", label = "p.signif")

print(p_organ_pos_final)







p_project_organ <- ggplot(final_df, aes(x = project, y = Levins_Mean, fill = organ)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.6, position = position_dodge(0.8)) +
  geom_point(aes(color = organ), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), alpha = 0.4, size = 1) +
  labs(title = "Niche fongique par projet et organe",
       x = "Projet",
       y = "Indice de Levins Moyen",
       fill = "Organe",
       color = "Organe") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 11, face = "bold"),
        title = element_text(size = 13, face = "bold"))

p_project_organ_final <- p_project_organ + 
  stat_compare_means(aes(group = organ), method = "kruskal.test", label = "p.signif")

print(p_project_organ_final)




final_df$plant_family <- metadata$plant_family[match(final_df$Sample, rownames(metadata))]
final_df$plant_family <- as.factor(final_df$plant_family)

p_plant_family <- ggplot(final_df, aes(x = plant_family, y = Levins_Mean, fill = plant_family)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 1, color = "black") +
  labs(title = "Niche fongique par famille de plante",
       x = "Famille de plante",
       y = "Indice de Levins Moyen") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 11, face = "bold"),
        title = element_text(size = 13, face = "bold"))

p_plant_family_final <- p_plant_family + 
  stat_compare_means(method = "kruskal.test", label.y = max(final_df$Levins_Mean, na.rm = TRUE) + 0.05) +
  stat_compare_means(label = "p.signif", method = "wilcox.test", 
                     ref.group = ".all.", hide.ns = FALSE)

print(p_plant_family_final)



library(patchwork)

combined_plot <- (p_project_final | p_organ_pos_final) / 
         (p_project_organ_final + p_plant_family_final ) +
  plot_annotation(tag_levels = 'A')

print(combined_plot)


