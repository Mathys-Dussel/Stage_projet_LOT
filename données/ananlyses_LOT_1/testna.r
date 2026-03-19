library(phyloseq)
library(ggplot2)
library(ggpubr)
library(patchwork)

ps <- readRDS("donnees_nettoyees.rds")
sample_names(ps) <- gsub("\\.", "-", sub("^X", "", sample_names(ps)))

library(vegan)
dist_bray <- distance(ps, method = "bray")


library(vegan)
ps_hellinger <- transform_sample_counts(ps, function(x) sqrt(x / sum(x)))

#nmds_bray_K2 <- ordinate(ps_hellinger, method = "NMDS", distance = "bray")

nmds_bray_K3 <- ordinate(ps_hellinger, method = "NMDS", distance = "bray", k = 3)




plot_ordination(ps_hellinger, nmds_bray_K3, color = "project", shape = "organ") +
  stat_ellipse(aes(group = project, fill = project), alpha = 0.1, geom = "polygon") +
  theme_bw() +
  labs(title = "NMDS - Bray-Curtis (K=3)")



  ps_lot1 <- subset_samples(ps_hellinger, project == "LOT1")
  nmds_lot1_K3 <- ordinate(ps_lot1, method = "NMDS", distance = "bray", k = 3)

  plot_ordination(ps_lot1, nmds_lot1_K3, color = "organ", shape = "position") +
    stat_ellipse(aes(group = organ, fill = organ), alpha = 0.1, geom = "polygon") +
    theme_bw() +
    labs(title = "NMDS - Bray-Curtis (K=3) - Uniquement LOT1")


