setwd("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses")

ps <- readRDS("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses/donnees/ps_final.rds")




library(fungaltraits)
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

sam_df <- data.frame(sample_data(ps), check.names = FALSE)
sam_df$Sample <- rownames(sam_df)

otu_df <- data.frame(otu_table(ps), check.names = FALSE)
if (!taxa_are_rows(ps)) otu_df <- data.frame(t(otu_table(ps)), check.names = FALSE)
otu_df$OTU <- rownames(otu_df)

tax_df <- data.frame(tax_table(ps), check.names = FALSE)
tax_df$OTU <- rownames(tax_df)

species_col <- grep("(?i)species", colnames(tax_df), value = TRUE)[1]
if (is.na(species_col)) species_col <- colnames(tax_df)[ncol(tax_df)]

tax_df$espece_clean <- tolower(gsub("_", " ", tax_df[[species_col]]))
especes_uniques <- unique(na.omit(tax_df$espece_clean))

ft_db <- fungal_traits()
col_species_ft <- grep("(?i)species", colnames(ft_db), value = TRUE)[1]
ft_db$espece_clean <- tolower(gsub("_", " ", ft_db[[col_species_ft]]))

df_especes <- data.frame(espece_clean = especes_uniques)
especes_traits <- df_especes %>%
  left_join(ft_db, by = "espece_clean") %>%
  filter(!is.na(.data[[col_species_ft]]))

taux_remplissage <- colSums(!is.na(especes_traits)) / nrow(especes_traits)
especes_identifiees_clean <- especes_traits[, taux_remplissage >= 0.5]

print(paste("Nombre d'espèces avec informations dans FungalTraits :", nrow(especes_identifiees_clean)))

otu_long <- otu_df %>%
  pivot_longer(col = -OTU, names_to = "Sample", values_to = "Abundance") %>%
  filter(Abundance > 0)

donnees_completes <- otu_long %>%
  left_join(sam_df, by = "Sample") %>%
  left_join(tax_df[, c("OTU", "espece_clean")], by = "OTU") %>%
  inner_join(especes_identifiees_clean, by = "espece_clean")

col_lifestyle <- grep("(?i)lifestyle|guild", colnames(donnees_completes), value = TRUE)[1]
col_organe    <- grep("(?i)organ|tissue|compartment", colnames(sam_df), value = TRUE)[1]
col_projet    <- grep("(?i)^(project|projet)$", colnames(sam_df), value = TRUE)[1]
col_famille   <- grep("(?i)family|famille", colnames(sam_df), value = TRUE)[1]

plot_lifestyle <- function(df, x_var, fill_var, title_suffix) {
  if (is.na(x_var) || is.na(fill_var)) return(NULL)
  
  df %>%
    filter(!is.na(.data[[x_var]]), !is.na(.data[[fill_var]])) %>%
    ggplot(aes(x = .data[[x_var]], fill = .data[[fill_var]])) +
    geom_bar(position = "fill") +
    theme_minimal() +
    labs(y = "Proportion", title = paste("Mode de vie par", title_suffix), fill = "Mode de vie") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

run_chi2_test <- function(df, var1, var2, label) {
  if (is.na(var1) || is.na(var2)) return(NULL)
  
  df_clean <- df %>% filter(!is.na(.data[[var1]]), !is.na(.data[[var2]]))
  tbl <- table(df_clean[[var1]], df_clean[[var2]])
  
  cat("\n---", label, "---\n")
  print(suppressWarnings(chisq.test(tbl, simulate.p.value = TRUE)))
}

if (!is.na(col_lifestyle)) {
  
  print(plot_lifestyle(donnees_completes, col_organe, col_lifestyle, "Organe"))
  run_chi2_test(donnees_completes, col_organe, col_lifestyle, "Test : Organe vs Mode de vie")
  
  print(plot_lifestyle(donnees_completes, col_projet, col_lifestyle, "Projet"))
  
  print(plot_lifestyle(donnees_completes, col_famille, col_lifestyle, "Famille de plante"))
  run_chi2_test(donnees_completes, col_famille, col_lifestyle, "Test : Famille vs Mode de vie")
  
  if (!is.na(col_famille)) {
    df_fam <- donnees_completes %>% filter(!is.na(.data[[col_famille]]), !is.na(.data[[col_lifestyle]]))
    familles <- unique(df_fam[[col_famille]])
    
    calc_pval <- function(f1, f2) {
      if (f1 == f2) return(1)
      sub_df <- df_fam[df_fam[[col_famille]] %in% c(f1, f2), ]
      tbl <- table(sub_df[[col_famille]], sub_df[[col_lifestyle]])
      if (nrow(tbl) > 1 && ncol(tbl) > 1) {
        suppressWarnings(chisq.test(tbl, simulate.p.value = TRUE)$p.value)
      } else {
        NA
      }
    }
    
    pval_df <- expand.grid(Famille1 = familles, Famille2 = familles)
    pval_df$p_value <- mapply(calc_pval, pval_df$Famille1, pval_df$Famille2)
    pval_df$signif <- ifelse(!is.na(pval_df$p_value) & pval_df$p_value < 0.05, "*", "")
    
    p_heatmap <- ggplot(pval_df, aes(x = Famille1, y = Famille2, fill = p_value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = signif), size = 6) +
      scale_fill_gradient2(low = "green", mid = "white", high = "red", midpoint = 0.1, 
                           limits = c(0, 1), name = "P-value", na.value = "grey80") +
      theme_minimal() +
      labs(title = "Significativité des différences de mode de vie entre familles",
           x = "Famille", y = "Famille") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    print(p_heatmap)
  }
}














df_fonctions <- donnees_completes %>%
  group_by(Sample, !!sym(col_lifestyle)) %>%
  summarise(Abundance = sum(Abundance), .groups = 'drop') %>%
  pivot_wider(names_from = !!sym(col_lifestyle), values_from = Abundance, values_fill = 0)

df_fonctions_mat <- as.data.frame(df_fonctions)
rownames(df_fonctions_mat) <- df_fonctions_mat$Sample
df_fonctions_mat$Sample <- NULL

metadata_sub <- sam_df[rownames(df_fonctions_mat), ]


library(vegan)

cca_mod <- cca(df_fonctions_mat ~ organ + project + plant_family, data = metadata_sub)


anova(cca_mod) 

plot(cca_mod, display = c("sp", "bp"), main = "CCA : Influence des facteurs sur les fonctions fongiques")


summary(cca_mod)$cont

anova(cca_mod, by = "term")


scores_cca <- scores(cca_mod, display = "bp")
print(scores_cca)
