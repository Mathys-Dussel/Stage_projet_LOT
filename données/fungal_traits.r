library(fungaltraits)
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(FactoMineR)
library(factoextra)

setwd("~/Documents/Etudes/Stage_projet_LOT/CRBE/données")

ps <- readRDS("donnees_nettoyees.rds")
sample_names(ps) <- gsub("\\.", "-", sub("^X", "", sample_names(ps)))

sam_df <- as.data.frame(as(sample_data(ps), "data.frame"))
sam_df$Sample <- rownames(sam_df)

otu_df <- as.data.frame(as(otu_table(ps), "matrix"))
if (!taxa_are_rows(ps)) otu_df <- as.data.frame(t(otu_df))
otu_df$OTU <- rownames(otu_df)

tax_df <- as.data.frame(as(tax_table(ps), "matrix"))
tax_df$OTU <- rownames(tax_df)

species_col <- grep("(?i)species", colnames(tax_df), value = TRUE)[1]
if (is.na(species_col)) species_col <- colnames(tax_df)[ncol(tax_df)]

tax_df$espece_clean <- tolower(gsub("_", " ", tax_df[[species_col]]))
especes_uniques <- unique(na.omit(tax_df$espece_clean))

ft_db <- fungal_traits()
col_species_ft <- grep("(?i)species", colnames(ft_db), value = TRUE)[1]
ft_db$espece_clean <- tolower(gsub("_", " ", ft_db[[col_species_ft]]))

df_especes <- data.frame(espece_clean = especes_uniques, stringsAsFactors = FALSE)
especes_traits <- left_join(df_especes, ft_db, by = "espece_clean") %>% 
    filter(!is.na(.data[[col_species_ft]]))

taux_remplissage <- colSums(!is.na(especes_traits)) / nrow(especes_traits)
especes_identifiees_clean <- especes_traits[, taux_remplissage >= 0.5]

print(paste("Nombre d'espèces identifiées :", nrow(especes_identifiees_clean)))

otu_long <- otu_df %>%
    pivot_longer(col = -OTU, names_to = "Sample", values_to = "Abundance") %>%
    filter(Abundance > 0)

donnees_completes <- otu_long %>%
    left_join(sam_df, by = "Sample") %>%
    left_join(tax_df[, c("OTU", "espece_clean")], by = "OTU") %>%
    inner_join(especes_identifiees_clean, by = "espece_clean")

col_lifestyle <- grep("(?i)lifestyle|guild", colnames(donnees_completes), value = TRUE)[1]
col_organe    <- grep("(?i)organ|tissue|compartment", colnames(sam_df), value = TRUE)[1]
col_projet    <- grep("(?i)project|projet", colnames(sam_df), value = TRUE)[1]
col_famille   <- grep("(?i)family|famille", colnames(sam_df), value = TRUE)[1]

plot_lifestyle <- function(df, x_col, fill_col, title_group) {
    if (is.na(x_col) || is.na(fill_col)) return(NULL)
    
    df %>%
        filter(!is.na(.data[[x_col]]), !is.na(.data[[fill_col]])) %>%
        ggplot(aes(x = .data[[x_col]], fill = .data[[fill_col]])) +
        geom_bar(position = "fill") +
        theme_minimal() +
        labs(y = "Proportion", title = paste("Mode de vie par", title_group)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

run_chi2_test <- function(df, var1, var2, label) {
    if (is.na(var1) || is.na(var2)) return(NULL)
    
    df_clean <- df %>% filter(!is.na(.data[[var1]]), !is.na(.data[[var2]]))
    tbl <- table(df_clean[[var1]], df_clean[[var2]])
    
    cat("\n--- Tableau de contingence :", label, "---\n")
    print(tbl)
    
    test <- suppressWarnings(chisq.test(tbl, simulate.p.value = TRUE))
    cat("\nRésultat Chi-2 :\n")
    print(test)
}

if (!is.na(col_lifestyle)) {
    
    print(plot_lifestyle(donnees_completes, col_organe, col_lifestyle, "Organe"))
    run_chi2_test(donnees_completes, col_organe, col_lifestyle, "Organe vs Mode de vie")
    
    print(plot_lifestyle(donnees_completes, col_projet, col_lifestyle, "Projet"))
    
    print(plot_lifestyle(donnees_completes, col_famille, col_lifestyle, "Famille de plante"))
    run_chi2_test(donnees_completes, col_famille, col_lifestyle, "Famille vs Mode de vie")
    
    if (!is.na(col_famille)) {
        df_fam <- filter(donnees_completes, !is.na(.data[[col_famille]]), !is.na(.data[[col_lifestyle]]))
        familles <- unique(df_fam[[col_famille]])
        
        calc_pairwise_pval <- function(f1, f2) {
            if (f1 == f2) return(1)
            sub_df <- df_fam[df_fam[[col_famille]] %in% c(f1, f2), ]
            tbl <- table(sub_df[[col_famille]], sub_df[[col_lifestyle]])
            if (nrow(tbl) > 1 && ncol(tbl) > 1) {
                suppressWarnings(chisq.test(tbl, simulate.p.value = TRUE)$p.value)
            } else NA
        }

        pval_df <- expand.grid(Famille1 = familles, Famille2 = familles, stringsAsFactors = FALSE)
        pval_df$p_value <- mapply(calc_pairwise_pval, pval_df$Famille1, pval_df$Famille2)
        
        pval_df$signif <- ifelse(!is.na(pval_df$p_value) & pval_df$p_value < 0.05, "*", "")
        
        p_heatmap <- ggplot(pval_df, aes(x = Famille1, y = Famille2, fill = p_value)) +
            geom_tile(color = "white") +
            geom_text(aes(label = signif), color = "black", size = 5) +
            scale_fill_gradientn(
                colors = c("red", "orange", "white", "lightblue", "blue"), 
                values = scales::rescale(c(0, 0.05, 0.1, 0.5, 1)),
                limits = c(0, 1), name = "P-value", na.value = "grey80"
            ) +
            theme_minimal() +
            labs(title = "Différences de composition de modes de vie (Heatmap P-values)",
                     x = "Famille 1", y = "Famille 2") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        print(p_heatmap)
    }
}