install.packages("devtools")
devtools::install_github("ropenscilabs/datastorr")
devtools::install_github("traitecoevo/fungaltraits")
library(fungaltraits)

fungal_traits()
head(fungal_traits())

setwd("~/Documents/Etudes/Stage_projet_LOT/CRBE/données")

library(phyloseq)
library(ggplot2)
library(ggpubr)
library(patchwork)

ps <- readRDS("donnees_nettoyees.rds")
sample_names(ps) <- gsub("\\.", "-", sub("^X", "", sample_names(ps)))

# Extraire les noms d'espèces uniques depuis la table taxonomique
species_col <- grep("(?i)species", colnames(tax_table(ps)), value = TRUE)
if (length(species_col) == 0) species_col <- colnames(tax_table(ps))[ncol(tax_table(ps))]
especes_uniques <- unique(na.omit(as.character(tax_table(ps)[, species_col[1]])))


library(dplyr)

# Charger la base de données FungalTraits
ft_db <- fungal_traits()

# Identifier la colonne contenant le nom des espèces dans FungalTraits (généralement "Species_name")
col_species_ft <- grep("(?i)species", colnames(ft_db), value = TRUE)[1]

# Créer un dataframe à partir des espèces uniques calculées précédemment
df_especes <- data.frame(espece = especes_uniques, stringsAsFactors = FALSE)

# Fusionner (Left join) pour ajouter les traits aux espèces trouvées (NA si non trouvé)
especes_traits <- merge(
    x = df_especes, 
    y = ft_db, 
    by.x = "espece", 
    by.y = col_species_ft, 
    all.x = TRUE
)

# Aperçu du résultat
head(especes_traits)
# Afficher uniquement les espèces dont les traits ont été trouvés dans la base de données
especes_identifiees <- especes_traits[!is.na(especes_traits[[2]]), ]

# Aperçu du résultat filtré
head(especes_identifiees)
# Vérifions le format des noms d'espèces dans les deux tables
head(df_especes$espece)
head(ft_db[[col_species_ft]])

# Standardisation des noms pour faciliter la correspondance (minuscules, remplacement des underscores par des espaces)
df_especes$espece_clean <- tolower(gsub("_", " ", df_especes$espece))

# Assurez-vous d'avoir bien identifié la colonne des espèces dans ft_db (généralement "Genus_species" ou "Name")
# S'il y a un doute, on peut forcer l'utilisation de la bonne colonne si connue, ex: "Name"
ft_db$espece_clean <- tolower(gsub("_", " ", ft_db[[col_species_ft]]))

# Refaire la fusion sur la colonne nettoyée
especes_traits <- merge(
    x = df_especes, 
    y = ft_db, 
    by = "espece_clean", 
    all.x = TRUE
)

# Filtrer et compter le nombre d'espèces identifiées
especes_identifiees <- especes_traits[!is.na(especes_traits[[col_species_ft]]), ]
nrow(especes_identifiees)
head(especes_identifiees)
# Calculer la proportion de valeurs non manquantes pour chaque colonne
taux_remplissage <- colSums(!is.na(especes_identifiees)) / nrow(especes_identifiees)

# Définir un seuil (par exemple, conserver les colonnes remplies au moins à 50%)
seuil <- 0.5

# Sélectionner les colonnes qui respectent ce seuil
colonnes_fournies <- names(taux_remplissage[taux_remplissage >= seuil])

# Filtrer le dataframe
especes_identifiees_clean <- especes_identifiees[, colonnes_fournies]

# Afficher les colonnes conservées et un aperçu
print(colonnes_fournies)
head(especes_identifiees_clean)


library(tidyr)
library(dplyr)
library(ggplot2)
library(FactoMineR)
library(factoextra)

# Extraire la table d'abondance, la taxonomie et les métadonnées de phyloseq
otu_df <- as.data.frame(as(otu_table(ps), "matrix"))
tax_df <- as.data.frame(as(tax_table(ps), "matrix"))
sam_df <- as.data.frame(as(sample_data(ps), "data.frame"))

# Gérer l'orientation de l'OTU table
if(taxa_are_rows(ps)) {
    otu_df$OTU <- rownames(otu_df)
} else {
    otu_df <- as.data.frame(t(otu_df))
    otu_df$OTU <- rownames(otu_df)
}

# Préparer la taxonomie
tax_df$OTU <- rownames(tax_df)
tax_df$espece_clean <- tolower(gsub("_", " ", tax_df[[species_col[1]]]))

# Format long pour lier les OTUs aux échantillons
otu_long <- pivot_longer(otu_df, cols = -OTU, names_to = "Sample", values_to = "Abundance") %>%
    filter(Abundance > 0) # Ne garder que les présences

# Identifier la colonne correspondant à l'organe/tissu dans les métadonnées
col_organe <- grep("(?i)organ|tissue|compartment", colnames(sam_df), value = TRUE)[1]
if(is.na(col_organe)) col_organe <- colnames(sam_df)[1] # Fallback sur la 1ère colonne si introuvable

sam_df$Sample <- rownames(sam_df)

# Croiser les abondances avec les métadonnées (échantillons) et la taxonomie
croisement <- left_join(otu_long, sam_df[, c("Sample", col_organe)], by = "Sample") %>%
    left_join(tax_df[, c("OTU", "espece_clean")], by = "OTU")

# Croiser avec les traits des espèces identifiées
donnees_traits_organes <- left_join(croisement, especes_identifiees_clean, by = "espece_clean") %>%
    filter(!is.na(espece)) # Retirer les espèces sans traits

# 1. Visualisation basique : Répartition du mode de vie (primary_lifestyle) par organe
col_lifestyle <- grep("(?i)lifestyle|guild", colnames(donnees_traits_organes), value = TRUE)[1]

if(!is.na(col_lifestyle)) {
    p_bars <- ggplot(donnees_traits_organes[!is.na(donnees_traits_organes[[col_lifestyle]]), ], 
                                    aes_string(x = col_organe, fill = col_lifestyle)) +
        geom_bar(position = "fill") +
        theme_minimal() +
        labs(y = "Proportion", title = paste("Répartition de", col_lifestyle, "par", col_organe)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    print(p_bars)
}

# Vérifier la présence de la colonne avant de faire les tests
if(!is.na(col_lifestyle)) {
    # Supprimer les NAs pour les deux variables
    df_test <- donnees_traits_organes %>%
        filter(!is.na(.data[[col_organe]]) & !is.na(.data[[col_lifestyle]]))
    
    # Créer une table de contingence (comptage de chaque catégorie)
    table_contingence <- table(df_test[[col_organe]], df_test[[col_lifestyle]])
    
    print("Tableau de contingence :")
    print(table_contingence)
    
    # Test du Chi-2 d'indépendance
    # simulate.p.value = TRUE est utile si les effectifs observés sont petits
    test_chi2 <- chisq.test(table_contingence, simulate.p.value = TRUE)
    
    print("Résultat du test statitisque (Chi-2) :")
    print(test_chi2)
}

# Identifier la colonne correspondant au projet dans les métadonnées
col_projet <- grep("(?i)project|projet", colnames(sam_df), value = TRUE)[1]

if(is.na(col_projet) && ncol(sam_df) > 1) {
    col_projet <- colnames(sam_df)[2] # Fallback
}

if(!is.na(col_projet) && !is.na(col_lifestyle)) {
    # Ajouter l'information de projet au croisement
    croisement_projet <- left_join(otu_long, sam_df[, c("Sample", col_projet)], by = "Sample") %>%
        left_join(tax_df[, c("OTU", "espece_clean")], by = "OTU")
    
    donnees_traits_projets <- left_join(croisement_projet, especes_identifiees_clean, by = "espece_clean") %>%
        filter(!is.na(espece))
    
    # Plot : Répartition du mode de vie par projet
    p_bars_projet <- ggplot(donnees_traits_projets[!is.na(donnees_traits_projets[[col_lifestyle]]), ], 
                            aes_string(x = col_projet, fill = col_lifestyle)) +
        geom_bar(position = "fill") +
        theme_minimal() +
        labs(y = "Proportion", title = paste("Répartition de", col_lifestyle, "par", col_projet)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    print(p_bars_projet)
}

# Identifier la colonne correspondant à la famille de la plante dans les métadonnées
col_famille <- grep("(?i)family|famille", colnames(sam_df), value = TRUE)[1]

if(!is.na(col_famille) && !is.na(col_lifestyle)) {
    # Ajouter l'information de famille au croisement
    croisement_famille <- left_join(otu_long, sam_df[, c("Sample", col_famille)], by = "Sample") %>%
        left_join(tax_df[, c("OTU", "espece_clean")], by = "OTU")
    
    donnees_traits_familles <- left_join(croisement_famille, especes_identifiees_clean, by = "espece_clean") %>%
        filter(!is.na(espece))
    
    # Plot : Répartition du mode de vie par famille de plante
    p_bars_famille <- ggplot(donnees_traits_familles[!is.na(donnees_traits_familles[[col_lifestyle]]), ], 
                            aes_string(x = col_famille, fill = col_lifestyle)) +
        geom_bar(position = "fill") +
        theme_minimal() +
        labs(y = "Proportion", title = paste("Répartition de", col_lifestyle, "par", col_famille)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    print(p_bars_famille)
}

# Vérifier la présence des colonnes avant de faire le test
if(!is.na(col_famille) && !is.na(col_lifestyle)) {
    # Supprimer les NAs pour les deux variables
    df_test_famille <- donnees_traits_familles %>%
        filter(!is.na(.data[[col_famille]]) & !is.na(.data[[col_lifestyle]]))
    
    # Créer une table de contingence (comptage de chaque catégorie par famille)
    table_contingence_famille <- table(df_test_famille[[col_famille]], df_test_famille[[col_lifestyle]])
    
    print("Tableau de contingence (Famille de la plante vs Mode de vie) :")
    print(table_contingence_famille)
    
    # Test du Chi-2 d'indépendance
    test_chi2_famille <- chisq.test(table_contingence_famille, simulate.p.value = TRUE)
    
    print("Résultat du test statistique (Chi-2) par famille de la plante :")
    print(test_chi2_famille)
}
# Création d'une table statistique croisée (Organe vs Mode de vie)
if(!is.na(col_organe) && !is.na(col_lifestyle)) {
    table_croisee <- donnees_traits_organes %>%
        # Retirer les valeurs manquantes pour les deux colonnes ciblées
        filter(!is.na(.data[[col_organe]]) & !is.na(.data[[col_lifestyle]])) %>%
        # Compter le nombre d'occurrences pour chaque paire
        count(.data[[col_organe]], .data[[col_lifestyle]]) %>%
        # Passer en format large pour faciliter la lecture stat
        pivot_wider(
            names_from = .data[[col_lifestyle]], 
            values_from = n, 
            values_fill = list(n = 0)
        )
    
    print("Table statistique croisée : Organes (lignes) x Mode de vie (colonnes)")
    print(as.data.frame(table_croisee))
}

if(!is.na(col_famille) && !is.na(col_lifestyle)) {
    # Extraire les familles uniques
    familles <- unique(df_test_famille[[col_famille]])
    n_fam <- length(familles)
    
    # Initialiser une matrice pour stocker les p-values
    pval_matrix <- matrix(NA, nrow = n_fam, ncol = n_fam, dimnames = list(familles, familles))
    
    # Calculer les p-values pour chaque paire de familles
    for(i in 1:n_fam) {
        for(j in 1:n_fam) {
            if(i == j) {
                pval_matrix[i, j] <- 1
            } else if (i < j) {
                # Sous-ensemble des deux familles
                sub_df <- df_test_famille %>% 
                    filter(.data[[col_famille]] %in% c(familles[i], familles[j]))
                
                # Table de contingence
                tbl <- table(sub_df[[col_famille]], sub_df[[col_lifestyle]])
                
                # Vérifier qu'il y a assez de données pour le test
                if(nrow(tbl) > 1 && ncol(tbl) > 1) {
                    test_pair <- suppressWarnings(chisq.test(tbl, simulate.p.value = TRUE))
                    pval_matrix[i, j] <- test_pair$p.value
                    pval_matrix[j, i] <- test_pair$p.value
                }
            }
        }
    }
    
    # Transformer la matrice en format long pour ggplot2
    pval_long <- as.data.frame(as.table(pval_matrix))
    colnames(pval_long) <- c("Famille1", "Famille2", "p_value")
    
    # Créer une variable pour marquer les p-values significatives (< 0.05)
    pval_long$signif <- ifelse(!is.na(pval_long$p_value) & pval_long$p_value < 0.05, "*", "")
    
    # Visualisation avec une heatmap ggraph
    p_heatmap_pval <- ggplot(pval_long, aes(x = Famille1, y = Famille2, fill = p_value)) +
        geom_tile(color = "white") +
        geom_text(aes(label = signif), color = "black", size = 5) +
        scale_fill_gradientn(colors = c("red", "orange", "white", "lightblue", "blue"), 
                              values = scales::rescale(c(0, 0.05, 0.1, 0.5, 1)),
                              limits = c(0, 1),
                              name = "P-value", na.value = "grey80") +
        theme_minimal() +
        labs(title = "Heatmap des P-values : Différence de composition en mode de vie (Famille vs Famille)",
             x = "Famille de la plante", y = "Famille de la plante") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    print(p_heatmap_pval)
}
