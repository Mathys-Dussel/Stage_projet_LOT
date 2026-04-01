setwd("~/Documents/Etudes/Stage_projet_LOT/CRBE/Analyses")

ps <- readRDS("donnees/ps_final.rds")

MC_LOT_1=read.csv("donnees/Macroclimat_LOT1.csv", header = TRUE, sep = ";", dec = ",")
head(MC_LOT_1)

MC_TEMP=MC_LOT_1[,c("T_ZONE1","T_ZONE2","T_ZONE3","T_ZONE4","T_ZONE5","T_ZONE6")]
MC_HUMI=MC_LOT_1[,c("H_ZONE1","H_ZONE2","H_ZONE3","H_ZONE4","H_ZONE5","H_ZONE6")]

# bottom (Z1)
# top of the trunk (Z2–Z3)
# the heart (Z4)
# the periphery of the canopy (Z5-Z6)

MC_TEMP$top_trunk <- rowMeans(MC_TEMP[, c("T_ZONE2", "T_ZONE3")], na.rm = TRUE)
MC_TEMP <- MC_TEMP[, !(names(MC_TEMP) %in% c("T_ZONE2", "T_ZONE3"))]

MC_HUMI$top_trunk <- rowMeans(MC_HUMI[, c("H_ZONE2", "H_ZONE3")], na.rm = TRUE)
MC_HUMI <- MC_HUMI[, !(names(MC_HUMI) %in% c("H_ZONE2", "H_ZONE3"))]

MC_TEMP$periph_canopy <- rowMeans(MC_TEMP[, c("T_ZONE5", "T_ZONE6")], na.rm = TRUE)
MC_TEMP <- MC_TEMP[, !(names(MC_TEMP) %in% c("T_ZONE5", "T_ZONE6"))]

MC_HUMI$periph_canopy <- rowMeans(MC_HUMI[, c("H_ZONE5", "H_ZONE6")], na.rm = TRUE)
MC_HUMI <- MC_HUMI[, !(names(MC_HUMI) %in% c("H_ZONE5", "H_ZONE6"))]

colnames(MC_TEMP)[colnames(MC_TEMP) == "T_ZONE4"] <- "Heart"
colnames(MC_HUMI)[colnames(MC_HUMI) == "H_ZONE4"] <- "Heart"

colnames(MC_TEMP)[colnames(MC_TEMP) == "T_ZONE1"] <- "Bottom"
colnames(MC_HUMI)[colnames(MC_HUMI) == "H_ZONE1"] <- "Bottom"

boxplot(MC_TEMP)
boxplot(MC_HUMI)

head(MC_TEMP)
head(MC_HUMI)

library(tidyr)
library(ggpubr)

MC_TEMP_long <- pivot_longer(MC_TEMP, cols = everything(), names_to = "Zone", values_to = "Temperature")
MC_HUMI_long <- pivot_longer(MC_HUMI, cols = everything(), names_to = "Zone", values_to = "Humidity")

my_comparisons <- list(
    c("Bottom", "Heart"), 
    c("Bottom", "top_trunk"), 
    c("Bottom", "periph_canopy"), 
    c("Heart", "top_trunk"), 
    c("Heart", "periph_canopy"), 
    c("top_trunk", "periph_canopy")
)

p_temp <- ggboxplot(MC_TEMP_long, x = "Zone", y = "Temperature", fill = "Zone") +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
    stat_compare_means(label.y = max(MC_TEMP_long$Temperature, na.rm = TRUE) + 5) +
    ggtitle("Différences de Température entre les zones")

p_humi <- ggboxplot(MC_HUMI_long, x = "Zone", y = "Humidity", fill = "Zone") +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
    stat_compare_means(label.y = max(MC_HUMI_long$Humidity, na.rm = TRUE) + 10) +
    ggtitle("Différences d'Humidité entre les zones")

print(p_temp)
print(p_humi)



ZONES_LOT_123=read.csv("donnees/LOT_zones.csv", header = TRUE, sep = ";", dec = ",")
head(ZONES_LOT_123)

ZONES_LOT_123$TreeZone[ZONES_LOT_123$TreeZone == 3] <- 2
ZONES_LOT_123$TreeZone[ZONES_LOT_123$TreeZone == 6] <- 5

ZONES_LOT_123$TreeZone <- as.character(ZONES_LOT_123$TreeZone)
ZONES_LOT_123$TreeZone[ZONES_LOT_123$TreeZone == "1"] <- "Bottom"
ZONES_LOT_123$TreeZone[ZONES_LOT_123$TreeZone == "2"] <- "Heart"
ZONES_LOT_123$TreeZone[ZONES_LOT_123$TreeZone == "4"] <- "top_trunk"
ZONES_LOT_123$TreeZone[ZONES_LOT_123$TreeZone == "5"] <- "periph_canopy"

head(ZONES_LOT_123)

ZONES_LOT_123$LOT.SampleCode  <- gsub("/.*", "", ZONES_LOT_123$LOT.SampleCode)












library(phyloseq)
library(ggplot2)
library(ggpubr)
library(patchwork)

ps <- readRDS("donnees/ps_final.rds")
sample_names(ps) <- gsub("\\.", "-", sub("^X", "", sample_names(ps)))

ps_LOT1 <- subset_samples(ps, startsWith(LOT_sampleID, "LOT01"))

ps_LOT1 <- prune_taxa(taxa_sums(ps_LOT1) > 0, ps_LOT1)



sample_data(ps_LOT1)$Zone <- "Unknown"
sample_data(ps_LOT1)$Zone[substr(sample_data(ps_LOT1)$LOT_sampleID, 9, 9) == "1"] <- "Bottom"
sample_data(ps_LOT1)$Zone[substr(sample_data(ps_LOT1)$LOT_sampleID, 9, 9) %in% c("2", "3")] <- "top_trunk"
sample_data(ps_LOT1)$Zone[substr(sample_data(ps_LOT1)$LOT_sampleID, 9, 9) == "4"] <- "Heart"
sample_data(ps_LOT1)$Zone[substr(sample_data(ps_LOT1)$LOT_sampleID, 9, 9) %in% c("5", "6")] <- "periph_canopy"

ps_LOT1_zones <- subset_samples(ps_LOT1, Zone != "Unknown")
richness <- estimate_richness(ps_LOT1_zones, measures = c("Observed", "Shannon"))
meta_data <- sample_data(ps_LOT1_zones)
meta_data$Observed <- richness$Observed
meta_data$ExpShannon <- exp(richness$Shannon)

df_richness <- data.frame(
    Sample = rownames(meta_data),
    Zone = meta_data$Zone,
    Observed = meta_data$Observed,
    ExpShannon = meta_data$ExpShannon
)

df_richness$Zone <- factor(df_richness$Zone, levels = c("Bottom", "top_trunk", "Heart", "periph_canopy"))

p_richness <- ggboxplot(df_richness, x = "Zone", y = "Observed", fill = "Zone") +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
    ggtitle("Différences de Richesse Spécifique (Observed) entre les zones")

p_expshannon <- ggboxplot(df_richness, x = "Zone", y = "ExpShannon", fill = "Zone") +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
    ggtitle("Différences de Diversité (ExpShannon) entre les zones")

print(p_richness)
print(p_expshannon)

kruskal.test(Observed ~ Zone, data = df_richness)
pairwise.wilcox.test(df_richness$Observed, df_richness$Zone, p.adjust.method = "fdr")

kruskal.test(ExpShannon ~ Zone, data = df_richness)
pairwise.wilcox.test(df_richness$ExpShannon, df_richness$Zone, p.adjust.method = "fdr")

library(dplyr)

mean_temp_zone <- MC_TEMP_long %>%
    group_by(Zone) %>%
    summarise(Mean_Temp = mean(Temperature, na.rm = TRUE))

mean_humi_zone <- MC_HUMI_long %>%
    group_by(Zone) %>%
    summarise(Mean_Humi = mean(Humidity, na.rm = TRUE))

df_richness <- merge(df_richness, mean_temp_zone, by = "Zone", all.x = TRUE)
df_richness <- merge(df_richness, mean_humi_zone, by = "Zone", all.x = TRUE)

p_corr_temp_obs <- ggscatter(df_richness, x = "Mean_Temp", y = "Observed",
                         add = "reg.line", conf.int = TRUE, 
                         cor.coef = TRUE, cor.method = "spearman",
                         xlab = "Température moyenne de la zone", 
                         ylab = "Richesse Spécifique (Observed)") +
    ggtitle("Relation entre la Température et la Richesse (Observed)")

p_corr_humi_obs <- ggscatter(df_richness, x = "Mean_Humi", y = "Observed",
                         add = "reg.line", conf.int = TRUE, 
                         cor.coef = TRUE, cor.method = "spearman",
                         xlab = "Humidité moyenne de la zone", 
                         ylab = "Richesse Spécifique (Observed)") +
    ggtitle("Relation entre l'Humidité et la Richesse (Observed)")

print(p_corr_temp_obs)
print(p_corr_humi_obs)

p_corr_temp_exp <- ggscatter(df_richness, x = "Mean_Temp", y = "ExpShannon",
                             add = "reg.line", conf.int = TRUE, 
                             cor.coef = TRUE, cor.method = "spearman",
                             xlab = "Température moyenne de la zone", 
                             ylab = "Diversité (ExpShannon)") +
    ggtitle("Relation entre la Température et la Diversité (ExpShannon)")

p_corr_humi_exp <- ggscatter(df_richness, x = "Mean_Humi", y = "ExpShannon",
                             add = "reg.line", conf.int = TRUE, 
                             cor.coef = TRUE, cor.method = "spearman",
                             xlab = "Humidité moyenne de la zone", 
                             ylab = "Diversité (ExpShannon)") +
    ggtitle("Relation entre l'Humidité et la Diversité (ExpShannon)")

print(p_corr_temp_exp)
print(p_corr_humi_exp)


library(MASS)

glm_humi_temp_obs <- glm.nb(Observed ~ Mean_Humi * Mean_Temp, data = df_richness)
summary(glm_humi_temp_obs)

glm_humi_temp_exp <- glm.nb(ExpShannon ~ Mean_Humi * Mean_Temp, data = df_richness)
summary(glm_humi_temp_exp)



library(DHARMa)

simulationOutput_humi_temp_obs <- simulateResiduals(fittedModel = glm_humi_temp_obs, plot = TRUE)
testUniformity(simulationOutput_humi_temp_obs)
testDispersion(simulationOutput_humi_temp_obs)
testZeroInflation(simulationOutput_humi_temp_obs)

simulationOutput_humi_temp_exp <- simulateResiduals(fittedModel = glm_humi_temp_exp, plot = TRUE)
testUniformity(simulationOutput_humi_temp_exp)
testDispersion(simulationOutput_humi_temp_exp)
testZeroInflation(simulationOutput_humi_temp_exp)


library(interactions)
interact_plot(glm_humi_temp_obs, pred = Mean_Humi, modx = Mean_Temp)
interact_plot(glm_humi_temp_exp, pred = Mean_Humi, modx = Mean_Temp)















library(phyloseq)
library(ggplot2)
library(ggpubr)
library(patchwork)

ps <- readRDS("donnees/ps_final.rds")

sample_data(ps)$LOT_sampleID <- gsub("/.*", "", sample_data(ps)$LOT_sampleID)

env_data <- read.csv("donnees/LOT_zones.csv", header = TRUE, sep = ";", dec = ",")
env_data$LOT.SampleCode <- gsub("/.*", "", env_data$LOT.SampleCode)

sam_data <- as(sample_data(ps), "data.frame")

sam_data$TreeZone <- env_data$TreeZone[match(sam_data$LOT_sampleID, env_data$LOT.SampleCode)]

sample_data(ps) <- sample_data(sam_data)

sum(sam_data$LOT_sampleID %in% env_data$LOT.SampleCode)

head(sample_data(ps)$TreeZone)


data_matche <- merge(sam_data, env_data, by.x = "LOT_sampleID", by.y = "LOT.SampleCode")

