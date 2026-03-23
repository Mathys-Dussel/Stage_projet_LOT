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


sample_data(ps)$LOT_sampleID <- gsub("/.*", "", sample_data(ps)$LOT_sampleID)



library(dplyr)

zones_clean <- ZONES_LOT_123 %>%
  select(LOT.SampleCode, TreeZone) %>%
  distinct(LOT.SampleCode, .keep_all = TRUE)

metadata <- data.frame(sample_data(ps))

new_metadata <- metadata %>%
  left_join(zones_clean, by = c("LOT_sampleID" = "LOT.SampleCode"))

rownames(new_metadata) <- rownames(metadata)

sample_data(ps) <- sample_data(new_metadata)



ZONES_LOT_123 %>% 
  count(LOT.SampleCode) %>% 
  filter(n > 1)

