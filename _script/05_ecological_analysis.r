# 05. ecological analysis

source("utils_ecology_myxo.R")
load_ecology_packages()

## 1. Correlation analysis for Figure S2
### 1-1. Data preparation
BioClim_data <- readRDS("../_data/_processed_data/BioClim_data.rds")
GeoCoord_details <- readRDS("../_data/_processed_data/GeoCoord_avg_11regions_detail.rds")
str(BioClim_data)

BioClim_data <- BioClim_data %>%
  rename(
    Northing = Northing_avg,
    BIO1 = v1_Mean_temp,
    BIO2 = v2_Mean_diurnal_range,
    BIO3 = v3_Isothermality,
    BIO4 = v4_Temp_Seasonality,
    BIO5 = v5_Max_Temp_Warmest_Month,
    BIO6 = v6_Min_Temp_Coldest_Month,
    BIO7 = v7_Temp_Annual_Range,
    BIO8 = v8_Mean_Temp_Wettest_Quarter,
    BIO9 = v9_Mean_Temp_Driest_Quarter,
    BIO10 = v10_Mean_Temp_Warmest_Quarter,
    BIO11 = v11_Mean_Temp_Coldest_Quarter,
    BIO12 = v12_Annual_Precipit,
    BIO13 = v13_Precipit_Wettest_Month,
    BIO14 = v14_Precipit_Driest_Month,
    BIO15 = v15_Precipit_Seasonality,
    BIO16 = v16_Precipit_Wettest_Quarter,
    BIO17 = v17_Precipit_Driest_Quarter,
    BIO18 = v18_Precipit_Warmest_Quarter,
    BIO19 = v19_Precipit_Coldest_Quarter
  ) %>% 
  mutate(Code_Subplot = paste(Code, Subplot_Nr, sep = "_")) %>% 
  dplyr::select(-Subplot, - Easting_avg) 

GeoCoord_elev <- GeoCoord_details %>% 
  dplyr::select(Code, Subplot_Nr, Elevation_avg) %>% 
  mutate(Code_Subplot = paste(Code, Subplot_Nr, sep = "_"))

colnames(GeoCoord_elev)[colnames(GeoCoord_elev) == "Elevation_avg"] <- "Elev"

BioClim_variables <- left_join(BioClim_data, GeoCoord_elev)
BioClim_variables <- dplyr::select(BioClim_variables, -Code_Subplot)
BioClim_variables$Subplot_Nr <- c("Khi2", "Khi", "Khi1", "Nor3", "Nor2", "Nor", "Nor1", "Kam2", "Kam", "Kam1",
                             "Tat", "B/V3", "B/V", "B/V2", "B/V1", "GAP", "Fr", "Cau", "Pyr", "SN",
                             "Jp5",  "Jp4",  "Jp3",  "Jp", "Jp2",  "Jp1")

### 1-2. Run PCA
# only bioclimatic + geographic variables
results_PCA <- prcomp(BioClim_variables[, -c(1:2)], scale = TRUE)
summary(results_PCA)

plot(results_PCA$x, type = "n")
text(results_PCA$x, labels = BioClim_variables$Subplot_Nr, col = "blue", cex = 1.5)
biplot(results_PCA)

### 1-3. Plot PCA
#### 1-3-1. convert as dataframe to plot with ggplot
pca_df <- data.frame(
  PC1 = results_PCA$x[, 1],
  PC2 = results_PCA$x[, 2],
  Subplot = BioClim_variables$Subplot_Nr
)

sdevs <- results_PCA$sdev
centers <- results_PCA$center   
contribRatio <- sdevs^2 / sum(sdevs^2) 
DF_loadings <- as.data.frame(results_PCA$rotation[, 1:2])
pca_scores <- results_PCA$x
pca_scores_df <- as.data.frame(pca_scores)
pca_scores_df$Subplot <- BioClim_variables$Subplot_Nr

xRatio <- (max(results_PCA$x[, "PC1"]) - min(results_PCA$x[, "PC1"])) / (max(DF_loadings[, "PC1"]) - min(DF_loadings[, "PC1"]))
yRatio <- (max(results_PCA$x[, "PC2"]) - min(results_PCA$x[, "PC2"])) / (max(DF_loadings[, "PC2"]) - min(DF_loadings[, "PC2"]))
vars <- c("Northing", "BIO1", "BIO2", "BIO3", 
          "BIO4", "BIO5", "BIO6", 
          "BIO7", "BIO8", "BIO9", 
          "BIO10", "BIO11", "BIO12", 
          "BIO13", "BIO14", "BIO15", 
          "BIO16", "BIO17", "BIO18", 
          "BIO19", "SC", "SC_lasting", "Elev")

#### 1-3-2. Plot
library(ggalt) # geom_encircle()

pca <- ggplot(data.frame(results_PCA$x, Subplot = as.factor(pca_scores_df$Subplot))) +　
  theme_bw(base_size = 20) +
  geom_segment(data = data.frame(PC1 = xRatio * DF_loadings$PC1, PC2 = yRatio * DF_loadings$PC2),
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.01, "npc")), lwd = 0.2) +
  geom_text(data = data.frame(PC1 = xRatio * DF_loadings$PC1, PC2 = yRatio * DF_loadings$PC2),
            aes(x = PC1, y = PC2 + 0.5, label = vars), size=5, color = "steelblue") +
  geom_encircle(data = subset(pca_df, 
                              Subplot == "Kam"|
                                Subplot == "Kam1"|
                                Subplot == "Kam2"), 
                aes(x = PC1, y = PC2), 
                color = "black", 
                alpha = 0.5, 
                size = 1) +    
  geom_encircle(data = subset(pca_df, Subplot %in% c("Kam", "Kam1", "Kam2")), 
                aes(x = PC1, y = PC2), 
                color = "black", 
                alpha = 0.5, 
                size = 1) +
  geom_encircle(data = subset(pca_df, Subplot %in% c("Khi", "Khi1", "Khi2")), 
                aes(x = PC1, y = PC2), 
                color = "black", 
                alpha = 0.5, 
                size = 1) +
  geom_encircle(data = subset(pca_df, Subplot %in% c("Nor", "Nor1", "Nor2", "Nor3")), 
                aes(x = PC1, y = PC2), 
                color = "black", 
                alpha = 0.5, 
                size = 1) +
  geom_encircle(data = subset(pca_df, Subplot %in% c("B/V", "B/V1", "B/V2", "B/V3")), 
                aes(x = PC1, y = PC2), 
                color = "black", 
                alpha = 0.5, 
                size = 1) +
  geom_encircle(data = subset(pca_df, Subplot %in% c("Jp", "Jp1", "Jp2", "Jp3", "Jp4", "Jp5")), 
                aes(x = PC1, y = PC2), 
                color = "black", 
                alpha = 0.5, 
                size = 1) +
  geom_text(data = subset(pca_df, Subplot %in% c("Nor", "GAP", "B/V", "Tat", "Kam", "Jp", "Khi", "Cau", "SN", "Fr", "Pyr")), 
            aes(x = PC1, y = PC2, label = Subplot), 
            size = 7, 
            fontface = "bold", 
            color = "black") +
  geom_text(data = subset(pca_df, !Subplot %in% c("Nor", "GAP", "B/V", "Tat", "Kam", "Jp", "Khi", "Cau", "SN", "Fr", "Pyr")), 
            aes(x = PC1, y = PC2, label = Subplot), 
            size = 4, 
            color = "black") +
  theme(
    panel.grid.major = element_blank(),   
    panel.grid.minor = element_blank(),   
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color = "black"), 
    axis.title = element_text(color = "black"), 
    axis.ticks = element_blank(),         
    legend.position = "none"  
  ) +
  annotate("text", x = min(pca_df$PC1), y = 5, label = "A", size = 15) +
  labs(title = "",
       x = "PC 1",
       y = "PC 2") 
pca

png("../_results/FigS2A_PCA.png", width = 800, height = 800)
print(pca)
dev.off()

### 1-4. Run nMDS
#### 1-4-1. Scaling for nMDS
BioClim_variables_scale <- as.data.frame(scale(BioClim_variables[, -c(1:2)]))

# distance calculation
dist_data <- vegdist(BioClim_variables_scale, method = "euclidian", diag = FALSE, upper = FALSE)

#### 1-4-2. Run nMDS
set.seed(123) 
result.mds <- metaMDS(dist_data, k = 2, trymax = 9999)
result.mds

plot(result.mds$points, type="n", main=paste("Euclidian/NMDS - Stress =", round(result.mds$stress, 3), ": with Coordinates"), cex.main=1.5, cex.lab=1.2)
text(result.mds$points, labels = BioClim_variables$Subplot_Nr, cex=1.5) #rownames(spe)

### 1-5. Plot nMDS
#### 1-5-1. convert as dataframe to plot with ggplot
nmds_df <- data.frame(
  MDS1 = result.mds$points[, 1],
  MDS2 = result.mds$points[, 2],
  Subplot = BioClim_variables$Subplot_Nr
)

p2 <- ggplot(nmds_df, aes(x = MDS1, y = MDS2, label = Subplot)) + 
  theme_bw(base_size = 20) + 
  geom_encircle(data = subset(nmds_df, 
                              Subplot == "Kam"|
                                Subplot == "Kam1"|
                                Subplot == "Kam2"), 
                aes(x = MDS1, y = MDS2),                 
                color = "black", 
                alpha = 0.5, 
                size = 1) +    # Regionを囲む
  geom_encircle(data = subset(nmds_df, Subplot %in% c("Kam", "Kam1", "Kam2")), 
                aes(x = MDS1, y = MDS2), 
                color = "black", 
                alpha = 0.5, 
                size = 1) +
  geom_encircle(data = subset(nmds_df, Subplot %in% c("Khi", "Khi1", "Khi2")), 
                aes(x = MDS1, y = MDS2),  
                color = "black", 
                alpha = 0.5, 
                size = 1) +
  geom_encircle(data = subset(nmds_df, Subplot %in% c("Nor", "Nor1", "Nor2", "Nor3")), 
                aes(x = MDS1, y = MDS2), 
                color = "black", 
                alpha = 0.5, 
                size = 1) +
  geom_encircle(data = subset(nmds_df, Subplot %in% c("B/V", "B/V1", "B/V2", "B/V3")), 
                aes(x = MDS1, y = MDS2), 
                color = "black", 
                alpha = 0.5, 
                size = 1) +
  geom_encircle(data = subset(nmds_df, Subplot %in% c("Jp", "Jp1", "Jp2", "Jp3", "Jp4", "Jp5")), 
                aes(x = MDS1, y = MDS2), 
                color = "black", 
                alpha = 0.5, 
                size = 1) +
  # Region with bold
  geom_text(data = subset(nmds_df, Subplot %in% c("Nor", "GAP", "B/V", "Tat", "Kam", "Jp", "Khi", "Cau", "SN", "Fr", "Pyr")), 
            aes(x = MDS1, y = MDS2, label = Subplot), 
            size = 7, 
            fontface = "bold", 
            color = "black") +
  # Region nomal size
  geom_text(data = subset(nmds_df, !Subplot %in% c("Nor", "GAP", "B/V", "Tat", "Kam", "Jp", "Khi", "Cau", "SN", "Fr", "Pyr")), 
            aes(x = MDS1, y = MDS2, label = Subplot), 
            size = 4, 
            color = "black") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"  
  ) +
  annotate("text", x = min(nmds_df$MDS1), y = 4.5, label = "B", size = 15, hjust = 0) +
  annotate("text", x = min(nmds_df$MDS1), y = -3, label = paste("Stress =", round(result.mds$stress, 3)), size = 6, hjust = 0) +
  labs(title = "",
       x = "MDS 1",
       y = "MDS 2")
p2

png("../_results/FigS2B_nMDS.png", width = 800, height = 800)
print(p2)
dev.off()

FigS2 <- pca|p2
FigS2
png("../_results/FigS2AB_pca_nMDS.png", width = 1200, height = 600)
print(FigS2)
dev.off()


ggsave(
  filename = "../_results/FigS2AB_pca_nMDS.pdf", 
  plot = FigS2,                    
  width = 20,                          
  height = 10,                         
  dpi = 300                           
)

## 2. Correlation plot for Figure S3A
### 2-1. Calculate peason r
#### 2-1-1. Data preparation
BioClim_variables_final  <-  filter (BioClim_variables, 
                                     Subplot_Nr == "Nor"| 
                                       Subplot_Nr == "GAP"|
                                       Subplot_Nr == "B/V" |
                                       Subplot_Nr == "Tat" |
                                       Subplot_Nr == "Kam1" |
                                       Subplot_Nr == "Kam2" |
                                       Subplot_Nr == "Jp" |
                                       Subplot_Nr == "Khi" |
                                       Subplot_Nr == "Cau" |
                                       Subplot_Nr == "SN" |
                                       Subplot_Nr == "Fr" |
                                       Subplot_Nr == "Pyr" )

BioClim_step1 <- BioClim_variables_final[, -c(1:2)]

##### 2-1-2. Plot with psch
psych::pairs.panels(BioClim_step1, cex.cor = 1.5)


pdf("../_results/FigS3A_correlation.pdf", width=16, height= 16)
psych::pairs.panels(BioClim_step1, cex.cor = 1.5)
dev.off()

##### 2-1-2. Plot with qraph
cor_matrix <- cor(BioClim_step1)

library(qgraph)
qgraph(cor_matrix,
       minimum = 0.8,
       labels = colnames(cor_matrix),layout = "groups",
       edge.labels = T,
       label.scale = F,
       label.cex = 1,
       edge.label.cex = 1.2)

pdf("../_results/FigS3B_correlation.pdf", width=10, height= 10)
qgraph(cor_matrix,
       minimum = 0.8,
       labels = colnames(cor_matrix),layout = "groups",
       edge.labels = T,
       label.scale = F,
       label.cex = 1,
       edge.label.cex = 1.2)
dev.off()

## 3. db-RDA for Figure 8 and S2, Table 2 and S4
### 3-1. Data preparation
#### 3-1-1. Bioclimatic variable selection for set I and II

# Set I

# Step I. Obtain correlated variable (r>0.7) and exclude automatically with carret
library(caret)
cor_matrix <- cor(BioClim_step1)
high_cor_vars <- findCorrelation(cor_matrix, cutoff = 0.7, verbose = TRUE)
print(high_cor_vars)

# Exclude the variable with r>0.7
BioClim_step2 <- BioClim_step1[, - high_cor_vars]
BioClim_step2

# Step II: Obtain VIF for the remaining variables (VIF<10) 
BioClim_SetI <- cbind(BioClim_variables_final[, 2], BioClim_step2)
vif_BioClim_1 <- diag(solve(cor(BioClim_SetI[, -1]))) # exclude region name
vif_BioClim_1

#       BIO2       BIO3       BIO5       BIO7       BIO8      BIO15      BIO18         SC 
# 665.987369 484.094660   4.234153 497.752697   9.622425   4.811592   4.105455   2.718498 

# stepwise exclude
# BIO2
vif_BioClim_1e <- diag(solve(cor(BioClim_SetI[, -c(1:2)])))
vif_BioClim_1e
#     BIO3     BIO5     BIO7     BIO8    BIO15    BIO18       SC 
# 2.905275 4.071619 3.113773 3.614981 2.240061 2.993443 2.538156 
# all are VIF<10

BioClim_SetI_final <- BioClim_SetI[, -2]

# StepIII: Standardization with scale()
BioClim_SetI_final_scale <- data.frame(scale(BioClim_SetI_final[, -1]))
rownames(BioClim_SetI_final_scale) <- BioClim_variables_final$Subplot_Nr

# Step IV: Obtain correlation r
psych::pairs.panels(BioClim_SetI_final_scale, cex.cor = 1)

# Set II (based on B. albescens paper)

# Step I: Set initial variables
# Elev
# v2_Mean_diurnal_range
# v3_Isothermality
# v5_Max_Temp_Warmest_Month
# v8_Mean_Temp_Wettest_Quarter
# v9_Mean_Temp_Driest_Quarter
# v15_Precipit_Seasonality
# v18_Precipit_Warmest_Quarter
# v19_Precipit_Coldest_Quarter
# SC
# SC_lasting
vif_BioClim_2 <- diag(solve(cor(BioClim_step1[, c("Elev", "BIO2", "BIO3", "BIO5", "BIO8", "BIO9",
                                                "BIO15", "BIO18", "BIO19", "SC", "SC_lasting")])))
vif_BioClim_2
#       Elev       BIO2       BIO3       BIO5       BIO8       BIO9      BIO15      BIO18      BIO19 
# 946.89384  578.90350 1341.72567  964.18516  221.07624 2450.19021  655.21135  424.90190 3397.20550 
# SC SC_lasting 
# 830.01066   55.65247 

# Step II: Obtain VIF for the remaining variables (VIF<10) 
# stepwise exclude, looking at VIF <10

# v19_Precipit_Coldest_Quarter
vif_BioClim_2e <- diag(solve(cor(BioClim_step1[, c("Elev", "BIO2", "BIO3", "BIO5", "BIO8", "BIO9",
                                                   "BIO15", "BIO18", "SC", "SC_lasting")])))
vif_BioClim_2e
#       Elev       BIO2       BIO3       BIO5       BIO8       BIO9      BIO15      BIO18         SC 
# 44.919477  15.522447  60.771276  41.660668  11.244557  69.036138   3.519183   3.175313  14.572273 
# SC_lasting 
# 9.796270

# v9_Mean_Temp_Driest_Quarter
vif_BioClim_2e <- diag(solve(cor(BioClim_step1[, c("Elev", "BIO2", "BIO3", "BIO5", "BIO8", 
                                                   "BIO15", "BIO18", "SC", "SC_lasting")])))
vif_BioClim_2e
#       Elev       BIO2       BIO3       BIO5       BIO8      BIO15      BIO18         SC SC_lasting 
# 14.133838   7.187687  15.389832   7.097831   4.011214   2.449482   3.111383  14.322547   9.382128

# v3_Isothermality
vif_BioClim_2e <- diag(solve(cor(BioClim_step1[, c("Elev", "BIO2", "BIO5", "BIO8",
                                                  "BIO15", "BIO18", "SC", "SC_lasting")])))
vif_BioClim_2e
#       Elev       BIO2       BIO5       BIO8      BIO15      BIO18         SC SC_lasting 
# 3.623953   5.014756   5.082743   3.230490   2.314420   3.108444   5.617571   4.286452 
# all are VIF<10

BioClim_step2_II <- BioClim_step1[, c("Elev", "BIO2", "BIO5", "BIO8",
                     "BIO15", "BIO18", "SC", "SC_lasting")]
BioClim_SetII <- cbind(BioClim_variables_final[, 2], 
                       BioClim_step2_II)

# Step III: Standardization with scale()
BioClim_SetII_scale <- data.frame(scale(BioClim_SetII[, -1]))
rownames(BioClim_SetII_scale) <- BioClim_variables_final$Subplot_Nr

# Step IV: Obtain correlation r
psych::pairs.panels(BioClim_SetII_scale, cex.cor = 1)

# Step V: Exclude, SC (vs SC_lasting) and BIO5 (vs Elev) (r>0.7)
vif_BioClim_2e <- diag(solve(cor(BioClim_step1[, c("Elev", "BIO2", "BIO8",
                                                   "BIO15", "BIO18", "SC_lasting")])))
vif_BioClim_2e
#       Elev       BIO2       BIO8      BIO15      BIO18 SC_lasting 
# 2.844455   4.047833   2.973797   2.310569   2.989746   1.578292 
# all are VIF<10

BioClim_SetII_final <- BioClim_step1[, c("Elev", "BIO2",  "BIO8",
                                      "BIO15", "BIO18", "SC_lasting")]
# repeat Step II to IV
BioClim_SetII_final <- cbind(BioClim_variables_final[, 2], 
                         BioClim_SetII_final)
BioClim_SetII_final_scale <- data.frame(scale(BioClim_SetII_final[, -1]))
rownames(BioClim_SetII_final_scale) <- BioClim_variables_final$Subplot_Nr
psych::pairs.panels(BioClim_SetII_final_scale, cex.cor = 1)

#### 3-1-2. Geographic distance data
Geo_dist_detail <- dplyr::select(GeoCoord_details, Code, Subplot_Nr, Easting_avg, Northing_avg) 
Geo_dist_detail$Subplot_Nr <- c("Khi2", "Khi", "Khi1", "Nor3", "Nor2", "Nor", "Nor1", "Kam2", "Kam", "Kam1",
                                  "Tat", "B/V3", "B/V", "B/V2", "B/V1", "GAP", "Fr", "Cau", "Pyr", "SN",
                                  "Jp5",  "Jp4",  "Jp3",  "Jp", "Jp2",  "Jp1")
Geo_dist_detail <- Geo_dist_detail %>% 
  filter(Subplot_Nr == "Nor"| 
                            Subplot_Nr == "GAP"|
                            Subplot_Nr == "B/V" |
                            Subplot_Nr == "Tat" |
                            Subplot_Nr == "Kam1" |
                            Subplot_Nr == "Kam2" |
                            Subplot_Nr == "Jp" |
                            Subplot_Nr == "Khi" |
                            Subplot_Nr == "Cau" |
                            Subplot_Nr == "SN" |
                            Subplot_Nr == "Fr" |
                            Subplot_Nr == "Pyr" ) %>% 
  dplyr::select(-Code)

Geo_dist_detail$Easting_avg <- round(Geo_dist_detail$Easting_avg, 5)
colnames(Geo_dist_detail)[colnames(Geo_dist_detail) == "Northing_avg"] <- "Northing"
colnames(Geo_dist_detail)[colnames(Geo_dist_detail) == "Easting_avg"] <- "Easting"

library(geosphere) #haversine distance
dis.hav_Coord_detail <- distm(Geo_dist_detail[, -1], fun = distHaversine) #haversine method
dimnames(dis.hav_Coord_detail) <- list(Geo_dist_detail$Subplot_Nr, Geo_dist_detail$Subplot_Nr)
#dis.hav_Coord_detail

# Standardization
dis.hav_Coord_stand <- scale(dis.hav_Coord_detail)
#dis.hav_Coord_stand

#### 3-1-3. Community data, Unifrac data
community_RT_rda <- readRDS("../_data/_processed_data/community_RT_rda.rds")
dist_d7a_rt_rda <- readRDS("../_results/beta_diversity/D7_weighted_rda.rds")


### 3-2. db-RDA: BioClim variables vs. geographic distance
#### 3-2-1. Calculation
# Set I
rda_resultI_geodis <- dbrda(dis.hav_Coord_stand ~ ., BioClim_SetI_final_scale)
rda_resultI_geodis

plot(rda_resultI_geodis, scaling =1, 
     main = "Triplot dbRDA Region: - Set I: Geographic Distance (standerized)", cex.main = 1)

# Set II
rda_resultII_geodis <- dbrda(dis.hav_Coord_stand ~ ., BioClim_SetII_final_scale)
rda_resultII_geodis

plot(rda_resultII_geodis, scaling =1, 
     main = "Triplot dbRDA Region: - Set II: Geographic Distance (standerized)", cex.main = 1)

#### 3-2-2. Permutation test: BioClim variables vs. geographic distance
set.seed(123) 
anova_rda_resultI_geodis <- anova(rda_resultI_geodis, by="terms", permutations = 999)
anova_rda_resultI_geodis

set.seed(123) 
anova_rda_resultII_geodis <- anova(rda_resultII_geodis, by="terms", permutations = 999)
anova_rda_resultII_geodis

#### 3-2-3. Significance Test: BioClim variables vs. geographic distance
# set I
reduced_modelI_geodis <- dbrda(dis.hav_Coord_stand ~ 1, data = BioClim_SetI_final_scale)
full_modelI_geodis <- dbrda(dis.hav_Coord_stand ~ ., data = BioClim_SetI_final_scale)
# Permutation Test
set.seed(123) 
anova_comparisonI_geodis <- anova(full_modelI_geodis, reduced_modelI_geodis, permutations = 999)
anova_comparisonI_geodis
# p-value
p_valueI_geodis <- anova_comparisonI_geodis$"Pr(>F)"[2]
p_valueI_geodis

# set II
reduced_modelII_geodis <- dbrda(dis.hav_Coord_stand ~ 1, data = BioClim_SetII_final_scale)
full_modelII_geodis <- dbrda(dis.hav_Coord_stand ~ ., data = BioClim_SetII_final_scale)
# Permutation Test
set.seed(123) 
anova_comparisonII_geodis <- anova(full_modelII_geodis, reduced_modelII_geodis, permutations = 999)
anova_comparisonII_geodis
# p-value
p_valueII_geodis <- anova_comparisonII_geodis$"Pr(>F)"[2]
p_valueII_geodis


### 3-3. db-RDA: BioClim variables vs. D4a, D5, D7a
#### 3-3-1. Calculation
# Set I
# D4a
rda_resultI_d4a <- dbrda(community_RT_rda ~ ., BioClim_SetI_final_scale, distance = "jaccard")
rda_resultI_d4a
plot(rda_resultI_d4a, scaling =1, 
     main = "Triplot dbRDA Region: - Set I: Weighted Jaccard Dissimilarity", cex.main = 1)

# D5
rda_resultI_d5 <- dbrda(community_RT_rda ~ ., BioClim_SetI_final_scale, distance = "horn")
rda_resultI_d5
plot(rda_resultI_d5, scaling =1, 
     main = "Triplot dbRDA Region: - Set I: Horn Dissimilarity", cex.main = 1)

# D7a
rda_resultI_d7a <- dbrda(dist_d7a_rt_rda~., BioClim_SetI_final_scale)
rda_resultI_d7a
plot(rda_resultI_d7a, scaling =1, 
     main = "Triplot dbRDA Region: - Set I: Weighted UniFrac Distance", cex.main = 1)

# Set II
rda_resultII_d4a <- dbrda(community_RT_rda ~ ., BioClim_SetII_final_scale, distance = "jaccard")
rda_resultII_d4a
plot(rda_resultII_d4a, scaling =1, 
     main = "Triplot dbRDA Region: - Set II: Weighted Jaccard Dissimilarity", cex.main = 1)

# D5
rda_resultII_d5 <- dbrda(community_RT_rda ~ ., BioClim_SetII_final_scale, distance = "horn")
rda_resultII_d5

plot(rda_resultII_d5, scaling =1, 
     main = "Triplot dbRDA Region: - Set II: Horn Dissimilarity", cex.main = 1)

# D7a
rda_resultII_d7a <- dbrda(dist_d7a_rt_rda~., BioClim_SetII_final_scale)
rda_resultII_d7a
plot(rda_resultII_d7a, scaling =1, 
     main = "Triplot dbRDA Region: - Set II: Weighted UniFrac Distance", cex.main = 1)

#### 3-3-2. Permutation test: BioClim variables vs. D4a, D5, D7a
#D4a
set.seed(456) 
anova_rda_resultI_d4a <- anova(rda_resultI_d4a, by="terms", permutations = 999)
anova_rda_resultI_d4a

set.seed(456) 
anova_rda_resultII_d4a <- anova(rda_resultII_d4a, by="terms", permutations = 999)
anova_rda_resultII_d4a

#d5
set.seed(456) 
anova_rda_resultI_d5 <- anova(rda_resultI_d5, by="terms", permutations = 999)
anova_rda_resultI_d5

set.seed(456) 
anova_rda_resultII_d5 <- anova(rda_resultII_d5, by="terms", permutations = 999)
anova_rda_resultII_d5

#d7a
set.seed(456) 
anova_rda_resultI_d7a <- anova(rda_resultI_d7a, by="terms", permutations = 999)
anova_rda_resultI_d7a

set.seed(456) 
anova_rda_resultII_d7a <- anova(rda_resultII_d7a, by="terms", permutations = 999)
anova_rda_resultII_d7a

#### 3-3-3. Significance Test: BioClim variables vs. D4a, D5, D7a
#D4a: set I
reduced_modelI_d4a <- dbrda(community_RT_rda ~ 1, BioClim_SetI_final_scale, distance = "jaccard")
full_modelI_d4a <- dbrda(community_RT_rda ~ ., data = BioClim_SetI_final_scale, distance = "jaccard")
# Permutation Test
set.seed(456) 
anova_comparisonI_d4a <- anova(full_modelI_d4a, reduced_modelI_d4a, permutations = 999)
anova_comparisonI_d4a
# p-value
p_valueI_d4a <- anova_comparisonI_d4a$"Pr(>F)"[2]
p_valueI_d4a

#D4a: set II
reduced_modelII_d4a <- dbrda(community_RT_rda ~ 1, BioClim_SetII_final_scale, distance = "jaccard")
full_modelII_d4a <- dbrda(community_RT_rda ~ ., BioClim_SetII_final_scale, distance = "jaccard")
# Permutation Test
set.seed(456) 
anova_comparisonII_d4a <- anova(full_modelII_d4a, reduced_modelII_d4a, permutations = 999)
anova_comparisonII_d4a
# p-value
p_valueII_d4a <- anova_comparisonII_d4a$"Pr(>F)"[2]
p_valueII_d4a

#D5
reduced_modelI_d5 <- dbrda(community_RT_rda ~ 1, BioClim_SetI_final_scale, distance = "horn")
full_modelI_d5 <- dbrda(community_RT_rda ~ ., data = BioClim_SetI_final_scale, distance = "horn")
# Permutation Test
set.seed(456) 
anova_comparisonI_d5 <- anova(full_modelI_d5, reduced_modelI_d5, permutations = 999)
anova_comparisonI_d5
# p-value
p_valueI_d5 <- anova_comparisonI_d5$"Pr(>F)"[2]
p_valueI_d5

reduced_modelII_d5 <- dbrda(community_RT_rda ~ 1, BioClim_SetII_final_scale, distance = "horn")
full_modelII_d5 <- dbrda(community_RT_rda ~ ., data = BioClim_SetII_final_scale, distance = "horn")
# Permutation Test
set.seed(456) 
anova_comparisonII_d5 <- anova(full_modelII_d5, reduced_modelII_d5, permutations = 999)
anova_comparisonII_d5
# p-value
p_valueII_d5 <- anova_comparisonII_d5$"Pr(>F)"[2]
p_valueII_d5

#d7a
reduced_modelI_d7a <- dbrda(dist_d7a_rt_rda ~ 1, data = BioClim_SetI_final_scale)
full_modelI_d7a <- dbrda(dist_d7a_rt_rda ~ ., data = BioClim_SetI_final_scale)
# Permutation Test
set.seed(456) 
anova_comparisonI_d7a <- anova(full_modelI_d7a, reduced_modelI_d7a, permutations = 999)
anova_comparisonI_d7a
# p-value
p_valueI_d7a <- anova_comparisonI_d7a$"Pr(>F)"[2]
p_valueI_d7a

reduced_modelII_d7a <- dbrda(dist_d7a_rt_rda ~ 1, data = BioClim_SetII_final_scale)
full_modelII_d7a <- dbrda(dist_d7a_rt_rda ~ ., data = BioClim_SetII_final_scale)
# Permutation Test
set.seed(456) 
anova_comparisonII_d7a <- anova(full_modelII_d7a, reduced_modelII_d7a, permutations = 999)
anova_comparisonII_d7a
# p-value
p_valueII_d7a <- anova_comparisonII_d7a$"Pr(>F)"[2]
p_valueII_d7a

#### 3-2-4. Plot
# Convert results to dataframe
# A. vs. Geographic distance (GD)
# Regional Score
site_scoresI_geodis <- scores(rda_resultI_geodis, display = "sites", scaling = 1)
site_scoresI_geodis_df <- as.data.frame(site_scoresI_geodis)
site_scoresI_geodis_df$Region <- BioClim_SetI$Subplot_Nr
site_scoresI_geodis_df
# Species Score
species_scoresI_geodis <- scores(rda_resultI_geodis, display = "species")
species_scoresI_geodis_df <- as.data.frame(species_scoresI_geodis)
# Enviromental Variable Score
biplot_scoresI_geodis <- scores(rda_resultI_geodis, display = "bp")
biplot_scoresI_geodis_df <- as.data.frame(biplot_scoresI_geodis)
biplot_scoresI_geodis_df
# extract summary
summary_rdaI_geodis <- summary(rda_resultI_geodis)
# Calculate 'Constrained Proportion' 
total_inertiaI_geodis <- summary_rdaI_geodis$tot.chi
constrained_inertiaI_geodis <- summary_rdaI_geodis$constr.chi
constrained_proportionI_geodis <- constrained_inertiaI_geodis / total_inertiaI_geodis

# Regional Score
site_scoresII_geodis <- scores(rda_resultII_geodis, display = "sites", scaling = 1)
site_scoresII_geodis_df <- as.data.frame(site_scoresII_geodis)
site_scoresII_geodis_df$Region <- BioClim_SetII$Subplot_Nr
site_scoresII_geodis_df
# Species Score
species_scoresII_geodis <- scores(rda_resultII_geodis, display = "species")
species_scoresII_geodis_df <- as.data.frame(species_scoresII_geodis)
# Enviromental Variable Score
biplot_scoresII_geodis <- scores(rda_resultII_geodis, display = "bp")
biplot_scoresII_geodis_df <- as.data.frame(biplot_scoresII_geodis)
biplot_scoresII_geodis_df
# extract summary
summary_rdaII_geodis <- summary(rda_resultII_geodis)
# Calculate 'Constrained Proportion' 
total_inertiaII_geodis <- summary_rdaII_geodis$tot.chi
constrained_inertiaII_geodis <- summary_rdaII_geodis$constr.chi
constrained_proportionII_geodis <- constrained_inertiaII_geodis / total_inertiaII_geodis

# B. vs. weighted jaccaard dissimilarity (D4a)
# Regional Score
site_scoresI_d4a <- scores(rda_resultI_d4a, display = "sites", scaling = 1)
site_scoresI_d4a_df <- as.data.frame(site_scoresI_d4a)
site_scoresI_d4a_df$Region <- BioClim_SetI$Subplot_Nr
# Species Score
species_scoresI_d4a <- scores(rda_resultI_d4a, display = "species")
species_scoresI_d4a_df <- as.data.frame(species_scoresI_d4a)
# Enviromental Variable Score
biplot_scoresI_d4a <- scores(rda_resultI_d4a, display = "bp")
biplot_scoresI_d4a_df <- as.data.frame(biplot_scoresI_d4a)
# extract summary
summary_rdaI_d4a <- summary(rda_resultI_d4a)
# Calculate 'Constrained Proportion' 
total_inertiaI_d4a <- summary_rdaI_d4a$tot.chi
constrained_inertiaI_d4a <- summary_rdaI_d4a$constr.chi
constrained_proportionI_d4a <- constrained_inertiaI_d4a / total_inertiaI_d4a

# Regional Score
site_scoresII_d4a <- scores(rda_resultII_d4a, display = "sites", scaling = 1)
site_scoresII_d4a_df <- as.data.frame(site_scoresII_d4a)
site_scoresII_d4a_df$Region <- BioClim_SetII$Subplot_Nr
# Species Score
species_scoresII_d4a <- scores(rda_resultII_d4a, display = "species")
species_scoresII_d4a_df <- as.data.frame(species_scoresII_d4a)
# Enviromental Variable Score
biplot_scoresII_d4a <- scores(rda_resultII_d4a, display = "bp")
biplot_scoresII_d4a_df <- as.data.frame(biplot_scoresII_d4a)
# extract summary
summary_rdaII_d4a <- summary(rda_resultII_d4a)
# Calculate 'Constrained Proportion' 
total_inertiaII_d4a <- summary_rdaII_d4a$tot.chi
constrained_inertiaII_d4a <- summary_rdaII_d4a$constr.chi
constrained_proportionII_d4a <- constrained_inertiaII_d4a / total_inertiaII_d4a

# C. vs. Horn dissimilarity (D5)
# Regional Score
site_scoresI_d5 <- scores(rda_resultI_d5, display = "sites", scaling = 1)
site_scoresI_d5_df <- as.data.frame(site_scoresI_d5)
site_scoresI_d5_df$Region <- BioClim_SetI$Subplot_Nr
# Species Score
species_scoresI_d5 <- scores(rda_resultI_d5, display = "species")
species_scoresI_d5_df <- as.data.frame(species_scoresI_d5)
# Enviromental Variable Score
biplot_scoresI_d5 <- scores(rda_resultI_d5, display = "bp")
biplot_scoresI_d5_df <- as.data.frame(biplot_scoresI_d5)
# extract summary
summary_rdaI_d5 <- summary(rda_resultI_d5)
# Calculate 'Constrained Proportion' 
total_inertiaI_d5 <- summary_rdaI_d5$tot.chi
constrained_inertiaI_d5 <- summary_rdaI_d5$constr.chi
constrained_proportionI_d5 <- constrained_inertiaI_d5 / total_inertiaI_d5

# Regional Score
site_scoresII_d5 <- scores(rda_resultII_d5, display = "sites", scaling = 1)
site_scoresII_d5_df <- as.data.frame(site_scoresII_d5)
site_scoresII_d5_df$Region <- BioClim_SetII$Subplot_Nr
# Species Score
species_scoresII_d5 <- scores(rda_resultII_d5, display = "species")
species_scoresII_d5_df <- as.data.frame(species_scoresII_d5)
# Enviromental Variable Score
biplot_scoresII_d5 <- scores(rda_resultII_d5, display = "bp")
biplot_scoresII_d5_df <- as.data.frame(biplot_scoresII_d5)
# extract summary
summary_rdaII_d5 <- summary(rda_resultII_d5)
# Calculate 'Constrained Proportion' 
total_inertiaII_d5 <- summary_rdaII_d5$tot.chi
constrained_inertiaII_d5 <- summary_rdaII_d5$constr.chi
constrained_proportionII_d5 <- constrained_inertiaII_d5 / total_inertiaII_d5

# D. vs. weighted UniFrac dissimilarity (D7a)
# Regional Score
site_scoresI_d7a <- scores(rda_resultI_d7a, display = "sites", scaling = 1)
site_scoresI_d7a_df <- as.data.frame(site_scoresI_d7a)
site_scoresI_d7a_df$Region <- BioClim_SetI$Subplot_Nr
# Species Score
species_scoresI_d7a <- scores(rda_resultI_d7a, display = "species")
species_scoresI_d7a_df <- as.data.frame(species_scoresI_d7a)
# Enviromental Variable Score
biplot_scoresI_d7a <- scores(rda_resultI_d7a, display = "bp")
biplot_scoresI_d7a_df <- as.data.frame(biplot_scoresI_d7a)
# extract summary
summary_rdaI_d7a <- summary(rda_resultI_d7a)
# Calculate 'Constrained Proportion' 
total_inertiaI_d7a <- summary_rdaI_d7a$tot.chi
constrained_inertiaI_d7a <- summary_rdaI_d7a$constr.chi
constrained_proportionI_d7a <- constrained_inertiaI_d7a / total_inertiaI_d7a

# Regional Score
site_scoresII_d7a <- scores(rda_resultII_d7a, display = "sites", scaling = 1)
site_scoresII_d7a_df <- as.data.frame(site_scoresII_d7a)
site_scoresII_d7a_df$Region <- BioClim_SetII$Subplot_Nr
# Species Score
species_scoresII_d7a <- scores(rda_resultII_d7a, display = "species")
species_scoresII_d7a_df <- as.data.frame(species_scoresII_d7a)
# Enviromental Variable Score
biplot_scoresII_d7a <- scores(rda_resultII_d7a, display = "bp")
biplot_scoresII_d7a_df <- as.data.frame(biplot_scoresII_d7a)
# extract summary
summary_rdaII_d7a <- summary(rda_resultII_d7a)
# Calculate 'Constrained Proportion' 
total_inertiaII_d7a <- summary_rdaII_d7a$tot.chi
constrained_inertiaII_d7a <- summary_rdaII_d7a$constr.chi
constrained_proportionII_d7a <- constrained_inertiaII_d7a / total_inertiaII_d7a

#### 2-4-2. Plot
library(ggrepel)
region_order <- c("Khi", "Nor", "Kam1", "Kam2", "Tat", "B/V", "GAP", "Fr", "Cau", "Pyr", "SN", "Jp")
site_scoresI_geodis_df$Region <- factor(site_scoresI_geodis_df$Region, levels = region_order)

A_geodis <- ggplot() + theme_bw(base_size =15) +
  geom_point(data = site_scoresI_geodis_df, aes(x = dbRDA1, y = dbRDA2, color = Region , shape = Region), size = 10, bg = "steelblue") +
  scale_shape_manual(values = c(15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 17)) + # symbol(16●, 17▲, 18◼)
  scale_color_manual(values = c("#f8766dff", "#db8e00ff", "#aea200ff", "#aea200ff", "#64b200ff", "#00bd5cff", "#00c1a7ff", "#00badeff", "#00a6ffff", "#b385ffff", "#ef67ebff", "#ff63b6ff")) +
  geom_segment(data = biplot_scoresI_geodis_df, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), 
               arrow = arrow(length = unit(0.01, "npc")), lwd = 0.5) +
  geom_text_repel(data = biplot_scoresI_geodis_df, 
                  aes(x = dbRDA1, y = dbRDA2, label = rownames(biplot_scoresI_geodis_df)), color = "steelblue", size = 7) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"  
  ) +
  annotate("text", x = min(site_scoresI_geodis_df$dbRDA1), 
           y = max(biplot_scoresI_geodis_df$dbRDA2)*0.95, 
           label = "A", size = 12, hjust = 0) +
  annotate("text", x = Inf, y = Inf,
           label = paste("italic (p) ==", format(p_valueI_geodis, scientific = FALSE)),
           hjust = 1.1, vjust = 2.1, parse = TRUE, size = 6, color = "black") +
  labs(title = paste("Constrained Proportion: ", round(constrained_proportionI_geodis[1], 3)), 
       x = "Set I - Geographical Distance (GD): dbRDA1", 
       y = "Set I - Geographical Distance (GD): dbRDA2")   

A_geodis

site_scoresI_d4a_df$Region <- factor(site_scoresI_d4a_df$Region, levels = region_order)

B_d4a <- ggplot() + theme_bw(base_size =15) +
  geom_point(data = site_scoresI_d4a_df, aes(x = dbRDA1, y = dbRDA2, color = Region , shape = Region), size = 10, bg = "steelblue") +
  scale_shape_manual(values = c(15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 17)) + # symbol(16●, 17▲, 18◼)
  scale_color_manual(values = c("#f8766dff", "#db8e00ff", "#aea200ff", "#aea200ff", "#64b200ff", "#00bd5cff", "#00c1a7ff", "#00badeff", "#00a6ffff", "#b385ffff", "#ef67ebff", "#ff63b6ff")) +
  geom_segment(data = biplot_scoresI_d4a_df, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), 
               arrow = arrow(length = unit(0.01, "npc")), lwd = 0.5) +
  geom_text_repel(data = biplot_scoresI_d4a_df, 
                  aes(x = dbRDA1, y = dbRDA2, label = rownames(biplot_scoresI_d4a_df)), color = "steelblue", size = 7) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"  
  ) +
  annotate("text", x = min(biplot_scoresI_d4a_df$dbRDA1), 
           y = max(biplot_scoresI_d4a_df$dbRDA2)*0.95, 
           label = "B", size = 12, hjust = 0) +
  annotate("text", x = Inf, y = Inf,
           label = paste("italic (p) ==", format(p_valueI_d4a, scientific = FALSE)),
           hjust = 1.1, vjust = 2.1, parse = TRUE, size = 6, color = "black") +
  labs(title = paste("Constrained Proportion: ", round(constrained_proportionI_d4a[1], 3)), 
       x = "Set I - Weighted Jaccard Distance (D4): dbRDA1", 
       y = "Set I - Weighted Jaccard Distance (D4): dbRDA2")   

B_d4a

site_scoresI_d5_df$Region <- factor(site_scoresI_d5_df$Region, levels = region_order)

C_d5 <- ggplot() + theme_bw(base_size =15) +
  geom_point(data = site_scoresI_d5_df, aes(x = dbRDA1, y = dbRDA2, color = Region , shape = Region), size = 10, bg = "steelblue") +
  scale_shape_manual(values = c(15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 17)) + # symbol(16●, 17▲, 18◼)
  scale_color_manual(values = c("#f8766dff", "#db8e00ff", "#aea200ff", "#aea200ff", "#64b200ff", "#00bd5cff", "#00c1a7ff", "#00badeff", "#00a6ffff", "#b385ffff", "#ef67ebff", "#ff63b6ff")) +
  geom_segment(data = biplot_scoresI_d5_df, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), 
               arrow = arrow(length = unit(0.01, "npc")), lwd = 0.5) +
  geom_text_repel(data = biplot_scoresI_d5_df, 
                  aes(x = dbRDA1, y = dbRDA2, label = rownames(biplot_scoresI_d5_df)), color = "steelblue", size = 7) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"  
  ) +
  annotate("text", x = min(site_scoresI_d5_df$dbRDA1), 
           y = max(biplot_scoresI_d5_df$dbRDA2)*0.95, 
           label = "C", size = 12, hjust = 0) +
  annotate("text", x = Inf, y = Inf,
           label = paste("italic (p) ==", format(p_valueI_d5, scientific = FALSE)),
           hjust = 1.1, vjust = 2.1, parse = TRUE, size = 6, color = "black") +
  labs(title = paste("Constrained Proportion: ", round(constrained_proportionI_d5[1], 3)), 
       x = "Set I - Horn Dissimilarity (D5): dbRDA1", 
       y = "Set I - Horn Dissimilarity (D5): dbRDA2")   

C_d5

site_scoresI_d7a_df$Region <- factor(site_scoresI_d7a_df$Region, levels = region_order)

D_d7a <- ggplot() + theme_bw(base_size =15) +
  geom_point(data = site_scoresI_d7a_df, aes(x = dbRDA1, y = dbRDA2, color = Region , shape = Region), size = 10, bg = "steelblue") +
  scale_shape_manual(values = c(15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 17)) + # symbol(16●, 17▲, 18◼)
  scale_color_manual(values = c("#f8766dff", "#db8e00ff", "#aea200ff", "#aea200ff", "#64b200ff", "#00bd5cff", "#00c1a7ff", "#00badeff", "#00a6ffff", "#b385ffff", "#ef67ebff", "#ff63b6ff")) +
  geom_segment(data = biplot_scoresI_d7a_df, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), 
               arrow = arrow(length = unit(0.01, "npc")), lwd = 0.5) +
  geom_text_repel(data = biplot_scoresI_d7a_df, 
                  aes(x = dbRDA1, y = dbRDA2, label = rownames(biplot_scoresI_d7a_df)), color = "steelblue", size = 7) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"  
  ) +
  annotate("text", x = min(biplot_scoresI_d7a_df$dbRDA1), 
           y = max(site_scoresI_d7a_df$dbRDA2)*0.95, 
           label = "D", size = 12, hjust = 0) +
  annotate("text", x = Inf, y = Inf,
           label = paste("italic (p) ==", format(p_valueI_d7a, scientific = FALSE)),
           hjust = 1.1, vjust = 2.1, parse = TRUE, size = 6, color = "black") +
  labs(title = paste("Constrained Proportion: ", round(constrained_proportionI_d7a[1], 3)), 
       x = "Set I - Weighted UniFrac (D7a): dbRDA1", 
       y = "Set I - Weighted UniFrac (D7a): dbRDA2")   

D_d7a


site_scoresII_geodis_df$Region <- factor(site_scoresII_geodis_df$Region, levels = region_order)

E_geodis <- ggplot() + theme_bw(base_size =15) +
  geom_point(data = site_scoresII_geodis_df, aes(x = dbRDA1, y = dbRDA2, color = Region , shape = Region), size = 10, bg = "steelblue") +
  scale_shape_manual(values = c(15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 17)) + # symbol(16●, 17▲, 18◼)
  scale_color_manual(values = c("#f8766dff", "#db8e00ff", "#aea200ff", "#aea200ff", "#64b200ff", "#00bd5cff", "#00c1a7ff", "#00badeff", "#00a6ffff", "#b385ffff", "#ef67ebff", "#ff63b6ff")) +
  geom_segment(data = biplot_scoresII_geodis_df, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), 
               arrow = arrow(length = unit(0.01, "npc")), lwd = 0.5) +
  geom_text_repel(data = biplot_scoresII_geodis_df, 
                  aes(x = dbRDA1, y = dbRDA2, label = rownames(biplot_scoresII_geodis_df)), color = "steelblue", size = 7) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"  
  ) +
  annotate("text", x = min(site_scoresII_geodis_df$dbRDA1), 
           y = max(biplot_scoresII_geodis_df$dbRDA2)*0.95, 
           label = "E", size = 12, hjust = 0) +
  annotate("text", x = Inf, y = Inf,
           label = paste("italic (p) ==", format(p_valueII_geodis, scientific = FALSE)),
           hjust = 1.1, vjust = 2.1, parse = TRUE, size = 6, color = "black") +
  labs(title = paste("Constrained Proportion: ", round(constrained_proportionII_geodis[1], 3)), 
       x = "Set II - Geographical Distance (GD): dbRDA1", 
       y = "Set II - Geographical Distance (GD): dbRDA2")   

E_geodis

site_scoresII_d4a_df$Region <- factor(site_scoresII_d4a_df$Region, levels = region_order)

F_d4a <- ggplot() + theme_bw(base_size =15) +
  geom_point(data = site_scoresII_d4a_df, aes(x = dbRDA1, y = dbRDA2, color = Region , shape = Region), size = 10) +
  scale_shape_manual(values = c(15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 17)) + # symbol(16●, 17▲, 18◼)
  scale_color_manual(values = c("#f8766dff", "#db8e00ff", "#aea200ff", "#aea200ff", "#64b200ff", "#00bd5cff", "#00c1a7ff", "#00badeff", "#00a6ffff", "#b385ffff", "#ef67ebff", "#ff63b6ff")) +
  geom_segment(data = biplot_scoresII_d4a_df, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), 
               arrow = arrow(length = unit(0.01, "npc")), lwd = 0.5) +
  geom_text_repel(data = biplot_scoresII_d4a_df, 
                  aes(x = dbRDA1, y = dbRDA2, label = rownames(biplot_scoresII_d4a_df)), color = "steelblue", size = 7) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"  
  ) +
  annotate("text", x = min(biplot_scoresII_d4a_df$dbRDA1), 
           y = max(biplot_scoresII_d4a_df$dbRDA2)*0.95, 
           label = "F", size = 12, hjust = 0) +
  annotate("text", x = Inf, y = Inf,
           label = paste("italic (p) ==", format(p_valueII_d4a, scientific = FALSE)),
           hjust = 1.1, vjust = 2.1, parse = TRUE, size = 6, color = "black") +
  labs(title = paste("Constrained Proportion: ", round(constrained_proportionII_d4a[1], 3)), 
       x = "Set II - Weighted Jaccard Distance (D4): dbRDA1", 
       y = "Set II - Weighted Jaccard Distance (D4): dbRDA2")   

F_d4a

site_scoresII_d5_df$Region <- factor(site_scoresII_d5_df$Region, levels = region_order)

G_d5 <- ggplot() + theme_bw(base_size =15) +
  geom_point(data = site_scoresII_d5_df, aes(x = dbRDA1, y = dbRDA2, color = Region , shape = Region), size = 10, bg = "steelblue") +
  scale_shape_manual(values = c(15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 17)) + # symbol(16●, 17▲, 18◼)
  scale_color_manual(values = c("#f8766dff", "#db8e00ff", "#aea200ff", "#aea200ff", "#64b200ff", "#00bd5cff", "#00c1a7ff", "#00badeff", "#00a6ffff", "#b385ffff", "#ef67ebff", "#ff63b6ff")) +
  geom_segment(data = biplot_scoresI_d5_df, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), 
               arrow = arrow(length = unit(0.01, "npc")), lwd = 0.5) +
  geom_text_repel(data = biplot_scoresII_d5_df, 
                  aes(x = dbRDA1, y = dbRDA2, label = rownames(biplot_scoresII_d5_df)), color = "steelblue", size = 7) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"  
  ) +
  annotate("text", x = min(site_scoresII_d5_df$dbRDA1), 
           y = max(biplot_scoresII_d5_df$dbRDA2)*0.95, 
           label = "G", size = 12, hjust = 0) +
  annotate("text", x = Inf, y = Inf,
           label = paste("italic (p) ==", format(p_valueII_d5, scientific = FALSE)),
           hjust = 1.1, vjust = 2.1, parse = TRUE, size = 6, color = "black") +
  labs(title = paste("Constrained Proportion: ", round(constrained_proportionII_d5[1], 3)), 
       x = "Set II - Horn Dissimilarity (D5): dbRDA1", 
       y = "Set II - Horn Dissimilarity (D5): dbRDA2")   

G_d5

site_scoresII_d7a_df$Region <- factor(site_scoresII_d7a_df$Region, levels = region_order)

H_d7a <- ggplot() + theme_bw(base_size =15) +
  geom_point(data = site_scoresII_d7a_df, aes(x = dbRDA1, y = dbRDA2, color = Region , shape = Region), size = 10) +
  scale_shape_manual(values = c(15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 17)) + # symbol(16●, 17▲, 18◼)
  scale_color_manual(values = c("#f8766dff", "#db8e00ff", "#aea200ff", "#aea200ff", "#64b200ff", "#00bd5cff", "#00c1a7ff", "#00badeff", "#00a6ffff", "#b385ffff", "#ef67ebff", "#ff63b6ff")) +
  geom_segment(data = biplot_scoresII_d7a_df, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), 
               arrow = arrow(length = unit(0.01, "npc")), lwd = 0.5) +
  geom_text_repel(data = biplot_scoresII_d7a_df, 
                  aes(x = dbRDA1, y = dbRDA2, label = rownames(biplot_scoresII_d7a_df)), color = "steelblue", size = 7) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"  
  ) +
  annotate("text", x = min(biplot_scoresII_d7a_df$dbRDA1), 
           y = max(site_scoresI_d7a_df$dbRDA2)*0.95, 
           label = "H", size = 12, hjust = 0) +
  annotate("text", x = Inf, y = Inf,
           label = paste("italic (p) ==", format(p_valueII_d7a, scientific = FALSE)),
           hjust = 1.1, vjust = 2.1, parse = TRUE, size = 6, color = "black") +
  labs(title = paste("Constrained Proportion: ", round(constrained_proportionII_d7a[1], 3)), 
       x = "Set I - Weighted UniFrac (D7a): dbRDA1", 
       y = "Set I - Weighted UniFrac (D7a): dbRDA2")   

H_d7a

{A_geodis|B_d4a|C_d5|D_d7a}/{E_geodis|F_d4a|G_d5|H_d7a}
Fig8 <- {A_geodis|B_d4a|C_d5|D_d7a}/{E_geodis|F_d4a|G_d5|H_d7a}+
  plot_layout(guides = "collect") &
  theme(legend.position = "right")  # または "right"
Fig8
  
pdf("../_results/Fig8_dbRDA.pdf", width = 20, height = 14)
png("../_results/Fig8_dbRDA.png", width = 2000, height = 1400) #for paper
print(Fig8)
dev.off()


# Set I
# Geographic Distance
df_geodis_I <- as.data.frame(anova_rda_resultI_geodis)
df_geodis_I$Model <- "Set I - Geographic Distance"

# D4a
df_d4a_I <- as.data.frame(anova_rda_resultI_d4a)
df_d4a_I$Model <- "Set I - D4a: Weighted Jaccard Dissimilarity"

# D5
df_d5_I <- as.data.frame(anova_rda_resultI_d5)
df_d5_I$Model <- "Set I - D5: Horn Dissimilarity"

# D7a
df_d7a_I <- as.data.frame(anova_rda_resultI_d7a)
df_d7a_I$Model <- "Set I - D7a: Weighted UniFrac"


# Set II
df_geodis_II <- as.data.frame(anova_rda_resultII_geodis)
df_geodis_II$Model <- "Set II - Geographic Distance"

# D4a
df_d4a_II <- as.data.frame(anova_rda_resultII_d4a)
df_d4a_II$Model <- "Set II - D4a: Weighted Jaccard Dissimilarity"

# D5
df_d5_II <- as.data.frame(anova_rda_resultII_d5)
df_d5_II$Model <- "Set II - D5: Horn Dissimilarity"

# D7a
df_d7a_II <- as.data.frame(anova_rda_resultII_d7a)
df_d7a_II$Model <- "Set II - D7a: Weighted UniFrac"


# すべての結果を結合
table_S6 <- dplyr::bind_rows(
  df_geodis_I,
  df_d4a_I,
  df_d5_I,
  df_d7a_I,
  df_geodis_II,
  df_d4a_II,
  df_d5_II,
  df_d7a_II
)

# Extract and sort only the required columns
table_S6 <- table_S6 %>%
  dplyr::select(Model, everything())

write.csv(table_S6, "../_results/TabS6_dbRDA_permutation_test.csv", row.names = FALSE)



# Creating a data frame using a function (preserving variable order)
create_results_table <- function(anova_result, model_name) {
  df <- as.data.frame(anova_result)
  df$Variable <- rownames(df)
  df$Model <- model_name
  df$Original_Order <- 1:nrow(df)  
  rownames(df) <- NULL
  return(df)
}

results_list <- list(
  # Set I
  create_results_table(anova_rda_resultI_geodis, "Set I - Geographic Distance"),
  create_results_table(anova_rda_resultI_d4a, "Set I - D4a: Weighted Jaccard"),
  create_results_table(anova_rda_resultI_d5, "Set I - D5: Horn"),
  create_results_table(anova_rda_resultI_d7a, "Set I - D7a: Weighted UniFrac"),
  
  # Set II
  create_results_table(anova_rda_resultII_geodis, "Set II - Geographic Distance"),
  create_results_table(anova_rda_resultII_d4a, "Set II - D4a: Weighted Jaccard"),
  create_results_table(anova_rda_resultII_d5, "Set II - D5: Horn"),
  create_results_table(anova_rda_resultII_d7a, "Set II - D7a: Weighted UniFrac")
)

# combine
table_S6 <- dplyr::bind_rows(results_list)

# data cleaning and format
table_S6_clean <- table_S6 %>%
  dplyr::select(Model,  Variable, Df, SumOfSqs, Variance,`F` = F, `Pr(>F)` = `Pr(>F)`, Original_Order) %>%
  dplyr::mutate(
    # clean variables
    Variable = dplyr::case_when(
      Variable == "(Intercept)" ~ "Intercept",
      TRUE ~ Variable
    ),
    # format SumOfSqs
    `SumOfSqs` = round(`SumOfSqs`, 3),
    # format Variance
    `Variance` = round(`Variance`, 3),
    # format F-value
    `F` = round(`F`, 3),
    # format p-value
    `Pr(>F)` = ifelse(`Pr(>F)` < 0.001, sprintf("%.3e", `Pr(>F)`),
                      sprintf("%.3f", `Pr(>F)`))
  ) %>%
  # sort by model
  dplyr::arrange(
    factor(Model, levels = c(
      "Set I - Geographic Distance",
      "Set I - D4a: Weighted Jaccard", 
      "Set I - D5: Horn",
      "Set I - D7a: Weighted UniFrac",
      "Set II - Geographic Distance",
      "Set II - D4a: Weighted Jaccard",
      "Set II - D5: Horn", 
      "Set II - D7a: Weighted UniFrac"
    )),
    Original_Order
  ) %>%
  dplyr::select(-Original_Order)  

write.csv(table_S6_clean, 
          "../_results/TabS6_dbRDA_permutation_test_formatted.csv", 
          row.names = FALSE, 
          na = "")

library(kableExtra)

table_S6_kable <- table_S6_clean %>%
  kbl(align = c("l", "l", "c", "c", "c"),
      caption = "Table S6: dbRDA Permutation Test Results") %>%
  kable_paper(full_width = FALSE) %>%
  column_spec(1, bold = TRUE) %>%  # 
  collapse_rows(columns = 1, valign = "top") %>%  # Combine models with the same name
  # significant p-vallue without astarisk 
  footnote(general = "Permutation test results for distance-based redundancy analysis. P-values are reported as exact values without significance symbols.",
           general_title = "Note:")

table_S6_kable

## 3. Partial Mantel Test for Table 2
### 3-1. Calculate distance
# set order "dist_geo" as default
dis.hav_Coord_stand
geo_order <- rownames(dis.hav_Coord_stand)

dis.hav_Coord_stand_ordered <- as.dist(as.matrix(dis.hav_Coord_stand)[geo_order, geo_order])
dis.hav_Coord_stand_ordered

#### 3-1-1. Bioclimatic distance
# Set I
env_distance_scaleI <- dist(BioClim_SetI_final_scale)  
env_distance_scaleI_ordered <- as.dist(as.matrix(env_distance_scaleI)[geo_order, geo_order])
env_distance_scaleI_ordered
# SetII
env_distance_scaleII <- dist(BioClim_SetII_final_scale) 
env_distance_scaleII_ordered <- as.dist(as.matrix(env_distance_scaleII)[geo_order, geo_order])
env_distance_scaleII_ordered

#### 3-1-2. Bioclimatic distance
dis_D4a <- vegdist(community_RT_rda, method = "jaccard")
dis_D4a_ordered <- as.dist(as.matrix(dis_D4a)[geo_order, geo_order])
dis_D4a_ordered

dis_D5 <- vegdist(community_RT_rda, method = "horn")
dis_D5_ordered <- as.dist(as.matrix(dis_D5)[geo_order, geo_order])
dis_D5_ordered

dist_d7a_rt_rda
dis_D7a_ordered <- as.dist(as.matrix(dist_d7a_rt_rda)[geo_order, geo_order])
dis_D7a_ordered


### 3-2. Run partial Mantel test
# Set I: D4a
# Partial Mantel（control Geographical Distance, see only influence by bioclimatic factors）
set.seed(123)
partial_mantelI_d4a_bio <- mantel.partial(dis_D4a_ordered, env_distance_scaleI_ordered, dis.hav_Coord_stand_ordered)

# Partial Mantel（Control Bioclimatic Factors, see only influence by geographical distance）
set.seed(123)
partial_mantelI_d4a_geo <- mantel.partial(dis_D4a_ordered, dis.hav_Coord_stand_ordered, env_distance_scaleI_ordered)

partial_mantelI_d4a_bio
partial_mantelI_d4a_geo

# Set I: D5
# Partial Mantel（control Geographical Distance, see only influence by bioclimatic factors）
set.seed(123)
partial_mantelI_d5_bio <- mantel.partial(dis_D5_ordered, env_distance_scaleI_ordered, dis.hav_Coord_stand_ordered)

# Partial Mantel（Control Bioclimatic Factors, see only influence by geographical distance）
set.seed(123)
partial_mantelI_d5_geo <- mantel.partial(dis_D5_ordered, dis.hav_Coord_stand_ordered, env_distance_scaleI_ordered)

partial_mantelI_d5_bio
partial_mantelI_d5_geo

# Set I: D7a
# Partial Mantel（control Geographical Distance, see only influence by bioclimatic factors）
set.seed(123)
partial_mantelI_d7a_bio <- mantel.partial(dis_D7a_ordered, env_distance_scaleI_ordered, dis.hav_Coord_stand_ordered)

# Partial Mantel（Control Bioclimatic Factors, see only influence by geographical distance）
set.seed(123)
partial_mantelI_d7a_geo <- mantel.partial(dis_D7a_ordered, dis.hav_Coord_stand_ordered, env_distance_scaleI_ordered)

partial_mantelI_d7a_bio
partial_mantelI_d7a_geo

# Set II: D4a
# Partial Mantel（control Geographical Distance, see only influence by bioclimatic factors）
set.seed(123)
partial_mantelII_d4a_bio <- mantel.partial(dis_D4a_ordered, env_distance_scaleII_ordered, dis.hav_Coord_stand_ordered)

# Partial Mantel（Control Bioclimatic Factors, see only influence by geographical distance）
set.seed(123)
partial_mantelII_d4a_geo <- mantel.partial(dis_D4a_ordered, dis.hav_Coord_stand_ordered, env_distance_scaleII_ordered)

partial_mantelII_d4a_bio
partial_mantelII_d4a_geo

# Set II: d5
# Partial Mantel（control Geographical Distance, see only influence by bioclimatic factors）
set.seed(123)
partial_mantelII_d5_bio <- mantel.partial(dis_D5_ordered, env_distance_scaleII_ordered, dis.hav_Coord_stand_ordered)

# Partial Mantel（Control Bioclimatic Factors, see only influence by geographical distance）
set.seed(123)
partial_mantelII_d5_geo <- mantel.partial(dis_D5_ordered, dis.hav_Coord_stand_ordered, env_distance_scaleII_ordered)

partial_mantelII_d5_bio
partial_mantelII_d5_geo

# Set II: d7a
# Partial Mantel（control Geographical Distance, see only influence by bioclimatic factors）
set.seed(123)
partial_mantelII_d7a_bio <- mantel.partial(dis_D7a_ordered, env_distance_scaleII_ordered, dis.hav_Coord_stand_ordered)

# Partial Mantel（Control Bioclimatic Factors, see only influence by geographical distance）
set.seed(123)
partial_mantelII_d7a_geo <- mantel.partial(dis_D7a_ordered, dis.hav_Coord_stand_ordered, env_distance_scaleII_ordered)

partial_mantelII_d7a_bio
partial_mantelII_d7a_geo


# partial Mantel test detail
partial_mantel_results <- data.frame(
  Set = c(rep("Set I", 6), rep("Set II", 6)),
  Community_Distance = rep(rep(c("D4: Weighted Jaccard", "D5: Horn", "D7a: Weighted UniFrac"), each = 2), times = 2),
  Effect_Tested = rep(c("Bioclimatic factors", "Geographical distance"), times = 6),
  Partial_Mantel_r = c(
    partial_mantelI_d4a_bio$statistic, partial_mantelI_d4a_geo$statistic,
    partial_mantelI_d5_bio$statistic, partial_mantelI_d5_geo$statistic,
    partial_mantelI_d7a_bio$statistic, partial_mantelI_d7a_geo$statistic,
    partial_mantelII_d4a_bio$statistic, partial_mantelII_d4a_geo$statistic,
    partial_mantelII_d5_bio$statistic, partial_mantelII_d5_geo$statistic,
    partial_mantelII_d7a_bio$statistic, partial_mantelII_d7a_geo$statistic
  ),
  p_value = c(
    partial_mantelI_d4a_bio$signif, partial_mantelI_d4a_geo$signif,
    partial_mantelI_d5_bio$signif, partial_mantelI_d5_geo$signif,
    partial_mantelI_d7a_bio$signif, partial_mantelI_d7a_geo$signif,
    partial_mantelII_d4a_bio$signif, partial_mantelII_d4a_geo$signif,
    partial_mantelII_d5_bio$signif, partial_mantelII_d5_geo$signif,
    partial_mantelII_d7a_bio$signif, partial_mantelII_d7a_geo$signif
  )
)

# Format adjustment
partial_mantel_results <- partial_mantel_results %>%
  dplyr::mutate(
    Partial_Mantel_r = round(Partial_Mantel_r, 3),
    p_value = ifelse(p_value < 0.001, "<0.001", 
                     ifelse(p_value < 0.01, sprintf("%.3f", p_value),
                            sprintf("%.3f", p_value)))
  ) %>%
  dplyr::arrange(Set, Community_Distance, Effect_Tested)

print(partial_mantel_results)

write.csv(partial_mantel_results, "../_results/Tab2_partial_mantel.csv", row.names = FALSE)

# with kable
library(kableExtra)

partial_mantel_results %>%
  kbl(align = c("l", "l", "l", "c", "c"),
      caption = "Partial Mantel Test Results: Testing Relationships Between Community Dissimilarity and Environmental/Geographical Factors",
      col.names = c("Variable Set", "Community Distance", "Effect Tested", "Mantel r", "p-value")) %>%
  kable_paper(full_width = FALSE) %>%
  collapse_rows(columns = 1:2, valign = "top") %>%
  footnote(general = "Partial Mantel tests examining the relationship between community dissimilarity and either bioclimatic factors (controlling for geographical distance) or geographical distance (controlling for bioclimatic factors).")
