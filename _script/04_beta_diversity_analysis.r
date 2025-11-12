# 04. beta diversity analysis

source("utils_ecology_myxo.R")
load_ecology_packages()

## 1. UniFrac heatmap for Figure S9
### 1-1. Load and data preparation
unweighted_unifrac <- readRDS("../_results/beta_diversity/D7_unweighted.rds")
weighted_unifrac <- readRDS("../_results/beta_diversity/D7_weighted.rds")

unifrac_matrix_unw <- as.matrix(unweighted_unifrac)
unifrac_matrix <- as.matrix(weighted_unifrac)


### 1-2. Plot
# weighted UniFrac
library(pheatmap)
library(gridExtra)

ph1 <- pheatmap(unifrac_matrix_unw, display_numbers = TRUE, main = "A", cex= 1.5)
ph2 <- pheatmap(unifrac_matrix, display_numbers = TRUE, main = "B", cex= 1.5)

grid.arrange(ph1$gtable, ph2$gtable, ncol = 2)

# save as pdf
pdf("../_results/FigS9_UniFrac_heatmap.pdf", width = 20, height = 10)
png("../_results/FigS9_UniFrac_heatmap.png", width = 1600, height = 800)
grid.arrange(ph1$gtable, ph2$gtable, ncol = 2)
dev.off()


## 2. Phylogenetic Diversity for Fugure 6
### 2-1. Load and data preparation
library(picante)
community_RT<- readRDS("../_data/_processed_data/community_RT.rds")
rooted_tree <- readRDS("../_data/_processed_data/rooted_tree.rds")
tree_dist_matrix <- readRDS("../_data/_processed_data/tree_dist_matrix.rds")

### 2-2. Calculate Phylogenetic Diversity
#### 2-2-1. Calculate PD
pd_result <- pd(community_RT, rooted_tree)
pd_result

#### 2-2-2. Calculate standardized effect size of PD
set.seed(123)
sd_pd_result<- ses.pd(community_RT, rooted_tree)
sd_pd_result

#### 2-2-3. Create dataframe from result
sd_pd_results_df <- data.frame(
  Region = rownames(sd_pd_result),
  PD_Observed = sd_pd_result$pd.obs,
  PD_Expected = sd_pd_result$pd.rand.mean,
  PD_SD = sd_pd_result$pd.rand.sd,  # standard deviation
  SR = sd_pd_result$ntaxa,
  p_value = sd_pd_result$pd.obs.p
)
sd_pd_results_df

### 2-3. Calculate wieighted Mean Pairwise Distances (MPD)
#### 2-3-1. Calculate Standardized effect size of MPD
set.seed(123)
mpd_result_ab <- ses.mpd(community_RT, tree_dist_matrix, abundance.weighted = TRUE,　null.model = "taxa.labels", runs = 999)
mpd_result_ab

#### 2-3-2. Create dataframe from result
mpd_results_df <- data.frame(
  Region = rownames(mpd_result_ab),
  MPD_Observed = mpd_result_ab$mpd.obs,
  MPD_Expected = mpd_result_ab$mpd.rand.mean,
  MPD_SD = mpd_result_ab$mpd.rand.sd,  # standard deviation
  SR = mpd_result_ab$ntaxa,
  p_value = mpd_result_ab$mpd.obs.p
)
mpd_results_df

### 2-4. Calculate Mean Nearest Taxon Distances (MNTD)
#### 2-4-1. Calculate Standardized effect size of MNTD
set.seed(123)
mntd_result_ab <- ses.mntd(community_RT, tree_dist_matrix, abundance.weighted = TRUE,　
                           null.model = "taxa.labels", runs = 999)
mntd_result_ab

#### 2-4-2. Create dataframe from result
mntd_results_df <- data.frame(
  Region = rownames(mntd_result_ab),
  MNTD_Observed = mntd_result_ab$mntd.obs,
  MNTD_Expected = mntd_result_ab$mntd.rand.mean,
  MNTD_SD = mntd_result_ab$mntd.rand.sd,  # standard deviation
  SR = mntd_result_ab$ntaxa,
  p_value = mntd_result_ab$mntd.obs.p
)
mntd_results_df

### 2-5. Plot
library(ggrepel)

PD_plot <- ggplot(sd_pd_results_df, aes(x = SR, y = PD_Observed)) +
  geom_point(size= 3) +
  geom_point(aes(y = PD_Expected), color = "blue", size = 3) +  # expected PD as a dot
  geom_errorbar(aes(ymin = PD_Expected - PD_SD, ymax = PD_Expected + PD_SD), color = "gray30", width = 1) +  # エラーバーを追加
  geom_segment(data = subset(sd_pd_results_df, p_value <= 0.05), 
               aes(x = SR, xend = SR, y = PD_Observed, yend = PD_Expected), color = "red", linetype = 3, linewidth = 1.5) +  # p-valueが0.05以下の場合のみ赤線を引く
  geom_text_repel(aes(label = Region), size = 4) +
  annotate("text", x = 60, y = 12, label = "A", size = 8, hjust = 0) +
  labs(x = "Number of RTs", y = "Phylogenetic Diversity (PD)", 
       title = "") +
  theme_classic(base_size = 10) +
  theme(title = element_text(size=11),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
 #       plot.margin = unit(c(1, 0.5, 0, 0), "cm")
        )
PD_plot

MPD_plot <- ggplot(mpd_results_df, aes(x = SR, y = MPD_Observed)) +
  geom_point(size = 3) +
  geom_point(aes(y = MPD_Expected), color = "blue", size = 3) +  # expected MNTD as a dot
  geom_errorbar(aes(ymin = MPD_Expected - MPD_SD, ymax = MPD_Expected + MPD_SD), color = "gray30", width = 1) +  # add errror bar
  geom_segment(data = subset(mpd_results_df, p_value <= 0.05), 
               aes(xend = SR, yend = MPD_Expected), color = "red", linetype =3, linewidth = 1.5) +
  geom_text_repel(aes(label = Region), size = 4) +
  annotate("text", x = 60, y = 1, label = "B", size = 8, hjust = 0) +
  labs(x = "Number of RTs", y = "Weighted Mean Pairwise Distance (MPD)", 
       title = "") +
  theme_classic(base_size = 10) +
  theme(title = element_text(size=11),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
#        plot.margin = unit(c(1, 0.5, 0, 0.5), "cm")
        )
MPD_plot

MNTD_plot <- ggplot(mntd_results_df, aes(x = SR, y = MNTD_Observed)) +
  geom_point(size = 3) +
  geom_point(aes(y = MNTD_Expected), color = "blue", size = 3) +  # expected MNTD as a dot
  geom_errorbar(aes(ymin = MNTD_Expected - MNTD_SD, ymax = MNTD_Expected + MNTD_SD), color = "gray30", width = 1) +  # add error bar
  geom_segment(data = subset(mntd_results_df, p_value <= 0.05), 
               aes(xend = SR, yend = MNTD_Expected), color = "red", linetype = 3, linewidth = 1.5) +
  geom_text_repel(aes(label = Region), size = 4) +
  annotate("text", x = 60, y = 0.17, label = "C", size = 8, hjust = 0) +
  labs(x = "Number of RTs", y = "Weighted Mean Nearest Taxon Distances (MNTD)", 
       title = "") +
  theme_classic(base_size = 10) +
  theme(title = element_text(size=11),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
#        plot.margin = unit(c(1, 0.5, 0, 0.5), "cm")
        )
MNTD_plot

PD_plot|MPD_plot|MNTD_plot
Fig6 <- PD_plot|MPD_plot|MNTD_plot

# save as pdf/png
pdf("../_results/Fig6_PhyloDiv.pdf", width = 20, height = 7)
png("../_results/Fig6_PhyloDiv.png", width = 1200, height = 420) #for paper
print(Fig6)
dev.off()

## 3. Regression analysis for Figure 7, Figure S13 and Table S4-5

source("utils_ecology_myxo.R")
load_ecology_packages()

### 3-1. load beta diversity disDissimilarity data 
dist_d1_rt <- readRDS("../_results/beta_diversity/D1_RT.rds")
dist_d2_rt <- readRDS("../_results/beta_diversity/D2_RT.rds")
dist_d3_rt <- readRDS("../_results/beta_diversity/D3_RT.rds")
dist_d456_rt <- readRDS("../_results/beta_diversity/beta_diversity_RT.rds")
dist_d4a_rt <- dist_d456_rt$D4a_jaccard_weighted
dist_d4b_rt <- dist_d456_rt$D4b_jaccard_unweighted
dist_d5_rt <- dist_d456_rt$D5_horn
dist_d6_rt <- dist_d456_rt$D6_chao

dist_d1_rg <- readRDS("../_results/beta_diversity/D1_RG.rds")
dist_d2_rg <- readRDS("../_results/beta_diversity/D2_RG.rds")
dist_d3_rg <- readRDS("../_results/beta_diversity/D3_RG.rds")
dist_d456_rg <- readRDS("../_results/beta_diversity/beta_diversity_RG.rds")
dist_d4a_rg <- dist_d456_rg$D4a_jaccard_weighted
dist_d4b_rg <- dist_d456_rg$D4b_jaccard_unweighted
dist_d5_rg <- dist_d456_rg$D5_horn
dist_d6_rg <- dist_d456_rg$D6_chao

dist_d7a <- readRDS("../_results/beta_diversity/D7_unweighted.rds")
dist_d7b <- readRDS("../_results/beta_diversity/D7_weighted.rds")

### 3-2. load beta diversity (Dissimilarity) and calculate distance
#### 3-2-1. Load detailed geographic coordinate data and prepare data
GeoCoord <- readRDS("../_data/_processed_data/GeoCoord_avg_11regions_detail.rds")
GeoCoord_11 <- GeoCoord %>% 
  filter(Subplot_Nr == "main") %>% 
  dplyr::select(Code, Easting_avg, Northing_avg)

#### 3-2-2. Calculate geographic distance matrix (Haversine Method) for analysis (deatails with digits)
library(geosphere) # Haversine Method
dist_geo <- distm(GeoCoord_11[ ,-1], fun = distHaversine) # Haversine method 
# Set region name
dimnames(dist_geo) <- list(GeoCoord_11$Code, GeoCoord_11$Code)
dist_geo

dist_geo_km <- dist_geo/1000
dist_geo_km


### 3-3. Data preparation for mantel test
# set order "dist_geo" as default
geo_order <- rownames(dist_geo_km)

# function to unite matrix
reorder_distance_matrix <- function(dist_geo_km, reference_order) {
  if (inherits(dist_geo_km, "dist")) {
    # convert dist object to matrix
    dist_geo_km <- as.matrix(dist_geo_km)
  }
  # order 
  dist_geo_km[reference_order, reference_order]
}

# convert to matrix and order for mantel test
dist_d1_rt_ordered <- reorder_distance_matrix(dist_d1_rt, geo_order)
dist_d2_rt_ordered <- reorder_distance_matrix(dist_d2_rt, geo_order)
dist_d3_rt_ordered <- reorder_distance_matrix(dist_d3_rt, geo_order)
dist_d4a_rt_ordered <- reorder_distance_matrix(dist_d4a_rt, geo_order)
dist_d4b_rt_ordered <- reorder_distance_matrix(dist_d4b_rt, geo_order)
dist_d5_rt_ordered <- reorder_distance_matrix(dist_d5_rt, geo_order)
dist_d6_rt_ordered <- reorder_distance_matrix(dist_d6_rt, geo_order)

dist_d1_rg_ordered <- reorder_distance_matrix(dist_d1_rg, geo_order)
dist_d2_rg_ordered <- reorder_distance_matrix(dist_d2_rg, geo_order)
dist_d3_rg_ordered <- reorder_distance_matrix(dist_d3_rg, geo_order)
dist_d4a_rg_ordered <- reorder_distance_matrix(dist_d4a_rg, geo_order)
dist_d4b_rg_ordered <- reorder_distance_matrix(dist_d4b_rg, geo_order)
dist_d5_rg_ordered <- reorder_distance_matrix(dist_d5_rg, geo_order)
dist_d6_rg_ordered <- reorder_distance_matrix(dist_d6_rg, geo_order)

dist_d7a_ordered <- reorder_distance_matrix(dist_d7a, geo_order)
dist_d7b_ordered <- reorder_distance_matrix(dist_d7b, geo_order)


### 3-4. Analysis
#### 3-4-1. D1
set.seed(123)
d1_rt_mantel_result <- mantel(log10(dist_geo_km), 
                           dist_d1_rt_ordered, method = "pearson", 
                           permutations = 9999, na.rm = TRUE)
d1_rt_mantel_result

set.seed(123)
d1_rg_mantel_result <- mantel(log10(dist_geo_km), 
                              dist_d1_rg_ordered, method = "pearson", 
                              permutations = 9999, na.rm = TRUE)
d1_rg_mantel_result 

# Convert Result to "dataframe"
distance_df <- as.data.frame(as.table(dist_geo_km))
dist_d1_RT_df <- as.data.frame(as.table(dist_d1_rt_ordered))
dist_d1_RG_df <- as.data.frame(as.table(dist_d1_rg_ordered))

# Merge dataframe: RT
dist_d1_RT_merged_df <- merge(distance_df, dist_d1_RT_df, by=c("Var1", "Var2"))
colnames(dist_d1_RT_merged_df) <- c("Site1", "Site2", "Distance", "Dissimilarity")
str(dist_d1_RT_merged_df)

dist_d1_RG_merged_df <- merge(distance_df, dist_d1_RG_df, by=c("Var1", "Var2"))
colnames(dist_d1_RG_merged_df) <- c("Site1", "Site2", "Distance", "Dissimilarity")
str(dist_d1_RG_merged_df)

dist_d1_RT_cleaned_df <- filter(dist_d1_RT_merged_df, Distance != 0)
dist_d1_RG_cleaned_df <- filter(dist_d1_RG_merged_df, Distance != 0)

# Linear Regression Model
dist_d1_LM_RT <- lm(log10(Distance) ~ Dissimilarity, data = dist_d1_RT_cleaned_df)
lm(log10(Distance) ~ Dissimilarity, data = dist_d1_RT_cleaned_df)
summary(dist_d1_LM_RT)

dist_d1_LM_RG <- lm(log10(Distance) ~ Dissimilarity, data = dist_d1_RG_cleaned_df)
lm(log10(Distance) ~ Dissimilarity, data = dist_d1_RG_cleaned_df)

# Function to convert Dissimilarity D1 to Similarity to visualize for Figures
convert_dissim_to_sim_D1 <- function(dissim_matrix, original_beta_results) {
  # obtain "max_d1" from beta_results
  max_d1 <- max(original_beta_results$D1_similarity)
  
  # convert Dissimilarity to Similarity: similarity = max_d1 - dissimilarity
  sim_matrix <- max_d1 - as.matrix(dissim_matrix)
  
  return(sim_matrix)
}

community_RT <- readRDS("../_data/_processed_data/community_RT.rds")
community_RG <- readRDS("../_data/_processed_data/community_RG.rds")

beta_results_rt <- calculate_beta_diversity1(community_RT)
beta_results_rg <- calculate_beta_diversity1(community_RG)

sim_d1_rt_matrix <- convert_dissim_to_sim_D1(dist_d1_rt, beta_results_rt)
sim_d1_rt_ordered <- reorder_distance_matrix(sim_d1_rt_matrix, geo_order)
sim_d1_RT_df <- as.data.frame(as.table(sim_d1_rt_ordered))
sim_d1_RT_merged_df <- merge(distance_df, sim_d1_RT_df, by=c("Var1", "Var2"))
colnames(sim_d1_RT_merged_df) <- c("Site1", "Site2", "Distance", "Similarity")
sim_d1_RT_cleaned_df <- filter(sim_d1_RT_merged_df, Distance != 0)

sim_d1_rg_matrix <- convert_dissim_to_sim_D1(dist_d1_rg, beta_results_rg)
sim_d1_rg_ordered <- reorder_distance_matrix(sim_d1_rg_matrix, geo_order)
sim_d1_RG_df <- as.data.frame(as.table(sim_d1_rg_ordered))
sim_d1_RG_merged_df <- merge(distance_df, sim_d1_RG_df, by=c("Var1", "Var2"))
colnames(sim_d1_RG_merged_df) <- c("Site1", "Site2", "Distance", "Similarity")
sim_d1_RG_cleaned_df <- filter(sim_d1_RG_merged_df, Distance != 0)

library(scales)
D1_RT <- ggplot(sim_d1_RT_cleaned_df, aes(x=Distance, y= Similarity)) + theme_classic(base_size = 12) +
  geom_point() + stat_smooth(method="lm", se=T) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 2),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(title= paste("r =", round(d1_rt_mantel_result$statistic, 2), ", p =", d1_rt_mantel_result$signif), 
       x="Distance (km, log10 scale)", y="Shared RT count: D1") +
  annotate("text", x = 250, y = 50, label = "D", size = 7, hjust = 0) +
  annotate("text", x = 400, y = 16, 
           label = paste("'Adj.' ~ italic(R) ^ 2 ~ ' = ' ~ ", round(summary(dist_d1_LM_RT)$adj.r.squared, 2)),
           parse = TRUE, size = 4)
D1_RT

D1_RG <- ggplot(sim_d1_RG_cleaned_df, aes(x=Distance, y= Similarity)) + theme_classic(base_size = 12) +
  geom_point() + stat_smooth(method="lm", se=T) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 2),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(title= paste("r =", round(d1_rg_mantel_result$statistic, 2), ", p =", d1_rg_mantel_result$signif), 
       x="Distance (km, log10 scale)", y="Shared RG count: D1") +
  annotate("text", x = 250, y = 30, label = "A", size = 7, hjust = 0) +
  annotate("text", x = 400, y = 15, 
           label = paste("'Adj.' ~ italic(R) ^ 2 ~ ' = ' ~ ", round(summary(dist_d1_LM_RG)$adj.r.squared, 2)),
           parse = TRUE, size = 4)

#### 3-4-2. D2
set.seed(123)
d2_rt_mantel_result <- mantel(log10(dist_geo_km), 
                              dist_d2_rt_ordered, method = "pearson", 
                              permutations = 9999, na.rm = TRUE)
d2_rt_mantel_result

set.seed(123)
d2_rg_mantel_result <- mantel(log10(dist_geo_km), 
                              dist_d2_rg_ordered, method = "pearson", 
                              permutations = 9999, na.rm = TRUE)
d2_rg_mantel_result 

#Convert Result to "dataframe"
dist_d2_RT_df <- as.data.frame(as.table(dist_d2_rt_ordered))
dist_d2_RG_df <- as.data.frame(as.table(dist_d2_rg_ordered))

# Merge dataframe: RT
dist_d2_RT_merged_df <- merge(distance_df, dist_d2_RT_df, by=c("Var1", "Var2"))
colnames(dist_d2_RT_merged_df) <- c("Site1", "Site2", "Distance", "Dissimilarity")
str(dist_d2_RT_merged_df)

dist_d2_RG_merged_df <- merge(distance_df, dist_d2_RG_df, by=c("Var1", "Var2"))
colnames(dist_d2_RG_merged_df) <- c("Site1", "Site2", "Distance", "Dissimilarity")
str(dist_d2_RG_merged_df)

dist_d2_RT_cleaned_df <- filter(dist_d2_RT_merged_df, Distance != 0)
dist_d2_RG_cleaned_df <- filter(dist_d2_RG_merged_df, Distance != 0)

# Linear Regression Model
dist_d2_LM_RT <- lm(log10(Distance) ~ Dissimilarity, data = dist_d2_RT_cleaned_df)
lm(log10(Distance) ~ Dissimilarity, data = dist_d2_RT_cleaned_df)

dist_d2_LM_RG <- lm(log10(Distance) ~ Dissimilarity, data = dist_d2_RG_cleaned_df)
lm(log10(Distance) ~ Dissimilarity, data = dist_d2_RG_cleaned_df)

D2_RT <- ggplot(dist_d2_RT_cleaned_df, aes(x=Distance, y= 100 - Dissimilarity)) + theme_classic(base_size = 12) +
  geom_point() + stat_smooth(method="lm", se=T) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 2),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(title= paste("r =", round(d2_rt_mantel_result$statistic, 2), ", p =", d2_rt_mantel_result$signif), 
       x="Distance (km, log10 scale)", y="Prop. of shared RTs (%): D2") +
  annotate("text", x = 250, y = 55, label = "E", size = 7, hjust = 0) +
  annotate("text", x = 400, y = 18, 
           label = paste("'Adj.' ~ italic(R) ^ 2 ~ ' = ' ~ ", round(summary(dist_d2_LM_RT)$adj.r.squared, 2)),
           parse = TRUE, size = 4)

D2_RG <- ggplot(dist_d2_RG_cleaned_df, aes(x=Distance, y=100-Dissimilarity)) + theme_classic(base_size = 12) +
  geom_point() + stat_smooth(method="lm", se=T) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 2),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(title= paste("r =", round(d2_rg_mantel_result$statistic, 2), ", p =", d2_rg_mantel_result$signif), 
       x="Distance (km, log10 scale)", y="Prop. of shared RGs (%): D2") +
  annotate("text", x = 250, y = 80, label = "B", size = 7, hjust = 0) +
  annotate("text", x = 400, y = 40, 
           label = paste("'Adj.' ~ italic(R) ^ 2 ~ ' = ' ~ ", round(summary(dist_d2_LM_RG)$adj.r.squared, 2)),
           parse = TRUE, size = 4)

#### 3-4-3. D3
set.seed(123)
d3_rt_mantel_result <- mantel(log10(dist_geo_km), 
                              dist_d3_rt_ordered, method = "pearson", 
                              permutations = 9999, na.rm = TRUE)
d3_rt_mantel_result

set.seed(123)
d3_rg_mantel_result <- mantel(log10(dist_geo_km), 
                              dist_d3_rg_ordered, method = "pearson", 
                              permutations = 9999, na.rm = TRUE)
d3_rg_mantel_result 

#Convert Result to "dataframe"
dist_d3_RT_df <- as.data.frame(as.table(dist_d3_rt_ordered))
dist_d3_RG_df <- as.data.frame(as.table(dist_d3_rg_ordered))

# Merge dataframe: RT
dist_d3_RT_merged_df <- merge(distance_df, dist_d3_RT_df, by=c("Var1", "Var2"))
colnames(dist_d3_RT_merged_df) <- c("Site1", "Site2", "Distance", "Dissimilarity")
str(dist_d3_RT_merged_df)

dist_d3_RG_merged_df <- merge(distance_df, dist_d3_RG_df, by=c("Var1", "Var2"))
colnames(dist_d3_RG_merged_df) <- c("Site1", "Site2", "Distance", "Dissimilarity")
str(dist_d3_RG_merged_df)

dist_d3_RT_cleaned_df <- filter(dist_d3_RT_merged_df, Distance != 0)
dist_d3_RG_cleaned_df <- filter(dist_d3_RG_merged_df, Distance != 0)

# Linear Regression Model
dist_d3_LM_RT <- lm(log10(Distance) ~ Dissimilarity, data = dist_d3_RT_cleaned_df)
lm(log10(Distance) ~ Dissimilarity, data = dist_d3_RT_cleaned_df)

dist_d3_LM_RG <- lm(log10(Distance) ~ Dissimilarity, data = dist_d3_RG_cleaned_df)
lm(log10(Distance) ~ Dissimilarity, data = dist_d3_RG_cleaned_df)

D3_RT <- ggplot(dist_d3_RT_cleaned_df, aes(x=Distance, y=100 -Dissimilarity)) + theme_classic(base_size = 12) +
  geom_point() + stat_smooth(method="lm", se=T) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 2),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(title= paste("r =", round(d3_rt_mantel_result$statistic, 2), ", p =", d3_rt_mantel_result$signif), 
       x="Distance (km, log10 scale)", y="Prop. of shared RTs (%) in specimen count: D3") +
  annotate("text", x = 250, y = 90, label = "D", size = 7, hjust = 0) +
  annotate("text", x = 400, y = 30, 
           label = paste("'Adj.' ~ italic(R) ^ 2 ~ ' = ' ~ ", round(summary(dist_d3_LM_RT)$adj.r.squared, 2)),
           parse = TRUE, size = 4)

D3_RG <- ggplot(dist_d3_RG_cleaned_df, aes(x=Distance, y= 100 -Dissimilarity)) + theme_classic(base_size = 12) +
  geom_point() + stat_smooth(method="lm", se=T) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 2),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(title= paste("r =", round(d3_rg_mantel_result$statistic, 2), ", p =", d3_rg_mantel_result$signif), 
       x="Distance (km, log10 scale)", y="Prop. of shared RGs (%) in specimen count: D3") +
  annotate("text", x = 250, y = 100, label = "A", size = 7, hjust = 0) +
  annotate("text", x = 400, y = 60, 
           label = paste("'Adj.' ~ italic(R) ^ 2 ~ ' = ' ~ ", round(summary(dist_d3_LM_RG)$adj.r.squared, 2)),
           parse = TRUE, size = 4)

#### 3-4-4. D4a
set.seed(123)
d4a_rt_mantel_result <- mantel(log10(dist_geo_km), 
                               dist_d4a_rt_ordered, method = "pearson", 
                               permutations = 9999, na.rm = TRUE)
d4a_rt_mantel_result

set.seed(123)
d4a_rg_mantel_result <- mantel(log10(dist_geo_km), 
                               dist_d4a_rg_ordered, method = "pearson", 
                               permutations = 9999, na.rm = TRUE)
d4a_rg_mantel_result 

#Convert Result to "dataframe"
dist_d4a_RT_df <- as.data.frame(as.table(dist_d4a_rt_ordered))
dist_d4a_RG_df <- as.data.frame(as.table(dist_d4a_rg_ordered))

# Merge dataframe: RT
dist_d4a_RT_merged_df <- merge(distance_df, dist_d4a_RT_df, by=c("Var1", "Var2"))
colnames(dist_d4a_RT_merged_df) <- c("Site1", "Site2", "Distance", "Dissimilarity")
str(dist_d4a_RT_merged_df)

dist_d4a_RG_merged_df <- merge(distance_df, dist_d4a_RG_df, by=c("Var1", "Var2"))
colnames(dist_d4a_RG_merged_df) <- c("Site1", "Site2", "Distance", "Dissimilarity")
str(dist_d4a_RG_merged_df)

dist_d4a_RT_cleaned_df <- filter(dist_d4a_RT_merged_df, Distance != 0)
dist_d4a_RG_cleaned_df <- filter(dist_d4a_RG_merged_df, Distance != 0)

# Linear Regression Model
dist_d4a_LM_RT <- lm(log10(Distance) ~ Dissimilarity, data = dist_d4a_RT_cleaned_df)
lm(log10(Distance) ~ Dissimilarity, data = dist_d4a_RT_cleaned_df)

dist_d4a_LM_RG <- lm(log10(Distance) ~ Dissimilarity, data = dist_d4a_RG_cleaned_df)
lm(log10(Distance) ~ Dissimilarity, data = dist_d4a_RG_cleaned_df)

D4_RT <- ggplot(dist_d4a_RT_cleaned_df, aes(x=Distance, y=1-Dissimilarity)) + theme_classic(base_size = 12) +
  geom_point() + stat_smooth(method="lm", se=T) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 2),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(title= paste("r =", round(d4a_rt_mantel_result$statistic, 2), ", p =", d4a_rt_mantel_result$signif), 
       x="Distance (km, log10 scale)", y="Weighted Jaccard Similarity in RT: D4") +
  annotate("text", x = 250, y = 0.35, label = "E", size = 7, hjust = 0) +
  annotate("text", x = 400, y = 0.1, 
           label = paste("'Adj.' ~ italic(R) ^ 2 ~ ' = ' ~ ", round(summary(dist_d4a_LM_RT)$adj.r.squared, 2)),
           parse = TRUE, size = 4)

D4_RG <- ggplot(dist_d4a_RG_cleaned_df, aes(x=Distance, y=1-Dissimilarity)) + theme_classic(base_size = 12) +
  geom_point() + stat_smooth(method="lm", se=T) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 2),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(title= paste("r =", round(d4a_rg_mantel_result$statistic, 2), ", p =", d4a_rg_mantel_result$signif), 
       x="Distance (km, log10 scale)", y="Weighted Jaccard Similarity in RG: D4") +
  annotate("text", x = 250, y = 0.45, label = "B", size = 7, hjust = 0) +
  annotate("text", x = 400, y = 0.2, 
           label = paste("'Adj.' ~ italic(R) ^ 2 ~ ' = ' ~ ", round(summary(dist_d4a_LM_RG)$adj.r.squared, 2)),
           parse = TRUE, size = 4)

#### 3-4-5. D5
set.seed(123)
d5_rt_mantel_result <- mantel(log10(dist_geo_km), 
                              dist_d5_rt_ordered, method = "pearson", 
                              permutations = 9999, na.rm = TRUE)
d5_rt_mantel_result

set.seed(123)
d5_rg_mantel_result <- mantel(log10(dist_geo_km), 
                              dist_d5_rg_ordered, method = "pearson", 
                              permutations = 9999, na.rm = TRUE)
d5_rg_mantel_result 

#Convert Result to "dataframe"
dist_d5_RT_df <- as.data.frame(as.table(dist_d5_rt_ordered))
dist_d5_RG_df <- as.data.frame(as.table(dist_d5_rg_ordered))

# Merge dataframe: RT
dist_d5_RT_merged_df <- merge(distance_df, dist_d5_RT_df, by=c("Var1", "Var2"))
colnames(dist_d5_RT_merged_df) <- c("Site1", "Site2", "Distance", "Dissimilarity")
str(dist_d5_RT_merged_df)

dist_d5_RG_merged_df <- merge(distance_df, dist_d5_RG_df, by=c("Var1", "Var2"))
colnames(dist_d5_RG_merged_df) <- c("Site1", "Site2", "Distance", "Dissimilarity")
str(dist_d5_RG_merged_df)

dist_d5_RT_cleaned_df <- filter(dist_d5_RT_merged_df, Distance != 0)
dist_d5_RG_cleaned_df <- filter(dist_d5_RG_merged_df, Distance != 0)

# Linear Regression Model
dist_d5_LM_RT <- lm(log10(Distance) ~ Dissimilarity, data = dist_d5_RT_cleaned_df)
lm(log10(Distance) ~ Dissimilarity, data = dist_d5_RT_cleaned_df)

dist_d5_LM_RG <- lm(log10(Distance) ~ Dissimilarity, data = dist_d5_RG_cleaned_df)
lm(log10(Distance) ~ Dissimilarity, data = dist_d5_RG_cleaned_df)

D5_RT <- ggplot(dist_d5_RT_cleaned_df, aes(x=Distance, y=1-Dissimilarity)) + theme_classic(base_size = 12) +
  geom_point() + stat_smooth(method="lm", se=T) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 2),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(title= paste("r =", round(d5_rt_mantel_result$statistic, 2), ", p =", d5_rt_mantel_result$signif), 
       x="Distance (km, log10 scale)", y="Horn Similarity in RT: D5") +
  annotate("text", x = 250, y = 0.8, label = "F", size = 7, hjust = 0) +
  annotate("text", x = 400, y = 0.2, 
           label = paste("'Adj.' ~ italic(R) ^ 2 ~ ' = ' ~ ", round(summary(dist_d5_LM_RT)$adj.r.squared, 2)),
           parse = TRUE, size = 4)

D5_RG <- ggplot(dist_d5_RG_cleaned_df, aes(x=Distance, y=1-Dissimilarity)) + theme_classic(base_size = 12) +
  geom_point() + stat_smooth(method="lm", se=T) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 2),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(title= paste("r =", round(d5_rg_mantel_result$statistic, 2), ", p =", d5_rg_mantel_result$signif), 
       x="Distance (km, log10 scale)", y="Horn Similarity in RG: D5") +
  annotate("text", x = 250, y = 0.8, label = "C", size = 7, hjust = 0) +
  annotate("text", x = 400, y = 0.35, 
           label = paste("'Adj.' ~ italic(R) ^ 2 ~ ' = ' ~ ", round(summary(dist_d5_LM_RG)$adj.r.squared, 2)),
           parse = TRUE, size = 4)

#### 3-4-6. D6
set.seed(123)
d6_rt_mantel_result <- mantel(log10(dist_geo_km), 
                              dist_d6_rt_ordered, method = "pearson", 
                              permutations = 9999, na.rm = TRUE)
d6_rt_mantel_result

set.seed(123)
d6_rg_mantel_result <- mantel(log10(dist_geo_km), 
                              dist_d6_rg_ordered, method = "pearson", 
                              permutations = 9999, na.rm = TRUE)
d6_rg_mantel_result 

#Convert Result to "dataframe"
dist_d6_RT_df <- as.data.frame(as.table(dist_d6_rt_ordered))
dist_d6_RG_df <- as.data.frame(as.table(dist_d6_rg_ordered))

# Merge dataframe: RT
dist_d6_RT_merged_df <- merge(distance_df, dist_d6_RT_df, by=c("Var1", "Var2"))
colnames(dist_d6_RT_merged_df) <- c("Site1", "Site2", "Distance", "Dissimilarity")
str(dist_d6_RT_merged_df)

dist_d6_RG_merged_df <- merge(distance_df, dist_d6_RG_df, by=c("Var1", "Var2"))
colnames(dist_d6_RG_merged_df) <- c("Site1", "Site2", "Distance", "Dissimilarity")
str(dist_d6_RG_merged_df)

dist_d6_RT_cleaned_df <- filter(dist_d6_RT_merged_df, Distance != 0)
dist_d6_RG_cleaned_df <- filter(dist_d6_RG_merged_df, Distance != 0)

# Linear Regression Model
dist_d6_LM_RT <- lm(log10(Distance) ~ Dissimilarity, data = dist_d6_RT_cleaned_df)
lm(log10(Distance) ~ Dissimilarity, data = dist_d6_RT_cleaned_df)
summary(dist_d6_LM_RT)

dist_d6_LM_RG <- lm(log10(Distance) ~ Dissimilarity, data = dist_d6_RG_cleaned_df)
lm(log10(Distance) ~ Dissimilarity, data = dist_d6_RG_cleaned_df)
summary(dist_d6_LM_RG)

D6_RT <- ggplot(dist_d6_RT_cleaned_df, aes(x=Distance, y=1-Dissimilarity)) + theme_classic(base_size = 12) +
  geom_point() + stat_smooth(method="lm", se=T) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 2),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(title= paste("r =", round(d6_rt_mantel_result$statistic, 2), ", p =", d6_rt_mantel_result$signif), 
       x="Distance (km, log10 scale)", y="Chao Similarity in RT: D6") +
  annotate("text", x = 250, y = 1, label = "F", size = 7, hjust = 0) +
  annotate("text", x = 400, y = 0.25, 
           label = paste("'Adj.' ~ italic(R) ^ 2 ~ ' = ' ~ ", round(summary(dist_d6_LM_RT)$adj.r.squared, 2)),
           parse = TRUE, size = 4)

D6_RG <- ggplot(dist_d6_RG_cleaned_df, aes(x=Distance, y=1-Dissimilarity)) + theme_classic(base_size = 12) +
  geom_point() + stat_smooth(method="lm", se=T) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 2),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(title= paste("r =", round(d6_rg_mantel_result$statistic, 2), ", p =", d6_rg_mantel_result$signif), 
       x="Distance (km, log10 scale)", y="Chao Similarity in RG: D6") +
  annotate("text", x = 250, y = 1.1, label = "C", size = 7, hjust = 0) +
  annotate("text", x = 400, y = 0.45, 
           label = paste("'Adj.' ~ italic(R) ^ 2 ~ ' = ' ~ ", round(summary(dist_d6_LM_RG)$adj.r.squared, 2)),
           parse = TRUE, size = 4)

#Figure 7
{D3_RG|D4_RG|D5_RG}/{D3_RT|D4_RT|D5_RT}
Fig7 <- {D3_RG|D4_RG|D5_RG}/{D3_RT|D4_RT|D5_RT}
Fig7
 
pdf("../_results/Fig7_regression.pdf", width = 15, height = 10)
png("../_results/Fig7_regression.png", width = 1000, height = 700) #for paper
print(Fig7)
dev.off()

 #Figure S13
{D1_RG|D2_RG|D6_RG}/{D1_RT|D2_RT|D6_RT}
FigS13 <- {D1_RG|D2_RG|D6_RG}/{D1_RT|D2_RT|D6_RT}
FigS13

pdf("../_results/FigS13_regression.pdf", width = 15, height = 10)
png("../_results/FigS13_regression.png", width = 1000, height = 700) #for paper
print(FigS13)
dev.off()


### 3-5. Statistics table for Table S4
library(knitr)
library(kableExtra) # table layout 

#### 3-5-1. Create dataframe
Table_S4_RG_raw <- data.frame(
  Index = c(
    "D1: shared RG count", 
    "D2: proportion of shared RG", 
    "D3: proportion of shared RG weighted by specimen count", 
    "D4a: weighted Jaccard disDissimilarity RG", 
    "D5: Horn disDissimilarity RG", 
    "D6: Horn disDissimilarity RG"
  ),
  r = c(
    d1_rg_mantel_result$statistic, 
    d2_rg_mantel_result$statistic, 
    d3_rg_mantel_result$statistic, 
    d4a_rg_mantel_result$statistic, 
    d5_rg_mantel_result$statistic, 
    d6_rg_mantel_result$statistic
  ),
  p = c(
    d1_rg_mantel_result$signif, 
    d2_rg_mantel_result$signif, 
    d3_rg_mantel_result$signif, 
    d4a_rg_mantel_result$signif, 
    d5_rg_mantel_result$signif, 
    d6_rg_mantel_result$signif
  )
  ) %>%
  mutate(
    r = round(r, 4),
    p = round(p, 4)
  )

Table_S4_RT_raw <- data.frame(
  Index = c(
    "D1: shared RT count", 
    "D2: proportion of shared RT", 
    "D3: proportion of shared RT weighted by specimen count", 
    "D4a: weighted Jaccard disDissimilarity RT", 
    "D5: Horn disDissimilarity RT", 
    "D6: Horn disDissimilarity RT"
  ),
  r = c(
    d1_rt_mantel_result$statistic, 
    d2_rt_mantel_result$statistic, 
    d3_rt_mantel_result$statistic, 
    d4a_rt_mantel_result$statistic, 
    d5_rt_mantel_result$statistic, 
    d6_rt_mantel_result$statistic
  ),
  p = c(
    d1_rt_mantel_result$signif, 
    d2_rt_mantel_result$signif, 
    d3_rt_mantel_result$signif, 
    d4a_rt_mantel_result$signif, 
    d5_rt_mantel_result$signif, 
    d6_rt_mantel_result$signif
  )
) %>%
  mutate(
    r = round(r, 4),
    p = round(p, 4)
  )

#### 3-5-2. tidy with kable (for a check) and save
Table_S4_RG_raw %>%
  kable(
    format = "html", # output format（HTML、latex、markdown etc.）
    digits = 4,      
    col.names = c(
      "Index",
      "Mantel statistics (*r*)", # Italic in Markdown
      "Significance (*p*)"      # Italic in Markdown
    ),
    caption = "Table S4. Mantel test results for disDissimilarity indices (RG)."
  ) %>%
  kable_styling(
    bootstrap_options = "striped", 
    full_width = F
  )

# save as CSV
write.csv(
  Table_S4_RG_raw, 
  file = "../_results/TabS4_RG_Mantel_results.csv", 
  row.names = FALSE, 
  quote = FALSE      
)


Table_S4_RT_raw %>%
  kable(
    format = "html", # output format（HTML、latex、markdown etc.）
    digits = 4,      
    col.names = c(
      "Index",
      "Mantel statistics (*r*)", # Italic in Markdown
      "Significance (*p*)"      # Italic in Markdown
    ),
    caption = "Table S5. Mantel test results for disDissimilarity indices (RT)."
  ) %>%
  kable_styling(
    bootstrap_options = "striped", 
    full_width = F
  )

# save as CSV
write.csv(
  Table_S4_RT_raw, 
  file = "../_results/TabS4_RT_Mantel_results.csv", 
  row.names = FALSE, 
  quote = FALSE      
)

### 3-5. Statistics table for Table S5
library(purrr) # map_dfc

#### 3-5-1. Create model list
model_list_RG <- list(
  dist_d1_LM_RG, dist_d2_LM_RG, dist_d3_LM_RG, 
  dist_d4a_LM_RG, dist_d5_LM_RG, dist_d6_LM_RG
)

model_list_RT <- list(
  dist_d1_LM_RT, dist_d2_LM_RT, dist_d3_LM_RT, 
  dist_d4a_LM_RT, dist_d5_LM_RT, dist_d6_LM_RT
)

#### 3-5-2. calculate model summary and list
summary_list_RG <- purrr::map(model_list_RG, summary)
summary_list_RT <- purrr::map(model_list_RT, summary)

#### 3-5-3. Function to extract statistics from the results
extract_stats <- function(s) {
  # s: summary(lm_object)
  
  # calculate p-value of F-statistics
  p_model <- pf(s$fstatistic["value"], s$fstatistic["numdf"], s$fstatistic["dendf"], lower.tail = FALSE)
  
  # coefficient table
  coefs <- s$coefficients
  
  tibble(
    "Adj. R^2" = s$adj.r.squared,
    "F-statistic" = s$fstatistic["value"],
    "Df (model)" = s$fstatistic["numdf"],
    "Df (residual)" = s$fstatistic["dendf"],
    "p-value (Model)" = p_model,
    "Intercept (beta0)" = coefs["(Intercept)", "Estimate"],　# Predictor estimate value
    "Slope (beta1) Estimate" = coefs[2, "Estimate"],
    "p-value (beta1)" = coefs[2, "Pr(>|t|)"]         # Predictor p-value
  )
}

#### 3-5-4. Apply for all models and create dataframe
stats_data_RG <- purrr::map_df(summary_list_RG, extract_stats)
stats_data_RT <- purrr::map_df(summary_list_RT, extract_stats)

#### 3-5-5. Combine with "Index" and create the table
Table_S5_RG_raw <- data.frame(
  Index = c(
    "D1: shared RG count", "D2: proportion of shared RG", 
    "D3: proportion of shared RG weighted by specimen count", 
    "D4a: weighted Jaccard disDissimilarity RG", 
    "D5: Horn disDissimilarity RG", "D6: Chao disDissimilarity RG"
  ),
  # convert call$formula to characteristic data
  "Model formula" = purrr::map_chr(model_list_RG, ~ as.character(.x$call$formula)[3]), # Portion excluding the response variable and the tilde
  stringsAsFactors = FALSE
) %>%
  bind_cols(stats_data_RG) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
Table_S5_RG_raw$`F-statistic` <- round(Table_S5_RG_raw$`F-statistic`, 2)
Table_S5_RG_raw$`p-value (Model)` <- c("< 0.001", "< 0.001", "< 0.001", "< 0.001", "< 0.001", "< 0.001")
Table_S5_RG_raw$`p-value (beta1)` <- c("< 0.001", "< 0.001", "< 0.001", "< 0.001", "< 0.001", "< 0.001")

# save as CSV
write.csv(
  Table_S5_RG_raw, 
  file = "../_results/TabS5_RG_statistics.csv", 
  row.names = FALSE, 
  quote = FALSE      
)

Table_S5_RT_raw <- data.frame(
  Index = c(
    "D1: shared RT count", "D2: proportion of shared RT", 
    "D3: proportion of shared RT weighted by specimen count", 
    "D4a: weighted Jaccard disDissimilarity RT", 
    "D5: Horn disDissimilarity RT", "D6: Chao disDissimilarity RT"
  ),
  # convert call$formula to characteristic data
  "Model formula" = purrr::map_chr(model_list_RT, ~ as.character(.x$call$formula)[3]), # Portion excluding the response variable and the tilde
  stringsAsFactors = FALSE
) %>%
  bind_cols(stats_data_RT) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

Table_S5_RT_raw$`F-statistic` <- round(Table_S5_RT_raw$`F-statistic`, 2)
Table_S5_RT_raw$`p-value (Model)` <- c("< 0.001", "< 0.001", "< 0.001", "< 0.001", "< 0.001", "< 0.001")
Table_S5_RT_raw$`p-value (beta1)` <- c("< 0.001", "< 0.001", "< 0.001", "< 0.001", "< 0.001", "< 0.001")

# save as CSV
write.csv(
  Table_S5_RT_raw, 
  file = "../_results/TabS5_RT_statistics.csv", 
  row.names = FALSE, 
  quote = FALSE      
)
