# 02. Basic statistics 

source("utils_ecology_myxo.R")
load_ecology_packages()

## 1. Basic statistics for Table 1: Records & sequence records overview
### 1-1. Load RDSdata and create dataframe: all record, all record (seq attempt), sequenced data
DF_all <- readRDS("../_data/_processed_data/DF_biogeo_dark_rec.rds") # all records
table(DF_all$Code)

DF_niv <- readRDS("../_data/_processed_data/DF_biogeo_dark.rds") # data with sequencing attempt
table(DF_niv$Code)

DF_niv_seq <- readRDS("../_data/_processed_data/DF_biogeo_seq.rds") # barcoded specimens
table(DF_niv_seq$Code)

all <- DF_all %>%
  count(Code, name = "Records") # group by Code, count the raw
  
niv_rec <- DF_niv %>%
  count(Code, name = "rec") # group by Code, count the raw

niv_seq <- DF_niv_seq %>%
  count(Code, name = "seq") # group by Code, count the raw

### 1-2. Calculate and join the basic data
basic_all <- left_join(all, niv_rec) 
basic_all <- left_join(basic_all, niv_seq)

basic_all <- basic_all  %>%
  mutate(
    # create new row 'Percent' 
    percentage = (seq / rec) * 100,
    percentage = round(percentage, 1)
  )


### 1-3. Calculate sequence data: RT
#### 1-3-1. Load community data
community_RT <- readRDS("../_data/_processed_data/community_RT.rds")


#### 1-3-2. Calculate number of RTs
community_RT_df <- as.data.frame(community_RT)

RT_counts_by_region <- community_RT_df %>%
  (function(x) { x > 0 }) %>%
  rowSums() %>%
  data.frame(
    Code = names(.),
    RT_Count = as.numeric(.),
    row.names = NULL
  ) 
RT_counts_by_region  <-  RT_counts_by_region[,-1]

#### 1-3-3. Clauculate percent of total (473)
col_seq_data_percentage_RT <- RT_counts_by_region %>%
  mutate(
    # create new row 'Percentage_RT' 
    Percentage_RT = (RT_Count / 473) * 100,
    Percentage_RT = round(Percentage_RT, 1)
  )

basic_all <- left_join(basic_all, col_seq_data_percentage_RT) #join the data

#### 1-3-4. Clauculate singleton RT
singleton_counts_RT_df <- community_RT_df %>%
  rownames_to_column(var = "Code") %>%
  rowwise() %>%
  mutate(Singleton_RT = sum(c_across(-Code) == 1)) %>%
  dplyr::select(Code, Singleton_RT) %>%
  ungroup()

basic_all <- left_join(basic_all, singleton_counts_RT_df)


#### 1-3-5. Calculate "percent(region)" and join the basic data
basic_all <- basic_all  %>%
  mutate(
    # create new row 'Percent' 
    percent_RT = (Singleton_RT / RT_Count) * 100,
    percent_RT = round(percent_RT, 1)
  )


### 1-4. Calculate sequence data: RG
#### 1-4-1. Load community data
community_RG <- readRDS("../_data/_processed_data/community_RG.rds")

#### 1-4-2. Calculate number of RGs
community_RG_df <- as.data.frame(community_RG)

RG_counts_by_region <- community_RG_df %>%
  (function(x) { x > 0 }) %>%
  rowSums() %>%
  data.frame(
    Code = names(.),
    RG_Count = as.numeric(.),
    row.names = NULL
  ) 
RG_counts_by_region  <-  RG_counts_by_region[,-1]

#### 1-4-3. Clauculate percent of total (120)
col_seq_data_percentage_RG <- RG_counts_by_region %>%
  mutate(
    # create new row 'Percentage_RG' 
    Percentage_RG = (RG_Count / 120) * 100,
    Percentage_RG = round(Percentage_RG, 1)
  )

basic_all <- left_join(basic_all, col_seq_data_percentage_RG) #join the data

#### 1-4-4. Clauculate singleton RG
singleton_counts_RG_df <- community_RG_df %>%
  rownames_to_column(var = "Code") %>%
  rowwise() %>%
  mutate(Singleton_RG = sum(c_across(-Code) == 1)) %>%
  dplyr::select(Code, Singleton_RG) %>%
  ungroup()

basic_all <- left_join(basic_all, singleton_counts_RG_df)

#### 1-4-5. Calculate "percent(region)" and join the basic data
 basic_all <- basic_all  %>%
  mutate(
    # create new row 'Percent' 
    percent_RG = (Singleton_RG / RG_Count) * 100,
    percent_RG = round(percent_RG, 1)
  )


#### 1-5. Calculate sequence data: Clade
#### 1-5-1. Load community data
community_clade <- readRDS("../_data/_processed_data/community_clade.rds")

#### 1-5-2. Calculate number of Clades
community_clade_df <- as.data.frame(community_clade)

clade_counts_by_region <- community_clade_df %>%
  (function(x) { x > 0 }) %>%
  rowSums() %>%
  data.frame(
    Code = names(.),
    clade_Count = as.numeric(.),
    row.names = NULL
  ) 
clade_counts_by_region  <-  clade_counts_by_region[,-1]

#### 1-5-3. Clauculate percent of total (54)
col_seq_data_percentage_clade <- clade_counts_by_region %>%
  mutate(
    # create new row 'Percentage_clade' 
    Percentage_clade = (clade_Count / 54) * 100,
    Percentage_clade = round(Percentage_clade, 1)
  )

basic_all <- left_join(basic_all, col_seq_data_percentage_clade) #join the data

#### 1-5-4. Clauculate singleton clade
singleton_counts_clade_df <- community_clade_df %>%
  rownames_to_column(var = "Code") %>%
  rowwise() %>%
  mutate(Singleton_clade = sum(c_across(-Code) == 1)) %>%
  dplyr::select(Code, Singleton_clade) %>%
  ungroup()

basic_all <- left_join(basic_all, singleton_counts_clade_df)

#### 1-5-5. Calculate "percent(region)" and join the basic data
basic_all <- basic_all  %>%
  mutate(
    # create new row 'Percent' 
    percent_clade = (Singleton_clade / clade_Count) * 100,
    percent_clade = round(percent_clade, 1)
  )

#### 1-6. Order the region (South to North)
Code_order <- c("Jp", "SN", "Pyr", "Cau", "Fr", "GAP", "B/V", "Tat", "Kam", "Nor", "Khi")

basic_all <- basic_all %>%
  mutate(
    Code = factor(Code, levels = Code_order)
  ) %>%
  arrange(Code)

## 2. Basic statistics for Table 1 & Figure 5 & S9 : SAC & Alpha diversity
### 2-1. Calculate iNEXT and create list containing iNEXT results, relevant values, and basic SAC plot
region_codes <- unique(DF_niv_seq$Code)
all_region_stats <- list()

for(code in region_codes) {
  cat("Calculating for:", code, "\n")
  region_data <- DF_niv_seq %>% filter(Code == code)
  all_region_stats[[code]] <- calculate_region_stats(region_data, code)
}

all_region_stats$GAP$plot_data #for test


### 2-2. alpha diversity index for Table 1
# Function: extract Shannon Diversity and convert to dataframe
extract_shannon_diversity_safe <- function(all_region_stats) {
  shannon_data <- list()
  
  for(code in names(all_region_stats)) {
    stats <- all_region_stats[[code]]
    
    # filter with both Assemblage and Diversity
    rg_shannon <- stats$inext_result$AsyEst %>%
      filter(Assemblage == "RG", Diversity == "Shannon diversity")
    
    rt_shannon <- stats$inext_result$AsyEst %>%
      filter(Assemblage == "RT", Diversity == "Shannon diversity")
    
    region_shannon <- data.frame(
      Code = code,
      Shannon_RT = round(rt_shannon$Observed, 1),
      Shannon_RG = round(rg_shannon$Observed, 1)
    )
    
    shannon_data[[code]] <- region_shannon
  }
  
  bind_rows(shannon_data)
}

# Function: extract Simpson Diversity and convert to dataframe
extract_simpson_diversity_safe <- function(all_region_stats) {
  simpson_data <- list()
  
  for(code in names(all_region_stats)) {
    stats <- all_region_stats[[code]]
    
    # filter with both Assemblage and Diversity
    rg_simpson <- stats$inext_result$AsyEst %>%
      filter(Assemblage == "RG", Diversity == "Simpson diversity")
    
    rt_simpson <- stats$inext_result$AsyEst %>%
      filter(Assemblage == "RT", Diversity == "Simpson diversity")
    
    region_simpson <- data.frame(
      Code = code,
      Simpson_RT = round(rt_simpson$Observed, 1),
      Simpson_RG = round(rg_simpson$Observed, 1)
    )
    
    simpson_data[[code]] <- region_simpson
  }
  
  bind_rows(simpson_data)
}

shannon_table <- extract_shannon_diversity_safe(all_region_stats)
print(shannon_table)

simpson_table <- extract_simpson_diversity_safe(all_region_stats)
print(simpson_table)

basic_all <- left_join(basic_all, shannon_table)


## 3. Calculate rarefied Nr. for Table 1
min_sample_size_rt <- min(rowSums(community_RT))
set.seed(123)
rarefied_data_rt <- round(rarefy(community_RT, sample = min_sample_size_rt), 1)
rarefied_data_rt

rarefied_rt_df <- data.frame(
  Code = names(rarefied_data_rt),
  Rarefied_RT = as.numeric(rarefied_data_rt)
)

basic_all <- left_join(basic_all, rarefied_rt_df)


min_sample_size_rg <- min(rowSums(community_RG))
set.seed(123)
rarefied_data_rg <- round(rarefy(community_RG, sample = min_sample_size_rg), 1)
rarefied_data_rg

rarefied_rg_df <- data.frame(
  Code = names(rarefied_data_rg),
  Rarefied_RG = as.numeric(rarefied_data_rg)
)

basic_all <- left_join(basic_all, rarefied_rg_df)

min_sample_size_clade <- min(rowSums(community_clade))
set.seed(123)
rarefied_data_clade <- round(rarefy(community_clade, sample = min_sample_size_clade), 1)
rarefied_data_clade

rarefied_clade_df <- data.frame(
  Code = names(rarefied_data_clade),
  Rarefied_clade = as.numeric(rarefied_data_clade)
)

basic_all <- left_join(basic_all, rarefied_clade_df)

## 4. Expected diversity RT/RG & percent for Table 1
# Function: extract expected diversity and convert to dataframe
extract_expect_cov_diversity <- function(all_region_stats) {
  expect_data <- list()
  
  for(code in names(all_region_stats)) {
    stats <- all_region_stats[[code]]
    # Extract from stats
    region_expect <- data.frame(
      Code = code,
      Expect_RT = round(stats$exp_RT, 0),
      percent_found_RT = round(stats$cov_RT, 1),
      Expect_RG = round(stats$exp_RG, 0),
      percent_found_RG = round(stats$cov_RG, 1)
    )
    
    expect_data[[code]] <- region_expect
  }
  
  bind_rows(expect_data)
}

expect_table <- extract_expect_cov_diversity(all_region_stats)
print(expect_table)

basic_all <- left_join(basic_all, expect_table)

print(basic_all)

long_data <- basic_all %>%
  pivot_longer(
    cols = -Code, 
    names_to = "Variable", 
    values_to = "Value"    
  )

transposed_df <- long_data %>%
  pivot_wider(
    names_from = Code,    
    values_from = Value,  
  ) %>%
  column_to_rownames(var = "Variable") 

print(transposed_df)

write.csv(
  x = transposed_df,
  file = "../_results/Tab1_basic_statistics.csv",
  row.names = TRUE # set variables as column name
)

saveRDS(basic_all, "../_data/_processed_data/basic_statistics.rds")

## 5. Calculate sequence coverage plot for Figure 1B
### 5-1.  join two data by code 
col_seq_data <- full_join(  
  niv_rec, 
  niv_seq, 
  by = "Code"
)

col_seq_data <- col_seq_data %>%
  mutate(fail = rec - seq,
         coverage = seq / rec * 100)  # Coverage in percentage
col_seq_data

### 5-2. Convert Dataframe to Long Data
#library(tidyr)
long_data <- col_seq_data %>%
  pivot_longer(cols = c(seq, fail), names_to = "Sequencing", values_to = "Count") %>%
  mutate(Sequencing = factor(Sequencing, levels = c( "fail", "seq"), labels = c( "Failed", "Successful")))
long_data

### 5-3. Create Dataframe for Label
text_data <- col_seq_data %>%
  mutate(label = paste(Code, "\nTotal:", rec, "\nCov.:", scales::percent(seq / rec, accuracy = 0.1))) %>%
  dplyr::select(Code, label)

### 5-4. Set Region Order
custom_order <- c("Jp", "SN", "Pyr", "Cau", "Fr", "GAP", "B/V", "Tat", "Kam", "Nor", "Khi" )

### 5-5. Plot
p <- ggplot(long_data, aes(x = factor(Code, levels = custom_order), y = Count, fill = Sequencing)) +
  geom_bar(stat = "identity", position = "stack", colour = "gray30") +
  scale_fill_manual(values = c("Failed" = "white", "Successful" = "gray30")) +
  labs(title = "B",
       x = "Region",
       y = "Record\n(Dark-Spored Nivi. Myxomycetes) ") +
  theme_minimal(base_size = 20) +
  theme(title = element_text(size = 25),  
        axis.title = element_text(size = 15), 
        axis.text.x = element_text(size = 13),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13)) +
  scale_x_discrete(labels = function(x) {
    # label below the bars
    labels <- text_data$label
    names(labels) <- text_data$Code
    return(labels[x])
  })

p

ggsave(
  filename = "../_results/Fig1B_seq_covererage.pdf", 
  plot = p,                    
  width = 20,                          
  height = 7,                         
  dpi = 300                           
)


## 6. Geographic distance for Fig 1C and Tab S1
### 6-1. Load and create data for distance calculation
source("utils_ecology_myxo.R")
load_ecology_packages()

GeoCoord <- readRDS("../_data/_processed_data/GeoCoord_avg_11regions_detail.rds")

GeoCoord_plot <- GeoCoord %>% 
  filter(Subplot_Nr == "main") %>% 
  dplyr::select(Code, Easting_avg, Northing_avg)

### 6-2.Calculate distance matrix (Haversine Method) for Tab S1
library(geosphere) # Haversine Method
dis.hav_Coord <- distm(GeoCoord_plot[ ,-1], fun = distHaversine) # haversine method 

#### 6-2-1. Set region name
dimnames(dis.hav_Coord) <- list(GeoCoord_plot$Code, GeoCoord_plot$Code)
dis.hav_Coord

#### 6-2-2. Convert to km and order for Tab 1 and save
dis.hav_Coord_km <- round(dis.hav_Coord / 1000, 0)
dis.hav_Coord_km

current_order <- dimnames(dis.hav_Coord_km)[[1]]
new_order <- rev(current_order)
Table_S1_raw <- dis.hav_Coord_km[new_order, new_order]
Table_S1_raw

Table_S1_raw[upper.tri(Table_S1_raw, diag = TRUE)] <- NA
diag(Table_S1_raw) <- 0 
Table_S1_raw

write.csv(
  x = Table_S1_raw,
  file = "../_results/TabS1_pairwise_distance.csv",
  row.names = TRUE # set variables as column name
)

saveRDS(dis.hav_Coord_km, "../_data/_processed_data/dis.hav_Coord_km.rds")

### 6-3. Calculate Dstance from GAP
#### 6-3-1. Set the Center: GAP
c_Northing <- GeoCoord_plot[[6, 3]] #extract as number using double [[]]
c_Easting <- GeoCoord_plot[[6, 2]]
c_Northing
c_Easting

#### 6-3-2. Calculate Distance from GAP in Vector
GeoCoord_plot_dist <- GeoCoord_plot %>% 
  mutate(vector_Northing = Northing_avg - c_Northing) %>% 
  mutate(vector_Easting = Easting_avg - c_Easting)

#### 6-3-3. Calculate Distance from GAP in Metre: Haversine Method
GeoCoord_plot_dist$distance <- distHaversine(
  matrix(c(c_Easting, c_Northing), ncol=2),
  matrix(c(GeoCoord_plot_dist$Easting_avg, GeoCoord_plot_dist$Northing_avg), ncol=2)
)

#### 6-3-4. Prepare Plot
dist_p <- ggplot() +
  geom_point(data = GeoCoord_plot_dist, aes(x = Easting_avg, y = Northing_avg), size = 3) +
  geom_segment(data = GeoCoord_plot_dist, 
               aes(x = c_Easting, y = c_Northing, 
                   xend = Easting_avg, yend = Northing_avg),
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_text(data = GeoCoord_plot_dist, 
            aes(x = Easting_avg, y = Northing_avg, label = paste0(Code, "\n (", round(distance / 1000), " km", ")")), 
            vjust = 0, hjust = 0.5, size = 3) +
  labs(title = "",
       x = "Easting", y = "Northing") +
  annotate("text", x = -5, y = 70, label = "C", col = "black", size = 8) +
  theme_classic(base_size =10)
dist_p

ggsave(
  filename = "../_results/Fig1C_biogeogr_dist.pdf", 
  plot = dist_p,                    
  width = 7,                          
  height = 4,                         
  dpi = 300                           
)

## 7. Elevational distribution for Figure 1D
### 7-1. Data preparation
DF_niv_seq_elev <- dplyr::select(DF_niv_seq, final_colno, Code, Clade, Elev)

### 7-2. Create GEN name from the first three letters of Clade and add as new column
DF_niv_seq_elev$GEN <- substr(DF_niv_seq_elev$Clade, 1, 3)
DF_niv_seq_elev <- DF_niv_seq_elev[, -3]
summary(DF_niv_seq_elev$Elev)
unique(DF_niv_seq_elev$GEN)


DF_niv_seq_elev <- DF_niv_seq_elev %>%
  mutate(GEN = ifelse(GEN %in% "COM", "Others", GEN)) #Group “COM” as “Others”
unique(DF_niv_seq_elev$GEN)


### 7-3. Plot
library(ggbeeswarm)
custom_order <- c("Jp", "SN", "Pyr", "Cau", "Fr", "GAP", "B/V", "Tat", "Kam", "Nor", "Khi")
GEN_levels <- c("MER", "LAM", "DID", "DDY", "PLS", "PHY", "BAD", "Others")
DF_niv_seq_elev$GEN <- factor(DF_niv_seq_elev$GEN, levels = GEN_levels)


elev_p <- ggplot(DF_niv_seq_elev) + 
  theme_bw(base_size = 10) +
  geom_quasirandom(aes(y = Elev, x = factor(Code, levels = custom_order), color = GEN), size = 0.5, alpha = 0.9) + # ggbeeswarmのgeom_quasirandomを追加
  geom_violin(aes(y = Elev, x = Code), scale = "count", alpha = 0) + 
  scale_color_manual(values = c("BAD" = "red4", "PHY" = "red4", 
                                "PLS" = "olivedrab", "DDY" = "tan1", 
                                "DID" = "snow4", "LAM" = "royalblue3", 
                                "MER" = "pink4", "Others" = "black")) +  
  labs(title = "", y = "Elevation (m a.s.l.)", x = "Region", color = "Genus") + # y-axe label "Elevation (m a.s.l.)" 
  annotate("text", x = "Jp", y = 2900, label = "D", col = "black", size = 8)  # annotation

elev_p

ggsave(
  filename = "../_results/Fig1D_elev_violin.pdf", 
  plot = elev_p,                    
  width = 4.5,                          
  height = 4,                         
  dpi = 300                           
)

### 7-4. Text in Results for Figure 1D (Correlation between latitude and elevation)
#### 7-4-1. Data preparation
region_lat <- DF_niv_seq %>%
  dplyr::group_by(Code) %>%
  dplyr::summarise(lat = mean(Northing, na.rm = TRUE))  

# Defining the “peak elevation” for each region at the 95th percentile
region_peak <- DF_niv_seq_elev %>%
  dplyr::group_by(Code) %>%
  dplyr::summarise(
    peak_elev_95 = quantile(Elev, probs = 0.95, na.rm = TRUE),
    peak_elev_max = max(Elev, na.rm = TRUE),
    n = dplyr::n()
  )

# Combining latitude and peak elevation
region_summary <- dplyr::left_join(region_peak, region_lat, by = "Code")

#### 7-4-2. Model fitting: linear regression model & Spearman's correlation test
lm_fit <- lm(peak_elev_95 ~ lat, data = region_summary)
lm_fit_summary <- summary(lm_fit)     
lm_fit_summary

# calculate p-value of F-statistics
p_model <- pf(lm_fit_summary$fstatistic["value"], lm_fit_summary$fstatistic["numdf"], lm_fit_summary$fstatistic["dendf"], lower.tail = FALSE)

cat("Correlation between latitude and elevation")
cat("linear regression: β =", round(lm_fit_summary$coefficients["lat", "Estimate"], 2), "m/°")
cat("Adj. R2 = ", round(lm_fit_summary$adj.r.squared, 3))
cat("P-value = ", p_model)
cat("Df (model)" = lm_fit_summary$fstatistic["numdf"])
cat("Df (residual)" = lm_fit_summary$fstatistic["dendf"])
cat("F-statistic = ", round(lm_fit_summary$fstatistic["value"], 3))

# Spearman's correlation (non-parametric)
cor_test <- cor.test(region_summary$peak_elev_95, region_summary$lat, method = "spearman")
cor_test

#### 7-4-3. Plot
# Scatter plot on the ligression model
library(ggplot2)
scatter_p <- ggplot(region_summary, aes(x = lat, y = peak_elev_95)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(x = "Latitude (°)", y = "Peak elevation (m a.s.l.)",
       title = paste0("Spearman ρ = ", round(cor_test$estimate, 3), ", p = ", round(cor_test$p.value, 3),
                      "; \nLinear regression: Adj. R^2 = ", round(lm_fit_summary$adj.r.squared, 3), 
                      ", Slope (β1) = ", round(lm_fit_summary$coefficients["lat", "Estimate"], 2), "m/°")) 
scatter_p

# Peak elevation on elev_p
elev_p_with_peaks <- elev_p +
  geom_point(data = region_summary, aes(x = Code, y = peak_elev_95),
             color = "black", size = 2, shape = 21, fill = "yellow") +
  geom_text(data = region_summary, aes(x = Code, y = peak_elev_95, label = round(peak_elev_95)),
            vjust = -0.7, size = 2.5) +
  labs(title = paste0("Spearman ρ = ", round(cor_test$estimate, 3), ", p = ", round(cor_test$p.value, 3),
                      "; \nLinear regression: Adj. R^2 = ", round(lm_fit_summary$adj.r.squared, 3), 
                      ", Slope (β1) = ", round(lm_fit_summary$coefficients["lat", "Estimate"], 2), "m/°"))

elev_p_with_peaks


## 8. nMDS calculation for Fig 3DE Table S2&3
### 8-1. Load data and preparation
source("utils_ecology_myxo.R")
load_ecology_packages()

community_RT <- readRDS("../_data/_processed_data/community_RT.rds")
community_RG <- readRDS("../_data/_processed_data/community_RG.rds")

### 8-2. Calculate Rarefaction: select the stable solution after multiple iterations
min_sample_size_RT <- min(rowSums(community_RT))
rarefied_data_RT <- rrarefy(community_RT, sample = min_sample_size_RT)

min_sample_size_RG <- min(rowSums(community_RG))
rarefied_data_RG <- rrarefy(community_RG, sample = min_sample_size_RG)

# RT
set.seed(123)
solutions_RT <- list()
stresses_RT <- numeric()

for(i in 1:10) {
  nmds_result_RT <- metaMDS(rarefied_data_RT,k =2, trymax = 1000)
  solutions_RT[[i]] <- nmds_result_RT
  stresses_RT[i] <- nmds_result_RT$stress
}

best_solution_RT <- solutions_RT[[which.min(stresses_RT)]]

# RG
set.seed(123)
solutions_RG <- list()
stresses_RG <- numeric()

for(i in 1:10) {
  nmds_result_RG <- metaMDS(rarefied_data_RG,k =2, trymax = 1000)
  solutions_RG[[i]] <- nmds_result_RG
  stresses_RG[i] <- nmds_result_RG$stress
}

best_solution_RG <- solutions_RG[[which.min(stresses_RG)]]


### 8-3. Plot with ggplot
#### 8-3-1. Extract ordinate data from the results as dataframe 
nmds_df_RT <- data.frame(
  MDS1 = best_solution_RT$points[, 1],
  MDS2 = best_solution_RT$points[, 2],
  Region = rownames(community_RT)
)
nmds_df_RT

nmds_df_RG <- data.frame(
  MDS1 = best_solution_RG$points[, 1],
  MDS2 = best_solution_RG$points[, 2],
  Region = rownames(community_RG)
)
nmds_df_RG

#### 8-3-2. Plot
region_order <- c("Khi", "Nor", "Kam", "Tat", "B/V", "GAP", "Fr", "Cau", "Pyr", "SN", "Jp")
nmds_df_RT$Region <- factor(nmds_df_RT$Region, levels = region_order)

NMDS_RT <- ggplot() + theme_bw(base_size = 10) +
  geom_point(data = nmds_df_RT, aes(x = MDS1, y = MDS2, color = Region, shape = Region), size = 6) +
  scale_shape_manual(values = c(15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 17)) +
  scale_color_manual(values = c("#f8766dff", "#db8e00ff", "#aea200ff", "#64b200ff", "#00bd5cff", "#00c1a7ff", "#00badeff", "#00a6ffff", "#b385ffff", "#ef67ebff", "#ff63b6ff")) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  ) +
  annotate("text", x = -1.8, y = 1.17, label = "E", size = 8, hjust = 0) +
  labs(title = paste("Stress =", round(best_solution_RT$stress, 3)), x = "RT: MDS1", y = "RT: MDS2")

NMDS_RT

nmds_df_RG$Region <- factor(nmds_df_RG$Region, levels = region_order)

NMDS_RG <- ggplot() + theme_bw(base_size =10) +
  geom_point(data = nmds_df_RG, aes(x = MDS1, y = MDS2, color = Region , shape = Region), size = 6) +
  scale_shape_manual(values = c(15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 17)) +
  scale_color_manual(values = c("#f8766dff", "#db8e00ff", "#aea200ff", "#64b200ff", "#00bd5cff", "#00c1a7ff", "#00badeff", "#00a6ffff", "#b385ffff", "#ef67ebff", "#ff63b6ff")) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"  # hide legend
  ) +
  annotate("text", x = -0.9, y = 0.65, label = "D", size = 8, hjust = 0) +
  labs(title = paste("Stress =", round(best_solution_RG$stress, 3)), x = "RG: MDS1", y = "RG: MDS2")   

NMDS_RG

NMDS_RG|NMDS_RT
p_nMDS <- NMDS_RG|NMDS_RT

pdf("../_results/Fig3DE_nMDS.pdf", width = 6, height = 4) #for paper
print(p_nMDS)
dev.off()

### 8-4. Pairwise nMDS points
desired_order <- c("Jp", "SN", "Pyr", "Cau", "Fr", "GAP", "B/V", "Tat", "Kam", "Nor", "Khi")
nmds_coords_RT <- nmds_result_RT$points[desired_order, ]  
pairwise_dist_RT <- round(dist(nmds_coords_RT), 2)
pairwise_dist_RT <- as.matrix(pairwise_dist_RT)
pairwise_dist_RT


write.csv(
  x = pairwise_dist_RT,
  file = "../_results/TabS2_nMDS_pairwise_distance_RT.csv",
  row.names = TRUE # set variables as column name
)


nmds_coords_RG <- nmds_result_RG$points[desired_order, ]  
pairwise_dist_RG <- round(dist(nmds_coords_RG), 2)
pairwise_dist_RG <- as.matrix(pairwise_dist_RG)
pairwise_dist_RG

write.csv(
  x = pairwise_dist_RG,
  file = "../_results/TabS3_nMDS_pairwise_distance_RG.csv",
  row.names = TRUE # set variables as column name
)

### 8-5. Mantel test (nMDS vs. GD) for discussion part 
#### 8-5-1. Data preparation for mantel test
pairwise_dist_RT <- read.csv("../_results/TabS2_nMDS_pairwise_distance_RT.csv", row.names = 1)
pairwise_dist_RG <- read.csv("../_results/TabS3_nMDS_pairwise_distance_RG.csv", row.names = 1)
dis.hav_Coord_km <- readRDS("../_data/_processed_data/dis.hav_Coord_km.rds")

colnames(pairwise_dist_RT) <- rownames(pairwise_dist_RT) #column name "B.V" -> "B/V"
colnames(pairwise_dist_RG) <- rownames(pairwise_dist_RG)


# set order "dist_geo" as default
geo_order <- rownames(dis.hav_Coord_km)

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
pairwise_dist_RT_ordered <- reorder_distance_matrix(pairwise_dist_RT, geo_order)
pairwise_dist_RG_ordered <- reorder_distance_matrix(pairwise_dist_RG, geo_order)

#### 8-5-2. Analysis
set.seed(123)
dist_RT_mantel_result <- mantel(log10(dis.hav_Coord_km), 
                              pairwise_dist_RT_ordered, method = "pearson", 
                              permutations = 9999, na.rm = TRUE)
dist_RT_mantel_result

cat("Mantel r (nMDS distance to Geographic distance RT):", round(dist_RT_mantel_result$statistic, 3))
cat("P-value (nMDS distance to Geographic distance RT):", round(dist_RT_mantel_result$signif, 5))
# Mantel r (nMDS distance to Geographic distance RT): 0.715
# P-value (nMDS distance to Geographic distance RT): 1e-04

set.seed(123)
dist_RG_mantel_result <- mantel(log10(dis.hav_Coord_km), 
                                pairwise_dist_RG_ordered, method = "pearson", 
                                permutations = 9999, na.rm = TRUE)
dist_RG_mantel_result

cat("Mantel r (nMDS distance to Geographic distance RG):", round(dist_RG_mantel_result$statistic, 3))
cat("P-value (nMDS distance to Geographic distance RG):", round(dist_RG_mantel_result$signif, 5))
# Mantel r (nMDS distance to Geographic distance RG): 0.616
# P-value (nMDS distance to Geographic distance RG): 0.0029

## 9. Relationship between the number of records for an RT and 
##    the number of regions where a RT was found for Figure S4A
### 9-1. Load data and preparation
source("utils_ecology_myxo.R")
load_ecology_packages()

DF_niv_seq <- readRDS("../_data/_processed_data/DF_biogeo_seq.rds") # sequenced data
table(DF_niv_seq$Code) # check data 

plot1_data <- DF_niv_seq %>%
  group_by(RT) %>%
  summarise(
    Total_Records = n(),  # Number of records for this RT
    Region_Count = n_distinct(Code)  # Number of Region where this RT occurs
  ) %>%
  ungroup() %>%
  arrange(desc(Total_Records))

### 9-3. Plot and save
p1 <- ggplot(plot1_data, aes(x = Total_Records, y = Region_Count)) +
  theme_bw(base_size = 12) +
  geom_point(alpha = 0.6, size = 2, color = "gray30") +
#  geom_smooth(method = "lm", color = "gray", linetype = "dashed", se = FALSE) +
  scale_x_log10(breaks = c(1, 10, 100, 1000), 
                expand = c(0, 0),
                limits = c(1, 500),
                labels = c("1", "10", "100", "1000")
                ) +  # X achse log scale
  labs(x = "Number of records per RT (log scale)", 
       y = "Number of regions with this RT") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 12), 
                     breaks = c(0, 2, 4, 6, 8, 10, 12)) +
  annotate("text", x = 350, y = 11.2, label = "A", col = "black", size = 6)  # annotation
p1


png("../_results/FigS4A_record_RT_detail1.png", width = 500, height = 400)
pdf("../_results/FigS4A_record_RT_detail1.pdf", width = 5, height = 4) #for paper
print(p1)
dev.off()


### 9-3. Cariculate Coefficient of Correlation （after logarithmic transformation）
correlation <- cor(
  log10(plot1_data$Total_Records), 
  plot1_data$Region_Count, 
  use = "complete.obs"
)

### 9-4. Statistics summary
cat("Coefficient of Correlation:", round(correlation, 3), "\n")
cat("Number of RTs:", nrow(plot1_data), "\n")
cat("Range: Number of Records:", min(plot1_data$Total_Records), "-", max(plot1_data$Total_Records), "\n")
cat("Range: Number of Region :", min(plot1_data$Region_Count), "-", max(plot1_data$Region_Count), "\n")

# Summarize the data distribution
summary_stats <- plot1_data %>%
  summarise(
    mean_records = mean(Total_Records),
    median_records = median(Total_Records),
    mean_regions = mean(Region_Count),
    median_regions = median(Region_Count)
  )

print(summary_stats)

## 10. Number of RT occurring in 1, 2, …, 11 regions for Figure S4B
### 10-1. Load data and preparation
plot2_data <- DF_niv_seq %>%
  distinct(RT, Code) %>%  # 重複を排除（出現/非出現のみ）
  count(RT, name = "Region_Count") %>%  # 各RTが出現する地域数
  count(Region_Count, name = "RT_Count") %>%  # 各地域数カテゴリーに含まれるRT数
  arrange(Region_Count)

### 10-2. Plot and save
p2 <- ggplot(plot2_data, aes(x = Region_Count, y = RT_Count)) +
  theme_bw(base_size = 12) +
  geom_col(fill = "gray60", alpha = 0.8, width = 0.8) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0.5, 12.5),
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 350)) +
  geom_text(aes(label = RT_Count), vjust = -0.5, cex = 3) +
  labs(x = "Number of regions", y = "Number of RTs") +
  annotate("text", x = 11.8, y = 325, label = "B", col = "black", size = 6)  # annotation

p2

png("../_results/FigS4B_record_RT_detail2.png", width = 500, height = 400)
pdf("../_results/FigS4B_record_RT_detail2.pdf", width = 5, height = 4) #for paper
print(p2)
dev.off()

## 11. Abundance (number of records) of the RTs in geometric classes for Figure S4C
### 11-1. Load data and preparation
plot3_data <- DF_niv_seq %>%
  count(RT, name = "Total_Count") %>%
  mutate(Frequency_Class = case_when(
    Total_Count == 1 ~ "1",
    Total_Count == 2 ~ "2",
    Total_Count <= 4 ~ "≤4",
    Total_Count <= 8 ~ "≤8", 
    Total_Count <= 16 ~ "≤16",
    Total_Count <= 32 ~ "≤32",
    Total_Count <= 64 ~ "≤64",
    Total_Count <= 128 ~ "≤128",
    TRUE ~ "128<"
  )) %>%
  # Set the class order
  mutate(Frequency_Class = factor(Frequency_Class,
                                  levels = c("1", "2", "≤4", "≤8", "≤16", 
                                             "≤32", "≤64", "≤128", "128<"))) %>%
  count(Frequency_Class, name = "RT_Count")

p3 <- ggplot(plot3_data, aes(x = Frequency_Class, y = RT_Count)) +
  theme_bw(base_size = 12) +
  geom_col(fill = "gray60", alpha = 0.8) +
  geom_text(aes(label = RT_Count), vjust = -0.5, cex = 3) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) +
  labs(x = "Number of records (geometric class)", y = "Number of RTs") +
  annotate("text", x = "128<", y = 185, label = "C", col = "black", size = 6)  # annotation

p3

png("../_results/FigS4C_record_RT_detail3.png", width = 500, height = 400)
pdf("../_results/FigS4C_record_RT_detail3.pdf", width = 5, height = 4) #for paper
print(p3)
dev.off()

 p1|p2|p3
FigS4ABC_record_RT_detail <- p1|p2|p3

pdf("../_results/FigS4ABCD_record_RT_detail.pdf", width = 15, height = 3.5) #for paper
print(FigS4ABC_record_RT_detail)
dev.off()

## 12. Accumulation curve by retion: Figure 5
### 12-1. Data preparation 
# load "DF_niv_seq" (see ### 1-1.)
# calculate "all_region_stats" (see ### 2-1.)
source("utils_ecology_myxo.R")

### 12-3. Adjust regional plot
plots <- list() # create list for storing the adjusted plots

# standard position
plots[["Nor"]] <- create_region_plot(all_region_stats[["Nor"]], "Nor",
                                     show_x_text = TRUE) 
plots$Nor

# adjust position for other regions and store into the list
#GAP
plots[["GAP"]] <- create_region_plot(
  all_region_stats[["GAP"]], "GAP", 
  x_pos_obs_RG = 1200,  
  y_pos_obs_RG = 35,  
  x_pos_obs_RT = 950,
  y_pos_obs_RT = 120,
  x_pos_exp = 1700,  
  y_pos_exp_RG = 60,
  y_pos_exp_RT = 170,
  show_y_axis = FALSE,   # show Y axis label 
  show_y_text = FALSE,   # show Y axis text
  show_x_axis = TRUE,  # show X axis label
  show_x_text = TRUE)   # show X axis text
plots$GAP #check the plot

#BF
plots[["B/V"]] <- create_region_plot(
  all_region_stats[["B/V"]], "B/V",
  x_pos_obs_RG = 500,  
  y_pos_obs_RG = 18,  
  x_pos_obs_RT = 250,
  y_pos_obs_RT = 70,
  x_pos_exp = 1700,  
  y_pos_exp_RG = 43,
  y_pos_exp_RT = 130,
  show_x_text = TRUE
)
plots$"B/V"

#Tat
plots[["Tat"]] <- create_region_plot(
  all_region_stats[["Tat"]], "Tat",
  x_pos_obs_RG = 500,  
  y_pos_obs_RG = 20,  
  x_pos_obs_RT = 250,
  y_pos_obs_RT = 65,
  x_pos_exp = 1700,  
  y_pos_exp_RG = 60,
  y_pos_exp_RT = 90,
  show_x_text = TRUE
)
plots$Tat

#Kam
plots[["Kam"]] <- create_region_plot(
  all_region_stats[["Kam"]], "Kam",
  x_pos_obs_RG = 500,  
  y_pos_obs_RG = 25,  
  x_pos_obs_RT = 250,
  y_pos_obs_RT = 80,
  x_pos_exp = 1700,  
  y_pos_exp_RG = 60,
  y_pos_exp_RT = 115,
  show_x_text = TRUE
)
plots$Kam

#Jp
plots[["Jp"]] <- create_region_plot(
  all_region_stats[["Jp"]], "Jp",
  x_pos_obs_RG = 500,  
  y_pos_obs_RG = 20,  
  x_pos_obs_RT = 250,
  y_pos_obs_RT = 95,
  x_pos_exp = 1700,  
  y_pos_exp_RG = 65,
  y_pos_exp_RT = 170,
  show_y_axis = TRUE,   # show Y axis label 
  show_y_text = TRUE,   # show Y axis text
  show_x_axis = FALSE,  # show X axis label
  show_x_text = TRUE)   # show X axis text
plots$Jp

#Khi
plots[["Khi"]] <- create_region_plot(
  all_region_stats[["Khi"]], "Khi",
  x_pos_obs_RG = 500,  
  y_pos_obs_RG = 25,  
  x_pos_obs_RT = 250,
  y_pos_obs_RT = 95,
  x_pos_exp = 1700,  
  y_pos_exp_RG = 55,
  y_pos_exp_RT = 150,
  show_x_text = TRUE
)
plots$Khi

#Cau
plots[["Cau"]] <- create_region_plot(
  all_region_stats[["Cau"]], "Cau",
  x_pos_obs_RG = 500,  
  y_pos_obs_RG = 28,  
  x_pos_obs_RT = 250,
  y_pos_obs_RT = 110,
  x_pos_exp = 1700,  
  y_pos_exp_RG = 60,
  y_pos_exp_RT = 140,
  show_x_text = TRUE
)
plots$Cau

#SN
plots[["SN"]] <- create_region_plot(
  all_region_stats[["SN"]], "SN",
  x_pos_obs_RG = 500,  
  y_pos_obs_RG = 30,  
  x_pos_obs_RT = 250,
  y_pos_obs_RT = 95,
  x_pos_exp = 1700,  
  y_pos_exp_RG = 120,
  y_pos_exp_RT = 170,
  show_x_text = TRUE
)
plots$SN

#Fr
plots[["Fr"]] <- create_region_plot(
  all_region_stats[["Fr"]], "Fr",
  x_pos_obs_RG = 500,  
  y_pos_obs_RG = 20,  
  x_pos_obs_RT = 270,
  y_pos_obs_RT = 120,
  x_pos_exp = 1700,  
  y_pos_exp_RG = 45,
  y_pos_exp_RT = 170,
  show_x_text = TRUE
)
plots$Fr

#Pyr
plots[["Pyr"]] <- create_region_plot(
  all_region_stats[["Pyr"]], "Pyr",
  x_pos_obs_RG = 500,  
  y_pos_obs_RG = 25,  
  x_pos_obs_RT = 400,
  y_pos_obs_RT = 100,
  x_pos_exp = 1700,  
  y_pos_exp_RG = 55,
  y_pos_exp_RT = 170,
  show_x_text = TRUE
)
plots$Pyr

## 2-4.Plot for Figure 5 with patchwork and save for Figure 5
plots$Jp|plots$SN|plots$Pyr|plots$Cau|plots$Fr|plots$GAP|plots$"B/V"|plots$Tat|plots$Kam|plots$Nor|plots$Khi

Fig5 <- plots$Jp|plots$SN|plots$Pyr|plots$Cau|plots$Fr|plots$GAP|plots$"B/V"|plots$Tat|plots$Kam|plots$Nor|plots$Khi
ggsave("../_results/Fig5_SAC.pdf", Fig5, width = 15, height = 5, dpi = 300)


## 13. Accumulation curve ratio for Figure S10
### 13-1. Function to process all region together 
library(purrr)

create_all_ratio_plots_correct <- function(all_region_stats) {
  map(names(all_region_stats), function(region_name) {
    stats <- all_region_stats[[region_name]]
    
    # Extract data including the Method column
    ratio_data <- inner_join(
      stats$inext_result$iNextEst$size_based %>% 
        filter(Assemblage == "RT") %>% dplyr::select(t, qD, Method),
      stats$inext_result$iNextEst$size_based %>% 
        filter(Assemblage == "RG") %>% dplyr::select(t, qD, Method),
      by = c("t", "Method"), suffix = c("_RT", "_RG")
    ) %>%
      mutate(SAC_Ratio = qD_RG / qD_RT)
    
    # Separate data by Method
    rarefaction_data <- ratio_data %>% filter(Method == "Rarefaction")
    observed_data <- ratio_data %>% filter(Method == "Observed")
    extrapolation_data <- ratio_data %>% filter(Method == "Extrapolation")
    
    # Maximum observation point (maximum t value in the Observed data)
    max_observed_t <- max(observed_data$t)
    
    # plot
    ggplot() +
      # Rarefaction: solid line
      geom_line(data = rarefaction_data, 
                aes(x = t, y = SAC_Ratio),
                color = "black", linewidth = 1.5) +
      # Extrapolation: dotted line
      geom_line(data = extrapolation_data, 
                aes(x = t, y = SAC_Ratio),
                color = "black", linetype = "dashed", linewidth = 1.5) +
      # Observed: point
      geom_point(data = observed_data, 
                 aes(x = t, y = SAC_Ratio),
                 color = "black", size = 5, shape = 16) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, 2000), 
                         labels = remove_last_three_digits) +
      scale_y_continuous(expand = c(0, 0), limits = c(0.2, 1.05)) +
      labs(title = region_name) +
      theme_classic(base_size = 15) +
      theme(legend.position = "none", 
            axis.title.x = element_blank(), 
            axis.title.y = element_blank(), 
            axis.text.y = element_blank())
  }) %>% setNames(names(all_region_stats))
}

# 実行
all_ratio_plots <- create_all_ratio_plots_correct(all_region_stats)
all_ratio_plots[["Nor"]]  # check Nor plot

all_ratio_plots[["Jp"]] <- all_ratio_plots[["Jp"]] +
  theme(axis.title.y = element_text(angle = 90),
        axis.text.y = element_text()) +
  ylab("RG/RT Ratio")

all_ratio_plots[["GAP"]] <- all_ratio_plots[["GAP"]] +
  theme(axis.title.x = element_text()) +
  xlab("Specimen (x100)")

pr10 <- all_ratio_plots[["Nor"]]  
pr6 <- all_ratio_plots[["GAP"]]
pr7 <- all_ratio_plots[["B/V"]]
pr8 <- all_ratio_plots[["Tat"]]
pr9 <- all_ratio_plots[["Kam"]]
pr1 <- all_ratio_plots[["Jp"]]
pr11 <- all_ratio_plots[["Khi"]]
pr4 <- all_ratio_plots[["Cau"]]
pr2 <- all_ratio_plots[["SN"]]
pr5 <- all_ratio_plots[["Fr"]]
pr3 <- all_ratio_plots[["Pyr"]]


pr1|pr2|pr3|pr4|pr5|pr6|pr7|pr8|pr9|pr10|pr11
FigS10 <- pr1|pr2|pr3|pr4|pr5|pr6|pr7|pr8|pr9|pr10|pr11

png("../_results/FigS10_SAC_ratio.png", width = 1000, height = 300)
print(FigS10)
dev.off()
ggsave("../_results/FigS10_SAC_ratio.pdf", FigS10, width = 15, height = 5, dpi = 300)


## 14. Basic statistics for Figure S7ABC: RAD
### 14-1. Data preparation: sort clade/RG/RT by count 
#RT
df_summary_rt <- DF_niv_seq %>%
  group_by(RT) %>%
  summarize(count = n()) %>%
  group_by(RT) %>%
  summarize(total_count = sum(count)) %>%
  arrange(desc(total_count))

abund_rt <- df_summary_rt$total_count #convert abundance row to vector

#RG
df_summary_rg <- DF_niv_seq %>%
  group_by(RG) %>%
  summarize(count = n()) %>%
  group_by(RG) %>%
  summarize(total_count = sum(count)) %>%
  arrange(desc(total_count))

abund_rg <- df_summary_rg$total_count #convert abundance row to vector

#Clade
df_summary_clade <- DF_niv_seq %>%
  group_by(Clade) %>%
  summarize(count = n()) %>%
  group_by(Clade) %>%
  summarize(total_count = sum(count)) %>%
  arrange(desc(total_count))

abund_clade <- df_summary_clade$total_count #convert abundance row to vector

### 14-2. Model-fitting and calculate AIC, BIC: with vegan
fit_all_rt <- radfit(abund_rt) 
fit_all_rg <- radfit(abund_rg) 
fit_all_clade <- radfit(abund_clade) 
plot(fit_all_rt, main = "")
plot(fit_all_rg, main = "")
plot(fit_all_clade, main = "")

par(mfrow = c(1, 1))

png("../_results/FigS7C_plot_RT.png", width = 600, height = 500)
pdf("../_results/FigS7C_plot_RT.pdf", width = 8, height = 7) #for paper
plot(fit_all_rt, main = "")
text(2, 1.3, "C", cex = 3)
dev.off()
cat("RT AIC:", round(AIC(fit_all_rt), 0), "\n")
# RT AIC: 5193 3681 1826 2414 1634 

png("../_results/FigS7B_plot_RG.png", width = 600, height = 500)
pdf("../_results/FigS7B_plot_RG.pdf", width = 8, height = 7) #for paper
plot(fit_all_rg, main = "")
text(2, 1.3, "B", cex = 3)
dev.off()
cat("RG AIC:", round(AIC(fit_all_rg), 0),  "\n")
# RG AIC: 3845 765 1420 2649 769 

png("../_results/FigS7A_plot_clade.png", width = 600, height = 500)
pdf("../_results/FigS7A_plot_clade.pdf", width = 8, height = 7) #for paper
plot(fit_all_clade, main = "")
text(2, 1.3, "A", cex = 3)
dev.off()
cat("Clade AIC:", round(AIC(fit_all_clade), 0), "\n")
# Clade AIC: 1663 547 1312 2467 551 



## 15. Basic statistics for Figure S7D: SAC & Alpha diversity
### 15-1. Create community data for clade, RG, and RT level for the complete data set
comm_list<- create_multi_level_communities(DF_niv_seq, "final_colno", c("RT", "RG", "Clade"))
res <- plot_combined_sac(comm_list, 
                         endpoint_inext = 10000,
                         colors = c("steelblue", "tomato3", "olivedrab", "goldenrod", "purple"),
                         title = "", 
                         xlab = "Number of Specimens (x100)", 
                         ylab = "RT/RG/Clade",
                         limits_y = c(0, 670),
                         limits_x = c(0, 10500), 
                         anotate_x = 400,
                         anotate_y = 630,
                         annotate_label = "D")

## 15-2. Check visulalized SAC plot and results of alpha diversity
res$plot
res$inext
richness_results <- filter(res$inext$AsyEst, Diversity == "Species richness")
cat("Estimated Clade richness:\n")
cat(round(richness_results[1, "Estimator"], 1), "±", round(richness_results[1, "s.e."], 1))
cat("Estimated RG richness:\n")
cat(round(richness_results[2, "Estimator"], 1), "±", round(richness_results[2, "s.e."], 1))
cat("Estimated RT richness:\n")
cat(round(richness_results[3, "Estimator"], 1), "±", round(richness_results[3, "s.e."], 1))

ggsave(
  filename = "../_results/FigS7D_total_sac_plot.pdf", 
  plot = res$plot,                    
  width = 8,                          
  height = 7,                         
  dpi = 300                           
)


data_to_save <- res$inext$AsyEst
write.csv(
  x = data_to_save,
  file = "../_results/combined_sac_inext_data.csv",
  row.names = TRUE
)
 

## 16. Basic statistics for Figure S11AB: Alpha diversity
### 16-1. Data preparation
alpha_div_rt <- rarefied_rt_df %>% 
  left_join(shannon_table[, -3], by = "Code") %>% 
  left_join(simpson_table[, -3], by = "Code")

alpha_div_rg <- rarefied_rg_df %>% 
  left_join(shannon_table[, -2], by = "Code") %>% 
  left_join(simpson_table[, -2], by = "Code")

# reshape the data                         
data_melted_rt <- melt(alpha_div_rt, 
                       id.vars = "Code", 
                       variable.name = "Index", 
                       value.name = "Value")

data_melted_rg <- melt(alpha_div_rg, 
                       id.vars = "Code", 
                       variable.name = "Index", 
                       value.name = "Value")

# 16-2. Plot
custom_order <- c("Jp", "SN", "Pyr", "Cau", "Fr", "GAP", "B/V", "Tat", "Kam", "Nor", "Khi") 
index_order_RT <- c("Rarefied_RT", "Shannon_RT", "Simpson_RT")
index_order_RG <- c("Rarefied_RG", "Shannon_RG", "Simpson_RG")

data_melted_rt$Index <- factor(data_melted_rt$Index, levels = index_order_RT)
data_melted_rg$Index <- factor(data_melted_rg$Index, levels = index_order_RG)

p1_RT <- ggplot(data_melted_rt, aes(x = factor(Code, levels = custom_order), y = Value, fill = Code)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_bw(base_size= 30) +
  labs(title = "B",
       x = "Code",
       y = "Value (RT)") +
  facet_wrap(~ Index, scales = "free_y") +
  scale_fill_manual(values = c(
    "Jp" = "#ff63b6ff",
    "SN" = "#ef67ebff",
    "Pyr" = "#b385ffff",
    "Cau" = "#00a6ffff",
    "Fr" = "#00badeff",
    "GAP" = "#00c1a7ff",
    "B/V" = "#00bd5cff",
    "Tat" = "#64b200ff",
    "Kam" = "#aea200ff",
    "Khi" = "#f8766dff",
    "Nor" = "#db8e00ff"
  ))+
  theme(legend.position="None", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, size=15), 
        axis.title.x = element_text(vjust = 0, size=15),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(vjust = 0, size=15),
        plot.margin = unit(c(1, 1, 1, 1), "cm"))  #上下左右の余白をXXcmに設定
p1_RT

p1_RG <- ggplot(data_melted_rg, aes(x = factor(Code, levels = custom_order), y = Value, fill = Code)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_bw(base_size= 30) +
  labs(title = "B",
       x = "Code",
       y = "Value (RG)") +
  facet_wrap(~ Index, scales = "free_y") +
  scale_fill_manual(values = c(
    "Jp" = "#ff63b6ff",
    "SN" = "#ef67ebff",
    "Pyr" = "#b385ffff",
    "Cau" = "#00a6ffff",
    "Fr" = "#00badeff",
    "GAP" = "#00c1a7ff",
    "B/V" = "#00bd5cff",
    "Tat" = "#64b200ff",
    "Kam" = "#aea200ff",
    "Khi" = "#f8766dff",
    "Nor" = "#db8e00ff"
  ))+
  theme(legend.position="None", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, size=15), 
        axis.title.x = element_text(vjust = 0, size=15),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(vjust = 0, size=15),
        plot.margin = unit(c(1, 1, 1, 1), "cm"))  #上下左右の余白をXXcmに設定
p1_RG

p1_RG|p1_RT

ggsave("../_results/FigS11_rarefy_alpha_new.pdf", width= 20, height= 10)

## 17. Michaelis-Menten kinetics curve for Figure S8
### 17-1.
source("utils_ecology_myxo.R")
DF_niv_seq <- readRDS("../_data/_processed_data/DF_biogeo_seq.rds") # sequenced data

### 17-2. Set labels to show on the plot
labels_to_show <- c("DIDmic", "DDYdub", "LAMovo", "LAMove", "LAMsau", "PLScha", "PLScar-gra",
                    "LAMarc-plv-spi/LAMcrs/DCHmet", "PHYver",
                    "DIDalp", "DIDmey", "DIDeur", "LAMesp", "PHYalb", "MERcar-cri-spi-ech")

plot_data <- prepare_mm_data(
  data = DF_niv_seq,
  group_var = "Clade",
  count_var = "RT",
  label_groups = labels_to_show,
  min_records = 1,    
  min_ribotypes = 1   
)

### 17-3. Search for optimal initial value
optimal_starts <- find_optimal_start(plot_data)

### 17-4. Plot
result_mm <- fit_mm_model(
  data = plot_data,
  Vmax_start = 104,
  Km_start = 25,
  x_var = "num_records",
  y_var = "num_ribotypes",
  label_var = "label",
  plot_title = "",
  x_lab = "Number of Samples",
  y_lab = "Number of Ribotypes"
)

# Check the used number of data
cat("Number of data:", nrow(plot_data), "\n")
cat("Missing values for variable x:", sum(is.na(plot_data$num_records)), "\n")
cat("Missing values for variable y:", sum(is.na(plot_data$num_ribotypes)), "\n")
cat("Number of data points used in the model:", nobs(result_mm$model), "\n")

# Adjust for final plot
final_plot <- result_mm$plot +
  scale_x_continuous(
    limits = c(0, 780),
    breaks = seq(0, 750, by = 200),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 55),
    breaks = seq(0, 50, by = 20),
    expand = c(0, 0)
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 12)
  )

print(final_plot)
png("../_results/FigS8_mm_clade_RT.png", width = 600, height = 400)
pdf("../_results/FigS8_mm_clade_RT.pdf", width = 8, height = 6) #for paper
print(final_plot)
dev.off()

### 17-5. Diagnose plot
diagnose_mm_model(result_mm, plot_data)
interpret_diagnostic_plots(result_mm, plot_data)

## 18. Latitudinal gradient for Figure S12
### 18-1. Load data and preparation
basic_statistics <- readRDS("../_data/_processed_data/basic_statistics.rds")
basic_rf <- dplyr::select(basic_statistics, Code, Rarefied_RT, Rarefied_RG, Rarefied_clade)
GeoCoord <- readRDS("../_data/_processed_data/GeoCoord_avg_11regions.rds")
GeoCoord_xy <- dplyr::select(GeoCoord, Code, Easting_avg, Northing_avg)

rf_grad <- left_join(GeoCoord_xy, basic_rf)

### 18-2. Simple peason r calculation
cor.test(rf_grad$Northing_avg, rf_grad$Rarefied_RT)
cor.test(rf_grad$Northing_avg, rf_grad$Rarefied_RG)
cor.test(rf_grad$Northing_avg, rf_grad$Rarefied_clade)

rf_grad_long <- rf_grad %>%
  pivot_longer(cols = starts_with("Rarefied_"),
               names_to = "Index_Type",
               values_to = "Rarefied_Value")

### 18-3. Plot
#### 18-3-1. Calculate statistics for all OTUs
stats_labels <- rf_grad_long %>%
  # Index_Type (Rarefied_RT, Rarefied_RG, ...) 
  group_by(Index_Type) %>%
  # lm by group and extract the summary
  do({
    model <- lm(Rarefied_Value ~ Northing_avg, data = .)
    r2_val <- broom::glance(model)$adj.r.squared
    p_beta1_val <- broom::tidy(model) %>% filter(term == "Northing_avg") %>% pull(p.value)
    
    r2_label <- paste0(
      "'Adj. ' * italic(R) ^ 2 ~ ' = ' * ", # combine 'Adj.'and '=' as character
      round(r2_val, 3)
    )
    
    p_label <- if (p_beta1_val < 0.001) {
      "'p < 0.001'" # italic "p"
    } else {
      paste0("italic(p) ~ '=' ~ ", round(p_beta1_val, 3))
    }
    
    tibble(
      r2_label_expr = r2_label,
      p_label_expr = p_label
    )
  })

#### 18-3-2. Prepare statistics data
# data for each plot R
rt_data <- rf_grad %>% 
  dplyr::select(Code, Northing_avg, Rarefied_RT) %>%
  rename(Rarefied_Value = Rarefied_RT)

rg_data <- rf_grad %>% 
  dplyr::select(Code, Northing_avg, Rarefied_RG) %>%
  rename(Rarefied_Value = Rarefied_RG)

clade_data <- rf_grad %>% 
  dplyr::select(Code, Northing_avg, Rarefied_clade) %>%
  rename(Rarefied_Value = Rarefied_clade)

# basic statistics for each plot Rarefied_RT
rt_stats <- stats_labels %>%
  filter(Index_Type == "Rarefied_RT")

rg_stats <- stats_labels %>%
  filter(Index_Type == "Rarefied_RG")

clade_stats <- stats_labels %>%
  filter(Index_Type == "Rarefied_clade")


#### 18-3-3. Plot
# RT
plot_RT <- ggplot(rt_data, aes(x = Northing_avg, y = Rarefied_Value)) +
  
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") + 
  geom_text(aes(label = Code), vjust = -0.8, hjust = 0.5, size = 3.5, check_overlap = TRUE) +
  
  geom_text(data = rt_stats,
            aes(x = 57, y = 62, label = r2_label_expr), 
            hjust = 0, vjust = 1, size = 4, parse = TRUE) +
  
  # p-value label
  geom_text(data = rt_stats,
            aes(x = 57, y = 56, label = p_label_expr), 
            hjust = 0, vjust = 1, size = 4, parse = TRUE) +
  annotate("text", 
           x = min(rt_data$Northing_avg), 
           y = max(rt_data$Rarefied_Value), label = "C", size = 7, hjust = 0) +
  labs(
    title = "",
    x = "Northing",
    y = "Rarefied RT value"
  ) +
  theme_bw(base_size = 12)

plot_RT

# RG
plot_RG <- ggplot(rg_data, aes(x = Northing_avg, y = Rarefied_Value)) +
  
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") + 
  geom_text(aes(label = Code), vjust = -0.8, hjust = 0.5, size = 3.5, check_overlap = TRUE) +
  
  geom_text(data = rg_stats,
            aes(x = 57, y = 30, label = r2_label_expr), 
            hjust = 0, vjust = 1, size = 4, parse = TRUE) +
  
  geom_text(data = rg_stats,
            aes(x = 57, y = 27, label = p_label_expr),
            hjust = 0, vjust = 1, size = 4, parse = TRUE) +
  annotate("text", 
           x = min(rg_data$Northing_avg), 
           y = max(rg_data$Rarefied_Value), label = "B", size = 7, hjust = 0) +
  labs(
    title = "",
    x = "Northing",
    y = "Rarefied RG value"
  ) +
  theme_bw(base_size = 12)

plot_RG


# Clade
plot_clade <- ggplot(clade_data, aes(x = Northing_avg, y = Rarefied_Value)) +
  
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") + 
  geom_text(aes(label = Code), vjust = -0.8, hjust = 0.5, size = 3.5, check_overlap = TRUE) +
  
  geom_text(data = clade_stats,
            aes(x = 57, y = 22.5, label = r2_label_expr), # 左上配置
            hjust = 0, vjust = 1, size = 4, parse = TRUE) +
  
  geom_text(data = clade_stats,
            aes(x = 57, y = 21, label = p_label_expr), # 少し下に配置
            hjust = 0, vjust = 1, size = 4, parse = TRUE) +
  annotate("text", 
           x = min(clade_data$Northing_avg), 
           y = max(clade_data$Rarefied_Value), label = "A", size = 7, hjust = 0) +
  labs(
    title = "",
    x = "Northing",
    y = "Rarefied Clade value"
  ) +
  theme_bw(base_size = 12)

plot_clade

plot_clade|plot_RG|plot_RT
FigS12 <- plot_clade|plot_RG|plot_RT

pdf("../_results/FigS12_latitude_rf.pdf", width = 15, height = 5)
png("../_results/FigS12_latitude_rf.png", width = 900, height = 300) #for paper
print(FigS12)
dev.off()
 