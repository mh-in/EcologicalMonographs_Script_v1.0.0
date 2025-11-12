# 01. data_cleaning_processing

## 1. Data cleaning
### 1-1. Load raw data and function from "utils_ecology_myxo.R"
DF <- read.csv("../_data/All_records.csv",
              skip = 3,
              header = TRUE)
str(DF)

source("utils_ecology_myxo.R")
load_ecology_packages()


### 1-2. Change column name
colnames(DF)[colnames(DF) == "rec_no"] <- "All_Rec_No"
colnames(DF)[colnames(DF) == "rec_no.1"] <- "Regio_Rec_No"
colnames(DF)[colnames(DF) == "F.M"] <- "ID_on"
colnames(DF)[colnames(DF) == "Morph.ID"] <- "Morph_ID"
colnames(DF)[colnames(DF) == "Mol.ID"] <- "Mol_ID"
colnames(DF)[colnames(DF) == "RG..p9."] <- "RG"
colnames(DF)[colnames(DF) == "Clade.1"] <- "Clade_ch"
colnames(DF)[colnames(DF) == "seq.attempt"] <- "seq_attempt"
colnames(DF)[colnames(DF) == "Myxo.type"] <- "Myxo_Typ"
colnames(DF)[colnames(DF) == "Morph.Remark"] <- "Morph_Remark"
colnames(DF)[colnames(DF) == "RT.Nr"] <- "RT"
colnames(DF)[colnames(DF) == "seq.quality.etc."] <- "Seq_Remark"
colnames(DF)[colnames(DF) == "Reference.RT"] <- "Reference_Seq"
colnames(DF)[colnames(DF) == "RT name"] <- "RT_Same"
colnames(DF)[colnames(DF) == "prec"] <- "Prec_x"
colnames(DF)[colnames(DF) == "prec1"] <- "Prec_y"

str(DF)

### 1-3. Convert regional code to regional name
DF <- convert_region_codes(DF, code_column = "Code")

### 1-4. Convert data type
DF$RT_Nr <- as.character(DF$RT)
DF$RG <- as.character(DF$RG)
DF$Easting <- as.numeric(DF$Easting)

### 1-5. Data preparation
DF_niv_rec <- select(DF, Code, final_colno, Morph_ID, CF,  # all records
                     Myxo_Typ, seq_attempt, Clade, RT, RG, 
                     Clade_ch, GenBank, Easting, Northing, Elev, Subplot, date, Country, Region, Municipality, Locality)
count(DF_niv_rec)
myxo_type <- count(DF_niv_rec, Myxo_Typ)
myxo_type
#   Myxo_Typ    n
# 1       --   19   <- perhaps myxomycete but inidentifiable (= did not attempt sequencing only included as records) 
# 2     nivi 6267   <- nivicolous including both dark- and bright spored myxomycetes (="TRI" must be excluded)
# 3 non-nivi   17   <- non-nivicolous (= exclude for analysis)
# 4 non/nivi    7   <- maybe nivicolous (= include for analysis)
count(DF, seq_attempt)
#   seq_attempt    n
# 1             6069   <- attempted sequencing = "nivi" and "non/nivi"
# 2          no  241   <- these may either not applicable/bright-spored nivicolous, non-nivicolous, or no myxos

# calculate Trichia spp.
target_string <- "TRI"
column_to_check <- DF_niv_rec$Morph_ID
is_match <- grepl(target_string, column_to_check)
count_result <- sum(is_match)


DF_niv <- filter(DF_niv_rec, Myxo_Typ != "non-nivi") # exclude non-nivicolous
DF_niv <- filter(DF_niv, seq_attempt != "no") # exclude scanty/sclerotized nivicolous, non-nivicolous, or no myxos
count(DF_niv)
# n
# 1 6069 <-  data with sequencing attempt

DF_seq_niv <- filter(DF_niv, Clade !="NA", Clade != "") # barcoded specimens = exclude specimens failed sequencing
count(DF_seq_niv)
#      n
# 1 5546 <-  barcoded specimens

### 1-6. Final check
cat("Total records:\n")
print(count(DF_niv_rec) - 19)  # 19 = no myxos (not attempt sequencing) 
cat("\nRegional distribution:\n")
print(table(DF_niv_rec$Code))
#  B/V  Cau   Fr  GAP   Jp  Kam  Khi  Nor  Pyr   SN  Tat 
# 658  308  538 1646  477  314  272  329  745  370  653 
print(paste("no myxos:", myxo_type[1, 2])) 
print(paste("non-nivi:", myxo_type[3, 2])) 
print(paste0("Trichia spp.: ", count_result))

cat("Specimens with sequencing attempt (all dark-spored Myxomycetes for analysis):\n")

cat("\nRegional distribution:\n")
print(table(DF_niv$Code))

cat("Barcoded specimens:\n")
count(DF_seq_niv)
cat("\nRegional distribution:\n")
print(table(DF_seq_niv$Code))

print(paste0("Total sequence coverage: ", count(DF_seq_niv)/count(DF_niv)*100))


### 1-7. Save relevant data as rds
saveRDS(DF_niv_rec, "../_data/_processed_data/DF_biogeo_dark_rec.rds") # all records
saveRDS(DF_niv, "../_data/_processed_data/DF_biogeo_dark.rds") # data with sequencing attempt
saveRDS(DF_seq_niv, "../_data/_processed_data/DF_biogeo_seq.rds") # barcoded specimens


## 2. Processing for analysis: Beta diversity calculation
### 2-1. Load RDS data and function from "utils_ecology_myxo.R"

DF <- readRDS("../_data/_processed_data/DF_biogeo_seq.rds")

source("utils_ecology_myxo.R")
load_ecology_packages()

### 2-2. Create community matrix for Clade, RG and RT
#### 2-2-1. for general analysis
community_clade <- create_community_matrix(DF, "Code", "Clade")
community_RG <- create_community_matrix(DF, "Code", "RG")
community_RT <- create_community_matrix(DF, "Code", "RT")
cat("Number of clade:", ncol(community_clade), ", Number of regions:", nrow(community_clade))
cat("Number of RG:", ncol(community_RG), ", Number of regions:", nrow(community_RG))
cat("Number of RT:", ncol(community_RT), ", Number of regions:", nrow(community_RT))

saveRDS(community_clade, "../_data/_processed_data/community_clade.rds")
saveRDS(community_RG, "../_data/_processed_data/community_RG.rds")
saveRDS(community_RT, "../_data/_processed_data/community_RT.rds")

#### 2-2-2. for dbRDA analysis
DF_seq_niv <- readRDS("../_data/_processed_data/DF_biogeo_seq.rds")
DF_seq_niv <- dplyr:: select(DF_seq_niv, final_colno, Code, RT, Subplot)

# Not Kam
DF_seq_niv_RT_1 <- DF_seq_niv %>% 
  filter(Code != "Kam") %>% 
  mutate(Subplot_code = Code) %>% 
  select(-Subplot, -Code)

# Kam
DF_seq_niv_RT_2 <- DF_seq_niv %>% 
  filter(Code == "Kam") %>% 
  mutate(Subplot_code = Subplot) %>% 
  select(-Subplot, -Code)

DF_seq_niv_rda<- rbind(DF_seq_niv_RT_1, DF_seq_niv_RT_2)

# Create community matrix for RT for 12 regions
community_RT_rda <- create_community_matrix(DF_seq_niv_rda, "Subplot_code", "RT")
cat("Number of RT:", ncol(community_RT_rda), ", Number of regions:", nrow(community_RT_rda))

saveRDS(community_RT_rda, "../_data/_processed_data/community_RT_rda.rds")


### 2-3. Beta diversity calculation: D1-D3 (as dissimilarity)
#### 2-3-1. Calculate diversity index 
beta_results_rt <- calculate_beta_diversity1(community_RT)
beta_results_rg <- calculate_beta_diversity1(community_RG)

#### 2-3-2. Convert to matrix
dist_d1_rt <- beta_to_distmatrix(beta_results_rt, "D1")
dist_d2_rt <- beta_to_distmatrix(beta_results_rt, "D2")
dist_d3_rt <- beta_to_distmatrix(beta_results_rt, "D3")
dist_d1_rt
dist_d2_rt
dist_d3_rt

dist_d1_rg <- beta_to_distmatrix(beta_results_rg, "D1")
dist_d2_rg <- beta_to_distmatrix(beta_results_rg, "D2")
dist_d3_rg <- beta_to_distmatrix(beta_results_rg, "D3")
dist_d1_rg
dist_d2_rg
dist_d3_rg


#### 2-3-3. Save matrix
saveRDS(dist_d1_rt, "../_results/beta_diversity/D1_RT.rds")
saveRDS(dist_d2_rt, "../_results/beta_diversity/D2_RT.rds")
saveRDS(dist_d3_rt, "../_results/beta_diversity/D3_RT.rds")
saveRDS(dist_d1_rg, "../_results/beta_diversity/D1_RG.rds")
saveRDS(dist_d2_rg, "../_results/beta_diversity/D2_RG.rds")
saveRDS(dist_d3_rg, "../_results/beta_diversity/D3_RG.rds")


### 2-4. Calculate beta diversity index: D4ab, D5, D6 (as dissimilarity)
rg_dissimilarity <- calculate_beta_diversity2(community_RG, "RG")
rt_dissimilarity <- calculate_beta_diversity2(community_RT, "RT")
rt_dissimilarity_rda <- calculate_beta_diversity2(community_RT_rda, "RT_rda")

rg_dissimilarity$D5_horn
rt_dissimilarity$D5_horn
rt_dissimilarity_rda$D5_horn

# Save results as 
save_beta_diversity(rg_dissimilarity, "../_results/beta_diversity", "beta_diversity")
save_beta_diversity(rt_dissimilarity, "../_results/beta_diversity", "beta_diversity")
save_beta_diversity(rt_dissimilarity_rda, "../_results/beta_diversity", "beta_diversity")


### 2-5. Check dataset
rg_dissimilarity <- load_beta_diversity("../_results/beta_diversity/beta_diversity_RG.rds")
rt_dissimilarity <- load_beta_diversity("../_results/beta_diversity/beta_diversity_RT.rds")
rt_dissimilarity_rda <- load_beta_diversity("../_results/beta_diversity/beta_diversity_RT_rda.rds")

rg_dissimilarity$D5_horn
rt_dissimilarity$D5_horn
rt_dissimilarity_rda$D5_horn

print_beta_diversity_summary(rg_dissimilarity)
print_beta_diversity_summary(rt_dissimilarity)
print_beta_diversity_summary(rt_dissimilarity_rda)

### 2-7. Calculate beta diversity index: D7ab
#### 2-7-1. load treefile

library(ape) #to read treefile
tree <- read.tree("../_data/473RT_mafft.treefile")

str(tree$tip.label) # Check tip label

#### 2-7-2. Modify tip label
tree$tip.label <- gsub("RT([0-9]+)_.*", "\\1", tree$tip.label) # extract the number after "RT" and create new label
head(tree$tip.label) # check new label

#### 2-7-3. Calculate Branch Length
tree_dist_matrix <- cophenetic(tree)  #ape::cophenetic()
str(tree_dist_matrix)

saveRDS(tree_dist_matrix, "../_data/_processed_data/tree_dist_matrix.rds")

#### 2-7-4. Create OTU Table with the Length Weighing the Specimen Count
community_RT_t <- t(community_RT) #transpose
head(community_RT_t)

community_RT_rda_t <- t(community_RT_rda) #transpose
head(community_RT_rda_t)

library(phyloseq)
otu_table <- otu_table(as.matrix(community_RT_t), taxa_are_rows = TRUE) # Create OTU table
head(otu_table)

otu_table_rda <- otu_table(as.matrix(community_RT_rda_t), taxa_are_rows = TRUE) # Create OTU table
head(otu_table_rda)

#### 2-7-4. Root teh tree at midpoint
library(phangorn) #midpoint()
rooted_tree <- midpoint(tree)
#head(rooted_tree)
saveRDS(rooted_tree, "../_data/_processed_data/rooted_tree.rds") 

#### 2-7-5. Create Phyloseq Object

physeq_object <- phyloseq(otu_table, phy_tree(rooted_tree))
physeq_object_rda <- phyloseq(otu_table_rda, phy_tree(rooted_tree))

#### 2-7-6. Calculate UniFrac Distance
weighted_unifrac <- UniFrac(physeq_object, weighted = TRUE) # D7a: weighted UniFrac
weighted_unifrac

weighted_unifrac_rda <- UniFrac(physeq_object_rda, weighted = TRUE) # D7a: weighted UniFrac
weighted_unifrac_rda

unweighted_unifrac <- UniFrac(physeq_object, weighted = FALSE) # D7b: unweighted UniFrac
unweighted_unifrac

unweighted_unifrac_rda <- UniFrac(physeq_object_rda, weighted = FALSE) # D7b: unweighted UniFrac
unweighted_unifrac_rda

saveRDS(weighted_unifrac, "../_results/beta_diversity/D7_weighted.rds")
saveRDS(unweighted_unifrac, "../_results/beta_diversity/D7_unweighted.rds")

saveRDS(weighted_unifrac_rda, "../_results/beta_diversity/D7_weighted_rda.rds")
saveRDS(unweighted_unifrac_rda, "../_results/beta_diversity/D7_unweighted_rda.rds")

### 3. Geographic Coordinate
#### 3-1. Average geographic coordinates for 11 mountains (rough coordinates for distance calculation)
DF_all <- readRDS("../_data/_processed_data/DF_biogeo_dark_rec.rds") # all record

GeoCoord <- DF_all %>% 
  select(Code, Easting, Northing, Elev) %>% 
  group_by(Code) %>% 
  summarise(Easting_avg = round(mean(Easting, na.rm = T), 1),
            Easting_max = round(max(Easting, na.rm = T), 1),
            Easting_min = round(min(Easting, na.rm = T), 1),
            Northing_avg = round(mean(Northing, na.rm = T), 1),
            Northing_max = round(max(Northing, na.rm = T), 1),
            Northing_min = round(min(Northing, na.rm = T), 1),
            Elevation_avg = round(mean(Elev, na.rm = T), 1),
            Elevation_min = round(min(Elev, na.rm = T), 1),
            Elevation_max = round(max(Elev, na.rm = T), 1)) %>% 
  arrange(desc(Northing_avg))
print(GeoCoord)
saveRDS(GeoCoord, "../_data/_processed_data/GeoCoord_avg_11regions.rds")

#### 3-2. Average geographic coordinates for region detail (detail coordinates for calculation distance)
load_ecology_packages()
# For Nor, Khi, Kam, 
GeoCoord_detail_1 <- DF_all %>% 
  filter(Code == "Nor" | Code == "Khi" | Code == "Kam") %>% 
  dplyr::select(Code, Region, Municipality, Easting, Northing, Elev) %>% 
  group_by(Code, Municipality) %>% 
  summarise(Easting_avg = mean(Easting, na.rm = T),
            Northing_avg = mean (Northing, na.rm = T),
            Elevation_avg = round(mean (Elev, na.rm = T), 1),
            Elevation_min = round(min(Elev, na.rm = T), 1),
            Elevation_max = round(max(Elev, na.rm = T), 1)) %>% 
  arrange(desc(Northing_avg))

colnames(GeoCoord_detail_1)[colnames(GeoCoord_detail_1) == "Municipality"] <- "Subplot"
GeoCoord_detail_1$Subplot_Nr <- c("Khi2", "Khi1", "Nor3", "Nor2", "Nor1", "Kam2", "Kam1")

#For Jp
GeoCoord_detail_2 <- DF_all %>% 
  filter(Code == "Jp") %>% 
  dplyr::select(Code, Region, Municipality, Easting, Northing, Elev) %>% 
  group_by(Code, Region) %>% 
  summarise(Easting_avg = mean(Easting, na.rm = T),
            Northing_avg = mean (Northing, na.rm = T),
            Elevation_avg = round(mean (Elev, na.rm = T), 1),
            Elevation_min = round(min(Elev, na.rm = T), 1),
            Elevation_max = round(max(Elev, na.rm = T), 1)) %>% 
  arrange(desc(Northing_avg))

colnames(GeoCoord_detail_2)[colnames(GeoCoord_detail_2) == "Region"] <- "Subplot"
subplot_mapping_jp <- c(
  "1" = "Myoko", "2" = "Shiga-H", "3" = "Yatsugatake", "4" = "Sugadaira", 
  "5" = "Norikura"
  )
GeoCoord_detail_2$Subplot <- subplot_mapping_jp[as.character(GeoCoord_detail_2$Subplot)]
GeoCoord_detail_2$Subplot_Nr <- c("Jp5", "Jp4", "Jp3", "Jp2", "Jp1")


# For BF1: Feldberg
GeoCoord_detail_BF1 <- DF_all %>% 
  filter(Locality == "Feldberg") %>%  #Feldberg in Germany
  dplyr::select(Code, Locality, Municipality, Easting, Northing, Elev) %>% 
  group_by(Code) %>% 
  summarise(Easting_avg = mean(Easting, na.rm = T),
            Northing_avg = mean (Northing, na.rm = T),
            Elevation_avg = round(mean (Elev, na.rm = T), 1),
            Elevation_min = round(min(Elev, na.rm = T), 1),
            Elevation_max = round(max(Elev, na.rm = T), 1)) %>% 
  arrange(desc(Northing_avg)) %>% 
  mutate(Subplot = "Feldberg in Germany")
GeoCoord_detail_BF1$Subplot_Nr <- "B/V1"

# For BF2: France
GeoCoord_detail_BF2 <- DF_all %>% 
  filter(Locality == "Col de la Schlucht" | Locality == "Vogesen") %>%  #two subplots in France
  dplyr::select(Code, Locality, Municipality, Easting, Northing, Elev) %>% 
  group_by(Code) %>% 
  summarise(Easting_avg = mean(Easting, na.rm = T),
            Northing_avg = mean (Northing, na.rm = T),
            Elevation_avg = round(mean (Elev, na.rm = T), 1),
            Elevation_min = round(min(Elev, na.rm = T), 1),
            Elevation_max = round(max(Elev, na.rm = T), 1)) %>% 
  arrange(desc(Northing_avg)) %>% 
  mutate(Subplot = "France")
GeoCoord_detail_BF2$Subplot_Nr <- "B/V2"

# For BF3: Others in Germany
GeoCoord_detail_BF3 <- DF_all %>% 
  filter(Locality == "Nationalpark Schwarzwald, Hornisgrinde" |
           Locality == "Nationalpark Schwarzwald, Seibleseckle" | 
           Locality == "Seibleseckle" |
           Locality == "Nationalpark Schwarzwald, Wilder See" |
           Locality == "Nationalpark Schwarzwald, Vogelskopf" |
           Locality == "Nationalpark Schwarzwald, Schliffkopf" |
           Locality == "Nationalpark Schwarzwald, Bärenteich" |
           Locality == "Nationalpark Schwarzwald, Rossecklestraße" |
           Locality == "Nationalpark Schwarzwald, Zuflucht" |
           Locality == "Nationalpark Schwarzwald, WDG" ) %>%  #other subplots in Germany
  dplyr::select(Code, Locality, Municipality, Easting, Northing, Elev) %>% 
  group_by(Code) %>% 
  summarise(Easting_avg = mean(Easting, na.rm = T),
            Northing_avg = mean (Northing, na.rm = T),
            Elevation_avg = round(mean (Elev, na.rm = T), 1),
            Elevation_min = round(min(Elev, na.rm = T), 1),
            Elevation_max = round(max(Elev, na.rm = T), 1)) %>% 
  arrange(desc(Northing_avg)) %>% 
  mutate(Subplot = "Others in Germany")
GeoCoord_detail_BF3$Subplot_Nr <- "B/V3"

# "Middle" for all region with two digits
GeoCoord <- DF_all %>% 
  dplyr::select(Code, Easting, Northing, Elev) %>% 
  group_by(Code) %>% 
  summarise(Easting_avg = mean(Easting, na.rm = T),
            Northing_avg = mean(Northing, na.rm = T),
            Elevation_avg = round(mean (Elev, na.rm = T), 1),
            Elevation_min = round(min(Elev, na.rm = T), 1),
            Elevation_max = round(max(Elev, na.rm = T), 1)) %>% 
  arrange(desc(Northing_avg))

GeoCoord <- mutate(GeoCoord, Subplot = "Middle")
GeoCoord$Subplot_Nr <- c("main", "main", "main", "main", "main", "main", "main", "main", "main", "main", "main")

GeoCoord_detail <- rbind(GeoCoord, 
                         GeoCoord_detail_1, GeoCoord_detail_2, 
                         GeoCoord_detail_BF1, GeoCoord_detail_BF2, GeoCoord_detail_BF3)
GeoCoord_detail <- arrange(GeoCoord_detail, desc(Northing_avg))

print(GeoCoord_detail)
saveRDS(GeoCoord_detail, "../_data/_processed_data/GeoCoord_avg_11regions_detail.rds")

## 4. Endemic/Shared RT/RG
source("utils_ecology_myxo.R")
load_ecology_packages()

### 4-1. load data
community_RT <- readRDS("../_data/_processed_data/community_RT.rds")
community_RG <- readRDS("../_data/_processed_data/community_RG.rds")

### 4-2. convert to dataframe
DF_community_RT <- community_RT %>%
  as.data.frame() %>%
  rownames_to_column(var = "Code")

DF_community_RG <- community_RG %>%
  as.data.frame() %>%
  rownames_to_column(var = "Code")

### 4-3. calculate
DF_long_RT <- DF_community_RT %>%
  pivot_longer(
    cols = -Code,           
    names_to = "RT",   
    values_to = "Count") %>%
  filter(Count > 0) 

DF_long_RG <- DF_community_RG %>%
  pivot_longer(
    cols = -Code,          
    names_to = "RG",   
    values_to = "Count") %>%
  filter(Count > 0) 


DF_RT_summary <- DF_long_RT %>%
  group_by(RT) %>%
  summarise(
    Num_Regions = n_distinct(Code), 
    .groups = 'drop')

DF_RG_summary <- DF_long_RG %>%
  group_by(RG) %>%
  summarise(
    Num_Regions = n_distinct(Code), 
    .groups = 'drop')

# combine with riginal data
DF_classified_RT <- DF_long_RT %>%
  left_join(DF_RT_summary, by = "RT") %>%
  mutate(
    Classification = case_when(
      Num_Regions == 1 ~ "Endemic",
      Num_Regions > 1  ~ "Shared",
      TRUE ~ NA_character_ 
    )
  )

DF_classified_RG <- DF_long_RG %>%
  left_join(DF_RG_summary, by = "RG") %>%
  mutate(
    Classification = case_when(
      Num_Regions == 1 ~ "Endemic",
      Num_Regions > 1  ~ "Shared",
      TRUE ~ NA_character_ 
    )
  )

### 4-4. Calculate endemic/shared RT/R by region
DF_result_RT <- DF_classified_RT %>%
  group_by(Code, Classification) %>%
  summarise(
    Total_Species_Count = n_distinct(RT), 
    Total_Occurrence_Count = sum(Count),      
    .groups = 'drop_last'
  ) %>%
  # convert to wide data
  pivot_wider(
    names_from = Classification,
    values_from = c(Total_Species_Count, Total_Occurrence_Count),
    values_fill = 0 
  )

DF_result_RG <- DF_classified_RG %>%
  group_by(Code, Classification) %>%
  summarise(
    Total_Species_Count = n_distinct(RG), 
    Total_Occurrence_Count = sum(Count),      
    .groups = 'drop_last'
  ) %>%
  # conveRG to wide data
  pivot_wider(
    names_from = Classification,
    values_from = c(Total_Species_Count, Total_Occurrence_Count),
    values_fill = 0 
  )
colnames(DF_result_RT)[colnames(DF_result_RT) == "Total_Occurrence_Count_Endemic"] <- "endemic"
colnames(DF_result_RT)[colnames(DF_result_RT) == "Total_Occurrence_Count_Shared"] <- "shared"
colnames(DF_result_RT)[colnames(DF_result_RT) == "Total_Species_Count_Endemic"] <- "endemic_RT"
colnames(DF_result_RT)[colnames(DF_result_RT) == "Total_Species_Count_Shared"] <- "shared_RT"


colnames(DF_result_RG)[colnames(DF_result_RG) == "Total_Occurrence_Count_Endemic"] <- "endemic"
colnames(DF_result_RG)[colnames(DF_result_RG) == "Total_Occurrence_Count_Shared"] <- "shared"
colnames(DF_result_RG)[colnames(DF_result_RG) == "Total_Species_Count_Endemic"] <- "endemic_RG"
colnames(DF_result_RG)[colnames(DF_result_RG) == "Total_Species_Count_Shared"] <- "shared_RG"


saveRDS(DF_result_RT, "../_data/_processed_data/endemic_shared_RT.rds")
saveRDS(DF_result_RG, "../_data/_processed_data/endemic_shared_RG.rds")

## 5. Environmental variables

source("utils_ecology_myxo.R")
load_ecology_packages()

library(raster)
library(ncdf4)
#library(rgdal)

### 5-1. Load data
bio1 <- raster("../_data/wc2_30s/wc2.1_30s_bio_1.tif")
bio2 <- raster("../_data/wc2_30s/wc2.1_30s_bio_2.tif")
bio3 <- raster("../_data/wc2_30s/wc2.1_30s_bio_3.tif")
bio4 <- raster("../_data/wc2_30s/wc2.1_30s_bio_4.tif")
bio5 <- raster("../_data/wc2_30s/wc2.1_30s_bio_5.tif")
bio6 <- raster("../_data/wc2_30s/wc2.1_30s_bio_6.tif")
bio7 <- raster("../_data/wc2_30s/wc2.1_30s_bio_7.tif")
bio8 <- raster("../_data/wc2_30s/wc2.1_30s_bio_8.tif")
bio9 <- raster("../_data/wc2_30s/wc2.1_30s_bio_9.tif")
bio10 <- raster("../_data/wc2_30s/wc2.1_30s_bio_10.tif")
bio11 <- raster("../_data/wc2_30s/wc2.1_30s_bio_11.tif")
bio12 <- raster("../_data/wc2_30s/wc2.1_30s_bio_12.tif")
bio13 <- raster("../_data/wc2_30s/wc2.1_30s_bio_13.tif")
bio14 <- raster("../_data/wc2_30s/wc2.1_30s_bio_14.tif")
bio15 <- raster("../_data/wc2_30s/wc2.1_30s_bio_15.tif")
bio16 <- raster("../_data/wc2_30s/wc2.1_30s_bio_16.tif")
bio17 <- raster("../_data/wc2_30s/wc2.1_30s_bio_17.tif")
bio18 <- raster("../_data/wc2_30s/wc2.1_30s_bio_18.tif")
bio19 <- raster("../_data/wc2_30s/wc2.1_30s_bio_19.tif")

GeoCoord_detail <- readRDS("../_data/_processed_data/GeoCoord_avg_11regions_detail.rds") # sequenced data
GeoCoord_bio <- dplyr::select(GeoCoord_detail, Code, Subplot_Nr, Subplot, Easting_avg, Northing_avg)

### 5-2. Extract BioClim data and add as a new column
v1_Mean_temp <- extract(bio1, cbind(GeoCoord_bio$Easting_avg, GeoCoord_bio$Northing_avg))
GeoCoord_bio$v1_Mean_temp <- v1_Mean_temp

v2_Mean_diurnal_range <- extract(bio2, cbind(GeoCoord_bio$Easting_avg, GeoCoord_bio$Northing_avg))
GeoCoord_bio$v2_Mean_diurnal_range <- v2_Mean_diurnal_range

v3_Isothermality <- extract(bio3, cbind(GeoCoord_bio$Easting_avg, GeoCoord_bio$Northing_avg))
GeoCoord_bio$v3_Isothermality <- v3_Isothermality

v4_Temp_Seasonality <- extract(bio4, cbind(GeoCoord_bio$Easting_avg, GeoCoord_bio$Northing_avg))
GeoCoord_bio$v4_Temp_Seasonality <- v4_Temp_Seasonality

v5_Max_Temp_Warmest_Month <- extract(bio5, cbind(GeoCoord_bio$Easting_avg, GeoCoord_bio$Northing_avg))
GeoCoord_bio$v5_Max_Temp_Warmest_Month <- v5_Max_Temp_Warmest_Month

v6_Min_Temp_Coldest_Month <- extract(bio6, cbind(GeoCoord_bio$Easting_avg, GeoCoord_bio$Northing_avg))
GeoCoord_bio$v6_Min_Temp_Coldest_Month <- v6_Min_Temp_Coldest_Month

v7_Temp_Annual_Range <- extract(bio7, cbind(GeoCoord_bio$Easting_avg, GeoCoord_bio$Northing_avg))
GeoCoord_bio$v7_Temp_Annual_Range <- v7_Temp_Annual_Range

v8_Mean_Temp_Wettest_Quarter <- extract(bio8, cbind(GeoCoord_bio$Easting_avg, GeoCoord_bio$Northing_avg))
GeoCoord_bio$v8_Mean_Temp_Wettest_Quarter <- v8_Mean_Temp_Wettest_Quarter

v9_Mean_Temp_Driest_Quarter <- extract(bio9, cbind(GeoCoord_bio$Easting_avg, GeoCoord_bio$Northing_avg))
GeoCoord_bio$v9_Mean_Temp_Driest_Quarter <- v9_Mean_Temp_Driest_Quarter

v10_Mean_Temp_Warmest_Quarter <- extract(bio10, cbind(GeoCoord_bio$Easting_avg, GeoCoord_bio$Northing_avg))
GeoCoord_bio$v10_Mean_Temp_Warmest_Quarter <- v10_Mean_Temp_Warmest_Quarter

v11_Mean_Temp_Coldest_Quarter <- extract(bio11, cbind(GeoCoord_bio$Easting_avg, GeoCoord_bio$Northing_avg))
GeoCoord_bio$v11_Mean_Temp_Coldest_Quarter <- v11_Mean_Temp_Coldest_Quarter

v12_Annual_Precipit <- extract(bio12, cbind(GeoCoord_bio$Easting_avg, GeoCoord_bio$Northing_avg))
GeoCoord_bio$v12_Annual_Precipit <- v12_Annual_Precipit

v13_Precipit_Wettest_Month <- extract(bio13, cbind(GeoCoord_bio$Easting_avg, GeoCoord_bio$Northing_avg))
GeoCoord_bio$v13_Precipit_Wettest_Month <- v13_Precipit_Wettest_Month

v14_Precipit_Driest_Month <- extract(bio14, cbind(GeoCoord_bio$Easting_avg, GeoCoord_bio$Northing_avg))
GeoCoord_bio$v14_Precipit_Driest_Month <- v14_Precipit_Driest_Month

v15_Precipit_Seasonality <- extract(bio15, cbind(GeoCoord_bio$Easting_avg, GeoCoord_bio$Northing_avg))
GeoCoord_bio$v15_Precipit_Seasonality <- v15_Precipit_Seasonality

v16_Precipit_Wettest_Quarter <- extract(bio16, cbind(GeoCoord_bio$Easting_avg, GeoCoord_bio$Northing_avg))
GeoCoord_bio$v16_Precipit_Wettest_Quarter <- v16_Precipit_Wettest_Quarter

v17_Precipit_Driest_Quarter <- extract(bio17, cbind(GeoCoord_bio$Easting_avg, GeoCoord_bio$Northing_avg))
GeoCoord_bio$v17_Precipit_Driest_Quarter <- v17_Precipit_Driest_Quarter

v18_Precipit_Warmest_Quarter <- extract(bio18, cbind(GeoCoord_bio$Easting_avg, GeoCoord_bio$Northing_avg))
GeoCoord_bio$v18_Precipit_Warmest_Quarter <- v18_Precipit_Warmest_Quarter

v19_Precipit_Coldest_Quarter <- extract(bio19, cbind(GeoCoord_bio$Easting_avg, GeoCoord_bio$Northing_avg))
GeoCoord_bio$v19_Precipit_Coldest_Quarter <- v19_Precipit_Coldest_Quarter



### 5-3. Extract snow data
#### 5-3-1. Open the NC file
Nor_nc <- nc_open("../_data/nc_data/Nor_data.nc")
GAP_nc <- nc_open("../_data/nc_data/GAP_data.nc")
BF_nc <- nc_open("../_data/nc_data/BF_data.nc")
Tat_nc <- nc_open("../_data/nc_data/Tat_data.nc")
Kam_nc <- nc_open("../_data/nc_data/Kam_data.nc")
Jp_nc <- nc_open("../_data/nc_data/Jp_data.nc")
Khi_nc <- nc_open("../_data/nc_data/Khi_data.nc")
Cau_nc <- nc_open("../_data/nc_data/Cau_data.nc")
SN_nc <- nc_open("../_data/nc_data/SN_data.nc")
Fr_nc <- nc_open("../_data/nc_data/Fr_data.nc")
Pyr_nc <- nc_open("../_data/nc_data/Pyr_data.nc")

#### 5-3-2. Obtain longitude/latitude variables in the  file
Nor_lat <- ncvar_get(Nor_nc, "latitude")
Nor_lon <- ncvar_get(Nor_nc, "longitude")
Nor_time <- ncvar_get(Nor_nc, "time")

BF_lat <- ncvar_get(BF_nc, "latitude")
BF_lon <- ncvar_get(BF_nc, "longitude")
BF_time <- ncvar_get(BF_nc, "time")

GAP_lat <- ncvar_get(GAP_nc, "latitude")
GAP_lon <- ncvar_get(GAP_nc, "longitude")
GAP_time <- ncvar_get(GAP_nc, "time")

Tat_lat <- ncvar_get(Tat_nc, "latitude")
Tat_lon <- ncvar_get(Tat_nc, "longitude")
Tat_time <- ncvar_get(Tat_nc, "time")

Kam_lat <- ncvar_get(Kam_nc, "latitude")
Kam_lon <- ncvar_get(Kam_nc, "longitude")
Kam_time <- ncvar_get(Kam_nc, "time")

Jp_lat <- ncvar_get(Jp_nc, "latitude")
Jp_lon <- ncvar_get(Jp_nc, "longitude")
Jp_time <- ncvar_get(Jp_nc, "time")

Khi_lat <- ncvar_get(Khi_nc, "latitude")
Khi_lon <- ncvar_get(Khi_nc, "longitude")
Khi_time <- ncvar_get(Khi_nc, "time")

Cau_lat <- ncvar_get(Cau_nc, "latitude")
Cau_lon <- ncvar_get(Cau_nc, "longitude")
Cau_time <- ncvar_get(Cau_nc, "time")

SN_lat <- ncvar_get(SN_nc, "latitude")
SN_lon <- ncvar_get(SN_nc, "longitude")
SN_time <- ncvar_get(SN_nc, "time")

Fr_lat <- ncvar_get(Fr_nc, "latitude")
Fr_lon <- ncvar_get(Fr_nc, "longitude")
Fr_time <- ncvar_get(Fr_nc, "time")

Pyr_lat <- ncvar_get(Pyr_nc, "latitude")
Pyr_lon <- ncvar_get(Pyr_nc, "longitude")
Pyr_time <- ncvar_get(Pyr_nc, "time")

#### 5-3-3. Extract the data
Nor_snow_cover <- ncvar_get(Nor_nc, "snowc")
GAP_snow_cover <- ncvar_get(GAP_nc, "snowc")
BF_snow_cover <- ncvar_get(BF_nc, "snowc")
Tat_snow_cover <- ncvar_get(Tat_nc, "snowc")
Kam_snow_cover <- ncvar_get(Kam_nc, "snowc")
Jp_snow_cover <- ncvar_get(Jp_nc, "snowc")
Khi_snow_cover <- ncvar_get(Khi_nc, "snowc")
Cau_snow_cover <- ncvar_get(Cau_nc, "snowc")
SN_snow_cover <- ncvar_get(SN_nc, "snowc")
Fr_snow_cover <- ncvar_get(Fr_nc, "snowc")
Pyr_snow_cover <- ncvar_get(Pyr_nc, "snowc")

#### 5-3-3. Close the NC file
nc_close(Nor_nc)
nc_close(GAP_nc)
nc_close(BF_nc)
nc_close(Tat_nc)
nc_close(Kam_nc)
nc_close(Jp_nc)
nc_close(Khi_nc)
nc_close(Cau_nc)
nc_close(SN_nc)
nc_close(Fr_nc)
nc_close(Pyr_nc)

### 5-4. Calculate average snow cover (SC) over the period
dim(Nor_snow_cover) # check dimension the data

# output
# [1]   7   5 120
# 7 x 5 grids x 12 monthes x 10 years
#
# from this result, extract average snow cover data within a grid
# -> reduce dimension 3 (lon, lat, time) to 2 (lon, lat)

Nor_snow_mean <- apply(Nor_snow_cover, c(1, 2), mean, na.rm = TRUE)
dim(Nor_snow_mean)
Nor_snow_mean
GAP_snow_mean <- apply(GAP_snow_cover, c(1, 2), mean, na.rm = TRUE)
BF_snow_mean <- apply(BF_snow_cover, c(1, 2), mean, na.rm = TRUE)
Tat_snow_mean <- apply(Tat_snow_cover, c(1, 2), mean, na.rm = TRUE)
Kam_snow_mean <- apply(Kam_snow_cover, c(1, 2), mean, na.rm = TRUE)
Jp_snow_mean <- apply(Jp_snow_cover, c(1, 2), mean, na.rm = TRUE)
Khi_snow_mean <- apply(Khi_snow_cover, c(1, 2), mean, na.rm = TRUE)
Cau_snow_mean <- apply(Cau_snow_cover, c(1, 2), mean, na.rm = TRUE)
SN_snow_mean <- apply(SN_snow_cover, c(1, 2), mean, na.rm = TRUE)
Fr_snow_mean <- apply(Fr_snow_cover, c(1, 2), mean, na.rm = TRUE)
Pyr_snow_mean <- apply(Pyr_snow_cover, c(1, 2), mean, na.rm = TRUE)

### 5-5. Create tiff file
#### 5-5-1. Create regional data
# Nor
Nor_lon_range <- c(min(Nor_lon), max(Nor_lon))
Nor_lat_range <- c(min(Nor_lat), max(Nor_lat))

# get dimension info
num_lat <- nrow(Nor_snow_mean)  # 7
num_lon <- ncol(Nor_snow_mean)  # 5

# create raster object
Nor_r <- raster(nrow = num_lat, ncol = num_lon) 

# set range
extent(Nor_r) <- c(Nor_lon_range[1], Nor_lon_range[2], Nor_lat_range[1], Nor_lat_range[2])

# set resolution
res_x <- diff(Nor_lon_range) / (num_lon )
res_y <- diff(Nor_lat_range) / (num_lat )
res(Nor_r) <- c(res_x, res_y)
# set value
values(Nor_r) <- as.vector(Nor_snow_mean)
Nor_r

# GAP
GAP_lon_range <- c(min(GAP_lon), max(GAP_lon))
GAP_lat_range <- c(min(GAP_lat), max(GAP_lat))

GAP_num_lat <- nrow(GAP_snow_mean)  # 3
GAP_num_lon <- ncol(GAP_snow_mean)  # 2

GAP_r <- raster(nrow = GAP_num_lat, ncol = GAP_num_lon)

extent(GAP_r) <- c(GAP_lon_range[1], GAP_lon_range[2], GAP_lat_range[1], GAP_lat_range[2])

res_x <- diff(GAP_lon_range) / (GAP_num_lon )
res_y <- diff(GAP_lat_range) / (GAP_num_lat )
res(GAP_r) <- c(res_x, res_y)

values(GAP_r) <- as.vector(GAP_snow_mean)
GAP_r

# B/V
BF_lon_range <- c(min(BF_lon), max(BF_lon))
BF_lat_range <- c(min(BF_lat), max(BF_lat))

BF_num_lat <- nrow(BF_snow_mean)  # 15
BF_num_lon <- ncol(BF_snow_mean)  # 9

BF_r <- raster(nrow = BF_num_lat, ncol = BF_num_lon)

extent(BF_r) <- c(BF_lon_range[1], BF_lon_range[2], BF_lat_range[1], BF_lat_range[2])

res_x <- diff(BF_lon_range) / (BF_num_lon )
res_y <- diff(BF_lat_range) / (BF_num_lat )
res(BF_r) <- c(res_x, res_y)

values(BF_r) <- as.vector(BF_snow_mean)
BF_r

# Tat
Tat_lon_range <- c(min(Tat_lon), max(Tat_lon))
Tat_lat_range <- c(min(Tat_lat), max(Tat_lat))

Tat_num_lat <- nrow(Tat_snow_mean)  # 6
Tat_num_lon <- ncol(Tat_snow_mean)  # 3

Tat_r <- raster(nrow = Tat_num_lat, ncol = Tat_num_lon)

extent(Tat_r) <- c(Tat_lon_range[1], Tat_lon_range[2], Tat_lat_range[1], Tat_lat_range[2])

res_x <- diff(Tat_lon_range) / (Tat_num_lon )
res_y <- diff(Tat_lat_range) / (Tat_num_lat )
res(Tat_r) <- c(res_x, res_y)

values(Tat_r) <- as.vector(Tat_snow_mean)
Tat_r

# Kam
Kam_lon_range <- c(min(Kam_lon), max(Kam_lon))
Kam_lat_range <- c(min(Kam_lat), max(Kam_lat))

Kam_num_lat <- nrow(Kam_snow_mean)  # 8
Kam_num_lon <- ncol(Kam_snow_mean)  # 31

Kam_r <- raster(nrow = Kam_num_lat, ncol = Kam_num_lon)

extent(Kam_r) <- c(Kam_lon_range[1], Kam_lon_range[2], Kam_lat_range[1], Kam_lat_range[2])

res_x <- diff(Kam_lon_range) / (Kam_num_lon )
res_y <- diff(Kam_lat_range) / (Kam_num_lat )
res(Kam_r) <- c(res_x, res_y)

values(Kam_r) <- as.vector(Kam_snow_mean)
Kam_r

# Jp
Jp_lon_range <- c(min(Jp_lon), max(Jp_lon))
Jp_lat_range <- c(min(Jp_lat), max(Jp_lat))

Jp_num_lat <- nrow(Jp_snow_mean)  # 12
Jp_num_lon <- ncol(Jp_snow_mean)  # 10

Jp_r <- raster(nrow = Jp_num_lat, ncol = Jp_num_lon)

extent(Jp_r) <- c(Jp_lon_range[1], Jp_lon_range[2], Jp_lat_range[1], Jp_lat_range[2])

res_x <- diff(Jp_lon_range) / (Jp_num_lon )
res_y <- diff(Jp_lat_range) / (Jp_num_lat )
res(Jp_r) <- c(res_x, res_y)

values(Jp_r) <- as.vector(Jp_snow_mean)
Jp_r

# Khi
Khi_lon_range <- c(min(Khi_lon), max(Khi_lon))
Khi_lat_range <- c(min(Khi_lat), max(Khi_lat))

Khi_num_lat <- nrow(Khi_snow_mean)  # 14
Khi_num_lon <- ncol(Khi_snow_mean)  # 4

Khi_r <- raster(nrow = Khi_num_lat, ncol = Khi_num_lon)

extent(Khi_r) <- c(Khi_lon_range[1], Khi_lon_range[2], Khi_lat_range[1], Khi_lat_range[2])

res_x <- diff(Khi_lon_range) / (Khi_num_lon )
res_y <- diff(Khi_lat_range) / (Khi_num_lat )
res(Khi_r) <- c(res_x, res_y)

values(Khi_r) <- as.vector(Khi_snow_mean)
Khi_r

# Cau
Cau_lon_range <- c(min(Cau_lon), max(Cau_lon))
Cau_lat_range <- c(min(Cau_lat), max(Cau_lat))

Cau_num_lat <- nrow(Cau_snow_mean)  # 4
Cau_num_lon <- ncol(Cau_snow_mean)  # 4

Cau_r <- raster(nrow = Cau_num_lat, ncol = Cau_num_lon)

extent(Cau_r) <- c(Cau_lon_range[1], Cau_lon_range[2], Cau_lat_range[1], Cau_lat_range[2])

res_x <- diff(Cau_lon_range) / (Cau_num_lon )
res_y <- diff(Cau_lat_range) / (Cau_num_lat )
res(Cau_r) <- c(res_x, res_y)

values(Cau_r) <- as.vector(Cau_snow_mean)
Cau_r

# SN
SN_lon_range <- c(min(SN_lon), max(SN_lon))
SN_lat_range <- c(min(SN_lat), max(SN_lat))

SN_num_lat <- nrow(SN_snow_mean)  # 4
SN_num_lon <- ncol(SN_snow_mean)  # 4

SN_r <- raster(nrow = SN_num_lat, ncol = SN_num_lon)

extent(SN_r) <- c(SN_lon_range[1], SN_lon_range[2], SN_lat_range[1], SN_lat_range[2])

res_x <- diff(SN_lon_range) / (SN_num_lon )
res_y <- diff(SN_lat_range) / (SN_num_lat )
res(SN_r) <- c(res_x, res_y)

values(SN_r) <- as.vector(SN_snow_mean)
SN_r

# Fr
Fr_lon_range <- c(min(Fr_lon), max(Fr_lon))
Fr_lat_range <- c(min(Fr_lat), max(Fr_lat))

Fr_num_lat <- nrow(Fr_snow_mean)  # 2
Fr_num_lon <- ncol(Fr_snow_mean)  # 2

Fr_r <- raster(nrow = Fr_num_lat, ncol = Fr_num_lon)

extent(Fr_r) <- c(Fr_lon_range[1], Fr_lon_range[2], Fr_lat_range[1], Fr_lat_range[2])

res_x <- diff(Fr_lon_range) / (Fr_num_lon )
res_y <- diff(Fr_lat_range) / (Fr_num_lat )
res(Fr_r) <- c(res_x, res_y)

values(Fr_r) <- as.vector(Fr_snow_mean)
Fr_r

# Pyr
Pyr_lon_range <- c(min(Pyr_lon), max(Pyr_lon))
Pyr_lat_range <- c(min(Pyr_lat), max(Pyr_lat))

Pyr_num_lat <- nrow(Pyr_snow_mean)  # 3
Pyr_num_lon <- ncol(Pyr_snow_mean)  # 2

Pyr_r <- raster(nrow = Pyr_num_lat, ncol = Pyr_num_lon)

extent(Pyr_r) <- c(Pyr_lon_range[1], Pyr_lon_range[2], Pyr_lat_range[1], Pyr_lat_range[2])

res_x <- diff(Pyr_lon_range) / (Pyr_num_lon )
res_y <- diff(Pyr_lat_range) / (Pyr_num_lat )
res(Pyr_r) <- c(res_x, res_y)

values(Pyr_r) <- as.vector(Pyr_snow_mean)
Pyr_r

#### 5-5-2. Save as tif
writeRaster(Nor_r, "../_data/_processed_data/Nor_snow.tiff", format = "GTiff", overwrite = TRUE)
writeRaster(GAP_r, "../_data/_processed_data/GAP_snow.tiff", format = "GTiff", overwrite = TRUE)
writeRaster(BF_r, "../_data/_processed_data/BF_snow.tiff", format = "GTiff", overwrite = TRUE)
writeRaster(Tat_r, "../_data/_processed_data/Tat_snow.tiff", format = "GTiff", overwrite = TRUE)
writeRaster(Kam_r, "../_data/_processed_data/Kam_snow.tiff", format = "GTiff", overwrite = TRUE)
writeRaster(Jp_r, "../_data/_processed_data/Jp_snow.tiff", format = "GTiff", overwrite = TRUE)
writeRaster(Khi_r, "../_data/_processed_data/Khi_snow.tiff", format = "GTiff", overwrite = TRUE)
writeRaster(Cau_r, "../_data/_processed_data/Cau_snow.tiff", format = "GTiff", overwrite = TRUE)
writeRaster(SN_r, "../_data/_processed_data/SN_snow.tiff", format = "GTiff", overwrite = TRUE)
writeRaster(Fr_r, "../_data/_processed_data/Fr_snow.tiff", format = "GTiff", overwrite = TRUE)
writeRaster(Pyr_r, "../_data/_processed_data/Pyr_snow.tiff", format = "GTiff", overwrite = TRUE)

### 5-6. Calculate snow cover data
#### 5-6-1. Load data
Nor_sc_mean <- raster("../_data/_processed_data/Nor_snow.tiff")
GAP_sc_mean <- raster("../_data/_processed_data/GAP_snow.tiff")
BF_sc_mean <- raster("../_data/_processed_data/BF_snow.tiff")
Tat_sc_mean <- raster("../_data/_processed_data/Tat_snow.tiff")
Kam_sc_mean <- raster("../_data/_processed_data/Kam_snow.tiff")
Jp_sc_mean <- raster("../_data/_processed_data/Jp_snow.tiff")
Khi_sc_mean <- raster("../_data/_processed_data/Khi_snow.tiff")
Cau_sc_mean <- raster("../_data/_processed_data/Cau_snow.tiff")
SN_sc_mean <- raster("../_data/_processed_data/SN_snow.tiff")
Fr_sc_mean <- raster("../_data/_processed_data/Fr_snow.tiff")
Pyr_sc_mean <- raster("../_data/_processed_data/Pyr_snow.tiff")

#### 5-6-2. Create regional coordinate data
Nor_coord <- filter(GeoCoord_bio, Code == "Nor")
GAP_coord <- filter(GeoCoord_bio, Code == "GAP")
BF_coord <- filter(GeoCoord_bio, Code == "B/V")
Tat_coord <- filter(GeoCoord_bio, Code == "Tat")
Kam_coord <- filter(GeoCoord_bio, Code == "Kam")
Jp_coord <- filter(GeoCoord_bio, Code == "Jp")
Khi_coord <- filter(GeoCoord_bio, Code == "Khi")
Cau_coord <- filter(GeoCoord_bio, Code == "Cau")
SN_coord <- filter(GeoCoord_bio, Code == "SN")
Fr_coord <- filter(GeoCoord_bio, Code == "Fr")
SN_coord$Easting360 <- abs(360+SN_coord$Easting_avg) # need to calculate with absolute value
Pyr_coord <- filter(GeoCoord_bio, Code == "Pyr")

#### 5-6-3.Extract SC data and add
Snow <- extract(Nor_sc_mean, cbind(Nor_coord$Easting_avg, Nor_coord$Northing_avg))
Nor_coord$SC <- Snow

Snow <- extract(GAP_sc_mean, cbind(GAP_coord$Easting_avg, GAP_coord$Northing_avg))
GAP_coord$SC <- Snow

Snow <- extract(BF_sc_mean, cbind(BF_coord$Easting_avg, BF_coord$Northing_avg))
BF_coord$SC <- Snow

Snow <- extract(Tat_sc_mean, cbind(Tat_coord$Easting_avg, Tat_coord$Northing_avg))
Tat_coord$SC <- Snow

Snow <- extract(Kam_sc_mean, cbind(Kam_coord$Easting_avg, Kam_coord$Northing_avg))
Kam_coord$SC <- Snow
# to avoid NA, replace by the coordinate within the range
Kam_coord[3, 25] <- extract(Kam_sc_mean, cbind(158.249, 52.966)) 

Snow <- extract(Jp_sc_mean, cbind(Jp_coord$Easting_avg, Jp_coord$Northing_avg))
Jp_coord$SC <- Snow

Snow <- extract(Khi_sc_mean, cbind(Khi_coord$Easting_avg, Khi_coord$Northing_avg))
Khi_coord$SC <- Snow

Snow <- extract(Cau_sc_mean, cbind(Cau_coord$Easting_avg, Cau_coord$Northing_avg))
Cau_coord$SC <- Snow

Snow <- extract(SN_sc_mean, cbind(SN_coord$Easting360, SN_coord$Northing_avg))
SN_coord$SC <- Snow

Snow <- extract(Fr_sc_mean, cbind(Fr_coord$Easting_avg, Fr_coord$Northing_avg))
Fr_coord$SC <- Snow

Snow <- extract(Pyr_sc_mean, cbind(Pyr_coord$Easting_avg, Pyr_coord$Northing_avg))
Pyr_coord$SC <- Snow


#### 5-3-6. Calculate and Combine No Snow Duration to Reional Geog. Coordinate
source("utils_ecology_myxo.R")

# Nor
# Vector for storing results
Nor_no_snow_periods <- numeric(nrow(Nor_coord))

for (i in 1:nrow(Nor_coord)) {
  coord <- Nor_coord[i, ]
  indices <- get_indices(Nor_lon, Nor_lat, coord$Easting_avg, coord$Northing_avg)
  lon_index <- indices[1]
  lat_index <- indices[2]
  
  # Extract snow cover data at each coordinate point
  single_point_data <- Nor_snow_cover[lon_index, lat_index, ]
  
  # Mask created with no snow cover
  single_point_mask <- single_point_data == 0
  
  # Calculate the period of time (in months) during which there is no snow cover
  no_snow_periods <- calculate_consecutive_no_snow_periods(
    array(single_point_mask, dim = c(1, 1, length(single_point_mask)))
  )
  
  # Calculate the average of the duration of snow-free conditions for each coordinate point
  Nor_no_snow_periods[i] <- no_snow_periods[1, 1]
}

Nor_coord$SC_lasting <- 12 - (Nor_no_snow_periods)


# GAP
GAP_no_snow_periods <- numeric(nrow(GAP_coord))

for (i in 1:nrow(GAP_coord)) {
  coord <- GAP_coord[i, ]
  indices <- get_indices(GAP_lon, GAP_lat, coord$Easting_avg, coord$Northing_avg)
  lon_index <- indices[1]
  lat_index <- indices[2]
  
  
  single_point_data <- GAP_snow_cover[lon_index, lat_index, ]
  
  single_point_mask <- single_point_data == 0
  
  no_snow_periods <- calculate_consecutive_no_snow_periods(
    array(single_point_mask, dim = c(1, 1, length(single_point_mask)))
  )
  
  GAP_no_snow_periods[i] <- no_snow_periods[1, 1]
}

GAP_coord$SC_lasting <- 12 - (GAP_no_snow_periods)

# BF
BF_no_snow_periods <- numeric(nrow(BF_coord))

for (i in 1:nrow(BF_coord)) {
  coord <- BF_coord[i, ]
  indices <- get_indices(BF_lon, BF_lat, coord$Easting_avg, coord$Northing_avg)
  lon_index <- indices[1]
  lat_index <- indices[2]
  
  single_point_data <- BF_snow_cover[lon_index, lat_index, ]
  
  single_point_mask <- single_point_data == 0
  
  no_snow_periods <- calculate_consecutive_no_snow_periods(
    array(single_point_mask, dim = c(1, 1, length(single_point_mask)))
  )
  
  BF_no_snow_periods[i] <- no_snow_periods[1, 1]
}

BF_coord$SC_lasting <- 12 - (BF_no_snow_periods)

# Tat
Tat_no_snow_periods <- numeric(nrow(Tat_coord))

for (i in 1:nrow(Tat_coord)) {
  coord <- Tat_coord[i, ]
  indices <- get_indices(Tat_lon, Tat_lat, coord$Easting_avg, coord$Northing_avg)
  lon_index <- indices[1]
  lat_index <- indices[2]
  
  single_point_data <- Tat_snow_cover[lon_index, lat_index, ]
  
  single_point_mask <- single_point_data == 0
  
  no_snow_periods <- calculate_consecutive_no_snow_periods(
    array(single_point_mask, dim = c(1, 1, length(single_point_mask)))
  )
  
  Tat_no_snow_periods[i] <- no_snow_periods[1, 1]
}

Tat_coord$SC_lasting <- 12 - (Tat_no_snow_periods)

# Kam
Kam_no_snow_periods <- numeric(nrow(Kam_coord))

for (i in 1:nrow(Kam_coord)) {
  coord <- Kam_coord[i, ]
  indices <- get_indices(Kam_lon, Kam_lat, coord$Easting_avg, coord$Northing_avg)
  lon_index <- indices[1]
  lat_index <- indices[2]
  
  single_point_data <- Kam_snow_cover[lon_index, lat_index, ]
  
  single_point_mask <- single_point_data == 0
  
  no_snow_periods <- calculate_consecutive_no_snow_periods(
    array(single_point_mask, dim = c(1, 1, length(single_point_mask)))
  )
  
  Kam_no_snow_periods[i] <- no_snow_periods[1, 1]
}

Kam_coord$SC_lasting <- 12 - (Kam_no_snow_periods)

# Jp
Jp_no_snow_periods <- numeric(nrow(Jp_coord))

for (i in 1:nrow(Jp_coord)) {
  coord <- Jp_coord[i, ]
  indices <- get_indices(Jp_lon, Jp_lat, coord$Easting_avg, coord$Northing_avg)
  lon_index <- indices[1]
  lat_index <- indices[2]
  
  single_point_data <- Jp_snow_cover[lon_index, lat_index, ]
  
  single_point_mask <- single_point_data == 0
  
  no_snow_periods <- calculate_consecutive_no_snow_periods(
    array(single_point_mask, dim = c(1, 1, length(single_point_mask)))
  )
  
  Jp_no_snow_periods[i] <- no_snow_periods[1, 1]
}

Jp_coord$SC_lasting <- 12 - (Jp_no_snow_periods)

# Khi
Khi_no_snow_periods <- numeric(nrow(Khi_coord))

for (i in 1:nrow(Khi_coord)) {
  coord <- Khi_coord[i, ]
  indices <- get_indices(Khi_lon, Khi_lat, coord$Easting_avg, coord$Northing_avg)
  lon_index <- indices[1]
  lat_index <- indices[2]
  
  single_point_data <- Khi_snow_cover[lon_index, lat_index, ]
  
  single_point_mask <- single_point_data == 0
  
  no_snow_periods <- calculate_consecutive_no_snow_periods(
    array(single_point_mask, dim = c(1, 1, length(single_point_mask)))
  )
  
  Khi_no_snow_periods[i] <- no_snow_periods[1, 1]
}

Khi_coord$SC_lasting <- 12 -(Khi_no_snow_periods)

# Cau
Cau_no_snow_periods <- numeric(nrow(Cau_coord))

for (i in 1:nrow(Cau_coord)) {
  coord <- Cau_coord[i, ]
  indices <- get_indices(Cau_lon, Cau_lat, coord$Easting_avg, coord$Northing_avg)
  lon_index <- indices[1]
  lat_index <- indices[2]
  
  single_point_data <- Cau_snow_cover[lon_index, lat_index, ]
  
  single_point_mask <- single_point_data == 0
  
  no_snow_periods <- calculate_consecutive_no_snow_periods(
    array(single_point_mask, dim = c(1, 1, length(single_point_mask)))
  )
  
  Cau_no_snow_periods[i] <- no_snow_periods[1, 1]
}

Cau_coord$SC_lasting <- 12 - (Cau_no_snow_periods)

#SN
SN_no_snow_periods <- numeric(nrow(SN_coord))

for (i in 1:nrow(SN_coord)) {
  coord <- SN_coord[i, ]
  indices <- get_indices(SN_lon, SN_lat, coord$Easting_avg, coord$Northing_avg)
  lon_index <- indices[1]
  lat_index <- indices[2]
  
  single_point_data <- SN_snow_cover[lon_index, lat_index, ]
  
  single_point_mask <- single_point_data == 0
  
  no_snow_periods <- calculate_consecutive_no_snow_periods(
    array(single_point_mask, dim = c(1, 1, length(single_point_mask)))
  )
  
  SN_no_snow_periods[i] <- no_snow_periods[1, 1]
}

SN_coord$SC_lasting <- 12 -(SN_no_snow_periods)

# Fr
Fr_no_snow_periods <- numeric(nrow(Fr_coord))

for (i in 1:nrow(Fr_coord)) {
  coord <- Fr_coord[i, ]
  indices <- get_indices(Fr_lon, Fr_lat, coord$Easting_avg, coord$Northing_avg)
  lon_index <- indices[1]
  lat_index <- indices[2]
  
  single_point_data <- Fr_snow_cover[lon_index, lat_index, ]
  
  single_point_mask <- single_point_data == 0
  
  no_snow_periods <- calculate_consecutive_no_snow_periods(
    array(single_point_mask, dim = c(1, 1, length(single_point_mask)))
  )
  
  Fr_no_snow_periods[i] <- no_snow_periods[1, 1]
}

Fr_coord$SC_lasting <- 12  # as the risults "-Inf"

#Pyr
Pyr_no_snow_periods <- numeric(nrow(Pyr_coord))

for (i in 1:nrow(Pyr_coord)) {
  coord <- Pyr_coord[i, ]
  indices <- get_indices(Pyr_lon, Pyr_lat, coord$Easting_avg, coord$Northing_avg)
  lon_index <- indices[1]
  lat_index <- indices[2]
  
  single_point_data <- Pyr_snow_cover[lon_index, lat_index, ]
  
  single_point_mask <- single_point_data == 0
  
  no_snow_periods <- calculate_consecutive_no_snow_periods(
    array(single_point_mask, dim = c(1, 1, length(single_point_mask)))
  )
  
  Pyr_no_snow_periods[i] <- no_snow_periods[1, 1]
}

Pyr_coord$SC_lasting <- 12 - (Pyr_no_snow_periods)


###
# Nor
row_to_replace_Nor <- 3 

Nor_coord_main_av <- Nor_coord %>%
  mutate(
    across(
      .cols = 6:26, 
      .fns = ~ {
        mean_val <- mean(.[-row_to_replace_Nor], na.rm = TRUE) 
        replace(., row_to_replace_Nor, mean_val)
      }
    )
  )

# BF
row_to_replace_BF <- 2 

BF_coord_main_av <- BF_coord %>%
  mutate(
    across(
      .cols = 6:26, 
      .fns = ~ {
        mean_val <- mean(.[-row_to_replace_BF], na.rm = TRUE) 
        replace(., row_to_replace_BF, mean_val)
      }
    )
  )

# Kam
row_to_replace_Kam <- 2 

Kam_coord_main_av <- Kam_coord %>%
  mutate(
    across(
      .cols = 6:26, 
      .fns = ~ {
        mean_val <- mean(.[-row_to_replace_Kam], na.rm = TRUE) 
        replace(., row_to_replace_Kam, mean_val)
      }
    )
  )

# Jp
row_to_replace_Jp <- 4 

Jp_coord_main_av <- Jp_coord %>%
  mutate(
    across(
      .cols = 6:26, 
      .fns = ~ {
        mean_val <- mean(.[-row_to_replace_Jp], na.rm = TRUE) 
        replace(., row_to_replace_Jp, mean_val)
      }
    )
  )

# Khi
row_to_replace_Khi <- 2 

Khi_coord_main_av <- Khi_coord %>%
  mutate(
    across(
      .cols = 6:26, 
      .fns = ~ {
        mean_val <- mean(.[-row_to_replace_Khi], na.rm = TRUE) 
        replace(., row_to_replace_Khi, mean_val)
      }
    )
  )


# SN
SN_coord_final <- dplyr:: select(SN_coord, -Easting360)


BioClim_data <- rbind(Khi_coord_main_av, Nor_coord_main_av, Kam_coord_main_av, 
                      Tat_coord, BF_coord_main_av, GAP_coord, Fr_coord, Cau_coord, 
                      Pyr_coord, SN_coord_final, Jp_coord_main_av)

saveRDS(BioClim_data, "../_data/_processed_data/BioClim_data.rds")
write.csv(
  x = BioClim_data,
  file = "../_data/_processed_data/BioClim_data.csv",
  row.names = TRUE # set variables as column name
)
