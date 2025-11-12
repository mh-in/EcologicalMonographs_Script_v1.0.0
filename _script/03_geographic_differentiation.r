# 03. geographic differentiation

source("utils_ecology_myxo.R")
load_ecology_packages()


## 1. Genetic differentiation within common clade: Fig 3ABC
### 1-1. Data preparation

DF_niv_seq <- readRDS("../_data/_processed_data/DF_biogeo_seq.rds") # sequenced data

DF_niv_seq_gen_dif <- DF_niv_seq %>% 
  dplyr::select(Code, Clade, RT) %>% 
  mutate(Clade_RT = paste(RT, Clade, sep = "_") ) 


DF_niv_seq_LAMesp <- DF_niv_seq_gen_dif %>% 
    filter(Clade == "LAMesp") %>% 
  dplyr::select(Code, Clade_RT) %>% 
  group_by(Clade_RT, Code, ) %>% 
  summarise(count = n()) 

DF_niv_seq_DIDalp <- DF_niv_seq_gen_dif %>% 
  filter(Clade == "DIDalp") %>% 
  dplyr::select(Code, Clade_RT) %>% 
  group_by(Clade_RT, Code, ) %>% 
  summarise(count = n()) 

DF_niv_seq_PLScha <- DF_niv_seq_gen_dif %>% 
  filter(Clade == "PLScha") %>% 
  dplyr::select(Code, Clade_RT) %>% 
  group_by(Clade_RT, Code, ) %>% 
  summarise(count = n()) 

### 1-2. Plot
# This plot needs to modify size, legend, location in inkscape on map
library(randomcoloR) # for random colour
library(gridExtra) 

all_clade_rts <- unique(c(DF_niv_seq_LAMesp$Clade_RT, DF_niv_seq_DIDalp$Clade_RT, DF_niv_seq_PLScha$Clade_RT))
set.seed(234)
custom_colors_map <- distinctColorPalette(length(all_clade_rts))
names(custom_colors_map) <- all_clade_rts

create_pie_chart <- function(df) {
  ggplot(df, aes(x = "", y = count, fill = Clade_RT)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar(theta = "y") +
    geom_text(aes(label = count), position = position_stack(vjust = 0.5), size = 5) +
    labs(title = unique(df$Code)) +
    theme_void(base_size = 15) +
    theme(legend.position = "none") + 
    scale_fill_manual(values = custom_colors_map)
}

# dammy data for Clade_RT of all LAMesp, DIDalp, PLScha
dummy_data <- tibble(
  Code = "DUMMY",
  Clade_RT = all_clade_rts,
  count = 1 
)

# create legend plot for dammy data
dummy_plot <- ggplot(dummy_data, aes(x = "", y = count, fill = Clade_RT)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = custom_colors_map) +
  labs(fill = "Clade_RT") + 
  theme(legend.position = "right") 

# extract legend from dammy data
get_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if (length(leg) > 0) return(tmp$grobs[[leg]]) else return(NULL)
}

# extract legend object
legend_grob <- get_legend(dummy_plot)


# LAMesp
plots_LAMesp_df <- DF_niv_seq_LAMesp %>%
  group_by(Code) %>%
  do(plot = create_pie_chart(.))

plots_list <- plots_LAMesp_df$plot

png("../_results/Fig3A_gen_dif_LAMesp.png", width = 600, height = 500)
pdf("../_results/Fig3A_gen_dif_LAMesp.pdf", width = 7, height = 8) #for paper
grid.arrange(
  grobs = c(lapply(plots_list, ggplotGrob), list(legend_grob)),
  ncol = 4 
)
dev.off()

# DIDalp
plots_DIDalp_df <- DF_niv_seq_DIDalp %>%
  group_by(Code) %>%
  do(plot = create_pie_chart(.))

plots_list <- plots_DIDalp_df$plot

png("../_results/Fig3B_gen_dif_DIDalp.png", width = 600, height = 500)
pdf("../_results/Fig3B_gen_dif_DIDalp.pdf", width = 7, height = 8) #for paper
grid.arrange(
  grobs = c(lapply(plots_list, ggplotGrob), list(legend_grob)),
  ncol = 4 
)
dev.off()

# PLScha
plots_PLScha_df <- DF_niv_seq_PLScha %>%
  group_by(Code) %>%
  do(plot = create_pie_chart(.))

plots_list <- plots_PLScha_df$plot

png("../_results/Fig3C_gen_dif_PLScha.png", width = 600, height = 500)
pdf("../_results/Fig3C_gen_dif_PLScha.pdf", width = 7, height = 8) #for paper
grid.arrange(
  grobs = c(lapply(plots_list, ggplotGrob), list(legend_grob)),
  ncol = 4 
)
dev.off()

## 2. Endemic/ Shared RT for Figure 4
### 2-1. Load data
source("utils_ecology_myxo.R")
load_ecology_packages()

### 2-2. Data preparation
DF_endemic_shared_RT <- readRDS("../_data/_processed_data/endemic_shared_RT.rds")
DF_endemic_shared_RG <- readRDS("../_data/_processed_data/endemic_shared_RG.rds")

DF_endemic_shared_RT <- DF_endemic_shared_RT[, -(2:3)]
DF_endemic_shared_RG <- DF_endemic_shared_RG[, -(2:3)]

endemic_shared_RT_long <- DF_endemic_shared_RT %>% 
  pivot_longer(cols = endemic:shared,
               names_to = "class",
               values_to = "count")

endemic_shared_RG_long <- DF_endemic_shared_RG %>% 
  pivot_longer(cols = endemic:shared,
               names_to = "class",
               values_to = "count")

### 2-4. Plot
# This plot needs to modify legend colour in inkscape

library(webr)

png("../_results/Fig4A_pie_endemic_RT.png", width = 500, height = 500)
pdf("../_results/Fig4A_pie_endemic_RT.pdf", width = 7, height = 7) #for paper
PieDonut(endemic_shared_RT_long, aes(class, Code, count=count), title = "",  
         titlesize = 8, 
         pieLabelSize = 8, 
         donutLabelSize = 7, 
         explode = 1, 
         explodePos= 0.1
)
dev.off()

png("../_results/Fig4B_pie_endemic_region_RT.png", width = 500, height = 500)
pdf("../_results/Fig4B_pie_endemic_region_RT.pdf", width = 7, height = 7) #for paper
PieDonut(endemic_shared_RT_long, aes(Code, class, count=count), title = "",  
         titlesize = 8, 
         pieLabelSize = 8, 
         donutLabelSize = 7, 
         explode = c(5, 10), 
         explodePos= 0.2
)
dev.off()

png("../_results/Fig4C_pie_endemic_RG.png", width = 500, height = 500)
pdf("../_results/Fig4C_pie_endemic_RG.pdf", width = 7, height = 7) #for paper
PieDonut(endemic_shared_RG_long, aes(class, Code, count=count), title = "",  
         titlesize = 8, 
         pieLabelSize = 8, 
         donutLabelSize = 7, 
         explode = 1, 
         explodePos= 0.1
)
dev.off()

png("../_results/Fig4D_pie_endemic_region_RG.png", width = 500, height = 500)
pdf("../_results/Fig4D_pie_endemic_region_RG.pdf", width = 7, height = 7) #for paper
PieDonut(endemic_shared_RG_long, aes(Code, class, count=count), title = "",  
         titlesize = 8, 
         pieLabelSize = 8, 
         donutLabelSize = 7, 
         explode = c(5, 10), 
         explodePos= 0.2
)
dev.off()
