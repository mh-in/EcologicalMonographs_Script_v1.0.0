# utils for ecological analysis myxomycetes

#' Required package
#'
#' @export
load_ecology_packages <- function() {
  required_packages <- c("dplyr", "tidyr", "reshape2", 
                         "proxy", 
                         "tibble", 
                         "rlang", "ggplot2", "vegan", "iNEXT", "patchwork")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    # Indicated library call (for renv detection)
    eval(parse(text = paste0("library(", pkg, ")")))
  }
}

#' Convert Region Codes to Names in Data Frame
#'
#' @param data Data frame containing region codes
#' @param code_column Name of the column containing region codes
#' @param mapping Named vector for code to name mapping
#' @param new_column_name Optional new column name (if different from original)
#'
#' @return Data frame with converted region names
#'
#' @export
convert_region_codes <- function(data, code_column = NULL, mapping = NULL, new_column_name = NULL) {
  
  # regional name mapping
  if (is.null(mapping)) {
    mapping <- c(
      "#1" = "Nor", "#2" = "GAP", "#3" = "B/V", "#4" = "Tat", 
      "#5" = "Kam", "#6" = "Jp", "#7" = "Khi", "#8" = "Cau",
      "#9" = "SN", "#10" = "Fr", "#11" = "Pyr"
    )
  }
  
  # regional code column auto-detection (when not indicated) 
  if (is.null(code_column)) {
    # search for numeric column
    potential_columns <- sapply(data, function(x) {
      if (is.numeric(x)) return(TRUE)
      if (is.character(x)) return(all(grepl("^[0-9]+$", na.omit(x))))
      return(FALSE)
    })
    
    if (any(potential_columns)) {
      code_column <- names(which(potential_columns))[1]
      message("Auto-detected code column: ", code_column)
    } else {
      stop("No suitable code column found. Please specify code_column parameter.")
    }
  }
  
  # converting
  original_values <- as.character(data[[code_column]])
  converted_values <- mapping[original_values]
  
  # new column name
  output_column <- if (!is.null(new_column_name)) new_column_name else code_column
  
  # add/update dataframe
  data[[output_column]] <- ifelse(is.na(converted_values), original_values, converted_values)
  
  # show summary after converting
  conversion_stats <- table(Original = original_values, Converted = data[[output_column]])
  message("Region conversion summary:")
  print(conversion_stats)
  
  return(data)
}

#' Create community matrix
#'
#' @description
#' Create species (OTU) x sampling site community matrix
#'
#' @param data raw data
#' @param site_col column for sampling site (character)
#' @param species_col column for OTU (character)
#' @param value_col column for abundance/presence data (character): default "abundance"
#' @param fun.aggregate Function; How to calculate duplicated species x sampling site: default = "sum"
#'
#' @return sampling site × species (OTU) community matrix
#'
#' @export
create_community_matrix <- function(data, site_col, species_col, 
                                    value_col = NULL, fun.aggregate = sum) {
  
  load_ecology_packages()
  
  # if "value_col" not indicated default "abundance"
  if (is.null(value_col)) {
    data <- data %>%
      dplyr::group_by(!!rlang::sym(site_col), !!rlang::sym(species_col)) %>%
      dplyr::summarise(abundance = dplyr::n(), .groups = "drop")
    value_col <- "abundance"
  } else {
    # if indicated sum
    data <- data %>%
      dplyr::group_by(!!rlang::sym(site_col), !!rlang::sym(species_col)) %>%
      dplyr::summarise(abundance = fun.aggregate(!!rlang::sym(value_col)), .groups = "drop")
  }
  
  # Create community matrix
  comm_matrix <- data %>%
    tidyr::pivot_wider(
      names_from = !!rlang::sym(species_col),
      values_from = abundance,
      values_fill = 0
    ) %>%
    tibble::column_to_rownames(var = site_col)
  
  return(as.matrix(comm_matrix))
}


# Beta Diversity Indices for Pairwise Regional Comparisons
#
# This script provides functions to calculate three beta diversity indices
# for pairwise comparisons between multiple regions.
#
# Functions:
#   - beta_d1(): Simple count of shared RTs/RGs between two regions
#   - beta_d2(): Adjusted proportional Jaccard-Dice Coefficient
#   - beta_d3(): Proportion of shared RTs/RGs weighted by specimen count
#   - calculate_beta_diversity(): Main function for pairwise comparisons
#
# Usage:
#   data_matrix <- matrix(c(...), nrow = species, ncol = regions)
#   results <- calculate_beta_diversity(data_matrix)
#

#' Calculate Beta Diversity D1 Index
#'
#' Computes the simple count of RTs/RGs shared between two regions.
#'
#' @param a Numeric vector of specimen counts for region A
#' @param b Numeric vector of specimen counts for region B
#' @return Integer value representing the number of shared RTs/RGs
#' @export
beta_d1 <- function(a, b) {
  if (length(a) != length(b)) {
    stop("Vectors a and b must have the same length")
  }
  sum(a > 0 & b > 0)
}

#' Calculate Beta Diversity D2 Index
#'
#' Computes the adjusted proportional Jaccard-Dice Coefficient.
#' Returns value as percentage (0-100%).
#'
#' @param a Numeric vector of specimen counts for region A
#' @param b Numeric vector of specimen counts for region B
#' @return Numeric value representing similarity percentage (0-100%)
#' @export
#'
#' @examples
#' region_A <- c(5, 0, 3, 2, 0, 1)
#' region_B <- c(2, 4, 0, 1, 0, 3)
#' beta_d2(region_A, region_B)
beta_d2 <- function(a, b) {
  if (length(a) != length(b)) {
    stop("Vectors a and b must have the same length")
  }
  
  shared <- sum(a > 0 & b > 0)
  total_a <- sum(a > 0)
  total_b <- sum(b > 0)
  
  if (total_a + total_b == 0) {
    return(0)
  } else {
    return((shared * 2) / (total_a + total_b) * 100)
  }
}

#' Calculate Beta Diversity D3 Index
#'
#' Computes the proportion of shared RTs/RGs weighted by specimen count.
#' Returns value as percentage (0-100%).
#'
#' @param a Numeric vector of specimen counts for region A
#' @param b Numeric vector of specimen counts for region B
#' @return Numeric value representing weighted similarity percentage (0-100%)
#' @export
#'
#' @examples
#' region_A <- c(5, 0, 3, 2, 0, 1)
#' region_B <- c(2, 4, 0, 1, 0, 3)
#' beta_d3(region_A, region_B)
beta_d3 <- function(a, b) {
  if (length(a) != length(b)) {
    stop("Vectors a and b must have the same length")
  }
  
  total_sum <- sum(a + b)
  shared_sum <- sum((a + b)[a > 0 & b > 0])
  
  if (total_sum == 0) {
    return(0)
  } else {
    return(shared_sum / total_sum * 100)
  }
}

#' Calculate pairwise beta diversity indices (similarity and dissimilarity: D1-D3)
#'
#' @description
#' Calculate three pairwise beta diversity indices (D1, D2, D3) between multiple regions/sites.
#' Returns both similarity and dissimilarity versions for use with vegan package.
#'
#' @param comm_matrix Community matrix with sites as rows and species as columns
#' @param method Which indices to calculate: "all" or specific indices ("D1", "D2", "D3")
#'
#' @return A data frame with pairwise comparisons and calculated indices (both similarity and dissimilarity)
#'
#' @export
calculate_beta_diversity1 <- function(comm_matrix, method = "all") {
  
  # Input validation
  if (!is.matrix(comm_matrix) && !is.data.frame(comm_matrix)) {
    stop("comm_matrix must be a matrix or data frame")
  }
  
  if (nrow(comm_matrix) < 2) {
    stop("comm_matrix must have at least 2 rows (sites)")
  }
  
  valid_methods <- c("all", "D1", "D2", "D3")
  if (!method %in% valid_methods) {
    stop("method must be one of: ", paste(valid_methods, collapse = ", "))
  }
  
  # Convert to matrix if data frame
  comm_matrix <- as.matrix(comm_matrix)
  
  # Get site names
  site_names <- rownames(comm_matrix)
  if (is.null(site_names)) {
    site_names <- paste0("Site_", 1:nrow(comm_matrix))
  }
  
  # Pre-calculate max D1 for dissimilarity conversion
  max_d1 <- .calculate_max_D1(comm_matrix)
  
  # Initialize results
  results <- list()
  pair_count <- 1
  
  # Calculate pairwise comparisons
  for (i in 1:(nrow(comm_matrix) - 1)) {
    for (j in (i + 1):nrow(comm_matrix)) {
      
      site_i <- comm_matrix[i, ]
      site_j <- comm_matrix[j, ]
      
      # Calculate indices based on method parameter
      if (method == "all" || method == "D1") {
        d1_sim <- .calculate_D1(site_i, site_j)
        d1_dissim <- max_d1 - d1_sim  # Convert to dissimilarity: max(D1) - D1
      }
      
      if (method == "all" || method == "D2") {
        d2_sim <- .calculate_D2(site_i, site_j)
        d2_dissim <- 100 - d2_sim  # Convert to dissimilarity: 100 - D2
      }
      
      if (method == "all" || method == "D3") {
        d3_sim <- .calculate_D3(site_i, site_j)
        d3_dissim <- 100 - d3_sim  # Convert to dissimilarity: 100 - D3
      }
      
      # Store results
      if (method == "all") {
        results[[pair_count]] <- data.frame(
          site1 = site_names[i],
          site2 = site_names[j],
          D1_similarity = d1_sim,
          D1_dissimilarity = d1_dissim,
          D2_similarity = d2_sim,
          D2_dissimilarity = d2_dissim,
          D3_similarity = d3_sim,
          D3_dissimilarity = d3_dissim,
          stringsAsFactors = FALSE
        )
      } else if (method == "D1") {
        results[[pair_count]] <- data.frame(
          site1 = site_names[i],
          site2 = site_names[j],
          D1_similarity = d1_sim,
          D1_dissimilarity = d1_dissim,
          stringsAsFactors = FALSE
        )
      } else if (method == "D2") {
        results[[pair_count]] <- data.frame(
          site1 = site_names[i],
          site2 = site_names[j],
          D2_similarity = d2_sim,
          D2_dissimilarity = d2_dissim,
          stringsAsFactors = FALSE
        )
      } else if (method == "D3") {
        results[[pair_count]] <- data.frame(
          site1 = site_names[i],
          site2 = site_names[j],
          D3_similarity = d3_sim,
          D3_dissimilarity = d3_dissim,
          stringsAsFactors = FALSE
        )
      }
      
      pair_count <- pair_count + 1
    }
  }
  
  # Combine all results
  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL
  
  return(result_df)
}

# Internal function for D1 index
.calculate_D1 <- function(a, b) {
  sum(a > 0 & b > 0)
}

# Internal function for D2 index
.calculate_D2 <- function(a, b) {
  shared <- sum(a > 0 & b > 0)
  total_a <- sum(a > 0)
  total_b <- sum(b > 0)
  
  if (total_a + total_b == 0) {
    return(0)
  } else {
    return((shared * 2) / (total_a + total_b) * 100)
  }
}

# Internal function for D3 index
.calculate_D3 <- function(a, b) {
  total_sum <- sum(a + b)
  shared_sum <- sum((a + b)[a > 0 & b > 0])
  
  if (total_sum == 0) {
    return(0)
  } else {
    return(shared_sum / total_sum * 100)
  }
}

# Internal function to calculate maximum D1 value
.calculate_max_D1 <- function(comm_matrix) {
  n_sites <- nrow(comm_matrix)
  max_d1 <- 0
  
  for (i in 1:(n_sites - 1)) {
    for (j in (i + 1):n_sites) {
      d1_value <- .calculate_D1(comm_matrix[i, ], comm_matrix[j, ])
      if (d1_value > max_d1) {
        max_d1 <- d1_value
      }
    }
  }
  
  return(max_d1)
}

#' Create beta diversity distance matrices for vegan package: D1-D3
#'
#' @description
#' Convert pairwise beta diversity results to distance matrix format for use with vegan package.
#' Returns dissimilarity matrices as required for Mantel tests and db-RDA.
#'
#' @param beta_results Output from calculate_beta_diversity1()
#' @param index Which index to use for distance matrix ("D1", "D2", "D3")
#'
#' @return A distance matrix object (dissimilarity)
#'
#' @export
beta_to_distmatrix <- function(beta_results, index = "D2") {
  
  valid_indices <- c("D1", "D2", "D3")
  if (!index %in% valid_indices) {
    stop("index must be one of: ", paste(valid_indices, collapse = ", "))
  }
  
  # Get all unique sites
  all_sites <- unique(c(beta_results$site1, beta_results$site2))
  n_sites <- length(all_sites)
  
  # Create empty distance matrix
  dist_mat <- matrix(0, nrow = n_sites, ncol = n_sites,
                     dimnames = list(all_sites, all_sites))
  
  # Fill the matrix with dissimilarity values
  for (i in 1:nrow(beta_results)) {
    site1 <- beta_results$site1[i]
    site2 <- beta_results$site2[i]
    
    # Use pre-calculated dissimilarity values
    if (index == "D1") {
      value <- beta_results$D1_dissimilarity[i]
    } else if (index == "D2") {
      value <- beta_results$D2_dissimilarity[i]
    } else if (index == "D3") {
      value <- beta_results$D3_dissimilarity[i]
    }
    
    dist_mat[site1, site2] <- value
    dist_mat[site2, site1] <- value
  }
  
  # Convert to dist object
  return(as.dist(dist_mat))
}

#' Calculate Multiple Beta Diversity Indices: D5-D6
#'
#' @description
#' Calculate multiple beta diversity indices for community data matrices
#' at different taxonomic levels (e.g., RG, RT levels).
#'
#' @param community_matrix A community data matrix (sites x species)
#' @param level_name Character string describing the level (e.g., "RG", "RT")
#'
#' @return A list containing all calculated distance matrices
#'
#' @examples
#' \dontrun{
#' # Calculate beta diversity for RG level
#' rg_dissimilarity <- calculate_beta_diversity(
#'   community_matrix = DF_RG_div,
#'   level_name = "RG"
#' )
#'
#' # Access individual distance matrices
#' horn_distances <- rg_dissimilarity$horn
#' bray_curtis_distances <- rg_dissimilarity$bray_curtis
#' }
#'
#' @export
#' @importFrom vegan vegdist
calculate_beta_diversity2 <- function(community_matrix, level_name = "") {
  
  # examine dataset
  if (!is.matrix(community_matrix) && !is.data.frame(community_matrix)) {
    stop("community_matrix must be a matrix or data frame")
  }
  
  if (any(community_matrix < 0)) {
    stop("community_matrix contains negative values")
  }
  
  # calculate diversity index
  dissimilarities <- list(
    D5_horn = vegan::vegdist(community_matrix, method = "horn"),
    bray_curtis = vegan::vegdist(community_matrix, method = "bray"),
    D6_chao = vegan::vegdist(community_matrix, method = "chao"),
    D4b_jaccard_unweighted = vegan::vegdist(community_matrix, binary = TRUE, method = "jaccard"),
    D4a_jaccard_weighted = vegan::vegdist(community_matrix, binary = FALSE, method = "jaccard"),
    morisita = vegan::vegdist(community_matrix, method = "morisita")
  )
  
  # calculate Dice distance separately as it not included in vegan
  dissimilarities$dice <- dist(community_matrix, method = "Dice")
  
  # add level name as attribution
  attr(dissimilarities, "level") <- level_name
  attr(dissimilarities, "n_sites") <- nrow(community_matrix)
  attr(dissimilarities, "n_species") <- ncol(community_matrix)
  attr(dissimilarities, "calculation_date") <- Sys.Date()
  
  return(dissimilarities)
}



#' Save Beta Diversity Results
#'
#' @param beta_results Output from calculate_beta_diversity()
#' @param output_dir Directory to save results
#' @param prefix Prefix for output files
#'
#' @export
save_beta_diversity <- function(beta_results, output_dir = "results", prefix = "beta_diversity") {
  
  # create output derectory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  level_name <- attr(beta_results, "level")
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # save results
  saveRDS(beta_results, 
          file.path(output_dir, 
                    paste0(prefix, "_", 
                          level_name, 
#                          "_", 
#                          timestamp, 
                           ".rds")))
  
  # save metadata
  metadata <- list(
    level = attr(beta_results, "level"),
    n_sites = attr(beta_results, "n_sites"),
    n_species = attr(beta_results, "n_species"),
    calculation_date = attr(beta_results, "calculation_date"),
    indices_calculated = names(beta_results)
  )
  
  yaml::write_yaml(metadata, 
                   file.path(output_dir, 
                             paste0(prefix, "_", 
                                    level_name, 
#                                   "_", 
#                                   timestamp, 
                                    "_metadata.yml")))
  
  message("Beta diversity results saved to: ", output_dir)
  return(invisible(TRUE))
}

#' Load Saved Beta Diversity Results
#'
#' @param file_path Path to the saved .rds file
#' @return Beta diversity results with metadata
#'
#' @export
load_beta_diversity <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  results <- readRDS(file_path)
  
  # also load metadata if exist
  metadata_file <- sub("\\.rds$", "_metadata.yml", file_path)
  if (file.exists(metadata_file)) {
    attr(results, "metadata") <- yaml::read_yaml(metadata_file)
  }
  
  return(results)
}

#' Print Summary of Beta Diversity Results
#'
#' @param beta_results Beta diversity results object
#'
#' @export
print_beta_diversity_summary <- function(beta_results) {
  cat("=== Beta Diversity Summary ===\n")
  cat("Level:", attr(beta_results, "level"), "\n")
  cat("Sites:", attr(beta_results, "n_sites"), "\n")
  cat("Species:", attr(beta_results, "n_species"), "\n")
  cat("Indices calculated:", paste(names(beta_results), collapse = ", "), "\n")
  cat("Date:", attr(beta_results, "calculation_date"), "\n\n")
  
  # show basic statictics of each distance matrix
  for (index_name in names(beta_results)) {
    dist_matrix <- beta_results[[index_name]]
    cat("•", index_name, ":\n")
    cat("  Range:", round(range(dist_matrix), 3), "\n")
    cat("  Mean:", round(mean(dist_matrix), 3), "\n")
    cat("  SD:", round(sd(dist_matrix), 3), "\n\n")
  }
}



#' Create community matrix for multiple OTU for iNEXT analysis (transposed community data)
#' 
#' @param data raw data
#' @param site_col column for sampling site
#' @param species_cols column for OTUs e.g., c(RT_Nr, RG, Clade) 
#' @param value_col 
#' @param fun.aggregate 
#' @return list of community
#' @export
create_multi_level_communities <- function(data, site_col, species_cols, 
                                           value_col = NULL, fun.aggregate = sum) {
  
  community_list <- list()
  
  for (species_col in species_cols) {
    cat("Creating community matrix for:", species_col, "\n")
    comm_matrix <- create_community_matrix(
      data = data,
      site_col = site_col,
      species_col = species_col,
      value_col = value_col,
      fun.aggregate = fun.aggregate
    )
    
    # transposed for iNEXT analysis (species x sampling site)
    comm_matrix_t <- t(comm_matrix)
    
    # show raw and column name after transpose
    cat("  Original dim:", dim(comm_matrix), "(sites x species)\n")
    cat("  Transposed dim:", dim(comm_matrix_t), "(species x sites)\n")
    cat("  Number of species:", nrow(comm_matrix_t), "\n")
    cat("  Number of sites:", ncol(comm_matrix_t), "\n")
    cat("  Total occurrences:", sum(comm_matrix_t > 0), "\n\n")
    
    community_list[[species_col]] <- comm_matrix_t
  }
  
  return(community_list)
}

#' Remove last two digits of specimens 
#' @description
#' Remove last two digits of specimens for SAC: Figure 5 
#' e.g., 2000 -> 20
#' Add legend "Specimen (x 100) "
#' 
#' @export
remove_last_three_digits <- function(x) {
  floor(x / 100)
}

#' Create list containing community multilevel community matrix and extract relevant data for all regions
#' @description
#' Create community multilevel community matrix by region and list of iNEXT results with extract relevant values for SAC: Figure 5
#' - observed number of RT
#' - observed number of RG
#' - expected number of RT
#' - expected number of RG
#' - sampling coverage: observed/expected number of RT x 100
#' - sampling coverage: observed/expected number of RG x 100
#' 
#' @export
calculate_region_stats <- function(region_data, region_code) {
  comm_list <- create_multi_level_communities(region_data, "final_colno", c("RT", "RG"))
  set.seed(123)
  res <- plot_combined_sac(comm_list)
  
  # extract basic statistics
  Ob <- res$inext$AsyEst %>% filter(Diversity == "Species richness")
  
  stats <- list(
    obs_RG = Ob[1, 3],
    obs_RT = Ob[2, 3],
    exp_RG = round(Ob[1, 4], 0),
    exp_RT = round(Ob[2, 4], 0),
    cov_RG = Ob[1, 3]/Ob[1, 4]*100,
    cov_RT = Ob[2, 3]/Ob[2, 4]*100,
    inext_result = res$inext,
    plot_data = res$plot
  )
  
  return(stats)
}

#' Plot accumulation curve with relevant values
#'
#' @description
#' Plot accumulation curve with relevant values for SAC: Figure 5
#' Default values are the postion for the curve of "Nor" 
#' 
#' @param stats: created list by function "calculate_region_stats()"
#' 
#' @param x_pos_obs_RG: text x position for observed number of RG
#' @param y_pos_obs_RG: text y position for observed number of RG
#' @param x_pos_obs_RT: text x position for observed number of RT
#' @param y_pos_obs_RT: text y position for observed number of RT
#' @param x_pos_exp: text x position for observed number of RT and RG (same)
#' @param y_pos_exp_RG:  text y position for observed number of RG
#' @param y_pos_exp_RT:  text y position for observed number of RT
#' @param x_pos_cov: text x sampling coverage of RT and RG (same)
#' @param y_pos_cov_RG:  text y position for ampling coverage of RG
#' @param y_pos_cov_RT:  text y position for ampling coverage of RT

create_region_plot <- function(stats, region_code, 
                               x_pos_obs_RG = 500, y_pos_obs_RG = 35,
                               x_pos_obs_RT = 250, y_pos_obs_RT = 110,
                               x_pos_exp = 1700, y_pos_exp_RG = 130, y_pos_exp_RT = 150,
                               x_pos_cov = 1400, y_pos_cov_RT = 200, y_pos_cov_RG = 185,
                               show_y_axis = FALSE, show_x_axis = FALSE, 
                               show_y_text = FALSE, show_x_text = FALSE)  {
  
  p <- ggiNEXT(stats$inext_result, type = 1, color.var="Assemblage") +
    theme_classic(base_size = 15) +
    theme(legend.position="None") +
    scale_y_continuous(expand = c(0, 0), limits=c(0,210)) +
    scale_x_continuous(expand = c(0, 0), limits=c(0,2000), labels = remove_last_three_digits) + 
    labs(x = "Specimen (x100)", y= "RT/RG") +
    annotate("text", x = x_pos_obs_RG, y = y_pos_obs_RG, label = stats$obs_RG, col = "tomato3") +
    annotate("text", x = x_pos_obs_RT, y = y_pos_obs_RT, label = stats$obs_RT, col = "steelblue") +
    annotate("text", x = x_pos_exp, y = y_pos_exp_RG, label = stats$exp_RG, col = "tomato3") +
    annotate("text", x = x_pos_exp, y = y_pos_exp_RT, label = stats$exp_RT, col = "steelblue") +
    annotate("text", x = x_pos_cov, y = y_pos_cov_RT, label = paste(round(stats$cov_RT,1), "%"), 
             col = "steelblue", fontface = "bold") +
    annotate("text", x = x_pos_cov, y = y_pos_cov_RG, label = paste(round(stats$cov_RG,1), "%"), 
             col = "tomato3", fontface = "bold") +
    ggtitle(region_code)
  
  # show or not axis label
  if (!show_y_axis) {
    p <- p + theme(axis.title.y = element_blank())
  }
  if (!show_x_axis) {
    p <- p + theme(axis.title.x = element_blank())
  }
  
  # show or not axis text
  if (!show_y_text) {
    p <- p + theme(axis.text.y = element_blank())
  }
  if (!show_x_text) {
    p <- p + theme(axis.text.x = element_blank())
  }
  
  return(p)
}


#' Plot accumulation curve in different OTU levels together with iNEXT
#' 
#' @param community_list list of community e.g., c(comm_RT, comm_RG, comm_Clade)
#' @param colors colour setting
#' @param labels label setting
#' @param title title
#' @param xlab label x-achse
#' @param ylab label y-achse
#' @return ggplot object
#' @export
plot_combined_sac <- function(community_list, 
                              endpoint_inext = 2000,
                              colors = c("steelblue", "tomato3", "olivedrab", "goldenrod", "purple"),
                              labels = names(community_list),
                              title = "", 
                              xlab = "Number of Specimens (x100)", 
                              ylab = "RT/RG/Clade",
                              limits_y = c(0, 650),
                              limits_x = c(0,2000), 
                              anotate_x = 20,
                              anotate_y = 650,
                              annotate_label = "D") {
  
  load_ecology_packages()
  
  # convert data to list
  data_list <- community_list
  for (i in seq_along(data_list)) {
    rownames(data_list[[i]]) <- rownames(data_list[[i]])
  }
  names(data_list) <- labels[1:length(data_list)]
  
  # create iNEXT object
  inext_obj <- iNEXT(data_list, q = 0, datatype = "incidence_raw", 
                     endpoint = endpoint_inext
                     )
  
  # colour setting
  color_values <- setNames(colors[1:length(data_list)], names(data_list))
  
  # plot
  p <- ggiNEXT(inext_obj, type = 1, facet.var = "None", color.var = "Assemblage") +
    theme_bw() +
    labs(title = title, x = xlab, y = ylab) +
    scale_color_manual(values = color_values) +
    scale_fill_manual(values = color_values) +
    scale_y_continuous(expand = c(0, 0), limits = limits_y) +
    scale_x_continuous(expand = c(0, 0), limits = limits_x, labels = remove_last_three_digits) +
    annotate("text", x = anotate_x, y = anotate_y, label = annotate_label, col = "black", size = 8) +
    theme(legend.position = "bottom", legend.title = element_blank())
  
  return(list(
    plot = p,
    inext = inext_obj
  ))
}



#' Michaelis-Menten kinetics-1
#' @description
#' Function for data prepatation for Michaelis-Menten kinetics
#' 
#' @export
prepare_mm_data <- function(
    data, 
    group_var = "Clade",
    count_var = "RT_Nr",
    x_name = "num_records",
    y_name = "num_ribotypes",
    label_groups = NULL,
    min_records = 1,  # minimum record filter
    min_ribotypes = 1 # minimum RT filter
) {
  
  # calculate basic statistics
  plot_data <- data %>%
    group_by(across(all_of(group_var))) %>%
    summarize(
      !!x_name := n(), # number of records
      !!y_name := n_distinct(!!sym(count_var)) # number of variation
    ) %>%
    ungroup() %>%
    # filter
    filter(!!sym(x_name) >= min_records & !!sym(y_name) >= min_ribotypes)
  
  # add label if indicated 
  if (!is.null(label_groups)) {
    plot_data <- plot_data %>%
      mutate(
        label = ifelse(!!sym(group_var) %in% label_groups, 
                       as.character(!!sym(group_var)), 
                       "")
      )
  } else {
    # add empty label column if not indicated
    plot_data$label <- ""
  }
  
  # show summary
  cat("=== Data Preparation Summary ===\n")
  cat("Total groups:", nrow(plot_data), "\n")
  cat("Records range:", range(plot_data[[x_name]]), "\n")
  cat("Ribotypes range:", range(plot_data[[y_name]]), "\n")
  if (!is.null(label_groups)) {
    cat("Labeled groups:", sum(plot_data$label != ""), "\n")
  }
  
  return(plot_data)
}

#' Michaelis-Menten kinetics-2
#' @description
#' Function for explore optimal initial values
#' 
#' @export
find_optimal_start <- function(plot_data, x_var = "num_records", y_var = "num_ribotypes") {
  
  x_vals <- plot_data[[x_var]]
  y_vals <- plot_data[[y_var]]
  
  # basic statistics
  x_summary <- summary(x_vals)
  y_summary <- summary(y_vals)
  
  cat("=== Data Summary ===\n")
  cat("X variable (records):", x_summary, "\n")
  cat("Y variable (ribotypes):", y_summary, "\n")
  
  # candidates Vmax_start
  Vmax_candidates <- c(
    max(y_vals) * 1.1,      # conservative
    max(y_vals) * 1.5,      # standard
    max(y_vals) * 2.0,      # optimistic
    quantile(y_vals, 0.95) * 1.5  # less influence by outlier value
  )
  
  # candidates Km_start
  Km_candidates <- c(
    median(x_vals),                    # median
    mean(x_vals),                      # average
    quantile(x_vals, 0.6),             # 60%
    quantile(x_vals, 0.75),            # 75%
    max(x_vals) * 0.3                  # 30% of maximum
  )
  
  # test candidate combinations
  results <- list()
  best_r2 <- -Inf
  best_combo <- NULL
  
  for (Vmax_cand in Vmax_candidates) {
    for (Km_cand in Km_candidates) {
      combo_name <- paste0("Vmax_", round(Vmax_cand,1), "_Km_", round(Km_cand,1))
      
      tryCatch({
        model <- nls(
          as.formula(paste(y_var, "~ (Vmax *", x_var, ") / (Km +", x_var, ")")),
          data = plot_data,
          start = list(Vmax = Vmax_cand, Km = Km_cand),
          control = nls.control(maxiter = 100, warnOnly = TRUE)
        )
        
        # calculate the goodness-of-fit
        RSS <- sum(residuals(model)^2)
        TSS <- sum((y_vals - mean(y_vals))^2)
        r2 <- 1 - (RSS / TSS)
        
        results[[combo_name]] <- list(
          Vmax_start = Vmax_cand,
          Km_start = Km_cand,
          r2 = r2,
          converged = TRUE
        )
        
        # update optimal combination
        if (r2 > best_r2) {
          best_r2 <- r2
          best_combo <- combo_name
        }
        
      }, error = function(e) {
        results[[combo_name]] <- list(
          Vmax_start = Vmax_cand,
          Km_start = Km_cand,
          r2 = NA,
          converged = FALSE,
          error = e$message
        )
      })
    }
  }
  
  # show results
  cat("\n=== Initial Value Analysis ===\n")
  successful_models <- sapply(results, function(x) x$converged)
  cat("Successful models:", sum(successful_models), "/", length(results), "\n")
  
  if (!is.null(best_combo)) {
    cat("Best combination:", best_combo, "\n")
    cat("Best R-squared:", round(best_r2, 4), "\n")
    cat("Recommended Vmax_start:", results[[best_combo]]$Vmax_start, "\n")
    cat("Recommended Km_start:", results[[best_combo]]$Km_start, "\n")
  }
  
  return(list(
    all_results = results,
    best_combo = best_combo,
    recommended_Vmax = if(!is.null(best_combo)) results[[best_combo]]$Vmax_start else NA,
    recommended_Km = if(!is.null(best_combo)) results[[best_combo]]$Km_start else NA
  ))
}

#' @examples
#' optimal_starts <- find_optimal_start(plot_data)



#' Michaelis-Menten kinetics-3
#' 
#' @description
#' Function for Michaelis-Menten Model Fit (without saving)
#' Fitting non-linear model
#' y = (Vmax * x) / (Km + x) 
#' Initial value (start) is important (use Michaelis-Menten kinetics-2)
#' - Initial Vmax value: slightly more than maximal num_ribotypes (e.g., max(plot_data$num_ribotypes) * 1.5)
#' - Initial Km value: infer x value from Vmax/2 
#' @export
fit_mm_model <- function(
    data, 
    x_var = "num_records", 
    y_var = "num_ribotypes",
    label_var = "label",
    Vmax_start = NULL,
    Km_start = NULL,
    auto_start = TRUE,
    max_iter = 100,
    make_plot = TRUE,
    plot_title = "",
    x_lab = "Number of records",
    y_lab = "Number of ribotypes",
    point_size = 3,
    line_color = "steelblue", 
    line_size = 1,
    hline_type = "dotted",
    hline_color = "gray50",
    text_adjust = TRUE,
    theme_style = "minimal",
    show_labels = TRUE,
    label_size = 3,
    label_alpha = 0.7) {

  # data prepatation
  plot_data <- data
  x_vals <- plot_data[[x_var]]
  y_vals <- plot_data[[y_var]]
  
  # automatic initial value
  if (auto_start) {
    if (is.null(Vmax_start)) Vmax_start <- max(y_vals, na.rm = TRUE) * 1.5
    if (is.null(Km_start)) Km_start <- median(x_vals, na.rm = TRUE)
  }
  
  # model fitting
  model_MM <- try(
    nls(
      as.formula(paste(y_var, "~ (Vmax *", x_var, ") / (Km +", x_var, ")")),
      data = plot_data,
      start = list(Vmax = Vmax_start, Km = Km_start),
      control = nls.control(maxiter = max_iter, warnOnly = TRUE)
    ),
    silent = TRUE
  )
  
  # initialize result object
  result <- list(
    success = FALSE,
    model = NULL,
    parameters = NULL,
    goodness_of_fit = NULL,
    plot = NULL
  )
  
  if (inherits(model_MM, "nls")) {
    result$success <- TRUE
    result$model <- model_MM
    
    # extract parameter
    coef_summary <- summary(model_MM)
    coef_vals <- coef(model_MM)
    coef_se <- coef_summary$coefficients[, "Std. Error"]
    coef_t <- coef_summary$coefficients[, "t value"]
    coef_p <- coef_summary$coefficients[, "Pr(>|t|)"]
    
    Vmax_est <- coef_vals["Vmax"]
    Km_est <- coef_vals["Km"]
    Vmax_se <- coef_se["Vmax"]
    Km_se <- coef_se["Km"]
    Vmax_p <- coef_p["Vmax"]
    Km_p <- coef_p["Km"]
    
    result$parameters <- list(
      Vmax = Vmax_est,
      Km = Km_est
    )
    
    result$parameters_stats <- list(
      Vmax = list(
        estimate = Vmax_est,
        std_error = Vmax_se,
        t_value = coef_t["Vmax"],
        p_value = Vmax_p
      ),
      Km = list(
        estimate = Km_est,
        std_error = Km_se,
        t_value = coef_t["Km"], 
        p_value = Km_p
      )
    )
    
    # calculate model fitting index
    RSS <- sum(residuals(model_MM)^2)
    TSS <- sum((y_vals - mean(y_vals, na.rm = TRUE))^2)
    Pseudo_R2 <- 1 - (RSS / TSS)
    
    # degree of freedom in Residual
    df_residual <- df.residual(model_MM)
    
    result$goodness_of_fit <- list(
      RSS = RSS,
      TSS = TSS,
      Pseudo_R2 = Pseudo_R2,
      AIC = AIC(model_MM),
      BIC = BIC(model_MM),
      df_residual = df_residual
    )
    
    # plot
    if (make_plot) {
      # create exploratory line
      x_range <- seq(min(x_vals), max(x_vals), length.out = 100)
      pred_y <- (Vmax_est * x_range) / (Km_est + x_range)
      pred_df <- data.frame(x = x_range, y = pred_y)
      names(pred_df) <- c(x_var, y_var)
      
      # base plot
      p <- ggplot(plot_data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
        geom_point(size = point_size) +
        geom_line(data = pred_df, color = line_color, linewidth = line_size) +
        geom_hline(yintercept = Vmax_est, linetype = hline_type, color = hline_color) +
        labs(
          x = x_lab,
          y = y_lab,
          title = plot_title
        )
      
      # add label if indicated
      if (!is.null(label_var) && label_var %in% names(plot_data)) {
        if (text_adjust) {
          p <- p + geom_text(aes(label = .data[[label_var]]), 
                             hjust = -0.1, vjust = -0.5, check_overlap = FALSE)
        } else {
          p <- p + geom_text(aes(label = .data[[label_var]]), check_overlap = FALSE)
        }
      }
      
      # apply theme
      p <- p + switch(theme_style,
                      "minimal" = theme_minimal(),
                      "bw" = theme_bw(),
                      "classic" = theme_classic(),
                      theme_minimal()
      )
      
      # basic theme
      p <- p + theme(
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid.major = element_line(color = "lightgray"),
        panel.grid.minor = element_blank()
      )
      
      result$plot <- p
      print(p)
    }
    
    # output
    cat("=== Michaelis-Menten Model Fit Results ===\n")
    cat("Model convergence: SUCCESS\n\n")
    
    cat("Parameter Estimates:\n")
    cat(sprintf("Vmax = %7.3f ± %6.3f (SE)", Vmax_est, Vmax_se))
    cat(sprintf(", t = %6.3f", coef_t["Vmax"]))
    cat(sprintf(", p = %7.4f", Vmax_p))
    cat(ifelse(Vmax_p < 0.05, " *\n", "\n"))
    
    cat(sprintf("Km   = %7.3f ± %6.3f (SE)", Km_est, Km_se))
    cat(sprintf(", t = %6.3f", coef_t["Km"]))
    cat(sprintf(", p = %7.4f", Km_p))
    cat(ifelse(Km_p < 0.05, " *\n", "\n"))
    cat("* p < 0.05\n\n")
    
    cat("Goodness of Fit:\n")
    cat("Pseudo R-squared:", round(Pseudo_R2, 4), "\n")
    cat("AIC:", round(AIC(model_MM), 2), "\n")
    cat("BIC:", round(BIC(model_MM), 2), "\n")
    cat("Residual df:", df_residual, "\n\n")
    
    # add confidence intervals
    conf_int <- try(confint(model_MM), silent = TRUE)
    if (!inherits(conf_int, "try-error")) {
      cat("95% Confidence Intervals:\n")
      cat(sprintf("Vmax: [%6.3f, %6.3f]\n", conf_int["Vmax", 1], conf_int["Vmax", 2]))
      cat(sprintf("Km:   [%6.3f, %6.3f]\n", conf_int["Km", 1], conf_int["Km", 2]))
    }
    
  } else {
    cat("=== Michaelis-Menten Model Fit Results ===\n")
    cat("Model convergence: FAILED\n")
    cat("Message:", as.character(model_MM), "\n")
    cat("Try adjusting starting values or checking data quality.\n")
  }
  
  # return results invisible
  invisible(result)
}

#' Save plot safely（for multiple format）
#' 
#' @export
save_plot_safely <- function(plot_obj, 
                             file_path,
                             width = 10,
                             height = 6,
                             dpi = 300,
                             format = NULL) {  
  
  # automatic format detection or specification
  if (is.null(format)) {
    file_ext <- tolower(tools::file_ext(file_path))
    format <- file_ext
  }
  
  cat("Saving as:", format, "format\n")
  cat("Target path:", file_path, "\n")
  
  # create directory
  dir_path <- dirname(file_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    cat("Created directory:", dir_path, "\n")
  }
  
  result <- list(success = FALSE, path = file_path, format = format)
  
  tryCatch({
    if (format %in% c("pdf", "svg")) {
      # vector format
      ggsave(file_path, plot = plot_obj, width = width, height = height)
      result$success <- TRUE
      cat("SUCCESS: Saved as", format, "\n")
      
    } else if (format == "png") {
      # PNG
      cat("Using cautious PNG saving...\n")
      
      # method 1: low resolution
      tryCatch({
        ggsave(file_path, plot = plot_obj, width = width, height = height, 
               dpi = 150,  # low resoution
               bg = "white")  #  background colour
        result$success <- TRUE
        cat("SUCCESS: PNG saved with lower DPI\n")
      }, error = function(e) {
        cat("Lower DPI failed:", e$message, "\n")
        
        # method 2: other png device
        png(file_path, width = width * 100, height = height * 100, 
            res = 150, type = "cairo")  
        print(plot_obj)
        dev.off()
        result$success <- TRUE
        cat("SUCCESS: PNG saved with Cairo device\n")
      })
      
    } else {
      # other format
      ggsave(file_path, plot = plot_obj, width = width, height = height, dpi = dpi)
      result$success <- TRUE
      cat("SUCCESS: Saved as", format, "\n")
    }
    
  }, error = function(e) {
    cat("Final save attempt failed:", e$message, "\n")
    result$error <- e$message
  })
  
  if (result$success) {
    return(invisible(result))
  } else {
    return(invisible(result))
  }
}

#' Michaelis-Menten kinetics-4
#' @description
#' Model diagnose
#' 
#' @export
diagnose_mm_model <- function(result, plot_data) {
  par(mfrow = c(2, 2))
  
  # 1. visulize model fitting
  plot(plot_data$num_records, plot_data$num_ribotypes,
       xlab = "Number of Records", ylab = "Number of Ribotypes",
       main = "Model Fit", pch = 16)
  x_range <- seq(min(plot_data$num_records), max(plot_data$num_records), length.out = 100)
  Vmax_est <- result$parameters$Vmax
  Km_est <- result$parameters$Km
  lines(x_range, (Vmax_est * x_range) / (Km_est + x_range), col = "red", lwd = 2)
  
  # 2. residual plot
  residuals <- resid(result$model)
  plot(fitted(result$model), residuals,
       xlab = "Fitted Values", ylab = "Residuals",
       main = "Residuals vs Fitted")
  abline(h = 0, col = "red")
  
  # 3. Q-Q plot
  qqnorm(residuals, main = "Q-Q Plot of Residuals")
  qqline(residuals)
  
  # 4. autocorrelation
  acf(residuals, main = "Autocorrelation of Residuals")
  
  par(mfrow = c(1, 1))
}

#' Michaelis-Menten kinetics-5
#' @description
#' Function for interpreting the diagnose (see Michaelis-Menten kinetics-4)
#' 
#' @export
interpret_diagnostic_plots <- function(result, plot_data) {
  cat("=== Diagnose Michaelis-Menten kinetics model fitting ===\n\n")
  
  # basic statistics
  Vmax_est <- result$parameters$Vmax
  Km_est <- result$parameters$Km
  r2 <- result$goodness_of_fit$Pseudo_R2
  
  cat("1. Expected parameter:\n")
  cat("   Vmax =", round(Vmax_est, 2), "(Theoretical maximum)\n")
  cat("   Km =", round(Km_est, 2), "(Half-saturation constant)\n")
  cat("   Pseudo-R² =", round(r2, 4), "\n\n")
  
  # Residual Analysis
  residuals <- resid(result$model)
  fitted_vals <- fitted(result$model)
  
  cat("2. Residual Analysis:\n")
  cat("   Mean residual:", round(mean(residuals), 4), "(better closer to 0 )\n")
  cat("   Standard deviation of residuals:", round(sd(residuals), 4), "\n")
  cat("   Number of standardised residual > |2| :", sum(abs(scale(residuals)) > 2), "\n\n")
  
  # Shapiro-Wilk test (normality)
  if (length(residuals) >= 3 && length(residuals) <= 5000) {
    norm_test <- shapiro.test(residuals)
    cat("3. Normality test for residuals:\n")
    cat("   Shapiro-Wilk test p-value =", round(norm_test$p.value, 4), "\n")
    if (norm_test$p.value > 0.05) {
      cat("   → The residuals follow a normal distribution (good)\n")
    } else {
      cat("   → Residuals deviate from normal distribution (Note)\n")
    }
  }
  
  cat("\n4. Evaluation of Model Fit:\n")
  if (r2 > 0.9) {
    cat("   ✓ Excellent fit (R² > 0.9)\n")
  } else if (r2 > 0.7) {
    cat("   ○ Good fit (R² > 0.7)\n")
  } else if (r2 > 0.5) {
    cat("   △ Acceptable (R² > 0.5)\n")
  } else {
    cat("   × Inadequate compliance (R² < 0.5)\n")
  }
  
  # Example of Calculating Prediction Intervals
  cat("\n5. Example of a prediction:\n")
  example_x <- c(10, 50, 100, 200)
  for (x in example_x) {
    pred_y <- (Vmax_est * x) / (Km_est + x)
    cat("   Number of specimens", x, "→ Expected ribotype count:", round(pred_y, 1), "\n")
  }
}

#' @examples
#' interpret_diagnostic_plots(result, plot_data)


#' Michaelis-Menten kinetics-6
#' 
#' @description
#' Function for alternative model consideration
#' 
#' @export
test_alternative_models <- function(plot_data) {
  models <- list()
  
  # 1. Linear model
  try({
    models$linear <- lm(num_ribotypes ~ num_records, data = plot_data)
  })
  
  # 2. Log-linear model
  try({
    models$log_linear <- lm(log(num_ribotypes + 1) ~ num_records, data = plot_data)
  })
  
  # 3. Power model
  try({
    models$power <- nls(num_ribotypes ~ a * num_records^b,
                        data = plot_data, start = list(a = 1, b = 0.5))
  })
  
  # 4. Logistic model
  try({
    models$logistic <- nls(num_ribotypes ~ Vmax / (1 + exp(-k * (num_records - Km))),
                           data = plot_data, start = list(Vmax = 20, k = 0.1, Km = 50))
  })
  
  # model comparison
  if (length(models) > 0) {
    cat("=== Model Comparison ===\n")
    for (name in names(models)) {
      if (!is.null(models[[name]])) {
        if (inherits(models[[name]], "lm")) {
          r2 <- summary(models[[name]])$r.squared
        } else if (inherits(models[[name]], "nls")) {
          RSS <- sum(residuals(models[[name]])^2)
          TSS <- sum((plot_data$num_ribotypes - mean(plot_data$num_ribotypes))^2)
          r2 <- 1 - (RSS / TSS)
        }
        cat(name, ": R² =", round(r2, 4), "\n")
      }
    }
  }
  
  return(models)
}

#' Snow Cover data-1
#' @description
#' Function for obtain snow cover data in the extracted rater from nc file
#' 
#' @param lon_vals: longitude variables (range) obtained from nc file
#' @param lat_vals: latitude variables (range) obtained from nc file
#' @export
get_indices <- function(lon_vals, lat_vals, easting, northing) {
  lon_index <- which.min(abs(lon_vals - easting))
  lat_index <- which.min(abs(lat_vals - northing))
  return(c(lon_index, lat_index))
}

#' Snow Cover data-2
#' @description
#' Function for obtain snow cover data with coordinate-based data
#' 
#' @export
calculate_consecutive_no_snow_periods <- function(mask_data) {
  num_lon <- dim(mask_data)[1]
  num_lat <- dim(mask_data)[2]
  num_time <- dim(mask_data)[3]
  
  no_snow_periods <- matrix(NA, nrow = num_lon, ncol = num_lat)
  
  for (i in 1:num_lon) {
    for (j in 1:num_lat) {
      # Retrieve snow-free conditions at each location
      mask <- mask_data[i, j, ]
      
      # Calculate the average of the duration of snow-free conditions for each coordinate point
      rle_vals <- rle(mask)
      no_snow_periods[i, j] <- max(rle_vals$lengths[rle_vals$values == TRUE], na.rm = TRUE)
    }
  }
  
  return(no_snow_periods)
}
