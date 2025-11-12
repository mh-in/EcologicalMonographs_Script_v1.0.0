# computational environment

sink("targeted_computational_environment.txt")

cat("=== Enhanced Targeted Computational Environment Report ===\n\n")
cat("Report generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Working directory:", getwd(), "\n\n")

cat("This report includes packages used in specified analysis scripts.\n\n")

target_scripts <- c(
  "01_data_cleaning_processing.R",
  "02_basic_statistics.r",
  "03_geographic_differentiation.r", 
  "04_beta_diversity_analysis.r",
  "05_ecological_analysis.r",
  "utils_ecology_myxo.R"
)

cat("TARGETED SCRIPTS:\n")
for(script in target_scripts) {
  exists <- file.exists(script)
  status <- ifelse(exists, "[FOUND]", "[NOT FOUND]")
  cat("-", script, status, "\n")
}
cat("\n")

cat("1. R VERSION INFORMATION:\n")
cat("Version:", R.version$version.string, "\n")
cat("Platform:", R.version$platform, "\n\n")

cat("2. PACKAGES FROM TARGETED SCRIPTS:\n\n")

all_target_packages <- character(0)

for(script in target_scripts) {
  if(file.exists(script)) {
    cat("Analyzing:", script, "\n")
    
    # A: detect with renv
    script_deps <- tryCatch({
      renv::dependencies(path = script, progress = FALSE)
    }, error = function(e) NULL)
    
    if(!is.null(script_deps) && nrow(script_deps) > 0) {
      packages_auto <- unique(script_deps$Package)
      all_target_packages <- c(all_target_packages, packages_auto)
      cat("  Auto-detected:", length(packages_auto), "packages\n")
    }
    
    # B: create package list manually
    content <- readLines(script, warn = FALSE)
    
    # detect list "required_packages"
    required_lines <- grep("required_packages.*c\\([^)]+\\)", content, value = TRUE)
    if(length(required_lines) > 0) {
      # estract package
      pkg_string <- gsub('.*required_packages.*c\\(([^)]+)\\).*', '\\1', required_lines[1])
      pkg_list <- unlist(strsplit(gsub('["\']', '', pkg_string), ","))
      pkg_list <- trimws(pkg_list)
      
      cat("  Manual list found:", length(pkg_list), "packages\n")
      all_target_packages <- c(all_target_packages, pkg_list)
    }
    
    # C: detect the use of the colon operator
    colon_matches <- regmatches(content, gregexpr("[a-zA-Z0-9.]+::[a-zA-Z0-9._]+", content))
    colon_packages <- unique(sapply(unlist(colon_matches), function(x) strsplit(x, "::")[[1]][1]))
    if(length(colon_packages) > 0) {
      cat("  Colon operator:", length(colon_packages), "packages\n")
      all_target_packages <- c(all_target_packages, colon_packages)
    }
    
    cat("\n")
  }
}

# remove duplicates
if(length(all_target_packages) > 0) {
  unique_packages <- unique(all_target_packages)
  cat("SUMMARY:\n")
  cat("Total unique packages found:", length(unique_packages), "\n\n")
  
  # obtain version info
  package_info <- data.frame(
    Package = character(0),
    Version = character(0),
    stringsAsFactors = FALSE
  )
  
  for(pkg in sort(unique_packages)) {
    version <- tryCatch({
      as.character(packageVersion(pkg))
    }, error = function(e) "Not installed")
    
    package_info <- rbind(package_info, data.frame(Package = pkg, Version = version))
  }
  
  print(package_info, row.names = FALSE)
  
  # verification of particularly important packages
  cat("\n4. KEY ECOLOGICAL PACKAGES VERIFICATION:\n")
  key_packages <- c("vegan", "iNEXT", "proxy", "phyloseq", "ape", "geosphere", "ggplot2", "patchwork")
  missing_key <- setdiff(key_packages, unique_packages)
  if(length(missing_key) > 0) {
    cat("WARNING: The following key packages were not detected:\n")
    cat(paste("-", missing_key, collapse = "\n"), "\n")
    cat("They are used in the analysis but may not be properly detected.\n")
  } else {
    cat("All key ecological packages were successfully detected.\n")
  }
  
} else {
  cat("No packages detected in any of the targeted scripts.\n")
}

cat("\n5. RENV STATUS:\n")
renv_summary <- capture.output(renv::status())
cat(paste(renv_summary, collapse = "\n"))

sink()

cat("Enhanced targeted computational environment report saved.\n")

