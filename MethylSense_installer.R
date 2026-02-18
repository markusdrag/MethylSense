#!/usr/bin/env Rscript

# ================================================================================
# MethylSense Complete Package Installer
# ================================================================================
# Version: 5.6.6
# Date: 2026-01-29
# Description: Installs all required R packages for the MethylSense workflow
#              (load_data, main analysis, and predict scripts)
#
# Usage: Rscript MethylSense_installer_v1.R
#
# Estimated Runtime: 15-45 minutes (depending on system and existing packages)
# ================================================================================

cat("\n")
cat("================================================================================\n")
cat("  METHYLSENSE COMPLETE PACKAGE INSTALLER v1.0.0\n")
cat("================================================================================\n")
cat("\n")
cat("This script will install all required R packages for:\n")
cat("  • MethylSense_load_data.R       - Data preprocessing\n")
cat("  • MethylSense_analysis.R        - DMR detection & ML training\n")
cat("  • MethylSense_predict.R         - Prediction on new samples\n")
cat("\n")
cat("Estimated installation time: 15-45 minutes\n")
cat("================================================================================\n\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))
options(timeout = 600) # Increase timeout for large packages

# ================================================================================
# PACKAGE LISTS
# ================================================================================

# Core packages (used across all scripts)
core_packages <- c(
  "optparse", # Command-line argument parsing
  "data.table", # Fast data manipulation
  "dplyr", # Data manipulation
  "qs" # Fast serialization for R objects
)

# Data loading packages (MethylSense_load_data.R)
data_packages <- c(
  "readxl", # Read Excel files (sample metadata)
  "parallel" # Parallel processing
)

# Machine learning packages (MethylSense_analysis.R + predict)
ml_packages <- c(
  "caret", # ML framework
  "randomForest", # Random Forest
  "ranger", # Fast Random Forest implementation
  "e1071", # SVM
  "glmnet", # Elastic Net (glmnet)
  "xgboost", # Gradient Boosting
  "naivebayes", # Naive Bayes
  "nnet", # Neural Networks
  "MASS", # LDA
  "class", # kNN
  "klaR", # PAM (Nearest Shrunken Centroids)
  "pls" # Partial Least Squares
)

# Model evaluation packages
eval_packages <- c(
  "pROC", # ROC curves and AUC
  "PRROC", # Precision-Recall curves
  "MLmetrics" # Additional ML metrics
)

# Visualization packages
viz_packages <- c(
  "ggplot2", # Publication-quality plots
  "pheatmap", # Heatmaps
  "RColorBrewer", # Color palettes
  "viridis", # Colorblind-friendly palettes
  "corrplot", # Correlation plots
  "factoextra", # PCA visualization
  "ggridges", # Ridge plots
  "gridExtra", # Arrange multiple plots (for predict script)
  "ggrepel", # Text labels for plots (reviewer script)
  "patchwork" # Combine multiple ggplots (CpG barplots)
)

# Tree/ensemble packages
tree_packages <- c(
  "tree", # Decision trees
  "rpart", # Recursive partitioning
  "rpart.plot" # Plot rpart trees
)

# Parallel processing packages
parallel_packages <- c(
  "parallel", # Base parallel
  "doParallel", # Parallel backend for caret
  "foreach" # Parallel foreach loops
)

# Additional utility packages
utility_packages <- c(
  "reshape2", # Data reshaping
  "cluster", # Clustering algorithms
  "jsonlite", # JSON I/O for nested CV hyperparameters
  "igraph", # Network analysis (reviewer script)
  "umap", # UMAP dimensionality reduction (reviewer script)
  "dendextend", # Dendrogram manipulation (reviewer script)
  "pagedown" # PDF generation for clinical reports (reviewer script)
)

# Statistical analysis packages (for reviewer script)
statistical_packages <- c(
  "lme4", # Linear mixed effects models (time series)
  "nlme", # Non-linear mixed effects (alternative)
  "MuMIn", # Multi-model inference (R² calculations)
  "ez", # ANOVA analysis
  "fossil" # Rand index for clustering metrics
)

# Combine all CRAN packages
all_cran_packages <- unique(c(
  core_packages,
  data_packages,
  ml_packages,
  eval_packages,
  viz_packages,
  tree_packages,
  parallel_packages,
  utility_packages,
  statistical_packages
))

# Bioconductor packages (methylation analysis)
bioc_packages <- c(
  "methylKit", # Methylation data analysis (CORE)
  "GenomicRanges", # Genomic coordinates manipulation
  "genomation", # Genomic region annotation
  "regioneR", # Region enrichment analysis
  "IRanges" # Integer ranges operations
)

# Optional Bioconductor packages (for reviewer script advanced features)
# These are optional and installed only if user wants advanced analyses
bioc_optional_packages <- c(
  "biomaRt", # Gene annotation (optional)
  "clusterProfiler", # Functional enrichment (optional)
  "WGCNA" # Co-expression network analysis (optional)
)

# ================================================================================
# INSTALLATION FUNCTIONS
# ================================================================================

# Function to check and install CRAN packages
install_cran_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("[INSTALL] %-20s ... ", pkg))
    tryCatch(
      {
        suppressMessages(
          install.packages(pkg, dependencies = TRUE, quiet = TRUE)
        )
        cat("✓ SUCCESS\n")
        return(TRUE)
      },
      error = function(e) {
        cat("✗ FAILED\n")
        cat(sprintf("          Error: %s\n", e$message))
        return(FALSE)
      }
    )
  } else {
    cat(sprintf("[SKIP]    %-20s ... already installed\n", pkg))
    return(TRUE)
  }
}

# Function to check and install Bioconductor packages
install_bioc_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("[INSTALL] %-20s ... ", pkg))
    tryCatch(
      {
        suppressMessages(
          BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
        )
        cat("✓ SUCCESS\n")
        return(TRUE)
      },
      error = function(e) {
        cat("✗ FAILED\n")
        cat(sprintf("          Error: %s\n", e$message))
        return(FALSE)
      }
    )
  } else {
    cat(sprintf("[SKIP]    %-20s ... already installed\n", pkg))
    return(TRUE)
  }
}

# ================================================================================
# INSTALL CRAN PACKAGES
# ================================================================================

cat("\n=== STEP 1/3: Installing CRAN packages ===\n\n")
cat(sprintf("Total packages to check: %d\n\n", length(all_cran_packages)))

cran_success <- sapply(all_cran_packages, install_cran_package)

cat(sprintf(
  "\nCRAN packages installed: %d/%d\n",
  sum(cran_success), length(all_cran_packages)
))

# ================================================================================
# INSTALL BIOCONDUCTOR PACKAGES
# ================================================================================

cat("\n=== STEP 2/3: Installing Bioconductor packages ===\n\n")

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("[INSTALL] BiocManager          ... ")
  tryCatch(
    {
      install.packages("BiocManager", quiet = TRUE)
      cat("✓ SUCCESS\n\n")
    },
    error = function(e) {
      cat("✗ FAILED\n")
      cat("ERROR: Cannot install BiocManager. Please install manually.\n")
      quit(status = 1)
    }
  )
} else {
  cat("[SKIP]    BiocManager          ... already installed\n\n")
}

cat(sprintf("Total Bioconductor packages to check: %d\n\n", length(bioc_packages)))

bioc_success <- sapply(bioc_packages, install_bioc_package)

cat(sprintf(
  "\nBioconductor packages installed: %d/%d\n",
  sum(bioc_success), length(bioc_packages)
))

# ================================================================================
# INSTALL OPTIONAL BIOCONDUCTOR PACKAGES
# ================================================================================

cat("\n=== STEP 2.5/3: Installing Optional Bioconductor packages ===\n\n")
cat("These packages enable advanced features in MethylSense_reviewer.R:\n")
cat("  - biomaRt: Gene annotation\n")
cat("  - clusterProfiler: Functional enrichment analysis\n")
cat("  - WGCNA: Co-expression network analysis\n\n")

cat(sprintf("Total optional packages to check: %d\n\n", length(bioc_optional_packages)))

bioc_optional_success <- sapply(bioc_optional_packages, install_bioc_package)

cat(sprintf(
  "\nOptional Bioconductor packages installed: %d/%d\n",
  sum(bioc_optional_success), length(bioc_optional_packages)
))

if (sum(bioc_optional_success) < length(bioc_optional_packages)) {
  cat("\nNOTE: Some optional packages failed to install.\n")
  cat("      MethylSense will work without them, but advanced reviewer features\n")
  cat("      (gene annotation, enrichment, WGCNA) will not be available.\n")
}

# ================================================================================
# VERIFICATION
# ================================================================================

cat("\n=== STEP 3/3: Verification ===\n\n")

# Check required packages (must have)
all_required_packages <- c(all_cran_packages, bioc_packages)
missing_required <- c()
installed_required <- c()

for (pkg in all_required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_required <- c(missing_required, pkg)
  } else {
    installed_required <- c(installed_required, pkg)
  }
}

# Check optional packages (nice to have)
missing_optional <- c()
installed_optional <- c()

for (pkg in bioc_optional_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_optional <- c(missing_optional, pkg)
  } else {
    installed_optional <- c(installed_optional, pkg)
  }
}

cat(sprintf("REQUIRED packages checked: %d\n", length(all_required_packages)))
cat(sprintf("  Successfully installed: %d\n", length(installed_required)))
cat(sprintf("  Failed installations:   %d\n", length(missing_required)))

cat(sprintf("\nOPTIONAL packages checked: %d\n", length(bioc_optional_packages)))
cat(sprintf("  Successfully installed: %d\n", length(installed_optional)))
cat(sprintf("  Failed installations:   %d\n", length(missing_optional)))

# Use required packages for final success/failure determination
all_packages <- all_required_packages
missing_packages <- missing_required
installed_packages <- installed_required

if (length(missing_packages) == 0) {
  cat("\n")
  cat("================================================================================\n")
  cat("  ✓ ALL PACKAGES INSTALLED SUCCESSFULLY!\n")
  cat("================================================================================\n\n")
  cat("You can now run the MethylSense workflow:\n\n")
  cat("  1. Load data:\n")
  cat("     Rscript MethylSense_load_data.R --bed_dir <dir> --output_file <file>.qs\n\n")
  cat("  2. Run analysis:\n")
  cat("     Rscript MethylSense_analysis.R --qs_file <file>.qs --window_file <bed>\n\n")
  cat("  3. Evaluate model (optional):\n")
  cat("     Rscript MethylSense_reviewer.R --model_dir <dir> --qs_file <file>.qs\n\n")
  cat("  4. Predict new samples:\n")
  cat("     Rscript MethylSense_predict.R --model_dir <dir> --new_data <file>.qs\n\n")
} else {
  cat("\n")
  cat("================================================================================\n")
  cat("  ✗ INSTALLATION INCOMPLETE\n")
  cat("================================================================================\n\n")
  cat(sprintf("The following %d package(s) failed to install:\n\n", length(missing_packages)))
  for (pkg in missing_packages) {
    cat(sprintf("  • %s\n", pkg))
  }
  cat("\n")
  cat("Please try installing them manually:\n\n")
  cat("For CRAN packages:\n")
  cat(sprintf("  install.packages(c('%s'))\n\n", paste(missing_packages, collapse = "', '")))
  cat("For Bioconductor packages:\n")
  cat(sprintf("  BiocManager::install(c('%s'))\n\n", paste(missing_packages, collapse = "', '")))
}

# ================================================================================
# SAVE SESSION INFO
# ================================================================================

cat("\n=== Saving package versions ===\n\n")
session_file <- "MethylSense_package_versions.txt"

tryCatch(
  {
    sink(session_file)
    cat("MethylSense Package Installation Summary\n")
    cat(paste0("Date: ", Sys.time(), "\n"))
    cat(paste0("R Version: ", R.version.string, "\n"))
    cat(paste0("Platform: ", R.version$platform, "\n\n"))
    cat(rep("=", 80), "\n\n", sep = "")
    print(sessionInfo())
    sink()
    cat(sprintf("✓ Package versions saved to: %s\n", session_file))
  },
  error = function(e) {
    sink() # Ensure sink is closed even on error
    cat("✗ Failed to save package versions\n")
  }
)

# ================================================================================
# FINAL STATUS
# ================================================================================

cat("\n")
cat("================================================================================\n")
if (length(missing_packages) == 0) {
  cat("  Installation Status: ✓ SUCCESS\n")
  quit(status = 0)
} else {
  cat("  Installation Status: ✗ INCOMPLETE\n")
  quit(status = 1)
}
cat("================================================================================\n\n")
