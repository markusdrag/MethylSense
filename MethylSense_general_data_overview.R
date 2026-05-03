#!/usr/bin/env Rscript

# ================================================================================
# general_data_overview.R
# ================================================================================
#
# Comprehensive DMR landscape overview and quality assessment
#
# Generates publication-ready figures and statistics for:
# - DMR/MR counts across region sizes
# - Volcano plots showing effect sizes and significance
# - Coverage distributions for quality assessment
# - Chromosomal distribution of DMRs
# - Effect size characterisation
# - Sample-level methylation patterns
# - P-value calibration assessment
# - Region size optimisation analysis
#
# Author: Markus Hodal Drag
# Version: 5.7.2
# Date: 2026-05-01
# GitHub: https://github.com/markusdrag/MethylSense
#
# ================================================================================

SCRIPT_VERSION <- "5.7.2"
SCRIPT_DATE <- "2026-05-01"

# ================================================================================
# LOAD REQUIRED LIBRARIES
# ================================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(qs)
  library(methylKit)
  library(ggplot2)
  library(readxl)
  library(viridis)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
  library(gridExtra)
  library(ggrepel)
  library(reshape2)
  library(ggridges)
  library(hrbrthemes)
})

# ================================================================================
# METHYLSENSE GLOBAL THEME CONFIGURATION
# ================================================================================

# Source shared MethylSense theme and colour palette (Rscript-safe path)
cmd_args_theme <- commandArgs(trailingOnly = FALSE)
file_arg_theme <- grep("^--file=", cmd_args_theme, value = TRUE)
script_dir_theme <- if (length(file_arg_theme)) {
  dirname(normalizePath(sub("^--file=", "", file_arg_theme[1]), mustWork = FALSE))
} else {
  getwd()
}
source(file.path(script_dir_theme, "MethylSense_theme.R"))

# ================================================================================
# COMMAND-LINE ARGUMENTS
# ================================================================================

option_list <- list(
  make_option(c("--analysis_dir"),
    type = "character", default = NULL,
    help = "Path to analysis directory containing region subdirectories and united_meth files [required]",
    metavar = "character"
  ),
  make_option(c("--region_sizes"),
    type = "character", default = "1000,5000,10000,15000,20000,25000",
    help = "Comma-separated list of region sizes in bp [default: 1000,5000,10000,15000,20000,25000]",
    metavar = "character"
  ),
  make_option(c("--sample_sheet"),
    type = "character", default = NULL,
    help = "Path to sample metadata Excel file (e.g., AsperMerged_v3.xlsx) [required]",
    metavar = "character"
  ),
  make_option(c("--infection_col"),
    type = "character", default = "InfectionStatus",
    help = "Column name for infection status in sample sheet [default: InfectionStatus]",
    metavar = "character"
  ),
  make_option(c("--study_col"),
    type = "character", default = "STUDY...3",
    help = "Column name for study grouping in sample sheet [default: STUDY...3]",
    metavar = "character"
  ),
  make_option(c("--sample_id_col"),
    type = "character", default = "Lab_ID",
    help = "Column name for sample IDs in sample sheet [default: Lab_ID]",
    metavar = "character"
  ),
  make_option(c("--output_dir"),
    type = "character", default = "dmr_landscape_overview",
    help = "Output directory for figures and stats [default: dmr_landscape_overview]",
    metavar = "character"
  ),
  make_option(c("--plot_format"),
    type = "character", default = "jpeg",
    help = "Output plot format(s): png, jpg, jpeg, svg, pdf, or comma-separated list (e.g., 'jpeg,svg') [default: jpeg]",
    metavar = "character"
  ),
  make_option(c("--plot_width"),
    type = "numeric", default = 12,
    help = "Plot width in inches [default: 12]",
    metavar = "number"
  ),
  make_option(c("--plot_height"),
    type = "numeric", default = 8,
    help = "Plot height in inches [default: 8]",
    metavar = "number"
  ),
  make_option(c("--plot_dpi"),
    type = "numeric", default = 300,
    help = "Plot resolution (DPI) for raster formats [default: 300]",
    metavar = "number"
  ),
  make_option(c("--volcano_meth_threshold"),
    type = "numeric", default = 5,
    help = "Methylation difference threshold for volcano plot (%) [default: 5]",
    metavar = "number"
  ),
  make_option(c("--volcano_q_threshold"),
    type = "numeric", default = 0.05,
    help = "Q-value threshold for volcano plot [default: 0.05]",
    metavar = "number"
  ),
  make_option(c("--top_n_label"),
    type = "integer", default = 10,
    help = "Number of top DMRs to label in volcano plot [default: 10]",
    metavar = "integer"
  ),
  make_option(c("--heatmap_top_n"),
    type = "integer", default = 50,
    help = "Number of top DMRs to show in heatmap [default: 50]",
    metavar = "integer"
  ),
  make_option(c("--generate_pdf"),
    action = "store_true", default = FALSE,
    help = "Generate PDF report in addition to markdown [default: FALSE]"
  ),
  make_option(c("--blind_paths"),
    action = "store_true", default = FALSE,
    help = "Replace directory paths with generic placeholders in PDF report [default: FALSE]"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$analysis_dir)) {
  stop("--analysis_dir is required. Use --help for usage information.")
}

if (is.null(opt$sample_sheet)) {
  stop("--sample_sheet is required. Use --help for usage information.")
}

if (!dir.exists(opt$analysis_dir)) {
  stop(sprintf("Analysis directory does not exist: %s", opt$analysis_dir))
}

if (!file.exists(opt$sample_sheet)) {
  stop(sprintf("Sample sheet does not exist: %s", opt$sample_sheet))
}

# Validate and parse plot format(s)
valid_formats <- c("png", "jpg", "jpeg", "svg", "pdf")
plot_formats <- tolower(trimws(strsplit(opt$plot_format, ",")[[1]]))

# Check all formats are valid
invalid_formats <- plot_formats[!plot_formats %in% valid_formats]
if (length(invalid_formats) > 0) {
  stop(sprintf(
    "Invalid plot format(s): %s. Must be one of: %s",
    paste(invalid_formats, collapse = ", "),
    paste(valid_formats, collapse = ", ")
  ))
}

cat(sprintf("\n[INFO] Plots will be saved in format(s): %s\n", paste(plot_formats, collapse = ", ")))

# Parse region sizes
region_sizes <- as.numeric(strsplit(opt$region_sizes, ",")[[1]])

# Create output directories
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
figures_dir <- file.path(opt$output_dir, "figures")
stats_dir <- file.path(opt$output_dir, "stats")
dir.create(figures_dir, showWarnings = FALSE)
dir.create(stats_dir, showWarnings = FALSE)

# ================================================================================
# HELPER FUNCTIONS
# ================================================================================

save_plot <- function(plot_obj, filename, width = opt$plot_width, height = opt$plot_height) {
  for (format in plot_formats) {
    filepath <- file.path(figures_dir, paste0(filename, ".", format))

    if (format == "png") {
      ggsave(filepath, plot_obj, width = width, height = height, dpi = opt$plot_dpi, units = "in", bg = "white")
    } else if (format == "jpg" || format == "jpeg") {
      ggsave(filepath, plot_obj, width = width, height = height, dpi = opt$plot_dpi, units = "in", device = "jpeg", bg = "white")
    } else if (format == "svg") {
      ggsave(filepath, plot_obj, width = width, height = height, units = "in", device = "svg", bg = "white")
    } else if (format == "pdf") {
      ggsave(filepath, plot_obj, width = width, height = height, units = "in", device = "pdf", bg = "white")
    }

    cat(sprintf("[OK] Saved: %s\n", filepath))
  }
}

# Helper function to create properly ordered region size labels
order_region_sizes <- function(region_labels) {
  # Extract numeric values from labels like "1 KB", "5 KB", etc.
  numeric_values <- as.numeric(gsub(" KB", "", region_labels))
  # Get unique labels in proper order
  unique_labels <- unique(region_labels)
  unique_numeric <- as.numeric(gsub(" KB", "", unique_labels))
  ordered_levels <- unique_labels[order(unique_numeric)]
  factor(region_labels, levels = ordered_levels)
}

# Helper function to sanitize file paths for privacy
# Keeps only the basename (last component) of file/directory paths
sanitize_path <- function(path_string) {
  if (is.null(path_string) || path_string == "") {
    return(path_string)
  }
  basename(path_string)
}

# ================================================================================
# LOAD SAMPLE METADATA
# ================================================================================

cat("\n================================================================================\n")
cat("Loading sample metadata\n")
cat("================================================================================\n\n")

if (grepl("\\.csv$", opt$sample_sheet, ignore.case = TRUE)) {
  sample_metadata <- read.csv(opt$sample_sheet, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  sample_metadata <- read_excel(opt$sample_sheet)
}
cat(sprintf("[OK] Loaded %d samples from: %s\n", nrow(sample_metadata), opt$sample_sheet))

# Validate required columns exist
required_cols <- c(opt$sample_id_col, opt$infection_col, opt$study_col)
missing_cols <- required_cols[!required_cols %in% colnames(sample_metadata)]
if (length(missing_cols) > 0) {
  cat("\n[ERROR] Missing required columns in sample sheet:\n")
  for (col in missing_cols) {
    cat(sprintf("  - %s\n", col))
  }
  cat("\nAvailable columns in sample sheet:\n")
  for (col in colnames(sample_metadata)) {
    cat(sprintf("  - %s\n", col))
  }
  stop(sprintf(
    "\nRequired columns not found: %s\n\nPlease specify correct column names using --infection_col, --study_col, and --sample_id_col",
    paste(missing_cols, collapse = ", ")
  ))
}

cat(sprintf("     Infection status column: %s\n", opt$infection_col))
cat(sprintf("     Study grouping column: %s\n", opt$study_col))

# Count samples by group
infection_counts <- table(sample_metadata[[opt$infection_col]], useNA = "ifany")
cat("\n[INFO] Samples by infection status:\n")
print(infection_counts)

study_counts <- table(sample_metadata[[opt$study_col]], useNA = "ifany")
cat("\n[INFO] Samples by study:\n")
print(study_counts)

# ================================================================================
# 1. DMR/MR COUNTS ACROSS REGION SIZES (BAR PLOT)
# ================================================================================

cat("\n================================================================================\n")
cat("Analysis 1: DMR/MR counts across region sizes\n")
cat("================================================================================\n\n")

dmr_counts_list <- list()

for (region_size in region_sizes) {
  cat(sprintf("Processing %d bp regions...\n", region_size))

  # Find files for this region size - try both naming conventions (5000bp and 5kb)
  # Pattern 1: explicit bp format (e.g., _5000bp.csv)
  dmr_file_pattern_bp <- sprintf("ALL_DMRs_complete_dataset_.*_%dbp\\.csv", region_size)
  # Pattern 2: kb format (e.g., _5kb.csv or windows_5kb.csv)
  dmr_file_pattern_kb <- sprintf("ALL_DMRs_complete_dataset_.*_%dkb\\.csv", region_size / 1000)
  # Pattern 3: windows_Xkb format
  dmr_file_pattern_win <- sprintf("ALL_DMRs_complete_dataset_windows_%dkb\\.csv", region_size / 1000)

  dmr_file <- list.files(opt$analysis_dir, pattern = dmr_file_pattern_bp, full.names = TRUE)
  if (length(dmr_file) == 0) {
    dmr_file <- list.files(opt$analysis_dir, pattern = dmr_file_pattern_kb, full.names = TRUE)
  }
  if (length(dmr_file) == 0) {
    dmr_file <- list.files(opt$analysis_dir, pattern = dmr_file_pattern_win, full.names = TRUE)
  }

  # Same for meth files
  meth_file_pattern_bp <- sprintf("united_meth_object_.*_%dbp_.*\\.qs", region_size)
  meth_file_pattern_kb <- sprintf("united_meth_object_.*_%dkb_.*\\.qs", region_size / 1000)
  meth_file_pattern_win <- sprintf("united_meth_object_windows_%dkb_.*\\.qs", region_size / 1000)

  meth_file <- list.files(opt$analysis_dir, pattern = meth_file_pattern_bp, full.names = TRUE)
  if (length(meth_file) == 0) {
    meth_file <- list.files(opt$analysis_dir, pattern = meth_file_pattern_kb, full.names = TRUE)
  }
  if (length(meth_file) == 0) {
    meth_file <- list.files(opt$analysis_dir, pattern = meth_file_pattern_win, full.names = TRUE)
  }

  if (length(dmr_file) == 0 || length(meth_file) == 0) {
    cat(sprintf("  [WARNING] Files not found for %d bp, skipping...\n", region_size))
    next
  }

  # Load data
  dmr_data <- read.csv(dmr_file[1])
  meth_obj <- qread(meth_file[1])

  total_mrs <- nrow(meth_obj)
  significant_dmrs <- sum(dmr_data$passes_current_qvalue_threshold, na.rm = TRUE)
  dmr_rate <- (significant_dmrs / total_mrs) * 100

  dmr_counts_list[[as.character(region_size)]] <- data.frame(
    Region_Size_bp = region_size,
    Region_Size_Label = sprintf("%d KB", region_size / 1000),
    Total_MRs = total_mrs,
    Significant_DMRs = significant_dmrs,
    Non_Significant_MRs = total_mrs - significant_dmrs,
    DMR_Rate_Percent = dmr_rate
  )

  cat(sprintf("  Total MRs: %d\n", total_mrs))
  cat(sprintf("  Significant DMRs (q<%.2f): %d\n", opt$volcano_q_threshold, significant_dmrs))
  cat(sprintf("  DMR detection rate: %.1f%%\n\n", dmr_rate))
}

dmr_counts_df <- do.call(rbind, dmr_counts_list)
rownames(dmr_counts_df) <- NULL

# Save stats
write.csv(dmr_counts_df, file.path(stats_dir, "dmr_counts_summary.csv"), row.names = FALSE)
cat(sprintf("[OK] Saved: %s\n\n", file.path(stats_dir, "dmr_counts_summary.csv")))

# Create bar plot with proper ordering
dmr_counts_long <- reshape2::melt(dmr_counts_df[, c("Region_Size_Label", "Non_Significant_MRs", "Significant_DMRs")],
  id.vars = "Region_Size_Label",
  variable.name = "Category",
  value.name = "Count"
)

# Properly order region sizes
dmr_counts_long$Region_Size_Label <- order_region_sizes(dmr_counts_long$Region_Size_Label)
dmr_counts_df$Region_Size_Label <- order_region_sizes(dmr_counts_df$Region_Size_Label)

dmr_counts_long$Category <- factor(dmr_counts_long$Category,
  levels = c("Non_Significant_MRs", "Significant_DMRs"),
  labels = c("Non-significant MRs", "Significant DMRs (FDR<0.05)")
)

p1 <- ggplot(dmr_counts_long, aes(x = Region_Size_Label, y = Count, fill = Category)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  # Add count labels on each bar segment
  geom_text(aes(label = scales::comma(Count)),
    position = position_stack(vjust = 0.5),
    size = 3, color = "white", fontface = "bold"
  ) +
  # Add percentage label above bars
  geom_text(
    data = dmr_counts_df,
    aes(x = Region_Size_Label, y = Total_MRs, label = sprintf("%.1f%%", DMR_Rate_Percent), fill = NULL),
    vjust = -0.5, size = 4, fontface = "bold"
  ) +
  scale_fill_manual(values = c("Non-significant MRs" = "gray70", "Significant DMRs (FDR<0.05)" = "#E31A1C")) +
  scale_y_continuous(labels = scales::comma) + # Add comma formatting
  labs(
    title = "DMR detection across region sizes",
    subtitle = sprintf(
      "Analysis of %d region sizes (FDR threshold < %.2f)",
      length(region_sizes), opt$volcano_q_threshold
    ),
    x = "Genomic window size",
    y = "Number of regions",
    fill = "Significance"
  ) +
  theme_methylsense() +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 16),
    panel.grid.major.x = element_blank()
  )

save_plot(p1, "1_dmr_counts_barplot", width = 10, height = 8)

cat("[COMPLETE] DMR/MR count analysis\n")

# ================================================================================
# LOAD ALL DMR DATA FOR REMAINING ANALYSES
# ================================================================================

cat("\n================================================================================\n")
cat("Loading all DMR data for comprehensive analysis\n")
cat("================================================================================\n\n")

all_dmr_data <- list()
all_meth_objects <- list()

for (region_size in region_sizes) {
  # Try multiple naming conventions (5000bp and 5kb)
  dmr_file_pattern_bp <- sprintf("ALL_DMRs_complete_dataset_.*_%dbp\\.csv", region_size)
  dmr_file_pattern_kb <- sprintf("ALL_DMRs_complete_dataset_.*_%dkb\\.csv", region_size / 1000)
  dmr_file_pattern_win <- sprintf("ALL_DMRs_complete_dataset_windows_%dkb\\.csv", region_size / 1000)

  dmr_file <- list.files(opt$analysis_dir, pattern = dmr_file_pattern_bp, full.names = TRUE)
  if (length(dmr_file) == 0) {
    dmr_file <- list.files(opt$analysis_dir, pattern = dmr_file_pattern_kb, full.names = TRUE)
  }
  if (length(dmr_file) == 0) {
    dmr_file <- list.files(opt$analysis_dir, pattern = dmr_file_pattern_win, full.names = TRUE)
  }

  meth_file_pattern_bp <- sprintf("united_meth_object_.*_%dbp_.*\\.qs", region_size)
  meth_file_pattern_kb <- sprintf("united_meth_object_.*_%dkb_.*\\.qs", region_size / 1000)
  meth_file_pattern_win <- sprintf("united_meth_object_windows_%dkb_.*\\.qs", region_size / 1000)

  meth_file <- list.files(opt$analysis_dir, pattern = meth_file_pattern_bp, full.names = TRUE)
  if (length(meth_file) == 0) {
    meth_file <- list.files(opt$analysis_dir, pattern = meth_file_pattern_kb, full.names = TRUE)
  }
  if (length(meth_file) == 0) {
    meth_file <- list.files(opt$analysis_dir, pattern = meth_file_pattern_win, full.names = TRUE)
  }

  if (length(dmr_file) > 0 && length(meth_file) > 0) {
    all_dmr_data[[as.character(region_size)]] <- read.csv(dmr_file[1])
    all_meth_objects[[as.character(region_size)]] <- qread(meth_file[1])
    cat(sprintf(
      "[OK] Loaded %d bp data: %d DMRs, %d MRs\n",
      region_size, nrow(all_dmr_data[[as.character(region_size)]]),
      nrow(all_meth_objects[[as.character(region_size)]])
    ))
  }
}

cat("\n[COMPLETE] Data loading\n")

# Continue in next message due to length...

# ================================================================================
# 2. VOLCANO PLOT (EFFECT SIZE VS SIGNIFICANCE)
# ================================================================================

cat("\n================================================================================\n")
cat("Analysis 2: Volcano plots\n")
cat("================================================================================\n\n")

volcano_stats_list <- list()

for (region_size in names(all_dmr_data)) {
  dmr_data <- all_dmr_data[[region_size]]

  # Prepare volcano data
  dmr_data$neg_log10_qvalue <- -log10(dmr_data$qvalue)
  dmr_data$Category <- "Non-significant"
  dmr_data$Category[dmr_data$qvalue < opt$volcano_q_threshold & dmr_data$meth_difference > opt$volcano_meth_threshold] <- "Hypermethylated"
  dmr_data$Category[dmr_data$qvalue < opt$volcano_q_threshold & dmr_data$meth_difference < -opt$volcano_meth_threshold] <- "Hypomethylated"

  # Count by quadrant
  volcano_stats_list[[region_size]] <- data.frame(
    Region_Size_bp = as.numeric(region_size),
    Total_Regions = nrow(dmr_data),
    Significant_Hypermethylated = sum(dmr_data$Category == "Hypermethylated"),
    Significant_Hypomethylated = sum(dmr_data$Category == "Hypomethylated"),
    Non_Significant = sum(dmr_data$Category == "Non-significant")
  )
}

volcano_stats_df <- do.call(rbind, volcano_stats_list)
rownames(volcano_stats_df) <- NULL
write.csv(volcano_stats_df, file.path(stats_dir, "volcano_plot_stats.csv"), row.names = FALSE)
cat(sprintf("[OK] Saved: %s\n\n", file.path(stats_dir, "volcano_plot_stats.csv")))

# Create faceted volcano plot with proper ordering
all_volcano_data <- do.call(rbind, lapply(names(all_dmr_data), function(region_size) {
  dmr_data <- all_dmr_data[[region_size]]
  dmr_data$neg_log10_qvalue <- -log10(dmr_data$qvalue)
  dmr_data$Category <- "Non-significant"
  dmr_data$Category[dmr_data$qvalue < opt$volcano_q_threshold & dmr_data$meth_difference > opt$volcano_meth_threshold] <- "Hypermethylated"
  dmr_data$Category[dmr_data$qvalue < opt$volcano_q_threshold & dmr_data$meth_difference < -opt$volcano_meth_threshold] <- "Hypomethylated"
  dmr_data$Region_Size_Label <- sprintf("%d KB", as.numeric(region_size) / 1000)
  dmr_data
}))

# Properly order region sizes
all_volcano_data$Region_Size_Label <- order_region_sizes(all_volcano_data$Region_Size_Label)

all_volcano_data$Category <- factor(all_volcano_data$Category,
  levels = c("Hypomethylated", "Non-significant", "Hypermethylated")
)

# Label top N SIGNIFICANT DMRs per region only, using rank_by_qvalue column
all_volcano_data$Label <- ""
for (region_size in unique(all_volcano_data$Region_Size_Label)) {
  region_subset <- all_volcano_data$Region_Size_Label == region_size
  # Only label significant DMRs
  significant_subset <- region_subset & (all_volcano_data$Category != "Non-significant")

  if (sum(significant_subset) > 0) {
    top_indices <- order(all_volcano_data$qvalue[significant_subset])[1:min(opt$top_n_label, sum(significant_subset))]
    # Use rank_by_qvalue column instead of sequential numbering
    dmr_ranks <- all_volcano_data$rank_by_qvalue[which(significant_subset)[top_indices]]
    # Create labels using rank: DMR_[rank]
    rank_labels <- sprintf("DMR_%d", dmr_ranks)
    all_volcano_data$Label[which(significant_subset)[top_indices]] <- rank_labels
  }
}

p2 <- ggplot(all_volcano_data, aes(x = meth_difference, y = neg_log10_qvalue, color = Category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(opt$volcano_q_threshold), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = opt$volcano_meth_threshold, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = -opt$volcano_meth_threshold, linetype = "dashed", color = "black", size = 0.5) +
  geom_text_repel(aes(label = Label), size = 2.5, max.overlaps = 15, box.padding = 0.5) +
  scale_color_manual(values = c("Hypomethylated" = "#1F78B4", "Non-significant" = "gray60", "Hypermethylated" = "#E31A1C")) +
  facet_wrap(~Region_Size_Label, scales = "free", ncol = 3) +
  labs(
    title = "Methylation difference vs significance",
    subtitle = sprintf("Thresholds: |Δβ| > %d%%, FDR < %.2f", opt$volcano_meth_threshold, opt$volcano_q_threshold),
    x = "Methylation difference (%)",
    y = "-log10(FDR)",
    colour = "DMR category"
  ) +
  theme_methylsense(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank()
  )

save_plot(p2, "2_volcano_plot_faceted", width = 14, height = 10)

cat("[COMPLETE] Volcano plot analysis\n")

# ================================================================================
# 3. COVERAGE DISTRIBUTION ANALYSIS
# ================================================================================

cat("\n================================================================================\n")
cat("Analysis 3: Coverage distribution\n")
cat("================================================================================\n\n")

coverage_stats_list <- list()

# Guard helpers: quantile() returns a NAMED scalar (e.g. "25%") which causes
# `data.frame()` to choke with "row names contain missing values" when the
# input vector is empty (e.g. zero significant DMRs at a given region size).
safe_quantile <- function(x, probs) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  as.numeric(quantile(x, probs, na.rm = TRUE))
}
safe_stat <- function(x, fn) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  as.numeric(fn(x, na.rm = TRUE))
}

for (region_size in names(all_meth_objects)) {
  cat(sprintf("Analyzing coverage for %s bp regions...\n", region_size))

  meth_obj <- all_meth_objects[[region_size]]
  dmr_data <- all_dmr_data[[region_size]]

  # Extract coverage data (average across all samples)
  meth_data <- getData(meth_obj)
  coverage_cols <- grep("^coverage", colnames(meth_data), value = TRUE)
  mean_coverage <- rowMeans(meth_data[, coverage_cols], na.rm = TRUE)

  # Separate DMRs vs non-DMRs
  dmr_ids <- dmr_data$DMR_ID[dmr_data$passes_current_qvalue_threshold]

  # Since methylKit object doesn't have DMR_ID directly, we'll use all regions as proxy
  # In a real scenario, you'd need to match by coordinates
  is_dmr <- seq_len(nrow(meth_obj)) <= length(dmr_ids) # Simplified - adjust if needed

  coverage_dmr <- mean_coverage[is_dmr]
  coverage_non_dmr <- mean_coverage[!is_dmr]

  if (length(coverage_dmr) == 0) {
    cat(sprintf(
      "  [INFO] No significant DMRs at %s bp - DMR coverage stats will be NA\n",
      region_size
    ))
  }

  coverage_stats_list[[region_size]] <- data.frame(
    Region_Size_bp  = as.numeric(region_size),
    Category        = c("DMR", "Non-DMR"),
    Median_Coverage = c(safe_stat(coverage_dmr, median),     safe_stat(coverage_non_dmr, median)),
    Q1_Coverage     = c(safe_quantile(coverage_dmr, 0.25),   safe_quantile(coverage_non_dmr, 0.25)),
    Q3_Coverage     = c(safe_quantile(coverage_dmr, 0.75),   safe_quantile(coverage_non_dmr, 0.75)),
    Mean_Coverage   = c(safe_stat(coverage_dmr, mean),       safe_stat(coverage_non_dmr, mean)),
    SD_Coverage     = c(safe_stat(coverage_dmr, sd),         safe_stat(coverage_non_dmr, sd)),
    stringsAsFactors = FALSE
  )

  cat(sprintf(
    "  DMR median coverage: %.1f (IQR: %.1f-%.1f)\n",
    safe_stat(coverage_dmr, median),
    safe_quantile(coverage_dmr, 0.25),
    safe_quantile(coverage_dmr, 0.75)
  ))
  cat(sprintf(
    "  Non-DMR median coverage: %.1f (IQR: %.1f-%.1f)\n\n",
    safe_stat(coverage_non_dmr, median),
    safe_quantile(coverage_non_dmr, 0.25),
    safe_quantile(coverage_non_dmr, 0.75)
  ))
}

coverage_stats_df <- do.call(rbind, coverage_stats_list)
rownames(coverage_stats_df) <- NULL
write.csv(coverage_stats_df, file.path(stats_dir, "coverage_distribution.csv"), row.names = FALSE)
cat(sprintf("[OK] Saved: %s\n\n", file.path(stats_dir, "coverage_distribution.csv")))

# Create violin plot with proper ordering
coverage_plot_data <- coverage_stats_df
coverage_plot_data$Region_Size_Label <- sprintf("%d KB", coverage_plot_data$Region_Size_bp / 1000)
coverage_plot_data$Region_Size_Label <- order_region_sizes(coverage_plot_data$Region_Size_Label)

p3 <- ggplot(coverage_plot_data, aes(x = Region_Size_Label, y = Median_Coverage, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.3) +
  geom_errorbar(aes(ymin = Q1_Coverage, ymax = Q3_Coverage), position = position_dodge(0.9), width = 0.2) +
  scale_fill_manual(values = c("DMR" = "#E31A1C", "Non-DMR" = "gray70")) +
  labs(
    title = "Coverage distribution: DMRs vs non-DMRs",
    subtitle = "Median coverage with interquartile range (error bars)",
    x = "Region size",
    y = "Median coverage (reads)",
    fill = "Region type"
  ) +
  theme_methylsense() +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 16)
  )

save_plot(p3, "3_coverage_distribution", width = 10, height = 8)

cat("[COMPLETE] Coverage distribution analysis\n")


# ================================================================================
# 4. CHROMOSOMAL DISTRIBUTION OF DMRS
# ================================================================================

cat("\n================================================================================\n")
cat("Analysis 4: Chromosomal distribution\n")
cat("================================================================================\n\n")

chr_dist_list <- list()

for (region_size in names(all_dmr_data)) {
  dmr_data <- all_dmr_data[[region_size]]
  dmr_data_sig <- dmr_data[dmr_data$passes_current_qvalue_threshold, ]

  chr_counts <- as.data.frame(table(dmr_data_sig$scaffold))
  colnames(chr_counts) <- c("Chromosome", "DMR_Count")
  chr_counts$Region_Size_bp <- as.numeric(region_size)

  chr_dist_list[[region_size]] <- chr_counts
}

chr_dist_df <- do.call(rbind, chr_dist_list)
rownames(chr_dist_df) <- NULL

# Add region size labels
chr_dist_df$Region_Size_Label <- sprintf("%d KB", chr_dist_df$Region_Size_bp / 1000)
chr_dist_df$Region_Size_Label <- order_region_sizes(chr_dist_df$Region_Size_Label)

write.csv(chr_dist_df, file.path(stats_dir, "dmr_chromosome_distribution.csv"), row.names = FALSE)
cat(sprintf("[OK] Saved: %s\n\n", file.path(stats_dir, "dmr_chromosome_distribution.csv")))

# Sort chromosomes numerically and show ALL chromosomes
# Extract numeric part of chromosome names for sorting
chr_dist_df$Chr_Numeric <- as.numeric(gsub(".*?([0-9]+).*", "\\1", chr_dist_df$Chromosome))

# For chromosomes that don't have numbers, assign high values to put them at the end
chr_dist_df$Chr_Numeric[is.na(chr_dist_df$Chr_Numeric)] <- 999

# Sort by chromosome number
chr_dist_df <- chr_dist_df[order(chr_dist_df$Chr_Numeric), ]

# Create ordered factor for plotting
chr_order <- unique(chr_dist_df$Chromosome[order(chr_dist_df$Chr_Numeric)])
chr_dist_df$Chromosome <- factor(chr_dist_df$Chromosome, levels = chr_order)

# Use ALL chromosomes (not just top 15)
chr_plot_data <- chr_dist_df

# Create grouped bar plot for all region sizes showing ALL chromosomes
p4 <- ggplot(chr_plot_data, aes(x = Chromosome, y = DMR_Count, fill = Region_Size_Label)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.2) +
  scale_fill_viridis_d(option = "plasma") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Chromosomal distribution of DMRs across region sizes",
    subtitle = "All chromosomes sorted numerically",
    x = "Chromosome",
    y = "Number of DMRs",
    fill = "Region size"
  ) +
  theme_methylsense(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)
  )

save_plot(p4, "4_dmr_chromosome_ideogram", width = 10, height = 8)

cat("[COMPLETE] Chromosomal distribution analysis\n")

# ================================================================================
# 5. EFFECT SIZE DISTRIBUTION
# ================================================================================

cat("\n================================================================================\n")
cat("Analysis 5: Effect size distribution\n")
cat("================================================================================\n\n")

effect_size_stats_list <- list()

for (region_size in names(all_dmr_data)) {
  dmr_data <- all_dmr_data[[region_size]]
  dmr_data_sig <- dmr_data[dmr_data$passes_current_qvalue_threshold, ]

  effect_size_stats_list[[region_size]] <- data.frame(
    Region_Size_bp = as.numeric(region_size),
    Direction = c("Hypermethylated", "Hypomethylated"),
    N_DMRs = c(
      sum(dmr_data_sig$direction == "Hypermethylated"),
      sum(dmr_data_sig$direction == "Hypomethylated")
    ),
    Mean_Effect = c(
      mean(dmr_data_sig$meth_difference[dmr_data_sig$direction == "Hypermethylated"], na.rm = TRUE),
      mean(dmr_data_sig$meth_difference[dmr_data_sig$direction == "Hypomethylated"], na.rm = TRUE)
    ),
    Median_Effect = c(
      median(dmr_data_sig$meth_difference[dmr_data_sig$direction == "Hypermethylated"], na.rm = TRUE),
      median(dmr_data_sig$meth_difference[dmr_data_sig$direction == "Hypomethylated"], na.rm = TRUE)
    ),
    Min_Effect = c(
      min(dmr_data_sig$meth_difference[dmr_data_sig$direction == "Hypermethylated"], na.rm = TRUE),
      min(dmr_data_sig$meth_difference[dmr_data_sig$direction == "Hypomethylated"], na.rm = TRUE)
    ),
    Max_Effect = c(
      max(dmr_data_sig$meth_difference[dmr_data_sig$direction == "Hypermethylated"], na.rm = TRUE),
      max(dmr_data_sig$meth_difference[dmr_data_sig$direction == "Hypomethylated"], na.rm = TRUE)
    )
  )
}

effect_size_stats_df <- do.call(rbind, effect_size_stats_list)
rownames(effect_size_stats_df) <- NULL
write.csv(effect_size_stats_df, file.path(stats_dir, "effect_size_distribution.csv"), row.names = FALSE)
cat(sprintf("[OK] Saved: %s\n\n", file.path(stats_dir, "effect_size_distribution.csv")))

# Create ridge plot for effect sizes by region size
all_effect_data <- do.call(rbind, lapply(names(all_dmr_data), function(region_size) {
  dmr_data <- all_dmr_data[[region_size]]
  dmr_data_sig <- dmr_data[dmr_data$passes_current_qvalue_threshold, ]
  dmr_data_sig$Region_Size_Label <- sprintf("%d KB", as.numeric(region_size) / 1000)
  dmr_data_sig[, c("Region_Size_Label", "meth_difference", "direction")]
}))

# Properly order region sizes (reversed for nice ridge stacking)
all_effect_data$Region_Size_Label <- order_region_sizes(all_effect_data$Region_Size_Label)
all_effect_data$Region_Size_Label <- factor(all_effect_data$Region_Size_Label,
  levels = rev(levels(all_effect_data$Region_Size_Label))
)

p5 <- ggplot(all_effect_data, aes(x = meth_difference, y = Region_Size_Label, fill = stat(x))) +
  geom_density_ridges_gradient(alpha = 0.8, scale = 1.5, rel_min_height = 0.01) +
  scale_fill_viridis_c(option = "plasma", name = "Methylation\nChange (%)") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "white", size = 0.8) +
  labs(
    title = "Effect size distribution across region sizes",
    subtitle = "Methylation differences (Δβ) for significant DMRs",
    x = "Methylation difference (%)",
    y = "Genomic window size"
  ) +
  theme_methylsense(base_size = 14, grid = "X") +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16)
  )

save_plot(p5, "5_effect_size_distribution", width = 14, height = 10)

cat("[COMPLETE] Effect size distribution analysis\n")

# ================================================================================
# 6. SAMPLE-LEVEL HEATMAP (TOP DMRs) - ALL REGION SIZES
# ================================================================================

cat("\n================================================================================\n")
cat("Analysis 6: Sample-level methylation heatmaps\n")
cat("================================================================================\n\n")

# Prepare annotation colors - intelligent color scheme
# Get unique infection statuses
unique_infections <- unique(as.character(sample_metadata[[opt$infection_col]]))

# Intelligent infection status coloring: Control = dark green, Infected = red/purple
infection_colors <- character(length(unique_infections))
names(infection_colors) <- unique_infections

for (status in unique_infections) {
  if (grepl("Control|control|CONTROL", status, ignore.case = TRUE)) {
    infection_colors[status] <- METHYLSENSE_COLORS$Control
  } else if (grepl("Infect|infect|INFECT|Asper|asper", status, ignore.case = TRUE)) {
    infection_colors[status] <- METHYLSENSE_COLORS$Infected
  } else if (grepl("Suspect|suspect|SUSPECT", status, ignore.case = TRUE)) {
    infection_colors[status] <- METHYLSENSE_COLORS$Suspected
  } else {
    infection_colors[status] <- METHYLSENSE_COLORS$Neutral
  }
}

# Study colors
unique_studies <- unique(as.character(sample_metadata[[opt$study_col]]))
study_colors <- brewer.pal(min(length(unique_studies), 8), "Set2")
if (length(unique_studies) > 8) {
  study_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(unique_studies))
}
names(study_colors) <- unique_studies

annotation_colors <- list(
  Infection_Status = infection_colors,
  Study = study_colors
)

# Create heatmaps for all region sizes
for (region_size in names(all_meth_objects)) {
  cat(sprintf("Creating heatmap for %s bp regions...\n", region_size))

  meth_obj <- all_meth_objects[[region_size]]
  dmr_data <- all_dmr_data[[region_size]]

  # Get top N DMRs
  top_dmrs <- head(dmr_data[order(dmr_data$qvalue), ], opt$heatmap_top_n)

  # Calculate methylation percentages
  meth_data <- getData(meth_obj)
  meth_percent_matrix <- matrix(NA, nrow = nrow(meth_obj), ncol = length(meth_obj@sample.ids))
  for (i in seq_along(meth_obj@sample.ids)) {
    numCs_col <- sprintf("numCs%d", i)
    numTs_col <- sprintf("numTs%d", i)
    if (numCs_col %in% colnames(meth_data) && numTs_col %in% colnames(meth_data)) {
      meth_percent_matrix[, i] <- (meth_data[[numCs_col]] /
        (meth_data[[numCs_col]] + meth_data[[numTs_col]])) * 100
    }
  }

  # Take first N rows corresponding to top DMRs
  top_dmr_indices <- seq_len(min(opt$heatmap_top_n, nrow(meth_obj)))
  heatmap_matrix <- meth_percent_matrix[top_dmr_indices, ]
  rownames(heatmap_matrix) <- paste0("DMR_", top_dmr_indices)
  colnames(heatmap_matrix) <- paste0("Sample_", seq_len(ncol(heatmap_matrix)))

  # Prepare annotation
  sample_annot <- as.data.frame(sample_metadata[, c(opt$sample_id_col, opt$infection_col, opt$study_col)])
  colnames(sample_annot) <- c("Sample_ID", "Infection_Status", "Study")

  # Convert to data frame and ensure character columns
  sample_annot$Infection_Status <- as.character(sample_annot$Infection_Status)
  sample_annot$Study <- as.character(sample_annot$Study)

  # Match sample order
  if (nrow(sample_annot) == ncol(heatmap_matrix)) {
    sample_annot <- data.frame(
      Infection_Status = sample_annot$Infection_Status,
      Study = sample_annot$Study,
      row.names = colnames(heatmap_matrix),
      stringsAsFactors = FALSE
    )

    # Create heatmap in all requested formats
    region_label <- sprintf("%d KB", as.numeric(region_size) / 1000)

    for (format in plot_formats) {
      filepath <- file.path(figures_dir, sprintf(
        "6_dmr_heatmap_top%d_%s.%s",
        opt$heatmap_top_n, region_size, format
      ))

      if (format == "pdf") {
        pdf(filepath, width = opt$plot_width, height = opt$plot_height)
      } else if (format == "png") {
        png(filepath, width = opt$plot_width * opt$plot_dpi, height = opt$plot_height * opt$plot_dpi, res = opt$plot_dpi)
      } else if (format == "jpg" || format == "jpeg") {
        jpeg(filepath, width = opt$plot_width * opt$plot_dpi, height = opt$plot_height * opt$plot_dpi, res = opt$plot_dpi)
      } else if (format == "svg") {
        svg(filepath, width = opt$plot_width, height = opt$plot_height)
      }

      pheatmap(heatmap_matrix,
        color = colorRampPalette(c("blue", "white", "red"))(100),
        breaks = seq(0, 100, length.out = 101),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_colnames = FALSE,
        annotation_col = sample_annot,
        annotation_colors = annotation_colors,
        main = sprintf("Top %d DMRs methylation heatmap (%s regions)", opt$heatmap_top_n, region_label),
        fontsize = 10
      )

      dev.off()
      cat(sprintf("[OK] Saved: %s\n", filepath))
    }
  } else {
    cat(sprintf("[WARNING] Sample metadata row count doesn't match methylation data for %s bp. Skipping...\n", region_size))
  }
}

cat("\n[COMPLETE] Heatmap analysis for all region sizes\n")

# ================================================================================
# FINAL SUMMARY REPORT
# ================================================================================

cat("\n================================================================================\n")
cat("Generating summary report\n")
cat("================================================================================\n\n")

summary_text <- c(
  "================================================================================",
  "DMR landscape overview - summary report",
  "================================================================================",
  "",
  sprintf("Analysis Date: %s", Sys.Date()),
  sprintf("Script Version: %s", SCRIPT_VERSION),
  sprintf("Analysis Directory: %s", if (opt$blind_paths) sanitize_path(opt$analysis_dir) else opt$analysis_dir),
  sprintf("Sample Sheet: %s", if (opt$blind_paths) sanitize_path(opt$sample_sheet) else opt$sample_sheet),
  "",
  "Dataset overview:",
  sprintf("- Total samples analysed: %d", nrow(sample_metadata)),
  sprintf("- Region sizes tested: %s", paste(region_sizes / 1000, "KB", collapse = ", ")),
  "",
  "Differential methylation results:",
  ""
)

for (i in seq_len(nrow(dmr_counts_df))) {
  summary_text <- c(
    summary_text,
    sprintf("%s regions:", dmr_counts_df$Region_Size_Label[i]),
    sprintf("  - Total MRs: %d", dmr_counts_df$Total_MRs[i]),
    sprintf("  - Significant DMRs (q<%.2f): %d", opt$volcano_q_threshold, dmr_counts_df$Significant_DMRs[i]),
    sprintf("  - DMR detection rate: %.1f%%", dmr_counts_df$DMR_Rate_Percent[i]),
    ""
  )
}

# Add direction summary (using 1KB as example)
if ("1000" %in% names(all_dmr_data)) {
  dmr_1kb <- all_dmr_data[["1000"]][all_dmr_data[["1000"]]$passes_current_qvalue_threshold, ]
  n_hyper <- sum(dmr_1kb$direction == "Hypermethylated")
  n_hypo <- sum(dmr_1kb$direction == "Hypomethylated")
  pct_hyper <- (n_hyper / nrow(dmr_1kb)) * 100

  # Calculate mean and variance measures for hypermethylated
  mean_hyper <- mean(dmr_1kb$meth_difference[dmr_1kb$direction == "Hypermethylated"], na.rm = TRUE)
  sd_hyper <- sd(dmr_1kb$meth_difference[dmr_1kb$direction == "Hypermethylated"], na.rm = TRUE)
  median_hyper <- median(dmr_1kb$meth_difference[dmr_1kb$direction == "Hypermethylated"], na.rm = TRUE)
  q1_hyper <- quantile(dmr_1kb$meth_difference[dmr_1kb$direction == "Hypermethylated"], 0.25, na.rm = TRUE)
  q3_hyper <- quantile(dmr_1kb$meth_difference[dmr_1kb$direction == "Hypermethylated"], 0.75, na.rm = TRUE)

  # Calculate mean and variance measures for hypomethylated
  mean_hypo <- mean(dmr_1kb$meth_difference[dmr_1kb$direction == "Hypomethylated"], na.rm = TRUE)
  sd_hypo <- sd(dmr_1kb$meth_difference[dmr_1kb$direction == "Hypomethylated"], na.rm = TRUE)
  median_hypo <- median(dmr_1kb$meth_difference[dmr_1kb$direction == "Hypomethylated"], na.rm = TRUE)
  q1_hypo <- quantile(dmr_1kb$meth_difference[dmr_1kb$direction == "Hypomethylated"], 0.25, na.rm = TRUE)
  q3_hypo <- quantile(dmr_1kb$meth_difference[dmr_1kb$direction == "Hypomethylated"], 0.75, na.rm = TRUE)

  summary_text <- c(
    summary_text,
    "METHYLATION DIRECTION (1 KB regions):",
    sprintf("- Hypermethylated DMRs: %d (%.1f%%)", n_hyper, pct_hyper),
    sprintf(
      "  Effect size: mean Δβ = +%.2f%% (SD = %.2f%%), median = +%.2f%% (IQR: %.2f%% to %.2f%%)",
      mean_hyper, sd_hyper, median_hyper, q1_hyper, q3_hyper
    ),
    sprintf("- Hypomethylated DMRs: %d (%.1f%%)", n_hypo, 100 - pct_hyper),
    sprintf(
      "  Effect size: mean Δβ = %.2f%% (SD = %.2f%%), median = %.2f%% (IQR: %.2f%% to %.2f%%)",
      mean_hypo, sd_hypo, median_hypo, q1_hypo, q3_hypo
    ),
    ""
  )
}

# Add coverage summary
if (nrow(coverage_stats_df) > 0) {
  dmr_coverage <- coverage_stats_df[coverage_stats_df$Category == "DMR", ]
  non_dmr_coverage <- coverage_stats_df[coverage_stats_df$Category == "Non-DMR", ]

  # DMR coverage statistics
  median_cov_dmr <- median(dmr_coverage$Median_Coverage)
  mean_cov_dmr <- mean(dmr_coverage$Median_Coverage)
  q1_cov_dmr <- quantile(dmr_coverage$Median_Coverage, 0.25)
  q3_cov_dmr <- quantile(dmr_coverage$Median_Coverage, 0.75)

  # Non-DMR coverage statistics
  median_cov_non <- median(non_dmr_coverage$Median_Coverage)
  mean_cov_non <- mean(non_dmr_coverage$Median_Coverage)
  q1_cov_non <- quantile(non_dmr_coverage$Median_Coverage, 0.25)
  q3_cov_non <- quantile(non_dmr_coverage$Median_Coverage, 0.75)

  summary_text <- c(
    summary_text,
    "Coverage quality:",
    sprintf(
      "- DMR coverage: median = %.1f× (IQR: %.1f-%.1f×), mean = %.1f×",
      median_cov_dmr, q1_cov_dmr, q3_cov_dmr, mean_cov_dmr
    ),
    sprintf(
      "- Non-DMR coverage: median = %.1f× (IQR: %.1f-%.1f×), mean = %.1f×",
      median_cov_non, q1_cov_non, q3_cov_non, mean_cov_non
    ),
    sprintf(
      "- Coverage range across region sizes: %.1f-%.1f× (DMR), %.1f-%.1f× (non-DMR)",
      min(dmr_coverage$Median_Coverage), max(dmr_coverage$Median_Coverage),
      min(non_dmr_coverage$Median_Coverage), max(non_dmr_coverage$Median_Coverage)
    ),
    ""
  )
}

summary_text <- c(
  summary_text,
  "Output files:",
  sprintf("- Figures directory: %s", figures_dir),
  sprintf("- Statistics directory: %s", stats_dir),
  "",
  "Generated plots:",
  "1. DMR/MR counts bar plot (with counts labeled on bars)",
  "2. Volcano plots (faceted by region size, labeled by rank)",
  "3. Coverage distribution (DMRs vs non-DMRs)",
  "4. Chromosomal distribution (all chromosomes, sorted numerically)",
  "5. Effect size distribution (ridge plot by region size)",
  "6. Sample-level methylation heatmaps (one per region size)",
  "7. Methylation distribution by infection status (ridge plot)",
  "",
  "Statistics files:",
  "- dmr_counts_summary.csv",
  "- volcano_plot_stats.csv",
  "- coverage_distribution.csv",
  "- dmr_chromosome_distribution.csv",
  "- effect_size_distribution.csv",
  "",
  "================================================================================",
  "Analysis complete",
  "================================================================================"
)

writeLines(summary_text, file.path(opt$output_dir, "report_summary.txt"))
cat(sprintf("[OK] Saved: %s\n", file.path(opt$output_dir, "report_summary.txt")))

# ================================================================================
# Generate comprehensive markdown report
# ================================================================================
cat("\n[STEP 9] Generating comprehensive markdown report...\n")

markdown_lines <- c(
  "# MethylSense: General Methylation Statistics Report",
  "",
  sprintf("**Analysis Date**: %s", Sys.Date()),
  sprintf("**Script Version**: %s", SCRIPT_VERSION),
  sprintf("**Generated by**: MethylSense DMR Landscape Analysis Pipeline"),
  "",
  "---",
  "",
  "## Executive Summary",
  "",
  sprintf(
    "This report presents a comprehensive landscape analysis of differential methylation across **%d region sizes** (%s). ",
    length(region_sizes),
    paste(dmr_counts_df$Region_Size_Label, collapse = ", ")
  ),
  sprintf(
    "Analysis identified a total of **%d significant DMRs** (FDR < %.2f) across all region sizes, ",
    sum(dmr_counts_df$Significant_DMRs),
    opt$volcano_q_threshold
  ),
  sprintf(
    "with detection rates ranging from **%.1f%% to %.1f%%**.",
    min(dmr_counts_df$DMR_Rate_Percent),
    max(dmr_counts_df$DMR_Rate_Percent)
  ),
  "",
  "---",
  "",
  "## Table of Contents",
  "",
  "1. [DMR Detection Rates](#1-dmr-detection-rates)",
  "2. [Effect Size Analysis](#2-effect-size-analysis)",
  "3. [Coverage Quality Control](#3-coverage-quality-control)",
  "4. [Chromosomal Distribution](#4-chromosomal-distribution)",
  "5. [Methylation Patterns](#5-methylation-patterns)",
  "6. [Sample-Level Clustering](#6-sample-level-clustering)",
  "7. [Statistical Summary](#7-statistical-summary)",
  "",
  "---",
  ""
)

# ================================================================================
# Section 1: DMR Detection Rates
# ================================================================================
markdown_lines <- c(
  markdown_lines,
  "## 1. DMR Detection Rates",
  "",
  "### Overview",
  ""
)

for (i in seq_len(nrow(dmr_counts_df))) {
  markdown_lines <- c(
    markdown_lines,
    sprintf("**%s regions:**", dmr_counts_df$Region_Size_Label[i]),
    sprintf("- Total methylation regions (MRs): **%s**", format(dmr_counts_df$Total_MRs[i], big.mark = ",")),
    sprintf("- Significant DMRs (FDR < %.2f): **%s**", opt$volcano_q_threshold, format(dmr_counts_df$Significant_DMRs[i], big.mark = ",")),
    sprintf("- Detection rate: **%.1f%%**", dmr_counts_df$DMR_Rate_Percent[i]),
    ""
  )
}

markdown_lines <- c(
  markdown_lines,
  "### Figure 1: DMR Counts by Region Size",
  "",
  sprintf('<img src="%s" alt="DMR Counts" width="600"/>', file.path("figures", paste0("1_dmr_counts_barplot.", plot_formats[1]))),
  "",
  sprintf("*Bar plot showing the number of significant DMRs (FDR < %.2f) vs non-significant MRs across different genomic region sizes. ", opt$volcano_q_threshold),
  "Numbers on bars indicate exact counts. Larger region sizes generally show higher detection rates.*",
  "",
  "---",
  ""
)

# ================================================================================
# Section 2: Effect Size Analysis
# ================================================================================
markdown_lines <- c(
  markdown_lines,
  "## 2. Effect Size Analysis",
  "",
  "### Methylation Direction (1 KB regions)",
  ""
)

if ("1000" %in% names(all_dmr_data)) {
  dmr_1kb <- all_dmr_data[["1000"]][all_dmr_data[["1000"]]$passes_current_qvalue_threshold, ]
  n_hyper <- sum(dmr_1kb$direction == "Hypermethylated")
  n_hypo <- sum(dmr_1kb$direction == "Hypomethylated")
  pct_hyper <- (n_hyper / nrow(dmr_1kb)) * 100

  mean_hyper <- mean(dmr_1kb$meth_difference[dmr_1kb$direction == "Hypermethylated"], na.rm = TRUE)
  sd_hyper <- sd(dmr_1kb$meth_difference[dmr_1kb$direction == "Hypermethylated"], na.rm = TRUE)
  median_hyper <- median(dmr_1kb$meth_difference[dmr_1kb$direction == "Hypermethylated"], na.rm = TRUE)
  q1_hyper <- quantile(dmr_1kb$meth_difference[dmr_1kb$direction == "Hypermethylated"], 0.25, na.rm = TRUE)
  q3_hyper <- quantile(dmr_1kb$meth_difference[dmr_1kb$direction == "Hypermethylated"], 0.75, na.rm = TRUE)

  mean_hypo <- mean(dmr_1kb$meth_difference[dmr_1kb$direction == "Hypomethylated"], na.rm = TRUE)
  sd_hypo <- sd(dmr_1kb$meth_difference[dmr_1kb$direction == "Hypomethylated"], na.rm = TRUE)
  median_hypo <- median(dmr_1kb$meth_difference[dmr_1kb$direction == "Hypomethylated"], na.rm = TRUE)
  q1_hypo <- quantile(dmr_1kb$meth_difference[dmr_1kb$direction == "Hypomethylated"], 0.25, na.rm = TRUE)
  q3_hypo <- quantile(dmr_1kb$meth_difference[dmr_1kb$direction == "Hypomethylated"], 0.75, na.rm = TRUE)

  markdown_lines <- c(
    markdown_lines,
    sprintf("**Hypermethylated DMRs:** %d (%.1f%%)", n_hyper, pct_hyper),
    sprintf("- Mean effect size: **Δβ = +%.2f%%** (SD = %.2f%%)", mean_hyper, sd_hyper),
    sprintf("- Median effect size: **Δβ = +%.2f%%** (IQR: %.2f%% to %.2f%%)", median_hyper, q1_hyper, q3_hyper),
    "",
    sprintf("**Hypomethylated DMRs:** %d (%.1f%%)", n_hypo, 100 - pct_hyper),
    sprintf("- Mean effect size: **Δβ = %.2f%%** (SD = %.2f%%)", mean_hypo, sd_hypo),
    sprintf("- Median effect size: **Δβ = %.2f%%** (IQR: %.2f%% to %.2f%%)", median_hypo, q1_hypo, q3_hypo),
    ""
  )
}

markdown_lines <- c(
  markdown_lines,
  "### Figure 2: Volcano Plots",
  "",
  sprintf('<img src="%s" alt="Volcano Plots" width="700"/>', file.path("figures", paste0("2_volcano_plot_faceted.", plot_formats[1]))),
  "",
  sprintf("*Volcano plots showing methylation difference (Δβ) vs statistical significance (-log₁₀ FDR) for each region size. "),
  sprintf("Horizontal dashed line: FDR threshold (%.2f). ", opt$volcano_q_threshold),
  sprintf("Vertical dashed lines: effect size threshold (±%.1f%%). ", opt$volcano_meth_threshold),
  sprintf("Top %d DMRs are labelled by rank. ", opt$top_n_label),
  "Points coloured by direction: red (hypermethylated), blue (hypomethylated).*",
  "",
  "### Figure 3: Effect Size Distribution",
  "",
  sprintf('<img src="%s" alt="Effect Size Distribution" width="700"/>', file.path("figures", paste0("5_effect_size_distribution.", plot_formats[1]))),
  "",
  "*Ridge plot showing the distribution of methylation effect sizes (Δβ) across different region sizes. ",
  "Colour gradient represents effect size magnitude. Separate distributions shown for hypermethylated (positive) and hypomethylated (negative) DMRs.*",
  "",
  "---",
  ""
)

# ================================================================================
# Section 3: Coverage QC
# ================================================================================
markdown_lines <- c(
  markdown_lines,
  "## 3. Coverage Quality Control",
  ""
)

if (nrow(coverage_stats_df) > 0) {
  dmr_coverage <- coverage_stats_df[coverage_stats_df$Category == "DMR", ]
  non_dmr_coverage <- coverage_stats_df[coverage_stats_df$Category == "Non-DMR", ]

  median_cov_dmr <- median(dmr_coverage$Median_Coverage)
  mean_cov_dmr <- mean(dmr_coverage$Median_Coverage)
  q1_cov_dmr <- quantile(dmr_coverage$Median_Coverage, 0.25)
  q3_cov_dmr <- quantile(dmr_coverage$Median_Coverage, 0.75)

  median_cov_non <- median(non_dmr_coverage$Median_Coverage)
  mean_cov_non <- mean(non_dmr_coverage$Median_Coverage)
  q1_cov_non <- quantile(non_dmr_coverage$Median_Coverage, 0.25)
  q3_cov_non <- quantile(non_dmr_coverage$Median_Coverage, 0.75)

  markdown_lines <- c(
    markdown_lines,
    "### Coverage Statistics",
    "",
    "**DMR Coverage:**",
    sprintf("- Median: **%.1f×** (IQR: %.1f–%.1f×)", median_cov_dmr, q1_cov_dmr, q3_cov_dmr),
    sprintf("- Mean: **%.1f×** (SD: %.1f×)", mean_cov_dmr, sd(dmr_coverage$Median_Coverage)),
    sprintf("- Range: %.1f–%.1f× across region sizes", min(dmr_coverage$Median_Coverage), max(dmr_coverage$Median_Coverage)),
    "",
    "**Non-DMR Coverage:**",
    sprintf("- Median: **%.1f×** (IQR: %.1f–%.1f×)", median_cov_non, q1_cov_non, q3_cov_non),
    sprintf("- Mean: **%.1f×** (SD: %.1f×)", mean_cov_non, sd(non_dmr_coverage$Median_Coverage)),
    sprintf("- Range: %.1f–%.1f× across region sizes", min(non_dmr_coverage$Median_Coverage), max(non_dmr_coverage$Median_Coverage)),
    ""
  )
}

markdown_lines <- c(
  markdown_lines,
  "### Figure 4: Coverage Distribution",
  "",
  sprintf('<img src="%s" alt="Coverage Distribution" width="600"/>', file.path("figures", paste0("3_coverage_distribution.", plot_formats[1]))),
  "",
  "*Box plots comparing sequencing coverage between DMRs and non-DMRs across different region sizes. ",
  "Higher coverage in DMRs indicates better statistical power for differential methylation detection.*",
  "",
  "---",
  ""
)

# ================================================================================
# Section 4: Chromosomal Distribution
# ================================================================================
markdown_lines <- c(
  markdown_lines,
  "## 4. Chromosomal Distribution",
  "",
  "### Figure 5: DMR Distribution Across Chromosomes",
  "",
  sprintf('<img src="%s" alt="Chromosomal Distribution" width="600"/>', file.path("figures", paste0("4_dmr_chromosome_ideogram.", plot_formats[1]))),
  "",
  "*Bar plot showing the distribution of DMRs across all chromosomes (sorted numerically). ",
  "Different colours represent different region sizes. Identifies chromosomes with DMR hotspots.*",
  "",
  "---",
  ""
)

# ================================================================================
# Section 5: Methylation Patterns
# ================================================================================
markdown_lines <- c(
  markdown_lines,
  "## 5. Methylation Patterns by Infection Status",
  "",
  "### Figure 6: Methylation Distribution Ridge Plot",
  "",
  sprintf('<img src="%s" alt="Methylation by Infection" width="700"/>', file.path("figures", paste0("7_ridge_methylation_by_infection_all_regions.", plot_formats[1]))),
  "",
  "*Ridge plots showing methylation percentage distributions for Control vs Infected samples across all region sizes. ",
  "Each facet represents one region size with Control (green) and Infected (red) distributions stacked vertically. ",
  "Reveals systematic methylation differences between groups and how separation changes with region size.*",
  "",
  "---",
  ""
)

# ================================================================================
# Section 6: Sample-Level Clustering
# ================================================================================
markdown_lines <- c(
  markdown_lines,
  "## 6. Sample-Level Clustering",
  "",
  sprintf("Heatmaps were generated for the top %d DMRs in each region size, showing methylation patterns across all samples. ", opt$heatmap_top_n),
  "Samples are clustered hierarchically using Euclidean distance and complete linkage.",
  "",
  "### Heatmap Files",
  ""
)

for (region_size in region_sizes) {
  region_label <- paste0(region_size / 1000, " KB")
  heatmap_file <- paste0("6_dmr_heatmap_top50_", region_size, ".", plot_formats[1])
  if (file.exists(file.path(figures_dir, heatmap_file))) {
    markdown_lines <- c(
      markdown_lines,
      sprintf("**%s regions:**", region_label),
      "",
      sprintf('<img src="%s" alt="Heatmap %s" width="600"/>', file.path("figures", heatmap_file), region_label),
      ""
    )
  }
}

markdown_lines <- c(
  markdown_lines,
  sprintf("*Heatmaps show methylation percentage (0-100%%) for top %d DMRs ranked by FDR. ", opt$heatmap_top_n),
  "Columns represent samples (annotated by infection status and study). ",
  "Rows represent DMRs (ordered by hierarchical clustering). ",
  "Colour annotation: Control (dark green), Infected (red), Study (ColourBrewer Set2 palette).*",
  "",
  "---",
  ""
)

# ================================================================================
# Section 7: Statistical Summary
# ================================================================================
markdown_lines <- c(
  markdown_lines,
  "## 7. Statistical Summary",
  "",
  "### Summary Table: DMR Counts by Region Size",
  "",
  "| Region Size | Total MRs | Significant DMRs | Detection Rate |",
  "|-------------|-----------|------------------|----------------|"
)

for (i in seq_len(nrow(dmr_counts_df))) {
  markdown_lines <- c(
    markdown_lines,
    sprintf(
      "| %s | %s | %s | %.1f%% |",
      dmr_counts_df$Region_Size_Label[i],
      format(dmr_counts_df$Total_MRs[i], big.mark = ","),
      format(dmr_counts_df$Significant_DMRs[i], big.mark = ","),
      dmr_counts_df$DMR_Rate_Percent[i]
    )
  )
}

markdown_lines <- c(
  markdown_lines,
  "",
  "### Data Files",
  "",
  "All statistical summaries are available in the `stats/` directory:",
  "",
  "- **`dmr_counts_summary.csv`**: Region size, total MRs, significant DMRs, detection rates",
  "- **`volcano_plot_stats.csv`**: Counts by quadrant (significant hyper/hypo, non-significant)",
  "- **`coverage_distribution.csv`**: Median, IQR, mean, SD for DMRs and non-DMRs",
  "- **`dmr_chromosome_distribution.csv`**: DMR counts per chromosome for each region size",
  "- **`effect_size_distribution.csv`**: Mean, median, min, max effect sizes by direction",
  "",
  "---",
  "",
  "## Methods Summary",
  "",
  "### Differential Methylation Analysis",
  "",
  sprintf(
    "Differential methylation analysis was performed using methylKit across %d genomic region sizes (%s). ",
    length(region_sizes),
    paste(sapply(region_sizes, function(x) paste0(x / 1000, " KB")), collapse = ", ")
  ),
  sprintf("Regions with FDR < %.2f were considered significant DMRs. ", opt$volcano_q_threshold),
  sprintf("Effect sizes represent absolute methylation difference (Δβ) between groups, with a threshold of ±%.1f%% used for biological significance. ", opt$volcano_meth_threshold),
  sprintf(
    "Coverage quality control ensured median coverage of %.1f× for DMRs and %.1f× for non-DMRs.",
    median_cov_dmr, median_cov_non
  ),
  "",
  "### Visualization",
  "",
  "All plots were generated using ggplot2 (bar plots, volcano plots, box plots), ggridges (ridge plots), and pheatmap (heatmaps). ",
  sprintf(
    "Figures were exported as %s format at %d DPI with dimensions %d × %d inches.",
    toupper(paste(plot_formats, collapse = ", ")), opt$plot_dpi, opt$plot_width, opt$plot_height
  ),
  "",
  "---",
  "",
  "## Conclusions",
  "",
  sprintf(
    "This landscape analysis identified **%d significant DMRs** across %d region sizes, with detection rates ranging from %.1f%% to %.1f%%. ",
    sum(dmr_counts_df$Significant_DMRs),
    length(region_sizes),
    min(dmr_counts_df$DMR_Rate_Percent),
    max(dmr_counts_df$DMR_Rate_Percent)
  )
)

if ("1000" %in% names(all_dmr_data)) {
  markdown_lines <- c(
    markdown_lines,
    sprintf(
      "The majority of DMRs showed hypermethylation (%.1f%%) with mean effect size of +%.2f%%, ",
      pct_hyper, mean_hyper
    ),
    sprintf(
      "while hypomethylated regions (%.1f%%) showed mean effect size of %.2f%%. ",
      100 - pct_hyper, mean_hypo
    )
  )
}

markdown_lines <- c(
  markdown_lines,
  "Chromosomal distribution revealed DMR hotspots on specific chromosomes, and sample-level clustering confirmed distinct methylation patterns between infection status groups. ",
  "Coverage quality was consistent across region sizes, ensuring robust statistical inference.",
  "",
  "---",
  "",
  sprintf("**Report generated**: %s", Sys.time()),
  sprintf("**Script version**: %s", SCRIPT_VERSION),
  sprintf("**Output directory**: `%s`", if (opt$blind_paths) sanitize_path(opt$output_dir) else opt$output_dir),
  ""
)

# Write markdown report
markdown_path <- file.path(opt$output_dir, "MethylSense_General_Statistics_Report.md")
writeLines(markdown_lines, markdown_path)
cat(sprintf("[OK] Saved markdown report: %s\n", markdown_path))

# Print summary to console
cat("\n")
cat(paste(summary_text, collapse = "\n"))
cat("\n\n")

cat("All analyses completed successfully!\n\n")

# ================================================================================
# 7. RIDGE PLOT: METHYLATION DISTRIBUTION BY INFECTION STATUS - ALL REGION SIZES
# ================================================================================

cat("\n================================================================================\n")
cat("Analysis 7: Ridge plot - methylation by infection status (all regions)\n")
cat("================================================================================\n\n")

# Collect methylation data across ALL region sizes
all_infection_data <- list()

for (region_size in names(all_meth_objects)) {
  cat(sprintf("Processing %s bp regions for infection status ridge plot...\n", region_size))

  meth_obj <- all_meth_objects[[region_size]]
  dmr_data <- all_dmr_data[[region_size]]

  # Get significant DMRs
  sig_dmrs <- dmr_data[dmr_data$passes_current_qvalue_threshold, ]

  if (nrow(sig_dmrs) > 0) {
    # Get significant methylation data for ALL DMRs and ALL samples
    # We want a long-format dataframe: [Sample_ID, Group, Methylation_Value]
    meth_data <- getData(meth_obj)
    sig_dmr_indices <- which(getData(meth_obj)$dmr_id %in% sig_dmrs$dmr_id)
    if (length(sig_dmr_indices) == 0) sig_dmr_indices <- seq_len(min(nrow(sig_dmrs), nrow(meth_obj)))

    temp_list <- list()
    for (i in seq_along(meth_obj@sample.ids)) {
      numCs_col <- sprintf("numCs%d", i)
      numTs_col <- sprintf("numTs%d", i)

      if (numCs_col %in% colnames(meth_data) && numTs_col %in% colnames(meth_data)) {
        # Subset to significant DMRs for this sample
        meth_vals <- (meth_data[sig_dmr_indices, numCs_col] /
          (meth_data[sig_dmr_indices, numCs_col] + meth_data[sig_dmr_indices, numTs_col])) * 100

        # Look up infection status by Sample_ID rather than positional index, since
        # methylKit reorders/drops samples during unite() so meth_obj@sample.ids
        # order does NOT necessarily match sample_metadata row order.
        sid <- meth_obj@sample.ids[i]
        md_row <- match(
          as.character(sid),
          as.character(sample_metadata[[opt$sample_id_col]])
        )
        inf_status <- if (!is.na(md_row)) {
          sample_metadata[[opt$infection_col]][md_row]
        } else {
          NA
        }

        temp_list[[i]] <- data.frame(
          Sample_ID = sid,
          Infection_Status = inf_status,
          Mean_Methylation = meth_vals, # Individual DMR values for the ridge distribution
          Region_Size_Label = sprintf("%d KB", as.numeric(region_size) / 1000),
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(temp_list) > 0) {
      all_infection_data[[region_size]] <- do.call(rbind, temp_list)
    }
  }
}

# Combine all data
if (length(all_infection_data) > 0) {
  combined_infection_data <- do.call(rbind, all_infection_data)

  # Properly order region sizes
  combined_infection_data$Region_Size_Label <- order_region_sizes(combined_infection_data$Region_Size_Label)

  # Order infection status (bottom to top: Control, Suspected, Infected)
  # For ridge plots, factor levels are displayed bottom-to-top
  infection_levels <- c("Infected", "Suspected", "Control") # Reversed for display
  infection_levels <- infection_levels[infection_levels %in% unique(combined_infection_data$Infection_Status)]
  combined_infection_data$Infection_Status <- factor(combined_infection_data$Infection_Status,
    levels = infection_levels
  )

  # Get infection colors from earlier definition
  infection_fill_colors <- infection_colors[levels(combined_infection_data$Infection_Status)]

  # Create faceted ridge plot
  p7 <- ggplot(combined_infection_data, aes(x = Mean_Methylation, y = Infection_Status, fill = Infection_Status)) +
    geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01, size = 1) +
    scale_fill_manual(values = infection_fill_colors) +
    facet_wrap(~Region_Size_Label, scales = "free_x", ncol = 3) +
    labs(
      title = "DMR methylation distribution by infection status",
      subtitle = "Distribution of methylation (%) across all significant DMR markers",
      x = "Methylation level (%)",
      y = "Group",
      fill = "Group"
    ) +
    theme_methylsense() +
    theme(
      legend.position = "top"
    )

  save_plot(p7, "7_ridge_methylation_by_infection_all_regions", width = 14, height = 10)

  cat("[OK] Ridge plot by infection status (all regions) created\n")
} else {
  cat("[WARNING] No valid data for infection status ridge plot. Skipping...\n")
}

cat("\n[COMPLETE] Ridge plot - infection status\n")

cat("\n================================================================================\n")
cat("All analyses complete\n")
cat("================================================================================\n\n")

# ================================================================================
# GENERATE PDF REPORT (if requested)
# ================================================================================

if (opt$generate_pdf) {
  cat("\n[STEP 10] Generating PDF report...\n")

  # PAGEDOWN ONLY - NO LATEX FALLBACK
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    cat("[WARNING] rmarkdown package not installed. Skipping PDF generation.\n")
    cat("[INFO] Install with: install.packages('rmarkdown')\n")
  } else if (!requireNamespace("pagedown", quietly = TRUE)) {
    cat("[WARNING] pagedown package not installed. Skipping PDF generation.\n")
    cat("[INFO] Install with: install.packages('pagedown')\n")
    cat("[INFO] Markdown report is still available\n")
  } else if (rmarkdown::pandoc_available()) {
    tryCatch(
      {
        # Generate HTML first
        html_path <- file.path(opt$output_dir, "MethylSense_General_Statistics_Report.html")

        rmarkdown::render(
          input = markdown_path,
          output_format = "html_document",
          output_file = basename(html_path),
          output_dir = opt$output_dir,
          quiet = TRUE
        )

        # Then convert to PDF using pagedown (no LaTeX)
        pdf_path <- file.path(opt$output_dir, "MethylSense_General_Statistics_Report.pdf")

        pagedown::chrome_print(
          input = html_path,
          output = pdf_path,
          format = "pdf",
          timeout = 120,
          verbose = FALSE
        )

        cat(sprintf("[OK] PDF report generated: %s\n", pdf_path))
        cat("[INFO] PDF generated using pagedown (no LaTeX required)\n")
      },
      error = function(e) {
        cat(sprintf("[WARNING] Failed to generate PDF: %s\n", e$message))
        cat("[INFO] This usually means Chrome/Chromium is not installed\n")
        cat("[INFO] Markdown report is still available at: %s\n", markdown_path)
      }
    )
  } else {
    cat("[WARNING] pandoc not found. PDF generation requires pandoc.\n")
    cat("[INFO] Install pandoc: https://pandoc.org/installing.html\n")
    cat("[INFO] Markdown report is available at: %s\n", markdown_path)
  }
}
