#!/usr/bin/env Rscript

# ================================================================================
#                    METHYLSENSE DIAGNOSTIC PREDICTION PIPELINE v3
# ================================================================================
# Purpose: Apply trained MethylSense models to predict disease status in new samples
#          v3: Enhanced with robust chromosome normalization and type handling
#
# Author: Markus Hodal Drag
# Version: 5.7.2
# Release Date: 2026-05-01
# GitHub: https://github.com/markusdrag/MethylSense
#
# Citation:
#   Drag MH, Hvilsom C, Poulsen LL, Jensen HE, Tahas SA, Leineweber C,
#   Cray C, Bertelsen MF, Bojesen AM.
#   MethylSense: high accuracy machine learning-based diagnostics for
#   Aspergillus fumigatus infection in chickens using host cell-free DNA
#   methylation and Nanopore sequencing. J Clin Microbiol 0:e01054-25.
#   https://doi.org/10.1128/jcm.01054-25
# ================================================================================

SCRIPT_VERSION <- "5.7.2"
SCRIPT_DATE <- "2026-05-01"

suppressPackageStartupMessages({
  library(optparse)
  library(methylKit)
  library(caret)
  library(qs)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(gridExtra)
  library(pROC)
})

# ================================================================================
# COMMAND LINE INTERFACE
# ================================================================================

option_list <- list(
  make_option(c("-q", "--qs_file"),
    type = "character",
    help = "Path to methylRawList .qs file (preprocessed samples)",
    metavar = "FILE"
  ),
  make_option(c("-m", "--model_dir"),
    type = "character",
    help = "Path to trained model directory",
    metavar = "DIRECTORY"
  ),
  make_option(c("-o", "--output_dir"),
    type = "character",
    help = "Output directory for prediction results",
    metavar = "DIRECTORY"
  ),
  make_option(c("-s", "--sample_ids"),
    type = "character",
    default = NULL,
    help = "Comma-separated sample IDs (default: all)",
    metavar = "STRING"
  ),
  make_option(c("-c", "--min_coverage"),
    type = "integer",
    default = 10,
    help = "Minimum coverage per CpG [default: %default]",
    metavar = "INTEGER"
  ),
  make_option(c("-t", "--confidence_threshold"),
    type = "double",
    default = 0.7,
    help = "Confidence threshold [default: %default]",
    metavar = "DOUBLE"
  ),
  make_option(c("-p", "--plots"),
    action = "store_true",
    default = FALSE,
    help = "Create diagnostic plots"
  ),
  make_option(c("-b", "--batch_mode"),
    action = "store_true",
    default = FALSE,
    help = "Batch mode (no individual plots)"
  ),
  make_option(c("-v", "--verbose"),
    action = "store_true",
    default = FALSE,
    help = "Verbose output"
  )
)

parser <- OptionParser(
  usage = "%prog [options]",
  option_list = option_list,
  description = "\nMethylSense Diagnostic Prediction Pipeline v3",
  epilogue = NULL
)

# Parse arguments - handle both command line and RStudio
opt <- tryCatch(
  {
    parse_args(parser, args = commandArgs(trailingOnly = TRUE))
  },
  error = function(e) {
    if (grepl("redundant long names", e$message)) {
      # Fallback: parse with positional = FALSE
      parse_args(parser, args = commandArgs(trailingOnly = TRUE), positional_arguments = FALSE)
    } else {
      stop(e)
    }
  }
)

# Validate required arguments
if (is.null(opt$qs_file) || is.null(opt$model_dir) || is.null(opt$output_dir)) {
  cat("\n[ERROR] Missing required arguments\n\n")
  cat("Required:\n")
  cat("  -q/--qs_file      : Path to .qs file\n")
  cat("  -m/--model_dir    : Path to trained model directory\n")
  cat("  -o/--output_dir   : Output directory\n\n")
  cat("Example:\n")
  cat("  Rscript script.R -q samples.qs -m model/ -o predictions/\n\n")
  stop("Missing required arguments", call. = FALSE)
}

# Handle plots flag
create_plots <- !is.null(opt$plots) && opt$plots

# Create output directory
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}

today <- format(Sys.Date(), "%Y%m%d")

# ================================================================================
# HEADER DISPLAY
# ================================================================================

cat("\n")
cat("================================================================================\n")
cat("                                                                                \n")
cat("        M   M  EEEEE  TTTTT  H   H  Y   Y  L      SSSS  EEEEE  N   N  SSSS  EEEEE\n")
cat("        MM MM  E        T    H   H   Y Y   L      S     E      NN  N  S     E    \n")
cat("        M M M  EEEE     T    HHHHH    Y    L       SSS  EEEE   N N N   SSS  EEEE \n")
cat("        M   M  E        T    H   H    Y    L          S E      N  NN      S E    \n")
cat("        M   M  EEEEE    T    H   H    Y    LLLLL  SSSS  EEEEE  N   N  SSSS  EEEEE\n")
cat("                                                                                \n")
cat("================================================================================\n")
cat("                     Diagnostic prediction pipeline v3                        \n")
cat("================================================================================\n")
cat(paste("Version:", SCRIPT_VERSION, "|", "Release Date:", SCRIPT_DATE, "\n"))
cat(paste("Start Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
cat(paste("Working Directory:", getwd(), "\n"))
cat(paste("R Version:", R.version.string, "\n"))
cat("\n")
cat("Please cite:\n")
cat("  Drag MH, Hvilsom C, Poulsen LL, Jensen HE, Tahas SA, Leineweber C,\n")
cat("  Cray C, Bertelsen MF, Bojesen AM.\n")
cat("  MethylSense: high accuracy machine learning-based diagnostics for\n")
cat("  Aspergillus fumigatus infection in chickens using host cell-free DNA\n")
cat("  methylation and Nanopore sequencing. J Clin Microbiol 0:e01054-25.\n")
cat("  https://doi.org/10.1128/jcm.01054-25\n")
cat("================================================================================\n\n")

# ================================================================================
# CONFIGURATION SUMMARY
# ================================================================================

cat("Configuration:\n")
cat("  Input File:        ", basename(opt$qs_file), "\n")
cat("  Model Directory:   ", opt$model_dir, "\n")
cat("  Output Directory:  ", opt$output_dir, "\n")
cat("  Min Coverage:      ", opt$min_coverage, "reads per CpG\n")
cat("  Confidence Thresh: ", opt$confidence_threshold, "\n")
cat("  Create Plots:      ", create_plots, "\n")
cat("  Batch Mode:        ", opt$batch_mode, "\n")
cat("\n")

# ================================================================================
# HELPER FUNCTIONS
# ================================================================================

log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  cat(timestamp, paste0("[", level, "]"), msg, "\n")
}

verbose_msg <- function(msg) {
  if (opt$verbose) {
    log_msg(paste("   [VERBOSE]", msg), level = "VERBOSE")
  }
}

# Robust chromosome normalization function
# Handles: numeric vs character, chr prefix, whitespace, case
normalize_chr <- function(chr) {
  # Convert to character if numeric
  chr <- as.character(chr)

  # Remove "chr" prefix if present (case-insensitive)
  chr <- sub("^chr", "", chr, ignore.case = TRUE)

  # Trim whitespace
  chr <- trimws(chr)

  # Handle special chromosomes - keep original case for numbers, uppercase for letters
  # This handles X, Y, M, MT properly
  # Use vectorized ifelse to handle vectors
  result <- ifelse(grepl("^[0-9]+$", chr),
    chr, # Numeric chromosome - keep as is
    toupper(chr)
  ) # Non-numeric (X, Y, M, MT) - uppercase

  return(result)
}

# ================================================================================
# STEP 1: LOAD AND VALIDATE MODEL
# ================================================================================

cat("================================================================================\n")
cat("Step 1: Loading trained model\n")
cat("================================================================================\n\n")

if (!dir.exists(opt$model_dir)) {
  stop(paste("[ERROR] Model directory not found:", opt$model_dir))
}

model_file <- file.path(opt$model_dir, "model.rds")
if (!file.exists(model_file)) {
  stop(paste("[ERROR] Model file not found:", model_file))
}

log_msg("Loading trained model...")
trained_model <- tryCatch(
  {
    readRDS(model_file)
  },
  error = function(e) {
    stop(paste("[ERROR] Failed to load model:", e$message))
  }
)

log_msg(paste("[OK] Loaded model:", trained_model$method))
log_msg(paste("     Model type:", trained_model$modelType))

dmr_file <- file.path(opt$model_dir, "dmr_details_with_coordinates.csv")
if (!file.exists(dmr_file)) {
  stop(paste("[ERROR] DMR details file not found:", dmr_file))
}

log_msg("Loading DMR coordinates...")
dmr_coords <- tryCatch(
  {
    read.csv(dmr_file, stringsAsFactors = FALSE)
  },
  error = function(e) {
    stop(paste("[ERROR] Failed to load DMR coordinates:", e$message))
  }
)

# Store original chromosome names for debugging
dmr_coords$scaffold_original <- dmr_coords$scaffold

# Normalize chromosome names in DMR coordinates
log_msg("Normalising DMR chromosome names...")
dmr_coords$scaffold <- normalize_chr(dmr_coords$scaffold)

n_markers <- nrow(dmr_coords)
log_msg(paste("[OK] Loaded", n_markers, "DMR markers"))
log_msg(paste("     Original format:", paste(head(unique(dmr_coords$scaffold_original), 5), collapse = ", ")))
log_msg(paste("     Normalised format:", paste(head(unique(dmr_coords$scaffold), 5), collapse = ", ")))

cat("\nDMR Summary:\n")
cat("  Total markers:     ", n_markers, "\n")
cat("  Hypermethylated:   ", sum(dmr_coords$direction == "Hypermethylated", na.rm = TRUE), "\n")
cat("  Hypomethylated:    ", sum(dmr_coords$direction == "Hypomethylated", na.rm = TRUE), "\n")
cat("  Mean width:        ", round(mean(dmr_coords$width_bp, na.rm = TRUE)), "bp\n")
cat("  Mean q-value:      ", signif(mean(dmr_coords$qvalue, na.rm = TRUE), 3), "\n")
cat("\n")

model_name <- trained_model$method
model_accuracy <- dmr_coords$accuracy[1]
model_auc <- dmr_coords$auc[1]

log_msg(paste("Model Performance (Training):"))
log_msg(paste("  Accuracy:", round(model_accuracy, 3)))
log_msg(paste("  AUC:     ", round(model_auc, 3)))

class_levels <- trained_model$levels
log_msg(paste("  Classes: ", paste(class_levels, collapse = " vs ")))

# ================================================================================
# STEP 2: LOAD NEW SAMPLES
# ================================================================================

cat("\n")
cat("================================================================================\n")
cat("Step 2: Loading new samples\n")
cat("================================================================================\n\n")

log_msg(paste("Loading samples from:", basename(opt$qs_file)))
start_load <- Sys.time()

new_samples <- tryCatch(
  {
    qread(opt$qs_file)
  },
  error = function(e) {
    stop(paste("[ERROR] Failed to load samples:", e$message))
  }
)

end_load <- Sys.time()
log_msg(paste(
  "[OK] Loaded", length(new_samples), "samples in",
  round(as.numeric(end_load - start_load, units = "secs")), "seconds"
))

verbose_msg("Normalising coverage using median method...")
new_samples <- normalizeCoverage(new_samples, method = "median")
log_msg("[OK] Coverage normalised")

# Normalize chromosome names in all samples
log_msg("Normalising chromosome names in samples...")
for (i in seq_along(new_samples)) {
  sample_data <- getData(new_samples[[i]])

  # Store original for debugging
  sample_data$chr_original <- sample_data$chr

  # Normalize chromosome names
  sample_data$chr <- normalize_chr(sample_data$chr)

  # Update the methylRaw object with normalized data
  new_samples[[i]]@.Data <- sample_data

  if (i == 1) {
    log_msg(paste("     Sample chr format (original):", paste(head(unique(sample_data$chr_original), 5), collapse = ", ")))
    log_msg(paste("     Sample chr format (normalised):", paste(head(unique(sample_data$chr), 5), collapse = ", ")))
  }
}
log_msg("[OK] Chromosome names normalised")

sample_ids <- sapply(new_samples, function(x) x@sample.id)
log_msg(paste(
  "Sample IDs:", paste(head(sample_ids, 3), collapse = ", "),
  if (length(sample_ids) > 3) paste("... and", length(sample_ids) - 3, "more") else ""
))

if (!is.null(opt$sample_ids)) {
  requested_ids <- trimws(strsplit(opt$sample_ids, ",")[[1]])
  log_msg(paste("Filtering to requested samples:", paste(requested_ids, collapse = ", ")))

  match_indices <- which(sample_ids %in% requested_ids)

  if (length(match_indices) == 0) {
    stop(paste(
      "[ERROR] None of the requested sample IDs found in data.\n",
      "Available IDs:", paste(sample_ids, collapse = ", ")
    ))
  }

  if (length(match_indices) < length(requested_ids)) {
    missing_ids <- setdiff(requested_ids, sample_ids[match_indices])
    log_msg(paste("[WARN] Could not find samples:", paste(missing_ids, collapse = ", ")))
  }

  new_samples <- new_samples[match_indices]
  sample_ids <- sample_ids[match_indices]

  log_msg(paste("[OK] Filtered to", length(new_samples), "samples"))
}

sample_sizes <- sapply(new_samples, nrow)
log_msg(paste(
  "Sample sizes: min =", format(min(sample_sizes), big.mark = ","),
  "| max =", format(max(sample_sizes), big.mark = ","),
  "| mean =", format(round(mean(sample_sizes)), big.mark = ",")
))

# Enhanced diagnostic output
cat("\n=== Diagnostic: data compatibility check ===\n")
test_sample <- getData(new_samples[[1]])

# Chromosome overlap check
sample_chrs <- unique(test_sample$chr)
dmr_chrs <- unique(dmr_coords$scaffold)
overlap_chrs <- intersect(sample_chrs, dmr_chrs)

cat("Chromosome Analysis:\n")
cat("  Sample chromosomes:", paste(head(sample_chrs, 10), collapse = ", "), "\n")
cat("  DMR chromosomes:   ", paste(head(dmr_chrs, 10), collapse = ", "), "\n")
cat("  Overlapping:       ", length(overlap_chrs), "of", length(dmr_chrs), "DMR chromosomes\n")

if (length(overlap_chrs) > 0) {
  cat("  Overlap list:      ", paste(head(overlap_chrs, 15), collapse = ", "), "\n")
} else {
  cat("  [WARNING] NO chromosome overlap detected!\n")
}

# Coordinate range check
cat("\nCoordinate Ranges:\n")
cat(
  "  DMR regions:       ", format(min(dmr_coords$start), big.mark = ","), "-",
  format(max(dmr_coords$end), big.mark = ","), "\n"
)
cat(
  "  Sample CpGs:       ", format(min(test_sample$start), big.mark = ","), "-",
  format(max(test_sample$end), big.mark = ","), "\n"
)

# Coverage check
cat("\nCoverage Analysis:\n")
cat("  Sample coverage:   ", min(test_sample$coverage), "-", max(test_sample$coverage), "\n")
cat("  Min threshold:     ", opt$min_coverage, "\n")
cat(
  "  CpGs >= threshold: ", sum(test_sample$coverage >= opt$min_coverage),
  "(", round(100 * sum(test_sample$coverage >= opt$min_coverage) / nrow(test_sample), 1), "%)\n"
)

cat("=========================================\n\n")

# ================================================================================
# STEP 3: EXTRACT DMR METHYLATION
# ================================================================================

cat("================================================================================\n")
cat("Step 3: Extracting DMR methylation\n")
cat("================================================================================\n\n")

log_msg("Creating genomic regions for DMRs...")

dmr_gr <- GRanges(
  seqnames = dmr_coords$scaffold,
  ranges = IRanges(start = dmr_coords$start, end = dmr_coords$end)
)

log_msg(paste("[OK] Created", length(dmr_gr), "genomic regions"))

log_msg("Extracting methylation values...")
extraction_start <- Sys.time()

methylation_matrix <- matrix(NA, nrow = length(new_samples), ncol = length(dmr_gr))
rownames(methylation_matrix) <- sample_ids
colnames(methylation_matrix) <- dmr_coords$DMR_ID

# Track extraction statistics
extraction_stats <- data.frame(
  sample_id = character(),
  cpgs_extracted = integer(),
  dmrs_with_data = integer(),
  stringsAsFactors = FALSE
)

for (i in seq_along(new_samples)) {
  sample_data <- new_samples[[i]]
  sample_df <- getData(sample_data)

  cpg_count <- 0
  dmr_count <- 0

  # Debug first sample
  if (i == 1) {
    verbose_msg("Processing first sample for detailed diagnostics...")
    verbose_msg(paste("  Total CpGs in sample:", nrow(sample_df)))
    verbose_msg(paste("  Chr data type:", class(sample_df$chr)))
    verbose_msg(paste("  First 3 chromosomes:", paste(head(sample_df$chr, 3), collapse = ", ")))
  }

  for (j in seq_along(dmr_gr)) {
    region <- dmr_gr[j]
    region_chr <- as.character(seqnames(region))
    region_start <- start(region)
    region_end <- end(region)

    # Normalize region chromosome (should already be normalized, but double-check)
    region_chr_norm <- normalize_chr(region_chr)

    # Extract CpGs in this region with sufficient coverage
    region_cpgs <- sample_df[
      sample_df$chr == region_chr_norm &
        sample_df$start >= region_start &
        sample_df$end <= region_end &
        sample_df$coverage >= opt$min_coverage,
    ]

    if (nrow(region_cpgs) > 0) {
      cpg_count <- cpg_count + nrow(region_cpgs)
      dmr_count <- dmr_count + 1

      total_cs <- sum(region_cpgs$numCs, na.rm = TRUE)
      total_coverage <- sum(region_cpgs$coverage, na.rm = TRUE)

      if (total_coverage > 0) {
        methylation_matrix[i, j] <- (total_cs / total_coverage) * 100
      }

      # Extra verbose output for first sample and first few DMRs
      if (i == 1 && j <= 3) {
        verbose_msg(paste(
          "  DMR", j, "(", region_chr_norm, ":", region_start, "-", region_end, "):",
          nrow(region_cpgs), "CpGs,",
          round((total_cs / total_coverage) * 100, 2), "% methylated"
        ))
      }
    }
  }

  # Store stats
  extraction_stats <- rbind(extraction_stats, data.frame(
    sample_id = sample_ids[i],
    cpgs_extracted = cpg_count,
    dmrs_with_data = dmr_count,
    stringsAsFactors = FALSE
  ))

  if (i %% 10 == 0 || i == length(new_samples)) {
    non_na <- sum(!is.na(methylation_matrix[i, ]))
    verbose_msg(paste(
      "Processed", i, "of", length(new_samples),
      "samples | DMRs with data:", non_na, "/", ncol(methylation_matrix),
      "| CpGs:", cpg_count
    ))
  }
}

extraction_end <- Sys.time()
log_msg(paste(
  "[OK] Extraction complete in",
  round(as.numeric(extraction_end - extraction_start, units = "secs")), "seconds"
))

# Enhanced diagnostic information
missing_count <- sum(is.na(methylation_matrix))
missing_pct <- round(100 * missing_count / length(methylation_matrix), 2)

cat("\nExtraction Diagnostics:\n")
cat("  Total values expected: ", length(methylation_matrix), "\n")
cat("  Missing values:        ", missing_count, " (", missing_pct, "%)\n", sep = "")
cat("  Successfully extracted:", length(methylation_matrix) - missing_count, "\n")
cat("  Mean CpGs per sample:  ", round(mean(extraction_stats$cpgs_extracted)), "\n")
cat("  Mean DMRs per sample:  ", round(mean(extraction_stats$dmrs_with_data), 1), "of", ncol(methylation_matrix), "\n")

# Check if ALL values are missing (critical error)
if (missing_count == length(methylation_matrix)) {
  cat("\n[ERROR] NO methylation values were extracted!\n")
  cat("\nDetailed Debugging Information:\n")

  sample_df <- getData(new_samples[[1]])

  cat("\nChromosome Format:\n")
  cat("  DMR (original):    ", paste(head(unique(dmr_coords$scaffold_original), 5), collapse = ", "), "\n")
  cat("  DMR (normalised):  ", paste(head(unique(dmr_coords$scaffold), 5), collapse = ", "), "\n")
  cat("  Sample (original): ", paste(head(unique(sample_df$chr_original), 5), collapse = ", "), "\n")
  cat("  Sample (normalised):", paste(head(unique(sample_df$chr), 5), collapse = ", "), "\n")

  # Check chromosome overlap
  dmr_chrs <- unique(dmr_coords$scaffold)
  sample_chrs <- unique(sample_df$chr)
  overlap <- intersect(dmr_chrs, sample_chrs)

  cat("\nChromosome Overlap:\n")
  cat("  DMR chromosomes:   ", length(dmr_chrs), "\n")
  cat("  Sample chromosomes:", length(sample_chrs), "\n")
  cat("  Overlapping:       ", length(overlap), "\n")

  if (length(overlap) == 0) {
    cat("\n  [CRITICAL] No chromosome overlap after normalisation!\n")
    cat("  This indicates incompatible genome assemblies.\n")
  } else {
    cat("  Overlapping chrs:  ", paste(overlap, collapse = ", "), "\n")

    # Check coordinate ranges on overlapping chromosomes
    cat("\nCoordinate Range Analysis:\n")
    cat(
      "  DMR start-end:     ", format(min(dmr_coords$start), big.mark = ","), "-",
      format(max(dmr_coords$end), big.mark = ","), "\n"
    )
    cat(
      "  Sample start-end:  ", format(min(sample_df$start), big.mark = ","), "-",
      format(max(sample_df$end), big.mark = ","), "\n"
    )

    # Check if ranges overlap at all
    for (chr in head(overlap, 3)) {
      dmr_on_chr <- dmr_coords[dmr_coords$scaffold == chr, ]
      sample_on_chr <- sample_df[sample_df$chr == chr, ]

      if (nrow(dmr_on_chr) > 0 && nrow(sample_on_chr) > 0) {
        cat(paste("\n  Chr", chr, ":\n"))
        cat(paste(
          "    DMR range:   ", format(min(dmr_on_chr$start), big.mark = ","), "-",
          format(max(dmr_on_chr$end), big.mark = ","), "\n"
        ))
        cat(paste(
          "    Sample range:", format(min(sample_on_chr$start), big.mark = ","), "-",
          format(max(sample_on_chr$end), big.mark = ","), "\n"
        ))

        if (max(sample_on_chr$end) < min(dmr_on_chr$start)) {
          cat("    [CRITICAL] No coordinate overlap! Sample doesn't reach DMR regions.\n")
        }
      }
    }
  }

  cat("\nCoverage Analysis:\n")
  cat("  Sample coverage range:", min(sample_df$coverage), "-", max(sample_df$coverage), "\n")
  cat("  Required threshold:   ", opt$min_coverage, "\n")
  cat(
    "  CpGs meeting threshold:", sum(sample_df$coverage >= opt$min_coverage),
    "(", round(100 * sum(sample_df$coverage >= opt$min_coverage) / nrow(sample_df), 1), "%)\n"
  )

  stop("\n[FATAL] Failed to extract methylation data. See detailed diagnostics above.")
}

# Only impute if there are some successful extractions
if (missing_count > 0 && missing_count < length(methylation_matrix)) {
  log_msg(paste("[WARN] Missing values:", missing_count, "(", missing_pct, "%)"))
  log_msg("        Imputing with column medians...")

  for (col in 1:ncol(methylation_matrix)) {
    col_data <- methylation_matrix[, col]
    if (any(is.na(col_data))) {
      median_val <- median(col_data, na.rm = TRUE)
      # Only use 50 as fallback if no valid data at all
      if (is.na(median_val)) {
        median_val <- 50
        verbose_msg(paste("DMR", col, "has no data, using default 50%"))
      }
      methylation_matrix[is.na(methylation_matrix[, col]), col] <- median_val
    }
  }

  log_msg("[OK] Missing value imputation complete")
}

methylation_df <- as.data.frame(methylation_matrix)

# Final sanity check
unique_values <- length(unique(as.vector(methylation_matrix)))
cat("  Unique methylation values: ", unique_values, "\n")
if (unique_values == 1) {
  cat("\n[WARNING] All methylation values are identical!\n")
  cat("This suggests a systematic extraction problem.\n\n")
}

log_msg(paste(
  "Feature matrix:", nrow(methylation_df), "samples x",
  ncol(methylation_df), "DMR features"
))

# ================================================================================
# STEP 4: MAKE PREDICTIONS
# ================================================================================

cat("\n")
cat("================================================================================\n")
cat("Step 4: Making predictions\n")
cat("================================================================================\n\n")

log_msg("Applying trained model to new samples...")
prediction_start <- Sys.time()

predicted_classes <- predict(trained_model, methylation_df)
predicted_probs <- predict(trained_model, methylation_df, type = "prob")

prediction_end <- Sys.time()
log_msg(paste(
  "[OK] Predictions complete in",
  round(as.numeric(prediction_end - prediction_start, units = "secs")), "seconds"
))

confidence_scores <- apply(predicted_probs, 1, max)

low_confidence_samples <- which(confidence_scores < opt$confidence_threshold)

if (length(low_confidence_samples) > 0) {
  log_msg(paste(
    "[WARN]", length(low_confidence_samples), "samples have confidence below",
    opt$confidence_threshold
  ))
}

results_df <- data.frame(
  Sample_ID = sample_ids,
  Predicted_Class = as.character(predicted_classes),
  Confidence = round(confidence_scores, 3),
  Confidence_Flag = ifelse(confidence_scores >= opt$confidence_threshold, "High", "Low"),
  stringsAsFactors = FALSE
)

for (i in 1:ncol(predicted_probs)) {
  class_name <- colnames(predicted_probs)[i]
  results_df[[paste0("Prob_", class_name)]] <- round(predicted_probs[, i], 3)
}

results_df <- results_df[order(-results_df$Confidence), ]

cat("\nPrediction Summary:\n")
prediction_table <- table(results_df$Predicted_Class)
for (i in 1:length(prediction_table)) {
  cat("  ", names(prediction_table)[i], ": ", prediction_table[i], " samples\n", sep = "")
}
cat("\n")

cat("Confidence Summary:\n")
cat("  Mean confidence:   ", round(mean(confidence_scores), 3), "\n")
cat("  Median confidence: ", round(median(confidence_scores), 3), "\n")
cat("  Min confidence:    ", round(min(confidence_scores), 3), "\n")
cat("  Max confidence:    ", round(max(confidence_scores), 3), "\n")
cat("  High confidence:   ", sum(results_df$Confidence_Flag == "High"), " samples\n")
cat("  Low confidence:    ", sum(results_df$Confidence_Flag == "Low"), " samples\n")
cat("\n")

# Display individual predictions
cat("Individual Predictions:\n")
cat("--------------------------------------------------------------------------------\n")
cat(sprintf("  %-25s %-15s %-12s %s\n", "Sample", "Prediction", "Confidence", "Flag"))
cat("--------------------------------------------------------------------------------\n")
for (i in 1:nrow(results_df)) {
  cat(sprintf(
    "  %-25s %-15s %-12.3f %s\n",
    results_df$Sample_ID[i],
    results_df$Predicted_Class[i],
    results_df$Confidence[i],
    results_df$Confidence_Flag[i]
  ))
}
cat("--------------------------------------------------------------------------------\n")
cat("\n")

# ================================================================================
# STEP 5: SAVE RESULTS
# ================================================================================

cat("================================================================================\n")
cat("Step 5: Saving results\n")
cat("================================================================================\n\n")

results_file <- file.path(opt$output_dir, paste0(today, "_predictions_detailed.csv"))
write.csv(results_df, results_file, row.names = FALSE)
log_msg(paste("[OK] Detailed results:", basename(results_file)))

clinical_report <- results_df[, c("Sample_ID", "Predicted_Class", "Confidence", "Confidence_Flag")]
clinical_file <- file.path(opt$output_dir, paste0(today, "_predictions_clinical_report.csv"))
write.csv(clinical_report, clinical_file, row.names = FALSE)
log_msg(paste("[OK] Clinical report:", basename(clinical_file)))

features_file <- file.path(opt$output_dir, paste0(today, "_methylation_features.csv"))
methylation_export <- cbind(Sample_ID = sample_ids, methylation_df)
write.csv(methylation_export, features_file, row.names = FALSE)
log_msg(paste("[OK] Feature matrix: ", basename(features_file)))

# Save extraction statistics
extraction_stats_file <- file.path(opt$output_dir, paste0(today, "_extraction_stats.csv"))
write.csv(extraction_stats, extraction_stats_file, row.names = FALSE)
log_msg(paste("[OK] Extraction stats:", basename(extraction_stats_file)))

summary_file <- file.path(opt$output_dir, paste0(today, "_prediction_summary.txt"))
writeLines(c(
  "================================================================================",
  "MethylSense diagnostic prediction summary v3",
  "================================================================================",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("Version:", SCRIPT_VERSION),
  "",
  "INPUT:",
  paste("  Samples file:", basename(opt$qs_file)),
  paste("  Model directory:", opt$model_dir),
  paste("  Total samples:", length(sample_ids)),
  "",
  "MODEL INFORMATION:",
  paste("  Model type:", model_name),
  paste("  Number of markers:", n_markers),
  paste("  Training accuracy:", round(model_accuracy, 3)),
  paste("  Training AUC:", round(model_auc, 3)),
  paste("  Classes:", paste(class_levels, collapse = " vs ")),
  "",
  "DATA EXTRACTION:",
  paste("  Mean CpGs extracted per sample:", round(mean(extraction_stats$cpgs_extracted))),
  paste("  Mean DMRs with data per sample:", round(mean(extraction_stats$dmrs_with_data), 1)),
  paste("  Missing values:", missing_count, "(", missing_pct, "%)"),
  "",
  "PREDICTION RESULTS:",
  paste("  Samples processed:", nrow(results_df)),
  sapply(1:length(prediction_table), function(i) {
    paste("  ", names(prediction_table)[i], ":", prediction_table[i], "samples")
  }),
  "",
  "CONFIDENCE METRICS:",
  paste("  Mean confidence:", round(mean(confidence_scores), 3)),
  paste("  Median confidence:", round(median(confidence_scores), 3)),
  paste("  Range:", round(min(confidence_scores), 3), "-", round(max(confidence_scores), 3)),
  paste(
    "  High confidence (>=", opt$confidence_threshold, "):",
    sum(results_df$Confidence_Flag == "High"), "samples"
  ),
  paste(
    "  Low confidence (<", opt$confidence_threshold, "):",
    sum(results_df$Confidence_Flag == "Low"), "samples"
  ),
  "",
  "OUTPUT FILES:",
  paste("  Detailed predictions:", basename(results_file)),
  paste("  Clinical report:", basename(clinical_file)),
  paste("  Feature matrix:", basename(features_file)),
  paste("  Extraction stats:", basename(extraction_stats_file)),
  "================================================================================"
), summary_file)

log_msg(paste("[OK] Summary report: ", basename(summary_file)))

# ================================================================================
# COMPLETION
# ================================================================================

cat("\n")
cat("================================================================================\n")
cat("Prediction pipeline complete\n")
cat("================================================================================\n\n")

cat("Summary:\n")
cat("  Samples processed: ", length(sample_ids), "\n")
cat("  Output directory:  ", opt$output_dir, "\n")
cat("\n")

cat("Key Output Files:\n")
cat("  1. Detailed predictions:       ", basename(results_file), "\n")
cat("  2. Clinical report:            ", basename(clinical_file), "\n")
cat("  3. Feature matrix:             ", basename(features_file), "\n")
cat("  4. Extraction statistics:      ", basename(extraction_stats_file), "\n")
cat("  5. Summary report:             ", basename(summary_file), "\n")
cat("\n")

if (length(low_confidence_samples) > 0) {
  cat("================================================================================\n")
  cat("Attention: low confidence predictions\n")
  cat("================================================================================\n\n")

  # Get the low confidence rows from the sorted results_df
  low_conf_rows <- results_df[results_df$Confidence_Flag == "Low", ]

  cat(paste(
    nrow(low_conf_rows),
    "sample(s) have prediction confidence below", opt$confidence_threshold, ":\n\n"
  ))

  for (i in 1:nrow(low_conf_rows)) {
    cat(sprintf(
      "  %s: %s (confidence = %.3f)\n",
      low_conf_rows$Sample_ID[i],
      low_conf_rows$Predicted_Class[i],
      low_conf_rows$Confidence[i]
    ))
  }

  cat("\nRecommendation: Review these samples carefully.\n")
  cat("\n")
}

cat("================================================================================\n")
cat("For questions or support, refer to the MethylSense documentation.\n")
cat("================================================================================\n\n")

log_msg("[SUCCESS] MethylSense prediction pipeline v3 completed successfully!")
log_msg(paste("[DIR] All results available in:", opt$output_dir))
