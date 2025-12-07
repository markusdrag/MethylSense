#!/usr/bin/env Rscript

# ================================================================================
#                    METHYLSENSE DATA PREPROCESSING PIPELINE
# ================================================================================
# Purpose: Convert raw bedFileOrg (13/18 col) → standardized bedFile (8 col)
#          → methylKit methylRaw objects for downstream DMR and ML analysis
#
# Author: Markus Hodal Drag
# Version: 5.13.7
# Release Date: 2025-12-03
# GitHub: https://github.com/markusdrag/MethylSense
#
# Citation:
#   Drag MH, Hvilsom C, Poulsen LL, Jensen HE, Tahas SA, Leineweber C,
#   Cray C, Bertelsen MF, Bojesen AM (2025)
#   New high accuracy diagnostics for avian Aspergillus fumigatus infection
#   using Nanopore methylation sequencing of host cell-free DNA and machine
#   learning prediction. bioRxiv 2025.04.11.648151
#   https://doi.org/10.1101/2025.04.11.648151
# ================================================================================

SCRIPT_VERSION <- "5.13.7"
SCRIPT_DATE <- "2025-12-03"

suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
  library(methylKit)
  library(data.table)
  library(qs)
  library(dplyr)
  library(parallel)
})

# ================================================================================
# COMMAND LINE INTERFACE
# ================================================================================

option_list <- list(
  make_option(c("-s", "--species"),
    type = "character",
    help = "Species name (must match Sample Sheet)",
    metavar = "CHARACTER"
  ),
  make_option(c("-i", "--sample_sheet"),
    type = "character",
    help = "Path to sample sheet (Excel format)",
    metavar = "FILE"
  ),
  make_option(c("-b", "--bed_dir"),
    type = "character",
    help = "Directory containing BED files",
    metavar = "DIRECTORY"
  ),
  make_option(c("-o", "--output_dir"),
    type = "character",
    help = "Output directory for processed objects",
    metavar = "DIRECTORY"
  ),
  make_option(c("-a", "--assembly"),
    type = "character",
    default = "",
    help = "Genome assembly name (default: auto-generated from species)",
    metavar = "CHARACTER"
  ),
  make_option(c("-c", "--cores"),
    type = "integer",
    default = 20,
    help = "Number of CPU cores [default: %default]",
    metavar = "INTEGER"
  ),
  make_option(c("-m", "--min_coverage"),
    type = "integer",
    default = 1,
    help = "Minimum coverage per CpG [default: %default]",
    metavar = "INTEGER"
  ),
  make_option("--force_convert",
    action = "store_true",
    default = FALSE,
    help = "Force re-conversion of existing bedFile files"
  ),
  make_option("--no_parallel",
    action = "store_true",
    default = FALSE,
    help = "Disable parallel processing"
  )
)

parser <- OptionParser(
  usage = "%prog [options]",
  option_list = option_list,
  description = "\nMethylSense Data Preprocessing Pipeline\nPrepares methylation data for DMR analysis and ML classification",
  epilogue = "Example:\n  Rscript methylsense_preprocess.R --species 'Minipig' --sample_sheet samples.xlsx --bed_dir ./beds --output_dir ./processed"
)

opt <- parse_args(parser)

# Validate required arguments
if (is.null(opt$species) || is.null(opt$sample_sheet) ||
  is.null(opt$bed_dir) || is.null(opt$output_dir)) {
  print_help(parser)
  stop("\nError: Missing required arguments.\n", call. = FALSE)
}

# Create output directory
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}

today <- format(Sys.Date(), "%Y%m%d")

# ================================================================================
# CREATE ASSEMBLY NAME - CRITICAL: DO THIS BEFORE ANY OUTPUT
# ================================================================================

if (is.null(opt$assembly) || opt$assembly == "") {
  # Auto-generate from species name
  species_clean <- gsub(" ", "_", opt$species)
  species_clean <- gsub("[^A-Za-z0-9_]", "", species_clean)
  assembly_name <- paste0(species_clean, "_assembly")
} else {
  # Use user-provided assembly name
  assembly_name <- opt$assembly
}

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
cat("                      DATA PREPROCESSING PIPELINE                              \n")
cat("================================================================================\n")
cat(paste("Version:", SCRIPT_VERSION, "|", "Release Date:", SCRIPT_DATE, "\n"))
cat(paste("Start Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
cat(paste("Working Directory:", getwd(), "\n"))
cat(paste("R Version:", R.version.string, "\n"))
cat("\n")
cat("Please cite:\n")
cat("  Drag MH, Hvilsom C, Poulsen LL, Jensen HE, Tahas SA, Leineweber C,\n")
cat("  Cray C, Bertelsen MF, Bojesen AM (2025)\n")
cat("  New high accuracy diagnostics for avian Aspergillus fumigatus infection\n")
cat("  using Nanopore methylation sequencing of host cell-free DNA and machine\n")
cat("  learning prediction. bioRxiv 2025.04.11.648151\n")
cat("  https://doi.org/10.1101/2025.04.11.648151\n")
cat("================================================================================\n\n")

# ================================================================================
# CONFIGURATION SUMMARY
# ================================================================================

cat("Configuration:\n")
cat("  Species:            ", opt$species, "\n")
cat("  Assembly:           ", assembly_name, "\n")
cat("  Sample Sheet:       ", basename(opt$sample_sheet), "\n")
cat("  BED Directory:      ", opt$bed_dir, "\n")
cat("  Output Directory:   ", opt$output_dir, "\n")
cat("  CPU Cores:          ", opt$cores, "\n")
cat("  Min Coverage:       ", opt$min_coverage, " reads per CpG\n")
cat("  Parallel Mode:      ", ifelse(opt$no_parallel, "Disabled", "Enabled"), "\n")
cat("  Force Convert:      ", ifelse(opt$force_convert, "Yes", "No"), "\n")
cat("\n")

# ================================================================================
# STEP 1: LOAD AND VALIDATE SAMPLE SHEET
# ================================================================================

cat("================================================================================\n")
cat("STEP 1: Loading Sample Sheet\n")
cat("================================================================================\n\n")

sheet_all <- tryCatch(
  {
    read_xlsx(opt$sample_sheet)
  },
  error = function(e) {
    stop("Failed to read sample sheet: ", e$message, call. = FALSE)
  }
)

sheet <- sheet_all %>% filter(Species == opt$species)

if (nrow(sheet) == 0) {
  stop("No samples found for species: ", opt$species, "\n",
    "Available species: ", paste(unique(sheet_all$Species), collapse = ", "),
    call. = FALSE
  )
}

cat("Sample Sheet Summary:\n")
cat("  Total samples in sheet: ", nrow(sheet_all), "\n")
cat("  Samples for", opt$species, ":", nrow(sheet), "\n")
cat("  Columns: bedFileOrg (input) -> bedFile (output)\n")
cat("\n")

original_bed_files <- file.path(opt$bed_dir, sheet$bedFileOrg)
converted_files <- file.path(opt$bed_dir, sheet$bedFile)
treatment_vector <- sheet$treatMethylkit

# Validate file existence
missing_originals <- original_bed_files[!file.exists(original_bed_files)]
existing_originals <- original_bed_files[file.exists(original_bed_files)]

if (length(missing_originals) > 0) {
  cat("Warning: Missing", length(missing_originals), "original bedFileOrg files:\n")
  for (missing in head(missing_originals, 3)) {
    cat("  -", basename(missing), "\n")
  }
  if (length(missing_originals) > 3) {
    cat("  - ... and", length(missing_originals) - 3, "more\n")
  }
  cat("\n")
}

cat("File Status:\n")
cat("  Found:   ", length(existing_originals), "bedFileOrg files\n")
cat("  Missing: ", length(missing_originals), "bedFileOrg files\n")
cat("\n")

# ================================================================================
# STEP 2: BED FILE CONVERSION
# ================================================================================

cat("================================================================================\n")
cat("STEP 2: Converting BED Files to Standard Format\n")
cat("================================================================================\n\n")

files_to_convert <- which(file.exists(original_bed_files) &
  (!file.exists(converted_files) | opt$force_convert))

if (length(files_to_convert) == 0) {
  cat("All bedFile files already exist. Skipping conversion.\n")
  cat("Use --force_convert to regenerate all files.\n\n")
} else {
  cat("Converting", length(files_to_convert), "files...\n")
  cat("Format: bedFileOrg (13 or 18 columns) -> bedFile (8 columns)\n\n")

  convert_function <- function(i) {
    infile <- original_bed_files[i]
    outfile <- converted_files[i]
    sample_id <- sheet$ID[i]

    tryCatch(
      {
        bed <- fread(infile, sep = "\t", header = FALSE, showProgress = FALSE)

        if (ncol(bed) < 10) {
          return(list(status = "ERROR", msg = "Insufficient columns", sample = sample_id, sites = 0))
        }

        if (ncol(bed) == 18) {
          bed$V10 <- as.numeric(bed$V10)
          bed$V12 <- as.numeric(bed$V12)
          bed$V13 <- as.numeric(bed$V13)
          unmeth_col <- bed$V12
          meth_col <- bed$V13
        } else if (ncol(bed) == 13) {
          bed$V10 <- as.numeric(bed$V10)
          bed$V11 <- as.numeric(bed$V11)
          bed$V13 <- as.numeric(bed$V13)
          unmeth_col <- bed$V11
          meth_col <- bed$V13
        } else {
          return(list(status = "ERROR", msg = paste("Unexpected columns:", ncol(bed)), sample = sample_id, sites = 0))
        }

        bed <- bed[!is.na(bed$V10) & !is.na(unmeth_col) & !is.na(meth_col) & bed$V10 >= opt$min_coverage]

        new <- data.table(
          chrBase = paste(bed$V1, bed$V2, sep = "."),
          chr = as.character(bed$V1),
          start = as.integer(bed$V2),
          stop = as.integer(bed$V2) + 1,
          strand = as.character(bed$V6),
          coverage = as.integer(bed$V10),
          freqC = round(pmin(100, (meth_col / bed$V10) * 100), 3),
          freqT = round(pmin(100, (unmeth_col / bed$V10) * 100), 3)
        )

        new <- new[coverage > 0 & freqC >= 0 & freqC <= 100 & freqT >= 0 & freqT <= 100]
        new <- new[!duplicated(new$chrBase)]

        fwrite(new, outfile, sep = "\t", col.names = TRUE, quote = FALSE)

        return(list(status = "SUCCESS", msg = "OK", sample = sample_id, sites = nrow(new)))
      },
      error = function(e) {
        return(list(status = "ERROR", msg = e$message, sample = sample_id, sites = 0))
      }
    )
  }

  if (opt$no_parallel) {
    results <- lapply(files_to_convert, convert_function)
  } else {
    results <- mclapply(files_to_convert, convert_function, mc.cores = opt$cores)
  }

  success_count <- sum(sapply(results, function(x) x$status == "SUCCESS"))
  error_count <- sum(sapply(results, function(x) x$status == "ERROR"))

  cat("\nConversion Summary:\n")
  cat("  Success: ", success_count, "/", length(results), " files\n")
  cat("  Failed:  ", error_count, "/", length(results), " files\n\n")

  if (success_count > 0) {
    cat("Sample Details (first 5):\n")
    for (result in head(results[sapply(results, function(x) x$status == "SUCCESS")], 5)) {
      cat(sprintf("  %-20s %15s CpG sites\n", result$sample, format(result$sites, big.mark = ",")))
    }
    if (success_count > 5) cat("  ... and", success_count - 5, "more\n")
    cat("\n")
  }

  if (error_count > 0) {
    cat("Errors:\n")
    for (result in results[sapply(results, function(x) x$status == "ERROR")]) {
      cat("  ", result$sample, ":", result$msg, "\n")
    }
    cat("\n")
  }
}

# ================================================================================
# STEP 3: CREATE METHYLRAW OBJECTS
# ================================================================================

cat("================================================================================\n")
cat("STEP 3: Creating methylKit Objects\n")
cat("================================================================================\n\n")

existing_converted <- converted_files[file.exists(converted_files)]
existing_indices <- which(file.exists(converted_files))
existing_treatment <- treatment_vector[existing_indices]
existing_sample_ids <- sheet$ID[existing_indices]

if (length(existing_converted) == 0) {
  stop("No converted bedFile files available for loading!", call. = FALSE)
}

cat("Loading", length(existing_converted), "samples into methylKit...\n")
cat("Assembly:", assembly_name, "\n")
cat("Context:  CpG\n")
cat("Min Cov: ", opt$min_coverage, "\n\n")

file_list <- as.list(existing_converted)
sample_list <- as.list(as.character(existing_sample_ids))

meth_obj <- tryCatch(
  {
    methylKit::methRead(
      file_list,
      sample.id = sample_list,
      assembly = assembly_name,
      treatment = as.integer(existing_treatment),
      mincov = opt$min_coverage,
      resolution = "base",
      pipeline = list(
        fraction = FALSE,
        chr.col = 2,
        start.col = 3,
        end.col = 4,
        coverage.col = 6,
        strand.col = 5,
        freqC.col = 7
      )
    )
  },
  error = function(e) {
    cat("Primary method failed, using fallback...\n")
    cat("Error:", e$message, "\n\n")

    meth_list <- list()
    successful_treatments <- c()
    successful_sample_ids <- c()

    for (i in seq_along(existing_converted)) {
      file_path <- existing_converted[i]
      sample_id <- existing_sample_ids[i]
      treatment_val <- existing_treatment[i]

      tryCatch(
        {
          dt <- fread(file_path, sep = "\t", header = TRUE, showProgress = FALSE)
          dt <- dt[dt$coverage >= opt$min_coverage]
          dt <- dt[complete.cases(dt)]
          dt <- dt[dt$coverage > 0 & dt$freqC >= 0 & dt$freqC <= 100 & dt$freqT >= 0 & dt$freqT <= 100]
          dt <- dt[!duplicated(dt$chrBase)]

          if (nrow(dt) == 0) stop("No valid CpG sites")

          methyl_data <- data.frame(
            chr = as.character(dt$chr),
            start = as.integer(dt$start),
            end = as.integer(dt$stop),
            strand = as.character(dt$strand),
            coverage = as.integer(dt$coverage),
            numCs = as.integer(round(dt$coverage * dt$freqC / 100)),
            numTs = as.integer(round(dt$coverage * dt$freqT / 100)),
            stringsAsFactors = FALSE
          )

          methyl_data <- methyl_data[methyl_data$coverage > 0 &
            methyl_data$numCs >= 0 &
            methyl_data$numTs >= 0 &
            (methyl_data$numCs + methyl_data$numTs) <= methyl_data$coverage * 1.1, ]

          if (nrow(methyl_data) == 0) stop("No valid data")

          meth_obj_single <- new("methylRaw",
            methyl_data,
            sample.id = as.character(sample_id),
            assembly = assembly_name,
            context = "CpG",
            resolution = "base"
          )

          meth_list <- append(meth_list, list(meth_obj_single))
          successful_treatments <- c(successful_treatments, treatment_val)
          successful_sample_ids <- c(successful_sample_ids, sample_id)

          cat("  ", sample_id, ":", format(nrow(meth_obj_single), big.mark = ","), "sites\n")
        },
        error = function(e) {
          cat("  ", sample_id, ": FAILED -", e$message, "\n")
        }
      )
    }

    if (length(meth_list) == 0) {
      stop("All creation methods failed!", call. = FALSE)
    }

    new("methylRawList", meth_list, treatment = as.integer(successful_treatments))
  }
)

cat("\nObject created successfully!\n\n")

# ================================================================================
# STEP 4: SUMMARY STATISTICS
# ================================================================================

cat("================================================================================\n")
cat("STEP 4: Summary Statistics\n")
cat("================================================================================\n\n")

actual_treatment <- getTreatment(meth_obj)
sample_sizes <- sapply(meth_obj, nrow)

cat("Dataset Overview:\n")
cat("  Total samples:     ", length(meth_obj), "\n")
cat("  Treatment groups:  ", paste(names(table(actual_treatment)), "=",
  table(actual_treatment),
  collapse = ", "
), "\n\n")

cat("CpG Sites per Sample:\n")
cat("  Minimum:           ", format(min(sample_sizes), big.mark = ","), "\n")
cat("  Maximum:           ", format(max(sample_sizes), big.mark = ","), "\n")
cat("  Mean:              ", format(round(mean(sample_sizes)), big.mark = ","), "\n")
cat("  Median:            ", format(round(median(sample_sizes)), big.mark = ","), "\n\n")

# ================================================================================
# STEP 5: SAVE OUTPUTS
# ================================================================================

cat("================================================================================\n")
cat("STEP 5: Saving Processed Data\n")
cat("================================================================================\n\n")

outfile <- file.path(
  opt$output_dir,
  paste0(
    today, "_", gsub(" ", "_", opt$species),
    "_n", length(meth_obj), "_methylRaw.qs"
  )
)
qsave(meth_obj, outfile)
cat("methylRaw object: ", basename(outfile), "\n")

summary_file <- file.path(
  opt$output_dir,
  paste0(
    today, "_", gsub(" ", "_", opt$species),
    "_preprocessing_summary.txt"
  )
)
writeLines(c(
  "================================================================================",
  "METHYLSENSE DATA PREPROCESSING SUMMARY",
  "================================================================================",
  paste0("Generated:              ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Version:                ", SCRIPT_VERSION),
  paste0("Species:                ", opt$species),
  paste0("Assembly:               ", assembly_name),
  paste0("Total samples:          ", length(meth_obj)),
  paste0("Treatment groups:       ", paste(names(table(actual_treatment)), "=",
    table(actual_treatment),
    collapse = ", "
  )),
  paste0("Min CpG sites:          ", format(min(sample_sizes), big.mark = ",")),
  paste0("Max CpG sites:          ", format(max(sample_sizes), big.mark = ",")),
  paste0("Mean CpG sites:         ", format(round(mean(sample_sizes)), big.mark = ",")),
  paste0("Median CpG sites:       ", format(round(median(sample_sizes)), big.mark = ",")),
  paste0("Min coverage filter:    ", opt$min_coverage),
  paste0("Output file:            ", basename(outfile)),
  "================================================================================",
  "",
  "Ready for downstream analysis:",
  "  - DMR identification with methylKit",
  "  - Machine learning classification",
  "  - Epigenetic biomarker discovery",
  "================================================================================"
), summary_file)
cat("Summary report:    ", basename(summary_file), "\n\n")

# ================================================================================
# COMPLETION
# ================================================================================

cat("================================================================================\n")
cat("PREPROCESSING COMPLETE\n")
cat("================================================================================\n\n")
cat("Status:            SUCCESS\n")
cat("Output directory:  ", opt$output_dir, "\n")
cat("Next step:         Run MethylSense analysis pipeline\n")
cat("\n")
cat("================================================================================\n\n")
