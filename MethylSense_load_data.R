#!/usr/bin/env Rscript

# ================================================================================
#                    METHYLSENSE DATA PREPROCESSING PIPELINE
# ================================================================================
# Purpose: Read modkit bedMethyl files (13 or 18 columns) directly into methylKit
#          methylRaw objects for downstream DMR and ML analysis, then verify the
#          saved object against the raw files site by site.
#
# Author: Markus Hodal Drag
# Version: 6.1.4
# Release Date: 2026-11-12
# GitHub: https://github.com/markusdrag/MethylSense
#
# Citation:
#   Drag MH, Hvilsom C, Poulsen LL, Jensen HE, Tahas SA, Leineweber C,
#   Cray C, Bertelsen MF, Bojesen AM.
#   MethylSense: high accuracy machine learning-based diagnostics for
#   Aspergillus fumigatus infection in chickens using host cell-free DNA
#   methylation and Nanopore sequencing. J Clin Microbiol 0:e01054-25.
#   https://doi.org/10.1128/jcm.01054-25
#
# Version history:
#   6.1.4 (2026-11-12)
#     - --traditional keeps the first count when a chromosome and start is
#       repeated, and also writes the legacy nonD files. The default load still
#       adds differing counts together.
#   6.1.3 (2026-11-12)
#     - --write_nond writes a legacy 8-column nonD file per sample. freqC is the
#       methylated percent from column 12. A repeated chromosome and start keeps
#       the first row. The files are not read back.
#   6.1.2 (2026-11-12)
#     - Repeated rows for the same chromosome, start and strand are collapsed.
#       Identical copies are kept once. Rows with different counts are added.
#   6.1.1 (2026-11-12)
#     - --min_coverage may be 0, which keeps CpG sites with no reads.
#   6.1.0 (2026-11-12)
#     - bedMethyl files are read directly. No intermediate _CpG_nonD.txt files are
#       written or read, so stale converted files can no longer enter a run.
#     - Methylated count is modkit column 12 (Nmod), unmethylated count is column
#       13 (Ncanonical). Each file must confirm this through its own column 11
#       (percent or fraction modified), or the run stops.
#     - Built-in QC: the saved .qs2 is read back from disk and compared with every
#       site of every raw bedMethyl file. Results go to _load_QC.csv and
#       _load_QC_report.txt. A failing object is renamed _QC_FAILED.qs2.
#     - Rows are filtered to one modification code (--mod_code, default auto).
#     - Duplicate sites, missing files and unreadable files stop the run.
#     - --force_convert is accepted and ignored.
#     - Script paths that contain spaces now work.
# ================================================================================

SCRIPT_VERSION <- "6.1.4"
SCRIPT_DATE <- "2026-11-12"

suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
  library(methylKit)
  library(data.table)
  library(qs2)
  library(parallel)
})

cmd_args_qs <- commandArgs(trailingOnly = FALSE)
file_arg_qs <- grep("^--file=", cmd_args_qs, value = TRUE)
script_dir_qs <- if (length(file_arg_qs)) {
  # Rscript writes spaces in the script path as "~+~".
  script_path_qs <- gsub("~+~", " ", sub("^--file=", "", file_arg_qs[1]), fixed = TRUE)
  dirname(normalizePath(script_path_qs, mustWork = FALSE))
} else {
  getwd()
}
qs_io_candidates <- c(
  file.path(script_dir_qs, "MethylSense_qs_io.R"),
  file.path(script_dir_qs, "..", "MethylSense_qs_io.R")
)
qs_io_hit <- qs_io_candidates[file.exists(qs_io_candidates)]
if (length(qs_io_hit) == 0L) {
  stop("Cannot find MethylSense_qs_io.R next to this script or in the parent directory.")
}
source(normalizePath(qs_io_hit[1]))

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
    help = "Path to sample sheet (Excel or CSV). Needs Species, ID, bedFileOrg and treatMethylkit columns",
    metavar = "FILE"
  ),
  make_option(c("-b", "--bed_dir"),
    type = "character",
    help = "Directory containing the modkit bedMethyl files named in bedFileOrg",
    metavar = "DIRECTORY"
  ),
  make_option(c("-o", "--output_dir"),
    type = "character",
    help = "Output directory for processed objects and QC files",
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
    help = "Minimum coverage (Nmod + Ncanonical) per CpG. 0 keeps sites with no reads [default: %default]",
    metavar = "INTEGER"
  ),
  make_option("--mod_code",
    type = "character",
    default = "auto",
    help = paste(
      "Modification code in bedMethyl column 4 to load. 'auto' uses the only code",
      "in the file, or 'm' (5mC) when several codes are present. Use 'C' for files",
      "made with --combine-mods. [default: %default]"
    ),
    metavar = "CHARACTER"
  ),
  make_option("--mapq_column",
    type = "character",
    default = NULL,
    help = paste(
      "Name of the sample-sheet column holding each sample's mean mapping",
      "quality (MAPQ). Combined with --min_mapq, samples below the threshold",
      "are dropped before any file is read. [default: no MAPQ filtering]"
    ),
    metavar = "CHARACTER"
  ),
  make_option("--min_mapq",
    type = "double",
    default = NULL,
    help = paste(
      "Minimum mean mapping quality to retain a sample (keeps samples with",
      "MAPQ >= this value, i.e. equal or better). Requires --mapq_column.",
      "[default: no MAPQ filtering]"
    ),
    metavar = "NUMERIC"
  ),
  make_option("--allow_missing",
    action = "store_true",
    default = FALSE,
    help = "Continue without samples whose bedMethyl file is missing (they are listed in the QC report)"
  ),
  make_option("--force_convert",
    action = "store_true",
    default = FALSE,
    help = "Ignored since 6.1.0, which writes no converted files. Kept so older job scripts still run"
  ),
  make_option("--write_nond",
    action = "store_true",
    default = FALSE,
    help = paste(
      "Also write one legacy 8-column nonD file per sample into --output_dir.",
      "freqC is the methylated percent (column 12) and freqT is the unmethylated",
      "percent (column 13). If the same chromosome and start appears more than once,",
      "the first row is kept and later rows are dropped. The qs2 object is unchanged.",
      "These nonD files are not read back."
    )
  ),
  make_option("--traditional",
    action = "store_true",
    default = FALSE,
    help = paste(
      "Build the qs2 object the way version 6.0.3 did when a site is repeated:",
      "keep the first count for a chromosome and start, and drop later rows.",
      "Also writes the legacy nonD files. Without this flag, differing counts are added."
    )
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
  description = "\nMethylSense Data Preprocessing Pipeline\nReads modkit bedMethyl files into a methylRawList and verifies it against the raw files",
  epilogue = "Example:\n  Rscript MethylSense_load_data.R --species 'Minipig' --sample_sheet samples.xlsx --bed_dir ./beds --output_dir ./processed"
)

opt <- parse_args(parser)
if (isTRUE(opt$traditional)) opt$write_nond <- TRUE

if (is.null(opt$species) || is.null(opt$sample_sheet) ||
  is.null(opt$bed_dir) || is.null(opt$output_dir)) {
  print_help(parser)
  stop("\nError: Missing required arguments.\n", call. = FALSE)
}

if (!is.null(opt$min_mapq) && is.null(opt$mapq_column)) {
  stop("\nError: --min_mapq requires --mapq_column so the script knows which sample-sheet column to read.\n", call. = FALSE)
}
if (!is.null(opt$mapq_column) && is.null(opt$min_mapq)) {
  stop("\nError: --mapq_column requires --min_mapq (the minimum mean mapping quality to keep).\n", call. = FALSE)
}
if (is.na(opt$min_coverage) || opt$min_coverage < 0) {
  stop("\nError: --min_coverage must be 0 or higher.\n", call. = FALSE)
}

# ================================================================================
# HELPER: SAMPLE-LEVEL MAPQ FILTERING
# ================================================================================
# Drops samples whose per-sample mean mapping quality (read from a user-named
# column in the sample sheet) is below `min_mapq`. Retains samples with
# MAPQ >= min_mapq (equal or better). Returns the sheet unchanged when either
# argument is not supplied.

filter_samples_by_mapq <- function(sheet, mapq_column = NULL, min_mapq = NULL) {
  if (is.null(mapq_column) || is.null(min_mapq)) {
    return(sheet)
  }

  if (!mapq_column %in% colnames(sheet)) {
    stop(
      "MAPQ column '", mapq_column, "' not found in sample sheet.\n",
      "  Available columns: ", paste(colnames(sheet), collapse = ", "),
      call. = FALSE
    )
  }

  mapq_values <- suppressWarnings(as.numeric(sheet[[mapq_column]]))
  n_before <- nrow(sheet)
  n_missing <- sum(is.na(mapq_values))

  keep <- !is.na(mapq_values) & mapq_values >= min_mapq
  dropped <- sheet[!keep, , drop = FALSE]
  sheet_filtered <- sheet[keep, , drop = FALSE]

  cat("MAPQ filtering:\n")
  cat("  Column:            ", mapq_column, "\n")
  cat("  Threshold:          MAPQ >=", min_mapq, "\n")
  cat("  Samples before:    ", n_before, "\n")
  cat("  Samples retained:  ", nrow(sheet_filtered), "\n")
  cat("  Samples dropped:   ", n_before - nrow(sheet_filtered), "\n")
  if (n_missing > 0) {
    cat("  Note: ", n_missing, "sample(s) had missing/non-numeric MAPQ and were dropped\n")
  }
  if (nrow(dropped) > 0) {
    drop_ids <- if ("ID" %in% colnames(dropped)) dropped$ID else as.character(seq_len(nrow(dropped)))
    cat("  Dropped samples:   ", paste(head(drop_ids, 10), collapse = ", "))
    if (length(drop_ids) > 10) cat(", ... and", length(drop_ids) - 10, "more")
    cat("\n")
  }
  cat("\n")

  if (nrow(sheet_filtered) == 0) {
    stop("No samples remain after MAPQ filtering (threshold MAPQ >= ", min_mapq, ").", call. = FALSE)
  }

  sheet_filtered
}

# ================================================================================
# HELPERS: BEDMETHYL READING
# ================================================================================
# modkit bedMethyl columns (https://github.com/nanoporetech/modkit#description-of-bedmethyl-output):
#   1 chrom, 2 start (0-based), 3 end, 4 modification code, 5 score, 6 strand,
#   7-9 display fields, 10 Nvalid_cov, 11 fraction modified (Nmod / Nvalid_cov,
#   written as 0-100 by older modkit and 0-1 by newer modkit), 12 Nmod,
#   13 Ncanonical, 14 Nother_mod, 15 Ndelete, 16 Nfail, 17 Ndiff, 18 Nnocall.

bed_stream_cmd <- function(path) {
  src <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    paste("gzip -dc", shQuote(path))
  } else {
    paste("cat", shQuote(path))
  }
  paste(src, "| tr ' ' '\\t'")
}

read_bedmethyl <- function(path) {
  probe <- fread(path, nrows = 1, header = FALSE, sep = "\t", showProgress = FALSE)
  space_separated <- FALSE
  # modkit before 0.2 separated columns 10 to 18 with spaces.
  if (ncol(probe) == 10L && is.character(probe[[10]]) && grepl(" ", probe[[10]][1])) {
    space_separated <- TRUE
    probe <- fread(cmd = paste(bed_stream_cmd(path), "| head -n 1"), header = FALSE, sep = "\t", showProgress = FALSE)
  }

  n_col <- ncol(probe)
  if (!n_col %in% c(13L, 18L)) {
    stop(sprintf("expected a modkit bedMethyl file with 13 or 18 columns, found %d", n_col))
  }
  if (!is.numeric(probe[[2]]) || !is.numeric(probe[[10]]) || !is.numeric(probe[[12]]) || !is.numeric(probe[[13]])) {
    stop("columns 2, 10, 12 and 13 are not numeric. The first row may be a header, which modkit does not write")
  }

  sel <- c(1, 2, 3, 4, 6, 10, 11, 12, 13, if (n_col == 18L) 14)
  nms <- c("chr", "start", "end", "code", "strand", "n_valid", "col11", "n_mod", "n_canonical", if (n_col == 18L) "n_other")

  bed <- if (space_separated) {
    fread(cmd = bed_stream_cmd(path), header = FALSE, sep = "\t", select = sel, col.names = nms, showProgress = FALSE)
  } else {
    fread(path, header = FALSE, sep = "\t", select = sel, col.names = nms, showProgress = FALSE)
  }
  if (n_col == 13L) bed[, n_other := NA_integer_]

  bed[, chr := as.character(chr)]
  bed[, code := as.character(code)]
  bed[, strand := as.character(strand)]
  for (col in c("start", "end", "n_valid", "n_mod", "n_canonical")) {
    if (anyNA(bed[[col]])) stop(sprintf("column %s has missing or non-numeric values", col))
  }

  list(bed = bed, n_col = n_col, space_separated = space_separated)
}

choose_mod_code <- function(codes, requested) {
  codes <- sort(unique(codes))
  if (!identical(requested, "auto")) {
    if (!requested %in% codes) {
      stop(sprintf("--mod_code %s is not in column 4. Codes in file: %s", requested, paste(codes, collapse = ", ")))
    }
    return(requested)
  }
  if (length(codes) == 1L) return(codes)
  if ("m" %in% codes) return("m")
  stop(sprintf(
    "column 4 holds several modification codes (%s). Set --mod_code, for example m for 5mC",
    paste(codes, collapse = ", ")
  ))
}

# Column 11 is modkit's own statement of the modified fraction. In modkit output it
# equals column 12 / column 10 exactly, which proves column 12 is the modified count.
# Simulated files may store a true level in column 11 instead. These are accepted
# with a warning when column 11 is still closer to the modified than to the
# unmodified fraction at nearly every site. A swapped file is closer at almost none.
DIRECTION_EXACT_MIN <- 0.99
DIRECTION_CLOSER_MIN <- 0.95

check_direction <- function(bed) {
  scale <- if (max(bed$col11, na.rm = TRUE) > 1 + 1e-9) 100 else 1
  informative <- bed[n_valid > 0 & n_mod != n_canonical]
  if (nrow(informative) == 0L) {
    return(list(scale = scale, agree_mod = NA_real_, agree_can = NA_real_, closer_mod = NA_real_, n_informative = 0L))
  }
  tol <- 0.0051 * scale
  dist_mod <- abs(informative$col11 - scale * informative$n_mod / informative$n_valid)
  dist_can <- abs(informative$col11 - scale * informative$n_canonical / informative$n_valid)
  list(
    scale = scale,
    agree_mod = mean(dist_mod <= tol),
    agree_can = mean(dist_can <= tol),
    closer_mod = mean(dist_mod < dist_can),
    n_informative = nrow(informative)
  )
}

direction_confirmed <- function(agree_mod, closer_mod) {
  is.na(agree_mod) || agree_mod >= DIRECTION_EXACT_MIN || closer_mod >= DIRECTION_CLOSER_MIN
}

methylkit_strand <- function(strand) {
  fifelse(strand == "-", "-", "+")
}

# A merged or concatenated bedMethyl file can list the same CpG more than once.
# An exact copy must not be counted twice. Two rows with different counts are
# separate pileups of the same site, so their counts are added.
collapse_bed_sites <- function(bed) {
  keys <- c("chr", "start", "strand")
  n_rows <- nrow(bed)
  bed <- unique(bed, by = c(keys, "n_mod", "n_canonical", "n_valid", "n_other"))
  rows_exact_duplicate <- n_rows - nrow(bed)
  sites_summed <- nrow(bed[, .N, by = keys][N > 1L])
  if (sites_summed > 0L) {
    bed <- bed[, .(
      n_valid = sum(as.numeric(n_valid)),
      n_mod = sum(as.numeric(n_mod)),
      n_canonical = sum(as.numeric(n_canonical)),
      n_other = if (all(is.na(n_other))) NA_real_ else sum(as.numeric(n_other), na.rm = TRUE)
    ), by = keys]
    too_big <- bed$n_mod > .Machine$integer.max |
      bed$n_canonical > .Machine$integer.max |
      bed$n_valid > .Machine$integer.max
    if (any(too_big, na.rm = TRUE)) stop("collapsed read counts exceed the integer limit")
  }
  list(bed = bed, rows_exact_duplicate = rows_exact_duplicate, sites_summed = sites_summed)
}

# Legacy 8-column file: chrBase, chr, start, stop, strand, coverage, freqC, freqT.
# freqC is the methylated percent from column 12. Coverage is column 10.
# chrBase is chromosome.start. The first row for that key is kept.
write_nond_file <- function(bed, path, min_coverage) {
  use <- bed[n_valid >= min_coverage & n_valid > 0]
  out <- data.table(
    chrBase = paste(use$chr, use$start, sep = "."),
    chr = use$chr,
    start = as.integer(use$start),
    stop = as.integer(use$start) + 1L,
    strand = use$strand,
    coverage = as.integer(use$n_valid),
    freqC = round(pmin(100, use$n_mod / use$n_valid * 100), 3),
    freqT = round(pmin(100, use$n_canonical / use$n_valid * 100), 3)
  )
  out <- out[freqC >= 0 & freqC <= 100 & freqT >= 0 & freqT <= 100]
  n_before <- nrow(out)
  out <- out[!duplicated(chrBase)]
  fwrite(out, path, sep = "\t", col.names = TRUE, quote = FALSE)
  list(
    sites = nrow(out),
    zero_coverage_dropped = sum(bed$n_valid == 0),
    chrbase_dropped = n_before - nrow(out)
  )
}

# Same rule as the legacy nonD file: one row per chromosome and start.
# The first row in the file is kept. Later rows, including the other strand, are dropped.
keep_first_chr_start <- function(bed, min_coverage) {
  use <- bed[n_valid >= min_coverage & n_valid > 0]
  n_before <- nrow(use)
  use <- use[!duplicated(use, by = c("chr", "start"))]
  list(bed = use, rows_dropped = n_before - nrow(use))
}

# ================================================================================
# SETUP
# ================================================================================

if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}

today <- format(Sys.Date(), "%Y%m%d")

if (is.null(opt$assembly) || opt$assembly == "") {
  species_clean <- gsub(" ", "_", opt$species)
  species_clean <- gsub("[^A-Za-z0-9_]", "", species_clean)
  assembly_name <- paste0(species_clean, "_assembly")
} else {
  assembly_name <- opt$assembly
}

use_cores <- if (opt$no_parallel) 1L else max(1L, opt$cores)
run_parallel <- function(x, fun) {
  if (use_cores == 1L) lapply(x, fun) else mclapply(x, fun, mc.cores = use_cores, mc.preschedule = FALSE)
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
cat("                      Data preprocessing pipeline                              \n")
cat("================================================================================\n")
cat(paste("Version:", SCRIPT_VERSION, "|", "Release Date:", SCRIPT_DATE, "\n"))
cat(paste("Start Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
cat(paste("Working Directory:", getwd(), "\n"))
cat(paste("R Version:", R.version.string, "\n"))
cat(paste("methylKit Version:", as.character(packageVersion("methylKit")), "\n"))
cat("\n")
cat("Please cite:\n")
cat("  Drag MH, Hvilsom C, Poulsen LL, Jensen HE, Tahas SA, Leineweber C,\n")
cat("  Cray C, Bertelsen MF, Bojesen AM.\n")
cat("  MethylSense: high accuracy machine learning-based diagnostics for\n")
cat("  Aspergillus fumigatus infection in chickens using host cell-free DNA\n")
cat("  methylation and Nanopore sequencing. J Clin Microbiol 0:e01054-25.\n")
cat("  https://doi.org/10.1128/jcm.01054-25\n")
cat("================================================================================\n\n")

cat("Configuration:\n")
cat("  Species:            ", opt$species, "\n")
cat("  Assembly:           ", assembly_name, "\n")
cat("  Sample Sheet:       ", basename(opt$sample_sheet), "\n")
cat("  BED Directory:      ", opt$bed_dir, "\n")
cat("  Output Directory:   ", opt$output_dir, "\n")
cat("  CPU Cores:          ", use_cores, "\n")
cat("  Min Coverage:       ", opt$min_coverage, " reads per CpG (Nmod + Ncanonical)\n")
cat("  Modification code:  ", opt$mod_code, "\n")
if (!is.null(opt$mapq_column) && !is.null(opt$min_mapq)) {
  cat("  MAPQ Filter:         keep", opt$mapq_column, ">=", opt$min_mapq, "\n")
} else {
  cat("  MAPQ Filter:         Disabled\n")
}
cat("  Missing files:      ", ifelse(opt$allow_missing, "Skipped and reported", "Stop the run"), "\n")
if (opt$force_convert) {
  cat("  Note: --force_convert is ignored. Version 6.1.0 reads bedMethyl files directly.\n")
}
cat("\n")

# ================================================================================
# STEP 1: LOAD AND VALIDATE SAMPLE SHEET
# ================================================================================

cat("================================================================================\n")
cat("Step 1: Loading sample sheet\n")
cat("================================================================================\n\n")

sheet_all <- tryCatch(
  {
    if (grepl("\\.csv$", opt$sample_sheet, ignore.case = TRUE)) {
      read.csv(opt$sample_sheet, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
    } else {
      as.data.frame(read_xlsx(opt$sample_sheet, col_types = "text"))
    }
  },
  error = function(e) {
    stop("Failed to read sample sheet: ", e$message, call. = FALSE)
  }
)

required_cols <- c("Species", "ID", "bedFileOrg", "treatMethylkit")
missing_cols <- setdiff(required_cols, colnames(sheet_all))
if (length(missing_cols) > 0) {
  stop("Sample sheet is missing column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

sheet <- sheet_all[!is.na(sheet_all$Species) & sheet_all$Species == opt$species, , drop = FALSE]

if (nrow(sheet) == 0) {
  stop("No samples found for species: ", opt$species, "\n",
    "Available species: ", paste(unique(sheet_all$Species), collapse = ", "),
    call. = FALSE
  )
}

cat("Sample Sheet Summary:\n")
cat("  Total samples in sheet: ", nrow(sheet_all), "\n")
cat("  Samples for", opt$species, ":", nrow(sheet), "\n")
cat("  Input column:            bedFileOrg (modkit bedMethyl)\n")
if ("bedFile" %in% colnames(sheet)) {
  cat("  Note: the bedFile column is not read. The object is built from bedFileOrg.\n")
}
if (isTRUE(opt$write_nond)) {
  cat("  Legacy nonD files:   one 8-column file per sample in the output directory\n")
}
cat("\n")

n_samples_pre_mapq <- nrow(sheet)
sheet <- filter_samples_by_mapq(sheet, opt$mapq_column, opt$min_mapq)

sheet$ID <- as.character(sheet$ID)
dup_ids <- unique(sheet$ID[duplicated(sheet$ID)])
if (length(dup_ids) > 0) {
  stop("Duplicate sample IDs in sample sheet: ", paste(dup_ids, collapse = ", "), call. = FALSE)
}
treat_raw <- trimws(as.character(sheet$treatMethylkit))
treat_blank <- is.na(treat_raw) | treat_raw == "" | toupper(treat_raw) == "NA"
treat_num <- suppressWarnings(as.numeric(treat_raw))
treat_bad <- !treat_blank & (is.na(treat_num) | treat_num != round(treat_num))
if (any(treat_bad)) {
  stop("treatMethylkit must be a whole number, or NA for unknown samples. Check: ",
    paste(sheet$ID[treat_bad], collapse = ", "),
    call. = FALSE
  )
}
treat_num[treat_blank] <- NA
sheet$treatMethylkit <- as.integer(treat_num)
if (any(treat_blank)) {
  cat("  Note:", sum(treat_blank), "sample(s) have no treatMethylkit value and are stored with NA treatment\n\n")
}

bed_paths <- file.path(opt$bed_dir, sheet$bedFileOrg)
file_found <- file.exists(bed_paths)
missing_samples <- sheet$ID[!file_found]

cat("File Status:\n")
cat("  Found:   ", sum(file_found), "bedMethyl files\n")
cat("  Missing: ", sum(!file_found), "bedMethyl files\n")
if (any(!file_found)) {
  for (k in head(which(!file_found), 10)) cat("    -", sheet$ID[k], ":", bed_paths[k], "\n")
  if (sum(!file_found) > 10) cat("    - ... and", sum(!file_found) - 10, "more\n")
  if (!opt$allow_missing) {
    stop("Missing bedMethyl file(s). Fix the sample sheet or rerun with --allow_missing.", call. = FALSE)
  }
}
cat("\n")

sheet <- sheet[file_found, , drop = FALSE]
bed_paths <- bed_paths[file_found]
if (nrow(sheet) == 0) stop("No bedMethyl files found.", call. = FALSE)

# ================================================================================
# STEP 2: READ BEDMETHYL FILES INTO METHYLRAW OBJECTS
# ================================================================================

cat("================================================================================\n")
cat("Step 2: Reading bedMethyl files into methylKit objects\n")
cat("================================================================================\n\n")
cat("Methylated reads = column 12 (Nmod), unmethylated reads = column 13 (Ncanonical)\n")
cat("Reading", nrow(sheet), "files...\n\n")

build_sample <- function(i) {
  path <- bed_paths[i]
  sample_id <- sheet$ID[i]
  tryCatch(
    {
      rd <- read_bedmethyl(path)
      bed <- rd$bed
      rows_in_file <- nrow(bed)

      code_used <- choose_mod_code(bed$code, opt$mod_code)
      rows_other_code <- sum(bed$code != code_used)
      bed <- bed[code == code_used]

      dirn <- check_direction(bed)
      if (!direction_confirmed(dirn$agree_mod, dirn$closer_mod)) {
        stop(sprintf(
          paste(
            "column 11 does not match Nmod / Nvalid_cov (exact on %.2f%% of informative sites,",
            "closer to Nmod than to Ncanonical on %.2f%%). Column 12 may not be the modified count,",
            "so the methylation direction cannot be confirmed"
          ),
          100 * dirn$agree_mod, 100 * dirn$closer_mod
        ))
      }

      if (rd$n_col == 18L) {
        count_identity <- mean(bed$n_mod + bed$n_canonical + bed$n_other == bed$n_valid)
        n_other_total <- sum(as.numeric(bed$n_other))
      } else {
        count_identity <- mean(bed$n_mod + bed$n_canonical <= bed$n_valid)
        n_other_total <- NA_real_
      }
      if (count_identity < 0.999) {
        stop(sprintf("Nmod + Ncanonical + Nother_mod equals Nvalid_cov on only %.2f%% of rows", 100 * count_identity))
      }

      rows_unstranded <- sum(!bed$strand %in% c("+", "-"))
      bed[, strand := methylkit_strand(strand)]

      nond_sites <- NA_integer_
      nond_chrbase_dropped <- NA_integer_
      nond_zero_dropped <- NA_integer_
      if (isTRUE(opt$write_nond)) {
        nond_path <- file.path(opt$output_dir, paste0(basename(path), "_CpG_nonD.txt"))
        nond <- write_nond_file(bed, nond_path, opt$min_coverage)
        nond_sites <- nond$sites
        nond_chrbase_dropped <- nond$chrbase_dropped
        nond_zero_dropped <- nond$zero_coverage_dropped
      }

      if (isTRUE(opt$traditional)) {
        kept <- keep_first_chr_start(bed, opt$min_coverage)
        bed <- kept$bed
        traditional_dropped <- kept$rows_dropped
        bed[, coverage := n_mod + n_canonical]
        rows_zero_coverage <- nond_zero_dropped
        if (is.na(rows_zero_coverage)) rows_zero_coverage <- sum(bed$coverage == 0)
        rows_below_min <- 0L
        collapsed <- list(rows_exact_duplicate = 0L, sites_summed = 0L)
      } else {
        traditional_dropped <- 0L
        collapsed <- collapse_bed_sites(bed)
        bed <- collapsed$bed
        bed[, coverage := n_mod + n_canonical]
        rows_zero_coverage <- sum(bed$coverage == 0)
        rows_below_min <- sum(bed$coverage > 0 & bed$coverage < opt$min_coverage)
        bed <- bed[coverage >= opt$min_coverage]
      }
      if (nrow(bed) == 0L) stop("no sites left after the coverage filter")

      df <- data.frame(
        chr = bed$chr,
        start = as.integer(bed$start),
        end = as.integer(bed$start) + 1L,
        strand = bed$strand,
        coverage = as.integer(bed$coverage),
        numCs = as.integer(bed$n_mod),
        numTs = as.integer(bed$n_canonical),
        stringsAsFactors = FALSE
      )
      df <- df[order(df$chr, df$start, df$end, df$strand), ]
      row.names(df) <- NULL

      obj <- new("methylRaw", df,
        sample.id = sample_id,
        assembly = assembly_name,
        context = "CpG",
        resolution = "base"
      )

      qc <- data.table(
        sample = sample_id,
        bed_file = basename(path),
        n_columns = rd$n_col,
        space_separated = rd$space_separated,
        rows_in_file = rows_in_file,
        mod_code = code_used,
        rows_other_code = rows_other_code,
        col11_scale = ifelse(dirn$scale == 100, "percent", "fraction"),
        informative_sites = dirn$n_informative,
        col11_agrees_with_Nmod = round(dirn$agree_mod, 6),
        col11_agrees_with_Ncanonical = round(dirn$agree_can, 6),
        col11_closer_to_Nmod = round(dirn$closer_mod, 6),
        count_identity = round(count_identity, 6),
        n_other_mod_total = n_other_total,
        rows_unstranded = rows_unstranded,
        rows_zero_coverage = rows_zero_coverage,
        rows_below_min_coverage = rows_below_min,
        rows_exact_duplicate = collapsed$rows_exact_duplicate,
        sites_summed = collapsed$sites_summed,
        traditional_dropped = traditional_dropped,
        nond_sites = nond_sites,
        nond_chrbase_dropped = nond_chrbase_dropped,
        nond_zero_dropped = nond_zero_dropped,
        sites_loaded = nrow(df),
        median_coverage = median(df$coverage),
        pct_sites_cov5 = round(100 * mean(df$coverage >= 5), 3),
        pct_plus_strand = round(100 * mean(df$strand == "+"), 3),
        mean_meth_pct = {
          cov_sum <- sum(as.numeric(df$coverage))
          if (cov_sum > 0) round(100 * sum(as.numeric(df$numCs)) / cov_sum, 3) else NA_real_
        }
      )
      list(status = "OK", obj = obj, qc = qc)
    },
    error = function(e) list(status = "ERROR", sample = sample_id, file = path, msg = conditionMessage(e))
  )
}

built <- run_parallel(seq_len(nrow(sheet)), build_sample)

failed <- vapply(built, function(x) !is.list(x) || !identical(x$status, "OK"), logical(1))
if (any(failed)) {
  cat("Failed to read", sum(failed), "file(s):\n")
  for (k in which(failed)) {
    x <- built[[k]]
    msg <- if (is.list(x)) x$msg else as.character(x)
    cat("  ", sheet$ID[k], ":", msg, "\n")
  }
  stop("bedMethyl reading failed. No object was written.", call. = FALSE)
}

build_qc <- rbindlist(lapply(built, `[[`, "qc"))
objs <- lapply(built, `[[`, "obj")
rm(built)
invisible(gc())

meth_obj <- new("methylRawList", objs, treatment = as.integer(sheet$treatMethylkit))
rm(objs)
invisible(gc())

for (k in seq_len(nrow(build_qc))) {
  cat(sprintf(
    "  %-24s %12s CpG sites   code %s   mean methylation %6.2f%%\n",
    build_qc$sample[k], format(build_qc$sites_loaded[k], big.mark = ","),
    build_qc$mod_code[k], build_qc$mean_meth_pct[k]
  ))
  if (build_qc$rows_exact_duplicate[k] > 0 || build_qc$sites_summed[k] > 0) {
    cat(sprintf(
      "    repeated rows collapsed: %s identical extra %s dropped, %s %s summed\n",
      format(build_qc$rows_exact_duplicate[k], big.mark = ","),
      ifelse(build_qc$rows_exact_duplicate[k] == 1, "row", "rows"),
      format(build_qc$sites_summed[k], big.mark = ","),
      ifelse(build_qc$sites_summed[k] == 1, "site", "sites")
    ))
  }
}
cat("\nObject created successfully!\n\n")

# ================================================================================
# STEP 3: SUMMARY STATISTICS
# ================================================================================

cat("================================================================================\n")
cat("Step 3: Summary statistics\n")
cat("================================================================================\n\n")

actual_treatment <- getTreatment(meth_obj)
sample_sizes <- sapply(meth_obj, nrow)
treatment_table <- table(actual_treatment, useNA = "ifany")
treatment_text <- paste(ifelse(is.na(names(treatment_table)), "NA", names(treatment_table)), "=",
  as.integer(treatment_table),
  collapse = ", "
)

cat("Dataset Overview:\n")
cat("  Total samples:     ", length(meth_obj), "\n")
cat("  Treatment groups:  ", treatment_text, "\n\n")

cat("CpG Sites per Sample:\n")
cat("  Minimum:           ", format(min(sample_sizes), big.mark = ","), "\n")
cat("  Maximum:           ", format(max(sample_sizes), big.mark = ","), "\n")
cat("  Mean:              ", format(round(mean(sample_sizes)), big.mark = ","), "\n")
cat("  Median:            ", format(round(median(sample_sizes)), big.mark = ","), "\n\n")

# ================================================================================
# STEP 4: SAVE OUTPUTS
# ================================================================================

cat("================================================================================\n")
cat("Step 4: Saving processed data\n")
cat("================================================================================\n\n")

file_prefix <- paste0(today, "_", gsub(" ", "_", opt$species))
outfile <- file.path(opt$output_dir, paste0(file_prefix, "_n", length(meth_obj), "_methylRaw.qs2"))
qs2::qs_save(meth_obj, outfile)
cat("methylRaw object: ", basename(outfile), "\n\n")

expected_ids <- sheet$ID
expected_treatment <- as.integer(sheet$treatMethylkit)
rm(meth_obj)
invisible(gc())

# ================================================================================
# STEP 5: QC OF THE SAVED OBJECT AGAINST THE RAW BEDMETHYL FILES
# ================================================================================

cat("================================================================================\n")
cat("Step 5: QC of the saved object against the raw bedMethyl files\n")
cat("================================================================================\n\n")

reloaded <- read_methylsense_qs(outfile)

object_checks <- c(
  is_methylRawList = methods::is(reloaded, "methylRawList"),
  sample_count = length(reloaded) == length(expected_ids),
  sample_ids_in_order = identical(as.character(getSampleID(reloaded)), expected_ids),
  treatment_in_order = identical(as.integer(getTreatment(reloaded)), expected_treatment)
)

verify_sample <- function(k) {
  tryCatch(
    {
      got <- as.data.table(getData(reloaded[[k]]))
      got[, chr := as.character(chr)]
      want_code <- build_qc$mod_code[k]

      raw <- read_bedmethyl(bed_paths[k])$bed
      raw <- raw[code == want_code]
      raw[, strand := methylkit_strand(strand)]
      if (isTRUE(opt$traditional)) {
        raw <- keep_first_chr_start(raw, opt$min_coverage)$bed
      } else {
        raw <- collapse_bed_sites(raw)$bed
        raw <- raw[(n_mod + n_canonical) >= opt$min_coverage]
      }
      raw <- raw[, .(chr, start, strand, n_mod, n_canonical)]

      m <- merge(got, raw, by = c("chr", "start", "strand"), all = TRUE)
      in_both <- m[!is.na(numCs) & !is.na(n_mod)]
      informative <- in_both[n_mod != n_canonical]

      data.table(
        sample = build_qc$sample[k],
        qc_sample_id_matches = identical(getSampleID(reloaded[[k]]), build_qc$sample[k]),
        qc_sites_in_object = nrow(got),
        qc_sites_expected_from_bed = nrow(raw),
        qc_sites_only_in_object = sum(is.na(m$n_mod)),
        qc_sites_only_in_bed = sum(is.na(m$numCs)),
        qc_end_is_start_plus_1 = mean(got$end == got$start + 1L),
        qc_coverage_is_numCs_plus_numTs = mean(got$coverage == got$numCs + got$numTs),
        qc_numCs_equals_Nmod = mean(in_both$numCs == in_both$n_mod),
        qc_numTs_equals_Ncanonical = mean(in_both$numTs == in_both$n_canonical),
        qc_informative_sites = nrow(informative),
        qc_informative_swapped = sum(informative$numCs == informative$n_canonical & informative$numTs == informative$n_mod),
        qc_mean_meth_bed_pct = {
          den <- sum(as.numeric(raw$n_mod + raw$n_canonical))
          if (den > 0) round(100 * sum(as.numeric(raw$n_mod)) / den, 3) else NA_real_
        },
        qc_mean_meth_object_pct = {
          den <- sum(as.numeric(got$coverage))
          if (den > 0) round(100 * sum(as.numeric(got$numCs)) / den, 3) else NA_real_
        },
        qc_error = NA_character_
      )
    },
    error = function(e) data.table(sample = build_qc$sample[k], qc_error = conditionMessage(e))
  )
}

verify_qc <- rbindlist(run_parallel(seq_len(nrow(build_qc)), verify_sample), fill = TRUE)
qc <- merge(build_qc, verify_qc, by = "sample", sort = FALSE)
qc <- qc[match(build_qc$sample, qc$sample)]

fail_reason <- function(r) {
  why <- character(0)
  if (!is.na(r$qc_error)) return(paste("QC error:", r$qc_error))
  if (!isTRUE(r$qc_sample_id_matches)) why <- c(why, "sample ID differs in saved object")
  if (!direction_confirmed(r$col11_agrees_with_Nmod, r$col11_closer_to_Nmod)) why <- c(why, "column 11 does not confirm column 12 as modified")
  if (r$count_identity < 0.999) why <- c(why, "count columns do not add up to Nvalid_cov")
  if (r$qc_sites_only_in_object > 0) why <- c(why, "sites in object that are not in the bedMethyl file")
  if (r$qc_sites_only_in_bed > 0) why <- c(why, "bedMethyl sites missing from object")
  if (r$qc_sites_in_object != r$sites_loaded) why <- c(why, "site count changed after saving")
  if (r$qc_end_is_start_plus_1 < 1) why <- c(why, "end is not start + 1")
  if (r$qc_coverage_is_numCs_plus_numTs < 1) why <- c(why, "coverage is not numCs + numTs")
  if (r$qc_numCs_equals_Nmod < 1) why <- c(why, "numCs differs from Nmod")
  if (r$qc_numTs_equals_Ncanonical < 1) why <- c(why, "numTs differs from Ncanonical")
  if (r$qc_informative_swapped > 0) why <- c(why, "methylated and unmethylated counts swapped")
  paste(why, collapse = "; ")
}

warn_reason <- function(r, cohort_median, cohort_mad) {
  why <- character(0)
  if (!r$mod_code %in% c("C", "m")) why <- c(why, sprintf("code %s is not 5mC (m) or combined C", r$mod_code))
  if (r$rows_other_code > 0) why <- c(why, sprintf("%s rows with other codes skipped", format(r$rows_other_code, big.mark = ",")))
  if (!is.na(r$n_other_mod_total) && r$n_other_mod_total > 0) why <- c(why, "Nother_mod > 0, so other modifications were left out of coverage")
  if (r$rows_unstranded > 0) why <- c(why, "unstranded rows were stored as + strand")
  if (r$rows_exact_duplicate > 0) {
    why <- c(why, sprintf(
      "%s identical extra %s dropped",
      format(r$rows_exact_duplicate, big.mark = ","),
      ifelse(r$rows_exact_duplicate == 1, "row was", "rows were")
    ))
  }
  if (r$sites_summed > 0) {
    why <- c(why, sprintf(
      "%s %s different counts and those counts were added",
      format(r$sites_summed, big.mark = ","),
      ifelse(r$sites_summed == 1, "site had", "sites had")
    ))
  }
  if (!is.na(r$traditional_dropped) && r$traditional_dropped > 0) {
    why <- c(why, sprintf(
      "%s later rows with the same chromosome and start were dropped. The first count was kept",
      format(r$traditional_dropped, big.mark = ",")
    ))
  }
  if (!is.na(r$nond_chrbase_dropped) && r$nond_chrbase_dropped > 0) {
    why <- c(why, sprintf(
      "%s later nonD rows with the same chromosome and start were dropped. The first row was kept",
      format(r$nond_chrbase_dropped, big.mark = ",")
    ))
  }
  if (is.na(r$col11_agrees_with_Nmod)) why <- c(why, "no informative sites to confirm direction")
  if (!is.na(r$col11_agrees_with_Nmod) && r$col11_agrees_with_Nmod < DIRECTION_EXACT_MIN) {
    why <- c(why, sprintf(
      "column 11 is not exactly Nmod / Nvalid_cov (exact on %.1f%% of sites) but is closer to Nmod on %.2f%%, as in simulated files",
      100 * r$col11_agrees_with_Nmod, 100 * r$col11_closer_to_Nmod
    ))
  }
  if (!is.na(cohort_mad) && cohort_mad > 0 && !is.na(r$mean_meth_pct) &&
    abs(r$mean_meth_pct - cohort_median) > 4 * cohort_mad) {
    why <- c(why, sprintf("mean methylation %.2f%% is far from the cohort median %.2f%%", r$mean_meth_pct, cohort_median))
  }
  paste(why, collapse = "; ")
}

cohort_median <- median(qc$mean_meth_pct)
cohort_mad <- if (nrow(qc) >= 5) mad(qc$mean_meth_pct) else NA_real_
qc[, qc_fail_reason := vapply(seq_len(.N), function(k) fail_reason(qc[k]), character(1))]
qc[, qc_warning := vapply(seq_len(.N), function(k) warn_reason(qc[k], cohort_median, cohort_mad), character(1))]
qc[, qc_status := fifelse(nzchar(qc_fail_reason), "FAIL", fifelse(nzchar(qc_warning), "WARN", "PASS"))]

all_pass <- all(object_checks) && !any(qc$qc_status == "FAIL")
qc_csv <- file.path(opt$output_dir, paste0(file_prefix, "_load_QC.csv"))
qc_txt <- file.path(opt$output_dir, paste0(file_prefix, "_load_QC_report.txt"))
fwrite(qc, qc_csv)

final_outfile <- outfile
if (!all_pass) {
  final_outfile <- sub("\\.qs2$", "_QC_FAILED.qs2", outfile)
  invisible(file.rename(outfile, final_outfile))
}

report <- c(
  "================================================================================",
  "MethylSense load QC report",
  "================================================================================",
  paste0("Generated:        ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Loader version:   ", SCRIPT_VERSION),
  paste0("Object:           ", basename(final_outfile)),
  paste0("Overall result:   ", ifelse(all_pass, "PASS", "FAIL")),
  "",
  "What is checked:",
  "  Column 11 of every bedMethyl file is modkit's own modified fraction. It must",
  "  equal column 12 / column 10, which confirms column 12 as the methylated count.",
  "  Files where column 11 is only closer to column 12 than to column 13 (simulated",
  "  data) pass with a warning.",
  "  Repeated rows for the same chromosome, start and strand are collapsed first.",
  "  Identical copies are kept once, and rows with different counts are added.",
  "  The saved .qs2 is read back from disk. Every site is matched to that collapsed",
  "  file by chromosome, start and strand. numCs must equal the collapsed column 12",
  "  and numTs must equal the collapsed column 13 at every site.",
  "",
  "Object checks:",
  paste0("  ", names(object_checks), ": ", ifelse(object_checks, "OK", "FAILED")),
  "",
  "Samples:",
  sprintf(
    "  %-24s %-4s %12s sites  meth %6.2f%%  col11 = Nmod on %7.3f%%  swapped %d",
    qc$sample, qc$qc_status, format(qc$sites_loaded, big.mark = ","),
    qc$mean_meth_pct, 100 * qc$col11_agrees_with_Nmod, as.integer(qc$qc_informative_swapped)
  ),
  ""
)
if (any(qc$qc_status == "FAIL")) {
  report <- c(report, "Failures:", paste0("  ", qc$sample[qc$qc_status == "FAIL"], ": ", qc$qc_fail_reason[qc$qc_status == "FAIL"]), "")
}
if (any(nzchar(qc$qc_warning))) {
  report <- c(report, "Warnings:", paste0("  ", qc$sample[nzchar(qc$qc_warning)], ": ", qc$qc_warning[nzchar(qc$qc_warning)]), "")
}
if (length(missing_samples) > 0) {
  report <- c(report, "Samples skipped because the bedMethyl file was missing:", paste0("  ", missing_samples), "")
}
report <- c(report, "================================================================================")
writeLines(report, qc_txt)

cat(paste(report[-(1:3)], collapse = "\n"), "\n\n")
cat("QC table:          ", basename(qc_csv), "\n")
cat("QC report:         ", basename(qc_txt), "\n\n")

if (!all_pass) {
  stop("Load QC failed. The object was renamed to ", basename(final_outfile),
    ". See ", basename(qc_txt), ".",
    call. = FALSE
  )
}

summary_file <- file.path(opt$output_dir, paste0(file_prefix, "_preprocessing_summary.txt"))
writeLines(c(
  "================================================================================",
  "MethylSense data preprocessing summary",
  "================================================================================",
  paste0("Generated:              ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Version:                ", SCRIPT_VERSION),
  paste0("Species:                ", opt$species),
  paste0("Assembly:               ", assembly_name),
  paste0("Total samples:          ", length(expected_ids)),
  paste0("Treatment groups:       ", treatment_text),
  paste0("Min CpG sites:          ", format(min(sample_sizes), big.mark = ",")),
  paste0("Max CpG sites:          ", format(max(sample_sizes), big.mark = ",")),
  paste0("Mean CpG sites:         ", format(round(mean(sample_sizes)), big.mark = ",")),
  paste0("Median CpG sites:       ", format(round(median(sample_sizes)), big.mark = ",")),
  paste0("Min coverage filter:    ", opt$min_coverage),
  paste0("Modification code:      ", paste(unique(qc$mod_code), collapse = ", ")),
  paste0(
    "MAPQ filter:            ",
    if (!is.null(opt$mapq_column) && !is.null(opt$min_mapq)) {
      paste0(
        opt$mapq_column, " >= ", opt$min_mapq,
        " (", n_samples_pre_mapq, " -> ", length(expected_ids), " samples)"
      )
    } else {
      "none"
    }
  ),
  paste0("Load QC:                ", ifelse(any(qc$qc_status == "WARN"), "PASS with warnings", "PASS")),
  paste0("Output file:            ", basename(final_outfile)),
  paste0("QC report:              ", basename(qc_txt)),
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
cat("Preprocessing complete\n")
cat("================================================================================\n\n")
cat("Status:            SUCCESS\n")
cat("Output directory:  ", opt$output_dir, "\n")
cat("Next step:         Run MethylSense analysis pipeline\n")
cat("\n")
cat("================================================================================\n\n")
