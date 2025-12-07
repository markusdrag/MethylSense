#!/usr/bin/env Rscript

# ================================================================================
# MethylSense: DMR Detection and Machine Learning Classification Pipeline
# ================================================================================
#
# Script Name: MethylSense_analysis.R
# Version: 5.13.6 (Public Release)
# Date: 2025-12-03
# GitHub: https://github.com/markusdrag/MethylSense
# Authors: Markus Hodal Drag, Christina Hvilsom, Louise Ladefoged Poulsen,
#          Henrik Elvang Jensen, Stamatios Alan Tahas, Christoph Leineweber,
#          Carolyn Cray, Mads Frost Bertelsen, Anders Miki Bojesen
#
# Description:
#   MethylSense is a comprehensive pipeline for detecting differentially
#   methylated regions (DMRs) from Nanopore whole-genome bisulfite sequencing
#   data and building machine learning classifiers to distinguish between
#   treatment groups. Originally developed for high-accuracy diagnostics of
#   avian Aspergillus fumigatus infection using host cell-free DNA methylation.
#
# Key Features:
#   - User-provided BED files for genomic window analysis
#   - Ten machine learning algorithms (RF, SVM, glmnet, nnet, kNN, LDA, NB, etc.)
#   - Cross-validation with 95% confidence intervals
#   - Nested cross-validation for unbiased performance estimation
#   - Predefined train/test splits for reproducibility
#   - Comprehensive visualization suite
#   - Caching system for faster re-runs
#   - Multi-class and binary classification support
#
# Citation:
#   If you use MethylSense, please cite:
#
#   Drag MH, Hvilsom C, Poulsen LL, Jensen HE, Tahas SA, Leineweber C,
#   Cray C, Bertelsen MF, Bojesen AM (2025)
#   New high accuracy diagnostics for avian Aspergillus fumigatus infection
#   using Nanopore methylation sequencing of host cell-free DNA and machine
#   learning prediction.
#   bioRxiv 2025.04.11.648151
#   https://doi.org/10.1101/2025.04.11.648151
#
# Contact: Markus Hodal Drag (first author)
# License: Academic Free License 3.0 (AFL-3.0)
# ================================================================================

# ================================================================================
# VERSION HISTORY
# ================================================================================
# v5.1.0 (2025-11-27) - PR-AUC BUG FIX (PUBLIC RELEASE)
#   FIXED: PR-AUC now correctly saved in analysis_summary.csv for multi-class designs
#   FIXED: Added pr_auc = 0 to error handler to prevent undefined variable errors
#   FIXED: Added pr_auc column to initial results_summary dataframe
#   IMPROVED: Multi-class PR-AUC values now consistently appear in all summary outputs
#
# v5.0.0 (2025-11-25) - AUTOMATED ANALYSIS REPORTS (PUBLIC RELEASE)
#   ADDED: generate_final_analysis_report() function for automated report generation
#   ADDED: Markdown and HTML analysis reports generated at pipeline completion
#   ADDED: Composite scoring system (40% sensitivity + 40% specificity + 20% accuracy)
#   ADDED: Best model recommendation based on nested CV metrics
#   ADDED: Top 10 models table with confidence intervals
#   IMPROVED: Professional UK English formatting in reports
#   IMPROVED: Complete model comparison with technical details
#
# v4.9.0 (2025-11-25) - PUBLIC RELEASE
#   NEW: Public version with traditional BED file analysis only
#   NOTE: Users must provide their own BED files for genomic window analysis
#
# v4.8.0 (2025-11-21) - FINAL FRIDAY EDITION
#   ADDED: PR curve visualization plots (binary and multi-class)
#   ADDED: plot_pr_curves() function with publication-quality formatting
#   ADDED: PR curve plots saved in multiple formats (png, pdf, svg, tiff)
#   IMPROVED: PR curves match ROC curve styling for consistency
#   IMPROVED: Multi-class PR curves with per-class legends and mean PR-AUC
#
# v4.7.0 (2025-11-21) - FRIDAY EDITION
#   FIXED: Critical sample name matching bug causing NA in treatment groups
#   FIXED: Positional indexing replaced with sample ID-based matching
#   ADDED: PR-AUC (Precision-Recall AUC) metric calculation for all models
#   ADDED: PR_AUC column in metrics.csv output
#   ADDED: PR_AUC statistics in summary CSV with mean and 95% CI
#   IMPROVED: Sample matching now uses getSampleID() from methylKit object
#   IMPROVED: Comprehensive sample matching diagnostics in logs
#
# v4.6.0 (2025-11-21)
#   FIXED: Group methylation summary returning NA for +2 group designs
#   FIXED: Heatmap annotation showing NA for +2 group designs
#   FIXED: Robust group name matching using as.character() conversion
#   ADDED: PRROC library for Precision-Recall curve computation
#   ADDED: Comprehensive debugging output for group matching diagnostics
#   IMPROVED: Placeholder entries for groups with no samples (prevents NA propagation)
#
# v4.5.0 (2025-11-21)
#   FIXED: Caret parallel processing limited by hardcoded 3-fold CV
#   ADDED: SLURM environment detection (SLURM_CPUS_PER_TASK, SLURM_JOB_CPUS_PER_NODE)
#   ADDED: Adaptive CV fold calculation based on available cores (5-10 folds)
#   IMPROVED: All trainControl instances now use adaptive cv_folds parameter
#   IMPROVED: Full CPU utilisation on SLURM clusters (8 cores = 8 parallel CV folds)
#   IMPROVED: Intelligent core allocation respects both --cores and SLURM limits
#
# v4.4.0 (2025-11-21)
#   FIXED: Parallel processing completely non-functional in v4.3.0
#   ADDED: doParallel and foreach libraries for proper parallel backend
#   ADDED: Parallel backend registration with makeCluster/registerDoParallel
#   FIXED: All allowParallel flags now conditional on --reproducible flag
#   ADDED: Proper cluster cleanup with stopCluster at script termination
#   ADDED: RAM usage warning for high core counts (>4 cores)
#   IMPROVED: 4-6x performance improvement in model training with parallel processing
#   IMPROVED: CPU utilisation now properly scales with --cores parameter
#
# v4.3.0 (2025-11-20)
#   FIXED: Methylation heatmap blank/white boxes in unsupervised clustering annotation
#   FIXED: Group methylation summary barplot missing groups in multi-group designs
#   IMPROVED: Consistent use of group_names_ordered variable throughout plotting functions
#   IMPROVED: Removed redundant group_names definition causing color mismatch in heatmap
#
# v4.2.0 (2025-11-17)
#   FIXED: DMR details duplication bug in 3-group designs
#   FIXED: negative_class vector handling - now properly collapses multiple classes
#   FIXED: neg_samples filtering to use %in% operator for multi-class support
#   IMPROVED: DMR details output now shows all negative classes in single row
#
# v4.1.0 (2025-11-14)
#   NEW: Sex chromosome exclusion flag (--exclude_sex_chroms)
#   NEW: Ability to exclude X, Y, Z, and W chromosomes from marker discovery
#   SAFETY: Prevents sex-linked markers in diagnostic panels
#   IMPROVED: Enhanced chromosome filtering with safety warnings
#
# v4.0.0 (2025-10-10) - MAJOR RELEASE
#   NEW: True nested cross-validation implementation (--nested_cv)
#   NEW: 95% confidence intervals for all CV metrics (AUC, sensitivity, specificity)
#   NEW: Complete train/test dataset export for validation
#   NEW: Cross-validation summary with publication-ready statistics
#   IMPROVED: Clean progress bars with ASCII characters
#   IMPROVED: Suppressed model training output for cleaner logs
#   IMPROVED: Fixed train/test split to prevent data leakage
#   IMPROVED: Enhanced plotting functions with better error handling
#   FIXED: Heatmap and bar plot generation
#   FIXED: Group methylation summary calculations
#
# v3.3.0 (2025-10-09)
#   NEW: Cross-validation repeats (--cv_repeats)
#   NEW: Organized output directory structure
#   NEW: Precomputed data caching system
#
# v3.2.0 (2025-10-08)
#   NEW: Predefined train/test splits from Excel
#   NEW: Treatment group mapping and filtering
#
# v3.1.0 (2025-10-07)
#   NEW: Multiple window sizes support
#   NEW: Region-specific caching
#
# v3.0.0 (2025-10-06)
#   Initial production release
#
# ================================================================================

# Script version for logging
SCRIPT_VERSION <- "5.13.6"
SCRIPT_DATE <- "2025-12-03"

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
cat("    High-Accuracy Epigenetic Diagnostics via Nanopore Methylation Sequencing   \n")
cat("================================================================================\n")
cat(paste("Version:", SCRIPT_VERSION, "|", "Release Date:", SCRIPT_DATE, "\n"))
cat(paste("Start Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
cat(paste("Working Directory:", getwd(), "\n"))
cat(paste("R Version:", R.version.string, "\n"))
cat("\n")
cat("Features:\n")
cat("  * 10 machine learning algorithms (RF, SVM, XGBoost, etc.)\n")
cat("  * Cross-validation with 95% confidence intervals\n")
cat("  * Nested CV for unbiased performance estimation\n")
cat("  * Custom group names and publication-ready colors\n")
cat("  * Comprehensive visualization suite\n")
cat("  * Species-agnostic methylation biomarker discovery\n")
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
# HELP AND DOCUMENTATION FUNCTIONS
# ================================================================================
show_help <- function() {
  cat("\n")
  cat("================================================================================\n")
  cat(paste0("                      METHYLSENSE ", SCRIPT_VERSION, " - HELP GUIDE                          \n"))
  cat("================================================================================\n")
  cat("\n")
  cat("DESCRIPTION:\n")
  cat("  MethylSense detects differentially methylated regions (DMRs) from Nanopore\n")
  cat("  bisulfite sequencing data and trains machine learning classifiers for\n")
  cat("  epigenetic diagnostics and biomarker discovery.\n")
  cat("\n")
  cat("USAGE:\n")
  cat("  Rscript MethylSense.R --qs_file <file> --output_dir <dir> [OPTIONS]\n")
  cat("\n")
  cat("================================================================================\n")
  cat("REQUIRED ARGUMENTS:\n")
  cat("================================================================================\n")
  cat("  --qs_file FILE        Path to methylRawList .qs file\n")
  cat("  --output_dir DIR      Output directory for results\n")
  cat("\n")
  cat("================================================================================\n")
  cat("CORE ANALYSIS OPTIONS:\n")
  cat("================================================================================\n")
  cat("  --bed_files FILE1,FILE2,...\n")
  cat("                        Comma-separated BED files for genomic windows\n")
  cat("  --models MODEL1,MODEL2,...\n")
  cat("                        ML models: rf,svm,glmnet,nnet,knn,lda,nb,ranger,xgboost\n")
  cat("                        Or use categories: all, fast, biomarker, interpretable\n")
  cat("                        [default: rf,svm,glmnet,nnet,knn,lda,nb]\n")
  cat("  --min_accuracy N      Minimum accuracy to save models [default: 0.7]\n")
  cat("\n")
  cat("================================================================================\n")
  cat("TREATMENT GROUP OPTIONS:\n")
  cat("================================================================================\n")
  cat("  --group_names NAME1,NAME2,...\n")
  cat("                        Custom group names (comma-separated, sets plot order)\n")
  cat("                        Example: --group_names Control,Infected\n")
  cat("                        Example: --group_names Healthy,Disease,Severe\n")
  cat("\n")
  cat("  --group_colors COLOR1,COLOR2,...\n")
  cat("                        Custom colors matching --group_names (comma-separated)\n")
  cat("                        Accepts hex codes (#RRGGBB) or R color names\n")
  cat("                        Example: --group_colors '#4575b4,#d73027'\n")
  cat("                        Example: --group_colors blue,red,orange\n")
  cat("                        Common palettes:\n")
  cat("                          Control vs Treatment: '#4575b4,#d73027' (blue,red)\n")
  cat("                          Healthy vs Disease:   '#1b9e77,#d95f02' (teal,orange)\n")
  cat("                          Three groups:         '#4575b4,#fee090,#d73027'\n")
  cat("\n")
  cat("  --positive_class NAME Positive class for binary ROC curves\n")
  cat("                        Must match one of --group_names\n")
  cat("                        Example: --positive_class Infected\n")
  cat("\n")
  cat("  --treatment_mapping CODE=NAME,...\n")
  cat("                        Map numeric treatment codes to group names\n")
  cat("                        Example: --treatment_mapping 0=Control,1=Infected\n")
  cat("\n")
  cat("  --min_group_size N    Minimum samples per group [default: 4]\n")
  cat("  --allow_small_groups  Allow groups with fewer samples (risky)\n")
  cat("\n")
  cat("================================================================================\n")
  cat("VALIDATION OPTIONS:\n")
  cat("================================================================================\n")
  cat("  --cv_repeats N        Number of CV repeats (0=disabled, 5-10 recommended)\n")
  cat("                        Always runs regular CV for quick validation\n")
  cat("                        Provides mean +/- 95% CI [default: 0]\n")
  cat("\n")
  cat("  --nested_cv           ADDITIONALLY run nested CV (gold standard)\n")
  cat("                        Requires --cv_repeats > 0\n")
  cat("                        Runs BOTH regular and nested CV\n")
  cat("                        Nested CV gives unbiased performance estimates\n")
  cat("                        Takes 5-10x longer than regular CV\n")
  cat("                        USE FOR PUBLICATION: Report nested CV results\n")
  cat("\n")
  cat("  Cross-Validation Types Explained:\n")
  cat("    Regular CV:  Fast, good for development, slightly optimistic\n")
  cat("    Nested CV:   Slow, gold standard, unbiased (required for papers)\n")
  cat("    Both saved:  Compare to estimate optimism bias\n")
  cat("\n")
  cat("================================================================================\n")
  cat("DMR FILTERING OPTIONS:\n")
  cat("================================================================================\n")
  cat("  --qvalue_threshold N  Q-value threshold for DMR significance [default: 0.05]\n")
  cat("  --min_coverage N      Minimum read coverage per CpG [default: 20]\n")
  cat("  --hyper_threshold N   Min % for hypermethylated DMRs [default: 5.0]\n")
  cat("  --hypo_threshold N    Max % for hypomethylated DMRs [default: -2.0]\n")
  cat("  --disable_meth_filter Disable methylation filtering (q-value only)\n")
  cat("  --min_group_size N    Min samples per group [default: 4]\n")
  cat("  --exclude_sex_chroms  Exclude sex chromosomes (X, Y, Z, W) from markers\n")
  cat("                        Recommended for safety to avoid sex-linked markers\n")
  cat("  --filter_standard_chroms\n")
  cat("                        Filter to standard chromosomes only (1-99 + sex)\n")
  cat("\n")
  cat("================================================================================\n")
  cat("FEATURE SELECTION OPTIONS:\n")
  cat("================================================================================\n")
  cat("  --max_markers N       Maximum DMRs as features [default: 100]\n")
  cat("  --min_markers N       Minimum DMRs to start training [default: 8]\n")
  cat("  --marker_step N       Step size for marker iteration [default: 1]\n")
  cat("\n")
  cat("================================================================================\n")
  cat("OUTPUT OPTIONS:\n")
  cat("================================================================================\n")
  cat("  --create_plots        Generate visualizations [default: TRUE]\n")
  cat("  --no_plots            Disable plot generation\n")
  cat("  --plot_format FORMAT  Output format: png, jpeg, svg, or 'all'\n")
  cat("                        Multiple formats: 'png,svg' [default: png]\n")
  cat("                        SVG recommended for Nature journals (vector graphics)\n")
  cat("  --plot_dpi N          Plot resolution in DPI [default: 300]\n")
  cat("                        Recommended: 300 (standard), 600 (high-res), 150 (draft)\n")
  cat("\n")
  cat("================================================================================\n")
  cat("OTHER OPTIONS:\n")
  cat("================================================================================\n")
  cat("  --reproducible        Fixed seeds for reproducibility [default: FALSE]\n")
  cat("  --force_ml            Force re-run even if cached results exist\n")
  cat("  --verbose             Enable verbose output\n")
  cat("  -h, --help            Show this help message\n")
  cat("  --show_output         Show detailed output directory structure\n")
  cat("  --version             Show version information\n")
  cat("  --list-models         Show available ML models and categories\n")
  cat("\n")
  cat("================================================================================\n")
  cat("QUICK START EXAMPLES:\n")
  cat("================================================================================\n")
  cat("\n")
  cat("1. Basic analysis:\n")
  cat("   Rscript MethylSense.R --qs_file data.qs --output_dir results/\n")
  cat("\n")
  cat("2. With custom group names and colors (publication-ready):\n")
  cat("   Rscript MethylSense.R --qs_file data.qs --output_dir results/ \\\n")
  cat("     --group_names Control,Infected \\\n")
  cat("     --group_colors '#4575b4,#d73027' \\\n")
  cat("     --positive_class Infected\n")
  cat("\n")
  cat("3. With 5-fold cross-validation:\n")
  cat("   Rscript MethylSense.R --qs_file data.qs --output_dir results/ \\\n")
  cat("     --cv_repeats 5 --group_names Control,Disease\n")
  cat("\n")
  cat("4. Nested CV for publication (gold standard):\n")
  cat("   Rscript MethylSense.R --qs_file data.qs --output_dir results/ \\\n")
  cat("     --cv_repeats 5 --nested_cv \\\n")
  cat("     --group_names Control,Infected \\\n")
  cat("     --group_colors blue,red\n")
  cat("\n")
  cat("5. With predefined train/test split:\n")
  cat("   Rscript MethylSense.R --qs_file data.qs --output_dir results/ \\\n")
  cat("     --train_test_file splits.xlsx --cv_repeats 5\n")
  cat("\n")
  cat("6. Three-group comparison with custom colors:\n")
  cat("   Rscript MethylSense.R --qs_file data.qs --output_dir results/ \\\n")
  cat("     --group_names Control,Mild,Severe \\\n")
  cat("     --group_colors '#4575b4,#fee090,#d73027'\n")
  cat("\n")
  cat("7. Fast biomarker models only:\n")
  cat("   Rscript MethylSense.R --qs_file data.qs --output_dir results/ \\\n")
  cat("     --models biomarker --cv_repeats 5\n")
  cat("\n")
  cat("For detailed output structure: Rscript MethylSense.R --show_output\n")
  cat("For available models: Rscript MethylSense.R --list-models\n")
  cat("\n")
  cat("8. High-resolution PNG for publication:\n")
  cat("   Rscript MethylSense.R --qs_file data.qs --output_dir results/ \\\n")
  cat("     --plot_format png --plot_dpi 600\n")
  cat("\n")
  cat("9. Vector graphics (SVG) for Nature journals:\n")
  cat("   Rscript MethylSense.R --qs_file data.qs --output_dir results/ \\\n")
  cat("     --plot_format svg --plot_dpi 300\n")
  cat("\n")
  cat("10. All formats (PNG + JPEG + SVG):\n")
  cat("    Rscript MethylSense.R --qs_file data.qs --output_dir results/ \\\n")
  cat("      --plot_format all --plot_dpi 300\n")
  cat("\n")
  cat("11. Multiple specific formats:\n")
  cat("    Rscript MethylSense.R --qs_file data.qs --output_dir results/ \\\n")
  cat("      --plot_format 'png,svg' --plot_dpi 600\n")
  cat("\n")
  cat("12. Excluding sex chromosomes (recommended for safety):\n")
  cat("    Rscript MethylSense.R --qs_file data.qs --output_dir results/ \\\n")
  cat("      --exclude_sex_chroms --group_names Control,Disease\n")
  cat("\n")
  cat("================================================================================\n")
  cat("COLOR RECOMMENDATIONS FOR PUBLICATIONS:\n")
  cat("================================================================================\n")
  cat("  Two groups (colorblind-safe):\n")
  cat("    Control vs Disease:  --group_colors '#4575b4,#d73027' (blue/red)\n")
  cat("    Healthy vs Infected: --group_colors '#1b9e77,#d95f02' (teal/orange)\n")
  cat("    Wild-type vs Mutant: --group_colors '#7570b3,#e7298a' (purple/pink)\n")
  cat("\n")
  cat("  Three groups (sequential):\n")
  cat("    Low/Med/High:        --group_colors '#4575b4,#fee090,#d73027'\n")
  cat("    Control/Mild/Severe: --group_colors '#1b9e77,#7570b3,#d95f02'\n")
  cat("\n")
  cat("  Color names (simple):\n")
  cat("    --group_colors blue,red\n")
  cat("    --group_colors steelblue,coral\n")
  cat("    --group_colors darkgreen,orange\n")
  cat("\n")
  cat("================================================================================\n")
  cat("CITATION:\n")
  cat("================================================================================\n")
  cat("  Drag MH, Hvilsom C, Poulsen LL, Jensen HE, Tahas SA, Leineweber C,\n")
  cat("  Cray C, Bertelsen MF, Bojesen AM (2025)\n")
  cat("  New high accuracy diagnostics for avian Aspergillus fumigatus infection\n")
  cat("  using Nanopore methylation sequencing of host cell-free DNA and machine\n")
  cat("  learning prediction. bioRxiv 2025.04.11.648151\n")
  cat("  https://doi.org/10.1101/2025.04.11.648151\n")
  cat("================================================================================\n\n")
}

show_output_structure <- function() {
  cat("================================================================================\n")
  cat("OUTPUT DIRECTORY STRUCTURE:\n")
  cat("================================================================================\n")
  cat("\n")
  cat("output_dir/\n")
  cat("|\n")
  cat("+-- [region_name]/                    # One folder per genomic window size\n")
  cat("|   |\n")
  cat("|   +-- [model_name]/                 # One folder per ML model (rf, svm, etc.)\n")
  cat("|       |\n")
  cat("|       +-- markers_[N]/              # Results for N markers\n")
  cat("|           |\n")
  cat("|           +-- model.rds             # Trained model object (for predictions)\n")
  cat("|           +-- predictions.csv       # Test set predictions\n")
  cat("|           +-- metrics.csv           # Performance metrics (Acc, AUC, Sens, Spec, F1)\n")
  cat("|           |\n")
  cat("|           +-- cross_validation/     # CV results (if --cv_repeats > 0)\n")
  cat("|           |   |\n")
  cat("|           |   +-- cv_summary_[model]_[N]markers.csv\n")
  cat("|           |   |   # Mean +/- 95% CI for all metrics\n")
  cat("|           |   |\n")
  cat("|           |   +-- cv_detailed_[model]_[N]markers.csv\n")
  cat("|           |   |   # Individual results for each CV repeat\n")
  cat("|           |   |\n")
  cat("|           |   +-- cv_predictions_[model]_[N]markers.csv\n")
  cat("|           |   |   # Predictions from all CV repeats\n")
  cat("|           |   |\n")
  cat("|           |   +-- nested_cv_summary_[model]_[N]markers.csv\n")
  cat("|           |   |   # Unbiased estimates (if --nested_cv enabled)\n")
  cat("|           |   |\n")
  cat("|           |   +-- nested_cv_hyperparameters_[model]_[N]markers.csv\n")
  cat("|           |       # Best hyperparameters per outer fold\n")
  cat("|           |\n")
  cat("|           +-- dmr_barplot_[N]markers.png\n")
  cat("|           |   # Bar plot of top DMRs with significance\n")
  cat("|           |\n")
  cat("|           +-- methylation_heatmap_[N]markers.png\n")
  cat("|           |   # Clustered heatmap of methylation patterns\n")
  cat("|           |\n")
  cat("|           +-- group_methylation_summary_[N]markers.png\n")
  cat("|           |   # Mean methylation per group with error bars\n")
  cat("|           |\n")
  cat("|           +-- group_methylation_stats_[N]markers.csv\n")
  cat("|           |   # Statistics for group summary plot\n")
  cat("|           |\n")
  cat("|           +-- roc_curves_multiclass.png\n")
  cat("|           |   # ROC curves (one-vs-rest for multi-class)\n")
  cat("|           |\n")
  cat("|           +-- confusion_matrix_plot.png\n")
  cat("|           |   # Prediction vs actual scatter plot\n")
  cat("|           |\n")
  cat("|           +-- metrics_plot.png\n")
  cat("|               # Bar chart of all performance metrics\n")
  cat("|\n")
  cat("+-- train_test_splits/\n")
  cat("|   |\n")
  cat("|   +-- TRAIN_dataset_[region]_[N]markers.csv\n")
  cat("|   |   # Complete training dataset with all features\n")
  cat("|   |\n")
  cat("|   +-- TEST_dataset_[region]_[N]markers.csv\n")
  cat("|   |   # Complete test dataset with all features\n")
  cat("|   |\n")
  cat("|   +-- VALIDATION_[region]_[N]markers.txt\n")
  cat("|   |   # Text summary of split with sample IDs\n")
  cat("|   |\n")
  cat("|   +-- split_[region]_[N]markers.csv\n")
  cat("|       # Simple train/test assignment per sample\n")
  cat("|\n")
  cat("+-- precomputed_data/\n")
  cat("|   |\n")
  cat("|   +-- precomputed_methylation_matrix_[region].rds\n")
  cat("|   |   # Cached methylation matrix for significant DMRs\n")
  cat("|   |\n")
  cat("|   +-- dmr_summary_[region].csv\n")
  cat("|       # Summary of all detected DMRs for this region\n")
  cat("|\n")
  cat("+-- complete_dmr_datasets/\n")
  cat("|   |\n")
  cat("|   +-- ALL_DMRs_complete_dataset_[region].csv\n")
  cat("|   |   # ALL detected DMRs (not filtered by significance)\n")
  cat("|   |\n")
  cat("|   +-- united_meth_object_[region]_[timestamp].qs\n")
  cat("|   |   # United methylation object (can be large)\n")
  cat("|   |\n")
  cat("|   +-- united_meth_metadata_[region]_[timestamp].rds\n")
  cat("|       # Metadata for united object\n")
  cat("|\n")
  cat("+-- clustering_analysis/\n")
  cat("|   |\n")
  cat("|   +-- [clustering plots and results]\n")
  cat("|\n")
  cat("+-- ml_training_logs/\n")
  cat("|   |\n")
  cat("|   +-- training_log_[timestamp].txt\n")
  cat("|       # Detailed log of ML training process\n")
  cat("|\n")
  cat("+-- unmatched_samples/\n")
  cat("    |\n")
  cat("    +-- unmatched_samples_[region].csv\n")
  cat("        # Samples that couldn't be matched (if using --train_test_file)\n")
  cat("\n")
  cat("KEY FILES FOR PUBLICATION:\n")
  cat("  - cv_summary_*.csv        : Mean performance +/- 95% CI\n")
  cat("  - nested_cv_summary_*.csv : Unbiased performance estimates\n")
  cat("  - metrics.csv             : Single train/test performance\n")
  cat("  - TRAIN/TEST_dataset_*.csv: Complete datasets for validation\n")
  cat("  - All .png files          : Figures for manuscript\n")
  cat("\n")
}

# Check if --help, --show_output, or --version was requested BEFORE full parsing
args <- commandArgs(trailingOnly = TRUE)

if ("--help" %in% args || "-h" %in% args) {
  suppressPackageStartupMessages({
    library(optparse)
  })
  show_help()
  quit(save = "no", status = 0)
}

if ("--show_output" %in% args) {
  suppressPackageStartupMessages({
    library(optparse)
  })
  show_output_structure()
  quit(save = "no", status = 0)
}

if ("--version" %in% args) {
  cat("\n")
  cat("================================================================================\n")
  cat("MethylSense: Interpretable ML for Epigenetic Diagnostics\n")
  cat("================================================================================\n")
  cat(paste("Version:", SCRIPT_VERSION, "\n"))
  cat(paste("Release Date:", SCRIPT_DATE, "\n"))
  cat("\n")
  cat("Citation:\n")
  cat("  Drag MH, Hvilsom C, Poulsen LL, Jensen HE, Tahas SA, Leineweber C,\n")
  cat("  Cray C, Bertelsen MF, Bojesen AM (2025)\n")
  cat("  bioRxiv 2025.04.11.648151 | doi: 10.1101/2025.04.11.648151\n")
  cat("\n")
  cat("For full help: Rscript MethylSense.R --help\n")
  cat("For output structure: Rscript MethylSense.R --show_output\n")
  cat("================================================================================\n")
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({
  library(jsonlite) # ADD THIS for nested CV hyperparameter storage
  library(optparse)
  library(methylKit)
  library(GenomicRanges)
  library(caret)
  library(randomForest)
  library(pROC)
  library(PRROC)
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(reshape2)
  library(qs)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(genomation)
  library(regioneR)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  library(xgboost)
  library(naivebayes)
  library(ranger)
  library(MASS)
  library(e1071)
  library(tree)
  library(rpart)
  library(rpart.plot)
  library(data.table)
  library(factoextra)
  library(cluster)
  library(corrplot)
  library(readxl)
})

cat("[OK] Libraries loaded successfully\n")

# NOTE: Reproducibility settings moved to after option parsing
# so they can be controlled by --reproducible flag

# ---------------------------
# INTEGRATED MODEL SUMMARY GENERATION
# ---------------------------

# Platform-specific null device
nullfile <- function() {
  if (.Platform$OS.type == "windows") {
    return("NUL")
  } else {
    return("/dev/null")
  }
}


# Function to generate model summary (integrated from your summary script)
generate_model_summary <- function(model_rds_path, output_dir) {
  tryCatch(
    {
      summary_file <- file.path(output_dir, "model_summary.txt")

      # Load model
      m <- readRDS(model_rds_path)

      # Helper functions
      safe <- function(expr, default = NULL) tryCatch(expr, error = function(e) default)
      line <- function(char = "-", n = 80) paste(rep(char, n), collapse = "")

      # Extract information
      method <- safe(m$method)
      metric <- safe(m$metric)
      maximize <- safe(m$maximize)
      outcome <- safe(names(m$trainingData)[ncol(m$trainingData)])
      n_obs <- safe(nrow(m$trainingData))
      n_feat <- safe(ncol(m$trainingData) - 1)

      tc <- safe(m$control)
      best_tune <- safe(m$bestTune)
      results_tbl <- safe(m$results)
      train_perf <- safe(caret::getTrainPerf(m))

      # Variable importance
      vi_tbl <- safe({
        vi <- caret::varImp(m)
        if (!is.null(vi$importance)) {
          imp <- vi$importance
          imp$Feature <- rownames(imp)
          rownames(imp) <- NULL
          ord <- order(-apply(
            imp[, setdiff(names(imp), "Feature"), drop = FALSE], 1,
            function(r) suppressWarnings(max(as.numeric(r), na.rm = TRUE))
          ))
          imp[ord, , drop = FALSE]
        } else {
          NULL
        }
      })

      # Coefficients extraction (handles glmnet properly)
      extract_coefs_df <- function(m) {
        cf <- safe(coef(m$finalModel))
        if (inherits(m$finalModel, "glmnet")) {
          lam <- safe(m$bestTune$lambda)
          cf <- safe(if (!is.null(lam)) coef(m$finalModel, s = lam) else coef(m$finalModel))
        }
        if (inherits(cf, "dgCMatrix") || is.matrix(cf)) {
          mat <- as.matrix(cf)
          if (ncol(mat) > 1) {
            sel <- 1L
            lam <- safe(as.numeric(m$bestTune$lambda))
            if (!is.null(lam) && !is.null(colnames(mat))) {
              diffs <- suppressWarnings(abs(as.numeric(colnames(mat)) - lam))
              if (all(is.finite(diffs))) sel <- which.min(diffs)
            }
            mat <- mat[, sel, drop = FALSE]
          }
          df <- data.frame(Term = rownames(mat), Estimate = as.numeric(mat[, 1]), row.names = NULL)
          df[!is.na(df$Estimate) & df$Estimate != 0, , drop = FALSE]
        } else {
          NULL
        }
      }

      coefs_df <- extract_coefs_df(m)

      # Timing
      times_list <- safe({
        t <- m$times
        if (inherits(t, "proc_time")) {
          v <- as.numeric(t)
          nms <- c("user", "system", "elapsed", "user.child", "system.child")[seq_along(v)]
          names(v) <- nms
          as.list(v)
        } else {
          t
        }
      })

      # Write summary
      con <- file(summary_file, open = "wt")
      sink(con)
      on.exit(
        {
          sink()
          close(con)
        },
        add = TRUE
      )

      cat("Caret Model Summary\n", line(), "\n\n")

      # Overview
      cat(line("="), "\nOVERVIEW\n", line("="), "\n")
      cat(sprintf("File:      %s\n", model_rds_path))
      cat(sprintf("Method:    %s\n", method))
      cat(sprintf("Metric:    %s (maximize: %s)\n", metric, maximize))
      cat(sprintf("Outcome:   %s\n", outcome))
      cat(sprintf("Data:      %s rows, %s predictors\n", n_obs, n_feat))

      # Resampling
      cat("\n", line("="), "\nRESAMPLING (TRAINCONTROL)\n", line("="), "\n")
      cat(sprintf("Method:        %s\n", safe(tc$method)))
      cat(sprintf("Number:        %s\n", safe(tc$number)))
      cat(sprintf("Repeats:       %s\n", safe(tc$repeats)))
      cat(sprintf("Class Probs:   %s\n", safe(tc$classProbs)))
      cat(sprintf("Save Preds:    %s\n", safe(tc$savePredictions)))
      resampling_summaryFn <- safe(if (!is.null(tc$summaryFunction)) deparse(body(tc$summaryFunction))[1])
      if (!is.null(resampling_summaryFn)) {
        cat("Summary Fn:    ")
        cat(resampling_summaryFn, "\n")
      }

      # Best parameters
      cat("\n", line("="), "\nBEST HYPERPARAMETERS\n", line("="), "\n")
      if (!is.null(best_tune)) {
        print(best_tune)
      } else {
        cat("(none recorded)\n")
      }

      # Tuning results
      cat("\n", line("="), "\nTUNING RESULTS (HEAD)\n", line("="), "\n")
      if (!is.null(results_tbl)) {
        print(utils::head(results_tbl, 20))
      } else {
        cat("(none recorded)\n")
      }

      # Training performance
      cat("\n", line("="), "\nTRAINING PERFORMANCE\n", line("="), "\n")
      if (!is.null(train_perf)) {
        print(train_perf)
      } else {
        cat("(not available)\n")
      }

      # Variable importance
      cat("\n", line("="), "\nVARIABLE IMPORTANCE (TOP 30)\n", line("="), "\n")
      if (!is.null(vi_tbl) && nrow(vi_tbl) > 0) {
        print(utils::head(vi_tbl, 30), row.names = FALSE)
      } else {
        cat("(not available for this model)\n")
      }

      # Final model details
      cat("\n", line("="), "\nFINAL MODEL DETAILS\n", line("="), "\n")
      cat(sprintf("Class: %s\n", paste(class(m$finalModel), collapse = ", ")))

      cat("\nCoefficients (non-zero, selected lambda if glmnet)\n", line("-", 45), "\n")
      if (!is.null(coefs_df) && nrow(coefs_df) > 0) {
        ord <- order(-abs(coefs_df$Estimate))
        print(utils::head(coefs_df[ord, , drop = FALSE], 50), row.names = FALSE)
        if (nrow(coefs_df) > 50) {
          cat(sprintf("... %d more coefficients omitted ...\n", nrow(coefs_df) - 50))
        }
      } else {
        cat("(not available for this model type)\n")
      }

      # Preprocessing
      cat("\n", line("="), "\nPREPROCESSING\n", line("="), "\n")
      pp_txt <- safe(paste(capture.output(str(m$preProcess)), collapse = "\n"))
      cat(ifelse(is.null(pp_txt), "(none or not recorded)\n", paste0(pp_txt, "\n")))

      # Timing
      cat("\n", line("="), "\nTIMING\n", line("="), "\n")
      if (!is.null(times_list)) {
        print(times_list)
      } else {
        cat("(not recorded)\n")
      }

      cat("\n", line(), "\nEnd of report\n")

      return(summary_file)
    },
    error = function(e) {
      log_warning(paste("Failed to generate model summary:", e$message), "Summary generation")
      return(NULL)
    }
  )
}

# ---------------------------
# ENHANCED LOGGING SYSTEM
# ---------------------------

# Global logging variables
log_file_path <- NULL
log_buffer <- character(0)
log_start_time <- Sys.time()

# Initialize logging system
initialize_logging <- function(output_dir) {
  log_file_path <<- file.path(output_dir, "run_log.txt")

  # Create log header
  header_lines <- c(
    paste(rep("=", 80), collapse = ""),
    "DMR + MACHINE LEARNING ANALYSIS - COMPLETE RUN LOG",
    paste("Script Version:", SCRIPT_VERSION),
    paste("Start Time:", format(log_start_time, "%Y-%m-%d %H:%M:%S")),
    paste("Working Directory:", getwd()),
    paste("Output Directory:", output_dir),
    paste("R Version:", R.version.string),
    paste(rep("=", 80), collapse = ""),
    ""
  )

  # Write header immediately
  writeLines(header_lines, log_file_path)

  if (!is.null(log_file_path) && is.character(log_file_path)) {
    cat("[LOG] Logging initialized:", basename(log_file_path), "\n")
  } else {
    cat("[LOG] Logging initialized: run_log.txt\n")
  }
}

# Enhanced logging function that writes to both console and log file
log_msg <- function(msg, level = "INFO", write_immediately = TRUE) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  elapsed_time <- as.numeric(Sys.time() - log_start_time, units = "secs")
  elapsed_str <- sprintf("[+%.1fs]", elapsed_time)

  # Format the complete log entry
  log_entry <- paste(timestamp, elapsed_str, paste0("[", level, "]"), msg)

  # Always print to console (with original timestamp format for consistency)
  console_msg <- paste(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg)
  cat(console_msg, "\n")

  # Add to log file
  if (!is.null(log_file_path)) {
    if (write_immediately) {
      # Write immediately for important messages
      cat(log_entry, "\n", file = log_file_path, append = TRUE)
    } else {
      # Buffer for less critical messages (flush periodically)
      log_buffer <<- c(log_buffer, log_entry)

      # Flush buffer if it gets too large
      if (length(log_buffer) >= 50) {
        flush_log_buffer()
      }
    }
  }
}

# Flush buffered log entries
flush_log_buffer <- function() {
  if (!is.null(log_file_path) && length(log_buffer) > 0) {
    cat(paste(log_buffer, collapse = "\n"), "\n", file = log_file_path, append = TRUE)
    log_buffer <<- character(0)
  }
}

# Enhanced verbose and debug messaging
verbose_msg <- function(msg) {
  if (exists("opt") && !is.null(opt$verbose) && opt$verbose) {
    log_msg(paste("   [VERBOSE]", msg), level = "VERBOSE", write_immediately = FALSE)
  }
}

debug_msg <- function(msg) {
  if (exists("opt") && !is.null(opt$debug) && opt$debug) {
    log_msg(paste("   [DEBUG]", msg), level = "DEBUG", write_immediately = FALSE)
  }
}

log_file_operation <- function(operation, file_path, details = NULL) {
  file_info <- if (!is.null(details)) paste0(" (", details, ")") else ""
  safe_basename <- if (!is.null(file_path) && is.character(file_path)) basename(file_path) else "unknown_file"
  log_msg(paste("[DIR]", operation, ":", safe_basename, file_info), level = "FILE")
}

# Function to log analysis steps
log_analysis_step <- function(step_name, details = NULL) {
  detail_info <- if (!is.null(details)) paste0(" | ", details) else ""
  log_msg(paste("[ANALYSIS]", step_name, detail_info), level = "ANALYSIS")
}

# Function to log errors and warnings
log_error <- function(error_msg, context = NULL) {
  context_info <- if (!is.null(context)) paste0(" [Context: ", context, "]") else ""
  log_msg(paste("[ERROR]", error_msg, context_info), level = "ERROR")
}

log_warning <- function(warning_msg, context = NULL) {
  context_info <- if (!is.null(context)) paste0(" [Context: ", context, "]") else ""
  log_msg(paste("[WARN]", warning_msg, context_info), level = "WARNING")
}
# ================================================================================
# UNIVERSAL PLOT SAVING FUNCTION
# ================================================================================
save_plot_universal <- function(plot_obj, file_path, width = 10, height = 8,
                                formats = NULL, dpi = NULL) {
  # Get formats and DPI from global opt if not specified
  if (is.null(formats)) {
    if (exists("opt", envir = .GlobalEnv)) {
      opt_global <- get("opt", envir = .GlobalEnv)
      formats <- opt_global$plot_format
    } else {
      formats <- "png"
    }
  }

  if (is.null(dpi)) {
    if (exists("opt", envir = .GlobalEnv)) {
      opt_global <- get("opt", envir = .GlobalEnv)
      dpi <- opt_global$plot_dpi
    } else {
      dpi <- 300
    }
  }

  # Parse formats
  if (formats == "all") {
    formats <- c("png", "jpeg", "svg")
  } else {
    formats <- unlist(strsplit(formats, ","))
    formats <- trimws(formats)
  }

  saved_files <- list()
  base_path <- tools::file_path_sans_ext(file_path)

  for (format in formats) {
    format <- tolower(format)
    output_file <- paste0(base_path, ".", format)

    tryCatch(
      {
        if (format == "svg") {
          ggsave(output_file, plot_obj,
            width = width, height = height,
            device = "svg", dpi = dpi
          )
          verbose_msg(paste("  Saved SVG (vector):", basename(output_file)))
        } else if (format == "png") {
          ggsave(output_file, plot_obj,
            width = width, height = height,
            device = "png", dpi = dpi
          )
          verbose_msg(paste("  Saved PNG", paste0(dpi, "dpi:"), basename(output_file)))
        } else if (format == "jpeg" || format == "jpg") {
          output_file <- paste0(base_path, ".jpeg") # Normalize to .jpeg
          ggsave(output_file, plot_obj,
            width = width, height = height,
            device = "jpeg", dpi = dpi, quality = 95
          )
          verbose_msg(paste("  Saved JPEG", paste0(dpi, "dpi:"), basename(output_file)))
        } else {
          log_warning(paste("Unknown plot format:", format, "- skipping"), "Plotting")
          next
        }
        saved_files[[format]] <- output_file
      },
      error = function(e) {
        log_warning(paste("Failed to save plot as", format, ":", e$message), "Plotting")
      }
    )
  }

  return(saved_files)
}

# ---------------------------
# MISSING HELPER FUNCTIONS - ADD THESE
# ---------------------------

# Function to check for completed analyses and filter region list
check_completed_analyses <- function(region_data, output_dir) {
  log_msg("[DEBUG] Checking for previously completed analyses...")

  completed_regions <- character(0)
  remaining_regions <- list()

  for (region_name in names(region_data)) {
    region_info <- region_data[[region_name]]
    reg_label <- region_info$label

    # Define expected output files for this region
    precomputed_matrix_file <- file.path(output_dir, "precomputed_data", paste0("precomputed_methylation_matrix_", reg_label, ".rds"))
    dmr_summary_file <- file.path(output_dir, "precomputed_data", paste0("dmr_summary_", reg_label, ".csv"))

    # Check if both critical files exist
    matrix_exists <- file.exists(precomputed_matrix_file)
    summary_exists <- file.exists(dmr_summary_file)

    if (matrix_exists && summary_exists) {
      completed_regions <- c(completed_regions, reg_label)
      verbose_msg(paste("[OK] Found completed analysis:", reg_label))
    } else {
      remaining_regions[[region_name]] <- region_info
      verbose_msg(paste("[NEXT] No previous analysis found for:", reg_label, "- will process"))
    }
  }

  return(list(
    completed = completed_regions,
    remaining = remaining_regions,
    total = length(region_data),
    completed_count = length(completed_regions),
    remaining_count = length(remaining_regions)
  ))
}

# Function to load and merge completed results
load_completed_results <- function(completed_regions, output_dir) {
  if (length(completed_regions) == 0) {
    return(data.frame())
  }

  log_msg("[DIR] Loading results from completed analyses...")
  all_results <- data.frame()

  # Look for existing results summary files
  possible_files <- c(
    file.path(output_dir, "analysis_summary.csv")
  )

  for (results_file in possible_files) {
    if (file.exists(results_file)) {
      tryCatch(
        {
          existing_results <- read.csv(results_file)
          region_results <- existing_results[existing_results$region %in% completed_regions, ]

          if (nrow(region_results) > 0) {
            all_results <- rbind(all_results, region_results)
          }
        },
        error = function(e) {
          verbose_msg(paste("[WARN] Could not load results from:", basename(results_file)))
        }
      )
    }
  }

  return(all_results)
}

# ---------------------------
# ORGANIZED OUTPUT DIRECTORY STRUCTURE
# ---------------------------

# Function to create organized directory structure
create_output_structure <- function(base_output_dir) {
  directories <- list(
    ml_logs = file.path(base_output_dir, "ml_training_logs"),
    train_test_splits = file.path(base_output_dir, "train_test_splits"),
    unmatched_samples = file.path(base_output_dir, "unmatched_samples"),
    precomputed_data = file.path(base_output_dir, "precomputed_data"),
    complete_dmr_datasets = file.path(base_output_dir, "complete_dmr_datasets"),
    cv_analysis_plots = file.path(base_output_dir, "cv_analysis_plots")
    # REMOVED: clustering_analysis - only create if clustering is enabled
  )

  log_msg("[DIR] Creating organized output directory structure...")

  for (dir_name in names(directories)) {
    dir_path <- directories[[dir_name]]
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    log_file_operation("Created directory", dir_path)
  }

  # Only create clustering directory if enabled
  if (!is.null(opt$perform_clustering) && opt$perform_clustering) {
    directories$clustering_analysis <- file.path(base_output_dir, "clustering_analysis")
    dir.create(directories$clustering_analysis, recursive = TRUE, showWarnings = FALSE)
    log_file_operation("Created directory", directories$clustering_analysis)
  }

  return(directories)
}

# ==============================================================================
# FINAL ANALYSIS REPORT GENERATION
# ==============================================================================

generate_final_analysis_report <- function(summary_csv_path, output_dir) {
  if (!file.exists(summary_csv_path)) {
    log_msg("[WARN] Cannot generate report: analysis_summary.csv not found")
    return(NULL)
  }

  log_msg("[REPORT] Generating final analysis report...")

  # Load data
  results <- read.csv(summary_csv_path, stringsAsFactors = FALSE)

  if (nrow(results) == 0) {
    log_msg("[WARN] Cannot generate report: no results in analysis_summary.csv")
    return(NULL)
  }

  # ===========================================================================
  # ANALYSIS AND MODEL SELECTION
  # ===========================================================================

  # Filter to models with nested CV results
  nested_results <- results[!is.na(results$nested_cv_mean_accuracy), ]

  if (nrow(nested_results) == 0) {
    log_msg("[WARN] No nested CV results found - using regular CV for report")
    nested_results <- results
    use_nested <- FALSE
  } else {
    use_nested <- TRUE
  }

  # Calculate composite score for ranking
  # Priority: Sensitivity (40%) + Specificity (40%) + Accuracy (20%)
  if (use_nested) {
    nested_results$composite_score <- (
      0.4 * nested_results$nested_cv_mean_sensitivity +
        0.4 * nested_results$nested_cv_mean_specificity +
        0.2 * nested_results$nested_cv_mean_accuracy
    )
    nested_results <- nested_results[order(-nested_results$composite_score), ]
  } else {
    nested_results$composite_score <- nested_results$accuracy
    nested_results <- nested_results[order(-nested_results$composite_score), ]
  }

  # Get top 10 and best model
  top10 <- head(nested_results, 10)
  best <- nested_results[1, ]

  # ===========================================================================
  # GENERATE MARKDOWN REPORT
  # ===========================================================================

  md_file <- file.path(output_dir, "ANALYSIS_REPORT.md")
  md <- c()

  # Header
  md <- c(md, "# MethylSense Analysis Report")
  md <- c(md, "")
  md <- c(md, paste0("**Generated:** ", format(Sys.time(), "%d %B %Y at %H:%M %Z")))
  md <- c(md, "")
  md <- c(md, "---")
  md <- c(md, "")

  # Executive Summary
  md <- c(md, "## Executive Summary")
  md <- c(md, "")
  md <- c(md, sprintf("- **Total models tested:** %d", nrow(results)))
  md <- c(md, sprintf("- **Genomic region:** %s", unique(results$region)[1]))
  if (use_nested) {
    md <- c(md, sprintf("- **Models with nested CV:** %d", sum(!is.na(results$nested_cv_mean_accuracy))))
  }
  md <- c(md, sprintf("- **Best performing model:** %s with %d markers", best$model, best$markers))
  md <- c(md, "")

  # Recommended Model
  md <- c(md, "## Recommended Model")
  md <- c(md, "")
  md <- c(md, sprintf("### %s (%d markers)", best$model, best$markers))
  md <- c(md, "")

  if (use_nested) {
    md <- c(md, "**Nested Cross-Validation Performance (Unbiased Estimates):**")
    md <- c(md, "")
    md <- c(md, sprintf(
      "- **Sensitivity:** %.1f%% (95%% CI: %.1f%%–%.1f%%)",
      best$nested_cv_mean_sensitivity * 100,
      (best$nested_cv_mean_sensitivity - 1.96 * best$nested_cv_sd_accuracy / sqrt(10)) * 100,
      (best$nested_cv_mean_sensitivity + 1.96 * best$nested_cv_sd_accuracy / sqrt(10)) * 100
    ))
    md <- c(md, sprintf(
      "- **Specificity:** %.1f%% (95%% CI: %.1f%%–%.1f%%)",
      best$nested_cv_mean_specificity * 100,
      (best$nested_cv_mean_specificity - 1.96 * best$nested_cv_sd_accuracy / sqrt(10)) * 100,
      (best$nested_cv_mean_specificity + 1.96 * best$nested_cv_sd_accuracy / sqrt(10)) * 100
    ))
    md <- c(md, sprintf(
      "- **Accuracy:** %.1f%% (95%% CI: %.1f%%–%.1f%%)",
      best$nested_cv_mean_accuracy * 100,
      best$nested_cv_95ci_lower * 100,
      best$nested_cv_95ci_upper * 100
    ))
    md <- c(md, sprintf("- **AUC:** %.3f", best$nested_cv_mean_auc))
  } else {
    md <- c(md, "**Cross-Validation Performance:**")
    md <- c(md, "")
    md <- c(md, sprintf("- **Accuracy:** %.1f%%", best$accuracy * 100))
    md <- c(md, sprintf("- **AUC:** %.3f", best$auc))
    md <- c(md, sprintf("- **Sensitivity:** %.1f%%", best$mean_sensitivity * 100))
    md <- c(md, sprintf("- **Specificity:** %.1f%%", best$mean_specificity * 100))
  }
  md <- c(md, "")

  # Justification
  md <- c(md, "**Justification:**")
  md <- c(md, "")
  md <- c(md, sprintf(
    "This model was selected based on optimal balance of sensitivity (%.1f%%) and specificity (%.1f%%), ",
    if (use_nested) best$nested_cv_mean_sensitivity * 100 else best$mean_sensitivity * 100,
    if (use_nested) best$nested_cv_mean_specificity * 100 else best$mean_specificity * 100
  ))
  md <- c(md, sprintf(
    "combined with strong overall accuracy (%.1f%%). ",
    if (use_nested) best$nested_cv_mean_accuracy * 100 else best$accuracy * 100
  ))

  if (use_nested && !is.na(best$cv_optimism_bias)) {
    if (abs(best$cv_optimism_bias) < 0.02) {
      md <- c(md, "The low optimism bias (<2%) indicates robust generalisation performance.")
    } else if (abs(best$cv_optimism_bias) < 0.05) {
      md <- c(md, "Moderate optimism bias detected; nested CV estimates are recommended for publication.")
    } else {
      md <- c(md, "Significant optimism bias (>5%) detected; nested CV provides unbiased performance estimates.")
    }
  }

  md <- c(md, "")
  md <- c(md, "---")
  md <- c(md, "")

  # Top 10 Models
  md <- c(md, "## Top 10 Models by Performance")
  md <- c(md, "")

  if (use_nested) {
    md <- c(md, "| Rank | Model | Markers | Sensitivity | Specificity | Accuracy | AUC |")
    md <- c(md, "|------|-------|---------|-------------|-------------|----------|-----|")

    for (i in 1:nrow(top10)) {
      r <- top10[i, ]
      md <- c(md, sprintf(
        "| %d | %s | %d | %.1f%% (%.1f–%.1f) | %.1f%% (%.1f–%.1f) | %.1f%% (%.1f–%.1f) | %.3f |",
        i,
        r$model,
        r$markers,
        r$nested_cv_mean_sensitivity * 100,
        (r$nested_cv_mean_sensitivity - 1.96 * r$nested_cv_sd_accuracy / sqrt(10)) * 100,
        (r$nested_cv_mean_sensitivity + 1.96 * r$nested_cv_sd_accuracy / sqrt(10)) * 100,
        r$nested_cv_mean_specificity * 100,
        (r$nested_cv_mean_specificity - 1.96 * r$nested_cv_sd_accuracy / sqrt(10)) * 100,
        (r$nested_cv_mean_specificity + 1.96 * r$nested_cv_sd_accuracy / sqrt(10)) * 100,
        r$nested_cv_mean_accuracy * 100,
        r$nested_cv_95ci_lower * 100,
        r$nested_cv_95ci_upper * 100,
        r$nested_cv_mean_auc
      ))
    }
  } else {
    md <- c(md, "| Rank | Model | Markers | Sensitivity | Specificity | Accuracy | AUC |")
    md <- c(md, "|------|-------|---------|-------------|-------------|----------|-----|")

    for (i in 1:nrow(top10)) {
      r <- top10[i, ]
      md <- c(md, sprintf(
        "| %d | %s | %d | %.1f%% | %.1f%% | %.1f%% | %.3f |",
        i, r$model, r$markers,
        r$mean_sensitivity * 100,
        r$mean_specificity * 100,
        r$accuracy * 100,
        r$auc
      ))
    }
  }

  md <- c(md, "")
  md <- c(md, "*Values shown with 95% confidence intervals where available.*")
  md <- c(md, "")
  md <- c(md, "---")
  md <- c(md, "")

  # Complete Results Table
  md <- c(md, "## Complete Model Comparison")
  md <- c(md, "")
  md <- c(md, sprintf("Full results for all %d models tested are available in `analysis_summary.csv`.", nrow(results)))
  md <- c(md, "")

  # Model Selection Rationale
  md <- c(md, "## Model Selection Rationale")
  md <- c(md, "")
  md <- c(md, "Models were ranked using a composite score prioritising:")
  md <- c(md, "")
  md <- c(md, "1. **Sensitivity** (40% weight) - Critical for detecting positive cases")
  md <- c(md, "2. **Specificity** (40% weight) - Important for minimising false positives")
  md <- c(md, "3. **Accuracy** (20% weight) - Overall correctness")
  md <- c(md, "")

  if (use_nested) {
    md <- c(md, "**Cross-Validation Strategy:**")
    md <- c(md, "")
    md <- c(md, "- **Nested CV** provides unbiased performance estimates suitable for publication")
    md <- c(md, "- **Regular CV** was used for model comparison and development")
    md <- c(md, "- Optimism bias quantified as the difference between regular and nested CV")
    md <- c(md, "")
  }

  md <- c(md, "For models with similar performance, preference was given to those with:")
  md <- c(md, "- Tighter confidence intervals (more stable)")
  md <- c(md, "- Fewer markers (greater parsimony)")
  md <- c(md, "")

  # Technical Details
  md <- c(md, "## Technical Details")
  md <- c(md, "")
  md <- c(md, sprintf("- **DMR Detection Method:** methylKit differential methylation"))
  md <- c(md, sprintf("- **Training Samples:** %d", best$train_samples))
  md <- c(md, sprintf("- **Test Samples:** %d", best$test_samples))
  if (use_nested) {
    md <- c(md, "- **Cross-Validation:** 10-fold nested CV with 5 repeats")
  } else {
    md <- c(md, "- **Cross-Validation:** 10-fold CV with 5 repeats")
  }
  md <- c(md, sprintf("- **Analysis Date:** %s", format(Sys.time(), "%d %B %Y")))
  md <- c(md, "")
  md <- c(md, "---")
  md <- c(md, "")
  md <- c(md, "*Report generated by MethylSense v5.1.0 (Public Release)*")
  md <- c(md, "")

  # Write Markdown
  writeLines(md, md_file)
  log_msg(paste("[REPORT] Markdown report saved to:", basename(md_file)))

  # ===========================================================================
  # GENERATE HTML REPORT
  # ===========================================================================

  html_file <- file.path(output_dir, "ANALYSIS_REPORT.html")

  # Convert markdown to HTML with professional styling
  html <- c("<!DOCTYPE html>")
  html <- c(html, '<html lang="en-GB">')
  html <- c(html, "<head>")
  html <- c(html, '  <meta charset="UTF-8">')
  html <- c(html, '  <meta name="viewport" content="width=device-width, initial-scale=1.0">')
  html <- c(html, "  <title>MethylSense Analysis Report</title>")
  html <- c(html, "  <style>")
  html <- c(html, '    body { font-family: "Helvetica Neue", Arial, sans-serif; line-height: 1.6; max-width: 1200px; margin: 40px auto; padding: 20px; background: #f5f5f5; }')
  html <- c(html, "    .container { background: white; padding: 40px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }")
  html <- c(html, "    h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }")
  html <- c(html, "    h2 { color: #34495e; margin-top: 30px; border-bottom: 2px solid #ecf0f1; padding-bottom: 8px; }")
  html <- c(html, "    h3 { color: #16a085; }")
  html <- c(html, "    table { border-collapse: collapse; width: 100%; margin: 20px 0; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }")
  html <- c(html, "    th { background: #34495e; color: white; padding: 12px; text-align: left; font-weight: 600; }")
  html <- c(html, "    td { padding: 10px; border-bottom: 1px solid #ecf0f1; }")
  html <- c(html, "    tr:hover { background: #f8f9fa; }")
  html <- c(html, "    tr:first-child td { background: #e8f5e9; font-weight: 600; }")
  html <- c(html, "    .metric { background: #ecf0f1; padding: 10px; border-radius: 4px; margin: 5px 0; }")
  html <- c(html, "    .highlight { background: #fff3cd; padding: 2px 6px; border-radius: 3px; }")
  html <- c(html, "    hr { border: none; border-top: 1px solid #dee2e6; margin: 30px 0; }")
  html <- c(html, "    ul { line-height: 1.8; }")
  html <- c(html, "    .timestamp { color: #6c757d; font-size: 0.9em; }")
  html <- c(html, "  </style>")
  html <- c(html, "</head>")
  html <- c(html, "<body>")
  html <- c(html, '<div class="container">')

  # Convert markdown to HTML (simple conversion)
  for (line in md) {
    if (grepl("^# ", line)) {
      html <- c(html, paste0("<h1>", gsub("^# ", "", line), "</h1>"))
    } else if (grepl("^## ", line)) {
      html <- c(html, paste0("<h2>", gsub("^## ", "", line), "</h2>"))
    } else if (grepl("^### ", line)) {
      html <- c(html, paste0("<h3>", gsub("^### ", "", line), "</h3>"))
    } else if (grepl("^\\*\\*.*:\\*\\*$", line)) {
      html <- c(html, paste0("<h4>", gsub("\\*\\*", "", line), "</h4>"))
    } else if (grepl("^- ", line)) {
      if (length(html) == 0 || !grepl("</ul>$", tail(html, 1))) {
        html <- c(html, "<ul>")
      }
      html <- c(html, paste0("<li>", gsub("^- ", "", line), "</li>"))
    } else if (grepl("^\\| ", line) && grepl("\\|$", line)) {
      # Table rows
      if (!grepl("<table>", tail(html, 3))) {
        html <- c(html, "<table>")
      }
      cells <- strsplit(line, "\\|")[[1]]
      cells <- cells[cells != ""]
      cells <- trimws(cells)

      if (grepl("^---", cells[1])) {
        next # Skip separator row
      } else if (grepl("Rank|Model", line)) {
        html <- c(html, "<tr>")
        for (cell in cells) {
          html <- c(html, paste0("<th>", cell, "</th>"))
        }
        html <- c(html, "</tr>")
      } else {
        html <- c(html, "<tr>")
        for (cell in cells) {
          html <- c(html, paste0("<td>", cell, "</td>"))
        }
        html <- c(html, "</tr>")
      }
    } else if (grepl("^---$", line)) {
      if (grepl("</ul>$", tail(html, 1))) {
        # Already closed
      } else if (grepl("<table>", tail(html, 20))) {
        html <- c(html, "</table>")
      }
      html <- c(html, "<hr>")
    } else if (line == "") {
      if (grepl("</ul>$", tail(html, 1))) {
        # Keep ul open
      } else if (grepl("<li>", tail(html, 1))) {
        html <- c(html, "</ul>")
      }
    } else {
      # Regular paragraph
      if (grepl("</ul>$", tail(html, 1))) {
        # ul already closed
      } else if (grepl("<li>", tail(html, 1))) {
        html <- c(html, "</ul>")
      }
      if (grepl("<table>", tail(html, 20)) && !grepl("</table>", tail(html, 1))) {
        html <- c(html, "</table>")
      }
      if (nchar(line) > 0) {
        # Handle bold text
        line <- gsub("\\*\\*(.*?)\\*\\*", "<strong>\\1</strong>", line)
        html <- c(html, paste0("<p>", line, "</p>"))
      }
    }
  }

  html <- c(html, "</div>")
  html <- c(html, "</body>")
  html <- c(html, "</html>")

  # Write HTML
  writeLines(html, html_file)
  log_msg(paste("[REPORT] HTML report saved to:", basename(html_file)))

  return(list(md_file = md_file, html_file = html_file))
}

# Enhanced option_list (keeping existing options)
option_list <- list(
  make_option(c("--qs_file"), type = "character", help = "Path to .qs methylRawList object"),

  # NEW: PREDEFINED TRAIN/TEST SPLIT OPTION
  make_option(c("--train_test_file"),
    type = "character", default = NULL,
    help = "Excel file with predefined train/test assignments (columns: ID, trainTestOverride)"
  ),

  # BED FILE INPUT OPTIONS
  make_option(c("--region_dir"),
    type = "character", default = NULL,
    help = "Directory with windowed BED files for analysis"
  ),
  make_option(c("--window_file"),
    type = "character", default = NULL,
    help = "Single BED file to use for analysis (alternative to --region_dir)"
  ),
  make_option(c("--reproducible"),
    action = "store_true", default = FALSE,
    help = "Enable strict reproducibility mode (disables parallelization) [default: %default]"
  ),
  make_option(c("--filter_standard_chroms"),
    action = "store_true", default = FALSE,
    help = "Filter to standard chromosomes only (1-99, Z, W, X, and Y) [default: %default]"
  ),
  make_option(c("--exclude_sex_chroms"),
    action = "store_true", default = FALSE,
    help = "Exclude sex chromosomes (X, Y, Z, W) from marker selection [default: %default]"
  ),

  # CLUSTERING OPTIONS
  make_option(c("--perform_clustering"),
    action = "store_true", default = TRUE,
    help = "Perform PCA and clustering analysis of DMRs [default: %default]"
  ),
  make_option(c("--n_clusters"),
    type = "integer", default = 0,
    help = "Number of clusters for k-means (0 = auto-determine) [default: %default]"
  ),

  # ANALYSIS PARAMETERS
  make_option(c("--min_accuracy"),
    type = "double", default = 0.7,
    help = "Minimum accuracy threshold to save results [default: %default]"
  ),
  make_option(c("--output_dir"),
    type = "character",
    help = "Output directory for models"
  ),
  make_option(c("--qval"),
    type = "double", default = 0.05,
    help = "Q-value threshold [default: %default]"
  ),
  make_option(c("--marker_range"),
    type = "character", default = "2:50",
    help = "Range of marker set sizes"
  ),
  make_option(c("--models"),
    type = "character", default = "rf,glmnet,xgboost",
    help = "Models: 'all', categories (fast,biomarker,interpretable,etc), or specific models [default: %default]"
  ),
  make_option(c("--list-models"),
    action = "store_true", default = FALSE,
    help = "Display available models and categories, then exit"
  ),
  make_option(c("--min_cpg"),
    type = "integer", default = 1000,
    help = "Min CpGs per sample [default: %default]"
  ),
  make_option(c("--cores"),
    type = "integer", default = 8,
    help = "Number of cores"
  ),
  make_option(c("--min_coverage"),
    type = "integer", default = 10,
    help = "Minimum coverage per region"
  ),
  make_option(c("--min_samples_per_region"),
    type = "double", default = 1.0,
    help = "Min fraction of samples per region"
  ),
  make_option(c("--test_ratio"),
    type = "double", default = 0.2,
    help = "Fraction of data for testing (ignored if --train_test_file specified) [default: %default]"
  ),
  make_option(c("--min_group_size"),
    type = "integer", default = 4,
    help = "Minimum samples required per group for ML analysis [default: %default]"
  ),

  # PLOTTING OPTIONS
  make_option(c("--create_plots"),
    action = "store_true", default = TRUE,
    help = "Create DMR bar plots and heatmaps [default: %default]"
  ),
  make_option(c("--plot_top_n"),
    type = "integer", default = 20,
    help = "Number of top DMRs to plot [default: %default]"
  ),
  make_option(c("--plot_format"),
    type = "character", default = "png",
    help = "Plot output format: png, jpeg, svg, or 'all' for all formats [default: png]"
  ),
  make_option(c("--plot_dpi"),
    type = "integer", default = 300,
    help = "Plot resolution in DPI (dots per inch) [default: 300]"
  ),

  # OTHER OPTIONS
  make_option(c("--dryrun"),
    action = "store_true", default = FALSE,
    help = "Preview only"
  ),
  make_option(c("--verbose"),
    action = "store_true", default = TRUE,
    help = "Verbose output"
  ),
  make_option(c("--debug"),
    action = "store_true", default = FALSE,
    help = "Enable debug output"
  ),
  make_option(c("--force_ml"),
    action = "store_true", default = FALSE,
    help = "Force machine learning analysis even if precomputed files exist"
  ),

  # CV OPTIONS
  make_option(c("--cv_repeats"),
    type = "integer", default = 0,
    help = "Number of cross-validation repeats (0 = disabled) [default: %default]"
  ),
  make_option(c("--nested_cv"),
    action = "store_true", default = FALSE,
    help = "Use nested cross-validation (more rigorous, slower) [default: %default]"
  ),

  # ADD THESE LINES AFTER THE EXISTING make_option() ENTRIES:
  make_option(c("--hyper_threshold"),
    type = "double", default = 5.0,
    help = "Minimum methylation difference for hypermethylation (%) [default: %default]"
  ),
  make_option(c("--hypo_threshold"),
    type = "double", default = -2.0,
    help = "Minimum methylation difference for hypomethylation (%) [default: %default]"
  ),
  make_option(c("--disable_meth_filter"),
    action = "store_true", default = FALSE,
    help = "Disable methylation difference filtering (use only q-value) [default: %default]"
  ),

  ## Group names, pos class and treat mapping!
  make_option(c("--group_names"),
    type = "character", default = NULL,
    help = "Comma-separated group names (e.g., 'Control,Infected' or 'Healthy,Disease')"
  ),
  make_option(c("--group_colors"),
    type = "character", default = NULL,
    help = "Comma-separated color codes matching group names (e.g., '#4575b4,#d73027' or 'blue,red')"
  ),
  make_option(c("--positive_class"),
    type = "character", default = NULL,
    help = "Which group to treat as positive class for ROC (must match one of --group_names)"
  ),
  make_option(c("--treatment_mapping"),
    type = "character", default = NULL,
    help = "Treatment code mapping (e.g., '0=Control,1=Infected' maps numeric codes to names)"
  ),
  # GROUP SIZE OPTIONS
  make_option(c("--allow_small_groups"),
    action = "store_true", default = FALSE,
    help = "Allow groups with fewer samples (risky for train/test split) [default: %default]"
  ),


  ## Software helper functions
  make_option(c("--version"),
    action = "store_true", default = FALSE,
    help = "Show version information and exit"
  ),
  make_option(c("--show_output"),
    action = "store_true", default = FALSE,
    help = "Show detailed output directory structure and exit"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Parse command line options
opt_parser <- OptionParser(
  option_list = option_list,
  usage = "usage: %prog [options]",
  description = "\nMethylSense: Interpretable ML for Epigenetic Diagnostics"
)
opt <- parse_args(opt_parser)

# Handle help flag
if (!is.null(opt$help) && opt$help) {
  show_help()
  quit(save = "no", status = 0)
}

# Handle version flag
if (!is.null(opt$version) && opt$version) {
  cat("\n")
  cat("================================================================================\n")
  cat("MethylSense: Interpretable ML for Epigenetic Diagnostics\n")
  cat("================================================================================\n")
  cat(paste("Version:", SCRIPT_VERSION, "\n"))
  cat(paste("Release Date:", SCRIPT_DATE, "\n"))
  cat("\n")
  cat("For full help: Rscript MethylSense.R --help\n")
  cat("================================================================================\n")
  quit(save = "no", status = 0)
}

# ===== IMMEDIATE CHECK FOR --list-models (must be first) =====
if (!is.null(opt$`list-models`) && opt$`list-models`) {
  cat("\n[MODEL] Available Machine Learning Models (10 total):\n\n")

  cat("Individual Models:\n")
  cat("  rf         : Random Forest - Ensemble method using decision trees (Speed: medium)\n")
  cat("  glmnet     : Elastic Net - Regularized linear regression (Speed: fast)\n")
  cat("  knn        : k-Nearest Neighbors - Instance-based learning (Speed: fast)\n")
  cat("  lda        : Linear Discriminant Analysis - Linear classification (Speed: fast)\n")
  cat("  svm        : Support Vector Machine - SVM with RBF kernel (Speed: slow)\n")
  cat("  xgboost    : XGBoost - Gradient boosting trees, often best performer (Speed: slow)\n")
  cat("  nb         : Naive Bayes - Probabilistic classifier, good for genomic data (Speed: fast)\n")
  cat("  ranger     : Ranger RF - Fast Random Forest implementation (Speed: medium)\n")
  cat("  logreg     : Logistic Regression - Linear logistic regression (Speed: fast)\n")
  cat("  nnet       : Neural Network - Feedforward neural network (Speed: medium)\n")

  cat("\n[DATA] Model Categories:\n")
  cat("  all                  : All 10 models\n")
  cat("  fast                 : knn, lda, glmnet, nb, logreg\n")
  cat("  biomarker           : xgboost, rf, glmnet, ranger, svm\n")
  cat("  interpretable       : glmnet, lda, nb, logreg\n")
  cat("  high_performance    : xgboost, rf, ranger\n")
  cat("  ensemble            : rf, xgboost, ranger\n")
  cat("  linear              : glmnet, lda, logreg\n")
  cat("  nonlinear           : rf, svm, knn, xgboost, ranger, nnet\n")

  cat("\n[INFO] Usage Examples:\n")
  cat("  --models all                  # All 10 models\n")
  cat("  --models biomarker           # Biomarker-optimized models\n")
  cat("  --models fast                # Fast models for quick testing\n")
  cat("  --models xgboost,rf,ranger   # Specific models\n")
  cat("\n")

  quit(save = "no", status = 0)
}
# ===== END immediate check =====
# ===== END --list-models check =====

# Set default strategy (only "traditional" available in public version)
opt$strategy <- "traditional"

# INSERT NEW PARSING LOGIC HERE:
# Parse the new group options

# Parse group names
if (!is.null(opt$group_names)) {
  group_names_list <- trimws(strsplit(opt$group_names, ",")[[1]])
  log_msg(paste("[DATA] Custom group names:", paste(group_names_list, collapse = " vs ")))
} else {
  group_names_list <- NULL
}

# NEW: Parse and validate group colors
group_colors_list <- NULL
if (!is.null(opt$group_colors)) {
  group_colors_list <- trimws(strsplit(opt$group_colors, ",")[[1]])

  # Validate color count matches group count
  if (!is.null(group_names_list) && length(group_colors_list) != length(group_names_list)) {
    stop(paste(
      "[ERROR] Number of colors (", length(group_colors_list),
      ") must match number of groups (", length(group_names_list), ")"
    ))
  }

  # Validate color codes
  invalid_colors <- character()
  for (i in seq_along(group_colors_list)) {
    color <- group_colors_list[i]
    # Try to validate color (R will error if invalid)
    test_color <- tryCatch(
      {
        col2rgb(color)
        TRUE
      },
      error = function(e) {
        FALSE
      }
    )

    if (!test_color) {
      invalid_colors <- c(invalid_colors, color)
    }
  }

  if (length(invalid_colors) > 0) {
    stop(paste(
      "[ERROR] Invalid color codes:", paste(invalid_colors, collapse = ", "),
      "\n[INFO] Use hex codes (#RRGGBB) or R color names (e.g., 'blue', 'red')"
    ))
  }

  # Create named color vector
  names(group_colors_list) <- group_names_list

  log_msg(paste("[COLORS] Custom group colors:"))
  for (i in seq_along(group_names_list)) {
    log_msg(paste("   ", group_names_list[i], "=", group_colors_list[i]))
  }
} else if (!is.null(group_names_list)) {
  log_msg("[COLORS] Using default color scheme (will auto-assign)")
}

if (!is.null(opt$positive_class)) {
  if (is.null(group_names_list)) {
    stop("[ERROR] --positive_class requires --group_names to be specified")
  }
  if (!opt$positive_class %in% group_names_list) {
    stop(paste("[ERROR] --positive_class must be one of:", paste(group_names_list, collapse = ", ")))
  }
  log_msg(paste("[TARGET] Positive class for ROC curves:", opt$positive_class))
} else {
  opt$positive_class <- NULL
}

# Validate methylation thresholds
if (opt$hyper_threshold <= 0) {
  stop("[ERROR] --hyper_threshold must be positive (e.g., 5.0 for 5% hypermethylation)")
}

if (opt$hypo_threshold >= 0) {
  stop("[ERROR] --hypo_threshold must be negative (e.g., -2.0 for 2% hypomethylation)")
}

if (abs(opt$hypo_threshold) > opt$hyper_threshold) {
  log_warning(
    "Hypomethylation threshold has larger absolute value than hypermethylation threshold",
    "This may result in more hypomethylated than hypermethylated DMRs"
  )
}

# Initialize the missing logical variable
if (is.null(opt$disable_meth_filter)) opt$disable_meth_filter <- FALSE

# Initialize any missing logical variables
if (is.null(opt$list_models)) opt$list_models <- FALSE
if (is.null(opt$dryrun)) opt$dryrun <- FALSE
if (is.null(opt$verbose)) opt$verbose <- TRUE
if (is.null(opt$debug)) opt$debug <- FALSE
if (is.null(opt$create_plots)) opt$create_plots <- TRUE
if (is.null(opt$perform_clustering)) opt$perform_clustering <- TRUE
if (is.null(opt$force_ml)) opt$force_ml <- FALSE

# ADD RIGHT AFTER:
if (is.null(opt$cv_repeats)) opt$cv_repeats <- 0
if (is.null(opt$nested_cv)) opt$nested_cv <- FALSE

if (opt$cv_repeats > 0) {
  log_msg(paste("[CV] Cross-validation enabled:", opt$cv_repeats, "repeats"))
  if (opt$nested_cv) {
    log_msg("[CV] Using nested cross-validation (gold standard)")
  }
} else {
  log_msg("[CV] Cross-validation disabled (use --cv_repeats N to enable)")
}

if (is.null(opt$reproducible)) opt$reproducible <- FALSE

# ================================================================================
# SLURM ENVIRONMENT DETECTION AND CORE ALLOCATION
# ================================================================================
# Detect SLURM allocation and adjust cores if running on cluster
if (Sys.getenv("SLURM_CPUS_PER_TASK") != "") {
  slurm_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
  log_msg(paste("[SLURM] Detected SLURM_CPUS_PER_TASK:", slurm_cores))
  if (opt$cores > slurm_cores) {
    log_msg(paste("[SLURM] Reducing --cores from", opt$cores, "to", slurm_cores, "to match SLURM allocation"))
    opt$cores <- slurm_cores
  }
} else if (Sys.getenv("SLURM_JOB_CPUS_PER_NODE") != "") {
  slurm_cores <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
  log_msg(paste("[SLURM] Detected SLURM_JOB_CPUS_PER_NODE:", slurm_cores))
  if (opt$cores > slurm_cores) {
    log_msg(paste("[SLURM] Reducing --cores from", opt$cores, "to", slurm_cores, "to match SLURM allocation"))
    opt$cores <- slurm_cores
  }
}

# REPRODUCIBILITY SETTINGS - CONDITIONAL BASED ON FLAG
if (opt$reproducible) {
  log_msg("[LOCKED] REPRODUCIBILITY MODE ENABLED - Disabling parallelization")
  Sys.setenv(OMP_NUM_THREADS = 1)
  options(mc.cores = 1)
  set.seed(123)
  RNGkind(sample.kind = "Rounding")
  log_msg("   * Single-core mode enforced")
  log_msg("   * RNG seed fixed at 123")
  log_msg("   * All parallel operations disabled")
} else {
  log_msg("[FAST] PERFORMANCE MODE ENABLED - Using parallel processing")
  log_msg(paste("   * Available cores:", opt$cores))

  # Register parallel backend for caret/foreach
  cl <- makeCluster(opt$cores)
  registerDoParallel(cl)

  # SECURITY FIX: Guarantee cluster cleanup even on errors
  on.exit(
    {
      if (exists("cl")) {
        try(stopCluster(cl), silent = TRUE)
        log_msg("[CLEANUP] Parallel cluster stopped")
      }
    },
    add = TRUE
  )

  Sys.setenv(OMP_NUM_THREADS = opt$cores)
  options(mc.cores = opt$cores)
  set.seed(123)
  log_msg("   * Parallel backend registered with doParallel")
  log_msg("   * Note: Results may vary slightly due to parallel execution order")

  # RAM usage warning for high core counts
  if (opt$cores > 4) {
    estimated_ram_gb <- opt$cores * 2
    log_msg(paste0(
      "   * WARNING: Using ", opt$cores, " cores may require approximately ",
      estimated_ram_gb, " GB RAM"
    ))
    log_msg("   * Parallel processing creates multiple dataset copies in memory")
    log_msg("   * Consider reducing --cores if system runs out of memory")
  }
}

# ================================================================================
# ADAPTIVE CV FOLD CALCULATION
# ================================================================================
# Calculate optimal CV folds based on available cores
# Caret parallelises across CV folds, so folds = parallel workers
# Range: 5-10 folds (standard practice, scaled to cores)
cv_folds <- min(10, max(5, opt$cores))

log_msg(paste("[CV] Adaptive CV folds set to:", cv_folds, "based on", opt$cores, "available cores"))
log_msg("   * Caret parallelises across CV folds (not repeats)")
log_msg(paste("   * With", cv_folds, "folds, up to", cv_folds, "cores will be utilised simultaneously"))

if (cv_folds < opt$cores) {
  log_msg(paste("   * NOTE:", opt$cores - cv_folds, "cores may be underutilised (capped at 10-fold CV)"))
}

# ---------------------------
# UPDATED TRAIN/TEST SPLIT FUNCTIONS WITH ORGANISED OUTPUT
# ---------------------------

# Function to load and validate train/test assignments from Excel
load_train_test_assignments <- function(file_path) {
  log_analysis_step("Loading train/test assignments", basename(file_path))

  # Check if file exists
  if (!file.exists(file_path)) {
    log_error(paste("Train/test file not found:", file_path))
    stop(paste("[ERROR] Train/test file not found:", file_path))
  }

  # Try to read Excel file
  tryCatch(
    {
      # Read Excel file
      assignments <- read_excel(file_path)
      log_file_operation("Loaded Excel file", file_path, paste(nrow(assignments), "rows"))

      # Check required columns
      required_cols <- c("ID", "trainTestOverride")
      missing_cols <- setdiff(required_cols, colnames(assignments))

      if (length(missing_cols) > 0) {
        error_msg <- paste("Missing required columns in Excel file:", paste(missing_cols, collapse = ", "))
        log_error(error_msg, "Excel file validation")
        stop(paste("[ERROR]", error_msg, "\n[INFO] Excel file must have columns: ID, trainTestOverride"))
      }

      # Clean and validate data
      assignments <- assignments %>%
        filter(!is.na(ID) & !is.na(trainTestOverride)) %>%
        mutate(
          ID = as.character(ID),
          trainTestOverride = tolower(trimws(as.character(trainTestOverride)))
        )

      # Validate trainTestOverride values
      valid_values <- c("train", "test")
      invalid_assignments <- assignments[!assignments$trainTestOverride %in% valid_values, ]

      if (nrow(invalid_assignments) > 0) {
        log_warning("Invalid trainTestOverride values found", "Data validation")
        for (i in 1:min(5, nrow(invalid_assignments))) {
          log_msg(paste("   Invalid entry - ID:", invalid_assignments$ID[i], "-> Value:", invalid_assignments$trainTestOverride[i]))
        }
        error_msg <- paste("trainTestOverride must be 'train' or 'test' only. Found invalid values in", nrow(invalid_assignments), "rows")
        log_error(error_msg, "Data validation")
        stop(paste("[ERROR]", error_msg))
      }

      # Check for duplicates
      duplicate_ids <- assignments$ID[duplicated(assignments$ID)]
      if (length(duplicate_ids) > 0) {
        error_msg <- paste("Duplicate sample IDs found:", paste(head(duplicate_ids, 5), collapse = ", "))
        log_error(error_msg, "Data validation")
        stop(paste("[ERROR]", error_msg))
      }

      # Summary
      assignment_table <- table(assignments$trainTestOverride)
      log_analysis_step(
        "Train/test assignments validated",
        paste("Train:", assignment_table["train"], "| Test:", assignment_table["test"])
      )

      return(assignments)
    },
    error = function(e) {
      log_error(paste("Error reading Excel file:", e$message), "File I/O")
      stop(paste(
        "[ERROR] Error reading Excel file:", e$message,
        "\n[INFO] Ensure file is a valid Excel (.xlsx) file with columns: ID, trainTestOverride"
      ))
    }
  )
}
# Function to apply predefined train/test split (FIXED)
apply_predefined_split <- function(input_data, assignments, sample_names, output_dirs) {
  log_analysis_step("Applying predefined train/test split")

  # Create sample ID mapping
  sample_mapping <- data.frame(
    methylkit_id = sample_names,
    sample_number = as.character(1:length(sample_names)),
    stringsAsFactors = FALSE
  )

  # Try to match assignments with different ID formats
  matched_assignments <- assignments
  matched_count <- 0

  # Strategy 1: Direct ID matching
  direct_matches <- assignments$ID %in% sample_mapping$methylkit_id
  if (sum(direct_matches) > 0) {
    log_msg(paste("[OK] Direct ID matches:", sum(direct_matches)))
    matched_count <- matched_count + sum(direct_matches)
  }

  # Strategy 2: Sample number matching
  number_matches <- assignments$ID %in% sample_mapping$sample_number
  if (sum(number_matches) > 0) {
    log_msg(paste("[OK] Sample number matches:", sum(number_matches)))
    matched_count <- matched_count + sum(number_matches)
  }

  # Strategy 3: Try with sample_ prefix
  sample_prefix_ids <- paste0("sample_", assignments$ID)
  prefix_matches <- sample_prefix_ids %in% sample_mapping$methylkit_id
  if (sum(prefix_matches) > 0) {
    log_msg(paste("[OK] Sample prefix matches:", sum(prefix_matches)))
    assignments$ID[prefix_matches] <- sample_prefix_ids[prefix_matches]
    matched_count <- matched_count + sum(prefix_matches)
  }

  if (matched_count == 0) {
    log_error("No sample ID matches found", "ID matching")
    debug_msg("Available methylKit sample IDs:")
    for (i in 1:min(10, length(sample_names))) {
      debug_msg(paste("   ", i, ":", sample_names[i]))
    }
    debug_msg("Assignment file sample IDs (first 10):")
    for (i in 1:min(10, nrow(assignments))) {
      debug_msg(paste("   ", i, ":", assignments$ID[i]))
    }
    stop(paste(
      "[ERROR] Could not match any sample IDs between methylKit data and assignment file.",
      "\n[INFO] Check that sample IDs in Excel file match methylKit sample names"
    ))
  }

  # Create final train/test indices
  train_indices <- c()
  test_indices <- c()
  unmatched_samples <- c()

  for (i in 1:length(sample_names)) {
    sample_id <- sample_names[i]
    sample_num <- as.character(i)

    # Find assignment for this sample
    assignment_row <- assignments[assignments$ID == sample_id |
      assignments$ID == sample_num |
      assignments$ID == paste0("sample_", sample_num), ]

    if (nrow(assignment_row) == 1) {
      if (assignment_row$trainTestOverride == "train") {
        train_indices <- c(train_indices, i)
      } else {
        test_indices <- c(test_indices, i)
      }
    } else {
      unmatched_samples <- c(unmatched_samples, sample_id)
    }
  }

  # Validate results
  if (length(unmatched_samples) > 0) {
    log_warning(paste("Unmatched samples (will be excluded):", length(unmatched_samples)), "Split validation")

    # SAVE UNMATCHED SAMPLES TO ORGANIZED DIRECTORY
    if (length(unmatched_samples) > 0) {
      unmatched_file <- file.path(output_dirs$unmatched_samples, "unmatched_samples.csv")
      write.csv(data.frame(unmatched_sample_id = unmatched_samples),
        unmatched_file,
        row.names = FALSE
      )
      log_file_operation("Saved unmatched samples", unmatched_file, paste(length(unmatched_samples), "samples"))
    }
  }

  if (length(train_indices) == 0 || length(test_indices) == 0) {
    log_error("Train/test split resulted in empty train or test set", "Split validation")
    stop("[ERROR] Train/test split resulted in empty train or test set")
  }

  # FIXED: Use input_data parameter instead of global ml_data
  train_labels <- input_data$label[train_indices]
  test_labels <- input_data$label[test_indices]

  train_table <- table(train_labels)
  test_table <- table(test_labels)

  log_analysis_step(
    "Predefined split results",
    paste("Train:", length(train_indices), "| Test:", length(test_indices))
  )

  # FIXED: Return data from input_data parameter, not global ml_data
  return(list(
    train_indices = train_indices,
    test_indices = test_indices,
    train_data = input_data[train_indices, ], #
    test_data = input_data[test_indices, ], #
    matched_samples = length(train_indices) + length(test_indices),
    unmatched_samples = unmatched_samples
  ))
}

# ---------------------------
# HELPER FUNCTIONS
# ---------------------------

# INSERT NEW HELPER FUNCTIONS HERE:

# Simple cache key generator (non-patented version)
generate_regioncounts_cache_key <- function(meth, gr, min_coverage) {
  # Generate a simple hash-based cache key for regionCounts operations
  # This is a basic implementation for the public version

  # Get sample IDs and region count
  n_samples <- length(meth)
  n_regions <- length(gr)

  # Create a simple key based on parameters
  key_string <- paste(
    "samples", n_samples,
    "regions", n_regions,
    "mincov", min_coverage,
    sep = "_"
  )

  # Add a hash of the genomic ranges
  if (requireNamespace("digest", quietly = TRUE)) {
    gr_hash <- digest::digest(paste(seqnames(gr), start(gr), end(gr), collapse = ""))
    key_string <- paste(key_string, substr(gr_hash, 1, 8), sep = "_")
  }

  return(key_string)
}

process_treatment_groups <- function(treatment_info, group_names = NULL, treatment_mapping = NULL) {
  log_analysis_step("Processing treatment groups")

  original_treatment <- treatment_info
  unique_treatments <- unique(original_treatment) # DON'T SORT - keep original order

  log_msg(paste("[DATA] Original treatment codes:", paste(unique_treatments, collapse = ", ")))

  # Apply treatment mapping if provided
  if (!is.null(treatment_mapping)) {
    log_msg("[PROCESS] Applying treatment code mapping...")

    # Parse mapping (format: "0=Control,1=Aspergillus_Infected")
    mapping_pairs <- strsplit(treatment_mapping, ",")[[1]]
    mapping_dict <- list()
    mapped_names_in_order <- character(length(mapping_pairs)) # NEW: Track order

    for (i in seq_along(mapping_pairs)) {
      pair <- mapping_pairs[i]
      parts <- strsplit(trimws(pair), "=")[[1]]
      if (length(parts) == 2) {
        code <- trimws(parts[1])
        name <- trimws(parts[2])
        mapping_dict[[code]] <- name
        mapped_names_in_order[i] <- name # NEW: Store in CLI order
        log_msg(paste("   Code", code, "->", name))
      }
    }

    # Apply mapping
    mapped_treatment <- sapply(original_treatment, function(x) {
      key <- as.character(x)
      if (key %in% names(mapping_dict)) {
        return(mapping_dict[[key]])
      } else {
        return(paste0("Group_", x)) # Fallback
      }
    })

    treatment_labels <- factor(mapped_treatment, levels = mapped_names_in_order) # NEW: Use CLI order
  } else if (!is.null(group_names)) {
    # Use provided group names
    if (length(group_names) != length(unique_treatments)) {
      stop(paste(
        "[ERROR] Number of group names (", length(group_names),
        ") doesn't match number of treatment groups (", length(unique_treatments), ")"
      ))
    }

    log_msg("[PROCESS] Applying custom group names...")
    name_mapping <- setNames(group_names, unique_treatments)

    for (i in seq_along(unique_treatments)) {
      log_msg(paste("   Code", unique_treatments[i], "->", group_names[i]))
    }

    mapped_treatment <- sapply(original_treatment, function(x) name_mapping[[as.character(x)]])
    treatment_labels <- factor(mapped_treatment, levels = group_names)
  } else {
    # Default: use generic Group_ labels
    log_msg("[DATA] Using default group labels...")
    treatment_labels <- factor(paste0("Group_", original_treatment))
  }

  # Log final group distribution
  group_table <- table(treatment_labels)
  log_msg("[DATA] Final treatment groups:")
  for (i in 1:length(group_table)) {
    log_msg(paste("   ", names(group_table)[i], ":", group_table[i], "samples"))
  }

  # NEW: Store group colors if provided
  group_colors <- NULL
  if (exists("group_colors_list", envir = .GlobalEnv) && !is.null(group_colors_list)) {
    group_colors <- group_colors_list
  }

  return(list(
    labels = treatment_labels,
    original_codes = original_treatment,
    group_names = levels(treatment_labels),
    group_colors = group_colors # NEW: Include colors
  ))
}

# ---------------------------
# WINDOW GENERATION FUNCTIONS - NOT INCLUDED
# ---------------------------
# Window generation is handled by external tools (e.g., GenomeToWindows).
# Users must provide their own BED files for genomic window analysis.
# ---------------------------

# Placeholder function for removed window generation
# ---------------------------
# DMR clustering function
# ---------------------------
perform_dmr_clustering <- function(methylation_matrix, dmrs, treatment_groups, output_dir, region_label) {
  log_msg("Performing PCA and clustering analysis of DMRs")

  tryCatch(
    {
      if (ncol(methylation_matrix) < 3 || nrow(methylation_matrix) < 3) {
        verbose_msg("Insufficient data for clustering analysis")
        return(NULL)
      }

      # Ensure dmrs matches the methylation matrix dimensions
      n_dmrs_to_use <- min(nrow(dmrs), ncol(methylation_matrix))
      dmrs_subset <- dmrs[1:n_dmrs_to_use, ]

      # Handle missing values
      meth_matrix_clean <- methylation_matrix[, 1:n_dmrs_to_use, drop = FALSE]

      for (col in 1:ncol(meth_matrix_clean)) {
        col_data <- meth_matrix_clean[, col]
        if (any(is.na(col_data))) {
          median_val <- median(col_data, na.rm = TRUE)
          if (is.na(median_val)) median_val <- 50
          meth_matrix_clean[is.na(meth_matrix_clean[, col]), col] <- median_val
        }
      }

      # Add cluster info
      dmrs_subset$cluster <- 1 # Simple default clustering

      return(list(
        dmrs_with_clusters = dmrs_subset,
        n_clusters = 1
      ))
    },
    error = function(e) {
      log_msg(paste("Error in clustering analysis:", e$message))
      return(NULL)
    }
  )
}

# ---------------------------
# Progress bar functions
# ---------------------------

create_progress_bar <- function(total, prefix = "", width = 50) {
  list(
    total = total,
    current = 0,
    prefix = prefix,
    width = width,
    start_time = Sys.time()
  )
}

update_progress_bar <- function(pb, current = NULL, suffix = "") {
  if (is.null(current)) {
    pb$current <- pb$current + 1
  } else {
    pb$current <- current
  }

  percent <- pb$current / pb$total
  filled_width <- round(pb$width * percent)

  # ASCII-SAFE progress bar
  bar <- paste0(
    "[",
    paste(rep("=", filled_width), collapse = ""),
    paste(rep("-", pb$width - filled_width), collapse = ""),
    "]"
  )

  # Calculate ETA
  elapsed <- as.numeric(Sys.time() - pb$start_time, units = "secs")
  if (pb$current > 0 && percent < 1) {
    eta_secs <- (elapsed / pb$current) * (pb$total - pb$current)
    eta_text <- if (eta_secs > 60) {
      paste0(" ETA: ", round(eta_secs / 60, 1), "m")
    } else {
      paste0(" ETA: ", round(eta_secs), "s")
    }
  } else {
    eta_text <- ""
  }

  # Print progress bar with carriage return
  cat("\r", pb$prefix, " ", bar, " ", round(percent * 100, 1), "% (",
    pb$current, "/", pb$total, ")", eta_text, " ", suffix,
    sep = ""
  )

  # CRITICAL: Force flush the output buffer
  flush.console()

  if (pb$current >= pb$total) {
    cat("\n")
    flush.console()
  }

  return(pb)
}

# ---------------------------
# CONFIDENCE INTERVAL CALCULATION
# ---------------------------
calculate_confidence_intervals <- function(values, confidence = 0.95) {
  # Calculate 95% CI using bootstrap method
  n <- length(values)
  if (n < 2) {
    return(list(mean = mean(values), lower = NA, upper = NA, sd = 0))
  }

  mean_val <- mean(values, na.rm = TRUE)
  sd_val <- sd(values, na.rm = TRUE)
  se_val <- sd_val / sqrt(n)

  # Use t-distribution for small samples
  t_crit <- qt((1 + confidence) / 2, df = n - 1)
  margin <- t_crit * se_val

  list(
    mean = mean_val,
    lower = mean_val - margin,
    upper = mean_val + margin,
    sd = sd_val
  )
}

# ---------------------------
# MODEL TRAINING HELPER FOR CV
# ---------------------------
perform_model_training <- function(model_name, train_data, ctrl, model_specs, tune_length = 3) {
  # Suppress all output
  sink(nullfile())

  model_fit <- NULL

  tryCatch(
    {
      if (model_name == "svm") {
        n_classes <- length(levels(train_data$label))

        if (n_classes == 2) {
          svm_ctrl <- trainControl(
            method = "cv", number = min(cv_folds, nrow(train_data)),
            classProbs = TRUE, summaryFunction = twoClassSummary,
            savePredictions = "final", allowParallel = !opt$reproducible, verboseIter = FALSE
          )
          svm_metric <- "ROC"
        } else {
          svm_ctrl <- trainControl(
            method = "cv", number = min(cv_folds, nrow(train_data)),
            classProbs = TRUE, summaryFunction = multiClassSummary,
            savePredictions = "final", allowParallel = !opt$reproducible, verboseIter = FALSE
          )
          svm_metric <- "Accuracy"
        }

        svm_grid <- expand.grid(sigma = c(0.01, 0.1), C = c(0.1, 1))
        model_fit <- train(label ~ .,
          data = train_data, method = model_specs[[model_name]],
          trControl = svm_ctrl, tuneGrid = svm_grid,
          preProcess = c("center", "scale"), metric = svm_metric
        )
      } else if (model_name == "nnet") {
        nnet_grid <- expand.grid(size = c(3, 5, 10), decay = c(0, 0.1, 0.01))
        model_fit <- train(label ~ .,
          data = train_data, method = "nnet",
          trControl = ctrl, tuneGrid = nnet_grid,
          preProcess = c("center", "scale"),
          metric = "Accuracy", trace = FALSE, linout = FALSE, maxit = 100
        )
      } else if (model_name == "glmnet") {
        model_fit <- train(label ~ .,
          data = train_data, method = "glmnet",
          trControl = ctrl, preProcess = c("center", "scale"),
          metric = "Accuracy", tuneLength = tune_length
        )
      } else if (model_name == "ranger") {
        n_features <- ncol(train_data) - 1
        mtry_values <- unique(c(floor(sqrt(n_features)), floor(n_features / 3), floor(n_features / 2)))
        mtry_values <- mtry_values[mtry_values <= n_features & mtry_values > 0]

        ranger_grid <- expand.grid(
          mtry = mtry_values,
          splitrule = c("gini", "extratrees"),
          min.node.size = c(1, 5)
        )

        model_fit <- train(label ~ .,
          data = train_data, method = "ranger",
          trControl = ctrl, tuneGrid = ranger_grid,
          metric = "Accuracy", importance = "impurity"
        )
      } else {
        # Default training
        model_fit <- train(label ~ .,
          data = train_data,
          method = model_specs[[model_name]],
          trControl = ctrl,
          preProcess = c("center", "scale"),
          metric = "Accuracy",
          tuneLength = tune_length
        )
      }
    },
    error = function(e) {
      sink()
      stop(paste("Model training failed:", e$message))
    }
  )

  sink()
  return(model_fit)
}

# ---------------------------
# CROSS-VALIDATION WITH NESTED CV SUPPORT
# ---------------------------
perform_cross_validation <- function(ml_data_subset, model_name, model_specs,
                                     ctrl, metric, n_repeats, output_dir,
                                     reg_label, n_markers, nested = FALSE) {
  cat("\n")
  log_msg("========== CV INFO START ==========")
  log_msg(paste("Region:", reg_label))
  log_msg(paste("Model:", model_name))
  log_msg(paste("Markers:", n_markers))
  log_msg(paste("Data dimensions:", nrow(ml_data_subset), "rows x", ncol(ml_data_subset), "cols"))
  log_msg(paste("Classes:", paste(levels(ml_data_subset$label), collapse = ", ")))
  log_msg(paste("Class counts:", paste(table(ml_data_subset$label), collapse = ", ")))
  log_msg(paste("N repeats:", n_repeats))
  log_msg(paste("Nested:", nested))
  log_msg("====================================")
  flush.console()

  # Initialize storage for nested CV ROC data
  fold_roc_data <- list()

  if (nested) {
    log_msg(paste("[NESTED CV] Starting", n_repeats, "outer folds for", model_name))
  } else {
    log_msg(paste("[CV] Starting", n_repeats, "repeats for", model_name))
  }

  # Initialize results storage
  cv_results <- data.frame(
    fold_num = integer(),
    accuracy = numeric(),
    auc = numeric(),
    mean_sensitivity = numeric(),
    mean_specificity = numeric(),
    mean_f1 = numeric(),
    train_samples = integer(),
    test_samples = integer(),
    stringsAsFactors = FALSE
  )

  # Additional storage for nested CV
  if (nested) {
    hyperparameter_results <- data.frame(
      fold_num = integer(),
      best_params = character(),
      inner_cv_score = numeric(),
      stringsAsFactors = FALSE
    )
  }

  all_predictions <- list()

  fold_roc_data <- list()


  # Main CV loop
  for (fold_idx in 1:n_repeats) {
    if (nested) {
      verbose_msg(paste("  Outer fold", fold_idx, "of", n_repeats))
    } else {
      verbose_msg(paste("  CV repeat", fold_idx, "of", n_repeats))
    }

    # Suppress warnings and output during training
    suppressWarnings({
      capture.output(
        {
          tryCatch(
            {
              # OUTER SPLIT (same for both regular and nested)
              set.seed(123 + fold_idx * 10)

              if (nested) {
                # For nested CV: Use k-fold instead of random split
                # Create stratified folds
                folds <- createFolds(ml_data_subset$label, k = n_repeats, list = TRUE)
                test_indices <- folds[[fold_idx]]
                train_indices <- setdiff(1:nrow(ml_data_subset), test_indices)
              } else {
                # For regular CV: Random split
                train_indices <- createDataPartition(
                  ml_data_subset$label,
                  p = 0.6,
                  list = FALSE,
                  times = 1
                )
                test_indices <- setdiff(1:nrow(ml_data_subset), train_indices)
              }

              train_data_outer <- ml_data_subset[train_indices, ]
              test_data_outer <- ml_data_subset[test_indices, ]

              # Alias for compatibility
              train_data_cv <- train_data_outer
              test_data_cv <- test_data_outer

              # NESTED CV: Inner loop for hyperparameter tuning
              if (nested) {
                verbose_msg(paste("    Inner CV for hyperparameter tuning..."))

                # Create inner CV control (adaptive folds on training data only)
                inner_ctrl <- trainControl(
                  method = "cv",
                  number = cv_folds,
                  classProbs = TRUE,
                  savePredictions = "none",
                  allowParallel = !opt$reproducible,
                  verboseIter = FALSE
                )

                # Train with hyperparameter search on INNER folds
                model_fit_cv <- perform_model_training(
                  model_name = model_name,
                  train_data = train_data_outer,
                  ctrl = inner_ctrl,
                  model_specs = model_specs,
                  tune_length = 5 # More extensive tuning for nested CV
                )

                # Extract best hyperparameters
                best_params <- as.character(jsonlite::toJSON(model_fit_cv$bestTune))
                inner_score <- max(model_fit_cv$results[[metric]], na.rm = TRUE)

                # Store hyperparameter info
                hyperparameter_results <- rbind(hyperparameter_results, data.frame(
                  fold_num = fold_idx,
                  best_params = best_params,
                  inner_cv_score = inner_score,
                  stringsAsFactors = FALSE
                ))

                verbose_msg(paste("    Best params:", best_params))
              } else {
                # REGULAR CV: Use default tuning
                model_fit_cv <- perform_model_training(
                  model_name = model_name,
                  train_data = train_data_outer,
                  ctrl = ctrl,
                  model_specs = model_specs,
                  tune_length = 3
                )
              }

              # Make predictions on outer test set
              pred_class_cv <- predict(model_fit_cv, test_data_outer)
              pred_prob_cv <- predict(model_fit_cv, test_data_outer, type = "prob")
              actual_labels_cv <- test_data_outer$label

              # Calculate metrics
              accuracy_cv <- sum(pred_class_cv == actual_labels_cv) / length(actual_labels_cv)

              classes <- levels(actual_labels_cv)
              n_classes <- length(classes)

              sensitivity <- numeric(n_classes)
              specificity <- numeric(n_classes)
              f1_score <- numeric(n_classes)

              for (i in seq_along(classes)) {
                class_name <- classes[i]
                tp <- sum(pred_class_cv == class_name & actual_labels_cv == class_name)
                tn <- sum(pred_class_cv != class_name & actual_labels_cv != class_name)
                fp <- sum(pred_class_cv == class_name & actual_labels_cv != class_name)
                fn <- sum(pred_class_cv != class_name & actual_labels_cv == class_name)

                sensitivity[i] <- if ((tp + fn) > 0) tp / (tp + fn) else 0
                specificity[i] <- if ((tn + fp) > 0) tn / (tn + fp) else 0
                precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0

                if (sensitivity[i] + precision > 0) {
                  f1_score[i] <- 2 * (precision * sensitivity[i]) / (precision + sensitivity[i])
                } else {
                  f1_score[i] <- 0
                }
              }

              # Calculate AUC with proper positive class handling (MATCHES MAIN MODEL)
              auc_cv <- 0
              if (n_classes == 2 && ncol(pred_prob_cv) >= 2) {
                tryCatch(
                  {
                    # Determine positive class (SAME logic as main model)
                    if (!is.null(opt$positive_class) && opt$positive_class %in% classes) {
                      pos_class <- opt$positive_class
                      neg_class <- setdiff(classes, opt$positive_class)
                    } else {
                      # Fallback: alphabetically last class
                      pos_class <- sort(classes)[2]
                      neg_class <- sort(classes)[1]
                    }

                    # Find the positive class column
                    if (pos_class %in% colnames(pred_prob_cv)) {
                      roc_obj <- pROC::roc(
                        response = actual_labels_cv,
                        predictor = pred_prob_cv[, pos_class], # <- FIXED: Use positive class
                        levels = c(neg_class, pos_class), # <- FIXED: Explicit ordering
                        direction = "<",
                        quiet = TRUE
                      )
                      auc_cv <- as.numeric(pROC::auc(roc_obj))

                      # Validate AUC
                      if (is.na(auc_cv) || !is.finite(auc_cv)) {
                        auc_cv <- 0
                      }

                      # NEW: Store ROC data for nested CV visualization
                      if (nested) {
                        fold_roc_data[[fold_idx]] <- data.frame(
                          sensitivity = roc_obj$sensitivities,
                          specificity = roc_obj$specificities,
                          stringsAsFactors = FALSE
                        )
                      }
                    }
                  },
                  error = function(e) {
                    auc_cv <- 0
                  }
                )
              } else if (n_classes > 2) {
                # Multi-class: calculate one-vs-rest AUC for each class
                auc_values <- numeric(n_classes)
                for (i in 1:n_classes) {
                  class_name <- classes[i]
                  if (class_name %in% colnames(pred_prob_cv)) {
                    tryCatch(
                      {
                        binary_labels <- ifelse(actual_labels_cv == class_name, 1, 0)
                        roc_obj <- pROC::roc(
                          response = binary_labels,
                          predictor = pred_prob_cv[, class_name],
                          quiet = TRUE
                        )
                        auc_values[i] <- as.numeric(pROC::auc(roc_obj))
                      },
                      error = function(e) {
                        auc_values[i] <- 0
                      }
                    )
                  }
                }
                auc_cv <- mean(auc_values, na.rm = TRUE)
              }


              # Store results
              cv_results <- rbind(cv_results, data.frame(
                fold_num = fold_idx,
                accuracy = accuracy_cv,
                auc = ifelse(is.na(auc_cv) || !is.finite(auc_cv), 0, auc_cv),
                mean_sensitivity = mean(sensitivity, na.rm = TRUE),
                mean_specificity = mean(specificity, na.rm = TRUE),
                mean_f1 = mean(f1_score, na.rm = TRUE),
                train_samples = nrow(train_data_outer),
                test_samples = nrow(test_data_outer),
                stringsAsFactors = FALSE
              ))

              all_predictions[[fold_idx]] <- data.frame(
                fold_num = fold_idx,
                sample_id = rownames(test_data_outer),
                actual = actual_labels_cv,
                predicted = pred_class_cv,
                stringsAsFactors = FALSE
              )
            },
            error = function(e) {
              if (!grepl("unused argument|verbose", e$message)) {
                log_warning(paste("Fold", fold_idx, "failed:", e$message), "CV")
              }
            }
          )
        },
        file = nullfile()
      )
    })
  }

  # ===== CALCULATE SUMMARY STATISTICS =====
  if (nrow(cv_results) == 0) {
    log_error("No CV results collected - all folds failed", "CV")
    return(NULL)
  }

  # Calculate mean and confidence intervals
  summary_stats <- data.frame(
    mean_accuracy = mean(cv_results$accuracy, na.rm = TRUE),
    sd_accuracy = sd(cv_results$accuracy, na.rm = TRUE),
    accuracy_95ci_lower = mean(cv_results$accuracy, na.rm = TRUE) -
      1.96 * sd(cv_results$accuracy, na.rm = TRUE) / sqrt(nrow(cv_results)),
    accuracy_95ci_upper = mean(cv_results$accuracy, na.rm = TRUE) +
      1.96 * sd(cv_results$accuracy, na.rm = TRUE) / sqrt(nrow(cv_results)),
    mean_auc = mean(cv_results$auc, na.rm = TRUE),
    sd_auc = sd(cv_results$auc, na.rm = TRUE),
    auc_95ci_lower = mean(cv_results$auc, na.rm = TRUE) -
      1.96 * sd(cv_results$auc, na.rm = TRUE) / sqrt(nrow(cv_results)),
    auc_95ci_upper = mean(cv_results$auc, na.rm = TRUE) +
      1.96 * sd(cv_results$auc, na.rm = TRUE) / sqrt(nrow(cv_results)),
    mean_sensitivity = mean(cv_results$mean_sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(cv_results$mean_sensitivity, na.rm = TRUE),
    mean_specificity = mean(cv_results$mean_specificity, na.rm = TRUE),
    sd_specificity = sd(cv_results$mean_specificity, na.rm = TRUE),
    mean_f1 = mean(cv_results$mean_f1, na.rm = TRUE),
    sd_f1 = sd(cv_results$mean_f1, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  log_msg(paste(
    "[CV]", model_name, "| Acc:", round(summary_stats$mean_accuracy, 3),
    "95% CI [", round(summary_stats$accuracy_95ci_lower, 3), "-",
    round(summary_stats$accuracy_95ci_upper, 3), "]"
  ))

  # ===== SAVE RESULTS =====
  cv_dir <- file.path(output_dir, "cross_validation")
  dir.create(cv_dir, showWarnings = FALSE, recursive = TRUE)

  # Calculate 95% CIs
  acc_ci <- calculate_confidence_intervals(cv_results$accuracy)
  auc_ci <- calculate_confidence_intervals(cv_results$auc)
  sens_ci <- calculate_confidence_intervals(cv_results$mean_sensitivity)
  spec_ci <- calculate_confidence_intervals(cv_results$mean_specificity)

  # Create summary with CIs
  cv_summary <- data.frame(
    model = model_name,
    region = reg_label,
    markers = n_markers,
    cv_type = ifelse(nested, "nested", "regular"),
    n_folds = n_repeats,
    mean_accuracy = acc_ci$mean,
    accuracy_95ci_lower = acc_ci$lower,
    accuracy_95ci_upper = acc_ci$upper,
    sd_accuracy = acc_ci$sd,
    mean_auc = auc_ci$mean,
    auc_95ci_lower = auc_ci$lower,
    auc_95ci_upper = auc_ci$upper,
    sd_auc = auc_ci$sd,
    mean_sensitivity = sens_ci$mean,
    sensitivity_95ci_lower = sens_ci$lower,
    sensitivity_95ci_upper = sens_ci$upper,
    mean_specificity = spec_ci$mean,
    specificity_95ci_lower = spec_ci$lower,
    specificity_95ci_upper = spec_ci$upper,
    mean_f1 = mean(cv_results$mean_f1, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  # Save results
  cv_dir <- file.path(output_dir, "cross_validation")
  dir.create(cv_dir, recursive = TRUE, showWarnings = FALSE)

  # File naming based on CV type
  prefix <- ifelse(nested, "nested_cv", "cv")

  write.csv(cv_results,
    file.path(cv_dir, paste0(prefix, "_detailed_", model_name, "_", n_markers, "markers.csv")),
    row.names = FALSE
  )

  write.csv(cv_summary,
    file.path(cv_dir, paste0(prefix, "_summary_", model_name, "_", n_markers, "markers.csv")),
    row.names = FALSE
  )

  # Save hyperparameters for nested CV
  if (nested && exists("hyperparameter_results") && nrow(hyperparameter_results) > 0) {
    write.csv(hyperparameter_results,
      file.path(cv_dir, paste0("nested_cv_hyperparameters_", model_name, "_", n_markers, "markers.csv")),
      row.names = FALSE
    )
  }

  if (length(all_predictions) > 0) {
    write.csv(do.call(rbind, all_predictions),
      file.path(cv_dir, paste0(prefix, "_predictions_", model_name, "_", n_markers, "markers.csv")),
      row.names = FALSE
    )
  }

  if (nested) {
    log_msg(paste(
      "[NESTED CV]", model_name, "| Acc:", round(acc_ci$mean, 3),
      "95% CI [", round(acc_ci$lower, 3), "-", round(acc_ci$upper, 3), "]"
    ))
  } else {
    log_msg(paste(
      "[CV]", model_name, "| Acc:", round(acc_ci$mean, 3),
      "95% CI [", round(acc_ci$lower, 3), "-", round(acc_ci$upper, 3), "]"
    ))
  }

  return(list(
    fold_results = cv_results,
    summary = summary_stats,
    fold_roc_data = if (nested) fold_roc_data else NULL # NEW: ROC data for nested CV plots
  ))
}


# ---------------------------
# Plotting functions
# ---------------------------
plot_dmr_barplot <- function(dmrs, treatment_groups, output_dir, region_label, n_markers) {
  tryCatch(
    {
      # ===== GET GROUP ORDER AND COLORS (for consistency, but not used in this plot) =====
      group_names_ordered <- NULL
      group_colors <- NULL

      # Priority 1: Check global environment for treatment_processed
      if (exists("treatment_processed", envir = .GlobalEnv)) {
        tp <- get("treatment_processed", envir = .GlobalEnv)
        if (!is.null(tp) && !is.null(tp$group_names)) {
          group_names_ordered <- tp$group_names
          group_colors <- tp$group_colors
        }
      }

      # Priority 2: Check global environment for opt
      if (is.null(group_names_ordered) && exists("opt", envir = .GlobalEnv)) {
        opt_global <- get("opt", envir = .GlobalEnv)
        if (!is.null(opt_global$group_names)) {
          group_names_ordered <- trimws(strsplit(opt_global$group_names, ",")[[1]])
          if (!is.null(opt_global$group_colors)) {
            group_colors <- trimws(strsplit(opt_global$group_colors, ",")[[1]])
            names(group_colors) <- group_names_ordered
          }
        }
      }

      # Priority 3: Fallback to factor levels
      if (is.null(group_names_ordered)) {
        group_names_ordered <- levels(treatment_groups)
      }

      # Set default colors if needed
      if (is.null(group_colors)) {
        n_groups <- length(group_names_ordered)
        group_colors <- brewer.pal(max(3, n_groups), "Set2")[1:n_groups]
        names(group_colors) <- group_names_ordered
      }

      log_msg(paste("[PLOT] DMR barplot (individual DMRs, not grouped)"))

      # Validate inputs
      if (is.null(dmrs) || nrow(dmrs) < 1) {
        log_warning("No DMRs provided for bar plot", "Plotting")
        return(NULL)
      }

      if (n_markers > nrow(dmrs)) {
        log_warning(paste("Requested", n_markers, "markers but only", nrow(dmrs), "available"), "Plotting")
        n_markers <- nrow(dmrs)
      }

      # Check required columns exist
      required_cols <- c("DMR_ID", "meth.diff", "qvalue")
      missing_cols <- setdiff(required_cols, colnames(dmrs))

      if (length(missing_cols) > 0) {
        log_warning(paste("Missing columns for bar plot:", paste(missing_cols, collapse = ", ")), "Plotting")
        return(NULL)
      }

      # Prepare data for plotting (use top N markers)
      plot_data <- dmrs[1:n_markers, ]

      dmr_plot_data <- data.frame(
        DMR_ID = plot_data$DMR_ID,
        meth_diff = plot_data$meth.diff,
        qvalue = plot_data$qvalue,
        significance = ifelse(plot_data$qvalue < 0.001, "***",
          ifelse(plot_data$qvalue < 0.01, "**",
            ifelse(plot_data$qvalue < 0.05, "*", "ns")
          )
        ),
        direction = ifelse(plot_data$meth.diff > 0, "Hypermethylated", "Hypomethylated"),
        stringsAsFactors = FALSE
      )

      # Create bar plot
      p1 <- ggplot(dmr_plot_data, aes(
        x = reorder(DMR_ID, abs(meth_diff)),
        y = meth_diff,
        fill = direction
      )) +
        geom_bar(stat = "identity", alpha = 0.8) +
        geom_text(aes(label = significance),
          hjust = ifelse(dmr_plot_data$meth_diff > 0, -0.1, 1.1),
          vjust = 0.5, size = 3
        ) +
        scale_fill_manual(values = c(
          "Hypermethylated" = "#d73027",
          "Hypomethylated" = "#4575b4"
        )) +
        # REMOVED: scale_x_discrete (not needed, x-axis is DMR_ID not Group)
        coord_flip() +
        labs(
          title = paste("Top", n_markers, "DMRs -", region_label),
          subtitle = paste("Significance: *** p<0.001, ** p<0.01, * p<0.05"),
          x = "DMR ID",
          y = "Methylation Difference (%)",
          fill = "Direction"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10),
          axis.text.y = element_text(size = 8),
          legend.position = "bottom"
        )

      # Save plot
      plot_file <- file.path(output_dir, paste0("dmr_barplot_", n_markers, "markers.png"))
      save_plot_universal(p1, plot_file, width = 10, height = max(6, n_markers * 0.3))

      log_msg(paste("[PLOT] DMR barplot saved:", basename(plot_file)))
      return(plot_file)
    },
    error = function(e) {
      log_error(paste("Bar plot failed:", e$message), "Plotting")
      log_error(paste("Call:", paste(deparse(e$call), collapse = " ")), "Plotting")
      return(NULL)
    }
  )
}

plot_group_methylation_summary <- function(methylation_matrix, treatment_groups,
                                           dmrs, output_dir, region_label, n_markers,
                                           group_names = NULL, group_colors = NULL) {
  tryCatch(
    {
      # Validate inputs
      if (length(treatment_groups) != nrow(methylation_matrix)) {
        log_warning("Treatment length mismatch in group summary", "Plotting")
        return(NULL)
      }

      # ===== GET GROUP ORDER AND COLORS =====
      group_names_ordered <- NULL
      group_colors_to_use <- NULL

      # Priority 1: Use explicitly passed parameters
      if (!is.null(group_names)) {
        group_names_ordered <- group_names
        group_colors_to_use <- group_colors
      }

      # Priority 2: Check global environment for treatment_processed
      if (is.null(group_names_ordered) && exists("treatment_processed", envir = .GlobalEnv)) {
        tp <- get("treatment_processed", envir = .GlobalEnv)
        if (!is.null(tp) && !is.null(tp$group_names)) {
          group_names_ordered <- tp$group_names
          group_colors_to_use <- tp$group_colors
        }
      }

      # Priority 3: Check global environment for opt
      if (is.null(group_names_ordered) && exists("opt", envir = .GlobalEnv)) {
        opt_global <- get("opt", envir = .GlobalEnv)
        if (!is.null(opt_global$group_names)) {
          group_names_ordered <- trimws(strsplit(opt_global$group_names, ",")[[1]])
          if (!is.null(opt_global$group_colors)) {
            group_colors_to_use <- trimws(strsplit(opt_global$group_colors, ",")[[1]])
            names(group_colors_to_use) <- group_names_ordered
          }
        }
      }

      # Priority 4: Fallback to factor levels
      if (is.null(group_names_ordered)) {
        group_names_ordered <- levels(treatment_groups)
      }

      # Set default colors if needed
      if (is.null(group_colors_to_use)) {
        n_groups <- length(group_names_ordered)
        if (n_groups <= 8) {
          group_colors_to_use <- brewer.pal(max(3, n_groups), "Set2")[1:n_groups]
        } else {
          group_colors_to_use <- rainbow(n_groups)
        }
        names(group_colors_to_use) <- group_names_ordered
      }

      log_msg(paste("[PLOT] Group methylation summary order:", paste(group_names_ordered, collapse = " -> ")))
      log_msg(paste("[PLOT] Using colors:", paste(names(group_colors_to_use), "=",
        group_colors_to_use,
        collapse = ", "
      )))

      # CRITICAL: Convert to character first, then to ordered factor
      treatment_groups <- factor(as.character(treatment_groups), levels = group_names_ordered)

      # DEBUG: Comprehensive group matching diagnostics
      log_msg(paste("[DEBUG] treatment_groups class:", class(treatment_groups)))
      log_msg(paste("[DEBUG] treatment_groups levels:", paste(levels(treatment_groups), collapse = ", ")))
      log_msg(paste("[DEBUG] treatment_groups first 10 values:", paste(head(as.character(treatment_groups), 10), collapse = ", ")))
      log_msg(paste("[DEBUG] group_names_ordered:", paste(group_names_ordered, collapse = ", ")))
      log_msg(paste("[DEBUG] methylation_matrix dimensions:", nrow(methylation_matrix), "rows x", ncol(methylation_matrix), "columns"))
      log_msg(paste("[DEBUG] treatment_groups length:", length(treatment_groups)))

      # Get colors
      if (!is.null(group_colors)) {
        # Use explicitly passed colors
        group_colors_to_use <- group_colors
      } else if (exists("treatment_processed", envir = .GlobalEnv) &&
        !is.null(treatment_processed) &&
        !is.null(treatment_processed$group_colors)) {
        # Use from treatment_processed
        group_colors_to_use <- treatment_processed$group_colors
      } else {
        # Generate default colors
        n_groups <- length(group_names_ordered)
        if (n_groups <= 8) {
          group_colors_to_use <- brewer.pal(max(3, n_groups), "Set2")[1:n_groups]
        } else {
          group_colors_to_use <- rainbow(n_groups)
        }
        names(group_colors_to_use) <- group_names_ordered
      }

      log_msg(paste("[PLOT] Using colors:", paste(names(group_colors_to_use), "=",
        group_colors_to_use,
        collapse = ", "
      )))

      # Rest of function...

      # Calculate mean and SD methylation for each group
      summary_data <- data.frame()

      for (group in group_names_ordered) {
        # ROBUST: Try as.character conversion for matching
        log_msg(paste("[DEBUG] Searching for group:", group))
        log_msg(paste("[DEBUG] treatment_groups unique values:", paste(unique(as.character(treatment_groups)), collapse = ", ")))
        log_msg(paste("[DEBUG] Exact comparison:", as.character(group), "==", paste(head(as.character(treatment_groups), 3), collapse = ", ")))

        group_indices <- which(as.character(treatment_groups) == as.character(group))

        log_msg(paste("[DEBUG] Group:", group, "| Matched indices:", length(group_indices), "| Sample IDs:", paste(head(group_indices, 5), collapse = ", ")))

        if (length(group_indices) == 0) {
          log_warning(paste("No samples found for group:", group, "- adding placeholder with NA"), "Plotting")
          # Add placeholder to maintain structure and prevent complete NA propagation
          summary_data <- rbind(summary_data, data.frame(
            Group = group,
            Mean_Methylation = NA_real_,
            SD = NA_real_,
            SE = NA_real_,
            N = 0,
            stringsAsFactors = FALSE
          ))
          next
        }

        group_matrix <- methylation_matrix[group_indices, , drop = FALSE]
        sample_means <- rowMeans(group_matrix, na.rm = TRUE)

        group_mean <- mean(sample_means, na.rm = TRUE)
        group_sd <- sd(sample_means, na.rm = TRUE)
        group_n <- length(sample_means)
        group_se <- group_sd / sqrt(group_n)

        log_msg(paste("[DEBUG] Group:", group, "| N:", group_n, "| Mean:", round(group_mean, 2), "| SD:", round(group_sd, 2)))

        summary_data <- rbind(summary_data, data.frame(
          Group = group,
          Mean_Methylation = group_mean,
          SD = group_sd,
          SE = group_se,
          N = group_n,
          stringsAsFactors = FALSE
        ))
      }

      # CRITICAL: Ensure Group column maintains order
      summary_data$Group <- factor(summary_data$Group, levels = group_names_ordered)

      # UPDATED: Get colors from treatment_processed or use defaults
      n_groups <- nrow(summary_data)

      if (exists("treatment_processed") && !is.null(treatment_processed) &&
        !is.null(treatment_processed$group_colors)) {
        # Use custom colors
        group_colors <- treatment_processed$group_colors
        log_msg(paste("[PLOT] Using custom colors:", paste(names(group_colors), "=", group_colors, collapse = ", ")))
      } else {
        # Default colors
        if (n_groups <= 8) {
          group_colors <- brewer.pal(max(3, n_groups), "Set2")[1:n_groups]
        } else {
          group_colors <- rainbow(n_groups)
        }
        names(group_colors) <- group_names_ordered
      }

      # Create plot with custom colors
      # CRITICAL: Ensure summary_data Group is an ordered factor
      summary_data$Group <- factor(summary_data$Group, levels = group_names_ordered)

      # DEBUG: Verify order
      log_msg(paste("[DEBUG] Summary data groups:", paste(levels(summary_data$Group), collapse = ", ")))
      log_msg(paste("[DEBUG] Summary data rows:", nrow(summary_data)))

      # Create plot with preserved order
      p_summary <- ggplot(summary_data, aes(x = Group, y = Mean_Methylation, fill = Group)) +
        geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
        geom_errorbar(aes(ymin = Mean_Methylation - SD, ymax = Mean_Methylation + SD),
          width = 0.2, linewidth = 0.8, color = "black"
        ) +
        geom_text(aes(label = paste0("n=", N)),
          vjust = -0.5, hjust = 0.5, size = 3.5,
          y = summary_data$Mean_Methylation + summary_data$SD + 2
        ) +
        scale_fill_manual(
          values = group_colors_to_use,
          limits = group_names_ordered, # CRITICAL
          breaks = group_names_ordered, # CRITICAL
          drop = FALSE
        ) +
        scale_x_discrete(
          limits = group_names_ordered, # CRITICAL: Force x-axis order
          breaks = group_names_ordered, # CRITICAL
          drop = FALSE
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        labs(
          title = paste("Mean Methylation Across", n_markers, "Significant DMRs"),
          subtitle = paste("Region:", region_label, "| Error bars: +/- 1 SD"),
          x = "Treatment Groups",
          y = "Mean Methylation Level (%)",
          caption = paste("Based on", n_markers, "top DMRs")
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 11),
          legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()
        )

      # Save plot
      summary_plot_file <- file.path(output_dir, paste0(
        "group_methylation_summary_",
        n_markers, "markers.png"
      ))

      save_plot_universal(p_summary, summary_plot_file, width = 10, height = 8)

      # Save stats
      summary_stats_file <- file.path(output_dir, paste0(
        "group_methylation_stats_",
        n_markers, "markers.csv"
      ))
      write.csv(summary_data, summary_stats_file, row.names = FALSE)

      log_msg(paste("[PLOT] Group summary saved:", basename(summary_plot_file)))
      return(list(plot = summary_plot_file, stats = summary_stats_file))
    },
    error = function(e) {
      log_error(paste("Group summary failed:", e$message), "Plotting")
      return(NULL)
    }
  )
}

# Enhanced ROC curve plotting function with proper binary/multiclass handling
plot_roc_curves <- function(pred_prob, actual_labels, model_name, output_dir,
                            region_label, n_markers, positive_class = NULL) {
  tryCatch(
    {
      # ===== STEP 1: GET GROUP ORDER AND COLORS =====
      # Try multiple sources in priority order
      group_names_ordered <- NULL
      color_palette <- NULL

      # Priority 1: Check global environment for treatment_processed
      if (exists("treatment_processed", envir = .GlobalEnv)) {
        tp <- get("treatment_processed", envir = .GlobalEnv)
        if (!is.null(tp) && !is.null(tp$group_names)) {
          group_names_ordered <- tp$group_names
          color_palette <- tp$group_colors
        }
      }

      # Priority 2: Check global environment for opt
      if (is.null(group_names_ordered) && exists("opt", envir = .GlobalEnv)) {
        opt_global <- get("opt", envir = .GlobalEnv)
        if (!is.null(opt_global$group_names)) {
          group_names_ordered <- trimws(strsplit(opt_global$group_names, ",")[[1]])
          if (!is.null(opt_global$group_colors)) {
            color_palette <- trimws(strsplit(opt_global$group_colors, ",")[[1]])
            names(color_palette) <- group_names_ordered
          }
        }
      }

      # Priority 3: Use actual_labels factor levels as fallback
      if (is.null(group_names_ordered)) {
        group_names_ordered <- levels(actual_labels)
      }

      # Set default colors if not found
      n_classes <- length(group_names_ordered)
      if (is.null(color_palette)) {
        color_palette <- c("#d73027", "#4575b4", "#1a9850", "#fee090", "#9970ab")[1:n_classes]
        names(color_palette) <- group_names_ordered
      }

      log_msg(paste("[ROC] Group order:", paste(group_names_ordered, collapse = " -> ")))
      log_msg(paste("[ROC] Colors:", paste(names(color_palette), "=", color_palette, collapse = ", ")))

      # ===== STEP 2: ENSURE ACTUAL_LABELS IS CORRECTLY ORDERED =====
      actual_labels <- factor(as.character(actual_labels), levels = group_names_ordered)

      classes <- levels(actual_labels)
      n_classes <- length(classes)

      if (n_classes == 2) {
        # ===== BINARY CLASSIFICATION =====
        log_msg(paste("[ROC] Binary classification:", paste(group_names_ordered, collapse = " vs ")))

        # Determine positive and negative classes
        if (!is.null(positive_class) && positive_class %in% group_names_ordered) {
          pos_class <- positive_class
          neg_class <- setdiff(group_names_ordered, positive_class)
        } else {
          # Default: SECOND group in ordered list is positive
          pos_class <- group_names_ordered[2]
          neg_class <- group_names_ordered[1]
        }

        log_msg(paste("[ROC] Reference (negative):", neg_class))
        log_msg(paste("[ROC] Target (positive):", pos_class))
        log_msg(paste("[ROC] Using positive class from CLI:", positive_class))

        # Verify pred_prob has the positive class column
        if (!pos_class %in% colnames(pred_prob)) {
          log_error(paste("Positive class", pos_class, "not found in prediction probabilities"), "ROC")
          log_msg(paste("[DEBUG] Available columns:", paste(colnames(pred_prob), collapse = ", ")))
          return(list(auc = 0))
        }

        # CRITICAL: Create ROC with POSITIVE CLASS as SECOND level
        # pROC expects: levels = c(control/negative, case/positive)
        roc_obj <- pROC::roc(
          response = actual_labels,
          predictor = pred_prob[, pos_class],
          levels = c(neg_class, pos_class), # CRITICAL: negative first, positive second
          direction = "<",
          quiet = TRUE
        )

        auc_value <- as.numeric(pROC::auc(roc_obj))

        # Validate AUC
        if (is.na(auc_value) || !is.finite(auc_value)) {
          log_warning("Invalid AUC value, setting to 0", "ROC")
          auc_value <- 0
        }

        log_msg(paste("[OK] Binary ROC AUC:", round(auc_value, 3)))

        # Get color for positive class
        roc_color <- if (pos_class %in% names(color_palette)) {
          color_palette[pos_class]
        } else {
          "#d73027" # Fallback red
        }

        # Create binary ROC plot
        roc_data <- data.frame(
          sensitivity = roc_obj$sensitivities,
          specificity = roc_obj$specificities
        )

        p_binary <- ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity)) +
          geom_line(color = roc_color, size = 1.2) +
          geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
          annotate("text",
            x = 0.7, y = 0.3,
            label = paste0("AUC = ", round(auc_value, 3)),
            size = 5, fontface = "bold"
          ) +
          annotate("text",
            x = 0.7, y = 0.2,
            label = paste0("Positive: ", pos_class),
            size = 4, color = roc_color, fontface = "bold"
          ) +
          annotate("text",
            x = 0.7, y = 0.15,
            label = paste0("Reference: ", neg_class),
            size = 3.5, color = "gray40"
          ) +
          labs(
            title = paste("Binary ROC Curve:", model_name),
            subtitle = paste(pos_class, "(target) vs", neg_class, "(reference) |", n_markers, "markers"),
            x = "1 - Specificity (False Positive Rate)",
            y = "Sensitivity (True Positive Rate)"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            legend.position = "none"
          ) +
          coord_equal() +
          xlim(0, 1) +
          ylim(0, 1)

        binary_roc_file <- file.path(output_dir, "roc_curve_binary.png")
        save_plot_universal(p_binary, binary_roc_file, width = 8, height = 8)

        log_msg(paste("[PLOT] Binary ROC curve saved:", basename(binary_roc_file)))

        return(list(
          binary_roc = binary_roc_file,
          auc = auc_value,
          roc_obj = roc_obj
        ))
      } else {
        # ===== MULTI-CLASS CLASSIFICATION =====
        log_msg(paste("[ROC] Multi-class classification:", n_classes, "classes"))
        log_msg(paste("[ROC] Multi-class order:", paste(group_names_ordered, collapse = " -> ")))

        all_roc_data <- list()
        auc_values <- numeric(n_classes)
        names(auc_values) <- group_names_ordered

        # Calculate one-vs-rest ROC for each class IN YOUR ORDER
        for (i in seq_along(group_names_ordered)) {
          class_name <- group_names_ordered[i]

          # Create binary labels: current class vs all others
          binary_labels <- factor(
            ifelse(actual_labels == class_name, class_name, "Other"),
            levels = c("Other", class_name)
          )

          # Get predicted probabilities for this class
          if (class_name %in% colnames(pred_prob)) {
            class_probs <- pred_prob[, class_name]

            roc_obj <- pROC::roc(
              response = binary_labels,
              predictor = class_probs,
              levels = c("Other", class_name),
              direction = "<",
              quiet = TRUE
            )

            auc_values[class_name] <- as.numeric(pROC::auc(roc_obj))

            all_roc_data[[class_name]] <- data.frame(
              sensitivity = roc_obj$sensitivities,
              specificity = roc_obj$specificities,
              class = class_name,
              auc = auc_values[class_name],
              stringsAsFactors = FALSE
            )
          }
        }

        # Combine all ROC data with ORDERED factor
        roc_data_combined <- do.call(rbind, all_roc_data)
        roc_data_combined$class <- factor(roc_data_combined$class, levels = group_names_ordered)

        # Create multi-class ROC plot
        p_multiclass <- ggplot(
          roc_data_combined,
          aes(
            x = 1 - specificity, y = sensitivity,
            color = class, linetype = class
          )
        ) +
          geom_line(size = 1.2) +
          geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
          scale_color_manual(
            values = color_palette,
            limits = group_names_ordered,
            breaks = group_names_ordered,
            labels = paste0(group_names_ordered, " (AUC=", round(auc_values[group_names_ordered], 3), ")"),
            drop = FALSE
          ) +
          scale_linetype_manual(
            values = rep("solid", n_classes),
            limits = group_names_ordered,
            breaks = group_names_ordered,
            labels = paste0(group_names_ordered, " (AUC=", round(auc_values[group_names_ordered], 3), ")"),
            drop = FALSE
          ) +
          labs(
            title = paste("Multi-class ROC Curves:", model_name),
            subtitle = paste(
              "One-vs-Rest |", n_markers, "markers | Mean AUC =",
              round(mean(auc_values), 3)
            ),
            x = "1 - Specificity (False Positive Rate)",
            y = "Sensitivity (True Positive Rate)",
            color = "Class",
            linetype = "Class"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            legend.position = "bottom"
          ) +
          coord_equal() +
          xlim(0, 1) +
          ylim(0, 1)

        # Remove the .png extension from multiclass_roc_file definition first:
        multiclass_roc_file <- file.path(output_dir, "roc_curves_multiclass") # No .png!
        save_plot_universal(p_multiclass, multiclass_roc_file, width = 10, height = 8)

        mean_auc <- mean(auc_values, na.rm = TRUE)

        log_msg(paste(
          "[OK] Multi-class AUC values (in order):",
          paste(group_names_ordered, "=", round(auc_values[group_names_ordered], 3), collapse = ", ")
        ))
        log_msg(paste("[PLOT] Multi-class ROC curves saved:", basename(multiclass_roc_file)))

        return(list(
          multiclass_roc = multiclass_roc_file,
          auc = mean_auc,
          auc_per_class = auc_values
        ))
      }
    },
    error = function(e) {
      log_error(paste("ROC curve creation failed:", e$message), "ROC plotting")
      return(list(auc = 0))
    }
  )
}

# ================================================================================
# PRECISION-RECALL CURVES - PUBLICATION QUALITY
# ================================================================================

plot_pr_curves <- function(pred_prob, actual_labels, model_name, output_dir,
                           region_label, n_markers, positive_class = NULL) {
  tryCatch(
    {
      # Get group order and colors (same as ROC)
      group_names_ordered <- NULL
      color_palette <- NULL

      if (exists("treatment_processed", envir = .GlobalEnv)) {
        tp <- get("treatment_processed", envir = .GlobalEnv)
        if (!is.null(tp) && !is.null(tp$group_names)) {
          group_names_ordered <- tp$group_names
          color_palette <- tp$group_colors
        }
      }

      if (is.null(group_names_ordered) && exists("opt", envir = .GlobalEnv)) {
        opt_global <- get("opt", envir = .GlobalEnv)
        if (!is.null(opt_global$group_names)) {
          group_names_ordered <- trimws(strsplit(opt_global$group_names, ",")[[1]])
          if (!is.null(opt_global$group_colors)) {
            color_palette <- trimws(strsplit(opt_global$group_colors, ",")[[1]])
            names(color_palette) <- group_names_ordered
          }
        }
      }

      if (is.null(group_names_ordered)) {
        group_names_ordered <- levels(actual_labels)
      }

      n_classes <- length(group_names_ordered)
      if (is.null(color_palette)) {
        color_palette <- c("#d73027", "#4575b4", "#1a9850", "#fee090", "#9970ab")[1:n_classes]
        names(color_palette) <- group_names_ordered
      }

      log_msg(paste("[PR] Group order:", paste(group_names_ordered, collapse = " -> ")))
      log_msg(paste("[PR] Colors:", paste(names(color_palette), "=", color_palette, collapse = ", ")))

      actual_labels <- factor(as.character(actual_labels), levels = group_names_ordered)

      if (n_classes == 2) {
        # Binary PR curve
        if (!is.null(positive_class) && positive_class %in% group_names_ordered) {
          pos_class <- positive_class
          neg_class <- setdiff(group_names_ordered, positive_class)
        } else {
          pos_class <- group_names_ordered[2]
          neg_class <- group_names_ordered[1]
        }

        log_msg(paste("[PR] Binary - Positive class:", pos_class))

        if (!pos_class %in% colnames(pred_prob)) {
          log_error(paste("Positive class", pos_class, "not found in predictions"), "PR")
          return(list(pr_auc = 0))
        }

        binary_labels <- ifelse(actual_labels == pos_class, 1, 0)
        pr_curve <- PRROC::pr.curve(
          scores.class0 = pred_prob[, pos_class],
          weights.class0 = binary_labels,
          curve = TRUE
        )

        pr_auc <- pr_curve$auc.integral
        log_msg(paste("[OK] Binary PR-AUC:", round(pr_auc, 3)))

        pr_color <- if (pos_class %in% names(color_palette)) {
          color_palette[pos_class]
        } else {
          "#d73027"
        }

        pr_data <- data.frame(
          recall = pr_curve$curve[, 1],
          precision = pr_curve$curve[, 2]
        )

        p_binary <- ggplot(pr_data, aes(x = recall, y = precision)) +
          geom_line(color = pr_color, size = 1.2) +
          geom_hline(
            yintercept = sum(binary_labels) / length(binary_labels),
            linetype = "dashed", color = "gray50"
          ) +
          annotate("text",
            x = 0.3, y = 0.3,
            label = paste0("PR-AUC = ", round(pr_auc, 3)),
            size = 5, fontface = "bold"
          ) +
          annotate("text",
            x = 0.3, y = 0.2,
            label = paste0("Positive: ", pos_class),
            size = 4, color = pr_color, fontface = "bold"
          ) +
          labs(
            title = paste("Precision-Recall Curve -", model_name),
            subtitle = paste(region_label, "|", n_markers, "markers |", pos_class, "vs", neg_class),
            x = "Recall (Sensitivity)",
            y = "Precision (PPV)"
          ) +
          scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
          scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 10),
            axis.title = element_text(size = 12, face = "bold"),
            axis.text = element_text(size = 10)
          )

        pr_file <- file.path(output_dir, paste0("pr_curve_", pos_class))
        save_plot_universal(p_binary, pr_file, width = 8, height = 7)
        log_msg(paste("[PLOT] Binary PR curve saved:", basename(pr_file)))

        return(list(binary_pr = pr_file, pr_auc = pr_auc))
      } else {
        # Multi-class PR curves
        log_msg(paste("[PR] Multi-class:", n_classes, "classes"))

        pr_auc_values <- numeric(n_classes)
        names(pr_auc_values) <- group_names_ordered
        pr_curves_list <- list()

        for (i in seq_along(group_names_ordered)) {
          class_name <- group_names_ordered[i]
          binary_labels <- ifelse(actual_labels == class_name, 1, 0)

          if (class_name %in% colnames(pred_prob)) {
            pr_curve <- PRROC::pr.curve(
              scores.class0 = pred_prob[, class_name],
              weights.class0 = binary_labels,
              curve = TRUE
            )
            pr_auc_values[i] <- pr_curve$auc.integral
            pr_curves_list[[class_name]] <- pr_curve
          } else {
            pr_auc_values[i] <- 0
          }
        }

        mean_pr_auc <- mean(pr_auc_values, na.rm = TRUE)
        log_msg(paste("[OK] Multi-class mean PR-AUC:", round(mean_pr_auc, 3)))

        # Create multi-class plot
        pr_data_list <- list()
        for (class_name in group_names_ordered) {
          if (class_name %in% names(pr_curves_list)) {
            pr_data_list[[class_name]] <- data.frame(
              recall = pr_curves_list[[class_name]]$curve[, 1],
              precision = pr_curves_list[[class_name]]$curve[, 2],
              class = class_name
            )
          }
        }

        pr_data_combined <- do.call(rbind, pr_data_list)
        pr_data_combined$class <- factor(pr_data_combined$class, levels = group_names_ordered)

        p_multi <- ggplot(pr_data_combined, aes(x = recall, y = precision, color = class)) +
          geom_line(size = 1.2) +
          scale_color_manual(
            values = color_palette,
            labels = paste0(group_names_ordered, " (PR-AUC=", round(pr_auc_values[group_names_ordered], 3), ")"),
            name = "Class (One-vs-Rest)"
          ) +
          labs(
            title = paste("Precision-Recall Curves -", model_name),
            subtitle = paste(
              "One-vs-Rest |", n_markers, "markers | Mean PR-AUC =",
              round(mean_pr_auc, 3)
            ),
            x = "Recall (Sensitivity)",
            y = "Precision (PPV)"
          ) +
          scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
          scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 10),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "right"
          )

        multiclass_pr_file <- file.path(output_dir, "pr_curves_multiclass")
        save_plot_universal(p_multi, multiclass_pr_file, width = 10, height = 7)
        log_msg(paste("[PLOT] Multi-class PR curves saved:", basename(multiclass_pr_file)))

        return(list(
          multiclass_pr = multiclass_pr_file,
          pr_auc = mean_pr_auc,
          pr_auc_per_class = pr_auc_values
        ))
      }
    },
    error = function(e) {
      log_error(paste("PR curve creation failed:", e$message), "PR plotting")
      return(list(pr_auc = 0))
    }
  )
}

# ================================================================================
# PREDICTION PROBABILITY RIDGE PLOTS - FIXED VERSION WITH ERROR HANDLING
# ================================================================================

plot_prediction_probabilities <- function(pred_prob, actual_labels, model_name,
                                          output_dir, region_label, n_markers) {
  tryCatch(
    {
      # ===== STEP 1: GET GROUP ORDER AND COLORS =====
      group_names_ordered <- NULL
      group_colors <- NULL

      # Priority 1: Check global environment for treatment_processed
      if (exists("treatment_processed", envir = .GlobalEnv)) {
        tp <- get("treatment_processed", envir = .GlobalEnv)
        if (!is.null(tp) && !is.null(tp$group_names)) {
          group_names_ordered <- tp$group_names
          group_colors <- tp$group_colors
        }
      }

      # Priority 2: Check global environment for opt
      if (is.null(group_names_ordered) && exists("opt", envir = .GlobalEnv)) {
        opt_global <- get("opt", envir = .GlobalEnv)
        if (!is.null(opt_global$group_names)) {
          group_names_ordered <- trimws(strsplit(opt_global$group_names, ",")[[1]])
          if (!is.null(opt_global$group_colors)) {
            group_colors <- trimws(strsplit(opt_global$group_colors, ",")[[1]])
            names(group_colors) <- group_names_ordered
          }
        }
      }

      # Priority 3: Fallback to factor levels
      if (is.null(group_names_ordered)) {
        group_names_ordered <- levels(actual_labels)
      }

      # Set default colors if needed
      if (is.null(group_colors)) {
        n_groups <- length(group_names_ordered)
        group_colors <- brewer.pal(max(3, n_groups), "Set2")[1:n_groups]
        names(group_colors) <- group_names_ordered
      }

      log_msg(paste("[PLOT] Ridge plot order:", paste(group_names_ordered, collapse = " -> ")))
      log_msg(paste(
        "[PLOT] Ridge plot colors:",
        paste(names(group_colors), "=", group_colors, collapse = ", ")
      ))

      # CRITICAL: Ensure actual_labels is correctly ordered
      actual_labels <- factor(as.character(actual_labels), levels = group_names_ordered)

      # ===== STEP 2: PREPARE DATA FOR PLOTTING =====
      prob_data <- data.frame(
        sample_id = rownames(pred_prob),
        actual = actual_labels,
        pred_prob,
        stringsAsFactors = FALSE
      )

      # Reshape for plotting
      prob_long <- reshape2::melt(prob_data,
        id.vars = c("sample_id", "actual"),
        variable.name = "predicted_class",
        value.name = "probability"
      )

      # CRITICAL: Convert to character first, then ordered factor
      prob_long$predicted_class <- as.character(prob_long$predicted_class)
      prob_long$actual <- as.character(prob_long$actual)

      prob_long$predicted_class <- factor(prob_long$predicted_class,
        levels = group_names_ordered
      )
      prob_long$actual <- factor(prob_long$actual, levels = group_names_ordered)

      # ===== STEP 3: DIAGNOSTIC CHECKS =====
      log_msg("[DEBUG] Ridge plot data diagnostics:")
      log_msg(paste("  Total rows:", format(nrow(prob_long), big.mark = ",")))
      log_msg(paste(
        "  Predicted class levels:",
        paste(levels(prob_long$predicted_class), collapse = ", ")
      ))
      log_msg(paste(
        "  Actual class levels:",
        paste(levels(prob_long$actual), collapse = ", ")
      ))

      # Check data availability per predicted class
      pred_class_table <- table(prob_long$predicted_class, useNA = "ifany")
      log_msg("  Predictions per class:")
      for (class_name in names(pred_class_table)) {
        log_msg(paste("    ", class_name, ":", pred_class_table[class_name], "rows"))
      }

      # Check for missing classes
      missing_classes <- character()
      for (class_name in group_names_ordered) {
        n_predictions <- sum(prob_long$predicted_class == class_name, na.rm = TRUE)
        if (n_predictions == 0) {
          missing_classes <- c(missing_classes, class_name)
          log_msg(paste("  [WARN] Zero predictions for class:", class_name))
        } else {
          class_data <- prob_long[prob_long$predicted_class == class_name, ]
          prob_range <- range(class_data$probability, na.rm = TRUE)
          log_msg(paste(
            "  [OK]", class_name, ":", n_predictions, "predictions |",
            "Prob range:", round(prob_range[1], 3), "-", round(prob_range[2], 3)
          ))
        }
      }

      if (length(missing_classes) > 0) {
        log_msg(paste(
          "[WARN] Model never predicts these classes:",
          paste(missing_classes, collapse = ", ")
        ))
        log_msg("       Ridge plot will show empty rows for these classes")
      }

      # Check actual class distribution in prob_long
      actual_class_table <- table(prob_long$actual, useNA = "ifany")
      log_msg("  Actual class distribution in plot data:")
      for (class_name in names(actual_class_table)) {
        log_msg(paste("    ", class_name, ":", actual_class_table[class_name], "rows"))
      }

      # ===== PLOT 1: SIMPLE RIDGE PLOT (FIXED VERSION) =====
      if (requireNamespace("ggridges", quietly = TRUE)) {
        tryCatch(
          {
            log_msg("[PLOT] Creating simple ridge plot...")

            # Check if we have any data at all
            if (nrow(prob_long) == 0) {
              log_msg("[ERROR] No data available for ridge plot")
              return(NULL)
            }

            # Create plot with defensive settings
            p_prob_ridge_simple <- ggplot(
              prob_long,
              aes(
                x = probability,
                y = predicted_class,
                fill = actual
              )
            ) +
              ggridges::geom_density_ridges(
                alpha = 0.7,
                scale = 0.9,
                rel_min_height = 0.001 # Show even very small densities
              ) +
              scale_fill_manual(
                values = group_colors,
                limits = group_names_ordered,
                breaks = group_names_ordered,
                drop = FALSE, # Keep all levels even if no data
                na.value = "gray50" # Color for any NA values
              ) +
              scale_y_discrete(
                limits = rev(group_names_ordered), # Reverse for top-to-bottom
                breaks = rev(group_names_ordered),
                drop = FALSE # Keep all levels even if no data
              ) +
              scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
              labs(
                title = paste("Prediction Probability Distributions:", model_name),
                subtitle = paste(
                  n_markers, "markers |", region_label,
                  if (length(missing_classes) > 0) {
                    paste(
                      "\nNote: Model never predicts:",
                      paste(missing_classes, collapse = ", ")
                    )
                  } else {
                    ""
                  }
                ),
                x = "Prediction Probability",
                y = "Predicted Class",
                fill = "Actual Class"
              ) +
              theme_minimal() +
              theme(
                plot.title = element_text(size = 14, face = "bold"),
                plot.subtitle = element_text(size = 10),
                legend.position = "bottom",
                legend.title = element_text(size = 11, face = "bold")
              )

            # Save plot
            prob_ridge_simple_file <- file.path(
              output_dir,
              "prediction_probabilities_ridge_simple"
            )
            save_plot_universal(p_prob_ridge_simple, prob_ridge_simple_file,
              width = 10, height = 6
            )
            log_msg(paste(
              "[OK] Simple ridge plot saved:",
              basename(prob_ridge_simple_file)
            ))
          },
          error = function(e) {
            log_msg(paste("[ERROR] Simple ridge plot failed:", e$message))
            log_msg(paste("       Error details:", toString(e)))
          }
        )
      } else {
        log_msg("[WARN] ggridges package not available - skipping ridge plots")
      }

      # ===== PLOT 2: RIDGE PLOT WITH JITTERED POINTS (FIXED VERSION) =====
      if (requireNamespace("ggridges", quietly = TRUE)) {
        tryCatch(
          {
            log_msg("[PLOT] Creating ridge plot with points...")

            # Check ggridges version for compatibility
            ggridges_version <- packageVersion("ggridges")
            has_jitter_position <- ggridges_version >= "0.5.0"

            if (has_jitter_position) {
              # NEWER VERSION: Use built-in jittered points
              p_prob_ridge_jitter <- ggplot(
                prob_long,
                aes(
                  x = probability,
                  y = predicted_class,
                  fill = actual
                )
              ) +
                ggridges::geom_density_ridges(
                  alpha = 0.6,
                  scale = 0.9,
                  jittered_points = TRUE,
                  point_shape = "|",
                  point_size = 2,
                  point_alpha = 0.5,
                  position = ggridges::position_points_jitter(width = 0.01, height = 0),
                  rel_min_height = 0.001
                ) +
                scale_fill_manual(
                  values = group_colors,
                  limits = group_names_ordered,
                  breaks = group_names_ordered,
                  drop = FALSE,
                  na.value = "gray50"
                ) +
                scale_y_discrete(
                  limits = rev(group_names_ordered),
                  breaks = rev(group_names_ordered),
                  drop = FALSE
                ) +
                scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
                labs(
                  title = paste(
                    "Prediction Probability Distributions (with points):",
                    model_name
                  ),
                  subtitle = paste(
                    n_markers, "markers |", region_label,
                    "| Vertical lines show individual samples",
                    if (length(missing_classes) > 0) {
                      paste(
                        "\nNote: Model never predicts:",
                        paste(missing_classes, collapse = ", ")
                      )
                    } else {
                      ""
                    }
                  ),
                  x = "Prediction Probability",
                  y = "Predicted Class",
                  fill = "Actual Class"
                ) +
                theme_minimal() +
                theme(
                  plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 9),
                  legend.position = "bottom"
                )
            } else {
              # OLDER VERSION: Manual overlay with geom_point
              # Create numeric y-positions for points
              y_positions <- setNames(
                1:length(group_names_ordered),
                rev(group_names_ordered)
              )

              prob_long$y_numeric <- y_positions[as.character(prob_long$predicted_class)]

              p_prob_ridge_jitter <- ggplot(
                prob_long,
                aes(
                  x = probability,
                  y = predicted_class,
                  fill = actual,
                  color = actual
                )
              ) +
                # First: density ridges
                ggridges::geom_density_ridges(
                  alpha = 0.5,
                  scale = 0.9,
                  rel_min_height = 0.001
                ) +
                # Second: points at the BASE of each ridge
                geom_point(
                  aes(y = y_numeric),
                  position = position_jitter(width = 0, height = 0.08),
                  alpha = 0.7,
                  size = 2.5,
                  shape = 16
                ) +
                scale_fill_manual(
                  values = group_colors,
                  limits = group_names_ordered,
                  breaks = group_names_ordered,
                  drop = FALSE,
                  na.value = "gray50"
                ) +
                scale_color_manual(
                  values = group_colors,
                  limits = group_names_ordered,
                  breaks = group_names_ordered,
                  drop = FALSE,
                  na.value = "gray50"
                ) +
                scale_y_continuous(
                  breaks = 1:length(group_names_ordered),
                  labels = rev(group_names_ordered),
                  limits = c(0.5, length(group_names_ordered) + 0.5)
                ) +
                scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
                labs(
                  title = paste(
                    "Prediction Probability Distributions (with points):",
                    model_name
                  ),
                  subtitle = paste(
                    n_markers, "markers |", region_label,
                    "| Points show individual samples"
                  ),
                  x = "Prediction Probability",
                  y = "Predicted Class",
                  fill = "Actual Class",
                  color = "Actual Class"
                ) +
                theme_minimal() +
                theme(
                  plot.title = element_text(size = 14, face = "bold"),
                  legend.position = "bottom"
                ) +
                guides(
                  fill = guide_legend(override.aes = list(alpha = 0.7)),
                  color = guide_legend(override.aes = list(alpha = 0.9, size = 4))
                )
            }

            prob_ridge_jitter_file <- file.path(
              output_dir,
              "prediction_probabilities_ridge_with_points"
            )
            save_plot_universal(p_prob_ridge_jitter, prob_ridge_jitter_file,
              width = 10, height = 6
            )
            log_msg(paste(
              "[OK] Ridge plot with points saved:",
              basename(prob_ridge_jitter_file)
            ))
          },
          error = function(e) {
            log_msg(paste("[ERROR] Ridge plot with points failed:", e$message))
            log_msg(paste("       Error details:", toString(e)))
          }
        )
      }

      # ===== PLOT 3: VIOLIN PLOT (YOUR WORKING VERSION) =====
      tryCatch(
        {
          log_msg("[PLOT] Creating violin plot...")

          p_prob_violin <- ggplot(
            prob_long,
            aes(
              x = predicted_class,
              y = probability,
              fill = actual
            )
          ) +
            geom_violin(alpha = 0.7, trim = FALSE, scale = "width") +
            geom_boxplot(width = 0.1, alpha = 0.5, outlier.alpha = 0.3) +
            scale_fill_manual(
              values = group_colors,
              limits = group_names_ordered,
              breaks = group_names_ordered,
              drop = FALSE
            ) +
            scale_x_discrete(
              limits = group_names_ordered,
              breaks = group_names_ordered,
              drop = FALSE
            ) +
            labs(
              title = paste("Prediction Probability Distributions:", model_name),
              subtitle = paste(n_markers, "markers |", region_label),
              x = "Predicted Class",
              y = "Prediction Probability",
              fill = "Actual Class"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = 14, face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom"
            ) +
            ylim(0, 1)

          prob_violin_file <- file.path(output_dir, "prediction_probabilities_violin")
          save_plot_universal(p_prob_violin, prob_violin_file, width = 10, height = 8)
          log_msg(paste("[OK] Violin plot saved:", basename(prob_violin_file)))
        },
        error = function(e) {
          log_msg(paste("[ERROR] Violin plot failed:", e$message))
        }
      )

      # ===== PLOT 3B: BAR PLOT WITH ERROR BARS (ALTERNATIVE) =====
      tryCatch(
        {
          log_msg("[PLOT] Creating bar plot with error bars...")

          # Calculate summary statistics
          bar_summary <- prob_long %>%
            group_by(predicted_class, actual) %>%
            summarise(
              mean_prob = mean(probability, na.rm = TRUE),
              sd_prob = sd(probability, na.rm = TRUE),
              n = n(),
              se_prob = sd_prob / sqrt(n),
              .groups = "drop"
            )

          # Ensure factors are ordered
          bar_summary$predicted_class <- factor(bar_summary$predicted_class,
            levels = group_names_ordered
          )
          bar_summary$actual <- factor(bar_summary$actual,
            levels = group_names_ordered
          )

          # Create bar plot with error bars
          p_prob_bar <- ggplot(
            bar_summary,
            aes(
              x = predicted_class,
              y = mean_prob,
              fill = actual
            )
          ) +
            geom_bar(
              stat = "identity",
              position = position_dodge(width = 0.9),
              alpha = 0.8,
              width = 0.8
            ) +
            geom_errorbar(
              aes(
                ymin = pmax(0, mean_prob - se_prob),
                ymax = pmin(1, mean_prob + se_prob)
              ),
              position = position_dodge(width = 0.9),
              width = 0.25,
              linewidth = 0.6
            ) +
            geom_text(aes(label = paste0("n=", n)),
              position = position_dodge(width = 0.9),
              vjust = -1.5,
              size = 3
            ) +
            scale_fill_manual(
              values = group_colors,
              limits = group_names_ordered,
              breaks = group_names_ordered,
              drop = FALSE
            ) +
            scale_x_discrete(
              limits = group_names_ordered,
              breaks = group_names_ordered,
              drop = FALSE
            ) +
            scale_y_continuous(
              limits = c(0, 1.2),
              expand = c(0, 0),
              breaks = seq(0, 1, 0.2)
            ) +
            geom_hline(
              yintercept = 0.5,
              linetype = "dashed",
              color = "gray40",
              alpha = 0.7
            ) +
            labs(
              title = paste("Mean Prediction Probability by Class:", model_name),
              subtitle = paste(
                n_markers, "markers |", region_label,
                "| Error bars show +/- SE | Dashed line at 0.5"
              ),
              x = "Predicted Class",
              y = "Mean Prediction Probability",
              fill = "Actual Class"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = 14, face = "bold"),
              plot.subtitle = element_text(size = 9),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
              axis.text.y = element_text(size = 11),
              legend.position = "bottom",
              legend.title = element_text(size = 11, face = "bold"),
              panel.grid.major.x = element_blank(),
              panel.grid.minor = element_blank()
            )

          prob_bar_file <- file.path(output_dir, "prediction_probabilities_barplot")
          save_plot_universal(p_prob_bar, prob_bar_file, width = 10, height = 8)
          log_msg(paste("[OK] Bar plot with error bars saved:", basename(prob_bar_file)))
        },
        error = function(e) {
          log_msg(paste("[ERROR] Bar plot failed:", e$message))
        }
      )

      # ===== PLOT 4: CONFIDENCE PLOT =====
      tryCatch(
        {
          log_msg("[PLOT] Creating confidence plot...")

          # Calculate max probability and predicted class for each sample
          prob_cols <- setdiff(colnames(prob_data), c("sample_id", "actual"))
          prob_data$max_prob <- apply(prob_data[, prob_cols], 1, max)
          prob_data$predicted <- prob_cols[apply(prob_data[, prob_cols], 1, which.max)]
          prob_data$correct <- (as.character(prob_data$actual) == prob_data$predicted)

          # Ensure actual is ordered factor
          prob_data$actual <- factor(as.character(prob_data$actual),
            levels = group_names_ordered
          )

          p_confidence <- ggplot(prob_data, aes(x = actual, y = max_prob, color = correct)) +
            geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
            geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
            scale_color_manual(
              values = c("FALSE" = "#d73027", "TRUE" = "#1a9850"),
              labels = c("FALSE" = "Incorrect", "TRUE" = "Correct")
            ) +
            scale_x_discrete(
              limits = group_names_ordered,
              breaks = group_names_ordered,
              drop = FALSE
            ) +
            labs(
              title = paste("Prediction Confidence:", model_name),
              subtitle = paste(n_markers, "markers |", region_label),
              x = "Actual Class",
              y = "Maximum Prediction Probability",
              color = "Prediction"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = 14, face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom"
            ) +
            ylim(0, 1)

          confidence_file <- file.path(output_dir, "prediction_confidence")
          save_plot_universal(p_confidence, confidence_file, width = 10, height = 8)
          log_msg(paste("[OK] Confidence plot saved:", basename(confidence_file)))
        },
        error = function(e) {
          log_msg(paste("[ERROR] Confidence plot failed:", e$message))
        }
      )

      # ===== RETURN SUMMARY =====
      return(list(
        ridge_simple = if (exists("prob_ridge_simple_file")) prob_ridge_simple_file else NULL,
        ridge_jitter = if (exists("prob_ridge_jitter_file")) prob_ridge_jitter_file else NULL,
        violin = if (exists("prob_violin_file")) prob_violin_file else NULL,
        confidence = if (exists("confidence_file")) confidence_file else NULL
      ))
    },
    error = function(e) {
      log_error(paste("Prediction probability plots failed:", e$message), "Plotting")
      log_error(paste("Error details:", toString(e)), "Plotting")
      return(NULL)
    }
  )
}

# ================================================================================
# NESTED CV VISUALIZATION SUITE
# ================================================================================

plot_nested_cv_results <- function(cv_result_regular, cv_result_nested,
                                   model_name, output_dir, region_label, n_markers) {
  if (is.null(cv_result_nested)) {
    log_msg("[PLOT] No nested CV results to plot")
    return(NULL)
  }

  log_msg("[PLOT] Creating nested CV visualization suite...")

  # Create nested CV output directory
  nested_dir <- file.path(output_dir, "nested_cv_visualizations")
  dir.create(nested_dir, showWarnings = FALSE, recursive = TRUE)

  plot_files <- list()

  # ===== PLOT 1: Comparison Forest Plot =====
  tryCatch(
    {
      if (!is.null(cv_result_regular)) {
        # Prepare data for all metrics
        comparison_data <- data.frame(
          Method = rep(c("Regular CV", "Nested CV"), 3),
          Metric = rep(c("Accuracy", "Sensitivity", "Specificity"), each = 2),
          Mean = c(
            cv_result_regular$summary$mean_accuracy, cv_result_nested$summary$mean_accuracy,
            cv_result_regular$summary$mean_sensitivity, cv_result_nested$summary$mean_sensitivity,
            cv_result_regular$summary$mean_specificity, cv_result_nested$summary$mean_specificity
          ),
          Lower = c(
            cv_result_regular$summary$accuracy_95ci_lower, cv_result_nested$summary$accuracy_95ci_lower,
            cv_result_regular$summary$mean_sensitivity - 1.96 * cv_result_regular$summary$sd_sensitivity / sqrt(nrow(cv_result_regular$fold_results)),
            cv_result_nested$summary$mean_sensitivity - 1.96 * cv_result_nested$summary$sd_sensitivity / sqrt(nrow(cv_result_nested$fold_results)),
            cv_result_regular$summary$mean_specificity - 1.96 * cv_result_regular$summary$sd_specificity / sqrt(nrow(cv_result_regular$fold_results)),
            cv_result_nested$summary$mean_specificity - 1.96 * cv_result_nested$summary$sd_specificity / sqrt(nrow(cv_result_nested$fold_results))
          ),
          Upper = c(
            cv_result_regular$summary$accuracy_95ci_upper, cv_result_nested$summary$accuracy_95ci_upper,
            cv_result_regular$summary$mean_sensitivity + 1.96 * cv_result_regular$summary$sd_sensitivity / sqrt(nrow(cv_result_regular$fold_results)),
            cv_result_nested$summary$mean_sensitivity + 1.96 * cv_result_nested$summary$sd_sensitivity / sqrt(nrow(cv_result_nested$fold_results)),
            cv_result_regular$summary$mean_specificity + 1.96 * cv_result_regular$summary$sd_specificity / sqrt(nrow(cv_result_regular$fold_results)),
            cv_result_nested$summary$mean_specificity + 1.96 * cv_result_nested$summary$sd_specificity / sqrt(nrow(cv_result_nested$fold_results))
          ),
          CV_Type = factor(rep(c("Regular CV", "Nested CV"), 3),
            levels = c("Regular CV", "Nested CV")
          )
        )

        # Order metrics
        comparison_data$Metric <- factor(comparison_data$Metric,
          levels = c("Accuracy", "Sensitivity", "Specificity")
        )

        # Calculate optimism bias
        bias_accuracy <- cv_result_regular$summary$mean_accuracy - cv_result_nested$summary$mean_accuracy
        bias_sensitivity <- cv_result_regular$summary$mean_sensitivity - cv_result_nested$summary$mean_sensitivity
        bias_specificity <- cv_result_regular$summary$mean_specificity - cv_result_nested$summary$mean_specificity

        p_comparison <- ggplot(comparison_data, aes(x = Mean, y = Metric, color = CV_Type, shape = CV_Type)) +
          geom_point(size = 4, position = position_dodge(width = 0.5)) +
          geom_errorbarh(aes(xmin = Lower, xmax = Upper),
            height = 0.3, size = 1,
            position = position_dodge(width = 0.5)
          ) +
          geom_vline(xintercept = 0.5, linetype = "dotted", color = "gray70", alpha = 0.5) +
          scale_color_manual(
            values = c("Regular CV" = "#4575b4", "Nested CV" = "#d73027"),
            labels = c(
              "Regular CV" = "Regular CV (may be optimistic)",
              "Nested CV" = "Nested CV (unbiased estimate)"
            )
          ) +
          scale_shape_manual(
            values = c("Regular CV" = 16, "Nested CV" = 17),
            labels = c(
              "Regular CV" = "Regular CV (may be optimistic)",
              "Nested CV" = "Nested CV (unbiased estimate)"
            )
          ) +
          labs(
            title = paste("Cross-Validation Performance Comparison:", model_name),
            subtitle = paste(
              n_markers, "markers |", region_label,
              "\nOptimism bias: Accuracy =", sprintf("%+.3f", bias_accuracy),
              "| Sensitivity =", sprintf("%+.3f", bias_sensitivity),
              "| Specificity =", sprintf("%+.3f", bias_specificity)
            ),
            x = "Performance Metric (95% CI)",
            y = "",
            color = "Method",
            shape = "Method"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 10, hjust = 0.5),
            axis.text.y = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.title = element_blank(),
            panel.grid.major.y = element_line(color = "gray90"),
            panel.grid.minor = element_blank()
          ) +
          xlim(0, 1)

        comparison_file <- file.path(nested_dir, "cv_comparison_forest.png")
        save_plot_universal(p_comparison, comparison_file, width = 10, height = 6)
        plot_files$comparison <- comparison_file

        log_msg("[PLOT] Forest plot created with Accuracy, Sensitivity, and Specificity")
      }
    },
    error = function(e) {
      log_error(paste("Forest plot failed:", e$message), "Plotting")
    }
  )
  # ===== PLOT 2: Nested CV ROC Curves (all outer folds) =====
  tryCatch(
    {
      if (!is.null(cv_result_nested$fold_roc_data)) {
        # Combine ROC data from all folds
        all_roc_data <- do.call(rbind, lapply(1:length(cv_result_nested$fold_roc_data), function(i) {
          fold_data <- cv_result_nested$fold_roc_data[[i]]
          if (!is.null(fold_data)) {
            fold_data$fold <- paste("Fold", i)
            fold_data$auc <- cv_result_nested$fold_results$auc[i]
            return(fold_data)
          }
          return(NULL)
        }))

        if (!is.null(all_roc_data) && nrow(all_roc_data) > 0) {
          p_nested_roc <- ggplot(
            all_roc_data,
            aes(
              x = 1 - specificity, y = sensitivity,
              group = fold, color = fold
            )
          ) +
            geom_line(alpha = 0.6, size = 1) +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
            annotate("text",
              x = 0.7, y = 0.3,
              label = paste0(
                "Mean AUC = ",
                round(cv_result_nested$summary$mean_auc, 3),
                "\n95% CI [",
                round(cv_result_nested$summary$auc_95ci_lower, 3), "-",
                round(cv_result_nested$summary$auc_95ci_upper, 3), "]"
              ),
              size = 4, fontface = "bold"
            ) +
            labs(
              title = paste("Nested CV ROC Curves:", model_name),
              subtitle = paste(
                n_markers, "markers |", region_label,
                "| Each curve is an independent outer fold"
              ),
              x = "1 - Specificity (False Positive Rate)",
              y = "Sensitivity (True Positive Rate)",
              color = "Outer Fold"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = 14, face = "bold"),
              legend.position = "right"
            ) +
            coord_equal() +
            xlim(0, 1) +
            ylim(0, 1)

          nested_roc_file <- file.path(nested_dir, "nested_cv_roc_all_folds.png")
          save_plot_universal(p_nested_roc, nested_roc_file, width = 10, height = 8)
          plot_files$nested_roc <- nested_roc_file

          log_msg("[PLOT] Nested ROC curves created")
        }
      }
    },
    error = function(e) {
      log_error(paste("Nested ROC plot failed:", e$message), "Plotting")
    }
  )
  # ===== PREPARE DATA FOR PLOTS 3 & 4 =====
  tryCatch(
    {
      fold_data <- cv_result_nested$fold_results
      fold_data$fold <- 1:nrow(fold_data)

      # DEBUG: Show what columns actually exist
      log_msg(paste("[DEBUG] Available columns in fold_results:", paste(colnames(fold_data), collapse = ", ")))
      log_msg(paste("[DEBUG] Looking for: accuracy, mean_sensitivity, mean_specificity, auc"))

      # Reshape for plotting - CRITICAL: use correct column names
      fold_long <- reshape2::melt(fold_data,
        id.vars = "fold",
        measure.vars = c("accuracy", "mean_sensitivity", "mean_specificity", "auc"),
        variable.name = "metric",
        value.name = "value"
      )

      fold_long$metric <- factor(fold_long$metric,
        levels = c("accuracy", "mean_sensitivity", "mean_specificity", "auc"),
        labels = c("Accuracy", "Sensitivity", "Specificity", "AUC")
      )

      # ===== PLOT 3: Fold-by-Fold Performance =====
      p_fold_performance <- ggplot(fold_long, aes(x = fold, y = value, color = metric, group = metric)) +
        geom_line(size = 1) +
        geom_point(size = 3) +
        geom_hline(data = data.frame(
          metric = c("Accuracy", "Sensitivity", "Specificity", "AUC"),
          mean_val = c(
            cv_result_nested$summary$mean_accuracy,
            cv_result_nested$summary$mean_sensitivity,
            cv_result_nested$summary$mean_specificity,
            cv_result_nested$summary$mean_auc
          )
        ), aes(yintercept = mean_val, color = metric), linetype = "dashed", alpha = 0.5) +
        facet_wrap(~metric, scales = "free_y", nrow = 2) +
        scale_color_brewer(palette = "Set1") +
        labs(
          title = paste("Nested CV Performance by Outer Fold:", model_name),
          subtitle = paste(
            n_markers, "markers |", region_label,
            "| Dashed lines show mean across folds"
          ),
          x = "Outer Fold",
          y = "Value",
          color = "Metric"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          legend.position = "none",
          strip.text = element_text(size = 12, face = "bold")
        ) +
        ylim(0, 1)

      fold_performance_file <- file.path(nested_dir, "nested_cv_fold_performance.png")
      save_plot_universal(p_fold_performance, fold_performance_file, width = 12, height = 8)
      plot_files$fold_performance <- fold_performance_file

      log_msg("[PLOT] Fold performance plot created")

      # ===== PLOT 4: Metric Distributions (Raincloud Plot) =====
      p_raincloud <- ggplot(fold_long, aes(x = metric, y = value, fill = metric)) +
        geom_violin(alpha = 0.5, trim = FALSE) +
        geom_boxplot(width = 0.2, alpha = 0.7, outlier.alpha = 0) +
        geom_jitter(width = 0.1, alpha = 0.6, size = 2) +
        stat_summary(
          fun = mean, geom = "point", shape = 23, size = 4,
          fill = "white", color = "black"
        ) +
        scale_fill_brewer(palette = "Set2") +
        labs(
          title = paste("Nested CV Performance Distribution:", model_name),
          subtitle = paste(
            n_markers, "markers |", region_label,
            "| Diamond = mean, box = IQR, dots = individual folds"
          ),
          x = "",
          y = "Value",
          fill = "Metric"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(size = 12, face = "bold"),
          legend.position = "none"
        ) +
        ylim(0, 1)

      raincloud_file <- file.path(nested_dir, "nested_cv_distributions.png")
      save_plot_universal(p_raincloud, raincloud_file, width = 10, height = 8)
      plot_files$distributions <- raincloud_file

      log_msg("[PLOT] Distribution plot created")
    },
    error = function(e) {
      log_error(paste("Fold plots failed:", e$message), "Plotting")
    }
  )
  # ===== PLOT 5: Optimism Bias Analysis =====
  tryCatch(
    {
      if (!is.null(cv_result_regular)) {
        bias_data <- data.frame(
          Metric = c("Accuracy", "Sensitivity", "Specificity", "AUC"),
          Regular = c(
            cv_result_regular$summary$mean_accuracy,
            cv_result_regular$summary$mean_sensitivity,
            cv_result_regular$summary$mean_specificity,
            cv_result_regular$summary$mean_auc
          ),
          Nested = c(
            cv_result_nested$summary$mean_accuracy,
            cv_result_nested$summary$mean_sensitivity,
            cv_result_nested$summary$mean_specificity,
            cv_result_nested$summary$mean_auc
          )
        )

        bias_data$Bias <- bias_data$Regular - bias_data$Nested
        bias_data$Metric <- factor(bias_data$Metric,
          levels = c("Accuracy", "Sensitivity", "Specificity", "AUC")
        )

        p_bias <- ggplot(bias_data, aes(x = Metric, y = Bias, fill = Metric)) +
          geom_bar(stat = "identity", alpha = 0.8) +
          geom_hline(yintercept = 0, linetype = "solid", color = "black") +
          geom_hline(
            yintercept = c(-0.05, 0.05), linetype = "dashed",
            color = "red", alpha = 0.5
          ) +
          geom_text(aes(label = sprintf("%.3f", Bias)),
            vjust = ifelse(bias_data$Bias > 0, -0.5, 1.5), size = 4
          ) +
          scale_fill_brewer(palette = "Set2") +
          labs(
            title = paste("Optimism Bias Analysis:", model_name),
            subtitle = paste(
              n_markers, "markers |", region_label,
              "| Positive = Regular CV overestimates | Red lines = +/-5% threshold"
            ),
            x = "",
            y = "Optimism Bias (Regular CV - Nested CV)",
            fill = "Metric"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            axis.text.x = element_text(size = 12, face = "bold"),
            legend.position = "none"
          ) +
          ylim(min(c(bias_data$Bias, -0.1)), max(c(bias_data$Bias, 0.1)))

        bias_file <- file.path(nested_dir, "optimism_bias_analysis.png")
        save_plot_universal(p_bias, bias_file, width = 10, height = 6)
        plot_files$bias_analysis <- bias_file

        log_msg("[PLOT] Bias analysis plot created")
      }
    },
    error = function(e) {
      log_error(paste("Bias plot failed:", e$message), "Plotting")
    }
  )

  log_msg(paste("[PLOT] Created", length(plot_files), "nested CV visualizations"))
  log_msg(paste("[DIR] Saved to:", basename(nested_dir)))

  return(plot_files)
}

# ================================================================================
# RANDOM FOREST VISUALIZATION SUITE
# ================================================================================
plot_random_forest_viz <- function(model_fit, train_data, dmrs, output_dir, region_label, n_markers) {
  tryCatch(
    {
      if (model_fit$method != "rf" && model_fit$method != "ranger") {
        return(NULL) # Only works for Random Forest
      }

      log_msg("[PLOT] Creating Random Forest visualization suite...")

      # Extract the actual randomForest object
      if (model_fit$method == "rf") {
        rf_model <- model_fit$finalModel
      } else if (model_fit$method == "ranger") {
        rf_model <- model_fit$finalModel
        if (!is.null(rf_model$variable.importance)) {
          importance_vals <- rf_model$variable.importance
        } else {
          log_msg("[WARN] Ranger model doesn't have importance scores")
          return(NULL)
        }
      } else {
        return(NULL)
      }

      saved_plots <- list()

      # ===== PLOT 1: Feature Importance =====
      tryCatch(
        {
          if (model_fit$method == "rf") {
            # Use randomForest package function explicitly
            if (!requireNamespace("randomForest", quietly = TRUE)) {
              log_warning("randomForest package not available", "Plotting")
              return(NULL)
            }
            importance_vals <- randomForest::importance(rf_model)
            if (is.matrix(importance_vals)) {
              importance_vals <- importance_vals[, 1]
            }
          } else if (model_fit$method == "ranger") {
            # Already extracted above - just verify it exists
            if (is.null(importance_vals)) {
              log_warning("Ranger importance values not found", "Plotting")
              return(NULL)
            }
          }

          importance_data <- data.frame(
            Feature = names(importance_vals),
            Importance = as.numeric(importance_vals),
            stringsAsFactors = FALSE
          )

          importance_data <- importance_data[order(-importance_data$Importance), ]

          if (nrow(importance_data) > 20) {
            importance_data <- head(importance_data, 20)
          }

          if (!is.null(dmrs) && nrow(importance_data) <= nrow(dmrs)) {
            for (i in 1:nrow(importance_data)) {
              feat_name <- importance_data$Feature[i]
              dmr_idx <- which(dmrs$DMR_ID == feat_name)
              if (length(dmr_idx) > 0) {
                importance_data$Meth_Diff[i] <- dmrs$meth.diff[dmr_idx[1]]
                importance_data$Direction[i] <- ifelse(dmrs$meth.diff[dmr_idx[1]] > 0, "Hyper", "Hypo")
              }
            }
          }

          if (!"Direction" %in% colnames(importance_data)) {
            importance_data$Direction <- "Unknown"
          }
          importance_data$Direction[is.na(importance_data$Direction)] <- "Unknown"

          p_importance <- ggplot(importance_data, aes(x = reorder(Feature, Importance), y = Importance)) +
            geom_bar(stat = "identity", aes(fill = Direction), alpha = 0.8) +
            coord_flip() +
            scale_fill_manual(values = c("Hyper" = "#d73027", "Hypo" = "#4575b4", "Unknown" = "gray70")) +
            labs(
              title = "Random Forest Feature Importance",
              subtitle = paste("Top DMRs | Region:", region_label),
              x = "DMR Features",
              y = "Mean Decrease in Gini",
              fill = "Methylation"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = 14, face = "bold"),
              axis.text.y = element_text(size = 8)
            )

          importance_file <- file.path(output_dir, "rf_feature_importance")
          save_plot_universal(p_importance, importance_file, width = 10, height = max(6, nrow(importance_data) * 0.3))
          saved_plots$importance <- importance_file
          log_msg("[PLOT] RF feature importance created")
        },
        error = function(e) {
          log_warning(paste("RF importance plot failed:", e$message), "Plotting")
        }
      )

      # ===== PLOT 2: Decision Tree Example =====
      tryCatch(
        {
          if (requireNamespace("rpart", quietly = TRUE) && requireNamespace("rpart.plot", quietly = TRUE)) {
            tree_model <- rpart::rpart(label ~ .,
              data = train_data, method = "class",
              control = rpart::rpart.control(maxdepth = 4, minsplit = 5)
            )

            # Get plot settings
            if (exists("opt", envir = .GlobalEnv)) {
              opt_global <- get("opt", envir = .GlobalEnv)
              plot_formats <- opt_global$plot_format
              plot_dpi <- opt_global$plot_dpi
            } else {
              plot_formats <- "png"
              plot_dpi <- 300
            }

            if (plot_formats == "all") {
              formats_to_save <- c("png", "jpeg", "svg")
            } else {
              formats_to_save <- unlist(strsplit(plot_formats, ","))
              formats_to_save <- trimws(formats_to_save)
            }

            tree_file_base <- file.path(output_dir, "rf_decision_tree_example")

            for (format in formats_to_save) {
              format <- tolower(format)
              tree_file <- paste0(tree_file_base, ".", format)

              if (format == "png") {
                png(tree_file, width = 12, height = 8, units = "in", res = plot_dpi)
              } else if (format == "jpeg") {
                jpeg(tree_file, width = 12, height = 8, units = "in", res = plot_dpi, quality = 95)
              } else if (format == "svg") {
                svg(tree_file, width = 12, height = 8)
              }

              rpart.plot::rpart.plot(tree_model,
                main = paste("Example Decision Tree\n", region_label, "-", n_markers, "markers"),
                type = 3,
                extra = 104,
                fallen.leaves = TRUE,
                box.palette = "RdBu",
                shadow.col = "gray",
                branch.lty = 3,
                split.font = 1,
                split.cex = 0.8
              )

              dev.off()
            }

            saved_plots$tree <- tree_file_base
            log_msg("[PLOT] RF decision tree created")
          }
        },
        error = function(e) {
          log_warning(paste("Decision tree plot failed:", e$message), "Plotting")
        }
      )

      # ===== PLOT 3: Feature Interactions =====
      tryCatch(
        {
          if (n_markers >= 4 && nrow(importance_data) >= 4) {
            top_features <- head(importance_data[order(-importance_data$Importance), ], 4)$Feature

            interaction_matrix <- matrix(0, nrow = length(top_features), ncol = length(top_features))
            rownames(interaction_matrix) <- top_features
            colnames(interaction_matrix) <- top_features

            for (i in 1:length(top_features)) {
              for (j in 1:length(top_features)) {
                if (i != j && top_features[i] %in% colnames(train_data) && top_features[j] %in% colnames(train_data)) {
                  feat1 <- train_data[[top_features[i]]]
                  feat2 <- train_data[[top_features[j]]]
                  if (is.numeric(feat1) && is.numeric(feat2)) {
                    interaction_matrix[i, j] <- abs(cor(feat1, feat2, use = "complete.obs"))
                  }
                }
              }
            }

            # Get plot settings
            if (exists("opt", envir = .GlobalEnv)) {
              opt_global <- get("opt", envir = .GlobalEnv)
              plot_formats <- opt_global$plot_format
              plot_dpi <- opt_global$plot_dpi
            } else {
              plot_formats <- "png"
              plot_dpi <- 300
            }

            if (plot_formats == "all") {
              formats_to_save <- c("png", "jpeg", "svg")
            } else {
              formats_to_save <- unlist(strsplit(plot_formats, ","))
              formats_to_save <- trimws(formats_to_save)
            }

            interaction_file_base <- file.path(output_dir, "rf_feature_interactions")

            for (format in formats_to_save) {
              format <- tolower(format)

              if (format == "png") {
                interaction_file <- paste0(interaction_file_base, ".png")
                pheatmap::pheatmap(
                  interaction_matrix,
                  main = paste("Top Feature Interactions\n", region_label),
                  color = colorRampPalette(c("white", "red"))(100),
                  display_numbers = TRUE,
                  number_format = "%.2f",
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  filename = interaction_file,
                  width = 8,
                  height = 6
                )
              } else if (format == "jpeg") {
                interaction_file <- paste0(interaction_file_base, ".jpeg")
                jpeg(interaction_file, width = 8, height = 6, units = "in", res = plot_dpi, quality = 95)
                print(pheatmap::pheatmap(
                  interaction_matrix,
                  main = paste("Top Feature Interactions\n", region_label),
                  color = colorRampPalette(c("white", "red"))(100),
                  display_numbers = TRUE,
                  number_format = "%.2f",
                  cluster_rows = TRUE,
                  cluster_cols = TRUE
                ))
                dev.off()
              } else if (format == "svg") {
                interaction_file <- paste0(interaction_file_base, ".svg")
                svg(interaction_file, width = 8, height = 6)
                print(pheatmap::pheatmap(
                  interaction_matrix,
                  main = paste("Top Feature Interactions\n", region_label),
                  color = colorRampPalette(c("white", "red"))(100),
                  display_numbers = TRUE,
                  number_format = "%.2f",
                  cluster_rows = TRUE,
                  cluster_cols = TRUE
                ))
                dev.off()
              }
            }

            saved_plots$interactions <- interaction_file_base
            log_msg("[PLOT] RF feature interactions created")
          }
        },
        error = function(e) {
          log_warning(paste("Interaction plot failed:", e$message), "Plotting")
        }
      )

      # ===== PLOT 4: Partial Dependence =====
      tryCatch(
        {
          if (nrow(importance_data) > 0) {
            top_feature <- importance_data[which.max(importance_data$Importance), "Feature"]

            if (top_feature %in% colnames(train_data)) {
              feature_values <- train_data[[top_feature]]
              feature_range <- seq(min(feature_values, na.rm = TRUE),
                max(feature_values, na.rm = TRUE),
                length.out = 20
              )

              partial_pred <- data.frame()
              for (val in feature_range) {
                temp_data <- train_data
                temp_data[[top_feature]] <- val

                pred_probs <- predict(rf_model, temp_data, type = "prob")
                avg_prob <- colMeans(pred_probs, na.rm = TRUE)

                partial_pred <- rbind(partial_pred, data.frame(
                  Feature_Value = val,
                  Class = names(avg_prob),
                  Probability = as.numeric(avg_prob),
                  stringsAsFactors = FALSE
                ))
              }

              p_partial <- ggplot(partial_pred, aes(x = Feature_Value, y = Probability, color = Class)) +
                geom_line(linewidth = 1.2) +
                geom_point(size = 2) +
                scale_color_brewer(type = "qual", palette = "Set1") +
                labs(
                  title = paste("Partial Dependence Plot"),
                  subtitle = paste("Most Important Feature:", top_feature),
                  x = paste("Methylation Level (%) -", top_feature),
                  y = "Predicted Probability",
                  color = "Class"
                ) +
                theme_minimal() +
                theme(
                  plot.title = element_text(size = 14, face = "bold"),
                  legend.position = "bottom"
                )

              partial_file <- file.path(output_dir, "rf_partial_dependence")
              save_plot_universal(p_partial, partial_file, width = 10, height = 6)
              saved_plots$partial <- partial_file
              log_msg("[PLOT] RF partial dependence created")
            }
          }
        },
        error = function(e) {
          log_warning(paste("Partial dependence plot failed:", e$message), "Plotting")
        }
      )

      log_msg(paste("[PLOT] Created", length(saved_plots), "Random Forest visualizations"))
      return(saved_plots)
    },
    error = function(e) {
      log_warning(paste("Random Forest visualization suite failed:", e$message), "Plotting")
      return(NULL)
    }
  )
}

plot_methylation_heatmap <- function(methylation_matrix, treatment_groups, dmrs,
                                     output_dir, region_label, n_markers) {
  tryCatch(
    {
      # ===== GET GROUP ORDER AND COLORS =====
      group_names_ordered <- NULL
      group_colors <- NULL

      # Priority 1: Check global environment for treatment_processed
      if (exists("treatment_processed", envir = .GlobalEnv)) {
        tp <- get("treatment_processed", envir = .GlobalEnv)
        if (!is.null(tp) && !is.null(tp$group_names)) {
          group_names_ordered <- tp$group_names
          group_colors <- tp$group_colors
        }
      }

      # Priority 2: Check global environment for opt
      if (is.null(group_names_ordered) && exists("opt", envir = .GlobalEnv)) {
        opt_global <- get("opt", envir = .GlobalEnv)
        if (!is.null(opt_global$group_names)) {
          group_names_ordered <- trimws(strsplit(opt_global$group_names, ",")[[1]])
          if (!is.null(opt_global$group_colors)) {
            group_colors <- trimws(strsplit(opt_global$group_colors, ",")[[1]])
            names(group_colors) <- group_names_ordered
          }
        }
      }

      # Priority 3: Fallback to factor levels
      if (is.null(group_names_ordered)) {
        group_names_ordered <- levels(treatment_groups)
      }

      # Set default colors if needed
      if (is.null(group_colors)) {
        n_groups <- length(group_names_ordered)
        group_colors <- brewer.pal(max(3, n_groups), "Set2")[1:n_groups]
        names(group_colors) <- group_names_ordered
      }

      log_msg(paste("[PLOT] Heatmap order:", paste(group_names_ordered, collapse = " -> ")))

      # CRITICAL: Ensure treatment_groups is correctly ordered
      treatment_groups <- factor(as.character(treatment_groups), levels = group_names_ordered)

      # DEBUG: Heatmap group matching diagnostics
      log_msg(paste("[HEATMAP DEBUG] treatment_groups levels:", paste(levels(treatment_groups), collapse = ", ")))
      log_msg(paste("[HEATMAP DEBUG] group_names_ordered:", paste(group_names_ordered, collapse = ", ")))
      log_msg(paste("[HEATMAP DEBUG] treatment_groups length:", length(treatment_groups)))
      log_msg(paste("[HEATMAP DEBUG] methylation_matrix rows:", nrow(methylation_matrix)))

      # Validate inputs
      if (is.null(methylation_matrix) || nrow(methylation_matrix) < 2) {
        log_warning("Insufficient samples for heatmap", "Plotting")
        return(NULL)
      }

      if (ncol(methylation_matrix) < 2) {
        log_warning("Insufficient markers for heatmap", "Plotting")
        return(NULL)
      }

      if (length(treatment_groups) != nrow(methylation_matrix)) {
        log_warning(paste(
          "Treatment length mismatch:", length(treatment_groups),
          "vs matrix rows:", nrow(methylation_matrix)
        ), "Plotting")
        return(NULL)
      }

      # Ensure rownames exist
      if (is.null(rownames(methylation_matrix))) {
        rownames(methylation_matrix) <- paste0("Sample_", 1:nrow(methylation_matrix))
      }

      # Prepare annotation for samples with ORDERED factor
      # CRITICAL: Ensure treatment_groups matches methylation_matrix dimensions
      if (length(treatment_groups) != nrow(methylation_matrix)) {
        log_warning(
          paste(
            "Treatment groups length mismatch:",
            length(treatment_groups), "vs", nrow(methylation_matrix), "samples"
          ),
          "Heatmap annotation"
        )
        # Truncate or pad as needed
        if (length(treatment_groups) > nrow(methylation_matrix)) {
          treatment_groups_matched <- treatment_groups[1:nrow(methylation_matrix)]
        } else {
          treatment_groups_matched <- rep(treatment_groups, length.out = nrow(methylation_matrix))
        }
      } else {
        treatment_groups_matched <- treatment_groups
      }

      # Create annotation with matched groups
      annotation_col <- data.frame(
        Group = factor(as.character(treatment_groups_matched), levels = group_names_ordered),
        row.names = rownames(methylation_matrix),
        stringsAsFactors = FALSE
      )

      # DEBUG: Verify annotation_col was created correctly
      log_msg(paste("[DEBUG] annotation_col created with", nrow(annotation_col), "rows"))
      log_msg(paste("[DEBUG] Group levels:", paste(levels(annotation_col$Group), collapse = ", ")))
      log_msg(paste("[DEBUG] First 3 groups:", paste(head(as.character(annotation_col$Group), 3), collapse = ", ")))

      # UPDATED: Get colors
      n_groups <- length(group_names_ordered)

      if (exists("treatment_processed") && !is.null(treatment_processed) &&
        !is.null(treatment_processed$group_colors)) {
        # Use custom colors
        group_colors <- treatment_processed$group_colors
        log_msg(paste("[PLOT] Heatmap using custom colors"))
      } else {
        # Default colors
        if (n_groups <= 8) {
          group_colors <- brewer.pal(max(3, n_groups), "Set2")[1:n_groups]
        } else {
          group_colors <- rainbow(n_groups)
        }
        names(group_colors) <- group_names_ordered
      }

      annotation_colors <- list(Group = group_colors)

      # Prepare row annotations (DMR info) - CORRECTED TO MATCH BAR PLOT
      if (!is.null(dmrs) && nrow(dmrs) >= ncol(methylation_matrix)) {
        # Get the corrected direction from dmrs (which should have been fixed earlier)
        if ("direction" %in% colnames(dmrs)) {
          # Use the already-corrected direction column
          dmr_directions <- dmrs$direction[1:ncol(methylation_matrix)]

          # Simplify the labels for heatmap (remove the "_in_ClassName" part)
          simplified_directions <- ifelse(grepl("Hypermethylated", dmr_directions),
            "Hyper", "Hypo"
          )

          annotation_row <- data.frame(
            Direction = simplified_directions,
            row.names = colnames(methylation_matrix)
          )

          log_msg(paste("[HEATMAP] Using corrected DMR directions from dmrs$direction"))
          log_msg(paste(
            "          Hyper:", sum(simplified_directions == "Hyper"),
            "| Hypo:", sum(simplified_directions == "Hypo")
          ))
        } else {
          # Fallback: use raw meth.diff (but warn user)
          log_warning("dmrs$direction column not found - using raw meth.diff", "Heatmap")
          annotation_row <- data.frame(
            Direction = ifelse(dmrs$meth.diff[1:ncol(methylation_matrix)] > 0,
              "Hyper", "Hypo"
            ),
            row.names = colnames(methylation_matrix)
          )
        }

        annotation_colors$Direction <- c("Hyper" = "#d73027", "Hypo" = "#4575b4")
      } else {
        annotation_row <- NULL
      }

      # Transpose matrix: samples as rows, DMRs as columns
      heatmap_matrix <- t(methylation_matrix)

      # Handle missing values
      if (any(is.na(heatmap_matrix))) {
        median_val <- median(heatmap_matrix, na.rm = TRUE)
        if (is.na(median_val)) median_val <- 50
        heatmap_matrix[is.na(heatmap_matrix)] <- median_val
      }

      # Create heatmap file path (base name without extension)
      plot_file_base <- file.path(output_dir, paste0("methylation_heatmap_", n_markers, "markers"))

      # Get plot settings from global opt
      if (exists("opt", envir = .GlobalEnv)) {
        opt_global <- get("opt", envir = .GlobalEnv)
        plot_formats <- opt_global$plot_format
        plot_dpi <- opt_global$plot_dpi
      } else {
        plot_formats <- "png"
        plot_dpi <- 300
      }

      # Parse formats
      if (plot_formats == "all") {
        formats_to_save <- c("png", "jpeg", "svg")
      } else {
        formats_to_save <- unlist(strsplit(plot_formats, ","))
        formats_to_save <- trimws(formats_to_save)
      }

      saved_files <- list()

      # Calculate dimensions
      plot_width <- max(8, ncol(heatmap_matrix) * 0.3)
      plot_height <- max(6, nrow(heatmap_matrix) * 0.2)

      # Generate heatmap for each format
      for (format in formats_to_save) {
        format <- tolower(format)

        tryCatch(
          {
            if (format == "png") {
              plot_file <- paste0(plot_file_base, ".png")
              pheatmap(
                heatmap_matrix,
                annotation_col = annotation_col,
                annotation_row = annotation_row,
                annotation_colors = annotation_colors,
                clustering_distance_rows = "euclidean",
                clustering_distance_cols = "euclidean",
                clustering_method = "ward.D2",
                color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
                scale = "row",
                show_rownames = TRUE,
                show_colnames = FALSE,
                fontsize = 10,
                fontsize_row = 8,
                fontsize_col = 8,
                main = paste("Methylation Heatmap -", region_label, "\n", n_markers, "Top DMRs"),
                filename = plot_file,
                width = plot_width,
                height = plot_height
              )
              saved_files$png <- plot_file
              verbose_msg(paste("  Saved heatmap PNG", paste0(plot_dpi, "dpi:"), basename(plot_file)))
            } else if (format == "jpeg" || format == "jpg") {
              plot_file <- paste0(plot_file_base, ".jpeg")
              jpeg(plot_file, width = plot_width, height = plot_height, units = "in", res = plot_dpi, quality = 95)
              print(pheatmap(
                heatmap_matrix,
                annotation_col = annotation_col,
                annotation_row = annotation_row,
                annotation_colors = annotation_colors,
                clustering_distance_rows = "euclidean",
                clustering_distance_cols = "euclidean",
                clustering_method = "ward.D2",
                color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
                scale = "row",
                show_rownames = TRUE,
                show_colnames = FALSE,
                fontsize = 10,
                fontsize_row = 8,
                fontsize_col = 8,
                main = paste("Methylation Heatmap -", region_label, "\n", n_markers, "Top DMRs")
              ))
              dev.off()
              saved_files$jpeg <- plot_file
              verbose_msg(paste("  Saved heatmap JPEG", paste0(plot_dpi, "dpi:"), basename(plot_file)))
            } else if (format == "svg") {
              plot_file <- paste0(plot_file_base, ".svg")
              svg(plot_file, width = plot_width, height = plot_height)
              print(pheatmap(
                heatmap_matrix,
                annotation_col = annotation_col,
                annotation_row = annotation_row,
                annotation_colors = annotation_colors,
                clustering_distance_rows = "euclidean",
                clustering_distance_cols = "euclidean",
                clustering_method = "ward.D2",
                color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
                scale = "row",
                show_rownames = TRUE,
                show_colnames = FALSE,
                fontsize = 10,
                fontsize_row = 8,
                fontsize_col = 8,
                main = paste("Methylation Heatmap -", region_label, "\n", n_markers, "Top DMRs")
              ))
              dev.off()
              saved_files$svg <- plot_file
              verbose_msg(paste("  Saved heatmap SVG (vector):", basename(plot_file)))
            }
          },
          error = function(e) {
            log_warning(paste("Failed to save heatmap as", format, ":", e$message), "Plotting")
          }
        )
      }

      # ===== NEW: EXPORT HEATMAP DATA TO CSV (FINAL FIX) =====
      tryCatch(
        {
          # Transpose back: samples as rows, DMRs as columns
          heatmap_export_data <- as.data.frame(t(heatmap_matrix))

          # Add sample metadata
          heatmap_export_data$Sample_ID <- rownames(heatmap_export_data)

          # DEBUG: Check dimensions and rownames
          log_msg(paste("[DEBUG] heatmap_export_data rows:", nrow(heatmap_export_data)))
          log_msg(paste("[DEBUG] annotation_col rows:", nrow(annotation_col)))
          log_msg(paste(
            "[DEBUG] heatmap_export_data first 3 rownames:",
            paste(head(rownames(heatmap_export_data), 3), collapse = ", ")
          ))
          log_msg(paste(
            "[DEBUG] annotation_col first 3 rownames:",
            paste(head(rownames(annotation_col), 3), collapse = ", ")
          ))
          log_msg(paste("[DEBUG] annotation_col$Group class:", class(annotation_col$Group)))
          log_msg(paste(
            "[DEBUG] annotation_col$Group first 3 values:",
            paste(head(annotation_col$Group, 3), collapse = ", ")
          ))

          # Find common rownames
          common_samples <- intersect(rownames(heatmap_export_data), rownames(annotation_col))
          log_msg(paste("[DEBUG] Common samples:", length(common_samples)))

          if (length(common_samples) > 0 && length(common_samples) == nrow(heatmap_export_data)) {
            # CRITICAL FIX: Convert factor to character FIRST, then subset
            group_vector <- as.character(annotation_col$Group)
            names(group_vector) <- rownames(annotation_col)

            # Now subset by rownames
            heatmap_export_data$Group <- group_vector[rownames(heatmap_export_data)]

            log_msg(paste(
              "[DEBUG] After assignment - first 3 groups:",
              paste(head(heatmap_export_data$Group, 3), collapse = ", ")
            ))
          } else if (nrow(heatmap_export_data) == nrow(annotation_col)) {
            # Fallback: match by position (assumes same order)
            log_msg("[DEBUG] Using position-based matching")
            heatmap_export_data$Group <- as.character(annotation_col$Group)
          } else {
            log_warning("Cannot assign groups - dimension mismatch", "Data export")
            heatmap_export_data$Group <- "Unknown"
          }

          # DEBUG: Check what we got
          log_msg(paste(
            "[DEBUG] Export - unique groups:",
            paste(unique(heatmap_export_data$Group), collapse = ", ")
          ))
          log_msg(paste("[DEBUG] Export - any NA?", any(is.na(heatmap_export_data$Group))))

          # Reorder columns: metadata first, then DMRs
          metadata_cols <- c("Sample_ID", "Group")
          dmr_cols <- setdiff(colnames(heatmap_export_data), metadata_cols)
          heatmap_export_data <- heatmap_export_data[, c(metadata_cols, dmr_cols)]

          # Save to CSV
          heatmap_data_file <- paste0(tools::file_path_sans_ext(plot_file_base), "_DATA.csv")
          write.csv(heatmap_export_data, heatmap_data_file, row.names = FALSE)

          log_msg(paste("[DATA] Heatmap data exported:", basename(heatmap_data_file)))
          log_msg(paste(
            "       Dimensions:", nrow(heatmap_export_data), "samples x",
            length(dmr_cols), "DMRs"
          ))

          # Show group distribution
          group_counts <- table(heatmap_export_data$Group, useNA = "always")
          log_msg(paste("       Groups:", paste(names(group_counts), "=",
            group_counts,
            collapse = ", "
          )))

          saved_files$data <- heatmap_data_file
        },
        error = function(e) {
          log_warning(paste("Failed to export heatmap data:", e$message), "Data export")
        }
      )
      # ===== END HEATMAP DATA EXPORT =====

      log_msg(paste("[PLOT] Heatmap saved in", length(saved_files), "format(s)"))
      return(saved_files)
    },
    error = function(e) {
      log_error(paste("Heatmap failed:", e$message), "Plotting")
      log_error(paste("Error class:", class(e)), "Plotting")
      return(NULL)
    }
  )
}

# Make sure this function exists in your script:

regionStats <- function(methylRawList_obj, regions_gr, column = "perc.meth") {
  regions_df <- data.frame(
    chr = as.character(seqnames(regions_gr)),
    start = start(regions_gr),
    end = end(regions_gr),
    stringsAsFactors = FALSE
  )

  n_samples <- length(methylRawList_obj)
  n_regions <- nrow(regions_df)
  result_matrix <- matrix(NA, nrow = n_samples, ncol = n_regions)

  for (i in 1:n_samples) {
    sample_data <- methylRawList_obj[[i]]
    for (j in 1:n_regions) {
      region <- regions_df[j, ]
      region_cpgs <- sample_data[
        sample_data$chr == region$chr &
          sample_data$start >= region$start &
          sample_data$end <= region$end,
      ]

      if (nrow(region_cpgs) > 0) {
        if (column == "perc.meth") {
          total_cs <- sum(region_cpgs$numCs, na.rm = TRUE)
          total_coverage <- sum(region_cpgs$coverage, na.rm = TRUE)
          if (total_coverage > 0) {
            result_matrix[i, j] <- (total_cs / total_coverage) * 100
          }
        }
      }
    }
  }

  rownames(result_matrix) <- paste0("sample_", 1:n_samples)
  colnames(result_matrix) <- paste0("region_", 1:n_regions)
  return(result_matrix)
}

# ---------------------------
# WINDOW GENERATION FUNCTIONS - NOT INCLUDED
# ---------------------------
# This public version ONLY supports traditional BED file input.
# Users must provide their own genomic regions via --window_file parameter.
# Use GenomeToWindows (github.com/markusdrag/GenomeToWindows) to generate regions.
# ---------------------------

# Stub functions that produce clear error messages
extract_cpg_positions <- function(meth_list) {
  stop(
    "ERROR: This function is not available in the public version.\n",
    "       Please provide genomic regions via --window_file or --region_dir.\n",
    "       Use GenomeToWindows (github.com/markusdrag/GenomeToWindows) to generate windowed regions."
  )
}

create_tiled_windows <- function(cpg_positions, window_size, min_cpgs) {
  stop(
    "ERROR: This function is not available in the public version.\n",
    "       Please provide genomic regions via --window_file or --region_dir.\n",
    "       Use GenomeToWindows (github.com/markusdrag/GenomeToWindows) to generate windowed regions."
  )
}

create_sliding_windows <- function(cpg_positions, window_size, min_cpgs, step_size = NULL) {
  stop(
    "ERROR: This function is not available in the public version.\n",
    "       Please provide genomic regions via --window_file or --region_dir.\n",
    "       Use GenomeToWindows (github.com/markusdrag/GenomeToWindows) to generate windowed regions."
  )
}

create_adaptive_windows <- function(cpg_positions, target_size, min_cpgs, output_dir = ".") {
  stop(
    "ERROR: This function is not available in the public version.\n",
    "       Please provide genomic regions via --window_file or --region_dir.\n",
    "       Use GenomeToWindows (github.com/markusdrag/GenomeToWindows) to generate windowed regions."
  )
}

generate_windows_from_methylation <- function(meth_list, window_sizes, strategy, min_cpgs, output_dir = ".") {
  stop(
    "ERROR: Automated window generation is not available in the public version.\n",
    "       Please provide genomic regions via --window_file or --region_dir.\n",
    "\n",
    "       Tools for generating BED files:\n",
    "         • GenomeToWindows: github.com/markusdrag/GenomeToWindows\n",
    "         • NanoporeToBED-Pipeline: github.com/markusdrag/NanoporeToBED-Pipeline"
  )
}


# Function to create detailed cache metadata
create_cache_metadata <- function(meth, gr, min_coverage, cache_key, rc_object = NULL) {
  metadata <- list(
    cache_key = cache_key,
    created_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    script_version = SCRIPT_VERSION,

    # Methylation data info
    n_samples = length(meth),
    treatment_groups = as.numeric(getTreatment(meth)),

    # Regions info
    n_regions = length(gr),
    unique_scaffolds = length(unique(as.character(seqnames(gr)))),
    total_bp_covered = sum(width(gr)),
    min_region_size = if (length(gr) > 0) min(width(gr)) else 0,
    max_region_size = if (length(gr) > 0) max(width(gr)) else 0,

    # Analysis parameters
    min_coverage = min_coverage,
    strand_aware = TRUE, # Fixed parameter in your script

    # System info
    r_version = R.version.string,
    methylkit_version = packageVersion("methylKit")
  )

  # Add post-processing info if rc_object is provided
  if (!is.null(rc_object)) {
    metadata$result_info <- list(
      n_result_samples = length(rc_object),
      total_regions_after_processing = sum(sapply(rc_object, nrow)),
      sample_region_counts = sapply(rc_object, nrow)
    )
  }

  return(metadata)
}

# Function to save regionCounts results with metadata
save_regioncounts_cache <- function(rc_object, meth, gr, min_coverage, cache_key, cache_dir) {
  log_analysis_step("Saving regionCounts to cache", cache_key)

  # Create cache files
  cache_file_rc <- file.path(cache_dir, paste0("regionCounts_", cache_key, ".rds"))
  cache_file_metadata <- file.path(cache_dir, paste0("regionCounts_", cache_key, "_metadata.json"))

  # Save the regionCounts object
  tryCatch(
    {
      saveRDS(rc_object, cache_file_rc)
      log_file_operation(
        "Saved regionCounts cache", cache_file_rc,
        paste(length(rc_object), "samples")
      )

      # Create and save metadata
      metadata <- create_cache_metadata(meth, gr, min_coverage, cache_key, rc_object)

      # Save metadata as JSON-like text file (since we don't have jsonlite)
      metadata_text <- capture.output({
        cat("{\n")
        for (i in seq_along(metadata)) {
          name <- names(metadata)[i]
          value <- metadata[[i]]

          if (is.list(value)) {
            cat(paste0('  "', name, '": {\n'))
            for (j in seq_along(value)) {
              subname <- names(value)[j]
              subvalue <- value[[j]]
              comma <- if (j < length(value)) "," else ""
              if (is.character(subvalue)) {
                cat(paste0('    "', subname, '": "', subvalue, '"', comma, "\n"))
              } else {
                cat(paste0('    "', subname, '": ', subvalue, comma, "\n"))
              }
            }
            comma <- if (i < length(metadata)) "," else ""
            cat(paste0("  }", comma, "\n"))
          } else {
            comma <- if (i < length(metadata)) "," else ""
            if (is.character(value)) {
              cat(paste0('  "', name, '": "', value, '"', comma, "\n"))
            } else if (is.numeric(value) && length(value) > 1) {
              cat(paste0('  "', name, '": [', paste(value, collapse = ", "), "]", comma, "\n"))
            } else {
              cat(paste0('  "', name, '": ', value, comma, "\n"))
            }
          }
        }
        cat("}\n")
      })

      writeLines(metadata_text, cache_file_metadata)
      log_file_operation("Saved cache metadata", cache_file_metadata)

      return(list(
        cache_file = cache_file_rc,
        metadata_file = cache_file_metadata,
        success = TRUE
      ))
    },
    error = function(e) {
      log_error(paste("Failed to save regionCounts cache:", e$message), "Cache saving")
      return(list(success = FALSE, error = e$message))
    }
  )
}

# Function to load regionCounts from cache with validation
load_regioncounts_cache <- function(cache_key, cache_dir, meth, gr, min_coverage) {
  log_analysis_step("Attempting to load regionCounts from cache", cache_key)

  cache_file_rc <- file.path(cache_dir, paste0("regionCounts_", cache_key, ".rds"))
  cache_file_metadata <- file.path(cache_dir, paste0("regionCounts_", cache_key, "_metadata.json"))

  # Check if cache files exist
  if (!file.exists(cache_file_rc)) {
    verbose_msg("Cache file not found")
    return(NULL)
  }

  if (!file.exists(cache_file_metadata)) {
    log_warning("Cache metadata not found - cache may be incomplete", "Cache validation")
    return(NULL)
  }

  # Validate cache against current parameters
  tryCatch(
    {
      # Read and parse metadata (simple parsing since we don't have jsonlite)
      metadata_lines <- readLines(cache_file_metadata)

      # Extract key validation parameters using simple regex
      cached_n_samples <- as.numeric(gsub(
        '.*"n_samples": ([0-9]+).*', "\\1",
        grep('"n_samples":', metadata_lines, value = TRUE)
      ))
      cached_n_regions <- as.numeric(gsub(
        '.*"n_regions": ([0-9]+).*', "\\1",
        grep('"n_regions":', metadata_lines, value = TRUE)
      ))
      cached_min_coverage <- as.numeric(gsub(
        '.*"min_coverage": ([0-9]+).*', "\\1",
        grep('"min_coverage":', metadata_lines, value = TRUE)
      ))

      # Validate parameters
      current_n_samples <- length(meth)
      current_n_regions <- length(gr)
      current_min_coverage <- min_coverage

      if (length(cached_n_samples) == 0 || cached_n_samples != current_n_samples) {
        log_warning("Cache invalid: sample count mismatch", "Cache validation")
        return(NULL)
      }

      if (length(cached_n_regions) == 0 || cached_n_regions != current_n_regions) {
        log_warning("Cache invalid: region count mismatch", "Cache validation")
        return(NULL)
      }

      if (length(cached_min_coverage) == 0 || cached_min_coverage != current_min_coverage) {
        log_warning("Cache invalid: coverage threshold mismatch", "Cache validation")
        return(NULL)
      }

      # Load the cached object
      log_analysis_step("Loading validated regionCounts cache")
      rc_cached <- readRDS(cache_file_rc)

      # Basic validation of loaded object
      if (!is.list(rc_cached) || length(rc_cached) != current_n_samples) {
        log_warning("Cache invalid: loaded object structure mismatch", "Cache validation")
        return(NULL)
      }

      log_analysis_step(
        "Successfully loaded regionCounts from cache",
        paste(length(rc_cached), "samples")
      )

      return(rc_cached)
    },
    error = function(e) {
      log_warning(paste("Failed to load/validate cache:", e$message), "Cache loading")
      return(NULL)
    }
  )
}

# Function to clean old cache files (optional maintenance)
clean_old_regioncounts_cache <- function(cache_dir, max_age_days = 30) {
  log_analysis_step(
    "Cleaning old regionCounts cache files",
    paste("older than", max_age_days, "days")
  )

  if (!dir.exists(cache_dir)) {
    return()
  }

  cache_files <- list.files(cache_dir, pattern = "^regionCounts_.*\\.rds$", full.names = TRUE)

  if (length(cache_files) == 0) {
    verbose_msg("No cache files found to clean")
    return()
  }

  current_time <- Sys.time()
  cleaned_count <- 0

  for (cache_file in cache_files) {
    file_age <- as.numeric(current_time - file.info(cache_file)$mtime, units = "days")

    if (file_age > max_age_days) {
      # Also remove corresponding metadata file
      metadata_file <- gsub("\\.rds$", "_metadata.json", cache_file)

      unlink(cache_file)
      if (file.exists(metadata_file)) unlink(metadata_file)

      cleaned_count <- cleaned_count + 1
      verbose_msg(paste("Removed old cache file:", basename(cache_file)))
    }
  }

  if (cleaned_count > 0) {
    log_analysis_step("Cache cleanup complete", paste("removed", cleaned_count, "old files"))
  } else {
    verbose_msg("No old cache files to remove")
  }
}

# Main function: Enhanced regionCounts with caching
regionCounts_with_cache <- function(meth, gr, min_coverage, cache_dir) {
  log_analysis_step("Starting regionCounts with caching system")

  # Create cache directory if it doesn't exist
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  # Generate cache key
  cache_key <- generate_regioncounts_cache_key(meth, gr, min_coverage)

  # Try to load from cache first
  rc_cached <- load_regioncounts_cache(cache_key, cache_dir, meth, gr, min_coverage)

  if (!is.null(rc_cached)) {
    log_analysis_step("Using cached regionCounts", "cache hit - saved 6-7 minutes!")
    return(rc_cached)
  }

  # Cache miss - compute regionCounts
  log_analysis_step(
    "Cache miss - computing regionCounts",
    "this will take 6-7 minutes but will be cached for future runs"
  )

  rc_start <- Sys.time()

  # Perform the actual regionCounts computation
  rc <- regionCounts(object = meth, regions = gr, cov.bases = min_coverage, strand.aware = TRUE)

  rc_end <- Sys.time()
  rc_time <- round(as.numeric(rc_end - rc_start, units = "mins"), 1)

  log_analysis_step(
    "regionCounts computation complete",
    paste(length(rc), "samples processed in", rc_time, "minutes")
  )

  # Save to cache for future use
  cache_result <- save_regioncounts_cache(rc, meth, gr, min_coverage, cache_key, cache_dir)

  if (cache_result$success) {
    log_analysis_step("regionCounts cached successfully", "future runs will be much faster")
  } else {
    log_warning("Failed to cache regionCounts - future runs will recompute", "Cache saving")
  }

  return(rc)
}

# ---------------------------
# Validation and setup
# ---------------------------
# ---------------------------
# Enhanced option validation with new train/test file checks
# ---------------------------

# Validate train/test file if provided
predefined_assignments <- NULL
if (!is.null(opt$train_test_file)) {
  log_msg("[TARGET] Predefined train/test split mode enabled")
  predefined_assignments <- load_train_test_assignments(opt$train_test_file)

  # Override test_ratio is ignored when predefined file is used
  if (opt$test_ratio != 0.2) { # 0.2 is the default
    log_msg("[WARN] --test_ratio is ignored when --train_test_file is specified")
  }

  log_msg("[OK] Train/test assignments loaded successfully")
} else {
  log_msg(paste("[RANDOM] Random train/test split mode: test_ratio =", opt$test_ratio))
}

# ============================================================================
# STRATEGY VALIDATION - PUBLIC VERSION
# ============================================================================
# This public version uses the "traditional" approach: analysis of user-provided
# genomic regions (BED format). This requires:
#   --window_file <path/to/regions.bed>  OR
#   --region_dir <path/to/bed_directory>
#
# Tools for generating input files:
#   • GenomeToWindows: Generate windowed genomic regions from reference genomes
#     https://github.com/markusdrag/GenomeToWindows
#
#   • NanoporeToBED-Pipeline: Convert Nanopore modBAM files to BED format
#     https://github.com/markusdrag/NanoporeToBED-Pipeline
# ============================================================================

# Only "traditional" strategy is available
if (opt$strategy != "traditional") {
  stop(
    "\n",
    "═══════════════════════════════════════════════════════════════════════════\n",
    "  ERROR: Invalid Strategy\n",
    "═══════════════════════════════════════════════════════════════════════════\n",
    "\n",
    "  This version of MethylSense only supports the 'traditional' strategy.\n",
    "\n",
    "  Strategy: traditional\n",
    "    Analyzes user-provided genomic regions in BED format.\n",
    "\n",
    "  Required input (choose one):\n",
    "    --window_file <path/to/regions.bed>\n",
    "        Single BED file with genomic regions to analyze\n",
    "\n",
    "    --region_dir <path/to/bed_directory>\n",
    "        Directory containing multiple BED files (one per region set)\n",
    "\n",
    "  Tools for generating BED files:\n",
    "    • GenomeToWindows - Generate windowed genomic regions\n",
    "      https://github.com/markusdrag/GenomeToWindows\n",
    "\n",
    "    • NanoporeToBED-Pipeline - Convert Nanopore modBAM to BED\n",
    "      https://github.com/markusdrag/NanoporeToBED-Pipeline\n",
    "\n",
    "  You specified: --strategy ", opt$strategy, "\n",
    "  Please use: --strategy traditional (or omit --strategy, as it defaults to traditional)\n",
    "═══════════════════════════════════════════════════════════════════════════\n"
  )
}

# Confirm strategy is traditional
opt$strategy <- "traditional"

# Set up parameters based on strategy
# REPLACE THE ENTIRE TRADITIONAL STRATEGY SECTION (around line 400) WITH THIS:

if (opt$strategy == "traditional") {
  # Traditional approach - use pre-made BED files
  region_files <- NULL

  # Check if single window file is specified
  if (!is.null(opt$window_file)) {
    # Use single specified file
    if (!file.exists(opt$window_file)) {
      stop(paste("[ERROR] Window file does not exist:", opt$window_file))
    }

    if (!grepl("\\.bed$", tolower(opt$window_file))) {
      stop(paste("[ERROR] Window file must be a .bed file:", opt$window_file))
    }

    region_files <- opt$window_file
    log_msg(paste("[DIR] Using single window file:", basename(opt$window_file)))
  } else if (!is.null(opt$region_dir)) {
    # Use region directory (existing code)
    if (!dir.exists(opt$region_dir)) {
      stop(paste("[ERROR] Region directory does not exist:", opt$region_dir))
    }

    bed_files <- list.files(opt$region_dir, pattern = "\\.bed$", full.names = TRUE)
    if (length(bed_files) == 0) {
      stop(paste("[ERROR] No .bed files found in:", opt$region_dir))
    }

    region_files <- bed_files
    log_msg(paste("[FOLDER] Using traditional strategy with", length(bed_files), "Window region files from:", opt$region_dir))
  } else {
    # Auto-detect region directory based on common locations
    possible_dirs <- c(
      "BED_files",
      "regions",
      "windows",
      paste0(dirname(opt$qs_file), "/BED_files"),
      paste0(dirname(opt$qs_file), "/regions"),
      paste0(dirname(opt$qs_file), "/windows")
    )

    found_dir <- NULL
    for (dir_path in possible_dirs) {
      if (dir.exists(dir_path) && length(list.files(dir_path, pattern = "\\.bed$")) > 0) {
        found_dir <- dir_path
        break
      }
    }

    if (is.null(found_dir)) {
      stop(paste(
        "[ERROR] Strategy 'traditional' requires --region_dir, --window_file, or BED files in standard locations.",
        "\n   Searched:", paste(possible_dirs, collapse = ", "),
        "\n   Please specify --region_dir /path/to/bed/files/ or --window_file /path/to/file.bed"
      ))
    } else {
      opt$region_dir <- found_dir
      region_files <- list.files(opt$region_dir, pattern = "\\.bed$", full.names = TRUE)
      log_msg(paste("[DIR] Auto-detected BED files in:", found_dir))
    }
  }

  # Validate that we have files to process
  if (is.null(region_files) || length(region_files) == 0) {
    stop("[ERROR] No BED files found to process")
  }

  log_msg(paste("[FOLDER] Found", length(region_files), "region file(s) to process"))
  for (rf in region_files) {
    verbose_msg(paste("Region file:", basename(rf)))
  }

  # Convert BED files to region list
  region_data <- list()
  for (i in seq_along(region_files)) {
    region_file <- region_files[i]
    reg_label <- tools::file_path_sans_ext(basename(region_file))

    # Read BED file
    bed_data <- fread(region_file, header = FALSE)
    colnames(bed_data)[1:3] <- c("chr", "start", "end")

    # Apply species-agnostic filtering if requested
    if (!is.null(opt$filter_standard_chroms) && opt$filter_standard_chroms) {
      # Universal filter: 1-99 + single-letter sex chromosomes
      valid_chroms <- c(
        as.character(1:99), # Covers all vertebrates
        "Z", "W", "z", "w", # Birds
        "X", "Y", "x", "y" # Mammals
      )

      bed_filtered <- bed_data[bed_data$chr %in% valid_chroms, ]
      log_msg(paste("[PACKAGE] Chromosome filtering:", nrow(bed_data), "->", nrow(bed_filtered), "regions"))
    } else {
      # No filtering - use all chromosomes from BED file
      bed_filtered <- bed_data
      log_msg(paste("[PACKAGE] Using all chromosomes:", nrow(bed_data), "regions"))
    }

    # Exclude sex chromosomes if requested (safety feature)
    if (!is.null(opt$exclude_sex_chroms) && opt$exclude_sex_chroms) {
      sex_chroms <- c("X", "Y", "Z", "W", "x", "y", "z", "w")
      bed_before <- bed_filtered
      bed_filtered <- bed_filtered[!(bed_filtered$chr %in% sex_chroms), ]
      log_msg(paste("[SAFETY] Sex chromosome exclusion:", nrow(bed_before), "->", nrow(bed_filtered), "regions"))
    }

    if (nrow(bed_filtered) == 0) {
      log_msg(paste("[WARN] No valid chromosomes found in", basename(region_file), "- skipping"))
      next
    }

    gr <- toGRanges(bed_filtered[, 1:3], format = "BED")
    region_data[[reg_label]] <- list(
      gr = gr,
      label = reg_label,
      file_path = region_file
    )
  }

  # Traditional strategy validated - no additional setup needed
}

# Validate test_ratio
if (opt$test_ratio <= 0 || opt$test_ratio >= 1) {
  stop("[ERROR] --test_ratio must be between 0 and 1 (exclusive)")
}

if (opt$min_accuracy < 0 || opt$min_accuracy > 1) {
  stop("[ERROR] --min_accuracy must be between 0 and 1")
}

train_ratio <- 1 - opt$test_ratio

# NOW continue with directory setup...

# REMOVE ALL OTHER DIRECTORY CREATION CODE AND USE ONLY THIS:

# Clean the base directory path - remove any trailing subdirectories
if (is.null(opt$output_dir)) {
  stop("[ERROR] --output_dir is required")
}
base_dir <- dirname(opt$output_dir)
if (basename(base_dir) == ".") {
  base_dir <- opt$output_dir
}

## REPLACE the entire directory creation section (around line 600-700) with this:

# Extract species name from output_dir FIRST
species_name <- basename(opt$output_dir) # Get "Chicken_" from the path
species_name <- gsub("_$", "", species_name) # Remove trailing underscore if present

# Clean the base directory path
base_dir <- dirname(opt$output_dir) # Get the parent directory
if (basename(base_dir) == ".") {
  base_dir <- opt$output_dir
}

# Build parameter pattern with methylation thresholds AND strategy
param_pattern <- paste0(
  "_", opt$strategy, # Strategy name
  "_q", gsub("\\.", "", sprintf("%.3f", opt$qval)),
  "_cov", opt$min_coverage,
  "_test", gsub("\\.", "", sprintf("%.2f", opt$test_ratio)),
  if (!opt$disable_meth_filter) {
    paste0(
      "_hyper", gsub("\\.", "", sprintf("%.1f", opt$hyper_threshold)),
      "_hypo", gsub("\\.", "", sprintf("%.1f", abs(opt$hypo_threshold)))
    )
  } else {
    "_nomethfilt"
  }
)

# Add strategy-specific info (traditional only)
if (opt$strategy == "traditional") {
  # For traditional, add cpg filter to output pattern
  param_pattern <- paste0(param_pattern, "_cpg", opt$min_cpg)
}

# Get only top-level directories from the clean base path
all_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

# Look for directories that contain our parameter pattern
matching_dirs <- all_dirs[grepl(param_pattern, basename(all_dirs), fixed = TRUE)]

log_msg(paste("Searching in base directory:", base_dir))
log_msg(paste("Parameter pattern:", param_pattern))
log_msg(paste("Species name:", species_name)) # DEBUG LINE
log_msg(paste("Found", length(all_dirs), "total directories"))

if (length(matching_dirs) > 0) {
  log_msg("Found matching directories:")
  for (dir in matching_dirs) {
    log_msg(paste("  ", basename(dir)))
  }

  # Use the most recent one
  dir_times <- file.info(matching_dirs)$mtime
  newest_idx <- which.max(dir_times)
  unique_output_dir <- matching_dirs[newest_idx]

  log_msg(paste("Using existing directory:", basename(unique_output_dir)))
} else {
  log_msg("No matching directories found, creating new one")

  # List what directories ARE available for debugging
  log_msg("Available directories:")
  for (dir in all_dirs) {
    log_msg(paste("  ", basename(dir)))
  }

  # Create new directory with species-specific naming
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
  new_dir_name <- paste0(species_name, param_pattern, "_", timestamp)
  unique_output_dir <- file.path(base_dir, new_dir_name)

  dir.create(unique_output_dir, showWarnings = FALSE, recursive = TRUE)
  log_msg(paste("Created new directory:", basename(unique_output_dir)))
}

# Set the final output directory
opt$output_dir <- unique_output_dir
# ONLY initialize these ONCE
if (!exists("logging_initialized") || !logging_initialized) {
  initialize_logging(unique_output_dir)
  output_dirs <- create_output_structure(unique_output_dir)

  # Save analysis settings for reproducibility
  settings_file <- file.path(unique_output_dir, "analysis_settings.txt")

  # Get original command line arguments
  args <- commandArgs(trailingOnly = TRUE)

  # Reconstruct the Rscript command
  script_name <- sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(script_name) == 0) script_name <- "MethylSense.R"

  rscript_command <- paste0(
    "Rscript ", basename(script_name), " \\\n"
  )

  # Add all arguments
  if (length(args) > 0) {
    arg_lines <- character()
    i <- 1
    while (i <= length(args)) {
      if (startsWith(args[i], "--")) {
        if (i < length(args) && !startsWith(args[i + 1], "--")) {
          # Flag with value
          arg_lines <- c(arg_lines, paste0("  ", args[i], " ", shQuote(args[i + 1], type = "sh")))
          i <- i + 2
        } else {
          # Flag without value (boolean)
          arg_lines <- c(arg_lines, paste0("  ", args[i]))
          i <- i + 1
        }
      } else {
        i <- i + 1
      }
    }
    rscript_command <- paste0(rscript_command, paste(arg_lines, collapse = " \\\n"))
  }

  writeLines(c(
    paste(rep("=", 80), collapse = ""),
    "METHYLSENSE ANALYSIS SETTINGS",
    paste(rep("=", 80), collapse = ""),
    paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste("Script Version:", SCRIPT_VERSION),
    paste("Output Directory:", unique_output_dir),
    "",
    "STRATEGY:",
    paste("  Strategy:", opt$strategy),
    paste("  Region directory:", opt$region_dir),
    "",
    "DATA FILTERING:",
    paste("  Q-value threshold:", opt$qval),
    paste("  Min coverage:", opt$min_coverage),
    paste("  Min CpG per sample:", opt$min_cpg),
    paste("  Min samples per region:", opt$min_samples_per_region),
    paste("  Methylation difference filter:", ifelse(opt$disable_meth_filter, "DISABLED", "ENABLED")),
    if (!opt$disable_meth_filter) {
      c(
        paste("    Hypermethylation threshold:", opt$hyper_threshold, "%"),
        paste("    Hypomethylation threshold:", opt$hypo_threshold, "%")
      )
    } else {
      NULL
    },
    "",
    "TREATMENT GROUPS:",
    if (!is.null(opt$group_names)) paste("  Group names:", opt$group_names) else "  Group names: (auto-detected)",
    if (!is.null(opt$group_colors)) paste("  Group colors:", opt$group_colors) else NULL,
    if (!is.null(opt$positive_class)) paste("  Positive class:", opt$positive_class) else NULL,
    if (!is.null(opt$treatment_mapping)) paste("  Treatment mapping:", opt$treatment_mapping) else NULL,
    paste("  Min group size:", opt$min_group_size),
    "",
    "MACHINE LEARNING:",
    paste("  Models:", opt$models),
    paste("  Marker range:", opt$marker_range),
    paste("  Test ratio:", opt$test_ratio),
    paste("  Min accuracy threshold:", opt$min_accuracy),
    paste("  Cross-validation repeats:", opt$cv_repeats),
    paste("  Nested CV:", opt$nested_cv),
    if (!is.null(opt$train_test_file)) paste("  Predefined split file:", opt$train_test_file) else NULL,
    "",
    "COMPUTATIONAL SETTINGS:",
    paste("  Cores:", opt$cores),
    paste("  Reproducible mode:", opt$reproducible),
    paste("  Force ML rerun:", opt$force_ml),
    "",
    "OUTPUT OPTIONS:",
    paste("  Create plots:", opt$create_plots),
    if (opt$create_plots) paste("  Plot format:", opt$plot_format) else NULL,
    if (opt$create_plots) paste("  Plot DPI:", opt$plot_dpi) else NULL,
    if (opt$create_plots) paste("  Plot top N DMRs:", opt$plot_top_n) else NULL,
    paste("  Perform clustering:", opt$perform_clustering),
    if (opt$perform_clustering && !is.null(opt$n_clusters)) paste("  Number of clusters:", opt$n_clusters) else NULL,
    "",
    paste(rep("=", 80), collapse = ""),
    "COMMAND TO REPRODUCE THIS ANALYSIS:",
    paste(rep("=", 80), collapse = ""),
    "",
    rscript_command,
    "",
    paste(rep("=", 80), collapse = ""),
    "END OF SETTINGS",
    paste(rep("=", 80), collapse = "")
  ), settings_file)

  log_file_operation("Saved analysis settings", settings_file)

  logging_initialized <- TRUE
}

# Print version info with strategy
cat(paste0("[INFO] DMR + MACHINE LEARNING ANALYSIS SCRIPT ", SCRIPT_VERSION, "\n"))
cat(paste0("[DATE] Release Date: ", SCRIPT_DATE, "\n"))
cat("[LOG] Expects: PROPER S4 methylRaw objects from preprocessing v2.0-FINAL\n")
cat("[SETUP] Reproducibility Mode:", ifelse(opt$reproducible, "ENABLED (single-core, deterministic)", "DISABLED (multi-core parallel)"), "\n")
cat("[TARGET] Strategy:", opt$strategy, "\n")

if (opt$strategy == "traditional") {
  cat("[DIR] BED files from:", opt$region_dir, "\n")
} else {
  stop("[ERROR] Invalid strategy. Only 'traditional' supported in public version.")
}

cat("[CONFIG] Cores:", opt$cores, "| Reproducible:", opt$reproducible, "| Models:", opt$models, "| Markers:", opt$marker_range, "\n")
cat("[DATA] Parameters: Q-val <=", opt$qval, "| Min coverage >=", opt$min_coverage, "| Min CpG >=", opt$min_cpg, "\n")

if (!opt$disable_meth_filter) {
  cat("[DNA] Methylation filter: Hyper >=", opt$hyper_threshold, "% | Hypo <=", opt$hypo_threshold, "%\n")
} else {
  cat("[DNA] Methylation filter: DISABLED (q-value only)\n")
}

cat("[TARGET] Train/Test split:", paste0(round(train_ratio * 100), "%/", round(opt$test_ratio * 100), "%\n"))

# NEW: Show group names and colors if specified
if (!is.null(opt$group_names)) {
  cat("[GROUPS] Treatment groups:", paste(trimws(strsplit(opt$group_names, ",")[[1]]), collapse = " -> "), "\n")

  if (!is.null(opt$group_colors)) {
    colors <- trimws(strsplit(opt$group_colors, ",")[[1]])
    names <- trimws(strsplit(opt$group_names, ",")[[1]])
    cat("[COLORS] Group colors:\n")
    for (i in seq_along(names)) {
      cat(sprintf("         %s = %s\n", names[i], colors[i]))
    }
  }
}

if (!is.null(opt$positive_class)) {
  cat("[TARGET] Positive class for ROC:", opt$positive_class, "\n")
}

cat("[PLOT] Plotting:", ifelse(opt$create_plots, paste("Enabled (top", opt$plot_top_n, "DMRs + CV analysis)"), "Disabled"), "\n")
if (opt$create_plots) {
  cat("[PLOT] Format:", opt$plot_format, "| DPI:", opt$plot_dpi, "\n")
}
cat("[TARGET] Accuracy threshold:", opt$min_accuracy, "(models below this won't be saved)\n")

if (opt$cv_repeats > 0) {
  cat("[CV] Cross-validation:", opt$cv_repeats, "repeats", ifelse(opt$nested_cv, "(nested CV enabled)", ""), "\n")
}

if (opt$debug) cat("[DEBUG] DEBUG MODE ENABLED\n")
if (opt$verbose) cat("[VERBOSE] VERBOSE MODE ENABLED\n")

if (opt$reproducible) {
  cat("[LOCKED] STRICT REPRODUCIBILITY: Fixed seeds, single-core, deterministic CV\n")
} else {
  cat("[FAST] PERFORMANCE MODE: Parallel processing enabled, results may vary slightly\n")
}

cat("\n")

# ---------------------------
# Load methylRawList - REVERTED TO USE NORMALIZED DATA CONSISTENTLY
# ---------------------------
log_msg(paste("[SYSTEM] Loading PROPER methylRawList from:", basename(opt$qs_file)))
verbose_msg("Reading .qs file into memory...")

start_load <- Sys.time()
meth_qread <- qread(opt$qs_file)
verbose_msg("[NORMALIZING] Normalising coverage using median method...")
meth <- normalizeCoverage(meth_qread, method = "median")
verbose_msg("[TARGET] Methylation data normalised!")
end_load <- Sys.time()

# =========================================================================
# INTELLIGENT SPECIES DETECTION FROM QS FILE
# =========================================================================
log_msg("[DATA] Detecting species from methylation data...")

# Extract species name from QS filename
qs_basename <- basename(opt$qs_file)
species_from_filename <- gsub("_n\\d+.*", "", gsub("^\\d{8}_", "", qs_basename))
log_msg(paste("[DATA] Species from filename:", species_from_filename))

# Get actual chromosomes/scaffolds present in the methylation data
meth_scaffolds <- unique(unlist(lapply(meth, function(x) {
  if (is(x, "methylRaw")) {
    return(unique(as.character(x$chr)))
  } else {
    return(character(0))
  }
})))

log_msg(paste("[DATA] Found", length(meth_scaffolds), "unique chromosomes/scaffolds in methylation data"))
log_msg(paste("[DATA] Example chromosomes:", paste(head(meth_scaffolds, 10), collapse = ", ")))

# Detect chromosome naming pattern
if (any(grepl("^chr", meth_scaffolds, ignore.case = TRUE))) {
  chr_pattern <- "chr_prefix"
} else if (any(grepl("^[0-9]+$", meth_scaffolds))) {
  chr_pattern <- "numeric"
} else if (any(grepl("^NC_", meth_scaffolds))) {
  chr_pattern <- "ncbi"
} else {
  chr_pattern <- "custom"
}

log_msg(paste("[DATA] Detected chromosome pattern:", chr_pattern))

# Store for later use
DETECTED_SPECIES <- species_from_filename
DETECTED_CHROMOSOMES <- meth_scaffolds
CHROMOSOME_PATTERN <- chr_pattern

log_msg(paste("[DATA] Loaded", length(meth), "samples in", round(as.numeric(end_load - start_load, units = "secs")), "seconds"))

if (length(meth) > 0) {
  first_obj <- meth[[1]]
  verbose_msg(paste("Object class:", paste(class(first_obj), collapse = ", ")))
  verbose_msg(paste("First sample columns:", paste(colnames(first_obj), collapse = ", ")))
  verbose_msg(paste("First sample rows:", nrow(first_obj)))
}

# Filter samples with enough CpGs
sample_sizes <- sapply(meth, nrow)
verbose_msg(paste("Sample sizes range:", min(sample_sizes), "to", max(sample_sizes), "CpGs"))

# valid_samples <- which(sample_sizes >= opt$min_cpg)

# if (length(valid_samples) < length(meth)) {
#  log_msg(paste("[DATA] Filtering:", length(meth) - length(valid_samples), "samples with <", opt$min_cpg, "CpGs"))
#
#  filtered_data <- meth[valid_samples]
#  filtered_treatment <- getTreatment(meth)[valid_samples]
#
#  verbose_msg("Recreating methylRawList with filtered samples...")
#  meth_new <- new("methylRawList", filtered_data)
#  meth_new@treatment <- as.integer(filtered_treatment)
#  meth <- meth_new
# }
log_msg(paste("[OK] Using", length(meth), "samples for analysis"))

# Filter chromosomes if requested
if (!is.null(opt$filter_standard_chroms) && opt$filter_standard_chroms) {
  log_msg("[DNA] Filtering to standard chromosomes (numeric + sex chromosomes)")

  # Universal filter: 1-2 digit numbers + single letter sex chromosomes
  # Matches: 1-99, Z, W, X, Y (and lowercase variants)
  valid_chroms <- c(
    as.character(1:99), # Covers all vertebrates (max ~50 chromosomes typically)
    "Z", "W", "z", "w", # Birds
    "X", "Y", "x", "y" # Mammals
  )

  meth_filtered <- lapply(meth, function(sample_data) {
    sample_data[sample_data$chr %in% valid_chroms, ]
  })

  # Recreate methylRawList
  meth_new <- new("methylRawList", meth_filtered)
  meth_new@treatment <- getTreatment(meth)
  meth <- meth_new

  log_msg(paste("[OK] Retained standard chromosomes only (numeric 1-99 + sex chromosomes)"))
}

# Exclude sex chromosomes if requested (safety feature)
if (!is.null(opt$exclude_sex_chroms) && opt$exclude_sex_chroms) {
  log_msg("[SAFETY] Excluding sex chromosomes (X, Y, Z, W) from analysis")

  sex_chroms <- c("X", "Y", "Z", "W", "x", "y", "z", "w")

  meth_filtered <- lapply(meth, function(sample_data) {
    sample_data[!(sample_data$chr %in% sex_chroms), ]
  })

  # Recreate methylRawList
  meth_new <- new("methylRawList", meth_filtered)
  meth_new@treatment <- getTreatment(meth)
  meth <- meth_new

  log_msg(paste("[OK] Sex chromosomes excluded from marker discovery"))
}

# =========================================================================
# TREATMENT GROUP VALIDATION AND FILTERING
# =========================================================================

# Extract treatment info BEFORE processing treatment groups
treatment_info <- getTreatment(meth)

# Enhanced treatment group processing with custom names/mapping
treatment_processed <- process_treatment_groups(
  treatment_info = treatment_info,
  group_names = group_names_list,
  treatment_mapping = opt$treatment_mapping
)

# Update treatment_info to use the processed labels
treatment_info <- treatment_processed$labels
group_names_final <- treatment_processed$group_names

# Show initial group distribution
treatment_table <- table(treatment_info)
log_msg("[DATA] Initial treatment groups:")
for (i in 1:length(treatment_table)) {
  log_msg(paste("   ", names(treatment_table)[i], ":", treatment_table[i], "samples"))
}

# Determine effective minimum based on flags
if (is.null(opt$min_group_size)) {
  opt$min_group_size <- 4
}

if (opt$allow_small_groups) {
  effective_min <- 2 # Absolute minimum for any ML
  log_msg(paste("[WARN] --allow_small_groups enabled: Using absolute minimum of", effective_min, "samples per group"))
  log_msg("[WARN] This may result in poor train/test splits and overfitting!")
} else {
  effective_min <- opt$min_group_size
  log_msg(paste("[INFO] Minimum group size:", effective_min, "samples (--min_group_size)"))
}

# Filter out groups with fewer than effective minimum
small_groups <- names(treatment_table[treatment_table < effective_min])

if (length(small_groups) > 0) {
  if (opt$allow_small_groups) {
    log_msg(paste(
      "[WARN] Removing", length(small_groups), "group(s) with <", effective_min,
      "samples:", paste(small_groups, collapse = ", ")
    ))
    log_msg("[WARN] Even with --allow_small_groups, groups need at least 2 samples")
  } else {
    log_msg(paste(
      "[FILTER] Removing", length(small_groups), "group(s) with <", effective_min,
      "samples:", paste(small_groups, collapse = ", ")
    ))
    log_msg(paste("[INFO] Reason: Minimum", effective_min, "samples needed for meaningful train/test split"))
    log_msg(paste("[TIP] Use --min_group_size 3 or --allow_small_groups to include smaller groups (risky)"))
  }

  # Create mask for samples to keep
  samples_to_keep <- !treatment_info %in% small_groups

  # Filter methylation data
  meth_filtered <- meth[samples_to_keep]
  treatment_info <- treatment_info[samples_to_keep]
  treatment_info <- droplevels(treatment_info)

  # Update group names
  group_names_final <- levels(treatment_info)

  # Recreate methylRawList
  meth_new <- new("methylRawList", meth_filtered)
  meth_new@treatment <- as.integer(treatment_info)
  meth <- meth_new

  # Update treatment_processed to match filtered data
  treatment_processed$labels <- treatment_info
  treatment_processed$group_names <- group_names_final

  log_msg(paste("[OK] Filtered to", length(meth), "samples in", length(group_names_final), "groups"))

  # Show updated distribution
  treatment_table <- table(treatment_info)
  log_msg("[DATA] Updated treatment groups:")
  for (i in 1:length(treatment_table)) {
    log_msg(paste("   ", names(treatment_table)[i], ":", treatment_table[i], "samples"))
  }
} else {
  log_msg(paste("[OK] All groups have >=", effective_min, "samples"))
}

# Recalculate treatment_table
treatment_table <- table(treatment_info)

# Validate positive class after filtering
if (!is.null(opt$positive_class)) {
  if (!opt$positive_class %in% levels(treatment_info)) {
    log_warning(paste("Positive class", opt$positive_class, "was filtered out or not found."))
    log_msg(paste("Available groups:", paste(levels(treatment_info), collapse = ", ")))

    # Auto-select first group as positive class
    opt$positive_class <- levels(treatment_info)[1]
    log_msg(paste("[PROCESS] Auto-selected positive class:", opt$positive_class))
  }
}

# Final validation
min_group_size <- min(treatment_table)

if (min_group_size < effective_min) {
  stop(paste(
    "[ERROR] After filtering, minimum group size is", min_group_size,
    "but effective minimum requires", effective_min,
    "\n[INFO] Either reduce --min_group_size or provide more balanced data"
  ))
}

log_msg(paste("[OK] Final dataset:", length(meth), "samples across", length(levels(treatment_info)), "groups"))
log_msg(paste("[OK] Minimum group size:", min_group_size, "samples"))

# Check if we have enough samples for the specified test ratio
min_test_samples <- ceiling(min_group_size * opt$test_ratio)
min_train_samples <- min_group_size - min_test_samples

if (min_train_samples < 2 || min_test_samples < 1) {
  # Calculate maximum safe test ratio
  max_safe_ratio <- (min_group_size - 2) / min_group_size
  stop(paste(
    "[ERROR] Test ratio", opt$test_ratio, "leaves insufficient samples for training/testing.",
    "\n   With minimum group size", min_group_size, ":",
    "\n   - Training samples per group:", min_train_samples,
    "\n   - Test samples per group:", min_test_samples,
    "\n[INFO] Solution: Use --test_ratio <=", round(max_safe_ratio, 2),
    "\n[INFO] Or use --min_group_size to increase minimum required samples"
  ))
}

log_msg(paste("[OK] Train/test validation passed"))
log_msg(paste(
  "   Expected split per group: ~", min_train_samples, "train /",
  min_test_samples, "test"
))

# Extract species name from filename or use strategy info
if (opt$strategy == "traditional") {
  species_string <- tolower(gsub(".*\\d{8}_(.*?)_n\\d+.*\\.qs$", "\\1", basename(opt$qs_file)))
  log_msg(paste("[DATA] Species detected:", species_string))
} else {
  species_string <- "generated_windows"
}


# ---------------------------
# Enhanced Model Selection System
# ---------------------------


# Define available models and parse selection
model_specs <- list(
  rf = "rf",
  glmnet = "glmnet",
  knn = "knn",
  lda = "lda",
  svm = "svmRadial",
  xgboost = "xgbTree",
  nb = "nb",
  ranger = "ranger",
  logreg = "glm",
  nnet = "nnet"
)

model_categories <- list(
  fast = c("knn", "lda", "glmnet", "nb", "logreg"),
  medium = c("rf", "ranger", "nnet"),
  slow = c("svm", "xgboost"),
  linear = c("glmnet", "lda", "logreg"),
  nonlinear = c("rf", "svm", "knn", "xgboost", "ranger", "nnet"),
  ensemble = c("rf", "xgboost", "ranger"),
  interpretable = c("glmnet", "lda", "nb", "logreg"),
  high_performance = c("xgboost", "rf", "ranger"),
  biomarker = c("xgboost", "rf", "glmnet", "ranger", "svm"),
  baseline = c("logreg", "nb"),
  # NEW: All models except xgboost (faster execution)
  all_except_xgb = c("rf", "glmnet", "knn", "lda", "svm", "nb", "ranger", "logreg", "nnet")
)

# Simple model parsing function
parse_model_selection <- function(models_string) {
  models_string <- tolower(trimws(models_string))
  requested_items <- strsplit(models_string, ",")[[1]]
  requested_items <- trimws(requested_items)

  selected_models <- character(0)

  for (item in requested_items) {
    if (item == "all") {
      selected_models <- c(selected_models, names(model_specs))
    } else if (item %in% names(model_categories)) {
      selected_models <- c(selected_models, model_categories[[item]])
    } else if (item %in% names(model_specs)) {
      selected_models <- c(selected_models, item)
    } else {
      log_msg(paste("[WARN] Unknown model/category ignored:", item))
    }
  }

  selected_models <- unique(selected_models)
  log_msg(paste("[MODEL] Selected models:", paste(selected_models, collapse = ", ")))
  return(selected_models)
}

# Parse models and setup - DEFINE model_list HERE
model_list <- parse_model_selection(opt$models)

# SECURITY FIX: Safe parsing of marker_range (no eval/parse)
marker_seq <- tryCatch(
  {
    parts <- strsplit(opt$marker_range, ":")[[1]]
    if (length(parts) != 2) stop("Invalid marker_range format")
    start <- as.integer(parts[1])
    end <- as.integer(parts[2])
    if (is.na(start) || is.na(end) || start > end) stop("Invalid marker range values")
    seq(start, end)
  },
  error = function(e) {
    stop(paste0(
      "[ERROR] Invalid --marker_range: ", e$message,
      ". Use format: 'start:end' (e.g., '2:50')"
    ))
  }
)

# NOW we can safely use model_list in verbose messages
verbose_msg(paste("Models to train:", paste(model_list, collapse = ", ")))
verbose_msg(paste("Marker range:", min(marker_seq), "to", max(marker_seq)))

# Initialize results summary data frame
results_summary <- data.frame(
  region = character(0),
  model = character(0),
  markers = integer(0),
  accuracy = numeric(0),
  auc = numeric(0),
  pr_auc = numeric(0),
  mean_sensitivity = numeric(0),
  mean_specificity = numeric(0),
  mean_f1 = numeric(0),
  train_time_sec = numeric(0),
  test_ratio = numeric(0),
  train_samples = integer(0),
  test_samples = integer(0),
  stringsAsFactors = FALSE
)

# Also initialize the total start time for timing
total_start_time <- Sys.time()

# Set up the current model list for use in the main loop
current_model_list <- model_list

verbose_msg(paste("Models to train:", paste(model_list, collapse = ", ")))
verbose_msg(paste("Marker range:", min(marker_seq), "to", max(marker_seq)))

# Generate windows from methylation data OR find existing BED files
if (opt$strategy == "traditional") {
  region_files <- NULL

  if (!is.null(opt$window_file)) {
    # Single file mode
    if (!file.exists(opt$window_file)) {
      stop(paste("[ERROR] Window file does not exist:", opt$window_file))
    }
    region_files <- opt$window_file
    log_msg(paste("[DIR] Using single window file:", basename(opt$window_file)))
  } else if (!is.null(opt$region_dir)) {
    # Directory mode with intelligent filtering
    if (!dir.exists(opt$region_dir)) {
      stop(paste("[ERROR] Region directory does not exist:", opt$region_dir))
    }

    all_bed_files <- list.files(opt$region_dir, pattern = "\\.bed$", full.names = TRUE)

    if (length(all_bed_files) == 0) {
      stop(paste("[ERROR] No .bed files found in:", opt$region_dir))
    }

    log_msg(paste("[FOLDER] Found", length(all_bed_files), "total BED files in directory"))

    # ===== INTELLIGENT SPECIES-BASED FILTERING =====
    if (exists("DETECTED_SPECIES") && !is.null(DETECTED_SPECIES)) {
      # Filter BED files by species name
      species_pattern <- gsub("_", "[_-]?", DETECTED_SPECIES) # Handle _ or - variations
      matching_files <- grep(species_pattern, basename(all_bed_files),
        ignore.case = TRUE, value = FALSE
      )

      if (length(matching_files) > 0) {
        region_files <- all_bed_files[matching_files]
        log_msg(paste(
          "[SMART FILTER] Matched", length(region_files),
          "BED files for species:", DETECTED_SPECIES
        ))
        log_msg(paste(
          "[FAST] Skipping", length(all_bed_files) - length(region_files),
          "BED files from other species"
        ))
      } else {
        log_msg(paste("[WARN] No BED files matched species pattern:", DETECTED_SPECIES))
        log_msg("[FALLBACK] Using all BED files in directory")
        region_files <- all_bed_files
      }
    } else {
      log_msg("[WARN] Species not detected - using all BED files")
      region_files <- all_bed_files
    }

    log_msg(paste("[FOLDER] Processing", length(region_files), "species-matched BED file(s)"))
    for (rf in region_files) {
      verbose_msg(paste("  - ", basename(rf)))
    }
  } else {
    # Auto-detect (existing code)
    possible_dirs <- c(
      "BED_files",
      "regions",
      "windows",
      paste0(dirname(opt$qs_file), "/BED_files"),
      paste0(dirname(opt$qs_file), "/regions"),
      paste0(dirname(opt$qs_file), "/windows")
    )

    found_dir <- NULL
    for (dir_path in possible_dirs) {
      if (dir.exists(dir_path) && length(list.files(dir_path, pattern = "\\.bed$")) > 0) {
        found_dir <- dir_path
        break
      }
    }

    if (is.null(found_dir)) {
      stop(paste(
        "[ERROR] Strategy 'traditional' requires --region_dir, --window_file, or BED files in standard locations.",
        "\n   Searched:", paste(possible_dirs, collapse = ", ")
      ))
    } else {
      opt$region_dir <- found_dir
      all_bed_files <- list.files(opt$region_dir, pattern = "\\.bed$", full.names = TRUE)

      # Apply species filtering even for auto-detected directories
      if (exists("DETECTED_SPECIES") && !is.null(DETECTED_SPECIES)) {
        species_pattern <- gsub("_", "[_-]?", DETECTED_SPECIES)
        matching_files <- grep(species_pattern, basename(all_bed_files),
          ignore.case = TRUE, value = FALSE
        )

        if (length(matching_files) > 0) {
          region_files <- all_bed_files[matching_files]
          log_msg(paste(
            "[AUTO-DETECT] Found", length(region_files),
            "matching BED files in:", found_dir
          ))
        } else {
          region_files <- all_bed_files
          log_msg(paste("[AUTO-DETECT] Using all BED files from:", found_dir))
        }
      } else {
        region_files <- all_bed_files
      }
    }
  }

  # Validate that we have files to process
  if (is.null(region_files) || length(region_files) == 0) {
    stop("[ERROR] No BED files found to process after species filtering")
  }

  # ===== LOAD AND VALIDATE BED FILES =====
  # CORRECTED TO MATCH v44 EXACTLY
  region_data <- list()

  for (i in seq_along(region_files)) {
    region_file <- region_files[i]
    reg_label <- tools::file_path_sans_ext(basename(region_file))

    log_msg(paste("[LOADING]", basename(region_file), "..."))

    # Read BED file
    bed_data <- fread(region_file, header = FALSE)
    colnames(bed_data)[1:3] <- c("chr", "start", "end")

    # Apply species-agnostic filtering if requested
    if (!is.null(opt$filter_standard_chroms) && opt$filter_standard_chroms) {
      # Universal filter: 1-99 + single-letter sex chromosomes
      valid_chroms <- c(
        as.character(1:99), # Covers all vertebrates
        "Z", "W", "z", "w", # Birds
        "X", "Y", "x", "y" # Mammals
      )
      bed_filtered <- bed_data[bed_data$chr %in% valid_chroms, ]
      log_msg(paste("[PACKAGE] Chromosome filtering:", nrow(bed_data), "->", nrow(bed_filtered), "regions"))
    } else {
      # No filtering - use all chromosomes from BED file
      bed_filtered <- bed_data
      log_msg(paste("[PACKAGE] Using all chromosomes:", nrow(bed_data), "regions"))
    }

    # Exclude sex chromosomes if requested (safety feature)
    if (!is.null(opt$exclude_sex_chroms) && opt$exclude_sex_chroms) {
      sex_chroms <- c("X", "Y", "Z", "W", "x", "y", "z", "w")
      bed_before <- bed_filtered
      bed_filtered <- bed_filtered[!(bed_filtered$chr %in% sex_chroms), ]
      log_msg(paste("[SAFETY] Sex chromosome exclusion:", nrow(bed_before), "->", nrow(bed_filtered), "regions"))
    }

    if (nrow(bed_filtered) == 0) {
      log_msg(paste("[WARN] No valid chromosomes found in", basename(region_file), "- skipping"))
      next
    }

    # Create GRanges and store
    gr <- toGRanges(bed_filtered[, 1:3], format = "BED")
    region_data[[reg_label]] <- list(
      gr = gr,
      label = reg_label,
      file_path = region_file
    )

    log_msg(paste("[OK] Loaded:", reg_label, "-", length(gr), "regions"))
  }

  if (length(region_data) == 0) {
    stop("[ERROR] No valid BED files loaded after species and chromosome filtering")
  }

  log_msg(paste("[SUCCESS] Loaded", length(region_data), "region sets for analysis"))
} else {
  stop(
    "[ERROR] Invalid strategy. This public version only supports strategy='traditional'.\n",
    "       Please provide genomic regions via --window_file or --region_dir."
  )
}

# ---------------------------
# Enhanced Model Selection System
# ---------------------------

# REQUIRED PACKAGES for all models:
# Standard: randomForest, glmnet, class (knn), MASS (lda), e1071 (svm), nnet
# Additional: xgboost, naivebayes, extraTrees
# Install missing packages with: install.packages(c("xgboost", "naivebayes", "extraTrees"))

# Updated model specs (9 models instead of 10):
model_specs <- list(
  rf = "rf",
  glmnet = "glmnet",
  knn = "knn",
  lda = "lda",
  svm = "svmRadial",
  xgboost = "xgbTree",
  nb = "nb",
  ranger = "ranger", # REPLACEMENT: Fast Random Forest
  logreg = "glm",
  nnet = "nnet"
)

# Model metadata for categorization and descriptions
model_info <- list(
  rf = list(
    method = "rf",
    name = "Random Forest",
    description = "Ensemble method using decision trees",
    speed = "medium",
    interpretability = "medium",
    categories = c("ensemble", "nonlinear", "robust")
  ),
  glmnet = list(
    method = "glmnet",
    name = "Elastic Net",
    description = "Regularized linear regression (Lasso/Ridge)",
    speed = "fast",
    interpretability = "high",
    categories = c("linear", "regularized", "fast", "interpretable")
  ),
  knn = list(
    method = "knn",
    name = "k-Nearest Neighbors",
    description = "Instance-based learning algorithm",
    speed = "fast",
    interpretability = "low",
    categories = c("fast", "nonlinear", "instance_based", "simple")
  ),
  lda = list(
    method = "lda",
    name = "Linear Discriminant Analysis",
    description = "Linear classification method",
    speed = "fast",
    interpretability = "high",
    categories = c("fast", "linear", "classical", "interpretable")
  ),
  svm = list(
    method = "svmRadial",
    name = "Support Vector Machine",
    description = "SVM with radial basis function kernel",
    speed = "slow",
    interpretability = "low",
    categories = c("nonlinear", "kernel_based", "robust")
  ),
  xgboost = list(
    method = "xgbTree",
    name = "XGBoost",
    description = "Gradient boosting trees - often best performer",
    speed = "slow",
    interpretability = "medium",
    categories = c("ensemble", "boosting", "nonlinear", "high_performance")
  ),
  nb = list(
    method = "nb",
    name = "Naive Bayes",
    description = "Probabilistic classifier, good for genomic data",
    speed = "fast",
    interpretability = "high",
    categories = c("fast", "probabilistic", "simple", "genomic_friendly")
  ),
  extratrees = list(
    method = "extraTrees",
    name = "Extra Trees",
    description = "Extremely randomized trees, often better than RF",
    speed = "medium",
    interpretability = "medium",
    categories = c("ensemble", "nonlinear", "fast_ensemble")
  ),
  logreg = list(
    method = "glm",
    name = "Logistic Regression",
    description = "Linear logistic regression - interpretable baseline",
    speed = "fast",
    interpretability = "high",
    categories = c("linear", "fast", "classical", "interpretable", "baseline")
  ),
  nnet = list(
    method = "nnet",
    name = "Neural Network",
    description = "Feedforward neural network for complex patterns",
    speed = "medium",
    interpretability = "low",
    categories = c("nonlinear", "neural", "complex_patterns")
  )
)

# Updated model categories:
model_categories <- list(
  # Speed-based categories
  fast = c("knn", "lda", "glmnet", "nb", "logreg"),
  medium = c("rf", "ranger", "nnet"), # CHANGED: ranger instead of extratrees
  slow = c("svm", "xgboost"),

  # Algorithm type categories
  linear = c("glmnet", "lda", "logreg"),
  nonlinear = c("rf", "svm", "knn", "xgboost", "ranger", "nnet"), # CHANGED

  # Ensemble methods
  ensemble = c("rf", "xgboost", "ranger"), # CHANGED: ranger instead of extratrees
  boosting = c("xgboost"),

  # Performance categories
  high_performance = c("xgboost", "rf", "ranger"), # CHANGED: ranger instead of extratrees
  baseline = c("logreg", "nb"),
  robust = c("rf", "svm", "ranger"), # CHANGED

  # Domain-specific categories
  biomarker = c("xgboost", "rf", "glmnet", "ranger", "svm"), # CHANGED: ranger instead of extratrees

  # Special groupings
  tree_based = c("rf", "xgboost", "ranger"), # CHANGED: ranger instead of extratrees
  fast_ensemble = c("ranger"), # CHANGED: ranger instead of extratrees
  all_except_xgb = c("rf", "glmnet", "knn", "lda", "svm", "nb", "ranger", "logreg", "nnet")
)

# Enhanced model selection parser
parse_model_selection <- function(models_string, available_models = model_specs, categories = model_categories) {
  # Clean up the input
  models_string <- tolower(trimws(models_string))

  # Parse comma-separated list
  requested_items <- strsplit(models_string, ",")[[1]]
  requested_items <- trimws(requested_items) # Remove whitespace

  selected_models <- character(0)
  processed_items <- character(0)
  invalid_items <- character(0)

  for (item in requested_items) {
    if (item == "all") {
      # Add all available models
      selected_models <- c(selected_models, names(available_models))
      processed_items <- c(processed_items, paste0(item, " (", length(names(available_models)), " models)"))
    } else if (item %in% names(categories)) {
      # Add models from category
      category_models <- categories[[item]]
      selected_models <- c(selected_models, category_models)
      processed_items <- c(processed_items, paste0(item, " (", paste(category_models, collapse = ","), ")"))
    } else if (item %in% names(available_models)) {
      # Add individual model
      selected_models <- c(selected_models, item)
      processed_items <- c(processed_items, item)
    } else {
      # Invalid item
      invalid_items <- c(invalid_items, item)
    }
  }

  # Remove duplicates while preserving order
  selected_models <- unique(selected_models)

  # Report results
  if (length(invalid_items) > 0) {
    log_msg(paste("[WARN] Invalid items ignored:", paste(invalid_items, collapse = ", ")))
    log_msg("[INFO] Use --list-models to see available options")
  }

  if (length(selected_models) == 0) {
    stop(paste(
      "[ERROR] No valid models specified.",
      "\n[INFO] Use --list-models to see available models and categories",
      "\n[INFO] Examples: --models all, --models fast, --models rf,svm"
    ))
  }

  log_msg(paste("[TARGET] Processed:", paste(processed_items, collapse = ", ")))
  log_msg(paste("[MODEL] Final model selection:", paste(selected_models, collapse = ", ")))

  return(selected_models)
}

# Parse models and setup
model_list <- parse_model_selection(opt$models)

# SECURITY FIX: Safe parsing of marker_range (no eval/parse)
marker_seq <- tryCatch(
  {
    parts <- strsplit(opt$marker_range, ":")[[1]]
    if (length(parts) != 2) stop("Invalid marker_range format")
    start <- as.integer(parts[1])
    end <- as.integer(parts[2])
    if (is.na(start) || is.na(end) || start > end) stop("Invalid marker range values")
    seq(start, end)
  },
  error = function(e) {
    stop(paste0(
      "[ERROR] Invalid --marker_range: ", e$message,
      ". Use format: 'start:end' (e.g., '2:50')"
    ))
  }
)

# ADD THESE MISSING VARIABLES:
current_model_list <- model_list

# ---------------------------
# ADD RESUME CHECK HERE - RIGHT BEFORE MAIN PROCESSING LOOP
# ---------------------------

region_data_original <- region_data

# ---------------------------
# ADD RESUME CHECK HERE - RIGHT BEFORE MAIN PROCESSING LOOP
# ---------------------------

region_data_original <- region_data

# Skip completed analysis check if using cached precomputed data
if (length(region_data) > 0 &&
  !is.null(region_data[[1]]$has_precomputed) &&
  region_data[[1]]$has_precomputed) {
  log_msg("[START] Using precomputed data - skipping completion check")
  analysis_status <- list(
    completed = character(0),
    remaining = region_data,
    total = length(region_data),
    completed_count = 0,
    remaining_count = length(region_data)
  )
} else {
  # Check for completed analyses and filter region list
  analysis_status <- check_completed_analyses(region_data, opt$output_dir)
}

if (opt$force_ml) {
  log_msg("[PROCESS] Force ML mode enabled - will use precomputed files but redo ML analysis")

  # Keep completed regions that have precomputed files but force ML rerun
  force_ml_regions <- list()
  for (region_name in names(region_data)) {
    region_info <- region_data[[region_name]]
    reg_label <- region_info$label

    # Check if precomputed files exist
    precomputed_matrix_file <- file.path(opt$output_dir, paste0("precomputed_methylation_matrix_", reg_label, ".rds"))
    dmr_summary_file <- file.path(opt$output_dir, paste0("dmr_summary_", reg_label, ".csv"))

    if (file.exists(precomputed_matrix_file) && file.exists(dmr_summary_file)) {
      log_msg(paste("[DIR] Will reuse precomputed files for:", reg_label))
      force_ml_regions[[region_name]] <- region_info
      force_ml_regions[[region_name]]$has_precomputed <- TRUE
      force_ml_regions[[region_name]]$precomputed_matrix_file <- precomputed_matrix_file
      force_ml_regions[[region_name]]$dmr_summary_file <- dmr_summary_file
    } else {
      log_msg(paste("[WARN] No precomputed files found for:", reg_label, "- will run full analysis"))
      force_ml_regions[[region_name]] <- region_info
      force_ml_regions[[region_name]]$has_precomputed <- FALSE
    }
  }

  region_data <- force_ml_regions
  existing_results <- data.frame() # Don't load existing ML results
} else {
  # Normal mode - load existing results and skip completed regions
  existing_results <- load_completed_results(analysis_status$completed, opt$output_dir)
  if (nrow(existing_results) > 0) {
    results_summary <- rbind(results_summary, existing_results)
  }
  region_data <- analysis_status$remaining
}

if (length(region_data) == 0 && !opt$force_ml) {
  log_msg("[SUCCESS] All analyses already completed!")

  # Still generate final summary
  if (nrow(results_summary) > 0) {
    summary_file <- file.path(opt$output_dir, "analysis_summary.csv")
    write.csv(results_summary, summary_file, row.names = FALSE)
    log_msg(paste("[DATA] Complete results summary saved to:", basename(summary_file)))

    best_results <- results_summary[order(-results_summary$accuracy), ]
    log_msg("[DATA] Best results from all analyses:")

    for (i in 1:min(5, nrow(best_results))) {
      result <- best_results[i, ]
      log_msg(sprintf(
        "   %s: %s | %d markers | Acc=%.3f | AUC=%.3f | Split=%d/%d",
        result$region, result$model, result$markers, result$accuracy, result$auc,
        result$train_samples, result$test_samples
      ))
    }
  } else if (length(region_data) == 0 && opt$force_ml) {
    # Restore region_data from original list since we want to force ML
    region_data <- region_data_original # We need to save this earlier
    log_msg("[PROCESS] Force ML mode: Restoring all regions for ML analysis")
  }


  total_end_time <- Sys.time()
  total_time <- round(as.numeric(total_end_time - total_start_time, units = "mins"), 1)
  log_msg(paste("[FAST] Analysis completed using cached results in", total_time, "minutes"))
  log_msg(paste("[DIR] All outputs available in:", opt$output_dir))

  cat("[COMPLETE] Analysis complete (resumed from cache)!\n")
  quit(save = "no", status = 0)
}

log_msg(paste("[START] Proceeding with", length(region_data), "remaining regions..."))


# ---------------------------
# Main processing loop
# ---------------------------

for (region_idx in seq_along(region_data)) {
  region_info <- region_data[[region_idx]]
  reg_label <- region_info$label

  # ADD THESE LINES (FIXED):
  cat("\n")
  log_msg("================================================================================")
  log_msg(paste("[REGION START]", reg_label, "-", region_idx, "of", length(region_data)))
  log_msg(paste("[REGION] Processing:", reg_label))
  flush.console()

  log_msg(paste("[FOLDER] Processing region", region_idx, "of", length(region_data), ":", reg_label))

  if (opt$dryrun) {
    log_msg("[TEST] Dry-run mode - skipping actual processing")
    next
  }
  region_start_time <- Sys.time()

  tryCatch(
    {
      # Check if we have precomputed data (cached or force_ml)
      if (!is.null(region_info$has_precomputed) && region_info$has_precomputed) {
        log_msg(paste("[DIR] Loading precomputed files for", reg_label))

        # Load precomputed methylation matrix
        precomputed_matrix_file <- region_info$precomputed_matrix_file
        dmr_summary_file <- region_info$dmr_summary_file

        log_msg(paste("[DATA] Loading precomputed matrix:", basename(precomputed_matrix_file)))
        X_all_dmrs <- readRDS(precomputed_matrix_file)

        log_msg(paste("[DATA] Loading DMR summary:", basename(dmr_summary_file)))
        dmrs <- read.csv(dmr_summary_file, stringsAsFactors = FALSE)

        # Rename columns to match expected format
        if ("scaffold" %in% colnames(dmrs)) {
          colnames(dmrs)[colnames(dmrs) == "scaffold"] <- "chr"
        }
        if ("meth_difference" %in% colnames(dmrs)) {
          colnames(dmrs)[colnames(dmrs) == "meth_difference"] <- "meth.diff"
        }
        if ("width_bp" %in% colnames(dmrs)) {
          colnames(dmrs)[colnames(dmrs) == "width_bp"] <- "width"
        }

        log_msg(paste("[OK] Loaded and standardized", nrow(dmrs), "DMR records"))

        # Validate and transpose if needed
        expected_samples <- length(meth)
        expected_dmrs <- nrow(dmrs)

        if (ncol(X_all_dmrs) == expected_samples && nrow(X_all_dmrs) == expected_dmrs) {
          X_all_dmrs <- t(X_all_dmrs)
          log_msg(paste("[OK] Transposed to:", nrow(X_all_dmrs), "samples x", ncol(X_all_dmrs), "DMRs"))
        }

        log_msg(paste("[OK] Loaded precomputed data:", nrow(X_all_dmrs), "samples x", ncol(X_all_dmrs), "DMRs"))
        # NEW CODE BLOCK - ADD THIS
        # Reconstruct ml_data from precomputed matrix
        log_msg("[PROCESS] Reconstructing ml_data from precomputed matrix...")

        # Get treatment labels
        if (exists("treatment_processed") && !is.null(treatment_processed)) {
          y <- factor(treatment_processed$labels[1:nrow(X_all_dmrs)])
        } else {
          y <- factor(treatment_info[1:nrow(X_all_dmrs)])
        }

        # Create ml_data from precomputed matrix
        ml_data <- bind_cols(X_all_dmrs, label = y)

        # Verify
        log_msg(paste("[OK] Reconstructed ml_data:", nrow(ml_data), "samples x", ncol(ml_data), "columns"))
        log_msg(paste("[DATA] Class distribution:", paste(table(ml_data$label), collapse = ", ")))
        # END NEW CODE BLOCK

        # Ensure DMR_ID column exists
        if (!"DMR_ID" %in% colnames(dmrs)) {
          dmrs$DMR_ID <- paste0("DMR_", 1:nrow(dmrs))
          # Traditional strategy - no merging
          if (!"n_merged" %in% colnames(dmrs)) dmrs$n_merged <- 1
          if (!"merge_span" %in% colnames(dmrs)) dmrs$merge_span <- dmrs$end - dmrs$start + 1
        }

        # Set marker range
        if (nrow(dmrs) < max(marker_seq)) {
          marker_seq_region <- marker_seq[marker_seq <= nrow(dmrs)]
        } else {
          marker_seq_region <- marker_seq
        }

        # Skip directly to ML training
        goto_ml_training <- TRUE
      } else {
        goto_ml_training <- FALSE

        # ONLY NOW load gr and do genomic processing
        gr <- region_info$gr
        log_msg(paste("[PACKAGE] Loaded", length(gr), "regions"))
      }

      # Windows are already loaded in gr
      log_msg(paste("[PACKAGE] Loaded", length(gr), "regions"))

      # Extract window size and filter by size (for generated windows, this might be mixed sizes)
      if (!goto_ml_training) {
        # Extract window size and filter by size (for generated windows, this might be mixed sizes)
        if (grepl("_(\\d+)bp", reg_label)) {
          window_size <- as.numeric(gsub(".*_(\\d+)bp.*", "\\1", reg_label))
          if (!is.na(window_size)) {
            verbose_msg(paste("Filtering regions by window size:", window_size, "bp"))
            # size_filtered <- gr[width(gr) >= window_size]
            size_filtered <- gr ## no filtering at all!!
            log_msg(paste("[PACKAGE] Size filtered:", length(gr), "->", length(size_filtered), "regions"))
            gr <- size_filtered
          }
        }

        if (length(gr) == 0) {
          log_msg("[WARN] No regions left after size filtering")
          next
        }

        # Check scaffold overlap
        verbose_msg("Checking scaffold overlap...")
        meth_scaffolds <- unique(unlist(lapply(meth, function(x) {
          if (is(x, "methylRaw")) {
            return(unique(as.character(x$chr)))
          } else {
            return(character(0))
          }
        })))
        bed_scaffolds <- unique(as.character(seqnames(gr)))
        common_scaffolds <- intersect(bed_scaffolds, meth_scaffolds)

        log_msg(paste("[PACKAGE] Scaffold overlap:", length(common_scaffolds), "common of", length(bed_scaffolds), "total"))
        if (opt$debug) {
          debug_msg(paste("Methylation scaffolds count:", length(meth_scaffolds)))
          debug_msg(paste("BED scaffolds count:", length(bed_scaffolds)))
          debug_msg(paste("First 5 methylation scaffolds:", paste(head(meth_scaffolds, 5), collapse = ", ")))
          debug_msg(paste("First 5 BED scaffolds:", paste(head(bed_scaffolds, 5), collapse = ", ")))
          if (length(common_scaffolds) > 0) {
            debug_msg(paste("First 5 common scaffolds:", paste(head(common_scaffolds, 5), collapse = ", ")))
          }
        }
        verbose_msg(paste("Methylation scaffolds:", length(meth_scaffolds)))
        verbose_msg(paste("BED scaffolds:", length(bed_scaffolds)))

        if (length(common_scaffolds) == 0) {
          log_msg("[ERROR] No scaffold overlap - skipping region")
          next
        }

        # Filter to common scaffolds
        gr <- gr[seqnames(gr) %in% common_scaffolds]
        log_msg(paste("[PACKAGE] After scaffold filtering:", length(gr), "regions"))


        # ---------------------------
        # CACHED REGIONCOUNTS WITH SMART CACHING
        # ---------------------------

        # Create regionCounts cache directory in organized structure
        regioncounts_cache_dir <- file.path(output_dirs$precomputed_data, "regioncounts_cache")
        dir.create(regioncounts_cache_dir, recursive = TRUE, showWarnings = FALSE)

        # Clean old cache files (optional - keeps cache directory manageable)
        if (region_idx == 1) { # Only clean on first region to avoid repeated cleaning
          clean_old_regioncounts_cache(regioncounts_cache_dir, max_age_days = 30)
        }

        # Use cached regionCounts - this will save/load automatically
        verbose_msg("Calculating region counts with coverage filtering (cached)...")
        rc <- regionCounts_with_cache(
          meth = meth,
          gr = gr,
          min_coverage = opt$min_coverage,
          cache_dir = regioncounts_cache_dir
        )

        # Log the results (same as before)
        total_regions_before <- sum(sapply(rc, nrow))
        log_msg(paste("[PACKAGE] Region counts completed:", length(rc), "samples,", total_regions_before, "total regions"))

        # IMPORTANT: Add this debug line to understand cache performance
        cache_files <- list.files(regioncounts_cache_dir, pattern = "regionCounts_.*\\.rds$")
        if (length(cache_files) > 0) {
          log_msg(paste("[DIR] regionCounts cache directory contains", length(cache_files), "cached files"))

          # Show cache file info for debugging
          if (opt$debug) {
            for (i in seq_along(cache_files)) {
              cache_file_path <- file.path(regioncounts_cache_dir, cache_files[i])
              cache_size_mb <- round(file.size(cache_file_path) / 1024^2, 1)
              cache_age_hours <- round(as.numeric(Sys.time() - file.info(cache_file_path)$mtime, units = "hours"), 1)
              debug_msg(paste("Cache", i, ":", cache_files[i], "(", cache_size_mb, "MB,", cache_age_hours, "hours old)"))
            }
          }
        }

        # No additional coverage filtering needed - already done in regionCounts
        rc_filt <- rc

        total_regions_before <- sum(sapply(rc, nrow))
        total_regions_after <- sum(sapply(rc_filt, nrow))

        log_msg(paste("[PACKAGE] After coverage filtering:", total_regions_after, "total regions across", length(rc_filt), "samples"))

        if (total_regions_after == 0) {
          log_msg("[WARN] No regions left after coverage filtering")
          next
        }

        # =========================================================================
        # SIMPLE UNITE THRESHOLD - DIRECT PERCENTAGE
        # =========================================================================
        #    log_msg(paste("[UNITE] Requiring data in", opt$min_samples_per_region * 100,
        #                  "% of samples per group"))
        #
        #    verbose_msg(paste("Uniting regions (destrand=FALSE, min.per.group =",
        #                     opt$min_samples_per_region, "as fraction)..."))

        # ===== DIAGNOSTIC =====
        log_msg("[DEBUG] Checking regionCounts before unite...")
        log_msg(paste("[DEBUG] Number of samples in rc_filt:", length(rc_filt)))

        # Check how many regions each sample has
        sample_region_counts <- sapply(rc_filt, nrow)
        log_msg(paste(
          "[DEBUG] Regions per sample - min:", min(sample_region_counts),
          "max:", max(sample_region_counts),
          "median:", median(sample_region_counts)
        ))

        # Check treatment distribution
        if (!is.null(rc_filt@treatment)) {
          treatment_dist <- table(rc_filt@treatment)
          log_msg(paste("[DEBUG] Treatment distribution:", paste(names(treatment_dist), "=", treatment_dist, collapse = ", ")))
        }

        # ===== CALCULATE min.per.group FROM FLAG =====
        # Convert fraction to integer samples per group (with R's L suffix requirement)
        treatment_table <- table(rc_filt@treatment)
        min_samples_per_group_raw <- floor(min(treatment_table) * opt$min_samples_per_region)
        min_samples_per_group <- as.integer(min_samples_per_group_raw)

        # Ensure it's a proper R integer
        if (!is.integer(min_samples_per_group)) {
          min_samples_per_group <- as.integer(min_samples_per_group)
        }

        # Safety check: ensure at least 1 sample per group
        if (min_samples_per_group < 1) {
          min_samples_per_group <- 1
          log_msg(paste("[WARN] Calculated min.per.group < 1, setting to 1"))
        }

        # CRITICAL: Determine if we should use NULL or the calculated value
        # When flag is 1.0, use NULL (methylKit's default behavior)
        # When flag is anything else, use the calculated value
        if (opt$min_samples_per_region == 1.0) {
          unite_min_per_group <- NULL
          log_msg(paste("[UNITE] Using min.per.group = NULL (flag = 1.0, methylKit default)"))
          log_msg(paste("  Group sizes:", paste(names(treatment_table), "=", treatment_table, collapse = ", ")))
          log_msg(paste("  Smallest group:", min(treatment_table), "samples"))
          log_msg(paste("  This will require regions present in ALL samples per group"))
        } else {
          unite_min_per_group <- min_samples_per_group
          log_msg(paste(
            "[UNITE] Using min.per.group =", min_samples_per_group,
            "samples (", round(opt$min_samples_per_region * 100, 1), "% of smallest group)"
          ))
          log_msg(paste("  Group sizes:", paste(names(treatment_table), "=", treatment_table, collapse = ", ")))
          log_msg(paste("  Smallest group:", min(treatment_table), "samples"))
          log_msg(paste(
            "  Requiring regions present in >=", min_samples_per_group,
            "samples per group"
          ))
        }
        # ===== END DIAGNOSTIC =====

        unite_start <- Sys.time()

        # Use the determined value (either NULL or calculated)
        united <- methylKit::unite(rc_filt, destrand = FALSE, min.per.group = unite_min_per_group)

        log_msg(paste("[UNITE] Complete -", nrow(united), "united regions"))
        flush.console()

        unite_end <- Sys.time()
        log_msg(paste(
          "[PACKAGE] United regions:", nrow(united), "in",
          round(as.numeric(unite_end - unite_start, units = "secs")), "seconds"
        ))

        if (nrow(united) < 1) {
          log_msg("[WARN] Too few united regions for meaningful analysis")
          log_msg(paste("   Try lowering --min_samples_per_region (currently", opt$min_samples_per_region, ")"))
          if (!is.null(unite_min_per_group)) {
            log_msg(paste("   Current requirement:", unite_min_per_group, "samples per group"))
          }
          next
        }

        # =========================================================================

        # Export the united methylation object to QS file (preserves raw count format)
        # Add this code after the united object is created (around line ~1400)

        # Create export filename based on region label and timestamp
        export_timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
        united_export_file <- file.path(opt$output_dir, paste0("united_meth_object_", reg_label, "_", export_timestamp, ".qs"))

        # Save the united object (this preserves the methylBase S4 object with raw counts)
        log_msg(paste("[SAVE] Exporting united methylation object to:", basename(united_export_file)))
        qsave(united, united_export_file, preset = "fast")

        # Verify the export and show column structure
        export_size_mb <- round(file.size(united_export_file) / 1024^2, 1)
        log_msg(paste("[OK] United object exported:", export_size_mb, "MB"))

        # Get the actual data to show structure
        united_data <- getData(united)
        log_msg(paste("[DATA] Object contains:", nrow(united_data), "regions x", length(united@treatment), "samples"))
        log_msg(paste("[DATA] Column structure: chr, start, end, strand, then coverage/numCs/numTs triplets for each sample"))
        log_msg(paste("[DATA] First few columns:", paste(colnames(united_data)[1:min(10, ncol(united_data))], collapse = ", ")))

        # Verify the expected column pattern (coverage, numCs, numTs repeating)
        sample_cols <- colnames(united_data)[5:ncol(united_data)] # Skip chr, start, end, strand
        coverage_cols <- sample_cols[seq(1, length(sample_cols), 3)]
        numCs_cols <- sample_cols[seq(2, length(sample_cols), 3)]
        numTs_cols <- sample_cols[seq(3, length(sample_cols), 3)]

        log_msg(paste("[DATA] Coverage columns:", length(coverage_cols), "samples"))
        log_msg(paste("[DATA] First coverage columns:", paste(head(coverage_cols, 3), collapse = ", ")))

        # Save metadata including treatment information
        united_metadata <- list(
          export_timestamp = export_timestamp,
          region_file = reg_label,
          n_regions = nrow(united_data),
          n_samples = length(united@treatment),
          treatment_groups = as.vector(table(united@treatment)),
          treatment_vector = united@treatment,
          column_structure = list(
            genomic_columns = c("chr", "start", "end", "strand"),
            sample_triplet_pattern = c("coverage", "numCs", "numTs"),
            total_columns = ncol(united_data)
          ),
          analysis_parameters = list(
            qval_threshold = opt$qval,
            min_coverage = opt$min_coverage,
            min_cpg = opt$min_cpg,
            strategy = opt$strategy
          )
        )

        # Save metadata as R object (more reliable than JSON for complex structures)
        metadata_file <- file.path(opt$output_dir, paste0("united_meth_metadata_", reg_label, "_", export_timestamp, ".rds"))
        saveRDS(united_metadata, metadata_file)
        log_msg(paste("[SUMMARY] Metadata saved to:", basename(metadata_file)))

        # ===== CRITICAL CORRECTIONS FOR DMR REPRODUCIBILITY =====
        # Replace the entire DMR calculation section in your new script with this:

        # Calculate differential methylation (EXACT REPRODUCTION OF OLD SCRIPT)
        verbose_msg("Calculating differential methylation (this will take time)...")
        diff_start <- Sys.time()

        # CRITICAL: Use the SAME parameters as old script
        diff <- calculateDiffMeth(united,
          overdispersion = "none", # OLD SCRIPT PARAMETER
          adjust = "fdr"
        ) # OLD SCRIPT PARAMETER

        diff_end <- Sys.time()
        log_msg(paste("[DATA] Differential methylation:", nrow(diff), "results in", round(as.numeric(diff_end - diff_start, units = "mins"), 1), "minutes"))

        # DEBUG: Check original DMR sizes from calculateDiffMeth
        if (nrow(diff) > 0) {
          original_sizes <- diff$end - diff$start + 1
          log_msg(paste(
            "[DATA] ORIGINAL DMR sizes from calculateDiffMeth:",
            "min=", min(original_sizes),
            "median=", median(original_sizes),
            "max=", max(original_sizes)
          ))
        }

        if (nrow(diff) > 0) {
          meth_diff_range <- range(diff$meth.diff, na.rm = TRUE)
          qval_range <- range(diff$qvalue, na.rm = TRUE)
          log_msg(paste("[DATA] Methylation difference range:", round(meth_diff_range[1], 2), "to", round(meth_diff_range[2], 2), "%"))
          log_msg(paste("[DATA] Q-value range:", signif(qval_range[1], 3), "to", signif(qval_range[2], 3)))
        }

        # Get DMRs
        # CRITICAL: Sort by qvalue first (EXACT reproduction of old script)
        diff_sorted <- diff[order(diff$qvalue), ]

        # CRITICAL: Filter by qvalue threshold FIRST
        dmrs_qval_filtered <- diff_sorted[diff_sorted$qvalue < opt$qval, ]
        log_msg(paste("[DATA] DMRs passing q-value filter (<=", opt$qval, "):", nrow(dmrs_qval_filtered)))

        # CRITICAL: Apply methylation difference filter AFTER q-value filter
        # Apply methylation difference filter AFTER q-value filter (if enabled)
        if (!opt$disable_meth_filter) {
          verbose_msg(paste("Applying methylation difference filter: >=", opt$hyper_threshold, "% hyper OR <=", opt$hypo_threshold, "% hypo"))

          dmrs_meth_filtered <- dmrs_qval_filtered[
            dmrs_qval_filtered$meth.diff >= opt$hyper_threshold | dmrs_qval_filtered$meth.diff <= opt$hypo_threshold,
          ]

          log_msg(paste("[DATA] DMRs after meth.diff filter (>=", opt$hyper_threshold, "% hyper OR <=", opt$hypo_threshold, "% hypo):", nrow(dmrs_meth_filtered)))
          log_msg(paste("[DATA] Hypermethylated DMRs (>=", opt$hyper_threshold, "%):", sum(dmrs_meth_filtered$meth.diff >= opt$hyper_threshold)))
          log_msg(paste("[DATA] Hypomethylated DMRs (<=", opt$hypo_threshold, "%):", sum(dmrs_meth_filtered$meth.diff <= opt$hypo_threshold)))

          # Handle case where filtering is too strict
          if (nrow(dmrs_meth_filtered) < 10) {
            log_msg("[WARN] Meth difference filter too strict, taking first 10 DMRs by q-value to avoid problems")
            dmrs <- dmrs_qval_filtered[1:min(10, nrow(dmrs_qval_filtered)), ]
          } else {
            log_msg("[OK] Methylation difference filtering successful")
            dmrs <- dmrs_meth_filtered
          }
        } else {
          log_msg("[DATA] Methylation difference filtering DISABLED - using only q-value threshold")
          dmrs <- dmrs_qval_filtered
        }

        log_msg(paste("[DATA] DMRs after meth.diff filter (>=5% hyper OR >=2% hypo):", nrow(dmrs_meth_filtered)))
        log_msg(paste("[DATA] Hypermethylated DMRs (>=5%):", sum(dmrs_meth_filtered$meth.diff >= 5)))
        log_msg(paste("[DATA] Hypomethylated DMRs (<=-2%):", sum(dmrs_meth_filtered$meth.diff <= -2)))

        # Handle case where filtering is too strict (matching old script logic)
        if (nrow(dmrs_meth_filtered) < 10) {
          log_msg("[WARN] Meth difference filter too strict, taking first 10 DMRs by q-value to avoid problems")
          dmrs <- dmrs_qval_filtered[1:min(10, nrow(dmrs_qval_filtered)), ]
        } else {
          log_msg("[OK] Methylation difference filtering successful")
          dmrs <- dmrs_meth_filtered
        }

        # Traditional strategy - no merging needed
        dmrs$n_merged <- 1
        dmrs$merge_span <- dmrs$end - dmrs$start + 1

        # Add DMR IDs
        dmrs <- dmrs[order(dmrs$qvalue, -abs(dmrs$meth.diff)), ]
        dmrs$DMR_ID <- paste0("DMR_", 1:nrow(dmrs))

        # Check DMR sizes
        log_msg(paste(
          "[DATA] DMR sizes:",
          "min=", min(dmrs$end - dmrs$start + 1),
          "median=", median(dmrs$end - dmrs$start + 1),
          "max=", max(dmrs$end - dmrs$start + 1)
        ))

        log_msg(paste("[DATA] Final DMRs for analysis:", nrow(dmrs)))

        if (nrow(dmrs) < max(marker_seq)) {
          log_msg(paste("[WARN] Only", nrow(dmrs), "DMRs found, adjusting marker range"))
          marker_seq_adj <- marker_seq[marker_seq <= nrow(dmrs)]
          if (length(marker_seq_adj) == 0) {
            log_msg("[WARN] No valid marker sets - skipping region")
            next
          }
          marker_seq_region <- marker_seq_adj
          log_msg(paste("[DATA] Adjusted marker range:", min(marker_seq_region), "to", max(marker_seq_region)))
        } else {
          marker_seq_region <- marker_seq
        }

        # Pre-compute methylation matrix ONLY for significant DMRs that pass q-value threshold
        log_analysis_step(
          "Pre-computing methylation matrix for significant DMRs",
          paste("q <=", opt$qval)
        )
        precompute_start <- Sys.time()

        # Determine maximum number of markers we might need
        max_markers_needed <- max(marker_seq_region)
        verbose_msg(paste("Maximum markers needed:", max_markers_needed, "of", nrow(dmrs), "available DMRs"))

        # Only process the top DMRs up to the maximum we'll actually use
        dmrs_to_process <- min(max_markers_needed, nrow(dmrs))
        verbose_msg(paste("Processing top", dmrs_to_process, "DMRs instead of all", nrow(dmrs)))

        # Create GRanges for only the DMRs we'll actually use
        if (is.data.frame(dmrs)) {
          # DMRs is already a data frame (from enhanced processing or regular processing)
          subset_dmrs <- dmrs[1:dmrs_to_process, ]
          all_dmrs_gr <- GRanges(
            seqnames = subset_dmrs$chr,
            ranges = IRanges(start = subset_dmrs$start, end = subset_dmrs$end)
          )
        } else {
          # DMRs is still a methylDiff object (shouldn't happen but safety check)
          subset_dmrs <- dmrs[1:dmrs_to_process, ]
          all_dmrs_gr <- as(subset_dmrs, "GRanges")
        }

        # Update dmrs to only contain the subset we're processing
        dmrs <- subset_dmrs

        # OPTIMIZED: Extract methylation ONLY for significant DMRs (not all 100k+ regions)
        verbose_msg("Extracting methylation for significant DMRs from united object...")

        # Get positions from united object
        united_data <- getData(united)
        united_positions <- paste(united_data$chr, united_data$start, united_data$end, sep = ".")

        # Get positions from significant DMRs (subset_dmrs)
        sig_positions <- paste(subset_dmrs$chr, subset_dmrs$start, subset_dmrs$end, sep = ".")

        # Find which rows in united match our significant DMRs
        match_indices <- match(sig_positions, united_positions)
        valid_matches <- !is.na(match_indices)

        if (sum(valid_matches) < length(sig_positions)) {
          log_msg(paste("[WARN]", sum(!valid_matches), "DMRs not found in united object"))
          subset_dmrs <- subset_dmrs[valid_matches, ]
          sig_positions <- sig_positions[valid_matches]
          match_indices <- match_indices[valid_matches]
        }

        log_msg(paste(
          "[DATA] Subsetting united object to significant DMRs only:",
          length(match_indices), "of", format(nrow(united_data), big.mark = ","),
          "total regions (",
          round(100 * length(match_indices) / nrow(united_data), 2), "%)"
        ))

        # Create a subset of united containing ONLY the significant DMRs
        # This avoids processing all 126k regions
        united_subset <- united[match_indices, ]

        # Extract methylation from the MUCH smaller subset
        perc_meth_subset <- percMethylation(united_subset, rowids = TRUE)

        # Transpose to get samples as rows, DMRs as columns
        perc_meth_all <- t(perc_meth_subset)

        log_msg(paste(
          "[OK] Fast extraction complete:", nrow(perc_meth_all), "samples x",
          ncol(perc_meth_all), "significant DMRs"
        ))

        # Verify sample count matches
        if (nrow(perc_meth_all) != length(meth)) {
          stop(paste("ERROR: Sample mismatch. Matrix:", nrow(perc_meth_all), "vs methylKit:", length(meth)))
        }

        # Extract methylation data for matched DMRs
        X_all_dmrs <- perc_meth_all[, 1:ncol(perc_meth_all), drop = FALSE]

        if (is.null(X_all_dmrs) || nrow(X_all_dmrs) == 0) {
          log_warning("Failed to extract methylation data", reg_label)
          next
        }

        # Update dmrs to match what we extracted
        dmrs <- subset_dmrs

        # Handle missing values
        for (col in 1:ncol(X_all_dmrs)) {
          col_data <- X_all_dmrs[, col]
          if (any(is.na(col_data))) {
            median_val <- median(col_data, na.rm = TRUE)
            if (is.na(median_val)) median_val <- 50
            X_all_dmrs[is.na(X_all_dmrs[, col]), col] <- median_val
          }
        }

        # Verify no NAs remain
        if (any(is.na(X_all_dmrs))) {
          log_msg("Warning: NAs still present, replacing with global median")
          X_all_dmrs[is.na(X_all_dmrs)] <- 50
        }

        # Convert to data frame and set proper column names
        X_all_dmrs <- as.data.frame(X_all_dmrs)
        colnames(X_all_dmrs) <- dmrs$DMR_ID

        # Final verification
        log_msg(paste("Final matrix dimensions:", nrow(X_all_dmrs), "samples x", ncol(X_all_dmrs), "DMRs"))
        if (nrow(X_all_dmrs) != length(meth) || ncol(X_all_dmrs) != nrow(dmrs)) {
          stop(paste("ERROR: Dimension mismatch. Expected:", length(meth), "samples x", nrow(dmrs), "DMRs"))
        }
        # Convert to data frame and set proper column names
        X_all_dmrs <- as.data.frame(X_all_dmrs)
        colnames(X_all_dmrs) <- dmrs$DMR_ID

        precompute_end <- Sys.time()
        precompute_time <- round(as.numeric(precompute_end - precompute_start, units = "secs"), 1)

        log_analysis_step(
          "Pre-computed methylation matrix",
          paste(nrow(X_all_dmrs), "samples x", ncol(X_all_dmrs), "DMRs in", precompute_time, "seconds")
        )

        # Save the pre-computed matrix to organized directory
        precomputed_matrix_file <- file.path(
          output_dirs$precomputed_data,
          paste0("precomputed_methylation_matrix_", reg_label, ".rds")
        )
        saveRDS(X_all_dmrs, precomputed_matrix_file)
        log_file_operation(
          "Saved precomputed matrix", precomputed_matrix_file,
          paste(nrow(X_all_dmrs), "samples x", ncol(X_all_dmrs), "DMRs")
        )

        # Create and save DMR summary
        dmr_summary <- data.frame(
          DMR_ID = dmrs$DMR_ID,
          rank = 1:nrow(dmrs),
          scaffold = dmrs$chr,
          start = dmrs$start,
          end = dmrs$end,
          strand = dmrs$strand,
          width_bp = dmrs$end - dmrs$start + 1,
          meth_difference = dmrs$meth.diff,
          qvalue = dmrs$qvalue,
          pvalue = dmrs$pvalue,
          direction = ifelse(dmrs$meth.diff > 0, "Hypermethylated", "Hypomethylated"),
          stringsAsFactors = FALSE
        )

        dmr_summary_file <- file.path(
          output_dirs$precomputed_data,
          paste0("dmr_summary_", reg_label, ".csv")
        )
        write.csv(dmr_summary, dmr_summary_file, row.names = FALSE)
        log_file_operation("Saved DMR summary", dmr_summary_file)

        # Create complete DMR dataset (all DMRs regardless of significance)
        log_analysis_step("Saving complete DMR dataset", "All detected DMRs")

        # Get all differential methylation results (before q-value filtering)
        if (is.data.frame(diff)) {
          all_dmrs_complete <- diff
        } else {
          all_dmrs_complete <- getData(diff) # This gets ALL results from calculateDiffMeth
        }

        all_dmrs_complete$DMR_ID <- paste0("ALL_DMR_", 1:nrow(all_dmrs_complete))
        all_dmrs_complete$rank_by_qvalue <- rank(all_dmrs_complete$qvalue, ties.method = "min")
        all_dmrs_complete$rank_by_pvalue <- rank(all_dmrs_complete$pvalue, ties.method = "min")
        all_dmrs_complete$rank_by_abs_meth_diff <- rank(-abs(all_dmrs_complete$meth.diff), ties.method = "min")

        # Add significance categories
        all_dmrs_complete$significance_category <- ifelse(
          all_dmrs_complete$qvalue <= 0.001, "Highly Significant (q<=0.001)",
          ifelse(all_dmrs_complete$qvalue <= 0.01, "Significant (q<=0.01)",
            ifelse(all_dmrs_complete$qvalue <= 0.05, "Moderately Significant (q<=0.05)",
              ifelse(all_dmrs_complete$qvalue <= 0.1, "Marginally Significant (q<=0.1)",
                "Not Significant (q>0.1)"
              )
            )
          )
        )

        all_dmrs_complete$passes_current_threshold <- all_dmrs_complete$qvalue <= opt$qval
        all_dmrs_complete$passes_meth_diff_threshold <- abs(all_dmrs_complete$meth.diff) >= 2

        # Create comprehensive summary
        all_dmrs_summary <- data.frame(
          region_file = reg_label,
          analysis_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          qvalue_threshold_used = opt$qval,
          DMR_ID = all_dmrs_complete$DMR_ID,
          rank_by_qvalue = all_dmrs_complete$rank_by_qvalue,
          rank_by_pvalue = all_dmrs_complete$rank_by_pvalue,
          rank_by_abs_meth_diff = all_dmrs_complete$rank_by_abs_meth_diff,
          scaffold = all_dmrs_complete$chr,
          start = all_dmrs_complete$start,
          end = all_dmrs_complete$end,
          strand = all_dmrs_complete$strand,
          width_bp = all_dmrs_complete$end - all_dmrs_complete$start + 1,
          meth_difference = all_dmrs_complete$meth.diff,
          pvalue = all_dmrs_complete$pvalue,
          qvalue = all_dmrs_complete$qvalue,
          significance_category = all_dmrs_complete$significance_category,
          passes_current_qvalue_threshold = all_dmrs_complete$passes_current_threshold,
          passes_meth_diff_threshold = all_dmrs_complete$passes_meth_diff_threshold,
          direction = ifelse(all_dmrs_complete$meth.diff > 0, "Hypermethylated", "Hypomethylated"),
          abs_meth_difference = abs(all_dmrs_complete$meth.diff),
          stringsAsFactors = FALSE
        )

        # Sort by q-value (most significant first)
        all_dmrs_summary <- all_dmrs_summary[order(all_dmrs_summary$qvalue, -all_dmrs_summary$abs_meth_difference), ]

        # Save complete DMR dataset to organized directory
        all_dmrs_file <- file.path(
          output_dirs$complete_dmr_datasets,
          paste0("ALL_DMRs_complete_dataset_", reg_label, ".csv")
        )
        write.csv(all_dmrs_summary, all_dmrs_file, row.names = FALSE)
        log_file_operation(
          "Saved complete DMR dataset", all_dmrs_file,
          paste(nrow(all_dmrs_summary), "total DMRs")
        )

        log_analysis_step(
          "DMR dataset statistics",
          paste(
            "Total:", nrow(all_dmrs_summary),
            "| Pass q-val:", sum(all_dmrs_summary$passes_current_qvalue_threshold),
            "| Pass meth-diff:", sum(all_dmrs_summary$passes_meth_diff_threshold)
          )
        )

        # Clustering analysis (if enabled)
        clustering_results <- NULL
        if (opt$perform_clustering && nrow(dmrs) >= 3) {
          clustering_results <- perform_dmr_clustering(
            X_all_dmrs, dmrs, treatment_info, output_dirs$clustering_analysis, reg_label
          )

          if (!is.null(clustering_results)) {
            dmrs <- clustering_results$dmrs_with_clusters
            log_analysis_step("DMR clustering complete", "Cluster assignments added")
          }
        }
        # Add this code AFTER the DMRs are calculated and sorted, BEFORE the machine learning loop
        # This should go right after this line: dmrs$DMR_ID <- paste0("DMR_", 1:nrow(dmrs))

        # Extract methylation data
        # verbose_msg("Extracting methylation features for all DMRs...")
        # X_all_dmrs <- regionStats(meth, all_dmrs_gr, column = "perc.meth")

        # Handle missing values
        # for (col in 1:ncol(X_all_dmrs)) {
        #  col_data <- X_all_dmrs[, col]
        #  if (any(is.na(col_data))) {
        #    median_val <- median(col_data, na.rm = TRUE)
        #    if (is.na(median_val)) median_val <- 50
        #    X_all_dmrs[is.na(X_all_dmrs[, col]), col] <- median_val
        #  }
        # }

        # Convert to data frame
        # X_all_dmrs <- as.data.frame(X_all_dmrs)
        # colnames(X_all_dmrs) <- dmrs$DMR_ID

        precompute_end <- Sys.time()
        log_msg(paste("Pre-computed matrix:", nrow(X_all_dmrs), "samples x", ncol(X_all_dmrs), "DMRs"))

        # Save the pre-computed matrix to disk for future use
        precomputed_matrix_file <- file.path(opt$output_dir, paste0("precomputed_methylation_matrix_", reg_label, ".rds"))
        saveRDS(X_all_dmrs, precomputed_matrix_file)

        precompute_end <- Sys.time()
        precompute_time <- round(as.numeric(precompute_end - precompute_start, units = "secs"), 1)

        log_msg(paste("[OK] Pre-computed methylation matrix:", nrow(X_all_dmrs), "samples x", ncol(X_all_dmrs), "DMRs in", precompute_time, "seconds"))
        log_msg(paste("[SAVE] Saved matrix to:", basename(precomputed_matrix_file)))

        # Also save a summary file with DMR information
        dmr_summary_file <- file.path(opt$output_dir, paste0("dmr_summary_", reg_label, ".csv"))

        dmr_summary <- data.frame(
          DMR_ID = dmrs$DMR_ID,
          rank = 1:nrow(dmrs),
          scaffold = dmrs$chr,
          start = dmrs$start,
          end = dmrs$end,
          width_bp = dmrs$end - dmrs$start + 1,
          meth_difference = dmrs$meth.diff,
          qvalue = dmrs$qvalue,
          pvalue = dmrs$pvalue,
          direction = ifelse(dmrs$meth.diff > 0, "Hypermethylated", "Hypomethylated"),
          stringsAsFactors = FALSE
        )
        write.csv(dmr_summary, dmr_summary_file, row.names = FALSE)
        log_msg(paste("[SUMMARY] DMR summary saved to:", basename(dmr_summary_file)))

        # NEW: Save ALL DMRs (not just significant ones) for record keeping
        log_msg("[SUMMARY] Saving complete DMR dataset (all detected DMRs regardless of significance)...")

        # Get all differential methylation results (before q-value filtering)
        if (is.data.frame(diff)) {
          all_dmrs_complete <- diff
        } else {
          all_dmrs_complete <- getData(diff) # This gets ALL results from calculateDiffMeth
        }

        all_dmrs_complete$DMR_ID <- paste0("ALL_DMR_", 1:nrow(all_dmrs_complete))
        all_dmrs_complete$rank_by_qvalue <- rank(all_dmrs_complete$qvalue, ties.method = "min")
        all_dmrs_complete$rank_by_pvalue <- rank(all_dmrs_complete$pvalue, ties.method = "min")
        all_dmrs_complete$rank_by_abs_meth_diff <- rank(-abs(all_dmrs_complete$meth.diff), ties.method = "min")

        # Add significance categories
        all_dmrs_complete$significance_category <- ifelse(
          all_dmrs_complete$qvalue <= 0.001, "Highly Significant (q<=0.001)",
          ifelse(all_dmrs_complete$qvalue <= 0.01, "Significant (q<=0.01)",
            ifelse(all_dmrs_complete$qvalue <= 0.05, "Moderately Significant (q<=0.05)",
              ifelse(all_dmrs_complete$qvalue <= 0.1, "Marginally Significant (q<=0.1)",
                "Not Significant (q>0.1)"
              )
            )
          )
        )

        all_dmrs_complete$passes_current_threshold <- all_dmrs_complete$qvalue <= opt$qval
        all_dmrs_complete$passes_meth_diff_threshold <- abs(all_dmrs_complete$meth.diff) >= 2

        # Create comprehensive summary
        all_dmrs_summary <- data.frame(
          region_file = reg_label,
          analysis_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          qvalue_threshold_used = opt$qval,
          DMR_ID = all_dmrs_complete$DMR_ID,
          rank_by_qvalue = all_dmrs_complete$rank_by_qvalue,
          rank_by_pvalue = all_dmrs_complete$rank_by_pvalue,
          rank_by_abs_meth_diff = all_dmrs_complete$rank_by_abs_meth_diff,
          scaffold = all_dmrs_complete$chr,
          start = all_dmrs_complete$start,
          end = all_dmrs_complete$end,
          width_bp = all_dmrs_complete$end - all_dmrs_complete$start + 1,
          meth_difference = all_dmrs_complete$meth.diff,
          pvalue = all_dmrs_complete$pvalue,
          qvalue = all_dmrs_complete$qvalue,
          significance_category = all_dmrs_complete$significance_category,
          passes_current_qvalue_threshold = all_dmrs_complete$passes_current_threshold,
          passes_meth_diff_threshold = all_dmrs_complete$passes_meth_diff_threshold,
          direction = ifelse(all_dmrs_complete$meth.diff > 0, "Hypermethylated", "Hypomethylated"),
          abs_meth_difference = abs(all_dmrs_complete$meth.diff),
          stringsAsFactors = FALSE
        )

        # Sort by q-value (most significant first)
        all_dmrs_summary <- all_dmrs_summary[order(all_dmrs_summary$qvalue, -all_dmrs_summary$abs_meth_difference), ]

        # Save complete DMR dataset
        all_dmrs_file <- file.path(opt$output_dir, paste0("ALL_DMRs_complete_dataset_", reg_label, ".csv"))
        write.csv(all_dmrs_summary, all_dmrs_file, row.names = FALSE)

        log_msg(paste("[SUMMARY] Complete DMR dataset saved to:", basename(all_dmrs_file)))
        log_msg(paste("[DATA] Total DMRs in complete dataset:", nrow(all_dmrs_summary)))
        log_msg(paste(
          "[DATA] DMRs passing q-value threshold (<=", opt$qval, "):",
          sum(all_dmrs_summary$passes_current_qvalue_threshold)
        ))
        log_msg(paste(
          "[DATA] DMRs passing meth diff threshold (>=2%):",
          sum(all_dmrs_summary$passes_meth_diff_threshold)
        ))
        log_msg(paste(
          "[DATA] DMRs passing both thresholds:",
          sum(all_dmrs_summary$passes_current_qvalue_threshold &
            all_dmrs_summary$passes_meth_diff_threshold)
        ))

        # ===== CLUSTERING CODE (only for significant DMRs, not all) =====
        # Perform clustering analysis if requested
        clustering_results <- NULL
        if (opt$perform_clustering && nrow(dmrs) >= 3) {
          log_msg("[ANALYSIS] Performing clustering on significant DMRs only...")

          clustering_results <- perform_dmr_clustering(
            X_all_dmrs, dmrs, treatment_info, opt$output_dir, reg_label
          )

          if (!is.null(clustering_results)) {
            dmrs <- clustering_results$dmrs_with_clusters
            log_msg(paste("[OK] DMRs now include cluster assignments (", nrow(dmrs), "significant DMRs only)"))
          }
        }
        # ===== END CLUSTERING CODE =====
      } # End of if (!goto_ml_training) block

      # Machine learning for each marker set - TRAIN FIRST, CHECK ACCURACY, THEN SAVE/PLOT

      # Either we just loaded precomputed data (force_ml) or we calculated it fresh
      # In both cases, we now have X_all_dmrs and dmrs ready for ML

      if (goto_ml_training) {
        log_msg("[START] Skipping to ML training using precomputed data")
      } else {
        log_msg("[START] Proceeding to ML training using freshly calculated data")
      }

      # ===== FIXED VERSION =====
      # CREATE ml_data BEFORE THE MARKER LOOP
      log_msg("[SETUP] Creating ml_data object for machine learning...")

      # Get treatment labels that match the sample count in X_all_dmrs
      if (exists("treatment_processed") && !is.null(treatment_processed)) {
        y <- factor(treatment_processed$labels[1:nrow(X_all_dmrs)])
      } else {
        y <- factor(treatment_info[1:nrow(X_all_dmrs)])
      }

      # ===== FILTER TO KEEP ONLY SELECTED GROUPS =====
      log_msg(paste("[FILTER] Initial sample count:", length(y)))
      log_msg(paste("[FILTER] Class distribution:", paste(table(y), collapse = ", ")))

      # If user specified specific groups via --group_names, filter to those
      if (!is.null(opt$group_names) && opt$group_names != "") {
        selected_groups <- trimws(strsplit(opt$group_names, ",")[[1]])
        log_msg(paste("[FILTER] Filtering to selected groups:", paste(selected_groups, collapse = ", ")))

        keep_indices <- y %in% selected_groups
        X_all_dmrs <- X_all_dmrs[keep_indices, , drop = FALSE]
        y <- y[keep_indices]
        y <- droplevels(y)

        log_msg(paste("[OK] After group filtering:", length(y), "samples"))
        log_msg(paste("[OK] Groups:", paste(levels(y), collapse = ", ")))
        log_msg(paste("[OK] Distribution:", paste(table(y), collapse = ", ")))
      }

      # CRITICAL: Remove filtered classes from ml_data BEFORE entering marker loop
      if (exists("insufficient_classes") && length(insufficient_classes) > 0) {
        log_msg(paste("[FIX] Removing filtered classes from ml_data:", paste(insufficient_classes, collapse = ", ")))

        # Filter X_all_dmrs to match
        keep_indices <- !y %in% insufficient_classes
        X_all_dmrs <- X_all_dmrs[keep_indices, , drop = FALSE]
        y <- y[keep_indices]
        y <- droplevels(y) # Remove unused factor levels

        log_msg(paste("[OK] Final ml_data after filtering:", length(y), "samples in", length(levels(y)), "classes"))
      }

      # Create the ml_data object (ALL DMRs, not just top N)
      ml_data <- bind_cols(X_all_dmrs, label = y)

      # Verify dimensions match
      if (length(y) != nrow(X_all_dmrs)) {
        log_error(paste(
          "Treatment length mismatch - X_all_dmrs has", nrow(X_all_dmrs),
          "samples but y has", length(y), "labels"
        ))
        next # Skip this region
      }

      # Create the ml_data object (ALL DMRs, not just top N)
      ml_data <- bind_cols(X_all_dmrs, label = y)

      # ===== ADD THIS DEBUG BLOCK =====
      log_msg("===== DEBUG: ml_data after creation =====")
      log_msg(paste("Dimensions:", nrow(ml_data), "x", ncol(ml_data)))
      log_msg(paste("Classes:", paste(table(ml_data$label), collapse = ", ")))

      # Check for NA rows
      na_per_row <- rowSums(is.na(ml_data[, -ncol(ml_data)]))
      log_msg(paste("Samples with ANY NAs:", sum(na_per_row > 0)))

      # Check NA distribution by class
      for (class_name in levels(ml_data$label)) {
        class_rows <- ml_data$label == class_name
        class_na_rows <- na_per_row[class_rows]
        log_msg(paste(
          "  Class", class_name, "- Total:", sum(class_rows),
          "| With NAs:", sum(class_na_rows > 0)
        ))
      }
      log_msg("=========================================")
      # ===== END DEBUG BLOCK =====

      # Verify it was created
      if (!exists("ml_data") || is.null(ml_data) || nrow(ml_data) == 0) {
        log_error("Failed to create ml_data object - skipping this region")
        next
      }

      log_msg(paste("[OK] Created ml_data:", nrow(ml_data), "samples x", ncol(ml_data), "columns"))
      log_msg(paste("[DATA] Class distribution:", paste(table(ml_data$label), collapse = ", ")))

      # Handle any missing values in the full dataset
      for (col in 1:(ncol(ml_data) - 1)) { # Skip the label column
        col_data <- ml_data[[col]]
        if (any(is.na(col_data))) {
          median_val <- median(col_data, na.rm = TRUE)
          if (is.na(median_val)) median_val <- 50
          ml_data[is.na(ml_data[[col]]), col] <- median_val
        }
      }

      # Remove rows with NA labels
      ml_data <- ml_data[!is.na(ml_data$label), ]
      ml_data$label <- droplevels(ml_data$label)

      log_msg(paste("[OK] ml_data ready for ML training"))
      # ===== END OF NEW BLOCK =====

      # Machine learning for each marker set - TRAIN FIRST, CHECK ACCURACY, THEN SAVE/PLOT
      model_results_for_region <- list()
      model_results_for_region <- list()
      # Create marker progress bar
      cat("\n")
      marker_pb <- create_progress_bar(length(marker_seq_region), prefix = "[DATA] Markers Progress:")

      # Check if we have enough samples per class for train/test split


      for (marker_idx in seq_along(marker_seq_region)) {
        n_markers <- marker_seq_region[marker_idx]

        # Update marker progress bar
        marker_pb <- update_progress_bar(
          marker_pb, marker_idx,
          paste("Processing", n_markers, "markers")
        )

        # LOG MARKER START
        log_msg(paste("[START] Processing", n_markers, "markers (", marker_idx, "of", length(marker_seq_region), ")"))
        flush.console()

        top_dmrs <- dmrs[1:n_markers, ]

        # Extract features for this marker set from the pre-existing ml_data
        verbose_msg("Extracting feature subset from ml_data...")
        dmr_indices <- 1:n_markers

        # Select only the DMR columns we need (ml_data has all DMRs + label column)
        feature_cols <- colnames(ml_data)[dmr_indices]
        ml_data_subset <- ml_data[, c(feature_cols, "label")]

        # ===== DEBUG: Check ml_data_subset =====
        log_msg(paste("[DEBUG] Marker", n_markers, "- ml_data_subset created"))
        log_msg(paste("  Dimensions:", nrow(ml_data_subset), "x", ncol(ml_data_subset)))
        log_msg(paste("  Classes:", paste(table(ml_data_subset$label), collapse = ", ")))

        # Check for NAs in this subset
        na_check <- sapply(1:nrow(ml_data_subset), function(i) {
          any(is.na(ml_data_subset[i, -ncol(ml_data_subset)]))
        })
        if (any(na_check)) {
          log_msg(paste("  [WARN] Samples with NAs:", sum(na_check)))
          # Show which classes have NAs
          for (class_name in levels(ml_data_subset$label)) {
            class_na <- sum(na_check[ml_data_subset$label == class_name])
            if (class_na > 0) {
              log_msg(paste("    Class", class_name, ":", class_na, "samples with NAs"))
            }
          }
        } else {
          log_msg("  No NAs detected")
        }
        # ===== END DEBUG =====

        # Handle NAs for THIS marker subset only
        na_rows <- apply(ml_data_subset[, -ncol(ml_data_subset)], 1, function(x) any(is.na(x)))
        if (any(na_rows)) {
          log_msg(paste("[WARN] Found", sum(na_rows), "samples with NA in", n_markers, "markers"))

          # Check class distribution before removal
          na_class_table <- table(ml_data_subset$label[na_rows])
          log_msg(paste("       NA samples by class:", paste(names(na_class_table), "=", na_class_table, collapse = ", ")))

          # Remove NA rows
          ml_data_subset <- ml_data_subset[!na_rows, ]

          # Check class distribution after removal
          remaining_classes <- table(ml_data_subset$label)
          log_msg(paste("       Remaining samples:", paste(names(remaining_classes), "=", remaining_classes, collapse = ", ")))

          # Verify we still have valid data for ML
          if (length(remaining_classes) < 2) {
            log_msg(paste("[ERROR] Marker set", n_markers, "- NA removal eliminated a class. Skipping."))
            next
          }

          if (any(remaining_classes < 2)) {
            log_msg(paste("[ERROR] Marker set", n_markers, "- Some classes have <2 samples. Skipping."))
            next
          }
        }

        verbose_msg(paste(
          "Feature subset:", nrow(ml_data_subset), "samples x",
          (ncol(ml_data_subset) - 1), "features"
        ))
        verbose_msg(paste("Class distribution:", paste(table(ml_data_subset$label), collapse = ", ")))

        # Create X matrix for this marker set
        X <- ml_data_subset[, -ncol(ml_data_subset)]

        verbose_msg(paste("Dataset created:", nrow(ml_data_subset), "samples x", ncol(ml_data_subset), "columns"))

        # Get current treatment labels
        current_treatment <- ml_data_subset$label

        # CRITICAL: Create train/test split HERE for this specific marker set
        verbose_msg(paste("Applying stratified train/test split: test_ratio =", opt$test_ratio))

        set.seed(123)

        # Check for very small classes (< 4 samples)
        class_counts <- table(ml_data_subset$label)

        # ============== INSERT THE NEW CODE HERE ==============
        insufficient_classes <- character(0)
        for (class_name in names(class_counts)) {
          n_samples <- class_counts[class_name]
          expected_train <- floor(n_samples * (1 - opt$test_ratio))
          expected_test <- ceiling(n_samples * opt$test_ratio)

          if (expected_train < 1 || expected_test < 1) {
            insufficient_classes <- c(insufficient_classes, class_name)
            log_msg(paste(
              "[WARN] Class", class_name, ":", n_samples, "samples ->",
              expected_train, "train +", expected_test, "test (insufficient)"
            ))
          } else {
            log_msg(paste(
              "[OK] Class", class_name, ":", n_samples, "samples ->",
              expected_train, "train +", expected_test, "test"
            ))
          }
        }

        if (length(insufficient_classes) > 0) {
          log_msg(paste(
            "[WARN] Removing classes with insufficient samples for train/test split:",
            paste(insufficient_classes, collapse = ", ")
          ))

          # CRITICAL FIX: Filter ml_data_subset, not ml_data!
          ml_data_subset <- ml_data_subset[!ml_data_subset$label %in% insufficient_classes, ]
          ml_data_subset$label <- droplevels(ml_data_subset$label)

          # Verify we still have enough classes
          remaining_classes <- levels(ml_data_subset$label)
          if (length(remaining_classes) < 2) {
            log_msg("[ERROR] Not enough classes remaining after filtering - skipping this marker count")
            next # Skip to next marker count
          }

          log_msg(paste(
            "[OK] Proceeding with", length(remaining_classes), "classes:",
            paste(remaining_classes, collapse = ", ")
          ))
          log_msg(paste("   Class distribution:", paste(table(ml_data_subset$label), collapse = " / ")))
        }
        # ============== END OF NEW CODE ==============

        small_classes <- names(class_counts[class_counts < 4])

        if (length(small_classes) > 0) {
          log_msg(paste("[WARN] Small classes detected:", paste(small_classes, "=", class_counts[small_classes], collapse = ", ")))
          log_msg("   Using manual stratified split to ensure all classes in both sets")

          # Manual stratified split for small classes
          train_indices <- c()
          test_indices <- c()

          for (class_name in levels(ml_data_subset$label)) {
            class_idx <- which(ml_data_subset$label == class_name)
            n_class <- length(class_idx)

            n_train <- max(1, floor(n_class * train_ratio))
            n_test <- n_class - n_train

            # Ensure at least 1 in each split
            if (n_train < 1) n_train <- 1
            if (n_test < 1) n_test <- 1

            # Randomly sample within class
            set.seed(123 + match(class_name, levels(ml_data_subset$label)))
            train_class_idx <- sample(class_idx, n_train)
            test_class_idx <- setdiff(class_idx, train_class_idx)

            train_indices <- c(train_indices, train_class_idx)
            test_indices <- c(test_indices, test_class_idx)

            log_msg(paste(
              "   ", class_name, ":", length(train_class_idx), "train +",
              length(test_class_idx), "test"
            ))
          }

          train_data <- ml_data_subset[train_indices, ]
          test_data <- ml_data_subset[test_indices, ]
        } else {
          # Standard createDataPartition for larger classes
          train_indices <- createDataPartition(
            ml_data_subset$label,
            p = train_ratio,
            list = FALSE,
            times = 1
          )

          train_data <- ml_data_subset[train_indices, ]
          test_data <- ml_data_subset[-train_indices, ]
        }


        verbose_msg(paste("Training set:", nrow(train_data), "samples (", round(train_ratio * 100), "%)", sep = ""))
        verbose_msg(paste("Test set:", nrow(test_data), "samples (", round(opt$test_ratio * 100), "%)", sep = ""))


        # =========================================================================
        # BULLETPROOF TRAIN/TEST DATASET EXPORT WITH VALIDATION
        # =========================================================================
        log_msg("[VALIDATION] Exporting and validating train/test datasets...")

        # Create export directory
        traintest_export_dir <- file.path(output_dirs$train_test_splits)
        dir.create(traintest_export_dir, recursive = TRUE, showWarnings = FALSE)

        # STEP 1: Add unique identifiers to track samples
        train_data$INTERNAL_ROW_ID <- paste0("TRAIN_", rownames(train_data))
        test_data$INTERNAL_ROW_ID <- paste0("TEST_", rownames(test_data))

        # STEP 2: Export COMPLETE datasets (all features + labels + IDs)
        train_export_file <- file.path(
          traintest_export_dir,
          paste0("TRAIN_dataset_", reg_label, "_", n_markers, "markers.csv")
        )
        test_export_file <- file.path(
          traintest_export_dir,
          paste0("TEST_dataset_", reg_label, "_", n_markers, "markers.csv")
        )

        # Export with all information
        write.csv(train_data, train_export_file, row.names = TRUE)
        write.csv(test_data, test_export_file, row.names = TRUE)

        log_file_operation(
          "Exported TRAIN dataset", train_export_file,
          paste(nrow(train_data), "samples x", ncol(train_data), "features")
        )
        log_file_operation(
          "Exported TEST dataset", test_export_file,
          paste(nrow(test_data), "samples x", ncol(test_data), "features")
        )

        # STEP 3: CRITICAL VALIDATION CHECKS
        validation_passed <- TRUE
        validation_log <- character()

        # Check 1: No sample overlap between train and test
        train_ids <- train_data$INTERNAL_ROW_ID
        test_ids <- test_data$INTERNAL_ROW_ID
        overlapping_ids <- intersect(train_ids, test_ids)

        if (length(overlapping_ids) > 0) {
          validation_passed <- FALSE
          error_msg <- paste(
            "CRITICAL ERROR: Found", length(overlapping_ids),
            "samples in BOTH train and test sets!"
          )
          log_error(error_msg, "DATA LEAKAGE DETECTED")
          validation_log <- c(validation_log, error_msg)
          validation_log <- c(validation_log, paste(
            "Overlapping IDs:",
            paste(head(overlapping_ids, 10), collapse = ", ")
          ))
        } else {
          success_msg <- paste("[PASS] VALIDATION PASSED: Zero overlap between train and test sets")
          log_msg(success_msg)
          validation_log <- c(validation_log, success_msg)
        }

        # Check 2: All samples accounted for
        total_samples <- nrow(ml_data_subset)
        train_count <- nrow(train_data)
        test_count <- nrow(test_data)

        if ((train_count + test_count) != total_samples) {
          validation_passed <- FALSE
          error_msg <- paste(
            "CRITICAL ERROR: Sample count mismatch!",
            "Original:", total_samples,
            "Train:", train_count,
            "Test:", test_count,
            "Total:", train_count + test_count
          )
          log_error(error_msg, "SAMPLE COUNT MISMATCH")
          validation_log <- c(validation_log, error_msg)
        } else {
          success_msg <- paste(
            "[PASS] VALIDATION PASSED: All", total_samples, "samples accounted for",
            "(", train_count, "train +", test_count, "test )"
          )
          log_msg(success_msg)
          validation_log <- c(validation_log, success_msg)
        }

        # Check 3: Feature consistency
        train_features <- setdiff(colnames(train_data), c("label", "INTERNAL_ROW_ID"))
        test_features <- setdiff(colnames(test_data), c("label", "INTERNAL_ROW_ID"))

        if (!identical(train_features, test_features)) {
          validation_passed <- FALSE
          error_msg <- "CRITICAL ERROR: Train and test have different features!"
          log_error(error_msg, "FEATURE MISMATCH")
          validation_log <- c(validation_log, error_msg)
        } else {
          success_msg <- paste(
            "[PASS] VALIDATION PASSED: Train and test have identical",
            length(train_features), "features"
          )
          log_msg(success_msg)
          validation_log <- c(validation_log, success_msg)
        }

        # Check 4: Class distribution
        train_class_dist <- table(train_data$label)
        test_class_dist <- table(test_data$label)

        class_check_msg <- paste("[PASS] Train classes:", paste(names(train_class_dist), "=",
          train_class_dist,
          collapse = ", "
        ))
        log_msg(class_check_msg)
        validation_log <- c(validation_log, class_check_msg)

        class_check_msg <- paste("[PASS] Test classes:", paste(names(test_class_dist), "=",
          test_class_dist,
          collapse = ", "
        ))
        log_msg(class_check_msg)
        validation_log <- c(validation_log, class_check_msg)

        # Check 5: No duplicate rows within each set
        train_duplicates <- sum(duplicated(train_data[, train_features]))
        test_duplicates <- sum(duplicated(test_data[, test_features]))

        if (train_duplicates > 0 || test_duplicates > 0) {
          warning_msg <- paste(
            "WARNING: Found duplicates - Train:", train_duplicates,
            "Test:", test_duplicates
          )
          log_warning(warning_msg, "DUPLICATE CHECK")
          validation_log <- c(validation_log, warning_msg)
        } else {
          success_msg <- "[PASS] VALIDATION PASSED: No duplicate samples in train or test"
          log_msg(success_msg)
          validation_log <- c(validation_log, success_msg)
        }

        # STEP 4: Save validation report
        validation_report_file <- file.path(
          traintest_export_dir,
          paste0("VALIDATION_", reg_label, "_", n_markers, "markers.txt")
        )

        validation_report <- c(
          paste(rep("=", 80), collapse = ""),
          "TRAIN/TEST SPLIT VALIDATION REPORT",
          paste(rep("=", 80), collapse = ""),
          paste("Region:", reg_label),
          paste("Markers:", n_markers),
          paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
          "",
          "DATASET SUMMARY:",
          paste("  Original samples:", total_samples),
          paste("  Training samples:", train_count, paste0("(", round(100 * train_count / total_samples, 1), "%)")),
          paste("  Test samples:", test_count, paste0("(", round(100 * test_count / total_samples, 1), "%)")),
          paste("  Features per sample:", length(train_features)),
          "",
          "SAMPLE IDENTIFIERS:",
          "  Training set IDs:",
          paste("   ", head(train_ids, 10)),
          if (length(train_ids) > 10) paste("    ... and", length(train_ids) - 10, "more"),
          "",
          "  Test set IDs:",
          paste("   ", head(test_ids, 10)),
          if (length(test_ids) > 10) paste("    ... and", length(test_ids) - 10, "more"),
          "",
          "VALIDATION RESULTS:",
          validation_log,
          "",
          if (validation_passed) {
            c(
              "", paste(rep("=", 80), collapse = ""),
              "[PASS][PASS][PASS] ALL VALIDATION CHECKS PASSED [PASS][PASS][PASS]",
              "This split is SAFE to use for model training.",
              paste(rep("=", 80), collapse = "")
            )
          } else {
            c(
              "", paste(rep("=", 80), collapse = ""),
              "[FAIL][FAIL][FAIL] VALIDATION FAILED - DATA LEAKAGE DETECTED [FAIL][FAIL][FAIL]",
              "DO NOT USE THIS SPLIT FOR MODEL TRAINING!",
              paste(rep("=", 80), collapse = "")
            )
          },
          "",
          "FILES EXPORTED:",
          paste("  Training data:", basename(train_export_file)),
          paste("  Test data:", basename(test_export_file)),
          paste("  Validation report:", basename(validation_report_file)),
          "",
          paste(rep("=", 80), collapse = "")
        )

        writeLines(validation_report, validation_report_file)
        log_file_operation("Saved validation report", validation_report_file)

        # STEP 5: STOP EXECUTION IF VALIDATION FAILED
        if (!validation_passed) {
          stop(paste(
            "CRITICAL ERROR: Train/test split validation FAILED for", reg_label,
            "with", n_markers, "markers. Check validation report:",
            validation_report_file
          ))
        }

        # STEP 6: Remove INTERNAL_ROW_ID before training (but keep in exported files)
        train_data$INTERNAL_ROW_ID <- NULL
        test_data$INTERNAL_ROW_ID <- NULL

        log_msg("[OK] Train/test datasets exported and validated successfully!")

        # Simple sanity check - just verify we have data
        if (nrow(train_data) < 2 || nrow(test_data) < 1) {
          log_msg("[ERROR] Invalid train/test split - insufficient samples")
          next
        }

        log_msg(paste(
          "[OK] Ready for ML training:",
          nrow(train_data), "train +", nrow(test_data), "test"
        ))

        # =========================================================================

        # CRITICAL: Verify positive class is valid for ROC curves
        if (!is.null(opt$positive_class)) {
          if (!opt$positive_class %in% levels(ml_data$label)) {
            log_warning(paste(
              "Positive class '", opt$positive_class, "' not in labels. Available:",
              paste(levels(ml_data$label), collapse = ", ")
            ))
            log_msg("Setting positive_class to NULL - will use alphabetical default")
            opt$positive_class <- NULL
          } else {
            log_msg(paste("[TARGET] Positive class for ROC curves:", opt$positive_class))
          }
        }

        # Check if we have enough samples for machine learning
        if (nrow(ml_data) < 10) {
          log_msg("[WARN] Too few samples for machine learning")
          next
        }

        # CRITICAL FIX: Remove any rows with NA values
        na_rows <- apply(ml_data, 1, function(x) any(is.na(x)))
        if (any(na_rows)) {
          # Check which classes will be affected
          na_class_table <- table(ml_data$label[na_rows])
          log_msg(paste("[WARN] Found", sum(na_rows), "samples with NA values"))
          log_msg(paste("       Affected classes:", paste(names(na_class_table), "=", na_class_table, collapse = ", ")))

          # Remove NA rows
          ml_data <- ml_data[!na_rows, ]

          # Verify we still have at least 2 classes
          remaining_class_table <- table(ml_data$label)
          log_msg(paste("       After removal:", paste(names(remaining_class_table), "=", remaining_class_table, collapse = ", ")))

          if (length(remaining_class_table) < 2) {
            log_msg("[ERROR] NA removal eliminated entire class(es) - cannot proceed with ML")
            log_msg("       This region needs better data quality or different filtering parameters")
            next # Skip this region
          }

          if (any(remaining_class_table < 2)) {
            log_msg("[ERROR] Some classes have <2 samples after NA removal - cannot train")
            next # Skip this region
          }
        }

        # CRITICAL FIX: Ensure label is a factor with no NAs
        ml_data$label <- droplevels(ml_data$label)
        if (any(is.na(ml_data$label))) {
          log_msg("[WARN] Removing samples with NA labels")
          ml_data <- ml_data[!is.na(ml_data$label), ]

          # Verify again
          if (length(levels(ml_data$label)) < 2) {
            log_msg("[ERROR] Lost classes after NA label removal")
            next
          }
        }

        # Verify we still have enough samples
        if (nrow(ml_data) < 10) {
          log_msg("[WARN] Too few samples after NA removal")
          next
        }

        class_table <- table(ml_data$label)
        # ADD THIS RIGHT AFTER: class_table <- table(ml_data$label)
        log_msg(paste("DEBUG: Class table details:"))
        for (i in 1:length(class_table)) {
          log_msg(paste("   Class", names(class_table)[i], ":", class_table[i], "samples"))
        }
        min_class_size <- min(class_table)

        log_msg(paste("DEBUG: Min class size calculated as:", min_class_size))
        log_msg(paste("DEBUG: ml_data dimensions:", nrow(ml_data), "rows x", ncol(ml_data), "columns"))
        log_msg(paste("DEBUG: Label levels:", paste(levels(ml_data$label), collapse = ", ")))
        min_class_size <- min(class_table)
        if (min_class_size < 2) {
          log_msg("[WARN] Some classes have fewer than 2 samples - skipping")
          next
        }

        # Check if we have enough samples for SVM (needs at least 4 per class)
        current_model_list <- model_list
        if ("svm" %in% current_model_list && min_class_size < 4) {
          log_msg("[WARN] Removing SVM from model list - needs at least 4 samples per class")
          current_model_list <- current_model_list[current_model_list != "svm"]
          if (length(current_model_list) == 0) {
            log_msg("[WARN] No suitable models for this dataset size")
            next
          }
        }

        # Enhanced train/test split with predefined option
        if (!is.null(predefined_assignments)) {
          # Use predefined assignments
          verbose_msg("Applying predefined train/test split from Excel file...")

          # Get sample names from CURRENT ml_data_subset (not the full ml_data!)
          sample_names <- rownames(ml_data_subset)
          if (is.null(sample_names)) {
            sample_names <- paste0("sample_", 1:nrow(ml_data_subset))
            rownames(ml_data_subset) <- sample_names
          }

          # Apply predefined split to CURRENT subset
          split_result <- apply_predefined_split(ml_data_subset, predefined_assignments, sample_names, output_dirs)

          train_data <- split_result$train_data
          test_data <- split_result$test_data

          # Log detailed split information
          log_msg("[DATA] Predefined split applied:")
          log_msg(paste("   Matched samples:", split_result$matched_samples, "of", nrow(ml_data_subset)))
          log_msg(paste("   Training set:", nrow(train_data), "samples"))
          log_msg(paste("   Test set:", nrow(test_data), "samples"))

          if (length(split_result$unmatched_samples) > 0) {
            log_msg(paste("   Unmatched samples:", length(split_result$unmatched_samples), "(excluded)"))
          }
        } else {
          # Use random split (existing code)


          verbose_msg(paste("Training set:", nrow(train_data), "samples (", round(train_ratio * 100), "%)", sep = ""))
          verbose_msg(paste("Test set:", nrow(test_data), "samples (", round(opt$test_ratio * 100), "%)", sep = ""))
        }


        # ===== ADD THIS BLOCK HERE =====
        # Save train/test split for validation
        split_record_file <- file.path(
          output_dirs$train_test_splits,
          paste0("split_", reg_label, "_", n_markers, "markers.csv")
        )

        split_record <- data.frame(
          sample_id = rownames(ml_data_subset),
          actual_class = ml_data_subset$label,
          split_assignment = ifelse(1:nrow(ml_data_subset) %in%
            match(rownames(train_data), rownames(ml_data_subset)),
          "train", "test"
          ),
          marker_count = n_markers,
          stringsAsFactors = FALSE
        )

        write.csv(split_record, split_record_file, row.names = FALSE)
        log_file_operation(
          "Saved train/test split", split_record_file,
          paste(nrow(train_data), "train /", nrow(test_data), "test")
        )
        # ===== END NEW BLOCK =====

        # Set up cross-validation with reproducibility control
        if (opt$reproducible) {
          ctrl <- trainControl(
            method = "repeatedcv",
            number = cv_folds,
            repeats = 3,
            verboseIter = FALSE,
            classProbs = TRUE,
            allowParallel = FALSE
          )
        } else {
          ctrl <- trainControl(
            method = "repeatedcv",
            number = cv_folds,
            repeats = 3,
            verboseIter = FALSE,
            classProbs = TRUE,
            allowParallel = TRUE
          )
        }
        metric <- "Accuracy"

        # STEP 1: Train all models with classifier progress bar
        cat("\n")
        classifier_pb <- create_progress_bar(length(current_model_list),
          prefix = paste("[MODEL] Classifiers (", n_markers, " markers):")
        )

        marker_results <- list()
        max_accuracy <- 0

        for (model_idx in seq_along(current_model_list)) {
          model_name <- current_model_list[model_idx]

          # Update classifier progress bar with current model info
          classifier_pb <- update_progress_bar(
            classifier_pb, model_idx,
            paste("Training", model_name)
          )

          # FORCE FLUSH OUTPUT
          flush.console()

          if (!(model_name %in% names(model_specs))) {
            log_msg(paste("[ERROR] Unknown model:", model_name))
            next
          }

          model_start <- Sys.time()

          # Show intermediate progress for slow models
          if (model_name %in% c("nnet", "svm", "glmnet", "xgboost")) {
            cat(paste0("  [", model_name, " training in progress...]"))
            flush.console()
          }

          tryCatch(
            {
              # SUPPRESS ALL MODEL OUTPUT to prevent progress bar breaks

              # Model-specific configurations
              if (model_name == "svm") {
                # SVM-specific setup - adaptive for 2 or 3+ classes
                n_classes <- length(levels(train_data$label))

                if (n_classes == 2) {
                  svm_ctrl <- trainControl(
                    method = "cv", number = min(cv_folds, nrow(train_data)),
                    classProbs = TRUE, summaryFunction = twoClassSummary,
                    savePredictions = "final", allowParallel = !opt$reproducible,
                    verboseIter = FALSE # CRITICAL: Suppress verbose output
                  )
                  svm_metric <- "ROC"
                } else {
                  svm_ctrl <- trainControl(
                    method = "cv", number = min(cv_folds, nrow(train_data)),
                    classProbs = TRUE, summaryFunction = multiClassSummary,
                    savePredictions = "final", allowParallel = !opt$reproducible,
                    verboseIter = FALSE # CRITICAL: Suppress verbose output
                  )
                  svm_metric <- "Accuracy"
                }

                svm_grid <- expand.grid(sigma = c(0.01, 0.1), C = c(0.1, 1))

                model_fit <- train(
                  label ~ .,
                  data = train_data,
                  method = model_specs[[model_name]],
                  trControl = svm_ctrl, tuneGrid = svm_grid,
                  preProcess = c("center", "scale"),
                  metric = svm_metric, verbose = FALSE, allowParallel = !opt$reproducible
                )
              } else if (model_name == "xgboost") {
                # XGBoost - suppress output
                xgb_ctrl <- trainControl(
                  method = "cv", number = cv_folds,
                  classProbs = TRUE, summaryFunction = multiClassSummary,
                  savePredictions = "final", allowParallel = !opt$reproducible,
                  verboseIter = FALSE # CRITICAL
                )

                xgb_grid <- expand.grid(
                  nrounds = c(50, 100), max_depth = c(3, 6),
                  eta = 0.1, gamma = 0, colsample_bytree = 1,
                  min_child_weight = 1, subsample = 1
                )

                model_fit <- train(
                  label ~ .,
                  data = train_data,
                  method = model_specs[[model_name]],
                  trControl = xgb_ctrl, tuneGrid = xgb_grid,
                  metric = "Accuracy", verbose = 0 # CRITICAL: verbose = 0
                )
              } else if (model_name == "nb") {
                # Naive Bayes
                set.seed(7)
                Grid <- expand.grid(usekernel = TRUE, adjust = 1, fL = c(0.2, 0.5, 0.8))
                model_fit <- train(
                  label ~ .,
                  data = train_data,
                  method = "nb", metric = metric,
                  trControl = ctrl, tuneGrid = Grid
                )
              } else if (model_name == "logreg") {
                # Logistic Regression
                model_fit <- train(
                  label ~ .,
                  data = train_data,
                  method = model_specs[[model_name]],
                  trControl = ctrl, family = "binomial",
                  metric = "Accuracy", verbose = FALSE
                )
              } else if (model_name == "nnet") {
                # Neural Network - special handling for trace parameter
                nnet_grid <- expand.grid(size = c(3, 5), decay = c(0, 0.1))

                # Create a special control that passes trace = FALSE through to nnet
                ctrl_nnet <- trainControl(
                  method = "cv",
                  number = cv_folds,
                  classProbs = TRUE,
                  allowParallel = !opt$reproducible,
                  verboseIter = FALSE
                )

                model_fit <- train( # <- FIXED: Use model_fit
                  label ~ .,
                  data = train_data, # <- FIXED: Use train_data
                  method = "nnet",
                  trControl = ctrl_nnet,
                  tuneGrid = nnet_grid,
                  preProcess = c("center", "scale"),
                  metric = "Accuracy",
                  trace = FALSE,
                  linout = FALSE,
                  maxit = 100
                )
              } else if (model_name == "ranger") {
                # Ranger (Fast Random Forest)
                n_features <- ncol(train_data) - 1
                max_mtry <- n_features
                mtry_values <- unique(sort(c(
                  max(1, floor(sqrt(n_features))),
                  max(1, floor(n_features / 3)),
                  max(1, floor(n_features / 2))
                )))
                mtry_values <- mtry_values[mtry_values <= max_mtry]

                ranger_grid <- expand.grid(
                  mtry = mtry_values,
                  splitrule = c("gini", "extratrees"),
                  min.node.size = c(1, 5)
                )

                model_fit <- train(
                  label ~ .,
                  data = train_data,
                  method = "ranger", trControl = ctrl,
                  tuneGrid = ranger_grid, metric = "Accuracy",
                  verbose = FALSE, importance = "impurity"
                )
              } else if (model_name == "glmnet") {
                # Glmnet
                if (opt$reproducible) {
                  set.seed(123)
                  seeds <- vector(mode = "list", length = 31)
                  for (i in 1:30) seeds[[i]] <- sample.int(1000, 100)
                  seeds[[31]] <- sample.int(1000, 1)

                  ctrl_reproducible <- trainControl(
                    method = "repeatedcv", number = cv_folds, repeats = 3,
                    classProbs = TRUE, summaryFunction = defaultSummary,
                    savePredictions = "none", seeds = seeds,
                    allowParallel = !opt$reproducible, verboseIter = FALSE # CRITICAL
                  )

                  set.seed(123)
                  model_fit <- train(
                    label ~ .,
                    data = train_data,
                    method = "glmnet", trControl = ctrl_reproducible,
                    preProcess = c("center", "scale"),
                    metric = "Accuracy", tuneLength = 10,
                    verbose = FALSE, allowParallel = !opt$reproducible
                  )
                } else {
                  set.seed(123)
                  ctrl_fast <- trainControl(
                    method = "repeatedcv", number = cv_folds, repeats = 3,
                    classProbs = TRUE, allowParallel = TRUE,
                    verboseIter = FALSE # CRITICAL
                  )

                  model_fit <- train(
                    label ~ .,
                    data = train_data,
                    method = "glmnet", trControl = ctrl_fast,
                    preProcess = c("center", "scale"),
                    metric = "Accuracy", tuneLength = 10,
                    verbose = FALSE
                  )
                }
              } else {
                # Default training for other models
                if (model_name %in% c("knn", "lda")) {
                  model_fit <- train(
                    label ~ .,
                    data = train_data,
                    method = model_specs[[model_name]],
                    trControl = ctrl,
                    preProcess = c("center", "scale", "zv"),
                    metric = "Accuracy", tuneLength = 3
                  )
                } else {
                  model_fit <- train(
                    label ~ .,
                    data = train_data,
                    method = model_specs[[model_name]],
                    trControl = ctrl,
                    preProcess = c("center", "scale", "zv"),
                    metric = "Accuracy", tuneLength = 3,
                    verbose = FALSE
                  )
                }
              }

              # Make predictions
              pred_class <- predict(model_fit, test_data)
              pred_prob <- predict(model_fit, test_data, type = "prob")
              actual_labels <- test_data$label

              # ADD THESE LINES:
              classes <- levels(actual_labels)
              n_classes <- length(classes)

              # Calculate basic accuracy
              accuracy <- sum(pred_class == actual_labels) / length(actual_labels)

              model_end <- Sys.time()
              model_time <- round(as.numeric(model_end - model_start, units = "secs"))

              # Store results
              marker_results[[model_name]] <- list(
                model_fit = model_fit,
                pred_class = pred_class,
                pred_prob = pred_prob,
                actual_labels = actual_labels,
                accuracy = accuracy,
                model_time = model_time,
                current_dmrs = top_dmrs,
                current_methylation_matrix = X
              )

              max_accuracy <- max(max_accuracy, accuracy)

              # LOG COMPLETION
              log_msg(paste("[COMPLETE]", model_name, "finished in", model_time, "seconds | Acc:", round(accuracy, 3)))

              # REPRINT progress bar after model completes
              classifier_pb <- update_progress_bar(
                classifier_pb, model_idx,
                paste(model_name, "- Acc:", round(accuracy, 3))
              )

              # FORCE OUTPUT
              flush.console()
            },
            error = function(e) {
              log_msg(paste("[ERROR] Model training failed:", e$message))
              classifier_pb <- update_progress_bar(
                classifier_pb, model_idx,
                paste(model_name, "- FAILED")
              )
            }
          )
        }

        # Complete classifier progress bar
        if (classifier_pb$current < classifier_pb$total) {
          classifier_pb <- update_progress_bar(classifier_pb, classifier_pb$total, "Complete")
        }
        cat("\n")

        # Show results summary for this marker set
        if (length(marker_results) > 0) {
          log_msg(paste("[DATA] Results for", n_markers, "markers:"))
          for (model_name in names(marker_results)) {
            result <- marker_results[[model_name]]
            status_icon <- if (result$accuracy >= opt$min_accuracy) "[OK]" else "[ERROR]"
            log_msg(sprintf(
              "   %s %s: Acc=%.3f (Time=%ds)",
              status_icon, model_name, result$accuracy, result$model_time
            ))
          }
        }

        # STEP 2: Check if ANY model meets accuracy threshold
        if (max_accuracy >= opt$min_accuracy) {
          log_msg(paste(
            "[SUCCESS] At least one model achieved accuracy >=", opt$min_accuracy,
            "(best:", round(max_accuracy, 3), ") - saving results and creating plots"
          ))

          # Save results for models that meet threshold
          for (model_name in names(marker_results)) {
            result <- marker_results[[model_name]]

            # Add this debug code before plotting functions

            if (result$accuracy >= opt$min_accuracy) {
              log_msg(paste("[SAVE] Saving", model_name, "model (accuracy:", round(result$accuracy, 3), ")"))

              # Create model directory first
              model_dir <- file.path(opt$output_dir, reg_label, model_name, paste0("markers_", n_markers))
              dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)


              # Safety check: ensure metrics vectors exist and have correct length
              if (!exists("sensitivity") || length(sensitivity) == 0) {
                sensitivity <- rep(0, n_classes)
              }
              if (!exists("specificity") || length(specificity) == 0) {
                specificity <- rep(0, n_classes)
              }
              if (!exists("f1_score") || length(f1_score) == 0) {
                f1_score <- rep(0, n_classes)
              }

              # REPLACE THE ENTIRE METRICS CALCULATION SECTION:
              # Create detailed DMR information FIRST (outside tryCatch)
              # Use actual column names from current_dmrs
              # Create detailed DMR information with error checking
              # ===== CALCULATE DIRECTION BASED ON POSITIVE CLASS (FIXED) =====
              # Determine which groups to compare
              group_names_ordered <- NULL
              if (exists("treatment_processed", envir = .GlobalEnv)) {
                tp <- get("treatment_processed", envir = .GlobalEnv)
                if (!is.null(tp) && !is.null(tp$group_names)) {
                  group_names_ordered <- tp$group_names
                }
              } else if (!is.null(opt$group_names)) {
                group_names_ordered <- trimws(strsplit(opt$group_names, ",")[[1]])
              } else {
                group_names_ordered <- levels(ml_data_subset$label)
              }

              # Determine positive and negative class (SAME logic as ROC curves)
              if (!is.null(opt$positive_class) && opt$positive_class %in% group_names_ordered) {
                pos_class <- opt$positive_class
                neg_class <- setdiff(group_names_ordered, opt$positive_class)
              } else {
                # Default: second group in order is positive
                pos_class <- group_names_ordered[2]
                neg_class <- group_names_ordered[1]
              }

              log_msg(paste("[DMR] Calculating direction relative to POSITIVE class:", pos_class))
              log_msg(paste("      Reference (negative) class:", neg_class))

              # Calculate ACTUAL methylation difference: positive_class - negative_class
              # Get mean methylation for each DMR in each group
              dmr_directions <- character(n_markers)
              actual_meth_diffs <- numeric(n_markers)

              for (i in 1:n_markers) {
                dmr_id <- result$current_dmrs$DMR_ID[i]

                # Get methylation values for this DMR
                if (dmr_id %in% colnames(ml_data_subset)) {
                  dmr_values <- ml_data_subset[[dmr_id]]
                  dmr_labels <- ml_data_subset$label

                  # Calculate mean for each group
                  pos_samples <- dmr_values[dmr_labels == pos_class]
                  neg_samples <- dmr_values[dmr_labels %in% neg_class]

                  pos_mean <- mean(pos_samples, na.rm = TRUE)
                  neg_mean <- mean(neg_samples, na.rm = TRUE)

                  # CORRECT direction: positive_class - negative_class
                  actual_diff <- pos_mean - neg_mean
                  actual_meth_diffs[i] <- actual_diff

                  # Direction from positive class perspective
                  if (actual_diff > 0) {
                    dmr_directions[i] <- paste0("Hypermethylated_in_", pos_class)
                  } else {
                    dmr_directions[i] <- paste0("Hypomethylated_in_", pos_class)
                  }
                } else {
                  dmr_directions[i] <- "Unknown"
                  actual_meth_diffs[i] <- NA
                }
              }

              # Log summary
              hyper_count <- sum(grepl("Hypermethylated", dmr_directions))
              hypo_count <- sum(grepl("Hypomethylated", dmr_directions))
              log_msg(paste("      Hypermethylated in", pos_class, ":", hyper_count, "DMRs"))
              log_msg(paste("      Hypomethylated in", pos_class, ":", hypo_count, "DMRs"))
              # ===== END DIRECTION CALCULATION =====

              # Create detailed DMR information FIRST (outside tryCatch)
              detailed_dmr_info <- tryCatch(
                {
                  data.frame(
                    region_file = rep(reg_label, n_markers),
                    model = rep(model_name, n_markers),
                    n_markers = rep(n_markers, n_markers),
                    accuracy = rep(result$accuracy, n_markers),
                    auc = rep(0, n_markers),
                    DMR_rank = 1:n_markers,
                    DMR_ID = result$current_dmrs$DMR_ID[1:n_markers],
                    scaffold = result$current_dmrs$chr[1:n_markers],
                    start = result$current_dmrs$start[1:n_markers],
                    end = result$current_dmrs$end[1:n_markers],
                    width_bp = result$current_dmrs$end[1:n_markers] - result$current_dmrs$start[1:n_markers] + 1,

                    # FIXED: Use calculated actual difference and direction
                    meth_difference_raw = result$current_dmrs$meth.diff[1:n_markers], # Original methylKit value
                    meth_difference_actual = actual_meth_diffs, # Corrected: pos_class - neg_class

                    qvalue = result$current_dmrs$qvalue[1:n_markers],
                    pvalue = result$current_dmrs$pvalue[1:n_markers],

                    # FIXED: Direction from positive class perspective
                    direction = dmr_directions,

                    # Add comparison info for transparency
                    positive_class = rep(pos_class, n_markers),
                    negative_class = rep(paste(neg_class, collapse = ", "), n_markers),
                    stringsAsFactors = FALSE
                  )
                },
                error = function(e) {
                  log_warning(paste("Could not create detailed DMR info:", e$message), "DMR info creation")
                  # Return minimal dataframe
                  data.frame(
                    region_file = reg_label,
                    model = model_name,
                    n_markers = n_markers,
                    accuracy = result$accuracy,
                    auc = 0,
                    stringsAsFactors = FALSE
                  )
                }
              )

              # THEN do metrics calculation in tryCatch
              tryCatch(
                {
                  # Calculate additional metrics - WRAPPED IN ERROR HANDLING
                  classes <- levels(result$actual_labels)
                  n_classes <- length(classes)

                  sensitivity <- rep(0, n_classes)
                  specificity <- rep(0, n_classes)
                  precision <- rep(0, n_classes)
                  f1_score <- rep(0, n_classes)
                  names(sensitivity) <- names(specificity) <- names(precision) <- names(f1_score) <- classes

                  for (i in seq_along(classes)) {
                    class_name <- classes[i]

                    tp <- sum(result$pred_class == class_name & result$actual_labels == class_name)
                    tn <- sum(result$pred_class != class_name & result$actual_labels != class_name)
                    fp <- sum(result$pred_class == class_name & result$actual_labels != class_name)
                    fn <- sum(result$pred_class != class_name & result$actual_labels == class_name)

                    sensitivity[i] <- if ((tp + fn) > 0) tp / (tp + fn) else 0
                    specificity[i] <- if ((tn + fp) > 0) tn / (tn + fp) else 0
                    precision[i] <- if ((tp + fp) > 0) tp / (tp + fp) else 0

                    if (sensitivity[i] + precision[i] > 0) {
                      f1_score[i] <- 2 * (precision[i] * sensitivity[i]) / (precision[i] + sensitivity[i])
                    } else {
                      f1_score[i] <- 0
                    }
                  }

                  # Calculate AUC with proper error handling and multi-class support
                  auc <- 0
                  tryCatch(
                    {
                      n_classes <- length(levels(result$actual_labels))

                      # Get ordered group names
                      group_names_ordered <- NULL
                      if (exists("treatment_processed", envir = .GlobalEnv)) {
                        tp <- get("treatment_processed", envir = .GlobalEnv)
                        group_names_ordered <- tp$group_names
                      } else if (!is.null(opt$group_names)) {
                        group_names_ordered <- trimws(strsplit(opt$group_names, ",")[[1]])
                      } else {
                        group_names_ordered <- levels(result$actual_labels)
                      }

                      if (n_classes == 2) {
                        # ===== BINARY CLASSIFICATION =====
                        # Determine positive and negative class
                        if (!is.null(opt$positive_class) && opt$positive_class %in% group_names_ordered) {
                          pos_class <- opt$positive_class
                          neg_class <- setdiff(group_names_ordered, opt$positive_class)
                        } else {
                          # Default: second group in order is positive
                          pos_class <- group_names_ordered[2]
                          neg_class <- group_names_ordered[1]
                        }

                        log_msg(paste("[AUC] Binary - Positive class:", pos_class, "| Negative class:", neg_class))

                        # Verify positive class column exists
                        if (pos_class %in% colnames(result$pred_prob)) {
                          roc_obj <- pROC::roc(
                            response = result$actual_labels,
                            predictor = result$pred_prob[, pos_class], # Use positive class column
                            levels = c(neg_class, pos_class), # CRITICAL: negative first, positive second
                            direction = "<",
                            quiet = TRUE
                          )
                          auc <- as.numeric(pROC::auc(roc_obj))
                          log_msg(paste("[OK] Binary AUC:", round(auc, 3)))
                        } else {
                          log_warning(paste("Positive class column", pos_class, "not found in predictions"))
                          auc <- 0
                        }
                      } else {
                        # ===== MULTI-CLASS CLASSIFICATION =====
                        # Calculate one-vs-rest AUC for each class IN ORDER
                        log_msg(paste("[AUC] Multi-class: calculating one-vs-rest AUC"))

                        auc_values <- numeric(n_classes)
                        names(auc_values) <- group_names_ordered

                        for (i in seq_along(group_names_ordered)) {
                          class_name <- group_names_ordered[i]
                          binary_labels <- ifelse(result$actual_labels == class_name, 1, 0)

                          if (class_name %in% colnames(result$pred_prob)) {
                            roc_obj <- pROC::roc(
                              response = binary_labels,
                              predictor = result$pred_prob[, class_name],
                              quiet = TRUE
                            )
                            auc_values[i] <- as.numeric(pROC::auc(roc_obj))
                          } else {
                            auc_values[i] <- 0
                          }
                        }

                        auc <- mean(auc_values, na.rm = TRUE)
                        log_msg(paste("[OK] Multi-class mean AUC:", round(auc, 3)))
                        log_msg(paste(
                          "     Per-class AUCs:",
                          paste(group_names_ordered, "=", round(auc_values, 3), collapse = ", ")
                        ))
                      }

                      # Sanity check
                      if (is.na(auc) || !is.finite(auc)) {
                        log_warning("AUC calculation returned invalid value, setting to 0")
                        auc <- 0
                      }
                    },
                    error = function(e) {
                      log_warning(paste("AUC calculation failed:", e$message), "AUC calculation")
                      auc <- 0
                    }
                  )

                  # ===== PR-AUC CALCULATION (FRIDAY EDITION!) =====
                  pr_auc <- 0
                  tryCatch(
                    {
                      n_classes <- length(levels(result$actual_labels))

                      if (n_classes == 2) {
                        # Binary classification PR-AUC
                        if (!is.null(opt$positive_class) && opt$positive_class %in% group_names_ordered) {
                          pos_class <- opt$positive_class
                          neg_class <- setdiff(group_names_ordered, opt$positive_class)
                        } else {
                          pos_class <- group_names_ordered[2]
                          neg_class <- group_names_ordered[1]
                        }

                        if (pos_class %in% colnames(result$pred_prob)) {
                          # Binary labels: 1 = positive, 0 = negative
                          binary_labels <- ifelse(result$actual_labels == pos_class, 1, 0)
                          pr_curve <- PRROC::pr.curve(
                            scores.class0 = result$pred_prob[, pos_class],
                            weights.class0 = binary_labels,
                            curve = FALSE
                          )
                          pr_auc <- pr_curve$auc.integral
                          log_msg(paste("[OK] Binary PR-AUC:", round(pr_auc, 3)))
                        }
                      } else {
                        # Multi-class: calculate one-vs-rest PR-AUC for each class
                        pr_auc_values <- numeric(n_classes)
                        names(pr_auc_values) <- group_names_ordered

                        for (i in seq_along(group_names_ordered)) {
                          class_name <- group_names_ordered[i]
                          binary_labels <- ifelse(result$actual_labels == class_name, 1, 0)

                          if (class_name %in% colnames(result$pred_prob)) {
                            pr_curve <- PRROC::pr.curve(
                              scores.class0 = result$pred_prob[, class_name],
                              weights.class0 = binary_labels,
                              curve = FALSE
                            )
                            pr_auc_values[i] <- pr_curve$auc.integral
                          } else {
                            pr_auc_values[i] <- 0
                          }
                        }

                        pr_auc <- mean(pr_auc_values, na.rm = TRUE)
                        log_msg(paste("[OK] Multi-class mean PR-AUC:", round(pr_auc, 3)))
                        log_msg(paste(
                          "     Per-class PR-AUCs:",
                          paste(group_names_ordered, "=", round(pr_auc_values, 3), collapse = ", ")
                        ))
                      }

                      if (is.na(pr_auc) || !is.finite(pr_auc)) {
                        log_warning("PR-AUC calculation returned invalid value, setting to 0")
                        pr_auc <- 0
                      }
                    },
                    error = function(e) {
                      log_warning(paste("PR-AUC calculation failed:", e$message), "PR-AUC calculation")
                      pr_auc <- 0
                    }
                  )

                  # Update detailed_dmr_info with calculated AUC and PR-AUC
                  detailed_dmr_info$auc <- auc
                  detailed_dmr_info$pr_auc <- pr_auc
                },
                error = function(e) {
                  log_msg(paste("Warning: Error in metrics calculation:", e$message))
                  log_msg("Continuing with basic results...")

                  # Set safe defaults
                  sensitivity <- c(0)
                  specificity <- c(0)
                  f1_score <- c(0)
                  auc <- 0
                  pr_auc <- 0
                }
              )

              dmr_details_file <- file.path(model_dir, "dmr_details_with_coordinates.csv")
              write.csv(detailed_dmr_info, dmr_details_file, row.names = FALSE)

              model_rds_file <- file.path(model_dir, "model.rds")
              saveRDS(result$model_fit, model_rds_file)

              # Generate automatic model summary
              summary_file <- generate_model_summary(model_rds_file, model_dir)
              # ... continue with the rest

              # Store results summary
              # Ensure all vectors have the same length

              # DEBUG SECTION
              if (length(sensitivity) > 0) {
                log_msg(paste("DEBUG: Sensitivity values:", paste(round(sensitivity, 3), collapse = ", ")))
              }
              if (length(specificity) > 0) {
                log_msg(paste("DEBUG: Specificity values:", paste(round(specificity, 3), collapse = ", ")))
              }

              # ===== PERFORM CROSS-VALIDATION - BOTH REGULAR AND NESTED =====
              cv_result_regular <- NULL
              cv_result_nested <- NULL

              if (opt$cv_repeats > 0) {
                # STEP 1: ALWAYS run regular CV (fast validation)
                log_msg(paste("[CV] Running REGULAR cross-validation:", opt$cv_repeats, "repeats for", model_name))

                cv_result_regular <- perform_cross_validation(
                  ml_data_subset = ml_data_subset,
                  model_name = model_name,
                  model_specs = model_specs,
                  ctrl = ctrl,
                  metric = metric,
                  n_repeats = opt$cv_repeats,
                  output_dir = model_dir,
                  reg_label = reg_label,
                  n_markers = n_markers,
                  nested = FALSE # Regular CV
                )

                log_msg(paste("[CV] Regular CV complete for", model_name))
                log_msg(paste(
                  "    Mean Accuracy:", round(cv_result_regular$summary$mean_accuracy, 3),
                  "95% CI [", round(cv_result_regular$summary$accuracy_95ci_lower, 3),
                  "-", round(cv_result_regular$summary$accuracy_95ci_upper, 3), "]"
                ))

                # STEP 2: Optionally run nested CV (gold standard, slow)
                if (opt$nested_cv) {
                  log_msg("[CV] Running NESTED cross-validation (unbiased, gold standard)...")
                  log_msg("    Note: This will take 5-10x longer...")

                  cv_result_nested <- perform_cross_validation(
                    ml_data_subset = ml_data_subset,
                    model_name = model_name,
                    model_specs = model_specs,
                    ctrl = ctrl,
                    metric = metric,
                    n_repeats = opt$cv_repeats,
                    output_dir = model_dir,
                    reg_label = reg_label,
                    n_markers = n_markers,
                    nested = TRUE # Nested CV
                  )

                  log_msg(paste("[CV] Nested CV complete for", model_name))
                  log_msg(paste(
                    "    Mean Accuracy (UNBIASED):", round(cv_result_nested$summary$mean_accuracy, 3),
                    "95% CI [", round(cv_result_nested$summary$accuracy_95ci_lower, 3),
                    "-", round(cv_result_nested$summary$accuracy_95ci_upper, 3), "]"
                  ))

                  # Calculate and report optimism bias
                  bias_estimate <- cv_result_regular$summary$mean_accuracy - cv_result_nested$summary$mean_accuracy
                  log_msg(paste("[ANALYSIS] Optimism bias estimate:", round(bias_estimate, 3)))

                  if (abs(bias_estimate) > 0.05) {
                    log_msg("    [WARN] Regular CV shows >5% optimism - use nested CV for publication")
                  } else {
                    log_msg("    [OK] Minimal optimism detected - both CV methods agree")
                  }
                }
              }
              # ===== END CV BLOCK =====

              # ===== CREATE NESTED CV VISUALIZATIONS =====
              if (!is.null(cv_result_nested) && opt$create_plots) {
                log_msg("[PLOT] Creating nested CV visualization suite...")

                nested_cv_plots <- plot_nested_cv_results(
                  cv_result_regular = cv_result_regular,
                  cv_result_nested = cv_result_nested,
                  model_name = model_name,
                  output_dir = model_dir,
                  region_label = reg_label,
                  n_markers = n_markers
                )

                if (!is.null(nested_cv_plots)) {
                  log_msg(paste("[OK] Created", length(nested_cv_plots), "nested CV plots"))
                  for (plot_name in names(nested_cv_plots)) {
                    log_msg(paste("    -", plot_name, ":", basename(nested_cv_plots[[plot_name]])))
                  }
                }
              }
              # ===== END NESTED CV VISUALIZATIONS =====

              # NOW CREATE results_row WITH BOTH REGULAR AND NESTED CV METRICS
              results_row <- data.frame(
                region = reg_label,
                model = model_name,
                markers = n_markers,
                accuracy = result$accuracy,
                auc = ifelse(is.na(auc), 0, auc),
                pr_auc = ifelse(exists("pr_auc") && !is.na(pr_auc), pr_auc, 0),
                mean_sensitivity = ifelse(length(sensitivity) > 0, mean(sensitivity, na.rm = TRUE), 0),
                mean_specificity = ifelse(length(specificity) > 0, mean(specificity, na.rm = TRUE), 0),
                mean_f1 = ifelse(length(f1_score) > 0, mean(f1_score, na.rm = TRUE), 0),
                train_time_sec = result$model_time,
                test_ratio = opt$test_ratio,
                train_samples = nrow(train_data),
                test_samples = nrow(test_data),

                # REGULAR CV METRICS (for development, slightly optimistic)
                cv_mean_accuracy = if (!is.null(cv_result_regular)) cv_result_regular$summary$mean_accuracy else NA,
                cv_mean_sensitivity = if (!is.null(cv_result_regular)) cv_result_regular$summary$mean_sensitivity else NA,
                cv_mean_specificity = if (!is.null(cv_result_regular)) cv_result_regular$summary$mean_specificity else NA,
                cv_mean_auc = if (!is.null(cv_result_regular)) cv_result_regular$summary$mean_auc else NA,
                cv_95ci_lower = if (!is.null(cv_result_regular)) cv_result_regular$summary$accuracy_95ci_lower else NA,
                cv_95ci_upper = if (!is.null(cv_result_regular)) cv_result_regular$summary$accuracy_95ci_upper else NA,
                cv_sd_accuracy = if (!is.null(cv_result_regular)) cv_result_regular$summary$sd_accuracy else NA,

                # NESTED CV METRICS (for publication, unbiased gold standard)
                nested_cv_mean_accuracy = if (!is.null(cv_result_nested)) cv_result_nested$summary$mean_accuracy else NA,
                nested_cv_mean_sensitivity = if (!is.null(cv_result_nested)) cv_result_nested$summary$mean_sensitivity else NA,
                nested_cv_mean_specificity = if (!is.null(cv_result_nested)) cv_result_nested$summary$mean_specificity else NA,
                nested_cv_mean_auc = if (!is.null(cv_result_nested)) cv_result_nested$summary$mean_auc else NA,
                nested_cv_95ci_lower = if (!is.null(cv_result_nested)) cv_result_nested$summary$accuracy_95ci_lower else NA,
                nested_cv_95ci_upper = if (!is.null(cv_result_nested)) cv_result_nested$summary$accuracy_95ci_upper else NA,
                nested_cv_sd_accuracy = if (!is.null(cv_result_nested)) cv_result_nested$summary$sd_accuracy else NA,

                # OPTIMISM BIAS ESTIMATE (difference between regular and nested)
                cv_optimism_bias = if (!is.null(cv_result_regular) && !is.null(cv_result_nested)) {
                  cv_result_regular$summary$mean_accuracy - cv_result_nested$summary$mean_accuracy
                } else {
                  NA
                },
                stringsAsFactors = FALSE
              )

              results_summary <- rbind(results_summary, results_row)

              # Save predictions and metrics to CSV files
              tryCatch(
                {
                  # Save predictions
                  predictions_df <- data.frame(
                    sample_id = rownames(test_data),
                    actual = result$actual_labels,
                    predicted = result$pred_class,
                    result$pred_prob,
                    stringsAsFactors = FALSE
                  )
                  predictions_file <- file.path(model_dir, "predictions.csv")
                  write.csv(predictions_df, predictions_file, row.names = FALSE)
                  log_file_operation("Saved predictions", predictions_file)

                  # Save metrics
                  metrics_df <- data.frame(
                    metric = c("Accuracy", "AUC", "PR_AUC", "Mean_Sensitivity", "Mean_Specificity", "Mean_F1"),
                    value = c(
                      result$accuracy,
                      ifelse(exists("auc") && !is.na(auc), auc, 0),
                      ifelse(exists("pr_auc") && !is.na(pr_auc), pr_auc, 0),
                      ifelse(exists("sensitivity") && length(sensitivity) > 0, mean(sensitivity, na.rm = TRUE), 0),
                      ifelse(exists("specificity") && length(specificity) > 0, mean(specificity, na.rm = TRUE), 0),
                      ifelse(exists("f1_score") && length(f1_score) > 0, mean(f1_score, na.rm = TRUE), 0)
                    ),
                    stringsAsFactors = FALSE
                  )

                  # Add per-class metrics if available
                  if (exists("sensitivity") && length(sensitivity) > 0) {
                    for (i in seq_along(sensitivity)) {
                      class_name <- names(sensitivity)[i]
                      metrics_df <- rbind(metrics_df, data.frame(
                        metric = c(
                          paste0("Sensitivity_", class_name),
                          paste0("Specificity_", class_name),
                          paste0("Precision_", class_name),
                          paste0("F1_", class_name)
                        ),
                        value = c(sensitivity[i], specificity[i], precision[i], f1_score[i]),
                        stringsAsFactors = FALSE
                      ))
                    }
                  }

                  metrics_file <- file.path(model_dir, "metrics.csv")
                  write.csv(metrics_df, metrics_file, row.names = FALSE)
                  log_file_operation("Saved metrics", metrics_file)
                },
                error = function(e) {
                  log_warning(paste("Failed to save predictions/metrics:", e$message))
                }
              )

              # NOW CREATE PLOTS (only for saved models)

              # Add this debug code before the plotting section:

              # ===== PLOTTING SECTION - CLEAN VERSION =====
              # DEBUG: Verify group order is preserved
              cat("\n=== PLOTTING DIAGNOSTIC ===\n")
              cat("treatment_processed exists:", exists("treatment_processed"), "\n")
              if (exists("treatment_processed")) {
                cat("Group names in treatment_processed:", paste(treatment_processed$group_names, collapse = ", "), "\n")
                cat("Group colors in treatment_processed:", paste(names(treatment_processed$group_colors), "=",
                  treatment_processed$group_colors,
                  collapse = ", "
                ), "\n")
              }
              cat("current_treatment levels:", paste(levels(current_treatment), collapse = ", "), "\n")
              cat("opt$group_names:", opt$group_names, "\n")
              cat("opt$group_colors:", opt$group_colors, "\n")
              cat("========================\n\n")


              if (opt$create_plots) {
                verbose_msg(paste("Creating plots for", n_markers, "markers in", model_name, "model"))

                # STEP 1: Prepare data for plotting
                current_methylation_matrix <- result$current_methylation_matrix
                n_samples_in_matrix <- nrow(current_methylation_matrix)

                # STEP 2: Get proper treatment labels WITH PRESERVED ORDER
                if (exists("treatment_processed") && !is.null(treatment_processed)) {
                  # FIXED: Match by rownames instead of positional indexing
                  matrix_sample_names <- rownames(current_methylation_matrix)

                  # Get full treatment labels with sample IDs as names
                  full_treatment_with_names <- treatment_processed$labels
                  # Use getSampleID from methylKit object to get sample names
                  if (exists("meth") && !is.null(meth)) {
                    meth_sample_ids <- getSampleID(meth)
                    names(full_treatment_with_names) <- meth_sample_ids

                    # Subset to match current matrix samples by name
                    current_treatment <- full_treatment_with_names[matrix_sample_names]

                    # Convert to factor with proper levels
                    current_treatment <- factor(current_treatment, levels = treatment_processed$group_names)

                    log_msg(paste("[DEBUG] Sample name matching:", sum(!is.na(current_treatment)), "of", length(current_treatment), "matched"))
                    log_msg(paste("[DEBUG] Matrix sample names (first 3):", paste(head(matrix_sample_names, 3), collapse = ", ")))
                    log_msg(paste("[DEBUG] Meth sample IDs (first 3):", paste(head(meth_sample_ids, 3), collapse = ", ")))
                  } else {
                    # Fallback to positional matching if meth object not available
                    log_warning("meth object not available - using positional matching", "Plotting")
                    current_treatment <- factor(
                      treatment_processed$labels[1:n_samples_in_matrix],
                      levels = treatment_processed$group_names
                    )
                  }
                  group_names_for_plots <- treatment_processed$group_names
                  group_colors_for_plots <- treatment_processed$group_colors

                  log_msg(paste("[PLOT] Using group order:", paste(group_names_for_plots, collapse = " -> ")))
                  if (!is.null(group_colors_for_plots)) {
                    log_msg(paste("[PLOT] Using custom colors:", paste(names(group_colors_for_plots), "=", group_colors_for_plots, collapse = ", ")))
                  }
                } else if (!is.null(opt$group_names)) {
                  # Use CLI-specified order
                  group_names_for_plots <- trimws(strsplit(opt$group_names, ",")[[1]])
                  current_treatment <- factor(treatment_info[1:n_samples_in_matrix], levels = group_names_for_plots)

                  if (!is.null(opt$group_colors)) {
                    group_colors_for_plots <- trimws(strsplit(opt$group_colors, ",")[[1]])
                    names(group_colors_for_plots) <- group_names_for_plots
                  } else {
                    group_colors_for_plots <- NULL
                  }
                } else {
                  # Fallback
                  current_treatment <- factor(treatment_info[1:n_samples_in_matrix])
                  group_names_for_plots <- levels(current_treatment)
                  group_colors_for_plots <- NULL
                }

                # STEP 3: Store these in a plotting context object for consistency
                plotting_context <- list(
                  treatment_groups = current_treatment,
                  group_names = group_names_for_plots,
                  group_colors = group_colors_for_plots
                )

                # STEP 3: Validate dimensions
                if (length(current_treatment) != n_samples_in_matrix) {
                  log_warning("Treatment length mismatch - using generic labels", "Plotting")
                  current_treatment <- factor(rep(paste0("Group_", 1:3), length.out = n_samples_in_matrix))
                }

                verbose_msg(paste("Plot treatment labels:", length(current_treatment), "samples"))
                log_msg(paste(
                  "Plotting:", nrow(current_methylation_matrix), "samples x",
                  ncol(current_methylation_matrix), "markers"
                ))

                # STEP 4: Handle missing values in methylation matrix
                for (col in 1:ncol(current_methylation_matrix)) {
                  if (any(is.na(current_methylation_matrix[, col]))) {
                    median_val <- median(current_methylation_matrix[, col], na.rm = TRUE)
                    if (is.na(median_val)) median_val <- 50
                    current_methylation_matrix[is.na(current_methylation_matrix[, col]), col] <- median_val
                  }
                }

                # STEP 5: Create DMR bar plot
                tryCatch(
                  {
                    log_msg("[PLOT] Creating DMR bar plot...")
                    bar_plot_file <- plot_dmr_barplot(
                      dmrs = result$current_dmrs,
                      treatment_groups = plotting_context$treatment_groups,
                      output_dir = model_dir,
                      region_label = paste(reg_label, model_name),
                      n_markers = n_markers
                    )

                    if (!is.null(bar_plot_file)) {
                      verbose_msg(paste("Bar plot saved:", basename(bar_plot_file)))
                    } else {
                      log_warning("Bar plot returned NULL", "Plotting")
                    }
                  },
                  error = function(e) {
                    log_error(paste("Bar plot failed:", e$message), "Plotting")
                  }
                )

                # STEP 6: Create group methylation summary WITH EXPLICIT PARAMS
                tryCatch(
                  {
                    log_msg("[PLOT] Creating group methylation summary...")
                    group_summary_files <- plot_group_methylation_summary(
                      methylation_matrix = current_methylation_matrix,
                      treatment_groups = plotting_context$treatment_groups,
                      dmrs = result$current_dmrs,
                      output_dir = model_dir,
                      region_label = paste(reg_label, model_name),
                      n_markers = n_markers,
                      group_names = plotting_context$group_names, # EXPLICIT
                      group_colors = plotting_context$group_colors # EXPLICIT
                    )

                    if (!is.null(group_summary_files) && "plot" %in% names(group_summary_files)) {
                      verbose_msg(paste("Group summary saved:", basename(group_summary_files$plot)))
                    } else {
                      log_warning("Group summary returned NULL", "Plotting")
                    }
                  },
                  error = function(e) {
                    log_error(paste("Group summary failed:", e$message), "Plotting")
                  }
                )


                # Create corrected DMRs for heatmap (with proper direction)
                dmrs_for_heatmap <- result$current_dmrs

                # Add corrected direction if we have detailed_dmr_info
                if (exists("detailed_dmr_info") && "direction" %in% colnames(detailed_dmr_info)) {
                  # Match DMR_IDs and update direction
                  for (i in 1:nrow(dmrs_for_heatmap)) {
                    dmr_id <- dmrs_for_heatmap$DMR_ID[i]
                    matching_row <- detailed_dmr_info[detailed_dmr_info$DMR_ID == dmr_id, ]
                    if (nrow(matching_row) > 0) {
                      dmrs_for_heatmap$direction[i] <- matching_row$direction[1]
                    }
                  }
                  log_msg(paste("[PLOT] Updated", nrow(dmrs_for_heatmap), "DMRs with corrected direction"))
                } else {
                  log_warning("detailed_dmr_info not available - heatmap will use raw direction", "Plotting")
                }

                # NOW call heatmap with corrected DMRs
                tryCatch(
                  {
                    log_msg("[PLOT] Creating methylation heatmap...")
                    heatmap_files <- plot_methylation_heatmap(
                      methylation_matrix = current_methylation_matrix,
                      treatment_groups = current_treatment,
                      dmrs = dmrs_for_heatmap, # <- USE CORRECTED VERSION
                      output_dir = model_dir,
                      region_label = paste(reg_label, model_name),
                      n_markers = n_markers
                    )

                    if (!is.null(heatmap_files)) {
                      if ("heatmap" %in% names(heatmap_files)) {
                        verbose_msg(paste("Heatmap saved:", basename(heatmap_files$heatmap)))
                      }
                      if ("kmeans_heatmap" %in% names(heatmap_files)) {
                        verbose_msg(paste("K-means heatmap saved:", basename(heatmap_files$kmeans_heatmap)))
                      }
                    } else {
                      log_warning("Heatmap returned NULL", "Plotting")
                    }
                  },
                  error = function(e) {
                    log_error(paste("Heatmap failed:", e$message), "Plotting")
                  }
                )

                # STEP 8: Create ROC curves
                tryCatch(
                  {
                    log_msg("[PLOT] Creating ROC curves...")
                    roc_files <- plot_roc_curves(
                      pred_prob = result$pred_prob,
                      actual_labels = result$actual_labels,
                      model_name = model_name,
                      output_dir = model_dir,
                      region_label = paste(reg_label, model_name),
                      n_markers = n_markers,
                      positive_class = opt$positive_class
                    )

                    if (!is.null(roc_files)) {
                      if (!is.null(roc_files$binary_roc)) {
                        verbose_msg(paste("Binary ROC saved:", basename(roc_files$binary_roc)))
                      }
                      if (!is.null(roc_files$multiclass_roc)) {
                        verbose_msg(paste("Multi-class ROC saved:", basename(roc_files$multiclass_roc)))
                      }
                    }
                  },
                  error = function(e) {
                    log_error(paste("ROC curves failed:", e$message), "Plotting")
                  }
                )

                # STEP 8.5: Create PR curves
                tryCatch(
                  {
                    log_msg("[PLOT] Creating PR curves...")
                    pr_files <- plot_pr_curves(
                      pred_prob = result$pred_prob,
                      actual_labels = result$actual_labels,
                      model_name = model_name,
                      output_dir = model_dir,
                      region_label = paste(reg_label, model_name),
                      n_markers = n_markers,
                      positive_class = opt$positive_class
                    )

                    if (!is.null(pr_files)) {
                      if (!is.null(pr_files$binary_pr)) {
                        verbose_msg(paste("Binary PR curve saved:", basename(pr_files$binary_pr)))
                      }
                      if (!is.null(pr_files$multiclass_pr)) {
                        verbose_msg(paste("Multi-class PR curves saved:", basename(pr_files$multiclass_pr)))
                      }
                    }
                  },
                  error = function(e) {
                    log_error(paste("PR curves failed:", e$message), "Plotting")
                  }
                )

                # STEP 9: Create prediction probability plots
                tryCatch(
                  {
                    log_msg("[PLOT] Creating prediction probability plots...")
                    prob_files <- plot_prediction_probabilities(
                      result$pred_prob,
                      result$actual_labels,
                      model_name,
                      model_dir,
                      paste(reg_label, model_name),
                      n_markers
                    )
                  },
                  error = function(e) {
                    log_error(paste("Probability plots failed:", e$message), "Plotting")
                  }
                )

                # STEP 10: Random Forest specific visualizations
                if (model_name == "rf" || model_name == "ranger") {
                  tryCatch(
                    {
                      log_msg("[PLOT] Creating Random Forest visualizations...")
                      rf_viz_files <- plot_random_forest_viz(
                        result$model_fit,
                        train_data,
                        result$current_dmrs,
                        model_dir,
                        paste(reg_label, model_name),
                        n_markers
                      )
                      if (!is.null(rf_viz_files)) {
                        verbose_msg(paste("Random Forest visualizations created:", length(rf_viz_files), "plots"))
                      }
                    },
                    error = function(e) {
                      log_error(paste("RF visualizations failed:", e$message), "Plotting")
                    }
                  )
                }

                # STEP 11: Confusion matrix plot
                tryCatch(
                  {
                    log_msg("[PLOT] Creating confusion matrix...")

                    # ===== GET GROUP ORDER AND COLORS =====
                    group_names_ordered <- NULL
                    point_colors <- NULL

                    # Priority 1: Check global environment for treatment_processed
                    if (exists("treatment_processed", envir = .GlobalEnv)) {
                      tp <- get("treatment_processed", envir = .GlobalEnv)
                      if (!is.null(tp) && !is.null(tp$group_names)) {
                        group_names_ordered <- tp$group_names
                        point_colors <- tp$group_colors
                      }
                    }

                    # Priority 2: Check global environment for opt
                    if (is.null(group_names_ordered) && exists("opt", envir = .GlobalEnv)) {
                      opt_global <- get("opt", envir = .GlobalEnv)
                      if (!is.null(opt_global$group_names)) {
                        group_names_ordered <- trimws(strsplit(opt_global$group_names, ",")[[1]])
                        if (!is.null(opt_global$group_colors)) {
                          point_colors <- trimws(strsplit(opt_global$group_colors, ",")[[1]])
                          names(point_colors) <- group_names_ordered
                        }
                      }
                    }

                    # Priority 3: Fallback to factor levels
                    if (is.null(group_names_ordered)) {
                      group_names_ordered <- levels(result$actual_labels)
                    }

                    # Set default colors if needed
                    if (is.null(point_colors)) {
                      n_groups <- length(group_names_ordered)
                      point_colors <- brewer.pal(max(3, n_groups), "Set2")[1:n_groups]
                      names(point_colors) <- group_names_ordered
                    }

                    log_msg(paste("[PLOT] Confusion matrix order:", paste(group_names_ordered, collapse = " -> ")))
                    log_msg(paste("[PLOT] Confusion matrix colors:", paste(names(point_colors), "=", point_colors, collapse = ", ")))

                    # CRITICAL FIX: Convert to character FIRST, then ordered factor
                    pred_results <- data.frame(
                      actual = as.character(result$actual_labels),
                      predicted = as.character(result$pred_class),
                      stringsAsFactors = FALSE
                    )

                    # NOW convert to ordered factors - CRITICAL ORDER
                    pred_results$actual <- factor(pred_results$actual, levels = group_names_ordered)
                    pred_results$predicted <- factor(pred_results$predicted, levels = group_names_ordered)

                    # DEBUG: Verify the order
                    log_msg(paste("[DEBUG] Confusion matrix actual levels:", paste(levels(pred_results$actual), collapse = ", ")))
                    log_msg(paste("[DEBUG] Confusion matrix predicted levels:", paste(levels(pred_results$predicted), collapse = ", ")))
                    log_msg(paste("[DEBUG] Sample data - actual:", paste(head(pred_results$actual), collapse = ", ")))
                    log_msg(paste("[DEBUG] Sample data - predicted:", paste(head(pred_results$predicted), collapse = ", ")))

                    # DEBUG: Verify factor levels
                    log_msg(paste("[DEBUG] Actual levels:", paste(levels(pred_results$actual), collapse = ", ")))
                    log_msg(paste("[DEBUG] Predicted levels:", paste(levels(pred_results$predicted), collapse = ", ")))

                    conf_matrix_plot <- ggplot(pred_results, aes(x = actual, y = predicted, color = actual)) +
                      geom_jitter(alpha = 0.6, width = 0.2, height = 0.2, size = 3) +
                      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
                      scale_color_manual(
                        values = point_colors,
                        limits = group_names_ordered,
                        breaks = group_names_ordered, # ADD THIS
                        drop = FALSE
                      ) +
                      scale_x_discrete(
                        limits = group_names_ordered,
                        breaks = group_names_ordered, # ADD THIS
                        drop = FALSE
                      ) +
                      scale_y_discrete(
                        limits = group_names_ordered,
                        breaks = group_names_ordered, # ADD THIS
                        drop = FALSE
                      ) +
                      labs(
                        title = paste("Predictions:", model_name, "with", n_markers, "markers"),
                        subtitle = paste("Accuracy:", round(result$accuracy, 3), "| Region:", reg_label),
                        x = "Actual Group",
                        y = "Predicted Group",
                        color = "Actual Group"
                      ) +
                      theme_minimal() +
                      theme(
                        plot.title = element_text(size = 14, face = "bold"),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
                        axis.text.y = element_text(size = 11),
                        legend.position = "right",
                        legend.title = element_text(size = 12, face = "bold")
                      )

                    conf_plot_file <- file.path(model_dir, "confusion_matrix_plot.png")
                    save_plot_universal(conf_matrix_plot, conf_plot_file, width = 10, height = 8)
                    log_msg(paste("[PLOT] Confusion matrix saved:", basename(conf_plot_file)))
                  },
                  error = function(e) {
                    log_error(paste("Confusion matrix failed:", e$message), "Plotting")
                    log_error(paste("Error details:", e$message), "Plotting")
                  }
                )

                # STEP 12: Metrics bar plot
                tryCatch(
                  {
                    log_msg("[PLOT] Creating metrics bar plot...")

                    metrics_plot_df <- data.frame(
                      metric = c("Accuracy", "AUC", "Mean_Sensitivity", "Mean_Specificity", "Mean_F1"),
                      value = c(
                        result$accuracy,
                        ifelse(exists("auc") && !is.na(auc), auc, 0),
                        ifelse(exists("sensitivity") && length(sensitivity) > 0, mean(sensitivity, na.rm = TRUE), 0),
                        ifelse(exists("specificity") && length(specificity) > 0, mean(specificity, na.rm = TRUE), 0),
                        ifelse(exists("f1_score") && length(f1_score) > 0, mean(f1_score, na.rm = TRUE), 0)
                      ),
                      stringsAsFactors = FALSE
                    )

                    metrics_plot <- ggplot(metrics_plot_df, aes(x = metric, y = value, fill = metric)) +
                      geom_bar(stat = "identity", alpha = 0.8) +
                      geom_text(aes(label = round(value, 3)), vjust = -0.3) +
                      ylim(0, 1.1) +
                      labs(
                        title = paste("Model Performance:", model_name),
                        subtitle = paste(n_markers, "markers | Region:", reg_label),
                        x = "Metric",
                        y = "Score"
                      ) +
                      theme_minimal() +
                      theme(legend.position = "none")

                    metrics_plot_file <- file.path(model_dir, "performance_metrics")
                    save_plot_universal(metrics_plot, metrics_plot_file, width = 8, height = 6)
                    verbose_msg(paste("Metrics plot saved:", basename(metrics_plot_file)))
                  },
                  error = function(e) {
                    log_error(paste("Metrics plot failed:", e$message), "Plotting")
                  }
                )

                log_msg(paste("[PLOT] Completed all plots for", model_name, "with", n_markers, "markers"))
              } # END if (opt$create_plots)
              # ===== END PLOTTING SECTION =====
            } else {
              verbose_msg(paste(
                "[SKIP] Skipping", model_name, "model (accuracy:", round(result$accuracy, 3),
                "< threshold:", opt$min_accuracy, ")"
              ))
            }
          }
        } else {
          log_msg(paste(
            "[ERROR] No models achieved accuracy >=", opt$min_accuracy,
            "(best:", round(max_accuracy, 3), ") - skipping plots and file saving for", n_markers, "markers"
          ))
        }

        cat("\n") # Add space between marker sets
      }

      # Complete marker progress bar
      if (marker_pb$current < marker_pb$total) {
        marker_pb <- update_progress_bar(marker_pb, marker_pb$total, "All markers complete")
      }
      cat("\n")

      region_end <- Sys.time()
      region_time <- round(as.numeric(region_end - region_start_time, units = "mins"), 1)
      log_msg(paste("[OK] Region", reg_label, "completed in", region_time, "minutes"))
    },
    error = function(e) {
      log_msg(paste("[ERROR] Error processing", reg_label, ":", e$message))
      log_msg(paste("Error details:", toString(e)))
      log_msg(paste("Call stack:", paste(deparse(e$call), collapse = " ")))
      traceback()
    }
  )
}

# Final summary
total_end_time <- Sys.time()
total_time <- round(as.numeric(total_end_time - total_start_time, units = "mins"), 1)

log_msg("[SUCCESS] DMR + ML ANALYSIS v2.3.2-DUAL-PROGRESS COMPLETE!")
log_msg(paste("[DATA] Total processing time:", total_time, "minutes"))

if (nrow(results_summary) > 0) {
  log_msg(paste("[DATA] Total models trained:", nrow(results_summary)))

  summary_file <- file.path(opt$output_dir, "analysis_summary.csv")
  write.csv(results_summary, summary_file, row.names = FALSE)
  log_msg(paste("[DATA] Results summary saved to:", basename(summary_file)))

  best_results <- results_summary[order(-results_summary$accuracy), ]
  log_msg("[DATA] Best results:")

  for (i in 1:min(5, nrow(best_results))) {
    result <- best_results[i, ]
    log_msg(sprintf(
      "   %s: %s | %d markers | Acc=%.3f | AUC=%.3f | Split=%d/%d",
      result$region, result$model, result$markers, result$accuracy, result$auc,
      result$train_samples, result$test_samples
    ))
  }
} else {
  log_msg("[WARN] No models were successfully trained")
}

# ===== CV COMPARISON SUMMARY =====
if (nrow(results_summary) > 0) {
  has_nested <- any(!is.na(results_summary$nested_cv_mean_accuracy))

  if (has_nested) {
    cat("\n")
    log_msg("================================================================================")
    log_msg("CROSS-VALIDATION COMPARISON SUMMARY")
    log_msg("================================================================================")

    cv_comparison <- results_summary[
      !is.na(results_summary$nested_cv_mean_accuracy),
      c(
        "region", "model", "markers",
        "cv_mean_accuracy", "nested_cv_mean_accuracy",
        "cv_optimism_bias"
      )
    ]

    if (nrow(cv_comparison) > 0) {
      for (i in 1:nrow(cv_comparison)) {
        log_msg(sprintf(
          "\n%s - %s (%d markers):",
          cv_comparison$region[i],
          cv_comparison$model[i],
          cv_comparison$markers[i]
        ))
        log_msg(sprintf(
          "  Regular CV:  %.3f (development/comparison)",
          cv_comparison$cv_mean_accuracy[i]
        ))
        log_msg(sprintf(
          "  Nested CV:   %.3f (unbiased/publication)",
          cv_comparison$nested_cv_mean_accuracy[i]
        ))
        log_msg(sprintf(
          "  Bias:        %+.3f",
          cv_comparison$cv_optimism_bias[i]
        ))

        if (abs(cv_comparison$cv_optimism_bias[i]) > 0.05) {
          log_msg("  [WARN] Optimism >5% detected - report nested CV in publications")
        } else if (abs(cv_comparison$cv_optimism_bias[i]) < 0.02) {
          log_msg("  [OK] Minimal bias - both methods agree well")
        } else {
          log_msg("  [INFO] Moderate bias - prefer nested CV for formal reporting")
        }
      }

      log_msg("\n================================================================================")
      log_msg("RECOMMENDATION:")
      log_msg("  * Use REGULAR CV for: Model comparison, feature selection, development")
      log_msg("  * Use NESTED CV for:  Publication, clinical validation, formal claims")
      log_msg("  * Report both in papers to demonstrate robustness")
      log_msg("================================================================================\n")
    }
  }
}
# ===== END CV COMPARISON =====

# ================================================================================
# GENERATE FINAL ANALYSIS REPORT
# ================================================================================
if (nrow(results_summary) > 0) {
  summary_file <- file.path(opt$output_dir, "analysis_summary.csv")
  if (file.exists(summary_file)) {
    tryCatch(
      {
        generate_final_analysis_report(summary_file, opt$output_dir)
      },
      error = function(e) {
        log_msg(paste("[WARN] Report generation failed:", e$message))
      }
    )
  }
}

log_msg(paste("[DIR] All outputs saved to:", unique_output_dir))

# ================================================================================
# CLEANUP
# ================================================================================
if (!opt$reproducible && exists("cl")) {
  stopCluster(cl)
  log_msg("[CLEANUP] Parallel cluster stopped")
}

cat("[COMPLETE] Analysis complete!\n")
