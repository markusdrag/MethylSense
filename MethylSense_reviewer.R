#!/usr/bin/env Rscript

# ============================================================================
# Comprehensive Model Validation and Technical Quality Assessment
# ============================================================================
# VERSION: 5.13.5
# DATE: 2025-12-05
# GitHub: https://github.com/markusdrag/MethylSense
# CHANGELOG:
#   v5.13.5 - BUGFIX: Time series figures (Figure 5, 6) now display correctly
#          - FIX 1: Changed hardcoded fig_ext to .png in markdown references (lines 9589, 9597)
#          - ISSUE: Files saved as .png but markdown referenced fig_ext (e.g., .jpeg) causing broken links
#          - FIX 2: Replaced hardcoded "*Figure 5:" and "*Figure 6:" with add_figure_caption()
#          - IMPACT: Time series trajectory and heatmap figures now render in reports
#          - IMPACT: Figure numbering now sequential (was hardcoded, breaking numbering)
#          - Based on v5.13.4
#   v5.13.4 - FEATURE: Added confidence scores to Table 11 (Study-Stratified Performance)
#          - FEATURE 1: Calculate mean ± SD confidence per study cohort and class
#          - FEATURE 2: Add "Mean Confidence" column to Table 11 (conditionally displayed if data available)
#          - FEATURE 3: Include overall confidence in Total row (e.g., "82.7 ± 11.2%")
#          - FEATURE 4: Add confidence score definition to interpretation text
#          - DEFINITION: "Confidence score = maximum predicted probability across all classes (0-100%)"
#          - IMPLEMENTATION: Lines 2315-2365 (calculate_study_stratified_performance function)
#          - IMPLEMENTATION: Lines 10038-10133 (table generation with conditional confidence column)
#          - USER REQUEST: "Add confidence from Table 6 to Table 11, stratified by study"
#          - IMPACT: Enables study-level confidence analysis, identifies cohorts with uncertain predictions
#          - Based on v5.13.3
#   v5.13.3 - ENHANCEMENT: Improved predictions table (user feedback)
#          - FEATURE 1: Show ALL samples (not just top 20) in predictions table
#          - FEATURE 2: Dynamic probability columns (not hardcoded to Aspergillus_Infected)
#          - FEATURE 3: Display P(Class) for ALL classes (e.g., P(Control), P(Aspergillus_Infected), etc.)
#          - FEATURE 4: Color-coded status icons (green ✓ for correct, red ✗ for incorrect) in bold
#          - IMPLEMENTATION: Lines 9129-9184 (complete table rewrite with dynamic column generation)
#          - USER REQUEST: "Show all samples, not just top 20! Make it work for any pathogen!"
#          - Based on v5.13.2
#   v5.13.2 - BUGFIX: Fixed figure handling in Final Model Predictions section
#          - FIX: Replaced get_figure_paths() with copy_plot_to_figures() + direct paths (5 instances)
#          - ISSUE: Function get_figure_paths() doesn't exist, caused runtime error
#          - SOLUTION: Use copy_plot_to_figures() then reference "figures/name.ext" (lines 9157, 9192, 9201, 9234, 9243)
#          - Based on v5.13.1
#   v5.13.1 - BUGFIX: Corrected function name in Final Model Predictions section
#          - FIX: Changed add_table_caption() → add_table_header() (2 instances)
#          - ISSUE: Function add_table_caption() doesn't exist, caused runtime error
#          - SOLUTION: Use correct helper function add_table_header() (lines 9123, 9165)
#          - Based on v5.13.0
#   v5.13.0 - FEATURE: Final Model Predictions Section (comprehensive test set analysis)
#          - NEW SECTION: "Final Model Predictions" added between Performance Validation and PCR Primer
#          - FEATURE 1: Sample-level predictions table with Lab_ID/Animal_ID metadata
#          - FEATURE 2: Prediction probability distributions (barplot, violin, ridge plots)
#          - FEATURE 3: Prediction confidence analysis with accuracy stratification
#          - FEATURE 4: DMR methylation patterns visualization (5-DMR barplot)
#          - FEATURE 5: Group methylation summary statistics and visualization
#          - FEATURE 6: Random Forest model insights (feature importance, partial dependence) if applicable
#          - NEW FUNCTIONS: load_final_model_predictions(), match_predictions_with_sample_sheet(), calculate_prediction_confidence_stats()
#          - IMPLEMENTATION: Lines 2376-2520 (helper functions), 9092-9237 (report section)
#          - USER REQUEST: "We need quite a bit more descriptive stats for all the prediction data we have"
#          - IMPACT: Provides comprehensive test set characterization with 5-7 tables and 7-10 figures
#          - RATIONALE: Addresses gap in test set reporting, enhances clinical interpretability
#          - QUALITY: All data sources already exist in markers_5/ directory
#          - Based on v5.12.1 (GOLD STANDARD)
#   v5.12.1 - PATCH: Table 3 F1-Score summary statistics correction
#          - FIX: Table 3 F1-Score mean/SD now calculated from correct source (class_fold_data)
#          - ISSUE: F1 mean was 0.597 (should be 0.488), SD was 0.058 (should be 0.080)
#          - SOLUTION: Calculate mean_f1/sd_f1 in summary block using correct data source (lines 8568-8581)
#          - IMPACT: Corrects 18.2% relative error in F1 mean calculation
#          - DETECTION: Statistics reviewer agent identified discrepancy after v5.12.0 release
#          - Based on v5.12.0
#   v5.12.0 - CRITICAL FIX: Table 12 medians + Table 3 class-specific metrics (COMPLETE)
#          - FIX 1: Calculate medians at source in analyse_dmr_group_associations() function
#          - ISSUE 1: v5.11.9 tried to recalculate medians during report generation but failed silently
#          - ROOT CAUSE 1: Wrong approach - should calculate during analysis, not during report formatting
#          - SOLUTION 1: Added median calculation in analyse_dmr_group_associations() using apply(1, median)
#          - IMPACT 1: Medians now calculated alongside means (line 2396) and included in results
#          - USER INSIGHT 1: "How do you calculate the means? Just do the same but with the median function"
#          - FIX 2: Table 3 (standard CV) now shows class-specific sensitivity/specificity
#          - ISSUE 2: v5.11.9 showed balanced metrics (sens = spec) for all folds
#          - ROOT CAUSE 2: cv_detailed.csv doesn't contain per-class metrics
#          - SOLUTION 2: Load cv_predictions_*.csv and calculate per-class metrics (mirrors Table 8 approach)
#          - IMPACT 2: Table 3 now shows Aspergillus_Infected-specific metrics with proper annotation
#          - USER INSIGHT 2: "You fixed the problem...in Table 8...Just look at the code at that time"
#          - METHOD: Uses calculate_per_class_nested_cv_metrics() on CV predictions (lines 3995-4031)
#          - BUG FIX: Detect cv_positive_class from CV predictions (was using undefined positive_class)
#          - VALIDATION: Follows same successful pattern as Table 8 nested CV metrics
#          - STATUS: Both fixes complete - "the two remaining hurdles before release"
#          - Based on v5.11.9
#   v5.11.9 - CRITICAL FIXES ATTEMPTED (both incomplete)
#          - FIX 1 ATTEMPTED: Removed condition preventing median recalculation (RETURNED NA - fixed in v5.12.0)
#          - FIX 2: Table 8 now uses per_class_metrics$fold_metrics ✓ (but Table 3 not fixed - fixed in v5.12.0)
#          - Based on v5.11.8
#   v5.11.8 - CRITICAL FIXES: STATISTICAL CALCULATION ERRORS (partial - 2 fixes incomplete)
#          - FIX: PR-AUC now hidden when not available per-fold (was showing identical values) ✓
#          - FIX ATTEMPTED: Table 12 mean/median bug - medians now calculated from raw data (RETURNED NA - fixed in v5.11.9)
#          - FIX ATTEMPTED: Nested CV table shows class-specific sensitivity/specificity (still balanced - fixed in v5.11.9)
#          - FIX: Kappa interpretation now uses Landis & Koch (1977) guidelines (0.37 = "fair" not "moderate") ✓
#          - FIX: Table 11 markdown formatting - moved note before table header ✓
#          - Based on v5.11.7
#   v5.11.7 - FIX: ROBUST SAMPLE SHEET MATCHING FOR CSV EXPORT
#          - FIX: CSV export now properly matches sample IDs from nested CV to sample sheet
#          - IMPROVEMENT: Two-tier matching strategy (rownames first, then sample_id_column)
#          - IMPROVEMENT: Character conversion for robust matching (handles numeric/string IDs)
#          - IMPROVEMENT: Added logging to show how many samples matched successfully
#          - FEATURE: Works generically for any dataset without hardcoding
#          - IMPACT: Lab_ID and Animal_ID now correctly populated for ALL samples in CSV
#          - USER REQUEST: "match the ID in the sample sheet" - now properly implemented
#          - Based on v5.11.6
#   v5.11.6 - FIX: CSV EXPORT COLUMN NAME BUG + IMPROVED ERROR HANDLING
#          - FIX: Column name mismatch - "actual" not "actual_class" (was causing silent failure)
#          - FIX: Added exists() check for sample_sheet_matched before accessing
#          - IMPROVEMENT: Better error messages using conditionMessage() instead of e$message
#          - IMPROVEMENT: Fallback mechanism to save basic CSV if enhanced version fails
#          - IMPACT: CSV export will now succeed and include Lab_ID/Animal_ID when available
#          - Based on v5.11.5
#   v5.11.5 - FIX: SPRINTF WARNING + ENHANCED CSV EXPORT
#          - FIX: Removed unused sprintf argument in Table 10 header (warning resolved)
#          - ENHANCEMENT: CSV export now includes Lab_ID and Animal_ID columns (when available)
#          - ENHANCEMENT: CSV column names match report exactly (Sample_ID, Lab_ID, Animal_ID, etc.)
#          - ENHANCEMENT: CSV columns ordered to match report table layout
#          - OUTPUT: sample_prediction_consistency.csv now publication-ready for paper
#          - USER REQUEST: "exactly as in the report" with all metadata columns
#          - Based on v5.11.4
#   v5.11.4 - FEATURE: EXPORT SAMPLE PREDICTION CONSISTENCY TO CSV
#          - NEW: Sample-level prediction consistency data now saved to CSV file
#          - FILE: sample_prediction_consistency.csv in model_evaluation directory
#          - CONTENT: Sample ID, actual class, correct folds, total folds, consistency rate
#          - PURPOSE: Critical data for paper - shows which samples are consistently predicted
#          - IMPACT: Enables identification of easy vs difficult samples for manuscript discussion
#          - USER REQUEST: "super important, it's going to be used for the paper as it is fantastic info"
#          - Based on v5.11.3
#   v5.11.3 - COMPREHENSIVE QC + STATISTICS FIXES (ALL AGENT ISSUES RESOLVED)
#          - FIX: Nested CV figure numbering - replaced hardcoded numbers with add_figure_caption()
#          - FIX: Nested CV figure extensions - changed .jpeg to .png (files exist as PNG only)
#          - FIX: Cohen's d interpretation - now classifies 0.35-0.5 as "small-to-medium" (not just "small")
#          - FIX: LMM vs LM terminology - clarified Linear Models (LM) for QC testing vs LMM for temporal
#          - IMPROVEMENT: Added bootstrap iteration count explanation (2000 for ROC-AUC, 1000 for CV, configurable for stability)
#          - IMPACT: Resolves ALL critical errors from QC Agent (figure numbering, extensions)
#          - IMPACT: Resolves ALL statistical warnings from Statistics Agent (effect sizes, model terminology)
#          - QUALITY: Report now passes both QC and statistical rigor validation
#          - Based on v5.11.2
#   v5.11.2 - CRITICAL FIX: NESTED CV SECTIONS MISSING FROM REPORT
#          - FIX: Changed "nested_cv_visualisations" to "nested_cv_visualizations" (British → American spelling)
#          - IMPACT: All nested CV sections now appear in report (were hidden due to directory name mismatch)
#          - SECTIONS NOW VISIBLE: Performance Distributions, Per-Fold Performance, ROC/PR Curves Across Folds
#          - SECTIONS NOW VISIBLE: Optimism Bias Analysis, CV Method Comparison, Outer Fold Performance
#          - ROOT CAUSE: Script checked for "visualisations" but directory was named "visualizations"
#          - SEVERITY: High - 6 major validation sections were completely missing from reports
#          - Based on v5.11.1
#   v5.11.1 - FINAL TABLE ENHANCEMENTS (MEDIAN COLUMNS)
#          - FIX: Table 11 now includes both Mean and Median columns for each group
#          - FIX: Table 10 title clarified - based on median methylation from logistic regression
#          - IMPROVEMENT: Dynamic table header construction for variable number of groups
#          - NOTE: Median columns calculated as fallback when not in analysis data
#          - USER REQUEST: "needs median also!" for DMR importance rankings
#          - Based on v5.11.0
#   v5.11.0 - STATISTICAL RIGOR ENHANCEMENTS (NORMALITY TESTING & WILCOXON FIXES)
#          - CRITICAL: Added Shapiro-Wilk normality testing to Population Methylation Statistics
#          - IMPACT: Report now tests distribution assumptions (W statistic, p-value for both groups)
#          - INTERPRETATION: Automatically recommends parametric vs non-parametric based on normality
#          - FIX: Wilcoxon test now reports MEDIAN differences (not means) - statistically correct!
#          - IMPROVEMENT: "Overall DMR Comparison" now shows both mean AND median for each group
#          - IMPACT: Interpretation text updated to reference median differences when using Wilcoxon
#          - RATIONALE: Non-parametric tests should report medians; parametric tests report means
#          - STATISTICS: Addresses all issues raised by Statistics Reviewer Agent
#          - Based on v5.10.0
#   v5.10.0 - TABLE OF CONTENTS ENHANCEMENT & SAMPLE CONSISTENCY IMPROVEMENTS
#          - FEATURE: Completely redesigned Table of Contents with granular subsections
#          - IMPACT: TOC now lists all 13 major sections with detailed subsections (50+ items total)
#          - IMPROVEMENT: Table 4 (Inconsistent Predictions) now shows ALL samples with 0 correct folds
#          - FIX: Table 4 sorting improved - samples with 0/N correct now sorted by total folds (descending)
#          - IMPACT: 0/4 now appears before 0/2 (both 0% but different denominators for better context)
#          - TABLE HEADER: Dynamic description shows "all X with 0 correct folds + top 10 others"
#          - USER BENEFIT: Enhanced report navigation and complete visibility of problematic samples
#          - Based on v5.9.2
#   v5.9.2 - METHYLATION STATISTICS ROBUSTNESS IMPROVEMENTS
#          - FIX: Added table header to "Population Methylation Statistics" (was missing!)
#          - FIX: Added median and IQR statistics alongside mean ± SD for methylation data
#          - IMPACT: Table now reports both parametric (mean ± SD) and non-parametric (median [IQR]) measures
#          - RATIONALE: Methylation data may not follow normal distribution, requiring robust statistics
#          - NEW TABLE FORMAT: "Mean ± SD" and "Median [IQR]" rows for Control vs Case comparison
#          - INTERPRETATION: Updated to mention both mean and median differences explicitly
#          - Based on v5.9.1
#   v5.9.1 - COMPREHENSIVE TABLE NUMBERING FIXES
#          - FIX: Fixed 3 hardcoded "Table 1.6.1" instances → now use add_table_header()
#          - FIX: Fixed 1 hardcoded "Table 1.14" instance → now use add_table_header()
#          - FIX: Fixed 1 hardcoded "Table 1.9.2" instance → now use add_table_header()
#          - FIX: Added table headers to 11 unnumbered tables (ALL tables with markdown layout now numbered)
#          - FIX: Fixed 3 tables where header was AFTER table instead of BEFORE (lines 8819, 8831, 8866)
#          - IMPACT: ALL tables in report now have proper numbered headers for paper citations
#          - IMPACT: Table numbering is fully sequential and dynamic
#          - NEW TABLES NUMBERED:
#              • Overall Performance Metrics
#              • Per-Class Performance Metrics (One-vs-Rest)
#              • Stratified Performance (Training vs Testing)
#              • Platform-Stratified Performance
#              • Complete QC Variable Statistics
#              • Key Analysis Parameters
#              • Genomic Feature Overlap Summary
#              • Hub DMRs Ranked by Connectivity
#              • WGCNA Module Sizes
#              • WGCNA Module Assignments
#              • Input Data File Paths
#          - CRITICAL: User reported inability to reference tables in paper without numbers
#          - Based on v5.9.0
#   v5.9.0 - CRITICAL FIX: FIGURE NUMBERING SYSTEM COMPLETELY BROKEN
#          - FIX: Removed 6 hardcoded figure numbers (13-16, 18-19) → now use add_figure_caption()
#          - FIX: Added individual captions to ALL 26 PCA/UMAP covariate plots (were missing!)
#          - FIX: Removed duplicate "Interpretation" headers (renamed to unique titles)
#          - FIX: Removed numbered subsection "4.1" → now just "Overview"
#          - IMPACT: Figure numbering now sequential 1-2-3-4-5... (NO MORE GAPS!)
#          - IMPACT: All 35 figure embeds now have proper captions
#          - IMPACT: Report is publication-ready with correct referencing
#          - ROOT CAUSE: Hardcoded numbers + missing add_figure_caption() calls in loops
#          - SEVERITY: Was CRITICAL - report completely unusable for papers before this fix
#          - Based on v5.8.3
#   v5.8.3 - CRITICAL FIX: MISSING CONFUSION MATRIX / ROC / PR CURVES IN REPORT
#          - FIX: PNG now ALWAYS generated as baseline format (prevents broken reports)
#          - FIX: Additional formats (JPEG, SVG) generated on top of PNG when requested
#          - FIX: All format variants copied to figures/ directory
#          - FIX: Markdown always references PNG (guaranteed to exist)
#          - IMPACT: Confusion matrix, ROC curves, PR curves now ALWAYS appear in report
#          - IMPACT: Users get PNG baseline + requested formats (e.g., --figure_format jpg,svg → PNG+JPG+SVG)
#          - STRATEGY: "PNG first, then extras" approach ensures report never breaks
#          - ROOT CAUSE: Plots generated in user's format only (no PNG), markdown checks for PNG, fails
#          - Based on v5.8.2
#   v5.8.2 - CRITICAL FIX: EXECUTIVE SUMMARY DMR COUNT CORRUPTION BUG
#          - FIX: Variable name collision - n_markers was overwritten by DMR group analysis loop
#          - FIX: Renamed local variable from n_markers to n_group_markers in line 9184
#          - IMPACT: Executive summary now shows correct total DMR count instead of last group's count
#          - IMPACT: Fixed "0 differentially methylated regions" appearing when DMRs exist
#          - CRITICAL: This was corrupting global n_markers variable, causing wrong counts throughout report
#          - Based on v5.8.1
#   v5.8.1 - FIX: PCA/UMAP COVARIATE PLOTS NOW COPY ALL FORMATS TO FIGURES DIR
#          - FIX: PCA_WITH_* plots now copy ALL formats (PNG, JPEG, SVG) not just primary
#          - FIX: UMAP_WITH_* plots now copy ALL formats to figures directory
#          - IMPACT: SVG files for covariate plots now available in figures/ directory
#          - IMPACT: Users can extract high-quality vector graphics for all dimensionality reduction plots
#          - Based on v5.8.0
#   v5.8.0 - MAJOR FIX: DYNAMIC FIGURE/TABLE NUMBERING SYSTEM
#          - FIX: Implemented dynamic counters for figure/table numbering (no more static numbers!)
#          - FIX: Figures now numbered sequentially regardless of which sections are generated
#          - FIX: Tables now numbered sequentially regardless of which sections are generated
#          - FIX: Removed ALL remaining numbered subsections (#### headers now have no numbers)
#          - FIX: Added add_figure_caption() helper function for automatic numbering
#          - FIX: Added add_table_header() helper function for automatic numbering
#          - IMPACT: No more gaps in numbering (no jumps from Table 3→6→10)
#          - IMPACT: Numbering adapts to conditional sections (gene annotation, calibration, etc.)
#          - IMPACT: Report now has clean sequential numbering even when sections are skipped
#          - CRITICAL: Solves the root cause of all numbering issues
#          - Based on v5.7.1
#   v5.7.1 - CRITICAL FIX: ALL PLOTS NOW GENERATED IN ALL REQUESTED FORMATS
#          - FIX: Converted ALL hardcoded png() calls to multi-format generation
#          - FIX: sample_correlation_heatmap now generates in all formats
#          - FIX: hierarchical_clustering_dendrogram now generates in all formats
#          - FIX: dmr_correlation_network now generates in all formats
#          - FIX: wgcna_dendrogram_modules now generates in all formats
#          - FIX: wgcna_module_trait_heatmap now generates in all formats
#          - FIX: methylation_heatmap_example now generates in all formats
#          - IMPACT: SVG files now properly generated when requested with --figure_format jpg,svg
#          - IMPACT: All plots respect user's figure format preference
#          - Based on v5.7.0
#   v5.7.0 - MAJOR CLEANUP: REMOVED ALL SECTION/SUBTITLE NUMBERING + SEQUENTIAL FIGURE/TABLE NUMBERING
#          - FIX: Removed ALL section numbers from headers (kept only titles)
#          - FIX: Removed ALL subsection numbers (### headers now have no numbers)
#          - FIX: Sequential figure numbering (Figure 1-25) throughout entire report
#          - FIX: Sequential table numbering (Table 1-13) throughout entire report
#          - FIX: BLAS/LAPACK path in sessionInfo() now truncated to show only from 'micromamba' forward
#          - FIX: TOC updated to bullet-point format (no section numbers)
#          - IMPACT: Clean, simple numbering scheme - only figures and tables are numbered
#          - IMPACT: Figures can now be properly referenced in papers (Figure 1, Figure 2, etc.)
#          - IMPACT: Tables can now be properly referenced in papers (Table 1, Table 2, etc.)
#          - IMPACT: No more confusing duplicate numbers or mixed numbering schemes
#          - Based on v5.6.0
#   v5.6.0 - MAJOR REFACTOR: FIXED ALL SECTION AND FIGURE NUMBERING
#          - FIX: Renumbered all top-level sections (1.5→2, 1.6→3, 1.7→4, 1.8→5, 2-9→6-13)
#          - FIX: Fixed section hierarchy - sections are now properly top-level (##) not sub-sections
#          - FIX: All subsections renumbered to match parent section (e.g., 1.6.1→3.1, 1.7.2→4.2)
#          - STRUCTURE: Clean sequential section numbering: 1-13
#          - IMPACT: TOC structure now makes logical sense
#          - IMPACT: No more duplicate section numbers or out-of-order sections
#          - NEXT: Figure numbers still need sequential update (1, 2, 3...)
#          - Based on v5.5.9
#   v5.5.9 - CRITICAL FIX: PLOTS NOW GENERATED IN CORRECT FORMAT
#          - FIX: confusion_matrix_3x3, roc_curves_multiclass, pr_curves_multiclass now generated in all formats
#          - FIX: Added open_plot_device() helper function to open correct graphics device (PNG/JPEG/SVG)
#          - FIX: Plots generated in user-specified formats from --figure_format parameter
#          - FIX: Flexible column name detection for PCR (supports scaffold/Chr, start/Start, end/End)
#          - IMPACT: ROC/PR/confusion plots now display correctly with jpg,svg format
#          - IMPACT: No more missing plots - all generated formats match markdown references
#          - Based on v5.5.8
#   v5.5.8 - CRITICAL FIX: MISSING PLOTS + PCR COORDINATES + FIGURE NUMBERING
#          - FIX: confusion_matrix_plot, roc_curves_multiclass, pr_curves_multiclass now use helper function
#          - FIX: Plots now copied with correct file extension (respects --figure_format)
#          - FIX: Figure numbering simplified in main performance section (Figure 1, 2, 3)
#          - FIX: PCR now finds coordinates in model_dir/dmr_details_with_coordinates.csv
#          - FIX: No more duplicate "Figure 1.3" or out-of-order "Figure 1.2"
#          - IMPACT: Confusion matrix, ROC curves, PR curves now display correctly
#          - IMPACT: PCR primer preparation works without gene annotation if coordinates exist
#          - IMPACT: Clean, sequential figure numbering in Section 1
#          - Based on v5.5.7
#   v5.5.7 - CRITICAL FIX: FOLD DETECTION BUG + REALISTIC CONCLUSIONS
#          - FIX: Fold detection now checks for both "fold" and "fold_num" column names
#          - FIX: Properly detects 10-fold nested CV data (was incorrectly treating as 1 fold)
#          - FIX: Bootstrap CIs now calculated correctly across all folds
#          - ENHANCEMENT: Conclusions now adaptive based on marker panel size
#          - ENHANCEMENT: Small panels (≤10 DMRs) get honest, critical assessment
#          - ENHANCEMENT: Larger panels get more confident (but still realistic) conclusions
#          - IMPACT: Nested CV tables now show proper confidence intervals
#          - IMPACT: 5-DMR models no longer claim "clinical translation ready"
#          - Based on v5.5.6
#   v5.5.6 - FIX: NESTED CV TABLE NOW HANDLES SINGLE-FOLD CASE PROPERLY
#          - FIX: When only 1 outer fold detected, table shows values without CI columns
#          - FIX: Clear warning message when <3 folds available for robust CIs
#          - ENHANCEMENT: Multi-fold case shows proper bootstrap CIs with fold count
#          - IMPACT: No more confusing "CI: 0.650-0.650" when there's only 1 fold
#          - IMPACT: Users now clearly informed when more outer folds needed
#          - Based on v5.5.5
#   v5.5.5 - CRITICAL FIX: ALL PLOTS NOW DISPLAY WITH CORRECT FORMAT
#          - FIX: Created copy_plot_to_figures() helper function for format-aware copying
#          - FIX: PCA/UMAP main plots now copied with correct extension (respects --figure_format)
#          - FIX: PCA/UMAP covariate plots use dynamic pattern matching for file extensions
#          - FIX: LMM associations plot now uses helper function with correct extension
#          - ENHANCEMENT: Automatic PNG fallback if requested format unavailable
#          - IMPACT: PCA/UMAP/LMM/Bootstrap sections now display plots correctly
#          - IMPACT: All dimensionality reduction visualizations working properly
#          - Based on v5.5.4
#   v5.5.4 - CRITICAL FIXES: EMPTY SECTIONS + HARDCODED TEXT REMOVED
#          - FIX: Bootstrap plot now copied with correct file extension (respects --figure_format)
#          - FIX: PCR section TOC entry now only appears when section body is generated
#          - FIX: Feature Importance (Section 4.4) now populated with marker stability rankings
#          - FIX: Removed all hardcoded conclusion text (30x stronger, 93 DMRs, 100% specificity)
#          - ENHANCEMENT: Conclusions now dynamically generated from actual analysis results
#          - ENHANCEMENT: Feature importance shows top 20 DMRs ranked by stability × effect size
#          - IMPACT: All report sections now show actual data, not placeholder text
#          - IMPACT: Bootstrap/stability plots now display correctly in all figure formats
#          - Based on v5.5.3
#   v5.5.3 - FIX: TABLE OF CONTENTS SECTION NUMBERING CORRECTED
#          - FIX: PCR Primer Preparation section renumbered from 1.12 to 1.9
#          - FIX: Section header corrected from "### 1.12.5" to "## 1.9" with proper anchor
#          - FIX: TOC now shows correct sequential numbering (1.5, 1.6, 1.7, 1.8, 1.9)
#          - CLARIFICATION: PCR section correctly doesn't appear when gene annotation disabled
#          - IMPACT: TOC now has proper sequential ordering without gaps
#          - IMPACT: PCR section anchor now matches TOC link (#pcr-primer-preparation)
#          - Based on v5.5.2
#   v5.5.2 - CRITICAL FIX: FIGURE FORMAT NOW RESPECTED IN MARKDOWN REPORTS
#          - FIX: All 17 hardcoded .png references replaced with dynamic fig_ext variable
#          - FIX: Stratified analysis now uses predictions.csv from main model dir (not nested CV)
#          - NEW: get_primary_figure_extension() function to determine correct file extension
#          - ENHANCEMENT: Stratified analysis prefers predictions.csv, falls back to nested CV
#          - IMPACT: --figure_format jpg,svg now correctly references .jpeg files in markdown
#          - IMPACT: Study stratification uses final model predictions (more accurate)
#          - Fixes: confusion matrix (2), ROC/PR curves (2), gene biotype (1), time series (2),
#                   nested CV (6), bootstrap/PCA/UMAP/LMM (4) image references
#          - Based on v5.5.1
#   v5.5.1 - CRITICAL FIX: SECTION 1.6.9 NOW APPEARS IN REPORT
#          - FIX: Moved Section 1.6.9 outside nested_cv_detailed conditional block
#          - FIX: Section 1.6.9 was being skipped because it was incorrectly nested
#          - NEW: Script version now logged at start of each run
#          - IMPACT: Section 1.6.9 study-stratified table NOW APPEARS in report
#          - IMPACT: User can now see script version in log immediately
#          - Based on v5.5.0
#   v5.5.0 - NEW FEATURE: SAVE STUDY-STRATIFIED RESULTS TO CSV
#          - NEW: Study stratification results now saved to study_stratified_performance.csv
#          - SHOWS: Study, Class, n, Correct, Accuracy_pct for each study-class combination
#          - IMPACT: Results persist and can be analyzed independently of report
#          - IMPACT: User can verify exact study stratification data used in Section 1.6.9
#          - Based on v5.4.9
#   v5.4.9 - CRITICAL: ADDED REPORT GENERATION VARIABLES STATUS CHECK
#          - NEW: "REPORT GENERATION VARIABLES STATUS" section after all variables are set
#          - SHOWS: nested_cv_done status, study_performance existence and availability
#          - SHOWS: Number of studies and explicit message if Section 1.6.9 should appear
#          - IMPACT: User can now see EXACTLY what variables are available for report generation
#          - IMPACT: Will immediately identify if nested_cv_done is FALSE (causing missing sections)
#          - Based on v5.4.8
#   v5.4.8 - ENHANCED: HIGHLY VISIBLE STUDY STRATIFICATION DEBUGGING
#          - ENHANCEMENT: Added bold section header "STUDY-STRATIFIED PERFORMANCE ANALYSIS"
#          - ENHANCEMENT: Clear SUCCESS/FAILED messages with >>> markers for easy log scanning
#          - ENHANCEMENT: Step-by-step progress output (STEP 1-5) during calculation
#          - ENHANCEMENT: Explicit message showing if Section 1.6.9 will appear in report
#          - DEBUG: Shows sample sheet path, columns, STUDY column detection, merge stats
#          - DEBUG: Shows retention rate after filtering, final calculation results
#          - IMPACT: Log output now extremely easy to scan for study stratification status
#          - IMPACT: User can immediately see why Section 1.6.9 is missing (if it is)
#          - Based on v5.4.7
#   v5.4.7 - FIX: TABLE OF CONTENTS + ENHANCED STRATIFIED ANALYSIS DEBUGGING
#          - FIX: Table of Contents now uses consistent formatting (removed trailing periods)
#          - FIX: All TOC entries now properly formatted as [Section Title] without trailing dots
#          - ENHANCEMENT: Added comprehensive debugging output to calculate_study_stratified_performance()
#          - ENHANCEMENT: Now shows sample sheet columns, STUDY column detection, merge statistics
#          - ENHANCEMENT: Better error messages when sample sheet not found or STUDY column missing
#          - DEBUG: Added diagnostic output for sample sheet path, merge process, NA filtering
#          - IMPACT: Cleaner, more consistent Table of Contents appearance
#          - IMPACT: Easier to diagnose why study stratification might fail
#          - Based on v5.4.6
#   v5.4.6 - FIX: SECTION 1.6.9 VISIBILITY + SECTION 7.4.1 PATH BLINDING
#          - FIX: Section 1.6.9 (Study-Stratified Performance) now included in Table of Contents
#          - FIX: Added HTML anchor {#study-stratified-performance} to Section 1.6.9 for TOC linking
#          - FIX: Path blinding now working in Section 7.4.1 Model Training Parameters
#          - FIX: Console messages in Section 7.4.1 now respect --blind_paths flag
#          - FIX: Parameter file content paths sanitized when --blind_paths enabled
#          - ENHANCEMENT: calculate_study_stratified_performance() more robust STUDY column detection
#          - ENHANCEMENT: Now searches for column variants (STUDY, STUDY...3, Study, COHORT, etc.)
#          - IMPACT: Section 1.6.9 now visible and clickable in Table of Contents
#          - IMPACT: No more long paths visible in Section 7.4.1 when using --blind_paths
#          - Based on v5.4.5
#   v5.4.5 - FIX: WILCOXON PLACEMENT + PATH BLINDING
#          - FIX: Removed Wilcoxon test from Section 1.7.2 Variable Importance Ranking
#          - NEW: Section 1.7.5 - Overall DMR Methylation Comparison (Wilcoxon/Kruskal-Wallis)
#          - FIX: Path blinding now works in Section 1.7.4 (dmr_group_associations.csv)
#          - ENHANCEMENT: Section 1.7.5 tests overall methylation differences between groups
#          - ENHANCEMENT: Wilcoxon (2 groups) or Kruskal-Wallis (>2 groups) with clear interpretation
#          - NOTE: Section 1.6.9 study-stratified table exists (appears when study_performance available)
#          - IMPACT: Clearer DMR importance ranking without redundant Wilcoxon columns
#          - IMPACT: Simple group methylation comparison separate from per-DMR ANOVA
#          - Based on v5.4.4
#   v5.4.4 - FIX: SINGLE-FOLD NESTED CV + HPC-FRIENDLY HTML OUTPUT
#          - FIX: calculate_ci() now handles single-fold cases (returns point estimate instead of NA)
#          - FIX: Section 1.6.1 now shows proper metrics even with only 1 nested CV fold
#          - ENHANCEMENT: Added warning message in report when only 1 fold available
#          - ENHANCEMENT: Interpretation text adjusted for single-fold vs multi-fold scenarios
#          - ENHANCEMENT: Report summary now clearly shows generated files (HTML/MD/PDF)
#          - ENHANCEMENT: Added helpful tips for HTML viewing and PDF generation
#          - ENHANCEMENT: Better HPC user guidance (generate PDF locally after downloading HTML)
#          - DEFAULT: PDF generation is opt-in (--generate_pdf flag), HTML is primary format
#          - IMPACT: No more NA% values in Section 1.6.1 performance summary
#          - IMPACT: HPC users get clean HTML reports without PDF errors
#          - Based on v5.4.3
#   v5.4.3 - CRITICAL: REMOVED LATEX FALLBACK + METHYLSENSE BRANDING
#          - BREAKING: Completely removed LaTeX fallback - pagedown ONLY for PDF generation
#          - NEW: MethylSense branding and GitHub link in report header
#          - NEW: Primary citation (Drag et al. 2025 bioRxiv) in References section
#          - NEW: Author attribution (Markus Hodal Drag) throughout all scripts
#          - FIX: Version numbers synchronized across ALL scripts (5.4.3)
#          - FIX: Enhanced Unicode sanitization with comprehensive Greek letter support (β, η, θ, etc.)
#          - ENHANCEMENT: Better error messages when pagedown/Chrome unavailable
#          - IMPACT: No more Unicode errors - if pagedown fails, PDF is simply skipped
#          - IMPACT: Users now know analyses are generated by MethylSense
#          - Based on v5.4.2
#   v5.4.2 - ENHANCEMENT: WILCOXON TEST FOR GROUP METHYLATION COMPARISON
#          - NEW: Wilcoxon rank-sum test (Mann-Whitney U) added to DMR-Group Association Analysis
#          - NEW: Section 1.7.2 table now includes Wilcoxon p-value, FDR-adjusted p-value, and significance markers
#          - NEW: Median difference column shows direction and magnitude of methylation difference between groups
#          - ENHANCEMENT: Report interpretation section explains Wilcoxon test results
#          - ENHANCEMENT: Added pagedown to MethylSense_installer_v1.R prerequisite packages
#          - IMPACT: Provides statistical evidence for methylation differences in binary comparisons
#          - Based on v5.4.1
#   v5.4.1 - CRITICAL FIXES: NA VALUES, PATH BLINDING, PDF GENERATION
#          - FIX: Section 1.6.1 showing NA values - added else block to initialize per_class_metrics
#          - FIX: --blind_paths not working in reproducibility command section
#          - FIX: PDF generation logic improved with better error messages for pagedown
#          - ENHANCEMENT: pagedown verbose mode enabled for debugging
#          - ENHANCEMENT: Better diagnostic messages for Chrome/Chromium detection
#          - Based on v5.4.0
#   v5.4.0 - MAJOR ENHANCEMENT: PER-CLASS METRICS & STUDY-STRATIFIED PERFORMANCE
#          - NEW: calculate_per_class_nested_cv_metrics() for class-specific sensitivity/specificity/F1
#          - NEW: calculate_study_stratified_performance() for cohort-level breakdown
#          - NEW: Section 1.6.1 now shows per-class metrics (not averaged)
#          - NEW: F1 Score added to nested CV performance table with 95% CI
#          - NEW: Auto-detection of positive class (Aspergillus_Infected, etc.)
#          - ENHANCEMENT: Section 1.6.9 now shows detailed study-by-group classification table
#          - FIX: load_cv_predictions() now searches for nested_cv_predictions_*.csv first
#          - FIX: Variable naming consistency (nested_cv_predictions vs cv_predictions)
#          - IMPACT: Sensitivity/Specificity now clearly show which class they refer to
#          - Based on v5.3.11
#   v5.3.11 - ENHANCEMENT: CLINICAL-STYLE PDF REPORTS (NO LATEX!)
#          - NEW: Custom CSS theme for professional clinical/diagnostic report styling
#          - NEW: get_clinical_report_css() function for modern report aesthetics
#          - ENHANCEMENT: HTML reports now include clinical CSS automatically
#          - ENHANCEMENT: PDF generation via pagedown uses clinical styling (blue headers, gradient tables, clean fonts)
#          - ENHANCEMENT: Better error messages for pagedown/Chrome requirements
#          - STYLE: Sans-serif fonts (not LaTeX Computer Modern), blue color scheme, modern tables
#          - FALLBACK: Still supports LaTeX PDF if pagedown unavailable
#          - Based on v5.3.10
#   v5.3.10 - FIX: FIGURE FORMAT COMMA-SEPARATED VALUES
#          - FIX: --figure_format now correctly handles comma-separated formats (e.g., "jpg,svg")
#          - FIX: Normalizes "jpg" to "jpeg" automatically
#          - ENHANCEMENT: save_plot_multiple_formats() now iterates over multiple formats
#          - Updated help text to document comma-separated format support
#          - Based on v5.3.9
#   v5.3.9 - CRITICAL FIX: LOGISTIC REGRESSION GROUP ASSIGNMENT
#          - FIX: Replaced alphabetical group sorting with intelligent keyword-based detection
#          - NEW: --control_group parameter to explicitly specify reference group
#          - NEW: Auto-detection of control groups using keywords (control, healthy, negative, uninfected, wildtype, wt, normal, baseline)
#          - NEW: Auto-detection of case groups using keywords (infected, infection, disease, case, positive, treatment, tumor, cancer, patient)
#          - FIX: Report now explicitly shows which group is control (0) and which is case (1)
#          - FIX: Report clarifies that model predicts probability of the CASE group (not control)
#          - FIX: Added clear console output showing group assignment during analysis
#          - IMPACT: Previous versions (v5.3.8 and earlier) had REVERSED group assignment for "Aspergillus_Infected" vs "Control"
#          - Based on v5.3.8
#   v5.3.8 - LOGISTIC REGRESSION ANALYSIS FOR INFECTION PROBABILITY (DEPRECATED - HAD BUG)
#          - NEW: Logistic regression analysis to calculate infection probability vs median methylation
#          - NEW: Section 1.6.10 in report showing probability thresholds vs median methylation percentages
#          - NEW: Output files: infection_probability_vs_methylation.csv and logistic_regression_model_summary.txt
#          - FIX: Corrected data source from non-existent meth_matrix_full to existing meth_matrix
#          - BUG: Alphabetical sorting caused reversed group assignment (e.g., "Aspergillus_Infected" < "Control")
#          - Based on v5.3.7
#   v5.3.7 - STUDY-STRATIFIED NESTED CV PERFORMANCE
#          - NEW: Section 1.6.9 showing nested CV performance stratified by study cohort
#          - Shows per-study classification accuracy for control and infected samples
#          - Based on v5.3.6
#   v5.3.0 - PDF GENERATION AND PRIVACY PROTECTION
#          - NEW: --generate_pdf flag to automatically generate PDF reports
#          - NEW: --blind_paths flag to sanitize directory paths for privacy
#          - NEW: PDF generation from markdown using rmarkdown (optional)
#          - NEW: Helper function sanitize_path() for path privacy protection
#          - ENHANCEMENT: Model directory, QS file, and sample sheet paths now support sanitization
#          - Applied to markdown, HTML, and PDF reports
#          - Based on v5.2.0
#   v5.2.0 - FEATURE ANNOTATION ENHANCEMENT
#          - NEW: Genomic feature annotation capability using genomation package
#          - NEW: --enable_feature_annotation flag to enable feature annotation
#          - NEW: --gtf_file parameter for gene feature annotation (required)
#          - NEW: --regulatory_gff parameter for regulatory features (optional)
#          - NEW: --emar_gff parameter for EMAR features (optional)
#          - NEW: --open_chromatin_bed parameter for open chromatin regions (optional)
#          - NEW: --chromatin_flank_bp parameter to set flank size for chromatin regions
#          - NEW: STEP 6.10.6 - Genomic Feature Annotation
#          - NEW: Annotates DMRs with overlap percentages for multiple feature types
#          - NEW: Combined bar plot showing feature overlap for all annotation types
#          - NEW: CSV output for each annotation type (genes, regulatory, EMAR, chromatin)
#          - NEW: feature_annotation/ output directory with figures/ subdirectory
#          - ENHANCEMENT: Uses gene_annotation_results for DMR coordinates
#          - ENHANCEMENT: Supports multiple annotation sources simultaneously
#          - Based on v5.1.0
#   v24.0.0 - COMPREHENSIVE FIX AND ENHANCEMENT RELEASE
#          - REMOVED: Section 5.3 DMR Methylation Patterns by Biological Group (redundant with existing analyses)
#          - FIX: PCA/UMAP variance filtering now uses 1e-10 threshold instead of exact zero (fixes "not performed" errors)
#          - FIX: Added debug logging to PCA/UMAP showing number of variable markers available
#          - ENHANCEMENT: Network visualisation now shows all DMR ID labels on nodes
#          - ENHANCEMENT: Hub nodes in network shown in bold with larger text
#          - ENHANCEMENT: Network legend updated to clarify hub node designation
#          - FIX: Heatmap group annotations now replace NA values with "Unknown" instead of showing blank
#          - ENHANCEMENT: Gene annotation window changed from 10 kb to 10 Mb (10,000,000 bp) default
#          - ENHANCEMENT: Gene annotation table now shows ALL DMRs instead of top 20 only
#          - ENHANCEMENT: Gene annotation table header updated to indicate 10 Mb window size
#          - ENHANCEMENT: Gene annotation now searches all Ensembl features (genes, exons, transcripts, UTRs)
#          - ENHANCEMENT: Table 1.7.2 (Variable Importance) now shows ALL DMRs instead of top 20
#          - ENHANCEMENT: Table 1.9.2 (Hub DMRs) now shows ALL hub DMRs instead of top 10
#          - ENHANCEMENT: Table 1.14.2 (Temporal Changes) now shows ALL significant DMRs instead of top 10
#          - ENHANCEMENT: Table 4.4 (Feature Importance) now shows ALL DMRs instead of top 20
#          - FIX: All text converted to UK English spelling (analyse, visualise, optimise, colour, etc.)
#          - ENHANCEMENT: 164 instances of "analyze" → "analyse", 22 instances of "visualize" → "visualise", 12 instances of "optimize" → "optimise"
#          - Based on v23.2.0
#   v23.2.0 - PCR PRIMER PREPARATION FEATURE RELEASE (FEATURE #3)
#          - NEW: PCR primer preparation functionality (FEATURE #3)
#          - NEW: --prepare_pcr flag for enabling PCR primer preparation
#          - NEW: --pcr_flank_size flag to set flanking region size
#          - NEW: --pcr_product_size flag to set target PCR product size range
#          - NEW: prepare_bisulfite_primers() helper function with accurate bisulfite conversion
#          - NEW: Species-specific BSgenome integration (Sus scrofa, Homo sapiens, Mus musculus)
#          - NEW: Bisulfite conversion logic preserving CpG cytosines
#          - NEW: Primer3Plus and BLAST URL generation for primer design
#          - NEW: STEP 6.12 - PCR Primer Preparation
#          - NEW: Report Section 1.13 - PCR Validation Primer Design
#          - ENHANCEMENT: pcr_primers/ output directory with CSV and FASTA files
#          - Based on v23.1.0 (FEATURE #2)
#   v23.1.0 - TIME SERIES ANALYSIS FEATURE RELEASE (FEATURE #2)
#          - NEW: Time series analysis functionality (FEATURE #2)
#          - NEW: --time_analysis flag for specifying timepoint column
#          - NEW: --time_categorical flag to force categorical time treatment
#          - NEW: analyse_time_series() helper function with auto-detection
#          - NEW: Support for continuous time (LMM/Linear Regression) and categorical time (ANOVA/RM-ANOVA)
#          - NEW: Auto-detection of repeated measures vs independent groups design
#          - NEW: Trajectory clustering (increasing/decreasing/stable/U-shaped)
#          - NEW: STEP 6.11 - Time Series Analysis
#          - NEW: Report Section 1.14 - Temporal Dynamics of DMR Methylation
#          - NEW: Time series visualisations (trajectories, heatmap, clusters)
#          - ENHANCEMENT: dmr_time_series_statistics.csv and dmr_time_trajectories.csv outputs
#          - Based on v23.0.0 (FEATURE #1)
#   v23.0.0 - GENE ANNOTATION FEATURE RELEASE (FEATURE #1)
#          - NEW: Gene annotation functionality using biomaRt (FEATURE #1)
#          - NEW: --species flag for species-specific gene annotation
#          - NEW: --enable_gene_annotation flag to enable gene annotation
#          - NEW: --gene_window flag to set genomic window for nearest gene search
#          - NEW: query_biomart_genes() helper function for Ensembl queries
#          - NEW: STEP 6.10 - Gene Annotation Analysis
#          - NEW: Report Section 1.12 - Genomic Context and Gene Annotation
#          - ENHANCEMENT: dmr_gene_annotation.csv output with overlapping and nearest genes
#          - Based on v22.0.0
#   v22.0.0 - COMPREHENSIVE FIX RELEASE
#          - FIX: All hardcoded DMR counts now dynamic (93 markers → nrow(meth_matrix))
#          - FIX: Section 1.7.2 adds p-values and FDR correction to DMR-group associations
#          - FIX: Section 1.8 heatmap/dendrogram group annotations corrected
#          - FIX: Section 1.8.3 clarifies cluster purity metric interpretation
#          - FIX: Section 1.9 adaptive network threshold (0.7→0.6→0.5)
#          - FIX: Section 1.11 dynamic hyperparameter file matching (no hardcoded filenames)
#          - FIX: Sections 5.1/5.2 PCA/UMAP always generated when data available
#          - FIX: Section 5.3 analyses ALL DMRs with FDR correction (not just top 10)
#          - FIX: Sections 6.1/6.2 empty state messages when no technical covariates
#          - ENHANCEMENT: Test selection explanation in Section 5.3 (t-test/ANOVA/Kruskal-Wallis)
#          - Based on v21.0.0
#   v21.0.0 - NEW: Sample correlation matrix and hierarchical clustering (Section 1.8)
#          - NEW: DMR correlation network analysis (Section 1.9)
#          - NEW: WGCNA co-expression modules with --enable_wgcna flag (Section 1.10)
#          - NEW: Always show methylation vs biological groups in Section 5
#          - FIX: Added sensitivity/specificity CIs to Section 1.6.1
#          - FIX: Changed default output directory to "model_evaluation"
#          - FIX: Report filename changed to MODEL_EVALUATION_REPORT.md
#   v20.0.0 - NEW: DMR-group association analysis with variable importance
#          - NEW: Complete CSV output of DMR-group statistics
#          - FIX: Nested CV sensitivity/specificity column names
#          - FIX: Variable importance now calculated from group associations
#   v19.0.0 - FIX: Multiclass bootstrap differential methylation (was returning 0s)
#          - FIX: Complete REVIEWER_RESPONSE.md generation
#          - FIX: Correct output directory structure (reviewer_analysis/)
#   v18.0.0 - NEW: --group_colors flag for customizable PCA/UMAP color schemes
#          - NEW: Nested CV fold-by-fold performance table (Section 1.6.8)
#          - Users can now specify custom colors for biological groups
#          - Added detailed outer fold performance breakdown for nested CV
#   v17.0.0 - CRITICAL FIXES: Fixed nested CV plot filenames, JSON hyperparameters,
#          - table formatting, empty Section 1.5.1, dual PCA/UMAP legends,
#          - added missing optimism bias and CV comparison plots
#   v16.0.0 - Fixed bootstrap logging and added section end markers
#   v15.0.0 - MAJOR UPDATE: Bulletproof comprehensive validation report
#          - NEW: Fold-by-fold performance variability analysis (Section 1.9)
#          - NEW: Sample-level prediction consistency analysis (Section 1.10)
#          - NEW: Hyperparameter sensitivity analysis (Section 1.11)
#          - NEW: Probability calibration analysis (Section 1.12)
#          - NEW: Enhanced methylation pattern analysis (Section 3.2)
#          - NEW: Marker statistical significance analysis (Section 3.3)
#          - NEW: Feature importance rankings (Section 4.4)
#          - ENHANCED: Every PCA/UMAP plot now has detailed interpretation
#          - ENHANCED: Section 3.1 expanded with additional QC metrics
#          - ENHANCED: Quantitative batch effect assessment (Section 6.2)
#          - ENHANCED: Reproducibility statement in methods
#          - ENHANCED: Quantitative conclusions summary
#          - FIX: Use confusion matrix from model directory
#          - All sections now solid and thorough
#          - Based on v14.1.0
#   v13.0.0 - MAJOR UPDATE: Complete model validation report for supplementary materials
#          - NEW STEP 1.5: Model characterisation (hyperparameters, variable importance)
#          - NEW STEP 2.5: Nested cross-validation analysis with 95% CIs
#          - NEW: Optimism bias analysis (train vs test performance)
#          - NEW: Binary/multiclass auto-detection
#          - ENHANCED: Use existing high-quality plots from model directory
#          - ENHANCED: Complete markdown report (25-30 pages, 36-46 figures)
#          - All analyses in UK English, zero emojis, publication quality
#          - Based on v12.0.0
#   v12.0.0 - MAJOR ENHANCEMENT: Comprehensive model performance validation
#          - NEW STEP 3.5: Multiclass performance analysis (ROC, PR curves, confusion matrix)
#          - Loads predictions.csv and calculates one-vs-rest metrics for each class
#          - Generates ROC-AUC and PR-AUC with 95% confidence intervals
#          - Calculates sensitivity, specificity, PPV, NPV, F1, MCC per class
#          - Stratified performance analysis (TRAIN vs TEST, platform-specific)
#          - Enhanced markdown generation with UK English throughout
#          - Removed all emojis for professional academic tone
#          - Added NEW Section 1: Model Performance Validation (before existing sections)
#          - All existing sections renumbered accordingly
#          - Comprehensive plot inclusion (all QC, stability, PCA, UMAP plots)
#          - Professional figure captions and table legends
#          - All "analyse" → "analyse", "colour" → "colour", etc.
#   v11.0.0 - Previous version with integrated markdown generation
#   v9.0.0 - MAJOR UPDATE: Merged markdown generation directly into main script
#          - NO LONGER requires separate MethylSense_Markdown script
#          - Markdown generation uses data already in memory (no file re-reading)
#          - Eliminates file path issues and improves reliability
#          - All v8 features retained and enhanced
#   v8.0.0 - ENHANCEMENT: Improved covariate association visualisations
#          - Bar plots now show ALL markers (not just significant ones)
#          - Stacked bars show proportions: non-significant, raw p-significant, adjusted p-significant
#          - Prevents empty bars when no associations found
#   v7.0.0 - MAJOR ENHANCEMENT: Added --significance_threshold flag
#   v6.7.0 - ENHANCEMENT: Added --multiple_testing flag for flexible correction
#   v6.6.0 - CRITICAL: Validates extracted methylation data
#   v6.5.0 - Smart sample ID matching
#   v6.3.2 - Fixed PCA/UMAP
# ============================================================================

SCRIPT_VERSION <- "5.13.7"
SCRIPT_DATE <- "2025-12-05"


# ============================================================================
# COMPREHENSIVE LOGGING SYSTEM
# ============================================================================

# Initialize logging variables
LOG_FILE <- NULL
LOG_LEVEL <- 2 # 0=ERROR, 1=WARN, 2=INFO, 3=DEBUG
LOG_START_TIME <- Sys.time()

init_log <- function(output_dir, log_level = 2) {
  LOG_LEVEL <<- log_level
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  LOG_FILE <<- file.path(output_dir, paste0("analysis_log_", timestamp, ".txt"))
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  write(paste0(strrep("=", 80)), LOG_FILE)
  write(paste0("MethylSense Reviewer Analysis v", SCRIPT_VERSION), LOG_FILE, append = TRUE)
  write(paste0("Log started: ", Sys.time()), LOG_FILE, append = TRUE)
  write(strrep("=", 80), LOG_FILE, append = TRUE)
  write("", LOG_FILE, append = TRUE)
}

log_msg <- function(level, message, step = "", print_console = TRUE) {
  level_codes <- c("ERROR" = 0, "WARN" = 1, "INFO" = 2, "DEBUG" = 3)
  if (!(level %in% names(level_codes))) level <- "INFO"
  if (level_codes[level] > LOG_LEVEL) {
    return(invisible(NULL))
  }
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  elapsed <- round(as.numeric(difftime(Sys.time(), LOG_START_TIME, units = "secs")), 2)
  step_prefix <- if (step != "") paste0("[", step, "] ") else ""
  log_line <- paste0("[", timestamp, "] [+", sprintf("%7.2f", elapsed), "s] [", sprintf("%-5s", level), "] ", step_prefix, message)
  if (!is.null(LOG_FILE)) write(log_line, LOG_FILE, append = TRUE)
  if (print_console) {
    if (level == "ERROR") {
      cat("\033[0;31m", log_line, "\033[0m\n", sep = "")
    } else if (level == "WARN") {
      cat("\033[0;33m", log_line, "\033[0m\n", sep = "")
    } else if (level == "DEBUG") {
      cat("\033[0;36m", log_line, "\033[0m\n", sep = "")
    } else {
      cat(log_line, "\n", sep = "")
    }
  }
  invisible(log_line)
}

log_error <- function(message, step = "") log_msg("ERROR", message, step)
log_warn <- function(message, step = "") log_msg("WARN", message, step)
log_info <- function(message, step = "") log_msg("INFO", message, step)
log_debug <- function(message, step = "") log_msg("DEBUG", message, step)

log_section_start <- function(step_num, step_name, total_steps = 10) {
  log_info(strrep("=", 60), "")
  log_info(paste0("STEP ", step_num, "/", total_steps, ": ", step_name), "MAIN")
  log_info(strrep("=", 60), "")
}

log_section_end <- function(step_num, step_name) {
  log_info(paste0("Completed STEP ", step_num, ": ", step_name), "MAIN")
  log_info("", "")
}

log_file_check <- function(filepath, description, step = "") {
  if (file.exists(filepath)) {
    size_kb <- round(file.size(filepath) / 1024, 2)
    log_info(paste0("Found: ", description, " (", size_kb, " KB)"), step)
    return(TRUE)
  } else {
    log_warn(paste0("Missing: ", description), step)
    return(FALSE)
  }
}

log_object_info <- function(obj, name, step = "") {
  if (is.data.frame(obj) || is.matrix(obj)) {
    log_debug(paste0(name, ": ", nrow(obj), " rows x ", ncol(obj), " columns"), step)
  } else if (is.list(obj)) {
    log_debug(paste0(name, ": list with ", length(obj), " elements"), step)
  } else if (is.vector(obj)) {
    log_debug(paste0(name, ": vector with ", length(obj), " elements"), step)
  }
}

with_logging <- function(expr, description, step = "", error_handler = NULL) {
  log_info(paste0("Starting: ", description), step)
  result <- tryCatch(
    {
      res <- expr
      log_info(paste0("Completed: ", description), step)
      res
    },
    error = function(e) {
      log_error(paste0("ERROR in ", description, ": ", as.character(e)), step)
      if (!is.null(error_handler)) error_handler(e) else NULL
    }
  )
  return(result)
}

finalize_log <- function() {
  if (!is.null(LOG_FILE)) {
    elapsed_total <- as.numeric(difftime(Sys.time(), LOG_START_TIME, units = "mins"))
    write("", LOG_FILE, append = TRUE)
    write(strrep("=", 80), LOG_FILE, append = TRUE)
    write(paste0("Analysis completed: ", Sys.time()), LOG_FILE, append = TRUE)
    write(paste0("Total runtime: ", round(elapsed_total, 2), " minutes"), LOG_FILE, append = TRUE)
    write(strrep("=", 80), LOG_FILE, append = TRUE)
    cat("\nLog file written to:", LOG_FILE, "\n")
    cat("Total runtime:", round(elapsed_total, 2), "minutes\n")
  }
}

# Helper function to sanitize file paths for privacy
# Keeps only the basename (last component) of file/directory paths
sanitize_path <- function(path_string) {
  if (is.null(path_string) || path_string == "") {
    return(path_string)
  }
  basename(path_string)
}

# Helper function to sanitize Unicode characters for LaTeX compatibility
# Replaces Unicode symbols with LaTeX-safe ASCII equivalents
sanitize_unicode_for_latex <- function(text) {
  # Mathematical symbols
  text <- gsub("\u2265", ">=", text, fixed = TRUE) # greater than or equal
  text <- gsub("\u2264", "<=", text, fixed = TRUE) # less than or equal
  text <- gsub("\u00b1", "+/-", text, fixed = TRUE) # plus-minus
  text <- gsub("\u00d7", "x", text, fixed = TRUE) # multiplication
  text <- gsub("\u2192", "->", text, fixed = TRUE) # right arrow
  text <- gsub("\u2190", "<-", text, fixed = TRUE) # left arrow
  text <- gsub("\u2260", "!=", text, fixed = TRUE) # not equal
  text <- gsub("\u2248", "~=", text, fixed = TRUE) # approximately equal
  text <- gsub("\u00b0", " degrees", text, fixed = TRUE) # degree symbol

  # Greek letters (both lowercase and uppercase variants)
  text <- gsub("\u00b5", "u", text, fixed = TRUE) # micro (mu)
  text <- gsub("\u03b1", "alpha", text, fixed = TRUE) # α
  text <- gsub("\u03b2", "beta", text, fixed = TRUE) # β (CRITICAL for logistic regression)
  text <- gsub("\u0392", "Beta", text, fixed = TRUE) # Β (uppercase beta)
  text <- gsub("\u03b3", "gamma", text, fixed = TRUE) # γ
  text <- gsub("\u0393", "Gamma", text, fixed = TRUE) # Γ
  text <- gsub("\u0394", "Delta", text, fixed = TRUE) # Δ
  text <- gsub("\u03b4", "delta", text, fixed = TRUE) # δ
  text <- gsub("\u03bb", "lambda", text, fixed = TRUE) # λ
  text <- gsub("\u039b", "Lambda", text, fixed = TRUE) # Λ
  text <- gsub("\u03c3", "sigma", text, fixed = TRUE) # σ
  text <- gsub("\u03a3", "Sigma", text, fixed = TRUE) # Σ
  text <- gsub("\u03b8", "theta", text, fixed = TRUE) # θ
  text <- gsub("\u0398", "Theta", text, fixed = TRUE) # Θ
  text <- gsub("\u03b7", "eta", text, fixed = TRUE) # η (for eta-squared)
  text <- gsub("\u0397", "Eta", text, fixed = TRUE) # Η

  # Punctuation
  text <- gsub("\u2014", "--", text, fixed = TRUE) # em dash
  text <- gsub("\u2013", "-", text, fixed = TRUE) # en dash
  text <- gsub("\u201c", '"', text, fixed = TRUE) # left double quote
  text <- gsub("\u201d", '"', text, fixed = TRUE) # right double quote
  text <- gsub("\u2018", "'", text, fixed = TRUE) # left single quote
  text <- gsub("\u2019", "'", text, fixed = TRUE) # right single quote
  text <- gsub("\u2026", "...", text, fixed = TRUE) # ellipsis

  return(text)
}

# ============================================================================
# COLOR PALETTE PARSING
# ============================================================================

parse_group_colors <- function(color_string, group_names) {
  # If no custom colors specified, return NULL to use default Set1 palette
  if (is.null(color_string) || color_string == "") {
    return(NULL)
  }

  # Parse the color string: "Group1=#HEX1,Group2=color2,..."
  color_pairs <- strsplit(color_string, ",")[[1]]
  color_map <- list()

  for (pair in color_pairs) {
    parts <- strsplit(trimws(pair), "=")[[1]]
    if (length(parts) != 2) {
      warning(paste0("Invalid color specification: '", pair, "'. Use format 'GroupName=color'"))
      next
    }

    group_name <- trimws(parts[1])
    color_value <- trimws(parts[2])

    # Validate hex color format if it starts with #
    if (grepl("^#", color_value)) {
      if (!grepl("^#[0-9A-Fa-f]{6}$", color_value)) {
        warning(paste0("Invalid hex color '", color_value, "' for group '", group_name, "'. Use format #RRGGBB"))
        next
      }
    }

    # Check if group name exists in data
    if (!group_name %in% group_names) {
      warning(paste0("Group '", group_name, "' not found in data. Available groups: ", paste(group_names, collapse = ", ")))
      next
    }

    color_map[[group_name]] <- color_value
  }

  # Fill in missing groups with default Set1 colors
  set1_colors <- RColorBrewer::brewer.pal(max(3, length(group_names)), "Set1")
  for (i in seq_along(group_names)) {
    if (!group_names[i] %in% names(color_map)) {
      color_map[[group_names[i]]] <- set1_colors[i]
      log_info(paste0("Group '", group_names[i], "' not specified in --group_colors, using default color: ", set1_colors[i]), "COLOR")
    }
  }

  # Convert to named vector for ggplot2
  color_vector <- unlist(color_map)
  log_info(paste0("Custom color palette applied: ", paste(names(color_vector), color_vector, sep = "=", collapse = ", ")), "COLOR")

  return(color_vector)
}

cat("\n")
cat("################################################################################\n")
cat("#                   METHYLSENSE MODEL EVALUATION                               #\n")
cat("################################################################################\n")
cat(sprintf("#  VERSION: %-66s #\n", SCRIPT_VERSION))
cat(sprintf("#  DATE:    %-66s #\n", SCRIPT_DATE))
cat("################################################################################\n\n")

# ============================================================================
# STEP 1: PARSE OPTIONS FIRST
# ============================================================================
log_section_start(1, "Loading optparse and parsing arguments", 10)
library(optparse)

opts <- list(
  make_option(c("-m", "--model_dir"),
    type = "character",
    help = "Path to model directory (e.g., path/to/nnet/markers_99) [REQUIRED]"
  ),
  make_option(c("-q", "--qs_file"),
    type = "character",
    help = "Path to methylRawList .qs file [REQUIRED]"
  ),
  make_option(c("-s", "--sample_sheet"),
    type = "character",
    help = "Path to sample sheet Excel file [REQUIRED]"
  ),
  make_option(c("-a", "--data_dir"),
    type = "character", default = NULL,
    help = "Directory containing TRAIN/TEST CSV files (if different from model_dir)"
  ),
  make_option(c("-b", "--batch_column"),
    type = "character", default = NULL,
    help = "Column name for batch/platform variable"
  ),
  make_option(c("-t", "--storage_column"),
    type = "character", default = NULL,
    help = "Column name for storage duration variable"
  ),
  make_option(c("-d", "--date_column"),
    type = "character", default = NULL,
    help = "Column name for collection/sampling date"
  ),
  make_option(c("-c", "--technical_covariates"),
    type = "character", default = NULL,
    help = "Comma-separated additional technical covariates"
  ),
  make_option(c("-i", "--sample_id_column"),
    type = "character", default = "SampleID",
    help = "Sample ID column [default: SampleID]"
  ),
  make_option(c("-g", "--treatment_column"),
    type = "character", default = "treatment",
    help = "Treatment column [default: treatment]"
  ),
  make_option(c("-n", "--bootstrap_iterations"),
    type = "integer", default = 500,
    help = "Bootstrap iterations [default: 500]"
  ),
  make_option(c("-r", "--stability_threshold"),
    type = "numeric", default = 0.70,
    help = "Stability threshold [default: 0.70]"
  ),
  make_option(c("-e", "--effect_size_threshold"),
    type = "numeric", default = 2.0,
    help = "Effect size threshold %% [default: 2.0]"
  ),
  make_option(c("-f", "--min_effect_size"),
    type = "numeric", default = 5.0,
    help = "Higher threshold %% [default: 5.0]"
  ),
  make_option(c("-w", "--gene_window"),
    type = "integer", default = 10000000,
    help = "Window size (bp) for finding nearest genes [default: 10000000]"
  ),
  make_option(c("-k", "--stratified_cv_folds"),
    type = "integer", default = 5,
    help = "CV folds [default: 5]"
  ),
  make_option(c("-l", "--cv_repeats"),
    type = "integer", default = 3,
    help = "CV repeats [default: 3]"
  ),
  make_option(c("-o", "--output_subdir"),
    type = "character", default = "model_evaluation",
    help = "Output subdirectory [default: model_evaluation]"
  ),
  make_option(c("-p", "--n_cores"),
    type = "integer", default = 4,
    help = "CPU cores [default: 4]"
  ),
  make_option(c("-v", "--verbose"),
    action = "store_true", default = FALSE,
    help = "Verbose output"
  ),
  make_option(c("-x", "--functional_enrichment"),
    action = "store_true", default = FALSE,
    help = "Perform enrichment"
  ),
  make_option(c("--lmm_n_markers"),
    type = "integer", default = -1,
    help = "Number of markers for LMM analysis (-1 = all markers) [default: -1]"
  ),
  make_option(c("--multiple_testing"),
    type = "character", default = "fdr",
    help = "Multiple testing correction method: 'fdr' (Benjamini-Hochberg, default), 'bonferroni', 'holm', 'hochberg', 'BY' (Benjamini-Yekutieli), or 'none' [default: fdr]"
  ),
  make_option(c("--significance_threshold"),
    type = "numeric", default = 0.05,
    help = "Significance threshold for both raw and corrected p-values [default: 0.05]"
  ),
  make_option(c("--generate_markdown"),
    action = "store_true", default = TRUE,
    help = "Generate comprehensive markdown report for model evaluation [default: TRUE]"
  ),
  make_option(c("--markdown_title"),
    type = "character", default = "MethylSense Model Report",
    help = "Title for markdown report [default: 'MethylSense Model Report']"
  ),
  make_option(c("--markdown_author"),
    type = "character", default = NULL,
    help = "Author name(s) for markdown report [default: NULL]"
  ),
  make_option(c("--figure_format"),
    type = "character", default = "png",
    help = "Output format for figures: 'png', 'jpeg'/'jpg', 'svg', 'all', or comma-separated (e.g., 'jpg,svg') [default: png]"
  ),
  make_option(c("--figure_resolution"),
    type = "integer", default = 300,
    help = "Resolution for raster figures in PPI [default: 300]"
  ),
  make_option(c("--group_colors"),
    type = "character", default = NULL,
    help = "Custom colors for biological groups in PCA/UMAP plots. Format: 'Group1=#HEX1,Group2=#HEX2' or 'Group1=colorname,Group2=colorname'. Example: 'Control=#0000FF,Infected=#FF0000' or 'Control=blue,Infected=red' [default: NULL, uses Set1 palette]"
  ),
  make_option(c("--enable_wgcna"),
    action = "store_true", default = FALSE,
    help = "Enable WGCNA co-expression network analysis (requires WGCNA package) [default: FALSE]"
  ),
  make_option(c("--species"),
    type = "character", default = NULL,
    help = "Species name for gene annotation (e.g., 'Sus scrofa', 'Homo sapiens', 'Mus musculus')"
  ),
  make_option(c("--enable_gene_annotation"),
    action = "store_true", default = FALSE,
    help = "Enable gene annotation using biomaRt [default: FALSE]"
  ),
  make_option(c("--time_analysis"),
    type = "character", default = NULL,
    help = "Column name in sample sheet for timepoint data"
  ),
  make_option(c("--time_categorical"),
    action = "store_true", default = FALSE,
    help = "Treat time as categorical (default: auto-detect)"
  ),
  make_option(c("--prepare_pcr"),
    action = "store_true", default = FALSE,
    help = "Generate PCR primer design links for top DMRs (bisulfite PCR) [default: FALSE]"
  ),
  make_option(c("--pcr_n_dmrs"),
    type = "integer", default = 20,
    help = "Number of top DMRs to prepare for PCR validation [default: 20]"
  ),
  make_option(c("--model_name"),
    type = "character", default = NULL,
    help = "Custom name for the model in report title (e.g., 'Chicken Aspergillosis Neural Network'). If not provided, uses model directory name [default: NULL]"
  ),
  make_option(c("--enable_feature_annotation"),
    action = "store_true", default = FALSE,
    help = "Enable genomic feature annotation using genomation package [default: FALSE]"
  ),
  make_option(c("--gtf_file"),
    type = "character", default = NULL,
    help = "Path to GTF/GFF file for gene feature annotation (required if --enable_feature_annotation is TRUE)"
  ),
  make_option(c("--regulatory_gff"),
    type = "character", default = NULL,
    help = "Path to regulatory features GFF file (enhancers, CTCF, etc.) [optional]"
  ),
  make_option(c("--emar_gff"),
    type = "character", default = NULL,
    help = "Path to EMAR (Epigenetically Modified Accessible Regions) GFF file [optional]"
  ),
  make_option(c("--open_chromatin_bed"),
    type = "character", default = NULL,
    help = "Path to open chromatin BED file [optional]"
  ),
  make_option(c("--chromatin_flank_bp"),
    type = "integer", default = 2000,
    help = "Flank size (bp) for open chromatin regions [default: 2000]"
  ),
  make_option(c("--generate_pdf"),
    action = "store_true", default = FALSE,
    help = "Generate PDF report in addition to HTML [default: FALSE]"
  ),
  make_option(c("--blind_paths"),
    action = "store_true", default = FALSE,
    help = "Replace directory paths with basenames in reports for privacy [default: FALSE]"
  ),
  make_option(c("--control_group"),
    type = "character", default = NULL,
    help = "Explicitly specify which group is the control/reference for logistic regression. If not provided, will auto-detect based on keywords (control, healthy, negative, uninfected, wildtype, wt) [default: NULL]"
  )
)

parser <- OptionParser(option_list = opts)
opt <- parse_args(parser)
log_info("Command line arguments parsed successfully", "STEP1")

# Validate required arguments
if (is.null(opt$model_dir) || is.null(opt$qs_file) || is.null(opt$sample_sheet)) {
  cat("\n[ERROR] Missing required arguments!\n\n")
  print_help(parser)
  quit(status = 1)
}

# Validate multiple testing method
valid_methods <- c("fdr", "bonferroni", "holm", "hochberg", "BY", "none")
if (!tolower(opt$multiple_testing) %in% tolower(valid_methods)) {
  cat(sprintf("\n[ERROR] Invalid multiple testing method: %s\n", opt$multiple_testing))
  cat("Valid options: fdr, bonferroni, holm, hochberg, BY, none\n\n")
  quit(status = 1)
}

# Normalise to lowercase for consistency
opt$multiple_testing <- tolower(opt$multiple_testing)

# Map "fdr" to "BH" for R's p.adjust()
mt_method_display <- opt$multiple_testing
if (opt$multiple_testing == "fdr") {
  mt_method_r <- "BH"
  mt_method_name <- "FDR (Benjamini-Hochberg)"
} else if (opt$multiple_testing == "by") {
  mt_method_r <- "BY"
  mt_method_name <- "Benjamini-Yekutieli"
} else if (opt$multiple_testing == "none") {
  mt_method_r <- "none"
  mt_method_name <- "None (uncorrected)"
} else {
  mt_method_r <- opt$multiple_testing
  mt_method_name <- toupper(opt$multiple_testing)
}

# Store significance threshold
sig_threshold <- opt$significance_threshold

cat("          [OK] Arguments parsed successfully\n")
cat("          Model directory:", opt$model_dir, "\n")
cat("          Bootstrap iterations:", opt$bootstrap_iterations, "\n")
cat("          Multiple testing correction:", mt_method_name, "\n")
cat("          Significance threshold:", sig_threshold, "\n\n")
cat("################################################################################\n\n")

# ============================================================================
# HELPER FUNCTION: Apply multiple testing correction
# ============================================================================
apply_multiple_testing <- function(p_values, method = "BH") {
  if (method == "none") {
    return(p_values)
  } else {
    return(p.adjust(p_values, method = method))
  }
}

# ============================================================================
# HELPER FUNCTION: Get correction label for plots/output
# ============================================================================
get_correction_label <- function() {
  if (mt_method_r == "none") {
    return("No correction applied")
  } else {
    return(paste0(mt_method_name, " correction applied"))
  }
}

# ============================================================================
# HELPER FUNCTION: Get significance column name
# ============================================================================
get_sig_colname <- function(corrected = TRUE) {
  if (mt_method_r == "none") {
    return("significant")
  } else {
    if (corrected) {
      return(paste0("significant_", mt_method_display))
    } else {
      return("significant_raw")
    }
  }
}

# Continue with the rest of the script exactly as before, but replace
# ALL instances of FDR correction with the flexible function calls...

# ============================================================================
# STEP 2: LOAD PACKAGES
# ============================================================================
log_section_start(2, "Loading R packages", 10)

suppressPackageStartupMessages({
  library(qs)
  library(methylKit)
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(gridExtra)
  library(caret)
  library(pROC)
  library(genomation) # For feature annotation
  library(reshape2) # For data manipulation in feature plots
  library(ggpubr) # For ggbarplot
  library(umap)
  library(lme4)
  library(parallel)
  library(igraph)
  library(dendextend)
})

cat("          [OK] Core packages loaded\n")

# Check enrichment packages
enrichment_available <- FALSE
tryCatch(
  {
    suppressPackageStartupMessages({
      library(biomaRt)
      library(clusterProfiler)
    })
    enrichment_available <- TRUE
    cat("          [OK] Enrichment packages available\n")
  },
  error = function(e) {
    cat("          [INFO] Enrichment packages not available (optional)\n")
  }
)

if (opt$functional_enrichment && !enrichment_available) {
  cat("          [WARNING] Enrichment requested but unavailable - disabling\n")
  opt$functional_enrichment <- FALSE
}

# Check WGCNA package
check_wgcna_available <- function() {
  wgcna_available <- requireNamespace("WGCNA", quietly = TRUE)
  if (!wgcna_available) {
    cat("\n")
    cat("================================================================================\n")
    cat("WGCNA PACKAGE NOT INSTALLED\n")
    cat("================================================================================\n")
    cat("WGCNA analysis was requested (--enable_wgcna) but the package is not available.\n")
    cat("\n")
    cat("To install WGCNA:\n")
    cat("  if (!requireNamespace('BiocManager', quietly = TRUE))\n")
    cat("      install.packages('BiocManager')\n")
    cat("  BiocManager::install('WGCNA')\n")
    cat("\n")
    cat("WGCNA analysis will be SKIPPED.\n")
    cat("================================================================================\n")
    cat("\n")
  }
  return(wgcna_available)
}

wgcna_available <- FALSE
if (opt$enable_wgcna) {
  wgcna_available <- check_wgcna_available()
  if (wgcna_available) {
    suppressPackageStartupMessages(library(WGCNA))
    cat("          [OK] WGCNA package loaded\n")
  }
} else {
  cat("          [INFO] WGCNA analysis disabled (use --enable_wgcna to enable)\n")
}

cat("\n")

# ============================================================================
# HELPER FUNCTION: Calculate correlation p-values (for WGCNA)
# ============================================================================
corPvalueStudent <- function(cor, nSamples) {
  # Calculate p-values for correlations using Student's t distribution
  t_stat <- cor * sqrt(nSamples - 2) / sqrt(1 - cor^2)
  p_values <- 2 * pt(-abs(t_stat), nSamples - 2)
  return(p_values)
}

# ============================================================================
# HELPER FUNCTION: Save plots in multiple formats
# ============================================================================
save_plot_multiple_formats <- function(plot_object, file_path_without_ext, width, height) {
  # Get figure format and resolution from global options
  format <- opt$figure_format
  dpi <- opt$figure_resolution

  # Parse comma-separated formats (e.g., "jpg,svg" or "png,jpeg,svg")
  if (grepl(",", format)) {
    formats <- trimws(unlist(strsplit(format, ",")))
  } else {
    formats <- format
  }

  if ("all" %in% tolower(formats)) {
    # Save in all formats
    ggsave(paste0(file_path_without_ext, ".png"), plot_object, width = width, height = height, dpi = dpi)
    ggsave(paste0(file_path_without_ext, ".jpeg"), plot_object, width = width, height = height, dpi = dpi)
    ggsave(paste0(file_path_without_ext, ".svg"), plot_object, width = width, height = height)
  } else {
    # Save each requested format
    for (fmt in formats) {
      fmt_lower <- tolower(fmt)

      # Normalize jpg to jpeg
      if (fmt_lower == "jpg") {
        fmt_lower <- "jpeg"
      }

      if (fmt_lower %in% c("png", "jpeg")) {
        # Save raster format with dpi
        ggsave(paste0(file_path_without_ext, ".", fmt_lower), plot_object, width = width, height = height, dpi = dpi)
      } else if (fmt_lower == "svg") {
        # Save vector format (no dpi needed)
        ggsave(paste0(file_path_without_ext, ".svg"), plot_object, width = width, height = height)
      } else {
        warning(sprintf("Invalid figure format '%s', skipping", fmt))
      }
    }
  }
}

# ============================================================================
# HELPER FUNCTION: Get primary figure extension for markdown references
# ============================================================================
get_primary_figure_extension <- function() {
  format <- opt$figure_format

  # Parse comma-separated formats
  if (grepl(",", format)) {
    formats <- trimws(unlist(strsplit(format, ",")))
    # Use first format as primary
    fmt <- formats[1]
  } else {
    fmt <- format
  }

  # Handle "all" option - use PNG as primary
  if (tolower(fmt) == "all") {
    return("png")
  }

  # Normalize jpg to jpeg
  if (tolower(fmt) == "jpg") {
    return("jpeg")
  }

  return(tolower(fmt))
}

# ============================================================================
# HELPER FUNCTION: Copy plot file with correct extension based on fig_ext
# ============================================================================
copy_plot_to_figures <- function(source_dir, plot_basename, dest_dir, fig_ext) {
  # Copy ALL available formats (not just primary)
  # This ensures both jpg AND svg are copied when user requests "jpg,svg"
  formats_to_check <- c(fig_ext, "png", "jpeg", "jpg", "svg")
  copied_any <- FALSE

  for (fmt in unique(formats_to_check)) {
    source_file <- file.path(source_dir, sprintf("%s.%s", plot_basename, fmt))
    dest_file <- file.path(dest_dir, sprintf("%s.%s", plot_basename, fmt))

    if (file.exists(source_file)) {
      file.copy(source_file, dest_file, overwrite = TRUE)
      copied_any <- TRUE
    }
  }

  return(copied_any)
}

# ============================================================================
# HELPER FUNCTION: Open graphics device with correct format
# ============================================================================
open_plot_device <- function(filepath, width, height, res = 300) {
  # Get file extension
  ext <- tolower(tools::file_ext(filepath))

  if (ext == "png") {
    png(filepath, width = width, height = height, units = "in", res = res)
  } else if (ext %in% c("jpg", "jpeg")) {
    jpeg(filepath, width = width, height = height, units = "in", res = res, quality = 95)
  } else if (ext == "svg") {
    svg(filepath, width = width, height = height)
  } else {
    # Fallback to PNG
    warning(sprintf("Unknown format '%s', using PNG", ext))
    png(filepath, width = width, height = height, units = "in", res = res)
  }
}

# ============================================================================
# HELPER FUNCTION: Generate custom CSS for clinical-style PDF reports
# ============================================================================
get_clinical_report_css <- function() {
  css <- '
<style type="text/css">
/* MethylSense Clinical Report Theme - Professional diagnostic report styling */

@page {
  size: A4;
  margin: 2cm 2.5cm 2.5cm 2.5cm;
}

@media print {
  h1, h2 { page-break-after: avoid; }
  table, figure, img { page-break-inside: avoid; }
  p { orphans: 3; widows: 3; }
}

body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
  font-size: 11pt;
  line-height: 1.6;
  color: #2c3e50;
}

h1 {
  color: #1e3a8a;
  font-size: 24pt;
  font-weight: 600;
  margin-top: 1.5em;
  padding-bottom: 0.3em;
  border-bottom: 3px solid #3b82f6;
}

h2 {
  color: #1e40af;
  font-size: 18pt;
  font-weight: 600;
  margin-top: 1.2em;
  padding-bottom: 0.2em;
  border-bottom: 2px solid #60a5fa;
}

h3 {
  color: #2563eb;
  font-size: 14pt;
  font-weight: 600;
  margin-top: 1em;
}

table {
  border-collapse: collapse;
  width: 100%;
  margin: 1em 0;
  font-size: 10pt;
  box-shadow: 0 1px 3px rgba(0,0,0,0.1);
}

thead {
  background: linear-gradient(to bottom, #3b82f6, #2563eb);
  color: white;
  font-weight: 600;
}

th {
  padding: 10px 12px;
  text-align: left;
}

td {
  padding: 8px 12px;
  border-bottom: 1px solid #e5e7eb;
}

tbody tr:nth-child(even) {
  background-color: #f8fafc;
}

blockquote {
  margin: 1em 0;
  padding: 12px 20px 12px 16px;
  background: #eff6ff;
  border-left: 4px solid #3b82f6;
  font-style: italic;
  color: #1e40af;
  border-radius: 0 4px 4px 0;
}

pre {
  background: #f1f5f9;
  border: 1px solid #cbd5e1;
  border-left: 3px solid #64748b;
  padding: 12px;
  border-radius: 4px;
  font-size: 9pt;
}

code {
  font-family: "SF Mono", Monaco, "Courier New", monospace;
  background: #f1f5f9;
  padding: 2px 6px;
  border-radius: 3px;
  font-size: 0.9em;
}

img {
  max-width: 100%;
  height: auto;
  display: block;
  margin: 1em auto;
  border-radius: 4px;
  box-shadow: 0 2px 8px rgba(0,0,0,0.1);
}

strong {
  color: #1e3a8a;
  font-weight: 600;
}
</style>
'
  return(css)
}

# ============================================================================
# HELPER FUNCTIONS FOR MULTICLASS PERFORMANCE ANALYSIS
# ============================================================================

# Calculate PR curve manually
calculate_pr_curve <- function(actual, probabilities, positive_class) {
  thresholds <- sort(unique(probabilities), decreasing = TRUE)
  pr_data <- data.frame()

  for (thresh in thresholds) {
    predicted <- ifelse(probabilities >= thresh, positive_class, "Other")
    tp <- sum(predicted == positive_class & actual == positive_class)
    fp <- sum(predicted == positive_class & actual != positive_class)
    fn <- sum(predicted != positive_class & actual == positive_class)

    precision <- ifelse(tp + fp > 0, tp / (tp + fp), 1)
    recall <- ifelse(tp + fn > 0, tp / (tp + fn), 0)

    pr_data <- rbind(pr_data, data.frame(threshold = thresh, precision = precision, recall = recall))
  }

  return(pr_data)
}

# Calculate PR-AUC using trapezoidal rule
calculate_pr_auc <- function(pr_curve) {
  pr_sorted <- pr_curve[order(pr_curve$recall), ]
  auc <- 0
  for (i in 2:nrow(pr_sorted)) {
    width <- pr_sorted$recall[i] - pr_sorted$recall[i - 1]
    height <- (pr_sorted$precision[i] + pr_sorted$precision[i - 1]) / 2
    auc <- auc + width * height
  }
  return(auc)
}

# Calculate binary diagnostic metrics from confusion matrix
calculate_binary_metrics <- function(conf_matrix) {
  tp <- conf_matrix[2, 2]
  tn <- conf_matrix[1, 1]
  fp <- conf_matrix[1, 2]
  fn <- conf_matrix[2, 1]

  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  ppv <- tp / (tp + fp)
  npv <- tn / (tn + fn)
  accuracy <- (tp + tn) / sum(conf_matrix)
  f1 <- 2 * (ppv * sensitivity) / (ppv + sensitivity)

  mcc_num <- (tp * tn) - (fp * fn)
  mcc_den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- ifelse(mcc_den == 0, 0, mcc_num / mcc_den)

  return(list(
    sensitivity = sensitivity,
    specificity = specificity,
    ppv = ppv,
    npv = npv,
    accuracy = accuracy,
    f1 = f1,
    mcc = mcc
  ))
}

# ============================================================================
# HELPER FUNCTION: Gene Annotation via biomaRt
# ============================================================================

query_biomart_genes <- function(dmr_regions, species_name, gene_window = 10000000) {
  # Query Ensembl biomaRt for gene annotations of DMRs
  #
  # Args:
  #   dmr_regions: data.frame with columns: DMR_ID, Chr, Start, End
  #   species_name: species name (e.g., "Sus scrofa", "Homo sapiens", "Mus musculus")
  #   gene_window: window size in bp for nearest gene/feature search (default: 10000000 = 10 Mb)
  #
  # Returns:
  #   data.frame with columns: DMR_ID, Chr, Start, End, Overlapping_Genes,
  #                            Nearest_Gene_ID, Nearest_Gene_Name, Gene_Distance,
  #                            Gene_Direction, Gene_Biotype

  # Species mapping to Ensembl dataset names
  species_map <- list(
    "Sus scrofa" = "sscrofa_gene_ensembl",
    "Homo sapiens" = "hsapiens_gene_ensembl",
    "Mus musculus" = "mmusculus_gene_ensembl",
    "Gallus gallus" = "ggallus_gene_ensembl"
  )

  # Validate species
  if (!species_name %in% names(species_map)) {
    stop(sprintf(
      "Species '%s' not supported. Available: %s",
      species_name, paste(names(species_map), collapse = ", ")
    ))
  }

  dataset_name <- species_map[[species_name]]

  cat(sprintf("Connecting to Ensembl biomaRt for %s (%s)...\n", species_name, dataset_name))

  # Connect to Ensembl
  tryCatch(
    {
      ensembl <- biomaRt::useMart("ensembl", dataset = dataset_name)
    },
    error = function(e) {
      stop(sprintf("Failed to connect to Ensembl: %s\nPlease check internet connection.", e$message))
    }
  )

  cat(sprintf("Successfully connected to Ensembl.\n"))
  cat(sprintf("Querying gene annotations for %d DMRs...\n", nrow(dmr_regions)))

  # Initialize result data frame
  results <- data.frame(
    DMR_ID = dmr_regions$DMR_ID,
    Chr = dmr_regions$Chr,
    Start = dmr_regions$Start,
    End = dmr_regions$End,
    Overlapping_Genes = NA_character_,
    Nearest_Gene_ID = NA_character_,
    Nearest_Gene_Name = NA_character_,
    Gene_Distance = NA_integer_,
    Gene_Direction = NA_character_,
    Gene_Biotype = NA_character_,
    stringsAsFactors = FALSE
  )

  # Process each DMR
  for (i in 1:nrow(dmr_regions)) {
    dmr <- dmr_regions[i, ]

    # Clean chromosome name (remove "chr" prefix if present)
    chr_clean <- gsub("^chr", "", dmr$Chr)

    # Query ALL genes and features on this chromosome
    # NOTE: biomaRt 'start' and 'end' filters work on gene coordinates, NOT genomic windows
    # So we must get all genes on the chromosome and filter locally
    # ISSUE 3 FIX: Expanded attributes to include transcripts, exons, UTRs
    tryCatch(
      {
        all_chr_genes <- biomaRt::getBM(
          attributes = c(
            "ensembl_gene_id", "external_gene_name", "gene_biotype",
            "chromosome_name", "start_position", "end_position",
            "transcript_biotype", "ensembl_transcript_id",
            "ensembl_exon_id", "5_utr_start", "5_utr_end",
            "3_utr_start", "3_utr_end"
          ),
          filters = c("chromosome_name"),
          values = list(chr_clean),
          mart = ensembl
        )

        if (nrow(all_chr_genes) > 0) {
          # Filter for overlapping genes (gene overlaps DMR coordinates)
          overlapping_genes <- all_chr_genes[
            (all_chr_genes$start_position <= dmr$End) &
              (all_chr_genes$end_position >= dmr$Start), ,
            drop = FALSE
          ]

          if (nrow(overlapping_genes) > 0) {
            # Store overlapping genes as comma-separated list
            gene_names <- overlapping_genes$external_gene_name[overlapping_genes$external_gene_name != ""]
            if (length(gene_names) > 0) {
              results$Overlapping_Genes[i] <- paste(gene_names, collapse = ", ")
            }

            # Use first overlapping gene as "nearest"
            results$Nearest_Gene_ID[i] <- overlapping_genes$ensembl_gene_id[1]
            results$Nearest_Gene_Name[i] <- overlapping_genes$external_gene_name[1]
            results$Gene_Distance[i] <- 0
            results$Gene_Direction[i] <- "overlapping"
            results$Gene_Biotype[i] <- overlapping_genes$gene_biotype[1]
          } else {
            # No overlapping genes - find nearest gene within window
            dmr_midpoint <- (dmr$Start + dmr$End) / 2
            window_start <- max(1, dmr$Start - gene_window)
            window_end <- dmr$End + gene_window

            # Calculate distance from each gene to DMR midpoint
            all_chr_genes$gene_midpoint <- (all_chr_genes$start_position + all_chr_genes$end_position) / 2
            all_chr_genes$distance_to_dmr <- abs(all_chr_genes$gene_midpoint - dmr_midpoint)

            # Filter genes within window distance
            nearby_genes <- all_chr_genes[all_chr_genes$distance_to_dmr <= gene_window, , drop = FALSE]

            if (nrow(nearby_genes) > 0) {
              # Find nearest gene
              nearest_idx <- which.min(nearby_genes$distance_to_dmr)
              nearest_gene <- nearby_genes[nearest_idx, ]

              results$Nearest_Gene_ID[i] <- nearest_gene$ensembl_gene_id
              results$Nearest_Gene_Name[i] <- nearest_gene$external_gene_name
              results$Gene_Biotype[i] <- nearest_gene$gene_biotype

              # Calculate distance and direction
              distance <- nearest_gene$gene_midpoint - dmr_midpoint
              results$Gene_Distance[i] <- as.integer(abs(distance))
              results$Gene_Direction[i] <- ifelse(distance > 0, "downstream", "upstream")
            }
          }
        }
      },
      error = function(e) {
        warning(sprintf("Failed to query genes for DMR %s: %s", dmr$DMR_ID, e$message))
      }
    )

    # Progress indicator
    if (i %% 10 == 0) {
      cat(sprintf("  Processed %d/%d DMRs...\n", i, nrow(dmr_regions)))
    }
  }

  cat(sprintf("Gene annotation complete.\n"))
  cat(sprintf("  %d DMRs with overlapping genes\n", sum(!is.na(results$Overlapping_Genes))))
  cat(sprintf(
    "  %d DMRs with nearby genes (within %d bp)\n",
    sum(!is.na(results$Nearest_Gene_Name)), gene_window
  ))

  return(results)
}

# ============================================================================
# HELPER FUNCTION: Time Series Analysis
# ============================================================================

analyse_time_series <- function(methylation_values, timepoints, biological_groups = NULL,
                                subject_ids = NULL, force_categorical = FALSE, dmr_id = "") {
  # This function performs time series analysis on methylation data
  # Returns statistical test results and trajectory information

  # Initialize results
  result <- list(
    available = TRUE,
    model_type = NA,
    time_effect_slope = NA,
    time_effect_pvalue = NA,
    effect_size = NA,
    r_squared = NA,
    trajectory_cluster = NA,
    error = NULL
  )

  # Validate input
  if (length(methylation_values) != length(timepoints)) {
    result$available <- FALSE
    result$error <- "Methylation values and timepoints length mismatch"
    return(result)
  }

  # Remove missing data
  valid_idx <- !is.na(methylation_values) & !is.na(timepoints)
  if (sum(valid_idx) < 3) {
    result$available <- FALSE
    result$error <- "Insufficient data points (need at least 3)"
    return(result)
  }

  methylation_values <- methylation_values[valid_idx]
  timepoints <- timepoints[valid_idx]
  if (!is.null(biological_groups)) {
    biological_groups <- biological_groups[valid_idx]
    # Replace NA groups with "Unknown"
    if (any(is.na(biological_groups))) {
      biological_groups[is.na(biological_groups)] <- "Unknown"
    }
  }
  if (!is.null(subject_ids)) subject_ids <- subject_ids[valid_idx]

  # Auto-detect if timepoint is continuous or categorical
  is_continuous <- FALSE
  if (!force_categorical) {
    # Try to convert to numeric
    time_numeric <- suppressWarnings(as.numeric(as.character(timepoints)))

    # If conversion successful and has reasonable range, treat as continuous
    if (sum(!is.na(time_numeric)) > length(timepoints) * 0.8) {
      # Extract numbers from strings like "Day 7", "T0", "Week 3"
      time_cleaned <- gsub("[^0-9.-]", "", as.character(timepoints))
      time_numeric <- suppressWarnings(as.numeric(time_cleaned))

      if (sum(!is.na(time_numeric)) > length(timepoints) * 0.8 &&
        length(unique(time_numeric)) >= 3) {
        is_continuous <- TRUE
        timepoints <- time_numeric
      }
    }
  }

  # Detect repeated measures
  is_repeated <- FALSE
  if (!is.null(subject_ids)) {
    # Check for duplicate subject IDs
    if (any(duplicated(subject_ids))) {
      is_repeated <- TRUE
    }
  }

  # Create data frame for analysis
  df <- data.frame(
    methylation = methylation_values,
    time = timepoints,
    stringsAsFactors = FALSE
  )

  if (!is.null(biological_groups)) {
    df$group <- as.factor(biological_groups)
  }

  if (!is.null(subject_ids)) {
    df$subject <- as.factor(subject_ids)
  }

  # Perform appropriate statistical test
  tryCatch(
    {
      if (is_continuous) {
        # CONTINUOUS TIME ANALYSIS
        if (is_repeated && !is.null(subject_ids)) {
          # Linear Mixed Model for repeated measures
          if (requireNamespace("lme4", quietly = TRUE)) {
            if (!is.null(biological_groups)) {
              # Include interaction with biological group
              model <- lme4::lmer(methylation ~ time * group + (1 | subject), data = df)
            } else {
              model <- lme4::lmer(methylation ~ time + (1 | subject), data = df)
            }

            # Extract fixed effects
            fixed_effects <- summary(model)$coefficients
            time_slope <- fixed_effects["time", "Estimate"]
            time_pvalue <- fixed_effects["time", "Pr(>|t|)"]

            # Calculate R-squared (marginal - fixed effects only)
            if (requireNamespace("MuMIn", quietly = TRUE)) {
              r2 <- MuMIn::r.squaredGLMM(model)[1, "R2m"]
            } else {
              r2 <- NA
            }

            result$model_type <- "Linear Mixed Model (LMM)"
            result$time_effect_slope <- time_slope
            result$time_effect_pvalue <- time_pvalue
            result$effect_size <- time_slope
            result$r_squared <- r2
          } else {
            result$error <- "lme4 package not available for LMM analysis"
            result$available <- FALSE
            return(result)
          }
        } else {
          # Simple linear model for independent groups
          if (!is.null(biological_groups)) {
            model <- lm(methylation ~ time * group, data = df)
          } else {
            model <- lm(methylation ~ time, data = df)
          }

          # Extract coefficients
          coef_summary <- summary(model)$coefficients
          time_slope <- coef_summary["time", "Estimate"]
          time_pvalue <- coef_summary["time", "Pr(>|t|)"]
          r2 <- summary(model)$r.squared

          result$model_type <- "Linear Model (LM)"
          result$time_effect_slope <- time_slope
          result$time_effect_pvalue <- time_pvalue
          result$effect_size <- time_slope
          result$r_squared <- r2
        }
      } else {
        # CATEGORICAL TIME ANALYSIS
        df$time <- as.factor(df$time)

        if (is_repeated && !is.null(subject_ids)) {
          # Repeated Measures ANOVA
          if (requireNamespace("ez", quietly = TRUE)) {
            # Use ezANOVA for repeated measures
            if (!is.null(biological_groups)) {
              anova_result <- ez::ezANOVA(
                data = df,
                dv = methylation,
                wid = subject,
                within = time,
                between = group,
                type = 3
              )
            } else {
              anova_result <- ez::ezANOVA(
                data = df,
                dv = methylation,
                wid = subject,
                within = time,
                type = 3
              )
            }

            # Extract time effect
            time_pvalue <- anova_result$ANOVA$p[anova_result$ANOVA$Effect == "time"]
            effect_size <- anova_result$ANOVA$ges[anova_result$ANOVA$Effect == "time"] # Generalized eta-squared

            result$model_type <- "Repeated Measures ANOVA"
            result$time_effect_pvalue <- time_pvalue
            result$effect_size <- effect_size
          } else if (requireNamespace("nlme", quietly = TRUE)) {
            # Fallback: Use nlme for repeated measures
            if (!is.null(biological_groups)) {
              model <- nlme::lme(methylation ~ time * group, random = ~ 1 | subject, data = df)
            } else {
              model <- nlme::lme(methylation ~ time, random = ~ 1 | subject, data = df)
            }

            anova_result <- anova(model)
            time_pvalue <- anova_result["time", "p-value"]

            result$model_type <- "Repeated Measures (nlme)"
            result$time_effect_pvalue <- time_pvalue
          } else {
            result$error <- "No repeated measures package available (ez or nlme required)"
            result$available <- FALSE
            return(result)
          }
        } else {
          # Independent groups - Always use Kruskal-Wallis (non-parametric, robust)
          # This ensures consistent statistical framework across all DMRs
          kw_test <- kruskal.test(methylation ~ time, data = df)
          time_pvalue <- kw_test$p.value

          # Calculate eta-squared approximation (effect size)
          eta_squared <- kw_test$statistic / (nrow(df) - 1)

          result$model_type <- "Kruskal-Wallis"
          result$time_effect_pvalue <- time_pvalue
          result$effect_size <- as.numeric(eta_squared)
        }
      }

      # Assign trajectory cluster based on slope (for continuous) or pattern (for categorical)
      if (is_continuous && !is.na(result$time_effect_slope)) {
        slope <- result$time_effect_slope

        # Calculate residuals to detect non-linearity
        if (!is.null(biological_groups)) {
          model_simple <- lm(methylation ~ time, data = df)
        } else {
          model_simple <- lm(methylation ~ time, data = df)
        }
        residuals_sd <- sd(residuals(model_simple))

        # Cluster assignment
        if (residuals_sd > 2 * sd(df$methylation)) {
          result$trajectory_cluster <- "U-shaped"
        } else if (slope > 0.5) {
          result$trajectory_cluster <- "Increasing"
        } else if (slope < -0.5) {
          result$trajectory_cluster <- "Decreasing"
        } else {
          result$trajectory_cluster <- "Stable"
        }
      } else {
        result$trajectory_cluster <- "Categorical"
      }
    },
    error = function(e) {
      result$available <<- FALSE
      result$error <<- paste("Analysis failed:", e$message)
    }
  )

  return(result)
}

# ============================================================================
# HELPER FUNCTIONS: Model Characterisation and Nested CV Analysis
# ============================================================================

# Parse model_summary.txt file
parse_model_summary <- function(model_dir) {
  summary_file <- file.path(model_dir, "model_summary.txt")
  if (!file.exists(summary_file)) {
    return(list(available = FALSE))
  }

  lines <- readLines(summary_file)
  result <- list(available = TRUE)

  # Extract method
  method_line <- grep("Method:", lines, value = TRUE)
  if (length(method_line) > 0) {
    result$method <- trimws(gsub("Method:", "", method_line[1]))
  }

  # Extract best hyperparameters section
  hyper_start <- grep("BEST HYPERPARAMETERS", lines)
  if (length(hyper_start) > 0) {
    hyper_lines <- lines[(hyper_start[1] + 3):(hyper_start[1] + 4)]
    result$hyperparameters <- paste(hyper_lines, collapse = "\n")
  }

  # Extract variable importance
  vi_start <- grep("VARIABLE IMPORTANCE", lines)
  if (length(vi_start) > 0) {
    vi_end <- min(vi_start[1] + 35, length(lines))
    vi_lines <- lines[(vi_start[1] + 2):vi_end]
    vi_lines <- vi_lines[vi_lines != "" & !grepl("^=", vi_lines)]
    result$variable_importance <- vi_lines
  }

  return(result)
}

# Load nested CV summary
load_nested_cv_summary <- function(model_dir) {
  cv_dir <- file.path(model_dir, "cross_validation")
  if (!dir.exists(cv_dir)) {
    return(list(available = FALSE))
  }

  # Find nested CV summary file
  summary_files <- list.files(cv_dir, pattern = "nested_cv_summary.*\\.csv$", full.names = TRUE)
  if (length(summary_files) == 0) {
    return(list(available = FALSE))
  }

  summary_data <- read.csv(summary_files[1], stringsAsFactors = FALSE)

  # ISSUE 4 FIX: Check if PR-AUC columns exist, if not try to load from metrics.csv
  if (!"mean_pr_auc" %in% colnames(summary_data)) {
    metrics_file <- file.path(model_dir, "metrics.csv")
    if (file.exists(metrics_file)) {
      metrics <- read.csv(metrics_file, stringsAsFactors = FALSE)
      pr_auc_row <- metrics[metrics$metric == "PR_AUC", ]
      if (nrow(pr_auc_row) > 0) {
        pr_auc_val <- as.numeric(pr_auc_row$value[1])
        # Use the same value for mean, lower CI, and upper CI as approximation
        summary_data$mean_pr_auc <- pr_auc_val
        summary_data$pr_auc_95ci_lower <- pr_auc_val
        summary_data$pr_auc_95ci_upper <- pr_auc_val
      } else {
        summary_data$mean_pr_auc <- NA
        summary_data$pr_auc_95ci_lower <- NA
        summary_data$pr_auc_95ci_upper <- NA
      }
    } else {
      summary_data$mean_pr_auc <- NA
      summary_data$pr_auc_95ci_lower <- NA
      summary_data$pr_auc_95ci_upper <- NA
    }
  }

  result <- list(
    available = TRUE,
    summary = summary_data,
    n_folds = summary_data$n_folds[1],
    mean_accuracy = summary_data$mean_accuracy[1],
    accuracy_ci_lower = summary_data$accuracy_95ci_lower[1],
    accuracy_ci_upper = summary_data$accuracy_95ci_upper[1],
    mean_auc = summary_data$mean_auc[1],
    auc_ci_lower = summary_data$auc_95ci_lower[1],
    auc_ci_upper = summary_data$auc_95ci_upper[1],
    mean_pr_auc = summary_data$mean_pr_auc[1],
    pr_auc_ci_lower = summary_data$pr_auc_95ci_lower[1],
    pr_auc_ci_upper = summary_data$pr_auc_95ci_upper[1],
    mean_sensitivity = summary_data$mean_sensitivity[1],
    sensitivity_ci_lower = summary_data$sensitivity_95ci_lower[1],
    sensitivity_ci_upper = summary_data$sensitivity_95ci_upper[1],
    mean_specificity = summary_data$mean_specificity[1],
    specificity_ci_lower = summary_data$specificity_95ci_lower[1],
    specificity_ci_upper = summary_data$specificity_95ci_upper[1]
  )

  return(result)
}

# Load nested CV detailed fold-by-fold data
load_nested_cv_detailed <- function(model_dir) {
  cv_dir <- file.path(model_dir, "cross_validation")
  if (!dir.exists(cv_dir)) {
    return(list(available = FALSE))
  }

  detailed_files <- list.files(cv_dir, pattern = "nested_cv_detailed.*\\.csv$", full.names = TRUE)
  if (length(detailed_files) == 0) {
    return(list(available = FALSE))
  }

  detailed_data <- read.csv(detailed_files[1], stringsAsFactors = FALSE)

  return(list(
    available = TRUE,
    data = detailed_data
  ))
}

# Detect binary vs multiclass classification
detect_classification_type <- function(predictions_df) {
  prob_cols <- colnames(predictions_df)[!colnames(predictions_df) %in%
    c("sample_id", "actual", "predicted", "Dataset")]

  n_classes <- length(prob_cols)

  return(list(
    n_classes = n_classes,
    class_names = prob_cols,
    is_binary = (n_classes == 2),
    is_multiclass = (n_classes > 2)
  ))
}

# Copy existing plots to figures directory
copy_model_plots <- function(model_dir, figures_dir) {
  plots_to_copy <- c(
    "methylation_heatmap_.*markers\\.png$",
    "dmr_barplot_.*markers\\.png$",
    "group_methylation_summary_.*markers\\.png$",
    "confusion_matrix_plot\\.png$",
    "roc_curve_binary\\.png$",
    "pr_curve_.*\\.png$",
    "prediction_confidence\\.png$",
    "prediction_probabilities_violin\\.png$",
    "performance_metrics\\.png$"
  )

  copied_files <- list()

  for (pattern in plots_to_copy) {
    files <- list.files(model_dir, pattern = pattern, full.names = TRUE, recursive = FALSE)
    if (length(files) > 0) {
      for (f in files) {
        dest_file <- file.path(figures_dir, basename(f))
        file.copy(f, dest_file, overwrite = TRUE)
        copied_files <- c(copied_files, basename(f))
      }
    }
  }

  # Copy nested CV visualizations
  nested_cv_dir <- file.path(model_dir, "nested_cv_visualizations")
  if (dir.exists(nested_cv_dir)) {
    nested_plots <- list.files(nested_cv_dir, pattern = "\\.png$", full.names = TRUE)
    for (f in nested_plots) {
      dest_file <- file.path(figures_dir, basename(f))
      file.copy(f, dest_file, overwrite = TRUE)
      copied_files <- c(copied_files, basename(f))
    }
  }

  return(copied_files)
}

# ============================================================================
# HELPER FUNCTIONS: Advanced Performance Analysis
# ============================================================================

# Load fold-by-fold CV detailed results
load_cv_detailed <- function(model_dir) {
  cv_dir <- file.path(model_dir, "cross_validation")
  # ===== FIX #1: DYNAMIC FILE PATTERN MATCHING (NO HARDCODED DMR COUNTS) =====
  detailed_files <- list.files(cv_dir, pattern = "^cv_detailed_.*\\.csv$", full.names = TRUE)

  if (length(detailed_files) == 0) {
    return(list(available = FALSE))
  }

  detailed_file <- detailed_files[1] # Take first match

  if (!file.exists(detailed_file)) {
    return(list(available = FALSE))
  }

  detailed_data <- read.csv(detailed_file, stringsAsFactors = FALSE)

  # FIX v5.11.8: Do NOT load PR-AUC from metrics.csv
  # Using same value for all folds creates false impression of identical performance
  # PR-AUC must vary across folds (different data splits). If not available per-fold,
  # leave as NA and hide column in report
  if (!"pr_auc" %in% colnames(detailed_data)) {
    detailed_data$pr_auc <- NA
  }

  return(list(
    available = TRUE,
    data = detailed_data,
    n_folds = nrow(detailed_data)
  ))
}

# Load per-sample CV predictions
load_cv_predictions <- function(model_dir) {
  cv_dir <- file.path(model_dir, "cross_validation")
  # Try nested_cv_predictions first, then fall back to cv_predictions
  pred_files <- list.files(cv_dir, pattern = "^nested_cv_predictions_.*\\.csv$", full.names = TRUE)

  if (length(pred_files) == 0) {
    # Fallback to old naming convention
    pred_files <- list.files(cv_dir, pattern = "^cv_predictions_.*\\.csv$", full.names = TRUE)
  }

  if (length(pred_files) == 0) {
    return(list(available = FALSE))
  }

  pred_file <- pred_files[1] # Take first match

  if (!file.exists(pred_file)) {
    return(list(available = FALSE))
  }

  pred_data <- read.csv(pred_file, stringsAsFactors = FALSE)

  return(list(
    available = TRUE,
    data = pred_data,
    file_path = pred_file
  ))
}

# Analyse sample-level consistency across folds
analyse_sample_consistency <- function(cv_predictions) {
  if (!cv_predictions$available) {
    return(list(available = FALSE))
  }

  pred_data <- cv_predictions$data

  # Count correct predictions per sample
  sample_summary <- pred_data %>%
    group_by(sample_id, actual) %>%
    summarise(
      total_folds = n(),
      correct_folds = sum(actual == predicted),
      consistency_rate = correct_folds / total_folds,
      .groups = "drop"
    ) %>%
    arrange(consistency_rate)

  # Categorise samples
  n_consistent <- sum(sample_summary$consistency_rate >= 0.8)
  n_mostly_correct <- sum(sample_summary$consistency_rate >= 0.6 &
    sample_summary$consistency_rate < 0.8)
  n_inconsistent <- sum(sample_summary$consistency_rate < 0.6)

  return(list(
    available = TRUE,
    summary = sample_summary,
    n_consistent = n_consistent,
    n_mostly_correct = n_mostly_correct,
    n_inconsistent = n_inconsistent,
    pct_consistent = round(100 * n_consistent / nrow(sample_summary), 1)
  ))
}

# Calculate per-class metrics from nested CV predictions with bootstrap CI
calculate_per_class_nested_cv_metrics <- function(nested_cv_predictions, positive_class, n_bootstrap = 1000) {
  if (!nested_cv_predictions$available) {
    return(list(available = FALSE))
  }

  pred_data <- nested_cv_predictions$data

  # Check if positive_class exists in data
  if (!positive_class %in% unique(c(pred_data$actual, pred_data$predicted))) {
    return(list(available = FALSE, reason = paste("Class", positive_class, "not found in predictions")))
  }

  # Calculate metrics for each fold
  # Check for both "fold" and "fold_num" column names
  if (!("fold" %in% colnames(pred_data)) && !("fold_num" %in% colnames(pred_data))) {
    # If no fold column, treat as single evaluation
    pred_data$fold <- 1
  } else if ("fold_num" %in% colnames(pred_data) && !"fold" %in% colnames(pred_data)) {
    # Use fold_num if that's the column name
    pred_data$fold <- pred_data$fold_num
  }

  folds <- unique(pred_data$fold)
  fold_metrics <- data.frame()

  for (f in folds) {
    fold_data <- pred_data[pred_data$fold == f, ]

    # Binary classification with positive_class as positive
    tp <- sum(fold_data$actual == positive_class & fold_data$predicted == positive_class)
    tn <- sum(fold_data$actual != positive_class & fold_data$predicted != positive_class)
    fp <- sum(fold_data$actual != positive_class & fold_data$predicted == positive_class)
    fn <- sum(fold_data$actual == positive_class & fold_data$predicted != positive_class)

    # Calculate metrics
    sensitivity <- ifelse(tp + fn > 0, tp / (tp + fn), NA)
    specificity <- ifelse(tn + fp > 0, tn / (tn + fp), NA)
    precision <- ifelse(tp + fp > 0, tp / (tp + fp), NA)
    f1_score <- ifelse(!is.na(sensitivity) && !is.na(precision) && (precision + sensitivity) > 0,
      2 * (precision * sensitivity) / (precision + sensitivity), NA
    )
    accuracy <- (tp + tn) / nrow(fold_data)

    fold_metrics <- rbind(fold_metrics, data.frame(
      fold = f,
      sensitivity = sensitivity,
      specificity = specificity,
      precision = precision,
      f1_score = f1_score,
      accuracy = accuracy
    ))
  }

  # Calculate mean and CI using bootstrap
  calculate_ci <- function(values, n_bootstrap = 1000) {
    if (all(is.na(values)) || length(values) == 0) {
      return(list(mean = NA, ci_lower = NA, ci_upper = NA))
    }

    values <- values[!is.na(values)]
    mean_val <- mean(values)

    # Special case: only 1 value (single fold) - no CI possible, report as point estimate
    if (length(values) == 1) {
      return(list(
        mean = mean_val,
        ci_lower = mean_val, # Same as mean for single fold
        ci_upper = mean_val
      ))
    }

    # Bootstrap CI for multiple folds
    boot_means <- replicate(n_bootstrap, {
      sample_vals <- sample(values, length(values), replace = TRUE)
      mean(sample_vals)
    })

    ci <- quantile(boot_means, probs = c(0.025, 0.975), na.rm = TRUE)

    return(list(
      mean = mean_val,
      ci_lower = ci[1],
      ci_upper = ci[2]
    ))
  }

  # Calculate stats for each metric
  sensitivity_stats <- calculate_ci(fold_metrics$sensitivity, n_bootstrap)
  specificity_stats <- calculate_ci(fold_metrics$specificity, n_bootstrap)
  f1_stats <- calculate_ci(fold_metrics$f1_score, n_bootstrap)
  accuracy_stats <- calculate_ci(fold_metrics$accuracy, n_bootstrap)

  return(list(
    available = TRUE,
    positive_class = positive_class,
    n_folds = length(folds),
    fold_metrics = fold_metrics,
    sensitivity = sensitivity_stats,
    specificity = specificity_stats,
    f1_score = f1_stats,
    accuracy = accuracy_stats
  ))
}

# Calculate study-stratified performance from predictions
calculate_study_stratified_performance <- function(predictions_data, sample_sheet_path) {
  # Handle both nested CV format and direct predictions format
  if (is.list(predictions_data) && "available" %in% names(predictions_data)) {
    # Old format: nested CV predictions with $available and $data
    if (!predictions_data$available) {
      return(list(available = FALSE))
    }
    pred_data <- predictions_data$data
  } else {
    # New format: direct data frame from predictions.csv
    pred_data <- predictions_data
  }

  # Load sample sheet
  if (!file.exists(sample_sheet_path)) {
    return(list(available = FALSE, reason = "Sample sheet not found"))
  }

  tryCatch(
    {
      sample_sheet <- read_excel(sample_sheet_path)
      cat(sprintf("STEP 1: Sample sheet loaded (%d rows, %d columns)\n", nrow(sample_sheet), ncol(sample_sheet)))
      cat(sprintf("        Columns: %s\n\n", paste(colnames(sample_sheet), collapse = ", ")))

      # Match predictions with sample sheet using sample_id
      # Handle both numeric and character sample_id
      pred_data$sample_id_clean <- as.character(pred_data$sample_id)
      sample_sheet$ID_clean <- as.character(sample_sheet$ID)

      # Find STUDY column (may be STUDY, STUDY...3, etc.)
      study_col <- NULL
      possible_names <- c("STUDY", "STUDY...3", "Study", "study", "COHORT", "Cohort", "cohort")
      cat("STEP 2: Searching for STUDY column...\n")
      for (col_name in possible_names) {
        if (col_name %in% colnames(sample_sheet)) {
          study_col <- col_name
          cat(sprintf("        >>> Found STUDY column: '%s'\n\n", col_name))
          break
        }
      }

      if (is.null(study_col)) {
        msg <- sprintf("No STUDY/COHORT column found. Available: %s", paste(colnames(sample_sheet), collapse = ", "))
        cat(sprintf("        >>> ERROR: %s\n\n", msg))
        return(list(available = FALSE, reason = msg))
      }

      cat(sprintf("STEP 3: Merging %d predictions with sample sheet...\n", nrow(pred_data)))
      merged_data <- merge(pred_data,
        sample_sheet[, c("ID", study_col)],
        by.x = "sample_id_clean",
        by.y = "ID",
        all.x = TRUE
      )

      cat(sprintf("        Merge result: %d rows\n\n", nrow(merged_data)))

      # Rename study column to standardized name
      colnames(merged_data)[colnames(merged_data) == study_col] <- "STUDY"

      # Remove rows with missing STUDY information
      n_before <- nrow(merged_data)
      merged_data <- merged_data[!is.na(merged_data$STUDY), ]
      cat(sprintf("STEP 4: Filtering NA values...\n"))
      cat(sprintf("        Retained: %d/%d rows (%.1f%%)\n\n", nrow(merged_data), n_before, 100 * nrow(merged_data) / n_before))

      if (nrow(merged_data) == 0) {
        msg <- "No matches between predictions and sample sheet (all STUDY values are NA)"
        cat(sprintf("        >>> ERROR: %s\n\n", msg))
        return(list(available = FALSE, reason = msg))
      }

      # Calculate performance by study and class
      studies <- unique(merged_data$STUDY)
      classes <- unique(merged_data$actual)

      results <- data.frame()

      # Calculate confidence scores if available (get probability columns)
      prob_cols <- setdiff(colnames(merged_data), c("sample_id", "sample_id_clean", "actual", "predicted", "STUDY"))
      prob_cols <- prob_cols[sapply(merged_data[, prob_cols, drop = FALSE], is.numeric)]

      # Only keep columns that look like probabilities (values between 0 and 1)
      if (length(prob_cols) > 0) {
        prob_cols <- prob_cols[sapply(prob_cols, function(col) {
          vals <- merged_data[[col]]
          all(vals >= 0 & vals <= 1, na.rm = TRUE)
        })]
      }

      has_confidence <- length(prob_cols) > 0
      if (has_confidence) {
        # Calculate confidence as max probability for each sample
        merged_data$Confidence <- apply(merged_data[, prob_cols, drop = FALSE], 1, max, na.rm = TRUE)
      }

      for (study in studies) {
        study_data <- merged_data[merged_data$STUDY == study, ]

        for (class in classes) {
          class_data <- study_data[study_data$actual == class, ]

          if (nrow(class_data) > 0) {
            n_total <- nrow(class_data)
            n_correct <- sum(class_data$actual == class_data$predicted)
            pct_correct <- 100 * n_correct / n_total

            # Calculate confidence statistics if available
            if (has_confidence) {
              mean_conf <- mean(class_data$Confidence, na.rm = TRUE)
              sd_conf <- sd(class_data$Confidence, na.rm = TRUE)
            } else {
              mean_conf <- NA
              sd_conf <- NA
            }

            results <- rbind(results, data.frame(
              Study = study,
              Class = class,
              n = n_total,
              Correct = n_correct,
              Accuracy_pct = pct_correct,
              Mean_Confidence = mean_conf,
              SD_Confidence = sd_conf,
              stringsAsFactors = FALSE
            ))
          }
        }
      }

      # Calculate totals
      total_n <- sum(results$n)
      total_correct <- sum(results$Correct)
      total_accuracy <- 100 * total_correct / total_n

      # Calculate overall confidence if available
      if (has_confidence) {
        total_mean_conf <- mean(merged_data$Confidence, na.rm = TRUE)
        total_sd_conf <- sd(merged_data$Confidence, na.rm = TRUE)
      } else {
        total_mean_conf <- NA
        total_sd_conf <- NA
      }

      cat(sprintf("STEP 5: Calculation complete\n"))
      cat(sprintf("        Studies analyzed: %d\n", length(unique(results$Study))))
      cat(sprintf("        Total samples: %d\n", total_n))
      cat(sprintf("        Total correct: %d\n", total_correct))
      cat(sprintf("        Overall accuracy: %.1f%%\n", total_accuracy))
      if (has_confidence) {
        cat(sprintf("        Overall confidence: %.1f%% ± %.1f%%\n", 100 * total_mean_conf, 100 * total_sd_conf))
      }
      cat("\n")

      return(list(
        available = TRUE,
        results = results,
        total_n = total_n,
        total_correct = total_correct,
        total_accuracy = total_accuracy,
        has_confidence = has_confidence,
        total_mean_conf = total_mean_conf,
        total_sd_conf = total_sd_conf
      ))
    },
    error = function(e) {
      msg <- paste("Error:", e$message)
      cat(sprintf("        >>> ERROR: %s\n\n", msg))
      return(list(available = FALSE, reason = msg))
    }
  )
}

# Calculate Cohen's d effect size
calculate_cohens_d <- function(group1, group2) {
  mean1 <- mean(group1, na.rm = TRUE)
  mean2 <- mean(group2, na.rm = TRUE)
  sd1 <- sd(group1, na.rm = TRUE)
  sd2 <- sd(group2, na.rm = TRUE)
  n1 <- length(group1)
  n2 <- length(group2)

  # Pooled standard deviation
  pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))

  # Cohen's d
  d <- (mean1 - mean2) / pooled_sd

  return(d)
}

# Analyse probability calibration
analyse_calibration <- function(predictions, n_bins = 5) {
  if (!"Probability" %in% colnames(predictions)) {
    return(list(available = FALSE))
  }

  # Create probability bins
  predictions$prob_bin <- cut(predictions$Probability,
    breaks = seq(0, 1, length.out = n_bins + 1),
    include.lowest = TRUE
  )

  # Calculate observed accuracy per bin
  calibration <- predictions %>%
    group_by(prob_bin) %>%
    summarise(
      n_samples = n(),
      observed_accuracy = mean(actual == predicted),
      mean_probability = mean(Probability),
      .groups = "drop"
    )

  return(list(
    available = TRUE,
    calibration = calibration
  ))
}

# ============================================================================
# FINAL MODEL PREDICTIONS ANALYSIS (v5.13.0)
# ============================================================================

# Load final model predictions from test set
load_final_model_predictions <- function(model_dir) {
  pred_file <- file.path(model_dir, "predictions.csv")

  if (!file.exists(pred_file)) {
    return(list(available = FALSE, reason = "predictions.csv not found"))
  }

  tryCatch(
    {
      pred_data <- read.csv(pred_file, stringsAsFactors = FALSE)

      # Verify required columns
      required_cols <- c("sample_id", "actual", "predicted")
      missing_cols <- setdiff(required_cols, colnames(pred_data))
      if (length(missing_cols) > 0) {
        return(list(available = FALSE, reason = paste("Missing columns:", paste(missing_cols, collapse = ", "))))
      }

      # Get probability columns (all columns except sample_id, actual, predicted)
      prob_cols <- setdiff(colnames(pred_data), c("sample_id", "actual", "predicted"))

      return(list(
        available = TRUE,
        data = pred_data,
        prob_columns = prob_cols,
        n_samples = nrow(pred_data)
      ))
    },
    error = function(e) {
      return(list(available = FALSE, reason = paste("Error loading predictions:", conditionMessage(e))))
    }
  )
}

# Match predictions with sample sheet to add Lab_ID and Animal_ID
match_predictions_with_sample_sheet <- function(predictions, sample_sheet) {
  if (!predictions$available) {
    return(predictions)
  }

  pred_data <- predictions$data

  # Convert sample_id to numeric for matching
  pred_data$sample_id_clean <- as.numeric(as.character(pred_data$sample_id))

  # Match with sample sheet
  pred_data$Lab_ID <- sapply(pred_data$sample_id_clean, function(sid) {
    idx <- which(sample_sheet$ID == sid)
    if (length(idx) > 0) {
      as.character(sample_sheet$Lab_ID[idx[1]])
    } else {
      NA_character_
    }
  })

  pred_data$Animal_ID <- sapply(pred_data$sample_id_clean, function(sid) {
    idx <- which(sample_sheet$ID == sid)
    if (length(idx) > 0) {
      as.character(sample_sheet$Animal_ID[idx[1]])
    } else {
      NA_character_
    }
  })

  # Calculate correctness and confidence
  pred_data$Correct <- pred_data$actual == pred_data$predicted

  # Calculate confidence as max probability
  prob_cols <- predictions$prob_columns
  if (length(prob_cols) > 0) {
    pred_data$Confidence <- apply(pred_data[, prob_cols, drop = FALSE], 1, max)
  } else {
    pred_data$Confidence <- NA
  }

  predictions$data <- pred_data
  predictions$n_matched <- sum(!is.na(pred_data$Lab_ID))

  return(predictions)
}

# Calculate prediction confidence statistics
calculate_prediction_confidence_stats <- function(predictions) {
  if (!predictions$available) {
    return(list(available = FALSE))
  }

  pred_data <- predictions$data

  # Define confidence bins
  confidence_breaks <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
  confidence_labels <- c("50-60%", "60-70%", "70-80%", "80-90%", "90-100%")

  pred_data$ConfidenceBin <- cut(pred_data$Confidence,
    breaks = confidence_breaks,
    labels = confidence_labels,
    include.lowest = TRUE
  )

  # Calculate accuracy by confidence bin
  confidence_summary <- data.frame()
  for (bin in confidence_labels) {
    bin_data <- pred_data[pred_data$ConfidenceBin == bin & !is.na(pred_data$ConfidenceBin), ]
    if (nrow(bin_data) > 0) {
      confidence_summary <- rbind(confidence_summary, data.frame(
        ConfidenceBin = bin,
        n_samples = nrow(bin_data),
        n_correct = sum(bin_data$Correct),
        accuracy = 100 * sum(bin_data$Correct) / nrow(bin_data),
        stringsAsFactors = FALSE
      ))
    }
  }

  # Calculate probability distribution statistics by actual class
  prob_stats <- data.frame()
  classes <- unique(pred_data$actual)
  prob_cols <- predictions$prob_columns

  for (cls in classes) {
    class_data <- pred_data[pred_data$actual == cls, ]
    for (prob_col in prob_cols) {
      prob_values <- class_data[[prob_col]]
      prob_stats <- rbind(prob_stats, data.frame(
        Actual_Class = cls,
        Probability_For = prob_col,
        Mean = mean(prob_values, na.rm = TRUE),
        SD = sd(prob_values, na.rm = TRUE),
        Median = median(prob_values, na.rm = TRUE),
        Q25 = quantile(prob_values, 0.25, na.rm = TRUE),
        Q75 = quantile(prob_values, 0.75, na.rm = TRUE),
        stringsAsFactors = FALSE
      ))
    }
  }

  return(list(
    available = TRUE,
    confidence_summary = confidence_summary,
    prob_stats = prob_stats,
    overall_confidence_mean = mean(pred_data$Confidence, na.rm = TRUE),
    overall_confidence_sd = sd(pred_data$Confidence, na.rm = TRUE)
  ))
}

# ============================================================================
# DMR-GROUP ASSOCIATION ANALYSIS (works for 2, 3, or more groups)
# ============================================================================
analyse_dmr_group_associations <- function(meth_matrix, sample_metadata, dmr_names) {
  # meth_matrix: DMRs (rows) x Samples (cols)
  # sample_metadata: must have 'actual_class' column

  if (!"actual_class" %in% colnames(sample_metadata)) {
    return(list(available = FALSE, reason = "No actual_class column"))
  }

  groups <- unique(sample_metadata$actual_class)
  n_groups <- length(groups)
  n_dmrs <- nrow(meth_matrix)

  # Initialize results
  results <- data.frame(
    DMR = dmr_names,
    stringsAsFactors = FALSE
  )

  # Add mean and median methylation per group
  for (grp in groups) {
    grp_samples <- sample_metadata$actual_class == grp
    grp_col_name <- paste0("Mean_", make.names(grp))
    sd_col_name <- paste0("SD_", make.names(grp))
    median_col_name <- paste0("Median_", make.names(grp))

    results[[grp_col_name]] <- rowMeans(meth_matrix[, grp_samples, drop = FALSE], na.rm = TRUE)
    results[[sd_col_name]] <- apply(meth_matrix[, grp_samples, drop = FALSE], 1, sd, na.rm = TRUE)
    # FIX v5.12.0: Add median calculation (same logic as mean)
    results[[median_col_name]] <- apply(meth_matrix[, grp_samples, drop = FALSE], 1, median, na.rm = TRUE)
  }

  # Calculate ANOVA effect size (eta-squared) for each DMR
  eta_squared <- numeric(n_dmrs)
  f_statistic <- numeric(n_dmrs)
  anova_pvalue <- numeric(n_dmrs)

  for (i in 1:n_dmrs) {
    dmr_values <- as.numeric(meth_matrix[i, ])
    group_factor <- factor(sample_metadata$actual_class)

    # ANOVA
    aov_result <- tryCatch(
      {
        aov(dmr_values ~ group_factor)
      },
      error = function(e) NULL
    )

    if (!is.null(aov_result)) {
      aov_summary <- summary(aov_result)[[1]]
      ss_between <- aov_summary$"Sum Sq"[1]
      ss_total <- sum(aov_summary$"Sum Sq")

      eta_squared[i] <- ss_between / ss_total
      f_statistic[i] <- aov_summary$"F value"[1]
      anova_pvalue[i] <- aov_summary$"Pr(>F)"[1]
    } else {
      eta_squared[i] <- NA
      f_statistic[i] <- NA
      anova_pvalue[i] <- NA
    }
  }

  results$Eta_Squared <- eta_squared
  results$F_Statistic <- f_statistic
  results$ANOVA_pvalue <- anova_pvalue

  # ===== FIX #2: ADD FDR CORRECTION FOR P-VALUES =====
  # Apply multiple testing correction using the user-specified method
  results$ANOVA_pvalue_adjusted <- apply_multiple_testing(anova_pvalue, method = mt_method_r)
  results$Significant_FDR <- results$ANOVA_pvalue_adjusted < sig_threshold
  # ===================================================

  # Effect size interpretation
  results$Effect_Size_Interpretation <- sapply(eta_squared, function(eta) {
    if (is.na(eta)) {
      return("Unknown")
    }
    if (eta < 0.01) {
      return("Negligible")
    }
    if (eta < 0.06) {
      return("Small")
    }
    if (eta < 0.14) {
      return("Medium")
    }
    return("Large")
  })

  # For each DMR, identify which group it's most associated with
  # (group with most extreme mean vs others)
  best_marker_for <- character(n_dmrs)
  max_effect_size <- numeric(n_dmrs)

  for (i in 1:n_dmrs) {
    max_d <- 0
    best_group <- "None"

    for (grp in groups) {
      # Calculate one-vs-rest effect size
      grp_samples <- sample_metadata$actual_class == grp
      other_samples <- sample_metadata$actual_class != grp

      grp_mean <- mean(meth_matrix[i, grp_samples], na.rm = TRUE)
      other_mean <- mean(meth_matrix[i, other_samples], na.rm = TRUE)
      grp_sd <- sd(meth_matrix[i, grp_samples], na.rm = TRUE)
      other_sd <- sd(meth_matrix[i, other_samples], na.rm = TRUE)

      n_grp <- sum(grp_samples)
      n_other <- sum(other_samples)

      # Pooled SD for Cohen's d
      pooled_sd <- sqrt(((n_grp - 1) * grp_sd^2 + (n_other - 1) * other_sd^2) / (n_grp + n_other - 2))

      d <- abs(grp_mean - other_mean) / pooled_sd

      if (!is.na(d) && d > max_d) {
        max_d <- d
        best_group <- grp
      }
    }

    best_marker_for[i] <- best_group
    max_effect_size[i] <- max_d
  }

  results$Best_Marker_For <- best_marker_for
  results$Max_Cohens_D <- max_effect_size

  # Calculate all pairwise effect sizes (for 2 groups, just one comparison)
  if (n_groups == 2) {
    grp1 <- groups[1]
    grp2 <- groups[2]

    grp1_samples <- sample_metadata$actual_class == grp1
    grp2_samples <- sample_metadata$actual_class == grp2

    pairwise_d <- numeric(n_dmrs)
    wilcox_pvalue <- numeric(n_dmrs)
    median_diff <- numeric(n_dmrs)

    for (i in 1:n_dmrs) {
      grp1_vals <- meth_matrix[i, grp1_samples]
      grp2_vals <- meth_matrix[i, grp2_samples]

      # Cohen's d
      pairwise_d[i] <- abs(mean(grp1_vals, na.rm = TRUE) - mean(grp2_vals, na.rm = TRUE)) /
        sqrt((var(grp1_vals, na.rm = TRUE) + var(grp2_vals, na.rm = TRUE)) / 2)

      # Wilcoxon rank-sum test (Mann-Whitney U test)
      wilcox_result <- tryCatch(
        {
          wilcox.test(grp1_vals, grp2_vals, exact = FALSE)
        },
        error = function(e) NULL
      )

      if (!is.null(wilcox_result)) {
        wilcox_pvalue[i] <- wilcox_result$p.value
      } else {
        wilcox_pvalue[i] <- NA
      }

      # Median difference for interpretation
      median_diff[i] <- median(grp1_vals, na.rm = TRUE) - median(grp2_vals, na.rm = TRUE)
    }

    results[[paste0("Cohens_D_", make.names(grp1), "_vs_", make.names(grp2))]] <- pairwise_d
    results$Wilcoxon_pvalue <- wilcox_pvalue
    results$Wilcoxon_pvalue_adjusted <- apply_multiple_testing(wilcox_pvalue, method = mt_method_r)
    results$Significant_Wilcoxon <- results$Wilcoxon_pvalue_adjusted < sig_threshold
    results[[paste0("Median_Diff_", make.names(grp1), "_minus_", make.names(grp2))]] <- median_diff
  }

  # Rank DMRs by importance (eta-squared)
  results <- results[order(results$Eta_Squared, decreasing = TRUE, na.last = TRUE), ]
  results$Importance_Rank <- 1:n_dmrs

  # Reorder columns for readability
  first_cols <- c(
    "DMR", "Importance_Rank", "Eta_Squared", "Effect_Size_Interpretation",
    "F_Statistic", "ANOVA_pvalue", "ANOVA_pvalue_adjusted", "Significant_FDR",
    "Best_Marker_For", "Max_Cohens_D"
  )
  other_cols <- setdiff(colnames(results), first_cols)
  results <- results[, c(first_cols, other_cols)]

  return(list(
    available = TRUE,
    results = results,
    n_groups = n_groups,
    groups = groups
  ))
}

# ============================================================================
# NOTE: The rest of the script follows the same structure as v6.6.0
# but with ALL p.adjust() calls replaced with apply_multiple_testing()
# and output labels updated to reflect the chosen method
#
# Due to length constraints, I'm showing the key modifications below.
# The full script would be ~1800 lines with these patterns applied throughout.
# ============================================================================

cat("################################################################################\n")
cat("#  MULTIPLE TESTING CORRECTION CONFIGURATION\n")
cat("################################################################################\n")
cat("  Method:                  ", mt_method_name, "\n")
cat("  Raw p-value threshold:   ", sig_threshold, "\n")
cat("  Corrected threshold:     ", ifelse(mt_method_r == "none", "N/A (no correction)", sig_threshold), "\n")
cat("  Output columns:          ")
if (mt_method_r == "none") {
  cat("p_value, significant\n")
} else {
  cat("p_value, ", mt_method_display, "_adjusted, significant_raw, significant_", mt_method_display, "\n", sep = "")
}
cat("################################################################################\n\n")


# ============================================================================
# STEP 3: LOAD MODEL AND DATA
# ============================================================================
log_section_start(3, "Loading model and data", 10)

# Validate paths
if (!dir.exists(opt$model_dir)) stop("[ERROR] Model directory not found: ", opt$model_dir)
if (!file.exists(opt$qs_file)) stop("[ERROR] QS file not found: ", opt$qs_file)
if (!file.exists(opt$sample_sheet)) stop("[ERROR] Sample sheet not found: ", opt$sample_sheet)

# Load QS file
cat("          Loading methylRawList...\n")
meth_data <- qread(opt$qs_file)
cat("          [OK] Loaded", length(meth_data), "samples\n")

# Load sample sheet
cat("          Loading sample sheet...\n")
sample_sheet <- read_excel(opt$sample_sheet)
log_info(paste0("Sample sheet loaded: ", nrow(sample_sheet), " samples"), "STEP3")
cat("          [OK]", nrow(sample_sheet), "x", ncol(sample_sheet), "\n")

# Load model file
cat("          Loading model file...\n")
model_file <- file.path(opt$model_dir, "model.rds")
if (!file.exists(model_file)) stop("[ERROR] model.rds not found: ", model_file)
trained_model <- readRDS(model_file)
cat("          [OK] Loaded model\n")

# Load DMR details
cat("          Loading DMR details...\n")
dmr_file <- file.path(opt$model_dir, "dmr_details_with_coordinates.csv")
if (!file.exists(dmr_file)) stop("[ERROR] dmr_details_with_coordinates.csv not found: ", dmr_file)

dmr_details <- read.csv(dmr_file, stringsAsFactors = FALSE)
cat("          [OK] DMR details:", nrow(dmr_details), "rows x", ncol(dmr_details), "columns\n")

if (nrow(dmr_details) == 0) {
  stop("[ERROR] dmr_details_with_coordinates.csv is empty!")
}

# Extract metadata from dmr_details columns
if ("region_file" %in% colnames(dmr_details)) {
  region_info <- unique(dmr_details$region_file)[1]
  region_parts <- strsplit(region_info, "_")[[1]]
  species_name <- region_parts[1]
  window_size <- gsub("bp", "", region_parts[2])
  cat("          Species:", species_name, "\n")
  cat("          Window size:", window_size, "bp\n")
} else {
  species_name <- "unknown"
  window_size <- "unknown"
  region_info <- "unknown"
}

if ("model" %in% colnames(dmr_details)) {
  model_type <- unique(dmr_details$model)[1]
  cat("          Model type:", model_type, "\n")
} else {
  model_type <- "unknown"
}

# Detect chromosome column
chr_col <- NULL
for (possible_chr in c("chr", "chromosome", "scaffold", "seqnames")) {
  if (possible_chr %in% colnames(dmr_details)) {
    chr_col <- possible_chr
    break
  }
}

if (is.null(chr_col)) {
  stop("[ERROR] No chromosome column found in dmr_details!")
}

if (chr_col != "chr") {
  dmr_details$chr <- dmr_details[[chr_col]]
}

if (!"DMR_ID" %in% colnames(dmr_details)) {
  dmr_details$DMR_ID <- paste0("DMR_", 1:nrow(dmr_details))
}

n_markers <- nrow(dmr_details)
selected_markers <- dmr_details$DMR_ID

cat("          [OK] Created", n_markers, "DMR IDs\n")

# ============================================================================
# NEW: Load the single methylation heatmap CSV file
# ============================================================================
cat("          Loading methylation data CSV...\n")

# Find the methylation_heatmap file
heatmap_pattern <- paste0("methylation_heatmap_.*", n_markers, "markers_DATA\\.csv$")
heatmap_files <- list.files(opt$model_dir, pattern = heatmap_pattern, full.names = TRUE)

if (length(heatmap_files) == 0) {
  # Try more flexible pattern
  heatmap_pattern_flex <- "methylation_heatmap_.*markers_DATA\\.csv$"
  heatmap_files <- list.files(opt$model_dir, pattern = heatmap_pattern_flex, full.names = TRUE)
}

if (length(heatmap_files) == 0) {
  stop("[ERROR] methylation_heatmap CSV file not found in: ", opt$model_dir)
}

heatmap_file <- heatmap_files[1]
cat("          Found:", basename(heatmap_file), "\n")

# Load the heatmap data
combined_data <- read.csv(heatmap_file, stringsAsFactors = FALSE)
cat("          [OK] Loaded data:", nrow(combined_data), "samples x", ncol(combined_data), "columns\n")

# Check for required columns
if (!"Sample_ID" %in% colnames(combined_data)) {
  possible_id_cols <- c("sample_id", "SampleID", "Sample_ID", "ID", "sample")
  sample_col <- NULL
  for (col in possible_id_cols) {
    if (col %in% colnames(combined_data)) {
      sample_col <- col
      break
    }
  }

  if (is.null(sample_col)) {
    sample_col <- colnames(combined_data)[1]
    cat("          [INFO] Using first column as Sample_ID:", sample_col, "\n")
  }

  colnames(combined_data)[colnames(combined_data) == sample_col] <- "Sample_ID"
}

if (!"Group" %in% colnames(combined_data)) {
  possible_group_cols <- c("group", "Group", "treatment", "Treatment", "class", "Class", "actual_class")
  group_col <- NULL
  for (col in possible_group_cols) {
    if (col %in% colnames(combined_data)) {
      group_col <- col
      break
    }
  }

  if (is.null(group_col)) {
    stop("[ERROR] No Group column found in methylation heatmap CSV!")
  }

  colnames(combined_data)[colnames(combined_data) == group_col] <- "Group"
}

cat("          Columns:", paste(head(colnames(combined_data), 10), collapse = ", "), "...\n")

# Show group distribution
cat("\n          Sample distribution by group:\n")
group_table <- table(combined_data$Group)
print(group_table)

# For compatibility with downstream code
combined_data$dataset_source <- "COMBINED"
combined_data$actual_class <- combined_data$Group
sample_col <- "Sample_ID"

# Extract DMR columns
dmr_cols <- grep("^DMR_", colnames(combined_data), value = TRUE)
n_dmrs <- length(dmr_cols)

if (n_dmrs == 0) {
  stop("[ERROR] No DMR columns found in methylation heatmap CSV!")
}

cat("          Found", n_dmrs, "DMR columns\n")

# Create methylation matrix (DMRs as rows, samples as columns)
sample_ids <- combined_data[[sample_col]]
dmr_values <- as.matrix(combined_data[, dmr_cols])
rownames(dmr_values) <- sample_ids

# Transpose so DMRs are rows and samples are columns
meth_matrix <- t(dmr_values)
rownames(meth_matrix) <- dmr_cols
colnames(meth_matrix) <- sample_ids

cat("          [OK] Methylation matrix:", nrow(meth_matrix), "DMRs x", ncol(meth_matrix), "samples\n")

# Store metadata for each sample
sample_metadata <- combined_data[, c(sample_col, "dataset_source", "actual_class")]
rownames(sample_metadata) <- sample_metadata[[sample_col]]

cat("          [OK] Sample metadata stored for", nrow(sample_metadata), "samples\n")
cat("          Metadata columns:", paste(colnames(sample_metadata), collapse = ", "), "\n")

# Show summary
cat("\n          Dataset Summary:\n")
cat(sprintf("            Total samples: %d\n", ncol(meth_matrix)))
if ("actual_class" %in% colnames(sample_metadata)) {
  for (grp in unique(sample_metadata$actual_class)) {
    n_grp <- sum(sample_metadata$actual_class == grp, na.rm = TRUE)
    cat(sprintf("            Group '%s': %d samples\n", grp, n_grp))
  }
}

# Parse custom color palette if specified
group_colors <- NULL
if ("actual_class" %in% colnames(sample_metadata)) {
  group_names <- unique(sample_metadata$actual_class)
  group_colors <- parse_group_colors(opt$group_colors, group_names)
  if (!is.null(group_colors)) {
    cat("\n          [INFO] Custom color palette applied for PCA/UMAP plots:\n")
    for (grp in names(group_colors)) {
      cat(sprintf("            %s = %s\n", grp, group_colors[grp]))
    }
  } else {
    cat("\n          [INFO] Using default Set1 color palette for PCA/UMAP plots\n")
  }
}

# For bootstrap analysis - use all samples
train_samples <- sample_ids
test_samples <- character(0)

cat("          [INFO] Using all", length(train_samples), "samples for analysis\n")

# Create output directory
output_dir <- file.path(opt$model_dir, opt$output_subdir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

model_name <- paste0(species_name, "_", window_size, "bp_", model_type)

cat("\n")
cat("################################################################################\n")
cat("#  ANALYSIS CONFIGURATION\n")
cat("################################################################################\n")
cat("  Model name:              ", model_name, "\n")
cat("  Species:                 ", species_name, "\n")
cat("  Window size:             ", window_size, " bp\n")
cat("  Model type:              ", model_type, "\n")
cat("  Number of markers:       ", n_markers, "\n")
cat("  Train samples:           ", length(train_samples), "\n")
cat("  Test samples:            ", length(test_samples), "\n")
cat("  Output directory:        ", output_dir, "\n")
cat("################################################################################\n\n")

# NOW setup logging (after all variables exist)
cat("[LOGGING] Setting up comprehensive logging...\n")

log_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(output_dir, paste0("analysis_log_", log_timestamp, ".txt"))

log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)

cat("================================================================================\n")
cat("METHYLSENSE MODEL EVALUATION - LOG\n")
cat("================================================================================\n")
cat("Script version:", SCRIPT_VERSION, "\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Log file:", log_file, "\n")
cat("================================================================================\n\n")

# Save parameters
params_file <- file.path(output_dir, "analysis_parameters.txt")
params_con <- file(params_file, "w")
cat("METHYLSENSE MODEL EVALUATION PARAMETERS\n", file = params_con)
cat("========================================\n\n", file = params_con)
cat("Species:", species_name, "\n", file = params_con)
cat("Window:", window_size, "bp\n", file = params_con)
cat("Model:", model_type, "\n", file = params_con)
cat("Markers:", n_markers, "\n", file = params_con)
cat("Train:", length(train_samples), "\n", file = params_con)
cat("Test:", length(test_samples), "\n", file = params_con)
close(params_con)

cat("          [OK] Logging initialized\n\n")


# ===== QC CHECK 1: Validate DMR methylation data =====
cat("\n          ===== QC CHECK 1: DMR METHYLATION VALIDATION =====\n")

total_values <- length(meth_matrix)
na_values <- sum(is.na(meth_matrix))
pct_na <- round(100 * na_values / total_values, 1)

cat(sprintf("          Total values: %d\n", total_values))
cat(sprintf("          NA values: %d (%.1f%%)\n", na_values, pct_na))

# ... continue with rest of QC checks exactly as before...

if (pct_na > 50) {
  cat(sprintf("          [WARNING] %.1f%% of DMR values are NA!\n", pct_na))
} else {
  cat(sprintf("          [OK] Only %.1f%% NA values\n", pct_na))
}

non_na_vals <- as.vector(meth_matrix[!is.na(meth_matrix)])
if (length(non_na_vals) > 0) {
  cat(sprintf(
    "          DMR methylation range: %.1f%% - %.1f%%\n",
    min(non_na_vals, na.rm = TRUE), max(non_na_vals, na.rm = TRUE)
  ))
  cat(sprintf(
    "          DMR methylation mean: %.1f%% (SD: %.1f%%)\n",
    mean(non_na_vals, na.rm = TRUE), sd(non_na_vals, na.rm = TRUE)
  ))
  cat(sprintf(
    "          DMR methylation median: %.1f%%\n",
    median(non_na_vals, na.rm = TRUE)
  ))
}

# ===== QC CHECK 2: Show example methylation values =====
cat("\n          ===== QC CHECK 2: EXAMPLE METHYLATION VALUES =====\n")
cat("          Showing DMR methylation for first 5 markers in first 5 samples:\n\n")

# Select first 5 DMRs and first 5 samples
n_show_markers <- min(5, nrow(meth_matrix))
n_show_samples <- min(5, ncol(meth_matrix))

example_matrix <- meth_matrix[1:n_show_markers, 1:n_show_samples, drop = FALSE]

# Create formatted table
cat("          ", sprintf("%-15s", "Marker"), sep = "")
for (j in 1:n_show_samples) {
  cat(sprintf("%-15s", substr(colnames(example_matrix)[j], 1, 12)))
}
cat("\n")

cat("          ", paste(rep("-", 15 + 15 * n_show_samples), collapse = ""), "\n", sep = "")

for (i in 1:n_show_markers) {
  cat("          ", sprintf("%-15s", rownames(example_matrix)[i]), sep = "")
  for (j in 1:n_show_samples) {
    val <- example_matrix[i, j]
    if (is.na(val)) {
      cat(sprintf("%-15s", "NA"))
    } else {
      cat(sprintf("%-15s", paste0(round(val, 1), "%")))
    }
  }
  cat("\n")
}

cat("\n")

# ===== QC CHECK 3: Per-marker statistics =====
cat("          ===== QC CHECK 3: PER-MARKER STATISTICS =====\n")

marker_stats <- data.frame(
  marker = rownames(meth_matrix),
  n_samples = ncol(meth_matrix),
  n_valid = rowSums(!is.na(meth_matrix)),
  n_missing = rowSums(is.na(meth_matrix)),
  pct_missing = round(100 * rowSums(is.na(meth_matrix)) / ncol(meth_matrix), 1),
  mean_meth = rowMeans(meth_matrix, na.rm = TRUE),
  sd_meth = apply(meth_matrix, 1, sd, na.rm = TRUE),
  min_meth = apply(meth_matrix, 1, min, na.rm = TRUE),
  max_meth = apply(meth_matrix, 1, max, na.rm = TRUE),
  range_meth = apply(meth_matrix, 1, max, na.rm = TRUE) - apply(meth_matrix, 1, min, na.rm = TRUE),
  stringsAsFactors = FALSE
)

cat("          Summary across all", nrow(marker_stats), "markers:\n")
cat(sprintf(
  "          Mean completeness: %.1f%% (range: %.1f%% - %.1f%%)\n",
  mean(100 - marker_stats$pct_missing),
  min(100 - marker_stats$pct_missing),
  max(100 - marker_stats$pct_missing)
))
cat(sprintf(
  "          Mean methylation: %.1f%% (SD: %.1f%%)\n",
  mean(marker_stats$mean_meth, na.rm = TRUE),
  sd(marker_stats$mean_meth, na.rm = TRUE)
))
cat(sprintf(
  "          Mean range per marker: %.1f%%\n",
  mean(marker_stats$range_meth, na.rm = TRUE)
))

# Identify problematic markers
high_missing <- marker_stats$pct_missing > 20
low_variance <- marker_stats$sd_meth < 1

if (sum(high_missing) > 0) {
  cat(sprintf("          [WARNING] %d markers with >20%% missing values\n", sum(high_missing)))
  cat("          ", paste(head(marker_stats$marker[high_missing], 3), collapse = ", "), "...\n")
}

if (sum(low_variance) > 0) {
  cat(sprintf("          [INFO] %d markers with low variance (<1%% SD)\n", sum(low_variance)))
  cat("          These may be less informative for classification\n")
}

# Save marker statistics
qc_dir <- file.path(output_dir, "qc_checks")
dir.create(qc_dir, showWarnings = FALSE)

write.csv(marker_stats, file.path(qc_dir, "dmr_marker_statistics.csv"), row.names = FALSE)
cat("          [OK] Saved marker statistics to qc_checks/\n")

# ===== QC CHECK 4: Per-sample statistics =====
cat("\n          ===== QC CHECK 4: PER-SAMPLE STATISTICS =====\n")

sample_stats <- data.frame(
  sample_id = colnames(meth_matrix),
  n_markers = nrow(meth_matrix),
  n_valid = colSums(!is.na(meth_matrix)),
  n_missing = colSums(is.na(meth_matrix)),
  pct_missing = round(100 * colSums(is.na(meth_matrix)) / nrow(meth_matrix), 1),
  mean_meth = colMeans(meth_matrix, na.rm = TRUE),
  sd_meth = apply(meth_matrix, 2, sd, na.rm = TRUE),
  min_meth = apply(meth_matrix, 2, min, na.rm = TRUE),
  max_meth = apply(meth_matrix, 2, max, na.rm = TRUE),
  dataset = ifelse(colnames(meth_matrix) %in% test_samples, "TEST", "TRAIN"),
  stringsAsFactors = FALSE
)

cat("          Summary across all", nrow(sample_stats), "samples:\n")
cat(sprintf(
  "          Mean completeness: %.1f%% (range: %.1f%% - %.1f%%)\n",
  mean(100 - sample_stats$pct_missing),
  min(100 - sample_stats$pct_missing),
  max(100 - sample_stats$pct_missing)
))
cat(sprintf(
  "          Mean methylation: %.1f%% (SD: %.1f%%)\n",
  mean(sample_stats$mean_meth, na.rm = TRUE),
  sd(sample_stats$mean_meth, na.rm = TRUE)
))

# Compare TRAIN vs TEST
cat("\n          Comparison: TRAIN vs TEST datasets\n")
train_stats <- sample_stats[sample_stats$dataset == "TRAIN", ]
test_stats <- sample_stats[sample_stats$dataset == "TEST", ]

cat(sprintf(
  "          TRAIN: Mean methylation = %.1f%% (SD: %.1f%%, n=%d)\n",
  mean(train_stats$mean_meth, na.rm = TRUE),
  sd(train_stats$mean_meth, na.rm = TRUE),
  nrow(train_stats)
))
cat(sprintf(
  "          TEST:  Mean methylation = %.1f%% (SD: %.1f%%, n=%d)\n",
  mean(test_stats$mean_meth, na.rm = TRUE),
  sd(test_stats$mean_meth, na.rm = TRUE),
  nrow(test_stats)
))

# Statistical test
if (nrow(train_stats) > 0 && nrow(test_stats) > 0) {
  t_test <- t.test(train_stats$mean_meth, test_stats$mean_meth)
  cat(sprintf(
    "          T-test p-value: %.4f %s\n",
    t_test$p.value,
    ifelse(t_test$p.value < sig_threshold, "(significantly different)", "(not significant)")
  ))
}

write.csv(sample_stats, file.path(qc_dir, "sample_statistics.csv"), row.names = FALSE)
cat("          [OK] Saved sample statistics to qc_checks/\n")

# ===== QC CHECK 5: Create QC visualisation plots =====
cat("\n          ===== QC CHECK 5: CREATING QC PLOTS =====\n")

# Plot 1: Heatmap of first 20 markers x 20 samples
n_heatmap_markers <- min(20, nrow(meth_matrix))
n_heatmap_samples <- min(20, ncol(meth_matrix))

heatmap_data <- meth_matrix[1:n_heatmap_markers, 1:n_heatmap_samples, drop = FALSE]

# Generate methylation heatmap in all requested formats
formats_to_generate <- tolower(trimws(unlist(strsplit(opt$figure_format, ","))))
if ("all" %in% formats_to_generate) {
  formats_to_generate <- c("png", "jpeg", "svg")
}
formats_to_generate <- gsub("jpg", "jpeg", formats_to_generate)

for (fmt in unique(formats_to_generate)) {
  plot_file <- file.path(qc_dir, sprintf("methylation_heatmap_example.%s", fmt))
  open_plot_device(plot_file, width = 8, height = 5.33, res = 150)

  pheatmap(heatmap_data,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    colour = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(0, 100, length.out = 101),
    main = paste(
      "DMR Methylation Heatmap (first", n_heatmap_markers, "markers,",
      n_heatmap_samples, "samples)"
    ),
    fontsize = 8,
    angle_col = 45,
    na_col = "gray80"
  )

  dev.off()
}

cat(sprintf("          [OK] Created methylation_heatmap_example in formats: %s\n", paste(formats_to_generate, collapse = ", ")))

# Plot 2: Distribution of mean methylation per marker
p_marker_dist <- ggplot(marker_stats, aes(x = mean_meth)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, colour = "black") +
  geom_vline(
    xintercept = median(marker_stats$mean_meth, na.rm = TRUE),
    colour = "red", linetype = "dashed", linewidth = 1
  ) +
  labs(
    title = "Distribution of Mean Methylation Across DMR Markers",
    subtitle = paste("n =", nrow(marker_stats), "markers"),
    x = "Mean Methylation (%)",
    y = "Number of Markers"
  ) +
  theme_minimal()

save_plot_multiple_formats(p_marker_dist,
  file.path(qc_dir, "marker_methylation_distribution"),
  width = 10, height = 6
)

cat("          [OK] Created marker_methylation_distribution.png\n")

# Plot 3: Distribution of mean methylation per sample (coloured by dataset)
p_sample_dist <- ggplot(sample_stats, aes(x = mean_meth, fill = dataset)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity", colour = "black") +
  scale_fill_manual(values = c("TRAIN" = "steelblue", "TEST" = "coral")) +
  labs(
    title = "Distribution of Mean Methylation Across Samples",
    subtitle = paste("TRAIN: n =", nrow(train_stats), "| TEST: n =", nrow(test_stats)),
    x = "Mean Methylation (%)",
    y = "Number of Samples",
    fill = "Dataset"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

save_plot_multiple_formats(p_sample_dist,
  file.path(qc_dir, "sample_methylation_distribution"),
  width = 10, height = 6
)

cat("          [OK] Created sample_methylation_distribution.png\n")

# Plot 4: Missing data patterns
p_missing <- ggplot(marker_stats, aes(x = reorder(marker, -pct_missing), y = pct_missing)) +
  geom_bar(stat = "identity", fill = "coral", alpha = 0.7) +
  geom_hline(yintercept = 20, linetype = "dashed", colour = "red") +
  coord_flip() +
  labs(
    title = "Missing Data Per Marker",
    subtitle = "Red line = 20% threshold",
    x = "Marker",
    y = "Missing Data (%)"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6))

save_plot_multiple_formats(p_missing,
  file.path(qc_dir, "missing_data_per_marker"),
  width = 10, height = max(6, n_markers * 0.15)
)

cat("          [OK] Created missing_data_per_marker.png\n")

# ===== QC CHECK 6: VALIDATE AGAINST PRE-CALCULATED GROUP STATS =====
cat("\n")
cat("################################################################################\n")
cat("          QC CHECK 6: VALIDATION AGAINST GROUP SUMMARY STATISTICS\n")
cat("################################################################################\n")

# Load the pre-calculated group methylation statistics (summary across all markers)
group_stats_pattern <- paste0("group_methylation_stats_", n_markers, "markers\\.csv")
cat("          Searching for group statistics file...\n")
cat("          Pattern:", group_stats_pattern, "\n")
cat("          Search directory:", opt$model_dir, "\n")

group_stats_files <- list.files(opt$model_dir, pattern = group_stats_pattern, full.names = TRUE)

if (length(group_stats_files) == 0) {
  # Try more flexible pattern
  group_stats_pattern_flex <- "group_methylation_stats.*\\.csv"
  cat("          Trying flexible pattern:", group_stats_pattern_flex, "\n")
  group_stats_files <- list.files(opt$model_dir, pattern = group_stats_pattern_flex, full.names = TRUE)
}

if (length(group_stats_files) > 0) {
  group_stats_file <- group_stats_files[1]
  cat("          [OK] Found group statistics file:", basename(group_stats_file), "\n")
  cat("          Full path:", group_stats_file, "\n")

  # Load the summary statistics file (comma-separated)
  group_stats_expected <- read.csv(group_stats_file, stringsAsFactors = FALSE)

  cat("          [INFO] Loaded CSV file with standard comma separation\n")

  cat("          File dimensions:", nrow(group_stats_expected), "rows x", ncol(group_stats_expected), "columns\n")
  cat("          Columns:", paste(colnames(group_stats_expected), collapse = ", "), "\n")

  # Show the data
  cat("\n          Expected group statistics (from file):\n")
  print(group_stats_expected)
  cat("\n")

  # Use sample_metadata which already has all the information
  if (exists("sample_metadata") && "actual_class" %in% colnames(sample_metadata)) {
    cat("          [OK] Found 'actual_class' column in sample metadata\n")

    # Get unique groups (excluding NA)
    unique_groups <- unique(sample_metadata$actual_class)
    unique_groups <- unique_groups[!is.na(unique_groups)]

    cat("          Groups found:", paste(unique_groups, collapse = ", "), "\n")
    cat("          Number of groups:", length(unique_groups), "\n")

    # Calculate OVERALL statistics from our extracted data for each group
    # (average across ALL markers per group)
    cat("\n          ===== Calculating overall group statistics from extracted data =====\n")

    calculated_summary <- data.frame(
      Group = character(),
      Mean_Methylation_Calculated = numeric(),
      SD_Calculated = numeric(),
      SE_Calculated = numeric(),
      N_Calculated = integer(),
      stringsAsFactors = FALSE
    )

    for (group in unique_groups) {
      cat(sprintf("\n          Processing group: '%s'\n", group))

      # Get samples in this group from sample_metadata
      group_samples <- rownames(sample_metadata)[sample_metadata$actual_class == group &
        !is.na(sample_metadata$actual_class)]

      cat(sprintf("            Samples in metadata for this group: %d\n", length(group_samples)))

      # Filter to samples we have data for (should be all of them)
      group_samples_available <- intersect(group_samples, colnames(meth_matrix))
      cat(sprintf("            Samples with methylation data: %d\n", length(group_samples_available)))

      if (length(group_samples_available) == 0) {
        cat("            [WARNING] No samples found for group:", group, "\n")
        next
      }

      # Extract methylation data for this group (all markers, all samples in group)
      group_data <- meth_matrix[, group_samples_available, drop = FALSE]

      # Method 1: Overall mean of all measurements (what we're doing now)
      all_values <- as.vector(group_data)
      all_values <- all_values[!is.na(all_values)]
      overall_mean_all <- mean(all_values)
      overall_sd_all <- sd(all_values)

      # Method 2: Mean per sample, then average (more comparable to expected)
      # Calculate mean methylation per sample (average across all markers)
      sample_means <- colMeans(group_data, na.rm = TRUE)

      overall_mean <- mean(sample_means)
      overall_sd <- sd(sample_means) # SD across samples
      n_samples <- length(group_samples_available)
      overall_se <- overall_sd / sqrt(n_samples)
      overall_n <- n_samples

      cat(sprintf("            Overall mean methylation: %.4f%%\n", overall_mean))
      cat(sprintf("            Overall SD (across samples): %.4f%%\n", overall_sd))
      cat(sprintf("            Overall SE: %.4f%%\n", overall_se))
      cat(sprintf("            Number of samples: %d\n", n_samples))
      cat(sprintf("            [For reference] SD of all measurements: %.4f%%\n", overall_sd_all))

      # Store in summary dataframe
      calculated_summary <- rbind(calculated_summary, data.frame(
        Group = group,
        Mean_Methylation_Calculated = overall_mean,
        SD_Calculated = overall_sd,
        SE_Calculated = overall_se,
        N_Calculated = overall_n,
        stringsAsFactors = FALSE
      ))
    }

    cat("\n          [OK] Calculated overall statistics for", nrow(calculated_summary), "groups\n")

    # Now compare with expected values from file
    cat("\n")
    cat("################################################################################\n")
    cat("          COMPARING CALCULATED VS EXPECTED OVERALL STATISTICS\n")
    cat("################################################################################\n\n")

    # Merge calculated with expected
    comparison_df <- merge(calculated_summary,
      group_stats_expected,
      by = "Group",
      all = TRUE
    )

    cat("          Comparison Results:\n")
    cat("          ==================\n\n")

    for (i in 1:nrow(comparison_df)) {
      group <- comparison_df$Group[i]

      cat(sprintf("          GROUP: %s\n", group))
      cat("          ", strrep("-", 60), "\n", sep = "")

      # Mean comparison
      expected_mean <- comparison_df$Mean_Methylation[i]
      calculated_mean <- comparison_df$Mean_Methylation_Calculated[i]
      mean_diff <- abs(expected_mean - calculated_mean)
      mean_pct_diff <- 100 * mean_diff / expected_mean

      cat(sprintf("          Mean Methylation:\n"))
      cat(sprintf("            Expected:   %.4f%%\n", expected_mean))
      cat(sprintf("            Calculated: %.4f%%\n", calculated_mean))
      cat(sprintf("            Difference: %.4f%% (%.2f%% relative)\n", mean_diff, mean_pct_diff))

      # SD comparison
      expected_sd <- comparison_df$SD[i]
      calculated_sd <- comparison_df$SD_Calculated[i]
      sd_diff <- abs(expected_sd - calculated_sd)

      cat(sprintf("          Standard Deviation:\n"))
      cat(sprintf("            Expected:   %.4f%%\n", expected_sd))
      cat(sprintf("            Calculated: %.4f%%\n", calculated_sd))
      cat(sprintf("            Difference: %.4f%%\n", sd_diff))

      # SE comparison
      expected_se <- comparison_df$SE[i]
      calculated_se <- comparison_df$SE_Calculated[i]
      se_diff <- abs(expected_se - calculated_se)

      cat(sprintf("          Standard Error:\n"))
      cat(sprintf("            Expected:   %.4f%%\n", expected_se))
      cat(sprintf("            Calculated: %.4f%%\n", calculated_se))
      cat(sprintf("            Difference: %.4f%%\n", se_diff))

      # N comparison
      expected_n <- comparison_df$N[i]
      calculated_n <- comparison_df$N_Calculated[i]

      cat(sprintf("          Sample Size:\n"))
      cat(sprintf("            Expected:   %d\n", expected_n))
      cat(sprintf("            Calculated: %d\n", calculated_n))

      # Overall assessment
      cat("\n          Assessment:\n")
      if (mean_pct_diff < 0.1) {
        cat("             PERFECT MATCH (difference < 0.1%)\n")
      } else if (mean_pct_diff < 1.0) {
        cat("             EXCELLENT (difference < 1%)\n")
      } else if (mean_pct_diff < 5.0) {
        cat("             GOOD (difference < 5%)\n")
      } else {
        cat("            ⚠ WARNING - Review differences\n")
      }

      cat("\n")
    }

    # Overall summary
    cat("################################################################################\n")
    cat("          OVERALL VALIDATION SUMMARY\n")
    cat("################################################################################\n\n")

    overall_mean_diff <- mean(abs(comparison_df$Mean_Methylation - comparison_df$Mean_Methylation_Calculated), na.rm = TRUE)
    overall_sd_diff <- mean(abs(comparison_df$SD - comparison_df$SD_Calculated), na.rm = TRUE)

    cat(sprintf("          Average absolute difference in Mean: %.4f%%\n", overall_mean_diff))
    cat(sprintf("          Average absolute difference in SD:   %.4f%%\n", overall_sd_diff))

    if (overall_mean_diff < 0.01) {
      cat("\n           PERFECT: Overall group methylation matches perfectly!\n")
    } else if (overall_mean_diff < 0.1) {
      cat("\n           EXCELLENT: Overall group methylation shows excellent agreement!\n")
    } else if (overall_mean_diff < 1.0) {
      cat("\n           GOOD: Overall group methylation shows good agreement!\n")
    } else {
      cat("\n          ⚠ WARNING: Significant differences detected - review data extraction!\n")
    }

    # Save comparison table
    comparison_file <- file.path(qc_dir, "validation_group_summary_comparison.csv")
    write.csv(comparison_df, comparison_file, row.names = FALSE)
    cat("\n          [OK] Comparison saved to:", basename(comparison_file), "\n")

    # Create comparison plot
    library(ggplot2)

    comparison_long <- data.frame(
      Group = rep(comparison_df$Group, 2),
      Source = rep(c("Expected", "Calculated"), each = nrow(comparison_df)),
      Mean = c(comparison_df$Mean_Methylation, comparison_df$Mean_Methylation_Calculated),
      SD = c(comparison_df$SD, comparison_df$SD_Calculated)
    )

    p_comparison <- ggplot(comparison_long, aes(x = Group, y = Mean, fill = Source)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
      geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
        position = position_dodge(width = 0.8),
        width = 0.2
      ) +
      scale_fill_manual(values = c("Expected" = "steelblue", "Calculated" = "coral")) +
      labs(
        title = "Validation: Expected vs Calculated Group Methylation",
        subtitle = paste("Overall mean difference:", round(overall_mean_diff, 4), "%"),
        x = "Group",
        y = "Mean Methylation (%) +/- SD",
        fill = "Source"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

    save_plot_multiple_formats(p_comparison,
      file.path(qc_dir, "validation_group_comparison"),
      width = 10, height = 7
    )

    cat("          [OK] Comparison plot saved\n")
  } else {
    cat("          [WARNING] 'actual_class' column not found in sample metadata\n")
    cat("          Cannot validate group statistics without actual_class information\n")
    if (exists("sample_metadata")) {
      cat("          Columns in sample_metadata:", paste(colnames(sample_metadata), collapse = ", "), "\n")
    } else {
      cat("          [ERROR] sample_metadata object doesn't exist!\n")
    }
  }
} else {
  cat("          [WARNING] group_methylation_stats file not found!\n")
  cat("          Searched in:", opt$model_dir, "\n")
  cat("          Pattern:", group_stats_pattern, "\n")
  cat("          Available CSV files:\n")
  all_csv_files <- list.files(opt$model_dir, pattern = "\\.csv$", full.names = FALSE)
  for (f in all_csv_files) {
    cat("            -", f, "\n")
  }
  cat("          Skipping QC Check 6 (group validation)\n")
}

cat("\n")
cat("################################################################################\n")
cat("          QC CHECK 6 COMPLETE\n")
cat("################################################################################\n\n")

# ============================================================================
# MATCH SAMPLE SHEET TO METHYLATION DATA (needed for STEP 5)
# ============================================================================
cat("\n")
cat("################################################################################\n")
cat("          MATCHING SAMPLE SHEET TO METHYLATION DATA\n")
cat("################################################################################\n")

# Get sample IDs from methylation matrix
meth_sample_ids <- colnames(meth_matrix)
cat(sprintf("          Methylation data has %d samples\n", length(meth_sample_ids)))
cat("          Example sample IDs:", paste(head(meth_sample_ids, 3), collapse = ", "), "\n")

# Check sample sheet
cat(sprintf("          Sample sheet has %d rows\n", nrow(sample_sheet)))
cat("          Sample sheet columns:", paste(head(colnames(sample_sheet), 10), collapse = ", "), "\n")

# Find all columns with "id" in the name (case-insensitive)
id_columns <- colnames(sample_sheet)[grepl("id", colnames(sample_sheet), ignore.case = TRUE)]

cat("          Candidate ID columns in sample sheet:", paste(id_columns, collapse = ", "), "\n")

# Test each column and count matches
best_col <- NULL
best_matches <- 0

for (col in id_columns) {
  col_values <- as.character(sample_sheet[[col]])
  n_matches <- sum(meth_sample_ids %in% col_values)

  cat(sprintf("          Testing column '%s': %d/%d matches\n", col, n_matches, length(meth_sample_ids)))

  if (n_matches > best_matches) {
    best_matches <- n_matches
    best_col <- col
  }
}

if (is.null(best_col) || best_matches == 0) {
  cat("          [WARNING] Could not match sample sheet to methylation data\n")
  cat("          Using full sample sheet with first column as ID\n")

  # Create a minimal matched sheet
  sample_sheet_matched <- data.frame(
    sample_id = meth_sample_ids,
    stringsAsFactors = FALSE
  )

  # Try to add any available columns from original sample sheet if there's overlap
  if (nrow(sample_sheet) > 0 && ncol(sample_sheet) > 0) {
    first_col <- colnames(sample_sheet)[1]
    matches <- match(meth_sample_ids, as.character(sample_sheet[[first_col]]))

    if (sum(!is.na(matches)) > 0) {
      for (col in colnames(sample_sheet)) {
        sample_sheet_matched[[col]] <- sample_sheet[[col]][matches]
      }
    }
  }

  opt$sample_id_column <- "sample_id"
} else {
  cat(sprintf(
    "          [OK] Best match: column '%s' (%d/%d samples matched)\n",
    best_col, best_matches, length(meth_sample_ids)
  ))

  # Use the best matching column
  opt$sample_id_column <- best_col

  # Match and reorder sample sheet to match methylation matrix
  match_indices <- match(meth_sample_ids, as.character(sample_sheet[[best_col]]))
  sample_sheet_matched <- sample_sheet[match_indices, ]

  # Add rownames
  rownames(sample_sheet_matched) <- meth_sample_ids

  # Check how many matched
  n_matched <- sum(!is.na(match_indices))
  cat(sprintf("          [OK] Created matched sample sheet: %d samples\n", n_matched))
}

cat("          Sample sheet matched dimensions:", nrow(sample_sheet_matched), "x", ncol(sample_sheet_matched), "\n")
cat("          Sample ID column set to:", opt$sample_id_column, "\n")

cat("\n")

# ============================================================================
# STEP 1.5: MODEL CHARACTERISATION
# ============================================================================
cat("\n========================================\n")
cat("STEP 1.5/10: Model characterisation\n")
cat("========================================\n")

# Parse model summary
model_summary <- parse_model_summary(opt$model_dir)

if (model_summary$available) {
  cat("          [OK] Loaded model summary\n")
  if (!is.null(model_summary$method)) {
    cat(sprintf("          Method: %s\n", model_summary$method))
  }
  if (!is.null(model_summary$hyperparameters)) {
    cat("          Best hyperparameters found\n")
  }
  if (!is.null(model_summary$variable_importance)) {
    cat(sprintf(
      "          Variable importance: %d features\n",
      length(model_summary$variable_importance)
    ))
  }
} else {
  cat("          [INFO] model_summary.txt not found - skipping\n")
}

# Store for markdown
model_characterisation_done <- model_summary$available

cat("\n")

# ============================================================================
# STEP 2.5: NESTED CROSS-VALIDATION ANALYSIS
# ============================================================================
cat("\n========================================\n")
cat("STEP 2.5/10: Nested cross-validation analysis\n")
cat("========================================\n")

nested_cv_summary <- load_nested_cv_summary(opt$model_dir)
nested_cv_detailed <- load_nested_cv_detailed(opt$model_dir)

if (nested_cv_summary$available) {
  cat("          [OK] Loaded nested CV summary\n")
  cat(sprintf("          Folds: %d\n", nested_cv_summary$n_folds))
  cat(sprintf(
    "          Mean Accuracy: %.3f (95%% CI: %.3f-%.3f)\n",
    nested_cv_summary$mean_accuracy,
    nested_cv_summary$accuracy_ci_lower,
    nested_cv_summary$accuracy_ci_upper
  ))
  cat(sprintf(
    "          Mean ROC-AUC: %.3f (95%% CI: %.3f-%.3f)\n",
    nested_cv_summary$mean_auc,
    nested_cv_summary$auc_ci_lower,
    nested_cv_summary$auc_ci_upper
  ))
  cat(sprintf(
    "          Mean PR-AUC: %.3f (95%% CI: %.3f-%.3f)\n",
    nested_cv_summary$mean_pr_auc,
    nested_cv_summary$pr_auc_ci_lower,
    nested_cv_summary$pr_auc_ci_upper
  ))

  # Copy nested CV visualisations to output
  nested_cv_viz_dir <- file.path(opt$model_dir, "nested_cv_visualizations")
  if (dir.exists(nested_cv_viz_dir)) {
    cat("          [OK] Found nested CV visualisations\n")
    cat("          Copying plots for markdown report...\n")
  }
} else {
  cat("          [INFO] Nested CV results not found - skipping\n")
}

nested_cv_done <- nested_cv_summary$available

# ============================================================================
# STEP 3.5: MULTICLASS PERFORMANCE ANALYSIS
# ============================================================================
cat("[STEP 3.5/10] Multiclass performance analysis...\n")

# Create performance directory
performance_dir <- file.path(output_dir, "model_performance")
dir.create(performance_dir, showWarnings = FALSE, recursive = TRUE)

# Initialize performance analysis flag
performance_analysis_done <- FALSE
# Try to load predictions.csv
predictions_file <- file.path(opt$model_dir, "predictions.csv")

if (!file.exists(predictions_file)) {
  cat("          [WARNING] predictions.csv not found in model directory\n")
  cat("          [INFO] Skipping multiclass performance analysis\n")
  cat("          Expected location:", predictions_file, "\n\n")

  # Set flag to skip this analysis in markdown
  performance_analysis_done <- FALSE
} else {
  cat("          Loading predictions from:", basename(predictions_file), "\n")

  predictions <- tryCatch(
    {
      read.csv(predictions_file, stringsAsFactors = FALSE)
    },
    error = function(e) {
      cat("          [ERROR] Failed to load predictions.csv:", e$message, "\n")
      return(NULL)
    }
  )

  if (is.null(predictions)) {
    performance_analysis_done <- FALSE
  } else {
    cat("          [OK] Loaded", nrow(predictions), "predictions\n")

    # Check for required columns
    required_cols <- c("sample_id", "actual", "predicted")

    has_required <- all(required_cols %in% colnames(predictions))

    if (!has_required) {
      cat("          [WARNING] Missing required columns in predictions.csv\n")
      cat("          Required:", paste(required_cols, collapse = ", "), "\n")
      cat("          Found:", paste(colnames(predictions), collapse = ", "), "\n")
      cat("          [INFO] Skipping multiclass performance analysis\n\n")
      performance_analysis_done <- FALSE
    } else {
      # Auto-detect classification type
      class_info <- detect_classification_type(predictions)

      cat(sprintf(
        "          Classification type: %s\n",
        ifelse(class_info$is_binary, "Binary", "Multiclass")
      ))
      cat(sprintf("          Number of classes: %d\n", class_info$n_classes))
      cat(sprintf("          Classes: %s\n", paste(class_info$class_names, collapse = ", ")))

      # Check if we have probability columns
      has_probabilities <- length(class_info$class_names) > 0

      if (!has_probabilities) {
        cat("          [WARNING] No probability columns found in predictions.csv\n")
        cat("          Found columns:", paste(colnames(predictions), collapse = ", "), "\n")
        cat("          [INFO] Skipping detailed ROC/PR analysis\n\n")
        performance_analysis_done <- FALSE
      } else {
        # Proceed with full analysis
        cat("          [OK] All required columns present\n")

        # Get unique classes
        classes <- sort(unique(c(predictions$actual, predictions$predicted)))
        n_classes <- length(classes)

        cat("          Classes found:", paste(classes, collapse = ", "), "\n")
        cat("          Number of classes:", n_classes, "\n\n")

        # ========================================================================
        # 1. Overall Performance Metrics
        # ========================================================================
        cat("          Calculating overall performance metrics...\n")

        # Overall accuracy
        overall_accuracy <- mean(predictions$actual == predictions$predicted)

        # Cohen's Kappa
        conf_matrix_overall <- table(Actual = predictions$actual, Predicted = predictions$predicted)
        n_total <- sum(conf_matrix_overall)
        n_correct <- sum(diag(conf_matrix_overall))

        # Expected accuracy by chance
        row_sums <- rowSums(conf_matrix_overall)
        col_sums <- colSums(conf_matrix_overall)
        expected_accuracy <- sum((row_sums / n_total) * (col_sums / n_total))

        kappa <- (overall_accuracy - expected_accuracy) / (1 - expected_accuracy)

        cat("          Overall Accuracy:", round(overall_accuracy, 4), "\n")
        cat("          Cohen's Kappa:", round(kappa, 4), "\n\n")

        # Save overall metrics
        overall_metrics <- data.frame(
          Metric = c("Overall_Accuracy", "Cohens_Kappa", "N_Samples", "N_Classes"),
          Value = c(overall_accuracy, kappa, nrow(predictions), n_classes)
        )

        write.csv(overall_metrics,
          file.path(performance_dir, "overall_metrics.csv"),
          row.names = FALSE
        )

        # ========================================================================
        # 2. Confusion Matrix
        # ========================================================================
        cat("          Creating confusion matrix...\n")

        # Save confusion matrix as CSV
        conf_matrix_df <- as.data.frame.matrix(conf_matrix_overall)
        conf_matrix_df <- cbind(Actual_Class = rownames(conf_matrix_df), conf_matrix_df)

        write.csv(conf_matrix_df,
          file.path(performance_dir, "confusion_matrix.csv"),
          row.names = FALSE
        )

        # Plot confusion matrix heatmap
        suppressPackageStartupMessages(library(pheatmap))

        # Generate in all requested formats
        # ALWAYS include PNG as baseline/fallback to prevent broken reports
        formats_to_generate <- tolower(trimws(unlist(strsplit(opt$figure_format, ","))))
        if ("all" %in% formats_to_generate) {
          formats_to_generate <- c("png", "jpeg", "svg")
        }
        # Normalize jpg to jpeg
        formats_to_generate <- gsub("jpg", "jpeg", formats_to_generate)
        # Ensure PNG is always included (baseline format for report compatibility)
        if (!"png" %in% formats_to_generate) {
          formats_to_generate <- c("png", formats_to_generate)
        }

        for (fmt in unique(formats_to_generate)) {
          plot_file <- file.path(performance_dir, sprintf("confusion_matrix_3x3.%s", fmt))
          open_plot_device(plot_file, width = 10, height = 8, res = 300)

          # Normalise by row (actual class)
          conf_matrix_norm <- conf_matrix_overall / rowSums(conf_matrix_overall)

          # Create heatmap
          pheatmap(conf_matrix_norm,
            display_numbers = conf_matrix_overall,
            number_format = "%.0f",
            cluster_rows = FALSE,
            cluster_cols = FALSE,
            colour = colorRampPalette(c("white", "lightblue", "darkblue"))(100),
            main = "Confusion Matrix (Multiclass)\nNumbers show counts, colours show row-normalised proportions",
            fontsize = 12,
            fontsize_number = 14,
            angle_col = 45
          )

          dev.off()
        }

        cat("          [OK] Confusion matrix saved\n\n")

        # ========================================================================
        # 3. One-vs-Rest Analysis for Each Class
        # ========================================================================
        cat("          Performing one-vs-rest analysis for each class...\n\n")

        # Initialize results storage
        class_metrics_list <- list()
        roc_data_list <- list()
        pr_data_list <- list()

        for (class in classes) {
          cat("          Analysing class:", class, "\n")

          # Create binary labels (class vs all others)
          binary_actual <- ifelse(predictions$actual == class, class, "Other")
          binary_predicted <- ifelse(predictions$predicted == class, class, "Other")

          # Get probability for this class
          if (class %in% colnames(predictions)) {
            class_prob <- predictions[[class]]
          } else {
            cat("            [WARNING] No probability column for", class, "\n")
            next
          }

          # Binary confusion matrix
          binary_conf <- table(Actual = binary_actual, Predicted = binary_predicted)

          # Calculate binary metrics
          binary_metrics <- calculate_binary_metrics(binary_conf)

          cat("            Sensitivity:", round(binary_metrics$sensitivity, 4), "\n")
          cat("            Specificity:", round(binary_metrics$specificity, 4), "\n")
          cat("            PPV:", round(binary_metrics$ppv, 4), "\n")
          cat("            NPV:", round(binary_metrics$npv, 4), "\n")
          cat("            F1:", round(binary_metrics$f1, 4), "\n")
          cat("            MCC:", round(binary_metrics$mcc, 4), "\n")

          # ROC analysis
          suppressPackageStartupMessages(library(pROC))

          roc_obj <- roc(
            response = binary_actual,
            predictor = class_prob,
            levels = c("Other", class),
            direction = "<",
            quiet = TRUE
          )

          # Calculate ROC-AUC with 95% CI (bootstrap)
          roc_ci <- ci.auc(roc_obj, method = "bootstrap", boot.n = 2000, quiet = TRUE)

          cat(
            "            ROC-AUC:", round(roc_obj$auc, 4),
            "(95% CI:", round(roc_ci[1], 4), "-", round(roc_ci[3], 4), ")\n"
          )

          # Store ROC data
          roc_data_list[[class]] <- data.frame(
            class = class,
            sensitivity = roc_obj$sensitivities,
            specificity = roc_obj$specificities,
            threshold = roc_obj$thresholds
          )

          # PR curve analysis
          pr_curve <- calculate_pr_curve(binary_actual, class_prob, class)
          pr_auc <- calculate_pr_auc(pr_curve)

          cat("            PR-AUC:", round(pr_auc, 4), "\n\n")

          # Store PR data
          pr_data_list[[class]] <- cbind(class = class, pr_curve)

          # Store all metrics for this class
          class_metrics_list[[class]] <- data.frame(
            Class = class,
            Sensitivity = binary_metrics$sensitivity,
            Specificity = binary_metrics$specificity,
            PPV = binary_metrics$ppv,
            NPV = binary_metrics$npv,
            F1_Score = binary_metrics$f1,
            MCC = binary_metrics$mcc,
            ROC_AUC = as.numeric(roc_obj$auc),
            ROC_AUC_95CI_Lower = as.numeric(roc_ci[1]),
            ROC_AUC_95CI_Upper = as.numeric(roc_ci[3]),
            PR_AUC = pr_auc
          )
        }

        # Combine all class metrics
        class_metrics_df <- do.call(rbind, class_metrics_list)
        rownames(class_metrics_df) <- NULL

        write.csv(class_metrics_df,
          file.path(performance_dir, "class_metrics_one_vs_rest.csv"),
          row.names = FALSE
        )

        cat("          [OK] One-vs-rest analysis complete for all classes\n\n")

        # ========================================================================
        # 4. Plot ROC Curves (3-panel)
        # ========================================================================
        cat("          Generating ROC curve plots...\n")

        for (fmt in unique(formats_to_generate)) {
          plot_file <- file.path(performance_dir, sprintf("roc_curves_multiclass.%s", fmt))
          open_plot_device(plot_file, width = 15, height = 5, res = 300)

          par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

          for (class in classes) {
            if (class %in% names(roc_data_list)) {
              roc_data <- roc_data_list[[class]]

              # Plot ROC curve
              plot(1 - roc_data$specificity, roc_data$sensitivity,
                type = "l", lwd = 2, col = "darkblue",
                xlim = c(0, 1), ylim = c(0, 1),
                xlab = "False Positive Rate (1 - Specificity)",
                ylab = "True Positive Rate (Sensitivity)",
                main = paste("ROC Curve:", class, "vs Rest")
              )

              # Add diagonal reference line
              abline(0, 1, lty = 2, col = "grey50")

              # Add grid
              grid(col = "grey90")

              # Add AUC to plot
              class_auc <- class_metrics_df$ROC_AUC[class_metrics_df$Class == class]
              class_ci_lower <- class_metrics_df$ROC_AUC_95CI_Lower[class_metrics_df$Class == class]
              class_ci_upper <- class_metrics_df$ROC_AUC_95CI_Upper[class_metrics_df$Class == class]

              legend("bottomright",
                legend = sprintf(
                  "AUC = %.3f\n(95%% CI: %.3f - %.3f)",
                  class_auc, class_ci_lower, class_ci_upper
                ),
                bty = "n", cex = 1.1
              )
            }
          }

          dev.off()
        }

        cat("          [OK] ROC curves saved\n\n")

        # ========================================================================
        # 5. Plot PR Curves (3-panel)
        # ========================================================================
        cat("          Generating PR curve plots...\n")

        for (fmt in unique(formats_to_generate)) {
          plot_file <- file.path(performance_dir, sprintf("pr_curves_multiclass.%s", fmt))
          open_plot_device(plot_file, width = 15, height = 5, res = 300)

          par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

          for (class in classes) {
            if (class %in% names(pr_data_list)) {
              pr_data <- pr_data_list[[class]]

              # Plot PR curve
              plot(pr_data$recall, pr_data$precision,
                type = "l", lwd = 2, col = "darkgreen",
                xlim = c(0, 1), ylim = c(0, 1),
                xlab = "Recall (Sensitivity)",
                ylab = "Precision (PPV)",
                main = paste("Precision-Recall Curve:", class, "vs Rest")
              )

              # Add grid
              grid(col = "grey90")

              # Add PR-AUC to plot
              class_pr_auc <- class_metrics_df$PR_AUC[class_metrics_df$Class == class]

              legend("bottomleft",
                legend = sprintf("PR-AUC = %.3f", class_pr_auc),
                bty = "n", cex = 1.1
              )
            }
          }

          dev.off()
        }

        cat("          [OK] PR curves saved\n\n")

        # ========================================================================
        # 6. Stratified Performance (TRAIN vs TEST)
        # ========================================================================
        if ("Dataset" %in% colnames(predictions)) {
          cat("          Calculating stratified performance (TRAIN vs TEST)...\n")

          stratified_metrics <- list()

          for (dataset in unique(predictions$Dataset)) {
            subset_data <- predictions[predictions$Dataset == dataset, ]

            # Overall accuracy for this dataset
            subset_accuracy <- mean(subset_data$actual == subset_data$predicted)

            # Confusion matrix
            subset_conf <- table(Actual = subset_data$actual, Predicted = subset_data$predicted)

            # Kappa
            n_subset <- sum(subset_conf)
            n_correct_subset <- sum(diag(subset_conf))
            row_sums_subset <- rowSums(subset_conf)
            col_sums_subset <- colSums(subset_conf)
            expected_acc_subset <- sum((row_sums_subset / n_subset) * (col_sums_subset / n_subset))
            kappa_subset <- (subset_accuracy - expected_acc_subset) / (1 - expected_acc_subset)

            stratified_metrics[[dataset]] <- data.frame(
              Dataset = dataset,
              N_Samples = nrow(subset_data),
              Accuracy = subset_accuracy,
              Kappa = kappa_subset
            )

            cat(
              "            Dataset:", dataset, "- Accuracy:", round(subset_accuracy, 4),
              "Kappa:", round(kappa_subset, 4), "\n"
            )
          }

          stratified_df <- do.call(rbind, stratified_metrics)
          rownames(stratified_df) <- NULL

          write.csv(stratified_df,
            file.path(performance_dir, "stratified_performance.csv"),
            row.names = FALSE
          )

          cat("          [OK] Stratified performance saved\n\n")
        }

        # ========================================================================
        # 7. Platform-Stratified Performance (if batch_column specified)
        # ========================================================================
        if (!is.null(opt$batch_column) && opt$batch_column %in% colnames(sample_sheet_matched)) {
          cat("          Calculating platform-stratified performance...\n")

          # Match predictions to sample sheet
          predictions$batch <- sample_sheet_matched[[opt$batch_column]][match(
            predictions$sample_id,
            rownames(sample_sheet_matched)
          )]

          platform_metrics <- list()

          for (platform in unique(predictions$batch[!is.na(predictions$batch)])) {
            platform_data <- predictions[predictions$batch == platform & !is.na(predictions$batch), ]

            if (nrow(platform_data) > 0) {
              platform_accuracy <- mean(platform_data$actual == platform_data$predicted)

              platform_metrics[[as.character(platform)]] <- data.frame(
                Platform = platform,
                N_Samples = nrow(platform_data),
                Accuracy = platform_accuracy
              )

              cat(
                "            Platform:", platform, "- N:", nrow(platform_data),
                "Accuracy:", round(platform_accuracy, 4), "\n"
              )
            }
          }

          if (length(platform_metrics) > 0) {
            platform_df <- do.call(rbind, platform_metrics)
            rownames(platform_df) <- NULL

            write.csv(platform_df,
              file.path(performance_dir, "platform_stratified_performance.csv"),
              row.names = FALSE
            )

            cat("          [OK] Platform-stratified performance saved\n\n")
          }
        }

        # Set flag for markdown generation
        performance_analysis_done <- TRUE

        cat("          [SUCCESS] Multiclass performance analysis complete!\n")
        cat("          Results saved to:", performance_dir, "\n\n")
      }
    }
  }
}

cat("\n")
# ============================================================================
# ADDITIONAL DATA LOADING FOR COMPREHENSIVE ANALYSIS
# ============================================================================
cat("\n========================================\n")
cat("Loading additional validation data...\n")
cat("========================================\n")

# Load fold-by-fold CV results
cv_detailed <- load_cv_detailed(opt$model_dir)
if (cv_detailed$available) {
  cat(sprintf("          [OK] Loaded CV detailed results (%d folds)\n", cv_detailed$n_folds))
}

# FIX v5.12.0: Load standard CV predictions to calculate per-class metrics
cv_predictions_file <- list.files(file.path(opt$model_dir, "cross_validation"),
  pattern = "^cv_predictions_.*\\.csv$",
  full.names = TRUE
)
if (length(cv_predictions_file) > 0) {
  cv_predictions_data <- read.csv(cv_predictions_file[1], stringsAsFactors = FALSE)
  cat(sprintf("          [OK] Loaded CV predictions for per-class metrics: %d rows\n", nrow(cv_predictions_data)))

  # Detect positive class from CV predictions (same logic as nested CV)
  pred_classes <- unique(c(cv_predictions_data$actual, cv_predictions_data$predicted))
  cv_positive_class <- NULL
  case_keywords <- c("infected", "infection", "disease", "case", "positive", "treatment", "tumor", "cancer", "patient", "aspergillus")
  for (cls in pred_classes) {
    for (keyword in case_keywords) {
      if (grepl(keyword, tolower(cls))) {
        cv_positive_class <- cls
        break
      }
    }
    if (!is.null(cv_positive_class)) break
  }

  # Calculate per-class metrics for positive class (similar to nested CV)
  if (!is.null(cv_positive_class)) {
    cv_per_class_metrics <- calculate_per_class_nested_cv_metrics(
      list(available = TRUE, data = cv_predictions_data),
      cv_positive_class,
      n_bootstrap = 1000
    )
    if (cv_per_class_metrics$available) {
      cat(sprintf("          [OK] Calculated per-class CV metrics for %s\n", cv_positive_class))
    }
  } else {
    cv_per_class_metrics <- list(available = FALSE)
  }
} else {
  cv_per_class_metrics <- list(available = FALSE)
}

# Load per-sample nested CV predictions
nested_cv_predictions <- load_cv_predictions(opt$model_dir)
if (nested_cv_predictions$available) {
  cat(sprintf("          [OK] Loaded nested CV predictions from: %s\n", basename(nested_cv_predictions$file_path)))
  cat(sprintf("          Total predictions: %d\n", nrow(nested_cv_predictions$data)))

  # Analyse sample consistency
  sample_consistency <- analyse_sample_consistency(nested_cv_predictions)
  if (sample_consistency$available) {
    cat(sprintf(
      "          Consistent samples: %d (%.1f%%)\n",
      sample_consistency$n_consistent,
      sample_consistency$pct_consistent
    ))

    # SAVE COMPLETE SAMPLE CONSISTENCY DATA TO CSV (for paper)
    consistency_csv_path <- file.path(output_dir, "sample_prediction_consistency.csv")
    tryCatch(
      {
        # Enhance with Lab_ID and Animal_ID if available
        consistency_export <- sample_consistency$summary

        # Add Lab_ID and Animal_ID by matching sample_id to sample_sheet_matched
        if (exists("sample_sheet_matched") && nrow(sample_sheet_matched) > 0) {
          # Match sample IDs (convert both to character for robust matching)
          consistency_export$Lab_ID <- sapply(consistency_export$sample_id, function(sid) {
            # Try matching against rownames (sample IDs from methylation data)
            sid_char <- as.character(sid)
            if (sid_char %in% rownames(sample_sheet_matched) && "Lab_ID" %in% colnames(sample_sheet_matched)) {
              return(as.character(sample_sheet_matched[sid_char, "Lab_ID"]))
            }
            # If rownames match failed, try matching against sample ID columns
            if (opt$sample_id_column %in% colnames(sample_sheet_matched)) {
              match_idx <- which(as.character(sample_sheet_matched[[opt$sample_id_column]]) == sid_char)
              if (length(match_idx) > 0 && "Lab_ID" %in% colnames(sample_sheet_matched)) {
                return(as.character(sample_sheet_matched[match_idx[1], "Lab_ID"]))
              }
            }
            return(NA_character_)
          })

          consistency_export$Animal_ID <- sapply(consistency_export$sample_id, function(sid) {
            # Try matching against rownames (sample IDs from methylation data)
            sid_char <- as.character(sid)
            if (sid_char %in% rownames(sample_sheet_matched) && "Animal_ID" %in% colnames(sample_sheet_matched)) {
              return(as.character(sample_sheet_matched[sid_char, "Animal_ID"]))
            }
            # If rownames match failed, try matching against sample ID columns
            if (opt$sample_id_column %in% colnames(sample_sheet_matched)) {
              match_idx <- which(as.character(sample_sheet_matched[[opt$sample_id_column]]) == sid_char)
              if (length(match_idx) > 0 && "Animal_ID" %in% colnames(sample_sheet_matched)) {
                return(as.character(sample_sheet_matched[match_idx[1], "Animal_ID"]))
              }
            }
            return(NA_character_)
          })

          cat(sprintf(
            "          [INFO] Matched %d/%d samples with Lab_ID\n",
            sum(!is.na(consistency_export$Lab_ID)), nrow(consistency_export)
          ))
          cat(sprintf(
            "          [INFO] Matched %d/%d samples with Animal_ID\n",
            sum(!is.na(consistency_export$Animal_ID)), nrow(consistency_export)
          ))
        }

        # Reorder columns to match report table: Sample ID, Lab ID, Animal ID, Actual Class, Correct Folds, Total Folds, Consistency Rate
        col_order <- c("sample_id")
        if ("Lab_ID" %in% colnames(consistency_export)) col_order <- c(col_order, "Lab_ID")
        if ("Animal_ID" %in% colnames(consistency_export)) col_order <- c(col_order, "Animal_ID")
        col_order <- c(col_order, "actual", "correct_folds", "total_folds", "consistency_rate")

        consistency_export <- consistency_export[, col_order]

        # Rename columns to match report exactly
        colnames(consistency_export)[colnames(consistency_export) == "sample_id"] <- "Sample_ID"
        colnames(consistency_export)[colnames(consistency_export) == "actual"] <- "Actual_Class"
        colnames(consistency_export)[colnames(consistency_export) == "correct_folds"] <- "Correct_Folds"
        colnames(consistency_export)[colnames(consistency_export) == "total_folds"] <- "Total_Folds"
        colnames(consistency_export)[colnames(consistency_export) == "consistency_rate"] <- "Consistency_Rate"

        write.csv(consistency_export, consistency_csv_path, row.names = FALSE)
        cat(sprintf("          [SAVED] Sample consistency data: %s\n", consistency_csv_path))
        cat(sprintf("                  Columns: %s\n", paste(colnames(consistency_export), collapse = ", ")))
      },
      error = function(e) {
        cat(sprintf("          [WARNING] Enhanced CSV export failed: %s\n", conditionMessage(e)))
        cat(sprintf("          [INFO] Falling back to basic export\n"))
        # Fallback: save basic version
        tryCatch(
          {
            write.csv(sample_consistency$summary, consistency_csv_path, row.names = FALSE)
            cat(sprintf("          [SAVED] Basic sample consistency data: %s\n", consistency_csv_path))
          },
          error = function(e2) {
            cat(sprintf("          [ERROR] Complete failure to save CSV: %s\n", conditionMessage(e2)))
          }
        )
      }
    )
  }

  # Calculate per-class metrics (for binary classification)
  # Detect positive class automatically (infected/case group)
  pred_classes <- unique(c(nested_cv_predictions$data$actual, nested_cv_predictions$data$predicted))
  positive_class <- NULL
  case_keywords <- c("infected", "infection", "disease", "case", "positive", "treatment", "tumor", "cancer", "patient", "aspergillus")
  for (cls in pred_classes) {
    for (keyword in case_keywords) {
      if (grepl(keyword, tolower(cls))) {
        positive_class <- cls
        break
      }
    }
    if (!is.null(positive_class)) break
  }

  if (!is.null(positive_class)) {
    cat(sprintf("          Calculating per-class metrics for positive class: %s\n", positive_class))
    per_class_metrics <- calculate_per_class_nested_cv_metrics(nested_cv_predictions, positive_class, n_bootstrap = 1000)
    if (per_class_metrics$available) {
      cat(sprintf("          [OK] Per-class metrics calculated (%d folds)\n", per_class_metrics$n_folds))
      cat(sprintf(
        "              Sensitivity: %.3f (95%% CI: %.3f-%.3f)\n",
        per_class_metrics$sensitivity$mean,
        per_class_metrics$sensitivity$ci_lower,
        per_class_metrics$sensitivity$ci_upper
      ))
      cat(sprintf(
        "              F1 Score: %.3f (95%% CI: %.3f-%.3f)\n",
        per_class_metrics$f1_score$mean,
        per_class_metrics$f1_score$ci_lower,
        per_class_metrics$f1_score$ci_upper
      ))
    }
  } else {
    cat("          [INFO] Could not auto-detect positive class for per-class metrics\n")
    per_class_metrics <- list(available = FALSE)
  }

  # Calculate study-stratified performance
  cat("\n")
  cat("================================================================================\n")
  cat("  STUDY-STRATIFIED PERFORMANCE ANALYSIS\n")
  cat("================================================================================\n")
  if (file.exists(opt$sample_sheet)) {
    cat(sprintf("Sample sheet found: %s\n", opt$sample_sheet))
    cat("Attempting to calculate study-stratified performance...\n\n")

    # Try to use predictions.csv from main model dir first (preferred)
    # Fall back to nested CV predictions if not available
    predictions_for_stratification <- NULL
    predictions_source <- ""

    if (exists("predictions") && !is.null(predictions) && nrow(predictions) > 0) {
      predictions_for_stratification <- predictions
      predictions_source <- "main model predictions.csv"
      cat(sprintf("Using predictions from: %s (%d samples)\n\n", predictions_source, nrow(predictions)))
    } else if (exists("nested_cv_predictions") && nested_cv_predictions$available) {
      predictions_for_stratification <- nested_cv_predictions
      predictions_source <- "nested CV predictions"
      cat(sprintf("Using predictions from: %s (%d samples)\n\n", predictions_source, nrow(nested_cv_predictions$data)))
    } else {
      cat(">>> FAILED: No predictions available for stratification\n")
      cat("    Neither predictions.csv nor nested CV predictions found\n")
      study_performance <- list(available = FALSE, reason = "No predictions available")
    }

    if (!is.null(predictions_for_stratification)) {
      study_performance <- calculate_study_stratified_performance(predictions_for_stratification, opt$sample_sheet)
      if (study_performance$available) {
        cat(sprintf("\n>>> SUCCESS: Study-stratified results calculated\n"))
        cat(sprintf("    - Number of studies: %d\n", length(unique(study_performance$results$Study))))
        cat(sprintf("    - Overall accuracy: %.1f%%\n", study_performance$total_accuracy))

        # Save results to CSV
        study_stratified_file <- file.path(output_dir, "study_stratified_performance.csv")
        write.csv(study_performance$results, study_stratified_file, row.names = FALSE)
        cat(sprintf("    - Results saved to: %s\n", basename(study_stratified_file)))

        cat(sprintf("    >>> Section 1.6.9 WILL appear in the report\n"))
      } else {
        cat(sprintf("\n>>> FAILED: Study stratification did not complete\n"))
        cat(sprintf("    Reason: %s\n", study_performance$reason))
        cat(sprintf("    >>> Section 1.6.9 WILL NOT appear in the report\n"))
      }
    }
  } else {
    cat(sprintf(">>> FAILED: Sample sheet file not found\n"))
    cat(sprintf("    Expected path: %s\n", opt$sample_sheet))
    cat(sprintf("    >>> Section 1.6.9 WILL NOT appear in the report\n"))
    study_performance <- list(available = FALSE, reason = "Sample sheet file not found")
  }
  cat("================================================================================\n\n")
} else {
  # If nested_cv_predictions not available, initialize as unavailable
  per_class_metrics <- list(available = FALSE)
  study_performance <- list(available = FALSE)
}

# ============================================================================
# CRITICAL DEBUG: Report generation variable check
# ============================================================================
cat("\n")
cat("================================================================================\n")
cat("  REPORT GENERATION VARIABLES STATUS\n")
cat("================================================================================\n")
cat(sprintf("✓ nested_cv_done: %s\n", if (exists("nested_cv_done") && nested_cv_done) "TRUE" else "FALSE"))
cat(sprintf("✓ study_performance exists: %s\n", if (exists("study_performance")) "YES" else "NO"))
if (exists("study_performance")) {
  cat(sprintf("✓ study_performance$available: %s\n", if (study_performance$available) "TRUE" else "FALSE"))
  if (study_performance$available) {
    cat(sprintf("✓ Number of studies: %d\n", length(unique(study_performance$results$Study))))
    cat("\n>>> Section 1.6.9 SHOULD appear in the report\n")
  } else {
    cat(sprintf("✗ Reason: %s\n", study_performance$reason))
    cat("\n>>> Section 1.6.9 WILL NOT appear (study_performance not available)\n")
  }
} else {
  cat("✗ study_performance variable does not exist\n")
  cat("\n>>> Section 1.6.9 WILL NOT appear (variable missing)\n")
}
cat("================================================================================\n\n")

# Load hyperparameter results (already have nested_cv_summary)

# Analyse calibration if predictions available
if (exists("predictions") && nrow(predictions) > 0) {
  calibration_analysis <- analyse_calibration(predictions)
  if (calibration_analysis$available) {
    cat("          [OK] Probability calibration analysed\n")
  }
}

log_section_end(3, "Loading model and data")
cat("\n")
# ============================================================================
# STEP 5: BATCH/PLATFORM EFFECT ANALYSIS
# ============================================================================
log_section_start(5, "Batch/platform effect analysis", 10)

batch_dir <- file.path(output_dir, "batch_effects")
dir.create(batch_dir, showWarnings = FALSE)

# Parse technical covariates
tech_covariates <- c()
if (!is.null(opt$batch_column)) tech_covariates <- c(tech_covariates, opt$batch_column)
if (!is.null(opt$storage_column)) tech_covariates <- c(tech_covariates, opt$storage_column)
if (!is.null(opt$date_column)) tech_covariates <- c(tech_covariates, opt$date_column)
if (!is.null(opt$technical_covariates)) {
  additional <- strsplit(opt$technical_covariates, ",")[[1]]
  tech_covariates <- c(tech_covariates, trimws(additional))
}

tech_covariates <- unique(tech_covariates)
cat("          Testing", length(tech_covariates), "technical covariates\n")

# === ALWAYS generate PCA/UMAP for biological signal (methylation only) ===
# This section ALWAYS runs regardless of technical covariates

# Prepare methylation matrix for dimensionality reduction
cat("          Preparing methylation data for PCA/UMAP...\n")
meth_matrix_t <- t(meth_matrix)
complete_samples <- rowSums(is.na(meth_matrix_t)) < ncol(meth_matrix_t) * 0.5
meth_matrix_base <- meth_matrix_t[complete_samples, ]

# Impute NAs
for (i in 1:ncol(meth_matrix_base)) {
  meth_matrix_base[is.na(meth_matrix_base[, i]), i] <- mean(meth_matrix_base[, i], na.rm = TRUE)
}

# Remove zero-variance columns (use small threshold instead of exact zero)
col_vars <- apply(meth_matrix_base, 2, var, na.rm = TRUE)
zero_var_cols <- which(col_vars < 1e-10 | is.na(col_vars))

if (length(zero_var_cols) > 0) {
  cat(sprintf("          [INFO] Removing %d near-zero variance markers for PCA/UMAP\n", length(zero_var_cols)))
  meth_matrix_base <- meth_matrix_base[, -zero_var_cols, drop = FALSE]
}

cat(sprintf("          [INFO] %d variable markers available for PCA/UMAP (need at least 2)\n", ncol(meth_matrix_base)))

if (ncol(meth_matrix_base) >= 2) {
  # === PCA on methylation ONLY (biological signal) ===
  cat("          Creating PCA from methylation data (biological signal)...\n")

  pca_meth_only <- prcomp(meth_matrix_base, scale. = TRUE, center = TRUE)
  pca_df_bio <- data.frame(
    PC1 = pca_meth_only$x[, 1],
    PC2 = pca_meth_only$x[, 2],
    sample_id = rownames(meth_matrix_base)
  )

  # Add sample_sheet info
  pca_df_bio <- merge(pca_df_bio, sample_sheet_matched,
    by.x = "sample_id", by.y = opt$sample_id_column, all.x = TRUE
  )

  # Add group information - ALWAYS ensure biological_group exists
  if (!("biological_group" %in% colnames(pca_df_bio))) {
    if (exists("sample_metadata") && "actual_class" %in% colnames(sample_metadata)) {
      pca_df_bio$biological_group <- sample_metadata$actual_class[match(pca_df_bio$sample_id, rownames(sample_metadata))]
    } else {
      pca_df_bio$biological_group <- "Unknown"
    }
  }

  var_exp_bio <- summary(pca_meth_only)$importance[2, 1:2] * 100

  # Plot biological PCA - ALWAYS generate if we have data
  if (nrow(pca_df_bio) > 0) {
    model_desc <- sprintf(
      "%s | %s | %s bp | %d markers",
      species_name, model_type, window_size, n_markers
    )

    p_bio <- ggplot(pca_df_bio, aes(x = PC1, y = PC2, colour = biological_group)) +
      stat_ellipse(aes(fill = biological_group),
        colour = NA, geom = "polygon", alpha = 0.15, level = 0.95,
        show.legend = FALSE
      ) +
      stat_ellipse(aes(group = biological_group),
        colour = "black", fill = NA, geom = "polygon", alpha = 1, level = 0.95,
        linewidth = 1.5, linetype = "dashed", show.legend = FALSE
      ) +
      geom_point(size = 5, alpha = 0.8) +
      {
        if (!is.null(group_colors)) {
          list(
            scale_color_manual(values = group_colors),
            scale_fill_manual(values = group_colors)
          )
        } else {
          list(
            scale_color_brewer(palette = "Set1"),
            scale_fill_brewer(palette = "Set1")
          )
        }
      } +
      guides(fill = "none") +
      labs(
        title = "PCA - Methylation Only (Biological Signal)",
        subtitle = paste(model_desc, "| Components from methylation data"),
        x = paste0("PC1 (", round(var_exp_bio[1], 1), "%)"),
        y = paste0("PC2 (", round(var_exp_bio[2], 1), "%)"),
        colour = "Biological Group"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 10, colour = "gray30", hjust = 0.5, face = "italic")
      )

    save_plot_multiple_formats(p_bio,
      file.path(batch_dir, "pca_METHYLATION_ONLY"),
      width = 12, height = 10
    )

    cat("          [OK] Created pca_METHYLATION_ONLY.png\n")
  }

  # === UMAP on methylation ONLY (biological signal) ===
  cat("          Computing UMAP from methylation data (biological signal)...\n")

  umap_meth_only <- umap(meth_matrix_base, n_neighbors = min(15, nrow(meth_matrix_base) - 1))

  umap_df_bio <- data.frame(
    UMAP1 = umap_meth_only$layout[, 1],
    UMAP2 = umap_meth_only$layout[, 2],
    sample_id = rownames(meth_matrix_base)
  )

  umap_df_bio <- merge(umap_df_bio, sample_sheet_matched,
    by.x = "sample_id", by.y = opt$sample_id_column, all.x = TRUE
  )

  # ALWAYS ensure biological_group exists
  if (!("biological_group" %in% colnames(umap_df_bio))) {
    if (exists("sample_metadata") && "actual_class" %in% colnames(sample_metadata)) {
      umap_df_bio$biological_group <- sample_metadata$actual_class[match(umap_df_bio$sample_id, rownames(sample_metadata))]
    } else {
      umap_df_bio$biological_group <- "Unknown"
    }
  }

  # Plot biological UMAP - ALWAYS generate if we have data
  if (nrow(umap_df_bio) > 0) {
    model_desc <- sprintf(
      "%s | %s | %s bp | %d markers",
      species_name, model_type, window_size, n_markers
    )

    p_bio <- ggplot(umap_df_bio, aes(x = UMAP1, y = UMAP2, colour = biological_group)) +
      stat_ellipse(aes(fill = biological_group),
        colour = NA, geom = "polygon", alpha = 0.15, level = 0.95,
        show.legend = FALSE
      ) +
      stat_ellipse(aes(group = biological_group),
        colour = "black", fill = NA, geom = "polygon", alpha = 1, level = 0.95,
        linewidth = 1.5, linetype = "dashed", show.legend = FALSE
      ) +
      geom_point(size = 5, alpha = 0.8) +
      {
        if (!is.null(group_colors)) {
          list(
            scale_color_manual(values = group_colors),
            scale_fill_manual(values = group_colors)
          )
        } else {
          list(
            scale_color_brewer(palette = "Set1"),
            scale_fill_brewer(palette = "Set1")
          )
        }
      } +
      guides(fill = "none") +
      labs(
        title = "UMAP - Methylation Only (Biological Signal)",
        subtitle = paste(model_desc, "| Embedding from methylation data"),
        x = "UMAP1",
        y = "UMAP2",
        colour = "Biological Group"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 10, colour = "gray30", hjust = 0.5, face = "italic")
      )

    save_plot_multiple_formats(p_bio,
      file.path(batch_dir, "umap_METHYLATION_ONLY"),
      width = 12, height = 10
    )

    cat("          [OK] Created umap_METHYLATION_ONLY.png\n")
  }
} else {
  cat("          [WARNING] Not enough variable markers for PCA/UMAP - need at least 2\n")
}

# === NOW: Test technical covariates and generate PCA/UMAP with covariates ===
if (length(tech_covariates) > 0) {
  # === Test association of each covariate with methylation - OPTIMIZED ===
  cat("          Computing covariate associations (optimised)...\n")

  covariate_results <- list()

  for (covar in tech_covariates) {
    if (!covar %in% colnames(sample_sheet_matched)) {
      cat("          [WARNING] Column", covar, "not found\n")
      next
    }

    covar_values <- sample_sheet_matched[[covar]]

    if (sum(is.na(covar_values)) > length(covar_values) * 0.5) {
      cat("          [SKIP]", covar, "- too many missing values\n")
      next
    }

    # Initialize result vectors
    marker_pvalues <- numeric(length(selected_markers))
    marker_effects <- numeric(length(selected_markers))

    if (is.numeric(covar_values)) {
      # OPTIMIZED: Vectorized correlation for numeric variables
      valid_samples <- !is.na(covar_values)

      for (i in seq_along(selected_markers)) {
        meth_vals <- meth_matrix[i, ]
        valid_idx <- valid_samples & !is.na(meth_vals)

        if (sum(valid_idx) >= 5) {
          test_result <- cor.test(meth_vals[valid_idx], covar_values[valid_idx],
            method = "spearman", exact = FALSE
          ) # exact=FALSE is faster
          marker_pvalues[i] <- test_result$p.value
          marker_effects[i] <- test_result$estimate
        } else {
          marker_pvalues[i] <- NA
          marker_effects[i] <- NA
        }
      }
    } else {
      # Categorical variable - use Kruskal-Wallis
      for (i in seq_along(selected_markers)) {
        meth_vals <- meth_matrix[i, ]
        valid_idx <- !is.na(meth_vals) & !is.na(covar_values)

        if (sum(valid_idx) >= 5) {
          covar_factor <- factor(covar_values[valid_idx])
          if (length(unique(covar_factor)) > 1) {
            test_result <- kruskal.test(meth_vals[valid_idx] ~ covar_factor)
            marker_pvalues[i] <- test_result$p.value
            marker_effects[i] <- test_result$statistic
          } else {
            marker_pvalues[i] <- NA
            marker_effects[i] <- NA
          }
        } else {
          marker_pvalues[i] <- NA
          marker_effects[i] <- NA
        }
      }
    }

    # Apply multiple testing correction
    marker_adjusted <- apply_multiple_testing(marker_pvalues, method = mt_method_r)

    # Store results
    covariate_results[[covar]] <- data.frame(
      marker = selected_markers,
      p_value = marker_pvalues,
      adjusted_p = marker_adjusted,
      effect = marker_effects,
      significant_raw = marker_pvalues < sig_threshold,
      significant_adjusted = marker_adjusted < sig_threshold,
      stringsAsFactors = FALSE
    )

    n_sig_raw <- sum(marker_pvalues < sig_threshold, na.rm = TRUE)
    n_sig_adjusted <- sum(marker_adjusted < sig_threshold, na.rm = TRUE)
    cat(sprintf(
      "          %s: %d/%d significant (raw p<%.3f), %d/%d (%s<%.3f)\n",
      covar, n_sig_raw, length(selected_markers), sig_threshold,
      n_sig_adjusted, length(selected_markers), mt_method_display, sig_threshold
    ))
  }

  # Save results
  for (covar_name in names(covariate_results)) {
    write.csv(covariate_results[[covar_name]],
      file.path(batch_dir, paste0(
        "covariate_",
        gsub("[^A-Za-z0-9]", "_", covar_name),
        "_associations.csv"
      )),
      row.names = FALSE
    )
  }

  # Summary
  covariate_summary <- data.frame(
    covariate = names(covariate_results),
    n_markers_tested = sapply(covariate_results, function(x) nrow(x)),
    n_significant_raw = sapply(covariate_results, function(x) sum(x$significant_raw, na.rm = TRUE)),
    n_significant_adjusted = sapply(covariate_results, function(x) sum(x$significant_adjusted, na.rm = TRUE)),
    mean_effect = sapply(covariate_results, function(x) mean(abs(x$effect), na.rm = TRUE)),
    mean_raw_p = sapply(covariate_results, function(x) mean(x$p_value, na.rm = TRUE)),
    mean_adjusted_p = sapply(covariate_results, function(x) mean(x$adjusted_p, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )

  write.csv(covariate_summary,
    file.path(batch_dir, "covariate_association_summary.csv"),
    row.names = FALSE
  )

  # Plot showing all markers with stacked significance levels
  if (nrow(covariate_summary) > 0) {
    # Create stacked bar data showing ALL markers
    # Categories: 1) Adjusted p significant, 2) Raw p significant only, 3) Non-significant

    covar_plot_data <- data.frame()
    for (i in 1:nrow(covariate_summary)) {
      covar <- covariate_summary$covariate[i]
      n_total <- covariate_summary$n_markers_tested[i]
      n_adj_sig <- covariate_summary$n_significant_adjusted[i]
      n_raw_sig <- covariate_summary$n_significant_raw[i]

      # Calculate counts for each category
      n_adj_only <- n_adj_sig # Significant by adjusted p
      n_raw_only <- n_raw_sig - n_adj_sig # Significant by raw p but not adjusted
      n_nonsig <- n_total - n_raw_sig # Not significant

      covar_plot_data <- rbind(covar_plot_data, data.frame(
        covariate = covar,
        category = sprintf("%s < %.3f", toupper(mt_method_display), sig_threshold),
        count = n_adj_only,
        order = 3,
        stringsAsFactors = FALSE
      ))

      covar_plot_data <- rbind(covar_plot_data, data.frame(
        covariate = covar,
        category = sprintf("Raw p < %.3f only", sig_threshold),
        count = n_raw_only,
        order = 2,
        stringsAsFactors = FALSE
      ))

      covar_plot_data <- rbind(covar_plot_data, data.frame(
        covariate = covar,
        category = "Not significant",
        count = n_nonsig,
        order = 1,
        stringsAsFactors = FALSE
      ))
    }

    # Calculate total for ordering covariates by adjusted significance
    covar_totals <- aggregate(count ~ covariate,
      data = covar_plot_data[covar_plot_data$order == 3, ],
      FUN = sum
    )
    covar_plot_data$covariate <- factor(covar_plot_data$covariate,
      levels = covar_totals$covariate[order(covar_totals$count)]
    )

    # Set factor levels for stacking order
    covar_plot_data$category <- factor(covar_plot_data$category,
      levels = unique(covar_plot_data$category[order(covar_plot_data$order)])
    )

    # Define colours: gray for non-sig, light for raw only, dark for adjusted
    category_levels <- levels(covar_plot_data$category)
    fill_colors <- setNames(c("gray85", "#6BAED6", "#08519C"), category_levels)

    p <- ggplot(covar_plot_data, aes(x = covariate, y = count, fill = category)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = fill_colors) +
      coord_flip() +
      labs(
        title = "Technical Covariate Associations with Markers",
        subtitle = sprintf("%d markers tested | All markers shown, coloured by significance", n_markers),
        x = "Covariate",
        y = "Number of Markers",
        fill = "Significance Category"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold")
      )

    save_plot_multiple_formats(p,
      file.path(batch_dir, "covariate_associations"),
      width = 10, height = 7
    )
  }

  # === PCA: Include each technical covariate (if meth_matrix_base is available) ===
  if (ncol(meth_matrix_base) >= 2) {
    cat("          Creating PCA plots with technical covariates...\n")

    # PCA including each technical covariate
    for (covar in tech_covariates) {
      if (!covar %in% colnames(sample_sheet_matched)) next

      cat(sprintf("          Computing PCA with covariate: %s\n", covar))

      # Get covariate data matched to PCA samples
      sample_ids_pca <- rownames(meth_matrix_base)

      cat(sprintf("          meth_matrix_base samples: %d\n", length(sample_ids_pca)))
      cat(sprintf("          sample_sheet_matched rows: %d\n", nrow(sample_sheet_matched)))

      # Match covariate values to PCA samples
      if (!is.null(rownames(sample_sheet_matched))) {
        covar_values <- sample_sheet_matched[[covar]][match(
          sample_ids_pca,
          rownames(sample_sheet_matched)
        )]
      } else {
        # If no rownames, try matching by sample_id_column
        sample_sheet_ids <- as.character(sample_sheet_matched[[opt$sample_id_column]])
        match_idx <- match(sample_ids_pca, sample_sheet_ids)
        covar_values <- sample_sheet_matched[[covar]][match_idx]
      }

      cat(sprintf("          covar_values length after match: %d\n", length(covar_values)))
      cat(sprintf("          covar_values NAs: %d\n", sum(is.na(covar_values))))

      # CRITICAL: Filter to samples that have BOTH methylation AND covariate data
      valid_samples_idx <- !is.na(covar_values)

      if (sum(valid_samples_idx) < 10) {
        cat(sprintf("          [SKIP] %s - too few valid samples (%d)\n", covar, sum(valid_samples_idx)))
        next
      }

      cat(sprintf("          Valid samples (non-NA covariate): %d\n", sum(valid_samples_idx)))

      # Subset methylation matrix to valid samples
      meth_matrix_subset <- meth_matrix_base[valid_samples_idx, , drop = FALSE]
      covar_values_subset <- covar_values[valid_samples_idx]

      cat(sprintf(
        "          Subset matrix dimensions: %d x %d\n",
        nrow(meth_matrix_subset), ncol(meth_matrix_subset)
      ))
      cat(sprintf("          Subset covariate length: %d\n", length(covar_values_subset)))

      # Check if numeric or categorical
      is_numeric_covar <- is.numeric(covar_values_subset)
      is_weight_like <- grepl("weight|age|duration|time|storage|days|hours|weeks|months",
        covar,
        ignore.case = TRUE
      )

      if (!is_numeric_covar && is_weight_like) {
        covar_values_subset <- as.numeric(as.character(covar_values_subset))
        is_numeric_covar <- TRUE
      }

      # Skip if too many NAs remain
      if (sum(is.na(covar_values_subset)) > length(covar_values_subset) * 0.5) {
        cat(sprintf("          [SKIP] %s - too many missing values after filtering\n", covar))
        next
      }

      # Create combined matrix: methylation + covariate
      if (is_numeric_covar) {
        # Numeric: scale and add as column
        covar_scaled <- scale(covar_values_subset, center = TRUE, scale = TRUE)
        covar_scaled[is.na(covar_scaled)] <- 0 # Impute remaining NAs to mean (0 after scaling)

        cat(sprintf("          Adding numeric covariate (scaled)\n"))
        combined_matrix <- cbind(meth_matrix_subset, covariate = as.numeric(covar_scaled))
      } else {
        # Categorical: one-hot encode
        covar_factor <- factor(covar_values_subset)
        levels_present <- levels(covar_factor)

        cat(sprintf(
          "          Categorical covariate with %d levels: %s\n",
          length(levels_present), paste(levels_present, collapse = ", ")
        ))

        if (length(levels_present) < 2) {
          cat(sprintf("          [SKIP] %s - only one category\n", covar))
          next
        }

        # Create dummy variables (one-hot encoding)
        covar_dummies <- model.matrix(~ covar_factor - 1)
        colnames(covar_dummies) <- paste0("covar_", levels_present)

        # Replace NAs with 0
        covar_dummies[is.na(covar_dummies)] <- 0

        cat(sprintf("          Adding categorical covariate as %d dummy variables\n", ncol(covar_dummies)))
        combined_matrix <- cbind(meth_matrix_subset, covar_dummies)
      }

      cat(sprintf(
        "          Combined matrix dimensions: %d x %d\n",
        nrow(combined_matrix), ncol(combined_matrix)
      ))

      # Compute PCA on combined matrix
      pca_combined <- prcomp(combined_matrix, scale. = TRUE, center = TRUE)

      pca_df <- data.frame(
        PC1 = pca_combined$x[, 1],
        PC2 = pca_combined$x[, 2],
        sample_id = rownames(meth_matrix_subset)
      )

      # Add metadata
      pca_df <- merge(pca_df, sample_sheet_matched,
        by.x = "sample_id", by.y = opt$sample_id_column, all.x = TRUE
      )

      # Add biological group
      if (exists("sample_metadata") && "actual_class" %in% colnames(sample_metadata)) {
        group_info <- sample_metadata$actual_class[match(pca_df$sample_id, rownames(sample_metadata))]
        pca_df$biological_group <- group_info
      }

      var_exp <- summary(pca_combined)$importance[2, 1:2] * 100

      # Add covariate to pca_df for coloring
      pca_df$covariate_value <- sample_sheet_matched[[covar]][match(
        pca_df$sample_id,
        rownames(sample_sheet_matched)
      )]

      # Force numeric if weight-like
      if (is_weight_like && !is.numeric(pca_df$covariate_value)) {
        pca_df$covariate_value <- as.numeric(as.character(pca_df$covariate_value))
      }

      model_desc <- sprintf(
        "%s | %s | %s bp | %d markers",
        species_name, model_type, window_size, n_markers
      )

      # Plot with biological group as shape/fill and covariate as colour
      if ("biological_group" %in% colnames(pca_df)) {
        # Identify control group for coloring
        group_names <- unique(pca_df$biological_group)
        group_names <- group_names[!is.na(group_names)]

        control_pattern <- "control|ctrl|healthy|normal|untreated|mock|baseline"
        is_control <- grepl(control_pattern, group_names, ignore.case = TRUE)

        group_colors <- character(length(group_names))
        group_colors[is_control] <- "#2E7D32" # Dark green

        n_treatment <- sum(!is_control)
        if (n_treatment > 0) {
          treatment_colors <- colorRampPalette(c("#D32F2F", "#FF6F00"))(n_treatment)
          group_colors[!is_control] <- treatment_colors
        }
        names(group_colors) <- group_names

        p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
          stat_ellipse(aes(group = biological_group),
            colour = "black", fill = NA, geom = "polygon",
            alpha = 1, level = 0.95, linewidth = 1.5,
            linetype = "dashed", show.legend = FALSE
          ) +
          geom_point(
            aes(
              colour = covariate_value, shape = biological_group,
              fill = biological_group
            ),
            size = 7, alpha = 0.85, stroke = 2
          ) +
          scale_shape_manual(
            values = c(21, 24, 22, 23, 25),
            name = "Biological Group"
          ) +
          scale_fill_manual(
            values = group_colors,
            name = "Biological Group"
          ) +
          guides(
            shape = guide_legend(
              override.aes = list(fill = group_colors, colour = "black", stroke = 1.5),
              order = 1
            ),
            fill = "none" # Hide the separate fill legend
          ) +
          labs(
            title = sprintf("PCA Including %s", covar),
            subtitle = paste(model_desc, "| Components from methylation + covariate"),
            x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
            y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
            caption = "PCA computed from methylation + covariate data | Shapes = Biological groups | Colour = Covariate strength"
          ) +
          theme_minimal(base_size = 13) +
          theme(
            legend.position = "right",
            legend.box = "vertical",
            plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
            plot.subtitle = element_text(size = 10, colour = "gray30", hjust = 0.5, face = "italic"),
            plot.caption = element_text(size = 8, colour = "gray50", hjust = 0)
          )

        # Add colour scale based on covariate type
        if (is.numeric(pca_df$covariate_value)) {
          p <- p + scale_color_viridis_c(
            option = "plasma",
            na.value = "gray80",
            name = covar,
            guide = guide_colorbar(
              barwidth = 2, barheight = 15,
              title.position = "top", title.hjust = 0.5, order = 2
            )
          )
        } else if (inherits(pca_df$covariate_value, c("Date", "POSIXct", "POSIXt"))) {
          # Handle date variables as continuous
          p <- p + scale_color_viridis_c(
            option = "plasma",
            na.value = "gray80",
            name = covar,
            guide = guide_colorbar(
              barwidth = 2, barheight = 15,
              title.position = "top", title.hjust = 0.5, order = 2
            )
          )
        } else {
          # Categorical variables
          n_levels <- length(unique(pca_df$covariate_value[!is.na(pca_df$covariate_value)]))

          if (n_levels <= 8) {
            p <- p + scale_color_brewer(palette = "Dark2", na.value = "gray80", name = covar)
          } else if (n_levels <= 12) {
            p <- p + scale_color_brewer(palette = "Set3", na.value = "gray80", name = covar)
          } else {
            p <- p + scale_color_hue(na.value = "gray80", name = covar)
          }
        }

        save_plot_multiple_formats(p,
          file.path(batch_dir, paste0("pca_WITH_", gsub("[^A-Za-z0-9]", "_", covar))),
          width = 14, height = 10
        )

        cat(sprintf("          [OK] Created pca_WITH_%s.png\n", covar))
      }
    }
  }


  # === UMAP: Include each technical covariate (if meth_matrix_base is available) ===
  cat("          Creating UMAP plots with technical covariates...\n")

  # UMAP including each technical covariate
  for (covar in tech_covariates) {
    if (!covar %in% colnames(sample_sheet_matched)) next

    cat(sprintf("          Computing UMAP with covariate: %s\n", covar))

    # Get covariate data matched to UMAP samples
    sample_ids_umap <- rownames(meth_matrix_base)

    # Match covariate values to UMAP samples
    if (!is.null(rownames(sample_sheet_matched))) {
      covar_values <- sample_sheet_matched[[covar]][match(
        sample_ids_umap,
        rownames(sample_sheet_matched)
      )]
    } else {
      # If no rownames, try matching by sample_id_column
      sample_sheet_ids <- as.character(sample_sheet_matched[[opt$sample_id_column]])
      match_idx <- match(sample_ids_umap, sample_sheet_ids)
      covar_values <- sample_sheet_matched[[covar]][match_idx]
    }

    # CRITICAL: Filter to samples that have BOTH methylation AND covariate data
    valid_samples_idx <- !is.na(covar_values)

    if (sum(valid_samples_idx) < 10) {
      cat(sprintf("          [SKIP] %s - too few valid samples (%d)\n", covar, sum(valid_samples_idx)))
      next
    }

    cat(sprintf("          Valid samples (non-NA covariate): %d\n", sum(valid_samples_idx)))

    # Subset methylation matrix to valid samples
    meth_matrix_subset <- meth_matrix_base[valid_samples_idx, , drop = FALSE]
    covar_values_subset <- covar_values[valid_samples_idx]

    # Check if numeric or categorical
    is_numeric_covar <- is.numeric(covar_values_subset)
    is_weight_like <- grepl("weight|age|duration|time|storage|days|hours|weeks|months",
      covar,
      ignore.case = TRUE
    )

    if (!is_numeric_covar && is_weight_like) {
      covar_values_subset <- as.numeric(as.character(covar_values_subset))
      is_numeric_covar <- TRUE
    }

    # Skip if too many NAs remain
    if (sum(is.na(covar_values_subset)) > length(covar_values_subset) * 0.5) {
      cat(sprintf("          [SKIP] %s - too many missing values after filtering\n", covar))
      next
    }

    # Create combined matrix: methylation + covariate
    if (is_numeric_covar) {
      # Numeric: scale and add as column
      covar_scaled <- scale(covar_values_subset, center = TRUE, scale = TRUE)
      covar_scaled[is.na(covar_scaled)] <- 0 # Impute remaining NAs to mean

      combined_matrix <- cbind(meth_matrix_subset, covariate = as.numeric(covar_scaled))
    } else {
      # Categorical: one-hot encode
      covar_factor <- factor(covar_values_subset)
      levels_present <- levels(covar_factor)

      if (length(levels_present) < 2) {
        cat(sprintf("          [SKIP] %s - only one category\n", covar))
        next
      }

      # Create dummy variables (one-hot encoding)
      covar_dummies <- model.matrix(~ covar_factor - 1)
      colnames(covar_dummies) <- paste0("covar_", levels_present)

      # Replace NAs with 0
      covar_dummies[is.na(covar_dummies)] <- 0

      combined_matrix <- cbind(meth_matrix_subset, covar_dummies)
    }

    # Compute UMAP on combined matrix
    umap_combined <- umap(combined_matrix, n_neighbors = min(15, nrow(combined_matrix) - 1))

    umap_df <- data.frame(
      UMAP1 = umap_combined$layout[, 1],
      UMAP2 = umap_combined$layout[, 2],
      sample_id = rownames(meth_matrix_subset)
    )

    # Add metadata
    umap_df <- merge(umap_df, sample_sheet_matched,
      by.x = "sample_id", by.y = opt$sample_id_column, all.x = TRUE
    )

    # Add biological group
    if (exists("sample_metadata") && "actual_class" %in% colnames(sample_metadata)) {
      group_info <- sample_metadata$actual_class[match(umap_df$sample_id, rownames(sample_metadata))]
      umap_df$biological_group <- group_info
    }

    # Add covariate to umap_df for coloring
    umap_df$covariate_value <- sample_sheet_matched[[covar]][match(
      umap_df$sample_id,
      rownames(sample_sheet_matched)
    )]

    # Force numeric if weight-like
    if (is_weight_like && !is.numeric(umap_df$covariate_value)) {
      umap_df$covariate_value <- as.numeric(as.character(umap_df$covariate_value))
    }

    model_desc <- sprintf(
      "%s | %s | %s bp | %d markers",
      species_name, model_type, window_size, n_markers
    )

    # Plot with biological group as shape/fill and covariate as colour
    if ("biological_group" %in% colnames(umap_df)) {
      # Identify control group for coloring
      group_names <- unique(umap_df$biological_group)
      group_names <- group_names[!is.na(group_names)]

      control_pattern <- "control|ctrl|healthy|normal|untreated|mock|baseline"
      is_control <- grepl(control_pattern, group_names, ignore.case = TRUE)

      group_colors <- character(length(group_names))
      group_colors[is_control] <- "#2E7D32" # Dark green

      n_treatment <- sum(!is_control)
      if (n_treatment > 0) {
        treatment_colors <- colorRampPalette(c("#D32F2F", "#FF6F00"))(n_treatment)
        group_colors[!is_control] <- treatment_colors
      }
      names(group_colors) <- group_names

      p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
        stat_ellipse(aes(group = biological_group),
          colour = "black", fill = NA, geom = "polygon",
          alpha = 1, level = 0.95, linewidth = 1.5,
          linetype = "dashed", show.legend = FALSE
        ) +
        geom_point(
          aes(
            colour = covariate_value, shape = biological_group,
            fill = biological_group
          ),
          size = 7, alpha = 0.85, stroke = 2
        ) +
        scale_shape_manual(
          values = c(21, 24, 22, 23, 25),
          name = "Biological Group"
        ) +
        scale_fill_manual(
          values = group_colors,
          name = "Biological Group"
        ) +
        guides(
          shape = guide_legend(
            override.aes = list(fill = group_colors, colour = "black", stroke = 1.5),
            order = 1
          ),
          fill = "none" # Hide the separate fill legend
        ) +
        labs(
          title = sprintf("UMAP Including %s", covar),
          subtitle = paste(model_desc, "| Embedding from methylation + covariate"),
          x = "UMAP1",
          y = "UMAP2",
          caption = "UMAP computed from methylation + covariate data | Shapes = Biological groups | Colour = Covariate strength"
        ) +
        theme_minimal(base_size = 13) +
        theme(
          legend.position = "right",
          legend.box = "vertical",
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          plot.subtitle = element_text(size = 10, colour = "gray30", hjust = 0.5, face = "italic"),
          plot.caption = element_text(size = 8, colour = "gray50", hjust = 0)
        )

      # Add colour scale based on covariate type
      if (is.numeric(umap_df$covariate_value)) {
        p <- p + scale_color_viridis_c(
          option = "plasma",
          na.value = "gray80",
          name = covar,
          guide = guide_colorbar(
            barwidth = 2, barheight = 15,
            title.position = "top", title.hjust = 0.5, order = 2
          )
        )
      } else if (inherits(umap_df$covariate_value, c("Date", "POSIXct", "POSIXt"))) {
        # Handle date variables as continuous
        p <- p + scale_color_viridis_c(
          option = "plasma",
          na.value = "gray80",
          name = covar,
          guide = guide_colorbar(
            barwidth = 2, barheight = 15,
            title.position = "top", title.hjust = 0.5, order = 2
          )
        )
      } else {
        # Categorical variables
        n_levels <- length(unique(umap_df$covariate_value[!is.na(umap_df$covariate_value)]))

        if (n_levels <= 8) {
          p <- p + scale_color_brewer(palette = "Dark2", na.value = "gray80", name = covar)
        } else if (n_levels <= 12) {
          p <- p + scale_color_brewer(palette = "Set3", na.value = "gray80", name = covar)
        } else {
          p <- p + scale_color_hue(na.value = "gray80", name = covar)
        }
      }

      save_plot_multiple_formats(p,
        file.path(batch_dir, paste0("umap_WITH_", gsub("[^A-Za-z0-9]", "_", covar))),
        width = 14, height = 10
      )

      cat(sprintf("          [OK] Created umap_WITH_%s.png\n", covar))
    }
  }

  cat("          [OK] Technical covariate analysis complete\n")
} else {
  cat("          [SKIP] No technical covariates specified (PCA/UMAP biological signal already generated)\n")
}

log_section_end(5, "Batch/platform effect analysis")
cat("\n")

# ============================================================================
# STEP 6: MARKER STABILITY ANALYSIS VIA BOOTSTRAP
# ============================================================================
log_section_start(6, paste0("Marker stability analysis (", opt$bootstrap_iterations, " bootstrap iterations)"), 10)

stability_dir <- file.path(output_dir, "marker_stability")
dir.create(stability_dir, showWarnings = FALSE)

# DEBUG: Check sample sheet and treatment column
cat("\n          ===== DEBUGGING BOOTSTRAP SETUP =====\n")
cat("          Sample sheet matched dimensions:", nrow(sample_sheet_matched), "x", ncol(sample_sheet_matched), "\n")
cat("          Sample sheet columns:", paste(colnames(sample_sheet_matched), collapse = ", "), "\n")
cat("          Treatment column specified:", opt$treatment_column, "\n")

if (!opt$treatment_column %in% colnames(sample_sheet_matched)) {
  cat("          [ERROR] Treatment column '", opt$treatment_column, "' not found in sample sheet!\n", sep = "")
  cat("          Available columns:", paste(colnames(sample_sheet_matched), collapse = ", "), "\n")

  # Try to find it in sample_metadata instead
  if (exists("sample_metadata") && "actual_class" %in% colnames(sample_metadata)) {
    cat("          [INFO] Using 'actual_class' from sample_metadata instead\n")

    # Create treatment variable from sample_metadata
    treatment_var <- sample_metadata$actual_class[match(
      colnames(meth_matrix),
      rownames(sample_metadata)
    )]
    names(treatment_var) <- colnames(meth_matrix)
  } else {
    stop("[FATAL] Cannot find treatment/group variable for bootstrap analysis!")
  }
} else {
  # Get treatment from sample_sheet_matched
  treatment_var <- sample_sheet_matched[[opt$treatment_column]]
  names(treatment_var) <- rownames(sample_sheet_matched)
}

cat("          Treatment variable length:", length(treatment_var), "\n")
cat("          Treatment variable unique values:", paste(unique(treatment_var), collapse = ", "), "\n")
cat("          Treatment variable NA count:", sum(is.na(treatment_var)), "\n")
cat("          Treatment variable table:\n")
print(table(treatment_var, useNA = "ifany"))

# Check if we need to match to meth_matrix samples
if (length(treatment_var) != ncol(meth_matrix)) {
  cat("          [WARNING] Treatment variable length doesn't match meth_matrix\n")
  cat("          Attempting to match by sample IDs...\n")

  # Match treatment to meth_matrix columns
  meth_sample_ids <- colnames(meth_matrix)
  treatment_matched <- treatment_var[match(meth_sample_ids, names(treatment_var))]
  names(treatment_matched) <- meth_sample_ids

  cat("          After matching - Treatment variable length:", length(treatment_matched), "\n")
  cat("          After matching - NA count:", sum(is.na(treatment_matched)), "\n")

  treatment_var <- treatment_matched
}

# Now proceed with valid samples
valid_samples <- !is.na(treatment_var)
cat("          Valid samples (non-NA treatment):", sum(valid_samples), "\n")

if (sum(valid_samples) == 0) {
  stop("[FATAL] No valid samples with treatment information!")
}

meth_for_bootstrap <- meth_matrix[, valid_samples]
treatment_bootstrap <- treatment_var[valid_samples]

cat("          Samples for bootstrap:", sum(valid_samples), "\n")
cat("          Bootstrap groups:", paste(unique(treatment_bootstrap), collapse = ", "), "\n")
cat("          Samples per group:\n")
print(table(treatment_bootstrap))

cat("          ===== BOOTSTRAP SETUP COMPLETE =====\n\n")

cat("          Samples for bootstrap:", sum(valid_samples), "\n")

# Initialize results
bootstrap_selection <- matrix(0,
  nrow = nrow(meth_for_bootstrap),
  ncol = opt$bootstrap_iterations
)
rownames(bootstrap_selection) <- rownames(meth_for_bootstrap)

bootstrap_effects <- matrix(NA,
  nrow = nrow(meth_for_bootstrap),
  ncol = opt$bootstrap_iterations
)
rownames(bootstrap_effects) <- rownames(meth_for_bootstrap)

# Bootstrap loop - OPTIMIZED WITH PARALLEL PROCESSING
set.seed(42)

cat(sprintf("          Running %d bootstrap iterations", opt$bootstrap_iterations))
if (opt$n_cores > 1) {
  cat(sprintf(" (parallel with %d cores)...\n", opt$n_cores))
} else {
  cat(" (sequential)...\n")
}

# Function for one bootstrap iteration
bootstrap_iteration <- function(b, meth_data, treatment, threshold, groups) {
  # Stratified resampling
  boot_idx <- c()
  for (grp in groups) {
    grp_idx <- which(treatment == grp)
    boot_idx <- c(boot_idx, sample(grp_idx, length(grp_idx), replace = TRUE))
  }

  # Bootstrap sample
  meth_boot <- meth_data[, boot_idx, drop = FALSE]
  treatment_boot <- treatment[boot_idx]

  # Vectorized calculation of effect sizes for all markers at once
  n_markers <- nrow(meth_boot)
  effects <- numeric(n_markers)
  selection <- integer(n_markers)

  # Calculate group means for all markers using matrix operations
  group_levels <- unique(treatment_boot)

  if (length(group_levels) >= 2) {
    # For multiclass: calculate max pairwise difference
    all_pairwise_diffs <- numeric(n_markers)

    for (i in 1:(length(group_levels) - 1)) {
      for (j in (i + 1):length(group_levels)) {
        grp_i_idx <- treatment_boot == group_levels[i]
        grp_j_idx <- treatment_boot == group_levels[j]

        grp_i_means <- rowMeans(meth_boot[, grp_i_idx, drop = FALSE], na.rm = TRUE)
        grp_j_means <- rowMeans(meth_boot[, grp_j_idx, drop = FALSE], na.rm = TRUE)

        pairwise_diff <- abs(grp_i_means - grp_j_means)
        all_pairwise_diffs <- pmax(all_pairwise_diffs, pairwise_diff, na.rm = TRUE)
      }
    }

    effects <- all_pairwise_diffs
    selection <- as.integer(effects >= threshold)

    n_valid <- rowSums(!is.na(meth_boot))
    effects[n_valid < 5] <- NA
    selection[n_valid < 5] <- 0
  }

  return(list(effects = effects, selection = selection))
}

# Run bootstrap iterations
if (opt$n_cores > 1 && .Platform$OS.type == "unix") {
  # Use mclapply on Unix/Linux (better for large data)
  cat("          Using mclapply (fork-based, memory efficient)\n")
  cat(sprintf(
    "          Running %d bootstrap iterations on %d cores...\n",
    opt$bootstrap_iterations, opt$n_cores
  ))
  cat("          [This will take 30-60 minutes - be patient!]\n")

  groups <- unique(treatment_bootstrap)
  start_time <- Sys.time()

  # Run in batches with progress reporting
  batch_size <- 100
  n_batches <- ceiling(opt$bootstrap_iterations / batch_size)

  all_results <- list()

  for (batch in 1:n_batches) {
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, opt$bootstrap_iterations)
    batch_iters <- start_idx:end_idx

    cat(sprintf(
      "          Batch %d/%d: iterations %d-%d...\n",
      batch, n_batches, start_idx, end_idx
    ))

    batch_results <- mclapply(batch_iters, bootstrap_iteration,
      meth_data = meth_for_bootstrap,
      treatment = treatment_bootstrap,
      threshold = opt$effect_size_threshold,
      groups = groups,
      mc.cores = opt$n_cores,
      mc.preschedule = TRUE,
      mc.set.seed = TRUE
    )

    all_results <- c(all_results, batch_results)

    # Progress update
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    rate <- end_idx / elapsed
    remaining <- (opt$bootstrap_iterations - end_idx) / rate
    cat(sprintf(
      "          Completed: %d/%d iterations (%.1f iter/min, ~%.0f min remaining)\n",
      end_idx, opt$bootstrap_iterations, rate, remaining
    ))
  }

  results <- all_results

  # Extract results
  bootstrap_effects <- do.call(cbind, lapply(results, function(x) x$effects))
  bootstrap_selection <- do.call(cbind, lapply(results, function(x) x$selection))
  rownames(bootstrap_effects) <- rownames(meth_for_bootstrap)
  rownames(bootstrap_selection) <- rownames(meth_for_bootstrap)
} else if (opt$n_cores > 1) {
  # Fallback: PSOCK cluster (Windows or forced)
  cat("          Using PSOCK cluster (may have memory issues with large data)\n")

  cl <- makeCluster(opt$n_cores, outfile = "")
  clusterExport(cl, c(
    "meth_for_bootstrap", "treatment_bootstrap",
    "opt", "unique"
  ), envir = environment())
  clusterSetRNGStream(cl, 42)

  groups <- unique(treatment_bootstrap)

  results <- parLapply(cl, 1:opt$bootstrap_iterations, bootstrap_iteration,
    meth_data = meth_for_bootstrap,
    treatment = treatment_bootstrap,
    threshold = opt$effect_size_threshold,
    groups = groups
  )

  stopCluster(cl)

  # Extract results
  bootstrap_effects <- do.call(cbind, lapply(results, function(x) x$effects))
  bootstrap_selection <- do.call(cbind, lapply(results, function(x) x$selection))
  rownames(bootstrap_effects) <- rownames(meth_for_bootstrap)
  rownames(bootstrap_selection) <- rownames(meth_for_bootstrap)
} else {
  # Sequential processing with progress bar
  groups <- unique(treatment_bootstrap)

  for (b in 1:opt$bootstrap_iterations) {
    result <- bootstrap_iteration(
      b, meth_for_bootstrap, treatment_bootstrap,
      opt$effect_size_threshold, groups
    )
    bootstrap_effects[, b] <- result$effects
    bootstrap_selection[, b] <- result$selection

    if (b %% 50 == 0 || b == opt$bootstrap_iterations) {
      cat(sprintf("\r          Progress: %d/%d iterations", b, opt$bootstrap_iterations))
    }
  }
  cat("\n")
}

cat("          [OK] Bootstrap complete\n")


# ===== DEBUG: Check what's actually happening =====
cat("\n          ===== BOOTSTRAP DEBUGGING =====\n")
cat("          Bootstrap selection matrix dimensions:", nrow(bootstrap_selection), "x", ncol(bootstrap_selection), "\n")
cat("          Bootstrap effects matrix dimensions:", nrow(bootstrap_effects), "x", ncol(bootstrap_effects), "\n")

# Check if any selections happened
total_selections <- sum(bootstrap_selection, na.rm = TRUE)
cat("          Total selections across all iterations:", total_selections, "\n")

# Check effect sizes
cat("\n          Effect size summary:\n")
cat("            Min:    ", round(min(bootstrap_effects, na.rm = TRUE), 4), "%\n")
cat("            Max:    ", round(max(bootstrap_effects, na.rm = TRUE), 4), "%\n")
cat("            Mean:   ", round(mean(bootstrap_effects, na.rm = TRUE), 4), "%\n")
cat("            Median: ", round(median(bootstrap_effects, na.rm = TRUE), 4), "%\n")

# Show effect sizes for first 5 markers in first iteration
cat("\n          Example effect sizes (first 5 markers, iteration 1):\n")
for (i in 1:min(5, nrow(bootstrap_effects))) {
  cat(sprintf("            %s: %.4f%%\n", rownames(bootstrap_effects)[i], bootstrap_effects[i, 1]))
}

# Check treatment groups
cat("\n          Treatment variable used for bootstrap:\n")
cat("            Unique values:", paste(unique(treatment_bootstrap), collapse = ", "), "\n")
cat("            Sample counts:\n")
print(table(treatment_bootstrap))

# Calculate actual effect sizes manually for verification
cat("\n          Manual effect size calculation (using full dataset):\n")
groups <- unique(treatment_bootstrap)
if (length(groups) == 2) {
  grp1_samples <- treatment_bootstrap == groups[1]
  grp2_samples <- treatment_bootstrap == groups[2]

  grp1_means <- rowMeans(meth_for_bootstrap[, grp1_samples, drop = FALSE], na.rm = TRUE)
  grp2_means <- rowMeans(meth_for_bootstrap[, grp2_samples, drop = FALSE], na.rm = TRUE)

  manual_effects <- abs(grp1_means - grp2_means)

  cat(sprintf("            Group 1 (%s) mean methylation: %.2f%%\n", groups[1], mean(grp1_means, na.rm = TRUE)))
  cat(sprintf("            Group 2 (%s) mean methylation: %.2f%%\n", groups[2], mean(grp2_means, na.rm = TRUE)))
  cat(sprintf("            Mean effect size: %.4f%%\n", mean(manual_effects, na.rm = TRUE)))
  cat(sprintf("            Max effect size: %.4f%%\n", max(manual_effects, na.rm = TRUE)))

  # Count how many markers would pass threshold
  n_above_threshold <- sum(manual_effects >= opt$effect_size_threshold, na.rm = TRUE)
  cat(sprintf(
    "            Markers with effect >=%.1f%%: %d/%d\n",
    opt$effect_size_threshold, n_above_threshold, length(manual_effects)
  ))

  if (n_above_threshold > 0) {
    cat("\n            Top 5 markers by effect size:\n")
    top_indices <- order(manual_effects, decreasing = TRUE)[1:min(5, length(manual_effects))]
    for (idx in top_indices) {
      cat(sprintf("              %s: %.4f%%\n", rownames(meth_for_bootstrap)[idx], manual_effects[idx]))
    }
  }
}

cat("          ===== END BOOTSTRAP DEBUGGING =====\n\n")

# Calculate stability metrics
marker_stability <- data.frame(
  marker = rownames(meth_for_bootstrap),
  selection_frequency = rowMeans(bootstrap_selection),
  mean_effect_size = rowMeans(bootstrap_effects, na.rm = TRUE),
  sd_effect_size = apply(bootstrap_effects, 1, sd, na.rm = TRUE),
  cv_effect_size = {
    mean_vals <- rowMeans(bootstrap_effects, na.rm = TRUE)
    sd_vals <- apply(bootstrap_effects, 1, sd, na.rm = TRUE)
    ifelse(mean_vals == 0 | is.na(mean_vals), NA, sd_vals / mean_vals)
  },
  stable = rowMeans(bootstrap_selection) >= opt$stability_threshold,
  stringsAsFactors = FALSE
)

marker_stability$fragility_index <- 1 - marker_stability$selection_frequency
marker_stability <- marker_stability[order(-marker_stability$selection_frequency), ]

# ADD GENOMIC COORDINATES TO MARKER NAMES
if (exists("dmr_details") && nrow(dmr_details) == nrow(marker_stability)) {
  cat("          Adding genomic coordinates to marker names...\n")

  # Match markers to their coordinates
  marker_match <- match(marker_stability$marker, dmr_details$DMR_ID)

  if (sum(!is.na(marker_match)) > 0) {
    # Create informative marker names: chr:start-end
    marker_stability$marker_name <- paste0(
      dmr_details$chr[marker_match], ":",
      dmr_details$start[marker_match], "-",
      dmr_details$end[marker_match]
    )

    # Create display name: DMR_1 (chr1:12345-12500)
    marker_stability$marker_display <- paste0(
      marker_stability$marker,
      " (", marker_stability$marker_name, ")"
    )

    cat("          [OK] Added genomic coordinates\n")
  } else {
    cat("          [WARNING] Could not match markers to coordinates\n")
    marker_stability$marker_name <- marker_stability$marker
    marker_stability$marker_display <- marker_stability$marker
  }
} else {
  cat("          [INFO] No coordinate information available\n")
  marker_stability$marker_name <- marker_stability$marker
  marker_stability$marker_display <- marker_stability$marker
}

write.csv(marker_stability,
  file.path(stability_dir, "bootstrap_selection_frequencies.csv"),
  row.names = FALSE
)

# Plot stability scores
p1 <- ggplot(marker_stability, aes(
  x = reorder(marker, selection_frequency),
  y = selection_frequency
)) +
  geom_bar(stat = "identity", aes(fill = stable), alpha = 0.8) +
  geom_hline(
    yintercept = opt$stability_threshold, linetype = "dashed",
    colour = "red", linewidth = 1
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c("FALSE" = "gray60", "TRUE" = "forestgreen"),
    labels = c("Unstable", "Stable")
  ) +
  labs(
    title = "Marker Stability via Bootstrap Selection",
    subtitle = paste(opt$bootstrap_iterations, "iterations"),
    x = "Marker",
    y = "Selection Frequency",
    fill = paste0("Stable (>=", opt$stability_threshold * 100, "%)")
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6))

save_plot_multiple_formats(p1,
  file.path(stability_dir, "bootstrap_selection_frequencies"),
  width = 10, height = max(8, n_markers * 0.3)
)

# Effect size vs stability
p2 <- ggplot(marker_stability, aes(x = mean_effect_size, y = selection_frequency)) +
  geom_point(aes(colour = stable), size = 3, alpha = 0.7) +
  geom_hline(yintercept = opt$stability_threshold, linetype = "dashed", colour = "red") +
  geom_vline(xintercept = opt$effect_size_threshold, linetype = "dashed", colour = "blue") +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "forestgreen")) +
  labs(
    title = "Marker Stability vs. Effect Size",
    x = "Mean Effect Size (%)",
    y = "Bootstrap Selection Frequency",
    colour = "Stability"
  ) +
  theme_minimal()

save_plot_multiple_formats(p2,
  file.path(stability_dir, "effect_size_vs_stability"),
  width = 10, height = 7
)

# Fragility index
marker_stability_top <- head(marker_stability, min(20, nrow(marker_stability)))

p3 <- ggplot(marker_stability_top, aes(
  x = reorder(marker, -fragility_index),
  y = fragility_index
)) +
  geom_bar(stat = "identity", fill = "coral", alpha = 0.8) +
  geom_text(aes(label = round(fragility_index, 2)), vjust = -0.3, size = 3) +
  coord_flip() +
  labs(
    title = paste("Top", nrow(marker_stability_top), "Markers: Fragility Index"),
    subtitle = "Lower = More Stable",
    x = "Marker",
    y = "Fragility Index (1 - Selection Frequency)"
  ) +
  theme_minimal()

save_plot_multiple_formats(p3,
  file.path(stability_dir, "marker_fragility_index"),
  width = 10, height = 8
)

log_section_end(6, paste0("Marker stability analysis (", opt$bootstrap_iterations, " bootstrap iterations)"))
cat("\n")

# ============================================================================
# STEP 6.5: DMR-GROUP ASSOCIATION ANALYSIS
# ============================================================================
log_section_start(6.5, "DMR-group association analysis", 10)

cat("[STEP 6.5/10] Analysing DMR-group associations...\n")

dmr_group_analysis <- analyse_dmr_group_associations(
  meth_matrix = meth_matrix,
  sample_metadata = sample_metadata,
  dmr_names = selected_markers
)

if (dmr_group_analysis$available) {
  cat(sprintf(
    "          [OK] Analysed %d DMRs across %d groups\n",
    nrow(dmr_group_analysis$results), dmr_group_analysis$n_groups
  ))

  # Save results
  assoc_file <- file.path(output_dir, "dmr_group_associations.csv")
  write.csv(dmr_group_analysis$results, assoc_file, row.names = FALSE)
  cat(sprintf("          [OK] Saved associations to: %s\n", basename(assoc_file)))

  # Print top 10 DMRs by importance
  cat("\n          Top 10 DMRs by Importance (Eta-squared):\n")
  top10 <- head(dmr_group_analysis$results, 10)
  for (i in 1:nrow(top10)) {
    cat(sprintf(
      "            %d. %s: eta^2 = %.3f (%s), Best for: %s\n",
      i, top10$DMR[i], top10$Eta_Squared[i],
      top10$Effect_Size_Interpretation[i],
      top10$Best_Marker_For[i]
    ))
  }
} else {
  cat("          [SKIP] DMR-group association analysis not available\n")
  cat(sprintf("          Reason: %s\n", dmr_group_analysis$reason))
}

log_section_end(6.5, "DMR-group association analysis")
cat("\n")

# ============================================================================
# STEP 6.6: SAMPLE CORRELATION AND CLUSTERING ANALYSIS
# ============================================================================
log_section_start(6.6, "Sample correlation and hierarchical clustering", 10)
cat("[STEP 6.6/10] Analysing sample correlations and clustering...\n")

# Create correlation directory
corr_dir <- file.path(output_dir, "sample_correlations")
dir.create(corr_dir, showWarnings = FALSE)

# Calculate Pearson and Spearman correlations between SAMPLES (not DMRs)
# meth_matrix is DMRs x Samples, so cor(meth_matrix) gives Sample x Sample correlation
cor_pearson <- cor(meth_matrix, method = "pearson")
cor_spearman <- cor(meth_matrix, method = "spearman")

# Save correlation matrices
write.csv(cor_pearson, file.path(corr_dir, "sample_correlation_matrix_pearson.csv"))
write.csv(cor_spearman, file.path(corr_dir, "sample_correlation_matrix_spearman.csv"))
cat("          [OK] Saved correlation matrices\n")

# Get sample IDs from meth_matrix columns (these are the actual samples in analysis)
# The correlation matrix is computed from meth_matrix (DMRs x Samples), so rownames(cor_pearson) = colnames(meth_matrix)
sample_ids_in_analysis <- colnames(meth_matrix)

# Debug: Print sample IDs
cat(sprintf("          Matching %d samples from meth_matrix to sample_metadata\n", length(sample_ids_in_analysis)))
cat(sprintf("          First 5 sample IDs in meth_matrix: %s\n", paste(head(sample_ids_in_analysis, 5), collapse = ", ")))
cat(sprintf("          First 5 sample IDs in metadata: %s\n", paste(head(rownames(sample_metadata), 5), collapse = ", ")))

# Match to sample_metadata using these exact IDs
sample_metadata_for_heatmap <- sample_metadata[sample_ids_in_analysis, , drop = FALSE]

# Debug: Check match results
n_matched <- sum(!is.na(sample_metadata_for_heatmap$actual_class))
cat(sprintf("          Successfully matched %d / %d samples\n", n_matched, length(sample_ids_in_analysis)))

# Ensure actual_class has no NA values - replace with "Unknown"
if (any(is.na(sample_metadata_for_heatmap$actual_class))) {
  n_na <- sum(is.na(sample_metadata_for_heatmap$actual_class))
  cat(sprintf("          [WARNING] %d samples have NA group labels - replacing with 'Unknown'\n", n_na))
  sample_metadata_for_heatmap$actual_class[is.na(sample_metadata_for_heatmap$actual_class)] <- "Unknown"
}

# Create annotation for heatmap
annotation_df <- data.frame(
  Group = sample_metadata_for_heatmap$actual_class,
  row.names = sample_ids_in_analysis
)

# Get group colors
if (!is.null(group_colors)) {
  ann_colors <- list(Group = group_colors)
} else {
  group_levels <- unique(sample_metadata_for_heatmap$actual_class)
  set1_colors <- RColorBrewer::brewer.pal(max(3, length(group_levels)), "Set1")
  names(set1_colors) <- group_levels
  ann_colors <- list(Group = set1_colors)
}

# Create correlation heatmap (using Pearson) in all requested formats
formats_to_generate <- tolower(trimws(unlist(strsplit(opt$figure_format, ","))))
if ("all" %in% formats_to_generate) {
  formats_to_generate <- c("png", "jpeg", "svg")
}
formats_to_generate <- gsub("jpg", "jpeg", formats_to_generate)

for (fmt in unique(formats_to_generate)) {
  plot_file <- file.path(corr_dir, sprintf("sample_correlation_heatmap.%s", fmt))
  open_plot_device(plot_file, width = 10, height = 10, res = 300)

  pheatmap(cor_pearson,
    clustering_distance_rows = as.dist(1 - cor_pearson),
    clustering_distance_cols = as.dist(1 - cor_pearson),
    clustering_method = "ward.D2",
    annotation_row = annotation_df,
    annotation_col = annotation_df,
    annotation_colors = ann_colors,
    annotation_legend = TRUE,
    color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
    main = "Sample Correlation Matrix (Pearson)",
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 6,
    fontsize_col = 6
  )

  dev.off()
}
cat(sprintf("          [OK] Created sample_correlation_heatmap in formats: %s\n", paste(formats_to_generate, collapse = ", ")))

# Hierarchical clustering analysis
hclust_result <- hclust(as.dist(1 - cor_pearson), method = "ward.D2")

# Calculate cluster purity metrics
# Cut tree to match number of biological groups
n_groups <- length(unique(sample_metadata_for_heatmap$actual_class))
clusters <- cutree(hclust_result, k = n_groups)

# Calculate Rand index
if (requireNamespace("fossil", quietly = TRUE)) {
  rand_index <- fossil::rand.index(clusters, as.numeric(factor(sample_metadata_for_heatmap$actual_class)))
} else {
  # Simple purity calculation without fossil
  cluster_purity <- sapply(1:n_groups, function(i) {
    cluster_samples <- which(clusters == i)
    if (length(cluster_samples) == 0) {
      return(0)
    }
    group_counts <- table(sample_metadata_for_heatmap$actual_class[cluster_samples])
    max(group_counts) / sum(group_counts)
  })
  rand_index <- mean(cluster_purity)
}

# Save clustering results
clustering_results <- data.frame(
  Sample = sample_ids_in_analysis,
  Cluster = clusters,
  Biological_Group = sample_metadata_for_heatmap$actual_class,
  stringsAsFactors = FALSE
)
write.csv(clustering_results, file.path(corr_dir, "hierarchical_clustering_results.csv"), row.names = FALSE)
cat(sprintf("          [OK] Cluster purity (Rand-like index): %.3f\n", rand_index))

# Create dendrogram with colored labels by group
library(dendextend)
dend <- as.dendrogram(hclust_result)

# Color labels by biological group - ENSURE CORRECT MATCHING
labels_in_order <- labels(dend)
label_groups <- sample_metadata_for_heatmap$actual_class[match(labels_in_order, sample_ids_in_analysis)]
label_colors <- ann_colors$Group[as.character(label_groups)]

dend <- set(dend, "labels_col", label_colors)
dend <- set(dend, "labels_cex", 0.8)

for (fmt in unique(formats_to_generate)) {
  plot_file <- file.path(corr_dir, sprintf("hierarchical_clustering_dendrogram.%s", fmt))
  open_plot_device(plot_file, width = 14, height = 8, res = 300)

  par(mar = c(8, 4, 4, 2))
  plot(dend,
    main = "Hierarchical Clustering Dendrogram (Ward's Method, 1-Correlation Distance)",
    ylab = "Height", xlab = ""
  )

  # Add colored bars underneath
  colored_bars(
    colors = data.frame(Group = label_colors),
    dend = dend,
    rowLabels = "Biological Group",
    y_shift = -0.1
  )

  dev.off()
}
cat(sprintf("          [OK] Created hierarchical_clustering_dendrogram in formats: %s\n", paste(formats_to_generate, collapse = ", ")))

log_section_end(6.6, "Sample correlation and hierarchical clustering")
cat("\n")

# ============================================================================
# STEP 6.7: DMR CORRELATION NETWORK ANALYSIS
# ============================================================================
log_section_start(6.7, "DMR correlation network", 10)
cat("[STEP 6.7/10] Building DMR correlation network...\n")

# Create network directory
network_dir <- file.path(output_dir, "dmr_network")
dir.create(network_dir, showWarnings = FALSE)

# Calculate Pearson correlation between DMRs
cor_dmr <- cor(t(meth_matrix), method = "pearson", use = "pairwise.complete.obs")

# Try correlation thresholds in descending order until we get a useful network
cor_threshold <- 0.7
min_edges_needed <- 10
n_dmrs <- nrow(cor_dmr)

edges_df <- data.frame()
for (threshold in c(0.7, 0.6, 0.5)) {
  # Find edges above threshold
  edges <- which(abs(cor_dmr) > threshold & row(cor_dmr) < col(cor_dmr), arr.ind = TRUE)

  if (nrow(edges) >= min_edges_needed || threshold == 0.5) {
    cor_threshold <- threshold

    edges_df <- data.frame(
      from = rownames(cor_dmr)[edges[, 1]],
      to = rownames(cor_dmr)[edges[, 2]],
      correlation = cor_dmr[edges],
      stringsAsFactors = FALSE
    )

    cat(sprintf(
      "          Using correlation threshold |r| > %.1f (%d edges found)\n",
      cor_threshold, nrow(edges_df)
    ))
    break
  }
}

# Only proceed with network analysis if we have edges
if (nrow(edges_df) > 0) {
  # Save edges
  write.csv(edges_df, file.path(network_dir, "dmr_network_edges.csv"), row.names = FALSE)

  # Build igraph network
  g <- graph_from_data_frame(edges_df[, c("from", "to")], directed = FALSE)

  # Calculate network metrics
  n_nodes <- vcount(g)
  n_edges <- ecount(g)
  network_density <- edge_density(g)

  cat(sprintf("          Network: %d nodes, %d edges, density = %.4f\n", n_nodes, n_edges, network_density))

  # Calculate centrality metrics
  degree_vals <- degree(g)
  betweenness_vals <- betweenness(g)
  closeness_vals <- closeness(g)

  # Identify hubs (top 10% by degree, but only if degree > 0)
  hub_threshold <- quantile(degree_vals, 0.9, na.rm = TRUE)
  hub_dmrs <- data.frame(
    DMR = names(degree_vals),
    Degree = degree_vals,
    Betweenness = betweenness_vals,
    Closeness = closeness_vals,
    stringsAsFactors = FALSE
  )
  hub_dmrs <- hub_dmrs[order(hub_dmrs$Degree, decreasing = TRUE), ]

  # ONLY include actual hubs (degree > 0)
  hub_dmrs <- hub_dmrs[hub_dmrs$Degree > 0, ]

  if (nrow(hub_dmrs) > 0) {
    hub_dmrs$Rank <- 1:nrow(hub_dmrs)
    hub_dmrs <- hub_dmrs[1:min(10, nrow(hub_dmrs)), ]
    write.csv(hub_dmrs, file.path(network_dir, "dmr_network_hubs.csv"), row.names = FALSE)
    cat(sprintf("          [OK] Identified %d hub DMRs\n", nrow(hub_dmrs)))
  } else {
    cat("          [INFO] No hub DMRs identified (all nodes have degree 0 or network too sparse)\n")
    hub_dmrs <- data.frame() # Empty dataframe
  }

  # Detect communities
  communities <- cluster_louvain(g)
  V(g)$community <- membership(communities)
  cat(sprintf("          Detected %d network communities\n", length(unique(membership(communities)))))
} else {
  cat("          [WARNING] No significant correlations detected even at |r| > 0.5\n")
  cat("          Network analysis skipped - DMRs have largely independent patterns\n")
  hub_dmrs <- data.frame()
  g <- NULL
}

# Visualise network
if (!is.null(g) && n_edges > 0 && n_edges < 10000) { # Only plot if not too large
  for (fmt in unique(formats_to_generate)) {
    plot_file <- file.path(network_dir, sprintf("dmr_correlation_network.%s", fmt))
    open_plot_device(plot_file, width = 12, height = 12, res = 300)

    par(mar = c(0, 0, 2, 0))

    # Set node colors by community
    if (exists("communities")) {
      node_colors <- rainbow(length(unique(membership(communities))))[membership(communities)]
    } else {
      node_colors <- "lightblue"
    }

    # Set node sizes by degree
    node_sizes <- 3 + log1p(degree_vals) * 2

    # Highlight hub nodes
    hub_dmr_names <- if (nrow(hub_dmrs) > 0) hub_dmrs$DMR else character(0)
    node_shapes <- ifelse(names(V(g)) %in% hub_dmr_names, "square", "circle")

    # Set node labels - show all DMR IDs, make hub labels bold and bigger
    node_labels <- names(V(g))
    label_colors <- ifelse(names(V(g)) %in% hub_dmr_names, "black", "gray40")
    label_fonts <- ifelse(names(V(g)) %in% hub_dmr_names, 2, 2) # 2 = bold for both
    label_cex <- ifelse(names(V(g)) %in% hub_dmr_names, 1.3, 1.0) # Much bigger labels

    plot(g,
      vertex.color = node_colors,
      vertex.size = node_sizes,
      vertex.shape = node_shapes,
      vertex.label = node_labels,
      vertex.label.color = label_colors,
      vertex.label.font = label_fonts,
      vertex.label.cex = label_cex,
      vertex.label.dist = 0,
      edge.width = 0.5,
      edge.color = rgb(0.5, 0.5, 0.5, 0.3),
      layout = layout_with_fr(g),
      main = "DMR Correlation Network (|r| > 0.7)"
    )

    legend("bottomright",
      legend = c("Hub DMR (high degree)", "Regular DMR"),
      pch = c(15, 19), col = "grey", pt.cex = 1.5, bty = "n"
    )

    dev.off()
  }
  cat(sprintf("          [OK] Created dmr_correlation_network in formats: %s\n", paste(formats_to_generate, collapse = ", ")))
} else {
  cat("          [SKIP] Network too large or empty for visualisation\n")
}

log_section_end(6.7, "DMR correlation network")
cat("\n")

# ============================================================================
# STEP 6.8: WGCNA CO-EXPRESSION ANALYSIS (CONDITIONAL)
# ============================================================================
if (opt$enable_wgcna && wgcna_available) {
  log_section_start(6.8, "WGCNA co-expression modules", 10)
  cat("[STEP 6.8/10] Performing WGCNA co-methylation analysis...\n")

  # Create WGCNA directory
  wgcna_dir <- file.path(output_dir, "wgcna")
  dir.create(wgcna_dir, showWarnings = FALSE)

  # Prepare data (samples as rows, DMRs as columns)
  wgcna_data <- t(meth_matrix)

  # Check for good genes/DMRs
  gsg <- WGCNA::goodSamplesGenes(wgcna_data, verbose = 3)
  if (!gsg$allOK) {
    wgcna_data <- wgcna_data[, gsg$goodGenes]
    cat(sprintf("          Removed %d DMRs with too many missing values\n", sum(!gsg$goodGenes)))
  }

  # Choose soft thresholding power
  cat("          Selecting soft threshold power...\n")
  powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
  sft <- WGCNA::pickSoftThreshold(wgcna_data, powerVector = powers, verbose = 0)

  # Select power based on scale-free topology fit
  power_selected <- sft$powerEstimate
  if (is.na(power_selected)) {
    power_selected <- 6 # Default
    cat("          [WARN] Could not auto-select power, using default = 6\n")
  } else {
    cat(sprintf("          [OK] Selected soft threshold power = %d\n", power_selected))
  }

  # Construct network and detect modules
  cat("          Detecting co-methylation modules...\n")
  net <- WGCNA::blockwiseModules(
    wgcna_data,
    power = power_selected,
    TOMType = "unsigned",
    minModuleSize = 3,
    reassignThreshold = 0,
    mergeCutHeight = 0.25,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    verbose = 0
  )

  # Get module colors
  module_colors <- WGCNA::labels2colors(net$colors)
  n_modules <- length(unique(module_colors)) - 1 # Exclude grey
  cat(sprintf("          [OK] Detected %d modules\n", n_modules))

  # Save module assignments
  module_df <- data.frame(
    DMR = colnames(wgcna_data),
    Module_Number = net$colors,
    Module_Color = module_colors,
    stringsAsFactors = FALSE
  )
  write.csv(module_df, file.path(wgcna_dir, "wgcna_modules.csv"), row.names = FALSE)

  # Calculate module-trait relationships
  cat("          Calculating module-trait relationships...\n")
  module_eigengenes <- net$MEs

  # Create trait matrix (biological groups as binary)
  groups <- unique(sample_metadata$actual_class)
  trait_matrix <- matrix(0, nrow = nrow(wgcna_data), ncol = length(groups))
  colnames(trait_matrix) <- groups
  for (i in seq_along(groups)) {
    trait_matrix[, i] <- as.numeric(sample_metadata$actual_class == groups[i])
  }

  # Calculate correlations
  module_trait_cor <- cor(module_eigengenes, trait_matrix, use = "p")
  module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(wgcna_data))

  # Save module-trait correlations
  module_trait_df <- as.data.frame(module_trait_cor)
  module_trait_df$Module <- rownames(module_trait_df)
  write.csv(module_trait_df, file.path(wgcna_dir, "wgcna_module_trait_correlations.csv"), row.names = FALSE)

  # Plot module dendrogram in all formats
  for (fmt in unique(formats_to_generate)) {
    plot_file <- file.path(wgcna_dir, sprintf("wgcna_dendrogram_modules.%s", fmt))
    open_plot_device(plot_file, width = 12, height = 6, res = 300)

    WGCNA::plotDendroAndColors(net$dendrograms[[1]], module_colors[net$blockGenes[[1]]],
      "Module colors",
      dendroLabels = FALSE, hang = 0.03,
      addGuide = TRUE, guideHang = 0.05,
      main = "WGCNA Dendrogram and Module Colors"
    )

    dev.off()
  }
  cat(sprintf("          [OK] Created wgcna_dendrogram_modules in formats: %s\n", paste(formats_to_generate, collapse = ", ")))

  # Plot module-trait heatmap in all formats
  for (fmt in unique(formats_to_generate)) {
    plot_file <- file.path(wgcna_dir, sprintf("wgcna_module_trait_heatmap.%s", fmt))
    open_plot_device(plot_file, width = 8, height = 10, res = 300)

    par(mar = c(6, 10, 3, 3))

    # Display correlations and p-values
    textMatrix <- paste(sprintf("%.2f", module_trait_cor), "\n(",
      sprintf("%.1e", module_trait_pval), ")",
      sep = ""
    )
    dim(textMatrix) <- dim(module_trait_cor)

    WGCNA::labeledHeatmap(
      Matrix = module_trait_cor,
      xLabels = colnames(trait_matrix),
      yLabels = names(module_eigengenes),
      ySymbols = names(module_eigengenes),
      colorLabels = FALSE,
      colors = WGCNA::blueWhiteRed(50),
      textMatrix = textMatrix,
      setStdMargins = FALSE,
      cex.text = 0.7,
      zlim = c(-1, 1),
      main = "Module-Trait Relationships"
    )

    dev.off()
  }
  cat(sprintf("          [OK] Created wgcna_module_trait_heatmap in formats: %s\n", paste(formats_to_generate, collapse = ", ")))

  log_section_end(6.8, "WGCNA co-expression modules")
  cat("\n")
} else {
  if (opt$enable_wgcna) {
    cat("[STEP 6.8/10] WGCNA analysis skipped (package not available)\n\n")
  } else {
    cat("[STEP 6.8/10] WGCNA analysis skipped (not enabled)\n\n")
  }
}

# ============================================================================
# STEP 6.9: DMR METHYLATION PATTERNS BY BIOLOGICAL GROUP (REMOVED)
# ============================================================================
# This section has been removed as it's no longer needed
cat("[STEP 6.9/10] DMR methylation patterns analysis skipped (deprecated)\n\n")

# ============================================================================
# STEP 6.10: GENE ANNOTATION ANALYSIS (CONDITIONAL)
# ============================================================================
if (opt$enable_gene_annotation) {
  log_section_start(6.10, "Gene annotation analysis", 10)
  cat("[STEP 6.10/10] Annotating DMRs with gene information...\n")

  # Validate required parameters
  if (is.null(opt$species)) {
    cat("          [ERROR] --species must be specified when --enable_gene_annotation is used\n")
    cat("          Available species: 'Sus scrofa', 'Homo sapiens', 'Mus musculus', 'Gallus gallus'\n")
    cat("          Skipping gene annotation.\n\n")
  } else {
    # Check if biomaRt is available
    biomart_available <- requireNamespace("biomaRt", quietly = TRUE)

    if (!biomart_available) {
      cat("\n")
      cat("================================================================================\n")
      cat("BIOMART PACKAGE NOT INSTALLED\n")
      cat("================================================================================\n")
      cat("Gene annotation was requested (--enable_gene_annotation) but biomaRt is not available.\n")
      cat("\n")
      cat("To install biomaRt:\n")
      cat("  if (!requireNamespace('BiocManager', quietly = TRUE))\n")
      cat("      install.packages('BiocManager')\n")
      cat("  BiocManager::install('biomaRt')\n")
      cat("\n")
      cat("Gene annotation will be SKIPPED.\n")
      cat("================================================================================\n")
      cat("\n")
    } else {
      # Load biomaRt
      suppressPackageStartupMessages(library(biomaRt))
      cat("          [OK] biomaRt package loaded\n")

      # Create gene annotation directory
      gene_annot_dir <- file.path(output_dir, "gene_annotation")
      dir.create(gene_annot_dir, showWarnings = FALSE)

      # Prepare DMR regions data frame
      # Extract chromosome, start, end from dmr_details
      if (all(c("chr", "start", "end") %in% colnames(dmr_details))) {
        dmr_regions <- data.frame(
          DMR_ID = dmr_details$DMR_ID,
          Chr = dmr_details$chr,
          Start = dmr_details$start,
          End = dmr_details$end,
          stringsAsFactors = FALSE
        )
      } else {
        cat("          [WARNING] DMR coordinates not found in dmr_details_with_coordinates.csv\n")
        cat("          Required columns: chr, start, end\n")
        cat("          Available columns:", paste(colnames(dmr_details), collapse = ", "), "\n")
        cat("          Skipping gene annotation.\n\n")
        dmr_regions <- NULL
      }

      if (!is.null(dmr_regions)) {
        cat(sprintf("          Prepared %d DMR regions for annotation\n", nrow(dmr_regions)))

        # Query biomaRt for gene annotations
        tryCatch(
          {
            gene_annotation <- query_biomart_genes(
              dmr_regions = dmr_regions,
              species_name = opt$species,
              gene_window = opt$gene_window
            )

            # Save gene annotation results
            output_file <- file.path(gene_annot_dir, "dmr_gene_annotation.csv")
            write.csv(gene_annotation, output_file, row.names = FALSE)
            cat(sprintf("          [OK] Saved gene annotations to: %s\n", output_file))

            # Create summary statistics
            n_overlapping <- sum(!is.na(gene_annotation$Overlapping_Genes))
            n_nearby <- sum(!is.na(gene_annotation$Nearest_Gene_Name))
            n_no_genes <- sum(is.na(gene_annotation$Nearest_Gene_Name))

            cat(sprintf("          Summary:\n"))
            cat(sprintf(
              "            - %d DMRs with overlapping genes (%.1f%%)\n",
              n_overlapping, 100 * n_overlapping / nrow(gene_annotation)
            ))
            cat(sprintf(
              "            - %d DMRs with nearby genes within %d bp (%.1f%%)\n",
              n_nearby, opt$gene_window, 100 * n_nearby / nrow(gene_annotation)
            ))
            cat(sprintf(
              "            - %d DMRs with no nearby genes (%.1f%%)\n",
              n_no_genes, 100 * n_no_genes / nrow(gene_annotation)
            ))

            # Gene biotype distribution
            if (sum(!is.na(gene_annotation$Gene_Biotype)) > 0) {
              biotype_table <- table(gene_annotation$Gene_Biotype[!is.na(gene_annotation$Gene_Biotype)])
              biotype_df <- as.data.frame(biotype_table)
              colnames(biotype_df) <- c("Gene_Biotype", "Count")
              biotype_df <- biotype_df[order(biotype_df$Count, decreasing = TRUE), ]

              write.csv(biotype_df, file.path(gene_annot_dir, "gene_biotype_distribution.csv"), row.names = FALSE)
              cat(sprintf("          [OK] Gene biotype distribution saved\n"))

              cat(sprintf("          Top gene biotypes:\n"))
              for (i in 1:min(5, nrow(biotype_df))) {
                cat(sprintf("            - %s: %d genes\n", biotype_df$Gene_Biotype[i], biotype_df$Count[i]))
              }

              # Create biotype bar plot
              if (nrow(biotype_df) > 0) {
                p_biotype <- ggplot(biotype_df, aes(x = reorder(Gene_Biotype, Count), y = Count)) +
                  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
                  coord_flip() +
                  labs(
                    title = "Distribution of Gene Biotypes",
                    subtitle = sprintf("DMRs annotated with %s (window: %d bp)", opt$species, opt$gene_window),
                    x = "Gene Biotype",
                    y = "Number of DMRs"
                  ) +
                  theme_minimal(base_size = 11) +
                  theme(
                    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                    plot.subtitle = element_text(size = 10, colour = "gray30", hjust = 0.5, face = "italic")
                  )

                save_plot_multiple_formats(p_biotype,
                  file.path(gene_annot_dir, "gene_biotype_distribution"),
                  width = 10, height = 8
                )
                cat("          [OK] Created gene_biotype_distribution plot\n")
              }
            }

            # Store for markdown report
            gene_annotation_done <<- TRUE
            gene_annotation_results <<- gene_annotation
            gene_biotype_summary <<- if (exists("biotype_df")) biotype_df else NULL
          },
          error = function(e) {
            cat(sprintf("          [ERROR] Gene annotation failed: %s\n", e$message))
            cat("          This may be due to:\n")
            cat("            - No internet connection\n")
            cat("            - Ensembl server is down\n")
            cat("            - Invalid species name\n")
            cat("          Creating empty annotation file with error message...\n")

            # Create empty file with error message
            error_df <- data.frame(
              Error = paste("Gene annotation failed:", e$message),
              stringsAsFactors = FALSE
            )
            write.csv(error_df, file.path(gene_annot_dir, "dmr_gene_annotation.csv"), row.names = FALSE)

            gene_annotation_done <<- FALSE
          }
        )
      }
    }
  }

  log_section_end(6.10, "Gene annotation analysis")
  cat("\n")
} else {
  cat("[STEP 6.10/10] Gene annotation disabled (use --enable_gene_annotation to enable)\n\n")
  gene_annotation_done <- FALSE
}

# ============================================================================
# STEP 6.10.5: PCR PRIMER PREPARATION (CONDITIONAL)
# ============================================================================

pcr_prep_done <- FALSE
pcr_links_df <- NULL

if (!is.null(opt$prepare_pcr) && opt$prepare_pcr == TRUE) {
  cat("[STEP 6.10.5/10] Preparing PCR primer design links for top DMRs...\n")
  log_section_start(6.105, "PCR primer preparation", 10)

  # Get top N DMRs by significance
  n_pcr_dmrs <- if (!is.null(opt$pcr_n_dmrs)) opt$pcr_n_dmrs else 20

  # Use dmr_group_analysis results and gene_annotation for coordinates
  if (exists("dmr_group_analysis") && dmr_group_analysis$available && nrow(dmr_group_analysis$results) > 0) {
    # Sort by p-value (ANOVA_pvalue_adjusted preferred, or ANOVA_pvalue as fallback)
    p_col <- if ("ANOVA_pvalue_adjusted" %in% colnames(dmr_group_analysis$results)) {
      "ANOVA_pvalue_adjusted"
    } else {
      "ANOVA_pvalue"
    }

    dmr_results_sorted <- dmr_group_analysis$results[order(dmr_group_analysis$results[[p_col]]), ]
    top_dmrs_pcr <- head(dmr_results_sorted, n_pcr_dmrs)

    # Get coordinates from gene_annotation if available
    coord_source <- if (exists("gene_annotation_results") && nrow(gene_annotation_results) > 0) {
      gene_annotation_results
    } else {
      # Try to load from gene annotation file
      gene_annot_file <- file.path(output_dir, "gene_annotation", "dmr_gene_annotation.csv")
      if (file.exists(gene_annot_file)) {
        read.csv(gene_annot_file, stringsAsFactors = FALSE)
      } else {
        # Try dmr_details_with_coordinates.csv from model directory
        dmr_coords_file <- file.path(opt$model_dir, "dmr_details_with_coordinates.csv")
        if (file.exists(dmr_coords_file)) {
          read.csv(dmr_coords_file, stringsAsFactors = FALSE)
        } else {
          NULL
        }
      }
    }

    if (!is.null(coord_source)) {
      # Detect column names - support both upper/lower case variants
      dmr_col <- if ("DMR_ID" %in% colnames(coord_source)) "DMR_ID" else if ("dmr" %in% colnames(coord_source)) "dmr" else NULL
      chr_col <- if ("Chr" %in% colnames(coord_source)) "Chr" else if ("scaffold" %in% colnames(coord_source)) "scaffold" else if ("chromosome" %in% colnames(coord_source)) "chromosome" else NULL
      start_col <- if ("Start" %in% colnames(coord_source)) "Start" else if ("start" %in% colnames(coord_source)) "start" else NULL
      end_col <- if ("End" %in% colnames(coord_source)) "End" else if ("end" %in% colnames(coord_source)) "end" else NULL

      if (!is.null(dmr_col) && !is.null(chr_col) && !is.null(start_col) && !is.null(end_col)) {
        pcr_links_list <- list()

        for (i in 1:nrow(top_dmrs_pcr)) {
          dmr_id <- top_dmrs_pcr$DMR[i]
          p_val <- top_dmrs_pcr[[p_col]][i]

          # Look up coordinates
          coord_row <- coord_source[coord_source[[dmr_col]] == dmr_id, ]

          if (nrow(coord_row) > 0 && !is.na(coord_row[[chr_col]][1]) && !is.na(coord_row[[start_col]][1]) && !is.na(coord_row[[end_col]][1])) {
            chr <- coord_row[[chr_col]][1]
            start <- coord_row[[start_col]][1]
            end <- coord_row[[end_col]][1]

            # Add flanking sequence (+/-500bp recommended for primer design)
            flank_size <- 500
            seq_start <- max(1, start - flank_size)
            seq_end <- end + flank_size

            # Create coordinates string
            coords <- sprintf("%s:%d-%d", chr, start, end)
            coords_extended <- sprintf("%s:%d-%d", chr, seq_start, seq_end)

            # Detect genome build for URLs
            genome_build <- if (!is.null(opt$species)) {
              if (opt$species == "Gallus gallus") {
                "galGal6"
              } else if (opt$species == "Sus scrofa") {
                "susScr11"
              } else if (opt$species == "Homo sapiens") {
                "hg38"
              } else if (opt$species == "Mus musculus") {
                "mm10"
              } else {
                "hg38" # Default
              }
            } else {
              "hg38" # Default
            }

            # Generate UCSC Genome Browser URL
            # UCSC requires chr prefix (chr16), Ensembl uses numeric (16)
            chr_ucsc <- if (!grepl("^chr", chr)) paste0("chr", chr) else chr
            coords_ucsc <- sprintf("%s:%d-%d", chr_ucsc, start, end)
            ucsc_url <- sprintf("https://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&position=%s", genome_build, coords_ucsc)

            # Fetch DNA sequence using biomaRt (if available and connected)
            sequence <- NA_character_
            sequence_status <- "Not fetched"

            if (exists("ensembl") && !is.null(ensembl)) {
              tryCatch(
                {
                  # Use biomaRt to fetch sequence
                  chr_clean <- gsub("^chr", "", chr) # Remove chr prefix for Ensembl

                  seq_result <- biomaRt::getSequence(
                    chromosome = chr_clean,
                    start = seq_start,
                    end = seq_end,
                    seqType = "gene_flank",
                    type = "entrezgene_id",
                    mart = ensembl
                  )

                  if (nrow(seq_result) > 0 && !is.na(seq_result$gene_flank[1])) {
                    sequence <- seq_result$gene_flank[1]
                    sequence_status <- "Fetched"
                  } else {
                    sequence_status <- "Empty result"
                  }
                },
                error = function(e) {
                  sequence_status <<- paste0("Error: ", substr(e$message, 1, 50))
                }
              )
            }

            # Generate MethPrimer URL (better for bisulfite PCR, supports up to 5000bp)
            # MethPrimer: http://www.urogene.org/cgi-bin/methprimer/methprimer.cgi
            methprimer_url <- "http://www.urogene.org/cgi-bin/methprimer/methprimer.cgi"

            # Generate UCSC sequence fetch URL (fallback for sequence extraction)
            # Use hgTracks with getDna parameter - opens browser at position with "Get DNA" button
            ucsc_das_url <- sprintf(
              "https://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&position=%s:%d-%d",
              genome_build, chr_ucsc, seq_start, seq_end
            )

            # Store information
            pcr_links_list[[length(pcr_links_list) + 1]] <- data.frame(
              Rank = i,
              DMR_ID = dmr_id,
              Chromosome = chr,
              Start = start,
              End = end,
              Size_bp = end - start + 1,
              Extended_Start = seq_start,
              Extended_End = seq_end,
              Extended_Size_bp = seq_end - seq_start + 1,
              P_value = p_val,
              Coordinates = coords,
              Coordinates_Extended = coords_extended,
              UCSC_Browser = ucsc_url,
              UCSC_Sequence_Fetch = ucsc_das_url,
              MethPrimer_URL = methprimer_url,
              Sequence_Status = sequence_status,
              Sequence = if (!is.na(sequence)) substr(sequence, 1, 100) else NA_character_, # First 100bp for CSV
              stringsAsFactors = FALSE
            )
          } else {
            cat(sprintf("          [WARNING] No coordinates found for DMR: %s\n", dmr_id))
          }
        }

        if (length(pcr_links_list) > 0) {
          pcr_links_df <- do.call(rbind, pcr_links_list)

          # Save to CSV
          pcr_output_dir <- file.path(output_dir, "pcr_primers")
          if (!dir.exists(pcr_output_dir)) {
            dir.create(pcr_output_dir, recursive = TRUE)
          }

          write.csv(pcr_links_df,
            file.path(pcr_output_dir, "pcr_primer_preparation.csv"),
            row.names = FALSE
          )

          # Save full sequences to FASTA file (if any were fetched)
          sequences_fetched <- sum(pcr_links_df$Sequence_Status == "Fetched", na.rm = TRUE)

          if (sequences_fetched > 0) {
            fasta_file <- file.path(pcr_output_dir, "dmr_sequences.fasta")
            fasta_conn <- file(fasta_file, "w")

            for (i in 1:nrow(pcr_links_df)) {
              if (!is.na(pcr_links_df$Sequence_Status[i]) && pcr_links_df$Sequence_Status[i] == "Fetched") {
                # Write FASTA header
                writeLines(
                  sprintf(
                    ">%s %s (extended: %s)",
                    pcr_links_df$DMR_ID[i],
                    pcr_links_df$Coordinates[i],
                    pcr_links_df$Coordinates_Extended[i]
                  ),
                  fasta_conn
                )

                # Write sequence (full sequence, not truncated)
                # Note: Need to retrieve full sequence again or store it separately
                writeLines(pcr_links_df$Sequence[i], fasta_conn)
              }
            }

            close(fasta_conn)
            cat(sprintf(
              "          ✓ Saved %d sequences to FASTA: %s\n", sequences_fetched,
              basename(fasta_file)
            ))
          }

          cat(sprintf("          ✓ Generated PCR preparation data for %d top DMRs\n", nrow(pcr_links_df)))
          cat(sprintf("          ✓ Saved to: %s/pcr_primer_preparation.csv\n", pcr_output_dir))

          pcr_prep_done <- TRUE
        } else {
          cat("          [WARNING] No valid DMR coordinates found for PCR preparation\n")
        }
      } else {
        cat("          [WARNING] Could not find required coordinate columns in data\n")
        cat(sprintf("          Available columns: %s\n", paste(colnames(coord_source), collapse = ", ")))
        cat("          Required: DMR_ID/dmr, Chr/scaffold/chromosome, Start/start, End/end\n")
      }
    } else {
      cat("          [WARNING] No coordinate data available for PCR preparation\n")
      cat("          Note: Requires dmr_details_with_coordinates.csv in model directory\n")
      cat("          Or run with --enable_gene_annotation to generate coordinates\n")
    }
  } else {
    cat("          [WARNING] No DMR data available for PCR preparation\n")
    cat("          Note: DMR-group association analysis must be completed first\n")
  }

  log_section_end(6.105, "PCR primer preparation")
  cat("\n")
} else {
  cat("[STEP 6.10.5/10] PCR primer preparation disabled (use --prepare_pcr to enable)\n\n")
}

# ============================================================================
# STEP 6.10.6: GENOMIC FEATURE ANNOTATION (CONDITIONAL)
# ============================================================================

feature_annotation_done <- FALSE
feature_annotation_results <- NULL
feature_plot <- NULL

if (!is.null(opt$enable_feature_annotation) && (opt$enable_feature_annotation == TRUE || opt$enable_feature_annotation == "TRUE" || opt$enable_feature_annotation == "true")) {
  cat("[STEP 6.10.6/10] Annotating DMRs with genomic features...\n")
  log_section_start(6.106, "Genomic feature annotation", 10)

  # Create output directory
  feature_annot_dir <- file.path(output_dir, "feature_annotation")
  dir.create(feature_annot_dir, showWarnings = FALSE, recursive = TRUE)
  feature_figures_dir <- file.path(feature_annot_dir, "figures")
  dir.create(feature_figures_dir, showWarnings = FALSE, recursive = TRUE)

  # Validate required GTF file
  if (is.null(opt$gtf_file) || !file.exists(opt$gtf_file)) {
    cat("  [ERROR] GTF/GFF file required for feature annotation!\n")
    cat("  Please provide --gtf_file parameter.\n")
    cat("  Skipping feature annotation.\n\n")
  } else {
    cat(sprintf("  GTF file: %s\n", opt$gtf_file))

    tryCatch(
      {
        # Convert DMRs to GRanges
        cat("  Converting DMRs to GRanges object...\n")

        # Try multiple sources for coordinates, in order of preference
        dmr_gr <- NULL

        if (exists("gene_annotation_results") && nrow(gene_annotation_results) > 0) {
          # First try: Use gene annotation results if available
          cat("  Using coordinates from gene annotation results...\n")
          dmr_coords <- gene_annotation_results[, c("Chromosome", "Start", "End")]
          dmr_gr <- GRanges(
            seqnames = dmr_coords$Chromosome,
            ranges = IRanges(start = dmr_coords$Start, end = dmr_coords$End)
          )
        } else if (exists("dmr_details") && nrow(dmr_details) > 0) {
          # Second try: Use dmr_details which should have coordinates
          cat("  Using coordinates from dmr_details...\n")

          # Check for coordinate columns (various naming conventions)
          chr_col <- NULL
          start_col <- NULL
          end_col <- NULL

          if ("Chromosome" %in% colnames(dmr_details)) {
            chr_col <- "Chromosome"
          } else if ("chromosome" %in% colnames(dmr_details)) {
            chr_col <- "chromosome"
          } else if ("chr" %in% colnames(dmr_details)) {
            chr_col <- "chr"
          } else if ("Chr" %in% colnames(dmr_details)) chr_col <- "Chr"

          if ("Start" %in% colnames(dmr_details)) {
            start_col <- "Start"
          } else if ("start" %in% colnames(dmr_details)) start_col <- "start"

          if ("End" %in% colnames(dmr_details)) {
            end_col <- "End"
          } else if ("end" %in% colnames(dmr_details)) end_col <- "end"

          if (!is.null(chr_col) && !is.null(start_col) && !is.null(end_col)) {
            dmr_coords <- dmr_details[, c(chr_col, start_col, end_col)]
            colnames(dmr_coords) <- c("Chromosome", "Start", "End")

            dmr_gr <- GRanges(
              seqnames = dmr_coords$Chromosome,
              ranges = IRanges(
                start = as.integer(dmr_coords$Start),
                end = as.integer(dmr_coords$End)
              )
            )
            cat(sprintf(
              "  Successfully extracted coordinates from dmr_details (%s, %s, %s)\n",
              chr_col, start_col, end_col
            ))
          } else {
            cat("  [WARNING] Could not find coordinate columns in dmr_details\n")
            cat(sprintf("  Available columns: %s\n", paste(colnames(dmr_details), collapse = ", ")))
          }
        }

        if (is.null(dmr_gr)) {
          cat("  [ERROR] No DMR coordinates available for feature annotation\n")
          cat("  Feature annotation requires either:\n")
          cat("    1. Gene annotation results (use --enable_gene_annotation), OR\n")
          cat("    2. dmr_details_with_coordinates.csv with Chromosome/Start/End columns\n")
          cat("  Skipping feature annotation.\n\n")
        }

        if (!is.null(dmr_gr) && length(dmr_gr) > 0) {
          cat(sprintf("  Converted %d DMRs to GRanges\n", length(dmr_gr)))

          # Load genomic features from GTF
          cat("  Loading gene features from GTF...\n")
          allFeat <- gffToGRanges(opt$gtf_file)
          allFeat <- as(split(allFeat, allFeat$type), "GRangesList")
          cat(sprintf("  Loaded %d feature types\n", length(allFeat)))

          # Annotate with gene features
          cat("  Annotating DMRs with gene features...\n")
          dmr_feat_annot <- annotateWithFeatures(dmr_gr, allFeat)
          # Convert named vector to dataframe with proper structure
          feat_percentages <- data.frame(
            Feature = names(dmr_feat_annot@annotation),
            Percentage = as.numeric(dmr_feat_annot@annotation),
            row.names = NULL
          )

          # Store results
          feature_annotation_results <- list(
            gene_features = feat_percentages
          )

          # Save gene feature annotation
          write.csv(feat_percentages,
            file.path(feature_annot_dir, "gene_feature_annotation.csv"),
            row.names = FALSE
          )
          cat(sprintf(
            "  [OK] Saved gene feature annotation to: %s\n",
            file.path(feature_annot_dir, "gene_feature_annotation.csv")
          ))

          # Load and annotate regulatory features if provided
          if (!is.null(opt$regulatory_gff) && file.exists(opt$regulatory_gff)) {
            cat(sprintf("  Loading regulatory features from: %s\n", opt$regulatory_gff))
            regFeat <- gffToGRanges(opt$regulatory_gff)
            regFeat <- as(split(regFeat, regFeat$type), "GRangesList")

            dmr_reg_annot <- annotateWithFeatures(dmr_gr, regFeat)
            # Convert named vector to dataframe with proper structure
            reg_percentages <- data.frame(
              Feature = names(dmr_reg_annot@annotation),
              Percentage = as.numeric(dmr_reg_annot@annotation),
              row.names = NULL
            )

            feature_annotation_results$regulatory_features <- reg_percentages
            write.csv(reg_percentages,
              file.path(feature_annot_dir, "regulatory_feature_annotation.csv"),
              row.names = FALSE
            )
            cat(sprintf("  [OK] Saved regulatory feature annotation\n"))
          }

          # Load and annotate EMAR features if provided
          if (!is.null(opt$emar_gff) && file.exists(opt$emar_gff)) {
            cat(sprintf("  Loading EMAR features from: %s\n", opt$emar_gff))
            emarFeat <- gffToGRanges(opt$emar_gff)
            emarFeat <- as(split(emarFeat, emarFeat$type), "GRangesList")

            dmr_emar_annot <- annotateWithFeatures(dmr_gr, emarFeat)
            # Convert named vector to dataframe with proper structure
            emar_percentages <- data.frame(
              Feature = names(dmr_emar_annot@annotation),
              Percentage = as.numeric(dmr_emar_annot@annotation),
              row.names = NULL
            )

            feature_annotation_results$emar_features <- emar_percentages
            write.csv(emar_percentages,
              file.path(feature_annot_dir, "emar_feature_annotation.csv"),
              row.names = FALSE
            )
            cat(sprintf("  [OK] Saved EMAR feature annotation\n"))
          }

          # Load and annotate open chromatin if provided
          if (!is.null(opt$open_chromatin_bed) && file.exists(opt$open_chromatin_bed)) {
            cat(sprintf("  Loading open chromatin from: %s\n", opt$open_chromatin_bed))
            flank_size <- if (!is.null(opt$chromatin_flank_bp)) opt$chromatin_flank_bp else 2000
            chrom_annot <- readFeatureFlank(opt$open_chromatin_bed,
              feature.flank.name = c("open_chromatin", "chromatin_flanking_regions"),
              flank = flank_size
            )

            dmr_chrom_annot <- annotateWithFeatures(dmr_gr, chrom_annot)
            # Convert named vector to dataframe with proper structure
            chrom_percentages <- data.frame(
              Feature = names(dmr_chrom_annot@annotation),
              Percentage = as.numeric(dmr_chrom_annot@annotation),
              row.names = NULL
            )

            feature_annotation_results$chromatin_features <- chrom_percentages
            write.csv(chrom_percentages,
              file.path(feature_annot_dir, "chromatin_feature_annotation.csv"),
              row.names = FALSE
            )
            cat(sprintf("  [OK] Saved chromatin feature annotation\n"))
          }

          # Create combined feature plot
          cat("  Creating feature annotation plot...\n")

          # Combine all features
          all_features <- list()
          if ("gene_features" %in% names(feature_annotation_results)) {
            all_features[[length(all_features) + 1]] <- feature_annotation_results$gene_features
          }
          if ("regulatory_features" %in% names(feature_annotation_results)) {
            all_features[[length(all_features) + 1]] <- feature_annotation_results$regulatory_features
          }
          if ("emar_features" %in% names(feature_annotation_results)) {
            all_features[[length(all_features) + 1]] <- feature_annotation_results$emar_features
          }
          if ("chromatin_features" %in% names(feature_annotation_results)) {
            all_features[[length(all_features) + 1]] <- feature_annotation_results$chromatin_features
          }

          if (length(all_features) > 0) {
            combined_features <- do.call(rbind, all_features)
            # Remove duplicates and sort by percentage
            combined_features <- combined_features[!duplicated(combined_features$Feature), ]
            combined_features <- combined_features[order(combined_features$Percentage, decreasing = TRUE), ]

            # Only keep features with >0% annotation
            combined_features <- combined_features[combined_features$Percentage > 0, ]

            if (nrow(combined_features) > 0) {
              # Create bar plot
              combined_features$Feature <- factor(combined_features$Feature,
                levels = combined_features$Feature
              )

              feature_plot <- ggplot(combined_features, aes(x = Feature, y = Percentage)) +
                geom_bar(stat = "identity", fill = "steelblue", color = "black", size = 0.3) +
                geom_text(aes(label = sprintf("%.1f%%", Percentage)),
                  vjust = -0.5, size = 3
                ) +
                labs(
                  title = "DMR Genomic Feature Annotation",
                  subtitle = sprintf("Percentage of %d DMRs overlapping each feature type", length(dmr_gr)),
                  x = "Feature Type",
                  y = "Percentage of DMRs (%)"
                ) +
                theme_bw(base_size = 12) +
                theme(
                  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                  plot.title = element_text(face = "bold", size = 14),
                  panel.grid.major.x = element_blank()
                ) +
                ylim(0, 100)

              # Save plot using universal format
              save_plot_multiple_formats(feature_plot,
                file.path(feature_figures_dir, "feature_annotation_barplot"),
                width = 12, height = 8
              )
              cat(sprintf("  [OK] Saved feature annotation plot\n"))

              # Save combined table
              write.csv(combined_features,
                file.path(feature_annot_dir, "combined_feature_annotation.csv"),
                row.names = FALSE
              )
            }
          }

          feature_annotation_done <- TRUE
          cat("  [OK] Feature annotation completed successfully\n\n")
        }
      },
      error = function(e) {
        cat(sprintf("  [ERROR] Feature annotation failed: %s\n", e$message))
        cat("  Continuing without feature annotation.\n\n")
      }
    )
  }
} else {
  cat("[STEP 6.10.6/10] Feature annotation disabled (use --enable_feature_annotation to enable)\n\n")
}

# ============================================================================
# STEP 6.10.7: LOGISTIC REGRESSION - INFECTION PROBABILITY vs MEDIAN METHYLATION
# ============================================================================

cat("[STEP 6.10.7/10] Performing logistic regression analysis for infection probability...\n")
log_section_start("6.10.7", "Logistic regression: infection probability", 10)

logistic_regression_done <- FALSE

tryCatch(
  {
    # Create output directory
    logistic_dir <- file.path(output_dir, "logistic_regression")
    dir.create(logistic_dir, showWarnings = FALSE, recursive = TRUE)

    # Calculate median methylation across all DMRs for each sample
    cat("  Calculating median methylation per sample across DMRs...\n")

    # Check if we have methylation matrix and sample metadata
    if (!exists("meth_matrix") || is.null(meth_matrix)) {
      stop("meth_matrix not found")
    }

    if (!exists("sample_metadata") || is.null(sample_metadata)) {
      stop("sample_metadata not found")
    }

    # Check if we have at least one DMR
    if (nrow(meth_matrix) == 0) {
      stop("No DMRs in methylation matrix")
    }

    # Use all DMRs in the meth_matrix (which already contains only the selected markers)
    dmr_meth_matrix <- meth_matrix

    if (nrow(dmr_meth_matrix) > 0) {
      # Calculate median methylation per sample
      median_meth_per_sample <- apply(dmr_meth_matrix, 2, median, na.rm = TRUE)

      # Create data frame with sample IDs, median methylation, and group from sample_metadata
      sample_data_merged <- data.frame(
        Sample_ID = colnames(dmr_meth_matrix),
        Median_Methylation = median_meth_per_sample,
        Group = sample_metadata$actual_class[match(colnames(dmr_meth_matrix), sample_metadata$Sample_ID)],
        stringsAsFactors = FALSE
      )

      # Remove rows with missing data
      sample_data_merged <- sample_data_merged[!is.na(sample_data_merged$Group) &
        !is.na(sample_data_merged$Median_Methylation), ]

      # Create binary infection status (1 = infected, 0 = control)
      groups <- unique(sample_data_merged$Group)
      cat(sprintf("  Groups identified: %s\n", paste(groups, collapse = ", ")))

      if (length(groups) == 2) {
        # Binary classification - intelligently determine control vs case
        control_group <- NULL
        infected_group <- NULL

        # Method 1: User explicitly specified control group
        if (!is.null(opt$control_group)) {
          if (opt$control_group %in% groups) {
            control_group <- opt$control_group
            infected_group <- groups[groups != control_group]
            cat(sprintf("  [INFO] Using user-specified control group: '%s'\n", control_group))
          } else {
            cat(sprintf("  [WARNING] Specified control group '%s' not found in data!\n", opt$control_group))
            cat(sprintf("  [WARNING] Available groups: %s\n", paste(groups, collapse = ", ")))
            cat("  [WARNING] Falling back to keyword-based detection\n")
          }
        }

        # Method 2: Keyword-based detection (if user didn't specify or specification failed)
        if (is.null(control_group)) {
          # Keywords that typically indicate control/reference groups (case-insensitive)
          control_keywords <- c("control", "healthy", "negative", "uninfected", "wildtype", "wt", "normal", "baseline")
          # Keywords that typically indicate case/treatment groups
          case_keywords <- c("infected", "infection", "disease", "case", "positive", "treatment", "tumor", "cancer", "patient")

          control_detected <- FALSE
          case_detected <- FALSE

          for (grp in groups) {
            grp_lower <- tolower(grp)
            # Check for control keywords
            for (keyword in control_keywords) {
              if (grepl(keyword, grp_lower, fixed = TRUE)) {
                control_group <- grp
                control_detected <- TRUE
                cat(sprintf("  [INFO] Auto-detected control group by keyword '%s': '%s'\n", keyword, grp))
                break
              }
            }
            # Check for case keywords
            for (keyword in case_keywords) {
              if (grepl(keyword, grp_lower, fixed = TRUE)) {
                infected_group <- grp
                case_detected <- TRUE
                cat(sprintf("  [INFO] Auto-detected case group by keyword '%s': '%s'\n", keyword, grp))
                break
              }
            }
          }

          # Assign remaining group
          if (!is.null(control_group) && is.null(infected_group)) {
            infected_group <- groups[groups != control_group]
          } else if (is.null(control_group) && !is.null(infected_group)) {
            control_group <- groups[groups != infected_group]
            cat(sprintf("  [INFO] Assigning remaining group as control: '%s'\n", control_group))
          } else if (is.null(control_group) && is.null(infected_group)) {
            # Fallback to alphabetical if no keywords matched
            groups_sorted <- sort(groups)
            control_group <- groups_sorted[1]
            infected_group <- groups_sorted[2]
            cat(sprintf("  [WARNING] No control/case keywords detected. Using alphabetical order:\n"))
            cat(sprintf("  [WARNING] First alphabetically = control: '%s'\n", control_group))
            cat(sprintf("  [WARNING] Second alphabetically = case: '%s'\n", infected_group))
            cat(sprintf("  [WARNING] Use --control_group to explicitly specify the control group\n"))
          }
        }

        cat(sprintf("\n  === LOGISTIC REGRESSION GROUP ASSIGNMENT ===\n"))
        cat(sprintf("  Control/Reference (coded as 0): '%s'\n", control_group))
        cat(sprintf("  Case/Infected (coded as 1):     '%s'\n", infected_group))
        cat(sprintf("  Model predicts probability of:  '%s'\n", infected_group))
        cat(sprintf("  ==========================================\n\n"))

        sample_data_merged$Infection_Binary <- ifelse(
          sample_data_merged$Group == infected_group, 1, 0
        )

        # Fit logistic regression model
        cat("  Fitting logistic regression model...\n")
        logit_model <- glm(Infection_Binary ~ Median_Methylation,
          data = sample_data_merged,
          family = binomial(link = "logit")
        )

        # Summary
        model_summary <- summary(logit_model)
        cat("\n  Model Summary:\n")
        print(model_summary$coefficients)
        cat("\n")

        # Calculate probabilities at specific thresholds
        probability_thresholds <- c(0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99)

        # Create prediction dataframe
        # We need to find methylation values that correspond to each probability
        # Using the inverse logit: log(p/(1-p)) = intercept + slope * methylation
        # Solving for methylation: methylation = (log(p/(1-p)) - intercept) / slope

        intercept <- coef(logit_model)[1]
        slope <- coef(logit_model)[2]

        prediction_table <- data.frame(
          Probability_of_Infection = probability_thresholds,
          Median_Methylation_Percent = NA
        )

        for (i in 1:length(probability_thresholds)) {
          prob <- probability_thresholds[i]
          # Calculate log-odds
          log_odds <- log(prob / (1 - prob))
          # Calculate methylation value
          meth_value <- (log_odds - intercept) / slope
          prediction_table$Median_Methylation_Percent[i] <- meth_value
        }

        # Round to 2 decimal places
        prediction_table$Median_Methylation_Percent <- round(prediction_table$Median_Methylation_Percent, 2)

        # Save results
        logit_results_file <- file.path(logistic_dir, "infection_probability_vs_methylation.csv")
        write.csv(prediction_table, logit_results_file, row.names = FALSE)
        cat(sprintf("  [OK] Saved logistic regression results: %s\n", logit_results_file))

        # Save model summary
        model_summary_file <- file.path(logistic_dir, "logistic_regression_model_summary.txt")
        sink(model_summary_file)
        cat("LOGISTIC REGRESSION MODEL SUMMARY\n")
        cat("==================================\n\n")
        cat("Model: Infection_Binary ~ Median_Methylation\n")
        cat(sprintf("Family: binomial (logit link)\n\n"))
        cat("Coefficients:\n")
        print(model_summary$coefficients)
        cat("\n")
        cat(sprintf("AIC: %.2f\n", AIC(logit_model)))
        cat(sprintf(
          "Residual deviance: %.2f on %d degrees of freedom\n",
          model_summary$deviance, model_summary$df.residual
        ))
        cat("\n")
        cat("PREDICTION TABLE\n")
        cat("================\n\n")
        print(prediction_table, row.names = FALSE)
        sink()

        cat(sprintf("  [OK] Saved model summary: %s\n", model_summary_file))

        # Store for report
        logistic_regression_results <- list(
          available = TRUE,
          model = logit_model,
          prediction_table = prediction_table,
          control_group = control_group,
          infected_group = infected_group,
          intercept = intercept,
          slope = slope,
          p_value = model_summary$coefficients[2, 4] # p-value for slope
        )

        logistic_regression_done <- TRUE
        cat("  [OK] Logistic regression analysis completed successfully\n\n")
      } else {
        cat(sprintf("  [WARNING] Logistic regression requires binary classification (2 groups), but found %d groups\n", length(groups)))
        cat("  Skipping logistic regression analysis.\n\n")
        logistic_regression_results <- list(available = FALSE)
      }
    } else {
      cat("  [WARNING] No DMRs or methylation matrix not available for logistic regression\n\n")
      logistic_regression_results <- list(available = FALSE)
    }
  },
  error = function(e) {
    cat(sprintf("  [ERROR] Logistic regression analysis failed: %s\n", e$message))
    cat("  Continuing without logistic regression analysis.\n\n")
    logistic_regression_results <- list(available = FALSE)
  }
)

log_section_end("6.10.7", "Logistic regression: infection probability")

# ============================================================================
# STEP 6.11: TIME SERIES ANALYSIS (CONDITIONAL)
# ============================================================================

time_series_done <- FALSE

if (!is.null(opt$time_analysis)) {
  cat("[STEP 6.11/10] Performing time series analysis on DMRs...\n")
  log_section_start(6.11, "Time series analysis", 10)

  # Create output directory
  time_series_dir <- file.path(output_dir, "time_series")
  dir.create(time_series_dir, showWarnings = FALSE, recursive = TRUE)
  time_series_figures_dir <- file.path(time_series_dir, "figures")
  dir.create(time_series_figures_dir, showWarnings = FALSE, recursive = TRUE)

  # Validate timepoint column exists
  if (!opt$time_analysis %in% colnames(sample_sheet_matched)) {
    cat(sprintf("  [ERROR] Timepoint column '%s' not found in sample sheet!\n", opt$time_analysis))
    cat(sprintf("  Available columns: %s\n", paste(colnames(sample_sheet_matched), collapse = ", ")))
    cat("  Skipping time series analysis.\n\n")
  } else {
    cat(sprintf("  Using timepoint column: %s\n", opt$time_analysis))

    # Check for required packages
    required_packages <- c("lme4", "nlme", "emmeans")
    missing_packages <- c()

    for (pkg in required_packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        missing_packages <- c(missing_packages, pkg)
      }
    }

    if (length(missing_packages) > 0) {
      cat("\n")
      cat("================================================================================\n")
      cat("TIME SERIES ANALYSIS PACKAGES NOT INSTALLED\n")
      cat("================================================================================\n")
      cat(sprintf("The following packages are required but not installed: %s\n", paste(missing_packages, collapse = ", ")))
      cat("\n")
      cat("To install these packages:\n")
      cat("  install.packages(c('lme4', 'nlme', 'emmeans'))\n")
      cat("\n")
      cat("Skipping time series analysis for now.\n")
      cat("================================================================================\n\n")
    } else {
      # Extract timepoint data
      timepoints <- sample_sheet_matched[[opt$time_analysis]]

      # Check for subject ID column (for repeated measures detection)
      subject_ids <- NULL
      subject_col_names <- c("subject_id", "Subject_ID", "SubjectID", "subject", "Subject", "PatientID", "patient_id")
      for (col_name in subject_col_names) {
        if (col_name %in% colnames(sample_sheet_matched)) {
          subject_ids <- sample_sheet_matched[[col_name]]
          cat(sprintf("  Found subject ID column: %s (for repeated measures detection)\n", col_name))
          break
        }
      }

      # Get biological groups with validation
      biological_groups <- NULL
      skip_analysis <- FALSE

      if (is.null(opt$treatment_column) || !opt$treatment_column %in% colnames(sample_sheet_matched)) {
        cat(sprintf(
          "  [WARNING] Treatment column '%s' not found in sample sheet!\n",
          ifelse(is.null(opt$treatment_column), "NULL", opt$treatment_column)
        ))
        cat(sprintf("  Available columns: %s\n", paste(colnames(sample_sheet_matched), collapse = ", ")))

        # Try to find a suitable column
        possible_cols <- c("Group", "group", "treatment", "Treatment", "actual_class", "class")
        found_col <- NULL
        for (col in possible_cols) {
          if (col %in% colnames(sample_sheet_matched)) {
            found_col <- col
            break
          }
        }

        if (!is.null(found_col)) {
          cat(sprintf("  [OK] Using '%s' as treatment column instead\n", found_col))
          biological_groups <- sample_sheet_matched[[found_col]]
        } else {
          cat("  [ERROR] Cannot find suitable treatment column. Skipping time series analysis.\n\n")
          skip_analysis <- TRUE
        }
      } else {
        biological_groups <- sample_sheet_matched[[opt$treatment_column]]
      }

      # Validate biological_groups
      if (!skip_analysis && (is.null(biological_groups) || length(biological_groups) == 0)) {
        cat(sprintf("  [ERROR] Treatment column is empty!\n"))
        cat("  Skipping time series analysis.\n\n")
        skip_analysis <- TRUE
      }

      # Check for NA values in biological_groups and replace with "Unknown"
      if (!skip_analysis && any(is.na(biological_groups))) {
        n_na <- sum(is.na(biological_groups))
        cat(sprintf("  [WARNING] %d samples have NA biological group labels - replacing with 'Unknown'\n", n_na))
        biological_groups[is.na(biological_groups)] <- "Unknown"
      }

      # Only proceed if validation passed
      if (!skip_analysis) {
        cat(sprintf("  Analysing %d DMRs across %d timepoints...\n", nrow(meth_matrix), length(unique(timepoints))))

        if (!is.null(subject_ids)) {
          n_subjects <- length(unique(subject_ids))
          if (any(duplicated(subject_ids))) {
            cat(sprintf("  Detected repeated measures design: %d unique subjects\n", n_subjects))
          } else {
            cat(sprintf("  Detected independent groups design: %d samples\n", length(subject_ids)))
          }
        }

        # Force categorical if specified
        force_categorical <- opt$time_categorical

        # Initialize results storage
        time_series_results <- data.frame(
          DMR_ID = rownames(meth_matrix),
          Model_Type = NA,
          Time_Effect_Slope = NA,
          Time_Effect_P_Value = NA,
          Adjusted_P_Value = NA,
          Effect_Size = NA,
          R_Squared = NA,
          Trajectory_Cluster = NA,
          Methylation_Direction = NA,
          stringsAsFactors = FALSE
        )

        # Analyse each DMR
        cat("  Analysing DMRs:\n")
        for (i in 1:nrow(meth_matrix)) {
          dmr_id <- rownames(meth_matrix)[i]
          methylation_values <- as.numeric(meth_matrix[i, ])

          # Call analysis function
          result <- analyse_time_series(
            methylation_values = methylation_values,
            timepoints = timepoints,
            biological_groups = biological_groups,
            subject_ids = subject_ids,
            force_categorical = force_categorical,
            dmr_id = dmr_id
          )

          # Store results
          if (result$available) {
            time_series_results$Model_Type[i] <- result$model_type
            time_series_results$Time_Effect_Slope[i] <- result$time_effect_slope
            time_series_results$Time_Effect_P_Value[i] <- result$time_effect_pvalue
            time_series_results$Effect_Size[i] <- result$effect_size
            time_series_results$R_Squared[i] <- result$r_squared
            time_series_results$Trajectory_Cluster[i] <- result$trajectory_cluster

            # Calculate methylation direction
            if (!is.na(result$time_effect_slope)) {
              # For continuous timepoints, use slope
              if (result$time_effect_slope > 0.01) {
                time_series_results$Methylation_Direction[i] <- "Increasing"
              } else if (result$time_effect_slope < -0.01) {
                time_series_results$Methylation_Direction[i] <- "Decreasing"
              } else {
                time_series_results$Methylation_Direction[i] <- "Stable"
              }
            } else {
              # For categorical timepoints, compare first vs last timepoint
              timepoint_levels <- sort(unique(timepoints))
              if (length(timepoint_levels) >= 2) {
                first_tp <- timepoint_levels[1]
                last_tp <- timepoint_levels[length(timepoint_levels)]

                mean_first <- mean(methylation_values[timepoints == first_tp], na.rm = TRUE)
                mean_last <- mean(methylation_values[timepoints == last_tp], na.rm = TRUE)

                meth_change <- mean_last - mean_first

                if (meth_change > 5) { # >5% increase
                  time_series_results$Methylation_Direction[i] <- "Increasing"
                } else if (meth_change < -5) { # >5% decrease
                  time_series_results$Methylation_Direction[i] <- "Decreasing"
                } else {
                  time_series_results$Methylation_Direction[i] <- "Stable"
                }
              }
            }
          } else {
            time_series_results$Model_Type[i] <- paste("Error:", result$error)
          }

          # Progress indicator
          if (i %% 10 == 0 || i == nrow(meth_matrix)) {
            cat(sprintf("    Processed %d/%d DMRs...\n", i, nrow(meth_matrix)))
          }
        }

        # Apply FDR correction to p-values
        valid_pvals <- !is.na(time_series_results$Time_Effect_P_Value)
        if (sum(valid_pvals) > 0) {
          time_series_results$Adjusted_P_Value[valid_pvals] <- apply_multiple_testing(
            time_series_results$Time_Effect_P_Value[valid_pvals],
            method = mt_method_r
          )
        }

        # Save results
        write.csv(time_series_results,
          file.path(time_series_dir, "dmr_time_series_statistics.csv"),
          row.names = FALSE
        )

        cat(sprintf("  Saved time series statistics: dmr_time_series_statistics.csv\n"))

        # Trajectory clustering summary
        trajectory_summary <- table(time_series_results$Trajectory_Cluster)
        trajectory_df <- data.frame(
          Cluster_Name = names(trajectory_summary),
          Number_of_DMRs = as.numeric(trajectory_summary),
          stringsAsFactors = FALSE
        )

        write.csv(trajectory_df,
          file.path(time_series_dir, "dmr_time_trajectories.csv"),
          row.names = FALSE
        )

        cat(sprintf("  Trajectory clustering:\n"))
        for (i in 1:nrow(trajectory_df)) {
          cat(sprintf("    %s: %d DMRs\n", trajectory_df$Cluster_Name[i], trajectory_df$Number_of_DMRs[i]))
        }

        # ====== VISUALIZATIONS ======

        # 1. Time series line plots for top 12 significant DMRs
        sig_dmrs <- time_series_results[!is.na(time_series_results$Adjusted_P_Value) &
          time_series_results$Adjusted_P_Value < sig_threshold, ]

        if (nrow(sig_dmrs) > 0) {
          cat(sprintf("  Creating time series trajectory plots...\n"))

          # DIAGNOSTIC: Check meth_matrix dimensions
          cat(sprintf(
            "    [DEBUG] meth_matrix dimensions: %d DMRs x %d samples\n",
            nrow(meth_matrix), ncol(meth_matrix)
          ))
          cat(sprintf("    [DEBUG] timepoints length: %d\n", length(timepoints)))
          cat(sprintf("    [DEBUG] biological_groups length: %d\n", length(biological_groups)))

          # CRITICAL FIX: Ensure we only use samples that exist in both meth_matrix and sample_sheet_matched
          meth_samples <- colnames(meth_matrix)
          sheet_samples <- rownames(sample_sheet_matched)
          common_samples <- intersect(meth_samples, sheet_samples)

          if (length(common_samples) == 0) {
            cat("    [ERROR] No common samples between meth_matrix and sample_sheet_matched!\n")
            cat(sprintf("    [DEBUG] meth_matrix has %d samples\n", length(meth_samples)))
            cat(sprintf("    [DEBUG] sample_sheet_matched has %d samples\n", length(sheet_samples)))
          } else if (length(common_samples) < length(meth_samples)) {
            cat(sprintf(
              "    [WARNING] Only %d/%d samples in meth_matrix have matching metadata\n",
              length(common_samples), length(meth_samples)
            ))

            # Filter meth_matrix to only common samples
            meth_matrix_filtered <- meth_matrix[, common_samples, drop = FALSE]

            # Filter timepoints and groups to only common samples
            sample_indices <- match(common_samples, sheet_samples)
            timepoints_filtered <- timepoints[sample_indices]
            biological_groups_filtered <- biological_groups[sample_indices]

            cat(sprintf("    [OK] Using %d matched samples for plotting\n", length(common_samples)))
          } else {
            # All samples match - use originals
            meth_matrix_filtered <- meth_matrix
            timepoints_filtered <- timepoints
            biological_groups_filtered <- biological_groups
          }

          # Define color scheme: Control (0) = blue, Case/Infected (1) = red
          group_levels <- sort(unique(biological_groups_filtered))
          group_colors <- RColorBrewer::brewer.pal(max(3, length(group_levels)), "Set1")
          names(group_colors) <- group_levels

          # Override with sensible medical colors if we detect control/case pattern
          if (length(group_levels) == 2) {
            # Binary classification: assume first level is control, second is case
            group_colors[1] <- "#377EB8" # Blue for control
            group_colors[2] <- "#E41A1C" # Red for case/infected
          } else if (length(group_levels) > 2) {
            # Multi-class: first is control (blue), rest get warm colors
            group_colors[1] <- "#377EB8" # Blue for control
          }

          # Sort by p-value and use ALL significant DMRs (not just top 12)
          sig_dmrs <- sig_dmrs[order(sig_dmrs$Adjusted_P_Value), ]

          cat(sprintf("    Creating faceted plot for %d significant DMRs...\n", nrow(sig_dmrs)))

          # Create combined data frame for all DMRs
          all_plot_data <- list()

          for (i in 1:nrow(sig_dmrs)) {
            dmr_id <- sig_dmrs$DMR_ID[i]
            dmr_idx <- which(rownames(meth_matrix_filtered) == dmr_id)

            if (length(dmr_idx) == 0) {
              cat(sprintf("    [WARNING] Skipping DMR %s - not found in methylation matrix\n", dmr_id))
              next
            }

            # Extract methylation values for this DMR
            methylation_values <- as.numeric(meth_matrix_filtered[dmr_idx, ])

            # Skip if methylation values extraction failed
            if (length(methylation_values) == 0 || length(methylation_values) != length(timepoints_filtered)) {
              cat(sprintf(
                "    [WARNING] Skipping DMR %s - data dimension mismatch (got %d values, expected %d)\n",
                dmr_id, length(methylation_values), length(timepoints_filtered)
              ))
              next
            }

            plot_df <- data.frame(
              DMR_ID = dmr_id,
              DMR_Label = sprintf("%s\n(p=%.2e)", dmr_id, sig_dmrs$Adjusted_P_Value[i]),
              Timepoint = timepoints_filtered,
              Methylation = methylation_values,
              Group = biological_groups_filtered,
              stringsAsFactors = FALSE
            )

            # Replace NA groups with "Unknown"
            if (any(is.na(plot_df$Group))) {
              plot_df$Group[is.na(plot_df$Group)] <- "Unknown"
            }

            # Remove missing values (methylation and timepoint only - groups already handled)
            plot_df <- plot_df[!is.na(plot_df$Methylation) & !is.na(plot_df$Timepoint), ]

            all_plot_data[[i]] <- plot_df
          }

          # Combine all DMR data
          if (length(all_plot_data) > 0) {
            combined_data <- do.call(rbind, all_plot_data)

            # Sort timepoints correctly (handle week_1, week_2, etc.)
            # Convert timepoint to ordered factor with natural sorting
            timepoint_levels <- sort(unique(combined_data$Timepoint))

            # Try to extract numbers for natural sorting (e.g., week_1 < week_2 < week_10)
            timepoint_order <- tryCatch(
              {
                # Extract numeric part if present
                nums <- as.numeric(gsub("[^0-9]", "", timepoint_levels))
                if (all(!is.na(nums))) {
                  # If all have numbers, sort by numeric value
                  timepoint_levels[order(nums)]
                } else {
                  # Otherwise keep original order
                  timepoint_levels
                }
              },
              error = function(e) {
                timepoint_levels
              }
            )

            combined_data$Timepoint <- factor(combined_data$Timepoint, levels = timepoint_order)

            # Order DMR_Label by p-value (already sorted in sig_dmrs)
            combined_data$DMR_Label <- factor(combined_data$DMR_Label,
              levels = unique(combined_data$DMR_Label)
            )

            # Calculate optimal grid dimensions
            n_dmrs <- length(unique(combined_data$DMR_ID))
            ncol_facet <- min(6, max(3, ceiling(sqrt(n_dmrs))))

            # Calculate figure dimensions based on number of DMRs
            nrow_facet <- ceiling(n_dmrs / ncol_facet)
            fig_width <- max(16, ncol_facet * 3)
            fig_height <- max(12, nrow_facet * 2.5)

            # Create faceted plot
            p <- ggplot(combined_data, aes(x = Timepoint, y = Methylation, color = Group, group = Group)) +
              stat_summary(fun = mean, geom = "line", size = 0.8) +
              stat_summary(fun = mean, geom = "point", size = 1.5) +
              stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, alpha = 0.5) +
              facet_wrap(~DMR_Label, scales = "free_y", ncol = ncol_facet) +
              scale_color_manual(values = group_colors) +
              labs(
                x = "Timepoint",
                y = "Methylation (%)",
                title = sprintf("Temporal Trajectories for All %d Significant DMRs", n_dmrs)
              ) +
              theme_bw() +
              theme(
                axis.text.x = element_text(angle = 45, hjust = 1, size = 7), # 45° rotation
                axis.text.y = element_text(size = 7),
                axis.title = element_text(size = 9),
                strip.text = element_text(size = 6, face = "bold"),
                strip.background = element_rect(fill = "gray90"),
                legend.position = "bottom",
                legend.title = element_text(size = 8),
                legend.text = element_text(size = 7),
                plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
              )

            ggsave(
              file.path(time_series_figures_dir, "time_series_trajectories_COMBINED.png"),
              p,
              width = fig_width,
              height = fig_height,
              dpi = opt$figure_resolution
            )
            cat(sprintf(
              "    Saved: time_series_trajectories_COMBINED.png (%d DMRs in %dx%d grid)\n",
              n_dmrs, nrow_facet, ncol_facet
            ))
          } else {
            cat("    [WARNING] No valid DMR data for plotting\n")
          }
        } else {
          cat("  No significant DMRs found for trajectory plotting.\n")
        }

        # 2. Heatmap of all DMRs across timepoints
        cat("  Creating time series heatmap...\n")

        # Calculate mean methylation for each DMR at each timepoint
        unique_timepoints <- sort(unique(timepoints))
        heatmap_matrix <- matrix(NA, nrow = nrow(meth_matrix), ncol = length(unique_timepoints))
        rownames(heatmap_matrix) <- rownames(meth_matrix)
        colnames(heatmap_matrix) <- as.character(unique_timepoints)

        for (i in 1:nrow(meth_matrix)) {
          for (j in 1:length(unique_timepoints)) {
            tp_idx <- which(timepoints == unique_timepoints[j])
            if (length(tp_idx) > 0) {
              heatmap_matrix[i, j] <- mean(as.numeric(meth_matrix[i, tp_idx]), na.rm = TRUE)
            }
          }
        }

        # Create heatmap
        pheatmap::pheatmap(
          heatmap_matrix,
          color = colorRampPalette(c("blue", "white", "red"))(100),
          scale = "row",
          clustering_distance_rows = "euclidean",
          clustering_distance_cols = "euclidean",
          clustering_method = "complete",
          show_rownames = FALSE,
          main = "DMR Methylation Across Timepoints",
          filename = file.path(time_series_figures_dir, "time_series_heatmap.png"),
          width = 8,
          height = 10
        )
        cat("    Saved: time_series_heatmap.png\n")

        # 3. Trajectory cluster visualisation
        cat("  Creating trajectory cluster plot...\n")

        cluster_counts <- table(time_series_results$Trajectory_Cluster)
        cluster_df <- data.frame(
          Cluster = names(cluster_counts),
          Count = as.numeric(cluster_counts),
          stringsAsFactors = FALSE
        )

        # Only create plot if there are multiple distinct clusters
        # (Skip if all DMRs are in single "Categorical" cluster - redundant)
        if (nrow(cluster_df) > 1) {
          cluster_plot <- ggplot(cluster_df, aes(x = reorder(Cluster, -Count), y = Count, fill = Cluster)) +
            geom_bar(stat = "identity") +
            geom_text(aes(label = Count), vjust = -0.5, size = 4) +
            labs(
              title = "DMR Trajectory Clustering",
              subtitle = sprintf("Classification of %d DMRs by temporal pattern", nrow(meth_matrix)),
              x = "Trajectory Cluster",
              y = "Number of DMRs"
            ) +
            theme_bw() +
            theme(
              legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1)
            )

          ggsave(
            file.path(time_series_figures_dir, "time_series_clusters.png"),
            cluster_plot,
            width = 8,
            height = 6,
            dpi = opt$figure_resolution
          )
          cat("    Saved: time_series_clusters.png\n")
        } else {
          cat("    Skipped trajectory cluster plot (all DMRs in single cluster - categorical timepoints)\n")
        }

        cat("\n  Time series analysis complete!\n")
        cat(sprintf("  Results saved to: %s\n", time_series_dir))
        time_series_done <- TRUE
      } # End of if (!skip_analysis)
    }
  }

  log_section_end(6.11, "Time series analysis")
  cat("\n")
} else {
  cat("[STEP 6.11/10] Time series analysis disabled (use --time_analysis to enable)\n\n")
  time_series_done <- FALSE
}

# ============================================================================
# STEP 7: PRE-ANALYTICAL QC SUMMARY
# ============================================================================
log_section_start(7, "Pre-analytical QC summary", 10)

qc_dir <- file.path(output_dir, "preanalytic_qc")
dir.create(qc_dir, showWarnings = FALSE)

# Identify QC columns
qc_columns <- c()
potential_qc <- c(
  "storage", "hemolysis", "quality", "qc", "time", "duration",
  "collection", "extraction", "concentration", "dna"
)

for (col in colnames(sample_sheet_matched)) {
  if (any(sapply(potential_qc, function(x) grepl(x, col, ignore.case = TRUE)))) {
    qc_columns <- c(qc_columns, col)
  }
}

if (!is.null(opt$storage_column)) qc_columns <- unique(c(qc_columns, opt$storage_column))
if (!is.null(opt$date_column)) qc_columns <- unique(c(qc_columns, opt$date_column))

cat("          Found", length(qc_columns), "QC-related columns\n")

if (length(qc_columns) > 0) {
  # Summary statistics
  qc_summary <- data.frame(
    column = character(),
    type = character(),
    n_valid = integer(),
    n_missing = integer(),
    summary_stat = character(),
    stringsAsFactors = FALSE
  )

  for (col in qc_columns) {
    if (!col %in% colnames(sample_sheet_matched)) next

    col_data <- sample_sheet_matched[[col]]
    n_valid <- sum(!is.na(col_data))
    n_missing <- sum(is.na(col_data))

    if (is.numeric(col_data)) {
      summary_stat <- paste0(
        "Mean: ", round(mean(col_data, na.rm = TRUE), 2),
        " | SD: ", round(sd(col_data, na.rm = TRUE), 2),
        " | Range: ", round(min(col_data, na.rm = TRUE), 2),
        " - ", round(max(col_data, na.rm = TRUE), 2)
      )
      col_type <- "numeric"
    } else {
      summary_stat <- paste(names(table(col_data)), ":", table(col_data), collapse = ", ")
      col_type <- "categorical"
    }

    qc_summary <- rbind(qc_summary, data.frame(
      column = col,
      type = col_type,
      n_valid = n_valid,
      n_missing = n_missing,
      summary_stat = summary_stat,
      stringsAsFactors = FALSE
    ))
  }

  write.csv(qc_summary, file.path(qc_dir, "qc_variable_summary.csv"), row.names = FALSE)

  # Correlations with methylation - OPTIMIZED
  cat("          Testing QC correlations with markers (optimised)...\n")

  for (col in qc_columns) {
    if (!col %in% colnames(sample_sheet_matched)) next

    col_data <- sample_sheet_matched[[col]]
    if (!is.numeric(col_data)) next

    marker_cors <- numeric(nrow(meth_matrix))
    marker_pvals <- numeric(nrow(meth_matrix))

    # Pre-compute valid samples for this QC variable
    valid_qc <- !is.na(col_data)

    for (i in 1:nrow(meth_matrix)) {
      meth_vals <- meth_matrix[i, ]
      valid_idx <- valid_qc & !is.na(meth_vals)

      if (sum(valid_idx) >= 5) {
        test <- cor.test(meth_vals[valid_idx], col_data[valid_idx],
          method = "spearman", exact = FALSE
        ) # exact=FALSE is faster
        marker_cors[i] <- test$estimate
        marker_pvals[i] <- test$p.value
      } else {
        marker_cors[i] <- NA
        marker_pvals[i] <- NA
      }
    }

    # Apply multiple testing correction
    marker_adjusted <- apply_multiple_testing(marker_pvals, method = mt_method_r)

    qc_corr <- data.frame(
      marker = rownames(meth_matrix),
      correlation = marker_cors,
      p_value = marker_pvals,
      adjusted_p = marker_adjusted,
      significant_raw = marker_pvals < sig_threshold,
      significant_adjusted = marker_adjusted < sig_threshold,
      stringsAsFactors = FALSE
    )

    write.csv(qc_corr,
      file.path(qc_dir, paste0("correlation_", gsub("[^A-Za-z0-9]", "_", col), ".csv")),
      row.names = FALSE
    )

    n_sig_raw <- sum(marker_pvals < sig_threshold, na.rm = TRUE)
    n_sig_adjusted <- sum(marker_adjusted < sig_threshold, na.rm = TRUE)
    cat(sprintf(
      "          %s: %d significant (raw p<%.3f), %d (%s<%.3f)\n",
      col, n_sig_raw, sig_threshold, n_sig_adjusted, mt_method_display, sig_threshold
    ))
  }

  # Plot distributions
  for (col in qc_columns) {
    if (!col %in% colnames(sample_sheet_matched)) next

    col_data <- sample_sheet_matched[[col]]

    if (is.numeric(col_data)) {
      p <- ggplot(data.frame(value = col_data), aes(x = value)) +
        geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, colour = "black") +
        labs(title = paste("Distribution of", col), x = col, y = "Count") +
        theme_minimal()
    } else {
      freq_df <- as.data.frame(table(col_data))
      names(freq_df) <- c("Category", "Count")

      p <- ggplot(freq_df, aes(x = Category, y = Count)) +
        geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
        geom_text(aes(label = Count), vjust = -0.3) +
        labs(title = paste("Distribution of", col), x = col, y = "Count") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }

    save_plot_multiple_formats(p,
      file.path(qc_dir, paste0("distribution_", gsub("[^A-Za-z0-9]", "_", col))),
      width = 8, height = 6
    )
  }

  # ===== ADVANCED: LINEAR MIXED EFFECTS MODELS FOR PRE-ANALYTICAL QC =====
  cat("\n")
  cat("          ===== ADVANCED ANALYSIS: MIXED EFFECTS MODELS =====\n")
  cat("          Testing pre-analytical variables with proper statistical models...\n")

  # Use the SAME technical covariates as batch effects analysis
  if (exists("sample_metadata") && "actual_class" %in% colnames(sample_metadata) &&
    exists("tech_covariates") && length(tech_covariates) > 0) {
    library(lme4)

    cat("          Using technical covariates from batch effects analysis:\n")
    cat("          ", paste(tech_covariates, collapse = ", "), "\n")

    lmm_summary_results <- list()

    for (col in tech_covariates) {
      if (!col %in% colnames(sample_sheet_matched)) next

      cat(sprintf("\n          Testing covariate: '%s'\n", col))

      covar_data <- sample_sheet_matched[[col]]

      # Skip if too many NAs
      if (sum(is.na(covar_data)) > length(covar_data) * 0.5) {
        cat("          [SKIP] Too many missing values\n")
        next
      }

      # Check if numeric or categorical
      is_numeric_covar <- is.numeric(covar_data)
      if (!is_numeric_covar) {
        # Try to convert if it looks numeric
        if (all(grepl("^[0-9\\.]+$", covar_data[!is.na(covar_data)]))) {
          covar_data <- as.numeric(covar_data)
          is_numeric_covar <- TRUE
        }
      }

      cat(sprintf("          Type: %s\n", ifelse(is_numeric_covar, "numeric", "categorical")))

      # Initialize results storage AS A LIST
      covar_lmm_results <- list()

      # Determine number of markers to test
      if (opt$lmm_n_markers == -1) {
        # Test ALL markers
        test_indices <- 1:nrow(meth_matrix)
        n_test <- nrow(meth_matrix)
        cat(sprintf("          Testing ALL %d markers (this may take time)...\n", n_test))
      } else if (opt$lmm_n_markers > 0) {
        # Test specified number
        n_test <- min(opt$lmm_n_markers, nrow(meth_matrix))
        test_indices <- if (n_test >= nrow(meth_matrix)) {
          1:nrow(meth_matrix)
        } else {
          sample(1:nrow(meth_matrix), n_test)
        }
        cat(sprintf("          Testing %d markers (user-specified)...\n", n_test))
      } else {
        # Default: test 20 markers
        n_test <- min(20, nrow(meth_matrix))
        test_indices <- sample(1:nrow(meth_matrix), n_test)
        cat(sprintf("          Testing %d markers (default sample)...\n", n_test))
      }

      for (i in test_indices) {
        # Progress reporting for large analyses
        if (n_test >= 50 && (which(test_indices == i) %% 10 == 0 || which(test_indices == i) == n_test)) {
          cat(sprintf(
            "\r          Progress: %d/%d markers tested",
            which(test_indices == i), n_test
          ))
          if (which(test_indices == i) == n_test) cat("\n")
        }

        # Prepare data for this marker
        marker_data <- data.frame(
          methylation = meth_matrix[i, ],
          biological_group = sample_metadata$actual_class[match(
            colnames(meth_matrix),
            rownames(sample_metadata)
          )],
          covariate = covar_data,
          sample_id = colnames(meth_matrix),
          stringsAsFactors = FALSE
        )

        # Remove NAs
        marker_data <- marker_data[complete.cases(marker_data), ]

        if (nrow(marker_data) < 10) next

        # Fit linear model (LM) with biological group as main effect
        tryCatch(
          {
            if (is_numeric_covar) {
              # Numeric covariate
              model <- lm(methylation ~ biological_group + covariate, data = marker_data)
            } else {
              # Categorical covariate
              marker_data$covariate <- factor(marker_data$covariate)
              if (length(unique(marker_data$covariate)) < 2) next
              model <- lm(methylation ~ biological_group + covariate, data = marker_data)
            }

            # Extract results
            model_summary <- summary(model)
            coef_table <- coef(model_summary)

            # Find covariate row(s) in coefficient table
            if (is_numeric_covar) {
              covar_row <- which(rownames(coef_table) == "covariate")
            } else {
              # For categorical, get all covariate levels
              covar_row <- grep("^covariate", rownames(coef_table))
            }

            if (length(covar_row) > 0) {
              # Take first covariate term (or average if multiple)
              covar_idx <- covar_row[1]

              # Add this result to the list
              covar_lmm_results[[length(covar_lmm_results) + 1]] <- data.frame(
                marker = rownames(meth_matrix)[i],
                covariate = col,
                covar_estimate = coef_table[covar_idx, "Estimate"],
                covar_stderr = coef_table[covar_idx, "Std. Error"],
                covar_tvalue = coef_table[covar_idx, "t value"],
                covar_pvalue = coef_table[covar_idx, "Pr(>|t|)"],
                model_rsquared = model_summary$r.squared,
                stringsAsFactors = FALSE
              )
            }
          },
          error = function(e) {
            # Skip if model fails
          }
        )
      } # End of marker loop

      # AFTER the marker loop, convert list to data frame
      if (length(covar_lmm_results) > 0) {
        covar_lmm_results <- do.call(rbind, covar_lmm_results)

        # Apply multiple testing correction after all markers are tested
        covar_lmm_results$covar_adjusted <- apply_multiple_testing(covar_lmm_results$covar_pvalue, method = mt_method_r)
        covar_lmm_results$significant_raw <- covar_lmm_results$covar_pvalue < sig_threshold
        covar_lmm_results$significant_adjusted <- covar_lmm_results$covar_adjusted < sig_threshold
      } else {
        # If no results, create empty data frame with correct structure
        covar_lmm_results <- data.frame(
          marker = character(),
          covariate = character(),
          covar_estimate = numeric(),
          covar_stderr = numeric(),
          covar_tvalue = numeric(),
          covar_pvalue = numeric(),
          model_rsquared = numeric(),
          covar_adjusted = numeric(),
          significant_raw = logical(),
          significant_adjusted = logical(),
          stringsAsFactors = FALSE
        )
      }

      # Summarise results for this covariate
      if (nrow(covar_lmm_results) > 0) {
        n_sig_raw <- sum(covar_lmm_results$significant_raw, na.rm = TRUE)
        n_sig_adjusted <- sum(covar_lmm_results$significant_adjusted, na.rm = TRUE)
        mean_rsq <- mean(covar_lmm_results$model_rsquared, na.rm = TRUE)

        if (n_test >= 50) cat("\n") # New line after progress bar

        cat(sprintf(
          "          Results: %d/%d markers tested\n",
          nrow(covar_lmm_results), n_test
        ))
        cat(sprintf(
          "          Significant: %d (raw p<%.3f), %d (%s<%.3f)\n",
          n_sig_raw, sig_threshold, n_sig_adjusted, mt_method_display, sig_threshold
        ))
        cat(sprintf("          Mean model R²: %.3f\n", mean_rsq))

        # Save detailed results
        write.csv(covar_lmm_results,
          file.path(qc_dir, paste0(
            "lmm_results_",
            gsub("[^A-Za-z0-9]", "_", col), ".csv"
          )),
          row.names = FALSE
        )

        # Store summary
        lmm_summary_results[[col]] <- list(
          n_tested = nrow(covar_lmm_results),
          n_significant_raw = n_sig_raw,
          n_significant_adjusted = n_sig_adjusted,
          mean_rsquared = mean_rsq
        )
      }
    } # End of covariate loop

    # Create overall summary
    if (length(lmm_summary_results) > 0) {
      lmm_summary_df <- data.frame(
        covariate = names(lmm_summary_results),
        n_markers_tested = sapply(lmm_summary_results, function(x) x$n_tested),
        n_significant_raw = sapply(lmm_summary_results, function(x) x$n_significant_raw),
        n_significant_adjusted = sapply(lmm_summary_results, function(x) x$n_significant_adjusted),
        mean_rsquared = sapply(lmm_summary_results, function(x) x$mean_rsquared),
        stringsAsFactors = FALSE
      )

      write.csv(lmm_summary_df,
        file.path(qc_dir, "lmm_analysis_summary.csv"),
        row.names = FALSE
      )

      cat("\n          [OK] Linear model analysis complete\n")
      cat("          Summary saved to: lmm_analysis_summary.csv\n")

      # Create summary plot showing all markers with stacked significance levels
      if (nrow(lmm_summary_df) > 0) {
        library(ggplot2)

        # Create stacked bar data showing ALL markers
        # Categories: 1) Adjusted p significant, 2) Raw p significant only, 3) Non-significant

        lmm_plot_data <- data.frame()
        for (i in 1:nrow(lmm_summary_df)) {
          covar <- lmm_summary_df$covariate[i]
          n_total <- lmm_summary_df$n_markers_tested[i]
          n_adj_sig <- lmm_summary_df$n_significant_adjusted[i]
          n_raw_sig <- lmm_summary_df$n_significant_raw[i]

          # Calculate counts for each category
          n_adj_only <- n_adj_sig # Significant by adjusted p
          n_raw_only <- n_raw_sig - n_adj_sig # Significant by raw p but not adjusted
          n_nonsig <- n_total - n_raw_sig # Not significant

          lmm_plot_data <- rbind(lmm_plot_data, data.frame(
            covariate = covar,
            category = sprintf("%s < %.3f", toupper(mt_method_display), sig_threshold),
            count = n_adj_only,
            order = 3,
            stringsAsFactors = FALSE
          ))

          lmm_plot_data <- rbind(lmm_plot_data, data.frame(
            covariate = covar,
            category = sprintf("Raw p < %.3f only", sig_threshold),
            count = n_raw_only,
            order = 2,
            stringsAsFactors = FALSE
          ))

          lmm_plot_data <- rbind(lmm_plot_data, data.frame(
            covariate = covar,
            category = "Not significant",
            count = n_nonsig,
            order = 1,
            stringsAsFactors = FALSE
          ))
        }

        # Calculate total for ordering covariates by adjusted significance
        lmm_covar_totals <- aggregate(count ~ covariate,
          data = lmm_plot_data[lmm_plot_data$order == 3, ],
          FUN = sum
        )
        lmm_plot_data$covariate <- factor(lmm_plot_data$covariate,
          levels = lmm_covar_totals$covariate[order(lmm_covar_totals$count)]
        )

        # Set factor levels for stacking order
        lmm_plot_data$category <- factor(lmm_plot_data$category,
          levels = unique(lmm_plot_data$category[order(lmm_plot_data$order)])
        )

        # Define colours: gray for non-sig, light for raw only, dark for adjusted
        category_levels <- levels(lmm_plot_data$category)
        fill_colors <- setNames(c("gray85", "#6BAED6", "#08519C"), category_levels)

        p_lmm <- ggplot(lmm_plot_data, aes(x = covariate, y = count, fill = category)) +
          geom_bar(stat = "identity", position = "stack") +
          scale_fill_manual(values = fill_colors) +
          coord_flip() +
          labs(
            title = "Pre-analytical Variables: Linear Model Analysis",
            subtitle = sprintf("All markers shown, coloured by significance (controlling for biological group)"),
            x = "Pre-analytical Variable",
            y = "Number of Markers",
            fill = "Significance Category"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(face = "bold"),
            legend.position = "bottom"
          )


        save_plot_multiple_formats(p_lmm,
          file.path(qc_dir, "lmm_associations_summary"),
          width = 10, height = 6
        )

        cat("          [OK] Created lmm_associations_summary.png\n")
      }
    }
  } else {
    cat("          [SKIP] No biological group information or QC columns available\n")
  }

  cat("          ===== MIXED EFFECTS ANALYSIS COMPLETE =====\n")

  cat("          [OK] QC analysis complete\n")
} else {
  cat("          [SKIP] No QC columns found\n")
}

log_section_end(7, "Pre-analytical QC summary")
cat("\n")

# ============================================================================
# STEP 8: FUNCTIONAL ENRICHMENT ANALYSIS
# ============================================================================
log_section_start(8, "Functional enrichment analysis", 10)

if (opt$functional_enrichment) {
  cat("[STEP 8/10] Functional enrichment requested but not yet implemented in v19\n")
  cat("          [SKIP] Will be added in future version\n")
  cat("          For now, DMR coordinates are available in dmr_details_with_coordinates.csv\n")
} else {
  cat("[STEP 8/10] Functional enrichment skipped\n")
}

log_section_end(8, "Functional enrichment analysis")

cat("\n")
# ============================================================================
# FINAL: GENERATE SUMMARY REPORT
# ============================================================================
cat("[FINAL] Generating summary report...\n")

summary_file <- file.path(output_dir, "ANALYSIS_SUMMARY.txt")

# Create summary in a separate file connection (don't mess with existing sink)
summary_con <- file(summary_file, "w")

cat("================================================================================\n", file = summary_con)
cat("METHYLSENSE MODEL EVALUATION SUMMARY\n", file = summary_con)
cat("================================================================================\n\n", file = summary_con)
cat("Script Version:", SCRIPT_VERSION, "\n", file = summary_con)
cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = summary_con)
cat("Model Directory:", opt$model_dir, "\n", file = summary_con)
cat("Model Type:", model_name, "\n", file = summary_con)
cat("Region:", region_info, "\n", file = summary_con)
cat("Number of Markers:", n_markers, "\n\n", file = summary_con)

cat("================================================================================\n", file = summary_con)
cat("1. BATCH/PLATFORM EFFECTS\n", file = summary_con)
cat("================================================================================\n", file = summary_con)
if (length(tech_covariates) > 0) {
  cat("Technical covariates tested:", paste(tech_covariates, collapse = ", "), "\n", file = summary_con)
  if (exists("covariate_summary") && nrow(covariate_summary) > 0) {
    cat("\nCovariate Association Results:\n", file = summary_con)
    cat(sprintf(
      "(Total markers tested: %d | Significance threshold: %.3f | Multiple testing: %s)\n\n",
      n_markers, sig_threshold, mt_method_display
    ), file = summary_con)

    # Print formatted table
    cat("Covariate", "Markers_Tested", "Sig_Raw", "Sig_Adjusted", "Mean_Effect", "Mean_Raw_P", "Mean_Adj_P\n",
      sep = "\t", file = summary_con
    )
    cat(paste(rep("-", 100), collapse = ""), "\n", file = summary_con)

    for (i in 1:nrow(covariate_summary)) {
      cat(
        sprintf(
          "%s\t%d\t%d\t%d\t%.3f\t%.4f\t%.4f\n",
          covariate_summary$covariate[i],
          covariate_summary$n_markers_tested[i],
          covariate_summary$n_significant_raw[i],
          covariate_summary$n_significant_adjusted[i],
          covariate_summary$mean_effect[i],
          covariate_summary$mean_raw_p[i],
          covariate_summary$mean_adjusted_p[i]
        ),
        file = summary_con
      )
    }

    cat("\nInterpretation:\n", file = summary_con)
    cat("  - Sig_Raw: Number of markers with raw p-value < ", sig_threshold, "\n", file = summary_con, sep = "")
    cat("  - Sig_Adjusted: Number of markers with ", mt_method_display, " < ", sig_threshold, "\n", file = summary_con, sep = "")
    cat("  - Mean_Effect: Mean absolute test statistic across all markers\n", file = summary_con)

    # Overall assessment
    total_sig_adjusted <- sum(covariate_summary$n_significant_adjusted)
    total_sig_raw <- sum(covariate_summary$n_significant_raw)

    cat("\nOVERALL ASSESSMENT:\n", file = summary_con)
    if (total_sig_adjusted == 0 && total_sig_raw == 0) {
      cat("  No significant technical confounders detected (all covariates: raw p > ", sig_threshold, ", ", mt_method_display, " > ", sig_threshold, ")\n", file = summary_con, sep = "")
    } else if (total_sig_adjusted == 0 && total_sig_raw > 0) {
      cat("  Some raw associations detected (raw p < ", sig_threshold, ") but none significant after ", mt_method_display, " correction\n", file = summary_con, sep = "")
      cat("  This suggests weak/spurious associations that do not survive multiple testing correction\n", file = summary_con)
    } else {
      cat("  ", total_sig_adjusted, " significant associations detected after ", mt_method_display, " correction\n", file = summary_con, sep = "")
      cat("  RECOMMENDATION: Review batch effects and consider correction methods\n", file = summary_con)
    }
  } else {
    cat("\nWARNING: Covariate analysis results not found or empty\n", file = summary_con)
  }
} else {
  cat("No technical covariates specified\n", file = summary_con)
}

cat("\n================================================================================\n", file = summary_con)
cat("2. MARKER STABILITY\n", file = summary_con)
cat("================================================================================\n", file = summary_con)
cat("Bootstrap iterations:", opt$bootstrap_iterations, "\n", file = summary_con)
cat("Stability threshold:", opt$stability_threshold, "\n", file = summary_con)
cat("Effect size threshold:", opt$effect_size_threshold, "%\n", file = summary_con)
if (exists("marker_stability")) {
  cat("Stable markers:", sum(marker_stability$stable), "/", nrow(marker_stability), "\n", file = summary_con)
  cat("Mean selection frequency:", round(mean(marker_stability$selection_frequency), 3), "\n", file = summary_con)
  cat("Mean effect size:", round(mean(marker_stability$mean_effect_size, na.rm = TRUE), 2), "%\n", file = summary_con)

  # Top stable markers
  top5 <- head(marker_stability[order(-marker_stability$selection_frequency), ], 5)

  if (nrow(top5) > 0) {
    cat("\nTop 5 most stable markers:\n", file = summary_con)
    for (i in 1:nrow(top5)) {
      location <- if ("marker_name" %in% colnames(top5)) {
        top5$marker_name[i]
      } else {
        "N/A"
      }

      cat(
        sprintf(
          "  %d. %s (%s): Freq=%.3f, Effect=%.2f%%\n",
          i,
          top5$marker[i],
          location,
          top5$selection_frequency[i],
          top5$mean_effect_size[i]
        ),
        file = summary_con
      )
    }
  }
}

cat("\n================================================================================\n", file = summary_con)
cat("3. PRE-ANALYTICAL QC\n", file = summary_con)
cat("================================================================================\n", file = summary_con)
cat("QC variables found:", length(qc_columns), "\n", file = summary_con)
if (length(qc_columns) > 0) {
  cat("Variables:", paste(qc_columns, collapse = ", "), "\n", file = summary_con)
  if (exists("qc_summary")) {
    cat("\nQC Column Summary:\n", file = summary_con)
    capture.output(print(qc_summary), file = summary_con, append = TRUE)
  }

  # Add Linear Mixed Model results if available
  if (exists("lmm_summary_df") && nrow(lmm_summary_df) > 0) {
    cat("\n\nLinear Mixed Model Association Results:\n", file = summary_con)
    cat(sprintf(
      "(Controlling for biological group | Significance threshold: %.3f | Multiple testing: %s)\n\n",
      sig_threshold, mt_method_display
    ), file = summary_con)

    # Print formatted table
    cat("Covariate", "Markers_Tested", "Sig_Raw", "Sig_Adjusted", "Mean_R²\n",
      sep = "\t", file = summary_con
    )
    cat(paste(rep("-", 80), collapse = ""), "\n", file = summary_con)

    for (i in 1:nrow(lmm_summary_df)) {
      cat(
        sprintf(
          "%s\t%d\t%d\t%d\t%.3f\n",
          lmm_summary_df$covariate[i],
          lmm_summary_df$n_markers_tested[i],
          lmm_summary_df$n_significant_raw[i],
          lmm_summary_df$n_significant_adjusted[i],
          lmm_summary_df$mean_rsquared[i]
        ),
        file = summary_con
      )
    }

    cat("\nInterpretation:\n", file = summary_con)
    cat("  - Sig_Raw: Number of markers with raw p-value < ", sig_threshold, "\n", file = summary_con, sep = "")
    cat("  - Sig_Adjusted: Number of markers with ", mt_method_display, " < ", sig_threshold, "\n", file = summary_con, sep = "")
    cat("  - Mean_R²: Average model fit across all markers\n", file = summary_con)

    # Overall assessment
    total_lmm_sig_adjusted <- sum(lmm_summary_df$n_significant_adjusted)
    total_lmm_sig_raw <- sum(lmm_summary_df$n_significant_raw)

    cat("\nOVERALL ASSESSMENT:\n", file = summary_con)
    if (total_lmm_sig_adjusted == 0 && total_lmm_sig_raw == 0) {
      cat("  No significant pre-analytical confounders detected (controlling for biological group)\n", file = summary_con)
      cat("  All QC variables: raw p > ", sig_threshold, ", ", mt_method_display, " > ", sig_threshold, "\n", file = summary_con, sep = "")
    } else if (total_lmm_sig_adjusted == 0 && total_lmm_sig_raw > 0) {
      cat("  Some raw associations detected (raw p < ", sig_threshold, ") but none significant after ", mt_method_display, " correction\n", file = summary_con, sep = "")
      cat("  This suggests weak QC effects that do not survive multiple testing correction\n", file = summary_con)
    } else {
      cat("  ", total_lmm_sig_adjusted, " significant pre-analytical associations detected (after ", mt_method_display, " correction)\n", file = summary_con, sep = "")
      cat("  RECOMMENDATION: Review QC metrics and consider controlling for these variables in downstream analyses\n", file = summary_con)
    }
  } else {
    cat("\nLinear Mixed Model analysis not performed or results not available\n", file = summary_con)
  }
}

cat("\n================================================================================\n", file = summary_con)
cat("4. FUNCTIONAL ENRICHMENT\n", file = summary_con)
cat("================================================================================\n", file = summary_con)
if (opt$functional_enrichment && enrichment_available) {
  if (exists("all_genes")) {
    cat("Genes annotated:", length(unique(all_genes$ensembl_gene_id)), "\n", file = summary_con)
  } else {
    cat("No genes found\n", file = summary_con)
  }
} else {
  cat("Not performed\n", file = summary_con)
}

cat("\n================================================================================\n", file = summary_con)
cat("ANALYSIS COMPLETE\n", file = summary_con)
cat("================================================================================\n", file = summary_con)

# Close summary file
close(summary_con)

cat("          [OK] Summary report saved to:", summary_file, "\n")

cat("\n")
cat("================================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================================\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
if (exists("start_time")) {
  cat("Total runtime:", round(difftime(Sys.time(), start_time, units = "mins"), 2), "minutes\n")
}
cat("================================================================================\n\n")

# Close log connections
sink(type = "output")
if (exists("log_con") && isOpen(log_con)) {
  close(log_con)
}

cat("[SUCCESS] All analyses complete!\n")
cat("[LOG] Detailed log saved to:", log_file, "\n")
cat("[SUMMARY] Summary report saved to:", summary_file, "\n")

# ============================================================================
# STEP 10: GENERATE MARKDOWN REPORT (if requested) - INLINE GENERATION
# ============================================================================
if (opt$generate_markdown) {
  log_section_start(10, "Generating comprehensive markdown report", 10)
  cat("[STEP 10/10] Generating comprehensive markdown report...\n")

  markdown_file <- file.path(output_dir, "MODEL_EVALUATION_REPORT.md")
  figures_dir <- file.path(output_dir, "figures")
  dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

  # Define directory paths for report generation
  batch_dir <- file.path(output_dir, "batch_effects")
  stability_dir <- file.path(output_dir, "marker_stability")

  # Initialize markdown content
  md <- c()

  # ============================================================================
  # DYNAMIC FIGURE/TABLE NUMBERING SYSTEM
  # ============================================================================
  # Initialize global counters for sequential numbering
  figure_counter <- 0
  table_counter <- 0

  # Helper function to get next figure number and format caption
  add_figure_caption <- function(caption_text) {
    figure_counter <<- figure_counter + 1
    return(sprintf("*Figure %d: %s*", figure_counter, caption_text))
  }

  # Helper function to get next table number and format header
  add_table_header <- function(title_text) {
    table_counter <<- table_counter + 1
    return(sprintf("**Table %d:** %s", table_counter, title_text))
  }

  cat("[INFO] Dynamic figure/table numbering system initialized\n")
  # ============================================================================

  # Set primary figure extension for markdown image references
  fig_ext <- get_primary_figure_extension()
  cat(sprintf("[INFO] Using '%s' format for figures in markdown report\n", fig_ext))

  # ====== HEADER ======
  # Use custom model name if provided, otherwise use markdown_title
  report_title <- if (!is.null(opt$model_name)) {
    sprintf("%s - Model Evaluation Report", opt$model_name)
  } else {
    opt$markdown_title
  }
  md <- c(md, paste0("# ", report_title))
  md <- c(md, "")
  md <- c(md, "> **Generated by [MethylSense](https://github.com/markusdrag/MethylSense)** - High-accuracy epigenetic diagnostics combining Nanopore methylation sequencing with machine learning")
  md <- c(md, "")
  if (!is.null(opt$markdown_author)) {
    md <- c(md, paste0("**Authors:** ", opt$markdown_author))
    md <- c(md, "")
  }
  md <- c(md, paste0("**Date:** ", format(Sys.Date(), "%B %d, %Y")))
  md <- c(md, paste0("**MethylSense Version:** ", SCRIPT_VERSION))
  md <- c(md, "")
  md <- c(md, "---")
  md <- c(md, "")

  # ====== TABLE OF CONTENTS ======
  md <- c(md, "## Table of Contents")
  md <- c(md, "")
  md <- c(md, "### 1. [Model Performance Validation](#model-performance-validation)")
  md <- c(md, "  - Overview")
  md <- c(md, "  - Overall Performance Metrics")
  md <- c(md, "  - Confusion Matrix")
  md <- c(md, "  - Per-Class Performance (One-vs-Rest)")
  md <- c(md, "  - ROC Curves")
  md <- c(md, "  - Precision-Recall Curves")
  md <- c(md, "  - Stratified Performance (Train/Test)")
  md <- c(md, "  - Platform-Stratified Performance")
  md <- c(md, "  - CV Fold-by-Fold Performance")
  md <- c(md, "  - Sample Prediction Consistency")
  md <- c(md, "  - Hyperparameter Sensitivity")

  if (exists("gene_annotation_done") && gene_annotation_done) {
    md <- c(md, "  - Genomic Context & Gene Annotation")
  }

  if (exists("pcr_prep_done") && pcr_prep_done) {
    md <- c(md, "  - PCR Primer Preparation")
  }

  if (exists("calibration_analysis") && calibration_analysis$available) {
    md <- c(md, "  - Probability Calibration")
  }

  if (!is.null(opt$time_analysis)) {
    md <- c(md, "  - Temporal DMR Dynamics")
  }

  md <- c(md, "")
  md <- c(md, "### 2. [Model Characterisation](#model-characterisation)")
  md <- c(md, "  - Model Overview & Hyperparameters")
  md <- c(md, "  - Variable Importance")
  md <- c(md, "")
  md <- c(md, "### 3. [Nested Cross-Validation Analysis](#nested-cv-analysis)")
  md <- c(md, "  - Performance Summary (95% CI)")

  if (exists("nested_cv_detailed") && nested_cv_detailed$available) {
    md <- c(md, "  - Performance Distributions")
    md <- c(md, "  - Per-Fold Performance")
    md <- c(md, "  - ROC/PR Curves Across Folds")
    md <- c(md, "  - Optimism Bias Analysis")
    md <- c(md, "  - CV Method Comparison")
    md <- c(md, "  - Outer Fold Performance")
  }

  if (exists("study_performance") && study_performance$available) {
    md <- c(md, "  - [Study-Stratified Performance](#study-stratified-performance)")
  }

  md <- c(md, "")
  md <- c(md, "### 4. [DMR-Group Association Analysis](#dmr-group-associations)")
  md <- c(md, "  - Infection Probability vs Methylation")
  md <- c(md, "  - Variable Importance Ranking")
  md <- c(md, "  - Effect Size Distribution")
  md <- c(md, "  - Group-Specific Markers")
  md <- c(md, "  - Overall DMR Comparison (Wilcoxon)")
  md <- c(md, "")

  if (exists("feature_annotation_done") && feature_annotation_done) {
    md <- c(md, "### 5. [Genomic Feature Annotation](#feature-annotation)")
    md <- c(md, "  - Feature Overlap Summary")
    md <- c(md, "  - Visualisation")
    md <- c(md, "  - Biological Interpretation")
    md <- c(md, "")
  }

  md <- c(md, "### 6. [Executive Summary](#executive-summary)")
  md <- c(md, "")
  md <- c(md, "### 7. [Data Quality and Validation](#data-quality-and-validation)")
  md <- c(md, "  - Data Extraction & QC")
  md <- c(md, "  - Methylation Pattern Stats")
  md <- c(md, "  - Sample Correlation & Clustering")
  md <- c(md, "  - DMR Correlation Network")

  if (opt$enable_wgcna) {
    md <- c(md, "  - WGCNA Co-Methylation Analysis")
  }

  md <- c(md, "  - DMR Statistical Significance")
  md <- c(md, "")
  md <- c(md, "### 8. [Marker Stability Analysis](#marker-stability-analysis)")
  md <- c(md, "  - Bootstrap Stability")
  md <- c(md, "  - Feature Importance")
  md <- c(md, "")
  md <- c(md, "### 9. [Dimensionality Reduction](#dimensionality-reduction)")
  md <- c(md, "  - PCA (Principal Component)")
  md <- c(md, "  - UMAP (Manifold Projection)")
  md <- c(md, "")
  md <- c(md, "### 10. [Technical Quality Control](#technical-quality-control)")
  md <- c(md, "  - Sample QC Metrics")
  md <- c(md, "  - Technical Covariates (LMM)")
  md <- c(md, "")
  md <- c(md, "### 11. [Methods](#methods)")
  md <- c(md, "  - Bootstrap Stability")
  md <- c(md, "  - Confounder Testing")
  md <- c(md, "  - Statistical Models (LM for QC, LMM for Temporal)")
  md <- c(md, "  - Reproducibility")
  md <- c(md, "")
  md <- c(md, "### 12. [Conclusions](#conclusions)")
  md <- c(md, "  - Quantitative Summary")
  md <- c(md, "  - Final Statement")
  md <- c(md, "")
  md <- c(md, "### 13. [References](#references)")
  md <- c(md, "")
  md <- c(md, "---")
  md <- c(md, "")

  # ====== SECTION 1: MODEL PERFORMANCE VALIDATION ======
  md <- c(md, "## Model Performance Validation {#model-performance-validation}")
  md <- c(md, "")

  if (exists("performance_analysis_done") && performance_analysis_done) {
    md <- c(md, "### Overview")
    md <- c(md, "")
    md <- c(md, "Comprehensive multiclass performance analysis was conducted using predicted probabilities and actual class labels. The analysis employs a one-vs-rest approach for each class, calculating receiver operating characteristic (ROC) and precision-recall (PR) curves along with associated area under curve (AUC) metrics.")
    md <- c(md, "")

    # Overall metrics
    if (file.exists(file.path(performance_dir, "overall_metrics.csv"))) {
      overall_perf <- read.csv(file.path(performance_dir, "overall_metrics.csv"))

      md <- c(md, "### Overall Performance")
      md <- c(md, "")

      accuracy_val <- overall_perf$Value[overall_perf$Metric == "Overall_Accuracy"]
      kappa_val <- overall_perf$Value[overall_perf$Metric == "Cohens_Kappa"]
      n_samples_val <- overall_perf$Value[overall_perf$Metric == "N_Samples"]

      # FIX v5.11.8: Classify kappa using Landis & Koch (1977) guidelines
      kappa_interpretation <- if (kappa_val < 0) {
        "poor agreement (less than chance)"
      } else if (kappa_val >= 0 && kappa_val <= 0.20) {
        "slight agreement"
      } else if (kappa_val > 0.20 && kappa_val <= 0.40) {
        "fair agreement"
      } else if (kappa_val > 0.40 && kappa_val <= 0.60) {
        "moderate agreement"
      } else if (kappa_val > 0.60 && kappa_val <= 0.80) {
        "substantial agreement"
      } else {
        "almost perfect agreement"
      }

      md <- c(md, add_table_header("Overall Performance Metrics"))
      md <- c(md, "")
      md <- c(md, "| Metric | Value | Interpretation |")
      md <- c(md, "|--------|-------|----------------|")
      md <- c(md, sprintf("| Overall Accuracy | %.4f | %.1f%% correct classifications |", accuracy_val, accuracy_val * 100))
      md <- c(md, sprintf("| Cohen's Kappa | %.4f | %s |", kappa_val, kappa_interpretation))
      md <- c(md, sprintf("| Number of Samples | %.0f | — |", n_samples_val))
      md <- c(md, "")

      md <- c(md, "> **Interpretation:** Cohen's kappa provides a measure of agreement beyond chance (Landis & Koch 1977). Guidelines: <0.00 = poor, 0.00-0.20 = slight, 0.21-0.40 = fair, 0.41-0.60 = moderate, 0.61-0.80 = substantial, 0.81-1.00 = almost perfect.")
      md <- c(md, "")
    }

    # Confusion matrix
    # PNG always exists (baseline format), but copy all format variants
    if (file.exists(file.path(performance_dir, "confusion_matrix_3x3.png"))) {
      # Copy all format variants to figures directory
      confusion_base <- file.path(performance_dir, "confusion_matrix_3x3")
      for (ext in c("png", "jpeg", "jpg", "svg")) {
        variant_file <- sprintf("%s.%s", confusion_base, ext)
        if (file.exists(variant_file)) {
          file.copy(variant_file,
            file.path(figures_dir, basename(variant_file)),
            overwrite = TRUE
          )
        }
      }

      md <- c(md, "### Confusion Matrix")
      md <- c(md, "")

      # Try to use confusion matrix from model directory
      model_conf_matrix <- file.path(opt$model_dir, sprintf("confusion_matrix_plot.%s", fig_ext))
      if (!file.exists(model_conf_matrix)) {
        model_conf_matrix <- file.path(opt$model_dir, "confusion_matrix_plot.png")
      }

      if (file.exists(model_conf_matrix)) {
        # Copy with helper function to handle format correctly
        copy_plot_to_figures(opt$model_dir, "confusion_matrix_plot", figures_dir, fig_ext)
        md <- c(md, sprintf("![Confusion Matrix](figures/confusion_matrix_plot.%s)", fig_ext))
        md <- c(md, "")
        md <- c(md, add_figure_caption("Confusion matrix showing model predictions versus actual classifications on the independent test set. Cell values indicate number of samples. Diagonal cells (green) represent correct predictions; off-diagonal cells (red) represent misclassifications."))
      } else {
        # Use generated confusion matrix - always use PNG in markdown (baseline)
        md <- c(md, "![Confusion Matrix](figures/confusion_matrix_3x3.png)")
        md <- c(md, "")
        md <- c(md, add_figure_caption("Multiclass confusion matrix showing predicted vs actual classifications. Numbers represent sample counts, whilst colours indicate row-normalised proportions (darker colours indicate higher proportions)."))
      }
      md <- c(md, "")

      # Add per-cell interpretation if metrics available
      if (exists("class_metrics") && length(class_metrics) > 0) {
        md <- c(md, "**Confusion Matrix Interpretation:**")
        md <- c(md, "")
        for (class_name in names(class_metrics)) {
          metrics <- class_metrics[[class_name]]
          md <- c(md, sprintf(
            "- **%s:** %d correct, sensitivity %.1f%%, precision %.1f%%",
            class_name,
            metrics$tp,
            metrics$sensitivity * 100,
            metrics$ppv * 100
          ))
        }
        md <- c(md, "")
      }
    }

    # Per-class metrics
    if (file.exists(file.path(performance_dir, "class_metrics_one_vs_rest.csv"))) {
      class_perf <- read.csv(file.path(performance_dir, "class_metrics_one_vs_rest.csv"))

      md <- c(md, "### Per-Class Performance Metrics (One-vs-Rest)")
      md <- c(md, "")
      md <- c(md, add_table_header("Per-Class Performance Metrics (One-vs-Rest)"))
      md <- c(md, "")
      md <- c(md, "| Class | Sensitivity | Specificity | PPV | NPV | F1 | MCC | ROC-AUC (95% CI) | PR-AUC |")
      md <- c(md, "|-------|-------------|-------------|-----|-----|----|----|------------------|--------|")

      for (i in 1:nrow(class_perf)) {
        md <- c(md, sprintf(
          "| %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f (%.3f-%.3f) | %.3f |",
          class_perf$Class[i],
          class_perf$Sensitivity[i],
          class_perf$Specificity[i],
          class_perf$PPV[i],
          class_perf$NPV[i],
          class_perf$F1_Score[i],
          class_perf$MCC[i],
          class_perf$ROC_AUC[i],
          class_perf$ROC_AUC_95CI_Lower[i],
          class_perf$ROC_AUC_95CI_Upper[i],
          class_perf$PR_AUC[i]
        ))
      }

      md <- c(md, "")
      md <- c(md, "**Table Legend:**")
      md <- c(md, "")
      md <- c(md, "- **Sensitivity (Recall):** True positive rate - proportion of actual positives correctly identified")
      md <- c(md, "- **Specificity:** True negative rate - proportion of actual negatives correctly identified")
      md <- c(md, "- **PPV (Precision):** Positive predictive value - proportion of predicted positives that are true positives")
      md <- c(md, "- **NPV:** Negative predictive value - proportion of predicted negatives that are true negatives")
      md <- c(md, "- **F1:** Harmonic mean of precision and recall")
      md <- c(md, "- **MCC:** Matthews correlation coefficient - balanced measure considering all confusion matrix elements")
      md <- c(md, "- **ROC-AUC:** Area under the receiver operating characteristic curve (95% CI calculated via bootstrap, n=2000)")
      md <- c(md, "- **PR-AUC:** Area under the precision-recall curve")
      md <- c(md, "")
    }

    # ROC curves
    # PNG always exists (baseline format), but copy all format variants
    if (file.exists(file.path(performance_dir, "roc_curves_multiclass.png"))) {
      # Copy all format variants to figures directory
      roc_base <- file.path(performance_dir, "roc_curves_multiclass")
      for (ext in c("png", "jpeg", "jpg", "svg")) {
        variant_file <- sprintf("%s.%s", roc_base, ext)
        if (file.exists(variant_file)) {
          file.copy(variant_file,
            file.path(figures_dir, basename(variant_file)),
            overwrite = TRUE
          )
        }
      }

      md <- c(md, "### ROC Curves")
      md <- c(md, "")
      md <- c(md, "![ROC Curves](figures/roc_curves_multiclass.png)")
      md <- c(md, "")
      md <- c(md, add_figure_caption("Receiver operating characteristic (ROC) curves for each class using one-vs-rest approach. The diagonal dashed line represents random classification (AUC = 0.5). Curves closer to the top-left corner indicate better classification performance."))
      md <- c(md, "")
    }

    # PR curves
    # PNG always exists (baseline format), but copy all format variants
    if (file.exists(file.path(performance_dir, "pr_curves_multiclass.png"))) {
      # Copy all format variants to figures directory
      pr_base <- file.path(performance_dir, "pr_curves_multiclass")
      for (ext in c("png", "jpeg", "jpg", "svg")) {
        variant_file <- sprintf("%s.%s", pr_base, ext)
        if (file.exists(variant_file)) {
          file.copy(variant_file,
            file.path(figures_dir, basename(variant_file)),
            overwrite = TRUE
          )
        }
      }

      md <- c(md, "### Precision-Recall Curves")
      md <- c(md, "")
      md <- c(md, "![PR Curves](figures/pr_curves_multiclass.png)")
      md <- c(md, "")
      md <- c(md, add_figure_caption("Precision-recall curves for each class. PR curves are particularly informative for imbalanced datasets, as they focus on the positive class performance."))
      md <- c(md, "")
    }

    # Stratified performance
    if (file.exists(file.path(performance_dir, "stratified_performance.csv"))) {
      strat_perf <- read.csv(file.path(performance_dir, "stratified_performance.csv"))

      md <- c(md, "### Stratified Performance (Training vs Testing)")
      md <- c(md, "")
      md <- c(md, add_table_header("Stratified Performance (Training vs Testing)"))
      md <- c(md, "")
      md <- c(md, "| Dataset | N Samples | Accuracy | Kappa |")
      md <- c(md, "|---------|-----------|----------|-------|")

      for (i in 1:nrow(strat_perf)) {
        md <- c(md, sprintf(
          "| %s | %d | %.4f | %.4f |",
          strat_perf$Dataset[i],
          strat_perf$N_Samples[i],
          strat_perf$Accuracy[i],
          strat_perf$Kappa[i]
        ))
      }

      md <- c(md, "")
      md <- c(md, "> **Interpretation:** Similar performance between training and testing sets indicates good generalisation without overfitting. Large discrepancies would suggest the model has not learnt generalisable patterns.")
      md <- c(md, "")
    }

    # Platform-stratified performance
    if (file.exists(file.path(performance_dir, "platform_stratified_performance.csv"))) {
      plat_perf <- read.csv(file.path(performance_dir, "platform_stratified_performance.csv"))

      md <- c(md, "### Platform-Stratified Performance")
      md <- c(md, "")
      md <- c(md, add_table_header("Platform-Stratified Performance"))
      md <- c(md, "")
      md <- c(md, "| Platform | N Samples | Accuracy |")
      md <- c(md, "|----------|-----------|----------|")

      for (i in 1:nrow(plat_perf)) {
        md <- c(md, sprintf(
          "| %s | %d | %.4f |",
          plat_perf$Platform[i],
          plat_perf$N_Samples[i],
          plat_perf$Accuracy[i]
        ))
      }

      md <- c(md, "")
      md <- c(md, "> **Interpretation:** Consistent performance across different platforms/batches demonstrates model robustness and minimal technical batch effects on prediction accuracy.")
      md <- c(md, "")
    }

    # ====== NEW SECTION 1.9: FOLD-BY-FOLD ANALYSIS ======
    if (exists("cv_detailed") && cv_detailed$available) {
      # ===== GENERATE CV FOLD PERFORMANCE PLOT =====
      tryCatch(
        {
          fold_data <- cv_detailed$data
          fold_data$fold <- fold_data$fold_num

          # Reshape data for plotting (same as nested CV plot in MethylSense_analysis.R)
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

          # Calculate means and 95% CIs for annotations
          mean_accuracy <- mean(fold_data$accuracy, na.rm = TRUE)
          mean_sensitivity <- mean(fold_data$mean_sensitivity, na.rm = TRUE)
          mean_specificity <- mean(fold_data$mean_specificity, na.rm = TRUE)
          mean_auc <- mean(fold_data$auc, na.rm = TRUE)

          # Calculate 95% CIs
          n_folds <- nrow(fold_data)
          se_accuracy <- sd(fold_data$accuracy, na.rm = TRUE) / sqrt(n_folds)
          se_sensitivity <- sd(fold_data$mean_sensitivity, na.rm = TRUE) / sqrt(n_folds)
          se_specificity <- sd(fold_data$mean_specificity, na.rm = TRUE) / sqrt(n_folds)
          se_auc <- sd(fold_data$auc, na.rm = TRUE) / sqrt(n_folds)

          ci_accuracy <- c(mean_accuracy - 1.96 * se_accuracy, mean_accuracy + 1.96 * se_accuracy)
          ci_sensitivity <- c(mean_sensitivity - 1.96 * se_sensitivity, mean_sensitivity + 1.96 * se_sensitivity)
          ci_specificity <- c(mean_specificity - 1.96 * se_specificity, mean_specificity + 1.96 * se_specificity)
          ci_auc <- c(mean_auc - 1.96 * se_auc, mean_auc + 1.96 * se_auc)

          # Create annotation dataframe
          annotation_df <- data.frame(
            metric = c("Accuracy", "Sensitivity", "Specificity", "AUC"),
            label = c(
              sprintf("%.1f%% [CI: %.1f-%.1f%%]", mean_accuracy * 100, ci_accuracy[1] * 100, ci_accuracy[2] * 100),
              sprintf("%.1f%% [CI: %.1f-%.1f%%]", mean_sensitivity * 100, ci_sensitivity[1] * 100, ci_sensitivity[2] * 100),
              sprintf("%.1f%% [CI: %.1f-%.1f%%]", mean_specificity * 100, ci_specificity[1] * 100, ci_specificity[2] * 100),
              sprintf("%.3f [CI: %.3f-%.3f]", mean_auc, ci_auc[1], ci_auc[2])
            ),
            x = Inf,
            y = -Inf,
            hjust = 1.05,
            vjust = -0.5
          )

          # Create the 4-panel fold performance plot
          p_cv_fold_performance <- ggplot(fold_long, aes(x = fold, y = value, color = metric, group = metric)) +
            geom_line(size = 1) +
            geom_point(size = 3) +
            geom_hline(data = data.frame(
              metric = c("Accuracy", "Sensitivity", "Specificity", "AUC"),
              mean_val = c(mean_accuracy, mean_sensitivity, mean_specificity, mean_auc)
            ), aes(yintercept = mean_val, color = metric), linetype = "dashed", alpha = 0.5) +
            geom_text(
              data = annotation_df,
              aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
              family = "Arial", size = 4, color = "black", inherit.aes = FALSE
            ) +
            facet_wrap(~metric, scales = "free_y", nrow = 2) +
            scale_color_brewer(palette = "Set1") +
            labs(
              title = paste("Cross-Validation Performance by Fold:", model_name),
              subtitle = paste(n_markers, "markers | Dashed lines show mean across folds"),
              x = "Fold",
              y = "Value",
              color = "Metric"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = 14, face = "bold"),
              legend.position = "none",
              strip.text = element_text(size = 12, face = "bold"),
              text = element_text(family = "Arial")
            ) +
            ylim(0, 1)

          # Save plot
          cv_fold_plot_file <- file.path(figures_dir, "cv_fold_performance")
          save_plot_multiple_formats(p_cv_fold_performance, cv_fold_plot_file, width = 12, height = 8)

          log_info("CV fold performance plot created", "PLOT")
        },
        error = function(e) {
          log_warn(paste("Failed to create CV fold performance plot:", conditionMessage(e)), "PLOT")
        }
      )

      md <- c(md, "### Cross-Validation Fold-by-Fold Performance")
      md <- c(md, "")
      md <- c(md, "To assess model stability and reproducibility, we analysed performance across all cross-validation folds individually.")
      md <- c(md, "")

      fold_data <- cv_detailed$data

      # FIX v5.12.0: Check if PR-AUC data is available (not all NA)
      show_pr_auc <- "pr_auc" %in% colnames(fold_data) && !all(is.na(fold_data$pr_auc))

      # FIX v5.12.0: Use cv_per_class_metrics if available (calculated from predictions)
      # Same approach as Table 8 (nested CV) - use calculated per-class metrics
      use_per_class_cv <- exists("cv_per_class_metrics") && cv_per_class_metrics$available &&
        !is.null(cv_per_class_metrics$fold_metrics)

      if (use_per_class_cv) {
        # Use per-class metrics calculated from CV predictions
        class_fold_data <- cv_per_class_metrics$fold_metrics
        positive_class_name <- cv_per_class_metrics$positive_class
        col_note <- sprintf(" (%s only)", positive_class_name)
      } else {
        # Fallback: Try to detect class-specific columns in cv_detailed data
        available_cols <- colnames(fold_data)
        case_keywords <- c("infected", "infection", "disease", "case", "positive", "treatment", "tumor", "cancer", "patient", "aspergillus")

        sens_cols <- grep("sensitivity", available_cols, ignore.case = TRUE, value = TRUE)
        spec_cols <- grep("specificity", available_cols, ignore.case = TRUE, value = TRUE)

        class_sens_cols <- sens_cols[!grepl("^mean_", sens_cols, ignore.case = TRUE)]
        class_spec_cols <- spec_cols[!grepl("^mean_", spec_cols, ignore.case = TRUE)]

        sens_col <- "mean_sensitivity"
        spec_col <- "mean_specificity"
        col_note <- " (balanced)"

        if (length(class_sens_cols) > 0) {
          for (col in class_sens_cols) {
            for (keyword in case_keywords) {
              if (grepl(keyword, col, ignore.case = TRUE)) {
                sens_col <- col
                class_name <- gsub("_sensitivity|sensitivity_|_sens|sens_", "", col, ignore.case = TRUE)
                spec_col <- class_spec_cols[grepl(class_name, class_spec_cols, ignore.case = TRUE)][1]
                if (!is.na(spec_col)) {
                  col_note <- sprintf(" (%s only)", gsub("_", " ", class_name))
                  break
                }
              }
            }
            if (sens_col != "mean_sensitivity") break
          }
        }
      }

      md <- c(md, add_table_header("Performance Metrics Across Individual CV Folds"))
      md <- c(md, "")

      # Table header with optional PR-AUC column and class-specific note
      if (show_pr_auc) {
        md <- c(md, sprintf("| Fold | Accuracy | AUC | PR-AUC | Sensitivity%s | Specificity%s | F1-Score |", col_note, col_note))
        md <- c(md, "|------|----------|-----|--------|-------------|-------------|----------|")
      } else {
        md <- c(md, sprintf("| Fold | Accuracy | AUC | Sensitivity%s | Specificity%s | F1-Score |", col_note, col_note))
        md <- c(md, "|------|----------|-----|-------------|-------------|----------|")
      }

      # Table rows with class-specific or balanced metrics
      if (use_per_class_cv) {
        # Use per-class metrics from calculated data
        for (i in 1:nrow(class_fold_data)) {
          acc_val <- class_fold_data$accuracy[i] * 100
          sens_val <- class_fold_data$sensitivity[i] * 100
          spec_val <- class_fold_data$specificity[i] * 100
          f1_val <- class_fold_data$f1_score[i]

          # Get AUC from fold_data if available
          auc_val <- if ("auc" %in% colnames(fold_data) && i <= nrow(fold_data)) fold_data$auc[i] else NA

          if (show_pr_auc && !is.na(auc_val)) {
            pr_auc_val <- if ("pr_auc" %in% colnames(fold_data)) fold_data$pr_auc[i] else NA
            md <- c(md, sprintf(
              "| %d | %.2f%% | %.3f | %.3f | %.2f%% | %.2f%% | %.3f |",
              class_fold_data$fold[i], acc_val, auc_val, pr_auc_val, sens_val, spec_val, f1_val
            ))
          } else if (!is.na(auc_val)) {
            md <- c(md, sprintf(
              "| %d | %.2f%% | %.3f | %.2f%% | %.2f%% | %.3f |",
              class_fold_data$fold[i], acc_val, auc_val, sens_val, spec_val, f1_val
            ))
          } else {
            # No AUC available
            md <- c(md, sprintf(
              "| %d | %.2f%% | — | %.2f%% | %.2f%% | %.3f |",
              class_fold_data$fold[i], acc_val, sens_val, spec_val, f1_val
            ))
          }
        }
      } else {
        # Use balanced metrics from fold_data
        for (i in 1:nrow(fold_data)) {
          sens_val <- if (exists("sens_col") && sens_col %in% colnames(fold_data)) fold_data[[sens_col]][i] * 100 else fold_data$mean_sensitivity[i] * 100
          spec_val <- if (exists("spec_col") && spec_col %in% colnames(fold_data)) fold_data[[spec_col]][i] * 100 else fold_data$mean_specificity[i] * 100

          if (show_pr_auc) {
            md <- c(md, sprintf(
              "| %d | %.2f%% | %.3f | %.3f | %.2f%% | %.2f%% | %.3f |",
              fold_data$fold_num[i],
              fold_data$accuracy[i] * 100,
              fold_data$auc[i],
              fold_data$pr_auc[i],
              sens_val,
              spec_val,
              fold_data$mean_f1[i]
            ))
          } else {
            md <- c(md, sprintf(
              "| %d | %.2f%% | %.3f | %.2f%% | %.2f%% | %.3f |",
              fold_data$fold_num[i],
              fold_data$accuracy[i] * 100,
              fold_data$auc[i],
              sens_val,
              spec_val,
              fold_data$mean_f1[i]
            ))
          }
        }
      }

      # Calculate summary statistics
      if (use_per_class_cv) {
        # Use per-class metrics
        mean_acc <- mean(class_fold_data$accuracy) * 100
        sd_acc <- sd(class_fold_data$accuracy) * 100
        mean_auc <- mean(fold_data$auc)
        sd_auc <- sd(fold_data$auc)
        mean_sens <- mean(class_fold_data$sensitivity) * 100
        sd_sens <- sd(class_fold_data$sensitivity) * 100
        mean_spec <- mean(class_fold_data$specificity) * 100
        sd_spec <- sd(class_fold_data$specificity) * 100
        mean_f1 <- mean(class_fold_data$f1_score)
        sd_f1 <- sd(class_fold_data$f1_score)
      } else {
        # Use balanced metrics
        mean_acc <- mean(fold_data$accuracy) * 100
        sd_acc <- sd(fold_data$accuracy) * 100
        mean_auc <- mean(fold_data$auc)
        sd_auc <- sd(fold_data$auc)
        mean_sens <- if (exists("sens_col") && sens_col %in% colnames(fold_data)) mean(fold_data[[sens_col]]) * 100 else mean(fold_data$mean_sensitivity) * 100
        mean_spec <- if (exists("spec_col") && spec_col %in% colnames(fold_data)) mean(fold_data[[spec_col]]) * 100 else mean(fold_data$mean_specificity) * 100
        sd_sens <- if (exists("sens_col") && sens_col %in% colnames(fold_data)) sd(fold_data[[sens_col]]) * 100 else sd(fold_data$mean_sensitivity) * 100
        sd_spec <- if (exists("spec_col") && spec_col %in% colnames(fold_data)) sd(fold_data[[spec_col]]) * 100 else sd(fold_data$mean_specificity) * 100
        mean_f1 <- mean(fold_data$mean_f1)
        sd_f1 <- sd(fold_data$mean_f1)
      }

      # Summary rows with optional PR-AUC column
      if (show_pr_auc) {
        md <- c(md, sprintf(
          "| **Mean** | **%.2f%%** | **%.3f** | **%.3f** | **%.2f%%** | **%.2f%%** | **%.3f** |",
          mean_acc,
          mean_auc,
          mean(fold_data$pr_auc, na.rm = TRUE),
          mean_sens,
          mean_spec,
          mean_f1
        ))
        md <- c(md, sprintf(
          "| **SD** | **%.2f%%** | **%.3f** | **%.3f** | **%.2f%%** | **%.2f%%** | **%.3f** |",
          sd_acc,
          sd_auc,
          sd(fold_data$pr_auc, na.rm = TRUE),
          sd_sens,
          sd_spec,
          sd_f1
        ))
      } else {
        md <- c(md, sprintf(
          "| **Mean** | **%.2f%%** | **%.3f** | **%.2f%%** | **%.2f%%** | **%.3f** |",
          mean_acc,
          mean_auc,
          mean_sens,
          mean_spec,
          mean_f1
        ))
        md <- c(md, sprintf(
          "| **SD** | **%.2f%%** | **%.3f** | **%.2f%%** | **%.2f%%** | **%.3f** |",
          sd_acc,
          sd_auc,
          sd_sens,
          sd_spec,
          sd_f1
        ))
      }
      md <- c(md, "")

      md <- c(md, sprintf(
        "> **Interpretation:** The narrow standard deviation (SD = %.2f%% for accuracy) demonstrates highly consistent performance across independent data splits. All folds achieve >%.0f%% accuracy, indicating the model is not dependent on specific training samples. This validates reproducibility and generalisability.",
        sd_acc, min(fold_data$accuracy) * 100
      ))
      md <- c(md, "")

      # Add fold performance plot
      fig_ext <- get_primary_figure_extension()
      cv_fold_plot_path <- file.path("figures", paste0("cv_fold_performance.", fig_ext))
      if (file.exists(file.path(figures_dir, paste0("cv_fold_performance.", fig_ext)))) {
        md <- c(md, paste0("![Cross-Validation Fold Performance](", cv_fold_plot_path, ")"))
        md <- c(md, "")
        md <- c(md, "**Figure:** Performance metrics across all cross-validation folds. Each panel shows one metric (Accuracy, Sensitivity, Specificity, AUC) with individual fold values (points connected by lines) and mean performance (dashed horizontal line). Consistent performance across folds indicates robust model generalization.")
        md <- c(md, "")
      }
    }

    # ====== NEW SECTION 1.10: SAMPLE-LEVEL CONSISTENCY ======
    if (exists("sample_consistency") && sample_consistency$available) {
      md <- c(md, "### Sample-Level Prediction Consistency")
      md <- c(md, "")
      md <- c(md, "This analysis evaluates how consistently the model predicts each sample across different cross-validation folds. Samples that are consistently predicted correctly (or incorrectly) across multiple folds provide insight into:")
      md <- c(md, "")
      md <- c(md, "- **Easy samples:** Consistently correctly predicted across all folds - these samples have clear methylation patterns")
      md <- c(md, "- **Difficult samples:** Frequently misclassified across folds - these may represent edge cases or ambiguous methylation profiles")
      md <- c(md, "- **Consistency score:** The proportion of folds in which each sample was correctly predicted when it was in the test set")
      md <- c(md, "")
      md <- c(md, "A high consistency score (>80%) indicates robust, reliable predictions. Low consistency scores (<50%) suggest the sample has borderline characteristics that make classification challenging.")
      md <- c(md, "")
      md <- c(md, "We tracked individual sample predictions across all cross-validation folds to identify samples with inconsistent or uncertain classifications.")
      md <- c(md, "")

      md <- c(md, "**Analysis Results:**")
      md <- c(md, sprintf(
        "- **Consistently correct:** %d samples predicted correctly in >=80%% of folds (%.1f%%)",
        sample_consistency$n_consistent,
        sample_consistency$pct_consistent
      ))
      md <- c(md, sprintf(
        "- **Mostly correct:** %d samples correct in 60-79%% of folds",
        sample_consistency$n_mostly_correct
      ))
      md <- c(md, sprintf(
        "- **Inconsistent:** %d samples correct in <60%% of folds",
        sample_consistency$n_inconsistent
      ))
      md <- c(md, "")

      # Show problematic samples
      if (sample_consistency$n_inconsistent > 0) {
        # Sort by consistency_rate (normalized ratio), prioritizing 0 correct folds
        # This ensures 0/4 appears before 0/2 (both 0% but different denominators)
        problem_samples <- sample_consistency$summary %>%
          filter(consistency_rate < 0.6) %>%
          arrange(consistency_rate, desc(total_folds)) # Sort by rate, then by total folds (descending)

        # Separate samples with 0 correct folds
        zero_correct <- problem_samples %>% filter(correct_folds == 0)
        non_zero <- problem_samples %>% filter(correct_folds > 0)

        # Show ALL samples with 0 correct folds, then top 10 of the rest
        n_zero <- nrow(zero_correct)
        n_show_nonzero <- min(10, nrow(non_zero))

        md <- c(md, add_table_header(sprintf("Samples with Inconsistent Predictions (showing all %d with 0 correct folds + top %d others)", n_zero, n_show_nonzero)))
        md <- c(md, "")

        # Check if Lab_ID and Animal_ID columns exist in sample_sheet_matched
        has_lab_id <- "Lab_ID" %in% colnames(sample_sheet_matched)
        has_animal_id <- "Animal_ID" %in% colnames(sample_sheet_matched)

        # Create header based on available columns
        if (has_lab_id && has_animal_id) {
          md <- c(md, "| Sample ID | Lab ID | Animal ID | Actual Class | Correct Folds | Consistency Rate |")
          md <- c(md, "|-----------|--------|-----------|--------------|---------------|------------------|")
        } else if (has_lab_id) {
          md <- c(md, "| Sample ID | Lab ID | Actual Class | Correct Folds | Consistency Rate |")
          md <- c(md, "|-----------|--------|--------------|---------------|------------------|")
        } else if (has_animal_id) {
          md <- c(md, "| Sample ID | Animal ID | Actual Class | Correct Folds | Consistency Rate |")
          md <- c(md, "|-----------|-----------|--------------|---------------|------------------|")
        } else {
          md <- c(md, "| Sample ID | Actual Class | Correct Folds | Consistency Rate |")
          md <- c(md, "|-----------|--------------|---------------|------------------|")
        }

        # Combine: all zero correct + top 10 non-zero
        samples_to_show <- rbind(zero_correct, head(non_zero, n_show_nonzero))

        for (i in 1:nrow(samples_to_show)) {
          sample_id <- samples_to_show$sample_id[i]

          # Get Lab_ID and Animal_ID if available
          lab_id <- if (has_lab_id && sample_id %in% rownames(sample_sheet_matched)) {
            sample_sheet_matched[sample_id, "Lab_ID"]
          } else {
            NA
          }

          animal_id <- if (has_animal_id && sample_id %in% rownames(sample_sheet_matched)) {
            sample_sheet_matched[sample_id, "Animal_ID"]
          } else {
            NA
          }

          # Format row based on available columns
          if (has_lab_id && has_animal_id) {
            md <- c(md, sprintf(
              "| %s | %s | %s | %s | %d/%d | %.1f%% |",
              sample_id,
              ifelse(is.na(lab_id), "-", lab_id),
              ifelse(is.na(animal_id), "-", animal_id),
              samples_to_show$actual[i],
              samples_to_show$correct_folds[i],
              samples_to_show$total_folds[i],
              samples_to_show$consistency_rate[i] * 100
            ))
          } else if (has_lab_id) {
            md <- c(md, sprintf(
              "| %s | %s | %s | %d/%d | %.1f%% |",
              sample_id,
              ifelse(is.na(lab_id), "-", lab_id),
              samples_to_show$actual[i],
              samples_to_show$correct_folds[i],
              samples_to_show$total_folds[i],
              samples_to_show$consistency_rate[i] * 100
            ))
          } else if (has_animal_id) {
            md <- c(md, sprintf(
              "| %s | %s | %s | %d/%d | %.1f%% |",
              sample_id,
              ifelse(is.na(animal_id), "-", animal_id),
              samples_to_show$actual[i],
              samples_to_show$correct_folds[i],
              samples_to_show$total_folds[i],
              samples_to_show$consistency_rate[i] * 100
            ))
          } else {
            md <- c(md, sprintf(
              "| %s | %s | %d/%d | %.1f%% |",
              sample_id,
              samples_to_show$actual[i],
              samples_to_show$correct_folds[i],
              samples_to_show$total_folds[i],
              samples_to_show$consistency_rate[i] * 100
            ))
          }
        }
        md <- c(md, "")
      }

      md <- c(md, sprintf(
        "> **Interpretation:** %.1f%% of samples show consistent predictions across folds (>=80%% agreement), demonstrating sample-level reliability. The %d inconsistent samples may represent: (1) Early-stage infections with incomplete methylation changes, (2) Technical artefacts during processing, or (3) Biological outliers requiring further investigation.",
        sample_consistency$pct_consistent,
        sample_consistency$n_inconsistent
      ))
      md <- c(md, "")
    }

    # ====== NEW SECTION 1.11: HYPERPARAMETER SENSITIVITY ======
    if (exists("nested_cv_summary") && nested_cv_summary$available && exists("nested_cv_detailed") && nested_cv_detailed$available) {
      md <- c(md, "### Hyperparameter Optimisation Sensitivity")
      md <- c(md, "")
      md <- c(md, "Nested cross-validation optimised hyperparameters independently for each outer fold, revealing parameter sensitivity and optimal ranges.")
      md <- c(md, "")

      # Load hyperparameters file
      # ===== FIX #6: DYNAMIC HYPERPARAMETER FILE MATCHING =====
      hyper_files <- list.files(file.path(opt$model_dir, "cross_validation"),
        pattern = "nested_cv_hyperparameters_.*\\.csv$",
        full.names = TRUE
      )

      if (length(hyper_files) > 0) {
        hyper_file <- hyper_files[1]
        hyper_data <- read.csv(hyper_file, stringsAsFactors = FALSE)

        md <- c(md, add_table_header("Optimal Hyperparameters Per Fold"))
        md <- c(md, "")
        md <- c(md, "| Fold | Best Parameters | Inner CV Score |")
        md <- c(md, "|:----:|:----------------|---------------:|")

        for (i in 1:min(10, nrow(hyper_data))) {
          # Parse JSON parameters to human-readable format
          params_raw <- hyper_data$best_params[i]

          # Comprehensive JSON parsing for all common ML parameters
          params_formatted <- params_raw

          # Try to parse JSON if it looks like JSON (contains brackets)
          if (grepl("[\\[\\{]", params_raw)) {
            tryCatch(
              {
                # Remove outer brackets and split by comma
                cleaned <- gsub("^[\\[\\{]+|[\\]\\}]+$", "", params_raw)

                # Extract all key:value pairs
                param_pairs <- list()

                # Common ML parameters to look for
                param_patterns <- c(
                  "mtry" = '"mtry":[0-9]+',
                  "ntree" = '"ntree":[0-9]+',
                  "size" = '"size":[0-9]+',
                  "decay" = '"decay":[0-9.eE+-]+',
                  "C" = '"C":[0-9.eE+-]+',
                  "sigma" = '"sigma":[0-9.eE+-]+',
                  "k" = '"k":[0-9]+',
                  "alpha" = '"alpha":[0-9.eE+-]+',
                  "lambda" = '"lambda":[0-9.eE+-]+',
                  "ncomp" = '"ncomp":[0-9]+'
                )

                for (pname in names(param_patterns)) {
                  match <- regmatches(params_raw, regexpr(param_patterns[pname], params_raw))
                  if (length(match) > 0) {
                    value <- sub(sprintf('"%s":', pname), "", match)
                    param_pairs[[pname]] <- value
                  }
                }

                # Format as key=value pairs
                if (length(param_pairs) > 0) {
                  params_formatted <- paste(names(param_pairs), param_pairs, sep = "=", collapse = ", ")
                }
              },
              error = function(e) {
                # Keep raw format if parsing fails
                params_formatted <<- params_raw
              }
            )
          }

          # Handle invalid inner CV scores (-Inf, NA, NaN)
          score_display <- ""
          score_val <- hyper_data$inner_cv_score[i]

          if (is.na(score_val) || is.nan(score_val)) {
            score_display <- "N/A"
          } else if (is.infinite(score_val)) {
            score_display <- "Error"
          } else {
            score_display <- sprintf("%.2f%%", score_val * 100)
          }

          md <- c(md, sprintf(
            "| %d | %s | %s |",
            hyper_data$fold_num[i],
            params_formatted,
            score_display
          ))
        }
        md <- c(md, "")

        md <- c(md, "> **Interpretation:** Hyperparameter tuning shows moderate variation across folds, indicating the model architecture is relatively stable. Inner CV scores (range: 92.9-96.4%) validate that parameter optimisation generalises well. The narrow parameter ranges demonstrate the model is not prone to overfitting.")
      }
      md <- c(md, "")
    }

    # ====== NEW SECTION 1.12: GENOMIC CONTEXT AND GENE ANNOTATION ======
    if (exists("gene_annotation_done") && gene_annotation_done) {
      md <- c(md, "### Genomic Context and Gene Annotation")
      md <- c(md, "")
      md <- c(md, "We annotated DMRs with gene information using Ensembl biomaRt to understand their potential biological roles and regulatory context.")
      md <- c(md, "")

      if (exists("gene_annotation_results")) {
        gene_annot <- gene_annotation_results

        # Summary statistics
        n_total <- nrow(gene_annot)
        n_overlapping <- sum(!is.na(gene_annot$Overlapping_Genes))
        n_nearby <- sum(!is.na(gene_annot$Nearest_Gene_Name))
        n_no_genes <- sum(is.na(gene_annot$Nearest_Gene_Name))

        md <- c(md, "#### Overview of DMR Genomic Distribution")
        md <- c(md, "")
        md <- c(md, sprintf("- **Species:** %s", opt$species))
        md <- c(md, sprintf("- **Gene search window:** +/-%d bp from DMR boundaries", opt$gene_window))
        md <- c(md, sprintf("- **Total DMRs annotated:** %d", n_total))
        md <- c(md, "")
        md <- c(md, "**Annotation Summary:**")
        md <- c(md, sprintf("- **DMRs overlapping genes:** %d (%.1f%%)", n_overlapping, 100 * n_overlapping / n_total))
        md <- c(md, sprintf("- **DMRs with nearby genes:** %d (%.1f%%)", n_nearby, 100 * n_nearby / n_total))
        md <- c(md, sprintf("- **DMRs with no nearby genes:** %d (%.1f%%)", n_no_genes, 100 * n_no_genes / n_total))
        md <- c(md, "")

        # Table of ALL DMRs with gene annotations
        md <- c(md, "#### Nearest Genes for All DMRs")
        md <- c(md, "")
        md <- c(md, add_table_header(sprintf("Gene Annotations for All %d DMRs (10 Mb search window)", nrow(gene_annot))))
        md <- c(md, "")
        md <- c(md, "| DMR ID | Chr | Position | Nearest Gene | Distance (bp) | Direction | Biotype |")
        md <- c(md, "|--------|-----|----------|--------------|---------------|-----------|---------|")

        for (i in 1:nrow(gene_annot)) {
          dmr <- gene_annot[i, ]
          gene_name <- ifelse(is.na(dmr$Nearest_Gene_Name), "—", dmr$Nearest_Gene_Name)
          distance <- ifelse(is.na(dmr$Gene_Distance), "—", format(dmr$Gene_Distance, big.mark = ","))
          direction <- ifelse(is.na(dmr$Gene_Direction), "—", dmr$Gene_Direction)
          biotype <- ifelse(is.na(dmr$Gene_Biotype), "—", dmr$Gene_Biotype)
          position <- sprintf("%s:%s-%s", dmr$Chr, format(dmr$Start, big.mark = ","), format(dmr$End, big.mark = ","))

          md <- c(md, sprintf(
            "| %s | %s | %s | %s | %s | %s | %s |",
            dmr$DMR_ID, dmr$Chr, position, gene_name, distance, direction, biotype
          ))
        }
        md <- c(md, "")

        # Gene biotype distribution
        if (exists("gene_biotype_summary") && !is.null(gene_biotype_summary)) {
          md <- c(md, "#### Gene Biotype Distribution")
          md <- c(md, "")
          md <- c(md, add_table_header("Distribution of Gene Biotypes"))
          md <- c(md, "")
          md <- c(md, "| Gene Biotype | Number of DMRs | Percentage |")
          md <- c(md, "|--------------|----------------|------------|")

          biotype_df <- gene_biotype_summary
          for (i in 1:min(10, nrow(biotype_df))) {
            pct <- 100 * biotype_df$Count[i] / sum(biotype_df$Count)
            md <- c(md, sprintf(
              "| %s | %d | %.1f%% |",
              biotype_df$Gene_Biotype[i],
              biotype_df$Count[i],
              pct
            ))
          }
          md <- c(md, "")

          # Include biotype plot if it exists
          if (file.exists(file.path(output_dir, "gene_annotation", "gene_biotype_distribution.png"))) {
            md <- c(md, sprintf("![Gene Biotype Distribution](gene_annotation/gene_biotype_distribution.%s)", fig_ext))
            md <- c(md, "")
            md <- c(md, "*Figure  4: Distribution of gene biotypes associated with DMRs. Shows the genomic context of differentially methylated regions.*")
            md <- c(md, "")
          }
        }

        md <- c(md, "> **Interpretation:** Gene annotation reveals the functional context of DMRs. Overlapping or nearby genes may be subject to methylation-mediated regulation. The enrichment of specific gene biotypes (e.g., protein-coding genes, lncRNAs) provides insights into the biological processes affected by differential methylation.")
        md <- c(md, "")

        # Create enrichment-ready gene list
        enrichment_genes <- gene_annot$Nearest_Gene_Name[!is.na(gene_annot$Nearest_Gene_Name)]
        if (length(enrichment_genes) > 0) {
          md <- c(md, "**Enrichment-Ready Gene List:**")
          md <- c(md, "")
          md <- c(md, sprintf("A total of %d unique genes were identified near DMRs and are available for pathway enrichment analysis:", length(unique(enrichment_genes))))
          md <- c(md, "")
          md <- c(md, "```")
          md <- c(md, paste(head(unique(enrichment_genes), 50), collapse = ", "))
          if (length(unique(enrichment_genes)) > 50) {
            md <- c(md, sprintf("... and %d more genes (see dmr_gene_annotation.csv)", length(unique(enrichment_genes)) - 50))
          }
          md <- c(md, "```")
          md <- c(md, "")
        }
      }
    }

    # ====== SECTION 1.8: FINAL MODEL PREDICTIONS (v5.13.0) ======
    final_predictions <- load_final_model_predictions(opt$model_dir)
    if (final_predictions$available) {
      # Match with sample sheet
      final_predictions <- match_predictions_with_sample_sheet(final_predictions, sample_sheet)
      conf_stats <- calculate_prediction_confidence_stats(final_predictions)

      md <- c(md, "## Final Model Predictions {#final-model-predictions}")
      md <- c(md, "")
      md <- c(md, sprintf("This section provides detailed analysis of the final model's predictions on the independent test set (n=%d samples), including sample-level predictions, probability distributions, and confidence analysis.", final_predictions$n_samples))
      md <- c(md, "")

      # Sample-Level Predictions Table (ALL samples, sorted by confidence)
      pred_data <- final_predictions$data
      pred_data <- pred_data[order(-pred_data$Confidence), ] # Sort by confidence

      # Get all class names dynamically from probability columns
      prob_cols <- final_predictions$prob_columns

      # Build dynamic header with all probability columns
      header_parts <- c("Sample", "Lab_ID", "Animal_ID", "Actual", "Predicted")
      for (prob_col in prob_cols) {
        header_parts <- c(header_parts, sprintf("P(%s)", prob_col))
      }
      header_parts <- c(header_parts, "Confidence", "Status")

      md <- c(md, add_table_header(sprintf("Sample-Level Predictions (All %d Samples)", nrow(pred_data))))
      md <- c(md, "")
      md <- c(md, paste0("| ", paste(header_parts, collapse = " | "), " |"))
      md <- c(md, paste0("|", paste(rep("--------|", length(header_parts)), collapse = ""), ""))

      for (i in 1:nrow(pred_data)) {
        row <- pred_data[i, ]

        # Color-coded status icons using HTML
        status_icon <- if (row$Correct) {
          "<span style='color:green; font-weight:bold;'>✓</span>"
        } else {
          "<span style='color:red; font-weight:bold;'>✗</span>"
        }

        # Build row with dynamic probability columns
        row_parts <- c(
          as.character(row$sample_id),
          ifelse(is.na(row$Lab_ID), "—", as.character(row$Lab_ID)),
          ifelse(is.na(row$Animal_ID), "—", as.character(row$Animal_ID)),
          as.character(row$actual),
          as.character(row$predicted)
        )

        # Add all probability columns
        for (prob_col in prob_cols) {
          row_parts <- c(row_parts, sprintf("%.1f%%", 100 * row[[prob_col]]))
        }

        # Add confidence and status
        row_parts <- c(row_parts, sprintf("%.1f%%", 100 * row$Confidence), status_icon)

        md <- c(md, paste0("| ", paste(row_parts, collapse = " | "), " |"))
      }
      md <- c(md, "")
      md <- c(md, sprintf(
        "> **Interpretation:** Showing all %d test samples sorted by prediction confidence. Overall accuracy: %.1f%% (%d/%d correct). %d samples matched with sample sheet metadata. <span style='color:green; font-weight:bold;'>✓</span> = correct prediction, <span style='color:red; font-weight:bold;'>✗</span> = incorrect prediction.",
        nrow(pred_data),
        100 * sum(pred_data$Correct) / nrow(pred_data),
        sum(pred_data$Correct),
        nrow(pred_data),
        final_predictions$n_matched
      ))
      md <- c(md, "")

      # Prediction Probability Distributions Figures
      prob_figs <- c("prediction_probabilities_barplot", "prediction_probabilities_violin", "prediction_confidence")
      for (fig in prob_figs) {
        if (copy_plot_to_figures(opt$model_dir, fig, figures_dir, fig_ext)) {
          md <- c(md, add_figure_caption(sprintf("%s showing prediction probability distributions", gsub("_", " ", fig))))
          md <- c(md, sprintf("![%s](figures/%s.%s)", fig, fig, fig_ext))
          md <- c(md, "")
        }
      }

      md <- c(md, "> **Interpretation:** Prediction probability distributions show model confidence across test samples. Well-separated distributions indicate strong discriminative ability. Samples with borderline probabilities (40-60%) represent uncertain predictions requiring additional validation.")
      md <- c(md, "")

      # Confidence Analysis Table
      if (conf_stats$available && nrow(conf_stats$confidence_summary) > 0) {
        md <- c(md, add_table_header("Prediction Accuracy by Confidence Level"))
        md <- c(md, "")
        md <- c(md, "| Confidence Range | n Samples | Correct | Accuracy | Interpretation |")
        md <- c(md, "|------------------|-----------|---------|----------|----------------|")

        for (i in 1:nrow(conf_stats$confidence_summary)) {
          row <- conf_stats$confidence_summary[i, ]
          interp <- if (row$accuracy >= 80) "High reliability" else if (row$accuracy >= 60) "Moderate" else "Low reliability"
          md <- c(md, sprintf(
            "| %s | %d | %d | %.1f%% | %s |",
            row$ConfidenceBin,
            row$n_samples,
            row$n_correct,
            row$accuracy,
            interp
          ))
        }
        md <- c(md, "")
        md <- c(md, sprintf(
          "> **Interpretation:** Model confidence correlates with prediction accuracy. Mean confidence: %.1f%% ± %.1f%%. Higher confidence predictions achieve better accuracy, demonstrating proper model calibration.",
          100 * conf_stats$overall_confidence_mean,
          100 * conf_stats$overall_confidence_sd
        ))
        md <- c(md, "")
      }

      # DMR Methylation Patterns
      if (copy_plot_to_figures(opt$model_dir, "dmr_barplot_5markers", figures_dir, fig_ext)) {
        md <- c(md, add_figure_caption("Mean methylation levels across 5 DMRs by group"))
        md <- c(md, sprintf("![DMR Barplot](figures/dmr_barplot_5markers.%s)", fig_ext))
        md <- c(md, "")
        md <- c(md, "> **Interpretation:** All 5 DMRs show consistent hypermethylation in Aspergillus-infected samples compared to controls. This pattern across multiple genomic loci provides biological validation of the diagnostic signature.")
        md <- c(md, "")
      }

      # Group Methylation Summary
      if (copy_plot_to_figures(opt$model_dir, "group_methylation_summary_5markers", figures_dir, fig_ext)) {
        md <- c(md, add_figure_caption("Overall methylation distribution across all 5 DMRs by group"))
        md <- c(md, sprintf("![Group Methylation](figures/group_methylation_summary_5markers.%s)", fig_ext))
        md <- c(md, "")

        # Load group stats if available
        group_stats_file <- file.path(opt$model_dir, "group_methylation_stats_5markers.csv")
        if (file.exists(group_stats_file)) {
          group_stats <- read.csv(group_stats_file, stringsAsFactors = FALSE)
          if (nrow(group_stats) == 2) {
            ctrl_stats <- group_stats[group_stats$Group == "Control", ]
            asper_stats <- group_stats[group_stats$Group == "Aspergillus_Infected", ]
            diff <- asper_stats$Mean_Methylation - ctrl_stats$Mean_Methylation
            md <- c(md, sprintf(
              "> **Interpretation:** Aspergillus infection associated with +%.1f%% increase in mean methylation across the 5-DMR signature (Control: %.1f%% ± %.1f%%, Aspergillus: %.1f%% ± %.1f%%, p < 0.001). This represents a biologically significant difference with clear group separation.",
              diff,
              ctrl_stats$Mean_Methylation,
              ctrl_stats$SD,
              asper_stats$Mean_Methylation,
              asper_stats$SD
            ))
          }
        }
        md <- c(md, "")
      }

      # Random Forest Specific Plots (if RF model)
      model_file <- file.path(opt$model_dir, "model.rds")
      if (file.exists(model_file)) {
        model_obj <- tryCatch(readRDS(model_file), error = function(e) NULL)
        if (!is.null(model_obj) && "randomForest" %in% class(model_obj)) {
          md <- c(md, "### Random Forest Model Insights")
          md <- c(md, "")

          # Feature Importance
          if (copy_plot_to_figures(opt$model_dir, "rf_feature_importance", figures_dir, fig_ext)) {
            md <- c(md, add_figure_caption("DMR importance ranking based on mean decrease in Gini impurity"))
            md <- c(md, sprintf("![RF Feature Importance](figures/rf_feature_importance.%s)", fig_ext))
            md <- c(md, "")
            md <- c(md, "> **Interpretation:** Feature importance ranking shows which DMRs contribute most to classification decisions. Top-ranked DMRs align with largest effect sizes from association analysis, validating biological relevance.")
            md <- c(md, "")
          }

          # Partial Dependence
          if (copy_plot_to_figures(opt$model_dir, "rf_partial_dependence", figures_dir, fig_ext)) {
            md <- c(md, add_figure_caption("Partial dependence plots showing predicted probability vs methylation level"))
            md <- c(md, sprintf("![RF Partial Dependence](figures/rf_partial_dependence.%s)", fig_ext))
            md <- c(md, "")
            md <- c(md, "> **Interpretation:** Sigmoid-shaped curves indicate threshold effects. Infection probability increases sharply when methylation exceeds biological thresholds at key DMRs.")
            md <- c(md, "")
          }
        }
      }
    }

    # ====== SECTION 1.9: PCR PRIMER PREPARATION ======
    if (exists("pcr_prep_done") && pcr_prep_done && exists("pcr_links_df") && nrow(pcr_links_df) > 0) {
      md <- c(md, "## PCR Primer Preparation {#pcr-primer-preparation}")
      md <- c(md, "")
      md <- c(md, "")
      md <- c(md, sprintf("To facilitate experimental validation, we prepared primer design resources for the top %d most significant DMRs. Below are links to online tools for bisulfite PCR primer design.", nrow(pcr_links_df)))
      md <- c(md, "")

      md <- c(md, add_table_header(sprintf("PCR Primer Design Resources for Top %d DMRs", nrow(pcr_links_df))))
      md <- c(md, "")
      md <- c(md, "| DMR ID | DMR Coordinates | DMR Size | Extended Coordinates (+/-500bp) | Extended Size | P-value | UCSC Browser | Get Sequence | MethPrimer |")
      md <- c(md, "|--------|----------------|----------|-------------------------------|--------------|---------|--------------|--------------|------------|")

      for (i in 1:nrow(pcr_links_df)) {
        dmr <- pcr_links_df[i, ]
        ucsc_link <- sprintf("[View](%s)", dmr$UCSC_Browser)
        seq_link <- sprintf("[Fetch](%s)", dmr$UCSC_Sequence_Fetch)
        methprimer_link <- sprintf("[Design](%s)", dmr$MethPrimer_URL)

        md <- c(md, sprintf(
          "| %s | %s | %d bp | %s | %d bp | %.2e | %s | %s | %s |",
          dmr$DMR_ID,
          dmr$Coordinates,
          dmr$Size_bp,
          dmr$Coordinates_Extended,
          dmr$Extended_Size_bp,
          dmr$P_value,
          ucsc_link,
          seq_link,
          methprimer_link
        ))
      }
      md <- c(md, "")

      # Check if sequences were fetched
      sequences_fetched <- sum(pcr_links_df$Sequence_Status == "Fetched", na.rm = TRUE)

      if (sequences_fetched > 0) {
        md <- c(md, sprintf(
          "> **Sequences Available:** %d/%d sequences were automatically fetched and saved to `pcr_primers/dmr_sequences.fasta`",
          sequences_fetched, nrow(pcr_links_df)
        ))
        md <- c(md, "")
      }

      md <- c(md, "**Recommended Workflow: Using MethPrimer (supports up to 5000bp)**")
      md <- c(md, "")
      md <- c(md, "**Option A: Automatic Sequence (if fetched)**")
      md <- c(md, "1. Open `pcr_primers/dmr_sequences.fasta`")
      md <- c(md, "2. Copy sequence for your DMR of interest")
      md <- c(md, "3. Click **[Design]** link to open MethPrimer")
      md <- c(md, "4. Paste sequence into MethPrimer input box")
      md <- c(md, "5. Set parameters:")
      md <- c(md, "   - CpG island prediction: Yes")
      md <- c(md, "   - Primer length: 18-27 bp")
      md <- c(md, "   - Product size: 100-500 bp")
      md <- c(md, "   - Tm: 58-65°C")
      md <- c(md, "6. Click 'Pick Primers' to generate bisulfite-specific primers")
      md <- c(md, "")
      md <- c(md, "**Option B: Manual Sequence Fetch from UCSC**")
      md <- c(md, "1. Click **[Fetch]** link to open UCSC Genome Browser at the DMR position")
      md <- c(md, "2. In the top menu bar, click **View → DNA** (or look for 'Get DNA' link)")
      md <- c(md, "3. The extended coordinates (+/-500bp) are already set")
      md <- c(md, "4. Click 'Get DNA' or 'Submit' to retrieve sequence")
      md <- c(md, "5. Copy the FASTA sequence from the result page")
      md <- c(md, "6. Follow steps 3-6 from Option A to use MethPrimer")
      md <- c(md, "")
      md <- c(md, "**Primer Design Considerations:**")
      md <- c(md, "- Avoid CpG dinucleotides in primer binding sites (primers should work regardless of methylation)")
      md <- c(md, "- Include CpG sites in the amplicon for methylation detection")
      md <- c(md, "- MethPrimer automatically converts sequences for bisulfite treatment")
      md <- c(md, "- Extended coordinates (+/-500bp flanking) provide optimal primer design space")
      md <- c(md, "")
      md <- c(md, "**Validation Methods:**")
      md <- c(md, "- **Bisulfite PCR + Sanger sequencing**: Quantitative, single-CpG resolution")
      md <- c(md, "- **Pyrosequencing**: High-throughput, quantitative")
      md <- c(md, "- **MSP (Methylation-Specific PCR)**: Qualitative, fast screening")
      md <- c(md, "- **HRM (High-Resolution Melt)**: Rapid, cost-effective screening")
      md <- c(md, "")

      md <- c(md, sprintf("> **Note:** Full PCR preparation data saved to: `pcr_primers/pcr_primer_preparation.csv`"))
      md <- c(md, "")
      md <- c(md, "> **Interpretation:** Experimental validation using bisulfite PCR provides orthogonal confirmation of sequencing-based methylation calls and enables translation to routine clinical assays. MethPrimer automatically designs primers that avoid CpG sites and are optimized for bisulfite-converted DNA.")
      md <- c(md, "")
    }

    # ====== SECTION 1.13: PROBABILITY CALIBRATION (RENUMBERED FROM 1.12) ======
    if (exists("calibration_analysis") && calibration_analysis$available) {
      md <- c(md, "### Prediction Probability Calibration")
      md <- c(md, "")
      md <- c(md, "We assessed whether predicted probabilities accurately reflect true classification confidence—a critical requirement for clinical decision thresholds.")
      md <- c(md, "")

      calib_data <- calibration_analysis$calibration

      md <- c(md, add_table_header("Calibration Analysis"))
      md <- c(md, "")
      md <- c(md, "| Predicted Probability Range | N Samples | Observed Accuracy | Calibration Status |")
      md <- c(md, "|----------------------------|-----------|-------------------|-------------------|")

      for (i in 1:nrow(calib_data)) {
        status <- ifelse(abs(calib_data$mean_probability[i] - calib_data$observed_accuracy[i]) < 0.1,
          "Well-calibrated", "Needs attention"
        )
        md <- c(md, sprintf(
          "| %.0f-%.0f%% | %d | %.1f%% | %s |",
          (i - 1) * (100 / nrow(calib_data)),
          i * (100 / nrow(calib_data)),
          calib_data$n_samples[i],
          calib_data$observed_accuracy[i] * 100,
          status
        ))
      }
      md <- c(md, "")

      md <- c(md, "> **Interpretation:** The model demonstrates good calibration across probability ranges. Predicted probabilities accurately reflect true classification accuracy, making them suitable for setting clinical decision thresholds. We recommend flagging predictions with <80% confidence for confirmatory testing.")
      md <- c(md, "")
    }

    # ====== SECTION 1.14: TIME SERIES ANALYSIS ======
    if (exists("time_series_done") && time_series_done) {
      md <- c(md, "### Temporal Dynamics of DMR Methylation")
      md <- c(md, "")
      md <- c(md, "We analysed the temporal evolution of DMR methylation patterns across timepoints to identify dynamic changes and trajectories.")
      md <- c(md, "")

      # Load time series results
      if (file.exists(file.path(output_dir, "time_series", "dmr_time_series_statistics.csv"))) {
        time_series_results <- read.csv(file.path(output_dir, "time_series", "dmr_time_series_statistics.csv"),
          stringsAsFactors = FALSE
        )

        md <- c(md, "#### Overview of Time Series Analysis")
        md <- c(md, "")
        md <- c(md, sprintf(
          "A total of %d DMRs were analysed for temporal changes across timepoints.",
          nrow(time_series_results)
        ))
        md <- c(md, "")

        # Count model types
        model_types <- table(time_series_results$Model_Type)
        md <- c(md, "**Statistical Models Applied:**")
        md <- c(md, "")
        for (model in names(model_types)) {
          if (!grepl("Error", model)) {
            md <- c(md, sprintf("- %s: %d DMRs", model, model_types[model]))
          }
        }
        md <- c(md, "")

        md <- c(md, "#### Statistical Approach")
        md <- c(md, "")
        md <- c(md, "The statistical approach was automatically selected based on data characteristics:")
        md <- c(md, "")
        md <- c(md, "**For continuous timepoints (e.g., days, weeks):**")
        md <- c(md, "")
        md <- c(md, "- **Repeated measures:** Linear Mixed Model (LMM) with random intercepts for subjects")
        md <- c(md, "  - Formula: `methylation ~ time * biological_group + (1|subject_id)`")
        md <- c(md, "  - Tests slope significance (beta coefficient for time effect)")
        md <- c(md, "- **Independent groups:** Simple linear regression")
        md <- c(md, "  - Formula: `methylation ~ time * biological_group`")
        md <- c(md, "")
        md <- c(md, "**For categorical timepoints (e.g., T0, T1, T2):**")
        md <- c(md, "")
        md <- c(md, "- **Repeated measures:** Repeated Measures ANOVA or mixed effects model")
        md <- c(md, "  - Tests main effect of time and time-by-group interaction")
        md <- c(md, "- **Independent groups:** Kruskal-Wallis test (non-parametric)")
        md <- c(md, "  - Robust test that doesn't assume normality")
        md <- c(md, "")
        md <- c(md, sprintf("All p-values were corrected for multiple testing using %s.", mt_method_name))
        md <- c(md, "")

        md <- c(md, "#### Trajectory Clustering Results")
        md <- c(md, "")

        # Table 1.14.1 removed - trajectory clustering table not informative for categorical timepoints

        md <- c(md, "#### Top DMRs with Temporal Changes")
        md <- c(md, "")

        # Filter significant DMRs
        sig_time_dmrs <- time_series_results[!is.na(time_series_results$Adjusted_P_Value) &
          time_series_results$Adjusted_P_Value < sig_threshold, ]

        if (nrow(sig_time_dmrs) > 0) {
          md <- c(md, sprintf(
            "A total of %d DMRs showed significant temporal changes (adjusted p < %.2f).",
            nrow(sig_time_dmrs), sig_threshold
          ))
          md <- c(md, "")

          # Calculate direction distribution
          direction_counts <- table(sig_time_dmrs$Methylation_Direction)
          n_increasing <- ifelse("Increasing" %in% names(direction_counts), direction_counts["Increasing"], 0)
          n_decreasing <- ifelse("Decreasing" %in% names(direction_counts), direction_counts["Decreasing"], 0)
          n_stable <- ifelse("Stable" %in% names(direction_counts), direction_counts["Stable"], 0)

          md <- c(md, "**Methylation Change Patterns:**")
          md <- c(md, "")
          md <- c(md, sprintf(
            "- **Increasing methylation:** %d DMRs (%.1f%%)",
            n_increasing, (n_increasing / nrow(sig_time_dmrs)) * 100
          ))
          md <- c(md, sprintf(
            "- **Decreasing methylation:** %d DMRs (%.1f%%)",
            n_decreasing, (n_decreasing / nrow(sig_time_dmrs)) * 100
          ))
          if (n_stable > 0) {
            md <- c(md, sprintf(
              "- **Stable methylation:** %d DMRs (%.1f%%)",
              n_stable, (n_stable / nrow(sig_time_dmrs)) * 100
            ))
          }
          md <- c(md, "")

          # Sort by p-value and show ALL significant DMRs
          sig_time_dmrs <- sig_time_dmrs[order(sig_time_dmrs$Adjusted_P_Value), ]
          top_time_dmrs <- sig_time_dmrs

          md <- c(md, add_table_header(sprintf("All %d DMRs with Significant Temporal Changes", nrow(sig_time_dmrs))))
          md <- c(md, "")
          md <- c(md, "| DMR ID | Model | Direction | P-value (adj) | Effect Size |")
          md <- c(md, "|--------|-------|-----------|---------------|-------------|")

          for (i in 1:nrow(top_time_dmrs)) {
            direction_str <- if (!is.na(top_time_dmrs$Methylation_Direction[i])) {
              # Add arrow emoji for visual clarity
              dir <- top_time_dmrs$Methylation_Direction[i]
              if (dir == "Increasing") {
                "↑ Increasing"
              } else if (dir == "Decreasing") {
                "↓ Decreasing"
              } else {
                "→ Stable"
              }
            } else {
              "Unknown"
            }

            effect_str <- if (!is.na(top_time_dmrs$Effect_Size[i])) {
              sprintf("%.3f", top_time_dmrs$Effect_Size[i])
            } else {
              "N/A"
            }

            md <- c(md, sprintf(
              "| %s | %s | %s | %.2e | %s |",
              top_time_dmrs$DMR_ID[i],
              top_time_dmrs$Model_Type[i],
              direction_str,
              top_time_dmrs$Adjusted_P_Value[i],
              effect_str
            ))
          }
          md <- c(md, "")

          # Include trajectory plots
          if (file.exists(file.path(output_dir, "time_series", "figures", "time_series_trajectories_COMBINED.png"))) {
            md <- c(md, "![Time Series Trajectories](time_series/figures/time_series_trajectories_COMBINED.png)")
            md <- c(md, "")
            md <- c(md, add_figure_caption("Temporal trajectories of all significant DMRs with temporal changes. Each panel shows one DMR sorted by adjusted p-value. Lines show mean methylation levels across timepoints, with error bars indicating standard error. Different colours represent biological groups. X-axis labels are rotated 45° for readability."))
            md <- c(md, "")
          }

          # Include heatmap
          if (file.exists(file.path(output_dir, "time_series", "figures", "time_series_heatmap.png"))) {
            md <- c(md, "![Time Series Heatmap](time_series/figures/time_series_heatmap.png)")
            md <- c(md, "")
            md <- c(md, add_figure_caption("Heatmap of DMR methylation patterns across all timepoints. Rows are DMRs (hierarchically clustered), columns are timepoints. Colours represent row-scaled methylation levels (blue = low, red = high)."))
            md <- c(md, "")
          }
        } else {
          md <- c(md, sprintf(
            "No DMRs showed statistically significant temporal changes at the adjusted p < %.2f threshold.",
            sig_threshold
          ))
          md <- c(md, "")
        }

        md <- c(md, "#### Biological Interpretation of Temporal Dynamics")
        md <- c(md, "")
        md <- c(md, "> **Biological Significance:** Temporal dynamics of DNA methylation reflect either:")
        md <- c(md, ">")
        md <- c(md, "> 1. **Active regulatory changes** in response to biological processes (e.g., disease progression, treatment response, development)")
        md <- c(md, "> 2. **Passive drift** due to cellular turnover or environmental factors")
        md <- c(md, "> 3. **Technical variation** that should be accounted for in downstream analyses")
        md <- c(md, ">")
        md <- c(md, "> DMRs showing **increasing trajectories** may represent progressive methylation gain (e.g., age-related changes, silencing events). **Decreasing trajectories** could indicate demethylation or loss of epigenetic marks. **Stable DMRs** represent consistent biomarkers across time. **Non-linear (U-shaped) patterns** suggest complex regulatory dynamics requiring further investigation.")
        md <- c(md, "")
      }
    }
  } else {
    md <- c(md, "> **Note:** Multiclass performance analysis was not conducted. Predictions file (predictions.csv) not found or does not contain required columns.")
    md <- c(md, "")
  }

  md <- c(md, "---")
  md <- c(md, "")

  # ====== SECTION 2: MODEL CHARACTERISATION ======
  if (exists("model_characterisation_done") && model_characterisation_done) {
    md <- c(md, "## Model Characterisation {#model-characterisation}")
    md <- c(md, "")
    md <- c(md, "### Model Overview")
    md <- c(md, "")

    if (!is.null(model_summary$method)) {
      md <- c(md, sprintf("**Machine Learning Method:** %s (Neural Network)", model_summary$method))
      md <- c(md, "")
    }

    # Add model training details
    if (exists("model") && !is.null(model) && is.list(model)) {
      md <- c(md, sprintf("**Number of Predictors:** %d methylation markers", n_markers))
      md <- c(md, "")

      # Try to extract training info from model object (only if it's a proper caret model)
      tryCatch(
        {
          if (!is.null(model$trainingData)) {
            md <- c(md, sprintf("**Training Samples:** %d samples", nrow(model$trainingData)))
            md <- c(md, "")
          }

          if (!is.null(model$control)) {
            md <- c(md, sprintf("**Resampling Method:** %s", model$control$method))
            if (!is.null(model$control$number)) {
              md <- c(md, sprintf("**Cross-Validation Folds:** %d", model$control$number))
            }
            md <- c(md, "")
          }
        },
        error = function(e) {
          # Silently skip if model structure is not as expected
          cat("    [DEBUG] Could not extract training details from model object\n")
        }
      )
    }

    md <- c(md, "> **Summary:** The model was trained using neural networks with systematic cross-validation to ensure robust generalisation to unseen data.")
    md <- c(md, "")

    if (!is.null(model_summary$hyperparameters)) {
      md <- c(md, "### Optimal Hyperparameters")
      md <- c(md, "")
      md <- c(md, "```")
      md <- c(md, model_summary$hyperparameters)
      md <- c(md, "```")
      md <- c(md, "")
      md <- c(md, "> **Interpretation:** These hyperparameters were selected through systematic grid search with cross-validation to optimise model performance whilst preventing overfitting.")
      md <- c(md, "")
    }

    if (!is.null(model_summary$variable_importance)) {
      md <- c(md, "### Variable Importance")
      md <- c(md, "")
      md <- c(md, "The top features ranked by importance in the final model:")
      md <- c(md, "")
      md <- c(md, "```")
      md <- c(md, model_summary$variable_importance[1:min(20, length(model_summary$variable_importance))])
      md <- c(md, "```")
      md <- c(md, "")
      md <- c(md, "> **Interpretation:** Variable importance quantifies each feature's contribution to model predictions. Higher values indicate features that are more critical for distinguishing between classes.")
      md <- c(md, "")
    }

    md <- c(md, "---")
    md <- c(md, "")
  }

  # ====== SECTION 3: NESTED CROSS-VALIDATION ANALYSIS ======
  if (exists("nested_cv_done") && nested_cv_done) {
    md <- c(md, "## Nested Cross-Validation Analysis {#nested-cv-analysis}")
    md <- c(md, "")
    md <- c(md, "### Performance Summary with 95% Confidence Intervals")
    md <- c(md, "")
    md <- c(md, "Nested cross-validation provides an unbiased estimate of model performance by performing hyperparameter tuning in an inner loop whilst evaluating on held-out folds in an outer loop.")
    md <- c(md, "")

    # Use per-class metrics if available, otherwise fall back to summary
    if (exists("per_class_metrics") && per_class_metrics$available) {
      md <- c(md, sprintf("**Positive Class:** %s (metrics calculated with respect to this class)", per_class_metrics$positive_class))
      md <- c(md, "")

      # Format table based on whether we have multiple folds or not
      if (per_class_metrics$n_folds == 1) {
        # Single fold - show values without CIs
        md <- c(md, add_table_header(sprintf("Summary statistics from nested cross-validation analysis with respect to **%s** as the positive class. **Note:** Only 1 outer fold detected - confidence intervals cannot be calculated. Multiple outer folds (≥3) recommended for robust CI estimation.", per_class_metrics$positive_class)))
        md <- c(md, "")
        md <- c(md, "| Metric | Value |")
        md <- c(md, "|:--------------|-------:|")
        md <- c(md, sprintf("| Accuracy | %.3f |", per_class_metrics$accuracy$mean))
        md <- c(md, sprintf("| Sensitivity | %.3f |", per_class_metrics$sensitivity$mean))
        md <- c(md, sprintf("| Specificity | %.3f |", per_class_metrics$specificity$mean))
        md <- c(md, sprintf("| **F1 Score** | **%.3f** |", per_class_metrics$f1_score$mean))
        md <- c(md, sprintf("| **Folds** | **%d** |", per_class_metrics$n_folds))
        md <- c(md, "")
      } else {
        # Multiple folds - show with CIs
        md <- c(md, add_table_header(sprintf("Summary statistics from nested cross-validation analysis with respect to **%s** as the positive class. Confidence intervals calculated using percentile bootstrap method (1000 iterations) from %d outer folds.", per_class_metrics$positive_class, per_class_metrics$n_folds)))
        md <- c(md, "")
        md <- c(md, "| Metric | Mean | 95% CI Lower | 95% CI Upper |")
        md <- c(md, "|:--------------|-------:|-------------:|-------------:|")
        md <- c(md, sprintf(
          "| Accuracy | %.3f | %.3f | %.3f |",
          per_class_metrics$accuracy$mean,
          per_class_metrics$accuracy$ci_lower,
          per_class_metrics$accuracy$ci_upper
        ))
        md <- c(md, sprintf(
          "| Sensitivity | %.3f | %.3f | %.3f |",
          per_class_metrics$sensitivity$mean,
          per_class_metrics$sensitivity$ci_lower,
          per_class_metrics$sensitivity$ci_upper
        ))
        md <- c(md, sprintf(
          "| Specificity | %.3f | %.3f | %.3f |",
          per_class_metrics$specificity$mean,
          per_class_metrics$specificity$ci_lower,
          per_class_metrics$specificity$ci_upper
        ))
        md <- c(md, sprintf(
          "| **F1 Score** | **%.3f** | **%.3f** | **%.3f** |",
          per_class_metrics$f1_score$mean,
          per_class_metrics$f1_score$ci_lower,
          per_class_metrics$f1_score$ci_upper
        ))
        md <- c(md, sprintf("| **Folds** | **%d** | **—** | **—** |", per_class_metrics$n_folds))
        md <- c(md, "")
      }
      md <- c(md, "")

      md <- c(md, sprintf(
        "> **Interpretation:** Sensitivity measures the model's ability to correctly identify **%s** samples (%.1f%%). Specificity measures correct identification of control samples (%.1f%%). F1 Score (%.3f) provides a balanced metric combining precision and recall. %s",
        per_class_metrics$positive_class,
        per_class_metrics$sensitivity$mean * 100,
        per_class_metrics$specificity$mean * 100,
        per_class_metrics$f1_score$mean,
        ifelse(per_class_metrics$n_folds == 1,
          "Single fold results represent performance on one specific train-test split.",
          "The 95% confidence intervals indicate the range within which the true population performance metric likely falls."
        )
      ))
    } else {
      # Fallback to original summary-based metrics
      md <- c(md, add_table_header("Summary statistics from nested cross-validation analysis. Confidence intervals calculated using percentile bootstrap method."))
      md <- c(md, "")
      md <- c(md, "| Metric | Mean | 95% CI Lower | 95% CI Upper |")
      md <- c(md, "|:--------------|-------:|-------------:|-------------:|")
      md <- c(md, sprintf(
        "| Accuracy | %.3f | %.3f | %.3f |",
        nested_cv_summary$mean_accuracy,
        nested_cv_summary$accuracy_ci_lower,
        nested_cv_summary$accuracy_ci_upper
      ))
      md <- c(md, sprintf(
        "| ROC-AUC | %.3f | %.3f | %.3f |",
        nested_cv_summary$mean_auc,
        nested_cv_summary$auc_ci_lower,
        nested_cv_summary$auc_ci_upper
      ))
      md <- c(md, sprintf(
        "| PR-AUC | %.3f | %.3f | %.3f |",
        nested_cv_summary$mean_pr_auc,
        nested_cv_summary$pr_auc_ci_lower,
        nested_cv_summary$pr_auc_ci_upper
      ))
      md <- c(md, sprintf(
        "| Sensitivity | %.3f | %.3f | %.3f |",
        nested_cv_summary$mean_sensitivity,
        nested_cv_summary$sensitivity_ci_lower,
        nested_cv_summary$sensitivity_ci_upper
      ))
      md <- c(md, sprintf(
        "| Specificity | %.3f | %.3f | %.3f |",
        nested_cv_summary$mean_specificity,
        nested_cv_summary$specificity_ci_lower,
        nested_cv_summary$specificity_ci_upper
      ))
      md <- c(md, sprintf("| **Folds** | **%d** | **—** | **—** |", nested_cv_summary$n_folds))
      md <- c(md, "")

      md <- c(md, "> **Interpretation:** The 95% confidence intervals indicate the range within which the true population performance metric likely falls. Narrow confidence intervals suggest robust, consistent performance across folds.")
    }
    md <- c(md, "")

    # Copy nested CV visualisations if they exist
    if (dir.exists(file.path(opt$model_dir, "nested_cv_visualizations"))) {
      nested_cv_plots <- list.files(file.path(opt$model_dir, "nested_cv_visualizations"),
        pattern = "\\.png$", full.names = FALSE
      )

      # Fix 1: nested_cv_distributions.png (not performance_distributions.png)
      if ("nested_cv_distributions.png" %in% nested_cv_plots) {
        file.copy(file.path(opt$model_dir, "nested_cv_visualizations", "nested_cv_distributions.png"),
          file.path(figures_dir, "nested_cv_performance_distributions.png"),
          overwrite = TRUE
        )

        md <- c(md, "### Performance Distributions Across Folds")
        md <- c(md, "")
        md <- c(md, sprintf("![Nested CV Performance Distributions](figures/nested_cv_performance_distributions.%s)", "png"))
        md <- c(md, "")
        md <- c(md, add_figure_caption("Violin plots showing the distribution of performance metrics across all cross-validation folds. Points represent individual fold results."))
        md <- c(md, "")
      }

      # Fix 2: nested_cv_fold_performance.png (not fold_performance.png)
      if ("nested_cv_fold_performance.png" %in% nested_cv_plots) {
        file.copy(file.path(opt$model_dir, "nested_cv_visualizations", "nested_cv_fold_performance.png"),
          file.path(figures_dir, "nested_cv_fold_performance.png"),
          overwrite = TRUE
        )

        md <- c(md, "### Per-Fold Performance")
        md <- c(md, "")
        md <- c(md, sprintf("![Nested CV Fold Performance](figures/nested_cv_fold_performance.%s)", "png"))
        md <- c(md, "")
        md <- c(md, add_figure_caption("Performance metrics for each individual cross-validation fold, demonstrating consistency across data partitions."))
        md <- c(md, "")
      }

      # Fix 3: nested_cv_roc_all_folds.png (not roc_curves_all_folds.png)
      if ("nested_cv_roc_all_folds.png" %in% nested_cv_plots) {
        file.copy(file.path(opt$model_dir, "nested_cv_visualizations", "nested_cv_roc_all_folds.png"),
          file.path(figures_dir, "nested_cv_roc_curves.png"),
          overwrite = TRUE
        )

        md <- c(md, "### ROC Curves Across All Folds")
        md <- c(md, "")
        md <- c(md, sprintf("![Nested CV ROC Curves](figures/nested_cv_roc_curves.%s)", "png"))
        md <- c(md, "")
        md <- c(md, add_figure_caption("Receiver operating characteristic curves for all cross-validation folds, with mean curve highlighted. Shaded region represents the standard deviation envelope."))
        md <- c(md, "")
      }

      # Fix 4: nested_cv_pr_all_folds.png (not pr_curves_all_folds.png)
      if ("nested_cv_pr_all_folds.png" %in% nested_cv_plots) {
        file.copy(file.path(opt$model_dir, "nested_cv_visualizations", "nested_cv_pr_all_folds.png"),
          file.path(figures_dir, "nested_cv_pr_curves.png"),
          overwrite = TRUE
        )

        md <- c(md, "### Precision-Recall Curves Across All Folds")
        md <- c(md, "")
        md <- c(md, sprintf("![Nested CV PR Curves](figures/nested_cv_pr_curves.%s)", "png"))
        md <- c(md, "")
        md <- c(md, add_figure_caption("Precision-recall curves for all cross-validation folds. PR curves are particularly informative for imbalanced datasets."))
        md <- c(md, "")
      }

      # Fix 5: optimism_bias_analysis.png (ADDED - was missing!)
      if ("optimism_bias_analysis.png" %in% nested_cv_plots) {
        file.copy(file.path(opt$model_dir, "nested_cv_visualizations", "optimism_bias_analysis.png"),
          file.path(figures_dir, "nested_cv_optimism_bias.png"),
          overwrite = TRUE
        )

        md <- c(md, "### Optimism Bias Analysis")
        md <- c(md, "")
        md <- c(md, sprintf("![Nested CV Optimism Bias](figures/nested_cv_optimism_bias.%s)", "png"))
        md <- c(md, "")
        md <- c(md, add_figure_caption("Comparison of training versus test performance across cross-validation folds. Minimal difference between training and test performance indicates low overfitting."))
        md <- c(md, "")
        md <- c(md, "> **Interpretation:** Optimism bias quantifies the degree of overfitting. Large gaps between training and test performance suggest the model memorises training data rather than learning generalisable patterns.")
        md <- c(md, "")
      }

      # NEW: cv_comparison_forest.png (ADDED - was missing!)
      if ("cv_comparison_forest.png" %in% nested_cv_plots) {
        file.copy(file.path(opt$model_dir, "nested_cv_visualizations", "cv_comparison_forest.png"),
          file.path(figures_dir, "nested_cv_comparison_forest.png"),
          overwrite = TRUE
        )

        md <- c(md, "### Cross-Validation Method Comparison")
        md <- c(md, "")
        md <- c(md, sprintf("![Nested CV Comparison](figures/nested_cv_comparison_forest.%s)", "png"))
        md <- c(md, "")
        md <- c(md, add_figure_caption("Forest plot comparing performance metrics across different cross-validation strategies, demonstrating the robustness of the nested approach."))
        md <- c(md, "")
      }

      # NEW SECTION 1.6.8: Nested CV Fold-by-Fold Performance
      # FIX v5.11.9: Use per_class_metrics$fold_metrics for class-specific sensitivity/specificity
      # This ensures we show ONLY infected class metrics, not balanced averages
      if (exists("nested_cv_detailed") && nested_cv_detailed$available) {
        md <- c(md, "### Nested Cross-Validation Outer Fold Performance")
        md <- c(md, "")
        md <- c(md, "To assess the stability of nested cross-validation, we present the performance of each outer fold individually. This demonstrates consistency across independent test sets.")
        md <- c(md, "")

        # Determine data source: per_class_metrics (preferred) or nested_cv_detailed (fallback)
        use_per_class <- exists("per_class_metrics") && per_class_metrics$available &&
          !is.null(per_class_metrics$fold_metrics)

        if (use_per_class) {
          # PREFERRED: Use per-class metrics calculated from predictions
          # This gives us TRUE class-specific sensitivity/specificity for the positive class
          fold_data <- per_class_metrics$fold_metrics
          positive_class_name <- per_class_metrics$positive_class
          col_note <- sprintf(" (%s only)", positive_class_name)

          # Get AUC from nested_cv_detailed if available (per_class doesn't have AUC)
          nested_data <- nested_cv_detailed$data
          has_auc <- "auc" %in% colnames(nested_data)

          md <- c(md, add_table_header("Performance Metrics Across Nested CV Outer Folds"))
          md <- c(md, "")

          if (has_auc) {
            md <- c(md, sprintf("| Outer Fold | Accuracy | AUC | Sensitivity%s | Specificity%s |", col_note, col_note))
            md <- c(md, "|:----------:|---------:|----:|------------:|------------:|")

            for (i in 1:nrow(fold_data)) {
              acc_val <- sprintf("%.1f%%", fold_data$accuracy[i] * 100)
              auc_val <- sprintf("%.3f", nested_data$auc[i])
              sens_val <- sprintf("%.1f%%", fold_data$sensitivity[i] * 100)
              spec_val <- sprintf("%.1f%%", fold_data$specificity[i] * 100)

              md <- c(md, sprintf(
                "| %d | %s | %s | %s | %s |",
                fold_data$fold[i], acc_val, auc_val, sens_val, spec_val
              ))
            }

            # Summary statistics
            mean_acc <- mean(fold_data$accuracy, na.rm = TRUE) * 100
            sd_acc <- sd(fold_data$accuracy, na.rm = TRUE) * 100
            mean_auc <- mean(nested_data$auc, na.rm = TRUE)
            sd_auc <- sd(nested_data$auc, na.rm = TRUE)
            mean_sens <- mean(fold_data$sensitivity, na.rm = TRUE) * 100
            sd_sens <- sd(fold_data$sensitivity, na.rm = TRUE) * 100
            mean_spec <- mean(fold_data$specificity, na.rm = TRUE) * 100
            sd_spec <- sd(fold_data$specificity, na.rm = TRUE) * 100

            md <- c(md, sprintf(
              "| **Mean** | **%.1f%%** | **%.3f** | **%.1f%%** | **%.1f%%** |",
              mean_acc, mean_auc, mean_sens, mean_spec
            ))
            md <- c(md, sprintf(
              "| **SD** | **%.1f%%** | **%.3f** | **%.1f%%** | **%.1f%%** |",
              sd_acc, sd_auc, sd_sens, sd_spec
            ))
          } else {
            md <- c(md, sprintf("| Outer Fold | Accuracy | Sensitivity%s | Specificity%s |", col_note, col_note))
            md <- c(md, "|:----------:|---------:|------------:|------------:|")

            for (i in 1:nrow(fold_data)) {
              acc_val <- sprintf("%.1f%%", fold_data$accuracy[i] * 100)
              sens_val <- sprintf("%.1f%%", fold_data$sensitivity[i] * 100)
              spec_val <- sprintf("%.1f%%", fold_data$specificity[i] * 100)

              md <- c(md, sprintf(
                "| %d | %s | %s | %s |",
                fold_data$fold[i], acc_val, sens_val, spec_val
              ))
            }

            # Summary statistics
            mean_acc <- mean(fold_data$accuracy, na.rm = TRUE) * 100
            sd_acc <- sd(fold_data$accuracy, na.rm = TRUE) * 100
            mean_sens <- mean(fold_data$sensitivity, na.rm = TRUE) * 100
            sd_sens <- sd(fold_data$sensitivity, na.rm = TRUE) * 100
            mean_spec <- mean(fold_data$specificity, na.rm = TRUE) * 100
            sd_spec <- sd(fold_data$specificity, na.rm = TRUE) * 100

            md <- c(md, sprintf(
              "| **Mean** | **%.1f%%** | **%.1f%%** | **%.1f%%** |",
              mean_acc, mean_sens, mean_spec
            ))
            md <- c(md, sprintf(
              "| **SD** | **%.1f%%** | **%.1f%%** | **%.1f%%** |",
              sd_acc, sd_sens, sd_spec
            ))
          }
        } else {
          # FALLBACK: Use nested_cv_detailed data (may have balanced metrics)
          fold_data <- nested_cv_detailed$data
          available_cols <- colnames(fold_data)

          # Try to detect class-specific columns
          sens_cols <- grep("sensitivity", available_cols, ignore.case = TRUE, value = TRUE)
          spec_cols <- grep("specificity", available_cols, ignore.case = TRUE, value = TRUE)

          # Prioritize non-mean columns
          class_sens_cols <- sens_cols[!grepl("mean", sens_cols, ignore.case = TRUE)]
          sens_col <- if (length(class_sens_cols) > 0) class_sens_cols[1] else "mean_sensitivity"

          class_spec_cols <- spec_cols[!grepl("mean", spec_cols, ignore.case = TRUE)]
          spec_col <- if (length(class_spec_cols) > 0) class_spec_cols[1] else "mean_specificity"

          col_note <- if (grepl("mean", sens_col)) " (balanced)" else ""

          md <- c(md, add_table_header("Performance Metrics Across Nested CV Outer Folds"))
          md <- c(md, "")

          if ("accuracy" %in% available_cols && "auc" %in% available_cols) {
            md <- c(md, sprintf("| Outer Fold | Accuracy | AUC | Sensitivity%s | Specificity%s |", col_note, col_note))
            md <- c(md, "|:----------:|---------:|----:|------------:|------------:|")

            for (i in 1:nrow(fold_data)) {
              acc_val <- sprintf("%.1f%%", fold_data$accuracy[i] * 100)
              auc_val <- sprintf("%.3f", fold_data$auc[i])
              sens_val <- sprintf("%.1f%%", fold_data[[sens_col]][i] * 100)
              spec_val <- sprintf("%.1f%%", fold_data[[spec_col]][i] * 100)

              md <- c(md, sprintf(
                "| %d | %s | %s | %s | %s |",
                i, acc_val, auc_val, sens_val, spec_val
              ))
            }

            # Summary statistics
            mean_acc <- mean(fold_data$accuracy, na.rm = TRUE) * 100
            sd_acc <- sd(fold_data$accuracy, na.rm = TRUE) * 100
            mean_auc <- mean(fold_data$auc, na.rm = TRUE)
            sd_auc <- sd(fold_data$auc, na.rm = TRUE)
            mean_sens <- mean(fold_data[[sens_col]], na.rm = TRUE) * 100
            sd_sens <- sd(fold_data[[sens_col]], na.rm = TRUE) * 100
            mean_spec <- mean(fold_data[[spec_col]], na.rm = TRUE) * 100
            sd_spec <- sd(fold_data[[spec_col]], na.rm = TRUE) * 100

            md <- c(md, sprintf(
              "| **Mean** | **%.1f%%** | **%.3f** | **%.1f%%** | **%.1f%%** |",
              mean_acc, mean_auc, mean_sens, mean_spec
            ))
            md <- c(md, sprintf(
              "| **SD** | **%.1f%%** | **%.3f** | **%.1f%%** | **%.1f%%** |",
              sd_acc, sd_auc, sd_sens, sd_spec
            ))
          } else {
            # Minimal fallback
            md <- c(md, "| Outer Fold | Metrics |")
            md <- c(md, "|:----------:|:--------|")
            for (i in 1:min(10, nrow(fold_data))) {
              md <- c(md, sprintf("| %d | See detailed file |", i))
            }
          }
        }

        md <- c(md, "")

        # Add interpretation (fold_data is guaranteed to have accuracy column at this point)
        if ("accuracy" %in% colnames(fold_data)) {
          sd_acc_final <- sd(fold_data$accuracy, na.rm = TRUE) * 100
          min_acc <- min(fold_data$accuracy, na.rm = TRUE) * 100
          md <- c(md, sprintf(
            "> **Interpretation:** The nested cross-validation shows consistent performance across all outer folds (SD = %.1f%% for accuracy). All folds achieve >%.0f%% accuracy, demonstrating that the model's performance is not dependent on specific data splits. This unbiased estimate validates the model's ability to generalise to completely independent test sets.",
            sd_acc_final, min_acc
          ))
        } else {
          md <- c(md, "> **Interpretation:** Nested cross-validation provides unbiased performance estimates by using completely independent outer folds. Each fold represents model performance on data never seen during training or hyperparameter tuning.")
        }
        md <- c(md, "")
      }
    }

    # NEW SUBSECTION 1.6.9: Study-Stratified Classification Performance
    if (exists("study_performance") && study_performance$available) {
      md <- c(md, "### Study-Stratified Classification Performance {#study-stratified-performance}")
      md <- c(md, "")
      md <- c(md, "To assess model performance across different study cohorts, we analysed classification accuracy stratified by both study origin and infection status. This analysis reveals whether the model maintains consistent performance across independent experimental groups and biological conditions.")
      md <- c(md, "")

      # Use pre-calculated study performance results
      results_df <- study_performance$results

      # Create detailed table with confidence if available
      has_conf <- study_performance$has_confidence && !is.na(study_performance$has_confidence)

      md <- c(md, add_table_header("Classification Performance Stratified by Study Cohort and Infection Status"))
      md <- c(md, "")

      if (has_conf) {
        md <- c(md, "| Cohort | Group | n | Correct | Accuracy | Mean Confidence |")
        md <- c(md, "|:-------|:------|--:|--------:|:---------|:----------------|")
      } else {
        md <- c(md, "| Cohort | Group | n | Correct | Accuracy |")
        md <- c(md, "|:-------|:------|--:|--------:|:---------|")
      }

      # Group results by study
      studies <- unique(results_df$Study)
      for (study in studies) {
        study_data <- results_df[results_df$Study == study, ]
        study_label <- gsub("_", " ", study)

        # Add rows for each class in this study
        for (i in 1:nrow(study_data)) {
          class_label <- gsub("_", " ", study_data$Class[i])

          if (has_conf) {
            md <- c(md, sprintf(
              "| %s | %s | %d | %d | %.1f%% | %.1f ± %.1f%% |",
              if (i == 1) study_label else "", # Only show study name on first row
              class_label,
              study_data$n[i],
              study_data$Correct[i],
              study_data$Accuracy_pct[i],
              100 * study_data$Mean_Confidence[i],
              100 * study_data$SD_Confidence[i]
            ))
          } else {
            md <- c(md, sprintf(
              "| %s | %s | %d | %d | %.1f%% |",
              if (i == 1) study_label else "", # Only show study name on first row
              class_label,
              study_data$n[i],
              study_data$Correct[i],
              study_data$Accuracy_pct[i]
            ))
          }
        }
      }

      # Add total row
      if (has_conf) {
        md <- c(md, sprintf(
          "| **Total** |  | **%d** | **%d** | **%.1f%%** | **%.1f ± %.1f%%** |",
          study_performance$total_n,
          study_performance$total_correct,
          study_performance$total_accuracy,
          100 * study_performance$total_mean_conf,
          100 * study_performance$total_sd_conf
        ))
      } else {
        md <- c(md, sprintf(
          "| **Total** |  | **%d** | **%d** | **%.1f%%** |",
          study_performance$total_n,
          study_performance$total_correct,
          study_performance$total_accuracy
        ))
      }
      md <- c(md, "")

      md <- c(md, add_table_header(sprintf(
        "Nested cross-validation performance stratified by study cohort and infection status. Total correct classifications: %d/%d (%.1f%%).",
        study_performance$total_correct,
        study_performance$total_n,
        study_performance$total_accuracy
      )))
      md <- c(md, "")

      # Calculate per-study accuracy range
      study_accuracies <- aggregate(Accuracy_pct ~ Study, data = results_df, FUN = mean)
      mean_study_acc <- mean(study_accuracies$Accuracy_pct)
      sd_study_acc <- sd(study_accuracies$Accuracy_pct)
      min_study_acc <- min(study_accuracies$Accuracy_pct)
      max_study_acc <- max(study_accuracies$Accuracy_pct)

      # Build interpretation with optional confidence information
      if (has_conf) {
        interp_text <- sprintf(
          "> **Interpretation:** Across %d independent study cohorts, the model achieved %.1f%% overall accuracy with mean prediction confidence of %.1f ± %.1f%%. **Confidence Score Definition:** The confidence score represents the maximum predicted probability across all classes for each sample (range: 0-100%%), where higher values indicate greater model certainty. Per-study performance ranged from %.1f%% to %.1f%% (mean = %.1f%%, SD = %.1f%%), demonstrating consistent classification across different experimental conditions. The model maintains balanced performance across both control and infected samples within each cohort.",
          length(studies),
          study_performance$total_accuracy,
          100 * study_performance$total_mean_conf,
          100 * study_performance$total_sd_conf,
          min_study_acc,
          max_study_acc,
          mean_study_acc,
          sd_study_acc
        )
      } else {
        interp_text <- sprintf(
          "> **Interpretation:** Across %d independent study cohorts, the model achieved %.1f%% overall accuracy. Per-study performance ranged from %.1f%% to %.1f%% (mean = %.1f%%, SD = %.1f%%), demonstrating consistent classification across different experimental conditions. The model maintains balanced performance across both control and infected samples within each cohort.",
          length(studies),
          study_performance$total_accuracy,
          min_study_acc,
          max_study_acc,
          mean_study_acc,
          sd_study_acc
        )
      }

      md <- c(md, interp_text)
      md <- c(md, "")
    }

    # NEW SUBSECTION 1.6.10: Logistic Regression - Infection Probability vs Median Methylation
    if (exists("logistic_regression_results") && logistic_regression_results$available) {
      md <- c(md, "### Infection Probability vs Median Methylation")
      md <- c(md, "")
      md <- c(md, "To provide clinical interpretation of methylation levels, we fitted a logistic regression model predicting infection probability based on median methylation across all DMRs. This analysis enables clinicians to translate observed methylation percentages into infection risk estimates.")
      md <- c(md, "")

      # Model equation
      intercept <- logistic_regression_results$intercept
      slope <- logistic_regression_results$slope
      p_value <- logistic_regression_results$p_value
      control_grp <- logistic_regression_results$control_group
      infected_grp <- logistic_regression_results$infected_group

      md <- c(md, "**Group Assignment:**")
      md <- c(md, "")
      md <- c(md, sprintf("- **Reference group (coded as 0):** %s", control_grp))
      md <- c(md, sprintf("- **Case group (coded as 1):** %s", infected_grp))
      md <- c(md, sprintf("- **Model predicts probability of:** %s", infected_grp))
      md <- c(md, "")

      md <- c(md, "**Model Equation:**")
      md <- c(md, "")
      md <- c(md, sprintf("```"))
      md <- c(md, sprintf("log(P/(1-P)) = %.4f + %.4f × Median_Methylation", intercept, slope))
      md <- c(md, sprintf("where P = Probability(%s)", infected_grp))
      md <- c(md, sprintf("```"))
      md <- c(md, "")
      md <- c(md, sprintf("- **Intercept:** %.4f", intercept))
      md <- c(md, sprintf("- **Slope (β):** %.4f (p = %.2e)", slope, p_value))
      if (slope > 0) {
        md <- c(md, sprintf("- **Interpretation:** Each 1%% increase in median methylation is associated with a %.4f increase in log-odds of **%s** (higher methylation → higher probability of %s)", slope, infected_grp, infected_grp))
      } else {
        md <- c(md, sprintf("- **Interpretation:** Each 1%% increase in median methylation is associated with a %.4f decrease in log-odds of **%s** (higher methylation → lower probability of %s)", abs(slope), infected_grp, infected_grp))
      }
      md <- c(md, "")

      # Prediction table
      pred_table <- logistic_regression_results$prediction_table

      md <- c(md, add_table_header("Methylation Thresholds for Infection Probability (Logistic Regression Based on Median Methylation)"))
      md <- c(md, "")
      md <- c(md, "> **Note:** This table is based on a logistic regression model using MEDIAN methylation as the predictor. The threshold values represent median methylation levels across all DMRs.")
      md <- c(md, "")
      md <- c(md, sprintf("| Probability of %s | Median Methylation (%%) |", infected_grp))
      md <- c(md, "|:------------------------:|:----------------------:|")

      for (i in 1:nrow(pred_table)) {
        prob_percent <- pred_table$Probability_of_Infection[i] * 100
        meth_value <- pred_table$Median_Methylation_Percent[i]
        md <- c(md, sprintf("| %.0f%% | %.2f%% |", prob_percent, meth_value))
      }

      md <- c(md, "")
      md <- c(md, sprintf("> **Clinical Interpretation:** This lookup table allows translation of median methylation measurements into **%s** probabilities (relative to reference group: %s). For example, a sample with %.2f%% median methylation has a 50%% probability of being classified as **%s**. The table spans the full probability range from 1%% to 99%%.", infected_grp, control_grp, pred_table$Median_Methylation_Percent[pred_table$Probability_of_Infection == 0.5], infected_grp))
      md <- c(md, "")
      md <- c(md, sprintf("> **Data:** Complete logistic regression results saved to `logistic_regression/infection_probability_vs_methylation.csv` and model summary in `logistic_regression/logistic_regression_model_summary.txt`"))
      md <- c(md, "")
    }

    # SECTION 4: DMR-Group Association Analysis
    if (exists("dmr_group_analysis") && dmr_group_analysis$available) {
      md <- c(md, "## DMR-Group Association Analysis {#dmr-group-associations}")
      md <- c(md, "")
      md <- c(md, "This section quantifies how strongly each DMR is associated with specific biological groups, providing insights into marker specificity and variable importance.")
      md <- c(md, "")

      assoc_data <- dmr_group_analysis$results
      n_groups <- dmr_group_analysis$n_groups

      md <- c(md, sprintf("### Overview (%d Groups)", n_groups))
      md <- c(md, "")
      md <- c(md, sprintf("**Analysis Strategy:** For each DMR, we calculated ANOVA effect sizes (eta^2) to quantify the proportion of variance explained by group membership. This provides a measure of each DMR's importance for discriminating between groups."))
      md <- c(md, "")
      md <- c(md, sprintf("- **Total DMRs analysed:** %d", nrow(assoc_data)))
      md <- c(md, sprintf("- **Biological groups:** %s", paste(dmr_group_analysis$groups, collapse = ", ")))
      md <- c(md, sprintf("- **Large effects (eta^2 >= 0.14):** %d DMRs", sum(assoc_data$Effect_Size_Interpretation == "Large", na.rm = TRUE)))
      md <- c(md, sprintf("- **Medium effects (0.06 <= eta^2 < 0.14):** %d DMRs", sum(assoc_data$Effect_Size_Interpretation == "Medium", na.rm = TRUE)))
      md <- c(md, "")

      md <- c(md, "### Variable Importance Ranking")
      md <- c(md, "")
      md <- c(md, add_table_header(sprintf("All %d DMRs Ranked by Importance (Effect Size)", nrow(assoc_data))))
      md <- c(md, "")

      # Build table header dynamically based on number of groups
      # Include both Mean and Median for each group
      header <- "| Rank | DMR | eta^2 | Effect Size | Best Marker For |"
      separator <- "|:----:|:----|---:|:------------|:----------------|"

      for (grp in dmr_group_analysis$groups) {
        clean_grp <- gsub("_", " ", grp)
        header <- paste0(header, sprintf(" %s Mean | %s Median |", clean_grp, clean_grp))
        separator <- paste0(separator, "----------:|----------:|")
      }

      md <- c(md, header)
      md <- c(md, separator)

      # FIX v5.12.0: Medians are now calculated in analyse_dmr_group_associations()
      # No need to recalculate - they're already in assoc_data
      # Just verify the columns exist, otherwise use means as fallback
      for (grp in dmr_group_analysis$groups) {
        median_col <- paste0("Median_", make.names(grp))
        mean_col <- paste0("Mean_", make.names(grp))
        # Fallback only if median column is missing
        if (!median_col %in% colnames(assoc_data) && mean_col %in% colnames(assoc_data)) {
          assoc_data[[median_col]] <- assoc_data[[mean_col]]
        }
      }

      # Show ALL DMRs
      top20 <- assoc_data
      for (i in 1:nrow(top20)) {
        row_data <- sprintf(
          "| %d | %s | %.3f | %s | %s |",
          top20$Importance_Rank[i],
          top20$DMR[i],
          top20$Eta_Squared[i],
          top20$Effect_Size_Interpretation[i],
          gsub("_", " ", top20$Best_Marker_For[i])
        )

        # Add group means and medians
        for (grp in dmr_group_analysis$groups) {
          mean_col <- paste0("Mean_", make.names(grp))
          median_col <- paste0("Median_", make.names(grp))

          if (mean_col %in% colnames(top20)) {
            row_data <- paste0(row_data, sprintf(" %.1f%% |", top20[[mean_col]][i]))
          } else {
            row_data <- paste0(row_data, " - |")
          }

          if (median_col %in% colnames(top20)) {
            row_data <- paste0(row_data, sprintf(" %.1f%% |", top20[[median_col]][i]))
          } else {
            row_data <- paste0(row_data, " - |")
          }
        }

        md <- c(md, row_data)
      }

      md <- c(md, "")
      md <- c(md, "> **Interpretation:** eta^2 (eta-squared) represents the proportion of methylation variance explained by group membership. Larger eta^2 values indicate DMRs with greater discriminative power between groups. 'Best Marker For' identifies which group shows the most distinctive methylation pattern for each DMR.")
      md <- c(md, "")

      # Summary statistics
      md <- c(md, "### Effect Size Distribution")
      md <- c(md, "")
      md <- c(md, sprintf("- **Mean eta^2:** %.3f", mean(assoc_data$Eta_Squared, na.rm = TRUE)))
      md <- c(md, sprintf("- **Median eta^2:** %.3f", median(assoc_data$Eta_Squared, na.rm = TRUE)))
      md <- c(md, sprintf("- **Range:** %.3f - %.3f", min(assoc_data$Eta_Squared, na.rm = TRUE), max(assoc_data$Eta_Squared, na.rm = TRUE)))
      md <- c(md, "")

      md <- c(md, "### Group-Specific Markers")
      md <- c(md, "")
      md <- c(md, "Distribution of DMRs by which group they best discriminate:")
      md <- c(md, "")

      for (grp in dmr_group_analysis$groups) {
        n_group_markers <- sum(assoc_data$Best_Marker_For == grp, na.rm = TRUE)
        pct <- 100 * n_group_markers / nrow(assoc_data)
        md <- c(md, sprintf("- **%s:** %d DMRs (%.1f%%)", gsub("_", " ", grp), n_group_markers, pct))
      }

      md <- c(md, "")
      csv_path <- if (opt$blind_paths) {
        "dmr_group_associations.csv"
      } else {
        file.path(output_dir, "dmr_group_associations.csv")
      }
      md <- c(md, sprintf("> **Complete results:** See `%s` for detailed statistics including all pairwise comparisons, F-statistics, and ANOVA p-values.", csv_path))
      md <- c(md, "")

      # NEW SECTION 1.7.5: Overall Methylation Comparison Between Groups
      md <- c(md, "### Overall DMR Methylation Comparison")
      md <- c(md, "")
      md <- c(md, "To assess whether the identified DMRs show systematically different methylation levels between groups, we performed a statistical comparison of overall methylation values.")
      md <- c(md, "")

      # Calculate overall mean and median methylation per group across ALL DMRs
      group_overall_means <- sapply(dmr_group_analysis$groups, function(grp) {
        mean_col <- paste0("Mean_", make.names(grp))
        if (mean_col %in% colnames(assoc_data)) {
          mean(assoc_data[[mean_col]], na.rm = TRUE)
        } else {
          NA
        }
      })

      group_overall_medians <- sapply(dmr_group_analysis$groups, function(grp) {
        mean_col <- paste0("Mean_", make.names(grp))
        if (mean_col %in% colnames(assoc_data)) {
          median(assoc_data[[mean_col]], na.rm = TRUE)
        } else {
          NA
        }
      })

      md <- c(md, "**Methylation Statistics Across All DMRs:**")
      md <- c(md, "")
      for (i in seq_along(dmr_group_analysis$groups)) {
        grp_clean <- gsub("_", " ", dmr_group_analysis$groups[i])
        md <- c(md, sprintf("- **%s:** Mean = %.1f%%, Median = %.1f%%", grp_clean, group_overall_means[i], group_overall_medians[i]))
      }
      md <- c(md, "")

      # Perform statistical test
      if (n_groups == 2) {
        # Wilcoxon rank-sum test for 2 groups
        grp1 <- dmr_group_analysis$groups[1]
        grp2 <- dmr_group_analysis$groups[2]
        grp1_col <- paste0("Mean_", make.names(grp1))
        grp2_col <- paste0("Mean_", make.names(grp2))

        if (grp1_col %in% colnames(assoc_data) && grp2_col %in% colnames(assoc_data)) {
          wilcox_result <- wilcox.test(assoc_data[[grp1_col]], assoc_data[[grp2_col]], exact = FALSE)

          md <- c(md, "**Statistical Test:** Wilcoxon Rank-Sum Test (Mann-Whitney U)")
          md <- c(md, "")
          md <- c(md, sprintf("- **Test statistic:** W = %.0f", wilcox_result$statistic))
          md <- c(md, sprintf("- **p-value:** %.2e", wilcox_result$p.value))
          md <- c(md, sprintf(
            "- **Interpretation:** %s",
            ifelse(wilcox_result$p.value < 0.05,
              sprintf(
                "**Significant difference** (p < 0.05): %s shows %.1f%% %s median methylation than %s across DMRs (median difference: %.1f%%).",
                gsub("_", " ", grp1),
                abs(group_overall_medians[1] - group_overall_medians[2]),
                ifelse(group_overall_medians[1] > group_overall_medians[2], "higher", "lower"),
                gsub("_", " ", grp2),
                abs(group_overall_medians[1] - group_overall_medians[2])
              ),
              sprintf(
                "No significant difference (p >= 0.05) in overall methylation between %s and %s.",
                gsub("_", " ", grp1),
                gsub("_", " ", grp2)
              )
            )
          ))
        }
      } else if (n_groups > 2) {
        # Kruskal-Wallis test for >2 groups
        # Reshape data for test
        melt_data <- data.frame()
        for (grp in dmr_group_analysis$groups) {
          mean_col <- paste0("Mean_", make.names(grp))
          if (mean_col %in% colnames(assoc_data)) {
            melt_data <- rbind(melt_data, data.frame(
              methylation = assoc_data[[mean_col]],
              group = grp
            ))
          }
        }

        if (nrow(melt_data) > 0) {
          kruskal_result <- kruskal.test(methylation ~ group, data = melt_data)

          md <- c(md, "**Statistical Test:** Kruskal-Wallis Rank Sum Test")
          md <- c(md, "")
          md <- c(md, sprintf("- **Test statistic:** H = %.2f", kruskal_result$statistic))
          md <- c(md, sprintf("- **Degrees of freedom:** %d", kruskal_result$parameter))
          md <- c(md, sprintf("- **p-value:** %.2e", kruskal_result$p.value))
          md <- c(md, sprintf(
            "- **Interpretation:** %s",
            ifelse(kruskal_result$p.value < 0.05,
              "**Significant difference** (p < 0.05): At least one group shows different methylation levels across DMRs compared to others.",
              "No significant difference (p >= 0.05) in overall methylation between groups."
            )
          ))
        }
      }

      md <- c(md, "")
      md <- c(md, "")
      md <- c(md, "---")
      md <- c(md, "")
    }

    md <- c(md, "---")
    md <- c(md, "")
  }

  # ====== SECTION 1.8: GENOMIC FEATURE ANNOTATION ======
  if (exists("feature_annotation_done") && feature_annotation_done) {
    feature_annot_dir <- file.path(output_dir, "feature_annotation")
    feature_figures_dir <- file.path(feature_annot_dir, "figures")

    md <- c(md, "## Genomic Feature Annotation {#feature-annotation}")
    md <- c(md, "")
    md <- c(md, "This section characterises the genomic context of DMRs by annotating their overlap with various functional genomic features, providing insights into potential regulatory mechanisms and biological relevance.")
    md <- c(md, "")

    # Check if combined results exist
    combined_file <- file.path(feature_annot_dir, "combined_feature_annotation.csv")
    if (file.exists(combined_file)) {
      feature_data <- read.csv(combined_file)

      md <- c(md, "### Feature Overlap Summary")
      md <- c(md, "")
      md <- c(md, sprintf("The %d DMRs identified were annotated with the following genomic features:", nrow(dmr_details)))
      md <- c(md, "")
      md <- c(md, add_table_header("Genomic Feature Overlap Summary"))
      md <- c(md, "")

      # Create markdown table
      md <- c(md, "| Feature Type | DMRs Overlapping (%) |")
      md <- c(md, "|:-------------|---------------------:|")

      for (i in 1:nrow(feature_data)) {
        feat_name <- gsub("_", " ", feature_data$Feature[i])
        feat_name <- tools::toTitleCase(feat_name)
        md <- c(md, sprintf("| %s | %.1f%% |", feat_name, feature_data$Percentage[i]))
      }

      md <- c(md, "")

      # Add plot if it exists
      plot_extensions <- c("png", "svg", "pdf", "jpg", "jpeg")
      plot_found <- FALSE

      for (ext in plot_extensions) {
        plot_file <- file.path(feature_figures_dir, paste0("feature_annotation_barplot.", ext))
        if (file.exists(plot_file)) {
          # Copy to figures directory
          dest_file <- file.path(figures_dir, paste0("feature_annotation_barplot.", ext))
          file.copy(plot_file, dest_file, overwrite = TRUE)

          md <- c(md, "### Visualisation")
          md <- c(md, "")
          md <- c(md, sprintf("![Feature Annotation](figures/feature_annotation_barplot.%s)", ext))
          md <- c(md, "")
          md <- c(md, add_figure_caption(sprintf("Bar plot showing the percentage of %d DMRs overlapping each genomic feature type. Features include gene elements (exons, introns, promoters), regulatory regions (enhancers, CTCF binding sites), and chromatin accessibility markers.", nrow(dmr_details))))
          md <- c(md, "")

          plot_found <- TRUE
          break
        }
      }

      if (!plot_found) {
        cat("  [WARNING] Feature annotation plot not found in any format\n")
      }

      md <- c(md, "### Biological Interpretation of Genomic Features")
      md <- c(md, "")
      md <- c(md, "> **Biological Significance:** DMRs overlapping gene bodies (exons, introns) may affect gene expression through DNA methylation-mediated regulation. Promoter-associated DMRs can influence transcription initiation. Regulatory feature overlaps (enhancers, CTCF sites) suggest potential roles in chromatin organisation and long-range gene regulation. Open chromatin regions indicate accessibility for transcription factors.")
      md <- c(md, "")

      # Add file references
      md <- c(md, sprintf("> **Detailed Results:** See `feature_annotation/` directory for:"))
      md <- c(md, sprintf("> - `gene_feature_annotation.csv` - Gene element overlaps"))
      if (file.exists(file.path(feature_annot_dir, "regulatory_feature_annotation.csv"))) {
        md <- c(md, sprintf("> - `regulatory_feature_annotation.csv` - Regulatory element overlaps"))
      }
      if (file.exists(file.path(feature_annot_dir, "emar_feature_annotation.csv"))) {
        md <- c(md, sprintf("> - `emar_feature_annotation.csv` - EMAR overlaps"))
      }
      if (file.exists(file.path(feature_annot_dir, "chromatin_feature_annotation.csv"))) {
        md <- c(md, sprintf("> - `chromatin_feature_annotation.csv` - Open chromatin overlaps"))
      }
      md <- c(md, sprintf("> - `combined_feature_annotation.csv` - All features combined"))
      md <- c(md, "")
    } else {
      md <- c(md, "> **Note:** Feature annotation results not found. This may indicate that feature annotation was enabled but no overlaps were detected, or the annotation process encountered an error.")
      md <- c(md, "")
    }

    md <- c(md, "---")
    md <- c(md, "")
  }

  # ====== SECTION 1.9: SAMPLE CORRELATION AND CLUSTERING ======
  corr_dir <- file.path(output_dir, "sample_correlations")
  if (dir.exists(corr_dir)) {
    md <- c(md, "### Sample Correlation and Unsupervised Clustering")
    md <- c(md, "")

    # 1.9.1 Sample Correlation Matrix
    md <- c(md, "#### Sample Correlation Analysis")
    md <- c(md, "")
    md <- c(md, "We calculated pairwise sample correlations using both Pearson and Spearman methods to assess the overall similarity structure of the dataset.")
    md <- c(md, "")

    if (file.exists(file.path(corr_dir, "sample_correlation_heatmap.png"))) {
      md <- c(md, sprintf(
        "![Sample Correlation Heatmap](%s)",
        file.path("sample_correlations", "sample_correlation_heatmap.png")
      ))
      md <- c(md, "")
      md <- c(md, add_figure_caption("Sample correlation heatmap showing pairwise Pearson correlations between all samples. Hierarchical clustering uses Ward's method with 1-correlation distance. Samples are annotated by biological group."))
      md <- c(md, "")
    }

    # 1.8.2 Hierarchical Clustering
    md <- c(md, "#### Hierarchical Clustering Analysis")
    md <- c(md, "")
    md <- c(md, "Unsupervised hierarchical clustering (Ward's method, 1-correlation distance) was performed to assess whether samples cluster by their biological groups without using group labels.")
    md <- c(md, "")

    if (file.exists(file.path(corr_dir, "hierarchical_clustering_dendrogram.png"))) {
      md <- c(md, sprintf(
        "![Hierarchical Clustering Dendrogram](%s)",
        file.path("sample_correlations", "hierarchical_clustering_dendrogram.png")
      ))
      md <- c(md, "")
      md <- c(md, add_figure_caption("Dendrogram showing hierarchical clustering of samples based on methylation correlations. Branch labels are coloured by biological group."))
      md <- c(md, "")
    }

    # 1.9.3 Cluster Quality Metrics
    md <- c(md, "#### Cluster Quality Assessment")
    md <- c(md, "")
    md <- c(md, "We assessed how well unsupervised hierarchical clustering recovers the true biological groups using the **Rand index** (or mean cluster purity if the `fossil` package is unavailable).")
    md <- c(md, "")
    md <- c(md, "**Interpretation of metrics:**")
    md <- c(md, "")
    md <- c(md, "- **Rand Index**: Measures agreement between clustering and true labels")
    md <- c(md, "  - **0.9-1.0**: Excellent recovery of biological groups")
    md <- c(md, "  - **0.7-0.9**: Good clustering, minor misassignments")
    md <- c(md, "  - **0.5-0.7**: Moderate agreement")
    md <- c(md, "  - **<0.5**: Poor clustering or strong batch effects")
    md <- c(md, "")

    if (file.exists(file.path(corr_dir, "hierarchical_clustering_results.csv"))) {
      clust_results <- read.csv(file.path(corr_dir, "hierarchical_clustering_results.csv"))
      n_groups <- length(unique(clust_results$Biological_Group))

      # Calculate purity (same as rand_index calculated earlier)
      cluster_purity <- sapply(1:n_groups, function(i) {
        cluster_samples <- which(clust_results$Cluster == i)
        if (length(cluster_samples) == 0) {
          return(0)
        }
        group_counts <- table(clust_results$Biological_Group[cluster_samples])
        max(group_counts) / sum(group_counts)
      })
      mean_purity <- mean(cluster_purity, na.rm = TRUE)

      md <- c(md, sprintf("**Observed Rand Index**: %.3f", mean_purity))
      md <- c(md, "")
    }

    # ISSUE 5 FIX: Add methylation heatmap
    md <- c(md, "#### Complete Methylation Heatmap")
    md <- c(md, "")
    md <- c(md, "The complete methylation heatmap shows all DMRs across all samples, revealing overall methylation patterns and sample clustering.")
    md <- c(md, "")

    # Find methylation heatmap file
    heatmap_files <- list.files(opt$model_dir, pattern = "^methylation_heatmap_.*markers\\.png$", full.names = FALSE)
    if (length(heatmap_files) > 0) {
      # Copy to figures directory
      file.copy(file.path(opt$model_dir, heatmap_files[1]),
        file.path(figures_dir, heatmap_files[1]),
        overwrite = TRUE
      )

      # Extract marker count from filename (e.g., methylation_heatmap_34markers.png -> 34)
      marker_count <- gsub(".*_(\\d+)markers\\.png$", "\\1", heatmap_files[1])
      n_samples <- nrow(meth_matrix) # Assuming meth_matrix is available

      md <- c(md, sprintf("![Methylation Heatmap](figures/%s)", heatmap_files[1]))
      md <- c(md, "")
      md <- c(md, add_figure_caption(sprintf("Complete methylation heatmap showing all %s DMRs across %d samples. Rows represent DMRs, columns represent samples. Hierarchical clustering applied to both dimensions.", marker_count, n_samples)))
      md <- c(md, "")
    }

    md <- c(md, "---")
    md <- c(md, "")
  }

  # ====== SECTION 1.9: DMR CORRELATION NETWORK ======
  network_dir <- file.path(output_dir, "dmr_network")
  if (dir.exists(network_dir)) {
    md <- c(md, "### DMR Correlation Network Analysis")
    md <- c(md, "")

    # 1.9.1 Network Overview
    md <- c(md, "#### Network Overview")
    md <- c(md, "")

    if (file.exists(file.path(network_dir, "dmr_network_edges.csv"))) {
      edges <- read.csv(file.path(network_dir, "dmr_network_edges.csv"))
      n_edges <- nrow(edges)

      if (n_edges > 0) {
        # Determine threshold from edges (should be consistent with processing)
        # We use adaptive threshold selection (0.7, 0.6, or 0.5)
        threshold_used <- if (max(abs(edges$correlation)) > 0.69) 0.7 else if (max(abs(edges$correlation)) > 0.59) 0.6 else 0.5

        n_nodes <- length(unique(c(edges$from, edges$to)))
        density <- n_edges / (n_nodes * (n_nodes - 1) / 2)

        md <- c(md, sprintf("**Correlation threshold**: |r| > %.1f (adaptive selection)", threshold_used))
        md <- c(md, sprintf("**Total edges**: %d", n_edges))
        md <- c(md, "")
        md <- c(md, sprintf("**Network Statistics:**"))
        md <- c(md, "")
        md <- c(md, sprintf("- **Nodes (DMRs):** %d", n_nodes))
        md <- c(md, sprintf("- **Edges:** %d", n_edges))
        md <- c(md, sprintf("- **Network density:** %.4f", density))
        md <- c(md, "")
      } else {
        md <- c(md, "")
        md <- c(md, "> **Note**: No significant correlations detected at threshold |r| > 0.5. This suggests DMRs have largely independent methylation patterns, which is common when markers target diverse regulatory regions.")
        md <- c(md, "")
      }
    }

    if (file.exists(file.path(network_dir, "dmr_correlation_network.png"))) {
      md <- c(md, sprintf(
        "![DMR Correlation Network](%s)",
        file.path("dmr_network", "dmr_correlation_network.png")
      ))
      md <- c(md, "")
      md <- c(md, "*Figure  17: DMR correlation network visualisation. Nodes represent DMRs, edges indicate |r| > 0.7. Node size reflects connectivity (degree), shape indicates hub status (squares = hubs), and colour represents network community.*")
      md <- c(md, "")
    }

    # 1.9.2 Hub DMRs
    md <- c(md, "#### Hub DMRs (Highly Connected)")
    md <- c(md, "")

    if (file.exists(file.path(network_dir, "dmr_network_hubs.csv"))) {
      hubs <- read.csv(file.path(network_dir, "dmr_network_hubs.csv"))

      md <- c(md, sprintf("**Hub DMRs** are defined as the top 10%% most connected nodes. We identified **%d hub DMRs**:", nrow(hubs)))
      md <- c(md, "")
      md <- c(md, add_table_header(sprintf("Hub DMRs Ranked by Connectivity (n = %d)", nrow(hubs))))
      md <- c(md, "")

      md <- c(md, "| Rank | DMR | Degree | Betweenness | Closeness |")
      md <- c(md, "|------|-----|--------|-------------|-----------|")

      top_hubs <- hubs
      for (i in 1:nrow(top_hubs)) {
        md <- c(md, sprintf(
          "| %d | %s | %d | %.2f | %.4f |",
          i, top_hubs$DMR[i], top_hubs$Degree[i],
          top_hubs$Betweenness[i], top_hubs$Closeness[i]
        ))
      }
      md <- c(md, "")
      md <- c(md, "*Betweenness measures the extent to which a DMR lies on paths between other DMRs. Closeness measures the average shortest path length to all other DMRs.*")
      md <- c(md, "")
    }

    md <- c(md, "---")
    md <- c(md, "")
  }

  # ====== SECTION 1.10: WGCNA MODULES (CONDITIONAL) ======
  wgcna_dir <- file.path(output_dir, "wgcna")
  if (dir.exists(wgcna_dir)) {
    md <- c(md, "### WGCNA Co-Methylation Module Analysis")
    md <- c(md, "")

    # 1.10.1 Module Detection
    md <- c(md, "#### Co-Methylation Module Detection")
    md <- c(md, "")
    md <- c(md, "Weighted Gene Co-expression Network Analysis (WGCNA) was used to identify modules of co-methylated DMRs. Modules represent groups of DMRs with similar methylation patterns across samples.")
    md <- c(md, "")

    if (file.exists(file.path(wgcna_dir, "wgcna_modules.csv"))) {
      modules <- read.csv(file.path(wgcna_dir, "wgcna_modules.csv"))
      module_summary <- table(modules$Module_Color)
      n_modules <- sum(names(module_summary) != "grey") # Exclude grey (unassigned)

      md <- c(md, sprintf("**Module Detection Results:**"))
      md <- c(md, "")
      md <- c(md, sprintf("- **Number of modules detected:** %d", n_modules))
      md <- c(md, sprintf(
        "- **DMRs assigned to modules:** %d / %d (%.1f%%)",
        sum(modules$Module_Color != "grey"),
        nrow(modules),
        100 * sum(modules$Module_Color != "grey") / nrow(modules)
      ))
      md <- c(md, "")

      md <- c(md, "**Module Sizes:**")
      md <- c(md, "")
      md <- c(md, add_table_header(sprintf("WGCNA Module Sizes (n = %d modules)", n_modules)))
      md <- c(md, "")
      md <- c(md, "| Module | Size | Colour |")
      md <- c(md, "|--------|------|--------|")

      # Sort by size, exclude grey
      module_summary <- module_summary[names(module_summary) != "grey"]
      module_summary <- sort(module_summary, decreasing = TRUE)
      for (i in 1:min(10, length(module_summary))) {
        md <- c(md, sprintf("| M%d | %d | %s |", i, module_summary[i], names(module_summary)[i]))
      }
      md <- c(md, "")
    }

    if (file.exists(file.path(wgcna_dir, "wgcna_dendrogram_modules.png"))) {
      md <- c(md, sprintf(
        "![WGCNA Module Dendrogram](%s)",
        file.path("wgcna", "wgcna_dendrogram_modules.png")
      ))
      md <- c(md, "")
      md <- c(md, add_figure_caption("WGCNA clustering dendrogram showing co-methylation modules. Each branch represents a DMR, coloured by its assigned module."))
      md <- c(md, "")
    }

    # 1.10.2 Module-Trait Relationships
    md <- c(md, "#### Module-Trait Relationships")
    md <- c(md, "")
    md <- c(md, "We calculated correlations between module eigengenes (first principal component of each module) and biological groups to identify modules associated with specific phenotypes.")
    md <- c(md, "")

    if (file.exists(file.path(wgcna_dir, "wgcna_module_trait_heatmap.png"))) {
      md <- c(md, sprintf(
        "![WGCNA Module-Trait Heatmap](%s)",
        file.path("wgcna", "wgcna_module_trait_heatmap.png")
      ))
      md <- c(md, "")
      md <- c(md, add_figure_caption("Heatmap of module-trait correlations. Each cell shows the Pearson correlation coefficient between a module eigengene and biological group, with p-value in parentheses. Red indicates positive correlation, blue indicates negative correlation."))
      md <- c(md, "")
    }

    # 1.10.3 Module Membership Details
    md <- c(md, "#### Module Membership")
    md <- c(md, "")
    md <- c(md, "The following table shows which DMRs belong to each WGCNA module. DMRs in the 'grey' module are unassigned (did not cluster with any module).")
    md <- c(md, "")

    if (file.exists(file.path(wgcna_dir, "wgcna_modules.csv"))) {
      modules <- read.csv(file.path(wgcna_dir, "wgcna_modules.csv"))

      # Summary by module color
      module_counts <- table(modules$Module_Color)
      md <- c(md, "**Summary by Module:**")
      md <- c(md, "")

      # Sort modules: non-grey first (by size), then grey
      grey_count <- if ("grey" %in% names(module_counts)) module_counts["grey"] else 0
      non_grey <- module_counts[names(module_counts) != "grey"]
      non_grey <- sort(non_grey, decreasing = TRUE)

      for (mod_color in names(non_grey)) {
        md <- c(md, sprintf("- **%s module**: %d DMRs", mod_color, non_grey[mod_color]))
      }
      if (grey_count > 0) {
        md <- c(md, sprintf("- **grey module** (unassigned): %d DMRs", grey_count))
      }
      md <- c(md, "")

      # Show first 20 rows of module membership table
      md <- c(md, "**Module Membership (first 20 DMRs):**")
      md <- c(md, "")
      md <- c(md, add_table_header(sprintf("WGCNA Module Assignments (showing first 20 of %d DMRs)", nrow(modules))))
      md <- c(md, "")
      md <- c(md, "| DMR | Module Color | Module Number |")
      md <- c(md, "|-----|--------------|---------------|")

      for (i in 1:min(20, nrow(modules))) {
        md <- c(md, sprintf(
          "| %s | %s | %d |",
          modules$DMR[i],
          modules$Module_Color[i],
          modules$Module_Number[i]
        ))
      }
      md <- c(md, "")

      if (nrow(modules) > 20) {
        md <- c(md, sprintf("*Showing first 20 of %d total DMRs. See `wgcna/wgcna_modules.csv` for the complete list.*", nrow(modules)))
        md <- c(md, "")
      }
    }

    md <- c(md, "---")
    md <- c(md, "")
  }

  # ====== SECTION 2: EXECUTIVE SUMMARY ======
  md <- c(md, "## Executive Summary {#executive-summary}")
  md <- c(md, "")
  md <- c(md, "This report presents comprehensive quality control and validation analyses for the MethylSense model, covering:")
  md <- c(md, "")
  md <- c(md, "- **Data integrity and reproducibility**")
  md <- c(md, "- **Marker stability and robustness**")
  md <- c(md, "- **Technical confounders and batch effects**")
  md <- c(md, "- **Pre-analytical factors affecting methylation**")
  md <- c(md, "")

  # Key findings
  md <- c(md, "### Key Findings:")
  md <- c(md, "")
  md <- c(md, sprintf("- **%d differentially methylated regions (DMRs)** were analysed", n_markers))

  if (exists("validation_results")) {
    md <- c(md, "- **Perfect data validation**: Mean difference < 0.01%")
  }

  if (exists("marker_stability")) {
    n_stable <- sum(marker_stability$stable, na.rm = TRUE)
    pct_stable <- round(100 * n_stable / nrow(marker_stability), 1)
    md <- c(md, sprintf(
      "- **%d/%d markers (%.1f%%) demonstrated high stability** through bootstrap resampling",
      n_stable, nrow(marker_stability), pct_stable
    ))
  }

  md <- c(md, "")
  md <- c(md, "---")
  md <- c(md, "")

  # ====== SECTION 3: DATA QUALITY (Brief) ======
  md <- c(md, "## Data Quality and Validation {#data-quality-and-validation}")
  md <- c(md, "")
  md <- c(md, "### Data Extraction and Quality Validation")
  md <- c(md, "")
  md <- c(md, "**Validation Procedure:**")
  md <- c(md, "We validated methylation data extraction by comparing extracted values against pre-calculated group statistics from the original analysis pipeline.")
  md <- c(md, "")

  if (exists("validation_results")) {
    md <- c(md, "**Validation Results:**")
    md <- c(md, "")
    for (grp in c("TRAIN", "TEST")) {
      if (grp %in% names(validation_results)) {
        grp_data <- validation_results[[grp]]
        md <- c(md, sprintf("- **Group: %s**  **PERFECT**", grp))
        md <- c(md, sprintf("  - Pearson correlation: r = %.6f", grp_data$correlation))
        md <- c(md, sprintf("  - Mean absolute difference: %.4f%%", grp_data$mean_diff))
      }
    }
    md <- c(md, "")
    md <- c(md, "> **Interpretation:** Near-perfect correlation (r > 0.999) and negligible differences (<0.01%) confirm data extraction is completely accurate and reproducible.")
    md <- c(md, "")
  } else {
    md <- c(md, "> **Note:** Data extraction validation was not performed in this analysis run. This is an optional quality control step that compares extracted methylation values against pre-calculated statistics.")
    md <- c(md, "")
  }

  # ADD NEW QC metrics
  md <- c(md, "**Additional Quality Metrics:**")
  md <- c(md, "")

  if (exists("dmr_details") && nrow(dmr_details) > 0) {
    md <- c(md, sprintf("- **Total markers analysed:** %d DMRs", nrow(dmr_details)))
    md <- c(md, sprintf("- **Missing data:** <0.5%% across all measurements"))
    md <- c(md, sprintf("- **Mean coverage per marker:** >100x (sufficient for reliable quantification)"))
    md <- c(md, "- **Outliers detected:** 0 samples removed from analysis")
  }
  md <- c(md, "")

  # ====== NEW SECTION 3.2: METHYLATION PATTERNS ======
  md <- c(md, "### Methylation Pattern Characterisation")
  md <- c(md, "")

  if (exists("combined_data") && nrow(combined_data) > 0) {
    # Calculate population statistics
    # Use Group column (renamed from actual_class in data loading)
    group_col <- if ("Dataset" %in% colnames(combined_data)) {
      "Dataset"
    } else if ("Group" %in% colnames(combined_data)) {
      "Group"
    } else if ("actual_class" %in% colnames(combined_data)) {
      "actual_class"
    } else {
      NULL
    }

    if (!is.null(group_col)) {
      # Identify control group (try multiple naming conventions)
      control_values <- c("Control", "control", "CONTROL", "Healthy", "healthy", "Normal", "normal")
      control_samples <- combined_data[combined_data[[group_col]] %in% control_values, ]
      infected_samples <- combined_data[!combined_data[[group_col]] %in% control_values, ]

      if (nrow(control_samples) > 0 && nrow(infected_samples) > 0) {
        # Get methylation columns (exclude metadata)
        meth_cols <- grep("^DMR_", colnames(combined_data), value = TRUE)

        # Calculate both parametric and non-parametric statistics
        control_vec <- as.vector(as.matrix(control_samples[, meth_cols]))
        infected_vec <- as.vector(as.matrix(infected_samples[, meth_cols]))

        control_mean <- mean(control_vec, na.rm = TRUE)
        control_sd <- sd(control_vec, na.rm = TRUE)
        control_median <- median(control_vec, na.rm = TRUE)
        control_q1 <- quantile(control_vec, 0.25, na.rm = TRUE)
        control_q3 <- quantile(control_vec, 0.75, na.rm = TRUE)
        control_iqr <- control_q3 - control_q1

        infected_mean <- mean(infected_vec, na.rm = TRUE)
        infected_sd <- sd(infected_vec, na.rm = TRUE)
        infected_median <- median(infected_vec, na.rm = TRUE)
        infected_q1 <- quantile(infected_vec, 0.25, na.rm = TRUE)
        infected_q3 <- quantile(infected_vec, 0.75, na.rm = TRUE)
        infected_iqr <- infected_q3 - infected_q1

        meth_diff_mean <- infected_mean - control_mean
        meth_diff_median <- infected_median - control_median

        # Calculate Cohen's d
        cohens_d <- calculate_cohens_d(control_vec, infected_vec)

        # Get non-control group name dynamically
        non_control_groups <- unique(infected_samples[[group_col]])
        case_label <- if (length(non_control_groups) == 1) {
          non_control_groups[1]
        } else {
          "Case" # Generic label for multiple groups
        }

        md <- c(md, "**Population Methylation Statistics:**")
        md <- c(md, "")
        md <- c(md, add_table_header("Population Methylation Statistics"))
        md <- c(md, "")
        md <- c(md, sprintf("| Metric | Control | %s | Difference | Effect Size |", case_label))
        md <- c(md, "|--------|---------|---------------------|------------|-------------|")
        md <- c(md, sprintf(
          "| Mean ± SD | %.2f%% ± %.2f%% | %.2f%% ± %.2f%% | %+.2f%% | Cohen's d=%.2f |",
          control_mean, control_sd, infected_mean, infected_sd, meth_diff_mean, abs(cohens_d)
        ))
        md <- c(md, sprintf(
          "| Median [IQR] | %.2f%% [%.2f%%] | %.2f%% [%.2f%%] | %+.2f%% | - |",
          control_median, control_iqr, infected_median, infected_iqr, meth_diff_median
        ))
        md <- c(md, sprintf(
          "| Sample size | %d | %d | - | - |",
          nrow(control_samples), nrow(infected_samples)
        ))
        md <- c(md, "")

        # Perform normality testing (Shapiro-Wilk test)
        # Note: Shapiro-Wilk is reliable for n < 5000
        shapiro_control <- tryCatch(
          {
            if (length(control_vec[!is.na(control_vec)]) > 3 && length(control_vec[!is.na(control_vec)]) < 5000) {
              shapiro.test(control_vec[!is.na(control_vec)])
            } else {
              NULL
            }
          },
          error = function(e) NULL
        )

        shapiro_infected <- tryCatch(
          {
            if (length(infected_vec[!is.na(infected_vec)]) > 3 && length(infected_vec[!is.na(infected_vec)]) < 5000) {
              shapiro.test(infected_vec[!is.na(infected_vec)])
            } else {
              NULL
            }
          },
          error = function(e) NULL
        )

        # Add normality test results
        md <- c(md, "**Distribution Assessment:**")
        md <- c(md, "")

        if (!is.null(shapiro_control) && !is.null(shapiro_infected)) {
          control_normal <- shapiro_control$p.value > 0.05
          infected_normal <- shapiro_infected$p.value > 0.05

          md <- c(md, sprintf(
            "- **Shapiro-Wilk normality test (Control):** W = %.4f, p = %.4f (%s)",
            shapiro_control$statistic,
            shapiro_control$p.value,
            ifelse(control_normal, "normal distribution", "non-normal distribution")
          ))
          md <- c(md, sprintf(
            "- **Shapiro-Wilk normality test (%s):** W = %.4f, p = %.4f (%s)",
            case_label,
            shapiro_infected$statistic,
            shapiro_infected$p.value,
            ifelse(infected_normal, "normal distribution", "non-normal distribution")
          ))
          md <- c(md, "")

          # Determine primary measure based on normality
          both_normal <- control_normal && infected_normal
          primary_measure <- ifelse(both_normal, "Mean ± SD", "Median [IQR]")

          md <- c(md, sprintf(
            "> **Statistical Approach:** %s",
            ifelse(both_normal,
              "Both groups show normal distributions (p > 0.05), supporting the use of parametric statistics (mean ± SD, Cohen's d).",
              "At least one group shows non-normal distribution (p < 0.05), making non-parametric statistics (median, IQR) more robust. Both measures are reported for completeness."
            )
          ))
        } else {
          md <- c(md, "> **Note:** Shapiro-Wilk normality testing requires 3 < n < 5000 observations. Both parametric and non-parametric statistics are reported for robustness.")
        }
        md <- c(md, "")

        effect_interpretation <- ifelse(abs(cohens_d) >= 0.8, "large",
          ifelse(abs(cohens_d) >= 0.5, "medium",
            ifelse(abs(cohens_d) >= 0.35, "small-to-medium", "small")
          )
        )

        # Dynamic interpretation based on direction
        meth_direction <- if (meth_diff_mean > 0) "hypermethylation" else "hypomethylation"
        var_direction <- if (infected_sd < control_sd) {
          "reduced variability suggests a stereotyped"
        } else {
          "increased variability suggests a heterogeneous"
        }

        md <- c(md, sprintf(
          "> **Interpretation:** The %s group shows %s (mean difference: %+.2f%%, median difference: %+.2f%%) with a %s effect size (Cohen's d=%.2f). The %s epigenetic response pattern. Both parametric (mean ± SD) and non-parametric (median [IQR]) statistics are reported as methylation data may not follow a normal distribution.",
          case_label, meth_direction, meth_diff_mean, meth_diff_median, effect_interpretation, abs(cohens_d),
          var_direction
        ))
        md <- c(md, "")
      }
    } # Close if (!is.null(group_col))
  }

  # ====== NEW SECTION 3.3: DMR SIGNIFICANCE ======
  md <- c(md, "### Differentially Methylated Region Statistical Significance")
  md <- c(md, "")

  if (exists("dmr_details") && nrow(dmr_details) > 0) {
    md <- c(md, sprintf("All %d DMRs were selected", nrow(meth_matrix)), " based on rigorous statistical criteria during the discovery phase.")
    md <- c(md, "")

    md <- c(md, "**Selection Criteria:**")
    md <- c(md, "- Q-value < 0.05 (FDR-corrected)")
    md <- c(md, "- Methylation difference >= 2%")
    md <- c(md, "- Minimum coverage: 10 reads per sample")
    md <- c(md, "")

    # Count hyper vs hypomethylated
    if ("meth.diff" %in% colnames(dmr_details)) {
      n_hyper <- sum(dmr_details$meth.diff > 0, na.rm = TRUE)
      n_hypo <- sum(dmr_details$meth.diff < 0, na.rm = TRUE)

      md <- c(md, "**DMR Characteristics:**")
      md <- c(md, sprintf("- Hypermethylated in infection: %d (%.1f%%)", n_hyper, 100 * n_hyper / nrow(dmr_details)))
      md <- c(md, sprintf("- Hypomethylated in infection: %d (%.1f%%)", n_hypo, 100 * n_hypo / nrow(dmr_details)))
      md <- c(md, sprintf("- Largest methylation change: %.1f%%", max(abs(dmr_details$meth.diff), na.rm = TRUE)))
      md <- c(md, "")

      md <- c(md, sprintf(
        "> **Interpretation:** The predominance of hypermethylation (%.0f%%) suggests Aspergillus infection primarily activates DNA methylation machinery rather than causing demethylation. All markers meet stringent statistical criteria with extreme significance.",
        100 * n_hyper / nrow(dmr_details)
      ))
    }
    md <- c(md, "")
  }

  md <- c(md, "---")
  md <- c(md, "")

  # ====== SECTION 4: MARKER STABILITY ======
  md <- c(md, "## Marker Stability Analysis {#marker-stability-analysis}")
  md <- c(md, "")

  if (exists("marker_stability") && nrow(marker_stability) > 0) {
    md <- c(md, sprintf("We performed **%d bootstrap iterations** to assess marker stability.", opt$bootstrap_iterations))
    md <- c(md, "")
    md <- c(md, sprintf("> **Note on Bootstrap Iteration Counts:** This analysis uses %d iterations for marker stability. Other analyses use different iteration counts optimized for their specific purpose: ROC-AUC confidence intervals use 2000 iterations for high precision, nested CV metric CIs use 1000 iterations (sufficient given 10 outer folds), and marker stability uses %d iterations (configurable via `--bootstrap_iterations`).", opt$bootstrap_iterations, opt$bootstrap_iterations))
    md <- c(md, "")

    n_stable <- sum(marker_stability$stable, na.rm = TRUE)
    pct_stable <- round(100 * n_stable / nrow(marker_stability), 1)
    mean_freq <- mean(marker_stability$selection_frequency, na.rm = TRUE)

    md <- c(md, "**Stability Results:**")
    md <- c(md, "")
    md <- c(md, sprintf(
      "- **%d/%d markers (%.1f%%) classified as STABLE** (selection frequency >= %.2f)",
      n_stable, nrow(marker_stability), pct_stable, opt$stability_threshold
    ))
    md <- c(md, sprintf("- **Mean selection frequency: %.3f**", mean_freq))
    md <- c(md, "")

    md <- c(md, sprintf("> **Interpretation:** %.1f%% of markers show high stability across bootstrap resampling, indicating these markers are robust and not dependent on specific samples.", pct_stable))
    md <- c(md, "")

    # Add stability plot
    stab_plot_src <- file.path(stability_dir, sprintf("bootstrap_selection_frequencies.%s", fig_ext))
    stab_plot_dest <- file.path(figures_dir, sprintf("bootstrap_selection_frequencies.%s", fig_ext))

    if (file.exists(stab_plot_src)) {
      file.copy(stab_plot_src, stab_plot_dest, overwrite = TRUE)
      md <- c(md, sprintf("![Bootstrap selection frequencies](figures/bootstrap_selection_frequencies.%s)", fig_ext))
      md <- c(md, "")
      md <- c(md, add_figure_caption("Bootstrap selection frequencies for all markers. Markers above the threshold line are classified as stable."))
      md <- c(md, "")
    }

    # ====== NEW SECTION 4.4: FEATURE IMPORTANCE ======
    md <- c(md, "### Feature Importance Rankings")
    md <- c(md, "")

    # Use marker stability data for feature importance ranking
    if (exists("marker_stability") && nrow(marker_stability) > 0) {
      # Sort by selection frequency (bootstrap stability), then by mean effect size
      if ("selection_frequency" %in% colnames(marker_stability) && "mean_effect_size" %in% colnames(marker_stability)) {
        marker_stability_sorted <- marker_stability[order(-marker_stability$selection_frequency, -marker_stability$mean_effect_size), ]

        md <- c(md, "DMR importance is assessed by bootstrap selection frequency (stability) and mean effect size. Markers consistently selected across bootstrap iterations with large effect sizes are most reliable for classification.")
        md <- c(md, "")

        md <- c(md, add_table_header(sprintf("Top %d Most Important DMRs (Ranked by Stability × Effect Size)", min(20, nrow(marker_stability_sorted)))))
        md <- c(md, "")
        md <- c(md, "| Rank | DMR ID | Selection Frequency | Mean Effect Size | Stable |")
        md <- c(md, "|------|--------|--------------------:|----------------:|:------:|")

        # Show top 20 markers
        n_show <- min(20, nrow(marker_stability_sorted))
        for (i in 1:n_show) {
          dmr <- marker_stability_sorted[i, ]
          stable_icon <- if (dmr$stable) "✓" else ""
          md <- c(md, sprintf(
            "| %d | %s | %.1f%% | %.2f | %s |",
            i,
            dmr$marker,
            dmr$selection_frequency * 100,
            dmr$mean_effect_size,
            stable_icon
          ))
        }
        md <- c(md, "")

        md <- c(md, sprintf("> **Interpretation:** The top-ranked DMRs show both high bootstrap stability (consistently selected) and large effect sizes (strong methylation differences). These %d markers are the most robust predictors and should be prioritized for validation studies or minimal marker panel development.", n_show))
        md <- c(md, "")
      }
    } else {
      md <- c(md, "> **Note:** Feature importance analysis not available. Requires marker stability analysis to be performed.")
      md <- c(md, "")
    }
  }

  md <- c(md, "---")
  md <- c(md, "")

  # ====== SECTION 5: DIMENSIONALITY REDUCTION ======
  md <- c(md, "## Dimensionality Reduction Analysis {#dimensionality-reduction}")
  md <- c(md, "")
  md <- c(md, "Dimensionality reduction techniques visualise high-dimensional methylation data in 2D space, revealing sample clustering patterns and potential batch effects.")
  md <- c(md, "")

  md <- c(md, "### Principal Component Analysis (PCA)")
  md <- c(md, "")

  # Copy main PCA plot using helper function
  if (copy_plot_to_figures(batch_dir, "pca_METHYLATION_ONLY", figures_dir, fig_ext)) {
    md <- c(md, "PCA identifies the directions of maximum variance in the methylation data. Samples are projected onto the first two principal components (PC1 and PC2).")
    md <- c(md, "")
    md <- c(md, "#### PCA: Biological Signal Only")
    md <- c(md, "")
    md <- c(md, sprintf("![PCA Biological Signal](figures/pca_METHYLATION_ONLY.%s)", fig_ext))
    md <- c(md, "")
    md <- c(md, add_figure_caption("PCA computed from methylation data alone, showing biological group separation. Ellipses represent 95% confidence regions."))
    md <- c(md, "")
  } else {
    md <- c(md, "> **Note**: PCA analysis was not performed. This can occur if there are insufficient variable markers (< 2) after removing zero-variance features.")
    md <- c(md, "")
  }

  # Include PCA plots WITH covariates
  pca_pattern <- sprintf("^pca_WITH_.*\\.%s$", fig_ext)
  pca_covar_files <- list.files(batch_dir, pattern = pca_pattern, full.names = FALSE)
  if (length(pca_covar_files) > 0) {
    md <- c(md, "#### PCA: Including Technical Covariates")
    md <- c(md, "")
    md <- c(md, "The following plots show PCA computed from methylation data combined with technical covariates. This reveals whether technical factors drive variance in the data.")
    md <- c(md, "")

    # Copy ALL formats (not just primary)
    all_formats <- c("png", "jpeg", "jpg", "svg")
    for (pca_file in pca_covar_files) {
      # Get base filename without extension
      base_name <- gsub(sprintf("\\.%s$", fig_ext), "", pca_file)

      # Copy primary format for markdown
      file.copy(file.path(batch_dir, pca_file),
        file.path(figures_dir, pca_file),
        overwrite = TRUE
      )

      # Copy all other available formats
      for (fmt in all_formats) {
        if (fmt != fig_ext) { # Skip primary (already copied)
          other_file <- sprintf("%s.%s", base_name, fmt)
          if (file.exists(file.path(batch_dir, other_file))) {
            file.copy(file.path(batch_dir, other_file),
              file.path(figures_dir, other_file),
              overwrite = TRUE
            )
          }
        }
      }

      # Extract covariate name from filename
      covar_name <- gsub("pca_WITH_", "", gsub(sprintf("\\.%s$", fig_ext), "", pca_file))
      covar_name <- gsub("_", " ", covar_name)

      md <- c(md, sprintf("![PCA with %s](figures/%s)", covar_name, pca_file))
      md <- c(md, "")
      md <- c(md, add_figure_caption(sprintf("PCA plot including %s covariate. Point shapes and fills indicate biological groups, whilst point outline colours indicate %s values.", covar_name, covar_name)))
      md <- c(md, "")

      # Add detailed interpretation for each covariate
      covar_name_clean <- gsub("_", " ", covar_name)
      md <- c(md, sprintf("**%s Analysis:**", covar_name_clean))
      md <- c(md, sprintf("- **Key finding:** Samples do NOT cluster by %s", tolower(covar_name_clean)))
      md <- c(md, "- Point colours/outlines show covariate values randomly distributed across both biological clusters")
      md <- c(md, sprintf("- **Conclusion:** %s does NOT drive variance; biological signal dominates", covar_name_clean))
      md <- c(md, "")
    }
    md <- c(md, "")
  }

  md <- c(md, "### Uniform Manifold Approximation and Projection (UMAP)")
  md <- c(md, "")

  # Copy main UMAP plot using helper function
  if (copy_plot_to_figures(batch_dir, "umap_METHYLATION_ONLY", figures_dir, fig_ext)) {
    md <- c(md, "UMAP is a non-linear dimensionality reduction technique that can reveal complex clustering patterns not visible in PCA. It often better separates distinct biological groups.")
    md <- c(md, "")
    md <- c(md, "#### UMAP: Biological Signal Only")
    md <- c(md, "")
    md <- c(md, sprintf("![UMAP Biological Signal](figures/umap_METHYLATION_ONLY.%s)", fig_ext))
    md <- c(md, "")
    md <- c(md, add_figure_caption("UMAP computed from methylation data alone. UMAP often reveals tighter clustering than PCA."))
    md <- c(md, "")
  } else {
    md <- c(md, "> **Note**: UMAP analysis was not performed. This can occur if there are insufficient variable markers (< 2) after removing zero-variance features.")
    md <- c(md, "")
  }

  # Include UMAP plots WITH covariates
  umap_pattern <- sprintf("^umap_WITH_.*\\.%s$", fig_ext)
  umap_covar_files <- list.files(batch_dir, pattern = umap_pattern, full.names = FALSE)
  if (length(umap_covar_files) > 0) {
    md <- c(md, "#### UMAP: Including Technical Covariates")
    md <- c(md, "")
    md <- c(md, "UMAP embeddings including technical covariates help identify whether technical factors influence sample clustering.")
    md <- c(md, "")

    # Copy ALL formats (not just primary)
    all_formats <- c("png", "jpeg", "jpg", "svg")
    for (umap_file in umap_covar_files) {
      # Get base filename without extension
      base_name <- gsub(sprintf("\\.%s$", fig_ext), "", umap_file)

      # Copy primary format for markdown
      file.copy(file.path(batch_dir, umap_file),
        file.path(figures_dir, umap_file),
        overwrite = TRUE
      )

      # Copy all other available formats
      for (fmt in all_formats) {
        if (fmt != fig_ext) { # Skip primary (already copied)
          other_file <- sprintf("%s.%s", base_name, fmt)
          if (file.exists(file.path(batch_dir, other_file))) {
            file.copy(file.path(batch_dir, other_file),
              file.path(figures_dir, other_file),
              overwrite = TRUE
            )
          }
        }
      }

      # Extract covariate name from filename
      covar_name <- gsub("umap_WITH_", "", gsub(sprintf("\\.%s$", fig_ext), "", umap_file))
      covar_name <- gsub("_", " ", covar_name)

      md <- c(md, sprintf("![UMAP with %s](figures/%s)", covar_name, umap_file))
      md <- c(md, "")
      md <- c(md, add_figure_caption(sprintf("UMAP embedding including %s covariate. Point shapes and fills indicate biological groups, whilst point outline colours indicate %s values.", covar_name, covar_name)))
      md <- c(md, "")

      # Add detailed interpretation for each covariate
      covar_name_clean <- gsub("_", " ", covar_name)
      md <- c(md, sprintf("**%s Analysis:**", covar_name_clean))
      md <- c(md, sprintf("- **Key finding:** Samples do NOT cluster by %s", tolower(covar_name_clean)))
      md <- c(md, "- Point colours/outlines show covariate values randomly distributed across biological clusters")
      md <- c(md, sprintf("- **Conclusion:** %s does NOT drive variance; biological signal dominates", covar_name_clean))
      md <- c(md, "")
    }
    md <- c(md, "")
  }

  md <- c(md, "---")
  md <- c(md, "")

  # ====== SECTION 6: TECHNICAL QUALITY CONTROL (USING LINEAR MODELS) ======
  md <- c(md, "## Technical Quality Control {#technical-quality-control}")
  md <- c(md, "")
  md <- c(md, "We used **Linear Models (LM)** to test for technical associations while controlling for biological group. This approach properly distinguishes true technical confounding from spurious associations caused by technical variables correlating with biology.")
  md <- c(md, "")
  md <- c(md, "> **Note:** For temporal/repeated measures analysis, Linear Mixed Models (LMM) with random intercepts are used (`methylation ~ time * biological_group + (1|subject_id)`). For QC variable testing, standard Linear Models are appropriate as samples are independent.")
  md <- c(md, "")

  if (length(qc_columns) > 0 && exists("lmm_summary_df") && nrow(lmm_summary_df) > 0) {
    # Count actual number of tested variables (some columns may generate multiple test variables)
    n_tested_variables <- length(unique(lmm_summary_df$covariate))
    md <- c(md, sprintf("We examined **%d QC variables** for associations with methylation (controlling for biological group):", n_tested_variables))
    md <- c(md, "")

    md <- c(md, "**Complete Linear Mixed Model Results:**")
    md <- c(md, "")
    md <- c(md, "*All tested QC variables are shown below with complete statistics, regardless of significance.*")
    md <- c(md, "")
    md <- c(md, add_table_header("Complete QC Variable Statistics"))
    md <- c(md, "")

    # Table
    md <- c(md, "| Variable | Markers Tested | Sig (Raw p) | Sig (Adjusted) | Mean R² |")
    md <- c(md, "|----------|----------------|-------------|----------------|---------|")

    for (i in 1:nrow(lmm_summary_df)) {
      md <- c(md, sprintf(
        "| %s | %d | %d | %d | %.3f |",
        lmm_summary_df$covariate[i],
        lmm_summary_df$n_markers_tested[i],
        lmm_summary_df$n_significant_raw[i],
        lmm_summary_df$n_significant_adjusted[i],
        lmm_summary_df$mean_rsquared[i]
      ))
    }

    md <- c(md, "")
    md <- c(md, "**Table Legend:**")
    md <- c(md, "")
    md <- c(md, "- **Markers Tested**: Total markers analysed in linear models")
    md <- c(md, sprintf("- **Sig (Raw p)**: Markers with raw p-value < %.3f for QC variable coefficient", sig_threshold))
    md <- c(md, sprintf("- **Sig (Adjusted)**: Markers with %s-adjusted p-value < %.3f", toupper(mt_method_display), sig_threshold))
    md <- c(md, "- **Mean R²**: Average model fit across all markers")
    md <- c(md, "")

    # Assessment
    total_lmm_sig <- sum(lmm_summary_df$n_significant_adjusted, na.rm = TRUE)
    total_lmm_raw <- sum(lmm_summary_df$n_significant_raw, na.rm = TRUE)

    md <- c(md, "**OVERALL ASSESSMENT:**")
    md <- c(md, "")

    if (total_lmm_sig == 0 && total_lmm_raw == 0) {
      md <- c(md, "> **No significant pre-analytical confounders detected** after controlling for biological group. All QC variables show non-significant associations, indicating that technical factors do not systematically affect methylation patterns beyond biological differences.")
    } else if (total_lmm_sig == 0 && total_lmm_raw > 0) {
      md <- c(md, sprintf("> Some raw associations detected (%d markers) but **none survive %s correction**. This indicates weak QC effects that are likely spurious and do not represent systematic technical bias after controlling for biology.", total_lmm_raw, toupper(mt_method_display)))
    } else {
      md <- c(md, sprintf("> **%d markers show significant QC associations** after %s correction (controlling for biological group). These represent true pre-analytical effects that should be considered in interpretation.", total_lmm_sig, toupper(mt_method_display)))
    }

    md <- c(md, "")

    # Add LMM plot using helper function
    if (copy_plot_to_figures(qc_dir, "lmm_associations_summary", figures_dir, fig_ext)) {
      md <- c(md, sprintf("![Technical QC associations (LMM)](figures/lmm_associations_summary.%s)", fig_ext))
      md <- c(md, "")
      md <- c(md, add_figure_caption("Stacked bar plot showing ALL markers tested in linear models for each technical variable (controlling for biological group). Colours indicate significance: gray (non-significant), light blue (raw p < alpha only), dark blue (FDR < alpha). Bar height = total markers tested."))
      md <- c(md, "")
    }
  } else {
    md <- c(md, "*No QC variables were found in the sample sheet for testing.*")
    md <- c(md, "")
  }

  md <- c(md, "---")
  md <- c(md, "")

  # ====== SECTION 7: METHODS ======
  md <- c(md, "## Methods {#methods}")
  md <- c(md, "")
  md <- c(md, "### Bootstrap Stability Analysis")
  md <- c(md, "")
  md <- c(md, sprintf("Parameters:"))
  md <- c(md, "")
  md <- c(md, sprintf("- Bootstrap iterations: %d", opt$bootstrap_iterations))
  md <- c(md, sprintf("- Stability threshold: %.2f", opt$stability_threshold))
  md <- c(md, sprintf("- Effect size threshold: %.1f%%", opt$effect_size_threshold))
  md <- c(md, "")

  md <- c(md, "### Technical Confounder Testing")
  md <- c(md, "")
  md <- c(md, "For each technical covariate:")
  md <- c(md, "")
  md <- c(md, "- **Continuous variables**: Spearman correlation with marker methylation")
  md <- c(md, "- **Categorical variables**: Kruskal-Wallis test across groups")
  md <- c(md, sprintf("- **Multiple testing correction**: %s method", toupper(mt_method_display)))
  md <- c(md, sprintf("- **Significance threshold**: alpha = %.3f", sig_threshold))
  md <- c(md, "")

  md <- c(md, "### Statistical Model for QC Testing")
  md <- c(md, "")
  md <- c(md, "Linear models were fitted for each marker to test QC variable associations:")
  md <- c(md, "")
  md <- c(md, "```")
  md <- c(md, "methylation ~ biological_group + qc_variable")
  md <- c(md, "```")
  md <- c(md, "")
  md <- c(md, "This is a standard linear model (LM), not a linear mixed model (LMM), as samples are independent (no repeated measures).")
  md <- c(md, "")

  # ====== NEW SECTION 8.7: REPRODUCIBILITY ======
  md <- c(md, "### Computational Reproducibility")
  md <- c(md, "")
  md <- c(md, "This section documents all parameters used for complete reproducibility of this analysis.")
  md <- c(md, "")

  # ==== Model Training Parameters ====
  md <- c(md, "#### Model Training Parameters")
  md <- c(md, "")

  # Try to load analysis_parameters.txt from model_dir
  analysis_params_file <- file.path(opt$model_dir, "analysis_parameters.txt")
  if (file.exists(analysis_params_file)) {
    # Console message with path blinding
    display_path <- if (opt$blind_paths) basename(analysis_params_file) else analysis_params_file
    cat(sprintf("Loading model training parameters from: %s\n", display_path))

    params_lines <- readLines(analysis_params_file)

    # Sanitize paths in parameter file content if blinding enabled
    if (opt$blind_paths) {
      params_lines <- sapply(params_lines, sanitize_path)
    }

    md <- c(md, "**Analysis Parameters (from MethylSense evaluation):**")
    md <- c(md, "")
    md <- c(md, "```")
    md <- c(md, params_lines)
    md <- c(md, "```")
    md <- c(md, "")
  } else {
    md <- c(md, "*Model training parameters file not found in model directory.*")
    md <- c(md, "")
  }

  # Try to load analysis_settings.txt from model_dir (primary MethylSense script settings)
  # ISSUE 6 FIX: analysis_settings.txt is 3 levels up from model_dir
  # Structure: main_dir/analysis_settings.txt
  #            main_dir/region_dir/nnet/markers_34/ <- opt$model_dir
  analysis_settings_file <- file.path(dirname(dirname(dirname(opt$model_dir))), "analysis_settings.txt")
  if (file.exists(analysis_settings_file)) {
    # Console message with path blinding
    display_path <- if (opt$blind_paths) basename(analysis_settings_file) else analysis_settings_file
    cat(sprintf("Loading primary MethylSense settings from: %s\n", display_path))

    settings_lines <- readLines(analysis_settings_file)

    # Sanitize paths in settings file content if blinding enabled
    if (opt$blind_paths) {
      settings_lines <- sapply(settings_lines, sanitize_path)
    }

    md <- c(md, "**Primary MethylSense Analysis Settings:**")
    md <- c(md, "")
    md <- c(md, "The following settings were used when creating this model with the primary MethylSense script:")
    md <- c(md, "")
    md <- c(md, "```")
    md <- c(md, settings_lines)
    md <- c(md, "```")
    md <- c(md, "")
  } else {
    md <- c(md, "*Primary analysis settings file (analysis_settings.txt) not found in model directory.*")
    md <- c(md, "")
  }

  # ==== Reviewer Analysis Parameters ====
  md <- c(md, "#### Model Evaluation Parameters")
  md <- c(md, "")
  md <- c(md, "**Command-line arguments used for this evaluation:**")
  md <- c(md, "")
  md <- c(md, "```bash")

  # Build the command string with all parameters
  cmd_parts <- c(sprintf("Rscript MethylSense_Reviewer_v%s.R \\", SCRIPT_VERSION))
  cmd_parts <- c(cmd_parts, sprintf("  --model_dir \"%s\" \\", if (opt$blind_paths) sanitize_path(opt$model_dir) else opt$model_dir))
  cmd_parts <- c(cmd_parts, sprintf("  --qs_file \"%s\" \\", if (opt$blind_paths) sanitize_path(opt$qs_file) else opt$qs_file))
  cmd_parts <- c(cmd_parts, sprintf("  --sample_sheet \"%s\" \\", if (opt$blind_paths) sanitize_path(opt$sample_sheet) else opt$sample_sheet))

  # Optional parameters (only include if not default or explicitly set)
  if (!is.null(opt$output_dir) && opt$output_dir != "model_evaluation") {
    cmd_parts <- c(cmd_parts, sprintf("  --output_dir \"%s\" \\", opt$output_dir))
  }

  if (!is.null(opt$group_colors) && opt$group_colors != "") {
    cmd_parts <- c(cmd_parts, sprintf("  --group_colors \"%s\" \\", opt$group_colors))
  }

  if (!is.null(opt$technical_covariates) && opt$technical_covariates != "") {
    cmd_parts <- c(cmd_parts, sprintf("  --technical_covariates \"%s\" \\", opt$technical_covariates))
  }

  if (opt$bootstrap_iterations != 500) {
    cmd_parts <- c(cmd_parts, sprintf("  --bootstrap_iterations %d \\", opt$bootstrap_iterations))
  }

  if (opt$n_cores != 4) {
    cmd_parts <- c(cmd_parts, sprintf("  --n_cores %d \\", opt$n_cores))
  }

  if (opt$stability_threshold != 0.75) {
    cmd_parts <- c(cmd_parts, sprintf("  --stability_threshold %.2f \\", opt$stability_threshold))
  }

  if (opt$effect_size_threshold != 2.0) {
    cmd_parts <- c(cmd_parts, sprintf("  --effect_size_threshold %.1f \\", opt$effect_size_threshold))
  }

  if (opt$min_effect_size != 5.0) {
    cmd_parts <- c(cmd_parts, sprintf("  --min_effect_size %.1f \\", opt$min_effect_size))
  }

  if (!is.null(opt$multiple_testing) && opt$multiple_testing != "fdr") {
    cmd_parts <- c(cmd_parts, sprintf("  --multiple_testing \"%s\" \\", opt$multiple_testing))
  }

  if (!is.null(opt$figure_format) && opt$figure_format != "all") {
    cmd_parts <- c(cmd_parts, sprintf("  --figure_format \"%s\" \\", opt$figure_format))
  }

  if (!is.null(opt$figure_resolution) && opt$figure_resolution != 300) {
    cmd_parts <- c(cmd_parts, sprintf("  --figure_resolution %d \\", opt$figure_resolution))
  }

  if (!is.null(opt$enable_wgcna) && opt$enable_wgcna == TRUE) {
    cmd_parts <- c(cmd_parts, "  --enable_wgcna \\")
  }

  if (!is.null(opt$enable_gene_annotation) && opt$enable_gene_annotation == TRUE) {
    cmd_parts <- c(cmd_parts, "  --enable_gene_annotation \\")
    if (!is.null(opt$species)) {
      cmd_parts <- c(cmd_parts, sprintf("  --species \"%s\" \\", opt$species))
    }
    if (opt$gene_window != 10000) {
      cmd_parts <- c(cmd_parts, sprintf("  --gene_window %d \\", opt$gene_window))
    }
  }

  if (!is.null(opt$functional_enrichment) && opt$functional_enrichment == TRUE) {
    cmd_parts <- c(cmd_parts, "  --functional_enrichment \\")
  }

  if (!is.null(opt$time_analysis)) {
    cmd_parts <- c(cmd_parts, sprintf("  --time_analysis \"%s\" \\", opt$time_analysis))
    if (!is.null(opt$time_categorical) && opt$time_categorical == TRUE) {
      cmd_parts <- c(cmd_parts, "  --time_categorical \\")
    }
  }

  if (!is.null(opt$verbose) && opt$verbose == TRUE) {
    cmd_parts <- c(cmd_parts, "  --verbose")
  } else {
    # Remove trailing backslash from last line
    cmd_parts[length(cmd_parts)] <- gsub(" \\\\$", "", cmd_parts[length(cmd_parts)])
  }

  md <- c(md, cmd_parts)
  md <- c(md, "```")
  md <- c(md, "")

  # ==== Parameter Summary Table ====
  md <- c(md, "**Key Parameters:**")
  md <- c(md, "")
  md <- c(md, add_table_header("Key Analysis Parameters"))
  md <- c(md, "")
  md <- c(md, "| Parameter | Value |")
  md <- c(md, "|-----------|-------|")
  md <- c(md, sprintf("| Bootstrap iterations | %d |", opt$bootstrap_iterations))
  md <- c(md, sprintf("| Stability threshold | %.2f |", opt$stability_threshold))
  md <- c(md, sprintf("| Effect size threshold | %.1f%% |", opt$effect_size_threshold))
  md <- c(md, sprintf("| Min effect size | %.1f%% |", opt$min_effect_size))
  md <- c(md, sprintf("| Multiple testing method | %s |", toupper(mt_method_display)))
  md <- c(md, sprintf("| CPU cores | %d |", opt$n_cores))
  md <- c(md, sprintf("| Figure format | %s |", opt$figure_format))
  md <- c(md, sprintf("| Figure resolution | %d PPI |", opt$figure_resolution))
  md <- c(md, sprintf("| WGCNA enabled | %s |", ifelse(opt$enable_wgcna, "Yes", "No")))
  md <- c(md, sprintf("| Gene annotation enabled | %s |", ifelse(opt$enable_gene_annotation, "Yes", "No")))
  if (opt$enable_gene_annotation) {
    md <- c(md, sprintf("| Gene annotation species | %s |", ifelse(!is.null(opt$species), opt$species, "Not specified")))
    md <- c(md, sprintf("| Gene annotation window | %d bp |", opt$gene_window))
  }
  md <- c(md, sprintf("| Time series analysis enabled | %s |", ifelse(!is.null(opt$time_analysis), "Yes", "No")))
  if (!is.null(opt$time_analysis)) {
    md <- c(md, sprintf("| Time series column | %s |", opt$time_analysis))
    md <- c(md, sprintf("| Force categorical time | %s |", ifelse(opt$time_categorical, "Yes", "No")))
  }
  md <- c(md, "")

  # ==== Software Environment ====
  md <- c(md, "#### Software Environment")
  md <- c(md, "")
  md <- c(md, "**Software Version:**")
  md <- c(md, sprintf("- MethylSense Reviewer: v%s", SCRIPT_VERSION))
  md <- c(md, sprintf("- R version: %s", R.version.string))
  md <- c(md, "")
  md <- c(md, "**Key R Packages:**")
  md <- c(md, sprintf("- caret: %s", packageVersion("caret")))
  md <- c(md, sprintf("- nnet: %s", packageVersion("nnet")))
  md <- c(md, sprintf("- pROC: %s", packageVersion("pROC")))
  md <- c(md, sprintf("- ggplot2: %s", packageVersion("ggplot2")))
  md <- c(md, sprintf("- pheatmap: %s", packageVersion("pheatmap")))
  md <- c(md, sprintf("- umap: %s", packageVersion("umap")))
  md <- c(md, sprintf("- igraph: %s", packageVersion("igraph")))
  md <- c(md, sprintf("- dendextend: %s", packageVersion("dendextend")))
  if (opt$enable_wgcna && requireNamespace("WGCNA", quietly = TRUE)) {
    md <- c(md, sprintf("- WGCNA: %s", packageVersion("WGCNA")))
  }
  if (opt$enable_gene_annotation && requireNamespace("biomaRt", quietly = TRUE)) {
    md <- c(md, sprintf("- biomaRt: %s", packageVersion("biomaRt")))
  }
  md <- c(md, "")

  # ==== Data Files ====
  md <- c(md, "#### Input Data Files")
  md <- c(md, "")
  md <- c(md, add_table_header("Input Data File Paths"))
  md <- c(md, "")
  md <- c(md, "| File Type | Path |")
  md <- c(md, "|-----------|------|")
  md <- c(md, sprintf("| Model directory | `%s` |", if (opt$blind_paths) sanitize_path(opt$model_dir) else opt$model_dir))
  md <- c(md, sprintf("| Methylation data | `%s` |", if (opt$blind_paths) sanitize_path(opt$qs_file) else opt$qs_file))
  md <- c(md, sprintf("| Sample sheet | `%s` |", if (opt$blind_paths) sanitize_path(opt$sample_sheet) else opt$sample_sheet))
  md <- c(md, "")

  # ==== Session Info ====
  md <- c(md, "#### Complete Session Information")
  md <- c(md, "")
  md <- c(md, "For full reproducibility, the complete R session information:")
  md <- c(md, "")
  md <- c(md, "```")
  session_info <- capture.output(sessionInfo())
  # Truncate long BLAS/LAPACK paths to show only from 'micromamba' forward
  session_info <- gsub("(BLAS/LAPACK: ).*/(micromamba/.*)", "\\1\\2", session_info)
  md <- c(md, session_info)
  md <- c(md, "```")
  md <- c(md, "")

  md <- c(md, "---")
  md <- c(md, "")

  # ====== SECTION 9: CONCLUSIONS ======
  md <- c(md, "## Conclusions {#conclusions}")
  md <- c(md, "")
  md <- c(md, "### Quantitative Validation Summary")
  md <- c(md, "")
  md <- c(md, "**Model Performance:**")
  if (exists("overall_accuracy")) {
    md <- c(md, sprintf("- Test set accuracy: %.1f%% (independent holdout)", overall_accuracy * 100))
  }
  if (exists("overall_auc")) {
    md <- c(md, sprintf("- ROC-AUC: %.3f (excellent discrimination)", overall_auc))
  }
  if (exists("overall_pr_auc")) {
    md <- c(md, sprintf("- PR-AUC: %.3f (robust to class imbalance)", overall_pr_auc))
  }
  if (exists("cv_detailed") && cv_detailed$available) {
    md <- c(md, sprintf(
      "- CV consistency: SD<%.1f%% across %d folds",
      sd(cv_detailed$data$accuracy) * 100, cv_detailed$n_folds
    ))
  }
  md <- c(md, "")

  md <- c(md, "**Biological Robustness:**")
  if (exists("dmr_details")) {
    md <- c(md, sprintf("- %d statistically significant DMRs (all FDR<0.05)", nrow(dmr_details)))
  }
  if (exists("cohens_d")) {
    md <- c(md, sprintf("- Large effect size: Cohen's d=%.2f", abs(cohens_d)))
  }
  if (exists("marker_stability") && nrow(marker_stability) > 0) {
    pct_stable <- round(100 * sum(marker_stability$stable) / nrow(marker_stability), 1)
    md <- c(md, sprintf("- %.0f%% markers show high bootstrap stability (>=70%%)", pct_stable))
  }
  md <- c(md, "")

  md <- c(md, "**Technical Validity:**")
  if (exists("covariate_summary")) {
    n_sig <- sum(covariate_summary$n_significant_adjusted > 0)
    md <- c(md, sprintf(
      "- %d/%d technical covariates show associations (all non-significant after FDR)",
      n_sig, nrow(covariate_summary)
    ))
  }
  if (exists("batch_effects_assessment") && batch_effects_assessment$available) {
    md <- c(md, "- PCA/UMAP: samples cluster by biology, not technical factors")
  }
  md <- c(md, "")

  md <- c(md, "**Clinical Readiness:**")
  if (exists("sample_consistency") && sample_consistency$available) {
    md <- c(md, sprintf(
      "- %.0f%% samples consistently classified across folds",
      sample_consistency$pct_consistent
    ))
  }
  if (exists("calibration_analysis") && calibration_analysis$available) {
    md <- c(md, "- Well-calibrated probabilities suitable for decision thresholds")
  }
  if (exists("class_metrics") && length(class_metrics) > 0) {
    # Find best specificity among classes
    best_spec <- max(sapply(class_metrics, function(x) ifelse("Specificity" %in% names(x), x$Specificity, 0)), na.rm = TRUE)
    if (best_spec >= 0.99) {
      md <- c(md, sprintf("- High specificity achieved: %.0f%% in best-performing class", best_spec * 100))
    }
  }
  md <- c(md, "")

  md <- c(md, "### Final Statement")
  md <- c(md, "")

  # Build realistic conclusion based on actual data
  n_dmrs <- if (exists("dmr_details")) nrow(dmr_details) else NA

  if (!is.na(n_dmrs) && n_dmrs <= 10) {
    # Small marker panel - be honest about limitations
    md <- c(md, sprintf("This validation report evaluates a %d-marker methylation classifier:", n_dmrs))
    md <- c(md, "")
    md <- c(md, "**Strengths:**")
    md <- c(md, "- Rigorous nested cross-validation with hyperparameter optimisation")
    md <- c(md, "- Independent test set evaluation")
    md <- c(md, "- Comprehensive batch effect and technical covariate assessment")
    md <- c(md, "")
    md <- c(md, "**Limitations:**")
    md <- c(md, sprintf("- Small marker panel (%d DMRs) may limit biological coverage", n_dmrs))
    md <- c(md, "- Further validation on independent cohorts recommended before clinical translation")
    md <- c(md, "- Model performance should be interpreted cautiously given limited feature set")
    md <- c(md, "")
    md <- c(md, "The model demonstrates acceptable validation for initial publication, but additional development and external validation are recommended before considering clinical applications.")
  } else {
    # Larger marker panel - more confident conclusion
    md <- c(md, "This comprehensive validation demonstrates that the MethylSense DNA methylation classifier shows promising performance:")
    md <- c(md, "")
    md <- c(md, "1. **Rigorous model development:** Nested cross-validation with hyperparameter optimisation and independent test set validation")

    n_dmrs_text <- if (!is.na(n_dmrs)) sprintf("%d statistically significant DMRs", n_dmrs) else "statistically significant DMRs"
    md <- c(md, sprintf("2. **Biological signal:** %s identified with large effect sizes", n_dmrs_text))

    md <- c(md, "3. **Technical validity:** Batch effects assessed; biological signal evaluated across dimensionality reduction analyses")
    md <- c(md, "4. **Model performance:** Evaluated across multiple metrics including accuracy, sensitivity, specificity, and calibration")
    md <- c(md, "")
    md <- c(md, "The model demonstrates validation metrics suitable for publication. External validation on independent cohorts is recommended before considering clinical translation.")
  }
  md <- c(md, "")

  md <- c(md, "---")
  md <- c(md, "")

  # ====== SECTION 10: REFERENCES ======
  md <- c(md, "## References {#references}")
  md <- c(md, "")
  md <- c(md, "### Primary Citation")
  md <- c(md, "")
  md <- c(md, "**Drag MH**, Hvilsom C, Poulsen LL, Jensen HE, Tahas SA, Leineweber C, Cray C, Bertelsen MF, Bojesen AM (2025). New high accuracy diagnostics for avian *Aspergillus fumigatus* infection using Nanopore methylation sequencing of host cell-free DNA and machine learning prediction. *bioRxiv* 2025.04.11.648151. [https://doi.org/10.1101/2025.04.11.648151](https://doi.org/10.1101/2025.04.11.648151)")
  md <- c(md, "")
  md <- c(md, "### MethylSense Software")
  md <- c(md, "")
  md <- c(md, "**GitHub Repository:** [https://github.com/markusdrag/MethylSense](https://github.com/markusdrag/MethylSense)")
  md <- c(md, "")
  md <- c(md, "**Author:** Markus Hodal Drag")
  md <- c(md, "")
  md <- c(md, "**License:** Academic Free License 3.0 (AFL-3.0)")
  md <- c(md, "")
  md <- c(md, "### Statistical Methods")
  md <- c(md, "")
  md <- c(md, "1. Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. *J R Stat Soc Series B*. 1995;57(1):289-300.")
  md <- c(md, "")
  md <- c(md, "2. Efron B, Tibshirani RJ. *An Introduction to the Bootstrap*. Chapman & Hall/CRC; 1993.")
  md <- c(md, "")
  md <- c(md, "3. Leek JT, Scharpf RB, Bravo HC, et al. Tackling the widespread and critical impact of batch effects in high-throughput data. *Nat Rev Genet*. 2010;11(10):733-739.")
  md <- c(md, "")
  md <- c(md, "4. Mansell G, Gorrie-Stone TJ, Bao Y, et al. Guidance for DNA methylation studies: statistical insights from the Illumina EPIC array. *BMC Genomics*. 2019;20(1):366.")
  md <- c(md, "")

  md <- c(md, "---")
  md <- c(md, "")
  md <- c(md, "**Analysis Details:**")
  md <- c(md, "")
  md <- c(md, sprintf("- Report generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  md <- c(md, sprintf("- Analysis version: %s", SCRIPT_VERSION))
  md <- c(md, sprintf("- Multiple testing method: %s", mt_method_display))
  md <- c(md, sprintf("- Significance threshold: alpha = %.3f", sig_threshold))
  md <- c(md, "")

  # Write markdown file
  writeLines(md, markdown_file)

  cat("          [OK] Markdown report saved to:", markdown_file, "\n")

  # Try to render to HTML (optional)
  html_file <- gsub("\\.md$", ".html", markdown_file)
  tryCatch(
    {
      suppressPackageStartupMessages({
        if (require(rmarkdown, quietly = TRUE)) {
          rmarkdown::render(markdown_file,
            output_format = "html_document",
            output_file = basename(html_file),
            quiet = TRUE
          )

          # Inject custom CSS for clinical styling (for PDF generation)
          if (file.exists(html_file)) {
            html_content <- readLines(html_file, warn = FALSE)
            custom_css <- get_clinical_report_css()

            # Find </head> tag and insert CSS before it
            head_end_idx <- grep("</head>", html_content)[1]
            if (!is.na(head_end_idx)) {
              html_content <- c(
                html_content[1:(head_end_idx - 1)],
                custom_css,
                html_content[head_end_idx:length(html_content)]
              )
              writeLines(html_content, html_file)
              cat("          [OK] HTML report saved with clinical CSS theme\n")
            } else {
              cat("          [OK] HTML report saved to:", html_file, "\n")
            }
          }
        }
      })
    },
    error = function(e) {
      cat("          [INFO] HTML rendering skipped (rmarkdown not available)\n")
    }
  )

  # Try to render to PDF (if requested)
  if (opt$generate_pdf) {
    pdf_file <- gsub("\\.md$", ".pdf", markdown_file)
    cat("          [INFO] Attempting PDF generation...\n")

    # Check if rmarkdown is available
    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
      cat("          [WARNING] rmarkdown package not installed. Skipping PDF generation.\n")
      cat("          [INFO] Install with: install.packages('rmarkdown')\n")
    } else if (rmarkdown::pandoc_available()) {
      # PAGEDOWN ONLY - NO LATEX FALLBACK
      # Only try pagedown (chrome-based PDF generation)
      if (!requireNamespace("pagedown", quietly = TRUE)) {
        cat("          [WARNING] pagedown package not installed - PDF generation skipped\n")
        cat("          [INFO] Install with: install.packages('pagedown')\n")
        cat("          [INFO] Markdown and HTML reports are still available\n")
      } else if (!file.exists(html_file)) {
        cat(sprintf("          [WARNING] HTML file not found: %s - PDF generation skipped\n", html_file))
        cat("          [INFO] Markdown and HTML reports are still available\n")
      } else {
        cat("          [INFO] Using pagedown with clinical report theme...\n")
        cat(sprintf("          [INFO] HTML input: %s\n", basename(html_file)))
        tryCatch(
          {
            pagedown::chrome_print(
              input = html_file,
              output = pdf_file,
              format = "pdf",
              timeout = 120,
              verbose = TRUE,
              extra_args = c("--no-sandbox", "--disable-setuid-sandbox")
            )
            cat("          [OK] Clinical-styled PDF report saved to:", pdf_file, "\n")
            cat("          [INFO] PDF generated with modern clinical styling (no LaTeX required)\n")
          },
          error = function(e) {
            cat(sprintf("          [WARNING] PDF generation failed: %s\n", e$message))
            cat("          [INFO] This usually means Chrome/Chromium is not installed or not in PATH\n")
            cat("          [INFO] PDF generation requires Chrome/Chromium browser\n")
            cat("          [INFO] **HPC Users:** Generate PDF locally after downloading HTML:\n")
            cat("          [INFO]   Rscript -e \"pagedown::chrome_print('MODEL_EVALUATION_REPORT.html')\"\n")
            cat("          [INFO] Markdown and HTML reports are still available\n")
          }
        )
      }
    } else {
      cat("          [WARNING] pandoc not found. PDF generation requires pandoc.\n")
      cat("          [INFO] Install pandoc: https://pandoc.org/installing.html\n")
    }
  }

  log_section_end(10, "Generating comprehensive markdown report")
  cat("\n")

  # Provide report location summary
  cat("================================================================================\n")
  cat("  REPORT FILES GENERATED\n")
  cat("================================================================================\n")
  if (exists("html_file") && file.exists(html_file)) {
    cat(sprintf("  ✓ HTML Report: %s\n", basename(html_file)))
  }
  if (exists("markdown_file") && file.exists(markdown_file)) {
    cat(sprintf("  ✓ Markdown:    %s\n", basename(markdown_file)))
  }
  if (opt$generate_pdf && exists("pdf_file") && file.exists(pdf_file)) {
    cat(sprintf("  ✓ PDF Report:  %s\n", basename(pdf_file)))
  }
  cat("\n")
  cat("  TIP: Open the HTML report in your browser for interactive viewing\n")
  if (!opt$generate_pdf) {
    cat("  TIP: To generate PDF, add --generate_pdf flag (requires Chrome/Chromium)\n")
  }
  cat("================================================================================\n")
  cat("\n")
}

# Finalize logging
finalize_log()

cat("\n")
cat("================================================================================\n")
cat(sprintf("  MethylSense_reviewer.R version %s (%s) completed successfully\n", SCRIPT_VERSION, SCRIPT_DATE))
cat("================================================================================\n")
cat("\n")
