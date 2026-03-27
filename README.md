<div align="center">

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="logos/methylsense_logo_dark.png?v=1">
  <source media="(prefers-color-scheme: light)" srcset="logos/methylsense_logo_white.png?v=1">
  <img alt="MethylSense Logo" src="logos/methylsense_logo_white.png?v=1" width="400">
</picture>

**High-accuracy epigenetic diagnostics via Nanopore methylation sequencing**

[![License: AFL-3.0](https://img.shields.io/badge/Licence-AFL--3.0-blue.svg)](https://opensource.org/licenses/AFL-3.0)
[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/version-5.7.0-brightgreen.svg)](https://github.com/markusdrag/MethylSense)

</div>

---

## Overview

MethylSense is a sort-of "semi-automated" R-based bioinformatics pipeline for training, testing, and cross-validation of machine learning models to predict the presence of clinical infections or inflammatory processes using (drum roll) epigenetics! After feeding Nanopore or Illumina methylation data and running the main analysis script, it will compute associated DMRs with your condition and find the best predictive ML models for diagnosing new cases and samples. You will also gain insights into the epigenome and generate publication-ready figures in short time - and lots more!



### So what does MethylSense actually do?

- **Detects differentially methylated regions (DMRs)**: Identifies custom sized genomic regions with statistically significant methylation differences between multiple groups (e.g., Control, Infected, Suspected) using the methylKit R package (Akalin et al., 2012) with false discovery rate (FDR) correction. Also supports covariates.

- **Trains machine learning classifiers**: Builds diagnostic models using nine machine learning algorithms via the caret framework (Kuhn, 2008), including Random Forest, SVM, XGBoost, and Elastic Net. Key features include Monte Carlo cross-validation and nested cross-validation for robust, unbiased model evaluation.

- **Diagnoses new samples**: Applies trained models to classify patient samples with probability scores and confidence estimates to support clinical decision-making.

- **Generates publication-ready reports**: Creates comprehensive evaluation reports with ROC curves, confusion matrices, Random Forest decision-tree graphics (and other model-specific plots), and performance metrics suitable for publication. You will have a lot of figures ready!

Originally developed for avian Aspergillus fumigatus infection diagnosis using host cell-free DNA methylation, achieving ~92% accuracy with Monte Carlo cross-validation.


---

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Installation](#installation)
3. [Sample Sheet Format](#sample-sheet-format)
4. [Quick Start](#quick-start)
5. [Example Output](#example-output)
6. [Pipeline Components](#pipeline-components)
7. [Output Structure](#output-structure)
8. [Command-Line Reference](#command-line-reference)
9. [Cross-Validation Methodology](#cross-validation-methodology)
10. [Clinical Applications](#clinical-applications)
11. [Troubleshooting](#troubleshooting)
12. [Citation](#citation)
13. [Acknowledgements](#acknowledgements)
14. [License and Authors](#license-and-authors)

---

## Prerequisites

### Upstream Pipelines

MethylSense requires processed Nanopore methylation data. You must run these pipelines first to prepare your sequencing data:

**NanoporeToBED-Pipeline** ([GitHub](https://github.com/markusdrag/NanoporeToBED-Pipeline))

This pipeline converts raw Nanopore sequencing output (modBAM files from Guppy or Dorado basecallers) into standardised CpG BED files. It extracts 5-methylcytosine (5mC) modification calls at CpG dinucleotide sites, preserving methylation frequencies and read coverage information. Run this pipeline first before MethylSense.

**GenomeToWindows** ([GitHub](https://github.com/markusdrag/GenomeToWindows))

This utility generates genomic window BED files (typically 1kb, 5kb, 10kb, or 25kb) from a reference genome. These windows define the regions for DMR detection. Aggregating CpG-level data into larger genomic windows reduces noise from individual CpG site variation and provides sufficient statistical power for differential methylation testing.

### System Requirements

- **R version:** 4.0 or higher (4.2+ recommended)
- **Operating system:** Linux, macOS, or Windows

**Running MethylSense from the terminal:**

| Platform | How to run |
|----------|------------|
| **Linux** | Use the native Terminal application |
| **macOS** | Use the native Terminal app (Applications → Utilities → Terminal) |
| **Windows** | Use [MobaXterm](https://mobaxterm.mobatek.net/) (recommended) or [WSL2](https://learn.microsoft.com/en-us/windows/wsl/) with Ubuntu |

> **Note for Windows users:** MobaXterm provides a Unix-like terminal environment with built-in X11 forwarding for R plots. Install R within MobaXterm's environment or use WSL2 with a native Linux R installation for best compatibility.

---

## Installation

### Quick Install

```bash
git clone https://github.com/markusdrag/MethylSense.git
cd MethylSense

Rscript MethylSense_installer.R

chmod +x MethylSense_*.R
```

The installer script automatically installs all required packages from both CRAN and Bioconductor. Installation typically takes 15-45 minutes depending on your system and existing packages. The script will report success or failure for each package and provide troubleshooting guidance if needed.

### Manual Installation

If you prefer manual installation, MethylSense depends on packages from both CRAN and Bioconductor:

```r
# CRAN packages
install.packages(c(
  "optparse", "qs", "caret", "pROC", "ggplot2", "hrbrthemes", "dplyr", "tidyr",
  "randomForest", "e1071", "glmnet", "nnet", "class", "MASS",
  "klaR", "ranger", "xgboost", "pheatmap", "RColorBrewer",
  "data.table", "readxl", "reshape2", "viridis", "corrplot",
  "factoextra", "ggridges", "gridExtra", "ggrepel", "naivebayes",
  "tree", "rpart", "rpart.plot", "cluster", "jsonlite", "MLmetrics",
  "doParallel", "foreach", "pls", "PRROC", "igraph", "umap",
  "dendextend", "lme4", "nlme", "MuMIn", "ez", "fossil"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("methylKit", "GenomicRanges", "genomation", "regioneR", "IRanges"))
```

### Verify Installation

```bash
Rscript -e "library(methylKit); library(caret); library(xgboost)"
Rscript MethylSense_load_data.R --help
```

---

## Sample Sheet Format

MethylSense requires a sample metadata file (**Excel `.xlsx` or comma-separated `.csv`**) that describes your samples, their group assignments, and file locations. This information is essential for the pipeline to correctly load, process, and analyse your data.

### Required columns

These columns must be present in your sample sheet:

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-top: 2px solid #333; border-bottom: 1px solid #333;">
<th style="padding: 8px; text-align: left;">Column</th>
<th style="padding: 8px; text-align: left;">Description</th>
<th style="padding: 8px; text-align: left;">Example</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">ID</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Unique sample identifier (must match BED filenames)</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">sample_001</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Species</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Species name used for filtering samples</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Gallus_gallus</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">bedFileOrg</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Original BED filename from NanoporeToBED</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">sample_001.CpG.bed</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">bedFile</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Output filename for converted BED (user-specified)</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">sample_001_8col.bed</td></tr>
<tr style="border-bottom: 2px solid #333;"><td style="padding: 8px;">treatMethylkit</td><td style="padding: 8px;">Numeric group code encoding each sample's treatment or infection status for methylKit. Each unique integer represents a distinct group (e.g., 0 = Control, 1 = Suspected, 2 = Infected). These codes are mapped to human-readable names via <code>--group_names</code> in ascending numeric order, or explicitly with <code>--treatment_mapping</code> (e.g., <code>0=Control,1=Suspected,2=Infected</code>). Plot colours are then auto-assigned from the MethylSense palette based on group name keywords.</td><td style="padding: 8px;">0, 1, 2</td></tr>
</tbody>
</table>


### Optional columns

These columns provide additional metadata for reporting and stratification:

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-top: 2px solid #333; border-bottom: 1px solid #333;">
<th style="padding: 8px; text-align: left;">Column</th>
<th style="padding: 8px; text-align: left;">Description</th>
<th style="padding: 8px; text-align: left;">Example</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Infection</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Human-readable group label</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Control, Infected, Suspected</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Study</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Study or cohort identifier for stratified analysis</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Study_A</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Batch</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Batch identifier for covariate adjustment (use with --covariate_cols)</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Batch_1, Batch_2</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Pre-analytic QC (numeric)</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Optional technical covariates: reviewer auto-detects columns whose names match keywords such as <code>storage</code>, <code>dna</code>, <code>hemolysis</code>, <code>time</code>, <code>concentration</code>, <code>extraction</code>, <code>quality</code>, <code>qc</code> (see script). Use several variables for meaningful histograms.</td><td style="padding: 8px; border-bottom: 1px solid #ddd;"><code>storage_time_days</code>, <code>cfDNA_concentration_ng_uL</code>, <code>hemolysis_index</code>, …</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Animal_ID</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Individual animal identifier</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Bird_42</td></tr>
<tr style="border-bottom: 2px solid #333;"><td style="padding: 8px;">Lab_ID</td><td style="padding: 8px;">Laboratory sample identifier</td><td style="padding: 8px;">LAB001</td></tr>
</tbody>
</table>


### Example sample sheet

```
ID            Species        bedFileOrg            bedFile                treatMethylkit  Infection
control_1     Gallus_gallus  control_1.CpG.bed     control_1_8col.bed     0               Control
control_2     Gallus_gallus  control_2.CpG.bed     control_2_8col.bed     0               Control
suspected_1   Gallus_gallus  suspected_1.CpG.bed   suspected_1_8col.bed   1               Suspected
suspected_2   Gallus_gallus  suspected_2.CpG.bed   suspected_2_8col.bed   1               Suspected
infected_1    Gallus_gallus  infected_1.CpG.bed    infected_1_8col.bed    2               Infected
infected_2    Gallus_gallus  infected_2.CpG.bed    infected_2_8col.bed    2               Infected
```

See the `example_data` folder for a complete working example: **`sample_metadata_30.csv`** (and matching BED files) provides 30 samples (10 per class). You can also use **`sample_metadata.xlsx`** with the same column layout.

---

## Quick Start

<img src="logos/methylsense_workflow_simple.svg" alt="MethylSense Workflow" style="max-width: 700px; width: 100%; height: auto;">

### Step 1: Load and preprocess data

Runtime: 5-10 minutes

This step converts raw CpG BED files into a unified methylKit object (qs format). The pipeline normalises coverage across samples, filters low-quality CpG sites, and prepares the data for downstream analysis.

```bash
Rscript MethylSense_load_data.R \
  --species "Gallus_gallus" \
  --sample_sheet ./example_data/sample_metadata_30.csv \
  --bed_dir ./example_data \
  --output_dir ./preprocessed
```

What happens: Each sample's CpG methylation data is loaded, filtered by minimum coverage thresholds, and combined into a single methylRawList object for efficient analysis.

Output: YYYYMMDD_Species_nN_methylRaw.qs

### Step 2: Train diagnostic models

Runtime: 1-2 hours (with 10-fold cross-validation)

This is the core analysis step. MethylSense identifies differentially methylated regions (DMRs) using the methylKit package, then trains machine learning classifiers to distinguish between your experimental groups.

```bash
Rscript MethylSense_analysis.R \
  --qs_file ./preprocessed/*_methylRaw.qs \
  --sample_sheet ./example_data/sample_metadata_30.csv \
  --output_dir ./training \
  --window_file ./example_data/windows/windows_5kb.bed \
  --models rf,svm,xgboost \
  --group_names Control,Suspected,Infected \
  --group_colors '#1B5E20,#006064,#4A148C' \
  --positive_class Infected \
  --cv_repeats 10 \
  --nested_cv \
  --nested_cv_repeats 5
```

> **Tip — Group colours:** The `--group_colors` parameter accepts comma-separated hex codes or colour names, one per group in the same order as `--group_names`. If omitted, colours are auto-assigned from the MethylSense palette based on group name keywords (e.g., "Control" → green, "Infected" → purple). Override with any valid R colour, for example `--group_colors 'darkgreen,orange,darkred'`.

What happens:
1. Methylation levels are calculated within genomic windows
2. Statistical testing identifies DMRs with FDR correction (Benjamini-Hochberg)
3. DMR methylation values become features for machine learning
4. Multiple models are trained with cross-validation for robust performance estimates
5. Nested cross-validation provides unbiased accuracy estimates suitable for publication

Output: Trained models (rds files), DMR coordinates, cross-validation results, performance metrics

### Optional: Data quality assessment and visualisation

Runtime: 5-15 minutes

**This step is optional but recommended** for exploring your DMR results and generating publication-ready visualisations. Run this after Step 2 to understand methylation patterns, assess data quality, and create comprehensive figures for presentations or publications.

```bash
Rscript MethylSense_general_data_overview.R \
  --analysis_dir ./training/training_* \
  --sample_sheet ./example_data/sample_metadata_30.csv \
  --region_sizes "5000" \
  --output_dir ./data_overview \
  --plot_format png,pdf \
  --infection_col Infection \
  --sample_id_col ID \
  --study_col Study
```

What happens: The script analyses your DMR results and generates publication-ready figures including:
- DMR counts across different window sizes
- Volcano plots showing effect sizes and significance
- Coverage distribution analysis (quality control)
- Chromosomal distribution of DMRs
- Effect size characterisation (hypermethylated vs hypomethylated)
- Sample-level methylation heatmaps for top DMRs

Output: Comprehensive report with figures and statistics in data_overview/

**When to use:** Run this after Step 2 to visualise DMR results, compare different window sizes, and generate exploratory figures. This is particularly useful for understanding your data before selecting the best model for clinical diagnosis.

### Step 3: Model evaluation and selection

Runtime: 5-10 minutes

Before using a model for clinical diagnosis, evaluate its performance to select the best model. This step generates comprehensive reports with publication-ready figures to help you choose the optimal diagnostic model.

```bash
# Find a trained model directory (pick one based on marker count)
Rscript MethylSense_reviewer.R \
  --model_dir ./training/training_*/windows_5kb/rf/markers_5 \
  --qs_file ./preprocessed/*_methylRaw.qs \
  --sample_sheet ./example_data/sample_metadata_30.csv \
  --sample_id_column ID \
  --treatment_column Infection \
  --output_subdir ./report
```

What happens: The reviewer analyses model performance on the training data, generates ROC curves, confusion matrices, and detailed reports. Compare performance across different models (rf, svm, xgboost) to select the best one for your diagnostic application.

Output: MODEL_EVALUATION_REPORT.md, MODEL_EVALUATION_REPORT.html, figures, tables

### Step 4: Diagnose new samples

Runtime: 5-10 minutes (including preprocessing)

Apply your trained model to classify new unknown samples. This is the real-world clinical workflow where you have new patient samples that need diagnosis.


**Step 4a: Generate unknown example samples (Tutorial only)**

For this tutorial, we provide a script to generate 3 "unknown" samples with hidden labels:

```bash
Rscript generate_unknown_samples.R
```

This creates example_data/unknown_samples/ with 3 BED files and a metadata sheet.


**Step 4b: Preprocess new samples**

Preprocess the unknown samples just like in Step 1. For unknown samples, the treatMethylkit and Infection columns should be NA:

```bash
Rscript MethylSense_load_data.R \
  --species "Gallus_gallus" \
  --sample_sheet ./example_data/unknown_samples/unknown_samples_metadata.xlsx \
  --bed_dir ./example_data/unknown_samples \
  --output_dir ./unknown_preprocessed
```


**Step 4c: Run predictions**

```bash
Rscript MethylSense_predict.R \
  --model_dir ./training/training_*/windows_5kb/rf/markers_5 \
  --qs_file ./unknown_preprocessed/*_methylRaw.qs \
  --output_dir ./unknown_predictions \
  --plots
```

What happens:
1. The prediction script loads your trained model and its DMR coordinates
2. It extracts methylation values from the new samples at those specific DMR regions
3. The model predicts each sample's infection status with a probability score
4. Samples with confidence ≥0.7 are classified; lower confidence samples are flagged for review

Output:
- *_predictions_detailed.csv - Full predictions with probabilities for each class
- *_predictions_clinical_report.csv - Simplified report for clinical use
- *_prediction_summary.txt - Overview of prediction results
- Diagnostic plots (if --plots flag used)

> **Tutorial Checkpoint:** Check your predictions! The 3 unknown samples have hidden true labels (Control, Infected, Suspected). Compare with your model's predictions to see accuracy.


**Example prediction output**


MethylSense produces detailed predictions with class probabilities for each sample:

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-top: 2px solid #333; border-bottom: 1px solid #333;">
<th style="padding: 8px; text-align: left;">Sample_ID</th>
<th style="padding: 8px; text-align: left;">True_Status</th>
<th style="padding: 8px; text-align: left;">Predicted_Class</th>
<th style="padding: 8px; text-align: left;">Confidence</th>
<th style="padding: 8px; text-align: right;">Prob_Control</th>
<th style="padding: 8px; text-align: right;">Prob_Suspected</th>
<th style="padding: 8px; text-align: right;">Prob_Infected</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">unknown_1</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Control</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Control</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">0.948 (High)</td><td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;"><strong>0.948</strong></td><td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">0.052</td><td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">0.000</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">unknown_2</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Infected</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Infected</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">0.850 (High)</td><td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">0.006</td><td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">0.144</td><td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;"><strong>0.850</strong></td></tr>
<tr style="border-bottom: 2px solid #333;"><td style="padding: 8px;">unknown_3</td><td style="padding: 8px;">Suspected</td><td style="padding: 8px;">Suspected</td><td style="padding: 8px;">0.820 (High)</td><td style="padding: 8px; text-align: right;">0.002</td><td style="padding: 8px; text-align: right;"><strong>0.820</strong></td><td style="padding: 8px; text-align: right;">0.178</td></tr>
</tbody>
</table>


*All 3 samples correctly classified with high confidence (≥0.7). The model correctly identified one sample from each class.*

---

## Example Output

MethylSense generates publication-ready visualisations across its pipeline scripts. The examples below were generated from a 3-class analysis (Control, Suspected, Infected) using 30 simulated samples with 10-fold nested cross-validation. **Many additional plots are generated by each script—see the output directories for the full suite.**

### From MethylSense_analysis.R

Plots generated during model training (use --cpg_report for CpG-level plots).

<table>
<tr>
<td width="33%">

**Confusion matrix**

<img src="screenshots/confusion_matrix.png" alt="Confusion Matrix" style="max-width: 100%; height: auto;">

</td>
<td width="33%">

**DMR barplot**

<img src="screenshots/dmr_barplot.png" alt="DMR Barplot" style="max-width: 100%; height: auto;">

</td>
<td width="33%">

**Group methylation summary**

<img src="screenshots/group_methylation_summary.png" alt="Group Methylation" style="max-width: 100%; height: auto;">

</td>
</tr>
<tr>
<td width="33%">

**Methylation heatmap**

<img src="screenshots/methylation_heatmap.png" alt="Heatmap" style="max-width: 100%; height: auto;">

</td>
<td width="33%">

**Prediction probabilities**

<img src="screenshots/prediction_probabilities_barplot.png" alt="Predictions" style="max-width: 100%; height: auto;">

</td>
<td width="33%">

**Random Forest — example decision tree**

<img src="screenshots/rf_decision_tree_example.png" alt="RF example decision tree" style="max-width: 100%; height: auto;">

*Illustrative classification tree for the selected DMR features*

</td>
</tr>
<tr>
<td colspan="3">

**Genome-wide DMR Manhattan plot**

<img src="screenshots/manhattan_plot.png" alt="Manhattan Plot" style="max-width: 100%; height: auto;">

</td>
</tr>
<tr>
<td colspan="3">

**CpG-level profile (DMR #4)**

<img src="screenshots/dmr_cpg_profile.png" alt="CpG Profile" style="max-width: 100%; height: auto;">

</td>
</tr>
</table>


### From MethylSense_reviewer.R

Model evaluation plots generated when reviewing a trained model.

<table>
<tr>
<td width="50%">

**UMAP dimensionality reduction**

<img src="screenshots/umap_methylation.png" alt="UMAP" style="max-width: 100%; height: auto;">

*Sample clustering based on methylation biomarkers*

</td>
<td width="50%">

**Multiclass ROC curves**

<img src="screenshots/roc_curves.png" alt="ROC Curves" style="max-width: 100%; height: auto;">

*One-vs-rest AUC for each class*

</td>
</tr>
<tr>
<td width="50%">

**DMR marker stability**

<img src="screenshots/bootstrap_selection_frequencies.png" alt="Bootstrap Stability" style="max-width: 100%; height: auto;">

*Bootstrap selection frequency for each DMR marker*

</td>
<td width="50%">

**Marker effect size vs stability**

<img src="screenshots/marker_stability.png" alt="Effect Size vs Stability" style="max-width: 100%; height: auto;">

*DMR selection consistency against methylation difference*

</td>
</tr>
<tr>
<td width="50%">

**Preanalytic QC - technical covariate associations**

<img src="screenshots/lmm_associations_summary.png" alt="LMM Associations Summary" style="max-width: 100%; height: auto;">

*Linear model analysis of technical covariates vs DMR methylation with FDR correction*

</td>
<td width="50%">

**CpG-level Manhattan plot**

<img src="screenshots/cpg_manhattan_plot.png" alt="CpG Manhattan" style="max-width: 100%; height: auto;">

*Genome-wide significance of individual CpG sites within DMRs*

</td>
</tr>
</table>


### From MethylSense_general_data_overview.R

Exploratory plots providing a landscape view of DMR results.

<table>
<tr>
<td width="50%">

**Volcano plot**

<img src="screenshots/volcano_plot.png" alt="Volcano Plot" style="max-width: 100%; height: auto;">

*Effect size vs significance with top DMRs labeled*

</td>
<td width="50%">

**Methylation by infection status**

<img src="screenshots/ridge_methylation.png" alt="Ridge Methylation" style="max-width: 100%; height: auto;">

*Ridge plot of methylation distribution by group*

</td>
</tr>
</table>


---

## Pipeline Components

<img src="logos/methylsense_architecture_simple.svg" alt="MethylSense Architecture" style="max-width: 650px; width: 100%; height: auto;">

### Core scripts

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-top: 2px solid #333; border-bottom: 1px solid #333;">
<th style="padding: 8px; text-align: left;">Script</th>
<th style="padding: 8px; text-align: left;">Purpose</th>
<th style="padding: 8px; text-align: left;">Key Output</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">MethylSense_load_data.R</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Preprocess BED files into methylKit format</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">qs methylRawList object</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">MethylSense_analysis.R</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">DMR detection and ML model training</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Trained models, DMR coordinates</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">MethylSense_predict.R</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Clinical diagnosis of new samples</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Predictions with confidence scores</td></tr>
<tr style="border-bottom: 2px solid #333;"><td style="padding: 8px;">MethylSense_reviewer.R</td><td style="padding: 8px;">Clinical interpretation and reporting</td><td style="padding: 8px;">Evaluation reports and figures</td></tr>
</tbody>
</table>


### Optional scripts

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-top: 2px solid #333; border-bottom: 1px solid #333;">
<th style="padding: 8px; text-align: left;">Script</th>
<th style="padding: 8px; text-align: left;">Purpose</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">MethylSense_general_data_overview.R</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Data quality assessment and exploratory visualisation</td></tr>
<tr style="border-bottom: 2px solid #333;"><td style="padding: 8px;">MethylSense_installer.R</td><td style="padding: 8px;">Complete package installer with verification</td></tr>
</tbody>
</table>


### Available machine learning models

MethylSense supports nine machine learning algorithms, implemented through the caret framework (Kuhn, 2008):

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-top: 2px solid #333; border-bottom: 1px solid #333;">
<th style="padding: 8px; text-align: left;">Model</th>
<th style="padding: 8px; text-align: left;">Code</th>
<th style="padding: 8px; text-align: left;">Speed</th>
<th style="padding: 8px; text-align: left;">Best Use Case</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Random Forest</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">rf</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Medium</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">General purpose, feature importance analysis</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Support Vector Machine</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">svm</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Medium</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">High accuracy classification</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">XGBoost</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">xgboost</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Medium</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Best overall accuracy</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Elastic Net</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">glmnet</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Fast</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Biomarker discovery, interpretable models</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Neural Network</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">nnet</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Medium</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Complex non-linear patterns</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">k-Nearest Neighbours</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">knn</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Fast</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Simple datasets, baseline comparisons</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Linear Discriminant Analysis</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">lda</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Fast</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Highly interpretable results</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Naive Bayes</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">nb</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Fast</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Baseline comparisons</td></tr>
<tr style="border-bottom: 2px solid #333;"><td style="padding: 8px;">Fast Random Forest</td><td style="padding: 8px;">ranger</td><td style="padding: 8px;">Fast</td><td style="padding: 8px;">Large datasets</td></tr>
</tbody>
</table>


Model shortcuts: --models all, --models fast, --models biomarker

---

## Output Structure

### After MethylSense_analysis.R

```
output_directory/
├── analysis_summary.csv          # Complete results for all models
├── analysis_settings.txt         # Reproducibility settings
├── ANALYSIS_REPORT.md            # Summary report (Markdown)
├── ANALYSIS_REPORT.html          # Summary report (HTML)
├── run_log.txt                   # Complete execution log
├── ml_training_logs/             # Detailed ML training diagnostics
├── train_test_splits/            # Data partition information
├── precomputed_data/             # Cached intermediate results
├── complete_dmr_datasets/        # DMR statistics per window size
├── cv_analysis_plots/            # Cross-validation visualisations
└── models/                       # Trained model files
    ├── rf_model/
    │   ├── model.rds             # Serialised model object
    │   ├── dmr_coords.csv        # DMR coordinates used
    │   └── model_summary.txt     # Model performance summary
    ├── svm_model/
    └── xgboost_model/
```

### After MethylSense_reviewer.R

```
report_directory/
├── MODEL_EVALUATION_REPORT.md    # Comprehensive clinical report
├── MODEL_EVALUATION_REPORT.html  # Interactive HTML version
├── MODEL_EVALUATION_REPORT.pdf   # Publication-ready PDF
├── figures/                      # All visualisations
│   ├── roc_curves/
│   ├── confusion_matrices/
│   ├── feature_importance/
│   └── methylation_heatmaps/
└── tables/                       # Summary statistics
```

---

## Command-Line Reference

All MethylSense scripts support `--help` to display full parameter lists. This section documents parameters for publication reporting and advanced usage.

### MethylSense_load_data.R

Preprocesses raw BED files into methylKit objects for analysis.

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-bottom: 2px solid #333;">
<th style="padding: 10px; text-align: left;">Parameter</th>
<th style="padding: 10px; text-align: left;">Type</th>
<th style="padding: 10px; text-align: left;">Default</th>
<th style="padding: 10px; text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr style="border-bottom: 1px solid #ddd;"><td style="padding: 10px;">--species</td><td style="padding: 10px;">string</td><td style="padding: 10px;">required</td><td style="padding: 10px;">Species name (must match sample sheet)</td></tr>
<tr style="border-bottom: 1px solid #ddd;"><td style="padding: 10px;">--sample_sheet</td><td style="padding: 10px;">file</td><td style="padding: 10px;">required</td><td style="padding: 10px;">Path to sample metadata (Excel or CSV)</td></tr>
<tr style="border-bottom: 1px solid #ddd;"><td style="padding: 10px;">--bed_dir</td><td style="padding: 10px;">dir</td><td style="padding: 10px;">required</td><td style="padding: 10px;">Directory containing CpG BED files</td></tr>
<tr style="border-bottom: 1px solid #ddd;"><td style="padding: 10px;">--output_dir</td><td style="padding: 10px;">dir</td><td style="padding: 10px;">required</td><td style="padding: 10px;">Output directory for processed .qs files</td></tr>
<tr style="border-bottom: 1px solid #ddd;"><td style="padding: 10px;">--assembly</td><td style="padding: 10px;">string</td><td style="padding: 10px;">auto</td><td style="padding: 10px;">Genome assembly name</td></tr>
<tr style="border-bottom: 1px solid #ddd;"><td style="padding: 10px;">--cores</td><td style="padding: 10px;">int</td><td style="padding: 10px;">20</td><td style="padding: 10px;">Number of CPU cores for parallel processing</td></tr>
<tr style="border-bottom: 1px solid #ddd;"><td style="padding: 10px;">--min_coverage</td><td style="padding: 10px;">int</td><td style="padding: 10px;">1</td><td style="padding: 10px;">Minimum read coverage per CpG site</td></tr>
<tr style="border-bottom: 1px solid #ddd;"><td style="padding: 10px;">--force_convert</td><td style="padding: 10px;">flag</td><td style="padding: 10px;">false</td><td style="padding: 10px;">Force re-conversion of existing files</td></tr>
<tr><td style="padding: 10px;">--no_parallel</td><td style="padding: 10px;">flag</td><td style="padding: 10px;">false</td><td style="padding: 10px;">Disable parallel processing</td></tr>
</tbody>
</table>


**Example:**

```bash
Rscript MethylSense_load_data.R \
  --species "Gallus_gallus" \
  --sample_sheet ./samples.xlsx \
  --bed_dir ./bed_files \
  --output_dir ./preprocessed \
  --cores 8
```

---

### MethylSense_analysis.R

Main analysis script for DMR detection and ML classifier training.


**Core parameters**

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-bottom: 2px solid #333;">
<th style="padding: 10px; text-align: left;">Parameter</th>
<th style="padding: 10px; text-align: left;">Type</th>
<th style="padding: 10px; text-align: left;">Default</th>
<th style="padding: 10px; text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--qs_file</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">file</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">required</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Path to preprocessed methylRawList .qs file</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--sample_sheet</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">file</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">—</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Sample metadata (Excel/CSV) for covariate adjustment</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--output_dir</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">dir</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">required</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Output directory for results</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--window_file</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">file</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">—</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Path to a single genomic window BED file (e.g., windows_5kb.bed)</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--bed_files</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">—</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Comma-separated BED files for genomic windows (alternative to --window_file)</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--region_dir</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">dir</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">—</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Directory containing multiple window BED files (alternative to --window_file)</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--models</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">rf,svm,glmnet,nnet,knn,lda,nb</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">ML models to train</td></tr>
<tr><td style="padding: 10px;">--cores</td><td style="padding: 10px;">int</td><td style="padding: 10px;">4</td><td style="padding: 10px;">CPU cores for parallel model training</td></tr>
</tbody>
</table>


**Treatment group options**

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-bottom: 2px solid #333;">
<th style="padding: 10px; text-align: left;">Parameter</th>
<th style="padding: 10px; text-align: left;">Type</th>
<th style="padding: 10px; text-align: left;">Default</th>
<th style="padding: 10px; text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--group_names</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">auto</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Custom group names, comma-separated</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--group_colors</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">auto</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Custom colours matching group_names</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--positive_class</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">—</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Name of positive class for binary metrics</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--treatment_mapping</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">—</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Map numeric codes to names (e.g., 0=Control,1=Infected)</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--min_group_size</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">int</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">4</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Minimum samples per group</td></tr>
<tr><td style="padding: 10px;">--allow_small_groups</td><td style="padding: 10px;">flag</td><td style="padding: 10px;">false</td><td style="padding: 10px;">Allow groups with fewer samples (risky)</td></tr>
</tbody>
</table>


**Cross-validation & performance**

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-bottom: 2px solid #333;">
<th style="padding: 10px; text-align: left;">Parameter</th>
<th style="padding: 10px; text-align: left;">Type</th>
<th style="padding: 10px; text-align: left;">Default</th>
<th style="padding: 10px; text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--cv_folds</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">int</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">10</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Cross-validation folds</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--cv_repeats</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">int</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">0</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Number of CV repetitions (5-10 recommended for publication)</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--nested_cv</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">flag</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">false</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Enable nested cross-validation (gold standard for publication)</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--nested_cv_repeats</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">int</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">1</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Number of times to repeat the entire nested k-fold CV with different seeds. For publication: use 5. Total evaluations = cv_repeats × nested_cv_repeats</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--test_ratio</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">float</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">0.3</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Proportion of data for test set (70:30 train:test split)</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--train_test_file</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">file</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">—</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Excel file with predefined train/test splits</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--min_accuracy</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">float</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">0.7</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Minimum accuracy to save trained models</td></tr>
<tr><td style="padding: 10px;">--reproducible</td><td style="padding: 10px;">flag</td><td style="padding: 10px;">false</td><td style="padding: 10px;">Fixed seeds for reproducibility</td></tr>
</tbody>
</table>


**DMR detection & filtering**

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-bottom: 2px solid #333;">
<th style="padding: 10px; text-align: left;">Parameter</th>
<th style="padding: 10px; text-align: left;">Type</th>
<th style="padding: 10px; text-align: left;">Default</th>
<th style="padding: 10px; text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--qvalue_threshold</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">float</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">0.05</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">FDR threshold for DMR significance</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--min_coverage</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">int</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">20</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Minimum read coverage per CpG</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--hyper_threshold</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">float</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">5.0</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Minimum % for hypermethylated DMRs</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--hypo_threshold</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">float</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">-2.0</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Maximum % for hypomethylated DMRs</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--disable_meth_filter</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">flag</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">false</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Disable methylation filtering (q-value only)</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--covariate_cols</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">—</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Comma-separated sample sheet columns for DMR covariate adjustment (e.g., 'Batch,Study')</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--exclude_sex_chroms</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">flag</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">false</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Exclude X, Y, Z, W chromosomes from markers</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--filter_standard_chroms</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">flag</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">false</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Filter to standard chromosomes only (1-99 + sex)</td></tr>
<tr><td style="padding: 10px;">--min_samples_per_region</td><td style="padding: 10px;">float</td><td style="padding: 10px;">1.0</td><td style="padding: 10px;">Minimum sample fraction per region</td></tr>
</tbody>
</table>


**Marker set size**

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-bottom: 2px solid #333;">
<th style="padding: 10px; text-align: left;">Parameter</th>
<th style="padding: 10px; text-align: left;">Type</th>
<th style="padding: 10px; text-align: left;">Default</th>
<th style="padding: 10px; text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--max_markers</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">int</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">100</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Maximum number of DMRs to use as features</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--min_markers</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">int</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">8</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Minimum number of DMRs to start training</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--marker_step</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">int</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">1</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Step size for marker iteration</td></tr>
<tr><td style="padding: 10px;">--marker_range</td><td style="padding: 10px;">string</td><td style="padding: 10px;">—</td><td style="padding: 10px;">Alternative: specify range as min:max (e.g., 2:50)</td></tr>
</tbody>
</table>


**Plot output**

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-bottom: 2px solid #333;">
<th style="padding: 10px; text-align: left;">Parameter</th>
<th style="padding: 10px; text-align: left;">Type</th>
<th style="padding: 10px; text-align: left;">Default</th>
<th style="padding: 10px; text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--create_plots</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">flag</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">true</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Generate visualisations</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--no_plots</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">flag</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">false</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Disable plot generation</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--plot_format</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">png</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Output format: png, jpg, svg, pdf, all</td></tr>
<tr><td style="padding: 10px;">--plot_dpi</td><td style="padding: 10px;">int</td><td style="padding: 10px;">300</td><td style="padding: 10px;">Resolution (300 standard, 600 high-res, 150 draft)</td></tr>
</tbody>
</table>


**Other options**

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-bottom: 2px solid #333;">
<th style="padding: 10px; text-align: left;">Parameter</th>
<th style="padding: 10px; text-align: left;">Type</th>
<th style="padding: 10px; text-align: left;">Default</th>
<th style="padding: 10px; text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--force_ml</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">flag</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">false</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Force re-run even with cached results</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--verbose</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">flag</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">false</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Enable verbose output</td></tr>
<tr><td style="padding: 10px;">--debug</td><td style="padding: 10px;">flag</td><td style="padding: 10px;">false</td><td style="padding: 10px;">Enable debug mode</td></tr>
</tbody>
</table>


**Model shortcuts:** `--models all` (all 9), `--models fast` (glmnet,ranger,lda,nb), `--models biomarker` (glmnet,rf)


**Example — Full research analysis (Drag et al. 2025 setup):**

```bash
Rscript MethylSense_analysis.R \
  --qs_file ./preprocessed/Chicken_n124.qs \
  --output_dir ./results/Chicken_Aspergillosis \
  --region_dir ./windowed_genomes \
  --group_names 'Control,Aspergillus_Infected' \
  --group_colors 'darkgreen,darkred' \
  --positive_class 'Aspergillus_Infected' \
  --treatment_mapping '0=Control,1=Aspergillus_Infected' \
  --min_accuracy 0.8 \
  --marker_range '2:50' \
  --hyper_threshold 5.0 --hypo_threshold -2.0 \
  --min_coverage 10 \
  --cv_repeats 10 --nested_cv --nested_cv_repeats 5 \
  --test_ratio 0.4 \
  --exclude_sex_chroms --filter_standard_chroms \
  --models 'nnet,glmnet,lda,knn,rf,nb,svm' \
  --cores 8 \
  --plot_format 'all' --plot_dpi 300 \
  --force_ml
```

**Example — Quick test run:**

```bash
Rscript MethylSense_analysis.R \
  --qs_file ./data.qs \
  --output_dir ./quick_test \
  --group_names Control,Disease \
  --models fast \
  --max_markers 20
```

---

### MethylSense_predict.R

Applies trained models to classify new samples.

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-bottom: 2px solid #333;">
<th style="padding: 10px; text-align: left;">Parameter</th>
<th style="padding: 10px; text-align: left;">Type</th>
<th style="padding: 10px; text-align: left;">Default</th>
<th style="padding: 10px; text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--qs_file</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">file</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">required</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Path to preprocessed samples (.qs)</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--model_dir</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">dir</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">required</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Directory containing trained model</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--output_dir</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">dir</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">required</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Output directory for predictions</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--sample_ids</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">all</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Comma-separated sample IDs to predict</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--min_coverage</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">int</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">10</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Minimum coverage per CpG</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--confidence_threshold</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">float</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">0.7</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Confidence threshold for classification</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--plots</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">flag</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">false</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Generate diagnostic plots</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--batch_mode</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">flag</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">false</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Batch mode (suppress individual plots)</td></tr>
<tr><td style="padding: 10px;">--verbose</td><td style="padding: 10px;">flag</td><td style="padding: 10px;">false</td><td style="padding: 10px;">Enable verbose output</td></tr>
</tbody>
</table>

**Example:**

```bash
Rscript MethylSense_predict.R \
  --qs_file ./new_samples/samples.qs \
  --model_dir ./results/rf/markers_10 \
  --output_dir ./predictions \
  --plots --confidence_threshold 0.8
```

---

### MethylSense_reviewer.R

Generates comprehensive model evaluation reports for publication.


**Core parameters**

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-bottom: 2px solid #333;">
<th style="padding: 10px; text-align: left;">Parameter</th>
<th style="padding: 10px; text-align: left;">Type</th>
<th style="padding: 10px; text-align: left;">Default</th>
<th style="padding: 10px; text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--model_dir</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">dir</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">required</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Path to trained model directory</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--qs_file</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">file</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">required</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Path to methylRawList .qs file</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--sample_sheet</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">file</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">required</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Path to sample metadata</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--output_subdir</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">model_review</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Output subdirectory name</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--model_name</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">auto</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Custom model name for report title</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--n_cores</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">int</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">4</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Number of CPU cores</td></tr>
<tr><td style="padding: 10px;">--verbose</td><td style="padding: 10px;">flag</td><td style="padding: 10px;">false</td><td style="padding: 10px;">Enable verbose output</td></tr>
</tbody>
</table>


**Statistical analysis**

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-bottom: 2px solid #333;">
<th style="padding: 10px; text-align: left;">Parameter</th>
<th style="padding: 10px; text-align: left;">Type</th>
<th style="padding: 10px; text-align: left;">Default</th>
<th style="padding: 10px; text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--bootstrap_iterations</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">int</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">500</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Bootstrap iterations for CI estimation</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--stability_threshold</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">float</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">0.70</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Feature stability threshold</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--effect_size_threshold</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">float</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">2.0</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Effect size threshold for significance</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--min_effect_size</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">float</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">5.0</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Minimum effect size (methylation %)</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--multiple_testing</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">fdr</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Correction: fdr, bonferroni, holm, BY, none</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--significance_threshold</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">float</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">0.05</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">P-value significance threshold</td></tr>
<tr><td style="padding: 10px;">--control_group</td><td style="padding: 10px;">string</td><td style="padding: 10px;">auto</td><td style="padding: 10px;">Reference group for logistic regression</td></tr>
</tbody>
</table>


**Batch & time effects**

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-bottom: 2px solid #333;">
<th style="padding: 10px; text-align: left;">Parameter</th>
<th style="padding: 10px; text-align: left;">Type</th>
<th style="padding: 10px; text-align: left;">Default</th>
<th style="padding: 10px; text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--batch_column</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">—</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Column for batch/platform variable</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--storage_column</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">—</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Column for storage duration variable</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--date_column</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">—</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Column for collection/sampling date</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--time_analysis</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">—</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Column for timepoint/longitudinal data</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--time_categorical</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">flag</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">false</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Treat time as categorical (vs continuous)</td></tr>
<tr><td style="padding: 10px;">--covariate_cols</td><td style="padding: 10px;">string</td><td style="padding: 10px;">—</td><td style="padding: 10px;">Comma-separated column names from sample sheet to adjust for confounders during DMR calling (e.g., Batch,Storage_Days)</td></tr>
</tbody>
</table>


**Network & annotation**

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-bottom: 2px solid #333;">
<th style="padding: 10px; text-align: left;">Parameter</th>
<th style="padding: 10px; text-align: left;">Type</th>
<th style="padding: 10px; text-align: left;">Default</th>
<th style="padding: 10px; text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--enable_wgcna</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">flag</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">false</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Enable WGCNA co-methylation network analysis</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--enable_gene_annotation</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">flag</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">false</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Enable biomaRt gene annotation</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--species</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">—</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Species for annotation (e.g., 'Gallus gallus')</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--enable_feature_annotation</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">flag</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">false</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Enable genomic feature annotation</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--gtf_file</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">file</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">—</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">GTF/GFF file for gene features</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--regulatory_gff</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">file</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">—</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Regulatory features GFF (enhancers, CTCF)</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--emar_gff</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">file</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">—</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">EMAR regions GFF file</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--open_chromatin_bed</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">file</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">—</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Open chromatin BED file</td></tr>
<tr><td style="padding: 10px;">--functional_enrichment</td><td style="padding: 10px;">flag</td><td style="padding: 10px;">false</td><td style="padding: 10px;">Enable functional enrichment analysis</td></tr>
</tbody>
</table>


**PCR validation**

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-bottom: 2px solid #333;">
<th style="padding: 10px; text-align: left;">Parameter</th>
<th style="padding: 10px; text-align: left;">Type</th>
<th style="padding: 10px; text-align: left;">Default</th>
<th style="padding: 10px; text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--prepare_pcr</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">flag</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">false</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Generate PCR primer design links</td></tr>
<tr><td style="padding: 10px;">--pcr_n_dmrs</td><td style="padding: 10px;">int</td><td style="padding: 10px;">20</td><td style="padding: 10px;">Number of top DMRs for PCR validation</td></tr>
</tbody>
</table>


**Report output**

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-bottom: 2px solid #333;">
<th style="padding: 10px; text-align: left;">Parameter</th>
<th style="padding: 10px; text-align: left;">Type</th>
<th style="padding: 10px; text-align: left;">Default</th>
<th style="padding: 10px; text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--generate_markdown</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">flag</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">true</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Generate markdown report</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--generate_pdf</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">flag</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">false</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Generate PDF report</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--markdown_title</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">auto</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Custom report title</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--markdown_author</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">—</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Author name(s) for report</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--figure_format</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">png</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Output format: png, jpg, svg, all</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--figure_resolution</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">int</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">300</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Figure resolution (DPI)</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--group_colors</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">auto</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Custom colors for PCA/UMAP plots</td></tr>
<tr><td style="padding: 10px;">--blind_paths</td><td style="padding: 10px;">flag</td><td style="padding: 10px;">false</td><td style="padding: 10px;">Hide full paths in reports for privacy</td></tr>
</tbody>
</table>


**Example — Full research evaluation (Drag et al. 2025 setup):**

```bash
Rscript MethylSense_reviewer.R \
  --model_dir ./results/svm/markers_35 \
  --qs_file ./preprocessed/Chicken_n124.qs \
  --sample_sheet ./AsperMerged_v3.xlsx \
  --model_name "Fast Model Diagnostic Test" \
  --covariate_cols "SEQUENCER,cfDNA_conc_ng_uL,number_of_reads,mapping_percentage,mean_coverage" \
  --time_analysis Sample_Week \
  --enable_wgcna \
  --bootstrap_iterations 1000 \
  --stability_threshold 0.75 \
  --effect_size_threshold 2.0 \
  --min_effect_size 5.0 \
  --multiple_testing fdr \
  --n_cores 8 \
  --prepare_pcr \
  --functional_enrichment \
  --enable_feature_annotation \
  --gtf_file ./genomes/Gallus_gallus.GRCg7b.108.gtf \
  --regulatory_gff ./genomes/Gallus_gallus.GRCg7b.regulatory_features.gff3 \
  --emar_gff ./genomes/Gallus_gallus.GRCg7b.EMARs.gff \
  --open_chromatin_bed ./genomes/openchrom-union.bed \
  --figure_format "jpg,svg" \
  --blind_paths \
  --generate_pdf \
  --verbose
```

**Example — Quick batch effect check:**

```bash
Rscript MethylSense_reviewer.R \
  --model_dir ./results/rf/markers_20 \
  --qs_file ./data.qs \
  --sample_sheet ./metadata.xlsx \
  --batch_column Platform \
  --storage_column StorageDays
```

---

### MethylSense_general_data_overview.R

Creates visualisations of the methylation data landscape.

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-bottom: 2px solid #333;">
<th style="padding: 10px; text-align: left;">Parameter</th>
<th style="padding: 10px; text-align: left;">Type</th>
<th style="padding: 10px; text-align: left;">Default</th>
<th style="padding: 10px; text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--analysis_dir</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">dir</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">required</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Path to analysis directory</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--sample_sheet</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">file</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">required</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Path to sample metadata</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--region_sizes</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">1000,5000,10000,15000,20000,25000</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Region sizes (bp)</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--plot_format</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">string</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">png</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Output format: png, jpg, svg, pdf</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--plot_width</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">float</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">12</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Plot width in inches</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--plot_height</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">float</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">8</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Plot height in inches</td></tr>
<tr><td style="padding: 10px; border-bottom: 1px solid #ddd;">--plot_dpi</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">int</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">300</td><td style="padding: 10px; border-bottom: 1px solid #ddd;">Resolution for raster formats</td></tr>
<tr><td style="padding: 10px;">--volcano_meth_threshold</td><td style="padding: 10px;">float</td><td style="padding: 10px;">5</td><td style="padding: 10px;">Methylation difference threshold (%)</td></tr>
</tbody>
</table>


**Example — High-resolution SVG for publication:**

```bash
Rscript MethylSense_general_data_overview.R \
  --analysis_dir ./results \
  --sample_sheet ./metadata.xlsx \
  --region_sizes "5000,10000,25000" \
  --plot_format svg --plot_dpi 600 \
  --plot_width 14 --plot_height 10
```

---

### Plot Output Specifications

For publication-quality figures:

| Script | Default format | Default DPI | Recommended for print |
|--------|---------------|-------------|----------------------|
| MethylSense_analysis.R | PNG | 300 | --plot_format svg --plot_dpi 600 |
| MethylSense_reviewer.R | PNG | 300 | --figure_format 'png,svg' --figure_resolution 600 --generate_pdf |
| MethylSense_general_data_overview.R | PNG | 300 | --plot_format svg --plot_dpi 600 |

**Common journal requirements:**
- Nature/Science: 300-600 DPI, prefer vector formats (SVG/PDF)
- PLoS/BMC: 300 DPI minimum, TIFF or PNG acceptable
- Elsevier: 300 DPI for halftones, 1000 DPI for line art

All MethylSense scripts generate figures at **300 DPI** by default, suitable for most journal requirements.

---

## Cross-Validation Methodology

MethylSense implements two complementary cross-validation strategies for model evaluation. Both are controlled through command-line parameters and generate separate output files.

### Standard Monte Carlo Cross-Validation (MCCV)

Enabled with `--cv_repeats N` (where N > 0). This is the default and faster method.

Each repeat draws a **new random 60/40 stratified split**. The model is trained on 60% of the data with internal k-fold CV for hyperparameter tuning, and evaluated on the held-out 40%. Metrics are averaged across all N repeats with 95% confidence intervals.

```bash
# Standard MCCV with 10 repeats
Rscript MethylSense_analysis.R \
  --qs_file ./data.qs \
  --output_dir ./results \
  --cv_repeats 10
```

**Output:** `cv_detailed_*.csv`, `cv_summary_*.csv`

### Nested Cross-Validation (Gold Standard)

Enabled by adding `--nested_cv`. This is the gold standard for unbiased performance estimation suitable for publication (Varma & Simon, 2006).

Nested CV uses two layers of cross-validation:
- **Outer loop**: K stratified, non-overlapping folds (set by `--cv_repeats`). Each fold is held out as a test set exactly once.
- **Inner loop**: Separate k-fold CV within the outer training set for hyperparameter tuning. The outer test fold is never seen during tuning.

The entire k-fold procedure can be **repeated N times** with different random seeds (set by `--nested_cv_repeats`), producing K × N total evaluations for more robust estimates.

```bash
# Publication-quality nested CV: 10 folds × 5 repetitions = 50 evaluations
Rscript MethylSense_analysis.R \
  --qs_file ./data.qs \
  --output_dir ./results \
  --cv_repeats 10 \
  --nested_cv \
  --nested_cv_repeats 5
```

**Output:** `nested_cv_detailed_*.csv`, `nested_cv_summary_*.csv`, `nested_cv_hyperparameters_*.csv`, plus nested CV visualisation plots

### Comparison

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-top: 2px solid #333; border-bottom: 1px solid #333;">
<th style="padding: 8px; text-align: left;">Aspect</th>
<th style="padding: 8px; text-align: left;">Standard MCCV</th>
<th style="padding: 8px; text-align: left;">Nested CV</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Outer split</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Random 60/40 (overlapping across repeats)</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">K-fold, non-overlapping</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Hyperparameter tuning</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Internal CV on train split</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Separate inner CV loop on outer train only</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Test data isolation</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Good (not used for tuning)</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Complete (fully isolated from all tuning)</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Bias</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Potentially mildly optimistic</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Unbiased (gold standard)</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Speed</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Fast</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">5–10× slower</td></tr>
<tr style="border-bottom: 2px solid #333;"><td style="padding: 8px;">Use case</td><td style="padding: 8px;">Development and exploration</td><td style="padding: 8px;">Publication-quality results</td></tr>
</tbody>
</table>

When `--nested_cv` is enabled, MethylSense runs **both** methods automatically and calculates the **optimism bias** (the difference between standard and nested CV accuracy). If the bias exceeds 5%, a warning is raised suggesting the nested CV results should be used for publication.

---

## Clinical Applications

This section provides comprehensive guidance for veterinarians, diagnosticians, and clinicians using MethylSense for disease diagnosis. The pipeline translates complex epigenetic data into actionable clinical information.

### Understanding cfDNA-Based Diagnosis

Cell-free DNA (cfDNA) is released into the bloodstream when cells die, carrying the methylation patterns of their tissue of origin. During infection, pathogen-induced host cell damage releases cfDNA with characteristic methylation signatures that differ from healthy individuals. MethylSense detects these disease-associated methylation patterns to classify samples as infected or healthy.

The diagnostic approach is fundamentally different from pathogen detection methods (PCR, culture): rather than detecting the pathogen directly, MethylSense detects the host's epigenetic response to infection. This provides complementary diagnostic information and may detect infection earlier or in cases where pathogen load is below detection limits of conventional methods.

### Sample Types and Collection

MethylSense can analyse cfDNA from both serum and plasma. Both sample types contain host cfDNA suitable for methylation analysis:

**Serum:** Blood allowed to clot before centrifugation. May contain slightly higher cfDNA concentrations due to cell lysis during clotting. Used in the original validation study (Drag et al., 2025).

**Plasma:** Blood collected in anticoagulant tubes (EDTA, citrate, or heparin) and centrifuged before clotting. May have lower background from lysed blood cells. Either sample type is suitable, but consistency within a study is recommended.

### Sample Processing Specifications

For avian aspergillosis diagnosis, the following specifications are based on Drag et al. (2025):

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-top: 2px solid #333; border-bottom: 1px solid #333;">
<th style="padding: 8px; text-align: left;">Parameter</th>
<th style="padding: 8px; text-align: left;">Specification</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Sample type</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Serum or plasma (cell-free DNA)</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Volume</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">200 µL minimum (serum or plasma)</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Pre-extraction</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Centrifugation 2 min at 400 × g to remove residual cells</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">DNA extraction</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Norgen Biotek cfDNA Micro Kit (cat. 55500) or equivalent cfDNA kit</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Quality control</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Bioanalyzer 2100, High Sensitivity DNA kit (confirm cfDNA fragment profile)</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Storage</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">-80°C (extracted cfDNA stable for months)</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Library preparation</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Native Barcoding Ligation Kit (SQK-NBD114.24), 11 µL cfDNA input</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Sequencer</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">PromethION P2 Solo or MinION (any ONT platform)</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Flow cell</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">R10.4.1 (version 14 chemistry)</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Basecalling</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Guppy v6.5.7 or Dorado, High Accuracy mode</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Methylation model</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">dna_r10.4.1_e8.2_400bps_modbases_5mc_cg_hac.cfg</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Methylation caller</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">ONT modkit v0.1.2 or later, pileup mode with cpg combine mods</td></tr>
<tr style="border-bottom: 2px solid #333;"><td style="padding: 8px;">Reference genome</td><td style="padding: 8px;">Gallus gallus (bGalGal1.mat.broiler.GRCg7b) or appropriate species</td></tr>
</tbody>
</table>


Note: During library preparation, bead clean-up steps should be multiplied by 2.5× to improve recovery of short cfDNA fragments, as recommended by Martignano et al. (2023).

For complete methodology details, see Drag et al. (2025).

### Clinical Workflow

The complete diagnostic workflow from sample collection to clinical decision:

1. **Sample collection** — Collect 200 µL serum or plasma from patient. For serum, allow blood to clot at room temperature for 30-60 minutes before centrifugation. For plasma, collect in anticoagulant tube and centrifuge promptly.

2. **Sample processing** — Centrifuge to separate serum/plasma from cells. Perform a second centrifugation at 400 × g for 2 minutes to remove residual cells that could contaminate cfDNA.

3. **DNA extraction** — Extract cell-free DNA using a kit designed for cfDNA recovery (low input, short fragments). Standard genomic DNA kits are not suitable as they lose short cfDNA fragments.

4. **Quality control** — Verify cfDNA profile on Bioanalyzer. Successful extraction shows characteristic cfDNA peak at approximately 160-180 bp (mononucleosomal fragment).

5. **Sequencing** — Prepare Oxford Nanopore library with methylation-aware protocol. Sequence on MinION, GridION, or PromethION with R10.4.1 flow cells.

6. **Data processing** — Run basecalling with methylation model enabled. Process with NanoporeToBED-Pipeline to generate CpG BED files suitable for MethylSense.

7. **Diagnosis** — Run MethylSense_predict.R with your trained model. The script outputs a diagnostic classification with confidence score for each sample.

8. **Interpretation** — Review prediction report. High confidence predictions (>80%) can inform clinical decisions. Medium or low confidence results warrant additional confirmatory testing.

### Interpreting Diagnostic Results

The prediction output includes a probability score for each sample indicating how strongly the methylation profile matches each class. The confidence score represents the maximum predicted probability across all classes:

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-top: 2px solid #333; border-bottom: 1px solid #333;">
<th style="padding: 8px; text-align: left;">Confidence Level</th>
<th style="padding: 8px; text-align: left;">Probability Range</th>
<th style="padding: 8px; text-align: left;">Clinical Interpretation</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">High</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">>80%</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Strong methylation signal consistent with predicted class. Prediction is reliable and can inform clinical decisions.</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">Medium</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">60-80%</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">Moderate methylation signal. Consider confirmatory testing (culture, PCR, endoscopy, imaging) before treatment decisions.</td></tr>
<tr style="border-bottom: 2px solid #333;"><td style="padding: 8px;">Low</td><td style="padding: 8px;"><60%</td><td style="padding: 8px;">Weak or ambiguous methylation signal. Do not use for clinical decisions. Investigate sample quality or consider resequencing.</td></tr>
</tbody>
</table>


### Causes of Low Confidence Predictions

Low confidence predictions warrant investigation before interpretation:

- **Sample quality issues** — Degraded cfDNA, insufficient DNA input, or contamination with genomic DNA from lysed blood cells
- **Technical factors** — Low sequencing coverage, library preparation failure, or basecalling errors
- **Biological factors** — Early-stage infection before host response develops, atypical disease presentation, or immunocompromised patients with altered epigenetic responses
- **Species mismatch** — Sample from a species different to the training data population

### Pre-trained Models

Pre-trained diagnostic models are available for download:

**Avian Aspergillosis (Chicken, n=124)**
- Accuracy: >95% with nested cross-validation
- Download: [Zenodo Repository](https://zenodo.org/records/15194046)
- DOI: [10.5281/zenodo.15194046](https://doi.org/10.5281/zenodo.15194046)
- Includes: Trained models, DMR coordinates, and validation results

**Additional species and disease models:** Contact for availability at markus.drag@sund.ku.dk

### Regulatory Notice

MethylSense is a research tool and has not been approved by FDA, EMA, or other regulatory bodies for clinical diagnostic use. Diagnostic decisions should be made in conjunction with other clinical findings and established diagnostic methods.

Models trained on one species may not generalise to other species without proper validation. Cross-species application requires separate validation studies before clinical use.

---

## Troubleshooting

### Common Issues

Package not found
```bash
Rscript MethylSense_installer.R
```

Not enough samples per group
```bash
Rscript MethylSense_analysis.R ... --allow_small_groups
```
Small groups (<10 samples) may produce unreliable or overfitted results. We recommend a minimum of 20 samples per group for robust classifier training.

Memory error
```bash
# Use larger windows (fewer features):
--window_file windows_10kb.bed

# Use fewer models:
--models rf,svm
```

Training too slow
```bash
--models fast
--cv_repeats 3
```

### Performance Benchmarks

<table style="border-collapse: collapse; width: 100%;">
<thead>
<tr style="border-top: 2px solid #333; border-bottom: 1px solid #333;">
<th style="padding: 8px; text-align: left;">Step</th>
<th style="padding: 8px; text-align: left;">Runtime</th>
<th style="padding: 8px; text-align: left;">Memory</th>
</tr>
</thead>
<tbody>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">load_data</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">5-10 min</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">2-4 GB</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">general_overview</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">5-15 min</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">2-6 GB</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">train (no CV)</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">20-30 min</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">4-8 GB</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">train (10× CV)</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">60-90 min</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">4-8 GB</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">train (nested CV)</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">2-3 hours</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">4-8 GB</td></tr>
<tr><td style="padding: 8px; border-bottom: 1px solid #ddd;">predict</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">1-2 min</td><td style="padding: 8px; border-bottom: 1px solid #ddd;">1-2 GB</td></tr>
<tr style="border-bottom: 2px solid #333;"><td style="padding: 8px;">reviewer</td><td style="padding: 8px;">5-10 min</td><td style="padding: 8px;">2-4 GB</td></tr>
</tbody>
</table>


Benchmarks for 30-50 samples, 2 groups, 8 CPU cores.

---

## Citation

If you use MethylSense in your research, please cite:

> Drag MH, Hvilsom C, Poulsen LL, Jensen HE, Tahas SA, Leineweber C, Cray C, Bertelsen MF, Bojesen AM (2025). New high accuracy diagnostics for avian Aspergillus fumigatus infection using Nanopore methylation sequencing of host cell-free DNA and machine learning prediction. bioRxiv 2025.04.11.648151. https://doi.org/10.1101/2025.04.11.648151

---

## Acknowledgements

MethylSense builds upon excellent work from the R and Bioconductor communities.

### Core Dependencies

methylKit (Akalin et al., 2012) — The foundation for DMR detection and methylation analysis in MethylSense. This Bioconductor package provides comprehensive tools for genome-wide DNA methylation analysis from bisulfite sequencing data.

> Akalin A, Kormaksson M, Li S, Garrett-Bakelman FE, Figueroa ME, Melnick A, Mason CE (2012). methylKit: a comprehensive R package for the analysis of genome-wide DNA methylation profiles. Genome Biology, 13, R87. https://doi.org/10.1186/gb-2012-13-10-r87

caret (Kuhn, 2008) — The machine learning framework powering model training, cross-validation, and hyperparameter tuning in MethylSense.

> Kuhn M (2008). Building Predictive Models in R Using the caret Package. Journal of Statistical Software, 28(5), 1-26. https://doi.org/10.18637/jss.v028.i05

pROC (Robin et al., 2011) — ROC curve analysis and AUC calculations.

> Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez JC, Müller M (2011). pROC: an open-source package for R and S+ to analyze and compare ROC curves. BMC Bioinformatics, 12, 77. https://doi.org/10.1186/1471-2105-12-77

### Additional R Packages

MethylSense also uses: GenomicRanges (Lawrence et al., 2013), ggplot2 (Wickham, 2016), randomForest (Liaw & Wiener, 2002), xgboost (Chen & Guestrin, 2016), glmnet (Friedman et al., 2010), e1071, ranger, pheatmap, and RColorBrewer.

### References

<details>
<summary>Expand full reference list</summary>

DNA Methylation
- Bird A (2002). DNA methylation patterns and epigenetic memory. Genes & Development, 16(1), 6-21.
- Jones PA (2012). Functions of DNA methylation: islands, start sites, gene bodies and beyond. Nature Reviews Genetics, 13(7), 484-492.

Nanopore Methylation
- Simpson JT, Workman RE, Zuzarte PC, David M, Dursi LJ, Timp W (2017). Detecting DNA cytosine methylation using nanopore sequencing. Nature Methods, 14(4), 407-410.

Statistical Methods
- Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society: Series B, 57(1), 289-300.
- Varma S, Simon R (2006). Bias in error estimation when using cross-validation for model selection. BMC Bioinformatics, 7, 91.

Related Pipelines
- [NanoporeToBED-Pipeline](https://github.com/markusdrag/NanoporeToBED-Pipeline)
- [GenomeToWindows](https://github.com/markusdrag/GenomeToWindows)

</details>

---

## License and Authors

License: Academic Free License 3.0 (AFL-3.0). See [LICENSE](LICENSE).

Lead Developer: Markus Hodal Drag

Co-authors: Christina Hvilsom, Louise Ladefoged Poulsen, Henrik Elvang Jensen, Stamatios Alan Tahas, Christoph Leineweber, Carolyn Cray, Mads Frost Bertelsen, Anders Miki Bojesen

Contact: markus.drag@sund.ku.dk | [GitHub Issues](https://github.com/markusdrag/MethylSense/issues)

---

<div align="center">

Developed for Oxford Nanopore methylation analysis

</div>
