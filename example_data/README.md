# MethylSense Example Data

This directory contains example Nanopore methylation data for testing and demonstrating the MethylSense pipeline.

## Dataset Overview

- **15 samples total** (Gallus gallus / chicken)
  - 5 control samples (control_1 to control_5)
  - 5 suspected infection samples (suspected_1 to suspected_5)
  - 5 confirmed infection samples (infected_1 to infected_5)

- **Progressive methylation patterns** reflecting infection progression:
  - Control: Baseline methylation (~25%)
  - Suspected: Moderate DMR effects (±12% in 100 DMR regions)
  - Infected: Strong DMR effects (±30% in 100 DMR regions)

- **Data characteristics**:
  - ~50,000 CpG sites per sample
  - Mean coverage: 20x
  - 100 true DMR regions (saved in true_dmrs.csv)
  - Chromosomes 1-5 only (reduced dataset for quick testing)

## Files

### BED Files (18-column Nanopore format)
- `control_*.CpG.bed` - Control samples (n=5)
- `suspected_*.CpG.bed` - Suspected infection samples (n=5)
- `infected_*.CpG.bed` - Confirmed infection samples (n=5)

### Metadata
- `sample_metadata.xlsx` - Sample sheet for MethylSense (Excel format)
- `sample_metadata.csv` - Sample sheet (CSV format)
- `true_dmrs.csv` - True DMR coordinates used for simulation

## Quick Start

```bash
# 1. Load and preprocess data
Rscript MethylSense_load_data.R \
  --species "Gallus_gallus" \
  --sample_sheet example_data/sample_metadata.xlsx \
  --bed_dir ./example_data \
  --output_dir ./example_output

# 2. Train models
Rscript MethylSense_analysis.R \
  --qs_file ./example_output/*_methylRaw.qs \
  --output_dir ./example_training \
  --window_file windows_5kb.bed \
  --models rf,svm,xgboost \
  --group_names Control,Suspected,Infected \
  --positive_class Infected \
  --cv_repeats 5

# 3. Evaluate models
Rscript MethylSense_reviewer.R \
  --model_dir ./example_training/models/rf_model/ \
  --qs_file ./example_output/*_methylRaw.qs \
  --sample_sheet example_data/sample_metadata.xlsx \
  --output_dir ./example_report
```

## Expected Results

With this example dataset, you should observe:

- **DMR Detection**: ~100 significant DMRs (FDR < 0.05)
- **Classification Accuracy**: >90% with cross-validation
- **Clear separation** between Control, Suspected, and Infected groups in PCA plots
- **Progressive methylation changes** in DMR heatmaps

## Data Generation

This synthetic dataset was generated using `generate_example_beds.R` to simulate realistic Nanopore methylation data with known DMRs. The progressive methylation patterns reflect observations from real avian aspergillosis infection studies.

To regenerate this dataset:

```bash
Rscript generate_example_beds.R
```
