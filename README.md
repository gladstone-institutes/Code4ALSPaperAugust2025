# Cell death analysis

This repository contains R functions and analysis scripts for changes in cell-death over time analysis.

## Authors
- Julia Kaye
- Reuben Thomas
- Stephanie Lam

## Repository Structure

```
├── functions/                          # Core R functions
│   ├── MRID_Logodd_Functions_RT_July1.R   # Log odds statistical functions
│   └── MRID_Ratio_Functions.R              # MRID ratio calculation functions
├── analysis/                           # Analysis scripts and notebooks
│   ├── combined_c9orf72_cellline/         # C9orf72 cell line comparisons
│   ├── combined_c9orf72_class/            # C9orf72 class comparisons
│   └── sals_analysis/                     # Sporadic ALS analysis
│       ├── sals_all_analysis.Rmd
│       ├── sals_high_islet1_analysis.Rmd
│       └── sals_low_islet1_analysis.Rmd
└── docs/                              # Generated HTML reports
    ├── combined_c9orf72_cellline_analysis.html
    ├── combined_c9orf72_class_analysis.html
    ├── sals_all_analysis.html
    ├── sals_high_islet1_analysis.html
    └── sals_low_islet1_analysis.html
```

## Description

This project provides tools for analyzing motor neuron death ratios in ALS research using MRID methodology. The analysis includes:

- Log odds statistical analysis for cell death ratios
- Comparison between control and ALS conditions
- Cell line and class-based comparisons
- Sporadic ALS analysis with different islet1 expression levels

## Usage

1. Load the required functions from the `functions/` directory
2. Use the R Markdown files in `analysis/` to run specific analyses
3. Generated HTML reports are saved in the `docs/` directory

## Requirements

- R with the following packages:
  - ggplot2
  - dplyr
  - xtable
  - lme4
  - plyr
  - gridExtra

## Input Data Structure

The analysis expects the following data structure:
```
Input path - Experiment1 - Experiment1_ratio_output.csv
                         - Experiment1_timepoint.csv
           - Experiment2 - Experiment2_ratio_output.csv
                         - Experiment2_timepoint.csv
```
