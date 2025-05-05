# QTL Mapping in Arabidopsis MAGIC Lines

This repository contains R code and supporting data for performing QTL mapping of root traits in Arabidopsis MAGIC lines using the `qtl2` and `atMAGIC` packages. The workflow includes genotype-phenotype integration, genome-wide scans, permutation-based thresholding, peak detection, and gene annotation overlay based on Arabidopsis gene coordinates.

## Overview

The analysis focuses on identifying QTL for five root traits measured in 139 MAGIC lines. The final output includes significant QTL peaks and a list of genes overlapping the QTL interval.

### Root Traits Analyzed:
- Root Penetration (RP)
- Number of Right Lateral roots (Num.RL)
- Length of Right Lateral roots (RL)
- Number of Right Adventitious roots (Num.RA)
- Length of Right Adventitious roots (RA)

## Files Included

- `qtl_mapping_final.R`: The full R script for data preparation, QTL analysis, and gene overlay.
- `raiz(mean).txt`: Phenotypic measurements for five root traits.
- `HSRIL to MAGIC.csv`: Mapping file linking experimental line IDs to MAGIC line names.

## Required R Packages

- `qtl2`
- `devtools`
- `atMAGIC` (installed via GitHub)
- `qtl2helper` (installed via GitHub)
- `BiocManager`
- `TxDb.Athaliana.BioMart.plantsmart12`
- `GenomicRanges`
- `IRanges`

The script includes install commands for all dependencies.

## How to Use

1. Clone or download this repository.
2. Open `qtl_mapping_final.R` in R or RStudio.
3. Place the input files in your working directory.
4. Run the script to execute the full workflow.

The script will:
- Load and clean phenotype and genotype data
- Perform genome-wide QTL scans
- Determine significance thresholds by permutation
- Identify QTL peaks and LOD intervals
- Overlay genes within QTL regions using Arabidopsis gene annotations

## Output Files

- `peaks_Num_RA.csv`: Significant QTL peaks for the Num.RA trait
- `lod_interval_Num_RA.csv`: LOD interval surrounding the peak QTL
- `candidate_genes_in_Num_RA_interval.csv`: List of genes overlapping the QTL interval

## Citation and License

This project is made available under the [MIT License](LICENSE). Please cite this repository if using it in your own research or derivative work.
