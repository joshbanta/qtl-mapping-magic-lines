# QTL Mapping in Arabidopsis MAGIC Lines

This repository contains R code and associated data files for mapping quantitative trait loci (QTL) related to root traits in Arabidopsis MAGIC lines using the [`qtl2`](https://kbroman.org/qtl2/) and [`atMAGIC`](https://github.com/tavareshugo/atMAGIC) packages.

## Overview

This project performs genome-wide QTL scans on five root traits across 139 MAGIC lines. It uses custom phenotype data and genotype probabilities provided in the `atMAGIC` R package.

The workflow includes:
- Formatting and merging phenotype data with MAGIC line identifiers
- Calculating genotype probabilities
- Running `scan1` for QTL detection
- Performing 10,000 permutations for empirical threshold estimation
- Identifying significant QTL peaks
- Extracting LOD intervals
- (Optional) Overlaying gene coordinates for candidate identification

## Repository Contents

- `qtl_mapping.R` – Main R script for executing the QTL workflow
- `raiz(mean).txt` – Raw phenotype data for five root traits
- `HSRIL to MAGIC.csv` – Mapping file for line ID conversion
- `57gene_cooridinate.csv` – List of candidate genes with genomic coordinates

## Usage Instructions

1. Clone the repository or download the ZIP.
2. Open R or RStudio and set the working directory to the repository folder.
3. Run the script:
   ```r
   source("qtl_mapping.R")
