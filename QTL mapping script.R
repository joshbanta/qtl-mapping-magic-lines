# -------------------------------
# Optional Setup: Install Required Packages
# -------------------------------

# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
library(devtools)

# Install atMAGIC and qtl2helper from GitHub if not already installed
if (!requireNamespace("atMAGIC", quietly = TRUE)) {
  devtools::install_github("tavareshugo/atMAGIC")
}
if (!requireNamespace("qtl2helper", quietly = TRUE)) {
  devtools::install_github("tavareshugo/qtl2helper")
}

# Load core packages
library(atMAGIC)        # Provides MAGIC genotype data (kover2009)
library(qtl2helper)     # Utilities for R/qtl2 workflows

# -------------------------------
# Step 1: Load Genotype Data
# -------------------------------

data("kover2009")        # Load pre-formatted genotype data object

# -------------------------------
# Step 2: Read and Clean Phenotype Data
# -------------------------------

# Read root trait phenotype data
pheno_raw <- read.table("raiz(mean).txt", header = TRUE)

# Rename column containing line IDs
colnames(pheno_raw)[which(colnames(pheno_raw) == "LM")] <- "SUBJECT.NAME"

# Format SUBJECT.NAME with MAGIC prefix for merging
pheno_raw$SUBJECT.NAME <- gsub("LM", "MAGIC.", pheno_raw$SUBJECT.NAME)

# Drop unnecessary columns (e.g., "IND")
pheno_clean <- pheno_raw[, -which(colnames(pheno_raw) == "IND")]

# Reformat ID column: remove "MAGIC." prefix and convert to numeric
pheno_clean[,1] <- as.numeric(gsub("MAGIC\\.", "", pheno_clean[,1]))

# -------------------------------
# Step 3: Merge with Line Label Reference
# -------------------------------

line_labels <- read.csv("HSRIL to MAGIC.csv")
colnames(pheno_clean)[1] <- "ours"

# Merge phenotype data with label reference
pheno_merged <- merge(pheno_clean, line_labels, by = "ours", all.x = TRUE)

# Reformat SUBJECT.NAME after merging
colnames(pheno_merged)[1] <- "SUBJECT.NAME"
pheno_merged$SUBJECT.NAME <- paste0("MAGIC.", pheno_merged$SUBJECT.NAME)

# Subset to keep phenotype columns only
pheno_final <- pheno_merged[, 1:ncol(pheno_clean)]

# -------------------------------
# Step 4: Attach Phenotype Data to Genotype Object
# -------------------------------

kover2009.updated <- add_pheno(kover2009, pheno_final, idcol = "SUBJECT.NAME")

# -------------------------------
# Step 5: Calculate Genotype Probabilities
# -------------------------------

kover2009_probs <- calc_genoprob(kover2009.updated)

# -------------------------------
# Step 6: Perform Genome-wide QTL Scan
# -------------------------------

kover2009_scan <- scan1(kover2009_probs, kover2009.updated$pheno)

# -------------------------------
# Step 7: Permutation Test for Significance Thresholds
# -------------------------------

kover2009_perm <- scan1perm(kover2009_probs, kover2009.updated$pheno,
                            n_perm = 10000, cores = 2)

kover2009_threshold <- summary(kover2009_perm)

# -------------------------------
# Step 8: Plot Genome Scans with Threshold Lines
# -------------------------------

traits <- c("RP", "Num.RL", "RL", "Num.RA", "RA")

for (trait in traits) {
  plot(kover2009_scan, kover2009.updated$pmap, lodcolumn = trait)
  abline(h = kover2009_threshold[, trait], col = "red")
}

# -------------------------------
# Step 9: Identify and Save Significant Peaks
# -------------------------------

trait_col <- which(colnames(pheno_final[, -1]) == "Num.RA")
threshold_val <- kover2009_threshold[, trait_col]

peaks_Num.RA <- find_peaks(kover2009_scan, kover2009.updated$pmap,
                           threshold = threshold_val)

peaks_Num.RA_filtered <- peaks_Num.RA[peaks_Num.RA$lodcolumn == "Num.RA", ]

write.csv(peaks_Num.RA_filtered, "peaks_Num_RA.csv", row.names = FALSE)

# -------------------------------
# Step 10: Compute and Save LOD Interval for Example QTL
# -------------------------------

lod_interval_Num.RA <- lod_int(kover2009_scan, kover2009.updated$pmap,
                               chr = 1, lodcolumn = 4)

write.csv(lod_interval_Num.RA, "lod_interval_Num_RA.csv", row.names = FALSE)

# -------------------------------
# Step 11: Retrieve Gene Coordinates and Overlay on LOD Interval
# -------------------------------

# Install Bioconductor annotation packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("TxDb.Athaliana.BioMart.plantsmart51")

# Load Arabidopsis gene annotation
library(TxDb.Athaliana.BioMart.plantsmart12)
txdb <- TxDb.Athaliana.BioMart.plantsmart12

transcripts <- transcriptsBy(txdb, by = "gene")
gene_ranges <- reduce(transcripts)
gene_df <- as.data.frame(gene_ranges)
rownames(gene_df) <- gene_df[,1]

# Install genomic tools if needed
if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
  BiocManager::install("GenomicRanges")
}
if (!requireNamespace("IRanges", quietly = TRUE)) {
  BiocManager::install("IRanges")
}
library(GenomicRanges)
library(IRanges)

# Create GRanges for LOD interval
lod_chr <- as.character(unique(lod_interval_Num.RA$chr))
lod_gr <- GRanges(
  seqnames = lod_chr,
  ranges = IRanges(start = min(lod_interval_Num.RA$start),
                   end = max(lod_interval_Num.RA$end))
)

# Create GRanges for all genes
gene_gr <- GRanges(
  seqnames = gene_df$seqnames,
  ranges = IRanges(start = gene_df$start, end = gene_df$end),
  gene_id = rownames(gene_df)
)

# Find and save genes overlapping the QTL interval
overlap_hits <- findOverlaps(gene_gr, lod_gr)
candidate_genes_in_interval <- gene_df[queryHits(overlap_hits), ]

write.csv(candidate_genes_in_interval, "candidate_genes_in_Num_RA_interval.csv", row.names = FALSE)
