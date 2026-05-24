# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Single-cell RNA sequencing (scRNA-seq) analysis pipeline for rice (Oryza sativa) seedlings. Generates a single-cell transcriptome atlas of leaf and root tissues.

## Project Structure

This is a sequential 4-step R pipeline:

1. **`1-Single-cell pepline.R`** - Preprocessing: reads raw data, removes protoplasting-sensitive genes, performs SCTransform normalization, integrates samples, runs PCA/UMAP, clusters cells (29 clusters)
2. **`2-Assignment of clusters to cell types by MICI.R`** - Cell type assignment using marker genes from `Table S5`
3. **`3-Determination of the organ origin in the Seedling samples.R`** - Assigns cells as "Aerial" or "Root" origin using differential expression
4. **`4-Reconstruction of developmental trajectory by monocle.R`** - Trajectory inference using Monocle (DDRTree)

## Running the Analysis

Execute scripts sequentially in R/RStudio:

```bash
Rscript "1-Single-cell pepline.R"
Rscript "2-Assignment of clusters to cell types by MICI.R"
Rscript "3-Determination of the organ origin in the Seedling samples.R"
Rscript "4-Reconstruction of developmental trajectory by monocle.R"
```

## Dependencies

### Core Libraries
- **Seurat** - Single-cell analysis framework (primary)
- **Monocle** - Trajectory inference
- **data.table** - Data I/O
- **tidyverse/dplyr** - Data manipulation
- **ggplot2/patchwork/ggpubr** - Visualization
- **Matrix** - Sparse matrix operations
- **parallel** - Multi-core processing (uses 28 cores in script 3)

## Input Files Required

- `Gene.length` - Gene length information
- `HTSeq_out_rep1_AP.txt`, `HTSeq_out_rep1_BP.txt` - Bulk RNA-seq counts (replicate 1)
- `HTSeq_out_rep2_AP.txt`, `HTSeq_out_rep2_BP.txt` - Bulk RNA-seq counts (replicate 2)
- `all.csv` - Single-cell expression matrix
- `Table S5` - Marker gene table for cell type assignment

## Key Parameters

- **Clustering resolution:** 0.75 (produces 29 clusters)
- **PCA dimensions:** 100 PCs
- **Protoplasting filter:** Fold-change cutoff of 3
- **Trajectory subsample:** 10,000 cells for Monocle analysis

## Data Structure

- **Samples:** Aerial (S_A1, S_A2) and Root (S_R1, S_R2)
- **Stress conditions:** HS (heat stress), ID (iron deficiency), LN (low nitrogen), Ctrl (control)
- **Cell types:** Assigned via MICI scoring system using marker genes
