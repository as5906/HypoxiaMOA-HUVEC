# MOA-seq Enrichment - TSS, Chromatin Accessibility, and cCRE Enrichment Analysis

This directory contains input files, scripts, and output files used to generate Figure 2, which illustrates MOA-seq signal distribution across transcription start sites (TSS), ENCODE cCREs, and DNase I hypersensitive sites (DHSs) across multiple cell types.

---

## Overview

This analysis includes:

- Signal profiling over TSS stratified by RNA-seq expression quartiles (Figure 2A–B)
- MOA-seq enrichment at cis-regulatory elements such as dELS, pELS, K4me3-marked promoters, and CTCF-bound regions (Figure 2C–D)
- Distribution of MOA-seq peak overlap with ENCODE DHSs across six cell types (Figure 2E–F)
- Empirical enrichment testing using base pair coverage distribution (`dist.py`)

---

## Contents

### Input Files

- `TSS_Q1/`, `TSS_Q2/`, `TSS_Q3/`, `TSS_Q4/`  
  These directories contain genes binned by RNA-seq expression quartiles.

- `UCSC_ENCODE_GM12878_PEAKS.bed`  
  `UCSC_ENCODE_HUVEC_PEAKS.bed`  
  `UCSC_ENCODE_HepG2_PEAKS.bed`  
  `UCSC_ENCODE_K562_PEAKS.bed`  
  `UCSC_ENCODE_NHEK_PEAKS.bed`  
  `UCSC_ENCODE_NHLF_PEAKS.bed`  
  DNase I hypersensitive peaks from ENCODE for six human cell lines.

- `hg38.chrom.sizes`  
  Genome size file used for sorting, shuffling, and plotting.

---

## Scripts

### `dist.py`

Generates a probability distribution of base pair overlaps between MOA-seq peaks and annotation features. Computes mean, standard deviation, and Z-score for observed overlap compared to a shuffled background distribution.

**Usage:**
```bash
python dist.py <input_file> <label_hour> <feature_name> <expected_bp>
```

**Arguments:**
- `<input_file>`: text file with one-column list of shuffled overlap base pair counts  
- `<label_hour>`: label for condition (e.g., `1hr`, `3hr`)  
- `<feature_name>`: feature being analyzed (e.g., `pELS`)  
- `<expected_bp>`: observed overlap base pairs for comparison  

**Output:**
- PNG file with histogram and KDE curve  
- Vertical lines for expected (observed) and average values  

**Console output:**
- Standard deviation  
- Mean  
- Z-score of observed relative to null  

**Example:**
```bash
python dist.py shuffled_lengths.txt 3hr pELS 147
```
This creates `3hr.pELS.png` and prints enrichment statistics.

---

### `cov_avg.py`

Computes average signal coverage over a set of intervals using MOA-seq signal files. This is typically used to quantify enrichment of MOA-seq reads over TSS or cCRE elements.

Note: See inline comments in the script for argument formatting and file requirements.

---

## Workflow Summary

### TSS Enrichment (Figure 2A–B)

- TSSs grouped into quartiles using RNA-seq data  
- Signal profile plots generated using tools like `computeMatrix` and `plotProfile`  

### cCRE Signal Overlap (Figure 2C–D)

- Intersect MOA-seq peaks with cCRE features (dELS, pELS, K4m3, CTCF)  
- Compute overlap statistics and visualize average enrichment  

### DNase Site Distribution (Figure 2E–F)

- Use `bedtools intersect` to find MOA-seq overlap with DHSs across six ENCODE cell lines  
- Summarize and visualize relative enrichment across cell types  

### Empirical Enrichment Analysis

- Use `bedtools shuffle` to generate null distributions of base pair overlaps  
- Analyze with `dist.py` to assess significance of observed enrichment  

---

## Dependencies

- Python 3.x with:
  - `pandas`, `numpy`, `matplotlib`, `seaborn`
- `bedtools`
- UCSC utilities (optional): `bigWigAverageOverBed`, `bigBedToBed`
- `deepTools` (optional for heatmap/profile plots)

---
