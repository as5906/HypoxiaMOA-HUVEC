# TF Coverage and Heatmap Analysis

## Overview

This directory contains the necessary input files, scripts, and commands to generate coverage plots and heatmaps of transcription factor (TF) occupancy relative to ChIP-seq summits and motif sites. These analyses support **Figure 3A–E** and **Figure S2**, highlighting baseline MOA-seq signal patterns and comparative analysis with publicly available ChIP-seq datasets.

---

## Contents

### Input Files

| Filename                                           | Description |
|----------------------------------------------------|-------------|
| `FLI1_motif_sites.bed`, `ETS1_motif_sites.bed`, etc. | BED files with motif-predicted binding sites for select TFs, provided in ChIP directory |
| `MOA0_merge_frenter_q20_chr_sort.bed`             | Baseline MOA-seq peak data (0h) |
| `MOA0_merge_frenter_q20_arcdn.bw`                 | BigWig file for MOA-seq 0h read density |

---

### Script

| Filename | Description |
|----------|-------------|
| `Script` | Bash script to generate heatmaps and coverage plots for Figures 3A–E. Includes summit overlap, per-TF signal computation, and ChIP-seq vs MOA-seq comparisons using `Cov_avg.py`. |
| `Cov_avg.py` | Python script to compute average coverage across motifs and calculate enrichment ratios. Used for generating CPM plots and peak-to-flank enrichment values. |

---

## Description of Script Functionality

### Figure 3A & 3B – Global Heatmap of MOA-seq Signal at ChIP Summits

1. **Download ReMap Peaks**  
   Downloads and decompresses `remap2022_HUVEC-C_nr_macs2_hg38_v1_0.bed`.

2. **Summit Extraction & Filtering**  
   - Extracts summit coordinates using columns 7–8.
   - Filters summits overlapping MOA-seq peaks using `bedtools intersect`.

3. **Matrix & Heatmap Generation**  
   - Computes read matrix using `computeMatrix reference-point`.
   - Generates heatmap using `plotHeatmap`.

---

### Figure 3C & Figure S2 – Per-TF Coverage Analysis

1. **Parse Summit BED with TF Names**  
   Extracts 4-column BED with TF names for TF-specific slicing.

2. **For Each TF (e.g., FLI1, ETS1, BRD4):**
   - Extracts TF-specific binding sites to separate BED files.
   - Runs `Cov_avg.py` using MOA-seq BED input, ±1000 bp window.
   - Stores results in TF-specific directories.

---

### Figure 3D & 3E – Direct ChIP vs MOA-seq Comparison at TF Motifs

1. **Download and Align ChIP-seq FASTQs**
   - Performs alignment using `bwa-mem2` to hg38.
   - Sorts, filters, and converts to BED format.

2. **Coverage Plot with Cov_avg.py**
   - Runs `Cov_avg.py` separately on ChIP BED and MOA BED.
   - Inputs are TF-specific motif BEDs (e.g., `FLI1_motif_sites.bed`).

---

## Output Files

| Filename / Folder                  | Description |
|------------------------------------|-------------|
| `heat.png`                         | Heatmap of MOA-seq signal at ReMap ChIP summits |
| `<TF>/Feature_cpm.png`             | CPM-normalized signal plot around TF binding sites |
| `<TF>_chip/Feature_cpm.png`        | Signal at motif sites using ChIP signal |
| `<TF>_moa/Feature_cpm.png`         | Signal at motif sites using MOA-seq signal |
| `<TF>/position_counts.csv`         | Raw position-level counts for TF binding sites |
| `<TF>/feature_enrichment_scores.csv` | Peak-to-flank enrichment ratios |

---

## Dependencies

### General (for `Script`)

- `bedtools`
- `awk`
- `bwa-mem2`
- `samtools`
- [`deepTools`](https://deeptools.readthedocs.io/en/develop/) (`computeMatrix`, `plotHeatmap`)
- `wget`
- `gzip`
- Unix shell (bash)

### Python (`Cov_avg.py`)

Install via pip:

```bash
pip install numpy pandas matplotlib plotly
```

| Python Module | Purpose |
|---------------|---------|
| `numpy`       | Array and matrix calculations |
| `pandas`      | Read/write and manipulate BED and signal data |
| `matplotlib`  | Generates PNG plots of signal around motifs |
| `plotly`      | Interactive visualization (optional, used if enabled) |

---

## Cov_avg.py – Summary

This Python script computes average read coverage ±N base pairs from feature centers and optionally generates normalized enrichment scores.

### Syntax

```bash
python Cov_avg.py <signal_BED> <feature_BED> <flank_bp> <strand_aware>
```

### Example

```bash
python Cov_avg.py MOA0_merge_frenter_q20_chr_sort.bed FLI1_motif_sites.bed 1000 False
```

### Output Files

- `Feature_cpm.png`: Normalized CPM signal around features
- `position_counts.csv`: Raw signal per position (relative to feature center)
- `feature_enrichment_scores.csv`: Peak-to-flank enrichment ratio

---

## Notes

- Enrichment is calculated as the signal at the feature center divided by the average signal in flanking regions.
- Output plots provide visual confirmation of TF-specific signal patterns in MOA-seq or ChIP-seq data.
- The analysis demonstrates how chromatin occupancy by key transcription factors aligns with MOA-seq signal at baseline.

---
