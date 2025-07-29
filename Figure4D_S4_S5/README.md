# ENCODE cCRE Intersection with Differential MOA-seq Peaks

## Overview

This directory contains peak files and scripts used to analyze the overlap between differential MOA-seq peaks and ENCODE candidate cis-regulatory elements (cCREs). The goal is to quantify how many MOA-seq peaks intersect known regulatory element classes and to determine the number of overlapping base pairs and peak counts for each class.
This directory supports Figure4D (Density calculation can be found in Supplementary Table 6), FigureS4, and FigureS5.

---

## Contents

### Peak Files

| Filename                           | Description                                                 |
|------------------------------------|-------------------------------------------------------------|
| `MOA0_P001_IDR_Nt_all_brr_crr.bed` | Baseline IDR peaks (normoxia) used as within-sample control |
| `MOA1_gain_diff.bed`               | Gain-of-signal peaks at 1h hypoxia                          |
| `MOA1_loss_diff.bed`               | Loss-of-signal peaks at 1h hypoxia                          |
| `MOA1_shared_diff.bed`            | Shared (unchanged) peaks at 1h hypoxia                      |
| `MOA3_gain_diff.bed`               | Gain-of-signal peaks at 3h hypoxia                          |
| `MOA3_loss_diff.bed`               | Loss-of-signal peaks at 3h hypoxia                          |
| `MOA3_shared_diff.bed`            | Shared (unchanged) peaks at 3h hypoxia                      |
| `MOA24_gain_diff.bed`              | Gain-of-signal peaks at 24h hypoxia                         |
| `MOA24_loss_diff.bed`              | Loss-of-signal peaks at 24h hypoxia                         |
| `MOA24_shared_diff.bed`           | Shared (unchanged) peaks at 24h hypoxia                     |

---

### Script

| Filename | Description |
|----------|-------------|
| `Script` | This script performs all intersection analyses between MOA-seq differential peak sets and ENCODE cCRE classes, recording both base pair and peak-level overlaps. |

---

## Description of Script Functionality

1. **Download and Convert ENCODE Regulatory Element Data**  
   - Downloads the `encodeCcreCombined.bb` file from UCSC (hg38).
   - Converts it to BED format using `bigBedToBed`.
   - Extracts subsets of cCREs into separate BED files:
     - `encodeCcre_dELS.bed`
     - `encodeCcre_pELS.bed`
     - `encodeCcre_PLS.bed`
     - `encodeCcre_CTCF.bed`
     - `encodeCcre_K4m3.bed`

2. **Intersect Each MOA Peak Set with Each cCRE Class**  
   For each cCRE type:
   - Computes number of overlapping peaks with each MOA-seq BED file.
   - Computes total number of overlapping base pairs.
   - Stores results in the following output files:
     - `peak_intersection_counts.txt`
     - `bp_intersection_counts.txt`

---

## Output Files

| Filename                     | Description                                                       |
|------------------------------|-------------------------------------------------------------------|
| `peak_intersection_counts.txt` | Line count of MOA-seq peaks overlapping each cCRE class          |
| `bp_intersection_counts.txt`   | Total number of base pairs overlapping cCRE classes for each set |

---

## Dependencies

- `wget`
- `bedtools`
- `awk`
- `bigBedToBed` (from UCSC Tools: https://hgdownload.soe.ucsc.edu/admin/exe/)

---

## Notes

- Each intersection is performed using both `-u` (unique peak count) and default mode (base pair overlap).
- This analysis helps determine whether specific classes of regulatory elements are more likely to gain or lose accessibility during hypoxia.
- Results support downstream visualization and interpretation of functional regulatory element enrichment in hypoxia-induced chromatin remodeling.
