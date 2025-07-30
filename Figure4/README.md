# MOA-seq Differential Peak Calling and Annotation Pipeline

This script automates the identification and comparison of differential MOA-seq peaks across hypoxia time points (1h, 3h, 24h) using `MACS3`, followed by sorting, visualization, and genomic annotation for enrichment analysis.

## Overview

This workflow includes:

- Differential peak calling between hypoxia and normoxia (0h) at each time point
- Heatmap generation for gain/loss peaks (Figure S3)
- Intersection with genomic annotations and cis-regulatory elements (Figure 4A–C)

## Contents

### Input
- MOA-seq BED files for each time point:
  - `MOA0_merge_frenter_q20_chr_sort.bed` (0h)
  - `MOA1_merge_frenter_q20_chr_sort.bed` (1h)
  - `MOA3_merge_frenter_q20_chr_sort.bed` (3h)
  - `MOA24_merge_frenter_q20_chr_sort.bed` (24h)

### Output
- `*_gain_diff.bed`, `*_loss_diff.bed`, `*_shared_diff.bed` for each time point
- Sorted BED files for plotting
- Heatmaps: `1.svg`, `3.svg`, `24.svg`
- BED files of differential peaks overlapping protein-coding genes or ENCODE cCREs:
  - `AlldiffMOA_ccre.bed`, `AlldiffMOA_nonccre.bed`
  - `AlldiffMOA_ccre_GENES.bed`, `AlldiffMOA_nonccre_GENES.bed`

## Workflow Summary

### 1. Differential Peak Calling (Figure 4A)

For each hypoxia time point:
- `macs3 callpeak`: generates signal tracks (`treat_pileup.bdg`, `control_lambda.bdg`)
- `macs3 bdgdiff`: compares to 0h to identify gain/loss peaks

Each comparison is placed in a separate directory (`1h/`, `3h/`, `24h/`).

### 2. Sorting for Heatmaps (Figure S3)

- Gain and loss BED files are sorted by column 5 (peak score).
- Sorted order is preserved in heatmap generation using `computeMatrix` and `plotHeatmap`.

### 3. Gene Annotation and ENRICHR Lists (Figure 4B, 4C)

- Downloads and processes GENCODE v45 annotations to extract ±200 bp around protein-coding genes.
- Downloads ENCODE cCREs and converts them from `.bb` to `.bed` format.
- Intersects MOA-seq peaks with gene annotations and cCREs using `bedtools`.
- Outputs gene lists suitable for submission to ENRICHR.

## Example Commands

To run the pipeline:

```bash
bash Script.sh
```

## Dependencies

- [MACS3](https://github.com/macs3-project/MACS)
- [bedtools](https://bedtools.readthedocs.io/)
- deepTools (`computeMatrix`, `plotHeatmap`)
- `wget`, `awk`, `sed`, `sort`, `shuf`
- `bigBedToBed` utility from UCSC tools: https://hgdownload.soe.ucsc.edu/admin/exe/

## Notes

- Genome size parameters (`-d1`, `-d2`, `-g`) are hard-coded and should reflect actual read depths.
- The `--nomodel` and `--extsize` options in `macs3 callpeak` reflect the fragment size used in MOA-seq.
- Ensure Internet access is available for downloading GENCODE and ENCODE files.
- Visualizations (SVG files) are time point-specific and reflect peak signal distributions relative to gain/loss peaks.

```bash
# Example: visualize gain/loss MOA-seq signal at 3h vs 0h
cd 3h
computeMatrix reference-point --referencePoint center -b 500 -a 500 \
  -R MOA3_gain_diff.bed MOA3_loss_diff.bed \
  -S "$Bed_3h" "$Bed_0h" --skipZeros --missingDataAsZero \
  --sortRegions keep -o 3hr_matrix.gz

plotHeatmap -m 3hr_matrix.gz --sortRegions no --colorMap Blues Greys \
  --legendLocation none -out 3.svg
```

## Outputs for Downstream Use

- `AlldiffMOA_ccre_GENES.bed` and `AlldiffMOA_nonccre_GENES.bed` can be submitted to:
  - [ENRICHR](https://maayanlab.cloud/Enrichr/)
  - [GSEA](https://www.gsea-msigdb.org/gsea/index.jsp)
- Figures 4 and S3 are generated as `.svg` files, suitable for publication.
