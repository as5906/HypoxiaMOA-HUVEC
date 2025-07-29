# MOA-seq Differential Peak and DEG Integration Pipeline

This directory contains scripts and data used to integrate differential MOA-seq peaks with differentially expressed genes (DEGs) across hypoxia time points (1h, 3h, 24h), and to generate gene lists for downstream pathway enrichment analysis (e.g., ENRICHR). It also includes the workflow to create randomized DEG subsets for statistical comparison.

## Contents

- `1h_DEG_q05.txt`, `3h_DEG_q05.txt`, `24h_DEG_q05.txt`:  
  Differentially expressed gene (DEG) lists at 1h, 3h, and 24h of hypoxia, filtered at q < 0.05.

- `Script`:  
  Main Bash script to:
  - Download and process GENCODE v45 gene annotations.
  - Identify protein-coding genes overlapping with MOA-seq gain/loss peaks at each time point.
  - Generate gene lists for each category (gain, loss, or both).
  - Create subsets of DEGs intersecting with MOA-seq peaks for enrichment analysis.
  - Create randomized DEG subsets of matching size to assess significance (Figure 5C).

- `Count.sh`:  
  Helper script to annotate each DEG as:
  - `gain`: overlaps a MOA-seq gain peak
  - `loss`: overlaps a MOA-seq loss peak
  - `gain,loss`: overlaps both
  - `na`: overlaps neither  
  Produces a 4-column output with gene ID and annotation.

---

## Workflow Summary

### Figure 5A — Identifying MOA-overlapping protein-coding genes
1. Downloads and processes GENCODE v45 annotations.
2. Extracts ±200 bp around gene bodies for protein-coding genes.
3. Uses `bedtools intersect` to find MOA gain/loss peaks overlapping with genes.
4. Outputs:
   - `MOA{1,3,24}_gain_GENES.txt`
   - `MOA{1,3,24}_loss_GENES.txt`
   - `MOA{1,3,24}_diff_GENES.bed` for use in ENRICHR.

### Figure 5B — Mapping DEGs to MOA-seq peaks
1. For each DEG file (`*_DEG_q05.txt`), genes are annotated using `Count.sh`:
   - Whether they overlap gain, loss, both, or neither MOA peaks.
2. Output files:
   - `1hr_MOA_DEG.txt`, `3hr_MOA_DEG.txt`, `24hr_MOA_DEG.txt`: Annotated DEG lists.
   - `*_MOA_DEG_subset.txt`: Only genes with MOA-seq overlap (used in enrichment).

### Figure 5C — Randomized controls
1. Generates 3 random subsets of DEGs (same size as MOA-overlapping set) per time point:
   - `*_MOA_DEG_rand_subset{1,2,3}.txt`
2. These are used as controls for comparison in pathway enrichment analysis.

---

## Notes

- `Count.sh` is internally called by `Script` and should not be run independently unless debugging.
- Gene lists for enrichment should be filtered to remove `na` (non-overlapping) entries.
- Gene identifiers are assumed to be consistent across DEGs and BED files (column 1 = gene name).

---

## Dependencies

- `bedtools`
- `awk`, `sed`, `shuf`
- Internet access to download GENCODE v45 from ENSEMBL

---

## Example Usage

```bash
bash Script
