# Motif Intersection Analysis with ReMap ChIP-seq

This directory contains BED files representing gain and loss motif regions from MOA-seq footprints across hypoxia time points (1h, 3h, 24h), along with a bash script used to calculate overlap percentages with ReMap 2022 transcription factor binding sites.

## Contents

### BED Files

| File Name                | Description                                      |
|--------------------------|--------------------------------------------------|
| `1hr_gain_motifs.bed`    | Motifs gained at 1-hour hypoxia                  |
| `1hr_loss_motifs.bed`    | Motifs lost at 1-hour hypoxia                   |
| `3hr_gain_motifs.bed`    | Motifs gained at 3-hour hypoxia                  |
| `3hr_loss_motifs.bed`    | Motifs lost at 3-hour hypoxia                   |
| `24hr_gain_motifs.bed`   | Motifs gained at 24-hour hypoxia                 |
| `24hr_loss_motifs.bed`   | Motifs lost at 24-hour hypoxia                  |

These files were generated as part of the motif analysis workflow and are compatible with downstream motif enrichment tools such as [MEME Suite's XStreme](https://meme-suite.org/meme/tools/xstreme).

---

### Script:

This script automates the download of selected ReMap 2022 ChIP-seq peak files for human TFs, decompresses them, and computes the percentage of motif sites (from the above BED files) that overlap with each TF binding site dataset using `bedtools`.

#### Features:

- Downloads `.bed.gz` files for selected TFs from ReMap 2022 (hg38, MACS2 narrow peaks)
- Converts motif names used in analysis to match standard ReMap TF names:
  - `ITF2` → TCF4
  - `NDF1` → NEUROD1
  - `PO2F1` → POU2F1
  - `ZN502` → ZNF502
  - `NF2L2` → NFE2L2
- Calculates overlap using `bedtools intersect -u`
- Outputs the percentage of intersecting sites per motif per timepoint
- Temporary files are cleaned up after execution

#### Usage

```bash
./Script
```

#### Dependencies

- wget
- gunzip
- bedtools
- awk, wc, bc

#### Example Output

1hr_gain_motifs.bed with SMAD3 intersection percentage: 22.17%
1hr_gain_motifs.bed with ITF2 intersection percentage: 15.04%
...

#### Notes
- This analysis complements Figure 6 of the manuscript by quantifying the consistency between predicted motifs and experimentally validated TF binding sites.
- BED files can also be used for input into motif enrichment or clustering tools for further regulatory characterization.

