# 🧬 Integration of HIF1A ChIP-seq with MOA-seq Dynamics and Motif Enrichment

This pipeline supports the analysis underlying Figure 8, which investigates HIF1A occupancy across MOA-seq-defined dynamic peaks, separating gain and loss regions and comparing HIF1A motif enrichment across clusters and timepoints.

---

## 🔧 Dependencies

**Required tools:**
- `bedtools` (v2.30+)
- `awk`, `cat`, `sort`, `shuf`, `wc`, `bc`
- MEME Suite (for downstream motif analysis)

---

## 📂 Input Files

- Differential MOA-seq peaks:
  - `MOA1_gain_diff.bed`, `MOA3_gain_diff.bed`, `MOA24_gain_diff.bed`
  - `MOA1_loss_diff.bed`, `MOA3_loss_diff.bed`, `MOA24_loss_diff.bed`
- ChIP-seq data:
  - `remap2022_HIF1A_nr_macs2_hg38_v1_0.bed`
- Cluster definitions:
  - `hif1_clusters` (HMCs 2, 3, 4, 5, 9)
  - `nonhif1_clusters` (HMCs 1, 6, 7, 8, 10)
- HIF1A motif calls:
  - `1h_gain_sites.tsv.HIF`, `3h_gain_sites.tsv.HIF`, `24h_gain_sites.tsv.HIF`
  - `1h_loss_sites.tsv.HIF`, `3h_loss_sites.tsv.HIF`, `24h_loss_sites.tsv.HIF`

---

## 🔬 Panel A–B: HIF1A ChIP Intersection with MOA-seq GAIN and LOSS Peaks

### Merge GAIN peaks and intersect with HIF1A ChIP
```bash
cat MOA1_gain_diff.bed MOA3_gain_diff.bed MOA24_gain_diff.bed > temp
bedtools sort -i temp > temp2
bedtools merge -i temp2 > all_gain_diff_merged.bed

bedtools intersect -u -a all_gain_diff_merged.bed -b remap2022_HIF1A_nr_macs2_hg38_v1_0.bed > all_gain_diff_merged_hif1a_chip.bed
```
Result: 439 GAIN peaks intersect HIF1A ChIP.

### Compare to HIF1A-associated and non-associated clusters
```bash
bedtools intersect -u -a all_gain_diff_merged_hif1a_chip.bed -b hif1_clusters > m  # 438
bedtools intersect -u -a all_gain_diff_merged_hif1a_chip.bed -b nonhif1_clusters > m  # 1
```

### Repeat for LOSS peaks
```bash
cat MOA1_loss_diff.bed MOA3_loss_diff.bed MOA24_loss_diff.bed > temp
bedtools sort -i temp > temp2
bedtools merge -i temp2 > all_loss_diff_merged.bed

bedtools intersect -u -a all_loss_diff_merged.bed -b remap2022_HIF1A_nr_macs2_hg38_v1_0.bed > all_loss_diff_merged_hif1a_chip.bed

bedtools intersect -u -a all_loss_diff_merged_hif1a_chip.bed -b hif1_clusters > m  # 0
bedtools intersect -u -a all_loss_diff_merged_hif1a_chip.bed -b nonhif1_clusters > m  # 197
```

---

## Motif Enrichment and Random Background Comparison

### Merge HIF1A motifs from GAIN and LOSS
```bash
# GAIN
cat 1h_gain_sites.tsv.HIF 3h_gain_sites.tsv.HIF 24h_gain_sites.tsv.HIF > temp
bedtools sort -i temp > temp2
bedtools merge -i temp2 > hif_motifs_gain_merged.bed

# LOSS
cat 1h_loss_sites.tsv.HIF 3h_loss_sites.tsv.HIF 24h_loss_sites.tsv.HIF > temp
bedtools sort -i temp > temp2
bedtools merge -i temp2 > hif_motifs_loss_merged.bed

rm temp temp2
```

### Intersect HIF1A motif with ChIP-bound peaks
```bash
bedtools intersect -u -a all_gain_diff_merged_hif1a_chip.bed -b hif_motifs_gain_merged.bed > temp  # 173
bedtools intersect -u -a all_loss_diff_merged_hif1a_chip.bed -b hif_motifs_loss_merged.bed > temp  # 29
```

### Random sampling to assess enrichment

100 iterations of random sampling were performed for both GAIN and LOSS peaks. In each iteration:

1. Random subset of N peaks drawn (`N = 438` for GAIN, `197` for LOSS).
2. Intersect with respective HIF1A motif set.
3. Count overlaps and store.
4. Compute mean, SD, and Z-score of observed motif count relative to random.

#### Example snippet (GAIN):
```bash
output_file="line_counts.txt"
total_lines=0
line_counts=()

for ((i=1; i<=100; i++))
do
    shuf -n 438 all_gain_diff_merged.bed > temp
    bedtools intersect -u -a temp -b hif_motifs_gain_merged.bed > temp2
    line_count=$(wc -l < temp2)
    echo "$line_count" >> "$output_file"
    ((total_lines += line_count))
    line_counts+=($line_count)
done

average_lines=$((total_lines / 100))
sum_squared_diff=0
for count in "${line_counts[@]}"; do
    diff=$((count - average_lines))
    squared_diff=$((diff * diff))
    sum_squared_diff=$((sum_squared_diff + squared_diff))
done
stdev=$(echo "scale=2; sqrt($sum_squared_diff / 100)" | bc)

observed_line_count=173
z_score=$(echo "scale=2; ($observed_line_count - $average_lines) / $stdev" | bc)
```

---

## 📈 Panel E: Generate 101 bp Sites for Motif Enrichment (MEME)

### From HIF1A-associated clusters (GAIN)
```bash
bedtools intersect -u -b hif_motifs_gain_merged.bed -a hif1_clusters > temp
awk '{mid=int(($2+$3)/2); print $1 "\t" mid-50 "\t" mid+50}' temp > hif1_clusters_hif1a_motif_101.bed
```

### From non-HIF1A clusters (LOSS)
```bash
bedtools intersect -u -b hif_motifs_loss_merged.bed -a nonhif1_clusters > temp
awk '{mid=int(($2+$3)/2); print $1 "\t" mid-50 "\t" mid+50}' temp > nonhif1_clusters_hif1a_motif_101.bed
```

These files are used as input to MEME Suite for co-motif analysis (not shown in this script).

---

## 📦 Output Files

| File | Description |
|------|-------------|
| `all_gain_diff_merged_hif1a_chip.bed` | GAIN peaks intersecting HIF1A ChIP |
| `all_loss_diff_merged_hif1a_chip.bed` | LOSS peaks intersecting HIF1A ChIP |
| `hif_motifs_gain_merged.bed` | All HIF1A motifs in GAIN peaks |
| `hif_motifs_loss_merged.bed` | All HIF1A motifs in LOSS peaks |
| `line_counts.txt` | Random motif counts for GAIN (100x) |
| `line_counts2.txt` | Random motif counts for LOSS (100x) |
| `hif1_clusters_hif1a_motif_101.bed` | Sites for MEME motif enrichment (HIF1 clusters) |
| `nonhif1_clusters_hif1a_motif_101.bed` | Sites for MEME motif enrichment (non-HIF1 clusters) |

---

## 📌 Summary

This workflow systematically assesses:
- Enrichment of HIF1A binding in GAIN vs. LOSS regions.
- Distribution of binding across HIF1-associated and non-associated modules.
- Motif-level evidence supporting HIF1A binding specificity.
- Statistical confidence in motif enrichment via permutation testing.
- Sequence-level extraction for co-motif discovery.
