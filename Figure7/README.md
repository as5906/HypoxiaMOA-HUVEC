# 🧬 Differential Peak Clustering and GSEA Pipeline

This folder contains the full pipeline to process differential MOA-seq peaks across hypoxia timepoints (1h, 3h, 24h), generate normalized coverage matrices, perform clustering, and run GSEA on each resulting cluster using gene-level differential expression statistics.

---
## Dependencies, Tools, and Libraries
Command-Line Tools
- BEDTools (v2.29 or later): Required for merging peak files and computing coverage. https://bedtools.readthedocs.io/
- awk, cat, paste, sort (Unix shell utilities): Used for data manipulation and matrix assembly.

📦 R Libraries
- Make sure R (v4.0 or later) is installed along with the following packages:
- ggplot2: for generating cluster profile plots
- preprocessCore: for quantile normalization (if added)
- tidyr and dplyr: for data reshaping and manipulation
- clusterProfiler: for Gene Set Enrichment Analysis (GSEA)
- enrichplot: for visualizing enrichment results
- org.Hs.eg.db: human gene annotation database
- svglite: for exporting plots to .svg

## 📁 Required Inputs

- Differential peaks:
  - `MOA1_gain_diff.bed`, `MOA3_gain_diff.bed`, `MOA24_gain_diff.bed`
  - `MOA1_loss_diff.bed`, `MOA3_loss_diff.bed`, `MOA24_loss_diff.bed`

- Fragment center BEDs:
  - `MOA0_merge_frenter_q20_chr_sort.bed`
  - `MOA1_merge_frenter_q20_chr_sort.bed`
  - `MOA3_merge_frenter_q20_chr_sort.bed`
  - `MOA24_merge_frenter_q20_chr_sort.bed`

---

## 🏗️ Step 1: Build Merged Differential Peak Set - Using "Script" File

```bash
# Define required file paths for gain/loss peaks
gain_1h=/home/ss19m/HYPOX/PEAKS/DIFF/MOA1_gain_diff.bed
gain_3h=/home/ss19m/HYPOX/PEAKS/DIFF/MOA3_gain_diff.bed
gain_24h=/home/ss19m/HYPOX/PEAKS/DIFF/MOA24_gain_diff.bed
loss_1h=/home/ss19m/HYPOX/PEAKS/DIFF/MOA1_loss_diff.bed
loss_3h=/home/ss19m/HYPOX/PEAKS/DIFF/MOA3_loss_diff.bed
loss_24h=/home/ss19m/HYPOX/PEAKS/DIFF/MOA24_loss_diff.bed

# Combine and merge
diff_out=all_diff_peaks.bed
cat "$gain_1h" "$gain_3h" "$gain_24h" "$loss_1h" "$loss_3h" "$loss_24h" > $diff_out
bedtools sort -i $diff_out > temp && mv temp $diff_out
bedtools merge -i $diff_out > all_diff_peaks_merge.bed
awk '{mid = int(($2 + $3) / 2); start = mid - 15; end = mid + 15; if (start < 0) start = 0; print $1, start, end}' OFS='\t' all_diff_peaks_merge.bed > all_diff_peaks_merge_win30.bed
```

---

## 📊 Compute Coverage at All Timepoints

```bash
# Define BEDs
Bed_0h=MOA0_merge_frenter_q20_chr_sort.bed
Bed_1h=MOA1_merge_frenter_q20_chr_sort.bed
Bed_3h=MOA3_merge_frenter_q20_chr_sort.bed
Bed_24h=MOA24_merge_frenter_q20_chr_sort.bed

# Coverage
bedtools coverage -a all_diff_peaks_merge_win30.bed -b $Bed_0h > coverage_0h.txt
bedtools coverage -a all_diff_peaks_merge_win30.bed -b $Bed_1h > coverage_1h.txt
bedtools coverage -a all_diff_peaks_merge_win30.bed -b $Bed_3h > coverage_3h.txt
bedtools coverage -a all_diff_peaks_merge_win30.bed -b $Bed_24h > coverage_24h.txt
```

---

## ⚖️ CPM Normalization

```bash
awk '{OFS="\t" ; print $1,$2,$3,$4*0.028393227}' coverage_0h.txt > temp && mv temp coverage_0h.txt
awk '{OFS="\t" ; print $1,$2,$3,$4*0.02100052331}' coverage_1h.txt > temp && mv temp coverage_1h.txt
awk '{OFS="\t" ; print $1,$2,$3,$4*0.01935495584}' coverage_3h.txt > temp && mv temp coverage_3h.txt
awk '{OFS="\t" ; print $1,$2,$3,$4*0.01205994235}' coverage_24h.txt > temp && mv temp coverage_24h.txt
```

---

## 🧱 Assemble Matrix for Clustering

```bash
awk '{print $4}' coverage_1h.txt > temp
awk '{print $4}' coverage_3h.txt > temp2
awk '{print $4}' coverage_24h.txt > temp3
paste -d "\t" temp temp2 temp3 > temp4
paste -d "\t" coverage_0h.txt temp4 > df
rm temp temp2 temp3 temp4
```

Alternatively, use `compiled_diffMOA_cpm_byHour.txt` for direct input to clustering.

---

## 🌐 Run Clustering in R - Using "Hierarchical_norm.R"

```r
# Hierarchical_norm.R
library(ggplot2)
library(preprocessCore)

# Read in file
gene_data_file <- readline(prompt = "Enter path to file: ")
gene_data <- read.delim(gene_data_file, header = TRUE, sep = "\t")
gene_column <- gene_data$GENE
gene_data <- gene_data[, -1]

# Standardize
standardized_data <- t(apply(gene_data, 1, scale))
gene_data <- cbind(GENE = gene_column, as.data.frame(standardized_data))
colnames(gene_data) <- c("GENE", "X0hr", "X1hr", "X3hr", "X24hr")
rownames(gene_data) <- gene_data[,1]
gene_data <- gene_data[,-1]

# Elbow method
wcss <- vector()
for (i in 1:30) {
  kmeans_model <- kmeans(gene_data, centers = i)
  wcss[i] <- kmeans_model$tot.withinss
}
plot(1:30, wcss, type = "b", pch = 19, xlab = "Clusters", ylab = "WCSS")

# Hierarchical clustering
hc <- hclust(dist(gene_data))
clusters <- cutree(hc, k = 10)

gene_data_df <- data.frame(Gene = rownames(gene_data), gene_data)
gene_data_df$Cluster <- as.factor(clusters)

library(tidyr)
library(dplyr)
gene_data_long <- gene_data_df %>%
  pivot_longer(cols = -c(Gene, Cluster), names_to = "Sample", values_to = "Expression")
gene_data_long$Sample <- factor(gene_data_long$Sample, levels = c("X0hr", "X1hr", "X3hr", "X24hr"))

combined_plot <- gene_data_long %>%
  filter(Cluster %in% unique(clusters)) %>%
  ggplot(aes(x = Sample, y = Expression, group = Gene, color = Gene)) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~Cluster, scales = "free") +
  ggtitle("MOA Profiles for Each Cluster") +
  theme_minimal()

print(combined_plot)
output_file <- readline(prompt = "Enter path to save SVG: ")
ggsave(output_file, plot = combined_plot, device = "svg")
```

---

## 🧬 GSEA for Cluster DEG Lists - Using "GSEA.R" File in GSEA Directory where the Log2FC files also are

Use `GSEA.r` script with the following structure for each cluster/timepoint log2FC file:

```r
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(svglite)
library(org.Hs.eg.db)

x <- as.integer(readline(prompt = "Enter cluster number: "))
input_file_24h <- readline(prompt = "Enter 24h log2fc path: ")
df_24h <- read.delim(input_file_24h, header = FALSE)
gene_list_24h <- sort(setNames(df_24h$V2, df_24h$V1), decreasing = TRUE)

gse_24h <- gseGO(geneList = gene_list_24h, ont = "BP", keyType = "SYMBOL", 
                 minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, 
                 verbose = TRUE, OrgDb = org.Hs.eg.db, pAdjustMethod = "none")

geneRatio_24h <- sapply(strsplit(gse_24h$core_enrichment, "/"), length) / gse_24h$setSize
selected_columns_24h <- gse_24h[, c("Description", "pvalue", "NES")]
selected_columns_24h$GeneRatio <- geneRatio_24h
selected_columns_24h$Sign <- ifelse(selected_columns_24h$NES > 0, "Positive", "Negative")

output_file_24h <- readline(prompt = "Enter output path for 24h: ")
write.csv(selected_columns_24h, file = output_file_24h, row.names = FALSE)
```

Repeat the above block for `3h` and `1h` using `clusterX_3log2fc` and `clusterX_1log2fc` inputs respectively.

---

## 📂 Input Files

```
cluster1_1log2fc ... cluster10_1log2fc
cluster1_3log2fc ... cluster10_3log2fc
cluster1_24log2fc ... cluster10_24log2fc
GSEA.r
```

---
