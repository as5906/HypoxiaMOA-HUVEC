library(DESeq2)

# Prompt user for counts file path
counts_file <- readline(prompt = "Enter the path to the counts file (e.g., '~/desktop/counts_reps_col.txt'): ")
countData <- read.table(counts_file, header = TRUE, row.names = 1)

# Prepare sample information
sample_names <- colnames(countData)
condition <- factor(rep(c("0h", "1h", "3h", "24h"), times = c(3, 3, 3, 3)))  # Update the times argument

colData <- data.frame(
  row.names = sample_names,
  condition = condition
)

# Filter out lowly expressed genes
# Keep genes with at least 10 reads in total across all samples
keep <- rowSums(countData) >= 10
countData <- countData[keep,]

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)

# Relevel to make "0h" the reference level
dds$condition <- relevel(dds$condition, ref = "0h")

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results
res_1h_vs_0h <- results(dds, contrast = c("condition", "1h", "0h"))
res_3h_vs_0h <- results(dds, contrast = c("condition", "3h", "0h"))
res_24h_vs_0h <- results(dds, contrast = c("condition", "24h", "0h"))

# Combine log2FC values into one data frame
log2fc_combined <- data.frame(
  Gene = rownames(res_1h_vs_0h),
  log2FC_1h = res_1h_vs_0h$log2FoldChange,
  log2FC_3h = res_3h_vs_0h$log2FoldChange,
  log2FC_24h = res_24h_vs_0h$log2FoldChange
)

# Prompt user for output file path to save results
output_file <- readline(prompt = "Enter the path to save the log2FC combined file (e.g., '~/desktop/log2FC_combined_test.tsv'): ")

# Write to tab-separated file
write.table(log2fc_combined, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

