# HypoxiaMOA-HUVEC

## Scripts for Analyzing MOA-seq Data in HUVECs under Hypoxia

### Description
This repository contains a collection of scripts developed to process MOA-seq data from HUVECs under hypoxic conditions. The pipeline covers initial processing of paired-end reads, CPM normalization, and downstream analyses, including clustering, differential expression, and Gene Set Enrichment Analysis (GSEA). The provided tools allow you to:
- Trim adapters
- Align reads
- Generate bedgraphs/bigwigs
- Perform hierarchical clustering of gene expression data
- Run DESeq2 for differential expression
- Explore enriched gene ontology terms across time points

## Installation & Dependencies
The MOA-seq processing pipeline requires standard Linux utilities and several bioinformatics tools:

### Bash Scripts Dependencies (MOAseq.GeneralPipeline.Linux.CPM.sh)
- `samtools`
- `cutadapt`
- `flash` (optional for paired-end merging)
- `bwa-mem2`
- `bedtools`
- `bedGraphToBigWig`
- `bc` and `gawk`

### R Scripts Dependencies
Requires `R` with the following packages:
- **For clustering:** `ggplot2`, `preprocessCore`, `tidyr`, `dplyr`
- **For differential expression:** `DESeq2`
- **For enrichment analysis:** `clusterProfiler`, `enrichplot`, `ggplot2`, `svglite`, `org.Hs.eg.db`

---

## Scripts & Usage

### **MOAseq.GeneralPipeline.Linux.CPM.sh**
This Bash script processes paired-end MOA-seq reads to generate CPM-normalized outputs. It performs the following steps:

1. Adapter removal using `cutadapt`
2. (Optional) Merging of paired-ends using `flash`
3. Alignment with `bwa-mem2` and conversion to a sorted BAM file
4. Generation of fragment BED files and calculation of fragment centers
5. Scaling for CPM normalization and creation of bedgraph and bigwig files

#### **Usage:**
```bash
./MOAseq.GeneralPipeline.Linux.CPM.sh -a <Paired-End_File1> -b <Paired-End_File2> -c <Reference_Genome_FASTA> -d <Genome_BWA_Index_Path/Prefix> [-e <Effective_Genome_Size>] [-f <MAPQ_Threshold>] [-g]
```
**Arguments:**
- `-a`, `-b`: Input paired-end FASTQ files
- `-c`: Reference genome in FASTA format
- `-d`: BWA index path and prefix for the reference genome
- `-e` (optional): Provide your own effective genome size
- `-f` (optional): Set a MAPQ threshold (default is 0)
- `-g` (optional): Flag to skip merging with `flash` if required

---

### **Hierarchical_norm.R**
This R script standardizes gene expression data and performs hierarchical clustering to explore MOA profiles over time. It:
- Reads a tab-delimited gene data file
- Standardizes the expression values
- Computes within-cluster sum of squares (WCSS) to help determine the optimal number of clusters via the elbow method
- Performs hierarchical clustering and assigns clusters
- Plots the gene expression profiles for each cluster and saves the combined plot as an SVG

#### **Usage:**
Run the script interactively in R or via `Rscript`.
```r
Rscript Hierarchical_norm.R
```
You will be prompted to input:
- The path to your gene data file (columns for 0h, 1h, 3h, and 24h)
- The output file path to save the combined cluster plot

---

### **Clean_Reps_for_DESEQ2.R**
This R script prepares and processes count data using `DESeq2` for differential expression analysis. It:
- Prompts for a counts file (tab-delimited with gene counts across replicates)
- Constructs sample information based on four time points (0h, 1h, 3h, and 24h)
- Filters out lowly expressed genes (keeping those with at least 10 reads across samples)
- Creates a `DESeq2` dataset, relevels conditions (using 0h as the reference), and runs the analysis
- Extracts log2 fold change values for 1h, 3h, and 24h (each versus 0h) and compiles them into a single table

#### **Usage:**
Run the script in R or via `Rscript`.
```r
Rscript Clean_Reps_for_DESEQ2.R
```
You will be prompted to specify:
- The path to your counts file
- The output file path where the combined log2FC table will be saved

---

### **GSEA.R**
This R script performs Gene Set Enrichment Analysis (GSEA) on differential expression results using the `clusterProfiler` package. It is designed to assess GO enrichment for a specified cluster across different time points. The script:
- Prompts for a cluster number
- Reads input files for log2 fold change data at 24h, 3h, and 1h (each expected in tab-delimited format with gene symbols and log2FC values)
- Sorts the gene lists and performs GSEA (using the BP ontology) for each time point
- Calculates gene ratios from core enrichment data and flags enrichment direction based on the Normalized Enrichment Score (NES)
- Writes CSV files with the selected GO terms (Description, p-value, NES, gene ratio, and enrichment sign) for each time point

#### **Usage:**
Run the script interactively in R or via `Rscript`.
```r
Rscript GSEA.R
```
The script will prompt for:
- A cluster number to analyze
- File paths for the 24h, 3h, and 1h log2FC input files
- Output file paths for the generated CSV results for each time point

---
