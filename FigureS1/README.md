# Peak Annotation Pipeline Using ChIPseeker

This repository provides a pipeline to annotate hypoxia-induced transcription factor binding peaks using `ChIPseeker` in R. The workflow includes downloading public peak data, filtering peaks based on overlap with regulatory elements, and generating genomic annotation plots.

---

## 🔧 Requirements

### Software
- `bedtools`
- `wget`
- `gzip`
- `bigBedToBed` (available from: https://hgdownload.soe.ucsc.edu/admin/exe/)

### R Packages
- `ChIPseeker`
- `TxDb.Hsapiens.UCSC.hg38.knownGene`
- `clusterProfiler`
- `AnnotationDbi`
- `org.Hs.eg.db`
- `rsvg`

---

## 📥 Step 1: Download and Prepare Input Files

### GEO NarrowPeak File (Hypoxia 0h Sample)

Download from GEO:  
[GSE289863 - NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289863)

```bash
# Download and unzip the narrowPeak file
FILE_PATH="GSE289863_0h_hpx_idr_allReps.narrowPeak.gz"
gzip -d "$FILE_PATH"
FILE_PATH="${FILE_PATH%.gz}"
ENCODE cCREs File
bash
Copy
Edit
# Download and convert ENCODE cCREs bigBed to BED format
wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb
bigBedToBed encodeCcreCombined.bb encodeCcreCombined.bed
rm encodeCcreCombined.bb
🧪 Step 2: Filter Peaks by Overlap
Remove Peaks Overlapping cCREs
bash
Copy
Edit
bedtools intersect -v -a GSE289863_0h_hpx_idr_allReps.narrowPeak -b encodeCcreCombined.bed > GSE289863_0h_hpx_idr_allReps_noCRE.narrowPeak
Remove Peaks Overlapping DHS
bash
Copy
Edit
# Requires UCSC_ENCODE_HUVEC_PEAKS.bed (DHS peaks in HUVEC)
bedtools intersect -v -a GSE289863_0h_hpx_idr_allReps.narrowPeak -b UCSC_ENCODE_HUVEC_PEAKS.bed > GSE289863_0h_hpx_idr_allReps_noDHS.narrowPeak
Remove Peaks Overlapping Both DHS and cCREs
bash
Copy
Edit
bedtools intersect -v -a GSE289863_0h_hpx_idr_allReps_noDHS.narrowPeak -b encodeCcreCombined.bed > GSE289863_0h_hpx_idr_allReps_noCRE_noDHS.narrowPeak
📊 Step 3: Annotate Peaks in R
Use the following R script to annotate the filtered and unfiltered peak sets and generate genomic distribution barplots.

r
Copy
Edit
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(rsvg)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

samplefiles <- list(
  "MOA0"    = "GSE289863_0h_hpx_idr_allReps.narrowPeak",
  "MOA0nc"  = "GSE289863_0h_hpx_idr_allReps_noCRE.narrowPeak",
  "MOA0nd"  = "GSE289863_0h_hpx_idr_allReps_noDHS.narrowPeak",
  "MOA0ncd" = "GSE289863_0h_hpx_idr_allReps_noCRE_noDHS.narrowPeak"
)

for (sample in names(samplefiles)) {
  peakAnno <- annotatePeak(samplefiles[[sample]], TxDb = txdb, tssRegion = c(-200, 200), verbose = FALSE)
  plot <- plotAnnoBar(peakAnno)
  pdf(file = paste0("plot_", sample, "_annoBar.pdf"), width = 10, height = 6)
  print(plot)
  dev.off()
}
📁 Output Files
This pipeline will generate the following barplot PDFs showing peak annotation distributions:

plot_MOA0_annoBar.pdf

plot_MOA0nc_annoBar.pdf

plot_MOA0nd_annoBar.pdf

plot_MOA0ncd_annoBar.pdf

📌 Notes
The file UCSC_ENCODE_HUVEC_PEAKS.bed must be manually downloaded or curated from ENCODE DHS datasets in HUVEC cells.

Ensure bigBedToBed is available in your PATH. You can download precompiled binaries from: UCSC Tools

📚 Citation
Please cite the following resources if you use this workflow:

GEO Dataset: GSE289863

ENCODE Project: ENCODE Candidate cis-Regulatory Elements (cCREs)

ChIPseeker: Yu et al., Bioinformatics (2015)

BEDTools: Quinlan and Hall, Bioinformatics (2010)
