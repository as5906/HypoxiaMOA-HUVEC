library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(rsvg)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
samplefiles <- list("GSE289863_0h_hpx_idr_allReps.narrowPeak")
# Convert to named list
samplefiles <- as.list(samplefiles)
names(samplefiles) <- "MOA0"
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, tssRegion=c(-200, 200), verbose=FALSE)
plotAnnoBar(peakAnnoList)
plot <- plotAnnoBar(peakAnnoList)
pdf(file = "plot_MOA0_P001_annoBar.pdf", width = 10, height = 6)
print(plot)
dev.off()


samplefiles <- list("GSE289863_0h_hpx_idr_allReps_noCRE_noDHS.narrowPeak")
# Convert to named list
samplefiles <- as.list(samplefiles)
names(samplefiles) <- "MOA0ncd"
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, tssRegion=c(-200, 200), verbose=FALSE)
plotAnnoBar(peakAnnoList)
plot <- plotAnnoBar(peakAnnoList)
pdf(file = "plot_MOA0_P001_noDHS_noCRE_annoBar.pdf", width = 10, height = 6)
print(plot)
dev.off()


samplefiles <- list("GSE289863_0h_hpx_idr_allReps_noDHS.narrowPeak")
# Convert to named list
samplefiles <- as.list(samplefiles)
names(samplefiles) <- "MOA0nd"
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, tssRegion=c(-200, 200), verbose=FALSE)
plotAnnoBar(peakAnnoList)
plot <- plotAnnoBar(peakAnnoList)
pdf(file = "plot_MOA0_P001_noDHS_annoBar.pdf", width = 10, height = 6)
print(plot)
dev.off()


samplefiles <- list("GSE289863_0h_hpx_idr_allReps_noCRE.narrowPeak")
# Convert to named list
samplefiles <- as.list(samplefiles)
names(samplefiles) <- "MOA0nc"
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, tssRegion=c(-200, 200), verbose=FALSE)
plotAnnoBar(peakAnnoList)
plot <- plotAnnoBar(peakAnnoList)
pdf(file = "plot_MOA0_P001_noCRE_annoBar.pdf", width = 10, height = 6)
print(plot)
dev.off()
