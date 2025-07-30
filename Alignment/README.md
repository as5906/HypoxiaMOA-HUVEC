# Genome hg38 Alignment 

This directory contains raw FASTQ files, alignment scripts, and output files for processing MOA-seq data across multiple time points and control conditions. The primary objective is to generate sorted BAM, BED, and bigWig files from paired-end sequencing data using the MOA-seq general pipeline.

---

## Contents

### `Script`
A bash script to automate the MOA-seq alignment pipeline across all samples. The script:
- Downloads and indexes the hg38 reference genome using `bwa-mem2`
- Runs `MOAseq.GeneralPipeline.Linux.CPM.sh` for each replicate at the following time points:
  - 0h
  - 1h
  - 3h
  - 24h
  - Control_1 through Control_4 (MNase controls)

Each replicate directory (e.g., `0h/Rep1`) contains paired-end FASTQ files and stores intermediate/output files from the pipeline.

---

## Required Tools

- `bwa-mem2`
- `samtools`
- `bedtools`
- `flash`
- `gawk`
- `bc`
- `pandas`, `numpy` (for `unique-kmers.py`)
- `bedGraphToBigWig` (UCSC)

Ensure these are accessible via your `$PATH`.

---

## Pipeline Script: `MOAseq.GeneralPipeline.Linux.CPM.sh`

This modular shell script aligns MOA-seq paired-end reads and generates coverage files.

### Usage

```bash
MOAseq.GeneralPipeline.Linux.CPM.sh \
  -a <paired_read_1.fastq.gz> \
  -b <paired_read_2.fastq.gz> \
  -c <reference_genome.fa.gz> \
  -d <bwa_index_prefix> \
  [-e <effective_genome_size>] \
  [-f <MAPQ_threshold>] \
  [-g] # skip FLASH merging
```

### Output

Each run produces:
- `output.sort.bam`: sorted BAM file
- `output.sort.bam.stats.txt`: BAM statistics
- `output.sort.bw`: bigWig of full fragments
- `output.shortenedReads.sort.bw`: bigWig of 20 bp fragment centers (frenters)
- `output.sort.bed` (raw file, full fragment), `output.shortenedReads.sort.bed` (frenter file, shortened fragment): BED files used for coverage tracks
- `tempFiles/`: directory storing intermediate files (SAM, unsorted BAM, fragment statistics)

---

## Key Steps Performed

1. **Optional FLASH merge** of paired-end reads
2. **Alignment** to the hg38 genome using `bwa-mem2`
3. **SAM to BAM** conversion with MAPQ filtering
4. **Fragment length calculation** via `samtools stats`
5. **Effective genome size estimation** using `unique-kmers.py` (unless provided)
6. **BAM to BED** and 20 bp center BED creation
7. **Coverage calculation** using `genomeCoverageBed`
8. **Track file generation** in bigWig format

---

## Notes

- Ensure FASTQ file names and directory structure match the expected format used in `Script`.
- All BED and bigWig files are sorted and ready for downstream analysis and visualization.
- Genome reference and index are assumed to be in `hg38.fa.gz` and `index*`.

---

## Example Call

```bash
MOAseq.GeneralPipeline.Linux.CPM.sh \
  -a Va_MOA3_R1.fq.gz \
  -b Va_MOA3_R2.fq.gz \
  -c hg38.fa.gz \
  -d index \
  -f 20
```

---

## Structure Summary

```
alignment/
│
├── 0h/
│   ├── Rep1/
│   ├── Rep2/
├── 1h/
│   ├── Rep1/
│   ├── Rep2/
├── 3h/
│   ├── Rep1/
│   ├── Rep2/
├── 24h/
│   ├── Rep1/
│   ├── Rep2/
├── Control_1/ to Control_4/
│   ├── Rep1/
│   ├── Rep2/
├── hg38.fa.gz
├── index.*
├── alignment_script.sh
├── MOAseq.GeneralPipeline.Linux.CPM.sh
```
