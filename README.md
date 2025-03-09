
# Prostate Cancer Transcriptomic Analysis: RNA-seq Workflow

## Introduction and Background

Prostate cancer is the second most common cancer in men worldwide, with significant heterogeneity in clinical outcomes. Molecular characterization of prostate tumors has identified key drivers such as gene fusions, somatic mutations, and dysregulated non-coding RNAs (e.g., PCA3). RNA sequencing (RNA-seq) enables comprehensive profiling of transcriptional aberrations, including alternative splicing, novel transcripts, and differential expression of coding and non-coding genes, which are critical for understanding tumor biology and identifying biomarkers.

This study analyzes a subset of the publicly available dataset **GSE22260** (NCBI GEO), which focuses on comparative transcriptomic analysis of prostate cancer tumors and matched normal tissues.

## Dataset

For this study, I selected a subset of 19 samples from the **GSE22260** dataset:
- 10 control samples representing normal tissues.
- 9 carcinoma samples representing tumor tissues.

The libraries in this dataset were:
- PolyA-selected, ensuring mRNA capture.
- Prepared as paired end reads to improve alignment accuracy and transcript assembly.
- Sequenced on an Illumina Genome Analyzer II, generating reads that are 36 bp long.

## Objectives

The primary objective of this study was to assess the prostate cancer-specific expression of the PCA3 gene using RNA-Seq data. The PCA3 gene is widely studied for its potential role as a biomarker in prostate cancer. The following steps were undertaken to achieve this:

1. Obtain data and references: RNA-Seq data was obtained from the **GSE22260** dataset.
2. Align the RNA-Seq FASTQ files:Reads were aligned to the human reference genome to generate transcriptomic profiles.
3. Quantify transcript expression:Expression levels of genes, including PCA3, were quantified.
4. Perform differential expression analysis:Differential expression analysis was conducted to identify significant genes associated with carcinoma compared to normal tissue samples.
