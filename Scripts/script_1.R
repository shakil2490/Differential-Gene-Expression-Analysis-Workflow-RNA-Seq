library(DESeq2)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(SummarizedExperiment)

count_data <- read.csv("C:/Project_2025/Diff_Genes_Expression_DESeq2/Data/count_matrix.csv", row.names = 1, check.names = FALSE)

# Clean column names to keep only sample IDs (e.g., "carcinoma_31")
colnames(count_data) <- gsub(
  "/scratch1/shakilma/gene_expression_pro/hisat2/|\\.bam", 
  "", 
colnames(count_data))

# Verify cleaned column names
colnames(count_data)

# Load metadata (with sample names as row names)
metadata <- read.csv("C:/Project_2025/Diff_Genes_Expression_DESeq2/Data/metadata.csv", row.names = 1)

# Verify alignment
all(colnames(count_data) == rownames(metadata))









