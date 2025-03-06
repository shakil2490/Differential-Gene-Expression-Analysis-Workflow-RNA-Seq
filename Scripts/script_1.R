library(DESeq2)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(SummarizedExperiment)
BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86) 
library(biomaRt)
library(org.Hs.eg.db)

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

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = metadata,
  design = ~ condition)

# Keep genes with ≥10 counts in at least 20% of samples
min_samples <- ceiling(0.2 * ncol(counts(dds)))
keep <- rowSums(counts(dds) >= 10) >= min_samples
dds <- dds[keep, ]

dds <- DESeq(dds)

# Extract results (carcinoma vs normal)
res <- results(dds,contrast = c("condition", "carcinoma", "normal"), 
                alpha = 0.05)  # Adjusted p-value threshold

# Map Ensembl IDs to gene symbols using EnsDb.Hsapiens.v86
symbols_ens <- mapIds(
  EnsDb.Hsapiens.v86,
  keys = rownames(res),       # Your Ensembl IDs
  column = "SYMBOL",          # Column for gene symbols
  keytype = "GENEID",         # Type of input IDs (Ensembl)
  multiVals = "first"         # Handle duplicates by taking the first match
)

# Map Ensembl IDs to gene symbols using org.Hs.eg.db
symbols_org <- mapIds(
  org.Hs.eg.db,
  keys = rownames(res),      
  column = "SYMBOL",          # Column for gene symbols
  keytype = "ENSEMBL",        # Type of input IDs (Ensembl)
  multiVals = "first"         # Handle duplicates by taking the first match
)

# Combine results (prioritize EnsDb, fall back to org.Hs.eg.db)
res$symbol <- ifelse(is.na(symbols_ens), symbols_org, symbols_ens)


# Check for remaining missing symbols
sum(is.na(res$symbol))

# Handle remaining missing symbols
res$symbol <- ifelse(is.na(res$symbol), rownames(res), res$symbol)

head(res$symbol)

# Identify non-standard gene symbols
non_standard_symbols <- grepl("^RP[0-9]+-|^AC[0-9]+\\.[0-9]+|^AL[0-9]+\\.[0-9]+", res$symbol)

# Count non-standard symbols
sum(non_standard_symbols)

# View non-standard symbols
head(res$symbol[non_standard_symbols])

res$symbol[non_standard_symbols] <- rownames(res)[non_standard_symbols]

# Convert to data frame and save
res_df <- as.data.frame(res)
write.csv(res_df, "deseq2_results_annotated.csv")


















