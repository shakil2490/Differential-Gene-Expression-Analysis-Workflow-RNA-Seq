#!/bin/bash
#SBATCH --job-name=deseq2_analysis
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=deseq2_%j.log

# Environment variables
export RNA_LAB="/scratch1/shakilma/gene_expression_pro"
export RNA_DATA_DIR="$RNA_LAB/fastq_file/fastq_QC_after"
export RNA_REFS_DIR="$RNA_LAB/ref"
export RNA_GTF="$RNA_LAB/gtf/Homo_sapiens.GRCh38.112.gtf"

# Reference files
export REF_FASTA="$RNA_REFS_DIR/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
export HISAT2_INDEX="$RNA_REFS_DIR/Homo_sapiens.GRCh38.dna_sm.primary_assembly"

# Output directories
export ALIGN_DIR="$RNA_LAB/hisat2"
export COUNTS_DIR="$RNA_LAB/counts"
export DE_DIR="$RNA_LAB/de/deseq2"

# Create directories
mkdir -p "$ALIGN_DIR" "$COUNTS_DIR" "$DE_DIR"

# --------------------------------------------------
# 1. Alignment with HISAT2 (same as before)
# --------------------------------------------------
if [ ! -f "${HISAT2_INDEX}.1.ht2" ]; then
  echo "Building HISAT2 index..."
  hisat2-build -p 8 "$REF_FASTA" "$HISAT2_INDEX"
fi

for R1 in "$RNA_DATA_DIR"/*_1.fastq.gz; do
  SAMPLE=$(basename "$R1" _1.fastq.gz)
  R2="${R1/_1.fastq.gz/_2.fastq.gz}"
  
  hisat2 -p 8 --rg-id="$SAMPLE" -x "$HISAT2_INDEX" \
    -1 "$R1" -2 "$R2" -S "$ALIGN_DIR/${SAMPLE}.sam"
  
  samtools sort -@ 8 -o "$ALIGN_DIR/${SAMPLE}.bam" "$ALIGN_DIR/${SAMPLE}.sam"
  rm "$ALIGN_DIR/${SAMPLE}.sam"
done

# --------------------------------------------------
# 2. Generate Count Matrix with featureCounts
# --------------------------------------------------
BAM_FILES=($ALIGN_DIR/*.bam)
featureCounts -T $SLURM_CPUS_PER_TASK -p -a "$RNA_GTF" \
  -o "$COUNTS_DIR/counts.txt" ${BAM_FILES[@]}

# Clean up featureCounts output
cut -f1,7- "$COUNTS_DIR/counts.txt" | grep -v '^#' > "$COUNTS_DIR/count_matrix.csv"

# --------------------------------------------------
# 3. Create Sample Metadata
# --------------------------------------------------
echo "sample,condition" > "$COUNTS_DIR/metadata.csv"
for bam in "${BAM_FILES[@]}"; do
  SAMPLE=$(basename "$bam" .bam)
  CONDITION=$(echo "$SAMPLE" | cut -d'_' -f1)  # Extract "carcinoma" or "normal"
  echo "$SAMPLE,$CONDITION" >> "$COUNTS_DIR/metadata.csv"
done
