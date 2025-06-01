# Mosquito Data

## Introduction

This project aims to determine the transcriptional changes associated with temperature-related feeding in _Anopheles coluzzii_ mosquitoes. The mosquitoes were fed at five temperatures (25, 30, 36, 40, 48). Mosquitoes that did not behaviourally respond to any of the temperature cues were collected and associated with a condition called "Unresponsive". Heads and bodies were collected individually.

## Project Pipeline Overview

Raw Read Quality Control
➤ Tool: FastQC
➤ Purpose: Assess the quality of raw sequencing reads.

Read Trimming
➤ Tool: Trimmomatic 
➤ Purpose: Remove adapter sequences and low-quality bases.

Post-trimming Quality Control
➤ Tool: FastQC
➤ Purpose: Confirm improvement in read quality after trimming.

Transcriptome Alignment / Quantification
➤ Tool: Salmon
➤ Purpose: Quasi-mapping reads to the transcriptome and quantifying transcript abundance.

Gene-level Aggregation and Differential Expression Analysis
➤ Tools: tximport → DESeq2
➤ Purpose: Import transcript-level estimates and perform differential gene expression analysis.

## Data Organisation: Directory Structure and Scripts

```mosquito_data/
├── FastQC_results/           # Quality reports for raw reads (FastQC output)
├── trimmed_reads/            # Trimmed FASTQ files (from Trimmomatic)
├── trimmed_fastqc_results/   # QC reports for trimmed reads (FastQC)
├── alignments/               # Quantification output (e.g., Salmon results)
├── logs/                     # Log files from trimming, alignment, etc.

git_repos/
├── Reference_genome/         # An. gambiae PEST genome/transcriptome (FASTA, GTF)
├── rnaseq_data/              # Raw data and RNA-seq processing scripts
│   ├── fastqc_array.sh             # Script to run FastQC on raw reads
│   ├── trimmomatic_array.sh        # Script for adapter trimming and quality filtering
│   └──  salmon_alignment.sh             # Script to quantify expression using Salmon 
  ```
  

## Metadata


 <pre> ```sampleName	replicate	temperature	tissue
UJ-3092-25-1B_quant	1	25	Body
UJ-3092-25-1H_quant	1	25	Head
UJ-3092-25-2B_quant	2	25	Body
UJ-3092-25-2H_quant	2	25	Head
UJ-3092-25-3B_quant	3	25	Body
UJ-3092-25-3H_quant	3	25	Head
UJ-3092-30-1B_quant	1	30	Body
UJ-3092-30-1H_quant	1	30	Head
UJ-3092-30-2B_quant	2	30	Body
UJ-3092-30-2H_quant	2	30	Head
UJ-3092-30-3B_quant	3	30	Body
UJ-3092-30-3H_quant	3	30	Head
UJ-3092-36-1B_quant	1	36	Body
UJ-3092-36-1H_quant	1	36	Head
UJ-3092-36-2B_quant	2	36	Body
UJ-3092-36-2H_quant	2	36	Head
UJ-3092-36-3B_quant	3	36	Body
UJ-3092-36-3H_quant	3	36	Head
UJ-3092-40-1B_quant	1	40	Body
UJ-3092-40-1H_quant	1	40	Head
UJ-3092-40-2B_quant	2	40	Body
UJ-3092-40-2H_quant	2	40	Head
UJ-3092-40-3B_quant	3	40	Body
UJ-3092-40-3H_quant	3	40	Head
UJ-3092-48-1B_quant	1	48	Body
UJ-3092-48-1H_quant	1	48	Head
UJ-3092-48-2B_quant	2	48	Body
UJ-3092-48-2H_quant	2	48	Head
UJ-3092-48-3B_quant	3	48	Body
UJ-3092-48-3H_quant	3	48	Head
UJ-3092-Unr-1B_quant	1	Unresponsive	Body
UJ-3092-Unr-1H_quant	1	Unresponsive	Head
UJ-3092-Unr-2B_quant	2	Unresponsive	Body
UJ-3092-Unr-2H_quant	2	Unresponsive	Head
UJ-3092-Unr-3B_quant	3	Unresponsive	Body
UJ-3092-Unr-3H_quant	3	Unresponsive	Head``` </pre>

## Installation / Requirements

- `FastQC/0.12.1-Java-17`
- `Java/17.0.6`
- `Trimmomatic/0.39-Java-17`
- `R` (v4.2.0)
  - Packages: `DESeq2`, `tximport`, `pheatmap`, `ggplot2

## How to run

Run quality control:

```bash
## Step 1: QC with FastQC

# Directory containing the raw data
ls "$raw_data_dir"

# Set Internal Field Separator (IFS) to newline to handle filenames with spaces
IFS=$'\n'

# Find all .fastq files in the raw_data directory and its subdirectories
fastq_files=($(find "$raw_data_dir" -type f \( -iname "*.fastq" -o -iname "*.fastq.gz" -o -iname "*.fq" -o -iname "*.fq.gz" \)))

# Write them to a file
printf "%s\n" "${fastq_files[@]}" > fastq_files

# Output the number of files found
echo "Found ${#fastq_files[@]} FASTQ files."

# Submit the SLURM job array
num_files=${#fastq_files[@]}
echo "Submitting job array for $num_files files..."
sbatch --array=0-$(($num_files - 1)) git_repos/rnaseq_data_scripts/fastqc_array.sh
```
## Assess the quality of reads and presence of adapters based on the output of FastQC

```bash
## Step 2: Trimming with Trimmomatic

# Set raw data directory
export raw_data_dir="temperature_samples/220211_A00181_0425_BHVMJNDSX2/"
ls "$raw_data_dir"

# Use newline separator
IFS=$'\n'

# Find all .fastq files in the raw_data directory and its subdirectories
fastq_files=($(find "$raw_data_dir" -type f \( -iname "*.fastq" -o -iname "*.fastq.gz" -o -iname "*.fq" -o -iname "*.fq.gz" \)))

# Output the number of files found
echo "Found ${#fastq_files[@]} FASTQ files."

# Find all paired .fastq files in the raw_data directory and its subdirectories

# Check for both R1 and R2 FASTQ files with updated pattern
R1_files=($(find "$raw_data_dir" -type f -name '*R1*.fastq.gz'))
R2_files=($(find "$raw_data_dir" -type f -name '*R2*.fastq.gz'))
num_files=${#R1_files[@]}

# Debugging: print R1 and R2 filenames found
echo "Number of R1 files: ${#R1_files[@]}"
echo "Number of R2 files: ${#R2_files[@]}"

# Submit the SLURM job array for paired files
num_files=${#R1_files[@]}
echo "Submitting job array for $num_files paired files..."
sbatch --array=0-$(($num_files - 1)) git_repos/rnaseq_data_scripts/trimmomatic_array.sh

```
## Rerun trimmed reads on FastQC to ensure adapters were removed successfully and read quality is good

```bash
## Step 3: QC with FastQC

# Directory containing the raw data
ls "$raw_data_dir"

# Set Internal Field Separator (IFS) to newline to handle filenames with spaces
IFS=$'\n'

# Find all .fastq files in the raw_data directory and its subdirectories
fastq_files=($(find "$raw_data_dir" -type f \( -iname "*.fastq" -o -iname "*.fastq.gz" -o -iname "*.fq" -o -iname "*.fq.gz" \)))

# Write them to a file
printf "%s\n" "${fastq_files[@]}" > fastq_files

# Output the number of files found
echo "Found ${#fastq_files[@]} FASTQ files."

# Submit the SLURM job array
num_files=${#fastq_files[@]}
echo "Submitting job array for $num_files files..."
sbatch --array=0-$(($num_files - 1)) git_repos/rnaseq_data_scripts/fastqc_array.sh
```

## Needed to trim PolyA/T and Poly G/C tails, standard Illumina recommends trimming 20.

```
sbatch --array=0-$((36 - 1)) git_repos/rnaseq_data_scripts/cutadapt_polya.sh
```
## QC Again to ensure they were removed successfully

## Then run Salmon

```bash
## Step 4: Salmon Alignment

## 4.1. Generating a decoy-aware transcriptome

# Directory containing the reference genome 
ref_genome_dir="/home/ar9416e/git_repos/Reference_genome"

# Directory containing the reference transcripts data
transcripts_dir="/home/ar9416e/git_repos/Reference_genome"

# Extract chromosome information from the genome file to make the decoys.txt file
grep "^>" "$ref_genome_dir/VectorBase-68_AgambiaePEST_Genome.fasta" | cut -d " " -f 1 > decoys.txt

# Check
head decoys.txt

# Remove the > symbol from the file
sed -i.bak -e 's/>//g' decoys.txt

# Check
head decoys.txt

# Combine the transcripts and genome file, in this order!
cat "$transcripts_dir/VectorBase-68_AgambiaePEST_AnnotatedTranscripts.fasta" "$ref_genome_dir/VectorBase-68_AgambiaePEST_Genome.fasta" > trans_and_gen.fa.gz

# Load required modules
module load Salmon/1.10.1-GCC-12.3.0
module load gzip

# Set desired output directory
index_dir="/home/ar9416e/mosquito_test/alignments/salmon_index"

# Create it if it doesn't exist
mkdir -p "$index_dir"

# Run salmon index with decoy-aware transcriptome
salmon index -t trans_and_gen.fa.gz -d decoys.txt -p 30 -i "$index_dir"

#Step 5: Run alignment with Salmon

# Set path to salmon index (make sure the index directory is used)
salmon_index="/home/ar9416e/mosquito_test/alignments/salmon_index"

# Trimmed fastq file folder
TRIM_DATA_DIR="/home/ar9416e/mosquito_test/trimmed_reads"

# Check for both R1 and R2 FASTQ files with updated pattern
R1_files=($(find "$TRIM_DATA_DIR" -type f -name '*R1*_paired.fastq.gz'))
R2_files=($(find "$TRIM_DATA_DIR" -type f -name '*R2*_paired.fastq.gz'))
num_files=${#R1_files[@]}

# Debug: Print the numbers of files found
echo "Number of R1 files: ${#R1_files[@]}"
echo "Number of R2 files: ${#R2_files[@]}"

# Export variables so that they are accessible by the job script
export R1_files
export R2_files
export salmon_index

# You can also export any other variables needed in salmon_alignment.sh, e.g.:
export output_dir="/home/ar9416e/mosquito_test/alignments"

# Submit the SLURM job array for paired files
echo "Submitting job array for $num_files paired files..."
sbatch --array=0-$(($num_files - 1)) /home/ar9416e/git_repos/rnaseq_data_scripts/salmon_alignment.sh

# Get % mapped reads
for d in *_quant; do
  log_file="$d/logs/salmon_quant.log"
  if [[ -f "$log_file" ]]; then
    rate=$(grep "Mapping rate" "$log_file" | awk -F'= ' '{print $2}')
    echo "$d: $rate"
  else
    echo "$d: log file not found"
  fi
done > percentage_mapped
# filepath: alignments/extract_mapping_rates.sh
```
## Assess mapping % to ensure good to use for DESeq

```
UJ-3092-25-1B_quant: 85.5997%
UJ-3092-25-1H_quant: 85.765%
UJ-3092-25-2B_quant: 86.6355%
UJ-3092-25-2H_quant: 87.2109%
UJ-3092-25-3B_quant: 79.1229%
UJ-3092-25-3H_quant: 86.2564%
UJ-3092-30-1B_quant: 87.1759%
UJ-3092-30-1H_quant: 85.6662%
UJ-3092-30-2B_quant: 88.8804%
UJ-3092-30-2H_quant: 84.6123%
UJ-3092-30-3B_quant: 86.185%
UJ-3092-30-3H_quant: 85.3965%
UJ-3092-36-1B_quant: 85.8969%
UJ-3092-36-1H_quant: 87.4188%
UJ-3092-36-2B_quant: 87.0317%
UJ-3092-36-2H_quant: 88.409%
UJ-3092-36-3B_quant: 86.476%
UJ-3092-36-3H_quant: 89.7034%
UJ-3092-40-1B_quant: 88.5589%
UJ-3092-40-1H_quant: 87.5719%
UJ-3092-40-2B_quant: 86.3222%
UJ-3092-40-2H_quant: 85.7973%
UJ-3092-40-3B_quant: 85.2677%
UJ-3092-40-3H_quant: 86.2961%
UJ-3092-48-1B_quant: 85.9177%
UJ-3092-48-1H_quant: 86.6066%
UJ-3092-48-2B_quant: 82.3994%
UJ-3092-48-2H_quant: 87.4921%
UJ-3092-48-3B_quant: 5.23747%
UJ-3092-48-3H_quant: 86.8481%
UJ-3092-Unr-1B_quant: 15.0958%
UJ-3092-Unr-1H_quant: 87.2495%
UJ-3092-Unr-2B_quant: 88.5369%
UJ-3092-Unr-2H_quant: 87.9834%
UJ-3092-Unr-3B_quant: 76.3759%
UJ-3092-Unr-3H_quant: 88.0948%
```
## Two samples below 75% to be removed for DESeq analysis

*Prepare RNA-seq data for differential expression analysis using DESeq2, starting from transcript-level quantification files (from Salmon) and a reference transcriptome in FASTA format.

```
setwd("/Users/ainhoarodriguezpereira/Documents/Mosquito_test")

library(tximport)
library(readr)
library(DESeq2)

# Prepare the dataframe with reference genome gene predictions
fasta <- "/Users/ainhoarodriguezpereira/Documents/Mosquito_test/Reference_genome/VectorBase-68_AgambiaePEST_AnnotatedTranscripts.fasta"
headers <- readLines(fasta)
headers <- headers[grepl("^>", headers)]

# Remove the ">" at the beginning
headers <- sub("^>", "", headers)

# Extract transcript and gene info
tx2gene <- do.call(rbind, lapply(headers, function(h) {
  tx <- strsplit(h, " ")[[1]][1]  # First part = transcript ID
  gene <- sub(".*gene=([^ ]+).*", "\\1", h)  # Extract gene=... value
  c(tx, gene)
}))

# Convert to data.frame
tx2gene_df <- as.data.frame(tx2gene, stringsAsFactors = FALSE)
colnames(tx2gene_df) <- c("TXNAME", "GENEID")

# Save
write.table(tx2gene_df, file = "tx2gene.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Prepare aligned samples for DESeq
samples <- list.files("/Users/ainhoarodriguezpereira/Documents/Mosquito_test/alignments", 
                      pattern = "quant\\.sf$", 
                      recursive = TRUE, 
                      full.names = TRUE)

names(samples) <- basename(dirname(samples))

tx2gene <- read_tsv("tx2gene.tsv", col_names = TRUE)
txi <- tximport(samples, type = "salmon", tx2gene = tx2gene)

# List all sample directories
sample_dirs <- list.files("/Users/ainhoarodriguezpereira/Documents/Mosquito_test/alignments", pattern = "_quant$", full.names = FALSE)

# Extract temperature and tissue
sample_info <- data.frame(
  sampleName = sample_dirs,
  temperature = sub("UJ-3092-(\\d+|Unr)-.*", "\\1", sample_dirs),
  tissue = sub(".*-(\\d+)([BH])_quant$", "\\2", sample_dirs)
)


# Then label B = Body, H = Head
sample_info$tissue <- ifelse(sample_info$tissue == "B", "Body", "Head")

# Remove bad samples
bad_samples <- c("UJ-3092-48-3B_quant", "UJ-3092-Unr-1B_quant")
sample_info <- sample_info[!sample_info$sampleName %in% bad_samples, ]

# Filter samples vector to match cleaned sample_info
samples <- samples[names(samples) %in% sample_info$sampleName]

# Ensure the order of sample_info matches samples
sample_info <- sample_info[match(names(samples), sample_info$sampleName), ]

# Save to CSV
write.csv(sample_info, "samples.csv", row.names = FALSE)
```
** Plot PCA

```

# List quant.sf files
samples <- list.files("alignments", pattern = "quant\\.sf$", recursive = TRUE, full.names = TRUE)
names(samples) <- basename(dirname(samples))

# Remove bad samples before tximport
bad_samples <- c("UJ-3092-48-3B_quant", "UJ-3092-Unr-1B_quant")
samples <- samples[!names(samples) %in% bad_samples]

# Re-run tximport with only the good samples
txi <- tximport(samples, type = "salmon", tx2gene = tx2gene)

colnames(txi$counts)

# Now create DESeq2 object
dds <- DESeqDataSetFromTximport(txi, colData = sample_info, design = ~ temperature * tissue)

colnames(dds)

vsd <- vst(dds, blind = FALSE)  # Use blind=FALSE if you already defined your design

plotPCA(vsd, intgroup = c("temperature", "tissue"))

library(ggplot2)

# Get PCA data
pcaData <- plotPCA(vsd, intgroup = c("temperature", "tissue"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

library(DESeq2)
library(ggplot2)
library(scales)

# Get PCA data from DESeq2 object
pcaData <- plotPCA(vsd, intgroup = c("temperature", "tissue"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Ensure temperature and tissue are factors with correct order and labels
pcaData$temperature <- factor(pcaData$temperature,
                              levels = c("Unr", "25", "30", "36", "40", "48"),
                              labels = c("Unresponsive", "25°C", "30°C", "36°C", "40°C", "48°C"))

pcaData$tissue <- factor(pcaData$tissue,
                         levels = c("Body", "Head"))

# Custom ggplot theme for publication
theme_pub <- function(base_size = 14, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15)),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12, color = "black"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.position = "right",
      panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray30", linewidth = 0.5),
      plot.margin = margin(15, 15, 15, 15)
    )
}

# PCA plot
plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = temperature, shape = tissue)) +
  geom_point(size = 4, alpha = 0.8, stroke = 0.8) +
  xlab(paste0("PC1 (", percentVar[1], "% variance)")) +
  ylab(paste0("PC2 (", percentVar[2], "% variance)")) +
  theme_pub() +
  scale_color_brewer(palette = "Set1", name = "Temperature") +
  scale_shape_manual(name = "Tissue", values = c("Body" = 16, "Head" = 17)) +
  labs(title = "PCA of gene expression by tissue and temperature") +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 4)),
    shape = guide_legend(order = 2, override.aes = list(size = 4))
  ) +
  coord_fixed(ratio = 1)

ggsave("PCA_plot.pdf", plot, width = 8, height = 6, dpi = 300)

# Optional: Scree plot of variance explained
pca_res <- prcomp(t(assay(vsd)))
variance_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)
scree_data <- data.frame(
  PC = factor(paste0("PC", seq_along(variance_explained)), levels = paste0("PC", seq_along(variance_explained))),
  Variance = variance_explained * 100
)


scree_plot <- ggplot(scree_data[1:10, ], aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("Principal Components") +
  ylab("Variance Explained (%)") +
  ggtitle("Scree Plot of PCA") +
  theme_pub()

# Save to file
ggsave("Scree_plot.pdf", scree_plot, width = 8, height = 6, dpi = 300)
```




