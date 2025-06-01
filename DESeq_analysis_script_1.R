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
