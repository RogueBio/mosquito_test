
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
