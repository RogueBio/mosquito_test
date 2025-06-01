
#Run DESeq -Unresponsive vs temperatures for bodies
# Subset only 'Body' samples and exclude bad ones
rownames(sample_info) <- colnames(txi$counts)
bad_samples <- c("UJ-3092-48-3B_quant", "UJ-3092-Unr-1B_quant")

body_samples <- sample_info$tissue == "Body" &
  !(rownames(sample_info) %in% bad_samples)

sample_info_body <- sample_info[body_samples, ]
sample_info_body$temperature <- factor(sample_info_body$temperature,
                                       levels = c("Unr", "25", "30", "36", "40", "48"))

txi_body <- list(
  counts = txi$counts[, body_samples],
  abundance = txi$abundance[, body_samples],
  length = txi$length[, body_samples]
)

if (!is.null(txi$countsFromAbundance)) {
  txi_body$countsFromAbundance <- txi$countsFromAbundance
}

stopifnot(all(colnames(txi_body$counts) == rownames(sample_info_body)))

# Run DESeq2
dds_body <- DESeqDataSetFromTximport(txi_body, colData = sample_info_body, design = ~ temperature)
dds_body <- DESeq(dds_body)

# Full results
res_body <- results(dds_body)
sig_genes_body <- res_body[which(res_body$padj < 0.05), ]
sig_genes_body <- sig_genes_body[order(sig_genes_body$padj), ]

cat("Number of significant DE genes in Body:", sum(res_body$padj < 0.05, na.rm = TRUE), "\n")

write.csv(as.data.frame(sig_genes_body), "significant_DE_genes_body.csv")
write.csv(as.data.frame(res_body), "all_DE_results_body.csv")

# Export each contrast vs Unr
temps <- c("25", "30", "36", "40", "48")

for (temp in temps) {
  contrast_name <- paste0("temperature_", temp, "_vs_Unr")
  
  res <- results(dds_body, name = contrast_name)
  
  sig <- res[which(res$padj < 0.05), ]
  sig <- sig[order(sig$padj), ]
  
  write.csv(as.data.frame(sig), paste0("sig_genes_", temp, "_vs_Unr.csv"))
  write.csv(as.data.frame(res), paste0("all_genes_", temp, "_vs_Unr.csv"))
  
  cat("Contrast", contrast_name, ":", nrow(sig), "significant genes (padj < 0.05)\n")
}

# PCA on vst-transformed data
vsd <- vst(dds_body, blind = TRUE)
pca_res <- prcomp(t(assay(vsd)))
percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2))[1:2], 1)
pcaData <- as.data.frame(pca_res$x[, 1:2])
pcaData$temperature <- colData(vsd)$temperature
pcaData$tissue <- colData(vsd)$tissue

# Custom ggplot theme
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

ggsave("PCA_plot_body.pdf", plot, width = 8, height = 6, dpi = 300)

# Scree plot
variance_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)

png("Scree_plot_body.png", width = 800, height = 600, res = 150)
barplot(
  variance_explained * 100,
  names.arg = paste0("PC", seq_along(variance_explained)),
  col = "steelblue",
  main = "Scree Plot",
  ylab = "Variance Explained (%)",
  xlab = "Principal Components",
  las = 2
)
dev.off()

# Create a new condition: "Unr" vs "Responsive"
sample_info_body$cond <- ifelse(sample_info_body$temperature == "Unr", "Unr", "Responsive")
sample_info_body$cond <- factor(sample_info_body$cond, levels = c("Unr", "Responsive"))

# Re-create DESeq2 object with new condition
dds_body_cond <- DESeqDataSetFromTximport(txi_body, colData = sample_info_body, design = ~ cond)
dds_body_cond <- DESeq(dds_body_cond)

# Run results for Responsive vs Unr
res_unr_vs_responsive <- results(dds_body_cond, contrast = c("cond", "Responsive", "Unr"))

# Filter and sort significant DE genes
sig_unr_vs_responsive <- res_unr_vs_responsive[which(res_unr_vs_responsive$padj < 0.05), ]
sig_unr_vs_responsive <- sig_unr_vs_responsive[order(sig_unr_vs_responsive$padj), ]

# Save results to CSV
write.csv(as.data.frame(res_unr_vs_responsive), "all_genes_Unr_vs_Responsive.csv")
write.csv(as.data.frame(sig_unr_vs_responsive), "sig_genes_Unr_vs_Responsive.csv")

# Print summary
cat("Unr vs Responsive: Found", nrow(sig_unr_vs_responsive), "significant genes (padj < 0.05)\n")
