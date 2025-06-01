#Run DESeq -Unresponsive vs temperatures for heads

# ------------------ Run DESeq2: Head samples only ------------------

# Setup
rownames(sample_info) <- colnames(txi$counts)
bad_samples <- c("UJ-3092-48-3B_quant", "UJ-3092-Unr-1B_quant")

# Subset to 'Head' samples only, excluding bad ones
head_samples <- sample_info$tissue == "Head" & !(rownames(sample_info) %in% bad_samples)
sample_info_head <- sample_info[head_samples, ]
sample_info_head$temperature <- factor(sample_info_head$temperature,
                                       levels = c("Unr", "25", "30", "36", "40", "48"))

# Subset txi object
txi_head <- list(
  counts = txi$counts[, head_samples],
  abundance = txi$abundance[, head_samples],
  length = txi$length[, head_samples]
)
if (!is.null(txi$countsFromAbundance)) {
  txi_head$countsFromAbundance <- txi$countsFromAbundance
}
stopifnot(all(colnames(txi_head$counts) == rownames(sample_info_head)))

# Create DESeq2 object and run DESeq
dds_head <- DESeqDataSetFromTximport(txi_head, colData = sample_info_head, design = ~ temperature)
dds_head <- DESeq(dds_head)

# Full DE results
res_head <- results(dds_head)
sig_genes_head <- res_head[which(res_head$padj < 0.05), ]
sig_genes_head <- sig_genes_head[order(sig_genes_head$padj), ]

cat("Number of significant DE genes in Head:", sum(res_head$padj < 0.05, na.rm = TRUE), "\n")

# Save DE results
write.csv(as.data.frame(sig_genes_head), "significant_DE_genes_head.csv")
write.csv(as.data.frame(res_head), "all_DE_results_head.csv")

# ------------------ Contrasts: Each temperature vs Unr ------------------

temps <- c("25", "30", "36", "40", "48")

for (temp in temps) {
  contrast_name <- paste0("temperature_", temp, "_vs_Unr")
  
  res <- results(dds_head, name = contrast_name)
  sig <- res[which(res$padj < 0.05), ]
  sig <- sig[order(sig$padj), ]
  
  write.csv(as.data.frame(sig), paste0("sig_genes_head_", temp, "_vs_Unr.csv"))
  write.csv(as.data.frame(res), paste0("all_genes_head_", temp, "_vs_Unr.csv"))
  
  cat("Contrast", contrast_name, ":", nrow(sig), "significant genes (padj < 0.05)\n")
}

# ------------------ PCA ------------------

vsd_head <- vst(dds_head, blind = TRUE)
pca_res <- prcomp(t(assay(vsd_head)))
percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2))[1:2], 1)

pcaData <- as.data.frame(pca_res$x[, 1:2])
pcaData$temperature <- colData(vsd_head)$temperature
pcaData$tissue <- colData(vsd_head)$tissue

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
  labs(title = "PCA of gene expression of Head Samples") +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 4)),
    shape = guide_legend(order = 2, override.aes = list(size = 4))
  ) +
  coord_fixed(ratio = 1)

ggsave("PCA_plot_head.pdf", plot, width = 8, height = 6, dpi = 300)

# Scree plot
variance_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)

png("Scree_plot_head.png", width = 800, height = 600, res = 150)
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

# ------------------ Collapse into "Unr" vs "Responsive" ------------------

sample_info_head$cond <- ifelse(sample_info_head$temperature == "Unr", "Unr", "Responsive")
sample_info_head$cond <- factor(sample_info_head$cond, levels = c("Unr", "Responsive"))

dds_head_cond <- DESeqDataSetFromTximport(txi_head, colData = sample_info_head, design = ~ cond)
dds_head_cond <- DESeq(dds_head_cond)

res_unr_vs_responsive_head <- results(dds_head_cond, contrast = c("cond", "Responsive", "Unr"))

sig_unr_vs_responsive_head <- res_unr_vs_responsive_head[which(res_unr_vs_responsive_head$padj < 0.05), ]
sig_unr_vs_responsive_head <- sig_unr_vs_responsive_head[order(sig_unr_vs_responsive_head$padj), ]

write.csv(as.data.frame(res_unr_vs_responsive_head), "all_genes_head_Unr_vs_Responsive.csv")
write.csv(as.data.frame(sig_unr_vs_responsive_head), "sig_genes_head_Unr_vs_Responsive.csv")

cat("Unr vs Responsive (Head): Found", nrow(sig_unr_vs_responsive_head), "significant genes (padj < 0.05)\n")


######### 36 comparisons ##########

# Define contrast pairs: 36°C vs other temps
compare_temps <- c("25", "30", "40", "48")

for (temp in compare_temps) {
  contrast_name <- paste0("36_vs_", temp)
  
  # Run contrast
  res <- results(dds_head, contrast = c("temperature", "36", temp))
  
  # Filter significant genes
  sig <- res[which(res$padj < 0.05), ]
  sig <- sig[order(sig$padj), ]
  
  # Save results
  write.csv(as.data.frame(sig), paste0("sig_genes_head_36_vs_", temp, ".csv"))
  write.csv(as.data.frame(res), paste0("all_genes_head_36_vs_", temp, ".csv"))
  
  cat("Contrast 36 vs", temp, ":", nrow(sig), "significant genes (padj < 0.05)\n")
}

# Create a new condition in sample_info_head: "36" vs "Other"
sample_info_head$cond <- ifelse(sample_info_head$temperature == "36", "36", "Other")
sample_info_head$cond <- factor(sample_info_head$cond, levels = c("36", "Other"))

# Re-create DESeq2 object with the new condition
dds_head_cond <- DESeqDataSetFromTximport(txi_head, colData = sample_info_head, design = ~ cond)
dds_head_cond <- DESeq(dds_head_cond)

# Extract results: 36 vs Other
res_36_vs_others <- results(dds_head_cond, contrast = c("cond", "36", "Other"))

# Filter and sort significant DE genes
sig_36_vs_others <- res_36_vs_others[which(res_36_vs_others$padj < 0.05), ]
sig_36_vs_others <- sig_36_vs_others[order(sig_36_vs_others$padj), ]

# Save results
write.csv(as.data.frame(res_36_vs_others), "all_genes_head_36_vs_Others.csv")
write.csv(as.data.frame(sig_36_vs_others), "sig_genes_head_36_vs_Others.csv")

# Print summary
cat("36 vs Other Temps (Heads):", nrow(sig_36_vs_others), "significant genes (padj < 0.05)\n")

