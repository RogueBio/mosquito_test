# Project Notebook

## FastQC and Trimming

### Adapter and Overrepresented Sequences
- **Adapter detected:** TruSeq Adapter, Index 6 (removed using Trimmomatic)
- **Most abundant overrepresented sequence:** 50-string A sequence (~1.5% of read count)
- **GC Content Before Trimming:** Peaks at ~35% and ~51%
- **GC Content After Trimming the overrepresented sequence:** Improved after removing the 50-string A sequence, but the decision to remove overrepresented sequences remains debatable.

**Note:** Removing overrepresented sequences corresponding to real biological transcripts (e.g., highly expressed genes) can bias results. Thus, these sequences were retained to avoid introducing unintended biases.

## Post-trimming FastQC

_(Include summary of FastQC results here, such as per-base sequence quality improvements, overall sequence length distributions, and any residual artifacts.)_

## Alignment and Mapping

### Low Mapping Samples Removed
- **UJ-3092-48-3B_quant:** 5.24% mapped reads
- **UJ-3092-Unr-1B_quant:** 15.10% mapped reads

Since these samples had extremely low mapping percentages, they were removed from further analysis to ensure reliable downstream results. This was justified as DESeq2 models gene expression dispersion across replicates. Removing samples with poor mapping rates might slightly shift dispersion estimates, but DESeq2 adjusts for this using empirical Bayes shrinkage, so I decided to remove these two reps, compromising on statistical power.

To assess whether low-mapping-quality samples distorted overall gene expression patterns, we performed a differential expression analysis and dimensionality reduction with and without those samples.

Step 1: Filtering Poor-Quality Samples
I identified and removed samples with abnormally low mapping rates. These samples could introduce bias or noise in downstream analyses, particularly in PCA and differential expression.

Step 2: Modeling Gene Expression
I used DESeq2 to fit a full interaction model of the form:
design=∼temperature∗tissue

Step 3: PCA for Expression Pattern Visualization
Using the variance-stabilized data from DESeq2 (RNA-seq gene count data is not normally distributed, thus it must be transformed for PCA).
I performed Principal Component Analysis (PCA). This unsupervised dimensionality reduction technique helped visualize global expression patterns and assess whether the removal of low-quality samples altered clustering or separation between conditions.

** Before removing samples:
<img width="846" alt="Screenshot 2025-05-30 at 10 28 11" src="https://github.com/user-attachments/assets/e8e380d8-33fb-40f8-abf7-d108cdac3b30" />
![Scree_plot](https://github.com/user-attachments/assets/c6fd7529-6577-43fd-be47-040cce0e6bb9)

**After removing samples:
<img width="902" alt="Screenshot 2025-05-30 at 10 29 56" src="https://github.com/user-attachments/assets/19c95cec-ac54-4261-9bbf-648831af8344" />
![Scree_plot](https://github.com/user-attachments/assets/5381e100-6e96-4463-9789-67e391df3c08)

Based on the PCA plots, removing those low-mapping samples seems to result in cleaner clustering and a more interpretable separation of conditions. Including them appears to introduce more variability, potentially masking true biological patterns. As mapping rates were low, I decided to remove them altogether. 



