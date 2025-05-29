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


