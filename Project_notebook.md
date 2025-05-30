## FastQC and Trimming

## Raw Read FastQC

Ran FastQC from raw read data, results:

```
FastQC Overview		Comments:
UJ-3092-25-1B_quant	R1	Presence adapters + Poly G
	R2	Presence adapters + Poly G
UJ-3092-25-1H_quant	R1	Presence adapters + Poly G
	R2	Presence adapters + Poly G
UJ-3092-25-2B_quant	R1	Presence adapters
	R2	Presence adapters + Poly G
UJ-3092-25-2H_quant	R1	Presence adapters
	R2	Presence adapters + Poly G
UJ-3092-25-3B_quant	R1	Presence adapters + Poly A
	R2	Presence adapters + Poly T
UJ-3092-25-3H_quant	R1	Presence adapters
	R2	Presence adapters
UJ-3092-30-1B_quant	R1	Presence adapters
	R2	Presence adapters + Poly G
UJ-3092-30-1H_quant	R1	Presence adapters
	R2	Presence adapters + Poly G
UJ-3092-30-2B_quant	R1	Presence adapters + Poly A
	R2	Presence adapters + Poly T
UJ-3092-30-2H_quant	R1	Presence adapters
	R2	Presence adapters
UJ-3092-30-3B_quant	R1	Presence adapters + Poly A
	R2	Presence adapters + Poly T
UJ-3092-30-3H_quant	R1	Presence adapters
	R2	Presence adapters + Poly G
UJ-3092-36-1B_quant	R1	Presence adapters
	R2	Presence adapters + Poly G
UJ-3092-36-1H_quant	R1	Presence adapters + Poly A
	R2	Presence adapters + Poly T
UJ-3092-36-2B_quant	R1	Presence adapters + Poly A
	R2	Presence adapters + Poly T
UJ-3092-36-2H_quant	R1	Presence adapters
	R2	Presence adapters
UJ-3092-36-3B_quant	R1	Presence adapters
	R2	Presence adapters + Poly T
UJ-3092-36-3H_quant	R1	Presence adapters
	R2	Presence adapters + Poly G + Poly T
UJ-3092-40-1B_quant	R1	Presence adapters
	R2	Presence adapters + Poly T
UJ-3092-40-1H_quant	R1	Presence adapters + Poly A
	R2	Presence adapters + Poly T
UJ-3092-40-2B_quant	R1	Presence adapters
	R2	Presence adapters + Poly G
UJ-3092-40-2H_quant	R1	Presence adapters
	R2	Presence adapters
UJ-3092-40-3B_quant	R1	Presence adapters + Poly A
	R2	Presence adapters + Poly T
UJ-3092-40-3H_quant	R1	Presence adapters + Poly A
	R2	Presence adapters + Poly T
UJ-3092-48-1B_quant	R1	Presence adapters + Poly G + Poly T
	R2	Presence adapters + Poly A
UJ-3092-48-1H_quant	R1	Presence adapters + Poly G + Poly T
	R2	Presence adapters
UJ-3092-48-2B_quant	R1	Presence adapters + Poly G + Poly T
	R2	Presence adapters + Poly A
UJ-3092-48-2H_quant	R1	Presence adapters + Poly G + Poly T
	R2	Presence adapters + Poly A
UJ-3092-48-3H_quant	R1	Presence adapters + Poly G + Poly T
	R2	Presence adapters + Poly A
UJ-3092-Unr-1H_quant	R1	Presence adapters, Coronavirus top hit most represented sequence, to trim
	R2	Presence adapters
UJ-3092-Unr-2B_quant	R1	Presence adapters + Poly G + Poly T
	R2	Presence adapters
UJ-3092-Unr-2H_quant	R1	Presence adapters
	R2	Presence adapters + Poly G + Poly T
UJ-3092-Unr-3B_quant	R1	Presence adapters
	R2	Presence adapters + Poly G + Poly T
UJ-3092-Unr-3H_quant	R1	Presence adapters + Poly A
	R2	Presence adapters + Poly G + Poly T
```
### Outcome:

 Adapter and Overrepresented Sequences
- **Adapter detected:** TruSeq Adapter, Index 6 (removed using Trimmomatic)
- **Most abundant overrepresented sequences in most reads:** Poly A/T and Poly C/G sequences found, as per Illumina guidelines 20 nt per string were removed with CutAdapt.
- **Removed Coronavirus contamination from one sample** -> 

![Screenshot 2025-05-30 at 22 24 25](https://github.com/user-attachments/assets/4aee29cb-fce0-444d-a24f-241d71e2813d)


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



