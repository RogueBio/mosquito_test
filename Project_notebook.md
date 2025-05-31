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
UJ-3092-Unr-1H_quant	R1	Presence adapters
	R2	Presence adapters
UJ-3092-Unr-2B_quant	R1	Presence adapters
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

**Note:** Removing overrepresented sequences corresponding to real biological transcripts (e.g., highly expressed genes) can bias results. Thus, these sequences were retained to avoid introducing unintended biases.

## Post-trimming FastQC

Samples still showed overrepresented sequences; they were blasted individually per sample, and none of the sequences blasted anything other than Anopheles sequences, so these were not further trimmed. Two samples were very small in size after trimming, these were further investigated.

- **UJ-3092-48-3B_quant:**
- **UJ-3092-Unr-1B_quant:**

Changed the parameters of Trimming to allow for more leniency, but the outcome of the trimming was still the same: Cutadapt log showed a high percentage of reads with adapter contamination

*** UJ-3092-48-3B**

Ran 
```
cutadapt \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -a A{20} -a T{20} -a G{20} -a C{20} \
  -A A{20} -A T{20} -A G{20} -A C{20} \
  -q 10,10 \
  --minimum-length 20 \
  --pair-filter=both \
  -o /home/ar9416e/mosquito_test/single_trim_reads_polyA/UJ-3092-48-3B_R1.fastq.gz \
  -p /home/ar9416e/mosquito_test/single_trim_reads_polyA/UJ-3092-48-3B_R2.fastq.gz \
  /home/ar9416e/temperature_samples/220211_A00181_0425_BHVMJNDSX2/Sample_UJ-3092-48-3B/UJ-3092-48-3B_S30_L001_R1_001.fastq.gz \
  /home/ar9416e/temperature_samples/220211_A00181_0425_BHVMJNDSX2/Sample_UJ-3092-48-3B/UJ-3092-48-3B_S30_L001_R2_001.fastq.gz \
  > /home/ar9416e/mosquito_test/single_trim_reads_polyA/UJ-3092-48-3B_cutadapt.log

```
And output of the log was:

```
This is cutadapt 4.2 with Python 3.10.4
Command line parameters: -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -a A{20} -a T{20} -a G{20} -a C{20} -A A{20} -A T{20} -A G{20} -A C{20} -q 10,10 --minimum-length 20 --pair-filter=both -o /home/ar9416e/mosquito_test/single_trim_reads_polyA/UJ-3092-48-3B_R1.fastq.gz -p /home/ar9416e/mosquito_test/single_trim_reads_polyA/UJ-3092-48-3B_R2.fastq.gz /home/ar9416e/temperature_samples/220211_A00181_0425_BHVMJNDSX2/Sample_UJ-3092-48-3B/UJ-3092-48-3B_S30_L001_R1_001.fastq.gz /home/ar9416e/temperature_samples/220211_A00181_0425_BHVMJNDSX2/Sample_UJ-3092-48-3B/UJ-3092-48-3B_S30_L001_R2_001.fastq.gz
Processing paired-end reads on 1 core ...
Finished in 51.101 s (35.712 µs/read; 1.68 M reads/minute).

=== Summary ===

Total read pairs processed:          1,430,926
  Read 1 with adapter:               1,414,427 (98.8%)
  Read 2 with adapter:               1,414,361 (98.8%)

== Read fate breakdown ==
Pairs that were too short:           1,357,961 (94.9%)
Pairs written (passing filters):        72,965 (5.1%)

Total basepairs processed:   432,139,652 bp
  Read 1:   216,069,826 bp
  Read 2:   216,069,826 bp
Quality-trimmed:                     993 bp (0.0%)
  Read 1:           194 bp
  Read 2:           799 bp
Total written (filtered):     14,309,872 bp (3.3%)
  Read 1:     6,391,006 bp
  Read 2:     7,918,866 bp
```

*** UJ-3092-Unr-1B**

Ran 
```
cutadapt \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -a A{20} -a T{20} -a G{20} -a C{20} \
  -A A{20} -A T{20} -A G{20} -A C{20} \
  -q 10,10 \
  --minimum-length 20 \
  --pair-filter=both \
  -o /home/ar9416e/mosquito_test/single_trim_reads_polyA/UJ-3092-Unr-1B_R1.fastq.gz \
  -p /home/ar9416e/mosquito_test/single_trim_reads_polyA/UJ-3092-Unr-1B_R2.fastq.gz \
  /home/ar9416e/temperature_samples/220211_A00181_0425_BHVMJNDSX2/Sample_UJ-3092-Unr-1B/UJ-3092-Unr-1B_S34_L001_R1_001.fastq.gz \
  /home/ar9416e/temperature_samples/220211_A00181_0425_BHVMJNDSX2/Sample_UJ-3092-Unr-1B/UJ-3092-Unr-1B_S34_L001_R2_001.fastq.gz \
  > /home/ar9416e/mosquito_test/single_trim_reads_polyA/UJ-3092-Unr-1B_cutadapt.log

```
And the log showed:

```
This is cutadapt 4.2 with Python 3.10.4
Command line parameters: -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -a A{20} -a T{20} -a G{20} -a C{20} -A A{20} -A T{20} -A G{20} -A C{20} -q 10,10 --minimum-length 20 --pair-filter=both -o /home/ar9416e/mosquito_test/single_trim_reads_polyA/UJ-3092-Unr-1B_R1.fastq.gz -p /home/ar9416e/mosquito_test/single_trim_reads_polyA/UJ-3092-Unr-1B_R2.fastq.gz /home/ar9416e/temperature_samples/220211_A00181_0425_BHVMJNDSX2/Sample_UJ-3092-Unr-1B/UJ-3092-Unr-1B_S34_L001_R1_001.fastq.gz /home/ar9416e/temperature_samples/220211_A00181_0425_BHVMJNDSX2/Sample_UJ-3092-Unr-1B/UJ-3092-Unr-1B_S34_L001_R2_001.fastq.gz
Processing paired-end reads on 1 core ...
Finished in 120.911 s (49.542 µs/read; 1.21 M reads/minute).

=== Summary ===

Total read pairs processed:          2,440,590
  Read 1 with adapter:               2,277,445 (93.3%)
  Read 2 with adapter:               2,278,184 (93.3%)

== Read fate breakdown ==
Pairs that were too short:           1,554,410 (63.7%)
Pairs written (passing filters):       886,180 (36.3%)

Total basepairs processed:   737,058,180 bp
  Read 1:   368,529,090 bp
  Read 2:   368,529,090 bp
Quality-trimmed:                   1,793 bp (0.0%)
  Read 1:           318 bp
  Read 2:         1,475 bp
Total written (filtered):    197,882,275 bp (26.8%)
  Read 1:    96,979,120 bp
  Read 2:   100,903,155 bp

```

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



