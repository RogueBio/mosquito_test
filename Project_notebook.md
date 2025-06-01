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
- Trimming strategies, removing Illumina adapters, polyA/T/G/C tails (≥20 nt), and low-quality bases (Q<10), ensuring a minimum length of 20 bp post-trim.

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
Here we can see that 98.8% of the reads had adapter contamination, and 94.9% were too short, only saving 26.8% of the reads for further alignment.

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
Here we can see that 93.3% of the reads had adapter contamination, and 63.7% were too short, only saving 26.8% of the reads for further alignment.

For contrast, an example of one of the samples that trimmed successfully:
```
Processing UJ-3092-48-3B_trimmed
This is cutadapt 4.2 with Python 3.10.4
Command line parameters: -a A{15} -A T{15} -o /home/ar9416e/mosquito_test/trimmed_reads_polyA/UJ-3092-48-3B_trimmed_R1_polyAT_trimmed.fastq.gz -p /home/ar9416e/mosquito_test/trimmed_reads_polyA/UJ-3092-48-3B_trimmed_R2_polyAT_trimmed.fastq.gz /home/ar9416e/mosquito_test/trimmed_reads/UJ-3092-48-3B_trimmed_R1_paired.fastq.gz /home/ar9416e/mosquito_test/trimmed_reads/UJ-3092-48-3B_trimmed_R2_paired.fastq.gz
Processing paired-end reads on 1 core ...
Finished in 1.456 s (27.869 µs/read; 2.15 M reads/minute).

=== Summary ===

Total read pairs processed:             52,239
  Read 1 with adapter:                   2,787 (5.3%)
  Read 2 with adapter:                   2,197 (4.2%)
Pairs written (passing filters):        52,239 (100.0%)

Total basepairs processed:    12,933,027 bp
  Read 1:     6,475,733 bp
  Read 2:     6,457,294 bp
Total written (filtered):     12,633,737 bp (97.7%)
  Read 1:     6,296,404 bp
  Read 2:     6,337,333 bp

```
Where only 4-5% of reads contained adapter sequences, thus it left 97.7% of the reads intact for alignment.

## Alignment and Mapping
Used validated Salmon index with --validateMappings, ensuring alignment rigor.

Ran Salmon Alignment for all samples, then generated a text file with the read mapping to assess successful or unsuccessful mapping.

```
UJ-3092-25-1H_trimmed_cut_R2_paired.fastq.gz_quant: 85.7753%
UJ-3092-25-2B_trimmed_cut_R2_paired.fastq.gz_quant: 86.6658%
UJ-3092-25-2H_trimmed_cut_R2_paired.fastq.gz_quant: 87.2251%
UJ-3092-25-3B_trimmed_cut_R2_paired.fastq.gz_quant: 79.1462%
UJ-3092-25-3H_trimmed_cut_R2_paired.fastq.gz_quant: 86.2708%
UJ-3092-30-1B_trimmed_cut_R2_paired.fastq.gz_quant: 87.1962%
UJ-3092-30-1H_trimmed_cut_R2_paired.fastq.gz_quant: 85.6791%
UJ-3092-30-2B_trimmed_cut_R2_paired.fastq.gz_quant: 88.8959%
UJ-3092-30-2H_trimmed_cut_R2_paired.fastq.gz_quant: 84.6324%
UJ-3092-30-3B_trimmed_cut_R2_paired.fastq.gz_quant: 86.206%
UJ-3092-30-3H_trimmed_cut_R2_paired.fastq.gz_quant: 85.4109%
UJ-3092-36-1B_trimmed_cut_R2_paired.fastq.gz_quant: 85.9194%
UJ-3092-36-1H_trimmed_cut_R2_paired.fastq.gz_quant: 87.4383%
UJ-3092-36-2B_trimmed_cut_R2_paired.fastq.gz_quant: 87.0504%
UJ-3092-36-2H_trimmed_cut_R2_paired.fastq.gz_quant: 88.4271%
UJ-3092-36-3B_trimmed_cut_R2_paired.fastq.gz_quant: 86.4982%
UJ-3092-36-3H_trimmed_cut_R2_paired.fastq.gz_quant: 89.7167%
UJ-3092-40-1B_trimmed_cut_R2_paired.fastq.gz_quant: 88.5846%
UJ-3092-40-1H_trimmed_cut_R2_paired.fastq.gz_quant: 87.5855%
UJ-3092-40-2B_trimmed_cut_R2_paired.fastq.gz_quant: 86.3431%
UJ-3092-40-2H_trimmed_cut_R2_paired.fastq.gz_quant: 85.8156%
UJ-3092-40-3B_trimmed_cut_R2_paired.fastq.gz_quant: 85.2924%
UJ-3092-40-3H_trimmed_cut_R2_paired.fastq.gz_quant: 86.3177%
UJ-3092-48-1B_trimmed_cut_R2_paired.fastq.gz_quant: 85.9451%
UJ-3092-48-1H_trimmed_cut_R2_paired.fastq.gz_quant: 86.6183%
UJ-3092-48-2B_trimmed_cut_R2_paired.fastq.gz_quant: 82.4108%
UJ-3092-48-2H_trimmed_cut_R2_paired.fastq.gz_quant: 87.5016%
UJ-3092-48-3B_trimmed_cut_R2_paired.fastq.gz_quant: 5.21449%
UJ-3092-48-3H_trimmed_cut_R2_paired.fastq.gz_quant: 86.8625%
UJ-3092-Unr-1B_trimmed_cut_R2_paired.fastq.gz_quant: 15.0974%
UJ-3092-Unr-1H_trimmed_cut_R2_paired.fastq.gz_quant: 87.2626%
UJ-3092-Unr-2B_trimmed_cut_R2_paired.fastq.gz_quant: 88.5424%
UJ-3092-Unr-2H_trimmed_cut_R2_paired.fastq.gz_quant: 87.9967%
UJ-3092-Unr-3B_trimmed_cut_R2_paired.fastq.gz_quant: 76.3877%
UJ-3092-Unr-3H_trimmed_cut_R2_paired.fastq.gz_quant: 88.1059%

```
### Two samples showed low mapping
- **UJ-3092-48-3B:** 5.24% mapped reads
- **UJ-3092-Unr-1B:** 15.10% mapped reads

Salmon logs for each of the low mapping samples:
**UJ-3092-48-3B:**
```
[2025-05-31 12:11:23.273] [jointLog] [info] setting maxHashResizeThreads to 6
[2025-05-31 12:11:23.273] [jointLog] [info] Fragment incompatibility prior below threshold.  Incompatible fragments will be ignored.
[2025-05-31 12:11:23.273] [jointLog] [info] Usage of --validateMappings implies use of minScoreFraction. Since not explicitly specified, it is being set to 0.65
[2025-05-31 12:11:23.273] [jointLog] [info] Setting consensusSlack to selective-alignment default of 0.35.
[2025-05-31 12:11:23.273] [jointLog] [info] parsing read library format
[2025-05-31 12:11:23.273] [jointLog] [info] There is 1 library.
[2025-05-31 12:11:23.274] [jointLog] [info] Loading pufferfish index
[2025-05-31 12:11:23.274] [jointLog] [info] Loading dense pufferfish index.
[2025-05-31 12:11:28.546] [jointLog] [info] done
[2025-05-31 12:11:28.570] [jointLog] [info] Index contained 15,859 targets
[2025-05-31 12:11:28.573] [jointLog] [info] Number of decoys : 153
[2025-05-31 12:11:28.573] [jointLog] [info] First decoy index : 15,706 
[2025-05-31 12:11:29.176] [jointLog] [info] Computed 227 rich equivalence classes for further processing
[2025-05-31 12:11:29.176] [jointLog] [info] Counted 2,724 total reads in the equivalence classes 
[2025-05-31 12:11:29.180] [jointLog] [warning] 2.72593% of fragments were shorter than the k used to build the index.
If this fraction is too large, consider re-building the index with a smaller k.
The minimum read size found was 0.


[2025-05-31 12:11:29.180] [jointLog] [info] Number of mappings discarded because of alignment score : 148,842
[2025-05-31 12:11:29.180] [jointLog] [info] Number of fragments entirely discarded because of alignment score : 4,439
[2025-05-31 12:11:29.180] [jointLog] [info] Number of fragments discarded because they are best-mapped to decoys : 104
[2025-05-31 12:11:29.180] [jointLog] [info] Number of fragments discarded because they have only dovetail (discordant) mappings to valid targets : 782
[2025-05-31 12:11:29.180] [jointLog] [warning] Only 2724 fragments were mapped, but the number of burn-in fragments was set to 5000000.
The effective lengths have been computed using the observed mappings.

[2025-05-31 12:11:29.180] [jointLog] [info] Mapping rate = 5.21449%

[2025-05-31 12:11:29.180] [jointLog] [info] finished quantifyLibrary()
[2025-05-31 12:11:29.182] [jointLog] [info] Starting optimizer
[2025-05-31 12:11:29.200] [jointLog] [info] Marked 0 weighted equivalence classes as degenerate
[2025-05-31 12:11:29.201] [jointLog] [info] iteration = 0 | max rel diff. = 1.2054
[2025-05-31 12:11:29.234] [jointLog] [info] iteration = 100 | max rel diff. = 0.051434
[2025-05-31 12:11:29.241] [jointLog] [info] iteration = 122 | max rel diff. = 0.00603909
[2025-05-31 12:11:29.241] [jointLog] [info] Finished optimizer
[2025-05-31 12:11:29.241] [jointLog] [info] writing output 

[2025-05-31 12:11:29.258] [jointLog] [warning] NOTE: Read Lib [[ /home/ar9416e/mosquito_test/trimmed_reads_polyA/UJ-3092-48-3B_trimmed_cut_R1_paired.fastq.gz, /home/ar9416e/mosquito_test/trimmed_reads_polyA/UJ-3092-48-3B_trimmed_cut_R2_paired.fastq.gz]] :

Detected a *potential* strand bias > 1% in an unstranded protocol check the file: /home/ar9416e/mosquito_test/alignments/UJ-3092-48-3B_trimmed_cut_R2_paired.fastq.gz_quant/lib_format_counts.json for details

[2025-05-31 12:11:29.177] [fileLog] [info] 
At end of round 0
==================
Observed 52239 total fragments (52239 in most recent round)

```
**UJ-3092-Unr-1B:**
```

[2025-05-31 12:11:23.317] [jointLog] [info] setting maxHashResizeThreads to 6
[2025-05-31 12:11:23.317] [jointLog] [info] Fragment incompatibility prior below threshold.  Incompatible fragments will be ignored.
[2025-05-31 12:11:23.317] [jointLog] [info] Usage of --validateMappings implies use of minScoreFraction. Since not explicitly specified, it is being set to 0.65
[2025-05-31 12:11:23.317] [jointLog] [info] Setting consensusSlack to selective-alignment default of 0.35.
[2025-05-31 12:11:23.317] [jointLog] [info] parsing read library format
[2025-05-31 12:11:23.317] [jointLog] [info] There is 1 library.
[2025-05-31 12:11:23.319] [jointLog] [info] Loading pufferfish index
[2025-05-31 12:11:23.319] [jointLog] [info] Loading dense pufferfish index.
[2025-05-31 12:11:28.778] [jointLog] [info] done
[2025-05-31 12:11:28.804] [jointLog] [info] Index contained 15,859 targets
[2025-05-31 12:11:28.806] [jointLog] [info] Number of decoys : 153
[2025-05-31 12:11:28.806] [jointLog] [info] First decoy index : 15,706 
[2025-05-31 12:11:29.484] [jointLog] [info] Automatically detected most likely library type as ISR

[2025-05-31 12:11:35.127] [jointLog] [info] Thread saw mini-batch with a maximum of 0.18% zero probability fragments
[2025-05-31 12:11:35.152] [jointLog] [info] Thread saw mini-batch with a maximum of 0.18% zero probability fragments
[2025-05-31 12:11:35.178] [jointLog] [info] Thread saw mini-batch with a maximum of 0.16% zero probability fragments
[2025-05-31 12:11:35.196] [jointLog] [info] Thread saw mini-batch with a maximum of 0.16% zero probability fragments
[2025-05-31 12:11:35.310] [jointLog] [info] Thread saw mini-batch with a maximum of 0.16% zero probability fragments
[2025-05-31 12:11:35.336] [jointLog] [info] Thread saw mini-batch with a maximum of 0.14% zero probability fragments
[2025-05-31 12:11:35.355] [jointLog] [info] Computed 689 rich equivalence classes for further processing
[2025-05-31 12:11:35.355] [jointLog] [info] Counted 119,590 total reads in the equivalence classes 
[2025-05-31 12:11:35.360] [jointLog] [warning] 0.141393% of fragments were shorter than the k used to build the index.
If this fraction is too large, consider re-building the index with a smaller k.
The minimum read size found was 0.


[2025-05-31 12:11:35.360] [jointLog] [info] Number of mappings discarded because of alignment score : 1,914,695
[2025-05-31 12:11:35.360] [jointLog] [info] Number of fragments entirely discarded because of alignment score : 47,750
[2025-05-31 12:11:35.360] [jointLog] [info] Number of fragments discarded because they are best-mapped to decoys : 6,424
[2025-05-31 12:11:35.360] [jointLog] [info] Number of fragments discarded because they have only dovetail (discordant) mappings to valid targets : 24,057
[2025-05-31 12:11:35.360] [jointLog] [warning] Only 119590 fragments were mapped, but the number of burn-in fragments was set to 5000000.
The effective lengths have been computed using the observed mappings.

[2025-05-31 12:11:35.360] [jointLog] [info] Mapping rate = 15.0974%

[2025-05-31 12:11:35.360] [jointLog] [info] finished quantifyLibrary()
[2025-05-31 12:11:35.363] [jointLog] [info] Starting optimizer
[2025-05-31 12:11:35.383] [jointLog] [info] Marked 0 weighted equivalence classes as degenerate
[2025-05-31 12:11:35.383] [jointLog] [info] iteration = 0 | max rel diff. = 74.9101
[2025-05-31 12:11:35.421] [jointLog] [info] iteration = 100 | max rel diff. = 0.497503
[2025-05-31 12:11:35.458] [jointLog] [info] iteration = 200 | max rel diff. = 0.0703318
[2025-05-31 12:11:35.466] [jointLog] [info] iteration = 222 | max rel diff. = 0.00895351
[2025-05-31 12:11:35.467] [jointLog] [info] Finished optimizer
[2025-05-31 12:11:35.467] [jointLog] [info] writing output 

[2025-05-31 12:11:35.356] [fileLog] [info] 
At end of round 0
==================
Observed 792121 total fragments (792121 in most recent round)

```
Trimming logs showed a very high percentage of reads too short after trimming, indicating degraded input or poor library quality.

Since these samples had extremely low mapping percentages, they were removed from further analysis to ensure reliable downstream results. This was justified as DESeq2 models gene expression dispersion across replicates. Removing samples with poor mapping rates might slightly shift dispersion estimates, but DESeq2 adjusts for this using empirical Bayes shrinkage, so I decided to remove these two reps, compromising on statistical power.

To assess whether low-mapping-quality samples distorted overall gene expression patterns, I performed a differential expression analysis and dimensionality reduction with and without those samples.

Step 1: Filtering Poor-Quality Samples
I identified and removed samples with abnormally low mapping rates. These samples could introduce bias or noise in downstream analyses, particularly in PCA and differential expression.

Step 2: Modeling Gene Expression
I used DESeq2 to fit a full interaction model of the form:
design=∼temperature∗tissue

Step 3: PCA for Expression Pattern Visualization
Using the variance-stabilized data from DESeq2 (RNA-seq gene count data is not normally distributed, thus it must be transformed for PCA).

<img width="941" alt="Screenshot 2025-06-01 at 12 41 06" src="https://github.com/user-attachments/assets/71e2d78b-7ea5-4b6e-8fd5-d29b65771322" />

<img width="941" alt="Screenshot 2025-06-01 at 12 52 53" src="https://github.com/user-attachments/assets/5d0f2004-2ca9-464c-91e3-98b4f8fb33aa" />

**Global PCA results** Tissue is the primary driver of variance (PC1, 78%), PC2 (9%) explains much less but has still meaningful variance. There’s a clear separation along PC1 between the two tissue types, tissue type dominates the expression variance. Temperature influences PC2, body samples show high temperatures such as 48°C higher on PC2. Head samples also show some vertical spreading, but it's less obvious. The Unresponsive variable shows heterogeneity, as it's scattered along the PC. 

<img width="941" alt="Screenshot 2025-06-01 at 13 07 41" src="https://github.com/user-attachments/assets/b547621e-b0a9-4dab-ab71-02da9a38f215" />

** Tissue-specific PCA results: Bodies ** The first principal component (PC1) captures the majority of the variance (72%), indicating a primary factor driving gene expression differences. The temperature values appear to be spread across PC1, suggesting temperature is a major contributor to gene expression variation. Unresponsive and 48C seem to disproportionately contribute to variance, and, specifically, Unresponsive seems to be heterogeneous. 


<img width="941" alt="Screenshot 2025-06-01 at 13 26 28" src="https://github.com/user-attachments/assets/edaa4337-7075-41e6-9163-6fcc4e934824" />


** Tissue-specific PCA results: Heads ** PC1 likely captures a dominant trend in temperature-driven gene expression differences. Samples with temperatures at 36°C and 48°C seem to be positioned differently compared to lower temperatures (25°C and 30°C), although 40°C seems to cluster closer to those at 25°C and 30°C, indicating a shift in gene expression patterns. The Unresponsive (Unr) samples and 48°C appear to be more dispersed, suggesting variability. Lower temperatures (25°C and 30°C) form more compact clusters, suggesting consistent gene expression profiles. Could this be suggesting there's a similar response between those samples at 48°C, 36°C and unresponsive compared to the response observed for those at 30°C, 25°C and 40°C?





