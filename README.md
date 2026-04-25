# S-locus genotyping pipeline for apricot

Hedhly A. et al. (2025)

This repository contains scripts and command-line examples used for S-locus genotyping of *Prunus armeniaca* accessions from whole-genome sequencing data.

The workflow maps filtered reads to a synthetic S-locus reference, calculates per-gene coverage, calls S-RNase and SFB alleles using Breadth of coverage, and includes scripts used to visualize the Breadth threshold and test reference-removal robustness.

## Reference version

In this study, `V31` is the label used for the synthetic S-locus reference sequence. If applying the workflow to another reference, update reference names, Bowtie2 index names, BED files, file names, and output labels accordingly.

## Repository structure

```text
Rscripts/
├── 01_automatic_S_genotype_calling.R
├── 02_breadth_threshold_distribution.R
└── 03_reference_removal_robustness_test.R
```

## Pipeline overview

1. K-mer based filtering of reads using `kmerRefFilter`
2. Read trimming using Trimmomatic
3. Mapping to the synthetic S-locus reference using Bowtie2
4. Filtering and sorting alignments with Samtools
5. Coverage calculation using Bedtools
6. Genotype calling based on Breadth of coverage
7. Breadth-threshold visualization
8. Reference-removal robustness analysis

## Requirements

- Python, for `kmerRefFilter`
- Trimmomatic 0.39
- Bowtie2
- Samtools
- Bedtools
- R
- R packages: `ggplot2`, `gridExtra`, `svglite`

## Example command-line workflow

The commands below show the structure of the analysis. Replace paths, file names, thread counts, and reference labels with those matching your system.

### 1. K-mer filtering

```bash
python kmerRefFilter.py \
  -P 12 \
  -r S_allele_database.fasta \
  -1 accession_id_1FWD.fastq.gz \
  -2 accession_id_2REV.fastq.gz \
  -ugzip \
  -mem 0 \
  -o filtered/
```

### 2. Read trimming

```bash
java -jar trimmomatic-0.39.jar PE \
  -threads 12 \
  -phred33 \
  filtered/accession_id_1FWD_filtered.fastq \
  filtered/accession_id_2REV_filtered.fastq \
  trimmed/accession_id_1FWD_filtered_trimmo_paired.fastq \
  trimmed/accession_id_1FWD_filtered_trimmo_unpaired.fastq \
  trimmed/accession_id_2REV_filtered_trimmo_paired.fastq \
  trimmed/accession_id_2REV_filtered_trimmo_unpaired.fastq \
  ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
  LEADING:3 \
  TRAILING:3 \
  SLIDINGWINDOW:4:15 \
  MINLEN:25
```

### 3. Mapping and BAM filtering

```bash
bowtie2 \
  --end-to-end \
  --very-sensitive \
  --maxins 550 \
  --no-unal \
  -x index/_Parm_GenmV31 \
  -1 trimmed/accession_id_1FWD_filtered_trimmo_paired.fastq \
  -2 trimmed/accession_id_2REV_filtered_trimmo_paired.fastq \
  --threads 12 \
| samtools view -b -F 0x4 -q 20 -e '[NM]<=2' - \
| samtools sort -o mapped/accession_id_F0x4q20eNM3.sort.bam
```

```bash
samtools index mapped/accession_id_F0x4q20eNM3.sort.bam
```

### 4. Coverage calculation

```bash
bedtools coverage \
  -a ref/_Parm_GenmV31_genes.bed \
  -b mapped/accession_id_F0x4q20eNM3.sort.bam \
  > coverage/accession_id_F0x4q20eNM3_BedT_coverage_hist.tsv
```

The R scripts expect coverage files ending in:

```text
_F0x4q20eNM3_BedT_coverage_hist.tsv
```

## R scripts

### `Rscripts/01_automatic_S_genotype_calling.R`

Calls S-locus genotypes from Bedtools coverage output.

For each accession, the script:

```text
reads Bedtools coverage output
extracts gene names and Breadth values
keeps genes with Breadth >= 0.50
collapses called genes per accession
writes Genotypes_table.csv
```

Output:

```text
Genotypes_table.csv
```

### `Rscripts/02_breadth_threshold_distribution.R`

Plots the distribution of all Breadth values and marks the 0.50 threshold used for allele calling.

The y-axis is manually split by stacking two y-axis ranges. This helps visualize both the large number of near-zero Breadth values and the lower-frequency high-Breadth values. The split range can be adjusted in the script with:

```r
y_break <- c(1000, 12000)
```

Outputs:

```text
Breadth_distribution_plot_data.csv
Breadth_distribution_y_axis_break.png
Breadth_distribution_y_axis_break.svg
```

### `Rscripts/03_reference_removal_robustness_test.R`

Generates the reference-removal robustness figure comparing mapping to the full reference versus modified references lacking target allele(s).

In this study, the tested cases were:

```text
S13 removal
S6 removal
S4 + S23 removal
```

The script checks whether reads from accessions carrying removed target alleles are reassigned above the 0.50 Breadth threshold to other alleles when the target allele is absent from the reference.

Expected modified-reference coverage file names follow:

```text
<Accession>_case_<RemapCaseName>_F0x4q20eNM3_BedT_coverage_hist.tsv
```

Outputs:

```text
remapping_breadth_plot_data.csv
remapping_all_nonzero_breadth.csv
remapping_breadth_robustness.png
remapping_breadth_robustness.pdf
remapping_breadth_robustness.svg
```

## Notes

- A minimum Breadth threshold of 0.50 was used for allele calling.
- The Breadth-threshold distribution script uses all Breadth values, including zero and low-coverage values. It does not filter for called genotypes.
- Mapping filters used in the example are: mapped reads only, mapping quality >= 20, and edit distance `NM <= 2`.
- The string `NM3` in file names reflects the historical naming convention of the workflow; the filtering command retains reads with `NM <= 2`, meaning reads with three or more edit-distance differences are excluded.
- Thread counts should be adjusted according to the user system.
- The scripts use generic placeholder paths and should be edited before running.
