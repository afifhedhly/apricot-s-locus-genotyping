# S-locus genotyping pipeline (Prunus armeniaca)

This repository contains the pipeline used for genotyping S-RNase and SFB alleles from whole-genome sequencing data, as described in:

Lora et al. (2025) – A comprehensive next-generation sequencing approach for S-locus genotyping in apricot

## Overview

The pipeline consists of the following steps:

1. K-mer–based filtering of reads using kmerRefFilter  
2. Read trimming using Trimmomatic  
3. Mapping to a synthetic S-locus reference sequence using Bowtie2  
4. Filtering and sorting alignments with Samtools  
5. Coverage calculation using Bedtools  
6. Genotype calling based on coverage breadth using a custom R script  

## Requirements

- Python (kmerRefFilter)
- Trimmomatic 0.39
- Bowtie2 (≥2.5)
- Samtools (≥1.22)
- Bedtools (≥2.31)
- R

## Example pipeline

### 1. K-mer filtering
```bash
python kmerRefFilter.py -P 12 -r S-allele_database.fasta \
-1 accession_id_1FWD.fastq.gz -2 accession_id_2REV.fastq.gz \
-ugzip -mem 0 -o outdir/
```

### 2. Trimming
```bash
java -jar trimmomatic-0.39.jar PE -threads 12 -phred33 accession_id_filered_1FWD.fastq accession_id_filered_2REV.fastq \
accession_id_filered_trimm_paired_1FWD accession_id_filered_trimm_unpaired_1FWD \
accession_id_filered_trimm_paired_2REV accession_id_filered_trimm_unpaired_2REV \
ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
```

### 3. Mapping
```bash
bowtie2 --end-to-end --very-sensitive --maxins 550 --no-unal \
-x index/_Parm_GenmV42 \
-1 accession_id.forward.trimmed.fastq \
-2 accession_id.reverse.trimmed.fastq \
--threads 12 | samtools view -b -F 0x4 -q 20 -e '[NM]<=2' - | samtools sort -o output_dir/accession_id.sort.bam
```

### 4. Coverage
```bash
bedtools coverage -a ref_genome.bed -b output_dir/accession_id.sort.bam > coverage.tsv
```

### 5. Genotype calling

See script: `genotype_calling.R`

---

### Notes

- A minimum breadth threshold of 50% was used for allele calling  
- Mapping parameters were optimized to minimize cross-mapping between closely related alleles  
- The pipeline can flag potential novel alleles through incomplete or ambiguous mapping patterns
