# Adaptation_RNA-seq_pipeline
A step-by-step data processing pipeline for RNA-seq data used for GSE255871 from the Greer Lab. 

## Introduction
This documentation will go through the pipeline used to process raw fastq RNA-sequencing files to a counts matrix and DGE table. This pipeline uses the following tools that require prior installation: 
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [CutAdapt 4.0](https://cutadapt.readthedocs.io/en/stable/installation.html)
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [SamTools](https://github.com/samtools/samtools)
- [FeatureCounts](https://subread.sourceforge.net/featureCounts.html)
- [WBcel235](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002985.6/)
    - Download the WBcel235 genome in GTF format

## Checking fastq quality using **FastQC**

```r
fastqc filename.fastq.gz
```
You can view FastQC Report from the generated .html file.

## Trimming adapters with **cutadapt** 
For each sample, there are two forward read files and two reverse read files. This code will trim bases with a Phred score lower than 30 and reads shorter than 20 base pairs are discarded. This will also trim the Illumina Truseq I5 (AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT) and I7 adapters (GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[i7]ATCTCGTATGCCGTCTTCTGCTTG)  

```r
#Sample cutadapt command for


```r
-q 30,30 -a AATGATACGGCGACCACCGAGATCTACACAATAACGTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A GATCGGAAGAGCACACGTCTGAACTCCAGTCACAATCGTTAATCTCGTATGCCGTCTTCTGCTTG -o firstN1one.1.fastq.gz -p firstN1one.2.fastq.gz LIB058795_TRA00258830_RNAseq1stgenN1_S22_L001_R1_001.fastq.gz LIB058795_TRA00258830_RNAseq1stgenN1_S22_L001_R2_001.fastq.gz

-q 30,30 -a AATGATACGGCGACCACCGAGATCTACACAATAACGTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A GATCGGAAGAGCACACGTCTGAACTCCAGTCACAATCGTTAATCTCGTATGCCGTCTTCTGCTTG -o firstN1two.1.fastq.gz -p firstN1two.2.fastq.gz LIB058795_TRA00258830_RNAseq1stgenN1_S22_L002_R1_001.fastq.gz LIB058795_TRA00258830_RNAseq1stgenN1_S22_L002_R2_001.fastq.gz


