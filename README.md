# Adaptation_RNA-seq_pipeline
A step-by-step data processing pipeline for RNA-seq data used for GSE255871 from the Greer Lab. 

## Installation

## Checking fastq quality using **FastQC**


## Trimming adapters with **cutadapt** 
For each sample, there are two forward read files and two reverse read files. This code will trim bases with a Phred score lower than 30 and reads shorter than 20 base pairs are discarded.  


```r
#Sample cutadapt command 
py -m cutadapt -q 30,30 -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o firstN1one.1.fastq.gz -p firstN1one.2.fastq.gz LIB058795_TRA00258830_RNAseq1stgenN1_S22_L001_R1_001.fastq.gz LIB058795_TRA00258830_RNAseq1stgenN1_S22_L001_R2_001.fastq.gz![image](https://github.com/user-attachments/assets/5ade6813-b106-4cb4-8107-4797fbf5e732)
"py -m cutadapt -q 30,30 -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o firstN1two.1.fastq.gz -p firstN1two.2.fastq.gz LIB058795_TRA00258830_RNAseq1stgenN1_S22_L002_R1_001.fastq.gz LIB058795_TRA00258830_RNAseq1stgenN1_S22_L002_R2_001.fastq.gz
"![image](https://github.com/user-attachments/assets/c2d993ce-8cf0-4563-b5fa-ee72b2973e45)


