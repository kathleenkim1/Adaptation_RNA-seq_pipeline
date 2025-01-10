# Adaptation_RNA-seq_pipeline
A step-by-step data processing pipeline for RNA-seq data used for GSE255871 from the Greer Lab. 

## Installation
This pipeline will use utadapt, Bowtie, FeatureCounts 


## Checking fastq quality using **FastQC**


## Trimming adapters with **cutadapt** 
For each sample, there are two forward read files and two reverse read files. This code will trim bases with a Phred score lower than 30 and reads shorter than 20 base pairs are discarded.  


-q 30,30 -a AATGATACGGCGACCACCGAGATCTACAC**XXXXXXXX**ACACTCTTTCCCTACACGACGCTCTTCCGATCT -A GATCGGAAGAGCACACGTCTGAACTCCAGTCAC**XXXXXXXX**ATCTCGTATGCCGTCTTCTGCTTG -o **Read1_output_filename**.fastq.gz -p **Read2_output_filname**.fastq.gz** Read1_filename**.fastq.gz** Read2_filename**.fastq.gz



```r
#Sample cutadapt command 
-q 30,30 -a AATGATACGGCGACCACCGAGATCTACACAATAACGTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A GATCGGAAGAGCACACGTCTGAACTCCAGTCACAATCGTTAATCTCGTATGCCGTCTTCTGCTTG -o firstN1one.1.fastq.gz -p firstN1one.2.fastq.gz LIB058795_TRA00258830_RNAseq1stgenN1_S22_L001_R1_001.fastq.gz LIB058795_TRA00258830_RNAseq1stgenN1_S22_L001_R2_001.fastq.gz

-q 30,30 -a AATGATACGGCGACCACCGAGATCTACACAATAACGTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A GATCGGAAGAGCACACGTCTGAACTCCAGTCACAATCGTTAATCTCGTATGCCGTCTTCTGCTTG -o firstN1two.1.fastq.gz -p firstN1two.2.fastq.gz LIB058795_TRA00258830_RNAseq1stgenN1_S22_L002_R1_001.fastq.gz LIB058795_TRA00258830_RNAseq1stgenN1_S22_L002_R2_001.fastq.gz


