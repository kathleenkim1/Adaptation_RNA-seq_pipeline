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

For the purpose of this tutorial, RNAseq1stgenN1 will be used as an example. 

## Checking fastq quality using **FastQC**
On command-line run fastqc on an individual fastq file or a folder containing all files. 
```bash
fastqc filename.fastq.gz
```
You can view FastQC Report from the generated .html file.

## Trimming adapters with **cutadapt** 
For each sample, there are two forward read files and two reverse read files from lane 1 (L1) and lane 2 (L2). For RNAseq1stgenN1 and all other samples will have these 4 types of files:   
- LIB058795_TRA00258830_RNAseq1stgenN1_S22_L001_R1_001.fastq.gz 
- LIB058795_TRA00258830_RNAseq1stgenN1_S22_L001_R2_001.fastq.gz
- LIB058795_TRA00258830_RNAseq1stgenN1_S22_L002_R1_001.fastq.gz 
- LIB058795_TRA00258830_RNAseq1stgenN1_S22_L002_R2_001.fastq.gz

On command-line run: 
```bash
#Sample cutadapt command
-q 30,30 -a AATGATACGGCGACCACCGAGATCTACACAATAACGTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A GATCGGAAGAGCACACGTCTGAACTCCAGTCACAATCGTTAATCTCGTATGCCGTCTTCTGCTTG -o firstN1one.1.fastq.gz -p firstN1one.2.fastq.gz LIB058795_TRA00258830_RNAseq1stgenN1_S22_L001_R1_001.fastq.gz LIB058795_TRA00258830_RNAseq1stgenN1_S22_L001_R2_001.fastq.gz 

-q 30,30 -a AATGATACGGCGACCACCGAGATCTACACAATAACGTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A GATCGGAAGAGCACACGTCTGAACTCCAGTCACAATCGTTAATCTCGTATGCCGTCTTCTGCTTG -o firstN1two.1.fastq.gz -p firstN1two.2.fastq.gz LIB058795_TRA00258830_RNAseq1stgenN1_S22_L002_R1_001.fastq.gz LIB058795_TRA00258830_RNAseq1stgenN1_S22_L002_R2_001.fastq.gz
```
This code will trim bases with a Phred score lower than 30 and reads shorter than 20 base pairs are discarded. This will also trim the Illumina Truseq I5 (AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT) and I7 adapters (GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[i7]ATCTCGTATGCCGTCTTCTGCTTG). i5 and i7 adapters for each sample can be found in [Hypoxia_Adaptation_Library_Indices.csv](https://github.com/kathleenkim1/Adaptation_RNA-seq_pipeline/blob/main/Hypoxia_Adaptation_Library_Indices.csv) file uploaded to this repository.  
After successfully trimming, you will have 4 output files per sample: 
- firstN1one.1.fastq.gz  
- firstN1one.2.fastq.gz 
- firstN1two.1.fastq.gz 
- firstN1two.2.fastq.gz

These files can know be concatenated by lanes:
```bash
cat firstN1one.1.fastq.gz firstN1two.1.fastq.gz > 1stgenN1read1combined.fastq.gz 
cat firstN1one.2.fastq.gz firstN1two.2.fastq.gz > 1stgenN1read2combined.fastq.gz
```
After concatenating, you will have 2 output files per sample: 
- 1stgenN1read1combined.fastq.gz 
- 1stgenN1read2combined.fastq.gz

  
## Aligning to reference sequence with **bowtie2** and converting to sorted bam file using **samtools**.
```bash
bowtie2 -x GCF_000002985.6_WBcel235_genomic  -1 1stgenN1read1combined.fastq.gz -2 1stgenN1read2combined.fastq.gz -S mapped/1stgenN1aligned.sam --threads 8
```
After alignment, you will have 1 file per sample:
- 1stgenN1aligned.sam

These files must now be converted from sam --> bam --> sorted bam using the following code. 

```bash
samtools view -b -o mapped/1stN1aligned.bam mapped/1stgenN1aligned.sam
samtools sort -@ 8 -o mapped/1stN1aligned.sorted.bam mapped/1stN1aligned.bam
samtools index mapped/1stN1aligned.sorted.bam
rm mapped/1stgenN1aligned.sam mapped/1stN1aligned.bam
```

## Generating counts matrix with **FeatureCounts**.
After generating .sorted.bam files for all 24 samples, FeatureCounts can be run like so:
```bash
featureCounts -p -O -T 8 -a
Caenorhabditis_elegans.WBcel235.109.gtf -o countsmatrix.txt
mapped/1stN1aligned.sorted.bam mapped/1stN2aligned.sorted.bam mapped/1stN3aligned.sorted.bam 
mapped/1stH1aligned.sorted.bam mapped/1stH2aligned.sorted.bam mapped/1stH3aligned.sorted.bam 
mapped/2ndN1aligned.sorted.bam mapped/2ndN2aligned.sorted.bam mapped/2ndN3aligned.sorted.bam 
mapped/2ndH1aligned.sorted.bam mapped/2ndH2aligned.sorted.bam mapped/2ndH3aligned.sorted.bam
mapped/3rdN1aligned.sorted.bam mapped/3rdN2aligned.sorted.bam mapped/3rdN3aligned.sorted.bam 
mapped/3rdH1aligned.sorted.bam mapped/3rdH2aligned.sorted.bam mapped/3rdH3aligned.sorted.bam 
mapped/4thN1aligned.sorted.bam mapped/4thN2aligned.sorted.bam mapped/4thN3aligned.sorted.bam 
mapped/4thH1aligned.sorted.bam mapped/4thH2aligned.sorted.bam mapped/4thH3aligned.sorted.bam
```
The `countsmatrix.txt` can be converted to a .csv for future analysis usage.  

## Creating DEG tables for each generation from counts matrix using R

```{r}
#Install R packages as needed: 
library(DESeq2)
library(AnnotationDbi)
library("org.Ce.eg.db")
library(tidyr)
library(dplyr)
library(ggplot2)

#from featurecounts count matrix, separate by generation
my_count_matrix <- read.delim("countsmatrix.csv", header = TRUE, row.names = 1, sep = ",")

first_exposure_count_matrix <- my_count_matrix[, c(1:6)]
second_exposure_count_matrix  <- my_count_matrix[, c(7:12)]
third_exposure_count_matrix  <- my_count_matrix[, c(13:18)]
fourth_exposure_count_matrix  <- my_count_matrix[, c(19:24)]


#function for running DE 
Get_DEGs = function(counts_matrix, dataframe_name){
counts <- counts_matrix

#filter any genes with counts less than 350
counts <- counts[which(rowSums(counts) > 350),]
condition <- factor(c("Normoxia", "Normoxia", "Normoxia","Hypoxia","Hypoxia","Hypoxia"))
coldata <- data.frame(row.names = colnames(counts), condition)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Hypoxia", "Normoxia"))
significant_genes <- na.omit(res)
significant_genes <- significant_genes[significant_genes$padj < 0.05,]
significant_genes_df <- as.data.frame(significant_genes)
significant_genes_df$symbol <- mapIds(org.Ce.eg.db, keys = rownames(significant_genes_df), keytype = "ENSEMBL", column = "SYMBOL")
print(summary(res))
print(nrow(significant_genes_df))
assign(dataframe_name, significant_genes_df, envir = .GlobalEnv)
}

#run function for each generation

Get_DEGs(first_exposure_count_matrix, "first_exposure_sig_df")
Get_DEGs(second_exposure_count_matrix, "second_exposure_sig_df")
Get_DEGs(third_exposure_count_matrix, "third_exposure_sig_df")
Get_DEGs(fourth_exposure_count_matrix, "fourth_exposure_sig_df")
```

