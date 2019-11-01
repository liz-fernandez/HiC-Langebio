---
title: "HiC alignment strategies"
output: 
  html_document:
    keep_md: true
---



## Learning objectives  
- Generate the genome index for mapping
- Generate a restriction fragment digested genome file 
- Run HiCUP truncator to truncate the 3' part of reads after the ligation motif
- Run HiCUP 


## Setting up the genome files

In order to run HiCUP to map and filter HiC reads, two files must be generated: the bowtie2 index, and the digested genome file.

First we generate a bowtie2 index for the genome.


```bash
bowtie2-build --threads 1 maize_mini2.fa maize_mini2
```

Then we use HiCUP digester to generate an *in silico* digested genome file. 
- --re1 is the restriction enzyme used in the procol. The restriction site, as well as the cut site (with ^) must be indicated. 
- --genome is the name of the genome for the output file (optional)


```bash
hicup_digester --re1 ^GATC,DpnII --genome maize_mini2 maize_mini2.fa
```

Let's inspect the Digested file.


```bash
head Digest_maize_mini2_DpnII_None_17-05-28_30-10-2019.txt
```

## HiCUP truncater 

The next step is to truncate the sequence downstream of a ligation site of the read.


```bash
hicup_truncater --re1 ^GATC,DpnII ZmEn_HiC_sub_1.fq.gz  ZmEn_HiC_sub_2.fq.gz
```

Let's inspect the truncation results.


```bash
less hicup_truncater_summary_wJYYjDJNrR_18-03-50_26-10-2019.txt
```

We expect a higher percentage of truncated reads with longer reads (~150 nts).
Another factor is the distribution of restriction fragment lengths.

After truncation, the average read length is ~80 which is still a reasonable read length. 

## HiCUP mapper

The next step is to map the read pairs to the reference genome. 

Forward and reverse reads are mapped independently, and then the resulting alignments are paired again to produce a paired end bam.


```bash
hicup_mapper --threads 1 --bowtie2 /data/software/bowtie2-2.3.5.1-linux-x86_64/bowtie2 --index maize_mini2 ZmEn_HiC_sub_1.trunc.fastq ZmEn_HiC_sub_2.trunc.fastq
```

We can inspect how many read pairs were correctly mapped in the hicup_mapper_summary file.


```bash
head hicup_mapper_summary_cKQMtNGDIb_18-27-31_26-10-2019.txt 
```

## HiCUP filter

After mapping, the resulting SAM file is parsed to filter out uninformative read pairs. 


```bash
hicup_filter --digest Digest_maize_mini2_DpnII_None_02-07-29_01-11-2019.txt ZmEn_HiC_sub_1_2.pair.sam --longest 800 --shortest 150
```

## HiCUP deduplicater

The final step of the workflow is to remove read pair duplicates. Duplicates may arise during the PCR protocol or in the sequencing step (optical duplicates). Removing duplicates is done by comparing the start and end coordinates of both reads of a read pair.


```bash
hicup_deduplicator --zip ZmEn_HiC_sub_1_2.filt.sam
```

## Run whole HiCUP pipeline 

An useful feature of HiCUP is that it can be run as a complete pipeline, producing a final html report. 

To do this we setup a configuration file:


```bash
#Example configuration file for the hicup Perl script - edit as required
########################################################################

#Directory to which output files should be written (optional parameter)
#Set to current working directory by default 
Outdir:

#Number of threads to use
Threads: 1

#Suppress progress updates (0: off, 1: on)
Quiet:0

#Retain intermediate pipeline files (0: off, 1: on)
Keep:0

#Compress outputfiles (0: off, 1: on)
Zip:1

#Path to the alignment program (Bowtie or Bowtie2)
#Remember to include the executable Bowtie/Bowtie2 filename.
#Note: ensure you specify the correct aligner i.e. Bowtie when 
#using Bowtie indices, or Bowtie2 when using Bowtie2 indices. 
#In the example below Bowtie2 is specified.
Bowtie2: /data/software/bowtie2-2.3.5.1-linux-x86_64/bowtie2

#Path to the reference genome indices
#Remember to include the basename of the genome indices
Index: maize_mini2

#Path to the genome digest file produced by hicup_digester
Digest: Digest_maize_mini2_DpnII_None_17-05-28_30-10-2019.txt

#FASTQ format (valid formats: 'Sanger', 'Solexa_Illumina_1.0', 'Illumina_1.3' or 'Illumina_1.5')
#If not specified, HiCUP will try to determine the format automatically by analysing
#one of the FASTQ files. All input FASTQ will assumed to be in this format
Format: 

#Maximum di-tag length (optional parameter)
Longest: 800

#Minimum di-tag length (optional parameter)
Shortest: 100

#FASTQ files to be analysed, placing paired files on adjacent lines
ZmEn_HiC_sub_1.fq.gz
ZmEn_HiC_sub_2.fq.gz
```


```bash
hicup --config hicup_config.txt
```

Finally let's inspect the html output.


```bash
open ZmMC_HiC_2.1.10_sub_1_2.HiCUP_summary_report.html
```


