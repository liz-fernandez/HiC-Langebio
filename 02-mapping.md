---
layout: page
title: Introduction to Next-generation Sequencing
subtitle: Read alignment
minutes: 5
---
> ## Learning Objectives {.objectives}
>
> * Align sequencing data to a reference (genomes or transcriptomes)
> * Understand how to interpret mass sequencing alignment data
> * First approach to the SAM and BAM coordinate formats

We will use the fastq files that we used in the previous practice, as well as a mock reference genome:

~~~ {.bash}
$ cd /usr/local/data
$ mkdir MAP
$ cd MAP
$ wget https://liz-fernandez.github.io/HiC-Langebio/datasets/genome/Sp_genome.fa
$ wget https://liz-fernandez.github.io/HiC-Langebio/datasets/Sp_ds.left.fq.gz
$ wget https://liz-fernandez.github.io/HiC-Langebio/datasets/Sp_ds.right.fq.gz
~~~

## Mapping the filtered reads to the genome

Once we verify that the reads are in the correct format, we will align the reads to the genome using Bowtie2 via TopHat.

You can find the manual in the following[link](https://ccb.jhu.edu/software/tophat/manual.shtml).

First we will generate a bowtie2 index for the genome:

~~~ {.bash}
$ bowtie2-build Sp_genome.fa Sp_genome 
~~~
 
~~~ {.bash}
$ ls *bt2
~~~ 

~~~ {.output}
Sp_genome.1.bt2        
Sp_genome.2.bt2              
Sp_genome.3.bt2       
Sp_genome.4.bt2           
Sp_genome.rev.1.bt2          
Sp_genome.rev.2.bt2
~~~

Usamos tophat2 para mapear las lecturas. Este programa nos permite dividir lecturas
que atraviesan sitios de splicing:

~~~ {.bash}
$ bowtie2 -x Sp_genome -1 Sp_ds.left.fq.gz -2 Sp_ds.right.fq.gz > ds_Sp_genome.sam
~~~ 

~~~ {.output}
101577 reads; of these:
  101577 (100.00%) were paired; of these:
    3927 (3.87%) aligned concordantly 0 times
    97412 (95.90%) aligned concordantly exactly 1 time
    238 (0.23%) aligned concordantly >1 times
    ----
    3927 pairs aligned concordantly 0 times; of these:
      193 (4.91%) aligned discordantly 1 time
    ----
    3734 pairs aligned 0 times concordantly or discordantly; of these:
      7468 mates make up the pairs; of these:
        5332 (71.40%) aligned 0 times
        2125 (28.45%) aligned exactly 1 time
        11 (0.15%) aligned >1 times
97.38% overall alignment rate
~~~

We explore the result, which is a file in SAM format.

~~~ {.bash}
$ head ds_Sp_genome.sam
~~~ 

### The SAM format

Let's see a smaller example of this format.
Suppose we have the following alignment:

~~~ {.output}
Coor	12345678901234 5678901234567890123456789012345
ref	AGCATGTTAGATAA**GATAGCTGTGCTAGTAGGCAGTCAGCGCCAT
+r001/1	      TTAGATAAAGGATA*CTG
+r002	     aaaAGATAA*GGATA
+r003	   gcctaAGCTAA
+r004	                 ATAGCT..............TCAGC
-r003	                        ttagctTAGGC
-r001/2	                                      CAGCGGCAT
~~~

The corresponding SAM format will be the following:

~~~ {.output}
@HD	VN:1.5	SO:coordinate
@SQ	SN:ref	LN:45
r001	99	ref	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*
r002	0	ref	9	30	3S6M1P1I4M	*	0	0	AAAAGATAAGGATA	*
r003	0	ref	9	30	5S6M	*	0	0	GCCTAAGCTAA	*	SA:Z:ref,29,-,6H5M,17,0;
r004	0	ref	16	30	6M14N5M	*	0	0	ATAGCTTCAGC	*
r003	2064	ref	29	17	6H5M	*	0	0	TAGGC	*	SA:Z:ref,9,+,5S6M,30,1;
r001	147	ref	37	30	9M	=	7	-39	CAGCGGCAT	*	NM:i:1
~~~

The SAM format is a plain text format that allows us to save sequencing data
in ASCII format delimited by tabs. 

It is made up of two core sections:

*  The headers
*  The alignment

The section of the **header** starts with the character `@` followed by one of the codes
of two letters that denote the characteristics of the alignments in this file.
Each line is delimited by tabs and, in addition to the lines that begin with
`@CO`, each data field has the format` TAG:VALUE`, where `TAG` is a string
Two characters that define the format and content of `VALUE`.

The header is not indispensable but contains information about the version of the
file as well as if it is ordered or not. Therefore it is advisable to include it.

The **alignment section** contains the following information:

1. **QNAME** Name of the reference, QNAME (SAM) / Name of the read (BAM).
It is used to group alignments that are together, as in the case of alignments
of paired end reads or a read that appears in multiple alignments.
2. **FLAG** Information set describing the alignment. Provides the following information:
	* Are there multiple fragments?
	* Are all the fragments well aligned?
	* Is this fragment aligned?
	* Has the following fragment not been aligned?
	* Is the reference the reverse string?
	* Is the following fragment the reverse chain?
	* Is this the first fragment?
	* Is this the last fragment?
	* Is this a secondary alignment?
	* Did this reading fail quality filters?
	* Is this reading a PCR or optical duplicate?
3. **RNAME** Name of the reference sequence.
3. **POS** Left alignment position (base 1).
3. **MAP** Alignment quality.
3. **CIGAR** CIGAR chain.
3. **RNEXT** Name of the paired end (mate) or the next read.
3. **PNEXT** Position of the pair (mate) or the next read.
3. **TLEN** Length of the alignment.
3. **SEQ** The test sequence of this alignment (in this case the sequence of the read).
3. **QUAL** The read quality.
3. **TAGs** Additional information.

> ## Cadenas CIGAR {.callout}
> The sequence aligned to the reference may have additional bases that are not in
> the reference or the read may be missing bases that are in the reference.
> The CIGAR string is a string that encodes each base and the characteristics of each in
> alignment.
>
> For example, the string CIGAR:
>
> ~~~ {.output}
> CIGAR: 3M1I3M1D5M
> ~~~
>
> indicates that the first 3 bases of the reading aligns with the reference (3M), the next base
> does not exist in the reference (1I), the following 3 bases align with the reference (3M), the
> next base does not exist in the reading (1D), and 5 more bases align with the reference (5M).

As you can see these files contain a lot of information that can be analyzed
using scripts that show alignment statistics. Programs like
[Picard](http://broadinstitute.github.io/picard/) perform
this type of analysis.



















