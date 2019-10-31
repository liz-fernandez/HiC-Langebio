---
title: "Generating a HiC matrix"
output: 
  html_document:
    keep_md: true
---



## Learning objectives {.objectives}
- Use juicer tools to obtain a multi resolution hic matrix (for juicebox visualization)
- Use cooler to obtain a multi resolution hic matrix (for HiCexplorer and HiGlass)

## HiC .pairs data format

In this section we will use our filtered hic read pairs in bam format, and convert it to a  data format commonly used by downstream HiC analysis tools. 

First, we will obtain a .pairs format. 

.pairs is a standard format based on plain text proposed by the [4DNucleome consortium](https://commonfund.nih.gov/4dnucleome). 

Here is an example: 


```bash
## pairs format v1.0
#sorted: chr1-chr2-pos1-pos2
#shape: upper triangle
#genome_assembly: hg38
#chromsize: chr1 249250621
#chromsize: chr2 243199373
#chromsize: chr3 198022430

#columns: readID chr1 pos1 chr2 pos2 strand1 strand2
EAS139:136:FC706VJ:2:2104:23462:197393 chr1 10000 chr1 20000 + +
EAS139:136:FC706VJ:2:8762:23765:128766 chr1 50000 chr1 70000 + +
EAS139:136:FC706VJ:2:2342:15343:9863 chr1 60000 chr2 10000 + +
EAS139:136:FC706VJ:2:1286:25:275154 chr1 30000 chr3 40000 + -

```

HiC data is often quite large, which means that it is not very efficient to work directly with plain text files. Tools for compression and indexing are available for the .pairs format. Some downstream analysis tools accept compressed and indexed files, which makes storage and access easier.

You can read more about the specification and supporting tools here:
https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md

We will use the bam2pairs command to obtain a .pairs file.


```bash
bam2pairs -l -c chr_file.txt ZmMC_HiC_1_1_2.hicup.bam ZmMC_HiC_1_1_2.hicup
```

Let's see how the .pairs format looks like:


```bash
head ZmMC_HiC_1_1_2.hicup.bsorted.pairs
```

## Binning and matrix formats 

The next step in the workflow is to aggregate the read level pairs into bins. Once the data is binned, other formats are used to store the matrix data. 

We will bin and store matrices using two common tools: juicer_tools, which stores the resulting matrix in the .hic format, and cooler, which uses the .cool format. 

Both formats are binary and compressed.

### .hic

The .hic is mainly supported by juicer, juicer_tools, and juicebox from the [Aiden lab](https://github.com/aidenlab/juicer). It is a binary format that stores a HiC matrix with multiple bin sizes and corrections in a single file.

Let's convert our .pairs file to a .hic file. 

The juicer_tools pre command bins and corrects the HiC matrix at several resolutions.


```bash
java -Xmx1G -jar juicer_tools_1.13.02.jar pre ZmMC_HiC_1_1_2.hicup.bsorted.pairs ZmMC_HiC_1_1_2.hicup.hic chr_file.txt
```

The resulting .hic file stores the matrix at various resolutions, as well as different corrections. This file can be directly uploaded to juicebox for visualization, and is ready to use with juicer_tools.  

### .cool

Another common format is .cool. This format is based on the hdf5 (hierachical data format 5) specification, which is also a binary format. There are multiple open source libraries for accessing and writing hdf5 files available for many programming languages. 

.cool files are compatible with HiCExplorer and cooler, common HiC data processing resources. 
Let's convert our .pairs file to a .cool file binned at 10Kb resolution.


```bash
cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 chr_file.txt:10000 ZmMC_HiC_1_1_2.hicup.bsorted.pairs ZmMC_HiC_1_1_2.hicup.10kb.cool
```

Unlike the .hic format, .cool files store a single matrix at a single resolution. However, the .mcool format can store multiresolution matrices. The HiGlass browser supports .mcool matrices.

Let's obtain a multi resolution matrix with cooler zoomify command.


```bash
cooler zoomify ZmMC_HiC_1_1_2.hicup.10kb.cool -r 10000,100000,500000 -o ZmMC_HiC_1_1_2.hicup.mcool
```


## Matrix correction

Matrix correction is necesary to remove biases like gc content or mappability.

In order to correct a matrix, it is assumed that if no biases were affecting the experiment, each bin should have equal "visibility" of contacts.

This translates to an intuitive solution to matrix correction: transforming the matrix in such a way that the total number of contacts of every row and every column is the same. 

Such a procedure is called "matrix balancing", and many algorithms have been described to achieve this for applications outside HiC data analysis.

For HiC data, the most common ones are called Knight-Ruiz (KR), and Iterative Correction (ICE). 

https://liorpachter.wordpress.com/2013/11/17/imakaev_explained/


Before matrix correction, bin level filtering is necesary to remove low count bins. 

An approach is the "MAD-max" filter.

To decide filtering values for this filter we can run a diagnostic plot of a histogram of counts per bin, together with the MAD statistic. 


```bash
hicCorrectMatrix diagnostic_plot --matrix ZmMC_HiC_1_1_2.hicup.mcool::/resolutions/10000 -o ZmMC_HiC_1_1_2_diagnostic_10kb.png
```

Let's explore the diagnostic plot and choose the cutoff values.


```bash
open ZmMC_HiC_2.1.10_sub_1_2_diagnostic_10kb.png
```

After deciding on minimum and maximum values, we can proceed to correct the matrix. 


```bash
hicCorrectMatrix correct --matrix ZmMC_HiC_2.1.10_sub_1_2.hicup.10kb.cool --correctionMethod ICE --outFileName ZmMC_HiC_2.1.10_sub_1_2.hicup.10kb.corrected.cool --filterThreshold -1 5

```

## Challenge:

- Obtain and correct matrices at 10kb and 100kb bin sizes for all other replicates and conditions.



## Matrix QC

When working with an experimental design that has multiple conditions and replicates, it is useful to assess how similar the replicates are, and how different the conditions are.

One approach for global comparison is to plot the linear distance versus the contact frequency. The interaction frequency should decay with distance. However, the shape of the decay curve is often characteristic of a single experimental condition. Replicates should have similar decay curves, while conditions should have different ones.


```bash
hicPlotDistVsCounts --matrices ZmMC_HiC_2.1.10_sub_1_2.hicup.1mb.corrected.cool -o plot_vs_counts.png 
```



## Sum matrices

A common practice in HiC data is to sum biological replicate matrices in order to increase sequencing depth, and thus matrix resolution. This can be done after checking that the biological replicates are indeed similar. It is advised to also conduct downstream analyses separately on each replicate to assess differences at those levels. 


```bash

```


