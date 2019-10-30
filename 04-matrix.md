---
title: "Generating the HiC matrix"
output: 
  html_document:
    keep_md: true
---



## Learning objectives  
- Use juicer tools to obtain a multi resolution hic matrix (for juicebox visualization)
- Use cooler to obtain a multi resolution hic matrix (for HiCexplorer and HiGlass)

## bam to pairs

pairs is a standard format based on plain text proposed by the 4DNucleome consortium. 


```bash
bam2pairs -l -c chr_file.txt ZmMC_HiC_2.1.10_sub_1_2.hicup.bam ZmMC_HiC_2.1.10_sub_1_2.hicup
```

Let's see how the pair format looks like:


```bash
head ZmMC_HiC_2.1.10_sub_1_2.hicup.bsorted.pairs
```

## pairs to .hic

HiC files often are quite large. This means that it is not very efficient to work directly with plain text files.

Different analysis tools support instead binary formats. 

For example,.hic is a binary and compressed format created and supported by juicer and juicer tools from the aiden lab. 

Let's convert our pairs file to a .hic file. 

The juicer_tools pre command bins and corrects the HiC matrix at several resolutions.


```bash
java -jar juicer_tools_1.13.02.jar pre ZmMC_HiC_2.1.10_sub_1_2.hicup.bsorted.pairs ZmMC_HiC_2.1.10_sub_1_2.hicup.hic chr_file.txt
```

The resulting .hic file stores the matrix at various resolutions, as well as different corrections. This file can be directly uploaded to juicebox for visualization, and is ready to use with juicer tools.  

## pairs to .cool

Another common format is .cool. This format is based on the hdf5 (hierachical data format 5) specification, which is also a binary format. There are multiple open source libraries for accessing and writing hdf5 files available for many programming languages. 

.cool files are compatible with HiCExplorer and cooler, common HiC data processing resources. 


```bash
cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 chr_file.txt:10000 ZmMC_HiC_2.1.10_sub_1_2.hicup.bsorted.pairs ZmMC_HiC_2.1.10_sub_1_2.hicup.10kb.cool
```

Unlike the .hic format, .cool files store a single matrix at a single resolution. However, the .mcool can store multiresolution matrices. The .mcool format is compatible with the HiGlass browser. 

Obtain a multi resolution matrix.


```bash
cooler zoomify ZmMC_HiC_2.1.10_sub_1_2.hicup.10kb.cool -r 50000,100000,500000,1000000 -o ZmMC_HiC_2.1.10_sub_1_2.hicup.mcool
```

## Matrix correction

Before matrix correction, bin level filtering is necesary to remove low count bins. 

An approach is the "MAD-max" filter.

To decide filtering values for this filter we can run a diagnostic plot of a histogram of counts per bin, together with the MAD statistic. 


```bash
hicCorrectMatrix diagnostic_plot --matrix ZmMC_HiC_2.1.10_sub_1_2.hicup.10kb.cool -o ZmMC_HiC_2.1.10_sub_1_2_diagnostic_10kb.png
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

- Obtain and correct matrices at 10kb, 50kb, 100kb and 500kb resolutions for all other replicates and conditions.
- Is every resolution appropriate for downstream analysis?

![https://liorpachter.wordpress.com/2013/11/17/imakaev_explained/]


```bash

hicCorrectMatrix diagnostic_plot --matrix ZmMC_HiC_2.1.10_sub_1_2.hicup.mcool::/resolutions/1000000 -o ZmMC_HiC_2.1.10_sub_1_2_diagnostic_1mb.png

hicCorrectMatrix correct --matrix ZmMC_HiC_2.1.10_sub_1_2.hicup.mcool::/resolutions/1000000 --correctionMethod ICE --outFileName ZmMC_HiC_2.1.10_sub_1_2.hicup.1mb.corrected.cool --filterThreshold -1 5

```

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


