---
title: "Feature annotation: Compartments, TADs and peaks"
output: 
  html_document:
    keep_md: true
---



## Learning objectives  
- Identify and visualize A/B compartments
- Identify and visualize TADs
- Identify and visualize interaction peaks

## Compartments

One of the coarsest levels of higher order chromatin organization is the spatial segregation of transcriptionally active and inactive genomic regions. In the earliest HiC maps, these regions were identified as a checker board pattern in the off diagonal. Two global compartments were proposed: A compartment of increased interactions among active regions, and B compartment of inactive regions.

The checker board pattern implies that active regions have a similar interaction profile, so their correlation should be higher than with inactive regions.

We can inspect the pearson correlation matrix to visualize this concept.


```bash
hicTransform --method pearson -m ZmMC_500k_corrected.cool --outFileName ZmMC_500k_pearson.cool

hicPlotMatrix -m ZmMC_500k_pearson.cool -o ZmMC_500k_pearson.png
```

Principal Component Analysis is a dimentionality reduction technique that finds linear combinations of the original matrix variables that a) maximize the variance, b) are independent, and c) are ordered by the percentage of explained variance. This means that applying PCA can be useful to reduce the global correlation patterns to a single vector that captures most of the variability.

The A / B compartments were originally defined in terms of the first component of the matrix. For example, regions with a positive value were assigned to the A compartment, while regions with a negative one to the B compartment.

We can use the hicPCA command to obtain this principal component.


```bash
hicPCA -noe 1 --matrix ZmMC_500k_corrected.cool  --format bigwig -o ZmMC_500k_pca1.bw
```

Now let's plot the first principal component beside the matrix.


```bash
hicPlotMatrix -m ZmMC_500k_corrected.cool -o ZmMC_500k_compartments.png --log1p --bigwig ZmMC_500k_pca1.bw --perChromosome
```

To emphasize compartmentalization, we can also plot the pearson correlation matrix. 


```bash

hicPlotMatrix -m ZmMC_500k_pearson.cool -o ZmMC_500k_compartments_pearson.png --bigwig ZmMC_500k_pca1.bw --perChromosome
```

The first component sign only tells us which regions are interacting more frequently with which others, so we can assign two compartments this way. However, from the interaction frequency alone we cannot tell wether a region is transcriptionally active or not. 

For this reason, additional information is often used to assign the A and B status to the compartments. RNA seq and histone modification ChIP tracks can be used for this purpose.

Let's make a new plot to add histone modification information (H3K27me3, which should mark the B compartment).


```bash
wget http://www.epigenome.cuhk.edu.hk/Hi-C_Data/bigwig/ZmMC_H3K4me3.bw
```



```bash
hicPlotMatrix --log1p -m ZmMC_500k_corrected.cool -o ZmMC_500k_histonemod.png --perChromosome --bigwig ZmMC_H3K4me3.bw
```

We can use the histone modification track to flip the sign if needed. In this case, it seems that the A compartment has the negative sign, so let's call the compartments again. 



## TADs

Topologically Associating Domains are defined as regions of increased self interaction in HiC maps. Several molecular mechanisms have been proposed for how they can arise, and many computational tools have been developed to identify them. 

To call TADs using HiCExplorer we can use the hicFindTADs tool. For each bin along the diagonal, hicFindTADs computes a TAD separation score. This score can be interpreted as the amount of contacts that each bin has with other upstream and downstream bins. 

In particular, hicFindTADs works in the following way:
- The corrected matrix is transformed to a Z score matrix, based on the average contact frequency at the same genomic distance.
- The average Z score within a window of upstream and downstream bins of length w is calculated for each bin. This results in a "diamond" submatrix. 
- The previous step is repeated for multiple values of w, and an average per bin is computed. This is the TAD separation score.
- The local minima of the resulting score track are putative TAD boundaries. Note that the TAD separation score can be thought of as the "boundary strength".
- To double check that these are likely boundaries, for each putative boundary bin, the distribution of Z scores within its diamond submatrix, is compared to the upstream and downstream diamond submatrices. 


```bash
hicFindTADs -m ZmMC_50k_corrected.cool --outPrefix ZmMC_50k_tads --correctForMultipleTesting fdr
```

Now we can use hicPlotTADs to visualize them.

First create a .ini file: 

```bash
[x-axis]
fontsize=16
where=top


[hic matrix]
file = ZmMC_50k_corrected.cool
title = Hi-C data
depth = 1000000
transform = log1p
file_type = hic_matrix

[spacer]

[tads]
file = ZmMC_50k_tads_domains.bed
file_type = domains
border color = white
color = gray
overlay previous = share-y
line width = 2
show data range = no

[spacer]

[tad score]
file = ZmMC_50k_tads_tad_score.bm
title = "TAD separation score"
height = 4
file_type = bedgraph_matrix
```

Now let's visualize the TADs.

```bash
hicPlotTADs --tracks hic_tads.ini -o ZmMC_50k_tads.png --region 2:122000000-126000000
```


Challenge: 
  - Call TADs for the ZmEn matrix, and make a plot with both matrices and tads.
  - Try other parameters and / or resolutions.


## Peaks

To detect peaks, we input the corrected matrix. 


```bash
hicDetectLoops -m ZmMC_50k_corrected.cool -o ZmMC_50k_loops.bedgraph --maxLoopDistance 20000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 --peakInteractionsThreshold 20 --maximumInteractionPercentageThreshold 0.1 --statisticalTest anderson-darling

```

We can plot the matrix with annotated loops.


```bash
hicPlotMatrix -m ZmMC_50k_corrected.cool -o ZmMC_50k_loops.png --log1p --region 2:110000000-114000000 --loops ZmMC_50k_loops.bedgraph
```

Challenge: 
  - Look at other loop regions 


In the maize Dong et al. paper, they call peaks using HiCCUPs. 

To use this juicer tool we need the .hic matrix. 

We generated a .hic matrix in a previous practical. Let's make a link to our current working directory.


```bash
ln -s /usr/local/data/hic_workshop_data/prac4_data/ZmEn_HiC_2.hic .
```



```bash
java -jar /usr/local/src/juicer/juicer_tools_1.13.02.jar hiccups --cpu --threads 2 -r 10000 ZmEn_HiC_2.hic -k KR maize_hiccups_loops
```

Now let's look at the aggregated peak plot. 


```bash
java -jar /usr/local/src/juicer/juicer_tools_1.13.02.jar apa -r 10000 ZmEn_HiC_2.hic maize_hiccups_loops maize_hiccups_apa 
```

* Challenge:
- Visualize the resulting peaks in juicebox.
