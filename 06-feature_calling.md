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
hicTransform --method pearson -m ZmMC_HiC_2.1.10_sub_1_2.hicup.1mb.corrected.cool --outFileName pearson.cool

hicPlotMatrix -m pearson.cool -o plot_pearson.png --region 2:120000000-132000000 
```

Principal Component Analysis is a dimentionality reduction technique that finds linear combinations of the original matrix variables that a) maximize the variance, b) are independent, and c) are ordered by the percentage of explained variance. This means that applying PCA can be useful to reduce the global correlation patterns to a single vector that captures most of the variability.

The A / B compartments were originally defined in terms of the first component of the matrix. For example, regions with a positive value were assigned to the A compartment, while regions with a negative one to the B compartment.

We can use the hicPCA command to obtain this principal component.


```bash
hicPCA -noe 1 –matrix ZmMC_HiC_2.1.10_sub_1_2.hicup.1mb.corrected.cool  --format bigwig -o pca1.bw
```

Now let's plot the first principal component beside the matrix.


```bash
hicPlotMatrix -m ZmMC_HiC_2.1.10_sub_1_2.hicup.1mb.corrected.cool -o plot_compartments.png --log1p --region 2:120000000-132000000 --bigwig pca1.bedgraph
```

To emphasize compartmentalization, we can also plot the pearson correlation matrix. 


```bash
hicTransform --method pearson -m ZmMC_HiC_2.1.10_sub_1_2.hicup.1mb.corrected.cool --outFileName pearson.cool


hicPlotMatrix -m pearson.cool -o plot_pearson.png --region 2:120000000-132000000 --bigwig pca1.bw
```

The first component sign only tells us which regions are interacting more frequently with which others, so we can assign two compartments this way. However, from the interaction frequency alone we cannot tell wether a region is transcriptionally active or not. 

For this reason, additional information is often used to assign the A and B status to the compartments. RNA seq and histone modification ChIP tracks can be used for this purpose.

Let's make a new plot to add histone modification information (H3K27me3, which should mark the B compartment).


```bash
hicPlotMatrix -m pearson.cool -o plot_pearson.png --region 2:120000000-132000000 --bigwig pca1.bw --bigwig ZmMC_H3K27me3.bw
```


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
hicFindTADs -m ZmMC_HiC_2.1.10_sub_1_2.hicup.1mb.corrected.cool --outPrefix hic_tads --numberOfProcessors 4
```

## Peaks

To detect peaks, we input the corrected matrix. 


```bash
hicDetectLoops -m ZmMC_HiC_2.1.10_sub_1_2.hicup.1mb.corrected.cool -o loops.bedgraph --maxLoopDistance 20000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 --peakInteractionsThreshold 20 --maximumInteractionPercentageThreshold 0.1 --statisticalTest anderson-darling

```

We can plot the matrix with annotated loops.


```bash
hicPlotMatrix -m ZmMC_HiC_2.1.10_sub_1_2.hicup.1mb.corrected.cool -o plot.png --log1p --region 1:18000000-22000000 --loops loops.bedgraph
```

To get an idea of the quality of the callset, besides plotting the annotated matrix we can also produce a plot showing the average of all the called peaks. 


```bash
hicAggregateContacts --matrix  ZmMC_HiC_2.1.10_sub_1_2.hicup.1mb.corrected.cool --BED loops.bedgraph --outFileName ZmMC_HiC_peaks --vMin 0.8 --vMax 2.2 --range 300000:1000000 --numberOfBins 30 --chromosomes 2 --avgType mean --transform obs/exp
```

In the maize Dong et al. paper, they call peaks using HiCCUPs. 

To use this juicer tool we need the .hic matrix. 


```bash
hiccups --cpu --threads 4 -c 2 -r 10000 ZmMC_HiC_2.1.10_sub_1_2.hicup.hic -k KR maize_hiccups_loops
```

Now let's look at the aggregated peak plot. 


```bash
apa -r 10000 -c 2 ZmMC_HiC_2.1.10_sub_1_2.hicup.hic maize_hiccups_loops maize_hiccups_apa 
```

* Challenge:
- Visualize the resulting peaks in juicebox.
