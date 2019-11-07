---
title: "Matrix Visualization"
output: 
  html_document:
    keep_md: true
---



## Learning objectives  
- Use juicebox to visualize a HiC matrix
- Use HiCExplorer to generate images of a matrix slice 

## Juicebox interactive visualization

A common visualization task for HiC data is to interactively explore the matrix for quality control purposes, to generate hypothesis or to confirm expected interaction patterns. 

For the first part of the practical, we are using the juicebox browser interactively. 

- Load a HiC matrix
- Explore multiple resolutions and corrections
- Lock resolutions and zoom
- Add compartment information.
- Explore a distance normalized matrix 
- Load a second HiC matrix as a control
- Add 1D tracks
- Add 2D tracks
- Export an interesting view

## HiCexplorer for command line visualization

Another visualization approach is to include a command in a script that automatically generates figures. This is useful for generating reproducible workflows.

We will use the HiCexplorer command line tools for this purpose.

First, let's plot a square matrix. This might be useful to explore large regions. 
Let's take a look at the whole chromosome.


```bash
hicPlotMatrix --perChromosome --matrix ZmMC_500k_corrected.cool --outFileName ZmMC_500k_corrected.png

```

Now let's try the same plot, but on a log scale. 


```bash
hicPlotMatrix --perChromosome --log1p --matrix ZmMC_500k_corrected.cool --outFileName ZmMC_500k_log.png 

```

Now open both plots and compare.


```bash
open full_matrix_log.png
open full_matrix.png
```

We can also plot a smaller region.


```bash
hicPlotMatrix --log1p --region 2:125000000-130000000 --matrix ZmMC_50k_corrected.cool --outFileName ZmMC_50k_region_log.png
```

We can compare differences between our conditions by plotting a difference matrix. 

```bash
hicCompareMatrices --operation log2ratio --matrices ZmMC_500k_corrected.cool ZmEn_500k_corrected.cool --outFileName ZmMC_ZmEn_500k_log2.cool
```


```bash
hicPlotMatrix --perChromosome --matrix ZmMC_ZmEn_500k_log2.cool --outFileName ZmMC_ZmEn_500k_log2.png
```


Another approach is to plot the upper triangle of the diagonal. This is particularly useful to plot the matrix alongside with other data tracks. 

The hicPlotTADs function requires a track configuration file like the following one:


```bash
[hic matrix]
file = ZmMC_10k_corrected.cool
title = Hi-C data
transform = log1p
file_type = hic_matrix
depth = 300000

```

Copy and paste the configuration file above and save it as "hic_track.ini". You can use the nano editor for this.


```bash
hicPlotTADs --tracks hic_track.ini --region 2:127000000-128000000 -out ZmMC_10k_horiz_mat.png
```




