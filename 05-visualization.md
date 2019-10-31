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


```bash
hicPlotMatrix --matrix ZmMC_HiC_1_1_2.hicup.mcool::/resolutions/100000 --region 2:120000000-132000000 --outFileName full_matrix.png
```

Let's try the same plot, but on a log scale. 


```bash
hicPlotMatrix --matrix ZmMC_HiC_1_1_2.hicup.mcool::/resolutions/100000 --region 2:120000000-132000000 --log1p --outFileName full_matrix_log.png
```

Now open both plots and compare.


```bash
open full_matrix_log.png
open full_matrix.png
```

Another approach is to plot the upper triangle of the diagonal. This is particularly useful to plot the matrix alongside with other data tracks. 

The hicPlotTADs function requires a track configuration file like the following one:


```bash
[x-axis]
where = top

[hic matrix]
file = ZmMC_HiC_2.1.10_sub_1_2.hicup.1mb.corrected.cool
title = Hi-C data
transform = log1p
file_type = hic_matrix

```



```bash
hicPlotTADs --tracks hic_track.ini --region 2:120000000-132000000 -out horiz_mat.png
```




