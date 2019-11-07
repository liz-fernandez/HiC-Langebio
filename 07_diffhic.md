---
title: "Calling genome-wide differential interactions"
output:
  html_document:
    keep_md: true
---

## Learning objectives
We will use diffHic to:  
- Count the read pairs into bins to measure interaction intensity  
- Filter our data to remove uninteresting bin pairs  
- Normalize our libraries to account for different biases  
- Estimate biological variability of the replicates  
- Obtain our list of differential interactions (DI)  


## 1. Getting HiCUP-mapped data into diffhic
  
Enter the docker container
```bash
docker run -i -t lizfernandez/hic-langebio:tag /bin/bash
#or
docker run 
```
  
Bam files need to be sorted:
```bash
for f in *.bam ; do file="${f%%.*}"; samtools sort -n $f >  ${file}.hicup.sorted.bam ; done
```
  
Copy these bam files as well as de digested genome in your host
```bash
#In your shell outside of docker:
docker cp <containerId>:/path/file /your/path/
```

### Working on RStudio:

Set your working directory to where your files are

```r
#example
setwd("/Users/americaramirezcolmenero/Documents/hic_workshop_langebio")
# or do: session > set working directory > choose directory
```

Install the required packages


```r
Packages <- c("diffHic","csaw", "edgeR", "GenomicRanges")
```

```r
BiocManager::install(Packages)
install.packages("statmod")
```

```r
lapply(Packages, library, character.only=T)
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: stats4
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind,
##     colnames, dirname, do.call, duplicated, eval, evalq, Filter,
##     Find, get, grep, grepl, intersect, is.unsorted, lapply, Map,
##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
##     pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
##     setdiff, sort, table, tapply, union, unique, unsplit, which,
##     which.max, which.min
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```
## Loading required package: IRanges
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: InteractionSet
```

```
## Loading required package: SummarizedExperiment
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Loading required package: DelayedArray
```

```
## Loading required package: matrixStats
```

```
## 
## Attaching package: 'matrixStats'
```

```
## The following objects are masked from 'package:Biobase':
## 
##     anyMissing, rowMedians
```

```
## Loading required package: BiocParallel
```

```
## 
## Attaching package: 'DelayedArray'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges
```

```
## The following objects are masked from 'package:base':
## 
##     aperm, apply, rowsum
```

```
## Loading required package: limma
```

```
## 
## Attaching package: 'limma'
```

```
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
```

```r
library(statmod)
```
  

Import HiCUP digest file into R and generate digest genome object

```r
digest <- read.csv("Digest_maize_mini2_DpnII_None_21-28-22_04-11-2019.txt", header=T, sep="\t", skip=1)
#The object zm.frag has the digested genome in the diffHic required format
zm.frag <- with(digest,GRanges(Chromosome, IRanges(Fragment_Start_Position,Fragment_End_Position)))
#lets look at it
zm.frag
```

```
## GRanges object with 42569 ranges and 0 metadata columns:
##                        seqnames            ranges strand
##                           <Rle>         <IRanges>  <Rle>
##       [1] 2:120000001-132000000              1-73      *
##       [2] 2:120000001-132000000             74-93      *
##       [3] 2:120000001-132000000            94-203      *
##       [4] 2:120000001-132000000           204-713      *
##       [5] 2:120000001-132000000           714-761      *
##       ...                   ...               ...    ...
##   [42565] 2:120000001-132000000 11999122-11999197      *
##   [42566] 2:120000001-132000000 11999198-11999474      *
##   [42567] 2:120000001-132000000 11999475-11999740      *
##   [42568] 2:120000001-132000000 11999741-11999814      *
##   [42569] 2:120000001-132000000 11999815-12000000      *
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

```r
#Generate pair.param, to store the restriction fragments.
#zm.param will hold other parameters later
zm.param <- pairParam(zm.frag)
zm.param
```

```
## Genome contains 42569 restriction fragments across 1 chromosome
## No discard regions are specified
## No limits on chromosomes for read extraction
## No cap on the read pairs per pair of restriction fragments
```
  
  
`preparePairs` creates the h5 files needed for counting data into bins, it matches the mapping location of each read with a restriction fragment:

```r
#Do this for every bam file

# Zea mays endosperm
preparePairs("ZmEn_HiC_1_1_2.hicup.sorted.bam", zm.param, file="En1_Zm.h5")
preparePairs("ZmEn_HiC_2_1_2.hicup.sorted.bam", zm.param, file="En2_Zm.h5")
# Zea mays mesophyll
preparePairs("ZmMC_HiC_1_1_2.hicup.sorted.bam", zm.param, file="MC1_Zm.h5")
preparePairs("ZmMC_HiC_2_1_2.hicup.sorted.bam", zm.param, file="MC2_Zm.h5")
```

Create **input**, an object that has the name of our files

```r
input <- c("En1_Zm.h5","En2_Zm.h5","MC1_Zm.h5","MC2_Zm.h5")
input
```

```
## [1] "En1_Zm.h5" "En2_Zm.h5" "MC1_Zm.h5" "MC2_Zm.h5"
```
  
##2. COUNTING READS INTO BINS

First we will indicate the size (in base pairs) of our bin

```r
bin.size <- 10000
```

`SquareCounts` is the function that counts our read pairs between the pairs of bins, for all of our libraries, using **input**, and will store this information in **zm_data**.


```r
  zm_data <- squareCounts(input, zm.param, width=bin.size, filter=1) 
```

Check the structure of *zm_data*:

```r
zm_data
```

```
## class: InteractionSet 
## dim: 515645 4 
## metadata(2): param width
## assays(1): counts
## rownames: NULL
## rowData names(0):
## colnames: NULL
## colData names(1): totals
## type: ReverseStrictGInteractions
## regions: 1200
```
*zm_data* contains 515,645 interactions (bin pairs) in all 4 libraries.


```r
head(anchors(zm_data))
```

```
## $first
## GRanges object with 515645 ranges and 1 metadata column:
##                         seqnames            ranges strand |    nfrags
##                            <Rle>         <IRanges>  <Rle> | <integer>
##        [1] 2:120000001-132000000           1-10505      * |        43
##        [2] 2:120000001-132000000       10506-21281      * |        28
##        [3] 2:120000001-132000000       10506-21281      * |        28
##        [4] 2:120000001-132000000       21282-30183      * |        34
##        [5] 2:120000001-132000000       21282-30183      * |        34
##        ...                   ...               ...    ... .       ...
##   [515641] 2:120000001-132000000 11989823-12000000      * |        39
##   [515642] 2:120000001-132000000 11989823-12000000      * |        39
##   [515643] 2:120000001-132000000 11989823-12000000      * |        39
##   [515644] 2:120000001-132000000 11989823-12000000      * |        39
##   [515645] 2:120000001-132000000 11989823-12000000      * |        39
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
## 
## $second
## GRanges object with 515645 ranges and 1 metadata column:
##                         seqnames            ranges strand |    nfrags
##                            <Rle>         <IRanges>  <Rle> | <integer>
##        [1] 2:120000001-132000000           1-10505      * |        43
##        [2] 2:120000001-132000000           1-10505      * |        43
##        [3] 2:120000001-132000000       10506-21281      * |        28
##        [4] 2:120000001-132000000           1-10505      * |        43
##        [5] 2:120000001-132000000       10506-21281      * |        28
##        ...                   ...               ...    ... .       ...
##   [515641] 2:120000001-132000000 11950034-11960081      * |        26
##   [515642] 2:120000001-132000000 11960082-11970341      * |        36
##   [515643] 2:120000001-132000000 11970342-11980481      * |        41
##   [515644] 2:120000001-132000000 11980482-11989822      * |        36
##   [515645] 2:120000001-132000000 11989823-12000000      * |        39
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
Each row is an interaction, each column is a library. 
Each bin has _nfrags_ fragments.
**Note:** The boundary of each bin is rounded to the closest restriction site.
  
This object also has a count matrix with the number of read pairs for each interaction, in each library

```r
dim(assay(zm_data))
```

```
## [1] 515645      4
```

```r
head(assay(zm_data))
```

```
##      [,1] [,2] [,3] [,4]
## [1,]   48   75   61   91
## [2,]   29   44   54   48
## [3,]   16   44   43   39
## [4,]   12   12   17   19
## [5,]   16   32   27   33
## [6,]   10   18   22   46
```

```r
zm_data$totals # total read pairs per library
```

```
## [1] 726265 606426 584113 485457
```
  
  
## 3. FILTERING BIN PAIRS

Visualize the distribution of the average abundance

```r
ave.ab <- aveLogCPM(asDGEList(zm_data))
hist(ave.ab, xlab="Average abundance", col="cadetblue3", main="Zea mays data before filtering")
```

![](07_diffhic_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

We filter bin pairs with very low absolute counts

```r
count.keep <- ave.ab >= aveLogCPM(2, lib.size=mean(zm_data$totals))
summary(count.keep)
```

```
##    Mode   FALSE    TRUE 
## logical  462732   52913
```

```r
zm_data_or <- zm_data # backup the non-filtered data
zm_data <- zm_data[count.keep,] #apply the filter
zm_data
```

```
## class: InteractionSet 
## dim: 52913 4 
## metadata(2): param width
## assays(1): counts
## rownames: NULL
## rowData names(0):
## colnames: NULL
## colData names(1): totals
## type: ReverseStrictGInteractions
## regions: 1200
```

Plot again to check changes

```r
ave.ab <- aveLogCPM(asDGEList(zm_data))
hist(ave.ab, xlab="Average abundance", col="blueviolet", main="Zea mays data after filtering")
```

![](07_diffhic_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

###### Other strategies (apply when having genome-wide data):

```r
#Retain only those bin pairs with abundances x-times higher than the median abundance across inter-chromosomal bin pairs.
direct <- filterDirect(zm_data)
keep <- direct$abundances > log2(3) + direct$threshold
data <- data[keep, ]
```

## 4. NORMALIZATION

Generate a MA plot comparing one library of each group

```r
ab <- aveLogCPM(asDGEList(zm_data))
o <- order(ab)
adj.counts <- cpm(asDGEList(zm_data), log=TRUE)
mval <- adj.counts[,3]-adj.counts[,2]
smoothScatter(ab, mval, xlab="A", ylab="M", main="En (1) vs MC (2)")
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")
```

![](07_diffhic_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

Perform non-linear normalization

```r
zm_data <- normOffsets(zm_data, se.out=TRUE)
```

The matrix of offsets has same dimensions as count matrix  and is stored as anelement of assays slot of the InteractionSet object.   
  
Create object for Matrix of offsets

```r
nb.off <- assay(zm_data, "offset")
head(nb.off)
```

```
##            [,1]      [,2]          [,3]        [,4]
## [1,] -0.4913710 0.4415066 -0.1288988384  0.17876321
## [2,] -0.3890223 0.3392252 -0.0821666812  0.13196379
## [3,] -0.3431208 0.2952738 -0.0625416628  0.11038860
## [4,] -0.1286706 0.1357166 -0.0009035591 -0.00614246
## [5,] -0.2798118 0.2380070 -0.0378370881  0.07964195
## [6,] -0.2577983 0.2197260 -0.0303318216  0.06840410
```

Plot after normalization

```r
ab <- aveLogCPM(asDGEList(zm_data))
o <- order(ab)
adj.counts <- cpm(asDGEList(zm_data), log=TRUE)
mval <- adj.counts[,3]-adj.counts[,2]
smoothScatter(ab, mval, xlab="A", ylab="M", main="En (1) vs MC (2) afterNLN")
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")
```

![](07_diffhic_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

## 5. MODELING BIOLOGICAL VARIABILITY

The NB model also consider extra-Poisson variability between biological replicates of the same conditon.  
Variability is modelled by estimating the dispersion parameter of the NB distribution  
  
Specify design matrix to describe the experimental setup

```r
design <- model.matrix(~factor(c("En", "En", "MC", "MC")))
colnames(design) <- c("Intercept", "MC")
design
```

```
##   Intercept MC
## 1         1  0
## 2         1  0
## 3         1  1
## 4         1  1
## attr(,"assign")
## [1] 0 1
## attr(,"contrasts")
## attr(,"contrasts")$`factor(c("En", "En", "MC", "MC"))`
## [1] "contr.treatment"
```

```r
y <- asDGEList(zm_data)
y
```

```
## An object of class "DGEList"
## $counts
##   Sample1 Sample2 Sample3 Sample4
## 1      48      75      61      91
## 2      29      44      54      48
## 3      16      44      43      39
## 4      12      12      17      19
## 5      16      32      27      33
## 52908 more rows ...
## 
## $samples
##         group lib.size norm.factors
## Sample1     1   726265            1
## Sample2     1   606426            1
## Sample3     1   584113            1
## Sample4     1   485457            1
## 
## $offset
##          [,1]     [,2]     [,3]     [,4]
## [1,] 12.80406 13.73693 13.16653 13.47419
## [2,] 12.90640 13.63465 13.21326 13.42739
## [3,] 12.95231 13.59070 13.23288 13.40581
## [4,] 13.16676 13.43114 13.29452 13.28928
## [5,] 13.01561 13.53343 13.25759 13.37507
## 52908 more rows ...
```
  
Estimate NB dispersion

```r
y <- estimateDisp(y, design)
y$common.dispersion
```

```
## [1] 0.03858963
```

```r
plotBCV(y)
```

![](07_diffhic_files/figure-html/unnamed-chunk-22-1.png)<!-- -->
  
Estimate QL dispersion  
Estimation of QL dispersion is performed to model variability of the dispersions.  
**Note:** robust=TRUE protect EB shrinkage against positive outliers (highly variable counts).

```r
library(statmod)
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)
summary(fit$df.prior)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.428  16.779  16.779  16.779  16.779  16.779
```

```r
plotQLDisp(fit)
```

![](07_diffhic_files/figure-html/unnamed-chunk-23-1.png)<!-- -->


## 6. TESTING FOR SIGNIFICANT INTERACTIONS

QL F-test

```r
result <- glmQLFTest(fit, coef=2)
topTags(result)
```

```
## Coefficient:  MC 
##           logFC   logCPM        F       PValue        FDR
## 30051 -2.360129 6.877610 46.52626 1.747771e-06 0.02758004
## 34718 -2.476954 6.375154 46.49172 1.756640e-06 0.02758004
## 46741 -2.708376 7.355397 45.67129 1.982645e-06 0.02758004
## 26579 -2.038484 8.219364 44.96066 2.204602e-06 0.02758004
## 43575 -1.918634 8.309786 41.99890 3.477439e-06 0.02758004
## 30050 -2.285957 6.503117 40.51552 4.407083e-06 0.02758004
## 19999 -2.085247 7.805370 40.13877 4.684964e-06 0.02758004
## 9866  -2.665358 7.312811 39.70436 5.029709e-06 0.02758004
## 29584 -2.654877 7.707343 39.61271 5.105982e-06 0.02758004
## 43640 -2.055887 7.719528 39.48740 5.212338e-06 0.02758004
```

Save significance statistics in rowData of InteractionSet object

```r
rowData(zm_data) <- cbind(rowData(zm_data), result$table)
```
  
Plot to visualize

```r
de <- decideTestsDGE(result, p.value=0.05, adjust.method="BH")
debins <- rownames(result)[as.logical(de)]
plotSmear(result, de.tags=debins)
```

![](07_diffhic_files/figure-html/unnamed-chunk-26-1.png)<!-- -->
  
Clustering based on significant bin pairs

```r
clustered.sig <- diClusters(zm_data, result$table, target=0.05, cluster.args=list(tol=1))
length(clustered.sig$interactions)
```

```
## [1] 8
```

```r
head(clustered.sig$interactions)
```

```
## ReverseStrictGInteractions object with 6 interactions and 0 metadata columns:
##                   seqnames1         ranges1                 seqnames2
##                       <Rle>       <IRanges>                     <Rle>
##   [1] 2:120000001-132000000 2549906-2559945 --- 2:120000001-132000000
##   [2] 2:120000001-132000000 4810045-4820002 --- 2:120000001-132000000
##   [3] 2:120000001-132000000 6259771-6269648 --- 2:120000001-132000000
##   [4] 2:120000001-132000000 6880275-6889819 --- 2:120000001-132000000
##   [5] 2:120000001-132000000 6980126-6989663 --- 2:120000001-132000000
##   [6] 2:120000001-132000000 8020044-8029946 --- 2:120000001-132000000
##               ranges2
##             <IRanges>
##   [1] 2549906-2559945
##   [2] 4810045-4820002
##   [3] 6259771-6269648
##   [4] 6880275-6889819
##   [5] 6970404-6989663
##   [6] 8020044-8029946
##   -------
##   regions: 9 ranges and 0 metadata columns
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

```r
clustered.sig$FDR
```

```
## [1] 0
```
  
Create objects for bin pairs identities
Combine Test combines p-values

```r
tabcomdata <- combineTests(clustered.sig$indices[[1]], result$table)
head(tabcomdata)
```

```
## DataFrame with 6 rows and 6 columns
##    nWindows  logFC.up logFC.down               PValue                 FDR
##   <integer> <integer>  <integer>            <numeric>           <numeric>
## 1         1         0          1  5.0297091356763e-06 5.2123380629601e-06
## 2         1         0          1 4.68496405602028e-06 5.2123380629601e-06
## 3         1         0          1 2.20460166356349e-06 5.2123380629601e-06
## 4         1         0          1 5.10598234433251e-06 5.2123380629601e-06
## 5         2         0          2 3.49554212929873e-06 5.2123380629601e-06
## 6         1         0          1 1.75664035890623e-06 5.2123380629601e-06
##     direction
##   <character>
## 1        down
## 2        down
## 3        down
## 4        down
## 5        down
## 6        down
```

```r
#getBestTest finds bin pair with lowest p-values, to find the strongest change in clusters
tabbestdata <- getBestTest(clustered.sig$indices[[1]], result$table)
head(tabbestdata)
```

```
## DataFrame with 6 rows and 6 columns
##        best             logFC           logCPM                F
##   <integer>         <numeric>        <numeric>        <numeric>
## 1      9866 -2.66535815179669 7.31281117912775 39.7043622176752
## 2     19999 -2.08524680381319 7.80536950485475 40.1387699719517
## 3     26579 -2.03848394291428 8.21936441174658 44.9606582587497
## 4     29584 -2.65487727057066 7.70734277439999 39.6127054827804
## 5     30051 -2.36012880009591 6.87761039510477 46.5262612155578
## 6     34718 -2.47695422935969  6.3751540718002  46.491715911073
##                 PValue                  FDR
##              <numeric>            <numeric>
## 1  5.0297091356763e-06 5.83540839352287e-06
## 2 4.68496405602028e-06 5.83540839352287e-06
## 3 2.20460166356349e-06 5.83540839352287e-06
## 4 5.10598234433251e-06 5.83540839352287e-06
## 5 3.49554212929873e-06 5.83540839352287e-06
## 6 1.75664035890623e-06 5.83540839352287e-06
```
  
Save statistics

```r
tabstat <- data.frame(tabcomdata[,1:4], logFC=tabbestdata$logFC, FDR=clustered.sig$FDR)
result.d <- as.data.frame(clustered.sig$interactions)[,,]
result.d <- cbind(result.d, tabstat)
o.d <- order(result.d$PValue)
write.table(result.d[o.d,], file="DIclustersData.tsv", sep="\t", quote=FALSE, row.names=FALSE)
```





