---
title: "Calling genome-wide differential interactions"
output:
  html_document:
    keep_md: true
---

## Learning objectives
We will use diffHic to:  
- Count the read pairs into bins and filter uninteresting bin pairs  
- Normalize our libraries to account for different biases  
- Estimate biological variability of the replicates  
- Obtain our list of differential interactions (DI)  


## 1. Getting HiCUP-mapped data into diffhic
  
Copy the sorted bam files and the digested genome in your host:
```bash
scp lbc@10.10.30.15:/Users/lbc/prac7_data.tar.gz .
tar xvzf prac7_data.tar.gz #uncompress
```
  
If you don´t have terminal you first need to enter the docker container:
```bash
docker run -i -t lizfernandez/hic-langebio:your_tag /bin/bash
scp lbc@10.10.30.15:/Users/lbc/prac7_data.tar.gz .
tar xvzf prac7_data.tar.gz#uncompress
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
# install biocmanager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

#bioconductor libraries
BiocManager::install(Packages)

#statmod
install.packages("statmod")
```




```r
lapply(Packages, library, character.only=T)
library(statmod)
```

Import HiCUP digest file into R and generate digest genome object

```r
digest <- read.csv("Digest_maize_mini2_DpnII_None_21-28-22_04-11-2019.txt", header=T, sep="\t", skip=1) #change name of the file

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
bin.size <- 50000
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
## dim: 28909 4 
## metadata(2): param width
## assays(1): counts
## rownames: NULL
## rowData names(0):
## colnames: NULL
## colData names(1): totals
## type: ReverseStrictGInteractions
## regions: 240
```
*zm_data* contains 515,645 interactions (bin pairs) in all 4 libraries.


```r
head(anchors(zm_data))
```

```
## $first
## GRanges object with 28909 ranges and 1 metadata column:
##                        seqnames            ranges strand |    nfrags
##                           <Rle>         <IRanges>  <Rle> | <integer>
##       [1] 2:120000001-132000000           1-50246      * |       163
##       [2] 2:120000001-132000000       50247-99992      * |       169
##       [3] 2:120000001-132000000       50247-99992      * |       169
##       [4] 2:120000001-132000000      99993-150558      * |       186
##       [5] 2:120000001-132000000      99993-150558      * |       186
##       ...                   ...               ...    ... .       ...
##   [28905] 2:120000001-132000000 11950034-12000000      * |       178
##   [28906] 2:120000001-132000000 11950034-12000000      * |       178
##   [28907] 2:120000001-132000000 11950034-12000000      * |       178
##   [28908] 2:120000001-132000000 11950034-12000000      * |       178
##   [28909] 2:120000001-132000000 11950034-12000000      * |       178
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
## 
## $second
## GRanges object with 28909 ranges and 1 metadata column:
##                        seqnames            ranges strand |    nfrags
##                           <Rle>         <IRanges>  <Rle> | <integer>
##       [1] 2:120000001-132000000           1-50246      * |       163
##       [2] 2:120000001-132000000           1-50246      * |       163
##       [3] 2:120000001-132000000       50247-99992      * |       169
##       [4] 2:120000001-132000000           1-50246      * |       163
##       [5] 2:120000001-132000000       50247-99992      * |       169
##       ...                   ...               ...    ... .       ...
##   [28905] 2:120000001-132000000 11750020-11799897      * |       195
##   [28906] 2:120000001-132000000 11799898-11850047      * |       175
##   [28907] 2:120000001-132000000 11850048-11900292      * |       191
##   [28908] 2:120000001-132000000 11900293-11950033      * |       180
##   [28909] 2:120000001-132000000 11950034-12000000      * |       178
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
## [1] 28909     4
```

```r
head(assay(zm_data))
```

```
##      [,1] [,2] [,3] [,4]
## [1,]  238  309  383  511
## [2,]  149  125  179  187
## [3,]  219  388  290  356
## [4,]   58   39   61   55
## [5,]  117  119  122  117
## [6,]  145  298  224  259
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

![](07_diffhic_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

We filter bin pairs with very low absolute counts

```r
count.keep <- ave.ab >= aveLogCPM(2, lib.size=mean(zm_data$totals))
summary(count.keep)
```

```
##    Mode   FALSE    TRUE 
## logical     551   28358
```

```r
zm_data_or <- zm_data # backup the non-filtered data
zm_data <- zm_data[count.keep,] # apply the filter
zm_data
```

```
## class: InteractionSet 
## dim: 28358 4 
## metadata(2): param width
## assays(1): counts
## rownames: NULL
## rowData names(0):
## colnames: NULL
## colData names(1): totals
## type: ReverseStrictGInteractions
## regions: 240
```

Plot again to check changes

```r
ave.ab <- aveLogCPM(asDGEList(zm_data))
hist(ave.ab, xlab="Average abundance", col="blueviolet", main="Zea mays data after filtering")
```

![](07_diffhic_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

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

![](07_diffhic_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

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
##              [,1]      [,2]        [,3]          [,4]
## [1,] -0.232805996 0.2620612 -0.05675096  0.0274957257
## [2,] -0.008384201 0.1884264 -0.06019441 -0.1198478186
## [3,] -0.190152329 0.2487886 -0.05786232 -0.0007739079
## [4,]  0.278478329 0.1172279 -0.08193442 -0.3137718484
## [5,]  0.079734145 0.1595868 -0.06153263 -0.1777883368
## [6,] -0.109757560 0.2230147 -0.05949599 -0.0537611676
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

![](07_diffhic_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

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
## 1     238     309     383     511
## 2     149     125     179     187
## 3     219     388     290     356
## 4      58      39      61      55
## 5     117     119     122     117
## 28353 more rows ...
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
## [1,] 13.06262 13.55749 13.23868 13.32292
## [2,] 13.28704 13.48385 13.23523 13.17558
## [3,] 13.10527 13.54421 13.23756 13.29465
## [4,] 13.57390 13.41265 13.21349 12.98165
## [5,] 13.37516 13.45501 13.23389 13.11764
## 28353 more rows ...
```
  
Estimate NB dispersion

```r
y <- estimateDisp(y, design)
y$common.dispersion
```

```
## [1] 0.02241369
```

```r
plotBCV(y)
```

![](07_diffhic_files/figure-html/unnamed-chunk-23-1.png)<!-- -->
  
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
##  0.3139 94.8057 94.8057 93.8711 94.8057 94.8057
```

```r
plotQLDisp(fit)
```

![](07_diffhic_files/figure-html/unnamed-chunk-24-1.png)<!-- -->


## 6. TESTING FOR SIGNIFICANT INTERACTIONS

QL F-test

```r
result <- glmQLFTest(fit, coef=2)
topTags(result)
```

```
## Coefficient:  MC 
##            logFC    logCPM        F       PValue          FDR
## 13346 -1.1231542 10.680540 40.05655 7.685422e-09 0.0001291583
## 9865  -1.1290963 10.190604 39.58786 9.109124e-09 0.0001291583
## 19439 -1.0759957 10.492930 36.64549 2.686483e-08 0.0002539443
## 7984  -1.7116681  6.387260 34.67566 5.622975e-08 0.0003986408
## 6493   1.9548499  5.787290 31.51096 2.013291e-07 0.0011418579
## 7997  -0.9821595 10.745390 30.43788 2.873087e-07 0.0013579166
## 16246 -1.9021883  5.909386 28.86006 5.355809e-07 0.0020249001
## 6268   2.0172868  5.567084 28.69803 5.712392e-07 0.0020249001
## 4749  -0.9134636  9.951511 26.16918 1.581572e-06 0.0049833576
## 14507 -1.5570337  6.166139 25.48952 2.088012e-06 0.0057935872
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

![](07_diffhic_files/figure-html/unnamed-chunk-27-1.png)<!-- -->
  
Clustering based on significant bin pairs

```r
clustered.sig <- diClusters(zm_data, result$table, target=0.05, cluster.args=list(tol=1))
length(clustered.sig$interactions)
```

```
## [1] 36
```

```r
head(clustered.sig$interactions)
```

```
## ReverseStrictGInteractions object with 6 interactions and 0 metadata columns:
##                   seqnames1         ranges1                 seqnames2
##                       <Rle>       <IRanges>                     <Rle>
##   [1] 2:120000001-132000000         1-50246 --- 2:120000001-132000000
##   [2] 2:120000001-132000000   249863-300231 --- 2:120000001-132000000
##   [3] 2:120000001-132000000 1549705-1599984 --- 2:120000001-132000000
##   [4] 2:120000001-132000000 2749835-2799978 --- 2:120000001-132000000
##   [5] 2:120000001-132000000 2900413-2949355 --- 2:120000001-132000000
##   [6] 2:120000001-132000000 4800039-4850036 --- 2:120000001-132000000
##               ranges2
##             <IRanges>
##   [1]         1-50246
##   [2]         1-50246
##   [3] 1450000-1500065
##   [4]   199996-249862
##   [5] 2900413-2949355
##   [6] 4800039-4850036
##   -------
##   regions: 41 ranges and 0 metadata columns
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

```r
clustered.sig$FDR
```

```
## [1] 0.02777778
```
  
Create objects for bin pairs identities
Combine Test combines p-values

```r
tabcomdata <- combineTests(clustered.sig$indices[[1]], result$table)
head(tabcomdata)
```

```
## DataFrame with 6 rows and 6 columns
##    nWindows  logFC.up logFC.down               PValue                  FDR
##   <integer> <integer>  <integer>            <numeric>            <numeric>
## 1         1         1          0 6.98300823925361e-05 6.98300823925361e-05
## 2         1         1          0 5.19655380695768e-05 5.85404213192927e-05
## 3         1         1          0 4.55319410551155e-06 1.09276658532277e-05
## 4         1         1          0 5.11079907655807e-05 5.85404213192927e-05
## 5         1         1          0 2.03040674036126e-05 3.17802794143501e-05
## 6         1         0          1 1.58157197900456e-06 6.32628791601824e-06
##     direction
##   <character>
## 1          up
## 2          up
## 3          up
## 4          up
## 5          up
## 6        down
```

```r
#getBestTest finds bin pair with lowest p-values, to find the strongest change in clusters
tabbestdata <- getBestTest(clustered.sig$indices[[1]], result$table)
head(tabbestdata)
```

```
## DataFrame with 6 rows and 6 columns
##        best              logFC           logCPM                F
##   <integer>          <numeric>        <numeric>        <numeric>
## 1         1  0.755133710812747 9.31224585499152 17.2767111629158
## 2        16   1.02581706935246 6.68036637000534 18.0278849330954
## 3       526   1.06287622995692 6.99529897295963 23.6052430478237
## 4      1545   1.67208692036701 5.43050462963205  17.978084998199
## 5      1770  0.813585444719294 9.66647056081663 20.1118927392774
## 6      4749 -0.913463593687237 9.95151138153677 26.1691841649811
##                 PValue                  FDR
##              <numeric>            <numeric>
## 1 6.98300823925361e-05 6.98300823925361e-05
## 2 5.19655380695768e-05 5.85404213192927e-05
## 3 4.55319410551155e-06 1.09276658532277e-05
## 4 5.11079907655807e-05 5.85404213192927e-05
## 5 2.03040674036126e-05 3.17802794143501e-05
## 6 1.58157197900456e-06 6.32628791601824e-06
```
  
Save statistics

```r
tabstat <- data.frame(tabcomdata[,2:6], logFC=tabbestdata$logFC, FDR=clustered.sig$FDR)
result.d <- as.data.frame(clustered.sig$interactions)[,c("seqnames1","start1","end1","seqnames2","start2","end2")]
result.d <- cbind(result.d, tabstat)
o.d <- order(result.d$PValue)
write.table(result.d[o.d,], file="DIclustersData.txt", sep="\t", quote=FALSE, row.names=FALSE)
```





