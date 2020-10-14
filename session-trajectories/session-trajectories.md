Trajectory inference
====================

Created by: Jules GILET (Institut Curie, France)

Overview
========

Transcriptional trajectories will be inferred from data by Nestorowa,
Hamey et al. ([Blood,
2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5305050/)) The
dataset consists in 1600 hematopoietic stem and progenitor cells from
mouse bone marrow (sequenced using the SMARTseq2 technology). Using flow
cytometry and index sorting, 12 HSPC of different phenotypes (about 10
cells each) have been included in the dataset, and will be used in this
lab as a biological prior for the identification of the root and the
branches in the transcriptional trajectory models.


## Datasets

You can find the datasets within `data.zip` in this directory. Unpack
it, and make sure that it creates a seperate directory named `data` that
includes the following files: 
* nestorowa_corrected_log2_transformed_counts.txt
* nestorowa_corrected_population_annotation.txt
* HTSeq_counts.txt

## Part I - Monocle2/DDRtree


```r
# Inference done with Monocle2/DDRtree available via Bioconductor

library(monocle)
library(biomaRt)
```

### Data loading

The authors provide an expression matrix that has been filtered (highly
expressed genes, high quality cells), scaled and log-normalized. An
annotation table is also provided, with each cell type labelled
according to the immunophenotyping done by flow cytometry.

```r
lognorm <- t(read.table('data/nestorowa_corrected_log2_transformed_counts.txt', sep=" ", header=TRUE))
anno_table <- read.table('data/nestorowa_corrected_population_annotation.txt')
```

To infer a trajectory with Monocle2/DDRtree, using non-normalized
UMI-based counts is highly recommended, as Monocle2 will scale and
normalize the data internally and is expecting data ditributed according
to a negative binomial.

The count matrix has been downloaded and will be used for Monocle2:

```r
counts <- read.table('data/HTSeq_counts.txt', sep="\t", header=TRUE, row.names='ID')

counts[1:5,1:5]
```

```
##                    HSPC_007 HSPC_013 HSPC_019 HSPC_025 HSPC_031
## ENSMUSG00000000001        0        7        1      185        2
## ENSMUSG00000000003        0        0        0        0        0
## ENSMUSG00000000028        4        1        2        4        3
## ENSMUSG00000000031        0        0        0        0        0
## ENSMUSG00000000037        0        0        0        0        0
```

```r
dim(counts)
```

```
## [1] 46175  1920
```

```r
lognorm[1:5,1:5]
```

```
##                HSPC_001 HSPC_002 HSPC_003 HSPC_004 HSPC_006
## X1110032F04Rik 0.000000 0.000000 0.000000 0.000000 0.000000
## X1110059E24Rik 0.000000 0.000000 2.795189 1.326478 7.348663
## X1300017J02Rik 0.000000 0.000000 0.000000 0.000000 0.000000
## X1600014C10Rik 0.000000 2.238601 0.000000 1.326478 4.946766
## X1700017B05Rik 1.225439 2.238601 1.989360 2.005685 0.000000
```

```r
dim(lognorm)
```

```
## [1] 3991 1645
```

Note that the count matrix is not filtered, and genes are labelled
according to ensembl gene IDs. We will first filter the matrix according
to the authors choices (ie. we keep the cells and genes present in the
lognorm matrix) and we will map the gene official symbols.

We filter the counts to keep only high quality cells:

```r
counts <- counts[ , colnames(lognorm) ]
dim(counts)
```

```
## [1] 46175  1645
```

We create a annotation data frame to label the cell types as defined by
the authors:

```r
pDat <- data.frame(cell=colnames(counts), celltype='undefined', stringsAsFactors=FALSE)
rownames(pDat) <- pDat$cell
pDat[ rownames(anno_table), 2] <- as.character(anno_table$celltype)
head(pDat)
```

```
##              cell  celltype
## HSPC_001 HSPC_001 undefined
## HSPC_002 HSPC_002 undefined
## HSPC_003 HSPC_003 undefined
## HSPC_004 HSPC_004 undefined
## HSPC_006 HSPC_006 undefined
## HSPC_008 HSPC_008 undefined
```

We create a feature annotation data frame that will contain gene
informations and matching symbols and IDs. The genes IDs in the counts
matrix are annotated using the biomaRt Bioconductor package:

```r
mart <- biomaRt::useDataset("mmusculus_gene_ensembl", biomaRt::useMart("ensembl"))
genes_table <- biomaRt::getBM(attributes=c("ensembl_gene_id", "external_gene_name"), values=rownames(counts), mart=mart)
rownames(genes_table) <- genes_table$ensembl_gene_id
head(genes_table)
```

```
##                       ensembl_gene_id external_gene_name
## ENSMUSG00000064372 ENSMUSG00000064372              mt-Tp
## ENSMUSG00000064371 ENSMUSG00000064371              mt-Tt
## ENSMUSG00000064370 ENSMUSG00000064370            mt-Cytb
## ENSMUSG00000064369 ENSMUSG00000064369              mt-Te
## ENSMUSG00000064368 ENSMUSG00000064368             mt-Nd6
## ENSMUSG00000064367 ENSMUSG00000064367             mt-Nd5
```

```r
fDat <- genes_table[ rownames(counts), ]
# to be consistent with Monocle naming conventions
colnames(fDat) <- c('ensembl_gene_id', 'gene_short_name')
head(fDat)
```

```
##                       ensembl_gene_id gene_short_name
## ENSMUSG00000000001 ENSMUSG00000000001           Gnai3
## ENSMUSG00000000003 ENSMUSG00000000003            Pbsn
## ENSMUSG00000000028 ENSMUSG00000000028           Cdc45
## ENSMUSG00000000031 ENSMUSG00000000031             H19
## ENSMUSG00000000037 ENSMUSG00000000037           Scml2
## ENSMUSG00000000049 ENSMUSG00000000049            Apoh
```

We can now use this table to filter the genes in the counts matrix that
are highly expressed according to the quality filters used by the
authors:

```r
fDat <- fDat[ fDat$gene_short_name %in% rownames(lognorm), ]
```

And we finally keep in the counts matrix only these genes:

```r
counts <- counts[ rownames(fDat), ]
dim(counts)
```

```
## [1] 3795 1645
```
```r
dim(fDat)
```

```
## [1] 3795    2
```

```r
dim(pDat)
```

```
## [1] 1645    2
```

We build a cell dataset object in an appropriate format for monocle.
Default method for modeling the expression values is
`VGAM::negbinomial.size()` and is adapted to counts.

```r
cds <- newCellDataSet(as.matrix(counts), phenoData=Biobase::AnnotatedDataFrame(pDat), featureData=Biobase::AnnotatedDataFrame(fDat))
cds
```

```
## CellDataSet (storageMode: environment)
## assayData: 3795 features, 1645 samples 
##   element names: exprs 
## protocolData: none
## phenoData
##   sampleNames: HSPC_001 HSPC_002 ... Prog_852 (1645 total)
##   varLabels: cell celltype Size_Factor
##   varMetadata: labelDescription
## featureData
##   featureNames: ENSMUSG00000000001 ENSMUSG00000000028 ...
##     ENSMUSG00000107235 (3795 total)
##   fvarLabels: ensembl_gene_id gene_short_name
##   fvarMetadata: labelDescription
## experimentData: use 'experimentData(object)'
## Annotation:
```

#### Trajectory inference

The monocle cds object is built and ready for trajectory inference.

```r
dir.create('monocle', showWarnings=FALSE)
saveRDS(cds, 'monocle/cds_hematopoiesis.rds')

# Monocle2 preprocess
# normalization and scaling
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
```

```
## Removing 18 outliers
```

```r
# No need for further filtering, the expression matrix has already be filtered
# We find the genes that are expressed:
cds <- detectGenes(cds, min_expr=0.1)
print(head(fData(cds)))
```

```
##                       ensembl_gene_id gene_short_name num_cells_expressed
## ENSMUSG00000000001 ENSMUSG00000000001           Gnai3                1613
## ENSMUSG00000000028 ENSMUSG00000000028           Cdc45                1438
## ENSMUSG00000000056 ENSMUSG00000000056            Narf                1333
## ENSMUSG00000000058 ENSMUSG00000000058            Cav2                 577
## ENSMUSG00000000078 ENSMUSG00000000078            Klf6                1560
## ENSMUSG00000000127 ENSMUSG00000000127             Fer                 578
```

```r
# We identify genes that are expressed in at least 10 cells
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
length(expressed_genes)
```

```
## [1] 3776
```

```r
print(head(pData(cds)))
```

```
##              cell  celltype Size_Factor num_genes_expressed
## HSPC_001 HSPC_001 undefined   3.0716832                2027
## HSPC_002 HSPC_002 undefined   0.3138974                1773
## HSPC_003 HSPC_003 undefined   1.3223730                2187
## HSPC_004 HSPC_004 undefined   0.6412046                2080
## HSPC_006 HSPC_006 undefined   0.6478616                2164
## HSPC_008 HSPC_008 undefined   0.7365603                2015
```

Identification of the ordering genes by differential testing (likelihood
ratio test) ie. genes that are presumed to be important in the
differentiation process captured in the sample. We used the cell types
identified by the authors to define the ordering genes by DE testing.
(Alternatively, a classical approach consist in clustering the cells,
then identify markers genes per clusters.)

```r
diff_test_res <- differentialGeneTest(cds[ expressed_genes, ], fullModelFormulaStr="~ celltype")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
length(ordering_genes)
```

```
## [1] 679
```

We mark the genes that will be used for the ordering :

```r
cds <- setOrderingFilter(cds, ordering_genes)
```

We use the DDRTree algorithm to infer a trajectory with potential
branching points.

```r
cds <- reduceDimension(cds, max_components = 2, method='DDRTree')
cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by="celltype")
```

![](session-trajectories_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```r
# Changing the cell color 
    cell_colors <-  c('lightblue','blue','red','black','orange','yellow','turquoise','lightgrey')
plot_cell_trajectory(cds, color_by="celltype") + scale_color_manual(values=cell_colors)
```

![](session-trajectories_files/figure-html/unnamed-chunk-2-2.png)<!-- -->

The most immature HSCs in the sample express E-Slam. We will define the
root of this model according to this subset of cells:

```r
table(pData(cds)$State, pData(cds)$celltype)[,"ESLAM"]
```

```
##  1  2  3 
## 10  0  0
```

State 1 defines the root in the model as it contains all 10 of the
E-Slam-expressing cells. Note that Monocle might return a different
state number containg these cells. Simply pass the correct state number
to the `orderCells` function:

```r
cds <- orderCells(cds, root_state = 1)
```

The pseudotime is now defined by the distance to the root:

```r
plot_cell_trajectory(cds, color_by = "Pseudotime")
```

#### Differential expression testing per branch

We look at the genes that are differentially expressed according to the
pseudotime model this time:

```r
diff_test_res <- differentialGeneTest(cds[ ordering_genes, ], fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(cds[ sig_gene_names[1:50], ], num_clusters = 3, cores=4, show_rownames=TRUE)
```

![](session-trajectories_files/figure-html/heatmap-1.png)<!-- -->

Differential expression per branch is done with a specific test:
Branched expression analysis modeling (BEAM). The test compares two
models with a likelihood ratio test for branch-dependent expression. The
full model is the product of smooth Pseudotime and the Branch a cell is
assigned to. The reduced model just includes Pseudotime. We look for
genes involved in the erythroid pathway

```r
BEAM_res <- BEAM(cds, branch_point = 1, cores = 4)
```

```r
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
```

```
##                    gene_short_name pval qval
## ENSMUSG00000004552            Ctse    0    0
## ENSMUSG00000004655            Aqp1    0    0
## ENSMUSG00000009350             Mpo    0    0
## ENSMUSG00000016494            Cd34    0    0
## ENSMUSG00000018819            Lsp1    0    0
## ENSMUSG00000020125           Elane    0    0
```

```r
plot_genes_branched_heatmap(cds[row.names(BEAM_res)[1:50]], branch_point = 1, num_clusters = 3, cores=4, use_gene_short_name=TRUE, show_rownames=TRUE)
```

![](session-trajectories_files/figure-html/unnamed-chunk-2-3.png)<!-- -->

There is a clear separation between genes that are involved in the
erythroid differentiation (eg. Gata1) on the left (cell fate1) with
genes involved in the leukocyte differentiation (eg. Sell, Ccl9).

```r
plot_genes_branched_pseudotime(cds[row.names(BEAM_res)[1:5]], branch_point = 1, color_by = "celltype", ncol = 1)  + scale_color_manual(values=cell_colors)
```

![](session-trajectories_files/figure-html/unnamed-chunk-2-4.png)<!-- -->

## Part II - Diffusion map

```r
# Analysis and inference done with the destiny package available via Bioconductor

# Trajectory inference by diffusion map an diffusion pseudotime

library(destiny)
library(ggplot2)
```

#### Data loading

We now will directly use the filtered, scaled, log-normalised expression
matrix provided by the authors of the article.

```r
lognorm <- t(read.table('data/nestorowa_corrected_log2_transformed_counts.txt', sep=" ", header=TRUE))
lognorm[1:5,1:5]
```

```
##                HSPC_001 HSPC_002 HSPC_003 HSPC_004 HSPC_006
## X1110032F04Rik 0,000000 0,000000 0,000000 0,000000 0,000000
## X1110059E24Rik 0,000000 0,000000 2,795189 1,326478 7,348663
## X1300017J02Rik 0,000000 0,000000 0,000000 0,000000 0,000000
## X1600014C10Rik 0,000000 2,238601 0,000000 1,326478 4,946766
## X1700017B05Rik 1,225439 2,238601 1,989360 2,005685 0,000000
```

We load the annotation of cell types that has been defined using flow
cytometry and index sorting. The cell subsets (final differentiation
stages) will be used to validate the trajectory model.

```r
anno_table <- read.table('data/nestorowa_corrected_population_annotation.txt')
pDat <- data.frame(cell=colnames(lognorm), celltype='undefined', stringsAsFactors=FALSE)
rownames(pDat) <- pDat$cell
pDat[ rownames(anno_table), 2] <- as.character(anno_table$celltype)
```

We build an expression set object for an easier integration with
destiny:

```r
eset <- Biobase::ExpressionSet(lognorm, phenoData=Biobase::AnnotatedDataFrame(pDat))
eset
```

```
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 3991 features, 1645 samples 
##   element names: exprs 
## protocolData: none
## phenoData
##   sampleNames: HSPC_001 HSPC_002 ... Prog_852 (1645 total)
##   varLabels: cell celltype
##   varMetadata: labelDescription
## featureData: none
## experimentData: use 'experimentData(object)'
## Annotation:
```

The expression set is ready for inference with destiny:

```r
dir.create('destiny', showWarnings=FALSE)
saveRDS(eset, 'destiny/eset_hematopoiesis.rds')
```

```r
# The process takes less than 60 seconds
dmap <- DiffusionMap(eset)

# We look at the global model
plot.DiffusionMap(dmap)
```

![](session-trajectories_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
plot.DiffusionMap(dmap, dims=c(1,2))
```

![](session-trajectories_files/figure-html/unnamed-chunk-3-2.png)<!-- -->

```r
plot.DiffusionMap(dmap, dims=c(2,3))
```

![](session-trajectories_files/figure-html/unnamed-chunk-3-3.png)<!-- -->

```r
plot.DiffusionMap(dmap, dims=c(1,3))
```

![](session-trajectories_files/figure-html/unnamed-chunk-3-4.png)<!-- -->

Components 1-2 describe well the branching process between erythroid
(red) and myeloid/lymphoid (white) lineages.

We use ggplot2 to have a better rendering and project the cell labels as
defined by flow cytometry experiment and index sorting.

```r
qplot(DC1, DC2, data=dmap, colour=celltype) + scale_color_manual(values=c('lightblue','brown','red','black','orange','yellow','blue','lightgrey')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
```

![](session-trajectories_files/figure-html/unnamed-chunk-3-5.png)<!-- -->

#### Pseudotime inference by diffusion 

The transcriptional distance between cells is estimated by random walk
along a neighborhood graph. The resulting “transcriptional” transition
probability between cells is used to infer a pseudo-time scale of the
differentiation process.

We first define a root cell (origin) for the model. We find the index of
a ESLAM positive cells:

```r
which(anno_table$celltype=="ESLAM")
```

```
##  [1] 19 20 21 22 23 24 25 26 27 28
```

We use this cell as a starting point

```r
dpt <- DPT(dmap, tips=19)
plot(dpt)
```

![](session-trajectories_files/figure-html/unnamed-chunk-3-6.png)<!-- -->

We can project the level of expression of known marker genes on the
trajectory model. Procr / Endothelial protein C is a marker of HSC
subsets:

```r
plot(dpt, col_by='Procr', pal=viridis::magma)
```

![](session-trajectories_files/figure-html/unnamed-chunk-3-7.png)<!-- -->

Gata1 is a key TF of the erythroid lineage

```r
plot(dpt, col_by='Gata1', pal=viridis::magma)
```

![](session-trajectories_files/figure-html/unnamed-chunk-3-8.png)<!-- -->

Cathepsin G is a marker of neutrophils

```r
plot(dpt, col_by='Ctsg', pal=viridis::magma)
```

![](session-trajectories_files/figure-html/unnamed-chunk-3-9.png)<!-- -->

```r
sessionInfo()
```

```
## R version 3.6.1 (2019-07-05)
## Platform: x86_64-conda_cos6-linux-gnu (64-bit)
## Running under: KDE neon User Edition 5.20
## 
## Matrix products: default
## BLAS/LAPACK: /home/mochar/anaconda3/envs/poep/lib/R/lib/libRblas.so
## 
## locale:
## [1] nl_NL.UTF-8
## 
## attached base packages:
##  [1] splines   stats4    parallel  stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] destiny_3.0.1       biomaRt_2.45.8      monocle_2.14.0     
##  [4] DDRTree_0.1.5       irlba_2.3.3         VGAM_1.1-3         
##  [7] ggplot2_3.3.2       Biobase_2.46.0      BiocGenerics_0.32.0
## [10] Matrix_1.2-17      
## 
## loaded via a namespace (and not attached):
##   [1] Rtsne_0.15                  colorspace_1.4-1           
##   [3] RcppEigen_0.3.3.7.0         class_7.3-17               
##   [5] rio_0.5.16                  XVector_0.26.0             
##   [7] RcppHNSW_0.3.0              GenomicRanges_1.38.0       
##   [9] proxy_0.4-24                hexbin_1.27.2              
##  [11] ggrepel_0.8.2               bit64_0.9-7                
##  [13] AnnotationDbi_1.48.0        RSpectra_0.16-0            
##  [15] ranger_0.12.1               xml2_1.2.0                 
##  [17] codetools_0.2-16            docopt_0.7.1               
##  [19] robustbase_0.93-6           knitr_1.22                 
##  [21] cluster_2.0.8               dbplyr_1.4.2               
##  [23] pheatmap_1.0.12             compiler_3.6.1             
##  [25] httr_1.4.0                  ggplot.multistats_1.0.0    
##  [27] assertthat_0.2.1            limma_3.42.2               
##  [29] htmltools_0.3.6             prettyunits_1.0.2          
##  [31] tools_3.6.1                 igraph_1.2.4.1             
##  [33] gtable_0.3.0                glue_1.3.1                 
##  [35] GenomeInfoDbData_1.2.2      RANN_2.6.1                 
##  [37] reshape2_1.4.3              dplyr_0.8.0.1              
##  [39] ggthemes_4.2.0              rappdirs_0.3.1             
##  [41] Rcpp_1.0.5                  carData_3.0-4              
##  [43] slam_0.1-47                 cellranger_1.1.0           
##  [45] vctrs_0.3.4                 lmtest_0.9-38              
##  [47] laeken_0.5.1                xfun_0.6                   
##  [49] stringr_1.4.0               openxlsx_4.2.2             
##  [51] XML_3.98-1.19               DEoptimR_1.0-8             
##  [53] zoo_1.8-6                   zlibbioc_1.32.0            
##  [55] MASS_7.3-51.3               scales_1.0.0               
##  [57] VIM_6.0.0                   pcaMethods_1.78.0          
##  [59] hms_0.4.2                   SummarizedExperiment_1.16.1
##  [61] RColorBrewer_1.1-2          SingleCellExperiment_1.8.0 
##  [63] yaml_2.2.0                  curl_3.3                   
##  [65] memoise_1.1.0               gridExtra_2.3              
##  [67] fastICA_1.2-2               stringi_1.4.3              
##  [69] RSQLite_2.1.1               S4Vectors_0.24.0           
##  [71] knn.covertree_1.0           e1071_1.7-4                
##  [73] TTR_0.24.2                  boot_1.3-25                
##  [75] zip_2.1.1                   densityClust_0.3           
##  [77] BiocParallel_1.20.1         GenomeInfoDb_1.22.1        
##  [79] rlang_0.4.7                 pkgconfig_2.0.2            
##  [81] matrixStats_0.57.0          bitops_1.0-6               
##  [83] qlcMatrix_0.9.7             evaluate_0.13              
##  [85] lattice_0.20-38             purrr_0.3.2                
##  [87] labeling_0.3                bit_1.1-14                 
##  [89] tidyselect_0.2.5            plyr_1.8.4                 
##  [91] magrittr_1.5                R6_2.4.0                   
##  [93] IRanges_2.20.0              combinat_0.0-8             
##  [95] DelayedArray_0.12.3         DBI_1.0.0                  
##  [97] pillar_1.3.1                haven_2.3.1                
##  [99] foreign_0.8-71              withr_2.1.2                
## [101] xts_0.12.1                  scatterplot3d_0.3-41       
## [103] abind_1.4-5                 RCurl_1.98-1.2             
## [105] sp_1.4-2                    nnet_7.3-14                
## [107] tibble_2.1.1                crayon_1.3.4               
## [109] car_3.0-10                  BiocFileCache_1.10.0       
## [111] rmarkdown_1.12              viridis_0.5.1              
## [113] progress_1.2.0              readxl_1.3.1               
## [115] grid_3.6.1                  data.table_1.12.2          
## [117] blob_1.1.1                  FNN_1.1.3                  
## [119] forcats_0.4.0               vcd_1.4-8                  
## [121] HSMMSingleCell_1.6.0        sparsesvd_0.2              
## [123] digest_0.6.18               tidyr_0.8.3                
## [125] openssl_1.3                 munsell_0.5.0              
## [127] viridisLite_0.3.0           smoother_1.1               
## [129] askpass_1.0
```