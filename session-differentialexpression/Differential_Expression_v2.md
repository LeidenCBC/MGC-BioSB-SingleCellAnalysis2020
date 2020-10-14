Differential Expression
================

Created by: Ahmed Mahfouz

# Overview

In this tutorial we will explore different methods to perform
differential expression analysis on scRNA-seq data.

Load required packages:

``` r
suppressMessages(require(Seurat))
```

    ## Warning: package 'Seurat' was built under R version 3.6.3

``` r
suppressMessages(require(scran))
```

    ## Warning: package 'S4Vectors' was built under R version 3.6.3

    ## Warning: package 'GenomeInfoDb' was built under R version 3.6.3

    ## Warning: package 'DelayedArray' was built under R version 3.6.3

    ## Warning: package 'matrixStats' was built under R version 3.6.3

``` r
suppressMessages(require(BiocFileCache))
```

    ## Warning: package 'dbplyr' was built under R version 3.6.3

``` r
suppressMessages(require(DropletUtils))
suppressMessages(require(scater))
```

    ## Warning: package 'ggplot2' was built under R version 3.6.3

``` r
suppressMessages(require(EnsDb.Hsapiens.v86))
suppressMessages(require(SingleCellExperiment))
suppressMessages(require(pheatmap))
```

## Data preprocessing

``` r
pbmc <- readRDS(file = "pbmc3k_final.rds")
```

## Default DE tests in Seurat

Seurat’s differential expression analysis can be performed using the
FindMarkers function. Differential expression is performed between
groups of cells. The function will automatically retrieve the cluster
identities from the seurat object using the `Idents()` function. To test
for differential expression between two specific groups of cells,
specify the ident.1 and ident.2 parameters.

As a default, Seurat performs differential expression based on the
non-parameteric Wilcoxon rank sum test.

List options for groups to perform differential expression
    on

``` r
levels(pbmc)
```

    ## [1] "Naive CD4 T"  "Memory CD4 T" "CD14+ Mono"   "B"            "CD8 T"       
    ## [6] "FCGR3A+ Mono" "NK"           "DC"           "Mk"

``` r
# Find differentially expressed features between CD14+ and FCGR3A+ Monocytes
monocyte.de.markers <- FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono")
# view results
head(monocyte.de.markers)
```

    ##                p_val avg_logFC pct.1 pct.2    p_val_adj
    ## FCGR3A 1.193617e-101 -2.617707 0.131 0.975 1.636926e-97
    ## LYZ     8.134552e-75  1.812078 1.000 0.988 1.115572e-70
    ## RHOC    4.479768e-68 -1.611576 0.162 0.864 6.143554e-64
    ## S100A8  7.471811e-65  2.610696 0.975 0.500 1.024684e-60
    ## S100A9  1.318422e-64  2.286734 0.996 0.870 1.808084e-60
    ## IFITM2  4.821669e-64 -1.445772 0.677 1.000 6.612437e-60

The results data frame has the following columns : p\_val : p\_val
(unadjusted) avg\_logFC : log fold-chage of the average expression
between the two groups. Positive values indicate that the feature is
more highly expressed in the first group. pct.1 : The percentage of
cells where the feature is detected in the first group pct.2 : The
percentage of cells where the feature is detected in the second group
p\_val\_adj : Adjusted p-value, based on bonferroni correction using all
features in the dataset.

If the ident.2 parameter is omitted or set to NULL, FindMarkers will
test for differentially expressed features between the group specified
by ident.1 and all other cells.

### Session info

``` r
sessionInfo()
```

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 18363)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.1252 
    ## [2] LC_CTYPE=English_United States.1252   
    ## [3] LC_MONETARY=English_United States.1252
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.1252    
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] pheatmap_1.0.12             EnsDb.Hsapiens.v86_2.99.0  
    ##  [3] ensembldb_2.10.2            AnnotationFilter_1.10.0    
    ##  [5] GenomicFeatures_1.38.2      AnnotationDbi_1.48.0       
    ##  [7] scater_1.14.6               ggplot2_3.3.2              
    ##  [9] DropletUtils_1.6.1          BiocFileCache_1.10.2       
    ## [11] dbplyr_1.4.4                scran_1.14.6               
    ## [13] SingleCellExperiment_1.8.0  SummarizedExperiment_1.16.1
    ## [15] DelayedArray_0.12.3         BiocParallel_1.20.1        
    ## [17] matrixStats_0.57.0          Biobase_2.46.0             
    ## [19] GenomicRanges_1.38.0        GenomeInfoDb_1.22.1        
    ## [21] IRanges_2.20.2              S4Vectors_0.24.4           
    ## [23] BiocGenerics_0.32.0         Seurat_3.2.2               
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] plyr_1.8.6               igraph_1.2.6             lazyeval_0.2.2          
    ##   [4] splines_3.6.2            listenv_0.8.0            digest_0.6.25           
    ##   [7] htmltools_0.5.0          viridis_0.5.1            memoise_1.1.0           
    ##  [10] magrittr_1.5             tensor_1.5               cluster_2.1.0           
    ##  [13] ROCR_1.0-11              limma_3.42.2             Biostrings_2.54.0       
    ##  [16] globals_0.13.1           R.utils_2.10.1           askpass_1.1             
    ##  [19] prettyunits_1.1.1        colorspace_1.4-1         blob_1.2.1              
    ##  [22] rappdirs_0.3.1           ggrepel_0.8.2            xfun_0.18               
    ##  [25] dplyr_1.0.2              crayon_1.3.4             RCurl_1.98-1.2          
    ##  [28] jsonlite_1.7.1           spatstat_1.64-1          spatstat.data_1.4-3     
    ##  [31] survival_3.2-7           zoo_1.8-8                glue_1.4.2              
    ##  [34] polyclip_1.10-0          gtable_0.3.0             zlibbioc_1.32.0         
    ##  [37] XVector_0.26.0           leiden_0.3.3             BiocSingular_1.2.2      
    ##  [40] Rhdf5lib_1.8.0           future.apply_1.6.0       HDF5Array_1.14.4        
    ##  [43] abind_1.4-5              scales_1.1.1             DBI_1.1.0               
    ##  [46] edgeR_3.28.1             miniUI_0.1.1.1           Rcpp_1.0.5              
    ##  [49] progress_1.2.2           viridisLite_0.3.0        xtable_1.8-4            
    ##  [52] reticulate_1.16          dqrng_0.2.1              bit_4.0.4               
    ##  [55] rsvd_1.0.3               htmlwidgets_1.5.2        httr_1.4.2              
    ##  [58] RColorBrewer_1.1-2       ellipsis_0.3.1           ica_1.0-2               
    ##  [61] XML_3.99-0.3             R.methodsS3_1.8.1        pkgconfig_2.0.3         
    ##  [64] uwot_0.1.8               deldir_0.1-29            locfit_1.5-9.4          
    ##  [67] tidyselect_1.1.0         rlang_0.4.8              reshape2_1.4.4          
    ##  [70] later_1.1.0.1            munsell_0.5.0            tools_3.6.2             
    ##  [73] generics_0.0.2           RSQLite_2.2.1            ggridges_0.5.2          
    ##  [76] evaluate_0.14            stringr_1.4.0            fastmap_1.0.1           
    ##  [79] yaml_2.2.1               goftest_1.2-2            bit64_4.0.5             
    ##  [82] knitr_1.30               fitdistrplus_1.1-1       purrr_0.3.4             
    ##  [85] RANN_2.6.1               pbapply_1.4-3            future_1.19.1           
    ##  [88] nlme_3.1-149             mime_0.9                 R.oo_1.24.0             
    ##  [91] biomaRt_2.42.1           compiler_3.6.2           curl_4.3                
    ##  [94] beeswarm_0.2.3           plotly_4.9.2.1           png_0.1-7               
    ##  [97] spatstat.utils_1.17-0    tibble_3.0.4             statmod_1.4.34          
    ## [100] stringi_1.5.3            lattice_0.20-41          ProtGenerics_1.18.0     
    ## [103] Matrix_1.2-18            vctrs_0.3.4              pillar_1.4.6            
    ## [106] lifecycle_0.2.0          lmtest_0.9-38            RcppAnnoy_0.0.16        
    ## [109] BiocNeighbors_1.4.2      data.table_1.13.0        cowplot_1.1.0           
    ## [112] bitops_1.0-6             irlba_2.3.3              rtracklayer_1.46.0      
    ## [115] httpuv_1.5.4             patchwork_1.0.1          R6_2.4.1                
    ## [118] promises_1.1.1           KernSmooth_2.23-17       gridExtra_2.3           
    ## [121] vipor_0.4.5              codetools_0.2-16         MASS_7.3-53             
    ## [124] assertthat_0.2.1         rhdf5_2.30.1             openssl_1.4.3           
    ## [127] withr_2.3.0              GenomicAlignments_1.22.1 Rsamtools_2.2.3         
    ## [130] sctransform_0.3.1        GenomeInfoDbData_1.2.2   hms_0.5.3               
    ## [133] mgcv_1.8-33              grid_3.6.2               rpart_4.1-15            
    ## [136] tidyr_1.1.2              rmarkdown_2.4            DelayedMatrixStats_1.8.0
    ## [139] Rtsne_0.15               shiny_1.5.0              ggbeeswarm_0.6.0
