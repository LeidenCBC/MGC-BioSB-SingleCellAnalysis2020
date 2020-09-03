## Required packages

### CRAN
ggplot2 (at least version 3.0.0)

tidyverse (version 1.2.1)

Seurat (at least version 3.0.0)

igraph (version 1.2.4.1)


cowplot (at least version 0.9.4)

Matrix (version 1.2-17)

### Bioconductor
monocle (version 2.10.1)

biomaRt (version 2.38.0)

destiny (version 2.12.0)

scater (version 1.10.1)

scran (version 1.10.2)

BiocParallel (version 1.16.6)

BiocNeighbors (version 1.0.0)


When using the install.packages() function to download a package, it could be that by default an older version was downloaded instead of the newer version that is needed. For example, just using install.packages(‘Seurat’) could install version 2 instead of version 3. If this happens, you can define the correct repository while installing the package (e.g.(install.pacakges(‘Seurat’, repo = “https://cran.r-project.org”)).

When installing or loading a package (e.g. library(Seurat)), this could fail because of errors in other packages. Using remove.packages(‘package’) this package that gives the error can be removed and installing  it again using install.packages() might work.
