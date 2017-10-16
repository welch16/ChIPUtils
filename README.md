
## ChIPUtils

The idea of this package is to have several of the functions used in a
exploratory analysis of ChIP-seq and ChIP-exo data all together in one
package.

We also included some commonly used quality control measures:

- PCR bottleneck coefficient

- Strand cross-correlation

### Installation

To install the package it is better to use:

```
## install.packages("devtools") ## to install devtools
devtools::install_github("welch16/ChIPUtils",ref = "devel")
## or 
devtools::install_github("welch16/ChIPUtils",ref = "devel",build_vignettes = TRUE) 
```
