
<!-- README.md is generated from README.Rmd. Please edit that file -->

# compassanalyzer

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install the development version of qNCIN like so:

``` r
# install.packages("devtools")
devtools::install_github("thereallda/qNCIN")
```

## Quick Start

### Load package

``` r
library(tidyverse)

# if you do not have `enONE` package installed, run the following code first: 
# devtools::install_github("thereallda/enONE")
if (!requireNamespace("enONE", quietly = TRUE)) { devtools::install_github("thereallda/enONE") }
library(enONE)
library(compassanalyzer)
```

### Load data

``` r
counts_df <- read.csv("data/Counts1.csv", row.names = 1)
meta <- read.csv("data/metadata1.csv")

# metadata for synthetic RNA
syn_id <- paste("syn",3:7, sep = "_")
syn_meta <- data.frame(
  id = syn_id,
  per = c(0.05,0.01,0.20,0,0.10)
)
```

### Filtering

``` r
counts_keep <- enONE::FilterLowExprGene(counts_df, 
                                 group = meta$condition,
                                 min.count = 20)
```

### Create object

``` r
Compass <- createCompass(counts_keep,
                         bio.group = meta$condition,
                         enrich.group = meta$enrich,
                         spike.in.prefix = "^FB",
                         input.id = "Input",
                         enrich.id = "Enrich",
                         synthetic.id = syn_meta$id)
```

### Calculate Ratio

``` r
Compass <- CompassAnalyze(Compass)
#> Global scaling...
#> Linear regression-based adjustment...
#> Computation of NCIN Ratio...
qratio_df <- getRatio(Compass, slot = "sample", filter = T)
head(qratio_df);dim(qratio_df)
#>                       Y1_YCNT    Y2_YCNT    Y3_YCNT
#> ENSMUSG00000100954 0.04605680 0.02392197 0.01336620
#> ENSMUSG00000051285 0.28605972 0.22196503 0.46417798
#> ENSMUSG00000048538 0.06555230 0.08062193 0.10834769
#> ENSMUSG00000057363 0.04013835 0.05196659 0.03195173
#> ENSMUSG00000033021 0.03507386 0.04809227 0.04212718
#> ENSMUSG00000061024 0.04004754 0.06341580 0.06380841
#> [1] 8807    3
```

### Synthetic RNA Calibration curve

``` r
synScatter(ratio.df = qratio_df, syn.meta = syn_meta)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />
