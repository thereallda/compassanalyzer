<!-- README.md is generated from README.Rmd. Please edit that file -->

# compassanalyzer <a href="https://github.com/thereallda/compassanalyzer"><img src="man/figures/logo.jpg" alt="CompassAnalyzer" align="right" height="138"/></a>

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install the development version of compassanalyzer like so:

```R
# install.packages("devtools")
devtools::install_github("thereallda/compassanalyzer")
```

## Quick Start

### Load package

```R
library(tidyverse)

# if you do not have `enONE` package installed, run the following code first: 
# devtools::install_github("thereallda/enONE")
if (!requireNamespace("enONE", quietly = TRUE)) { devtools::install_github("thereallda/enONE") }
library(enONE)
library(compassanalyzer)
```

### Load data

```R
counts_df <- read.csv("data/Counts_demo.csv", row.names = 1)
meta <- read.csv("data/meta_demo.csv")
head(meta)
#>        id   condition replicate enrich biology
#> 1 Y1.CRNT Young.Input         1  Input   Young
#> 2 Y2.CRNT Young.Input         2  Input   Young
#> 3 Y3.CRNT Young.Input         3  Input   Young
#> 4 O1.CRNT   Old.Input         1  Input     Old
#> 5 O2.CRNT   Old.Input         2  Input     Old
#> 6 O3.CRNT   Old.Input         3  Input     Old
# rownames of metadata should be consistent with the colnames of counts_mat
rownames(meta) <- meta$id

# metadata for synthetic RNA
syn_id <- paste("syn",1:4, sep = "_")
syn_meta <- data.frame(
  id = syn_id,
  per = c(0.05,0.01,0.20,0)
)
```

### Filtering

```R
counts_keep <- enONE::FilterLowExprGene(counts_df,
                                        group = meta$condition,
                                        min.count = 20)
```

### Create object

```R
Compass <- createCompass(counts_keep,
                         col.data = meta,
                         spike.in.prefix = "^FB",
                         input.id = "Input",
                         enrich.id = "Enrich",
                         synthetic.id = syn_meta$id)
```

### Use enONE to perform RUV normalization

```R
Enone <- createEnone(data = counts_keep,
                     col.data = meta,
                     spike.in.prefix = "FB",
                     input.id = "Input",
                     enrich.id = "Enrich"
)
Enone <- enONE(Enone,
               scaling.method = c("TMM"),
               ruv.norm = TRUE, ruv.k = 3,
               eval.pam.k = 2:4, eval.pc.n = 3,
               return.norm = TRUE
)
#> Gene set selection for normalization and assessment...
#> - The number of negative control genes for normalization: 1000 
#> - Estimate dispersion & Fit GLM... 
#> - Testing differential genes... 
#> - The number of positive evaluation genes: 500 
#> - Estimate dispersion & Fit GLM... 
#> - Testing differential genes... 
#> - The number of negative evaluation genes: 500 
#> - Estimate dispersion & Fit GLM... 
#> - Testing differential genes... 
#> Apply normalization...
#> - Scaling... 
#> - Regression-based normalization... 
#> Perform assessment...
# use RUVs_k3
norm.factors <- enONE::getFactor(Enone, slot="sample", method="TMM_RUVs_k3")
names(norm.factors)
#> [1] "normFactor"   "adjustFactor" "alpha"
```

### Calculate Ratio

RUV factors from enONE can be passed into CompassAnalyze for ratio calculation.

To obtain further accurate NCIN capping ratio, `ratio.shrinkage = TRUE` can be applied.

```R
Compass <- CompassAnalyze(Compass,
                           adjust = TRUE,
                           prop.top.enrich = 0.8,
                           decreasing = TRUE,
                           pseudo.count = 1,
                           enone.ruv.factor = norm.factors,
                           ratio.shrinkage = TRUE)
#> Global scaling...
#> Adjustment...
#> Computation of NCIN Ratio...
# prior ratio
ratio_pri_df <- getRatio(Compass, slot = "sample", ratio.shrinkage = F, filter = T)
head(ratio_pri_df)
#>                       Y1.YCNT    Y2.YCNT   Y3.YCNT    O1.YCNT    O2.YCNT
#> ENSMUSG00000097426 0.48710389 0.99000476 0.4936661 0.73750361 0.24354681
#> ENSMUSG00000102095 0.82912948 0.85157477 0.4795359 0.35906067 0.77368374
#> ENSMUSG00000100954 0.04565622 0.05748047 0.1249243 0.02888247 0.02975081
#> ENSMUSG00000051285 0.89163156 0.61455594 0.6860022 0.88205010 0.91585482
#> ENSMUSG00000048538 0.28660934 0.15923379 0.1990585 0.24455298 0.16976447
#> ENSMUSG00000057363 0.32945563 0.33476174 0.3458781 0.10474213 0.19740588
#>                       O3.YCNT
#> ENSMUSG00000097426 0.43702484
#> ENSMUSG00000102095 0.36310817
#> ENSMUSG00000100954 0.02574135
#> ENSMUSG00000051285 0.39694064
#> ENSMUSG00000048538 0.10497186
#> ENSMUSG00000057363 0.26909517

# To get shrunk ratio
qratio_ls <- getRatio(Compass, slot = "sample", ratio.shrinkage = T, filter = T)
names(qratio_ls)
#> [1] "Old"   "Young"
head(qratio_ls$Young)
#>               GeneID    Y1.YCNT    Y2.YCNT    Y3.YCNT ratio.shrunk
#> 2 ENSMUSG00000000049 0.04597514 0.03833789 0.04767986   0.04399781
#> 3 ENSMUSG00000000056 0.62025719 0.45506273 0.54518925   0.54031978
#> 4 ENSMUSG00000000078 0.31871490 0.19851584 0.33836185   0.28532273
#> 5 ENSMUSG00000000085 0.79234184 0.62394300 0.78914356   0.73589163
#> 6 ENSMUSG00000000088 0.05789702 0.05725308 0.05573969   0.05696327
#> 7 ENSMUSG00000000120 0.28749458 0.28231165 0.23467812   0.26817968
#>   ratio.shrunk.sd
#> 2    0.0028723388
#> 3    0.0477495990
#> 4    0.0437030429
#> 5    0.0557679713
#> 6    0.0006394031
#> 7    0.0168079741
```

### Synthetic RNA Calibration curve

```R
syn_ratio <- synRatio(Compass, ratio.shrinkage=TRUE, syn.meta = syn_meta)
head(syn_ratio)
#>   GeneID      ratio    ratio.sd     ratio.se  per
#> 1  syn_1 0.06052904 0.006964224 0.0028431327 0.05
#> 2  syn_2 0.03576192 0.002397542 0.0009787922 0.01
#> 3  syn_3 0.20577485 0.026321464 0.0107456927 0.20
#> 4  syn_4 0.01825970 0.001318250 0.0005381733 0.00
synScatter(ratio.df = syn_ratio, syn.meta = syn_meta)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />
