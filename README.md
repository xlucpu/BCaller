# BCaller

<!-- badges: start -->
<!-- badges: end -->

The goal of BCaller is to calculate an immune-related genes pair index (IRGPI) from single-sample perspective for bladder cancer only. The IRGPI risk score could be used to stratify prognosis of bladder cancer.

## Citation

For now, you can cite the following bioRxiv preprint

## Installation

You may install this package with:

``` r
if (!require("devtools")) 
    install.packages("devtools")
devtools::install_github("xlucpu/BCaller")
```

## Example
``` r
# load R package and internal data set
library(BCaller)
load(system.file("extdata", "demo.RData", package = "BCaller", mustWork = TRUE)) # load example data

# calculate IRGPI
IRGPI  <- calIRGPI(expr = demo) 

# print
head(IRGPI)
#        SampleID     IRGPI
# 1 BLCA-A9KO-01A 1.4756139
# 2 BLCA-A9KP-01A 1.5510460
# 3 BLCA-A9KQ-01A 0.5828771
# 4 BLCA-A9KR-01A 1.9368114
# 5 BLCA-A9KT-01A 1.6784781
# 6 BLCA-A9KW-01A 2.5174033
```

