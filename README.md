# ccc: differential Cell-Cell Communication analysis

![Under Development](https://img.shields.io/badge/status-under%20development-orange)

## Installation

To install `ccc` from GitHub, use the following command:

```r
devtools::install_github("CHENGHAO-WANG/ccc")
```

## Examples

```r
## Run with example simulated data
data(data.sim, package = "ccc")
data(lr.sim, package = "ccc")
expression_matrix <- log_normalize(data.sim$counts)
metadata <- data.sim$metadata

# Run sequentially
a <- ccc_diff(
  expression_matrix = expression_matrix, metadata = metadata,
  id_col = "sample", lr = lr.sim, sender = "CT1", receiver = c("CT2", "CT3"),
  contrast = c(grp2 = 1, grp1 = -1), lmm_re = TRUE, logmm_re = TRUE
)
head(a$summary)
head(a$test)

# Run in parallel with 4 cores
library(future)
oplan <- plan(multisession, workers = 4L)
a <- ccc_diff(
  expression_matrix = expression_matrix, metadata = metadata,
  id_col = "sample", lr = lr.sim, sender = "CT1", receiver = c("CT2", "CT3"),
  contrast = c(grp2 = 1, grp1 = -1), lmm_re = TRUE, logmm_re = TRUE
)
plan(oplan)
head(a$summary)
head(a$test)
```