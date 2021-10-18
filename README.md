# MR-Squared
This R package performs multiply robust Mendelian randomization analysis with many invalid instruments.

## Setup
Use the following command in R to install the package:
```
install.package("MRSquared")
```

## Arguments
- `Y` length n continuous outcome vector
- `A` length n continuous exposure vector
- `G.mat` n by K binary SNP matrix
- `k` 2<=k<=K, the number of valid IVs

## Return
Point estimate of causal effect of A on Y and standard error.

## Example
```
library(MRSquared)
```
