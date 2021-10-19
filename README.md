# MRSquared
This R package performs multiply robust Mendelian randomization analysis with many invalid instruments.

## Setup
Use the following command in R to install the package:
```
library(devtools)
devtools::install_github("https://github.com/zhonghualiu/MRSquared")
```

## Arguments
- `Y` a continuous outcome vector of length n
- `A` a continuous exposure vector of length n
- `G.mat` a n by K binary SNP matrix
- `k` 1<=k<=K, the number of valid IVs

## Return
The point estimate of the causal effect of A on Y and standard error.

## Example
```
library(MRSquared)
```
