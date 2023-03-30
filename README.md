# Ensemble Fractional Imputation

Ensemble Fractional Imputation is a nonparametric imputation method that
can handle high dimensional categorical data. It uses the idea of random
forest to fit multiple reduced models and then combine multiple models
using model weights.

## Installation

``` r
library(devtools)
install_github("yonghyun-K/EFI", dependencies = T)

library(EFI)
```

## Example commands

``` r
# Import data and generate missingness.
data(HairEyeColor)
p = 3
Y = do.call("rbind", apply(as.data.frame.table(HairEyeColor), 1, 
  function(x) matrix(rep(x[1:p], each = x[p+1]), nc = p)))
Y = as.data.frame(Y)
for(k in 1:p){
  Y[[k]] = factor(Y[[k]])
}
names(Y) <- names(dimnames(HairEyeColor))
(n = nrow(Y)); sum(HairEyeColor)
delta = matrix(rbinom(n * p, 1, 0.5), nr = n, nc = p)
Y[delta == 0] = NA

# Ensemble Fractional Imputation.
cand.edges = as.list(data.frame(combn(p, 2)))
dp = doublep(Y, cand.edges, freq = F)
EFI = efi(Y, dp, freq = F)
estimate(EFI, "(Hair == \"Black\") & (Eye == \"Brown\")")
estimate(EFI, "(Hair == \"Black\") & (Sex == \"Male\")")
```

## Externel Links
- [CRAN Task View: Missing Data](https://cran.r-project.org/web/views/MissingData.html)

- [FHDI](https://github.com/cran/FHDI)

- [mice](https://github.com/amices/mice)
- https://stefvanbuuren.name/fimd/

- [missForest](https://github.com/stekhoven/missForest)

- [GAIN](https://github.com/jsyoon0823/GAIN)
