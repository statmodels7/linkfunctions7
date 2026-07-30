# Evaluate Inverse Link Function

Evaluate Inverse Link Function

## Usage

``` r
linkinv(x, eta)
```

## Arguments

- x:

  An object of class `link`.

- eta:

  A numeric vector of linear predictors.

## Value

A numeric vector of probabilities/means.

## Examples

``` r
linkinv(logit_link(), c(-1, 0, 1))
#> [1] 0.2689414 0.5000000 0.7310586
linkinv(log_link(), c(0, 1))
#> [1] 1.000000 2.718282
```
