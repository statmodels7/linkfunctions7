# Evaluate Forward Link Function

Evaluate Forward Link Function

## Usage

``` r
linkfun(x, theta)
```

## Arguments

- x:

  An object of class `link`.

- theta:

  A numeric vector of parameters.

## Value

A numeric vector of the linear predictor.

## Examples

``` r
linkfun(logit_link(), c(0.25, 0.5, 0.75))
#> [1] -1.098612  0.000000  1.098612
linkfun(log_link(), c(1, exp(1)))
#> [1] 0 1
```
