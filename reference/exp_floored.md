# A Floored Exponential

`exp(eta)`, bounded below by
[`exp_floor`](https://statmodels7.github.io/linkfunctions7/reference/exp_floor.md).

## Usage

``` r
exp_floored(eta)
```

## Arguments

- eta:

  A numeric vector of linear predictors.

## Value

A numeric vector, never smaller than
[`exp_floor`](https://statmodels7.github.io/linkfunctions7/reference/exp_floor.md).
