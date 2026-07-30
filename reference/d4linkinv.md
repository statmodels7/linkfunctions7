# 4th Derivative of Inverse Link Function

The fourth derivative of \\g^{-1}(\eta)\\ with respect to \\\eta\\.

## Usage

``` r
d4linkinv(x, eta)
```

## Arguments

- x:

  An object of class `link`.

- eta:

  A numeric vector.

## Value

A numeric vector of the same length as `eta`, missing wherever `eta` is.

## See also

[`linkinvderiv`](https://statmodels7.github.io/linkfunctions7/reference/linkinvderiv.md),
which routes to this generic by order.

## Examples

``` r
d4linkinv(logit_link(), 0)
#> [1] 0
d4linkinv(probit_link(), 1)
#> [1] 0.4839414
```
