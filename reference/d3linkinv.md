# 3rd Derivative of Inverse Link Function

The third derivative of \\g^{-1}(\eta)\\ with respect to \\\eta\\.

## Usage

``` r
d3linkinv(x, eta)
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
d3linkinv(logit_link(), 0)
#> [1] -0.125
d3linkinv(probit_link(), 1)
#> [1] 0
```
