# 2nd Derivative of Inverse Link Function

The second derivative of \\g^{-1}(\eta)\\ with respect to \\\eta\\.

## Usage

``` r
d2linkinv(x, eta)
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
d2linkinv(logit_link(), 0)     # zero, by symmetry about eta = 0
#> [1] 0
d2linkinv(probit_link(), 1)
#> [1] -0.2419707
```
