# 4th Derivative of Link Function

The fourth derivative of \\g(\theta)\\ with respect to \\\theta\\.

## Usage

``` r
d4linkfun(x, theta)
```

## Arguments

- x:

  An object of class `link`.

- theta:

  A numeric vector.

## Value

A numeric vector of the same length as `theta`, missing wherever `theta`
is.

## See also

[`linkderiv`](https://statmodels7.github.io/linkfunctions7/reference/linkderiv.md),
which routes to this generic by order.

## Examples

``` r
d4linkfun(logit_link(), 0.5)
#> [1] 0
d4linkfun(log_link(), c(1, 2))
#> [1] -6.000 -0.375
```
