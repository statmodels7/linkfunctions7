# 2nd Derivative of Link Function

The second derivative of \\g(\theta)\\ with respect to \\\theta\\.

## Usage

``` r
d2linkfun(x, theta)
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
d2linkfun(logit_link(), 0.5)
#> [1] 0
d2linkfun(log_link(), c(1, 2))
#> [1] -1.00 -0.25
```
