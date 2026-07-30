# 3rd Derivative of Link Function

The third derivative of \\g(\theta)\\ with respect to \\\theta\\.

## Usage

``` r
d3linkfun(x, theta)
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
d3linkfun(logit_link(), 0.5)
#> [1] 32
d3linkfun(log_link(), c(1, 2))
#> [1] 2.00 0.25
```
