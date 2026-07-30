# Evaluate Derivative of Link Function by Order

A convenience router over
[`linkfun`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md)
and the four `d*linkfun` generics.

Routes to the correct forward derivative generic based on order.

## Usage

``` r
linkderiv(x, theta, order = 1)

linkderiv.link(x, theta, order = 1)
```

## Arguments

- x:

  An object of class `link`.

- theta:

  A numeric vector.

- order:

  An integer (0 to 4).

## Value

A numeric vector of the same length as `theta`: the requested derivative
of \\g\\ evaluated there.

A numeric vector of the same length as `theta`.

## Details

This dispatches twice, once on itself and once on the order-specific
generic. Where that matters, call
[`dlinkfun`](https://statmodels7.github.io/linkfunctions7/reference/dlinkfun.md)
and its siblings directly.

## Examples

``` r
lk <- logit_link()
linkderiv(lk, 0.5, order = 0)   # the link itself
#> [1] 0
linkderiv(lk, 0.5, order = 1)
#> [1] 4

# every order at once
vapply(0:4, function(k) linkderiv(lk, 0.25, order = k), numeric(1))
#> [1]    -1.098612     5.333333   -14.222222   132.740741 -1517.037037
```
