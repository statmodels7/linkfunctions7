# Evaluate Derivative of Inverse Link Function by Order

A convenience router over
[`linkinv`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)
and the four `d*linkinv` generics.

Routes to the correct inverse derivative generic based on order.

## Usage

``` r
linkinvderiv(x, eta, order = 1)

linkinvderiv.link(x, eta, order = 1)
```

## Arguments

- x:

  An object of class `link`.

- eta:

  A numeric vector.

- order:

  An integer (0 to 4).

## Value

A numeric vector of the same length as `eta`: the requested derivative
of \\g^{-1}\\ evaluated there.

A numeric vector of the same length as `eta`.

## Details

This dispatches twice, once on itself and once on the order-specific
generic. Where that matters, call
[`dlinkinv`](https://statmodels7.github.io/linkfunctions7/reference/dlinkinv.md)
and its siblings directly.

## Examples

``` r
lk <- logit_link()
linkinvderiv(lk, 0, order = 0)  # the inverse link itself
#> [1] 0.5
linkinvderiv(lk, 0, order = 1)
#> [1] 0.25

vapply(0:4, function(k) linkinvderiv(lk, 0.5, order = k), numeric(1))
#> [1]  0.62245933  0.23500371 -0.05755679 -0.09635676  0.10475593
```
