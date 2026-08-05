# 1st Derivative of Inverse Link Function

The first derivative of \\g^{-1}(\eta)\\ with respect to \\\eta\\.

## Usage

``` r
dlinkinv(x, eta)
```

## Arguments

- x:

  An object of class `link`.

- eta:

  A numeric vector.

## Value

A numeric vector of the same length as `eta`, missing wherever `eta` is.

## Details

This and its higher-order siblings are the generics a modeling routine
working on the unconstrained scale actually wants. Call them directly
rather than through
[`linkinvderiv`](https://statmodels7.github.io/linkfunctions7/reference/linkinvderiv.md)
in a hot loop: the router dispatches once on itself and then again on
the order-specific generic, which is about a third of the cost of the
call.

## See also

[`linkinvderiv`](https://statmodels7.github.io/linkfunctions7/reference/linkinvderiv.md),
which routes to this generic by order.

## Examples

``` r
dlinkinv(logit_link(), 0)      # p(1 - p) at p = 0.5
#> [1] 0.25
dlinkinv(log_link(), c(0, 1))
#> [1] 1.000000 2.718282
```
