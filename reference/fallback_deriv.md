# The Body Shared by Every Numerical Fallback

Computes the order-`order` derivative of a link, in either direction, by
differentiating once the highest order the link supplies analytically.

## Usage

``` r
fallback_deriv(x, v, order, inverse)
```

## Arguments

- x:

  An object of class `link`.

- v:

  A numeric vector: \\\theta\\ going forward, \\\eta\\ coming back.

- order:

  The derivative order wanted, 1 to 4.

- inverse:

  Logical; `TRUE` for the inverse-link direction.

## Value

A numeric vector of the same length as `v`.

## Details

The whole design is in the two lines that pick `m` and `gap`: never
differentiate numerically more than the number of orders actually
missing. A link analytic to the second order asks for a first difference
to reach the third, not three; a link supplying nothing but
[`linkfun`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md)
is the only case in which a fourth-order stencil is applied to the
function itself.

Note the recursion is only apparent. The base function is fetched
through
[`linkderiv`](https://statmodels7.github.io/linkfunctions7/reference/linkderiv.md),
which dispatches to the link's own method for an order it implements —
so the chain always terminates on analytic code, and never on another
fallback.
