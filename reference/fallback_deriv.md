# The Body Shared by Every Numerical Fallback

Computes the order-`order` derivative of a link, in either direction, by
one numerical differentiation of the highest order that the link
supplies analytically.

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

  The derivative order wanted, 1 to 5.

- inverse:

  Logical; `TRUE` for the inverse-link direction.

## Value

A numeric vector of the same length as `v`.

## Details

The function takes the highest analytic order \\m\\ below the requested
one and applies a single central stencil of order `order - m` to it,
never a chain of lower-order stencils. For a link analytic to the second
order, the third derivative is one first difference of the second; only
for a link that supplies nothing but
[`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md)
is a stencil of the full order applied to the function itself.

The recursion is only apparent. The base function is obtained through
[`linkderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkderiv.md)
or
[`linkinvderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinvderiv.md),
which dispatch to the link's own method for an order that the link
implements, so the chain always ends on analytic code and never on
another fallback.

## Methods

The ten registrations on the base class
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
have this function as their whole body:
[`dlinkfun()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkfun.md)
through
[`d5linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/d5linkfun.md)
going out and
[`dlinkinv()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkinv.md)
through
[`d5linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d5linkinv.md)
coming back, each passing its order and its direction. A link inherits
them for the orders it does not implement itself, so a link defined with
nothing but
[`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md)
and
[`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)
still has a method for every derivative generic. S7 requires the formals
of a method to match those of the generic, and the two directions name
their argument differently, so the ten wrappers are written out instead
of generated.
