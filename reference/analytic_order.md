# Highest Analytically Implemented Derivative Order

The largest \\k \le 4\\ for which the link's own class registers a
method for the order-\\k\\ derivative generic; `0` when it registers
none, so that only
[`linkfun`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md)
or
[`linkinv`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)
is available.

## Usage

``` r
analytic_order(x, inverse = FALSE)
```

## Arguments

- x:

  An object of class `link`.

- inverse:

  Logical; `TRUE` to ask about the inverse-link generics.

## Value

An integer between 0 and 4.

## Details

Detection uses the documented S7 property that a method records the
class it was registered on in its `signature` attribute: a method
inherited from the base
[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
class is a fallback, anything else is the link's own. Comparing method
objects with [`identical()`](https://rdrr.io/r/base/identical.html) does
not work for this, because S7 wraps them.

Comparing the recorded *class* with
[`identical()`](https://rdrr.io/r/base/identical.html) does not work
either, and the reason is worth recording because the failure is silent
and spectacular. [`identical()`](https://rdrr.io/r/base/identical.html)
on two S7 class objects is object identity, so it returns `FALSE` for a
class re-created from the same definition — and that is exactly what
happens whenever the package's code is re-evaluated rather than loaded,
as under test coverage instrumentation. The base fallback is then
mistaken for the link's own method, this function returns 4 for a link
that implements nothing, and
[`fallback_deriv`](https://statmodels7.github.io/linkfunctions7/reference/fallback_deriv.md)
walks down through
[`linkderiv`](https://statmodels7.github.io/linkfunctions7/reference/linkderiv.md)
into precisely the chain of nested first differences the design exists
to avoid. It is not a small error: the fourth derivative of the log link
came back wrong by a factor of 900. Classes are therefore compared by
name and package, with identity kept only as a fast path.

The search stops at the first missing order rather than continuing,
since a link that implements the first and the third but not the second
is not a case worth optimising for, and stopping keeps the answer
meaning "everything up to here is exact".

## See also

[`link_fallback_orders`](https://statmodels7.github.io/linkfunctions7/reference/link_fallback_orders.md),
which reports this to the user.
