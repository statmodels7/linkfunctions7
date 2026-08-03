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

The recorded class is compared by name and package rather than with
[`identical()`](https://rdrr.io/r/base/identical.html), since
[`identical()`](https://rdrr.io/r/base/identical.html) on S7 class
objects tests object identity and returns `FALSE` for a class re-created
from the same definition, as happens when the package's code is
re-evaluated under coverage instrumentation; identity is kept only as a
fast path.

The search stops at the first missing order, so the answer always means
that every order up to it is analytic.

## See also

[`link_fallback_orders`](https://statmodels7.github.io/linkfunctions7/reference/link_fallback_orders.md),
which reports this to the user.
