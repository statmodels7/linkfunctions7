# 2nd Derivative of a Link Function

The second derivative of the link \\g(\theta)\\ with respect to the
parameter, on the parameter scale.

## Usage

``` r
d2linkfun(x, theta)
```

## Arguments

- x:

  An object of class `link`.

- theta:

  A numeric vector of parameter values, inside `x@link_bounds`. The
  domain is not checked: outside it the formula is evaluated as written,
  so the result may be `NaN`, `NA` or an ordinary number, and no error
  is signaled.

## Value

A numeric vector of the same length as `theta`, missing wherever `theta`
is.

## Details

Every link has a method for this generic. If the class of a link
registers no method for it, the numerical method of the base class is
used: it applies a single central stencil to the highest order that the
link supplies analytically, never a chain of lower-order differences.
[`link_fallback_orders()`](https://statmodels7.github.io/linkfunctions7/reference/link_fallback_orders.md)
reports which orders of a given link are exact, and
[`check_link()`](https://statmodels7.github.io/linkfunctions7/reference/check_link.md)
leaves a fallback order unchecked, because comparing it with a finite
difference of the order below would repeat the same computation and
could not detect an error.

In code where speed matters, call this generic directly:
[`linkderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkderiv.md)
and
[`linkinvderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinvderiv.md)
route by order and so dispatch twice, once on themselves and once here,
and on a short vector the routed call takes two to three times as long
as the direct one.

## See also

[`linkderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkderiv.md),
which routes to this generic by order, and
[`d2linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d2linkinv.md)
for the same order in the other direction.

## Examples

``` r
# The log link's forward derivatives are -1 / t^2, so at theta = 2:
d2linkfun(log_link(), 2) - (-1 / 2^2)
#> [1] 0

# Missingness propagates instead of being filled in.
d2linkfun(logit_link(), c(0.5, NA))
#> [1]  0 NA
```
