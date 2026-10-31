# 6th Derivative of a Link Function

Computes the sixth derivative of the link \\g(\theta)\\ with respect to
the parameter, on the parameter scale.

## Usage

``` r
d6linkfun(x, theta)
```

## Arguments

- x:

  An object of class `link`.

- theta:

  A numeric vector of parameter values, inside `x@link_bounds`. The
  domain is not checked: outside it the formula is evaluated as written,
  so the result may be `NaN`, `NA` or an ordinary number, and the
  function does not signal an error.

## Value

A numeric vector of the same length as `theta`, missing wherever `theta`
is.

## Details

The sixth order is reached by the link-scale sixth derivative of a
log-likelihood, which combines the derivatives of the family with the
partial Bell polynomials in the derivatives of the inverse link up to
the sixth.

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

The result is a signed infinity only where the derivative itself exceeds
the largest double, and the powers and products are arranged so that no
intermediate quantity overflows or underflows first. For the log,
lower-bounded and upper-bounded links this happens for \\\theta\\ (or
the distance to the bound) below about \\9.3 \times 10^{-52}\\, and for
the logit, probit, cloglog, loglog and cauchit links at the
corresponding probabilities close to 0.

## See also

[`linkderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkderiv.md),
which routes to this generic by order, and
[`d6linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d6linkinv.md)
for the same order in the other direction.

## Examples

``` r
# The log link's sixth derivative is -120 / t^6, so at theta = 2:
d6linkfun(log_link(), 2) - (-120 / 2^6)
#> [1] 0

# Missingness propagates instead of being filled in.
d6linkfun(logit_link(), c(0.5, NA))
#> [1]  0 NA
```
