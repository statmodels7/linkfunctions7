# 2nd Derivative of an Inverse Link Function

The second derivative of the inverse link \\g^{-1}(\eta)\\ with respect
to the linear predictor. It enters the chain rule (the formula of Faa di
Bruno) that carries the higher derivatives of the log-likelihood from
\\\theta\\ onto \\\eta\\.

## Usage

``` r
d2linkinv(x, eta)
```

## Arguments

- x:

  An object of class `link`.

- eta:

  A numeric vector of linear predictors, inside the range that
  [`eta_bounds()`](https://statmodels7.github.io/linkfunctions7/reference/eta_bounds.md)
  returns for the link. The derivative is computed from `eta` directly,
  without the clamp that
  [`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)
  applies.

## Value

A numeric vector of the same length as `eta`, missing wherever `eta` is.

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

[`linkinvderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinvderiv.md),
which routes to this generic by order, and
[`d2linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/d2linkfun.md)
for the same order in the other direction.

## Examples

``` r
# Every derivative of exp is exp, so the log link's inverse gives the
# same number at every order.
d2linkinv(log_link(), 1) - exp(1)
#> [1] 0

# The logit's inverse derivatives are polynomials in theta. At eta = 0 the
# logistic is symmetric about 1/2, so its even-order derivatives vanish
# there while the first is the Bernoulli variance, 1/4.
d2linkinv(logit_link(), 0)
#> [1] 0
```
