# The Range of Predictors a Link Admits

The image of the link's parameter bounds under
[`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md),
which is the set of predictors on which the inverse link is defined.
Where a link maps onto the whole real line the result is `c(-Inf, Inf)`;
otherwise its finite end or ends bound the values that a caller may pass
to
[`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md).

## Usage

``` r
eta_bounds(x)
```

## Arguments

- x:

  A
  [`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
  object.

## Value

A numeric vector of length two, sorted, with `-Inf` or `Inf` in either
position where that end is unbounded.

## Details

A link need not map onto the whole real line: the square root reaches
only the positive half, and so do
[`inverse_link()`](https://statmodels7.github.io/linkfunctions7/reference/inverse_link.md),
[`inverse_sq_link()`](https://statmodels7.github.io/linkfunctions7/reference/inverse_sq_link.md)
and
[`power_link()`](https://statmodels7.github.io/linkfunctions7/reference/power_link.md)
at any non-zero exponent. The bounds are returned sorted, since a
decreasing link reverses them, and are infinite in the directions where
they cannot be established.

## Use outside the package

Inside the package the function keeps a finite-difference grid within
the set on which the inverse link is defined, since a stencil point
outside it returns `NaN` and makes a numerical derivative missing.
Outside the package it serves a related purpose:
[link_bounds()](https://statmodels7.github.io/linkfunctions7/reference/link.md)
gives the set that a link maps **onto**, and a consumer that carries an
unconstrained vector needs the set that it maps **from**.

A family that reads a free vector in \\\mathbb{R}^d\\ and applies an
inverse link to each coordinate needs the map to be defined and
injective on the whole real line, which fails for a link with finite eta
bounds:
[`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)
of a square root link is even, so `-2` and `2` both give 4 and the round
trip returns the absolute value. Such a family can reject the link at
construction by testing `all(is.infinite(eta_bounds(link)))`, and its
error message can then name the link.

## See also

[`link_bounds_clamp()`](https://statmodels7.github.io/linkfunctions7/reference/link_bounds_clamp.md),
which keeps
[`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)'s
result strictly inside the parameter bounds, and
[`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md),
the map whose image this is.

## Examples

``` r
# A link from the whole real line onto the positive half of the theta axis.
eta_bounds(log_link())
#> [1] -Inf  Inf

# The square root reaches only the positive half of the eta axis, so its
# inverse is even there and the round trip returns the absolute value.
eta_bounds(sqrt_link())
#> [1]   0 Inf
linkfun(sqrt_link(), linkinv(sqrt_link(), -2))
#> [1] 2

# Which is the question a consumer carrying an unconstrained vector asks.
links <- list(log_link(), softplus_link(), logit_link(), sqrt_link(),
              inverse_link(), inverse_sq_link(), power_link(0.5))
vapply(links, function(l) all(is.infinite(eta_bounds(l))), logical(1))
#> [1]  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE
```
