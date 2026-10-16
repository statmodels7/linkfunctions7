# S7 Class for a Lower Bounded Link

Carries the shifted log on \\(l, \infty)\\: \\\eta = \log(\theta - l)\\,
with inverse \\\theta = l + e^{\eta}\\.

Every derivative is the log link's, because the shift is a constant that
differentiates away, and the exponential is floored at
[`exp_floor()`](https://statmodels7.github.io/linkfunctions7/reference/exp_floor.md)
so the parameter never reaches the bound exactly.

## Usage

``` r
LowerBoundedLink(
  link_name = character(0),
  link_bounds = integer(0),
  link_params = NULL,
  lwr = integer(0)
)
```

## Arguments

- link_name:

  A character string naming the link, set by the constructor and shown
  by [`print()`](https://rdrr.io/r/base/print.html).

- link_bounds:

  A length-two numeric vector, the open interval in which the parameter
  lies. The constructor sets it, and the Value section gives its value
  for this link.

- link_params:

  A list of the link's own parameters, empty where it has none. The
  constructor sets it.

- lwr:

  The lower endpoint.

## Value

An S7 object of class `LowerBoundedLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
and carrying its three properties `link_name`, `link_bounds` and
`link_params`. Its `link_bounds` are `c(lwr, Inf)`, and its
`link_params` holds `lwr`.

## Methods

Twelve methods are registered on this class:
[`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md)
and
[`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md),
the five derivatives of the link,
[`dlinkfun()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkfun.md)
through
[`d5linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/d5linkfun.md),
and the five derivatives of the inverse link,
[`dlinkinv()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkinv.md)
through
[`d5linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d5linkinv.md).
The forward derivatives are the log link's read at \\\theta -
\mathrm{lwr}\\, because the shift is a constant that differentiates
away. Every inverse derivative is
[`exp_floored()`](https://statmodels7.github.io/linkfunctions7/reference/exp_floored.md)
of \\\eta\\ itself, because the exponential is its own derivative to
every order.

## See also

[`bounded_link()`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.md),
the constructor that users call.
