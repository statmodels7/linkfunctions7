# S7 Class for the Inverse Link

Carries the reciprocal transformation \\\eta = 1/\theta\\ on \\(0,
\infty)\\, which is its own inverse.

It is the canonical link of the Gamma family. Its image is \\(0,
\infty)\\, so the inverse link is defined only for a positive linear
predictor, and the map is decreasing: a larger parameter gives a smaller
\\\eta\\.

## Usage

``` r
InverseLink(
  link_name = character(0),
  link_bounds = integer(0),
  link_params = NULL
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

## Value

An S7 object of class `InverseLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
and carrying its three properties `link_name`, `link_bounds` and
`link_params`. Its `link_bounds` are `c(0, Inf)` and it carries no link
parameters.

## Methods

Fourteen methods are registered on this class:
[`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md)
and
[`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md),
the six derivatives of the link,
[`dlinkfun()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkfun.md)
through
[`d6linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/d6linkfun.md),
and the six derivatives of the inverse link,
[`dlinkinv()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkinv.md)
through
[`d6linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d6linkinv.md).
The link is its own inverse, so the two directions carry the same
expressions, \\1/\theta\\ and \\1/\eta\\, and the six derivatives
\\(-1)^k k!\\z^{-(k+1)}\\ are written out in both.

## See also

[`inverse_link()`](https://statmodels7.github.io/linkfunctions7/reference/inverse_link.md),
the constructor that users call.
