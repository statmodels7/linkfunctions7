# S7 Class for the InverseSq Link

Carries the inverse-square transformation \\\eta = 1/\theta^2\\ on \\(0,
\infty)\\, with inverse \\\theta = 1/\sqrt{\eta}\\.

It is the canonical link of the inverse Gaussian family. Its image is
\\(0, \infty)\\, so the inverse link is defined only for a positive
linear predictor, and like
[`inverse_link()`](https://statmodels7.github.io/linkfunctions7/reference/inverse_link.md)
the map is decreasing.

## Usage

``` r
InverseSqLink(
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

An S7 object of class `InverseSqLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
and carrying its three properties `link_name`, `link_bounds` and
`link_params`. Its `link_bounds` are `c(0, Inf)` and it carries no link
parameters.

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
The two directions are not the same map (\\1/\theta^2\\ in the forward
direction and \\\eta^{-1/2}\\ in the inverse one), so the two sets of
derivatives are written out separately, with integer falling factorials
in the forward direction and half-integer powers in the inverse one.

## See also

[`inverse_sq_link()`](https://statmodels7.github.io/linkfunctions7/reference/inverse_sq_link.md),
the constructor that users call.
