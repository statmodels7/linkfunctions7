# S7 Class for the Sqrt Link

Carries the square-root transformation \\\eta = \sqrt{\theta}\\ on \\(0,
\infty)\\, with inverse \\\theta = \eta^2\\.

Its image is \\(0, \infty)\\. On the whole real line the inverse link is
not one-to-one, since \\\eta\\ and \\-\eta\\ give the same \\\theta\\;
[`eta_bounds()`](https://statmodels7.github.io/linkfunctions7/reference/eta_bounds.md)
returns the range on which it is.

## Usage

``` r
SqrtLink(
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

An S7 object of class `SqrtLink`, inheriting from
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
The forward derivatives are half-integer falling factorials. The inverse
map is \\\eta^2\\, so its derivatives terminate: the second is the
constant two and the third to sixth are exactly zero.
[`const_like()`](https://statmodels7.github.io/linkfunctions7/reference/const_like.md)
builds all of them, so that a missing value still propagates.

## See also

[`sqrt_link()`](https://statmodels7.github.io/linkfunctions7/reference/sqrt_link.md),
the constructor that users call.
