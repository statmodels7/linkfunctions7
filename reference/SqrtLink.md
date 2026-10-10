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
  lies. Set by the constructor; see Value for this link's.

- link_params:

  A list of the link's own parameters, empty where it has none. Set by
  the constructor.

## Value

An S7 object of class `SqrtLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
and carrying its three properties `link_name`, `link_bounds` and
`link_params`. Its `link_bounds` are `c(0, Inf)` and it carries no link
parameters.

## Methods

Twelve methods are registered on this class:
[`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md)
and
[`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md),
and the five derivative orders in each direction,
[`dlinkfun()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkfun.md)
through
[`d5linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/d5linkfun.md)
going out and
[`dlinkinv()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkinv.md)
through
[`d5linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d5linkinv.md)
coming back. The forward derivatives are half-integer falling
factorials. The inverse map is \\\eta^2\\, so its derivatives terminate:
the second is the constant two and the third to fifth are exactly zero,
all built by
[`const_like()`](https://statmodels7.github.io/linkfunctions7/reference/const_like.md)
so that a missing value still propagates.

## See also

[`sqrt_link()`](https://statmodels7.github.io/linkfunctions7/reference/sqrt_link.md),
the constructor users call.
