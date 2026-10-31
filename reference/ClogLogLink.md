# S7 Class for the Complementary Log-Log Link

Carries the complementary log-log transformation \\\eta =
\log(-\log(1-\theta))\\ on \\(0, 1)\\, with inverse \\\theta = 1 -
\exp(-e^{\eta})\\.

Unlike the logit and the probit it is asymmetric about \\\theta = 1/2\\,
approaching one faster than zero, so it is the link of a
proportional-hazards model for a binary outcome.

## Usage

``` r
ClogLogLink(
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

An S7 object of class `ClogLogLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
and carrying its three properties `link_name`, `link_bounds` and
`link_params`. Its `link_bounds` are `c(0, 1)` and it carries no link
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
[`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md)
is written through [`log1p()`](https://rdrr.io/r/base/Log.html) rather
than as \\\log(-\log(1 - \theta))\\, which rounds to \\-\infty\\ for a
small \\\theta\\ where the true value is finite and representable. The
twelve derivatives come from a compiled kernel.

## See also

[`cloglog_link()`](https://statmodels7.github.io/linkfunctions7/reference/cloglog_link.md),
the constructor that users call.
