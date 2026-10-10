# S7 Class for the Power Link

Carries the power transformation \\\eta = \theta^{\lambda}\\ on \\(0,
\infty)\\ for a non-zero exponent, with inverse \\\theta =
\eta^{1/\lambda}\\. The exponent is stored in `link_params`, so one
class serves every \\\lambda\\.

At \\\lambda = 0\\ the constructor returns a
[LogLink](https://statmodels7.github.io/linkfunctions7/reference/LogLink.md)
instead, the limit of the Box-Cox transformation \\(\theta^\lambda -
1)/\lambda\\.

## Usage

``` r
PowerLink(
  link_name = character(0),
  link_bounds = integer(0),
  link_params = NULL,
  lambda = integer(0)
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

- lambda:

  The exponent of the transformation.

## Value

An S7 object of class `PowerLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
and carrying its three properties `link_name`, `link_bounds` and
`link_params`. Its `link_bounds` are `c(0, Inf)` and its `link_params`
holds `lambda`.

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
coming back. Both directions are a power, so all ten derivatives are
falling factorials in the exponent. Each is wrapped in
[`na_from()`](https://statmodels7.github.io/linkfunctions7/reference/na_from.md)
because `NA^0` is one in R, which would turn a missing parameter into a
number as soon as an exponent reached zero.

## See also

[`power_link()`](https://statmodels7.github.io/linkfunctions7/reference/power_link.md),
the constructor users call.
