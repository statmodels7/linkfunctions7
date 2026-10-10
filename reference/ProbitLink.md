# S7 Class for the Probit Link

Carries the probit transformation \\\eta = \Phi^{-1}(\theta)\\ on \\(0,
1)\\, with \\\Phi\\ the standard normal distribution function and the
inverse \\\theta = \Phi(\eta)\\.

It is symmetric about \\\theta = 1/2\\, like the logit, and it reaches
its bounds faster, because the tails of the normal distribution are
lighter than those of the logistic.

## Usage

``` r
ProbitLink(
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

An S7 object of class `ProbitLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
and carrying its three properties `link_name`, `link_bounds` and
`link_params`. Its `link_bounds` are `c(0, 1)` and it carries no link
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
coming back.
[`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md)
and
[`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)
delegate to [`stats::qnorm()`](https://rdrr.io/r/stats/Normal.html) and
[`stats::pnorm()`](https://rdrr.io/r/stats/Normal.html). The ten
derivatives come from a compiled kernel, one call per order and
direction.

## See also

[`probit_link()`](https://statmodels7.github.io/linkfunctions7/reference/probit_link.md),
the constructor users call.
