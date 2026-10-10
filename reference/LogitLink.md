# S7 Class for the Logit Link

Carries the logit transformation \\\eta = \log(\theta/(1-\theta))\\ on
\\(0, 1)\\, with inverse \\\theta = 1/(1+e^{-\eta})\\. It is the
canonical link for a probability, and its linear predictor is the
log-odds.

Every inverse derivative is a polynomial in \\\theta\\ itself, the first
being \\\theta(1-\theta)\\, the variance of a Bernoulli.

## Usage

``` r
LogitLink(
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

An S7 object of class `LogitLink`, inheriting from
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
delegate to [`stats::qlogis()`](https://rdrr.io/r/stats/Logistic.html)
and [`stats::plogis()`](https://rdrr.io/r/stats/Logistic.html), which
remain accurate near both boundaries. The ten derivatives come from a
compiled kernel, one call per order and direction.

## See also

[`logit_link()`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.md),
the constructor users call.
