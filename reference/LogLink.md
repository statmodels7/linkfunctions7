# S7 Class for the Logarithmic Link

Carries the log transformation \\\eta = \log\theta\\ on \\(0, \infty)\\,
with inverse \\\theta = e^{\eta}\\. It is the standard link for a
positive parameter such as a scale or a rate.

The inverse is floored at
[`exp_floor()`](https://statmodels7.github.io/linkfunctions7/reference/exp_floor.md),
so a parameter reported by this link is never exactly zero and can be
divided into.

## Usage

``` r
LogLink(link_name = character(0), link_bounds = integer(0), link_params = NULL)
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

An S7 object of class `LogLink`, inheriting from
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
coming back. The forward derivatives are
\\(-1)^{k-1}(k-1)!\\\theta^{-k}\\, written out. Every inverse derivative
is
[`exp_floored()`](https://statmodels7.github.io/linkfunctions7/reference/exp_floored.md)
of \\\eta\\, the exponential being its own derivative to every order,
and the floor is what keeps \\\theta\\ strictly inside \\(0, \infty)\\.

## See also

[`log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.md),
the constructor users call.
