# S7 Class for the Softplus Link

Carries the softplus transformation on \\(0, \infty)\\, whose inverse
\\\theta = \log(1 + e^{a\eta})/a\\ is a smooth approximation of
\\\max(0, \eta)\\ that sharpens as the scale \\a\\ grows.

It is an alternative to the log link for a positive parameter. For a
large positive \\\eta\\ the softplus grows linearly, \\\theta \approx
\eta\\, where the inverse of the log link grows exponentially; for a
large negative \\\eta\\ both approach zero exponentially. The scale is
stored in `link_params`, so one class serves every \\a\\.

## Usage

``` r
SoftplusLink(
  link_name = character(0),
  link_bounds = integer(0),
  link_params = NULL,
  a = integer(0)
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

- a:

  The scale parameter, strictly positive.

## Value

An S7 object of class `SoftplusLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
and carrying its three properties `link_name`, `link_bounds` and
`link_params`. Its `link_bounds` are `c(0, Inf)` and its `link_params`
holds `a`.

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
The softplus is an antiderivative of the logistic, so its inverse-link
derivatives are the logistic ones shifted a place, \\h^{(k+1)} = a^k
\sigma^{(k)}\\, and they reuse the same polynomials as the logit. The
forward set comes from a compiled kernel written in \\u =
-\mathrm{expm1}(-a\theta)\\, which is finite where the
\\\mathrm{expm1}(a\theta)\\ form overflows.

## See also

[`softplus_link()`](https://statmodels7.github.io/linkfunctions7/reference/softplus_link.md),
the constructor that users call.
