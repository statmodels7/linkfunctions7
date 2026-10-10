# S7 Class for an Upper Bounded Link

Carries the reflected log on \\(-\infty, u)\\: \\\eta = \log(u -
\theta)\\, with inverse \\\theta = u - e^{\eta}\\.

It is the mirror image of
[`LowerBoundedLink()`](https://statmodels7.github.io/linkfunctions7/reference/LowerBoundedLink.md),
the log of the distance below `upr`. The reflection makes the map
decreasing, so the odd-order forward derivatives change sign relative to
those of
[LowerBoundedLink](https://statmodels7.github.io/linkfunctions7/reference/LowerBoundedLink.md)
while the even ones do not, and every inverse derivative changes sign.

## Usage

``` r
UpperBoundedLink(
  link_name = character(0),
  link_bounds = integer(0),
  link_params = NULL,
  upr = integer(0)
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

- upr:

  The upper endpoint.

## Value

An S7 object of class `UpperBoundedLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
and carrying its three properties `link_name`, `link_bounds` and
`link_params`. Its `link_bounds` are `c(-Inf, upr)`, and its
`link_params` holds `upr`.

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
coming back. The reflection makes the map decreasing, so every inverse
derivative is the negative of
[`exp_floored()`](https://statmodels7.github.io/linkfunctions7/reference/exp_floored.md)
of \\\eta\\, and the forward ones are the log's read at \\\mathrm{upr} -
\theta\\, with the change of sign that the reflection introduces at the
odd orders.

## See also

[`bounded_link()`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.md),
the constructor users call.
