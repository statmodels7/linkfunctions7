# S7 Class for a Doubly Bounded Link

Carries the scaled logit on a finite interval \\(l, u)\\: the parameter
is mapped to its position \\p = (\theta - l)/(u - l)\\ within the
interval and that position is carried by the logit, so \\\eta =
\log(p/(1-p))\\.

Its derivatives are the logistic polynomials scaled by the width \\u -
l\\, so a bounded link costs no more than a logit.

## Usage

``` r
DoublyBoundedLink(
  link_name = character(0),
  link_bounds = integer(0),
  link_params = NULL,
  lwr = integer(0),
  upr = integer(0),
  width = integer(0)
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

- lwr, upr:

  The interval endpoints.

- width:

  The interval width, `upr - lwr`, set by the constructor.

## Value

An S7 object of class `DoublyBoundedLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
and carrying its three properties `link_name`, `link_bounds` and
`link_params`. Its `link_bounds` are the endpoints given, and its
`link_params` holds `lwr` and `upr`.

## Details

The interval width \\W = \mathrm{upr} - \mathrm{lwr}\\ is stored as its
own property, set by the constructor, so that each method reads one
property instead of two.

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
coming back. The forward derivatives are the logit's divided by \\W^k\\
and the inverse ones the logit's multiplied by \\W\\, because \\p\\ is
\\\theta\\ rescaled by the width. The inverse set therefore calls the
compiled kernel of the logit, and the forward set is the logit's
expressions written in \\p\\.

## See also

[`bounded_link()`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.md),
the constructor users call.
