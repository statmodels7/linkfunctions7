# S7 Class for the Identity Link

Carries the identity transformation \\\eta = \theta\\ on the whole real
line, for a parameter that needs no chart because it is already
unconstrained. Every derivative is a constant: the first is one and the
higher ones are zero, in both directions.

[`bounded_link()`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.md)
returns an object of this class when neither endpoint is given.

## Usage

``` r
IdentityLink(
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

An S7 object of class `IdentityLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
and carrying its three properties `link_name`, `link_bounds` and
`link_params`. Its `link_bounds` are `c(-Inf, Inf)` and it carries no
link parameters.

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
return their argument. The first derivative is one and every higher one
is zero, each built by
[`const_like()`](https://statmodels7.github.io/linkfunctions7/reference/const_like.md)
so that the result is `NA` wherever the argument is `NA`, although the
derivative does not depend on it.

## See also

[`identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.md),
the constructor users call.
