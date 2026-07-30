# S7 Class for a Doubly Bounded Link

The class
[`bounded_link`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.md)
instantiates when given both endpoints. It is the logit of \\p =
(\theta - \mathrm{lwr})/W\\, the position of \\\theta\\ within the
interval.

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

- lwr, upr:

  The interval endpoints.

- width:

  The interval width, `upr - lwr`, set by the constructor.

## Value

An S7 object of class `DoublyBoundedLink`, inheriting from
[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md).

## Details

The interval width \\W = \mathrm{upr} - \mathrm{lwr}\\ is stored as its
own property rather than recomputed. Every method needs it, and reading
two S7 properties and subtracting cost about a third of a call to
[`dlinkinv()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkinv.md);
the constructor is the only place it can change.

## See also

[`bounded_link`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.md),
the constructor users call.
