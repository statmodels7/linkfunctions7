# S7 Class for an Upper Bounded Link

The class
[`bounded_link`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.md)
instantiates when given only an upper endpoint. It is the mirror image
of
[`LowerBoundedLink`](https://statmodels7.github.io/linkfunctions7/reference/LowerBoundedLink.md),
the log of the distance below `upr`.

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

- upr:

  The upper endpoint.

## Value

An S7 object of class `UpperBoundedLink`, inheriting from
[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md).

## See also

[`bounded_link`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.md),
the constructor users call.
