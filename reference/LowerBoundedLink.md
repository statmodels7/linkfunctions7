# S7 Class for a Lower Bounded Link

The class
[`bounded_link`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.md)
instantiates when given only a lower endpoint. It is the log link
shifted to start at `lwr`.

## Usage

``` r
LowerBoundedLink(
  link_name = character(0),
  link_bounds = integer(0),
  link_params = NULL,
  lwr = integer(0)
)
```

## Arguments

- lwr:

  The lower endpoint.

## Value

An S7 object of class `LowerBoundedLink`, inheriting from
[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md).

## See also

[`bounded_link`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.md),
the constructor users call.
