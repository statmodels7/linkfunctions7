# S7 Class for the Softplus Link

The class
[`softplus_link`](https://statmodels7.github.io/linkfunctions7/reference/softplus_link.md)
instantiates. Carries the scale parameter `a` in a dedicated property.

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

- a:

  The scale parameter, strictly positive.

## Value

An S7 object of class `SoftplusLink`, inheriting from
[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md).

## See also

[`softplus_link`](https://statmodels7.github.io/linkfunctions7/reference/softplus_link.md),
the constructor users call.
