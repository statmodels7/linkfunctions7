# S7 Class for the Power Link

The class
[`power_link`](https://statmodels7.github.io/linkfunctions7/reference/power_link.md)
instantiates for a non-zero exponent. At `lambda = 0` the power link is
the log link by continuity, and
[`power_link`](https://statmodels7.github.io/linkfunctions7/reference/power_link.md)
returns a
[`LogLink`](https://statmodels7.github.io/linkfunctions7/reference/LogLink.md)
instead.

## Usage

``` r
PowerLink(
  link_name = character(0),
  link_bounds = integer(0),
  link_params = NULL,
  lambda = integer(0)
)
```

## Arguments

- lambda:

  The exponent of the transformation.

## Value

An S7 object of class `PowerLink`, inheriting from
[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md).

## See also

[`power_link`](https://statmodels7.github.io/linkfunctions7/reference/power_link.md),
the constructor users call.
