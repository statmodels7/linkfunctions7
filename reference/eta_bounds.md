# The Range a Stencil May Evaluate the Inverse Link On

The image of the link's parameter bounds under
[`linkfun`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md),
used to keep a finite-difference grid inside the set the inverse link is
defined on.

## Usage

``` r
eta_bounds(x)
```

## Arguments

- x:

  A
  [`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
  object.

## Value

A numeric vector of length two.

## Details

A link need not map onto the whole real line: the square root reaches
only the positive half, and a stencil straying outside returns `NaN`,
which would make a numerical derivative missing rather than inaccurate.
The bounds are returned sorted, since a decreasing link reverses them,
and are infinite in the directions where they cannot be established.

## See also

[`link_bounds_clamp`](https://statmodels7.github.io/linkfunctions7/reference/link_bounds_clamp.md)
